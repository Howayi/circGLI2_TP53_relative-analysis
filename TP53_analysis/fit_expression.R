options(stringsAsFactors = FALSE)
Sys.setenv(LANGUAGE = 'en')
suppressPackageStartupMessages({
  library(data.table)
  library(DESeq2)
  library(edgeR)
  library(ggplot2)
  library(jsonlite)
})
script_arg <- sub('^--file=', '', commandArgs(FALSE)[grepl('^--file=', commandArgs(FALSE))])
root <- normalizePath(dirname(script_arg), winslash = '/', mustWork = TRUE)
input_dir <- Sys.getenv('TP53_INPUT_DIR', unset = file.path(dirname(root), 'input_data', 'counts'))
out <- file.path(root, 'results')
fig <- file.path(root, 'figures')
dir.create(out, showWarnings = FALSE)
dir.create(fig, showWarnings = FALSE)
set.seed(530930)
x <- fread(file.path(input_dir, 'gene_counts.matrix.tsv'))
gene_version <- x$gene_id
counts <- as.matrix(x[, -1])
rownames(counts) <- sub('\\.[0-9]+$', '', gene_version)
storage.mode(counts) <- 'integer'
meta <- as.data.frame(fread(file.path(root, 'sample_metadata_qc.tsv')))
rownames(meta) <- meta$sample
stopifnot(identical(colnames(counts), meta$sample), !anyDuplicated(rownames(counts)))
meta$condition <- factor(meta$condition, levels = c('N', 'T'))
meta$patient <- factor(meta$patient)
keep <- rowSums(counts >= 10) >= 4
fwrite(data.frame(gene_id = rownames(counts), gene_id_version = gene_version,
                  kept = keep, samples_count_ge10 = rowSums(counts >= 10)),
       file.path(out, 'expression_filter.tsv'), sep = '\t')
message('Count matrix: ', nrow(counts), ' genes, ', ncol(counts), ' samples; retained: ', sum(keep))
dge <- calcNormFactors(DGEList(counts = counts[keep, ]), method = 'TMM')
logcpm <- cpm(dge, log = TRUE, prior.count = 2)
fwrite(data.frame(gene_id = rownames(logcpm), logcpm, check.names = FALSE),
       file.path(out, 'TMM_logCPM_ensembl.tsv'), sep = '\t')
meta$TMM_factor <- dge$samples$norm.factors
meta$effective_library_size <- dge$samples$lib.size * dge$samples$norm.factors
fwrite(meta, file.path(out, 'sample_metadata_normalized.tsv'), sep = '\t')

# PCA is an exploratory library/composition check, not a sample exclusion rule.
variable <- order(apply(logcpm, 1, var), decreasing = TRUE)[seq_len(min(3000, nrow(logcpm)))]
pc <- prcomp(t(logcpm[variable, ]), center = TRUE, scale. = FALSE)
pd <- data.frame(meta, PC1 = pc$x[, 1], PC2 = pc$x[, 2])
fwrite(pd, file.path(out, 'PCA_coordinates.tsv'), sep = '\t')
p <- ggplot(pd, aes(PC1, PC2, color = condition, shape = assigned_fraction < 0.08)) +
  geom_point(size = 3) + geom_text(aes(label = sample), nudge_y = 1.4, size = 3) +
  scale_color_manual(values = c(N = '#2878B5', T = '#D9534F')) + theme_bw(base_size = 11) +
  labs(title = 'RNA-seq PCA: 19 samples', shape = 'Assigned <8%',
       x = sprintf('PC1 (%.1f%%)', 100 * pc$sdev[1]^2 / sum(pc$sdev^2)),
       y = sprintf('PC2 (%.1f%%)', 100 * pc$sdev[2]^2 / sum(pc$sdev^2)))
ggsave(file.path(fig, 'PCA_all_samples.png'), p, width = 8, height = 5, dpi = 180)
qcplot <- ggplot(meta, aes(factor(sample, levels = sample), 100 * assigned_fraction, fill = condition)) +
  geom_col() + scale_fill_manual(values = c(N = '#2878B5', T = '#D9534F')) +
  theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = 'Sample', y = 'Assigned / featureCounts summary total (%)', title = 'Fragment assignment quality')
ggsave(file.path(fig, 'featurecounts_assignment.png'), qcplot, width = 9, height = 4, dpi = 180)

fit_one <- function(excluded, label) {
  ids <- meta$sample[meta$paired & !as.character(meta$patient) %in% excluded]
  md <- droplevels(meta[ids, , drop = FALSE])
  cm <- counts[keep, ids, drop = FALSE]
  cm <- cm[rowSums(cm) > 0, , drop = FALSE]
  design_matrix <- model.matrix(~ patient + condition, md)
  stopifnot(qr(design_matrix)$rank == ncol(design_matrix))
  message('Fitting ', label, ': ', length(ids) / 2, ' matched pairs; design ~ patient + condition; T vs N')
  dds <- DESeqDataSetFromMatrix(cm, colData = md, design = ~ patient + condition)
  dds <- DESeq(dds, quiet = TRUE, parallel = FALSE)
  # Automatic replacement cannot trigger with one sample per patient-condition cell.
  res <- as.data.frame(results(dds, contrast = c('condition', 'T', 'N'), alpha = 0.05))
  res$gene_id <- rownames(res)
  res <- res[, c('gene_id', 'baseMean', 'log2FoldChange', 'lfcSE', 'stat', 'pvalue', 'padj')]
  fwrite(res, file.path(out, paste0('DESeq2_', label, '_T_vs_N.tsv')), sep = '\t')
  sf <- data.frame(sample = colnames(dds), size_factor = sizeFactors(dds), analysis = label)
  fwrite(sf, file.path(out, paste0('DESeq2_', label, '_size_factors.tsv')), sep = '\t')
  png(file.path(fig, paste0('DESeq2_', label, '_dispersion.png')), width = 1000, height = 700, res = 130)
  plotDispEsts(dds)
  dev.off()
  nc <- counts(dds, normalized = TRUE)
  idx <- match('ENSG00000141510', rownames(nc))
  stopifnot(!is.na(idx))
  tp <- data.frame(sample = ids, patient = md$patient, condition = md$condition,
                   raw_count = cm[idx, ], normalized_count = nc[idx, ])
  fwrite(tp, file.path(out, paste0('TP53_', label, '_counts.tsv')), sep = '\t')
  fwrite(res[res$gene_id == 'ENSG00000141510', ],
         file.path(out, paste0('TP53_', label, '_DESeq2.tsv')), sep = '\t')
  list(results = res, ids = ids, size_factors = sf, TP53 = tp)
}
main <- fit_one(character(), 'paired8')
sens6 <- fit_one(c('3', '34'), 'exclude3_34_paired6')
sens5 <- fit_one(c('3', '34', '57'), 'exclude_lowAssigned_paired5')
bundle <- list(counts = counts, gene_version = gene_version, meta = meta, keep = keep,
               logcpm = logcpm, main = main, sens6 = sens6, sens5 = sens5,
               PCA_variance_explained = pc$sdev^2 / sum(pc$sdev^2))
saveRDS(bundle, file.path(out, 'expression_analysis.rds'))
writeLines(capture.output(sessionInfo()), file.path(root, 'session_expression.txt'))
write_json(list(filter = 'count >=10 in >=4 of 19 samples', normalization = 'edgeR TMM logCPM; prior.count=2',
                main = list(design = '~ patient + condition', contrast = 'T vs N', n_pairs = 8),
                sensitivity6 = list(excluded_patients = c('3', '34'), n_pairs = 6),
                sensitivity5 = list(rule = 'exclude entire pair if either sample assigned fraction <8%',
                                    excluded_patients = c('3', '34', '57'), n_pairs = 5),
                unpaired_samples = c('41T', '44N', '56N'), seed = 530930),
           file.path(root, 'analysis_config.json'), pretty = TRUE, auto_unbox = TRUE)
message('Expression fits complete')
