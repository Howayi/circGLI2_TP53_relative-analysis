options(stringsAsFactors = FALSE)
script_arg <- sub('^--file=', '', commandArgs(FALSE)[grepl('^--file=', commandArgs(FALSE))])
root <- normalizePath(dirname(script_arg), winslash = '/', mustWork = TRUE)
input_dir <- Sys.getenv('TP53_INPUT_DIR', unset = file.path(dirname(root), 'input_data', 'counts'))
suppressPackageStartupMessages({library(data.table); library(ggplot2); library(jsonlite)})
out <- file.path(root, 'results')
fig <- file.path(root, 'figures')
b <- readRDS(file.path(out, 'expression_analysis.rds'))
annotation <- fread(file.path(out, 'input_gene_annotation.tsv'))
gene <- annotation[symbol == 'MDM2']
stopifnot(nrow(gene) == 1L)
gene_id <- gene$gene_id
original <- fread(file.path(input_dir, 'gene_counts.matrix.tsv'))
stopifnot(identical(as.character(original$gene_id), b$gene_version),
          identical(colnames(original)[-1], colnames(b$counts)),
          all(as.matrix(original[, -1]) == b$counts))

samples <- data.frame(sample = b$meta$sample, patient = b$meta$patient,
                      condition = b$meta$condition, paired = b$meta$paired,
                      raw_count = as.numeric(b$counts[gene_id, ]),
                      TMM_logCPM = as.numeric(b$logcpm[gene_id, ]),
                      assigned_fraction = b$meta$assigned_fraction)
samples$DESeq2_normalized_count_paired8 <- NA_real_
summaries <- list()
for (label in c('main', 'sens6', 'sens5')) {
  fit <- b[[label]]
  de <- fit$results[fit$results$gene_id == gene_id, ]
  stopifnot(nrow(de) == 1L, is.finite(de$log2FoldChange), is.finite(de$lfcSE),
            identical(fit$ids, fit$size_factors$sample))
  nc <- b$counts[gene_id, fit$ids] / fit$size_factors$size_factor
  row <- data.frame(analysis = label, n_pairs = length(fit$ids) / 2,
                     gene_id = gene_id, symbol = 'MDM2', de[, -1],
                     fold_T_over_N = 2^de$log2FoldChange,
                     log2FC_CI95_low = de$log2FoldChange - qnorm(0.975) * de$lfcSE,
                     log2FC_CI95_high = de$log2FoldChange + qnorm(0.975) * de$lfcSE)
  row$fold_CI95_low <- 2^row$log2FC_CI95_low
  row$fold_CI95_high <- 2^row$log2FC_CI95_high
  summaries[[label]] <- row
  sd <- samples[match(fit$ids, samples$sample), ]
  sd$DESeq2_normalized_count <- as.numeric(nc)
  fwrite(sd, file.path(out, paste0('MDM2_', label, '_sample_expression.tsv')), sep = '\t')
  if (label == 'main') samples$DESeq2_normalized_count_paired8[match(fit$ids, samples$sample)] <- nc
}
summary <- as.data.frame(rbindlist(summaries))
fwrite(summary, file.path(out, 'MDM2_DESeq2_comparisons.tsv'), sep = '\t')
fwrite(samples, file.path(out, 'MDM2_all19_sample_expression.tsv'), sep = '\t')
paired <- samples[samples$paired, ]
patients <- unique(as.character(paired$patient))
pairs <- do.call(rbind, lapply(patients, function(pat) {
  n <- paired[paired$sample == paste0(pat, 'N'), ]
  t <- paired[paired$sample == paste0(pat, 'T'), ]
  data.frame(patient = pat, N_normalized_count = n$DESeq2_normalized_count_paired8,
              T_normalized_count = t$DESeq2_normalized_count_paired8,
              descriptive_fold_T_over_N = t$DESeq2_normalized_count_paired8 / n$DESeq2_normalized_count_paired8,
              N_raw_count = n$raw_count, T_raw_count = t$raw_count,
              N_assigned_fraction = n$assigned_fraction, T_assigned_fraction = t$assigned_fraction)
}))
pairs$descriptive_log2FC <- log2(pairs$descriptive_fold_T_over_N)
fwrite(pairs, file.path(out, 'MDM2_patient_pair_ratios.tsv'), sep = '\t')

colors <- setNames(c('#2878B5', '#D9534F'), c('N', 'T'))
p <- ggplot(paired, aes(condition, DESeq2_normalized_count_paired8, group = patient)) +
  geom_line(color = '#999999', linewidth = 0.65) +
  geom_point(aes(color = condition), size = 3) + scale_color_manual(values = colors) +
  ggrepel::geom_text_repel(data = paired[paired$condition == 'N', ], aes(label = patient),
                          nudge_x = -0.14, direction = 'y', hjust = 1,
                          color = '#205E8C', size = 3.6, seed = 530930,
                          box.padding = 0.3, point.padding = 0.2,
                          min.segment.length = 0, segment.color = '#AAAAAA', max.overlaps = Inf) +
  ggrepel::geom_text_repel(data = paired[paired$condition == 'T', ], aes(label = patient),
                          nudge_x = 0.14, direction = 'y', hjust = 0,
                          color = '#A63835', size = 3.6, seed = 530930,
                          box.padding = 0.3, point.padding = 0.2,
                          min.segment.length = 0, segment.color = '#AAAAAA', max.overlaps = Inf) +
  scale_x_discrete(expand = expansion(add = 0.45)) +
  scale_y_continuous(expand = expansion(mult = 0.13)) +
  theme_bw(base_size = 12) +
  labs(title = 'MDM2 mRNA: 8 matched tumor-adjacent pairs',
       subtitle = sprintf('DESeq2 T/N = %.2f (95%% CI %.2f-%.2f); P = %.3f; FDR = %.3f',
                           summary$fold_T_over_N[1], summary$fold_CI95_low[1], summary$fold_CI95_high[1],
                           summary$pvalue[1], summary$padj[1]),
       x = 'N = adjacent tissue; T = tumor', y = 'DESeq2 normalized fragment count',
       caption = 'Numbers identify patients; each line connects a matched pair.\nRNA abundance does not establish protein abundance or DNA amplification.')
ggsave(file.path(fig, 'MDM2_paired_expression.png'), p, width = 8, height = 5, dpi = 180)
ggsave(file.path(fig, 'MDM2_paired_expression.pdf'), p, width = 8, height = 5)
q <- ggplot(pairs, aes(factor(patient, levels = patients), descriptive_log2FC)) +
  geom_hline(yintercept = 0, color = '#888888') +
  geom_col(aes(fill = descriptive_log2FC > 0), width = 0.65, show.legend = FALSE) +
  geom_text(aes(label = sprintf('%.2fx', descriptive_fold_T_over_N)),
            vjust = ifelse(pairs$descriptive_log2FC > 0, -0.3, 1.2), size = 3.6) +
  scale_fill_manual(values = c('FALSE' = '#2878B5', 'TRUE' = '#D9534F')) + theme_bw(base_size = 12) +
  scale_y_continuous(expand = expansion(mult = 0.18)) +
  labs(title = 'MDM2 expression change in each patient', x = 'Patient', y = 'Descriptive log2(T/N)',
       caption = 'Ratios use paired8 DESeq2 normalized counts; each patient has one N and one T sample, so these are descriptive changes.')
ggsave(file.path(fig, 'MDM2_patient_fold_changes.png'), q, width = 8, height = 4.5, dpi = 180)
write_json(list(gene_id_version = gene$gene_id_version, cached_models_reused = TRUE,
                unchanged_raw_matrix_verified = TRUE, design = '~ patient + condition', contrast = 'T vs N',
                primary_pairs = patients, main_fold_change = summary$fold_T_over_N[1],
                main_pvalue = summary$pvalue[1], main_padj = summary$padj[1],
                patients_increased = sum(pairs$descriptive_fold_T_over_N > 1),
                patients_decreased = sum(pairs$descriptive_fold_T_over_N < 1),
                interpretation = 'No significant group-level MDM2 mRNA overexpression; protein and DNA copy number were not measured'),
           file.path(root, 'MDM2_validation.json'), pretty = TRUE, auto_unbox = TRUE)
writeLines(capture.output(sessionInfo()), file.path(root, 'session_MDM2.txt'))
print(summary)
print(pairs)
