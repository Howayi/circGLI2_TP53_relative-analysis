options(stringsAsFactors = FALSE)
Sys.setenv(LANGUAGE = 'en')
script_arg <- sub('^--file=', '', commandArgs(FALSE)[grepl('^--file=', commandArgs(FALSE))])
root <- normalizePath(dirname(script_arg), winslash = '/', mustWork = TRUE)
input_dir <- Sys.getenv('TP53_INPUT_DIR', unset = file.path(dirname(root), 'input_data', 'counts'))
.libPaths(c(file.path(root, 'R_library'), .libPaths()))
suppressPackageStartupMessages({
  library(data.table)
  library(decoupleR)
  library(edgeR)
  library(limma)
  library(GSVA)
  library(fgsea)
  library(ggplot2)
  library(pheatmap)
  library(jsonlite)
})
out <- file.path(root, 'results')
fig <- file.path(root, 'figures')
resdir <- file.path(root, 'resources')
set.seed(530930)
b <- readRDS(file.path(out, 'expression_analysis.rds'))
mapping <- as.data.frame(fread(file.path(resdir, 'gencode_v48_gene_mapping.tsv')))
ann <- mapping[match(b$gene_version, mapping$gene_id_version), ]
stopifnot(all(ann$gene_id_version == b$gene_version), !anyNA(ann$symbol))
fwrite(ann, file.path(out, 'input_gene_annotation.tsv'), sep = '\t')
la <- ann[match(rownames(b$logcpm), ann$gene_id), ]
eligible <- which(la$gene_type == 'protein_coding' & la$symbol != '')
# Resolve repeated symbols without looking at condition or differential statistics.
ordered <- eligible[order(rowMeans(b$logcpm)[eligible], decreasing = TRUE)]
selected <- ordered[!duplicated(la$symbol[ordered])]
mat <- b$logcpm[selected, , drop = FALSE]
rownames(mat) <- la$symbol[selected]
mat <- mat[order(rownames(mat)), , drop = FALSE]
symbol_to_gene <- setNames(rownames(b$logcpm)[selected], la$symbol[selected])
centered <- sweep(mat, 1, rowMeans(mat), '-')
fwrite(data.frame(symbol = rownames(mat), mat, check.names = FALSE),
       file.path(out, 'TMM_logCPM_symbols.tsv'), sep = '\t')
fwrite(data.frame(symbol = rownames(centered), centered, check.names = FALSE),
       file.path(out, 'TMM_logCPM_gene_centered_symbols.tsv'), sep = '\t')
fwrite(data.frame(gene_id = rownames(b$logcpm)[selected], symbol = la$symbol[selected],
                  mean_logCPM = rowMeans(b$logcpm)[selected]),
       file.path(out, 'activity_gene_selection.tsv'), sep = '\t')
message('Exact GENCODE v48 mapping verified. Activity background: ', nrow(mat), ' expressed protein-coding symbols')
annotate_de <- function(x) {
  i <- match(x$gene_id, ann$gene_id)
  data.frame(x, symbol = ann$symbol[i], gene_type = ann$gene_type[i])
}
for (label in c('main', 'sens6', 'sens5')) {
  de <- annotate_de(b[[label]]$results)
  fwrite(de, file.path(out, paste0('annotated_DESeq2_', label, '.tsv')), sep = '\t')
}

ct <- fread(file.path(resdir, 'collectri_regulons.csv'))
ct <- as.data.frame(ct[, .(source, target, mor = weight)])
load(file.path(resdir, 'dorothea_hs.rda'))
da <- as.data.frame(dorothea_hs)
names(da)[names(da) == 'tf'] <- 'source'
nets <- list(CollecTRI = ct, DoRothEA_AB = da[da$confidence %in% c('A', 'B'), ])
gmt <- strsplit(readLines(file.path(resdir, 'hallmark_p53.gmt'))[1], '\t')[[1]]
stopifnot(gmt[1] == 'HALLMARK_P53_PATHWAY')
hallmark <- gmt[-c(1, 2)]
stopifnot(length(hallmark) == 200, !anyDuplicated(hallmark))
fwrite(data.frame(symbol = hallmark, expressed = hallmark %in% rownames(mat)),
       file.path(out, 'Hallmark_p53_gene_coverage.tsv'), sep = '\t')

test_pairs <- function(values, ids) {
  md <- droplevels(b$meta[ids, , drop = FALSE])
  pts <- unique(as.character(md$patient))
  d <- vapply(pts, function(pt) values[paste0(pt, 'T')] - values[paste0(pt, 'N')], numeric(1))
  tt <- t.test(d, mu = 0)
  signs <- as.matrix(expand.grid(rep(list(c(-1, 1)), length(d))))
  p_exact <- mean(abs(drop(signs %*% d / length(d))) >= abs(mean(d)) - 1e-12)
  c(n_pairs = length(d), mean_T_minus_N = mean(d), median_T_minus_N = median(d),
    ci95_low = tt$conf.int[1], ci95_high = tt$conf.int[2],
    paired_t_p = tt$p.value, exact_signflip_p = p_exact,
    pairs_increased = sum(d > 0), pairs_decreased = sum(d < 0))
}
all_test <- list()
tp_scores <- data.frame(sample = b$meta$sample, patient = b$meta$patient,
                        condition = b$meta$condition, paired = b$meta$paired,
                        TP53_raw_count = as.numeric(b$counts['ENSG00000141510', ]),
                        TP53_logCPM = as.numeric(mat['TP53', ]))
coverage <- list()
tp_contrast <- list()
classics <- c('TP53', 'CDKN1A', 'MDM2', 'BBC3', 'BAX', 'GADD45A', 'PMAIP1', 'DDB2', 'RRM2B', 'SESN1', 'SESN2', 'BTG2', 'CCNG1', 'TIGAR')

for (nm in names(nets)) {
  net <- unique(nets[[nm]][, c('source', 'target', 'mor')])
  # Do not let TP53's own RNA expression contribute to its inferred TF activity.
  net <- net[!(net$source == 'TP53' & net$target == 'TP53'), ]
  stopifnot(!anyDuplicated(paste(net$source, net$target)))
  known <- net[net$source == 'TP53', ]
  known$expressed <- known$target %in% rownames(mat)
  fwrite(known, file.path(out, paste0(nm, '_TP53_regulon_coverage.tsv')), sep = '\t')
  coverage[[nm]] <- data.frame(network = nm, TP53_total_targets = nrow(known),
                               TP53_expressed_targets = sum(known$expressed),
                               positive = sum(known$expressed & known$mor > 0),
                               negative = sum(known$expressed & known$mor < 0))
  fwrite(net, file.path(resdir, paste0(nm, '_used_network.tsv')), sep = '\t')
  acts <- decoupleR::run_ulm(mat = centered, network = net, .source = 'source',
                            .target = 'target', .mor = 'mor', minsize = 5)
  # Per-sample ULM enrichment p-values are not the patient-level group-test p-values.
  fwrite(acts, file.path(out, paste0(nm, '_all_TF_sample_ULM.tsv')), sep = '\t')
  a <- dcast(as.data.table(acts), source ~ condition, value.var = 'score')
  score <- as.matrix(a[, -1])
  rownames(score) <- a$source
  score <- score[, b$meta$sample, drop = FALSE]
  tp_scores[[paste0(nm, '_TP53_ULM')]] <- as.numeric(score['TP53', ])
  for (label in c('main', 'sens6', 'sens5')) {
    ids <- b[[label]]$ids
    md <- droplevels(b$meta[ids, , drop = FALSE])
    use_score <- score[, ids, drop = FALSE]
    if (label != 'main') {
      # Recompute normalization and gene centering after excluding low-quality pairs.
      sdge <- calcNormFactors(DGEList(counts = b$counts[b$keep, ids, drop = FALSE]), method = 'TMM')
      smat <- cpm(sdge, log = TRUE, prior.count = 2)[symbol_to_gene[rownames(mat)], , drop = FALSE]
      rownames(smat) <- rownames(mat)
      smat <- sweep(smat, 1, rowMeans(smat), '-')
      sa <- decoupleR::run_ulm(mat = smat, network = net, .source = 'source',
                               .target = 'target', .mor = 'mor', minsize = 5)
      sw <- dcast(as.data.table(sa), source ~ condition, value.var = 'score')
      use_score <- as.matrix(sw[, -1])
      rownames(use_score) <- sw$source
      use_score <- use_score[, ids, drop = FALSE]
      fwrite(sa, file.path(out, paste0(nm, '_', label, '_recomputed_sample_ULM.tsv')), sep = '\t')
    }
    design <- model.matrix(~ patient + condition, md)
    fit <- eBayes(lmFit(use_score, design), trend = FALSE)
    tt <- topTable(fit, coef = 'conditionT', number = Inf, sort.by = 'none')
    tt$TF <- rownames(tt)
    fwrite(tt, file.path(out, paste0(nm, '_', label, '_TF_activity_difference.tsv')), sep = '\t')
    tp <- tt[tt$TF == 'TP53', ]
    tests <- test_pairs(setNames(use_score['TP53', ], colnames(use_score)), ids)
    all_test[[paste(nm, label)]] <- data.frame(network = nm, analysis = label,
                                              as.list(tests), limma_p = tp$P.Value,
                                              limma_FDR_all_TFs = tp$adj.P.Val,
                                              tested_TFs = nrow(tt))

    # Infer contrast-level enrichment from DESeq2 Wald statistics; separate from sample tests.
    de <- annotate_de(b[[label]]$results)
    rk <- de$stat[match(symbol_to_gene[rownames(mat)], de$gene_id)]
    names(rk) <- rownames(mat)
    rk <- rk[is.finite(rk)]
    ca <- decoupleR::run_ulm(mat = matrix(rk, ncol = 1, dimnames = list(names(rk), 'T_vs_N')),
                             network = net, .source = 'source', .target = 'target', .mor = 'mor', minsize = 5)
    ca$FDR_all_TFs <- p.adjust(ca$p_value, method = 'BH')
    fwrite(ca, file.path(out, paste0(nm, '_', label, '_contrast_ULM.tsv')), sep = '\t')
    tp_contrast[[paste(nm, label)]] <- data.frame(network = nm, analysis = label,
                                                ca[ca$source == 'TP53', c('score', 'p_value', 'FDR_all_TFs')])
  }
  de <- annotate_de(b$main$results)
  target_de <- merge(known, de, by.x = 'target', by.y = 'symbol', all.x = TRUE)
  target_de$signed_Wald_stat <- target_de$mor * target_de$stat
  target_de <- target_de[order(abs(target_de$signed_Wald_stat), decreasing = TRUE, na.last = TRUE), ]
  fwrite(target_de, file.path(out, paste0(nm, '_TP53_target_DESeq2.tsv')), sep = '\t')
}

# ssGSEA is a pathway signature score; it is not a TP53-specific regulon score.
ssparam <- ssgseaParam(mat, geneSets = list(HALLMARK_P53_PATHWAY = hallmark),
                       minSize = 10, alpha = 0.25, normalize = FALSE, verbose = FALSE)
ss <- gsva(ssparam, verbose = FALSE, BPPARAM = BiocParallel::SerialParam())
tp_scores$Hallmark_p53_ssGSEA <- as.numeric(ss[1, b$meta$sample])
for (label in c('main', 'sens6', 'sens5')) {
  all_test[[paste('Hallmark', label)]] <- data.frame(network = 'Hallmark_p53_ssGSEA', analysis = label,
                                                    as.list(test_pairs(setNames(ss[1, ], colnames(ss)), b[[label]]$ids)),
                                                    limma_p = NA_real_, limma_FDR_all_TFs = NA_real_, tested_TFs = NA_integer_)
}
fwrite(rbindlist(all_test), file.path(out, 'TP53_activity_group_tests.tsv'), sep = '\t')
fwrite(rbindlist(tp_contrast), file.path(out, 'TP53_contrast_regulon_enrichment.tsv'), sep = '\t')
fwrite(rbindlist(coverage), file.path(out, 'TP53_regulon_coverage_summary.tsv'), sep = '\t')
fwrite(tp_scores, file.path(out, 'TP53_per_sample_scores.tsv'), sep = '\t')

fg_results <- list()
for (label in c('main', 'sens6', 'sens5')) {
  de <- annotate_de(b[[label]]$results)
  rank_stats <- de$stat[match(symbol_to_gene[rownames(mat)], de$gene_id)]
  names(rank_stats) <- rownames(mat)
  rank_stats <- sort(rank_stats[is.finite(rank_stats)], decreasing = TRUE)
  fg <- fgseaMultilevel(pathways = list(HALLMARK_P53_PATHWAY = hallmark), stats = rank_stats,
                        minSize = 10, maxSize = 1000, eps = 0, nproc = 1)
  fg$analysis <- label
  # One prespecified pathway is tested; fgsea padj equals the within-one-pathway adjustment.
  fg$leadingEdge <- vapply(fg$leadingEdge, paste, collapse = ';', character(1))
  fg_results[[label]] <- fg
  if (label == 'main') {
    png(file.path(fig, 'Hallmark_p53_GSEA.png'), width = 1050, height = 650, res = 140)
    print(plotEnrichment(hallmark, rank_stats) + ggtitle('HALLMARK_P53_PATHWAY: paired T vs N'))
    dev.off()
  }
}
fwrite(rbindlist(fg_results), file.path(out, 'Hallmark_p53_GSEA.tsv'), sep = '\t')

de <- annotate_de(b$main$results)
fwrite(de[de$symbol %in% classics, ], file.path(out, 'classic_p53_targets_DESeq2.tsv'), sep = '\t')
tp_meta <- tp_scores
tp_meta$Assigned_fraction <- b$meta$assigned_fraction
correlations <- do.call(rbind, lapply(c('CollecTRI_TP53_ULM', 'DoRothEA_AB_TP53_ULM', 'Hallmark_p53_ssGSEA'), function(nm) {
  cr <- cor.test(tp_meta$TP53_logCPM, tp_meta[[nm]], method = 'spearman', exact = FALSE)
  qc <- cor.test(tp_meta$Assigned_fraction, tp_meta[[nm]], method = 'spearman', exact = FALSE)
  data.frame(score = nm, rho_TP53_mRNA = unname(cr$estimate), p_TP53_mRNA = cr$p.value,
             rho_assigned = unname(qc$estimate), p_assigned = qc$p.value)
}))
fwrite(correlations, file.path(out, 'score_correlations_exploratory.tsv'), sep = '\t')

long <- melt(as.data.table(tp_scores), id.vars = c('sample', 'patient', 'condition', 'paired'),
             measure.vars = c('TP53_logCPM', 'CollecTRI_TP53_ULM', 'DoRothEA_AB_TP53_ULM', 'Hallmark_p53_ssGSEA'),
             variable.name = 'measure', value.name = 'score')
long <- long[paired == TRUE]
p <- ggplot(long, aes(condition, score, group = patient)) +
  geom_line(color = '#999999', linewidth = 0.5) + geom_point(aes(color = condition), size = 2.5) +
  facet_wrap(~ measure, scales = 'free_y', ncol = 2) +
  scale_color_manual(values = c(N = '#2878B5', T = '#D9534F')) + theme_bw(base_size = 11) +
  labs(title = 'TP53 expression and inferred p53 activity: 8 matched pairs', x = 'N = adjacent; T = tumor', y = 'Score')
ggsave(file.path(fig, 'TP53_paired_expression_activity.png'), p, width = 10, height = 7, dpi = 180)
ggsave(file.path(fig, 'TP53_paired_expression_activity.pdf'), p, width = 10, height = 7)
heatgenes <- intersect(classics[-1], rownames(mat))
z <- t(scale(t(mat[heatgenes, b$main$ids, drop = FALSE])))
md <- data.frame(condition = b$meta[b$main$ids, 'condition'])
rownames(md) <- b$main$ids
pheatmap(z, cluster_cols = FALSE, annotation_col = md, fontsize = 9,
         main = 'Classic p53 target expression (row Z-score)',
         filename = file.path(fig, 'classic_p53_targets_heatmap.png'), width = 9, height = 5)
writeLines(capture.output(sessionInfo()), file.path(root, 'session_activity.txt'))
versions <- data.frame(package = c('DESeq2', 'edgeR', 'limma', 'decoupleR', 'GSVA', 'fgsea', 'data.table', 'ggplot2', 'pheatmap', 'jsonlite'),
                       version = vapply(c('DESeq2', 'edgeR', 'limma', 'decoupleR', 'GSVA', 'fgsea', 'data.table', 'ggplot2', 'pheatmap', 'jsonlite'),
                                        function(p) as.character(packageVersion(p)), character(1)))
fwrite(versions, file.path(root, 'package_versions.tsv'), sep = '\t')
write_json(list(activity_universe = nrow(mat), gene_symbols = 'exact GENCODE v48 gene annotation',
                duplicate_symbol_rule = 'retain highest mean logCPM over all samples',
                primary = 'CollecTRI + decoupleR ULM on gene-centered TMM logCPM',
                robustness = 'DoRothEA A/B + ULM',
                sample_center = 'gene mean across all 19 samples',
                ULM_score = 'regression slope t statistic; relative to the study cohort',
                activity_test = 'limma ~ patient + condition; BH correction over all TFs per network',
                paired_confirmation = 'paired t-test and all 2^n within-pair sign permutations',
                pathway_score = 'GSVA ssGSEA alpha=0.25 normalize=FALSE',
                pathway_enrichment = 'fgseaMultilevel on DESeq2 Wald statistics; one prespecified Hallmark pathway',
                no_TP53_autoregulation = TRUE), file.path(root, 'activity_config.json'), pretty = TRUE, auto_unbox = TRUE)
message('Activity inference complete')
