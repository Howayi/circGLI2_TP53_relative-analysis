script_arg <- sub('^--file=', '', commandArgs(FALSE)[grepl('^--file=', commandArgs(FALSE))])
root <- normalizePath(dirname(script_arg), winslash = '/', mustWork = TRUE)
input_dir <- Sys.getenv('TP53_INPUT_DIR', unset = file.path(dirname(root), 'input_data', 'counts'))
suppressPackageStartupMessages({library(data.table); library(ggplot2)})
d <- fread(file.path(root, 'results/TP53_per_sample_scores.tsv'))
cts <- fread(file.path(input_dir, 'gene_counts.matrix.tsv'))
d[, TP53_raw_count := as.numeric(cts[grepl('^ENSG00000141510\\.', gene_id), -1])]
fwrite(d, file.path(root, 'results/TP53_per_sample_scores.tsv'), sep = '\t')
l <- melt(d, id.vars = c('sample', 'patient', 'condition', 'paired'),
          measure.vars = c('TP53_logCPM', 'CollecTRI_TP53_ULM', 'DoRothEA_AB_TP53_ULM', 'Hallmark_p53_ssGSEA'),
          variable.name = 'measure', value.name = 'score')
l <- l[paired == TRUE]
l[, condition := factor(condition, levels = c('N', 'T'))]
l[, measure := factor(measure, levels = c('TP53_logCPM', 'CollecTRI_TP53_ULM', 'DoRothEA_AB_TP53_ULM', 'Hallmark_p53_ssGSEA'),
                      labels = c('TP53 mRNA (logCPM)', 'CollecTRI + ULM', 'DoRothEA A/B + ULM', 'p53 pathway (ssGSEA)'))]
p <- ggplot(l, aes(condition, score, group = patient)) +
  geom_line(color = '#999999', linewidth = 0.5) + geom_point(aes(color = condition), size = 2.5) +
  ggrepel::geom_text_repel(data = l[condition == 'N'], aes(label = patient),
                          nudge_x = -0.16, direction = 'y', hjust = 1,
                          color = '#205E8C', size = 3, seed = 530930,
                          box.padding = 0.28, point.padding = 0.2,
                          min.segment.length = 0, segment.color = '#AAAAAA', max.overlaps = Inf) +
  ggrepel::geom_text_repel(data = l[condition == 'T'], aes(label = patient),
                          nudge_x = 0.16, direction = 'y', hjust = 0,
                          color = '#A63835', size = 3, seed = 530930,
                          box.padding = 0.28, point.padding = 0.2,
                          min.segment.length = 0, segment.color = '#AAAAAA', max.overlaps = Inf) +
  scale_x_discrete(expand = expansion(add = 0.5)) +
  scale_y_continuous(expand = expansion(mult = 0.14)) +
  facet_wrap(~ measure, scales = 'free_y', ncol = 2) +
  scale_color_manual(values = c(N = '#2878B5', T = '#D9534F')) + theme_bw(base_size = 11) +
  labs(title = 'TP53 expression and inferred p53 activity: 8 matched pairs',
       x = 'N = adjacent tissue; T = tumor', y = 'Score',
       caption = 'Numbers identify patients; each line connects a matched pair.\nHigher TF scores indicate higher relative inferred activity.')
ggsave(file.path(root, 'figures/TP53_paired_expression_activity.png'), p, width = 11, height = 8, dpi = 180)
ggsave(file.path(root, 'figures/TP53_paired_expression_activity.pdf'), p, width = 11, height = 8)
