args <- commandArgs(trailingOnly = FALSE)
script_arg <- sub('^--file=', '', args[grepl('^--file=', args)])
root <- normalizePath(dirname(script_arg), winslash = '/', mustWork = TRUE)
.libPaths(c(file.path(root, 'R_library'), .libPaths()))
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(ggrepel)
  library(jsonlite)
})
options(encoding = 'UTF-8')
out <- file.path(root, 'figures')
dir.create(out, showWarnings = FALSE)
read_data <- function(name) fread(file.path(root, 'source_data', name), na.strings = 'NA', encoding = 'UTF-8')
a <- read_data('allele_evidence_annotated.tsv')
pair <- read_data('paired_evidence_annotated.tsv')
candidates <- read_data('tumor_only_candidates_review.tsv')
coverage <- read_data('coverage_recomputed.tsv')
ann <- read_data('variant_annotations.tsv')
metadata <- read_data('sample_metadata.tsv')
blocks <- read_data('canonical_transcript_blocks.tsv')
stopifnot(nrow(ann) == 628, nrow(pair) == 727, nrow(candidates) == 34)

font <- if ('Microsoft YaHei' %in% systemfonts::system_fonts()$family) 'Microsoft YaHei' else 'sans'
pal <- c(normal = '#447EA7', tumor = '#C87739', neutral = '#909AA3', accent = '#896591')
theme_set(theme_classic(base_size = 8, base_family = font) +
  theme(axis.line = element_line(linewidth = 0.25), axis.ticks = element_line(linewidth = 0.25),
        axis.text = element_text(colour = '#26333D'), plot.title = element_text(face = 'bold', size = 12),
        plot.subtitle = element_text(size = 8.5, margin = margin(b = 8)),
        plot.caption = element_text(size = 6.6, hjust = 0, colour = '#4F5C66', lineheight = 1.25),
        strip.background = element_rect(fill = '#EEF2F5', colour = NA),
        strip.text = element_text(face = 'bold', size = 8),
        legend.position = 'bottom', legend.title = element_text(size = 7.5),
        legend.text = element_text(size = 7), plot.margin = margin(10, 12, 10, 10)))

exports <- list()
save_figure <- function(plot, name, width_mm, height_mm) {
  filename <- file.path(out, name)
  width <- width_mm / 25.4
  height <- height_mm / 25.4
  for (kind in c('png', 'pdf', 'svg', 'tiff')) {
    path <- paste0(filename, '.', kind)
    if (kind == 'png') ragg::agg_png(path, width = width, height = height, units = 'in', res = 300, background = 'white')
    if (kind == 'pdf') grDevices::cairo_pdf(path, width = width, height = height, family = font, bg = 'white')
    if (kind == 'svg') svglite::svglite(path, width = width, height = height, bg = 'white')
    if (kind == 'tiff') ragg::agg_tiff(path, width = width, height = height, units = 'in', res = 600, compression = 'lzw', background = 'white')
    print(plot)
    dev.off()
    stopifnot(file.exists(path), file.info(path)$size > 1000)
  }
  exports[[length(exports) + 1]] <<- data.table(figure = name, width_mm = width_mm, height_mm = height_mm)
}

patients <- c(1, 3, 6, 28, 33, 34, 42, 57)
paired_samples <- unlist(lapply(patients, function(p) paste0(p, c('N', 'T'))))
sample_order <- c(paired_samples, '41T', '44N', '56N')
genes <- c('CDKN2A', 'CDKN2B', 'CTNNB1')

# The denominator is the union of padded target exons, not the complete gene or only its CDS.
cov_long <- melt(coverage, id.vars = c('target', 'sample', 'patient', 'condition', 'pair_status', 'target_bases'),
                 measure.vars = c('frac_depth_ge10', 'frac_depth_ge30'), variable.name = 'threshold', value.name = 'fraction')
cov_long[, target := factor(target, levels = genes)]
cov_long[, sample_axis := factor(ifelse(pair_status == 'paired', sample, paste0(sample, ' *')),
                                levels = rev(ifelse(sample_order %in% c('41T', '44N', '56N'), paste0(sample_order, ' *'), sample_order)))]
cov_long[, threshold := factor(threshold, levels = c('frac_depth_ge10', 'frac_depth_ge30'), labels = c('Depth ≥ 10', 'Depth ≥ 30'))]
cov_long[, label_text := sprintf('%.1f', 100 * fraction)]
p_cov <- ggplot(cov_long, aes(target, sample_axis, fill = fraction)) +
  geom_tile(colour = 'white', linewidth = 0.4) +
  geom_text(aes(label = label_text, colour = fraction > 0.55), size = 2.4) +
  scale_colour_manual(values = c('FALSE' = '#26333D', 'TRUE' = 'white'), guide = 'none') +
  scale_fill_gradient(low = '#F1F5F8', high = '#1B628D', limits = c(0, 1), breaks = seq(0, 1, 0.25), labels = scales::label_percent()) +
  facet_wrap(~threshold, nrow = 1) +
  labs(title = 'Coverage limits the interpretation of negative results', subtitle = '19 samples; cells show target bases reaching the specified depth (%)',
       x = NULL, y = NULL, fill = 'Target bases covered',
       caption = 'N: normal; T: tumor; *: unpaired. Targets: merged exons plus 10-bp flanks, including UTRs and isoforms.\nLow coverage does not prove absence of mutation; RNA depth cannot directly establish DNA deletion.') +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank())
save_figure(p_cov, '01_coverage_heatmap', 183, 190)

status_names <- c(candidate_tumor_only = 'Tumor-only candidate', alt_in_both = 'ALT in both', normal_alt_only = 'Normal ALT only',
                  tumor_alt_normal_low_alt = 'Low normal ALT count', tumor_alt_normal_underpowered = 'Normal underpowered',
                  tumor_alt_normal_no_coverage = 'No normal coverage')
status_cols <- c(candidate_tumor_only = '#C87739', alt_in_both = '#447EA7', normal_alt_only = '#829FA5',
                 tumor_alt_normal_low_alt = '#AFA0BC', tumor_alt_normal_underpowered = '#CCD2D8',
                 tumor_alt_normal_no_coverage = '#626F7A')
counts <- pair[, .N, by = .(target, patient, label)]
counts[, target := factor(target, levels = genes)]
counts[, patient := factor(patient, levels = patients)]
counts[, label := factor(label, levels = names(status_names))]
totals <- pair[, .N, by = .(target, patient)]
totals[, target := factor(target, levels = genes)]
totals[, patient := factor(patient, levels = patients)]
p_status <- ggplot(counts, aes(patient, N, fill = label)) + geom_col(width = 0.72) +
  geom_text(data = totals, aes(patient, N, label = N), inherit.aes = FALSE, vjust = -0.45, size = 2.5) +
  facet_wrap(~target, nrow = 1, scales = 'free_y') +
  scale_fill_manual(values = status_cols, labels = status_names, drop = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.14))) +
  labs(title = 'Candidate categories are not confirmed diagnoses', subtitle = '8 paired patients; 727 patient × site × ALT comparison records',
       x = 'Patient ID', y = 'Comparison records', fill = NULL,
       caption = 'All 34 tumor-only records are in CTNNB1 (31 distinct site/ALT pairs); none are in CDKN2A or CDKN2B.\nExact-count P0 corrects 7 labels: 6 candidates removed, 1 added. Y-axis scales differ by gene.') +
  theme(legend.position = 'bottom', axis.text.x = element_text(size = 7)) +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE))
save_figure(p_status, '02_candidate_status_overview', 183, 125)

hero_positions <- c(41239336, 41239897, 41240109, 41239922, 21968200, 41224549)
hero <- a[pos %in% hero_positions & !(pos == 41239897 & alt != 'G')]
hero <- merge(hero, metadata, by = 'sample', all.x = TRUE)
hero[, patient := as.integer(patient)]
hero <- hero[patient %in% patients]
hero[, condition := factor(condition, levels = c('N', 'T'))]
hero[, patient := factor(patient, levels = patients)]
hero[, patient_x := as.numeric(patient) + ifelse(condition == 'N', -0.12, 0.12)]
hero_labels <- c('41239336' = 'CTNNB1  chr3:41239336 C>T\np.Asp780= (synonymous)',
                 '41239897' = 'CTNNB1  chr3:41239897 T>G\nCanonical transcript 3′UTR',
                 '41240109' = 'CTNNB1  chr3:41240109 C>CTAAT\nCanonical 3′UTR (4-bp insertion)',
                 '41239922' = 'CTNNB1  chr3:41239922 T>A\nCanonical transcript 3′UTR',
                 '21968200' = 'CDKN2A  chr9:21968200 C>G\n3′UTR in both p16 and ARF',
                 '41224549' = 'CTNNB1  chr3:41224549 G>A\np.Ala13Thr (low allele fraction)')
hero[, facet := factor(as.character(pos), levels = as.character(hero_positions), labels = unname(hero_labels[as.character(hero_positions)]))]
hero[, evidence := factor(ifelse(depth < 20, 'Depth < 20', 'Depth ≥ 20'), levels = c('Depth ≥ 20', 'Depth < 20'))]
hero[, count_text := ifelse(depth > 0, paste0(alt_reads, '/', depth), 'No coverage')]
line_data <- hero[, if (all(depth > 0)) .SD else NULL, by = .(facet, patient)]
p_hero <- ggplot(hero, aes(patient_x, vaf_exact, colour = condition)) +
  geom_line(data = line_data, aes(group = patient), colour = '#CCD2D8', linewidth = 0.45, na.rm = TRUE) +
  geom_point(aes(shape = evidence), size = 2.0, stroke = 0.7, na.rm = TRUE) +
  geom_text(aes(label = ifelse(depth == 0, 'NA', '')), y = -0.045, colour = '#85909A', size = 2.1) +
  facet_wrap(~facet, ncol = 2) +
  scale_colour_manual(values = c(N = unname(pal['normal']), T = unname(pal['tumor'])), labels = c(N = 'Normal N', T = 'Tumor T')) +
  scale_shape_manual(values = c('Depth ≥ 20' = 16, 'Depth < 20' = 1)) +
  scale_x_continuous(breaks = seq_along(patients), labels = patients) +
  scale_y_continuous(limits = c(-0.065, 1.06), breaks = c(0, .25, .5, .75, 1), labels = scales::label_percent(accuracy = 1)) +
  labs(title = 'Paired allele fractions require functional context', subtitle = 'One sample per point; grey lines connect paired N/T samples; GRCh38 coordinates',
       x = 'Patient ID', y = 'ALT / original total AD', colour = NULL, shape = NULL,
       caption = 'p.Asp780= preserves the canonical amino acid; canonical 3′UTR signals do not directly imply protein changes.\nSome strong signals also occur in other normals. Open: depth < 20; NA: no coverage. RNA VAF is not a DNA genotype.') +
  theme(strip.text = element_text(size = 7.7), panel.spacing = grid::unit(6, 'mm'))
fwrite(hero, file.path(root, 'source_data', 'hero_figure_data.tsv'), sep = '\t', na = 'NA')
save_figure(p_hero, '03_paired_VAF_evidence', 183, 215)

# The heatmap contains every candidate allele after exact-count P0 relabelling.
candidate_ids <- unique(candidates$variant_id)
row_info <- ann[variant_id %in% candidate_ids]
row_info[, short_context := fifelse(canonical_context == 'CDS', protein_change,
                                    fifelse(canonical_context == '3_prime_UTR', '3′UTR',
                                            fifelse(canonical_context == '5_prime_UTR', '5′UTR', 'Outside canonical exons')))]
row_info <- row_info[order(pos, alt)]
row_info[, row_label := paste0(format(pos, scientific = FALSE, trim = TRUE), ' ', ref, '>', alt, '  ', short_context)]
hm <- merge(a[variant_id %in% candidate_ids & sample %in% paired_samples], row_info[, .(variant_id, row_label)], by = 'variant_id')
hm[, row_label := factor(row_label, levels = rev(row_info$row_label))]
hm[, sample := factor(sample, levels = paired_samples)]
mark <- merge(candidates[, .(variant_id, sample = paste0(patient, 'T'), one_direction_ALT, p0_exact_above_threshold)], row_info[, .(variant_id, row_label)], by = 'variant_id')
mark[, row_label := factor(row_label, levels = rev(row_info$row_label))]
mark[, sample := factor(sample, levels = paired_samples)]
p_hm <- ggplot(hm, aes(sample, row_label, fill = vaf_exact)) +
  geom_tile(colour = 'white', linewidth = 0.18) +
  geom_tile(data = mark, fill = NA, colour = '#C87739', linewidth = 0.5) +
  geom_point(data = hm[depth > 0 & depth < 20], shape = 1, size = 1.1, stroke = 0.35, colour = '#667581') +
  geom_point(data = mark[one_direction_ALT == TRUE], aes(sample, row_label), inherit.aes = FALSE, shape = 4, size = 1.35, stroke = 0.4, colour = '#A35225') +
  scale_fill_gradient(low = '#F1F6FA', high = '#175D88', na.value = '#C5CDD3', trans = 'sqrt', limits = c(0, 1),
                      breaks = c(0, .01, .05, .1, .5, 1), labels = scales::label_percent()) +
  labs(title = 'Evidence matrix for all CTNNB1 tumor-only candidates', subtitle = '31 site/ALT pairs × 16 paired samples; orange borders mark corrected candidates',
       x = NULL, y = 'chr3 position, REF>ALT and canonical consequence', fill = 'ALT fraction\nSquare-root scale',
       caption = '×: one ALT orientation; ○: depth <20; grey: no coverage. Candidate labels use exact-count P0 ≤0.05.\nThese are evidence flags, not validated filters; RNA-library strand bias alone does not establish false positives.') +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text.y = element_text(size = 5.9),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7), legend.position = 'bottom') +
  guides(fill = guide_colourbar(barwidth = grid::unit(95, 'mm'), barheight = grid::unit(4, 'mm'), title.position = 'top', title.hjust = 0))
fwrite(hm, file.path(root, 'source_data', 'candidate_heatmap_data.tsv'), sep = '\t', na = 'NA')
save_figure(p_hm, '04_all_tumor_only_candidate_heatmap', 183, 245)

# Missense candidates are deliberately shown as weak read evidence, not a mutation burden.
miss <- candidates[canonical_consequence == 'missense_variant']
miss[, evidence_label := paste0(patient, 'T  ', protein_change)]
miss <- miss[order(as.integer(patient), pos)]
miss[, evidence_label := factor(evidence_label, levels = rev(evidence_label))]
miss[, qc_group := factor(ifelse(one_direction_ALT, 'One ALT orientation', 'Both ALT orientations'),
                         levels = c('Both ALT orientations', 'One ALT orientation'))]
miss[, detail := paste0(T_alt, '/', T_depth, '; F/R ', T_alt_fwd, '/', T_alt_rev)]
p_miss <- ggplot(miss, aes(T_vaf_exact, evidence_label)) +
  geom_segment(aes(x = 0, xend = T_vaf_exact, yend = evidence_label), colour = '#CCD2D8', linewidth = 0.6) +
  geom_point(aes(colour = qc_group), size = 2.8) +
  geom_text(aes(x = .023, label = detail), hjust = 0, size = 2.35, colour = '#43515D') +
  scale_colour_manual(values = c('Both ALT orientations' = '#447EA7', 'One ALT orientation' = '#C87739'), drop = FALSE) +
  scale_x_continuous(limits = c(0, .042), breaks = c(0, .005, .01, .015, .02), labels = scales::label_percent(accuracy = .1)) +
  labs(title = 'Missense candidates still have limited read support', subtitle = '9 canonical missense records after exact-count P0 correction; ALT support: 3–8 each',
       x = 'Tumor ALT fraction', y = NULL, colour = NULL,
       caption = 'Right: ALT/total AD and forward/reverse ALT counts. None meets the high-support review criteria.\nPredicted amino acid changes alone do not prove pathogenicity, driver activity or somatic DNA origin.') +
  theme(legend.position = 'bottom')
fwrite(miss, file.path(root, 'source_data', 'missense_figure_data.tsv'), sep = '\t', na = 'NA')
save_figure(p_miss, '05_missense_candidate_evidence', 183, 140)

ink <- merge(a[target %in% c('CDKN2A', 'CDKN2B') & grepl('T$', sample) & alt_reads >= 3],
             ann[, .(variant_id, ARF_context, ARF_protein_change)], by = 'variant_id')
ink <- ink[canonical_consequence == 'missense_variant' | (!is.na(ARF_protein_change) & ARF_protein_change != '')]
ink[, patient := as.integer(sub('T$', '', sample))]
normal_counts <- copy(a[grepl('N$', sample)])
normal_counts[, patient := as.integer(sub('N$', '', sample))]
ink <- merge(ink, normal_counts[, .(variant_id, patient, normal_depth = depth, normal_alt = alt_reads)],
             by = c('variant_id', 'patient'), all.x = TRUE)
ink[, change_label := ifelse(is.na(protein_change) | protein_change == '', paste0('ARF ', ARF_protein_change),
                            ifelse(is.na(ARF_protein_change) | ARF_protein_change == '', protein_change,
                                   paste0(protein_change, ' / ARF ', ARF_protein_change)))]
ink <- ink[order(target, patient, pos)]
ink[, row_label := factor(paste(sample, change_label), levels = rev(paste(sample, change_label)))]
ink[, detail := paste0(alt_reads, '/', depth, '; N ', ifelse(is.na(normal_depth), 'unpaired', paste0(normal_alt, '/', normal_depth)))]
p_ink <- ggplot(ink, aes(vaf_exact, row_label)) +
  geom_segment(aes(x = 0, xend = vaf_exact, yend = row_label), colour = '#CCD2D8', linewidth = 0.6) +
  geom_point(colour = pal['tumor'], size = 2.7) +
  geom_text(aes(x = .16, label = detail), hjust = 0, size = 2.25, colour = '#43515D') +
  facet_grid(target ~ ., scales = 'free_y', space = 'free_y') +
  scale_x_continuous(limits = c(0, .28), breaks = c(0, .05, .1, .15), labels = scales::label_percent(accuracy = 1)) +
  labs(title = 'INK4 protein signals: limited coverage and pairing', subtitle = 'p16, ARF or p15 protein-altering candidates with ALT ≥3 in tumor samples',
       x = 'Tumor ALT fraction', y = NULL,
       caption = 'All 10: ALT 3–4, one orientation. 41T is unpaired; other normals have low depth.\nARF uses a separate transcript; p.Cys123Ter alone does not prove loss. Right: T and N ALT/total AD.') +
  theme(strip.text.y = element_text(angle = 0), axis.text.y = element_text(size = 6.5))
fwrite(ink, file.path(root, 'source_data', 'ink4_protein_candidate_evidence.tsv'), sep = '\t', na = 'NA')
save_figure(p_ink, '06_ink4_protein_candidate_evidence', 183, 150)

status_build <- ggplot_build(p_status)$data[[1]]
heatmap_build <- ggplot_build(p_hm)$data
missense_build <- ggplot_build(p_miss)$data
comparison_records <- sum(status_build$ymax - status_build$ymin)
orange_records <- sum((status_build$ymax - status_build$ymin)[status_build$fill == status_cols[['candidate_tumor_only']]])
stopifnot(comparison_records == 727, orange_records == 34,
          nrow(heatmap_build[[1]]) == 496, nrow(heatmap_build[[2]]) == 34,
          nrow(missense_build[[2]]) == 9)
write_json(list(comparison_records = comparison_records, tumor_only_records = orange_records,
                heatmap_cells = nrow(heatmap_build[[1]]), heatmap_candidate_borders = nrow(heatmap_build[[2]]),
                missense_points = nrow(missense_build[[2]])),
           file.path(root, 'plot_data_QA.json'), auto_unbox = TRUE, pretty = TRUE)
fwrite(rbindlist(exports), file.path(root, 'figure_exports.tsv'), sep = '\t')
capture.output(sessionInfo(), file = file.path(root, 'session_R.txt'))
cat('Updated figures 02, 04, and 05; retained figures 01, 03, and 06 in all four formats.\n')
