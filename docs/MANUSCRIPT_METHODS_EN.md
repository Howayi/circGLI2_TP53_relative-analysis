# Manuscript-ready Methods: TP53 mutation status and p53 pathway footprint

> Archive note (not manuscript text): This file retains the existing English section and English appendices of the supplied bilingual methods. Dataset counts, input snapshots, and validation records describe the original analysis. Referenced data and companion scripts are not bundled with this code archive; original project-layout paths are preserved. Replace the repository URL/DOI placeholder before submission. The analysis concerns TCGA PanCancer Atlas mutation-stratified cohorts; GTEx was not included in MUT/WT grouping, pathway scoring, or statistical testing.

### Data sources and study cohorts

This retrospective pan-cancer study used publicly available somatic-mutation and RNA-expression data from 32 TCGA PanCancer Atlas 2018 cohorts, accessed through the cBioPortal REST API on 21 September 2026. The included cohorts were ACC, BLCA, BRCA, CESC, CHOL, COADREAD, DLBC, ESCA, GBM, HNSC, KICH, KIRC, KIRP, LAML, LGG, LIHC, LUAD, LUSC, MESO, OV, PAAD, PCPG, PRAD, SARC, SKCM, STAD, TGCT, THCA, THYM, UCEC, UCS, and UVM. Fully expanded cBioPortal study, mutation- and expression-profile, mutation- and RNA-sample-list identifiers, together with sequenced-sample and unique-patient denominators, are provided in Supplementary Table S1 (`results/tables/supplementary_table_s1_data_ids.tsv`); the source study table and per-profile manifest are provided in `inputs/study_ids.tsv` and `results/tables/data_manifest.tsv`. Only patients with both mutation-sequencing status and available RNA-expression data were retained. The final matched data set contained 9,793 patients, including 3,480 classified as TP53-mutant (MUT) and 6,313 classified as TP53 wild type (WT). Supplementary Data S1 (`results/tables/patient_activity_scores.tsv`; SHA-256 `25622d984c6473c2c2677b826a6c1f9936f4614c50943aac6286700698e69da6`) contains 9,793 unique TCGA patient identifiers and the 9,793 corresponding selected RNA sample identifiers, with one selected sample per patient, together with sample-type code, TP53 status, expression-profile identifier, model-gene coverage, and pathway score.

### TP53 somatic-mutation status and sample selection

TP53 variants were obtained from the `{study_id}_mutations` profile for each cohort, and the analysis denominator was restricted to the corresponding `{study_id}_sequenced` sample list. Retained protein-altering consequences were `Missense_Mutation`, `Nonsense_Mutation`, `Nonstop_Mutation`, `Frame_Shift_Del`, `Frame_Shift_Ins`, `In_Frame_Del`, `In_Frame_Ins`, `Splice_Site`, and `Translation_Start_Site`; synonymous variants were not included. No `Nonstop_Mutation` event was present in the frozen data, although this consequence was part of the prespecified filter. Genomic coordinates and annotations were used as supplied in the GRCh37/hg19 source data. When a patient had more than one RNA sample, one sample was selected according to the prespecified TCGA sample-type priority `01, 03, 09, 02, 04, 05, 06, 07, 08`; ties were resolved by choosing the sample with greater signature-gene coverage and then the lexicographically first sample identifier. A selected RNA sample was classified as TP53-MUT when at least one retained TP53 event was present in that same sample. Otherwise, it was classified as TP53-WT only if it belonged to the mutation-sequenced roster. Ten patients had multiple sequenced samples, and one had a patient-level label discordant with that of the selected RNA sample; all primary analyses used the selected-sample label. WT was therefore an operational definition indicating that no retained event was detected by this workflow and did not exclude TP53-pathway impairment through copy-number change, epigenetic silencing, viral oncogenesis, or upstream regulatory mechanisms. Retained sample-level mutation events are provided as Supplementary Data S2 (`inputs/tp53_mutations.tsv`). Mutation retrieval and curation were performed by the companion script `../TP53_pan_cancer_analysis/tp53_pancancer_analysis.py`; the pathway script reads the frozen mutation tables and patient roster rather than re-downloading mutation data.

For provenance, the canonical PanCancer Atlas MC3 file is `mc3.v0.2.8.PUBLIC.maf.gz` (GDC UUID `1c8cfe5f-e52d-41ba-94da-f15ea1337efc`; MD5 `639ad8f8386e98dacc22e439188aa8fa`; size 753,339,089 bytes). This complete MAF was not downloaded or reprocessed by the workflow. The analyzed mutation calls came from the cBioPortal profiles and the frozen API cache at `../TP53_pan_cancer_analysis/results/raw/cbioportal_tp53_api_cache.json.gz`; the MC3 file identifier is therefore reported for provenance only.

### RNA-expression data and normalization

For each cohort, expression values were retrieved from `{study_id}_rna_seq_v2_mrna_median_all_sample_Zscores` using the `{study_id}_rna_seq_v2_mrna` sample list. This cBioPortal profile contains gene-wise expression z-scores derived from log RNA-Seq V2 RSEM values, with all expression-profiled samples in the same study serving as the reference population. The underlying TCGA PanCancer Atlas resource was `EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.tsv` (GDC UUID `3586c0da-64d0-4b74-a449-5ff4d9136611`; official MD5 `02e72c33071307ff6570621480d3c90b`; size 1,882,540,959 bytes). The underlying matrix was generated from RSEM quantification, upper-quartile normalized, and adjusted with EB++ for sequencing-platform, sequencing-center, and PRAD plate effects. Thus, the workflow contained three distinct scaling stages: (i) RSEM quantification, upper-quartile normalization, and EB++ correction of the underlying expression matrix; (ii) study-wise, gene-wise all-sample z standardization by cBioPortal; and (iii) cancer-wise z standardization of the weighted PROGENy score in the present analysis. No further log transformation or joint cross-cancer standardization was performed.

Expression data were intersected with the mutation-sequenced roster by exact sample identifier. Patients without RNA data were excluded and were neither assigned to the WT group nor assigned an expression value of zero. Cohort-specific matching and signature-coverage metrics are reported in `results/tables/cohort_expression_coverage.tsv`. The 32 cache files in `results/raw/cbioportal_expression/*.json.gz` contain API records retrieved in five-gene chunks, merged and reserialized together with the request parameters and their SHA-256 fingerprints without altering the returned expression values; the matched wide expression matrix is provided in `results/raw/progeny_p53_expression_zscores_wide.tsv.gz`. These files are frozen API-record input snapshots, but the expression values are processed TCGA/cBioPortal measurements rather than raw FASTQ/BAM sequencing data.

### PROGENy-derived p53/DDR transcriptional footprint

The p53 model from `model_human_full` in Bioconductor `progeny` version 1.34.0 (Bioconductor 3.23) was used, retaining the 100 response genes with the smallest model `p.value`. The source package (`inputs/progeny_1.34.0.tar.gz`; SHA-256 `3e8fe926eef24e587e4a1089e4efe0b81a8aa8bd385c91018edeffa9757791d1`) and frozen weight table (`inputs/progeny_p53_top100_weights.tsv`; SHA-256 `c2b9aba4a7793d1ab7b5ebccda4001dd66a6e654b39a722f4bc8ec0b38c945d9`) accompany the analysis. Ninety-nine of the 100 model genes mapped to cBioPortal Entrez Gene identifiers; seven historical symbols were mapped explicitly through stable Entrez identifiers, whereas the withdrawn locus `LOC727916` could not be mapped. Full mapping and coverage information is provided in `results/tables/signature_gene_coverage.tsv`.

For sample \(s\), the raw transcriptional-footprint score was calculated as

\[
R_s=\sum_{g=1}^{G} z_{sg}w_g,
\]

where \(G=99\) is the number of successfully mapped model genes, \(z_{sg}\) is the cBioPortal expression z-score for gene \(g\) in sample \(s\), and \(w_g\) is the signed PROGENy weight. A missing gene-level z-score contributed zero, corresponding to the reference mean for that gene; samples with observations for fewer than 80 of the 99 mapped model genes were excluded. Included samples contained 85–94 observed signature genes (median, 94). The raw score was subsequently standardized within each cancer type as

\[
A_s=\frac{R_s-\overline{R}_{c}}{SD(R_c)},
\]

where \(\overline{R}_{c}\) and \(SD(R_c)\) denote the mean and sample standard deviation of the raw scores among included samples in cancer type \(c\). All figures and inferential analyses used \(A_s\). This quantity is referred to as a “PROGENy-derived p53/DNA-damage-response transcriptional footprint score”; it is not TP53 mRNA abundance, a direct measurement of TP53 protein activity, or a clinical score. The original PROGENy TCGA application used study-wise DESeq2 variance-stabilizing transformation before application of the model weights. The present analysis is a transparent gene-z-score adaptation, and its numerical values are not interchangeable with scores obtained using the original preprocessing. Moreover, because the scores were standardized within cancer type, absolute score levels were not compared across cancer types.

### Statistical analysis and visualization

Within each cancer type, sample size, median, and interquartile range were summarized separately for the TP53-MUT and TP53-WT groups. Inferential testing was performed only when both groups contained at least 10 patients. Distributions were compared using a script-implemented, two-sided asymptotic Mann–Whitney U test with average ranks, tie-adjusted variance, and a 0.5 continuity correction applied to \(|U-n_{MUT}n_{WT}/2|\). Effect size was expressed as Cliff's delta, \(P(A_{MUT}>A_{WT})-P(A_{MUT}<A_{WT})\), such that a negative value indicated a lower distribution in the MUT group. Samples were resampled with replacement separately within the MUT and WT groups for 5,000 stratified bootstrap replicates, and the 95% confidence interval was defined by the 2.5th and 97.5th percentiles of the resulting Cliff's-delta estimates. The two-sided Mann–Whitney P values from the 23 eligible cancer types were adjusted by the Benjamini–Hochberg procedure, and BH-FDR \(q<0.05\) was considered statistically significant. The base random seed was fixed at `20260921`, with a deterministic cancer-code offset used for each cohort. Exact test and bootstrap implementations are archived in `tp53_progeny_activity.py`. Visualizations comprised cancer-specific MUT/WT distribution plots, a forest plot of Cliff's delta with 95% bootstrap confidence intervals, and a heatmap of group medians, generated using Matplotlib 3.10.6.

### Treatment of GTEx data

GTEx data were not included in the primary TP53 mutation-stratified analysis. GTEx comprises non-tumor tissues obtained from donors independent of TCGA and does not provide somatic TP53 mutation calls directly comparable with those generated for TCGA tumors. GTEx samples were therefore not assigned a TP53-WT status and were not used as the WT denominator; all MUT-versus-WT comparisons were performed within individual TCGA cancer types. Directly combining GTEx normal tissues with this analysis would confound tissue source, donor, tumor status, and data-processing workflow. A TCGA–GTEx tumor-versus-normal comparison would constitute a separate analysis requiring a prespecified GTEx release, tissue mapping, complete GTEx sample identifiers, harmonized expression processing, and explicit control of cross-project batch effects; no such analysis is represented in the present results.

### Code, sample identifiers, and data availability

Mutation retrieval and curation are implemented in `../TP53_pan_cancer_analysis/tp53_pancancer_analysis.py`. The pathway-analysis script `tp53_progeny_activity.py` reads the frozen mutation inputs and implements expression API retrieval, request-fingerprint validation, sample matching, PROGENy scoring, statistical analysis, automated validation, and figure generation. `extract_progeny_weights.py` re-extracts and verifies the frozen p53 top-100 weights from the accompanying Bioconductor source package. Python 3.12 is the recommended reproduction environment, and dependency versions are fixed in `requirements.txt`: NumPy 2.3.5, pandas 3.0.1, Matplotlib 3.10.6, and pyreadr 0.5.3. The pathway analysis can be reproduced offline using `python tp53_progeny_activity.py --offline`; expression slices can be retrieved again from cBioPortal using `python tp53_progeny_activity.py --force-download`.

No new wet-laboratory or primary sequencing data were generated in this study. Primary TCGA FASTQ/BAM files are maintained by the NCI Genomic Data Commons and were neither reprocessed nor redistributed here. The publicly available processed data used in this analysis comprised TCGA PanCancer Atlas TP53 mutation profiles, sequenced sample lists, RSEM-derived expression profiles, and the EB++-adjusted PanCancer Atlas expression resource. Barcodes for all 9,793 included samples are provided in Supplementary Data S1. API input snapshots, frozen model weights, derived matrices, result tables, file hashes, and retrieval metadata are provided under `results/raw/`, `inputs/`, and `results/tables/data_manifest.tsv`. Code and redistributable data will be deposited at `[repository name; URL/DOI]` before publication. Because GTEx was not used in the present analysis, no GTEx sample identifiers or GTEx-derived matrices are reported.

## Supplementary Table S1. Exact TCGA PanCancer Atlas identifiers

The machine-readable version of this table is `results/tables/supplementary_table_s1_data_ids.tsv` (SHA-256 `b55b06f19f3ddc6e9597652a52a2411db5867b2fa1da6a296503267e221734e6`). It explicitly contains the fully expanded expression-profile and RNA sample-list IDs, distinguishes sequenced samples from unique sequenced patients, and gives final matched, TP53-MUT, and TP53-WT counts. For compact display below, expression identifiers can also be obtained by appending `_rna_seq_v2_mrna_median_all_sample_Zscores` and `_rna_seq_v2_mrna`, respectively, to each listed `study_id`.

| Cancer | `study_id` | Mutation profile ID | Mutation-sequenced sample-list ID | Sequenced patients |
|---|---|---|---|---:|
| ACC | `acc_tcga_pan_can_atlas_2018` | `acc_tcga_pan_can_atlas_2018_mutations` | `acc_tcga_pan_can_atlas_2018_sequenced` | 91 |
| BLCA | `blca_tcga_pan_can_atlas_2018` | `blca_tcga_pan_can_atlas_2018_mutations` | `blca_tcga_pan_can_atlas_2018_sequenced` | 410 |
| BRCA | `brca_tcga_pan_can_atlas_2018` | `brca_tcga_pan_can_atlas_2018_mutations` | `brca_tcga_pan_can_atlas_2018_sequenced` | 1,066 |
| CESC | `cesc_tcga_pan_can_atlas_2018` | `cesc_tcga_pan_can_atlas_2018_mutations` | `cesc_tcga_pan_can_atlas_2018_sequenced` | 291 |
| CHOL | `chol_tcga_pan_can_atlas_2018` | `chol_tcga_pan_can_atlas_2018_mutations` | `chol_tcga_pan_can_atlas_2018_sequenced` | 36 |
| COADREAD | `coadread_tcga_pan_can_atlas_2018` | `coadread_tcga_pan_can_atlas_2018_mutations` | `coadread_tcga_pan_can_atlas_2018_sequenced` | 534 |
| DLBC | `dlbc_tcga_pan_can_atlas_2018` | `dlbc_tcga_pan_can_atlas_2018_mutations` | `dlbc_tcga_pan_can_atlas_2018_sequenced` | 41 |
| ESCA | `esca_tcga_pan_can_atlas_2018` | `esca_tcga_pan_can_atlas_2018_mutations` | `esca_tcga_pan_can_atlas_2018_sequenced` | 182 |
| GBM | `gbm_tcga_pan_can_atlas_2018` | `gbm_tcga_pan_can_atlas_2018_mutations` | `gbm_tcga_pan_can_atlas_2018_sequenced` | 390 |
| HNSC | `hnsc_tcga_pan_can_atlas_2018` | `hnsc_tcga_pan_can_atlas_2018_mutations` | `hnsc_tcga_pan_can_atlas_2018_sequenced` | 515 |
| KICH | `kich_tcga_pan_can_atlas_2018` | `kich_tcga_pan_can_atlas_2018_mutations` | `kich_tcga_pan_can_atlas_2018_sequenced` | 65 |
| KIRC | `kirc_tcga_pan_can_atlas_2018` | `kirc_tcga_pan_can_atlas_2018_mutations` | `kirc_tcga_pan_can_atlas_2018_sequenced` | 402 |
| KIRP | `kirp_tcga_pan_can_atlas_2018` | `kirp_tcga_pan_can_atlas_2018_mutations` | `kirp_tcga_pan_can_atlas_2018_sequenced` | 276 |
| LAML | `laml_tcga_pan_can_atlas_2018` | `laml_tcga_pan_can_atlas_2018_mutations` | `laml_tcga_pan_can_atlas_2018_sequenced` | 200 |
| LGG | `lgg_tcga_pan_can_atlas_2018` | `lgg_tcga_pan_can_atlas_2018_mutations` | `lgg_tcga_pan_can_atlas_2018_sequenced` | 514 |
| LIHC | `lihc_tcga_pan_can_atlas_2018` | `lihc_tcga_pan_can_atlas_2018_mutations` | `lihc_tcga_pan_can_atlas_2018_sequenced` | 366 |
| LUAD | `luad_tcga_pan_can_atlas_2018` | `luad_tcga_pan_can_atlas_2018_mutations` | `luad_tcga_pan_can_atlas_2018_sequenced` | 566 |
| LUSC | `lusc_tcga_pan_can_atlas_2018` | `lusc_tcga_pan_can_atlas_2018_mutations` | `lusc_tcga_pan_can_atlas_2018_sequenced` | 484 |
| MESO | `meso_tcga_pan_can_atlas_2018` | `meso_tcga_pan_can_atlas_2018_mutations` | `meso_tcga_pan_can_atlas_2018_sequenced` | 86 |
| OV | `ov_tcga_pan_can_atlas_2018` | `ov_tcga_pan_can_atlas_2018_mutations` | `ov_tcga_pan_can_atlas_2018_sequenced` | 523 |
| PAAD | `paad_tcga_pan_can_atlas_2018` | `paad_tcga_pan_can_atlas_2018_mutations` | `paad_tcga_pan_can_atlas_2018_sequenced` | 179 |
| PCPG | `pcpg_tcga_pan_can_atlas_2018` | `pcpg_tcga_pan_can_atlas_2018_mutations` | `pcpg_tcga_pan_can_atlas_2018_sequenced` | 178 |
| PRAD | `prad_tcga_pan_can_atlas_2018` | `prad_tcga_pan_can_atlas_2018_mutations` | `prad_tcga_pan_can_atlas_2018_sequenced` | 494 |
| SARC | `sarc_tcga_pan_can_atlas_2018` | `sarc_tcga_pan_can_atlas_2018_mutations` | `sarc_tcga_pan_can_atlas_2018_sequenced` | 255 |
| SKCM | `skcm_tcga_pan_can_atlas_2018` | `skcm_tcga_pan_can_atlas_2018_mutations` | `skcm_tcga_pan_can_atlas_2018_sequenced` | 438 |
| STAD | `stad_tcga_pan_can_atlas_2018` | `stad_tcga_pan_can_atlas_2018_mutations` | `stad_tcga_pan_can_atlas_2018_sequenced` | 436 |
| TGCT | `tgct_tcga_pan_can_atlas_2018` | `tgct_tcga_pan_can_atlas_2018_mutations` | `tgct_tcga_pan_can_atlas_2018_sequenced` | 149 |
| THCA | `thca_tcga_pan_can_atlas_2018` | `thca_tcga_pan_can_atlas_2018_mutations` | `thca_tcga_pan_can_atlas_2018_sequenced` | 489 |
| THYM | `thym_tcga_pan_can_atlas_2018` | `thym_tcga_pan_can_atlas_2018_mutations` | `thym_tcga_pan_can_atlas_2018_sequenced` | 123 |
| UCEC | `ucec_tcga_pan_can_atlas_2018` | `ucec_tcga_pan_can_atlas_2018_mutations` | `ucec_tcga_pan_can_atlas_2018_sequenced` | 517 |
| UCS | `ucs_tcga_pan_can_atlas_2018` | `ucs_tcga_pan_can_atlas_2018_mutations` | `ucs_tcga_pan_can_atlas_2018_sequenced` | 57 |
| UVM | `uvm_tcga_pan_can_atlas_2018` | `uvm_tcga_pan_can_atlas_2018_mutations` | `uvm_tcga_pan_can_atlas_2018_sequenced` | 80 |

## Supplementary files and reporting level

| Proposed item | File | Contents |
|---|---|---|
| Supplementary Table S1 | `results/tables/supplementary_table_s1_data_ids.tsv` | Fully expanded TCGA study, mutation/expression-profile and mutation/RNA sample-list IDs; sequenced-sample, unique-patient, matched, MUT and WT counts |
| Supplementary Data S1 | `results/tables/patient_activity_scores.tsv` | All 9,793 patient/sample IDs, TP53 status, sample type, coverage, profile ID and activity score |
| Supplementary Data S2 | `inputs/tp53_mutations.tsv` | Retained sample-level TP53 protein-altering events |
| Supplementary Data S3 | `results/tables/cancer_mut_vs_wt_statistics.tsv` | Group sizes, descriptive statistics, U, P, Cliff's delta, bootstrap CI and BH q |
| Supplementary Data S4 | `results/tables/cohort_expression_coverage.tsv`; `results/tables/signature_gene_coverage.tsv` | Cohort matching and model-gene mapping/coverage |
| Source-data manifest | `results/tables/data_manifest.tsv` | Dataset/profile IDs, normalization, URLs, retrieval time, local files and SHA-256 hashes |
| Frozen expression inputs | `results/raw/cbioportal_expression/*.json.gz`; `results/raw/progeny_p53_expression_zscores_wide.tsv.gz` | API response snapshots and matched expression matrix |
| Analysis code | `../TP53_pan_cancer_analysis/tp53_pancancer_analysis.py`; `tp53_progeny_activity.py`; `extract_progeny_weights.py`; `requirements.txt` | Mutation retrieval/curation, pathway analysis, model extraction and fixed dependencies |
| Validation | `results/analysis_metadata.json`; `results/validation_report.txt` | Machine-readable method summary and automated checks |

## Reproduction commands

```powershell
python -m venv .venv
.venv\Scripts\python -m pip install -r requirements.txt

# Fully offline reproduction from frozen inputs
.venv\Scripts\python tp53_progeny_activity.py --offline

# Optional re-download of the cBioPortal expression slices
.venv\Scripts\python tp53_progeny_activity.py --force-download

# Optional re-extraction and hash verification of PROGENy p53 weights
.venv\Scripts\python extract_progeny_weights.py
```

## Suggested references

1. Cerami E, et al. The cBio Cancer Genomics Portal: an open platform for exploring multidimensional cancer genomics data. *Cancer Discovery*. 2012;2:401–404. doi:10.1158/2159-8290.CD-12-0095.
2. Gao J, et al. Integrative analysis of complex cancer genomics and clinical profiles using the cBioPortal. *Science Signaling*. 2013;6:pl1. doi:10.1126/scisignal.2004088.
3. Hoadley KA, et al. Cell-of-origin patterns dominate the molecular classification of 10,000 tumors from 33 types of cancer. *Cell*. 2018;173:291–304.e6. doi:10.1016/j.cell.2018.03.022.
4. Schubert M, et al. Perturbation-response genes reveal signaling footprints in cancer gene expression. *Nature Communications*. 2018;9:20. doi:10.1038/s41467-017-02391-6.
5. GTEx Consortium. The GTEx Consortium atlas of genetic regulatory effects across human tissues. *Science*. 2020;369:1318–1330. doi:10.1126/science.aaz1776.

Authoritative resource pages: [cBioPortal API](https://docs.cbioportal.org/web-api-and-clients/), [cBioPortal RNA expression normalization FAQ](https://docs.cbioportal.org/user-guide/faq/), [GDC PanCancer Atlas](https://gdc.cancer.gov/about-data/publications/pancanatlas), [Bioconductor PROGENy](https://bioconductor.org/packages/release/bioc/html/progeny.html), and [GTEx documentation](https://gtexportal.org/home/documentationPage).
