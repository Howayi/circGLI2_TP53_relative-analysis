# TP53/MDM2 expression and p53 activity analysis

## Paired RNA-seq module

| Script | Purpose |
| --- | --- |
| `prepare_resources.py` | Verify agreement between the count matrix and featureCounts source; create sample metadata; download GENCODE v48, CollecTRI, DoRothEA, Hallmark p53, and decoupleR resources |
| `fit_expression.R` | Low-expression filtering, TMM/logCPM, PCA, paired DESeq2 analysis, and sensitivity analyses |
| `infer_activity.R` | Exact GENCODE mapping, CollecTRI/DoRothEA A/B with ULM, paired limma comparisons, sign-flip tests, ssGSEA, and fgsea |
| `compare_MDM2.R` | MDM2 expression, paired comparisons, and relationships with TP53 measures |
| `render_summary.R` | Summary plots of TP53 expression and pathway/activity measures |
| `validate_outputs.py` | Independently recalculate ULM and paired statistics; verify resource hashes and original-cohort output assertions |
| `install_local.R` | Historical Windows installer for a downloaded decoupleR ZIP; run from the repository root |

Provide `gene_counts.matrix.tsv`, `gene_counts.featureCounts.txt`, and `gene_counts.featureCounts.txt.summary` in one input directory. The matrix must have `gene_id` as its first column and integer fragment counts in the remaining columns. Sample names follow the `<patient number>N` / `<patient number>T` convention.

```bash
# Run from the repository root, using upstream 04_counts or an equivalent count directory.
export TP53_INPUT_DIR=/path/to/rnaseq_upstream/04_counts
python TP53_analysis/prepare_resources.py
# Confirm that required resources downloaded successfully and install the R dependencies.
Rscript TP53_analysis/fit_expression.R
Rscript TP53_analysis/infer_activity.R
Rscript TP53_analysis/compare_MDM2.R
Rscript TP53_analysis/render_summary.R
python TP53_analysis/validate_outputs.py
```

In Windows PowerShell, set `$env:TP53_INPUT_DIR = 'D:/path/to/04_counts'`. Outputs are written to this module's `results/`, `figures/`, and `resources/` directories and are excluded by `.gitignore`.

The original analysis used 19 samples, including eight matched pairs. Filtering retained genes with `count >=10 in >=4 samples`. DESeq2 used `~ patient + condition` with T versus N. Sensitivity analyses excluded patients `3/34` or `3/34/57`; the random seed was `530930`. These settings, some plot titles, and validation assertions are embedded in the scripts and require review for a new cohort. ULM excludes the TP53-to-TP53 autoregulatory edge. Activity and gene-expression statistics are exported separately.

## TCGA pan-cancer PROGENy module

`pan_cancer/tp53_progeny_activity.py` uses matched TCGA/cBioPortal mutation status and expression profiles to calculate a PROGENy p53/DDR footprint, compare MUT/WT groups, and estimate Cliff's delta, bootstrap intervals, and BH-adjusted P values.

```bash
python TP53_analysis/pan_cancer/tp53_progeny_activity.py \
  --inputs /path/to/frozen_inputs --out /path/to/progeny_results

# Offline execution requires all frozen inputs and the API caches under the same --out directory.
python TP53_analysis/pan_cancer/tp53_progeny_activity.py \
  --inputs /path/to/frozen_inputs --out /path/to/progeny_results --offline
```

Required inputs are `tp53_patient_status.tsv`, `tp53_mutations.tsv`, `study_ids.tsv`, and `progeny_p53_top100_weights.tsv`. Final validation also reads `progeny_1.34.0.tar.gz`. The script verifies fixed SHA-256 hashes for the model package and weight file; a regenerated or reformatted weight table may fail these checks.

These inputs, the frozen API caches, and the companion scripts `tp53_pancancer_analysis.py` and `extract_progeny_weights.py` were not among the supplied files. Complete offline reproduction therefore requires additional resources. The supplied [`MANUSCRIPT_METHODS_EN.md`](../docs/MANUSCRIPT_METHODS_EN.md) describes the original methods and retains paths from the original project layout.
