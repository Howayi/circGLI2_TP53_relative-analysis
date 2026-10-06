# INK4/CTNNB1 exploratory RNA variant evidence

This directory contains the `corrected_exact_P0` revision. `prepare_analysis.py` calculates `P0=(1-T_alt/T_depth)^N_depth` from integer counts, corrects eligible `candidate_tumor_only` / `tumor_alt_normal_underpowered` labels, and retains the original labels and P0 values.

## Required inputs

Provide the following files from the output directory of [`ink4_ctnnb1_check.sh`](../upstream/ink4_ctnnb1_check.sh):

- `allele_counts.long.tsv`, `tumor_normal_pairs.tsv`, `sample_pairing.tsv`, and `sample_names.txt`
- `target_exons.bed`, `depth_per_site.tsv.gz`, and `depth_summary.tsv`
- `candidates.vcf.gz` and `pileup.split.bcf`

Also place the original public Ensembl GRCh38 JSON snapshots in this module's `references/` directory: `CTNNB1.lookup.json`, `CTNNB1.cds.json`, `CDKN2A.lookup.json`, `CDKN2A.cds.json`, `CDKN2A_ARF.cds.json`, `CDKN2B.lookup.json`, and `CDKN2B.cds.json`. These snapshots are not bundled with this code archive. Their original hashes are recorded in [`reference_checksums.json`](reference_checksums.json). Updated public annotations may differ from the original snapshots.

The original reference transcripts were CTNNB1 `ENST00000349496.11`, CDKN2A `ENST00000304494.10`, and CDKN2B `ENST00000276925.7`, with ARF `ENST00000579755.2` additionally examined. Annotation is local to the specified transcripts; it does not implement a full VEP or clinical pathogenicity assessment.

## Execution

```bash
python INK4_CTNNB1_analysis/prepare_analysis.py --input /path/to/07_ink4_ctnnb1
Rscript INK4_CTNNB1_analysis/run_figures.R
python INK4_CTNNB1_analysis/build_report.py
python INK4_CTNNB1_analysis/validate_and_bundle.py --input /path/to/07_ink4_ctnnb1
```

`prepare_analysis.py` checks counts, VCF records, strand sums, coverage, and reference bases, then creates `source_data/` and `validation.json`. `render_figures.R` / `run_figures.R` regenerate six figures from these tables. `build_report.py` creates an interpretation report, and `validate_and_bundle.py` validates local outputs and packages the results. Generated files are excluded by Git.

The original drawing script regenerated only selected revised figures and required the other images to exist. The publication copy generates all six figures from code. Its validator removes byte comparisons against the old image directory while retaining integer-P0, source-hash, count, and export-format checks. Regenerated figures still require manual visual review.

The renderer, report, and validator retain patient lists and assertions for the original cohort: 628 variant records, 727 paired-evidence records, and 34 revised tumor-side candidate records. Review these assertions and report text when applying the code to another dataset. Exploratory RNA ALT evidence requires further review and DNA confirmation before establishing a somatic DNA mutation.
