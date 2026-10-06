# circGLI2_TP53_relative analysis

Analysis code for circGLI2, TP53/p53 signaling, and INK4–CTNNB1 research. The project includes RNA-seq/circRNA processing, exploratory RNA variant screening, downstream expression and enrichment analysis, paired TP53/MDM2 analysis, and TCGA pan-cancer PROGENy analysis.

This repository contains analysis scripts and documentation. Generated figures, sequencing data, sample result tables, dependency libraries, download caches, and result archives are excluded. Plotting scripts are included so that figures can be generated locally.

## Repository structure

| Directory | Contents |
| --- | --- |
| [`upstream/`](upstream/) | `rnaseq_upstream922.sh`: QC, trimming, alignment, featureCounts, and optional circRNA analysis; `CLEAR_array.sh`: CLEAR SGE array jobs; `ink4_ctnnb1_check.sh`: exploratory RNA allele screening for CDKN2A, CDKN2B, and CTNNB1 |
| [`downstream/`](downstream/README.md) | MLY R scripts for count-to-TPM conversion, differential expression, volcano plots, heatmaps, GO/KEGG analysis, chord diagrams, and RNA/proteomics comparisons |
| [`TP53_analysis/`](TP53_analysis/README.md) | Paired RNA-seq analysis of TP53 expression, ULM transcription factor activity, p53 gene sets, and MDM2; the TCGA PROGENy script is in `pan_cancer/` |
| [`INK4_CTNNB1_analysis/`](INK4_CTNNB1_analysis/README.md) | The `corrected_exact_P0` revision: local transcript annotation, P0 recalculation from integer counts, evidence preparation, plotting, and validation |
| [`docs/`](docs/) | Dependencies, reproduction requirements, the supplied English methods, provenance, and validation notes |
| [`environment/`](environment/) | Python dependencies and R package versions recorded for the original TP53 analysis |

## Getting started

1. Clone or download the repository and prepare the tools and R/Python packages described in [`docs/REPRODUCIBILITY.md`](docs/REPRODUCIBILITY.md).
2. Prepare FASTQ files and matching FASTA/GTF/index resources, or existing count matrices and checkpoint outputs. These data are not bundled with the code.
3. Run upstream scripts on Linux/HPC. CLEAR requires SGE and the separate `clear_py27` environment. R/Python modules can run in a suitable local environment.
4. Review the study design before using another dataset. Some TP53/INK4 thresholds, patient lists, and assertions are specific to the original cohort. The MLY files are interactive research scripts whose paths and preceding objects must be prepared before execution.

## Example workflow

```bash
# Run from the repository root; set paths and strandedness for the actual experiment.
RAW=/path/to/01.RawData OUT=/path/to/rnaseq_upstream \
  REF=/path/to/hg38 LIBRARY_TYPE=reverse MODE=both \
  bash upstream/rnaseq_upstream922.sh

# Check upstream inputs and tools after preparing the references and environment.
RAW=/path/to/01.RawData REF=/path/to/hg38 PREFLIGHT_ONLY=1 \
  bash upstream/rnaseq_upstream922.sh

# Generate exploratory RNA variant evidence.
BAM_DIR=/path/to/rnaseq_upstream/03_bam REF=/path/to/hg38 \
  bash upstream/ink4_ctnnb1_check.sh
```

Gene-mode counts are written to `rnaseq_upstream/04_counts/`. This directory can be used directly as `TP53_INPUT_DIR`. The TP53 and INK4 module READMEs provide their execution sequences. Before submitting CLEAR jobs, edit the cluster paths at the top of the script and create the log directory.

## Interpretation and reproduction scope

- TP53 expression, ULM activity, ssGSEA, GSEA, and PROGENy footprint scores measure different aspects of the data and are reported separately.
- RNA ALT evidence and the `candidate_tumor_only` label are exploratory. The screening script calls for subsequent manual review and DNA confirmation.
- Pan-cancer PROGENy analysis requires frozen mutation tables, sample rosters, and model weights. These inputs and the companion scripts mentioned in the supplied methods were not provided. See [`TP53_analysis/README.md`](TP53_analysis/README.md).
- [`MANUSCRIPT_METHODS_EN.md`](docs/MANUSCRIPT_METHODS_EN.md) retains the English section of the supplied historical methods. Its cohort counts, input snapshots, and validation records describe the original analysis; the biological analyses were not rerun during packaging.
- The publication copy retains the statistical calculations and study thresholds. Directory organization, encoding, input paths, and plotting/validation portability adjustments are recorded in [`docs/source_manifest.json`](docs/source_manifest.json).

Static publication checks are summarized in [`docs/VALIDATION.md`](docs/VALIDATION.md). Complete reproduction requires the corresponding inputs and a compatible software environment.
