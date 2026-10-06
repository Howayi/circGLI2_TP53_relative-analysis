# Environment, inputs, and reproduction scope

## Upstream and checkpoint scripts

Linux Bash >=4 is required. `rnaseq_upstream922.sh` uses FastQC, fastp, HISAT2, samtools, MultiQC, featureCounts, and standard GNU tools. Its circRNA/both modes additionally require BWA, StringTie, CIRIquant, Java, Perl, and the resources checked by the script. Reference paths (`REF/FASTA/GTF/HISAT_INDEX`), strandedness, and conda settings can be supplied through environment variables; see the script header.

`CLEAR_array.sh` requires SGE `qsub`, `NSLOTS/SGE_TASK_ID`, and `clear_quant` in the original `clear_py27` environment. Edit its `BASE_DIR`, cleaned-input paths, references/indexes, and SGE log directives before submission. The log directory must already exist.

```bash
# After editing the CLEAR paths, run from the actual BASE_DIR.
mkdir -p logs
qsub -t 1-$(wc -l < samples.txt) /path/to/repository/upstream/CLEAR_array.sh
```

CLEAR requires `samples.txt` with one sample per line. The main upstream script creates its own sample manifest; prepare CLEAR's sample list separately. Cleaned paired inputs follow `<sample>_1.clean.fq.gz` / `<sample>_2.clean.fq.gz` naming.

`upstream/ink4_ctnnb1_check.sh` requires bcftools, indexed `*.sorted.bam` files, matching hg38 FASTA/FAI files, and a GENCODE GTF. Defaults are CDKN2A/CDKN2B/CTNNB1, at least three ALT reads, MAPQ/BASEQ 20, and MAX_P0 0.05. Record actual software versions and confirm strandedness and reference consistency when running the workflow.

## Python

The supplied pan-cancer methods record Python 3.12 as the reproduction environment.

```bash
python -m venv .venv
# Activate the virtual environment before installing dependencies.
python -m pip install -r environment/requirements.txt
```

The core dependency list contains NumPy, pandas, requests, SciPy, Matplotlib, and pypdf as directly used by the scripts. It is not a complete environment lockfile. NumPy/pandas/Matplotlib versions recorded in the supplied pan-cancer methods are listed separately in `environment/requirements_pan_cancer.txt`; this combination was not installed or revalidated during packaging. pyreadr is needed by the unavailable companion model-extraction script and is not included among the current scripts' core dependencies.

## R

The original TP53 session recorded R 4.6.1. Original package versions are listed in [`../environment/TP53_R_package_versions.tsv`](../environment/TP53_R_package_versions.tsv). This TSV is a historical version record, not an installation lockfile.

- TP53: data.table, DESeq2, edgeR, limma, decoupleR, GSVA, fgsea, ggplot2, pheatmap, jsonlite, and ggrepel.
- INK4 plotting: data.table, ggplot2, patchwork, ggrepel, jsonlite, ragg, svglite, systemfonts, and scales.
- MLY also uses GenomicFeatures, GenomicRanges, DGEobj.utils, readr, dplyr, tidyr, tidyverse, statmod, tinyarray, EnhancedVolcano, clusterProfiler, enrichplot, GOplot, DOSE, ggnewscale, topGO, circlize, ComplexHeatmap, org.Hs.eg.db, org.Mm.eg.db, ggridges, stringr, biomaRt, GO.db, homologene, and VennDiagram. Install packages according to the sections being run.

Install CRAN and Bioconductor packages for the chosen R/Bioconductor release. The TP53 code uses the original GSVA and decoupleR APIs; check compatibility when upgrading. Installing newer versions alone does not reproduce the original results.

## External inputs and supplied methods

FASTQ/BAM/VCF files, sample tables, expression matrices, result tables, figures, and public reference/model snapshots are excluded. The paired TP53 module can download its public resources; download errors are recorded in `resource_manifest.json` and must be resolved before continuing.

The pan-cancer script reads four main TSV inputs and validates `progeny_1.34.0.tar.gz` at the end of the workflow. The companion scripts `tp53_pancancer_analysis.py` and `extract_progeny_weights.py`, frozen API caches, and original `requirements.txt` were not found among the supplied files. The supplied pan-cancer script is retained in full, while complete reproduction requires those additional resources.

Prepare the INK4 reference JSON snapshots separately and verify their hashes. Generated reports, session logs, analysis outputs, dependency libraries, reference sequences, caches, and result archives are excluded. This is a traceable code archive; complete reproduction requires the corresponding inputs and compatible environments.
