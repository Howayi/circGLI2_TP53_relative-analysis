# Downstream MLY R research scripts

The six supplied files retain their original filenames and use UTF-8/LF. The two TPM scripts contain different implementations and are both preserved for provenance.

| File | Contents |
| --- | --- |
| `MLY_count_to_TPM.R` | GTF-derived exon lengths, gene-ID handling, TPM calculation, and merging with differential-expression results |
| `MLY count转化为TPM.R` | Earlier TPM workflow using `DGEobj.utils::convertCounts`; the original filename is retained |
| `MLY_DEG_ANALYS.R` | Differential-expression research fragments, including edgeR, heatmaps, and volcano plots |
| `MLY_ggplot.R` | Labeled volcano plots and differential-gene heatmaps |
| `MLY_GO_KEGG.R` | RNA/proteomics GO, KEGG, GSEA, and visualization fragments |
| `MLY_chord.R` | Shared RNA/proteomics GO terms, mouse-to-human mapping, chord diagrams, and Venn diagrams |

These are interactive research records containing multiple workflow fragments, Windows paths, and references to objects created earlier in a session. Select and execute the relevant sections in R/RStudio after preparing their inputs and objects. The files do not form a directly executable sequential pipeline.

## Items to review before use

- Replace `D:/Bioinfo/...` and `E:/...` with actual input/output paths and prepare preceding objects. The `getwd(...)` call in `MLY_DEG_ANALYS.R` does not set the working directory.
- `MLY_DEG_ANALYS.R` includes a historical two-column normal/tumor example with fixed `bcv=0.4`. Confirm biological replication, paired design, and filtering for the actual study. Loading multiple statistical packages does not establish that every model is fully implemented.
- Review `exp1 <- exp[, -3]` and the later reassignment `effLen = gfe$length` in `MLY_count_to_TPM.R`; they may remove a sample or overwrite lengths aligned by gene ID. The earlier TPM script also requires explicit agreement between the length-vector names and count-matrix rows. TPM calculation requires correctly matched gene lengths; DESeq2/edgeR count models require raw counts.
- The GO/KEGG files mix human/mouse and RNA/proteomics sections. Use consistent species, gene IDs, background sets, and gene objects. The filter `change != 'up'` includes stable genes as well as downregulated genes.
- Scripts include package installation, `View()`, plot output, and example-result statements. Adjust the selected sections for the execution environment.

The statistical workflows were archived without being rewritten or rerun. Generated images are excluded; plotting code is retained.
