# Source provenance and publication adjustments

All 11 supplied files are represented in the archive: three shell scripts in `upstream/`, six MLY R scripts in `downstream/`, the pan-cancer Python script in `TP53_analysis/pan_cancer/`, and the English section of the supplied methods in `docs/MANUSCRIPT_METHODS_EN.md`. The checkpoint script follows the current local layout in `upstream/`.

`TP53_analysis/` includes resource preparation, expression modeling, activity inference, MDM2 comparison, summary plotting, independent validation, and the historical Windows installation helper. `INK4_CTNNB1_analysis/` includes the preparation, plotting, UTF-8 entry point, reporting, and validation scripts from `corrected_exact_P0`. The earlier duplicate INK4 version is excluded.

Text uses UTF-8 without a BOM and LF line endings. The supplied source files were left intact. Trailing whitespace and excess final blank lines were removed from the MLY publication copies; their parsed R expressions were verified to be identical before and after this formatting change. TP53 input paths use `TP53_INPUT_DIR`, and R modules locate themselves from the script directory. INK4 requires an explicit `--input`. Its renderer generates all six figures, and its validator omits comparisons against old images and replaces the historical visual-review claim with a requirement for manual review. Statistical models, filtering/candidate thresholds, exact-P0 logic, and cohort assertions were retained.

Documentation is provided in English. The Chinese section of the supplied bilingual methods was omitted from the publication copy, and the existing English section and English appendices were preserved. Original project paths and historical analytical claims in that methods document remain unchanged.

Figures, outputs, caches, third-party libraries, public reference sequences, and result archives are excluded. Historical R package versions and public reference checksums are retained as reproduction metadata. Both MLY TPM implementations remain available, with their required adaptations described in the module README.

[`source_manifest.json`](source_manifest.json) records each selected source, its original SHA-256, the publication-copy SHA-256, and the adjustments made. Source names are logical labels rather than personal computer or messaging-app paths.
