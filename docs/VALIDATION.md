# Publication checks

Validation date: 2026-10-06.

| Check | Result |
| --- | --- |
| Python syntax/AST | Six analysis scripts passed |
| Bash `-n` | Three shell scripts passed |
| R `parse(..., encoding='UTF-8')` | 13 R scripts passed |
| MLY formatting | Parsed R expressions were identical before and after formatting for all six files |
| Text encoding | UTF-8 without a BOM; LF line endings |
| JSON and source hashes | JSON files parse; SHA-256 hashes of the 24 selected publication copies agree with the manifest |
| Archive contents | 22 code files plus READMEs, dependencies/version records, methods, and provenance; no generated figures, sequencing data, result tables, caches, or third-party package libraries |
| Local links in authored Markdown | Targets verified; the supplied historical methods retain the original project-layout paths |
| Common credential-format scan | No GitHub tokens, AWS access keys, private keys, or Slack token formats found |
| English documentation update | All README and supporting-document prose is in English; the supplied methods retain only their existing English section and English appendices |

These are static packaging checks. RNA-seq, DESeq2/ULM, variant-evidence analysis, and TCGA retrieval/scoring were not rerun. Syntax checks do not establish that data and dependencies are available or that the original study results have been reproduced. Cohort-specific validation code remains available for execution with the corresponding inputs.

The Windows R startup reported that the host's `C.UTF-8` locale was unavailable. Setting `English_United States.utf8` explicitly allowed all 13 files to parse. INK4's `run_figures.R` also handles UTF-8 explicitly. No packages were installed or upgraded during packaging.

The documentation update preserves all 22 analysis code files byte for byte. Its Markdown links, manifest paths/hashes, English-only prose, and rebuilt archive were checked separately.
