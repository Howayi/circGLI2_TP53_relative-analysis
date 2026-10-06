dir.create('TP53_analysis/R_library', recursive = TRUE, showWarnings = FALSE)
install.packages('TP53_analysis/resources/decoupleR_2.17.0.zip',
                 repos = NULL, type = 'win.binary', lib = 'TP53_analysis/R_library')
.libPaths(c(normalizePath('TP53_analysis/R_library'), .libPaths()))
library(decoupleR)
print(packageVersion('decoupleR'))
print(formals(run_ulm))
