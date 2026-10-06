# Read the drawing script as UTF-8 explicitly, independently of Windows locale.
args <- commandArgs(trailingOnly = FALSE)
script_arg <- sub('^--file=', '', args[grepl('^--file=', args)])
root <- normalizePath(dirname(script_arg), winslash = '/', mustWork = TRUE)
if (.Platform$OS.type == 'windows') {
  suppressWarnings(Sys.setlocale('LC_CTYPE', 'Chinese (Simplified)_China.utf8'))
}
source(file.path(root, 'render_figures.R'), encoding = 'UTF-8')
