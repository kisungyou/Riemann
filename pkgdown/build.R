#!/usr/bin/env Rscript
# Run from the repository root after installing this checkout and its suggested
# packages. Select an isolated installation with R_LIBS_USER when appropriate.
stopifnot(file.exists("DESCRIPTION"), dir.exists("pkgdown"))
library(Riemann)
expected <- read.dcf("DESCRIPTION", fields = "Version")[1L, 1L]
if (as.character(packageVersion("Riemann")) != expected) {
  stop("Install this checkout before rebuilding its website.")
}
rmarkdown::render("README.Rmd", quiet = TRUE)
readme <- readLines("README.md", warn = FALSE)
writeLines(sub("[[:blank:]]+$", "", readme), "README.md")
pkgdown::build_site(".", examples = TRUE, run_dont_run = FALSE, lazy = FALSE,
                   preview = FALSE, install = FALSE, new_process = FALSE)
source("pkgdown/normalize-html.R")
