#!/usr/bin/env Rscript
# Usage: Rscript tools/maintenance/verify-package.R PACKAGE_SOURCE LIBRARY OUTPUT
# Always select an explicit installed library; a global older Riemann is unsafe.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) stop("Supply package source, installed library, and output directory.")
source_root <- normalizePath(args[1L], mustWork = TRUE)
library_root <- normalizePath(args[2L], mustWork = TRUE)
output <- args[3L]
dir.create(output, recursive = TRUE, showWarnings = FALSE)
output <- normalizePath(output, mustWork = TRUE)
.libPaths(c(library_root, .libPaths()))
library(Riemann, lib.loc = library_root)
library(testthat)
source_version <- read.dcf(file.path(source_root, "DESCRIPTION"))[1L, "Version"]
stopifnot(as.character(packageVersion("Riemann")) == source_version,
          normalizePath(find.package("Riemann")) == file.path(library_root, "Riemann"))
writeLines(capture.output(sessionInfo()), file.path(output, "session-info.txt"))
result <- test_dir(file.path(source_root, "tests", "testthat"),
                   package = "Riemann", load_package = "installed",
                   reporter = "summary", stop_on_failure = FALSE,
                   stop_on_warning = FALSE)
saveRDS(result, file.path(output, "results.rds"))
table <- as.data.frame(result)
table$result <- NULL
write.csv(table, file.path(output, "results.csv"), row.names = FALSE)
summary <- c(test_blocks = nrow(table), assertions = sum(table$nb),
             passed = sum(table$passed), failed = sum(table$failed),
             errors = sum(table$error), warnings = sum(table$warning),
             skipped = sum(table$skipped))
writeLines(capture.output(print(summary)), file.path(output, "summary.txt"))
print(summary)
if (any(summary[c("failed", "errors", "warnings", "skipped")] != 0L)) {
  stop("The complete suite did not finish cleanly; inspect results.csv and captured output.")
}
