## Export the supplied ERP array without loading or installing a package.
## Arguments: explicit ERP.rda path; output directory.
## Supplied under the Riemann package's GPL-3 terms.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop('Usage: Rscript export-reference.R ERP.rda output-directory')
input <- normalizePath(args[1L], mustWork = TRUE)
dir.create(args[2L], recursive = TRUE, showWarnings = FALSE)
output <- normalizePath(args[2L], mustWork = TRUE)
env <- new.env(parent = emptyenv())
loaded <- load(input, envir = env)
stopifnot('ERP' %in% loaded, identical(dim(env$ERP$covariance), c(32L, 32L, 216L)),
          length(env$ERP$label) == 216L, all(is.finite(env$ERP$covariance)))
array_file <- file.path(output, 'package-erp.float64')
label_file <- file.path(output, 'package-erp-labels.txt')
writeBin(as.double(env$ERP$covariance), array_file, size = 8L, endian = 'little')
writeLines(as.character(env$ERP$label), label_file)
write.csv(data.frame(file = c(input, array_file, label_file),
                    md5 = unname(tools::md5sum(c(input, array_file, label_file)))),
          file.path(output, 'reference-export.csv'), row.names = FALSE)
cat('Exported unchanged 32 x 32 x 216 array in R column-major order.\n')
print(table(env$ERP$label))
