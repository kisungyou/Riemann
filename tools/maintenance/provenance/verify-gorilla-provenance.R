# Reads package data and prints evidence; never modifies datasets.
# Rscript tools/maintenance/provenance/verify-gorilla-provenance.R [package-root]
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 1L) stop('Supply at most one package-root argument.')
if (length(args)) {
  root <- normalizePath(args[1L], mustWork = TRUE)
} else {
  invoked <- grep('^--file=', commandArgs(trailingOnly = FALSE), value = TRUE)
  if (!length(invoked)) stop('Run with Rscript or supply a package-root argument.')
  script <- normalizePath(sub('^--file=', '', invoked[1L]), mustWork = TRUE)
  root <- normalizePath(file.path(dirname(script), '..', '..', '..'), mustWork = TRUE)
}
stopifnot(file.exists(file.path(root, 'DESCRIPTION')))
packaged <- new.env(parent = emptyenv())
upstream <- new.env(parent = emptyenv())
load(file.path(root, 'data', 'gorilla.rda'), envir = packaged)
utils::data('gorm.dat', package = 'shapes', envir = upstream)
utils::data('gorf.dat', package = 'shapes', envir = upstream)
landmark_order <- c(1L, 5L, 4L, 3L, 2L, 8L, 7L, 6L)
stopifnot(identical(packaged$gorilla$male, upstream$gorm.dat[landmark_order, , , drop = FALSE]),
          identical(packaged$gorilla$female, upstream$gorf.dat[landmark_order, , , drop = FALSE]))
cat('R:', as.character(getRversion()), '\n')
cat('shapes:', as.character(utils::packageVersion('shapes')), '\n')
cat('Landmark order:', paste(landmark_order, collapse = ','), '\n')
cat('Male: exact identical array after reordering, 8 x 2 x 29.\n')
cat('Female: exact identical array after reordering, 8 x 2 x 30.\n')
cat('No coordinate transform or subject permutation needed.\n')
print(tools::md5sum(file.path(root, 'data', 'gorilla.rda')))
