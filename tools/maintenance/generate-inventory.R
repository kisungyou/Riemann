#!/usr/bin/env Rscript
# Regenerate installed method contracts from static source and test references.
# Usage: Rscript tools/maintenance/generate-inventory.R [package-root]
# No package functions are evaluated, and test references do not prove validity.
maintenance_script <- local({
  sourced <- unlist(lapply(sys.frames(), function(frame) frame$ofile), use.names = FALSE)
  if (length(sourced)) return(tail(sourced, 1L))
  invoked <- grep('^--file=', commandArgs(trailingOnly = FALSE), value = TRUE)
  if (!length(invoked)) stop('Run this entry point with Rscript or source().')
  sub('^--file=', '', invoked[1L])
})
source(file.path(dirname(normalizePath(maintenance_script, mustWork = TRUE)), 'paths.R'),
       local = TRUE)
args <- commandArgs(trailingOnly = TRUE)
maintenance_context <- riemann_maintenance_paths(maintenance_script,
  if (length(args)) args[1L] else NULL, require_package = TRUE)
root <- maintenance_context$package_root
relative <- function(path) substring(normalizePath(path), nchar(root) + 2L)

functions <- list()
for (file in list.files(file.path(root, 'R'), pattern = '\\.R$', full.names = TRUE)) {
  lines <- readLines(file, warn = FALSE)
  expressions <- parse(file, keep.source = TRUE)
  references <- attr(expressions, 'srcref')
  for (i in seq_along(expressions)) {
    expr <- expressions[[i]]
    if (!is.call(expr) || !is.symbol(expr[[1L]]) ||
        !as.character(expr[[1L]]) %in% c('<-', '=') || length(expr) != 3L ||
        !is.symbol(expr[[2L]]) || !is.call(expr[[3L]]) ||
        !identical(expr[[3L]][[1L]], as.name('function'))) next
    name <- as.character(expr[[2L]])
    first <- as.integer(references[[i]])[1L]
    j <- first - 1L
    doc <- character()
    while (j >= 1L && grepl('^\\s*(#|$)', lines[j])) {
      if (grepl("^\\s*#'", lines[j])) doc <- c(sub("^\\s*#' ?", '', lines[j]), doc)
      j <- j - 1L
    }
    functions[[name]] <- list(source = paste0(relative(file), ':', first),
      export = any(grepl('^@export(?:\\s|$)', doc, perl = TRUE)))
  }
}

exports <- methods <- character()
for (expr in parse(file.path(root, 'NAMESPACE'))) {
  if (!is.call(expr)) next
  name <- as.character(expr[[1L]])
  if (name == 'export') exports <- c(exports,
    vapply(as.list(expr)[-1L], as.character, ''))
  if (name == 'S3method') methods <- c(methods,
    paste(as.character(expr[[2L]]), as.character(expr[[3L]]), sep = '.'))
}
public <- sort(unique(c(exports, methods,
  names(functions)[vapply(functions, `[[`, logical(1), 'export')])))
unresolved <- setdiff(public, names(functions))
if (length(unresolved)) stop('Resolve public functions before regenerating: ',
  paste(unresolved, collapse = ', '))
test_files <- list.files(file.path(root, 'tests'), pattern = '\\.R$',
                         recursive = TRUE, full.names = TRUE)
test_names <- lapply(test_files, function(file) all.names(parse(file), unique = TRUE))
inventory <- do.call(rbind, lapply(public, function(name) {
  tests <- test_files[vapply(test_names, function(x) name %in% x, logical(1))]
  data.frame(name = name, source = functions[[name]]$source,
    direct_test_files = paste(sort(unique(vapply(tests, relative, ''))), collapse = '; '),
    stringsAsFactors = FALSE)
}))
contract_environment <- new.env(parent = globalenv())
maintenance_context$inventory <- inventory
contract_environment$riemann_contract_context <- maintenance_context
sys.source(file.path(maintenance_context$maintenance_root, 'generate-method-contracts.R'),
           envir = contract_environment)
cat(nrow(inventory), 'public exports/S3 methods; inst/method-contracts.csv regenerated.\n')
