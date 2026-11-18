# Resolve maintenance inputs relative to this checkout, independent of getwd().
riemann_maintenance_paths <- function(script_file, package_root = NULL,
                                      require_package = FALSE) {
  maintenance_root <- dirname(normalizePath(script_file, mustWork = TRUE))
  if (is.null(package_root))
    package_root <- file.path(maintenance_root, '..', '..')
  package_root <- normalizePath(package_root, mustWork = require_package)
  if (require_package) {
    description <- file.path(package_root, 'DESCRIPTION')
    if (!file.exists(description)) stop('Supply the Riemann package root.')
    metadata <- read.dcf(description)
    if (!'Package' %in% colnames(metadata) || metadata[1L, 'Package'] != 'Riemann')
      stop('The package root must contain the Riemann DESCRIPTION.')
  }
  list(maintenance_root = maintenance_root, package_root = package_root)
}
