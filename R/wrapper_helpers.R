# Common ingestion keeps singleton matrix dimensions and feature order intact.
riem_matrix_input <- function(input, square = FALSE) {
  if (check_3darray(input, symmcheck = square)) {
    dims <- dim(input)[1:2]
    labels <- dimnames(input)[1:2]
    if (is.null(dimnames(input))) labels <- NULL
    data <- lapply(seq_len(dim(input)[3L]), function(i) {
      matrix(input[, , i], dims[1L], dims[2L], dimnames = labels)
    })
  } else if (is.list(input)) data <- input
  else stop("Input must be a nonempty three-dimensional array or list of matrices.", call. = FALSE)
  if (!check_list_eqsize(data, check.square = square) ||
      !all(vapply(data, is.matrix, logical(1)))) {
    stop("Observations must be nonempty matrices with matching dimensions.", call. = FALSE)
  }
  for (i in seq_along(data)) {
    x <- data[[i]]
    if (!is.numeric(x) || is.complex(x) || any(!is.finite(x))) {
      stop("Observation ", i, " must be a finite real matrix.", call. = FALSE)
    }
    if (!identical(dimnames(x), dimnames(data[[1L]]))) {
      stop("Observations must have consistent feature and landmark ordering.", call. = FALSE)
    }
  }
  data
}

riem_wrap_matrices <- function(data, manifold) {
  structure(list(data = data, size = dim(data[[1L]]), name = manifold,
                 dimnames = dimnames(data[[1L]])), class = "riemdata")
}

riem_full_column_rank <- function(x) {
  if (nrow(x) < ncol(x) || max(abs(x)) == 0) return(FALSE)
  values <- svd(x / max(abs(x)), nu = 0L, nv = 0L)$d
  min(values) > 64 * .Machine$double.eps * max(dim(x)) * max(values)
}
