#' Inspect Wrapped Observations and Geometric Summaries
#'
#' Standard methods report the observation space, selected geometry, and actual
#' numerical termination. A finite estimate is not by itself proof of convergence.
#' Summary plots display matrix entries or ambient coordinates; they do not imply
#' an isometric view of the manifold.
#' @param x,object A wrapped dataset or fitted geometric summary.
#' @param ... Additional arguments passed to the plotting function where relevant.
#' @return Print methods return their object invisibly. Summary methods return
#'   structured metadata and numerical diagnostics. Plot methods return invisibly.
#' @name riem-inspect
NULL

#' @rdname riem-inspect
#' @method print riemdata
#' @export
print.riemdata <- function(x, ...) {
  cat(length(x$data), x$name, "observations; representation",
      paste(x$size, collapse = " x "), "\n")
  if (!is.null(x$geometry)) print(x$geometry)
  invisible(x)
}

#' @rdname riem-inspect
#' @method summary riemdata
#' @export
summary.riemdata <- function(object, ...) {
  riem_validate_data(object)
  out <- list(manifold = object$name, observations = length(object$data),
              representation_dim = object$size,
              geometry = riem.geometry(object))
  if (object$name == "spd") {
    out$eigenvalue_ratio <- vapply(object$data, function(x) {
      scale <- max(abs(x))
      values <- eigen(x / scale, symmetric = TRUE, only.values = TRUE)$values
      min(values) / max(values)
    }, numeric(1))
  }
  out
}

#' @rdname riem-inspect
#' @method print riem_summary
#' @export
print.riem_summary <- function(x, ...) {
  cat("Geometric", x$task, "under", x$geometry$geometry_id, "on",
      x$geometry$manifold_id, "\n")
  cat("Termination:", x$termination, "; converged:", isTRUE(x$converged),
      "; iterations:", x$iterations, "\n")
  cat("Objective:", format(x$objective, digits = 8), "\n")
  if (!is.null(x$estimand)) cat("Estimand:", x$estimand, "\n")
  print(if (x$task == "mean") x$mean else x$median)
  invisible(x)
}

#' @rdname riem-inspect
#' @method summary riem_summary
#' @export
summary.riem_summary <- function(object, ...) {
  if (!identical(object$schema_version, 1L)) {
    stop("Unsupported summary schema; refit the summary with this package version.", call. = FALSE)
  }
  object
}

#' @rdname riem-inspect
#' @method plot riem_summary
#' @export
plot.riem_summary <- function(x, ...) {
  estimate <- if (x$task == "mean") x$mean else x$median
  if (ncol(estimate) > 1L) {
    graphics::image(estimate, main = paste(x$task, "matrix entries"), ...)
  } else {
    graphics::plot(seq_along(estimate), as.numeric(estimate), type = "h",
                   xlab = "Ambient coordinate", ylab = "Value",
                   main = paste(x$geometry$manifold_id, x$task), ...)
  }
  invisible(x)
}
