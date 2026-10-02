# Internal controls shared by the public mean and median interfaces.
riem_summary_controls <- function(parameters, riemobj) {
  allowed <- c("maxiter", "eps", "max_backtrack", "trace", "init")
  if (length(parameters)) {
    supplied <- names(parameters)
    if (is.null(supplied) || anyNA(supplied) || any(!nzchar(supplied)) ||
        anyDuplicated(supplied)) {
      stop("Summary controls must have unique, nonempty names.", call. = FALSE)
    }
    unknown <- setdiff(supplied, allowed)
    if (length(unknown)) {
      stop("Unknown summary control: ", paste(unknown, collapse = ", "), call. = FALSE)
    }
  }
  controls <- list(maxiter = 50L, eps = 1e-5, max_backtrack = 50L,
                   trace = FALSE, init = NULL)
  controls[names(parameters)] <- parameters
  for (name in c("maxiter", "max_backtrack")) {
    value <- controls[[name]]
    if (!is.numeric(value) || is.complex(value) || length(value) != 1L ||
        !is.finite(value) || value < 1 || value != floor(value) ||
        value > .Machine$integer.max) {
      stop(name, " must be a finite positive integer.", call. = FALSE)
    }
    controls[[name]] <- as.integer(value)
  }
  if (controls$max_backtrack > 1024L) {
    stop("max_backtrack must not exceed 1024.", call. = FALSE)
  }
  value <- controls$eps
  if (!is.numeric(value) || is.complex(value) || length(value) != 1L ||
      !is.finite(value) || value <= 0) {
    stop("eps must be a finite positive number.", call. = FALSE)
  }
  value <- controls$trace
  if (!is.logical(value) || length(value) != 1L || is.na(value)) {
    stop("trace must be TRUE or FALSE.", call. = FALSE)
  }
  initial <- controls$init
  if (!is.null(initial) &&
      (!is.matrix(initial) || !is.numeric(initial) || is.complex(initial) ||
       any(!is.finite(initial)) || !identical(dim(initial), dim(riemobj$data[[1L]])))) {
    stop("init must be a finite real matrix with the observation dimensions.", call. = FALSE)
  }
  controls
}

riem_fit_summary <- function(riemobj, weight, geometry, parameters, statistic, call) {
  riem_validate_data(riemobj)
  resolved <- riem_resolve_geometry(riemobj, geometry, capability = statistic)
  controls <- riem_summary_controls(parameters, riemobj)
  n <- length(riemobj$data)
  weights <- if (is.null(weight)) rep(1 / n, n) else
    check_weight(weight, n, paste0("riem.", statistic))
  output <- inference_summary(
    riemobj$name, riemobj$data, weights, statistic, resolved$backend,
    controls$maxiter, controls$eps, controls$init,
    controls$max_backtrack, controls$trace)
  output$distvec <- NULL
  dimnames(output[[statistic]]) <- dimnames(riemobj$data[[1L]])
  output$task <- statistic
  output$geometry <- resolved
  output$call <- call
  output$weights <- weights
  output$nobs <- n
  output$input_dimensions <- dim(riemobj$data[[1L]])
  output$controls <- controls
  output$schema_version <- 1L
  output$package_version <- as.character(utils::packageVersion("Riemann"))
  output$objective_scale <- if (statistic == "mean")
    "normalized_weighted_sum_of_squared_distances" else "normalized_weighted_sum_of_distances"
  class(output) <- c("riem_summary", "riemfit")
  if (!isTRUE(output$converged)) {
    warning("The ", statistic, " solver did not converge (", output$termination,
            "); the last accepted estimate is returned.", call. = FALSE)
  }
  output
}
