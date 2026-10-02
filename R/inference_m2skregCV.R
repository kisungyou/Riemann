#' Manifold-to-Scalar Kernel Regression with K-Fold Cross Validation
#'
#' Selects a Gaussian kernel bandwidth by the sum of squared held-out prediction
#' errors over all folds, then fits the selected smoother to all observations.
#' A candidate is eligible only if every fold has a finite loss. Failed candidates
#' remain in the returned tables with infinite total loss and a recorded reason.
#' If no candidate succeeds, fitting stops with an error. Equal losses are resolved
#' by choosing the first candidate in the supplied order.
#'
#' @inheritParams riem.m2skreg
#' @param bandwidths A nonempty vector of finite, strictly positive bandwidths.
#' @param kfold Integer number of folds between two and the sample size. Random
#'   balanced folds use R's current random-number state. When \code{foldid} is
#'   supplied, an explicitly supplied \code{kfold} must agree with its number
#'   of distinct folds.
#' @param foldid Optional vector of fold identifiers, one per observation, with
#'   at least two distinct nonmissing identifiers. Numeric identifiers must be
#'   finite; character and factor identifiers must be nonempty. Observations with
#'   the same identifier are held out together. Use this to encode grouped or
#'   otherwise scientifically appropriate validation splits.
#'
#' @details Distances are computed once under the resolved geometry. This is
#'   appropriate only when preprocessing defining that geometry was fixed
#'   independently of the held-out observations. This function does not fit
#'   data-dependent alignment, scaling, or other preprocessing inside each fold;
#'   such pipelines require an external resampling loop. The selected CV loss is
#'   a tuning criterion, not an unbiased estimate of final predictive performance.
#'
#' @return An \code{m2skreg} fit retaining \code{ypred}, \code{bandwidth}, and
#'   \code{inputs}, with the geometry metadata of \code{riem.m2skreg}.
#'   \code{ypred} contains full-training fitted values. \code{errors} is a
#'   two-column matrix of all candidate bandwidths and total CV SSE;
#'   \code{fold_errors} and \code{fold_failure} have candidates in rows and folds
#'   in columns. \code{candidate_status} records success or failure,
#'   \code{foldid} stores the supplied/generated fold identifiers, and
#'   \code{cv_prediction} contains held-out predictions for the selected candidate.
#'
#' @examples
#' X <- wrap.euclidean(matrix(0:5, ncol = 1))
#' fit <- riem.m2skregCV(X, c(0, 3, 0, 0, 0, 0),
#'   bandwidths = c(0.1, 0.5, 1, 2, 10), foldid = rep(1:3, each = 2))
#' fit$errors
#' fit$cv_prediction
#'
#' @concept inference
#' @export
riem.m2skregCV <- function(riemobj, y,
                          bandwidths = seq(0.01, 1, length.out = 10),
                          geometry = NULL, kfold = 5, foldid = NULL) {
  riem_validate_data(riemobj)
  n <- length(riemobj$data)
  if (n < 2L) stop("Cross-validation requires at least two observations.", call. = FALSE)
  y <- riem_regression_response(y, n)
  bandwidths <- riem_regression_bandwidth(bandwidths, multiple = TRUE)
  geometry <- riem_resolve_geometry(riemobj, geometry, capability = "distance")
  if (is.null(foldid)) {
    kfold <- riem_regression_integer(kfold, "kfold", 2L, n)
    foldid <- sample(rep(seq_len(kfold), length.out = n))
  } else {
    if (length(foldid) != n || !is.null(dim(foldid)) || anyNA(foldid) ||
        !(is.numeric(foldid) || is.character(foldid) || is.factor(foldid)) ||
        is.complex(foldid) ||
        (is.numeric(foldid) && any(!is.finite(foldid))) ||
        ((is.character(foldid) || is.factor(foldid)) &&
         any(!nzchar(trimws(as.character(foldid)))))) {
      stop("'foldid' must give one finite, nonmissing, nonempty identifier per observation.",
           call. = FALSE)
    }
    nfold <- length(unique(foldid))
    if (nfold < 2L) stop("'foldid' must contain at least two distinct folds.", call. = FALSE)
    if (!missing(kfold) && riem_regression_integer(kfold, "kfold", 2L, n) != nfold) {
      stop("'kfold' does not match the number of distinct 'foldid' values.", call. = FALSE)
    }
    kfold <- nfold
  }
  fold_labels <- unique(foldid)
  splitgp <- split(seq_len(n), factor(match(foldid, fold_labels), levels = seq_len(kfold)))
  distances <- basic_pdist(riemobj$name, riemobj$data, geometry$backend)
  ncandidate <- length(bandwidths)
  fold_errors <- matrix(Inf, ncandidate, kfold,
                        dimnames = list(NULL, as.character(fold_labels)))
  fold_failure <- matrix(NA_character_, ncandidate, kfold,
                         dimnames = dimnames(fold_errors))
  best_sse <- Inf
  best_index <- NA_integer_
  best_prediction <- NULL
  for (i in seq_along(bandwidths)) {
    heldout_prediction <- rep(NA_real_, n)
    for (j in seq_along(splitgp)) {
      ids <- splitgp[[j]]
      trial <- tryCatch({
        result <- riem_kernel_predict_distances(distances[-ids, ids, drop = FALSE],
                                                y[-ids], bandwidths[i])
        loss <- sum((result$prediction - y[ids])^2)
        if (!is.finite(loss)) stop("The held-out squared-error loss is nonfinite.")
        list(prediction = result$prediction, loss = loss)
      }, error = function(e) e)
      if (inherits(trial, "error")) {
        fold_failure[i, j] <- conditionMessage(trial)
      } else {
        fold_errors[i, j] <- trial$loss
        heldout_prediction[ids] <- trial$prediction
      }
    }
    total <- sum(fold_errors[i, ])
    if (is.finite(total) && total < best_sse) {
      best_sse <- total
      best_index <- i
      best_prediction <- heldout_prediction
    }
  }
  totals <- rowSums(fold_errors)
  if (is.na(best_index)) {
    reason <- unique(stats::na.omit(as.vector(fold_failure)))
    if (!length(reason)) reason <- "The summed squared-error loss is nonfinite."
    stop(paste0("No bandwidth has finite errors in all folds. ", reason[1L]), call. = FALSE)
  }
  selected_bandwidth <- bandwidths[best_index]
  result <- riem_kernel_predict_distances(distances, y, selected_bandwidth)
  output <- riem_regression_object(riemobj, y, selected_bandwidth, geometry,
                                   result, match.call())
  output$errors <- cbind(bandwidth = bandwidths, SSE = totals)
  output$fold_errors <- fold_errors
  output$fold_failure <- fold_failure
  output$candidate_status <- ifelse(is.finite(totals), "ok", "failed")
  output$foldid <- foldid
  output$cv_prediction <- best_prediction
  output$cv <- list(loss = "SSE", selected_index = best_index,
                    selected_sse = best_sse, fold_labels = fold_labels,
                    tie_break = "first_candidate", preprocessing = "fixed")
  output
}

# Internal compatibility helper, also useful for direct SSE reference checks.
riem.m2skregCV.each <- function(pdistmat, y, id.now, bandwidth) {
  y <- riem_regression_response(y, nrow(pdistmat))
  bandwidth <- riem_regression_bandwidth(bandwidth)
  if (!is.numeric(id.now) || !length(id.now) || any(!is.finite(id.now)) ||
      any(id.now != floor(id.now)) || any(id.now < 1 | id.now > length(y)) ||
      anyDuplicated(id.now) || length(id.now) >= length(y)) {
    stop("Held-out indices must be a nonempty proper subset of the observations.", call. = FALSE)
  }
  result <- riem_kernel_predict_distances(pdistmat[-id.now, id.now, drop = FALSE],
                                          y[-id.now], bandwidth)
  loss <- sum((result$prediction - y[id.now])^2)
  if (!is.finite(loss)) stop("The held-out squared-error loss is nonfinite.", call. = FALSE)
  loss
}
