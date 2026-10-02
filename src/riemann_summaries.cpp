#include "riemann_src.h"
#include <algorithm>
#include <cmath>
#include <limits>

namespace {

const double roundoff = 64.0 * std::numeric_limits<double>::epsilon();

bool same_point(const arma::mat& x, const arma::mat& y) {
  return x.n_rows == y.n_rows && x.n_cols == y.n_cols &&
    arma::all(arma::vectorise(x) == arma::vectorise(y));
}

bool valid_point(const std::string& mfd, const arma::mat& x,
                 const arma::mat& exemplar) {
  if (x.n_rows != exemplar.n_rows || x.n_cols != exemplar.n_cols ||
      x.n_elem == 0 || !x.is_finite()) return false;
  const double tol = 1e-7;
  if (mfd == "spd" || mfd == "correlation") {
    if (x.n_rows != x.n_cols) return false;
    const double scale = arma::abs(x).max();
    if (scale == 0.0) return false;
    const arma::mat scaled = x / scale;
    if (arma::abs(scaled - scaled.t()).max() > 2.0 * roundoff) return false;
    arma::mat factor;
    if (!arma::chol(factor, scaled / 2.0 + scaled.t() / 2.0)) return false;
    if (mfd == "correlation" && arma::max(arma::abs(x.diag() - 1.0)) > tol)
      return false;
  } else if (mfd == "sphere" || mfd == "landmark") {
    if (std::abs(arma::norm(x, "fro") - 1.0) > tol) return false;
    if (mfd == "landmark" && arma::norm(arma::mean(x, 0), 2) > tol)
      return false;
  } else if (mfd == "grassmann" || mfd == "stiefel" || mfd == "rotation") {
    if (arma::norm(x.t() * x - arma::eye<arma::mat>(x.n_cols, x.n_cols), "fro") >
        tol * std::max(1.0, std::sqrt(static_cast<double>(x.n_cols)))) return false;
    if (mfd == "rotation" && (x.n_rows != x.n_cols || arma::det(x) <= 0.0))
      return false;
  } else if (mfd == "multinomial") {
    if (arma::any(arma::vectorise(x) < 0.0) ||
        std::abs(arma::accu(x) - 1.0) > tol) return false;
  }
  return true;
}

arma::vec validate(const std::string& mfd, const arma::field<arma::mat>& data,
                   arma::vec weights, const RiemannSummaryControl& control) {
  if (data.n_elem == 0 || weights.n_elem != data.n_elem || !weights.is_finite() ||
      arma::any(weights < 0.0) || weights.max() <= 0.0)
    Rcpp::stop("Summary weights must be finite, nonnegative, and have positive sum.");
  weights /= weights.max();
  weights /= arma::accu(weights);
  if (control.maxiter < 1 || !std::isfinite(control.tolerance) ||
      control.tolerance <= 0.0 || control.max_backtrack < 1)
    Rcpp::stop("Summary controls require positive iterations, tolerance, and backtracking limit.");
  for (arma::uword i = 0; i < data.n_elem; ++i) {
    if (!valid_point(mfd, data(i), data(0)))
      Rcpp::stop("Summary observation %d is invalid or has incompatible dimensions.",
                 static_cast<int>(i + 1));
  }
  return weights;
}

double distance(const std::string& mfd, const arma::mat& x, const arma::mat& y,
                bool extrinsic = false) {
  if (same_point(x, y)) return 0.0;
  const double value = extrinsic ? riem_distext(mfd, x, y) : riem_dist(mfd, x, y);
  if (!std::isfinite(value) || value < 0.0)
    Rcpp::stop("A geometry operation returned a nonfinite or negative distance.");
  return value;
}

double objective(const std::string& mfd, const arma::mat& x,
                 const arma::field<arma::mat>& data, const arma::vec& weights,
                 bool squared, bool extrinsic = false, arma::vec* distances = NULL) {
  double total = 0.0;
  if (distances != NULL) {
    distances->set_size(data.n_elem);
    distances->fill(NA_REAL);
  }
  for (arma::uword i = 0; i < data.n_elem; ++i) {
    if (weights(i) == 0.0) continue;
    const double d = distance(mfd, x, data(i), extrinsic);
    if (distances != NULL) (*distances)(i) = d;
    total += weights(i) * (squared ? d * d : d);
  }
  if (!std::isfinite(total)) Rcpp::stop("The summary objective is not finite.");
  return total;
}

double tangent_norm(const std::string& mfd, const arma::mat& x, const arma::mat& u) {
  if (!u.is_finite()) Rcpp::stop("A geometry operation returned a nonfinite tangent.");
  const double squared = riem_metric(mfd, x, u, u);
  if (!std::isfinite(squared) || squared < 0.0)
    Rcpp::stop("A geometry operation returned an invalid tangent norm.");
  return std::sqrt(squared);
}

void trace_add(RiemannSummaryResult& fit, const RiemannSummaryControl& control,
               double stationarity = NA_REAL, double alpha = 0.0) {
  if (!control.keep_trace) return;
  fit.trace_objective.push_back(fit.objective);
  fit.trace_stationarity.push_back(stationarity);
  fit.trace_step.push_back(fit.step_norm);
  fit.trace_alpha.push_back(alpha);
}

void trace_stationarity(RiemannSummaryResult& fit,
                        const RiemannSummaryControl& control, double value) {
  if (control.keep_trace && !fit.trace_stationarity.empty())
    fit.trace_stationarity.back() = value;
}

arma::mat initialize(const std::string& mfd, const arma::field<arma::mat>& data,
                     const arma::vec& weights, const arma::mat* initial,
                     RiemannSummaryResult& fit) {
  arma::mat x;
  if (initial != NULL) {
    x = *initial;
    fit.initialization = "supplied";
  } else {
    // Do not send zero-weight observations to an initializer that might form
    // undefined intermediate values before multiplying by their weights.
    const arma::uvec support = arma::find(weights > 0.0);
    arma::field<arma::mat> active(support.n_elem);
    arma::vec active_weights(support.n_elem);
    for (arma::uword j = 0; j < support.n_elem; ++j) {
      active(j) = data(support(j));
      active_weights(j) = weights(support(j));
    }
    fit.initialization = "weighted_geometry_initializer";
    try {
      x = riem_initialize(mfd, active, active_weights);
    } catch (const std::exception&) {
      x = active(0);
      fit.initialization = "first_positive_weight_observation";
    }
    if (!valid_point(mfd, x, data(0))) {
      x = active(0);
      fit.initialization = "first_positive_weight_observation";
    }
  }
  if (!valid_point(mfd, x, data(0)))
    Rcpp::stop("The supplied summary initializer is not a valid compatible observation.");
  return x;
}

bool all_identical(const arma::field<arma::mat>& data, const arma::vec& weights,
                   arma::uword& first) {
  const arma::uvec support = arma::find(weights > 0.0);
  first = support(0);
  for (arma::uword j = 1; j < support.n_elem; ++j)
    if (!same_point(data(first), data(support(j)))) return false;
  return true;
}

struct MedianDirection {
  arma::mat residual;
  double residual_norm;
  double coincident_weight;
  double inverse_distance_sum;
  double stationarity;
  arma::uword nearest;
};

MedianDirection median_direction(const std::string& mfd, const arma::mat& x,
                                  const arma::field<arma::mat>& data,
                                  const arma::vec& weights) {
  MedianDirection out;
  out.residual.zeros(x.n_rows, x.n_cols);
  out.coincident_weight = 0.0;
  out.inverse_distance_sum = 0.0;
  double nearest_distance = std::numeric_limits<double>::infinity();
  out.nearest = 0;
  for (arma::uword i = 0; i < data.n_elem; ++i) {
    if (weights(i) == 0.0) continue;
    const double d = distance(mfd, x, data(i));
    if (d < nearest_distance) {
      nearest_distance = d;
      out.nearest = i;
    }
    // Only genuine coincidences contribute subgradient mass. Never discard
    // an observation merely because its positive distance is small.
    if (d == 0.0) {
      out.coincident_weight += weights(i);
    } else {
      const arma::mat log = riem_log(mfd, x, data(i));
      if (!log.is_finite()) Rcpp::stop("The median logarithm is not finite.");
      out.residual += weights(i) * (log / d);
      out.inverse_distance_sum += weights(i) / d;
    }
  }
  out.residual = riem_project_tangent(mfd, x, out.residual);
  out.residual_norm = tangent_norm(mfd, x, out.residual);
  out.stationarity = std::max(0.0, out.residual_norm - out.coincident_weight);
  return out;
}

bool backtrack(const std::string& mfd, const arma::field<arma::mat>& data,
                const arma::vec& weights, const RiemannSummaryControl& control,
                const arma::mat& direction, double decrease,
                bool squared, RiemannSummaryResult& fit) {
  const double comparison_tolerance = roundoff * std::abs(fit.objective);
  double alpha = 1.0;
  bool found = false;
  arma::mat best;
  double best_objective = fit.objective;
  double best_step = 0.0;
  double best_alpha = 0.0;
  for (int back = 0; back < control.max_backtrack; ++back, alpha *= 0.5) {
    try {
      arma::mat candidate = riem_exp(mfd, fit.estimate, direction, alpha);
      if (!valid_point(mfd, candidate, data(0))) continue;
      const double candidate_objective = objective(mfd, candidate, data, weights, squared);
      if (candidate_objective <= fit.objective - 1e-4 * alpha * decrease + comparison_tolerance) {
        const double step = distance(mfd, fit.estimate, candidate);
        if (found && candidate_objective >= best_objective - comparison_tolerance) break;
        found = true;
        best = candidate;
        best_objective = candidate_objective;
        best_step = step;
        best_alpha = alpha;
        // An admissible full step can oscillate with barely decreasing cost.
        // Check finer steps while they improve the cost; every retained point
        // still satisfies the same Armijo rule and bounded trial budget.
      } else if (found) {
        break;
      }
    } catch (const std::exception&) {
      // Keep the last accepted point when a trial leaves the supported domain.
    }
  }
  if (found) {
    fit.estimate = best;
    fit.objective = best_objective;
    fit.step_norm = best_step;
    ++fit.iterations;
    trace_add(fit, control, NA_REAL, best_alpha);
    return true;
  }
  fit.termination = "line_search_failed";
  return false;
}

} // namespace

RiemannSummaryResult riem_summary_mean(const std::string& mfd,
  const arma::field<arma::mat>& data, arma::vec weights,
  const RiemannSummaryControl& control, const arma::mat* initial) {
  if (mfd == "stiefel" || mfd == "correlation")
    Rcpp::stop("Intrinsic summaries are unavailable for %s pending a compatible geometry implementation.", mfd.c_str());
  weights = validate(mfd, data, weights, control);
  RiemannSummaryResult fit;
  fit.estimand = "intrinsic_frechet_mean";
  arma::uword first = 0;
  if (mfd == "euclidean" || all_identical(data, weights, first)) {
    fit.estimate = (mfd == "euclidean") ?
      riem_initialize(mfd, data, weights) : data(first);
    fit.objective = objective(mfd, fit.estimate, data, weights, true, false, &fit.distances);
    fit.converged = true;
    fit.termination = "closed_form";
    fit.initialization = "closed_form";
    fit.gradient_norm = 0.0;
    trace_add(fit, control, 0.0);
    return fit;
  }
  fit.estimate = initialize(mfd, data, weights, initial, fit);
  fit.objective = objective(mfd, fit.estimate, data, weights, true);
  trace_add(fit, control);
  while (true) {
    Rcpp::checkUserInterrupt();
    try {
      arma::mat direction(fit.estimate.n_rows, fit.estimate.n_cols, arma::fill::zeros);
      for (arma::uword i = 0; i < data.n_elem; ++i) {
        if (weights(i) == 0.0 || same_point(fit.estimate, data(i))) continue;
        direction += weights(i) * riem_log(mfd, fit.estimate, data(i));
      }
      direction = riem_project_tangent(mfd, fit.estimate, direction);
      const double norm = tangent_norm(mfd, fit.estimate, direction);
      fit.gradient_norm = 2.0 * norm;
      trace_stationarity(fit, control, fit.gradient_norm);
      if (fit.gradient_norm <= control.tolerance) {
        fit.converged = true;
        fit.termination = "stationary";
        break;
      }
      if (fit.iterations >= control.maxiter) break;
      if (fit.iterations > 0 && fit.step_norm <= roundoff * std::sqrt(fit.objective)) {
        fit.termination = "stagnation";
        break;
      }
      if (!backtrack(mfd, data, weights, control, direction, 2.0 * norm * norm, true, fit))
        break;
    } catch (const std::exception& error) {
      fit.termination = "invalid_geometry_operation";
      fit.message = error.what();
      break;
    }
  }
  fit.objective = objective(mfd, fit.estimate, data, weights, true, false, &fit.distances);
  return fit;
}

RiemannSummaryResult riem_summary_median(const std::string& mfd,
  const arma::field<arma::mat>& data, arma::vec weights,
  const RiemannSummaryControl& control, const arma::mat* initial) {
  if (mfd == "stiefel" || mfd == "correlation")
    Rcpp::stop("Intrinsic summaries are unavailable for %s pending a compatible geometry implementation.", mfd.c_str());
  weights = validate(mfd, data, weights, control);
  RiemannSummaryResult fit;
  fit.estimand = "intrinsic_frechet_median";
  arma::uword first = 0;
  if (all_identical(data, weights, first) || weights.max() >= 0.5) {
    // Triangle inequality certifies an observation with at least half the
    // weight as a minimizer, including cases where a logarithm is nonunique.
    if (weights.max() >= 0.5) first = weights.index_max();
    fit.estimate = data(first);
    fit.objective = objective(mfd, fit.estimate, data, weights, false, false, &fit.distances);
    fit.converged = true;
    fit.termination = "closed_form";
    fit.initialization = "weighted_majority_or_identical";
    trace_add(fit, control);
    return fit;
  }
  fit.estimate = initialize(mfd, data, weights, initial, fit);
  fit.objective = objective(mfd, fit.estimate, data, weights, false);
  trace_add(fit, control);
  while (true) {
    Rcpp::checkUserInterrupt();
    try {
      MedianDirection direction = median_direction(mfd, fit.estimate, data, weights);
      fit.subgradient_residual = direction.stationarity;
      trace_stationarity(fit, control, fit.subgradient_residual);
      if (fit.subgradient_residual <= control.tolerance) {
        fit.converged = true;
        fit.termination = "stationary";
        break;
      }
      if (fit.iterations >= control.maxiter) break;
      // A smooth iteration may approach a nonsmooth optimum without ever
      // landing exactly on it. Test the nearest observation itself; accept it
      // only with its own subgradient certificate and a nonincreasing cost.
      if (!same_point(fit.estimate, data(direction.nearest))) {
        try {
          const arma::mat& candidate = data(direction.nearest);
          const MedianDirection at_observation = median_direction(mfd, candidate, data, weights);
          if (at_observation.stationarity <= control.tolerance) {
            const double candidate_objective = objective(mfd, candidate, data, weights, false);
            if (candidate_objective <= fit.objective + roundoff * fit.objective) {
              fit.step_norm = distance(mfd, fit.estimate, candidate);
              fit.estimate = candidate;
              fit.objective = candidate_objective;
              fit.subgradient_residual = at_observation.stationarity;
              ++fit.iterations;
              fit.converged = true;
              fit.termination = "stationary_at_observation";
              trace_add(fit, control, fit.subgradient_residual, 1.0);
              break;
            }
          }
        } catch (const std::exception&) {
          // A candidate at an observation can have a cut-locus ambiguity even
          // when the current iterate admits a valid direction.
        }
      }
      if (fit.iterations > 0 && fit.step_norm <= roundoff * fit.objective) {
        fit.termination = "stagnation";
        break;
      }
      if (!std::isfinite(direction.inverse_distance_sum) || direction.inverse_distance_sum <= 0.0)
        Rcpp::stop("The median update has an invalid inverse-distance weight sum.");
      const double multiplier = (1.0 - direction.coincident_weight / direction.residual_norm) /
        direction.inverse_distance_sum;
      const arma::mat step_direction = multiplier * direction.residual;
      const double norm = tangent_norm(mfd, fit.estimate, step_direction);
      if (!backtrack(mfd, data, weights, control, step_direction,
                     direction.stationarity * norm, false, fit)) break;
    } catch (const std::exception& error) {
      fit.termination = "invalid_geometry_operation";
      fit.message = error.what();
      break;
    }
  }
  fit.objective = objective(mfd, fit.estimate, data, weights, false, false, &fit.distances);
  return fit;
}

RiemannSummaryResult riem_summary_extrinsic(const std::string& mfd,
  const arma::field<arma::mat>& data, arma::vec weights,
  const RiemannSummaryControl& control, bool median, const arma::mat* initial) {
  weights = validate(mfd, data, weights, control);
  if (mfd == "landmark")
    Rcpp::stop("Extrinsic landmark summaries are unavailable: raw embedding and Procrustes distance disagree.");
  const arma::uvec support = arma::find(weights > 0.0);
  const int rows = data(0).n_rows;
  const int cols = data(0).n_cols;
  const arma::vec exemplar = riem_equiv(mfd, data(support(0)), rows, cols);
  if (exemplar.n_elem == 0 || !exemplar.is_finite())
    Rcpp::stop("The equivariant embedding is empty or nonfinite.");
  arma::field<arma::mat> embedded(data.n_elem);
  for (arma::uword i = 0; i < data.n_elem; ++i) {
    if (weights(i) == 0.0) {
      embedded(i).zeros(exemplar.n_elem, 1);
      continue;
    }
    const arma::vec point = riem_equiv(mfd, data(i), rows, cols);
    if (point.n_elem != exemplar.n_elem || !point.is_finite())
      Rcpp::stop("Equivariant embeddings have incompatible sizes or nonfinite entries.");
    embedded(i) = point;
  }
  arma::mat embedded_initial;
  const arma::mat* initial_pointer = NULL;
  if (initial != NULL) {
    if (!valid_point(mfd, *initial, data(0))) Rcpp::stop("The supplied initializer is invalid.");
    embedded_initial = riem_equiv(mfd, *initial, rows, cols);
    initial_pointer = &embedded_initial;
  }
  RiemannSummaryResult fit = median ?
    riem_summary_median("euclidean", embedded, weights, control, initial_pointer) :
    riem_summary_mean("euclidean", embedded, weights, control, initial_pointer);
  fit.ambient_objective = fit.objective;
  fit.ambient_subgradient_residual = fit.subgradient_residual;
  fit.estimate = riem_invequiv(mfd, arma::vectorise(fit.estimate), rows, cols);
  if (!valid_point(mfd, fit.estimate, data(0)))
    Rcpp::stop("The inverse embedding did not return a valid manifold point.");
  fit.objective = objective(mfd, fit.estimate, data, weights, !median, true, &fit.distances);
  const bool full_chart = mfd == "spd" || mfd == "euclidean";
  fit.estimand = median ? (full_chart ? "extrinsic_frechet_median" : "projected_ambient_geometric_median") :
    "extrinsic_frechet_mean";
  fit.diagnostic_scope = full_chart ? "embedding_chart" : "ambient_embedding";
  if (!full_chart) {
    fit.gradient_norm = NA_REAL;
    fit.subgradient_residual = NA_REAL;
  }
  return fit;
}
