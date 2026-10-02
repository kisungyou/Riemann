#include <RcppArmadillo.h>
#include "riemann_src.h"

namespace {

arma::field<arma::mat> summary_data(const Rcpp::List& data) {
  arma::field<arma::mat> out(data.size());
  for (int i = 0; i < data.size(); ++i) out(i) = Rcpp::as<arma::mat>(data[i]);
  return out;
}

Rcpp::List summary_list(const RiemannSummaryResult& fit, bool median, bool keep_trace) {
  Rcpp::List out;
  out[median ? "median" : "mean"] = fit.estimate;
  out["variation"] = fit.objective;
  out["objective"] = fit.objective;
  out["distvec"] = fit.distances;
  out["converged"] = fit.converged;
  out["termination"] = fit.termination;
  out["iterations"] = fit.iterations;
  out["gradient_norm"] = fit.gradient_norm;
  out["subgradient_residual"] = fit.subgradient_residual;
  out["step_norm"] = fit.step_norm;
  out["initialization"] = fit.initialization;
  out["estimand"] = fit.estimand;
  out["diagnostic_scope"] = fit.diagnostic_scope;
  out["message"] = fit.message;
  out["ambient_objective"] = fit.ambient_objective;
  out["ambient_subgradient_residual"] = fit.ambient_subgradient_residual;
  if (keep_trace) {
    const int n = fit.trace_objective.size();
    Rcpp::IntegerVector iteration(n);
    for (int i = 0; i < n; ++i) iteration[i] = i;
    out["trace"] = Rcpp::DataFrame::create(
      Rcpp::Named("iteration") = iteration,
      Rcpp::Named("objective") = fit.trace_objective,
      Rcpp::Named("stationarity") = fit.trace_stationarity,
      Rcpp::Named("step_norm") = fit.trace_step,
      Rcpp::Named("step_size") = fit.trace_alpha);
  } else {
    out["trace"] = R_NilValue;
  }
  return out;
}

} // namespace

// [[Rcpp::export]]
Rcpp::List inference_summary(std::string mfdname, Rcpp::List data,
                            arma::vec myweight, std::string statistic,
                            std::string geometry, int myiter, double myeps,
                            Rcpp::Nullable<Rcpp::NumericMatrix> initial = R_NilValue,
                            int max_backtrack = 50, bool keep_trace = false) {
  if (statistic != "mean" && statistic != "median")
    Rcpp::stop("Unknown summary statistic.");
  if (geometry != "intrinsic" && geometry != "extrinsic")
    Rcpp::stop("Unknown summary geometry backend.");
  const bool median = statistic == "median";
  const RiemannSummaryControl control(myiter, myeps, max_backtrack, keep_trace);
  const arma::field<arma::mat> observations = summary_data(data);
  arma::mat initial_matrix;
  const arma::mat* initial_pointer = NULL;
  if (initial.isNotNull()) {
    initial_matrix = Rcpp::as<arma::mat>(initial.get());
    initial_pointer = &initial_matrix;
  }
  const RiemannSummaryResult fit = geometry == "extrinsic" ?
    riem_summary_extrinsic(mfdname, observations, myweight, control, median, initial_pointer) :
    (median ? riem_summary_median(mfdname, observations, myweight, control, initial_pointer) :
              riem_summary_mean(mfdname, observations, myweight, control, initial_pointer));
  return summary_list(fit, median, keep_trace);
}

// Preserve the existing native entry points for internal callers and older R
// interfaces; all of them delegate to the same summary implementation.
// [[Rcpp::export]]
Rcpp::List inference_mean_intrinsic(std::string mfdname, Rcpp::List& data,
                                  arma::vec myweight, int myiter, double myeps) {
  return inference_summary(mfdname, data, myweight, "mean", "intrinsic", myiter, myeps);
}

// [[Rcpp::export]]
Rcpp::List inference_mean_extrinsic(std::string mfdname, Rcpp::List& data,
                                  arma::vec myweight, int myiter, double myeps) {
  return inference_summary(mfdname, data, myweight, "mean", "extrinsic", myiter, myeps);
}

// [[Rcpp::export]]
Rcpp::List inference_median_intrinsic(std::string mfdname, Rcpp::List& data,
                                    arma::vec myweight, int myiter, double myeps) {
  return inference_summary(mfdname, data, myweight, "median", "intrinsic", myiter, myeps);
}

// [[Rcpp::export]]
Rcpp::List inference_median_extrinsic(std::string mfdname, Rcpp::List& data,
                                    arma::vec myweight, int myiter, double myeps) {
  return inference_summary(mfdname, data, myweight, "median", "extrinsic", myiter, myeps);
}
