#ifndef RIEMANN_SRC_H_
#define RIEMANN_SRC_H_

#include "RcppArmadillo.h"
#include <string>
#include <vector>

struct RiemannSummaryControl {
  int maxiter;
  double tolerance;
  int max_backtrack;
  bool keep_trace;
  RiemannSummaryControl(int iterations = 50, double eps = 1e-5,
                        int backtrack = 50, bool trace = false) :
    maxiter(iterations), tolerance(eps), max_backtrack(backtrack),
    keep_trace(trace) {}
};

struct RiemannSummaryResult {
  arma::mat estimate;
  arma::vec distances;
  double objective;
  double gradient_norm;
  double subgradient_residual;
  double step_norm;
  double ambient_objective;
  double ambient_subgradient_residual;
  bool converged;
  int iterations;
  std::string termination;
  std::string initialization;
  std::string estimand;
  std::string diagnostic_scope;
  std::string message;
  std::vector<double> trace_objective;
  std::vector<double> trace_stationarity;
  std::vector<double> trace_step;
  std::vector<double> trace_alpha;
  RiemannSummaryResult() : objective(NA_REAL), gradient_norm(NA_REAL),
    subgradient_residual(NA_REAL), step_norm(0.0), ambient_objective(NA_REAL),
    ambient_subgradient_residual(NA_REAL), converged(false), iterations(0),
    termination("maxiter"), diagnostic_scope("manifold") {}
};

RiemannSummaryResult riem_summary_mean(const std::string& mfd,
  const arma::field<arma::mat>& data, arma::vec weights,
  const RiemannSummaryControl& control, const arma::mat* initial = NULL);
RiemannSummaryResult riem_summary_median(const std::string& mfd,
  const arma::field<arma::mat>& data, arma::vec weights,
  const RiemannSummaryControl& control, const arma::mat* initial = NULL);
RiemannSummaryResult riem_summary_extrinsic(const std::string& mfd,
  const arma::field<arma::mat>& data, arma::vec weights,
  const RiemannSummaryControl& control, bool median,
  const arma::mat* initial = NULL);

// OPERATIONS ==================================================================
arma::mat riem_initialize(std::string mfd, arma::field<arma::mat> data, arma::vec weight);
arma::mat riem_initialize_cube(std::string mfd, arma::cube mydata, arma::vec weight);
arma::mat riem_exp(std::string mfd, arma::mat x, arma::mat d, double t);
arma::mat riem_log(std::string mfd, arma::mat x, arma::mat y);
double    riem_dist(std::string mfd, arma::mat x, arma::mat y);
double    riem_distext(std::string mfd, arma::mat x, arma::mat y);
arma::vec riem_equiv(std::string mfd, arma::mat x, int m, int n);
arma::mat riem_invequiv(std::string mfd, arma::vec x, int m, int n);
double    riem_metric(std::string mfd, arma::mat x, arma::mat d1, arma::mat d2);
arma::mat riem_project_tangent(const std::string& mfd, const arma::mat& x,
                              const arma::mat& direction);

// OTHER FUNCTIONS TO BE USED IN OTHER CPP MODULES =============================
arma::mat internal_mean(std::string mfd, std::string dtype, arma::cube data, int iter, double eps);
arma::mat internal_mean_init(std::string mfd, std::string dtype, arma::cube data, int iter, double eps, arma::mat Sinit);
arma::mat internal_logvectors(std::string mfd, arma::cube data); // row-stacked vectors
arma::uvec helper_sample(int N, int m, arma::vec prob, bool replace);
arma::uvec helper_setdiff(arma::uvec& x, arma::uvec& y);


#endif
