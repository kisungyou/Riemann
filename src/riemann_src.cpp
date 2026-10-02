#include <RcppArmadillo.h>
#include "riemann_src.h"
#include "riemann_manifolds.h"
#include <RcppArmadilloExtensions/sample.h>

using namespace Rcpp;
using namespace arma;
using namespace std;

// OPERATIONS ==================================================================
// (01) riem_initialize & riem_initialize_cube
// (02) riem_exp
// (03) riem_log
// (04) riem_dist
// (05) riem_distext
// (06) riem_equiv
// (07) riem_invequiv
// (08) riem_metric

// INTERNAL FUNCTIONS TO BE USED ===============================================
// (01) internal_mean

// (01) riem_initialize ========================================================
arma::mat riem_initialize(std::string mfd, arma::field<arma::mat> data, arma::vec weight){
  arma::mat output;
  if (mfd=="sphere"){
    output = sphere_initialize(data, weight);
  } else if (mfd=="landmark"){
    output = landmark_initialize(data, weight);
  } else if (mfd=="spdk"){
    output = spdk_initialize(data, weight);
  } else if (mfd=="multinomial"){
    output = multinomial_initialize(data, weight);
  } else if (mfd=="stiefel"){
    output = stiefel_initialize(data, weight);
  } else if (mfd=="spd"){
    output = spd_initialize(data, weight);
  } else if (mfd=="correlation"){
    output = correlation_initialize(data, weight);
  } else if (mfd=="grassmann"){
    output = grassmann_initialize(data, weight);
  } else if (mfd=="euclidean"){
    output = euclidean_initialize(data, weight);
  } else if (mfd=="rotation"){
    output = rotation_initialize(data, weight);
  } else {
    std::string err = "* Riemann : 'initialization' is not implemented for " + mfd + " manifold.";
    Rcpp::stop(err);
  }
  return(output);
}
arma::mat riem_initialize_cube(std::string mfd, arma::cube mydata, arma::vec weight){
  int N = mydata.n_slices;
  arma::field<arma::mat> data(N);
  for (int n=0; n<N; n++){
    data(n) = mydata.slice(n);
  }
  return(riem_initialize(mfd, data, weight));
}
// (02) riem_exp ===============================================================
arma::mat riem_exp(std::string mfd, arma::mat x, arma::mat d, double t){
  arma::mat output;
  if (mfd=="sphere"){
    output = sphere_exp(x, d, t);
  } else if (mfd=="landmark"){
    output = landmark_exp(x, d, t);
  } else if (mfd=="spdk"){
    output = spdk_exp(x,d,t);
  } else if (mfd=="multinomial"){
    output = multinomial_exp(x, d, t);
  } else if (mfd=="stiefel"){
    output = stiefel_exp(x, d, t);
  } else if (mfd=="grassmann"){
    output = grassmann_exp(x, d, t);
  } else if (mfd=="rotation"){
    output = rotation_exp(x, d, t);
  } else if (mfd=="spd"){
    output = spd_exp(x, d, t);
  } else if (mfd=="euclidean"){
    output = euclidean_exp(x, d, t);
  } else if (mfd=="correlation"){
    output = correlation_exp(x, d, t);
  } else {
    std::string err = "* Riemann : 'exponential map' is not implemented for " + mfd + " manifold.";
    Rcpp::stop(err);
  }
  return(output);
}


// (03) riem_log ===============================================================
arma::mat riem_log(std::string mfd, arma::mat x, arma::mat y){
  arma::mat output;
  if (mfd=="sphere"){
    output = sphere_log(x, y);
  } else if (mfd=="landmark"){
    output = landmark_log(x, y);
  } else if (mfd=="spdk"){
    output = spdk_log(x,y);
  } else if (mfd=="multinomial"){
    output = multinomial_log(x, y);
  } else if (mfd=="grassmann"){
    output = grassmann_log(x, y); 
  } else if (mfd=="rotation"){
    output = rotation_log(x, y);
  } else if (mfd=="stiefel"){
    output = stiefel_log(x, y);
  } else if (mfd=="spd"){
    output = spd_log(x, y);
  } else if (mfd=="euclidean"){
    output = euclidean_log(x, y);
  } else if (mfd=="correlation"){
    output = correlation_log(x, y);
  } else {
    std::string err = "* Riemann : 'logarithm map' is not implemented for " + mfd + " manifold.";
    Rcpp::stop(err);
  }
  return(output);
}


// (04) riem_dist ==============================================================
double riem_dist(std::string mfd, arma::mat x, arma::mat y){
  double output;
  if (mfd=="sphere"){
    output = sphere_dist(x, y);
  } else if (mfd=="landmark"){
    output = landmark_dist(x, y);
  } else if (mfd=="spdk"){
    output = spdk_dist(x, y);
  } else if (mfd=="multinomial"){
    output = multinomial_dist(x, y);
  } else if (mfd=="grassmann"){
    output = grassmann_dist(x, y); 
  } else if (mfd=="stiefel"){
    output = stiefel_dist(x, y); 
  } else if (mfd=="rotation"){
    output = rotation_dist(x, y);
  } else if (mfd=="spd"){
    output = spd_dist(x, y);
  } else if (mfd=="euclidean"){
    output = euclidean_dist(x, y);
  } else if (mfd=="correlation"){
    output = correlation_dist(x, y);
  } else {
    std::string err = "* Riemann : 'geodesic distance' is not implemented for " + mfd + " manifold.";
    Rcpp::stop(err);
  }
  return(output); 
}

// (05) riem_distext ===========================================================
double riem_distext(std::string mfd, arma::mat x, arma::mat y){
  double output;
  if (mfd=="sphere"){
    output = sphere_distext(x, y);
  } else if (mfd=="landmark"){
    output = landmark_distext(x, y);
  } else if (mfd=="multinomial"){
    output = multinomial_distext(x, y);
  } else if (mfd=="grassmann"){
    output = grassmann_distext(x, y);
  } else if (mfd=="stiefel"){
    output = stiefel_distext(x, y);
  } else if (mfd=="rotation"){
    output = rotation_distext(x, y);
  } else if (mfd=="spd"){
    output = spd_distext(x, y);
  } else if (mfd=="euclidean"){
    output = euclidean_distext(x, y);
  } else {
    std::string err = "* Riemann : 'extrinsic distance' is not implemented for " + mfd + " manifold.";
    Rcpp::stop(err);
  }
  return(output); 
}

// (06) riem_equiv =============================================================
arma::vec riem_equiv(std::string mfd, arma::mat x, int m, int n){
  arma::vec output;
  if (mfd=="sphere"){
    output = sphere_equiv(x, m, n);
  } else if (mfd=="landmark"){
    output = landmark_equiv(x, m, n);
  } else if (mfd=="multinomial"){
    output = multinomial_equiv(x, m, n);
  } else if (mfd=="grassmann"){
    output = grassmann_equiv(x, m, n);
  } else if (mfd=="stiefel"){
    output = stiefel_equiv(x, m, n);
  } else if (mfd=="spd"){
    output = spd_equiv(x, m, n);
  } else if (mfd=="euclidean"){
    output = euclidean_equiv(x, m, n);
  } else if (mfd=="rotation"){
    output = rotation_equiv(x, m, n);
  } else {
    std::string err = "* Riemann : 'equivariant embedding' is not implemented for " + mfd + " manifold.";
    Rcpp::stop(err);
  }
  return(output); 
}
// (07) riem_invequiv ==========================================================
arma::mat riem_invequiv(std::string mfd, arma::vec x, int m, int n){
  arma::mat output;
  if (mfd=="sphere"){
    output = sphere_invequiv(x, m, n);
  } else if (mfd=="landmark"){
    output = landmark_invequiv(x, m, n);
  } else if (mfd=="multinomial"){
    output = multinomial_invequiv(x, m, n);
  } else if (mfd=="grassmann"){
    output = grassmann_invequiv(x, m, n);
  } else if (mfd=="stiefel"){
    output = stiefel_invequiv(x, m, n);
  } else if (mfd=="spd"){
    output = spd_invequiv(x, m, n);
  } else if (mfd=="euclidean"){
    output = euclidean_invequiv(x, m, n);
  } else if (mfd=="rotation"){
    output = rotation_invequiv(x, m, n);
  } else {
    std::string err = "* Riemann : 'inverse equivariant embedding' is not implemented for " + mfd + " manifold.";
    Rcpp::stop(err);
  }
  return(output);  
}
// (08) riem_metric ============================================================
double riem_metric(std::string mfd, arma::mat x, arma::mat d1, arma::mat d2){
  double output;
  if (mfd=="sphere"){
    output = sphere_metric(x,d1,d2);
  } else if (mfd=="landmark"){
    output = landmark_metric(x,d1,d2);
  } else if (mfd=="spdk"){
    output = spdk_metric(x,d1,d2);
  } else if (mfd=="grassmann"){
    output = grassmann_metric(x,d1,d2);
  } else if (mfd=="multinomial"){
    output = multinomial_metric(x,d1,d2);
  } else if (mfd=="stiefel"){
    output = stiefel_metric(x,d1,d2);
  } else if (mfd=="rotation"){
    output = rotation_metric(x,d1,d2);
  } else if (mfd=="spd"){
    output = spd_metric(x,d1,d2);
  } else if (mfd=="euclidean"){
    output = euclidean_metric(x,d1,d2);
  } else if (mfd=="correlation"){
    output = correlation_metric(x,d1,d2);
  } else {
    std::string err = "* Riemann : 'Riemannian metric' is not implemented for " + mfd + " manifold.";
    Rcpp::stop(err);
  }
  return(output);  
}



// OTHER FUNCTIONS TO BE USED IN OTHER CPP MODULES =============================
arma::mat riem_project_tangent(const std::string& mfd, const arma::mat& x,
                              const arma::mat& direction) {
  if (mfd == "sphere") return sphere_proj(x, direction);
  if (mfd == "spd") return direction / 2.0 + direction.t() / 2.0;
  if (mfd == "grassmann") return grassmann_proj(x, direction);
  if (mfd == "rotation") return direction / 2.0 - direction.t() / 2.0;
  if (mfd == "multinomial") return direction - arma::accu(direction) * x;
  if (mfd == "spdk") return legacy_horizontal_projection(x, direction);
  if (mfd == "landmark") {
    arma::mat out = direction;
    out.each_row() -= arma::mean(out, 0);
    out -= arma::accu(x % out) * x;
    return legacy_horizontal_projection(x, out);
  }
  return direction;
}

namespace {
arma::mat checked_internal_mean(std::string mfd, std::string dtype,
                                const arma::cube& data, int iter, double eps,
                                const arma::mat* initial) {
  if (data.n_slices == 0) Rcpp::stop("An internal mean requires at least one observation.");
  arma::field<arma::mat> observations(data.n_slices);
  for (arma::uword i = 0; i < data.n_slices; ++i) observations(i) = data.slice(i);
  arma::vec weights(data.n_slices, arma::fill::ones);
  RiemannSummaryControl control(iter, eps);
  RiemannSummaryResult fit;
  if (dtype == "intrinsic") {
    fit = riem_summary_mean(mfd, observations, weights, control, initial);
  } else if (dtype == "extrinsic") {
    fit = riem_summary_extrinsic(mfd, observations, weights, control, false, initial);
  } else {
    Rcpp::stop("Unknown internal mean geometry backend.");
  }
  if (!fit.converged)
    Rcpp::stop("An internal mean did not converge (%s).", fit.termination.c_str());
  return fit.estimate;
}
} // namespace

arma::mat internal_mean(std::string mfd, std::string dtype, arma::cube data, int iter, double eps) {
  return checked_internal_mean(mfd, dtype, data, iter, eps, NULL);
}

arma::mat internal_mean_init(std::string mfd, std::string dtype, arma::cube data,
                             int iter, double eps, arma::mat Sinit) {
  return checked_internal_mean(mfd, dtype, data, iter, eps, &Sinit);
}

arma::mat internal_logvectors(std::string mfd, arma::cube data){
  // PARAMETERS
  int nrow = data.n_rows;
  int ncol = data.n_cols;
  int N    = data.n_slices;
  
  // COMPUTE MEAN
  arma::mat X = internal_mean(mfd, "intrinsic", data, 50, 1e-6);
  
  // EXEMPLARY VECTOR
  arma::vec logvec = arma::vectorise(riem_log(mfd, X, data.slice(0)));
  int P = logvec.n_elem;
  
  // ITERATION
  arma::mat output(N,P,fill::zeros);
  for (int n=0; n<N; n++){
    output.row(n) = arma::trans(arma::vectorise(riem_log(mfd, X, data.slice(n))));
  }
  return(output);
}
arma::uvec helper_sample(int N, int m, arma::vec prob, bool replace){
  arma::uvec x     = arma::linspace<arma::uvec>(0L, N-1L, N);
  arma::vec myprob = prob/arma::accu(prob);
  arma::uvec output = Rcpp::RcppArmadillo::sample(x, m, replace, myprob);
  return(output);
}
arma::uvec helper_setdiff(arma::uvec& x, arma::uvec& y){
  std::vector<int> a = arma::conv_to< std::vector<int> >::from(arma::sort(x));
  std::vector<int> b = arma::conv_to< std::vector<int> >::from(arma::sort(y));
  std::vector<int> out;
  
  std::set_difference(a.begin(), a.end(), b.begin(), b.end(),
                      std::inserter(out, out.end()));
  
  return arma::conv_to<arma::uvec>::from(out);
}
