#include <RcppArmadillo.h>
#include "riemann_src.h"

// Internal bridge shared by geometry validation and native-path regression tests.
// [[Rcpp::export]]
Rcpp::List geometry_operations(std::string manifold, arma::mat x,
                               arma::mat y, arma::mat u, arma::mat v) {
  return Rcpp::List::create(
    Rcpp::Named("distance") = riem_dist(manifold,x,y),
    Rcpp::Named("log") = riem_log(manifold,x,y),
    Rcpp::Named("exp") = riem_exp(manifold,x,u,1.0),
    Rcpp::Named("metric") = riem_metric(manifold,x,u,v));
}
