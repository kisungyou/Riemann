#ifndef RIEMANN_SRC_H_
#define RIEMANN_SRC_H_

#include "RcppArmadillo.h"

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

// OTHER FUNCTIONS TO BE USED IN OTHER CPP MODULES =============================
arma::mat internal_mean(std::string mfd, std::string dtype, arma::cube data, int iter, double eps);
arma::mat internal_mean_init(std::string mfd, std::string dtype, arma::cube data, int iter, double eps, arma::mat Sinit);
arma::mat internal_logvectors(std::string mfd, arma::cube data); // row-stacked vectors
arma::uvec helper_sample(int N, int m, arma::vec prob, bool replace);
arma::uvec helper_setdiff(arma::uvec& x, arma::uvec& y);


#endif