#ifndef RIEMANN_GENERAL_H_
#define RIEMANN_GENERAL_H_

#include "RcppArmadillo.h"

arma::uword mat_rank(arma::mat A);            // compute matrix rank
arma::mat   mat_symm(arma::mat A, bool diag); // matrix symmetrization
arma::mat   mat_diaghalf(arma::mat D);        // compute D^{1/2}
arma::mat   mat_diaginvhalf(arma::mat D);     // compute D^{-1/2}
arma::mat   mat_cov2cor(arma::mat A);         // convert cov -> cor

arma::mat cpp_rmvnorm(int n, arma::vec mu, arma::mat sigma); // rmvnorm


#endif

