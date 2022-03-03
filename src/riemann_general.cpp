#include <RcppArmadillo.h>
#include "riemann_general.h"

using namespace Rcpp;
using namespace arma;
using namespace std;

// [[Rcpp::export]]
arma::uword mat_rank(arma::mat A){
  arma::uword r = arma::rank(A);
  return(r);  
}
// [[Rcpp::export]]
arma::mat mat_symm(arma::mat A, bool diag){
  arma::mat B = (A+A.t())/2.0;
  if (diag==false){
    return(B);
  } else {
    arma::mat C = arma::diagmat(B);
    return(C);
  }
}
// [[Rcpp::export]]
arma::mat mat_diaghalf(arma::mat D){
  arma::vec dsqrt = arma::sqrt(arma::diagvec(D));
  return(arma::diagmat(dsqrt));
}
// [[Rcpp::export]]
arma::mat mat_diaginvhalf(arma::mat D){
  arma::vec dsqrt = arma::sqrt(arma::diagvec(D));
  return(arma::diagmat(1.0/dsqrt));
}
// [[Rcpp::export]]
arma::mat mat_cov2cor(arma::mat A){
  arma::vec d = arma::sqrt(arma::diagvec(A));
  arma::mat D = arma::diagmat(1.0/d);
  return(D*A*D);
}


// [[Rcpp::export]]
arma::mat cpp_rmvnorm(int n, arma::vec mu, arma::mat sigma){
  int p = sigma.n_rows;
  arma::mat Y(n,p,fill::randn);
  arma::mat output = (Y*arma::chol(sigma)) + arma::repmat(mu, 1, n).t();
  return(output);
}

