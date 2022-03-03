#include <RcppArmadillo.h>
#include "riemann_general.h"

using namespace Rcpp;
using namespace arma;
using namespace std;

// SPECIAL FUNCTIONS ===========================================================
// runif_sphere
// runif_stiefel


// runif_sphere ----------------------------------------------------------------
// [[Rcpp::export]]
arma::mat runif_sphere(int n, int p){
  arma::mat output(n,p,fill::randn);
  arma::rowvec outvec(p,fill::zeros);
  for (int i=0; i<n; i++){
    outvec = output.row(i);
    output.row(i) = outvec/arma::norm(outvec,2);
  }
  return(output);
}

// runif_stiefel ---------------------------------------------------------------
// [[Rcpp::export]]
arma::cube runif_stiefel(int p, int k, int N){
  arma::cube output(p,k,N,fill::randn);
  arma::mat X(p,k,fill::zeros);
  arma::mat H(k,k,fill::zeros);
  for (int n=0; n<N; n++){
    X = output.slice(n);
    H = arma::real(arma::powmat(X.t()*X, -0.5));
    output.slice(n) = X*H;
  }
  return(output);
}

