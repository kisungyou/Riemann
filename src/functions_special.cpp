#include <RcppArmadillo.h>
#include "riemann_general.h"

using namespace Rcpp;
using namespace arma;
using namespace std;

// SPECIAL FUNCTIONS ===========================================================
// 01. sphere_runif


// 01. sphere_runif ------------------------------------------------------------
// [[Rcpp::export]]
arma::mat sphere_runif(int n, int p){
  arma::mat output(n,p,fill::randn);
  arma::rowvec outvec(p,fill::zeros);
  for (int i=0; i<n; i++){
    outvec = output.row(i);
    output.row(i) = outvec/arma::norm(outvec,2);
  }
  return(output);
}
  