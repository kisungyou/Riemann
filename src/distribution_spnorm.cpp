/* Everything from RiemSphere
 * (01) cppdist_int_1toN : intrinsic distance from a vector to row-vecs
 * (02) cppdist_ext_1toN : extrinsic distance
 */

#include <RcppArmadillo.h>
#include "riemann_src.h"
#include <math.h>

using namespace Rcpp;
using namespace arma;
using namespace std;

// (01) cppdist_int_1toN : intrinsic distance from a vector to row-vecs ========
//[[Rcpp::export]]
arma::vec cppdist_int_1toN(arma::vec x, arma::mat &Y){
  // parameters
  int N = Y.n_rows;
  int p = x.n_elem;
  
  // preliminary
  arma::vec output(N,fill::zeros);
  arma::rowvec xvec = x.t();
  arma::rowvec yvec(p,fill::zeros);
  
  // iteration
  for (int n=0;n<N;n++){
    yvec = Y.row(n);
    if (arma::norm(xvec-yvec,2) > 1e-10){
      output(n) = std::acos(arma::dot(xvec, yvec));  
    } else {
      output(n) = 0;
    }
  }
  
  // return
  return(output);
}
// (02) cppdist_ext_1toN : extrinsic distance ==================================
//[[Rcpp::export]]
arma::vec cppdist_ext_1toN(arma::vec x, arma::mat &Y){
  // parameters
  int N = Y.n_rows;
  int p = x.n_elem;
  
  // preliminary
  arma::vec output(N,fill::zeros);
  arma::rowvec xvec = x.t();
  arma::rowvec yvec(p,fill::zeros);
  
  // iteration
  for (int n=0;n<N;n++){
    yvec      = Y.row(n);
    output(n) = arma::norm(xvec-yvec, 2);
  }
  
  // return
  return(output);
}

