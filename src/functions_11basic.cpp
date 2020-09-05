#include <RcppArmadillo.h>
#include "riemann_src.h"

using namespace Rcpp;
using namespace arma;
using namespace std;

// 1. basic_pdist
// 2. basic_pdist2
// 3. basic_interpolate


// 1. basic_pdist ==============================================================
// [[Rcpp::export]]
arma::mat basic_pdist(std::string mfdname, Rcpp::List& data, std::string dtype){
  // PREPARE
  int N = data.size();
  arma::field<arma::mat> mydata(N);
  for (int n=0; n<N; n++){
    mydata(n) = Rcpp::as<arma::mat>(data[n]);
  }
  
  // COMPUTE
  arma::mat output(N,N,fill::zeros);
  for (int i=0; i<(N-1); i++){
    for (int j=(i+1); j<N; j++){
      if (dtype=="intrinsic"){
        output(i,j) = riem_dist(mfdname, mydata(i), mydata(j));
      } else {
        output(i,j) = riem_distext(mfdname, mydata(i), mydata(j));
      }
      output(j,i) = output(i,j);
    }
  }
  return(output);
}

// 2. basic_pdist2 =============================================================
// [[Rcpp::export]]
arma::mat basic_pdist2(std::string mfdname, Rcpp::List& data1, Rcpp::List& data2, std::string dtype){
  // PREPARE
  int M = data1.size();
  int N = data2.size();
  
  // COMPUTE
  arma::mat mat1;
  arma::mat mat2;
  arma::mat output(M,N,fill::zeros);
  for (int m=0; m<M; m++){
    mat1 = Rcpp::as<arma::mat>(data1[m]);
    for (int n=0; n<N; n++){
      mat2 = Rcpp::as<arma::mat>(data2[n]);
      if (dtype=="intrinsic"){
        output(m,n) = riem_dist(mfdname, mat1, mat2);
      } else {
        output(m,n) = riem_distext(mfdname, mat1, mat2);
      }
    }
  }
  return(output);
}

// 3. basic_interpolate ========================================================
// [[Rcpp::export]]
arma::cube basic_interpolate(std::string mfdname, std::string dtype, arma::mat mat1, arma::mat mat2, arma::vec vect){
  // PREPARE
  int p = mat1.n_rows;
  int k = mat1.n_cols;
  int tt = vect.n_elem;
  arma::cube output(p,k,tt,fill::zeros);
  
  if (dtype=="intrinsic"){
    arma::mat logxy = riem_log(mfdname, mat1, mat2);
    for (int i=0; i<tt; i++){
      output.slice(i) = riem_exp(mfdname, mat1, logxy, vect(i));
    }
  } else if (dtype=="extrinsic"){
    arma::vec vec1 = riem_equiv(mfdname, mat1, p, k);
    arma::vec vec2 = riem_equiv(mfdname, mat2, p, k);
    arma::vec veci(vec1.n_elem, fill::zeros);
    for (int i=0; i<tt; i++){
      veci = (1.0-vect(i))*vec1 + vect(i)*vec2;
      output.slice(i) = riem_invequiv(mfdname, veci, p, k);
    }
  }
  return(output);
}