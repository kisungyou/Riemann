#include <RcppArmadillo.h>
#include "riemann_src.h"

using namespace Rcpp;
using namespace arma;
using namespace std;

// 1. basic_pdist
// 2. basic_pdist2


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