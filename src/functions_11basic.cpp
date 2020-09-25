#include <RcppArmadillo.h>
#include "riemann_src.h"

using namespace Rcpp;
using namespace arma;
using namespace std;

// 1. basic_pdist
// 2. basic_pdist2
// 3. basic_interpolate
// 4. basic_curvedist_lp


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

// 4. basic_curvedist_lp =======================================================
// [[Rcpp::export]]
double basic_curvedist_lp(std::string mfd, std::string geo, Rcpp::List& data1, Rcpp::List& data2, arma::vec vect, double myp){
  // PARAMETER AND DATA PREP
  int N = vect.n_elem;
  arma::mat exemplar = Rcpp::as<arma::mat>(data1[0]);
  int nrow = exemplar.n_rows;
  int ncol = exemplar.n_cols;  
  
  arma::cube mydata1(nrow,ncol,N,fill::zeros);
  arma::cube mydata2(nrow,ncol,N,fill::zeros);
  for (int n=0; n<N; n++){
    mydata1.slice(n) = Rcpp::as<arma::mat>(data1[n]);
    mydata2.slice(n) = Rcpp::as<arma::mat>(data2[n]);
  }
  
  // COMPUTE DISTANCE ACROSS POINTS P-TH POWER
  double    dval = 0.0;
  arma::vec distp(N,fill::zeros);
  for (int n=0; n<N; n++){
    if (geo=="intrinsic"){
      dval = riem_dist(mfd, mydata1.slice(n), mydata2.slice(n));
    } else {
      dval = riem_distext(mfd, mydata1.slice(n), mydata2.slice(n));
    }
    distp(n) = std::pow(dval, myp);
  }
  
  // AGGREGATE
  double output = 0.0;
  for (int n=0; n<(N-1); n++){
    output += (distp(n) + distp(n+1))*std::abs(vect(n+1)-vect(n))/2.0;
  }
  double outval = std::pow(output, (1.0/myp));
  return(outval);
}