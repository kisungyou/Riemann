#include <RcppArmadillo.h>
#include "riemann_src.h"

using namespace Rcpp;
using namespace arma;
using namespace std;

/* CURVES
 * (01) curvedist_lp
 * (02) curvedist_dtwbasic : Wikipedia's work [dtw::symmetric type 1]
 */

// (01). curvedist_lp ==========================================================
// [[Rcpp::export]]
double curvedist_lp(std::string mfd, std::string geo, Rcpp::List& data1, Rcpp::List& data2, arma::vec vect, double myp){
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

// (02) curvedist_dtwbasic : Wikipedia's work ==================================
// [[Rcpp::export]]
double curvedist_dtwbasic(std::string mfd, std::string geo, Rcpp::List& data1, Rcpp::List& data2){
  // PARAMETER AND DATA PREP
  int N = data1.size();
  int M = data2.size();
  arma::mat exemplar = Rcpp::as<arma::mat>(data1[0]);
  int nrow = exemplar.n_rows;
  int ncol = exemplar.n_cols;  
  
  arma::cube mydata1(nrow,ncol,N,fill::zeros);
  arma::cube mydata2(nrow,ncol,M,fill::zeros);
  for (int n=0; n<N; n++){
    mydata1.slice(n) = Rcpp::as<arma::mat>(data1[n]);
  }
  for (int m=0; m<M; m++){
    mydata2.slice(m) = Rcpp::as<arma::mat>(data2[m]);
  }
  
  arma::mat DTW(N+1, M+1, fill::none); DTW.fill(arma::datum::inf);
  DTW(0,0) = 0.0;
  
  double cost = 0.0;
  arma::vec ctriple(3,fill::zeros);
  for (int i=1; i<=N; i++){
    for (int j=1; j<=M; j++){
      if (geo=="intrinsic"){
        cost = riem_dist(mfd, mydata1.slice(i-1), mydata2.slice(j-1));
      } else {
        cost = riem_distext(mfd, mydata1.slice(i-1), mydata2.slice(j-1));
      }
      ctriple(0) = DTW(i-1, j);  // insertion
      ctriple(1) = DTW(i, j-1);  // deletion
      ctriple(2) = DTW(i-1,j-1); // match
      DTW(i,j)   = cost + ctriple.min();
    }
  }

  double output = DTW(N, M);
  return(output);
}

