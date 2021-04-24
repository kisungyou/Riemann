#include <RcppArmadillo.h>
#include "riemann_src.h"

using namespace Rcpp;
using namespace arma;
using namespace std;

// DISTRIBUTION ON SPHERE
// (01) acg_density : (Tyler, 1986-7)
// (02) acg_mle     : MLE


// (02) acg_mle ----------------------------------------------------------------
// [[Rcpp::export]]
arma::mat acg_mle(Rcpp::List& data, int maxiter, double abstol){
  // PREPARE
  int N = data.size();
  arma::field<arma::mat> mydata(N);
  for (int n=0; n<N; n++){
    mydata(n) = Rcpp::as<arma::mat>(data[n]);
  }
  int p = mydata(0).n_rows;
  double pp = static_cast<double>(p);
  
  // INITIALIZE
  arma::mat Aold(p,p,fill::eye);      // tr(A)=p
  arma::mat Aoldinv(p,p,fill::zeros);
  arma::mat Anew(p,p,fill::zeros); 
  double    Ainc = 0.0;
  
  arma::vec form2(N,fill::zeros);
  double    form2invsum = 0.0;

  // FIXED-POINT ITERATION
  for (int it=0; it<maxiter; it++){
    // A inverse
    Aoldinv = arma::inv(Aold);
    Anew.fill(0.0);
    
    // main computation
    for (int n=0; n<N; n++){
      form2(n) = arma::as_scalar(arma::trans(mydata(n))*Aoldinv*mydata(n));
    }
    form2invsum = 0.0;
    for (int n=0; n<N; n++){
      form2invsum += 1.0/form2(n);
    }
    for (int n=0; n<N; n++){
      Anew += (mydata(n)*arma::trans(mydata(n)))/form2(n);
    }
    Anew *= (pp/form2invsum);
    
    // updater
    Ainc = arma::norm(Aold-Anew,"fro");
    Aold = Anew;
    
    if (Ainc < abstol){
      break;
    }
  }
  return(Aold);
}

// (01) acg_density ------------------------------------------------------------
// [[Rcpp::export]]
arma::vec acg_density(Rcpp::List& data, arma::mat A){
  // PREPARE
  int N = data.size();
  arma::field<arma::mat> mydata(N);
  for (int n=0; n<N; n++){
    mydata(n) = Rcpp::as<arma::mat>(data[n]);
  }
  int p = A.n_rows;
  double pp = static_cast<double>(p);
  
  // PRELIMINARY COMPUTATION
  arma::mat Ainv  = arma::inv_sympd(A);
  double    Acoef = 1.0/std::sqrt(arma::det(A));
  
  // ITERATION
  arma::vec output(N,fill::zeros);
  for (int n=0; n<N; n++){
     output(n) = std::pow(arma::as_scalar(arma::trans(mydata(n))*Ainv*mydata(n)), -pp/2.0)*Acoef;
  }
  return(output);
}