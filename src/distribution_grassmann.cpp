#include <RcppArmadillo.h>
#include "riemann_src.h"
#include "riemann_general.h"

using namespace Rcpp;
using namespace arma;
using namespace std;

// DISTRIBUTION ON GRASSMANN MANIFOLD
// (01) macg_density : (Chikuse, 1990) density evaluation for MACG
// (02) macg_sample  : sampling
// (03) macg_mle     : (Chikuse, 1990) MLE 

// (03) macg_mle ---------------------------------------------------------------
// [[Rcpp::export]]
arma::mat macg_mle(Rcpp::List& data, int maxiter, double abstol){
  // PREPARE
  int N = data.size();
  arma::field<arma::mat> mydata(N);
  for (int n=0; n<N; n++){
    mydata(n) = Rcpp::as<arma::mat>(data[n]);
  }
  int p = mydata(0).n_rows; // notation convention from Mardia & Jupp
  int r = mydata(0).n_cols; 
  
  double NN = static_cast<double>(N);
  double pp = static_cast<double>(p);
  double rr = static_cast<double>(r);
  double cc = pp/(NN*rr);
  
  // INITIALIZE
  arma::mat Sold(p,p,fill::eye);      // tr(S)=p
  arma::mat Soldinv(p,p,fill::zeros);
  arma::mat Stmp(p,p,fill::zeros);
  arma::mat Snew(p,p,fill::zeros); 
  double    Sinc = 0.0;
  
  arma::mat Xi(p,r,fill::zeros);
  
  // ITERATION
  for (int it=0; it<maxiter; it++){
    // S inverse
    Soldinv = arma::inv(Sold);
    Stmp.fill(0.0);

    // main iteration part + update
    for (int n=0; n<N; n++){
      Xi = mydata(n);
      Stmp += (Xi*arma::inv(Xi.t()*Soldinv*Xi)*Xi.t())*cc;
    }
    Snew = (Stmp/arma::trace(Stmp))*pp;
    
    // updater
    Sinc = arma::norm(Sold-Snew,"fro");
    Sold = Snew;
    if (Sinc < abstol){
      break;
    }
  }
  
  // RETURN
  return(Sold);
}

// (02) macg_sample ------------------------------------------------------------
// [[Rcpp::export]]
arma::cube macg_sample(int n, int r, arma::mat sigma){
  int p = sigma.n_rows;
  
  arma::vec mu(p,fill::zeros);
  arma::mat Xi(p,r,fill::zeros);
  arma::mat Xsq(r,r,fill::zeros);
  
  arma::cube output(p,r,n,fill::zeros);

  for (int i=0; i<n; i++){
    Xi  = arma::trans(cpp_rmvnorm(r, mu, sigma));    // (p x r)
    Xsq = arma::real(arma::powmat(Xi.t()*Xi, -0.5)); // (r x r)
    output.slice(i) = Xi*Xsq;
  }
  return(output);
}

// (01) macg_density : density evaluation for MACG -----------------------------
// [[Rcpp::export]]
arma::vec macg_density(Rcpp::List& data, arma::mat sigma){
  // PREPARE
  int N = data.size();
  arma::field<arma::mat> mydata(N);
  for (int n=0; n<N; n++){
    mydata(n) = Rcpp::as<arma::mat>(data[n]);
  }
  int p = mydata(0).n_rows; // notation convention from Mardia & Jupp
  int r = mydata(0).n_cols; 
  
  double pp = static_cast<double>(p);
  double rr = static_cast<double>(r);
  
  // COMPUTE
  arma::vec output(N,fill::zeros);
  arma::mat tmpmat(r,r,fill::zeros);
  double    tmpdet = 0.0;
  double    comval = std::pow(arma::as_scalar(arma::det(sigma)), -rr/2.0);
  arma::mat siginv = arma::inv_sympd(sigma);
  
  for (int n=0; n<N; n++){
    tmpmat = arma::trans(mydata(n))*siginv*mydata(n); // X'*Sigma^{-1}*X
    tmpdet = arma::as_scalar(arma::det(tmpmat));      // det()
    output(n) = std::pow(tmpdet, -(pp/2.0))*comval;   // det()^{-p/2}
  }
  return(output);
}
