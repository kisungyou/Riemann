#include <RcppArmadillo.h>
#include "riemann_src.h"
#include <map>

using namespace Rcpp;
using namespace arma;
using namespace std;

// =============================================================================
// WASSERSTEIN GEOMETRY
// =============================================================================
// (01) spdwass_sylvester : solve the sylvester equation L_A[X]; A-SPD, X-Symm
// (02) spdwass_log       : logarithmic map
// (03) spdwass_exp       : exponential map
// (04) spdwass_metric    : riemannian metric
// (05) spdwass_baryRU02  : wasserstein barycenter by Ruschendorf & Uckelman
// (06) spdwass_baryAE16  :                        by Alvarez-Esteban

// =============================================================================
// SPECIAL FUNCTIONS ON SPD MANIFOLD
// =============================================================================
// (01) src_spd_dist    : compute distance of two SPD matrices : ADD ON HERE!
//      src_spd_pdist   : pairwise distances



// =============================================================================
// WASSERSTEIN GEOMETRY
// =============================================================================
// (01) spdwass_sylvester
// [[Rcpp::export]]
arma::mat spdwass_sylvester(arma::mat A, arma::mat X){
  // eigen-decomposition
  arma::vec Lambda;
  arma::mat Q;
  arma::eig_sym(Lambda, Q, A);
  
  int N = A.n_rows;
  arma::mat C = Q.t()*X*Q;
  arma::mat E(N,N,fill::zeros);
  for (int n=0; n<N; n++){
    E(n,n) = C(n,n)/(2.0*Lambda(n));
  }
  for (int i=0; i<(N-1); i++){
    for (int j=(i+1); j<N; j++){
      E(i,j) = C(i,j)/(Lambda(i)+Lambda(j));
      E(j,i) = E(i,j);
    }
  }
  arma::mat output=Q*E*Q.t();
  return(output);
}
// (02) spdwass_log
// [[Rcpp::export]]
arma::mat spdwass_log(arma::mat C, arma::mat X){
  return(arma::sqrtmat_sympd(C*X) + arma::sqrtmat_sympd(X*C) - (2.0*C));
}
// (03) spdwass_exp
// [[Rcpp::export]]
arma::mat spdwass_exp(arma::mat C, arma::mat V, double t=1.0){
  arma::mat LCV = spdwass_sylvester(C,V);
  arma::mat output = C + (t*V) + (t*t)*(LCV*C*LCV);
  return(output);
}
// (04) spdwass_metric
// [[Rcpp::export]]
double spdwass_metric(arma::mat S, arma::mat X, arma::mat Y){
  arma::mat LSX = spdwass_sylvester(S, X);
  return(arma::trace(LSX.t()*Y)/2.0);
}
// (05) spdwass_baryRU02
// [[Rcpp::export]]
arma::mat spdwass_baryRU02(arma::field<arma::mat> spdlist, arma::vec weight, int maxiter, double abstol){
  // preparation
  int N = spdlist.n_elem;
  int p = spdlist(0).n_rows;
  arma::vec proportion = weight/arma::accu(weight);
  
  
  arma::cube Stower(p,p,N,fill::zeros);
  for (int n=0; n<N; n++){
    Stower.slice(n) = arma::logmat_sympd(spdlist(n));
  }
  arma::mat Stmean = arma::mean(Stower, 2);
  arma::mat Sold   = arma::expmat_sym(Stmean);
  
  arma::mat Snew(p,p,fill::zeros);
  arma::mat Shalf(p,p,fill::zeros);
  
  double Sinc = 1000.0;
  
  // iteration
  for (int it=0; it<maxiter; it++){
    // compute the half
    Shalf = arma::sqrtmat_sympd(Sold);
    // compute a target
    Snew.fill(0.0);
    for (int n=0; n<N; n++){
      Snew += proportion(n)*arma::sqrtmat_sympd(Shalf*spdlist(n)*Shalf);
    }
    // updater
    Sinc = arma::norm(Sold-Snew,"fro");
    Sold = Snew;
    if (Sinc < abstol){
      break;
    }
  }
  return(Sold);
}
// (06) spdwass_baryAE16
// [[Rcpp::export]]
arma::mat spdwass_baryAE16(arma::field<arma::mat> spdlist, arma::vec weight, int maxiter, double abstol){
  // preparation
  int N = spdlist.n_elem;
  int p = spdlist(0).n_rows;
  arma::vec proportion = weight/arma::accu(weight);
  
  arma::cube Stower(p,p,N,fill::zeros);
  for (int n=0; n<N; n++){
    Stower.slice(n) = arma::logmat_sympd(spdlist(n));
  }
  arma::mat Stmean = arma::mean(Stower, 2);
  arma::mat Sold   = arma::expmat_sym(Stmean);
  
  arma::mat Stmp(p,p,fill::zeros);
  arma::mat Snew(p,p,fill::zeros);
  arma::mat Shalf(p,p,fill::zeros);
  arma::mat Shinv(p,p,fill::zeros);
  
  double Sinc = 1000.0;
  
  // iteration
  for (int it=0; it<maxiter; it++){
    // compute the half
    Shalf = arma::sqrtmat_sympd(Sold);
    Shinv = arma::inv_sympd(Shalf);
    
    // compute a target
    Stmp.fill(0.0);
    for (int n=0; n<N; n++){
      Stmp += proportion(n)*arma::sqrtmat_sympd(Shalf*spdlist(n)*Shalf);
    }
    Snew = Shinv*Stmp*Stmp*Shinv;
    
    // updater
    Sinc = arma::norm(Sold-Snew,"fro");
    Sold = Snew;
    if (Sinc < abstol){
      break;
    }
  }
  return(Sold);
}


// =============================================================================
// SPECIAL FUNCTIONS ON SPD MANIFOLD
// =============================================================================
// (01) spd_dist  : compute distance of two SPD matrices -----------------------
double src_spd_dist(arma::mat X, arma::mat Y, std::string geometry){
  double output = 0.0;
  if (geometry=="airm"){                                              // 1. AIRM
    output = riem_dist("spd",X,Y);
  } else if (geometry=="lerm"){                                       // 2. LERM
    output = riem_distext("spd",X,Y);
  } else if (geometry=="jeffrey"){                                    // 3. Jeffrey
    double term1 = arma::trace(arma::solve(X,Y))/2.0;
    double term2 = arma::trace(arma::solve(Y,X))/2.0;
    double term3 = static_cast<double>(X.n_rows);
    output = term1 + term2 - term3;
  } else if (geometry=="stein"){                                      // 4. Stein
    output = std::sqrt(std::log(arma::det((X+Y)/2.0)) - 0.5*std::log(arma::det(X*Y)));
  } else if (geometry=="wasserstein"){                                // 5. Wasserstein
    arma::mat Xsqrt = arma::sqrtmat_sympd(X);
    output = std::sqrt(arma::trace(X+Y-2.0*arma::sqrtmat_sympd((Xsqrt*Y*Xsqrt))));
  }
  return(output);
}
// [[Rcpp::export]]
arma::mat src_spd_pdist(arma::cube &data, std::string geometry){
  // PRELIMINARY
  int N = data.n_slices;
  
  //arma::mat exmat  = Rcpp::as<arma::mat>(data[0]);
  //int p = exmat.n_rows;
  //int N = data.size();
  //double pp = static_cast<double>(p);
  
  // COMPUTE
  arma::mat distance(N,N,fill::zeros);
  for (int i=0; i<(N-1); i++){
    for (int j=(i+1); j<N; j++){
      distance(i,j) = src_spd_dist(data.slice(i), data.slice(j), geometry);
      distance(j,i) = distance(i,j);
    }
  }
  
  // RETURN
  return(distance);
}
