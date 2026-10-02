#include <RcppArmadillo.h>
#include "riemann_src.h"
#include <map>

using namespace Rcpp;
using namespace arma;
using namespace std;

// =============================================================================
// SPECIAL FUNCTIONS ON SPD MANIFOLD : ADD ON HERE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// =============================================================================
// (01) src_spd_dist    : compute distance of two SPD matrices : ADD ON HERE!
//      src_spd_pdist   : pairwise distances

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
// OTHERS
// =============================================================================
// (01) src_spd_variation : given a 3d array and a frechet mean, compute var


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
// Normalize before factorization so a change of physical units cannot overflow
// a determinant or a matrix product. No regularization is applied here.
static arma::mat spd_scaled_cholesky(const arma::mat& X, double& scale) {
  scale = arma::abs(X).max();
  arma::mat L;
  if (!(scale > 0.0) || !std::isfinite(scale) ||
      !arma::chol(L, X / scale, "lower")) {
    Rcpp::stop("SPD distance requires a numerically positive-definite matrix.");
  }
  return L;
}

static double spd_stein_stable(const arma::mat& X, const arma::mat& Y) {
  double sx, sy;
  arma::mat L = spd_scaled_cholesky(X, sx);
  sy = arma::abs(Y).max();
  arma::mat left = arma::solve(arma::trimatl(L), Y / sy);
  arma::mat relative = arma::solve(arma::trimatl(L), left.t()).t();
  relative = 0.5 * (relative + relative.t());
  arma::vec values;
  if (!arma::eig_sym(values, relative) || !values.is_finite() ||
      arma::any(values <= 0.0)) {
    Rcpp::stop("The relative SPD spectrum is not numerically positive definite.");
  }
  double ratio = sy / sx;
  double logscale = (ratio > 0.0 && std::isfinite(ratio)) ?
    std::log(ratio) : std::log(sy) - std::log(sx);
  double squared = 0.0;
  for (arma::uword i = 0; i < values.n_elem; ++i) {
    // log cosh(log(lambda)/2) is the scalar Stein divergence.
    // log1p preserves its quadratic behavior near lambda = 1.
    double t = 0.5 * std::abs(std::log(values(i)) + logscale);
    double sh = (t < 20.0) ? std::sinh(t / 2.0) : 0.0;
    squared += (t < 20.0) ? std::log1p(2.0 * sh * sh) :
      t + std::log1p(std::exp(-2.0 * t)) - std::log(2.0);
  }
  return std::sqrt(squared);
}

static double spd_wasserstein_stable(const arma::mat& X, const arma::mat& Y) {
  double sx, sy;
  arma::mat A = spd_scaled_cholesky(X, sx);
  arma::mat B = spd_scaled_cholesky(Y, sy);
  double root_scale = std::sqrt(std::max(sx, sy));
  A *= std::sqrt(sx) / root_scale;
  B *= std::sqrt(sy) / root_scale;
  arma::mat U, V;
  arma::vec singular;
  if (!arma::svd(U, singular, V, B.t() * A)) {
    Rcpp::stop("The Wasserstein orthogonal alignment failed.");
  }
  // The Procrustes residual equals the Bures/Wasserstein distance, without
  // subtracting nearly equal traces for identical or neighboring matrices.
  return root_scale * arma::norm(A - B * U * V.t(), "fro");
}

// (01) spd_dist  : compute distance of two SPD matrices -----------------------
double src_spd_dist(arma::mat X, arma::mat Y, std::string geometry){
  if (arma::all(arma::vectorise(X) == arma::vectorise(Y))) return 0.0;
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
    output = spd_stein_stable(X, Y);
  } else if (geometry=="wasserstein"){                                // 5. Wasserstein
    output = spd_wasserstein_stable(X, Y);
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




// =============================================================================
// OTHERS
// =============================================================================
// (01) src_spd_variation : given a 3d array and a frechet mean, compute var
// [[Rcpp::export]]
double src_spd_variation(arma::cube &data3d, arma::mat &fmean){
  int n = data3d.n_slices;
  double output = 0.0;
  double tmpval = 0.0;
  for (int i=0; i<n; i++){
    tmpval = riem_dist("spd",data3d.slice(i),fmean);
    output += (tmpval*tmpval);
  }
  return(output);
}
