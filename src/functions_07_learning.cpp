#include <RcppArmadillo.h>
#include "riemann_src.h"
#include <algorithm>

using namespace Rcpp;
using namespace arma;
using namespace std;

// 1. learning_seb        : smallest enclosing ball
// 2. learning_rmml       : riemannian manifold metric learning
// 3. learning_coreset18B : lightweight coreset

// 1. learning_seb : smallest enclosing ball ===================================
arma::mat learning_seb_aa2013(std::string mfdname, arma::field<arma::mat> mydata, int myiter, double myeps){
  // PREPARE
  int N = mydata.n_elem;
  int p = mydata(0).n_rows;
  int k = mydata(0).n_cols;
  
  arma::vec initweight(N,fill::ones);
  arma::mat cold = riem_initialize(mfdname, mydata, initweight);
  arma::mat clog(p,k,fill::zeros);
  arma::mat cnew(p,k,fill::zeros);
  arma::mat fnow(p,k,fill::zeros);
  double    cinc = 0.0;
  arma::vec cdists(N,fill::zeros);
  
  // MAIN ITERATION
  for (int it=0; it<myiter; it++){
    // 1. compute distances and find the target
    for (int n=0; n<N; n++){
      cdists(n) = riem_dist(mfdname, cold, mydata(n));
    }
    fnow = mydata(cdists.index_max());
    // 2. compute using geodesic
    clog = riem_log(mfdname, cold, fnow);
    cnew = riem_exp(mfdname, cold, clog, (1.0/static_cast<double>(it+1)));
    // 3. Update
    cinc = arma::norm(cold-cnew, 2);
    cold = cnew;
    if (cinc < myeps){
      break;
    }
  }
  
  // RETURN
  return(cold);
}
// [[Rcpp::export]]
Rcpp::List learning_seb(std::string mfdname, Rcpp::List& data, int myiter, double myeps, std::string method){
  // PREPARE
  int N = data.size();
  arma::field<arma::mat> mydata(N);
  for (int n=0; n<N; n++){
    mydata(n) = Rcpp::as<arma::mat>(data[n]);
  }
  
  // CENTER : BRANCHING THE METHOD
  arma::mat center;
  if (method=="aa2013"){
    center = learning_seb_aa2013(mfdname, mydata, myiter, myeps);
  }
  // RADIUS : COMPUTE THE DISTANCE AND RETURN THE LARGEST
  arma::vec vec_radius(N,fill::zeros);
  for (int n=0; n<N; n++){
    vec_radius(n) = riem_dist(mfdname, center, mydata(n));
  }
  double radius = vec_radius.max();
  
  // RETURN
  Rcpp::List result;
  result["center"] = center;
  result["radius"] = radius;
  return(result);
}

// 2. learning_rmml : riemannian manifold metric learning ======================
arma::mat gcurve(arma::mat A, arma::mat B, double t=0.5){
  arma::mat Asq    = arma::sqrtmat_sympd(A);
  arma::mat Asqinv = arma::inv_sympd(Asq);
  arma::mat C   = Asqinv*B*Asqinv;
  arma::mat Ct  = arma::real(arma::powmat(C, t));
  
  arma::mat output = Asq*Ct*Asq;
  return(output);
}
arma::cube helper_scatter(arma::mat X, arma::uvec label){
  int n = X.n_rows;
  int p = X.n_cols;
  
  arma::mat S(p,p,fill::zeros);
  arma::mat D(p,p,fill::zeros);
  arma::rowvec vecdiff(p,fill::zeros);
  for (int i=0; i<(n-1); i++){
    for (int j=(i+1); j<n; j++){
      vecdiff = X.row(i)-X.row(j);
      if (label(i)==label(j)){
        S += vecdiff.t()*vecdiff;
      } else {
        D += vecdiff.t()*vecdiff;
      }
    }
  }
  arma::cube output = arma::join_slices(S,D);
  return(output);
}
arma::mat helper_mahadist1(arma::mat X, arma::mat A){
  int n = X.n_rows;
  int p = X.n_cols;
  
  arma::mat output(n,n,fill::zeros);
  arma::rowvec vecdiff(p,fill::zeros);
  for (int i=0; i<(n-1); i++){
    for (int j=(i+1); j<n; j++){
      vecdiff = X.row(i)-X.row(j);
      output(i,j) = std::sqrt(arma::accu(vecdiff*A*vecdiff.t()));
      output(j,i) = output(i,j);
    }
  }
  return(output);
} 
arma::mat alg_GMMLreg(arma::mat X, arma::uvec label, double lambda){
  // PARAMETERS
  int n = X.n_rows;
  int p = X.n_cols;
  arma::mat A0(p,p,fill::eye);
  
  
  // COMPUTE S AND D
  arma::cube tmpscatter = helper_scatter(X, label);
  arma::mat S = tmpscatter.slice(0);
  arma::mat D = tmpscatter.slice(1);
  
  
  arma::mat check_Slbd = S + lambda*A0;
  if (arma::rank(check_Slbd) < p){
    Rcpp::stop("* riem.rmml : scatter matrix is rank deficient. We recommend to increase a regularization parameter 'lambda'");
  }
  
  arma::mat LHS = arma::inv_sympd(S + lambda*arma::inv_sympd(A0));
  arma::mat RHS = D + lambda*A0;
  
  // COMPUTE : bilinear form
  arma::mat A    = gcurve(LHS, RHS, 0.5);
  
  // COMPUTE : pairwise distance
  arma::mat MAHA = helper_mahadist1(X, A);
  
  // WRAP AND RETURN
  return(MAHA);
}
// [[Rcpp::export]]
arma::mat learning_rmml(std::string mfdname, Rcpp::List& data, double lambda, arma::uvec label){
  // PREPARE
  arma::mat tmpmat = Rcpp::as<arma::mat>(data[0]);
  int nrow = tmpmat.n_rows;
  int ncol = tmpmat.n_cols;
  int nrcs = nrow*ncol;
  int N    = data.size();
  arma::cube mydata(nrow,ncol,N,fill::zeros);
  for (int n=0; n<N; n++){
    mydata.slice(n) = Rcpp::as<arma::mat>(data[n]);
  }
  
  // EQUIVARIANT EMBEDDING
  arma::mat eqmat(N, nrcs, fill::zeros);
  for (int n=0; n<N; n++){
    eqmat.row(n) = arma::trans(riem_equiv(mfdname, mydata.slice(n), nrow, ncol));
  }

  // COMPUTE GMML
  arma::mat output = alg_GMMLreg(eqmat, label, lambda);
  return(output);
}

// 3. learning_coreset18B : lightweight coreset ================================
// [[Rcpp::export]]
Rcpp::List learning_coreset18B(std::string mfdname, std::string geoname, Rcpp::List& data, int M, int myiter, double myeps){
  // PARAMETER AND DATA PREP
  int N = data.size(); double NN = static_cast<double>(N);
  arma::mat exemplar = Rcpp::as<arma::mat>(data[0]);
  int nrow = exemplar.n_rows;
  int ncol = exemplar.n_cols;  
  arma::cube mydata(nrow,ncol,N,fill::zeros);
  for (int n=0; n<N; n++){
    mydata.slice(n) = Rcpp::as<arma::mat>(data[n]);
  }
  
  // STEP 1. COMPUTE MEAN AND DISTANCE
  arma::mat Xmean = internal_mean(mfdname, geoname, mydata, myiter, myeps);
  arma::vec distsq(N,fill::zeros);
  double    dval = 0.0;
  for (int n=0; n<N; n++){
    if (geoname=="intrinsic"){
      dval = riem_dist(mfdname, Xmean, mydata.slice(n));
    } else {
      dval = riem_distext(mfdname, Xmean, mydata.slice(n));
    }
    distsq(n) = dval*dval;
  }
  double distsqsum = arma::accu(distsq);
  
  // STEP 2. COMPUTE PROBABILITY
  arma::vec probability(N,fill::zeros);
  for (int n=0; n<N; n++){
    probability(n) = (0.5/NN) + (0.5*distsq(n)/distsqsum);
  }
  
  // STEP 3. DRAW INDEX FOR CORESET
  arma::uvec coreid = helper_sample(N, M, probability, false);
  
  // RETURN
  Rcpp::List output;
  output["qx"] = probability;
  output["id"] = coreid;
  return(output);
}