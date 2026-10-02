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
  const arma::uword N = data.size();
  if (N < 2 || label.n_elem != N || !std::isfinite(lambda) || lambda < 0.0)
    Rcpp::stop("RMML requires at least two observations, compatible labels, and finite nonnegative regularization.");
  const arma::mat exemplar = Rcpp::as<arma::mat>(data[0]);
  const arma::vec first = riem_equiv(mfdname, exemplar, exemplar.n_rows, exemplar.n_cols);
  if (first.n_elem == 0 || !first.is_finite())
    Rcpp::stop("RMML encountered an empty or nonfinite equivariant embedding.");
  // Projector embeddings can be wider than the original matrix representation.
  arma::mat eqmat(N, first.n_elem, fill::zeros);
  eqmat.row(0) = first.t();
  for (arma::uword n = 1; n < N; ++n) {
    const arma::mat point = Rcpp::as<arma::mat>(data[n]);
    if (point.n_rows != exemplar.n_rows || point.n_cols != exemplar.n_cols)
      Rcpp::stop("RMML observations have incompatible dimensions.");
    const arma::vec embedded = riem_equiv(mfdname, point, point.n_rows, point.n_cols);
    if (embedded.n_elem != first.n_elem || !embedded.is_finite())
      Rcpp::stop("RMML encountered inconsistent or nonfinite equivariant embeddings.");
    eqmat.row(n) = embedded.t();
  }
  return alg_GMMLreg(eqmat, label, lambda);
}

// 3. learning_coreset18B : lightweight coreset ================================
// [[Rcpp::export]]
Rcpp::List learning_coreset18B(std::string mfdname, std::string geoname, Rcpp::List& data, int M, int myiter, double myeps){
  const int N = data.size();
  if (N < 1 || M < 1 || myiter < 1 || !std::isfinite(myeps) || myeps <= 0.0)
    Rcpp::stop("Coreset construction requires nonempty data and positive sampling and mean controls.");
  if (geoname != "intrinsic" && geoname != "extrinsic")
    Rcpp::stop("Unknown coreset geometry.");
  const arma::mat exemplar = Rcpp::as<arma::mat>(data[0]);
  arma::cube observations(exemplar.n_rows, exemplar.n_cols, N);
  for (int n = 0; n < N; ++n) observations.slice(n) = Rcpp::as<arma::mat>(data[n]);
  const arma::mat mean = internal_mean(mfdname, geoname, observations, myiter, myeps);
  arma::vec distances(N);
  for (int n = 0; n < N; ++n) {
    distances(n) = arma::approx_equal(mean, observations.slice(n), "absdiff", 0.0) ? 0.0 :
      (geoname == "intrinsic" ? riem_dist(mfdname, mean, observations.slice(n)) :
       riem_distext(mfdname, mean, observations.slice(n)));
  }
  if (!distances.is_finite() || arma::any(distances < 0.0))
    Rcpp::stop("Coreset construction encountered invalid distances.");
  arma::vec probability(N, fill::ones);
  probability /= static_cast<double>(N);
  const double scale = distances.max();
  if (scale > 0.0) {
    // Normalize before squaring: q is invariant to a common distance scale.
    const arma::vec squared = arma::square(distances / scale);
    probability = 0.5 * probability + 0.5 * squared / arma::accu(squared);
  }
  probability /= arma::accu(probability);
  // The weight 1 / (M q_i) is valid for independent draws with replacement.
  const arma::uvec indices = helper_sample(N, M, probability, true);
  return Rcpp::List::create(Rcpp::Named("qx") = probability,
                            Rcpp::Named("id") = indices);
}
