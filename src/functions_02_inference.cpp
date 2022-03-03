#include <RcppArmadillo.h>
#include "riemann_src.h"

using namespace Rcpp;
using namespace arma;
using namespace std;

// 1. inference_mean_intrinsic
//    inference_mean_extrinsic
// 2. inference_median_intrinsic
//    inference_median_extrinsic

// 1. inference_mean_intrinsic/extrinsic =======================================
// [[Rcpp::export]]
Rcpp::List inference_mean_intrinsic(std::string mfdname, Rcpp::List& data, arma::vec myweight, int myiter, double myeps){
  // PREPARE
  int N = data.size();
  arma::field<arma::mat> mydata(N);
  for (int n=0; n<N; n++){
    mydata(n) = Rcpp::as<arma::mat>(data[n]);
  }
  int prow = mydata(0).n_rows;
  int pcol = mydata(0).n_cols;
  
  // COMPUTE : MEAN
  // initialization
  arma::mat Sold= riem_initialize(mfdname, mydata, myweight);
  arma::mat Stmp(prow, pcol, fill::zeros);
  arma::mat Snew(prow, pcol, fill::zeros);
  double    Sinc = 0.0;
  
  // iteration
  for (int it=0; it<myiter; it++){
    Stmp.fill(0.0);          // reset the temporary gradient matrix
    for (int n=0; n<N; n++){ // inner iteration
      Stmp += 2.0*myweight(n)*riem_log(mfdname, Sold, mydata(n));
    }
    Snew = riem_exp(mfdname, Sold, Stmp, 1.0);
    Sinc = arma::norm(Sold-Snew,"fro");
    Sold = Snew;
    if (Sinc < myeps){
      break;
    }
  }
  
  // COMPUTE : DISTANCE AND VARIATION
  double variation = 0.0;
  double tmpdval = 0.0;
  arma::vec distvec(N,fill::zeros);
  for (int n=0; n<N; n++){
    tmpdval    = riem_dist(mfdname, Sold, mydata(n));
    distvec(n) = tmpdval;
    variation += myweight(n)*tmpdval*tmpdval;
  }
  
  // RETURN
  Rcpp::List output;
  output["mean"] = Sold;
  output["variation"] = variation;
  output["distvec"]   = distvec;
  return(output);
}
// [[Rcpp::export]]
Rcpp::List inference_mean_extrinsic(std::string mfdname, Rcpp::List& data, arma::vec myweight, int myiter, double myeps){
  // PREPARE
  int N = data.size();
  arma::field<arma::mat> mydata(N);
  for (int n=0; n<N; n++){
    mydata(n) = Rcpp::as<arma::mat>(data[n]);
  }
  int prow = mydata(0).n_rows;
  int pcol = mydata(0).n_cols;
  
  // COMPUTE : MEAN
  // initialization
  arma::vec exemplar = riem_equiv(mfdname, mydata(0), prow, pcol);
  int eqdim = exemplar.n_elem;
  arma::vec Soldvec(eqdim, fill::zeros);
  
  // iteration
  for (int n=0; n<N; n++){
    Soldvec += myweight(n)*riem_equiv(mfdname, mydata(n), prow, pcol);
  }
  arma::mat Sold = riem_invequiv(mfdname, Soldvec, prow, pcol);
  
  // COMPUTE : DISTANCE AND VARIATION
  double variation = 0.0;
  double tmpdval = 0.0;
  arma::vec distvec(N,fill::zeros);
  for (int n=0; n<N; n++){
    tmpdval    = riem_distext(mfdname, Sold, mydata(n));
    distvec(n) = tmpdval;
    variation += myweight(n)*tmpdval*tmpdval;
  }
  
  // RETURN
  Rcpp::List output;
  output["mean"] = Sold;
  output["variation"] = variation;
  output["distvec"]   = distvec;
  return(output);
}

// 2. inference_median_intrinsic/extrinsic =====================================
arma::vec run_weiszfeld(arma::mat X, arma::vec weight, int myiter, double myeps){
  // parameter
  int n = X.n_cols;
  int p = X.n_rows;
  
  // initialize
  arma::vec mold = arma::mean(X, 1);
  arma::vec mnew(p,fill::zeros);
  double    minc = 10000.0;
  
  arma::vec tmp1(p,fill::zeros);
  double    tmp2 = 0.0;
  
  arma::vec vec_norm(n,fill::zeros);
  arma::vec vec_weight = weight/arma::accu(weight);
  
  // let's iterate
  int M = 0;
  arma::uvec nonsingular;
  for (int it=0; it<myiter; it++){
    // 1. compute log-pulled vectors and norm
    for (int i=0; i<n; i++){
      vec_norm(i) = arma::norm(mold-X.col(i), 2);
    }
    // 2. find the one with non-singular distance;
    nonsingular = arma::find(vec_norm > 1e-10);
    M = nonsingular.n_elem;
    if (M < 1){
      break;
    }
    // 3. update numerator & denominator
    tmp1.fill(0.0);
    tmp2 = 0.0;
    for (int j=0; j<M; j++){
      tmp1 += vec_weight(nonsingular(j))*X.col(nonsingular(j))/vec_norm(nonsingular(j));
      tmp2 += vec_weight(nonsingular(j))/vec_norm(nonsingular(j));
    }
    // 4. compute the updated solution
    mnew = tmp1/tmp2;
    // 5. update
    minc = arma::norm(mold-mnew,2);
    mold = mnew;
    if (minc < myeps){
      break;
    }
  }
  
  // update
  return(mold);
}
// [[Rcpp::export]]
Rcpp::List inference_median_intrinsic(std::string mfdname, Rcpp::List& data, arma::vec myweight, int myiter, double myeps){
  // PREPARE
  int N = data.size();
  arma::field<arma::mat> mydata(N);
  for (int n=0; n<N; n++){
    mydata(n) = Rcpp::as<arma::mat>(data[n]);
  }
  int prow = mydata(0).n_rows;
  int pcol = mydata(0).n_cols;
  
  // COMPUTE : MEDIAN
  // initialization
  arma::mat  Sold= riem_initialize(mfdname, mydata, myweight);
  arma::mat  Snew(prow,pcol,fill::zeros);
  arma::mat  Stmp(prow,pcol,fill::zeros);
  arma::cube Slogs(prow,pcol,N,fill::zeros);
  arma::vec  Sdist(N,fill::zeros);
  double     Sinc = 0.0;
  
  // iteration
  arma::uvec nonsingular;
  arma::mat  tmp1(prow,pcol,fill::zeros);
  double     tmp2 = 0.0;
  for (int it=0; it<myiter; it++){
    // 1. compute log-pulled vectors and norm
    for (int n=0; n<N; n++){
      Stmp = riem_log(mfdname, Sold, mydata(n));
      Slogs.slice(n) = Stmp;
      Sdist(n) = std::sqrt(riem_metric(mfdname, Sold, Stmp, Stmp));
    }
    // 2. find the one with singular-distance
    nonsingular = arma::find(Sdist > 1e-10);
    int M = nonsingular.n_elem;
    if (M < 1){
      break;
    }
    // 3. update numerator & denominator
    tmp1.fill(0.0);
    tmp2 = 0.0;
    for (int j=0; j<M; j++){
      tmp1 += myweight(nonsingular(j))*Slogs.slice(nonsingular(j))/Sdist(nonsingular(j));
      tmp2 += myweight(nonsingular(j))/Sdist(nonsingular(j));
    }
    // 4. update to the new one
    Stmp = tmp1/tmp2;
    Snew = riem_exp(mfdname, Sold, Stmp, 1.0);
    Sinc = arma::norm(Sold-Snew,"fro");
    Sold = Snew;
    // 5. update information
    if (Sinc < myeps){
      break;
    }
  }
  
  // COMPUTE : DISTANCE AND VARIATION
  double variation = 0.0;
  double tmpdval   = 0.0;
  arma::vec distvec(N,fill::zeros);
  for (int n=0; n<N; n++){
    tmpdval    = riem_dist(mfdname, Sold, mydata(n));
    distvec(n) = tmpdval;
    variation += myweight(n)*tmpdval;
  }
  
  // WRAP AND RETURN  
  Rcpp::List result;
  result["median"]    = Sold;
  result["variation"] = variation;
  result["distvec"]   = distvec;
  return(result);
}
// [[Rcpp::export]]
Rcpp::List inference_median_extrinsic(std::string mfdname, Rcpp::List& data, arma::vec myweight, int myiter, double myeps){
  // PREPARE
  int N = data.size();
  arma::field<arma::mat> mydata(N);
  for (int n=0; n<N; n++){
    mydata(n) = Rcpp::as<arma::mat>(data[n]);
  }
  int prow = mydata(0).n_rows;
  int pcol = mydata(0).n_cols;
  
  // EXTRINSIC COMPUTATION
  // 1. prepare
  arma::mat objcols(prow*pcol, N, fill::zeros);
  for (int n=0; n<N; n++){
    objcols.col(n) = riem_equiv(mfdname, mydata(n), prow, pcol);
  }
  // 2. compute
  arma::vec extmed = run_weiszfeld(objcols, myweight, myiter, myeps);
  arma::mat outmat = riem_invequiv(mfdname, extmed, prow, pcol);
  
  // DISTANCE AND VARIATION
  double variation = 0.0;
  double tmpdval   = 0.0;
  arma::vec distvec(N,fill::zeros);
  for (int n=0; n<N; n++){
    tmpdval    = riem_distext(mfdname, outmat, mydata(n));
    distvec[n] = tmpdval;
    variation += myweight(n)*tmpdval;
  }
  
  // WRAP AND RETURN  ----------------------------------------------------------
  Rcpp::List result;
  result["median"]    = outmat;
  result["variation"] = variation;
  result["distvec"]   = distvec;
  return(result);
}

