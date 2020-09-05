#include <RcppArmadillo.h>
#include "riemann_src.h"

using namespace Rcpp;
using namespace arma;
using namespace std;

// 1. visualize_pga : (intrinsic) method by Fletcher

// 1. visualize_pga : (intrinsic) method by Fletcher ===========================
// [[Rcpp::export]]
Rcpp::List visualize_pga(std::string mfdname, Rcpp::List& data){
  // PREPARE
  arma::mat tmpdata = Rcpp::as<arma::mat>(data[0]);
  int N = data.size();
  int p = tmpdata.n_rows;
  int k = tmpdata.n_cols;
  
  arma::cube mydata(p,k,N,fill::zeros);
  for (int n=0; n<N; n++){
    mydata.slice(n) = Rcpp::as<arma::mat>(data[n]);
  }
  
  // 1. COMPUTE CENTER
  arma::mat mycenter = internal_mean(mfdname, "intrinsic", mydata, 100, 1e-5);
  // 2. COMPUTE LOGARITHM
  arma::mat  singlelog(p,k,fill::zeros);
  arma::mat  rowlogs(N,p*k,fill::zeros);
  for (int n=0; n<N; n++){
    singlelog      = riem_log(mfdname, mycenter, mydata.slice(n));
    rowlogs.row(n) = arma::trans(arma::vectorise(singlelog,0)); 
  }
  // 3. EIGEN-DECOMPOSITION
  arma::vec eigval;
  arma::mat eigvec;
  arma::eig_sym(eigval, eigvec, rowlogs.t()*rowlogs); // increasing order
  
  // 4. DO THE EMBEDDING - If necessary, those tail columns are principal geodesics
  arma::mat cppembed = rowlogs*(eigvec.tail_cols(2));
  
  // WRAP AND RETURN
  Rcpp::List result;
  result["center"] = mycenter;
  result["embed"]  = cppembed;
  return(result);
}
  