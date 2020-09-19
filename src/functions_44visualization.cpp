#include <RcppArmadillo.h>
#include "riemann_src.h"

using namespace Rcpp;
using namespace arma;
using namespace std;

// 1. visualize_pga    : (intrinsic) method by Fletcher
// 2. visualize_kpca   : kernel principal component analysis
// 3. visualize_isomap : weighted distance function

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

// 2. visualize_kpca : kernel principal component analysis =====================
// [[Rcpp::export]]
Rcpp::List visualize_kpca(std::string mfdname, Rcpp::List& data, double sigma, int ndim){
  // PREPARE
  arma::mat tmpdata = Rcpp::as<arma::mat>(data[0]);
  int N = data.size();
  int p = tmpdata.n_rows;
  int k = tmpdata.n_cols;
  
  arma::cube mydata(p,k,N,fill::zeros);
  for (int n=0; n<N; n++){
    mydata.slice(n) = Rcpp::as<arma::mat>(data[n]);
  }
  
  // COMPUTE KERNEL MATRIX : exp(-d^2/(2*sigma^2))
  double dval = 0.0;
  arma::mat mat_kernel(N,N,fill::ones);
  for (int i=0; i<(N-1); i++){
    for (int j=(i+1); j<N; j++){
      dval = riem_distext(mfdname, mydata.slice(i), mydata.slice(j));
      mat_kernel(i,j) = std::exp(-(dval*dval)/(2*sigma*sigma));
      mat_kernel(j,i) = mat_kernel(i,j);
    }
  }
  
  // COMPUTE CENTERED KERNEL MATRIX
  double    NN  = static_cast<double>(N);
  double    NN2 = NN*NN;
  double    sum_all     = arma::accu(mat_kernel);
  arma::vec vec_rowsums = arma::sum(mat_kernel, 1);
  arma::mat mat_centered(N,N,fill::zeros);
  for (int i=0; i<N; i++){
    for (int j=i; j<N; j++){
      if (i==j){
        mat_centered(i,i) = mat_kernel(i,i) - (2.0/NN)*vec_rowsums(i) + (1.0/NN2)*sum_all;
      } else {
        mat_centered(i,j) = mat_kernel(i,j) - (1.0/NN)*vec_rowsums(i) - (1.0/NN)*vec_rowsums(j) + (1.0/NN2)*sum_all;
        mat_centered(j,i) = mat_centered(i,j);
      }
    }
  }
  
  // EIGENDECOMPOSITION FOR KPCA : ASCENDING ORDER
  arma::vec eigval;
  arma::mat eigvec;
  arma::eig_sym(eigval, eigvec, mat_centered);

  // FINALIZE
  Rcpp::List output;
  output["embed"] = mat_centered*eigvec.tail_cols(ndim);
  output["vars"]  = arma::reverse(eigval); // change to descending order
  return(output);
}

// 3. visualize_isomap : weighted distance function ============================
// [[Rcpp::export]]
arma::mat visualize_isomap(std::string mfdname, Rcpp::List& data, std::string geometry, int nnbd){
  // PREPARE
  arma::mat tmpdata = Rcpp::as<arma::mat>(data[0]);
  int N = data.size();
  int p = tmpdata.n_rows;
  int k = tmpdata.n_cols;
  
  arma::cube mydata(p,k,N,fill::zeros);
  for (int n=0; n<N; n++){
    mydata.slice(n) = Rcpp::as<arma::mat>(data[n]);
  }
  
  // COMPUTE PAIRWISE DISTANCE
  arma::mat mat_dist(N,N,fill::zeros);
  for (int i=0; i<(N-1); i++){
    for (int j=(i+1); j<N; j++){
      if (geometry=="intrinsic"){
        mat_dist(i,j) = riem_dist(mfdname, mydata.slice(i), mydata.slice(j));
      } else {
        mat_dist(i,j) = riem_distext(mfdname, mydata.slice(i), mydata.slice(j));
      }
      mat_dist(j,i) = mat_dist(i,j);
    }
  }
  
  // COMPUTE NEAREST NEIGHBOR : K-NN WITH INTERSECTION TYPE
  arma::uvec tmpidx;
  arma::field<arma::uvec> record_minimal(N);
  for (int n=0; n<N; n++){
    tmpidx = arma::sort_index(mat_dist.col(n));
    record_minimal(n) = tmpidx.head(nnbd+1);
  }
  arma::mat mat_index(N,N,fill::zeros);
  arma::uvec uvec1;
  arma::uvec uvec2;
  arma::uvec commons;
  for (int i=0; i<N; i++){
    uvec1 = record_minimal(i);
    for (int j=(i+1); j<N; j++){
      uvec2 = record_minimal(j);
      commons = arma::intersect(uvec1, uvec2);
      if (commons.n_elem > 0){
        mat_index(i,j) = 1.0;
        mat_index(j,i) = 1.0;
      }
    }
  }
  
  // RETURN THE WEIGHTED
  arma::mat output = mat_dist%mat_index;
  return(output);
}