#include <RcppArmadillo.h>
#include "riemann_src.h"
#include <algorithm>

using namespace Rcpp;
using namespace arma;
using namespace std;

// 1. spatial_moran_global    
// 2. spatial_geary_global
// 3. spatial_moran_local

// 1. spatial_moran_globa : Moran's I ==========================================
// [[Rcpp::export]]
Rcpp::List spatial_moran_global(std::string mfdname, 
                                std::string geometry,
                                arma::field<arma::mat> &mydata,
                                arma::mat &W,
                                int ntest){
  // PREP ----------------------------------------------------------------------
  // parameters (m x p x n) dimensionality
  int N = mydata.n_elem;
  int M = mydata(0).n_rows;
  int P = mydata(0).n_cols;
  
  // put the data as a 3d array
  arma::cube data3d(M,P,N,fill::zeros);
  for (int n=0; n<N; n++){
    data3d.slice(n) = mydata(n); 
  }
  
  // COMPUTE -------------------------------------------------------------------
  // frechet mean
  arma::mat frechet_mean = internal_mean(mfdname, geometry, data3d, 100, 1e-8);
  
  // common denoinator
  double denominator = 0.0;
  double dval = 0.0;
  for (int n=0; n<N; n++){
    dval = riem_dist(mfdname, frechet_mean, data3d.slice(n));
    denominator += dval*dval;
  }
  
  // Logarithmic matrices
  arma::cube fixed_logs(M,P,N,fill::zeros);
  for (int n=0; n<N; n++){
    fixed_logs.slice(n) = riem_log(mfdname, frechet_mean, data3d.slice(n));
  }
  
  // pairwise metric computation
  arma::mat fixed_pdmetric(N,N,fill::zeros);
  for (int i=0; i<(N-1); i++){
    for (int j=(i+1); j<N; j++){
      fixed_pdmetric(i,j) = riem_metric(mfdname, frechet_mean, fixed_logs.slice(i), fixed_logs.slice(j));
      fixed_pdmetric(j,i) = fixed_pdmetric(i,j);
    }
  }
  
  // compute the statistic
  double out_statistic = 0.0;
  for (int i=0; i<N; i++){
    for (int j=0; j<N; j++){
      if (W(i,j) > arma::datum::eps){
        out_statistic += W(i,j)*fixed_pdmetric(i,j)/denominator;
      }
    }
  }
  
  // iterate over permutations
  double tmp_statistic = -1;
  arma::uvec perm_vec;
  arma::vec out_permuted(ntest,fill::zeros);
  for (int it=0; it<ntest; it++){
    // get a permutation
    perm_vec.reset();
    perm_vec = arma::randperm(N);
    
    // compute a statistic
    tmp_statistic = 0.0;
    for (int i=0; i<N; i++){
      for (int j=0; j<N; j++){
        if (W(i,j) > arma::datum::eps){
          tmp_statistic += W(i,j)*fixed_pdmetric(perm_vec(i), perm_vec(j))/denominator;
        }
      }
    }
    out_permuted(it) = tmp_statistic;
  }
  
    // RETURN
    Rcpp::List output;
    output["statistic"] = out_statistic;
    output["permuted"] = out_permuted;
    return(output);
}




// 2. spatial_geary_global : Geary's C =========================================
// [[Rcpp::export]]
Rcpp::List spatial_geary_global(std::string mfdname, 
                                std::string geometry,
                                arma::field<arma::mat> &mydata,
                                arma::mat &W,
                                int ntest){
  // PREP ----------------------------------------------------------------------
  // parameters (m x p x n) dimensionality
  int N = mydata.n_elem;
  double NN = static_cast<double>(N);
  
  int M = mydata(0).n_rows;
  int P = mydata(0).n_cols;
  
  // put the data as a 3d array
  arma::cube data3d(M,P,N,fill::zeros);
  for (int n=0; n<N; n++){
    data3d.slice(n) = mydata(n); 
  }
  
  // COMPUTE -------------------------------------------------------------------
  // frechet mean
  arma::mat frechet_mean = internal_mean(mfdname, geometry, data3d, 100, 1e-8);
  
  // common denoinator
  double denominator = 0.0;
  double dval = 0.0;
  for (int n=0; n<N; n++){
    dval = riem_dist(mfdname, frechet_mean, data3d.slice(n));
    denominator += dval*dval*2.0*NN;
  }
  
  // pairwise distance matrix
  arma::mat fixed_pdmat(N,N,fill::zeros);
  for (int i=0; i<(N-1); i++){
    for (int j=(i+1); j<N; j++){
      fixed_pdmat(i,j) = riem_dist(mfdname, data3d.slice(i), data3d.slice(j));
      fixed_pdmat(j,i) = fixed_pdmat(i,j);
    }
  }
  
  // compute the statistic
  double out_statistic = 0.0;
  for (int i=0; i<N; i++){
    for (int j=0; j<N; j++){
      if (W(i,j) > arma::datum::eps){
        dval = fixed_pdmat(i,j);
        out_statistic += (NN-1.0)*W(i,j)*dval*dval/denominator;
      }
    }
  }
  
  // iterate over permutations
  double tmp_statistic = -1;
  arma::uvec perm_vec;
  arma::vec out_permuted(ntest,fill::zeros);
  for (int it=0; it<ntest; it++){
    // get a permutation
    perm_vec.reset();
    perm_vec = arma::randperm(N);
    
    // compute a statistic
    tmp_statistic = 0.0;
    for (int i=0; i<N; i++){
      for (int j=0; j<N; j++){
        if (W(i,j) > arma::datum::eps){
          dval = fixed_pdmat(perm_vec(i), perm_vec(j));
          tmp_statistic += (NN-1.0)*W(i,j)*dval*dval/denominator;
        }
      }
    }
    out_permuted(it) = tmp_statistic;
  }
  
  // RETURN
  Rcpp::List output;
  output["statistic"] = out_statistic;
  output["permuted"] = out_permuted;
  return(output);
}





// 3. spatial_moran_local ======================================================
// [[Rcpp::export]]
Rcpp::List spatial_moran_local(std::string mfdname, 
                               std::string geometry,
                               arma::field<arma::mat> &mydata,
                               arma::mat &W,
                               int ntest){
  // PREP ----------------------------------------------------------------------
  // parameters (m x p x n) dimensionality
  int N = mydata.n_elem; double NN = static_cast<double>(N);
  int M = mydata(0).n_rows;
  int P = mydata(0).n_cols;
  
  // put the data as a 3d array
  arma::cube data3d(M,P,N,fill::zeros);
  for (int n=0; n<N; n++){
    data3d.slice(n) = mydata(n); 
  }
  
  // COMPUTE -------------------------------------------------------------------
  // frechet mean
  arma::mat frechet_mean = internal_mean(mfdname, geometry, data3d, 100, 1e-8);
  
  // common denoinator
  double denominator = 0.0;
  double dval = 0.0;
  for (int n=0; n<N; n++){
    dval = riem_dist(mfdname, frechet_mean, data3d.slice(n));
    denominator += dval*dval/NN;
  }
  
  // Logarithmic matrices
  arma::cube fixed_logs(M,P,N,fill::zeros);
  for (int n=0; n<N; n++){
    fixed_logs.slice(n) = riem_log(mfdname, frechet_mean, data3d.slice(n));
  }
  
  // pairwise metric computation
  arma::mat fixed_pdmetric(N,N,fill::zeros);
  for (int i=0; i<(N-1); i++){
    for (int j=(i+1); j<N; j++){
      fixed_pdmetric(i,j) = riem_metric(mfdname, frechet_mean, fixed_logs.slice(i), fixed_logs.slice(j));
      fixed_pdmetric(j,i) = fixed_pdmetric(i,j);
    }
  }
  
  // compute the statistic
  double val_statistic = 0.0;
  arma::vec out_statistic(N,fill::zeros);
  for (int i=0; i<N; i++){
    val_statistic = 0.0;
    for (int j=0; j<N; j++){
      if (W(i,j) > arma::datum::eps){
        val_statistic += W(i,j)*fixed_pdmetric(i,j)/denominator;
      }
    }
    out_statistic(i) = val_statistic;
  }
  
  
  // iterate over permutation
  double tmp_stat = -1.0;
  arma::uvec perm_vec;
  arma::mat out_permuted(N,ntest,fill::zeros); // stack as rows
  
  for (int it=0; it<ntest; it++){
    // get a permutation
    perm_vec.reset();
    perm_vec = arma::randperm(N);
    
    // run over each indices
    for (int i=0; i<N; i++){
      tmp_stat = 0.0;
      for (int j=0; j<N; j++){
        if (W(i,j) > arma::datum::eps){
          tmp_stat+= W(i,j)*fixed_pdmetric(perm_vec(i), perm_vec(j))/denominator;
        }
      }
      out_permuted(i,it) = tmp_stat;
    }
  }
  
  // RETURN
  Rcpp::List output;
  output["statistic"] = out_statistic; // (N) vector
  output["permuted"] = out_permuted;   // (N x nsim) matrix
  return(output);
}
