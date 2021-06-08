#include <RcppArmadillo.h>
#include "riemann_src.h"

using namespace Rcpp;
using namespace arma;
using namespace std;

// HELPER FUNCTIONS FOR CLUSTERING
// helper_nunique          : number of unique elements
// helper_centers          : given uvec & datacube, return the centers of cubes
// helper_kmeans_cost      : compute the cost of k-means objective
// helper_assign_centroids : assign each individual to the closest centroid

// MAIN ALGORITHMS
// 1. clustering_nmshift      : nonlinear mean shift; POTENTIAL TO SPEED UP BY OPENMP
// 2. clustering_kmeans_lloyd : stop if empty cluster is generated
// 3. clustering_kmeans_macqueen 

int helper_nunique(arma::uvec x){
  arma::uvec y = arma::unique(x);
  int numy = y.n_elem;
  return(numy);
}
// helper_kmeans_cost      : compute the cost of k-means objective
double helper_kmeans_cost(std::string mfd, std::string dtype, arma::cube data, arma::cube centroids, arma::uvec label){
  int K  = centroids.n_slices;
  double output = 0.0;
  arma::uvec idx;
  arma::cube subdata;
  double tmpval;
  for (int k=0; k<K; k++){
    idx = arma::find(label==k);
    if (idx.n_elem > 1){
      subdata = data.slices(idx);
      int M = subdata.n_slices;
      for (int m=0; m<M; m++){
        if (dtype=="intrinsic"){
          tmpval = riem_dist(mfd, subdata.slice(m), centroids.slice(k));
        } else {
          tmpval = riem_distext(mfd, subdata.slice(m), centroids.slice(k));
        }
        output += tmpval*tmpval;
      }
    }
  }
  return(output);
}
arma::cube helper_centers(std::string mfd, std::string dtype, arma::cube data, arma::uvec label){
  int K = helper_nunique(label);
  int nrow = data.n_rows;
  int ncol = data.n_cols;
  
  arma::uvec idx;
  arma::cube output(nrow,ncol,K,fill::zeros);
  for (int k=0; k<K; k++){
    idx = arma::find(label==k);
    if (idx.n_elem==1){
      output.slice(k) = data.slice(idx(0));
    } else if (idx.n_elem > 1){
      output.slice(k) = internal_mean(mfd,dtype,data.slices(idx),50,1e-5);
    } else if (idx.n_elem < 1){
      output.slice(k).fill(arma::datum::nan);
    }
  }
  return(output);
}
arma::uvec helper_assign_centroids(std::string mfd, std::string dtype, arma::cube data, arma::cube centroids){
  int N = data.n_slices;
  int K = centroids.n_slices;
  
  arma::mat distmat(N,K,fill::zeros);
  for (int n=0; n<N; n++){
    for (int k=0; k<K; k++){
      if (dtype=="intrinsic"){
        distmat(n,k) = riem_dist(mfd, data.slice(n), centroids.slice(k));
      } else {
        distmat(n,k) = riem_distext(mfd, data.slice(n), centroids.slice(k));
      }
    }
  }
  
  arma::uvec output(N,fill::zeros);
  for (int n=0; n<N; n++){
    output(n) = index_min(distmat.row(n));
  }
  return(output);
}


// main 1. clustering_nmshift : nonlinear mean shift by Subbarao  ==============
arma::mat clustering_nmshift_single(std::string mfdname, int id, arma::field<arma::mat> mydata, double myh, int myiter, double myeps){
  // PREPARE
  int N = mydata.n_elem;
  int nrow = mydata(0).n_rows;
  int ncol = mydata(0).n_cols;
  
  arma::mat Yold = mydata(id);
  arma::mat Ytmp(nrow,ncol,fill::zeros);
  arma::mat Ynew(nrow,ncol,fill::zeros);
  double    Yinc = 0.0;
  arma::vec Ydists(N,fill::zeros);
  
  arma::mat term1(nrow, ncol, fill::zeros);
  double    term2 = 0.0;
  double    gval = 0.0;
  double    h2   = myh*myh;
  
  // ITERATION
  for (int it=0; it<myiter; it++){
    // 1. compute distances
    for (int n=0; n<N; n++){
      Ydists(n) = riem_dist(mfdname, Yold, mydata(n));
    }
    // 2. compute the updater
    term1.fill(0.0);
    term2 = 0.0;
    for (int n=0; n<N; n++){
      gval = std::exp(-(Ydists(n)*Ydists(n))/h2);
      term1 += gval*riem_log(mfdname, Yold, mydata(n));
      term2 += gval;
    }
    Ytmp = term1/term2;
    Ynew = riem_exp(mfdname, Yold, Ytmp, 1.0);
    // 3. update and stopping criterion
    Yinc = arma::norm(Yold-Ynew,"fro");
    Yold = Ynew;
    if (Yinc < myeps){
      break;
    }
  }
  
  // RETURN
  return(Yold);
}
// [[Rcpp::export]]
Rcpp::List clustering_nmshift(std::string mfdname, Rcpp::List& data, double h, int iter, double eps){
  // PREPARE
  int N = data.size();
  arma::field<arma::mat> mydata(N);
  for (int n=0; n<N; n++){
    mydata(n) = Rcpp::as<arma::mat>(data[n]);
  }
  int nrow = mydata(0).n_rows;
  int ncol = mydata(0).n_cols;
  
  // COMPUTE ---------------------- POSSIBLY, PARALLEL LATER ----------------------------------------------------------------------------------
  arma::cube limpts(nrow,ncol,N,fill::zeros);
  for (int n=0; n<N; n++){
    limpts.slice(n) = clustering_nmshift_single(mfdname, n, mydata, h, iter, eps);
  }
  
  // PAIRWISE DISTANCES
  arma::mat pdmat(N,N,fill::zeros);
  for (int i=0; i<(N-1); i++){
    for (int j=(i+1); j<N; j++){
      pdmat(i,j) = riem_dist(mfdname, limpts.slice(i), limpts.slice(j));
      pdmat(j,i) = pdmat(i,j);
    }
  }
  
  // RETURN
  Rcpp::List output;
  output["limits"]   = limpts;
  output["distance"] = pdmat;
  return(output);
}

// 2. clustering_kmeans_lloyd : stop if empty cluster is generated =============
// [[Rcpp::export]]
Rcpp::List clustering_kmeans_lloyd(std::string mfdname, std::string geotype, Rcpp::List& data, int iter, double eps, arma::uvec initlabel){
  // PREPARE
  // data : for clustering, cube is a better option.
  int N = data.size();
  arma::mat exemplar = Rcpp::as<arma::mat>(data[0]);
  int nrow = exemplar.n_rows;
  int ncol = exemplar.n_cols;
  
  arma::cube mydata(nrow,ncol,N,fill::zeros);
  for (int n=0; n<N; n++){
    mydata.slice(n) = Rcpp::as<arma::mat>(data[n]);
  }
  
  // labels 
  arma::uvec oldlabel = initlabel - initlabel.min(); 
  arma::cube oldmeans = helper_centers(mfdname, geotype, mydata, oldlabel);
  if (oldmeans.has_nan()){
    Rcpp::stop("Lloyd's Algorithm Terminated at Initialization.");
  }
  int K     = oldmeans.n_slices;
  double KK = static_cast<double>(K);
  
  arma::uvec newlabel(N,fill::zeros);
  arma::cube newmeans(nrow,ncol,K,fill::zeros);
  double     meanincs = 0.0;
  
  // MAIN ITERATION
  for (int it=0; it<iter; it++){
    // 1. Assignment Step
    //    If empty cluster appears, stop.
    newlabel = helper_assign_centroids(mfdname, geotype, mydata, oldmeans);
    if (helper_nunique(newlabel) < K){
      break;
    }
    // 2. Update Step
    newmeans = helper_centers(mfdname, geotype, mydata, newlabel);
    if (newmeans.has_nan()){
      break;
    }
    // 3. Termination : average increment of means is very small
    meanincs = 0.0;
    for (int k=0; k<K; k++){
      meanincs += arma::norm(oldmeans.slice(k)-newmeans.slice(k),"fro")/KK;
    }
    oldlabel = newlabel;
    oldmeans = newmeans;
    if (meanincs < eps){
      break;
    }
  }
  
  // SSE
  double SSE = helper_kmeans_cost(mfdname, geotype, mydata, oldmeans, oldlabel);
  
  // RETURN
  Rcpp::List output;
  output["label"] = oldlabel;
  output["means"] = oldmeans;
  output["WCSS"]   = SSE;
  return(output);
}

//3. clustering_kmeans_macqueen 
// [[Rcpp::export]]
Rcpp::List clustering_kmeans_macqueen(std::string mfdname, std::string geotype, Rcpp::List& data, int iter, double eps, arma::uvec initlabel){
  // PREPARE
  // data : for clustering, cube is a better option.
  int N = data.size();
  arma::mat exemplar = Rcpp::as<arma::mat>(data[0]);
  int nrow = exemplar.n_rows;
  int ncol = exemplar.n_cols;
  
  arma::cube mydata(nrow,ncol,N,fill::zeros);
  for (int n=0; n<N; n++){
    mydata.slice(n) = Rcpp::as<arma::mat>(data[n]);
  }
  
  // labels 
  arma::uvec oldlabel = initlabel - initlabel.min(); 
  arma::cube oldmeans = helper_centers(mfdname, geotype, mydata, oldlabel);
  if (oldmeans.has_nan()){
    Rcpp::stop("MacQueen's Algorithm Terminated at Initialization.");
  }
  int K     = oldmeans.n_slices;
  double KK = static_cast<double>(K);
  
  arma::uvec newlabel(N,fill::zeros);
  arma::cube newmeans(nrow,ncol,K,fill::zeros);
  double     meanincs = 0.0;
  
  // MAIN ITERATION
  arma::uword idnow;
  arma::uvec update_order;
  arma::vec  distctd(K,fill::zeros);
  arma::uword newclass;
  arma::uword oldclass;
  for (int it=0; it<iter; it++){
    // Random Permutation
    update_order = arma::randperm(N);
    newlabel = oldlabel;
    newmeans = oldmeans;
    for (int n=0; n<N; n++){
      // 1. compute distance to the centroids
      idnow = update_order(n);
      for (int k=0; k<K; k++){
        if (geotype=="intrinsic"){
          distctd(k) = riem_dist(mfdname, mydata.slice(idnow), newmeans.slice(k));  
        } else {
          distctd(k) = riem_distext(mfdname, mydata.slice(idnow), newmeans.slice(k));
        }
      }
      // 2. re-assign to the nearest and re-compute
      newclass = distctd.index_min();
      oldclass = newlabel(idnow);
      if (oldclass!=newclass){
        newlabel(idnow) = newclass;
        newmeans.slice(oldclass) = internal_mean_init(mfdname, geotype, mydata.slices(arma::find(newlabel==oldclass)), 50, 1e-5, newmeans.slice(oldclass));
        newmeans.slice(newclass) = internal_mean_init(mfdname, geotype, mydata.slices(arma::find(newlabel==newclass)), 50, 1e-5, newmeans.slice(newclass));
      }
    }
    if (helper_nunique(newlabel) < K){ // if there is any empty cluster, stop.
      break;
    }
    
    // Update & Termination
    meanincs = 0.0;
    for (int k=0; k<K; k++){
      meanincs += arma::norm(oldmeans.slice(k)-newmeans.slice(k),"fro")/KK;
    }
    oldlabel = newlabel;
    oldmeans = newmeans;
    if (meanincs < eps){
      break;
    }
  }
  
  // SSE
  double SSE = helper_kmeans_cost(mfdname, geotype, mydata, oldmeans, oldlabel);
  
  // RETURN
  Rcpp::List output;
  output["label"] = oldlabel;
  output["means"] = oldmeans;
  output["WCSS"]   = SSE;
  return(output);
}