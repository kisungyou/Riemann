#include <RcppArmadillo.h>
#include "riemann_src.h"
#include <algorithm>

using namespace Rcpp;
using namespace arma;
using namespace std;

// HELPER FUNCTIONS FOR CLUSTERING
// helper_nunique          : number of unique elements
// helper_centers          : given uvec & datacube, return the centers of cubes
// helper_kmeans_cost      : compute the cost of k-means objective
// helper_assign_centroids : assign each individual to the closest centroid
// helper_kmeans_initlabel : compute the approximate cluster label

// MAIN ALGORITHMS
// 1. clustering_nmshift       : nonlinear mean shift; POTENTIAL TO SPEED UP BY OPENMP
// 2. clustering_kmeans_lloyd  : stop if empty cluster is generated
// 3. clustering_kmeans_macqueen 
// 4. clustering_clrq          : competitive learning riemannian quantization
// 5. clustering_sup_intrinsic : self-updating process
// 6. clustering_kmeans18B     : lightweight coreset

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
arma::uvec helper_kmeans_initlabel(std::string mfdname, arma::cube data, int K){
  int nrow = data.n_rows;
  int ncol = data.n_cols;
  int N    = data.n_slices;
  
  // apply kmeans from Armadillo
  arma::mat logvecs = internal_logvectors(mfdname, data); 
  int P    = logvecs.n_cols;
  arma::mat means(P,K);
  bool status = arma::kmeans(means, arma::trans(logvecs), K, random_subset, 50, false);
  arma::mat centers = means.t();
  
  // compute the distance
  arma::mat distmat(N,K,fill::zeros);
  for (int n=0; n<N; n++){
    for (int k=0; k<K; k++){
      distmat(n,k) = arma::norm(logvecs.row(n)-centers.row(k),2);
    }
  }
  
  // compute the label
  arma::uvec output(N,fill::zeros);
  for (int n=0; n<N; n++){
    output(n) = arma::index_min(distmat.row(n));
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

// 4. clustering_clrq         : competitive learning riemannian quantization
// [[Rcpp::export]]
Rcpp::List clustering_clrq(std::string mfdname, Rcpp::List& data, arma::uvec init_label, double par_a, double par_b){
  // PREPARE DATA : cube is a better option
  int N = data.size();
  arma::mat exemplar = Rcpp::as<arma::mat>(data[0]);
  int nrow = exemplar.n_rows;
  int ncol = exemplar.n_cols;  
  
  arma::cube my_data(nrow,ncol,N,fill::zeros);
  for (int n=0; n<N; n++){
    my_data.slice(n) = Rcpp::as<arma::mat>(data[n]);
  }
  int K = init_label.n_elem;
  
  // PREPARE CENTROIDS : initial label
  arma::cube my_centroids(nrow,ncol,K,fill::zeros);
  for (int k=0; k<K; k++){
    my_centroids.slice(k) = my_data.slice(init_label(k)-1);
  }
  
  // STOCHASTIC GRADIENT DESCENT
  double    tmp_gain = 0.0;
  arma::mat tmp_log(nrow,ncol,fill::zeros);
  arma::vec my_dists(K,fill::zeros);
  arma::uword optid;
  for (int n=0; n<N; n++){
    // 1. compute distance from new observation to centroids
    for (int k=0; k<K; k++){
      my_dists(k) = riem_dist(mfdname, my_centroids.slice(k), my_data.slice(n));
    }
    // 2. pick the centroid with the minimal distance
    optid = my_dists.index_min();
    // 3. update the corresponding centroid
    tmp_gain = par_a/(1.0 + par_b*std::sqrt(static_cast<double>(n+1)));
    tmp_log  = riem_log(mfdname, my_centroids.slice(optid), my_data.slice(n));
    my_centroids.slice(optid) = riem_exp(mfdname, my_centroids.slice(optid), tmp_log, tmp_gain);
  }
  
  // COMPUTE PAIRWISE DISTANCES
  arma::mat pairwise_distance(N,K,fill::zeros);
  for (int n=0; n<N; n++){
    for (int k=0; k<K; k++){
      pairwise_distance(n,k) = riem_dist(mfdname, my_data.slice(n), my_centroids.slice(k));
    }
  }
  
  // RETURN
  Rcpp::List output;
  output["centers"] = my_centroids;
  output["pdist2"]  = pairwise_distance;
  return(output);
}

// 5. clustering_sup_intrinsic : self-updating process =========================
arma::mat clustering_sup_intrinsic_singlemean(std::string mfdname, arma::cube input_data, arma::vec input_weight, arma::mat input_init){
  // PREPARE PARAMETERS with standard choices
  int maxiter = 50;   
  double eps  = 1e-5;
  int nrow = input_data.n_rows;
  int ncol = input_data.n_cols;
  int N    = input_data.n_slices;
  
  arma::mat Xold = input_init;
  arma::mat Xtmp(nrow,ncol,fill::zeros);
  arma::mat Xnew(nrow,ncol,fill::zeros);
  arma::vec Xweight = input_weight/arma::accu(input_weight);
  double    Xinc = 100.0;
  
  for (int it=0; it<maxiter; it++){
    // 1. compute the gradient
    Xtmp.fill(0.0);
    for (int n=0; n<N; n++){
      Xtmp += 2.0*Xweight(n)*riem_log(mfdname, Xold, input_data.slice(n));
    }
    // 2. compute the target
    Xnew = riem_exp(mfdname, Xold, Xtmp, 1.0);
    Xinc = arma::norm(Xold-Xnew,"fro");
    // 3. update
    Xold = Xnew;
    if (Xinc < eps){
      break;
    }
  }
  return(Xold);
}
// [[Rcpp::export]]
Rcpp::List clustering_sup_intrinsic(std::string mfdname, Rcpp::List& data, arma::vec weight, double multiplier, int maxiter, double eps){
  // PREPARE DATA : cube is a better option
  int N = data.size();
  arma::mat exemplar = Rcpp::as<arma::mat>(data[0]);
  int nrow = exemplar.n_rows;
  int ncol = exemplar.n_cols;  
  
  arma::cube my_data_old(nrow,ncol,N,fill::zeros);
  arma::cube my_data_new(nrow,ncol,N,fill::zeros);
  for (int n=0; n<N; n++){
    my_data_old.slice(n) = Rcpp::as<arma::mat>(data[n]);
  }
  
  // PREPARE GAMMA & LAMBDA
  arma::mat my_data_dist(N,N,fill::zeros);
  double gamma = 0.0;
  for (int i=0; i<(N-1); i++){
    for (int j=(i+1); j<N; j++){
      my_data_dist(i,j) = riem_dist(mfdname, my_data_old.slice(i), my_data_old.slice(j));
      my_data_dist(j,i) = my_data_dist(i,j);
      gamma += my_data_dist(i,j)*2.0/(static_cast<double>(N*(N-1)));
    }
  }
  double lambda = multiplier*gamma;
  
  // ITERATION
  arma::mat  now_init(nrow,ncol,fill::zeros);
  arma::vec  now_dist(N,fill::zeros);
  arma::uvec now_smaller(N,fill::zeros);
  arma::vec  now_weight(N,fill::zeros);
  arma::cube now_data(nrow,ncol,1,fill::zeros);
  arma::vec  compare_distance(N,fill::zeros);
  for (int it=0; it<maxiter; it++){
    // 1. compute pairwise distance [except for iteration=0]
    if (it > 0){
      for (int i=0; i<(N-1); i++){
        for (int j=(i+1); j<N; j++){
          my_data_dist(i,j) = riem_dist(mfdname, my_data_old.slice(i), my_data_old.slice(j));
          my_data_dist(j,i) = my_data_dist(i,j);
        }
      }
    }
    // 2. iteratively update for each element
    for (int n=0; n<N; n++){
      // 2-1. reset the current ones
      now_smaller.reset(); 
      now_weight.reset();
      now_data.reset();
      
      // 2-2. update
      now_dist    = my_data_dist.col(n);
      now_smaller = arma::find(now_dist <= gamma);
      if (now_smaller.n_elem < 1){
        my_data_new.slice(n) = my_data_old.slice(n);
      } else {
        // 2-2-1. find the corresponding information
        now_weight = arma::exp(-now_dist.elem(now_smaller)/lambda)%(weight.elem(now_smaller));
        now_data   = my_data_old.slices(now_smaller);
        now_init   = my_data_old.slice(n);  
        
        // 2-2-2. compute the local update
        my_data_new.slice(n) = clustering_sup_intrinsic_singlemean(mfdname, now_data, now_weight, now_init);
      }
    }
    // 3. compute the maximum distance for stopping criterion
    for (int n=0; n<N; n++){
      compare_distance(n) = riem_dist(mfdname, my_data_old.slice(n), my_data_new.slice(n));
    }
    my_data_old = my_data_new;
    if (compare_distance.max() < eps){
      break;
    }
  }
  
  // COMPUTE PAIRWISE DISTANCE
  arma::mat my_solution_pdist(N,N,fill::zeros);
  for (int i=0; i<(N-1); i++){
    for (int j=(i+1); j<N; j++){
      my_solution_pdist(i,j) = riem_dist(mfdname, my_data_old.slice(i), my_data_old.slice(j));
      my_solution_pdist(j,i) = my_solution_pdist(i,j);
    }
  }
  
  // WRAP AND RETURN
  Rcpp::List output;
  output["limits"]   = my_data_old;
  output["distance"] = my_solution_pdist;
  return(output);
}

// 6. clustering_kmeans18B     : weighted lightweight coreset =================
// Shared construction keeps the sampling design identical in both public APIs.
Rcpp::List learning_coreset18B(std::string mfdname, std::string geoname,
  Rcpp::List& data, int M, int myiter, double myeps);

namespace {
double coreset_distance(const std::string& mfd, const std::string& geometry,
                        const arma::mat& x, const arma::mat& y) {
  if (arma::approx_equal(x, y, "absdiff", 0.0)) return 0.0;
  const double value = geometry == "intrinsic" ? riem_dist(mfd, x, y) : riem_distext(mfd, x, y);
  if (!std::isfinite(value) || value < 0.0)
    Rcpp::stop("Coreset clustering encountered invalid distances.");
  return value;
}

arma::mat coreset_distances(const std::string& mfd, const std::string& geometry,
                           const arma::cube& data, const arma::cube& centers) {
  arma::mat result(data.n_slices, centers.n_slices);
  for (arma::uword i = 0; i < data.n_slices; ++i)
    for (arma::uword j = 0; j < centers.n_slices; ++j)
      result(i,j) = coreset_distance(mfd, geometry, data.slice(i), centers.slice(j));
  return result;
}

arma::uvec coreset_assign(const arma::mat& distances) {
  arma::uvec labels(distances.n_rows);
  for (arma::uword i = 0; i < distances.n_rows; ++i) labels(i) = distances.row(i).index_min();
  return labels;
}

double coreset_cost(const arma::mat& distances, const arma::uvec& labels,
                    const arma::vec& weights) {
  double value = 0.0;
  for (arma::uword i = 0; i < labels.n_elem; ++i) {
    const double d = distances(i, labels(i));
    value += weights(i) * d * d;
  }
  if (!std::isfinite(value)) Rcpp::stop("Coreset clustering objective overflowed.");
  return value;
}

arma::cube coreset_initialize(const std::string& mfd, const std::string& geometry,
                             const arma::cube& data, const arma::vec& weights, int K) {
  const int M = data.n_slices;
  arma::cube centers(data.n_rows, data.n_cols, K);
  arma::uvec draw = helper_sample(M, 1, weights / arma::accu(weights), true);
  centers.slice(0) = data.slice(draw(0));
  arma::vec nearest(M);
  for (int i = 0; i < M; ++i)
    nearest(i) = coreset_distance(mfd, geometry, data.slice(i), centers.slice(0));
  for (int k = 1; k < K; ++k) {
    const double scale = nearest.max();
    if (scale == 0.0) {
      // Preserve the original iid draw even if its support is smaller than K.
      centers.slice(k) = centers.slice(0);
    } else {
      arma::vec probability = weights % arma::square(nearest / scale);
      draw = helper_sample(M, 1, probability / arma::accu(probability), true);
      centers.slice(k) = data.slice(draw(0));
      for (int i = 0; i < M; ++i)
        nearest(i) = std::min(nearest(i), coreset_distance(mfd, geometry,
          data.slice(i), centers.slice(k)));
    }
  }
  return centers;
}

arma::mat coreset_weighted_mean(const std::string& mfd, const std::string& geometry,
                               const arma::cube& data, const arma::vec& weights,
                               const arma::uvec& members, const arma::mat& initial) {
  arma::field<arma::mat> subset(members.n_elem);
  for (arma::uword i = 0; i < members.n_elem; ++i) subset(i) = data.slice(members(i));
  RiemannSummaryControl control(200, 1e-8);
  RiemannSummaryResult fit = geometry == "intrinsic" ?
    riem_summary_mean(mfd, subset, weights.elem(members), control, &initial) :
    riem_summary_extrinsic(mfd, subset, weights.elem(members), control, false, &initial);
  if (!fit.converged) Rcpp::stop("A weighted coreset mean did not converge (%s).", fit.termination.c_str());
  return fit.estimate;
}
} // namespace

// [[Rcpp::export]]
Rcpp::List clustering_kmeans18B(std::string mfdname, std::string geotype, Rcpp::List& data, int K, int M, int maxiter){
  const int N = data.size();
  if (N < 1 || K < 1 || K > N || M < K || maxiter < 1)
    Rcpp::stop("Coreset clustering requires 1 <= k <= min(N, M) and positive maxiter.");
  Rcpp::List coreset = learning_coreset18B(mfdname, geotype, data, M, 200, 1e-8);
  const arma::uvec indices = Rcpp::as<arma::uvec>(coreset["id"]);
  const arma::vec probabilities = Rcpp::as<arma::vec>(coreset["qx"]);
  const arma::vec weights = 1.0 / (static_cast<double>(M) * probabilities.elem(indices));
  const arma::mat exemplar = Rcpp::as<arma::mat>(data[0]);
  arma::cube full(exemplar.n_rows, exemplar.n_cols, N);
  for (int i = 0; i < N; ++i) full.slice(i) = Rcpp::as<arma::mat>(data[i]);
  const arma::cube sampled = full.slices(indices);
  arma::cube centers = coreset_initialize(mfdname, geotype, sampled, weights, K);
  arma::mat distances = coreset_distances(mfdname, geotype, sampled, centers);
  arma::uvec labels = coreset_assign(distances);
  double objective = coreset_cost(distances, labels, weights);
  std::vector<double> history(1, objective);
  bool converged = false;
  int iterations = 0;
  std::string termination = "iteration_limit";
  for (int iteration = 0; iteration < maxiter; ++iteration) {
    ++iterations;
    arma::cube updated = centers;
    for (int k = 0; k < K; ++k) {
      const arma::uvec members = arma::find(labels == static_cast<arma::uword>(k));
      if (members.n_elem > 0) updated.slice(k) = coreset_weighted_mean(mfdname,
        geotype, sampled, weights, members, centers.slice(k));
    }
    arma::mat next_distances = coreset_distances(mfdname, geotype, sampled, updated);
    arma::uvec next_labels = coreset_assign(next_distances);
    const double next_objective = coreset_cost(next_distances, next_labels, weights);
    if (next_objective > objective + 1e-10 * std::abs(objective)) {
      termination = "objective_increase_rejected";
      break;
    }
    const bool stable = arma::all(next_labels == labels);
    centers = updated;
    distances = next_distances;
    labels = next_labels;
    objective = next_objective;
    history.push_back(objective);
    if (stable) {
      converged = true;
      termination = "stable_assignment";
      break;
    }
  }
  const arma::mat full_distances = coreset_distances(mfdname, geotype, full, centers);
  const arma::uvec full_labels = coreset_assign(full_distances);
  const double wcss = coreset_cost(full_distances, full_labels, arma::ones<arma::vec>(N));
  return Rcpp::List::create(Rcpp::Named("means") = centers,
    Rcpp::Named("cluster") = full_labels, Rcpp::Named("wcss") = wcss,
    Rcpp::Named("coreid") = indices, Rcpp::Named("weight") = weights,
    Rcpp::Named("iterations") = iterations, Rcpp::Named("converged") = converged,
    Rcpp::Named("termination") = termination, Rcpp::Named("objective_history") = history);
}
