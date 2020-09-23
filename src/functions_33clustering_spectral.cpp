#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;
using namespace std;





/* Spectral Clustering Routines
 * (01) cpp_scNJW : Ng, Jordan, Weiss (2002)
 * (02) cpp_scUL  : unnormalized laplacian
 * (03) cpp_scSM  : Shi and Malik (2000)
 * (04) cpp_sc05Z : Zelnik-Manor and Perona (2005)
 */

// AUXILIARY ROUTINES
arma::urowvec label_kmeans(arma::mat data, int K, int maxiter){
  // parameters 
  int N = data.n_rows;
  
  // run k-means
  arma::mat means;
  bool status = arma::kmeans(means, arma::trans(data), K, random_subset, maxiter, false); // it returns (K x K) column means
  if (status == false){
    Rcpp::Rcout << "* k-means failed" << std::endl;
  }
  // need to compute pairwise distance matrix
  arma::mat kdist(K,N,fill::zeros);
  arma::colvec dcoli;
  for (int i=0; i<N; i++){
    dcoli = arma::trans(data.row(i));
    for (int j=0; j<K; j++){
      kdist(j,i) = arma::norm(means.col(j)-dcoli,2);
    }
  }
  urowvec gaus_ids = arma::index_min(kdist, 0);
  return(gaus_ids);
}
arma::urowvec label_gmm(arma::mat data, int K, int maxiter){
  arma::gmm_full model;
  bool status = model.learn(data.t(), K, maha_dist, random_subset, maxiter, maxiter, 1e-10, false);
  if (status == false){
    Rcpp::Rcout << "* GMM failed" << std::endl;
  }
  urowvec gaus_ids = model.assign(data.t(), prob_dist);
  return(gaus_ids);
}
Rcpp::List sc_unnormalized(arma::mat W, int K, bool usekmeans, int maxiter){
  // build laplacian
  arma::mat A = W; A.diag().fill(0.0);
  int  N = A.n_rows;
  arma::vec Dvec = arma::sum(A, 1);
  arma::mat Dmat = arma::diagmat(Dvec);
  arma::mat L = Dmat - A;
  
  arma::vec eigval;
  arma::mat eigvec;
  arma::eig_sym(eigval, eigvec, L);
  
  arma::mat dat = eigvec.head_cols(K); // (N x K) columns are smallest eigenvectors
  arma::urowvec output;
  if (usekmeans==true){
    output = label_kmeans(dat, K, maxiter);
  } else {
    output = label_gmm(dat, K, maxiter);
  }
  
  return Rcpp::List::create(Rcpp::Named("values")=eigval,
                            Rcpp::Named("embeds")=dat,
                            Rcpp::Named("labels")=output);
}
Rcpp::List sc_normalNJW(arma::mat W, int K, bool usekmeans, int maxiter){
  // build laplacian
  arma::mat A = W; A.diag().fill(0.0);
  int N = A.n_rows;
  arma::vec Dvec = arma::sum(A, 1);
  arma::vec Dhalfinv(N,fill::zeros);
  double Dvalue = 0.0;
  for (int n=0;n<N;n++){
    Dvalue = Dvec(n);
    if (Dvalue > arma::datum::eps){
      Dhalfinv(n) = 1.0/std::sqrt(static_cast<float>(Dvalue));
    }
  }
  arma::mat Dhalfmat = arma::diagmat(Dhalfinv);
  arma::mat L = arma::eye(N,N) - Dhalfmat*A*Dhalfmat;
  
  arma::vec eigval;
  arma::mat eigvec;
  arma::eig_sym(eigval, eigvec, L);
  
  arma::mat dat = eigvec.head_cols(K); // (N x K) columns are smallest eigenvectors
  for (int n=0;n<N;n++){
    dat.row(n) = dat.row(n)/arma::norm(dat.row(n),2);
  }
  arma::urowvec output;
  if (usekmeans==true){
    output = label_kmeans(dat, K, maxiter);
  } else {
    output = label_gmm(dat, K, maxiter);
  }
  return Rcpp::List::create(Rcpp::Named("values")=eigval,
                            Rcpp::Named("embeds")=dat,
                            Rcpp::Named("labels")=output);
}
Rcpp::List sc_normalSM(arma::mat W, int K, bool usekmeans, int maxiter){
  // build laplacian
  arma::mat A = W; A.diag().fill(0.0);
  int N = A.n_rows;
  arma::vec Dvec = arma::sum(A, 1);
  arma::vec Dinv(N,fill::zeros);
  double Dvalue = 0.0;
  for (int n=0;n<N;n++){
    Dvalue = Dvec(n);
    if (Dvalue > arma::datum::eps){
      Dinv(n) = 1.0/Dvalue;
    }
  }
  arma::mat Dinvmat = arma::diagmat(Dinv);
  arma::mat L = arma::eye(N,N) - Dinvmat*A;
  
  arma::cx_vec cxval;
  arma::cx_mat cxvec;
  arma::eig_gen(cxval, cxvec, L);
  arma::vec eigval = arma::real(cxval);
  arma::mat eigvec = arma::real(cxvec);
  
  arma::mat dat = eigvec.head_cols(K); // (N x K) columns are smallest eigenvectors
  arma::urowvec output;
  if (usekmeans==true){
    output = label_kmeans(dat, K, maxiter);
  } else {
    output = label_gmm(dat, K, maxiter);
  }
  
  return Rcpp::List::create(Rcpp::Named("values")=eigval,
                            Rcpp::Named("embeds")=dat,
                            Rcpp::Named("labels")=output);
}

// (01) cpp_scNJW : Ng, Jordan, Weiss (2002) ===================================
// [[Rcpp::export]]
Rcpp::List cpp_scNJW(arma::mat& D, int K, double sigma, bool usekmeans, int maxiter){
  // prepare
  int N = D.n_rows;
  
  // compute laplacian with zero diagonals
  arma::mat A = arma::exp(-(D%D)/(sigma*sigma));
  A.diag().fill(0.0);
  
  // run the clustering
  return(sc_normalNJW(A, K, usekmeans, maxiter));
}

// (02) cpp_scUL  : unnormalized laplacian =====================================
// [[Rcpp::export]]
Rcpp::List cpp_scUL(arma::mat& D, int K, double sigma, bool usekmeans, int maxiter){
  // prepare
  int N = D.n_rows;
  
  // compute laplacian with zero diagonals
  arma::mat A = arma::exp(-(D%D)/(sigma*sigma));
  A.diag().fill(0.0); 
  
  // step 3. run the clustering 
  return(sc_unnormalized(A, K, usekmeans, maxiter));
}

// (03) cpp_scSM  : Shi and Malik (2000) =======================================
// [[Rcpp::export]]
Rcpp::List cpp_scSM(arma::mat& D, int K, double sigma, bool usekmeans, int maxiter){
  // prepare
  int N = D.n_rows;
  
  // compute laplacian with zero diagonals
  arma::mat A = arma::exp(-(D%D)/(sigma*sigma));
  A.diag().fill(0.0); 
  
  // run the clustering
  return(sc_normalSM(A, K, usekmeans, maxiter));
}

// (04) cpp_sc05Z : Zelnik-Manor and Perona (2005) =============================
// [[Rcpp::export]]
Rcpp::List cpp_sc05Z(arma::mat& D, int K, int nnbd, bool usekmeans, int maxiter){
  // prepare
  int N = D.n_rows;
  
  // find nnbd-th nearest distance
  arma::vec dist_nnbd(N,fill::zeros);
  arma::vec dist_tmp(N,fill::zeros);
  for (int n=0; n<N; n++){
    dist_tmp     = arma::sort(D.col(n), "ascend");
    dist_nnbd(n) = dist_tmp(nnbd);
  }
  
  // compute laplacian
  arma::mat A(N,N,fill::zeros);
  for (int i=0; i<(N-1); i++){
    for (int j=(i+1); j<N; j++){
      A(i,j) = std::exp(-(D(i,j)*D(i,j))/(dist_nnbd(i)*dist_nnbd(j)));
      A(j,i) = A(i,j);
    }
  }
  
  // compute with NJW
  return(sc_normalNJW(A, K, usekmeans, maxiter));
}