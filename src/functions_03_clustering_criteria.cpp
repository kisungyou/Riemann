#include <RcppArmadillo.h>
#include "riemann_src.h"
#include <algorithm>

using namespace Rcpp;
using namespace arma;
using namespace std;

// AUXILIARY FUNCTIONS
// cvi_helper_classindex : per-class index
// cvi_helper_classmean  : per-class mean
// cvi_helper_pdist      : compute pairwise distance
// cvi_helper_distance   : two matrices distance without ifelse
// cvi_helper_nw         : the number of pairs in all classes


double cvi_helper_distance(std::string mfd, std::string dtype, arma::mat X, arma::mat Y){
  double output = 0.0;
  if (dtype=="intrinsic"){
    output = riem_dist(mfd, X, Y);
  } else {
    output = riem_distext(mfd, X, Y);
  }
  return(output);
}
arma::field<arma::uvec> cvi_helper_classindex(arma::uvec label){
  int N = label.n_elem;
  int K = label.max() + 1;
  
  arma::field<arma::uvec> output(K);
  for (int k=0; k<K; k++){
    output(k) = arma::find(label==k);
  }
  return(output);
}
arma::cube cvi_helper_classmean(std::string mfd, std::string dtype, arma::cube X, arma::uvec label){
  int nrow = X.n_rows;
  int ncol = X.n_cols;
  int N    = X.n_slices;
  int K    = label.max() + 1;
  
  arma::uvec cid;
  arma::cube output(nrow,ncol,K,fill::zeros);
  for (int k=0; k<K; k++){
    cid.reset();
    cid = arma::find(label==k);
    if (cid.n_elem < 2){
      output.slice(k) = X.slice(cid(0));
    } else {
      output.slice(k) = internal_mean(mfd, dtype, X.slices(cid), 100, 1e-5);
    }
  }
  return(output);
}
arma::mat cvi_helper_pdist(std::string mfd, std::string dtype, arma::cube X){
  int N    = X.n_slices;
  int nrow = X.n_rows;
  int ncol = X.n_cols;
  
  
  arma::mat output(N,N,fill::zeros);
  for (int i=0; i<(N-1); i++){
    for (int j=(i+1); j<N; j++){
      output(i,j) = cvi_helper_distance(mfd, dtype, X.slice(i), X.slice(j));
      output(j,i) = output(i,j);
    }
  }
  return(output);
}
int cvi_helper_nw(arma::uvec label){
  arma::field<arma::uvec> classindex = cvi_helper_classindex(label);
  int K = classindex.n_elem;
  
  int tmp = 0;
  int output = 0;
  for (int k=0; k<K; k++){
    tmp = classindex(k).n_elem;
    output += tmp*(tmp-1)/2;
  }
  return(output);
}

// INTERNAL CRITERIA FUNCTION
// (01; +) Dunn index
// (02; +) Calinski-Harabasz
// (03; -) Gamma : Memory Loss : Implement in R
// (04; -) C Index
// (05; -) Davies-Bouldin
// (06; +) Silhouette
// (07; +) Generalized Dunn
// (08; +) Score function

// (08; +) Score function ------------------------------------------------------
// [[Rcpp::export]]
double cvi_internal_score(std::string mfd, std::string dtype,  Rcpp::List& data, arma::uvec mylabel){
  // PRELIMINARY
  arma::uvec label = mylabel - mylabel.min(); // starting from 0 to (k-1)
  arma::mat exmat  = Rcpp::as<arma::mat>(data[0]);
  int nrow  = exmat.n_rows;
  int ncol  = exmat.n_cols;
  int N     = data.size();
  int K     = label.max() + 1;
  double NK = static_cast<double>(N*K);
  arma::cube X(nrow,ncol,N,fill::zeros);
  for (int n=0; n<N; n++){
    X.slice(n) = Rcpp::as<arma::mat>(data[n]);
  }
  
  // COMPUTE : DISTANCE MATRIX & CLASS INDEX
  arma::mat distmat = cvi_helper_pdist(mfd, dtype, X);
  arma::field<arma::uvec> classindex = cvi_helper_classindex(label);
  
  arma::cube mean_class  = cvi_helper_classmean(mfd,dtype,X,label);
  arma::mat  mean_global = internal_mean(mfd,dtype,X,100,1e-5);
  
  // compute BCD
  double bcd = 0.0;
  for (int k=0; k<K; k++){
    bcd += static_cast<double>(classindex(k).n_elem)*cvi_helper_distance(mfd,dtype,mean_class.slice(k),mean_global)/NK;
  }
  // compute WCD
  arma::vec vec_wcd(K,fill::zeros);
  for (int k=0; k<K; k++){
    arma::uvec idk = classindex(k);
    double ck  = static_cast<double>(idk.n_elem);
    double tmp = 0.0;
    for (int it=0; it<ck; it++){
      tmp += cvi_helper_distance(mfd,dtype,X.slice(idk(it)),mean_class.slice(k));
    }
    vec_wcd(k) = tmp/ck;
  }
  double wcd = arma::accu(vec_wcd);
  
  // compute output
  return(1 - (1.0/std::exp(std::exp(wcd+bcd))));
}

// (07; +) Generalized Dunn ----------------------------------------------------
// [[Rcpp::export]]
double cvi_internal_gdxx(std::string mfd, std::string dtype,  Rcpp::List& data, arma::uvec mylabel, int delta, int Delta){
  // PRELIMINARY
  arma::uvec label = mylabel - mylabel.min(); // starting from 0 to (k-1)
  arma::mat exmat  = Rcpp::as<arma::mat>(data[0]);
  int nrow = exmat.n_rows;
  int ncol = exmat.n_cols;
  int N    = data.size();
  int K    = label.max() + 1;
  arma::cube X(nrow,ncol,N,fill::zeros);
  for (int n=0; n<N; n++){
    X.slice(n) = Rcpp::as<arma::mat>(data[n]);
  }
  
  // COMPUTE : DISTANCE MATRIX & CLASS INDEX
  arma::mat distmat = cvi_helper_pdist(mfd, dtype, X);
  arma::field<arma::uvec> classindex = cvi_helper_classindex(label);
  
  // auxiliary computation
  arma::cube mean_class = cvi_helper_classmean(mfd,dtype,X,label);
  
  // denominator
  arma::vec vec_den(K,fill::zeros);
  arma::vec vec_tmp(K,fill::zeros);
  for (int k=0; k<K; k++){
    vec_tmp.fill(0.0);
    if (delta==4){
      for (int j=0; j<K; j++){
        if (k==j){
          vec_tmp(j) = arma::datum::inf;
        } else {
          vec_tmp(j) = cvi_helper_distance(mfd,dtype,mean_class.slice(j),mean_class.slice(k));
        }
      }
      vec_den(k) = vec_tmp.min();
    } else if (delta==3){
      for (int j=0; j<K; j++){
        if (j==k){
          vec_tmp(j) = arma::datum::inf;
        } else {
          vec_tmp(j) = arma::accu(distmat(classindex(k),classindex(j)))/static_cast<double>(classindex(k).n_elem*classindex(j).n_elem);
        }
      }
      vec_den(k) = vec_tmp.min();
    } else if (delta==5){
      for (int j=0; j<K; j++){
        if (j==k){
          vec_tmp(j) = arma::datum::inf;
        } else {
          arma::uvec idck = classindex(k); int nk = idck.n_elem;
          arma::uvec idcj = classindex(j); int nj = idcj.n_elem;
          
          double sumk = 0.0;
          double sumj = 0.0;
          
          for (int it=0; it<nk; it++){
            sumk += cvi_helper_distance(mfd,dtype,mean_class.slice(k),X.slice(idck(it)));
          }
          for (int it=0; it<nj; it++){
            sumj += cvi_helper_distance(mfd,dtype,mean_class.slice(j),X.slice(idcj(it)));
          }
          vec_tmp(j) = (sumk+sumj)/static_cast<double>(nk+nj);
        }
      }
      vec_den(k) = vec_tmp.min();
    }
  }
  double val_den = vec_den.min();
  
  // compute numerator
  arma::vec vec_num(K,fill::zeros);
  if (Delta==1){
    for (int k=0; k<K; k++){
      arma::mat partmat = distmat(classindex(k), classindex(k));
      vec_num(k) = partmat.max();
    }
  } else if (Delta==3){
    for (int k=0; k<K; k++){
      double sumk = 0.0;
      arma::uvec idk = classindex(k);
      int nk = idk.n_elem;
      for (int it=0; it<nk; it++){
        sumk += cvi_helper_distance(mfd,dtype,mean_class.slice(k),X.slice(idk(it)))*2.0/static_cast<double>(nk);
      }
      vec_num(k) = sumk;
    }
  }
  double val_num = vec_num.max();
  
  double output = val_den/val_num;
  return(output);
}

// (05; -) Davies-Bouldin ------------------------------------------------------
// [[Rcpp::export]]
double cvi_internal_db(std::string mfd, std::string dtype, Rcpp::List& data, arma::uvec mylabel){
  // PRELIMINARY
  arma::uvec label = mylabel - mylabel.min(); // starting from 0 to (k-1)
  arma::mat exmat  = Rcpp::as<arma::mat>(data[0]);
  int nrow = exmat.n_rows;
  int ncol = exmat.n_cols;
  int N    = data.size();
  int K    = label.max() + 1;
  arma::cube X(nrow,ncol,N,fill::zeros);
  for (int n=0; n<N; n++){
    X.slice(n) = Rcpp::as<arma::mat>(data[n]);
  }
  int Nw = cvi_helper_nw(label);
  arma::cube mean_class = cvi_helper_classmean(mfd, dtype, X, label);
  arma::field<arma::uvec> classindex = cvi_helper_classindex(label);
  
  // preliminary compute
  arma::vec vec_S(K,fill::zeros);
  arma::uvec tgtidx(2,fill::zeros);
  int ntgt = 0;
  double tmpsum = 0.0;
  for (int k=0; k<K; k++){
    tgtidx.reset();
    tgtidx = classindex(k);
    ntgt   = tgtidx.n_elem;
    
    tmpsum = 0.0;
    for (int i=0; i<ntgt; i++){
      tmpsum += cvi_helper_distance(mfd,dtype,mean_class.slice(k),X.slice(tgtidx(i)));
    }
    vec_S(k) = tmpsum/static_cast<double>(ntgt);
  }
  arma::mat mean_dists = cvi_helper_pdist(mfd, dtype, mean_class);
  mean_dists.diag().fill(arma::datum::inf);
  
  // main computation
  double output = 0.0;
  arma::rowvec row_scores(K,fill::zeros);
  for (int k=0; k<K; k++){
    row_scores.fill(0.0);
    for (int j=0; j<K; j++){
      if (j==k){
        row_scores(j) = 0.0;
      } else {
        row_scores(j) = (vec_S(k)+vec_S(j))/mean_dists(k,j);
      }
    }
    output += row_scores.max();
  }
  
  // return output
  return(output/static_cast<double>(K));
}



// (04; -) C Index -------------------------------------------------------------
// [[Rcpp::export]]
double cvi_internal_ci(std::string mfd, std::string dtype, Rcpp::List& data, arma::uvec mylabel){
  // PRELIMINARY
  arma::uvec label = mylabel - mylabel.min(); // starting from 0 to (k-1)
  arma::mat exmat  = Rcpp::as<arma::mat>(data[0]);
  int nrow = exmat.n_rows;
  int ncol = exmat.n_cols;
  int N    = data.size();
  int K    = label.max() + 1;
  arma::cube X(nrow,ncol,N,fill::zeros);
  for (int n=0; n<N; n++){
    X.slice(n) = Rcpp::as<arma::mat>(data[n]);
  }
  int Nw = cvi_helper_nw(label);
  int Nwhalf = static_cast<int>(std::ceil(static_cast<double>(Nw/2.0)));
  
  // COMPUTE : Prerequisites
  arma::mat distmat = cvi_helper_pdist(mfd, dtype, X);
  arma::vec disttri = arma::vectorise(arma::trimatu(distmat,1));
  arma::field<arma::uvec> classindex = cvi_helper_classindex(label);
  
  // COMPUTE SCORES
  double SC = 0.0;
  for (int k=0; k<K; k++){
    SC += arma::accu(distmat(classindex(k), classindex(k)));
  }
  arma::vec dtri_sorted = arma::sort(disttri, "ascend");
  double SminC = 2.0*arma::accu(dtri_sorted.head(Nwhalf));
  double SmaxC = 2.0*arma::accu(dtri_sorted.tail(Nwhalf));
  
  // RETURN OUTPUT
  double output = (SC-SminC)/(SmaxC-SminC);
  return(output);
}

// (03; -) Gamma ---------------------------------------------------------------
// (02; +) Calinski-Harabasz ---------------------------------------------------
// [[Rcpp::export]]
double cvi_internal_ch(std::string mfd, std::string dtype, Rcpp::List& data, arma::uvec mylabel){
  // PRELIMINARY
  arma::uvec label = mylabel - mylabel.min(); // starting from 0 to (k-1)
  arma::mat exmat  = Rcpp::as<arma::mat>(data[0]);
  int nrow = exmat.n_rows;
  int ncol = exmat.n_cols;
  int N    = data.size();
  int K    = label.max() + 1;
  arma::cube X(nrow,ncol,N,fill::zeros);
  for (int n=0; n<N; n++){
    X.slice(n) = Rcpp::as<arma::mat>(data[n]);
  }
  
  // COMPUTE : CLASS INDEX, CLASS MEAN, AND GLOBAL MEAN
  arma::mat distmat = cvi_helper_pdist(mfd, dtype, X);
  arma::field<arma::uvec> classindex = cvi_helper_classindex(label);
  arma::cube mean_class  = cvi_helper_classmean(mfd, dtype, X, label);
  arma::mat  mean_global = internal_mean(mfd, dtype, X, 100, 1e-5);
    
  // COMPUTE : DENOMINATOR
  double val_den = 0.0;
  for (int k=0; k<K; k++){
    val_den += (static_cast<double>(classindex(k).n_elem))*cvi_helper_distance(mfd, dtype, mean_class.slice(k), mean_global);
  }
  
  // COMPUTE : NUMERATOR
  arma::vec vec_numdist(K,fill::zeros);
  double tmpval  = 0.0; 
  int tmpsize = 0;
  arma::uvec tmpclass(2,fill::zeros);
  for (int k=0; k<K; k++){
    tmpclass.reset();
    tmpclass = classindex(k);
    tmpval   = 0.0;
    tmpsize  = tmpclass.n_elem;
    for (int i=0; i<tmpsize; i++){
      tmpval += cvi_helper_distance(mfd, dtype, X.slice(tmpclass(i)), mean_class.slice(k));
    }
    vec_numdist(k) = tmpval;
  }
  double val_num = arma::accu(vec_numdist);
  
  double output = val_den/val_num;
  return(output);  
}


// (01; +) Dunn index ----------------------------------------------------------
// [[Rcpp::export]]
double cvi_internal_dunn(std::string mfd, std::string dtype,  Rcpp::List& data, arma::uvec mylabel){
  // PRELIMINARY
  arma::uvec label = mylabel - mylabel.min(); // starting from 0 to (k-1)
  arma::mat exmat  = Rcpp::as<arma::mat>(data[0]);
  int nrow = exmat.n_rows;
  int ncol = exmat.n_cols;
  int N    = data.size();
  int K    = label.max() + 1;
  arma::cube X(nrow,ncol,N,fill::zeros);
  for (int n=0; n<N; n++){
    X.slice(n) = Rcpp::as<arma::mat>(data[n]);
  }
  
  // COMPUTE : DISTANCE MATRIX & CLASS INDEX
  arma::mat distmat = cvi_helper_pdist(mfd, dtype, X);
  arma::field<arma::uvec> classindex = cvi_helper_classindex(label);    

  // COMPUTE : DENOMINATOR
  arma::uvec ci(2,fill::zeros);
  arma::uvec cj(2,fill::zeros);
  arma::uvec ck(2,fill::zeros);
  
  arma::mat mat_den(K,K,fill::zeros); 
  mat_den.diag().fill(arma::datum::inf);
  for (int i=0; i<(K-1); i++){
    ci.reset();
    ci = classindex(i);
    for (int j=(i+1); j<K; j++){
      cj.reset();
      cj = classindex(j);
      mat_den(i,j) = distmat(ci, cj).min();
      mat_den(j,i) = mat_den(i,j);
    }
  }
  double val_den = mat_den.min();

  // COMPUTE : NUMERATOR
  arma::vec vec_num(K,fill::zeros);
  for (int k=0; k<K; k++){
    ck.reset();
    ck = classindex(k);
    vec_num(k) = distmat(ck,ck).max();
  }
  double val_num = vec_num.max();

  // RETURN
  double output = val_den/val_num;
  return(output);
}

