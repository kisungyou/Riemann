#include <RcppArmadillo.h>
#include "riemann_src.h"

using namespace Rcpp;
using namespace arma;
using namespace std;

// 1. visualize_pga    : (intrinsic) method by Fletcher
// 2. visualize_kpca   : kernel principal component analysis
// 3. visualize_isomap : weighted distance function
// 4. visualize_cmds   : classical muldimensional scaling
// 5. visualize_sammon : sammon mapping adapted from Rdimtools

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

// 4. visualize_cmds ===========================================================
arma::mat engine_cmds(arma::mat pdist, int ndim){ // given distance matrix, return (n x ndim)
  int N = pdist.n_rows;
  arma::mat D2 = arma::pow(pdist, 2.0);
  arma::mat J  = arma::eye<arma::mat>(N,N) - (arma::ones<arma::mat>(N,N)/(static_cast<double>(N)));
  arma::mat B  = -0.5*J*D2*J;
  
  arma::vec eigval;
  arma::mat eigvec;
  
  arma::eig_sym(eigval, eigvec, B);
  arma::mat Y = eigvec.tail_cols(ndim)*arma::diagmat(arma::sqrt(eigval.tail(ndim)));
  return(Y);
}
double engine_stress(arma::mat D, arma::mat Dhat){
  int N = D.n_rows;
  
  double tobesq = 0.0;
  double term1  = 0.0; // numerator
  double term2  = 0.0; // denominator
  for (int i=0;i<(N-1);i++){
    for (int j=(i+1);j<N;j++){
      tobesq = D(i,j)-Dhat(i,j);
      term1 += (tobesq*tobesq);
      term2 += D(i,j)*D(i,j);
    }
  }
  return(sqrt(term1/term2));  
}
// [[Rcpp::export]]
Rcpp::List visualize_cmds(std::string mfd, std::string geo, Rcpp::List& data, int ndim){
  // PREPARE
  arma::mat tmpdata = Rcpp::as<arma::mat>(data[0]);
  int N = data.size();
  int nrow = tmpdata.n_rows;
  int ncol = tmpdata.n_cols;
  
  arma::cube mydata(nrow,ncol,N,fill::zeros);
  for (int n=0; n<N; n++){
    mydata.slice(n) = Rcpp::as<arma::mat>(data[n]);
  }
  
  // COMPUTE PAIRWISE DISTANCE MATRIX
  arma::mat pdist(N,N,fill::zeros);
  for (int i=0; i<(N-1); i++){
    for (int j=(i+1); j<N; j++){
      if (geo=="intrinsic"){
        pdist(i,j) = riem_dist(mfd, mydata.slice(i), mydata.slice(j));
      } else {
        pdist(i,j) = riem_distext(mfd, mydata.slice(i), mydata.slice(j));
      }
      pdist(j,i) = pdist(i,j);
    }
  }
  
  // COMPUTE EMBEDDING & ITS DISTANCE
  arma::mat Y = engine_cmds(pdist, ndim);
  arma::mat DY(N,N,fill::zeros);
  for (int i=0; i<(N-1); i++){
    for (int j=(i+1); j<N; j++){
      DY(i,j) = arma::norm(Y.row(i)-Y.row(j),2);
      DY(j,i) = DY(i,j);
    }
  }
  double stress = engine_stress(pdist, DY);
  
  // WRAP AND RETURN
  Rcpp::List output;
  output["embed"] = Y;
  output["stress"] = stress;
  return(output);
}

// 5. visualize_sammon =========================================================
// [[Rcpp::export]]
Rcpp::List visualize_sammon(std::string mfd, std::string geo, Rcpp::List& data, int ndim, int maxiter, double abstol){
  // PREPARE
  arma::mat tmpdata = Rcpp::as<arma::mat>(data[0]);
  int N = data.size();
  int nrow = tmpdata.n_rows;
  int ncol = tmpdata.n_cols;
  
  arma::cube mydata(nrow,ncol,N,fill::zeros);
  for (int n=0; n<N; n++){
    mydata.slice(n) = Rcpp::as<arma::mat>(data[n]);
  }
  
  // COMPUTE PAIRWISE DISTANCE MATRIX
  arma::mat DPJX(N,N,fill::zeros);
  double c = 0.0;
  for (int i=0; i<(N-1); i++){
    for (int j=(i+1); j<N; j++){
      if (geo=="intrinsic"){
        DPJX(i,j) = riem_dist(mfd, mydata.slice(i), mydata.slice(j));
      } else {
        DPJX(i,j) = riem_distext(mfd, mydata.slice(i), mydata.slice(j));
      }
      DPJX(j,i) = DPJX(i,j);
      c += DPJX(i,j);
    }
  }
  
  // INITIALIZE VIA MDS
  arma::mat Yold = engine_cmds(DPJX, ndim);
  arma::mat Ynew(N,ndim,fill::zeros);
  double    Yinc = 0.0;
  arma::mat DPJY(N,N,fill::zeros);
  
  // MAIN ITERATION FOR SAMMON MAPPING
  double der1, der2, dpjstar, dpj;
  for (int it=0; it<maxiter; it++){
    // 1. compute DPJ for current iterate
    for (int i=0; i<(N-1); i++){
      for (int j=(i+1); j<N; j++){
        DPJY(i,j) = arma::norm(Yold.row(i)-Yold.row(j),2);
        DPJY(j,i) = DPJY(i,j);
      }
    }
    // 2. update for Ynew
    for (int p=0;p<N;p++){
      for (int q=0;q<ndim;q++){
        der1 = 0;
        der2 = 0;
        
        for (int j=0;j<N;j++){
          if (p!=j){
            dpjstar = DPJX(p,j);
            dpj     = DPJY(p,j);
            der1 += ((dpjstar-dpj)/(dpjstar*dpj))*(Yold(p,q)-Yold(j,q));
            der2 += (dpjstar-dpj)-(pow((Yold(p,q)-Yold(j,q)),2)/dpj)*(1+((dpjstar-dpj)/dpj));
          }
        }
        
        der1 *= (-2/c);
        der2 *= (-2/c);
        if (der2 < 0){
          der2 *= -1;
        }
        
        Ynew(p,q) = Yold(p,q) - 0.3*(der1/der2);
      }
    }
    // 3. updating information
    Yinc = arma::norm(Yold-Ynew,"fro");
    Yold = Ynew;
    if (Yinc < abstol){
      break;
    }
  }
  
  // COMPUTE THE STRESS
  for (int i=0; i<(N-1); i++){
    for (int j=(i+1); j<N; j++){
      DPJY(i,j) = arma::norm(Yold.row(i)-Yold.row(j),2);
      DPJY(j,i) = DPJY(i,j);
    }
  }
  double stress = engine_stress(DPJX, DPJY);
  
  // WRAP AND RETURN
  Rcpp::List output;
  output["embed"]  = Yold;
  output["stress"] = stress;
  return(output);
}
