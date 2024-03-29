#ifndef RIEMANN_MANIFOLDS_H
#define RIEMANN_MANIFOLDS_H

#define ARMA_NO_DEBUG

// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include "riemann_general.h"

using namespace Rcpp;
using namespace arma;
using namespace std;

// OPERATION
// (mat) initialize
// (mat) exp
// (mat) proj
// (mat) log
// (double) dist
// (double) metric
// (double) distext
// (vec) equiv
// (mat) invequiv


// 01. SPHERE ==================================================================
arma::mat sphere_initialize(arma::field<arma::mat> data, arma::vec weight){
  int N       = data.n_elem;
  double wsum = arma::accu(weight);
  
  arma::mat outmat(data(0).n_rows, data(0).n_cols, fill::zeros);
  for (int n=0; n<N; n++){
    outmat += (weight(n)/wsum)*data(n);
  }
  outmat /= arma::norm(outmat, "fro");
  return(outmat);
}
arma::mat sphere_exp(arma::mat x, arma::mat d, double t){
  double nrm_td = arma::norm(t*d, "fro"); // theta
  arma::mat out;
  if (nrm_td < 1e-15){ // very close
    out = x;
  } else {
    out = cos(nrm_td)*x + ((sin(nrm_td))/nrm_td)*t*d;
    out /= arma::norm(out, "fro");
  }
  return(out);
}
arma::mat sphere_proj(arma::mat x, arma::mat u){
  return(u-x*(arma::dot(x,u)));
}
double sphere_dist(arma::mat x, arma::mat y){
  arma::vec vecx = arma::vectorise(x);
  arma::vec vecy = arma::vectorise(y);
  arma::vec vecxy = vecx-vecy;
  double dotxy = arma::dot(vecx, vecy);
  
  if (arma::norm(vecxy, 2) < arma::datum::eps){
    return(0.0);
  } else if (std::sqrt(dotxy*dotxy) >= (1.0-arma::datum::eps)){
    return(arma::datum::pi);
  } else {
    return(std::acos(arma::dot(vecx, vecy)));  
  }
}
arma::mat sphere_log(arma::mat x, arma::mat y){
  arma::mat v = sphere_proj(x,y-x);
  double di = sphere_dist(x,y);
  if (di > 1e-6){
    double nv = arma::norm(v, "fro");
    v = v*(di/nv);
  }
  return(v);
}
double sphere_metric(arma::mat x, arma::mat d1, arma::mat d2){
  return(arma::dot(arma::vectorise(d1),arma::vectorise(d2)));
}
double sphere_distext(arma::mat x, arma::mat y){
  return(arma::norm(arma::vectorise(x)-arma::vectorise(y),2));
}
arma::vec sphere_equiv(arma::mat x, int m, int n){
  arma::vec out = arma::vectorise(x,0);
  return(out);
}
arma::mat sphere_invequiv(arma::vec x, int m, int n){
  arma::mat out  = arma::reshape(x,m,n);
  double outsize = arma::norm(out,"fro");
  return((out/outsize));
}

// 02. SPD =====================================================================
arma::mat spd_initialize(arma::field<arma::mat> data, arma::vec weight){
  int N       = data.n_elem;
  double wsum = arma::accu(weight);

  arma::mat outmat(data(0).n_rows, data(0).n_cols, fill::zeros);
  for (int n=0; n<N; n++){
    outmat += (weight(n)/wsum)*data(n);
  }
  return(outmat);
}
arma::mat spdaux_symm(arma::mat x){
  return((x+x.t())/2.0);
}
arma::mat spd_exp(arma::mat x, arma::mat eta, double t){
  arma::mat teta = t*eta;
  arma::mat tmp  = arma::expmat(arma::solve(x, teta));
  arma::mat yy   = x*tmp;
  return((yy+yy.t())/2.0);
}
double spd_dist(arma::mat x, arma::mat y){
  arma::mat sol = arma::solve(x,y);
  arma::cx_mat cxXY = arma::logmat(sol);
  arma::mat XY = arma::real(cxXY);
  return(std::sqrt(arma::as_scalar(arma::trace(XY*XY))));
}
arma::mat spd_proj(arma::mat x, arma::mat u){
  return(spdaux_symm(u));
}
arma::mat spd_log(arma::mat x, arma::mat y){
  arma::mat tmp  = arma::real(arma::logmat(arma::solve(x,y)));
  arma::mat yy   = x*tmp;
  return((yy+yy.t())/2.0);
}
double spd_metric(arma::mat x, arma::mat u, arma::mat v){
  arma::mat sol1 = arma::solve(x,u);
  arma::mat sol2 = arma::solve(x,v);
  // use of trinner
  return(arma::as_scalar(arma::trace(sol1.t()*sol2)));
}
arma::vec spd_equiv(arma::mat x, int m, int n){
  arma::cx_mat logx = arma::logmat(x);
  arma::mat realogx = arma::real(logx);
  arma::vec out = arma::vectorise(realogx,0);
  return(out);
}
arma::mat spd_invequiv(arma::vec x, int m, int n){
  arma::mat tmpx = arma::reshape(x,m,n);
  arma::mat output = arma::expmat(((tmpx + tmpx.t())/2.0));
  return(output);
}
double spd_distext(arma::mat x, arma::mat y){
  int m = x.n_rows;
  int n = x.n_cols;
  
  arma::vec xext = spd_equiv(x, m, n);
  arma::vec yext = spd_equiv(y, m, n);
  
  return(arma::as_scalar(arma::norm(xext-yext,"fro")));
}

// 03. CORRELATION =============================================================
// NOTE : find geodesic minimizing D for d(A, DBD) with some default parameters
arma::mat corr_airm_findD(arma::mat C1, arma::mat C2){
  // parameters and prepare
  int p = C1.n_rows;
  int maxiter    = 50;  
  double stopeps = 1e-10;
  
  arma::mat Dold(p,p,fill::eye);
  arma::mat Dnew(p,p,fill::zeros);
  arma::mat Dtmp(p,p,fill::zeros);
  arma::mat Ctmp(p,p,fill::zeros);
  
  arma::mat Dhalf(p,p,fill::zeros);
  arma::mat Dhalfinv(p,p,fill::zeros);
  
  arma::mat Ddel1(p,p,fill::zeros);
  arma::mat Ddel2(p,p,fill::zeros);
  arma::mat Ddel3(p,p,fill::zeros);
  
  double costold = spd_dist(C1,C2);
  double costnew = 0.0;
  double Dincval = 0.0;
  
  arma::vec vec_delta = arma::linspace(0.0001, 1.9999, 20);
  arma::vec vec_cost(20,fill::zeros);
  
  // main iteration
  for (int it=0; it<maxiter; it++){
    // compute intermediate values
    Ddel1 = Dold*arma::real(arma::logmat(C2*Dold*arma::solve(C1, Dold)));
    Ddel2 = 2.0*mat_symm(Ddel1, true);
    
    Dhalf    = mat_diaghalf(Dold);
    Dhalfinv = mat_diaginvhalf(Dold);
    Ddel3    = Dhalfinv*Ddel2*Dhalfinv;
    
    // iterate over the values
    for (int i=0;i<20;i++){
      Dtmp = Dhalf*arma::expmat_sym(-vec_delta(i)*Ddel3)*Dhalf;
      Ctmp = Dtmp*C2*Dtmp;
      vec_cost(i) = spd_dist(C1, Ctmp);
    }
    
    // optimal one
    arma::uword minid = arma::index_min(vec_cost);
    Dnew    = Dhalf*arma::expmat_sym(-vec_delta(minid)*Ddel3)*Dhalf;
    costnew = vec_cost(minid);
    
    // updating information
    Dincval = arma::norm(Dold - Dnew, "fro")/arma::norm(Dold, "fro");
    Dold    = Dnew;
    if ((it>0)&&((costnew > costold)||(Dincval < stopeps))){
      break;
    }
    costold = costnew;
  }
  return(Dold);
}
arma::mat correlation_initialize(arma::field<arma::mat> data, arma::vec weight){
  int N       = data.n_elem;
  double wsum = arma::accu(weight);

  arma::mat outmat(data(0).n_rows, data(0).n_cols, fill::zeros);
  for (int n=0; n<N; n++){
    outmat += (weight(n)/wsum)*data(n);
  }
  return(outmat);
}
double correlation_dist(arma::mat A, arma::mat B){
  arma::mat D = corr_airm_findD(A,B);
  arma::mat C = D*B*D;
  double dval = spd_dist(A,C);
  return(dval);
}
arma::mat correlation_exp(arma::mat X, arma::mat eta, double t){
  arma::mat Ztmp = spd_exp(X, eta, t);
  arma::mat Zout = mat_cov2cor(Ztmp);
  return(Zout);
}
arma::mat correlation_log(arma::mat X, arma::mat Y){
  int p = X.n_rows;
  arma::mat D = corr_airm_findD(X, Y); // diagonal scalar
  arma::mat Z = D*Y*D;                 // optimal position
  arma::mat logXY = spd_log(X, Z);
  return(logXY);
}
double correlation_metric(arma::mat X, arma::mat eta1, arma::mat eta2){
  double outval = spd_metric(X, eta1, eta2);
  return(outval);
}

// 04. STIEFEL =================================================================
double stiefel_metric(arma::mat x, arma::mat d1, arma::mat d2){
  return(arma::as_scalar(arma::dot(arma::vectorise(d1), arma::vectorise(d2))));
}
arma::mat stiefel_proj(arma::mat x, arma::mat u){
  arma::mat A = (x.t()*u);
  return(u - x*((A+A.t())/2.0));
}
arma::mat stiefel_nearest(arma::mat x){
  arma::mat xtx = (x.t()*x);
  arma::mat output = x*arma::real(arma::powmat(xtx, -0.5));
  return(output);
}
arma::mat stiefel_exp(arma::mat x, arma::mat u, double t){
  const int n = x.n_rows;
  const int p = x.n_cols;
  
  arma::mat Ip(p,p,fill::eye);
  arma::mat Zp(p,p,fill::zeros);
  
  arma::mat tu = t*u;
  arma::mat term1 = arma::join_horiz(x, tu);
  
  arma::mat term21 = arma::join_horiz((x.t()*tu), -((tu.t())*tu));
  arma::mat term22 = arma::join_horiz(Ip, (x.t()*tu));
  arma::mat term2  = arma::expmat(arma::join_vert(term21, term22));
  
  arma::mat term3  = arma::join_vert(arma::expmat(-(x.t()*tu)), Zp);
  
  arma::mat output = term1*term2*term3;
  return(output);
}
arma::mat stiefel_log(arma::mat U0, arma::mat U1){
  const int n = U0.n_rows;
  const int p = U0.n_cols;
  const double tau = 1e-6;   // default convergence threshold
  const int maxiter = 12345; // default maximum number of iterations
  
  // 1.
  arma::mat M = U0.t()*U1;   
  // 2. thin QR of normal component of U1
  arma::mat Q,N;
  arma::qr_econ(Q,N,U1-(U0*M));
  // 3. orthogonal completion + procrustes preprocessing ------------- no QR_ECON ?
  arma::mat V, Vaway;
  arma::qr(V, Vaway, arma::join_vert(M,N)); 
  
  arma::mat D, R; 
  arma::vec vecS;
  arma::svd(D,vecS,R,V.submat(p,p,(2*p)-1,(2*p)-1));
  arma::mat S = arma::diagmat(vecS);
  V.cols(p,(2*p)-1) = V.cols(p,(2*p)-1)*(R*D.t());
  V = arma::join_horiz(arma::join_vert(M,N),V.cols(p,(2*p)-1)); 
  
  // 4. for looping
  arma::cx_mat LVcx;
  arma::mat LV, C, Phi;
  double normC;
  for (int k=0;k<maxiter;k++){
    LVcx = arma::logmat(V);
    LV   = arma::real(LVcx);
    
    // lower (pxp) diagonal block
    C = LV.submat(p,p,(2*p)-1,(2*p)-1);
    // convergence check
    normC = arma::norm(C, 2);
    if (normC < tau){
      break;
    }
    // matrix exponential
    Phi = arma::expmat(-C);
    // update last p columns
    V.cols(p,(2*p)-1) = V.cols(p,(2*p)-1)*Phi;
  }
  
  // 5. prepare output
  arma::mat Delta = (U0*LV.submat(0,0,(p-1),(p-1))) + (Q*LV.submat(p,0,(2*p)-1,p-1));
  return(Delta);
}
double stiefel_dist(arma::mat x, arma::mat y){
  arma::mat delta = stiefel_log(x,y);
  double output = std::sqrt(stiefel_metric(x,delta,delta));
  return(output);
}
arma::vec stiefel_equiv(arma::mat x, int m, int n){
  arma::vec output = arma::vectorise(x,0);
  return(output);
}
arma::mat stiefel_invequiv(arma::vec x, int m, int n){
  arma::mat mu = arma::reshape(x,m,n);
  arma::mat rhs = arma::pinv(arma::real(arma::sqrtmat(mu.t()*mu)));
  arma::mat output = mu*rhs;
  return(output);
}
double stiefel_distext(arma::mat x, arma::mat y){
  int m = x.n_rows;
  int n = x.n_cols;
  
  arma::vec xext = stiefel_equiv(x, m, n);
  arma::vec yext = stiefel_equiv(y, m, n);
  
  return(arma::as_scalar(arma::norm(xext-yext,"fro")));
}
arma::mat stiefel_initialize(arma::field<arma::mat> data, arma::vec weight){
  int N = data.n_elem;
  double wsum = arma::accu(weight);
  
  arma::mat outmat(data(0).n_rows, data(0).n_cols, fill::zeros);
  for (int n=0; n<N; n++){
    outmat += (weight(n)/wsum)*data(n);
  }
  arma::mat output = stiefel_nearest(outmat);
  return(output);
}

// 05. GRASSMANN ===============================================================
double grassmann_metric(arma::mat x, arma::mat d1, arma::mat d2){
  return(arma::as_scalar(arma::dot(arma::vectorise(d1), arma::vectorise(d2))));
}
double grassmann_dist(arma::mat X, arma::mat Y){
  arma::mat XY = X.t()*Y;
  arma::vec s  = arma::svd(XY);
  const int N  = s.n_elem;
  
  arma::vec theta(N,fill::zeros);
  for (int i=0;i<s.n_elem;i++){
    if (s(i) > 1){
      s(i) = 1.0;
    }
    theta(i) = std::acos(static_cast<float>(s(i)));
  }
  
  double output = 0.0;
  for (int i=0;i<s.n_elem;i++){
    output += theta(i)*theta(i);
  }
  return(std::sqrt(output));
}
arma::mat grassmann_proj(arma::mat x, arma::mat u){
  return(u-x*((x.t()*u)));
}
arma::mat grassmann_nearest(arma::mat x){
  return(stiefel_nearest(x));
}
arma::mat grassmann_exp(arma::mat x, arma::mat d, double t){
  const int n = x.n_rows;
  const int p = x.n_cols;
  
  arma::mat u, v, sin_s, cos_s;
  arma::vec s;
  
  arma::mat tu = t*d;
  arma::svd_econ(u,s,v,tu);
  cos_s = arma::diagmat(arma::cos(s));
  sin_s = arma::diagmat(arma::sin(s));
  
  arma::mat Y = x*v*cos_s*v.t() + u*sin_s*v.t();
  arma::mat Q,R;
  arma::qr_econ(Q,R,Y);
  return(Q);
}
arma::mat grassmann_log(arma::mat x, arma::mat y){
  int n = x.n_rows;
  int p = x.n_cols;
  
  arma::mat ytx = y.t()*x;             // (p,p)
  arma::mat At = (y.t())-(ytx*x.t());  // (p,n)
  arma::mat Bt = arma::solve(ytx, At); // (p,n)
  
  arma::mat u,v;
  arma::vec s;
  
  arma::svd_econ(u,s,v,(Bt.t())); // (n,p) svd : (n,n) (n,p) (p,p) for [u,s,v]
  
  arma::mat U = u.cols(0,(p-1));
  arma::mat S = arma::diagmat(arma::atan(s));
  arma::mat V = v.cols(0,(p-1));
  
  arma::mat output = U*S*V.t();
  return(output);
}
arma::vec grassmann_equiv(arma::mat x, int n, int p){
  arma::vec output = arma::vectorise((x*x.t()),0);
  return(output);
}
arma::mat grassmann_invequiv(arma::vec x, int n, int p){
  arma::mat tmpx = arma::reshape(x,n,n);  // (n x n) in equivariants
  arma::mat symx = (tmpx + tmpx.t())/2.0; // symmetrization for sake
  
  arma::vec eigval;
  arma::mat eigvec;
  eig_sym(eigval, eigvec, symx); // in an ascending order
  
  arma::mat output = arma::fliplr(eigvec.tail_cols(p)); // yes.
  return(output);
}
double grassmann_distext(arma::mat x, arma::mat y){
  int m = x.n_rows;
  int n = x.n_cols;
  
  arma::vec xext = grassmann_equiv(x, m, n);
  arma::vec yext = grassmann_equiv(y, m, n);
  
  return(arma::norm(xext-yext,2));
}
arma::mat grassmann_initialize(arma::field<arma::mat> data, arma::vec weight){
  int N       = data.n_elem;
  double wsum = arma::accu(weight);
  
  arma::mat outmat(data(0).n_rows, data(0).n_cols, fill::zeros);
  for (int n=0; n<N; n++){
    outmat += (weight(n)/wsum)*data(n);
  }
  return(grassmann_nearest(outmat));  
}


// 06. ROTATION ================================================================
// https://www.manopt.org/reference/manopt/manifolds/rotations/rotationsfactory.html
arma::mat rotation_initialize(arma::field<arma::mat> data, arma::vec weight){
  int N       = data.n_elem;
  double wsum = arma::accu(weight);
  
  arma::mat outmat(data(0).n_rows, data(0).n_cols, fill::zeros);
  for (int n=0; n<N; n++){
    outmat += (weight(n)/wsum)*data(n);
  }
  arma::mat Q;
  arma::mat R;
  arma::qr(Q,R,outmat);
  return(Q);
}
double rotation_metric(arma::mat X, arma::mat eta1, arma::mat eta2){
  return(arma::dot(arma::vectorise(eta1), arma::vectorise(eta2)));
}
arma::mat rotation_log(arma::mat X, arma::mat Y){
  arma::mat Utmp = arma::real(arma::logmat(arma::trans(X)*Y));
  arma::mat Uout = (Utmp - arma::trans(Utmp))/2.0; // skew-symmetric part
  return(Uout);
}
arma::mat rotation_exp(arma::mat X, arma::mat eta, double t){
  arma::mat exptU  = arma::real(arma::expmat(t*eta));
  arma::mat output = X*exptU;
  return(output);
}
double rotation_dist(arma::mat X, arma::mat Y){
  arma::mat logXY = rotation_log(X, Y);
  double output   = arma::norm(arma::vectorise(logXY), 2);
  return(output);
}
arma::vec rotation_equiv(arma::mat x, int m, int n){
  arma::vec output = arma::vectorise(x,0);
  return(output);
}
arma::mat rotation_nearest(arma::mat x){
  return(stiefel_nearest(x));
}
arma::mat rotation_invequiv(arma::vec x, int m, int n){
  arma::mat resmat = arma::reshape(x, m, n);
  arma::mat identity(m,m,fill::eye);
  if (arma::norm((resmat.t()*resmat - identity), "fro")/std::sqrt(static_cast<double>(m)) > 1e-10){
    return(rotation_nearest(resmat));
  } else {
    return(resmat);
  }
}
double rotation_distext(arma::mat x, arma::mat y){
  int m = x.n_rows;
  int n = x.n_cols;
  
  arma::vec xext = rotation_equiv(x,m,n);
  arma::vec yext = rotation_equiv(y,m,n);
  
  return(arma::norm(xext-yext,2));
}

// 07. MULTINOMIAL =============================================================
// Astrom (2017) "Image Labeling by Assignment"
// https://www.manopt.org/reference/manopt/manifolds/multinomial/multinomialfactory.html
arma::mat multinomial_initialize(arma::field<arma::mat> data, arma::vec weight){
  int N       = data.n_elem;
  double wsum = arma::accu(weight);
  
  arma::mat outmat(data(0).n_rows, data(0).n_cols, fill::zeros);
  for (int n=0; n<N; n++){
    outmat += (weight(n)/wsum)*data(n);
  }
  outmat /= arma::accu(arma::abs(outmat));
  return(outmat);
}
double multinomial_metric(arma::mat x, arma::mat d1, arma::mat d2){
  arma::vec X    = arma::vectorise(x);
  arma::vec eta  = arma::vectorise(d1);
  arma::vec zeta = arma::vectorise(d2);
  double output = arma::accu((eta%zeta)/X);
  return(output); 
}
arma::mat multinomial_log(arma::mat x, arma::mat y){
  arma::mat a = arma::sqrt(x%y);
  double    s = arma::accu(a);
  arma::mat u = (2.0*std::acos(s)/std::sqrt(1.0-(s*s)))*(a - (s*x));
  return(u);
}
arma::mat multinomial_exp(arma::mat x, arma::mat u, double t){
  int n = x.n_rows;
  arma::mat y(n,1,fill::zeros);
  arma::mat z(n,1,fill::zeros);
  arma::mat tU = t*u;
  arma::mat s  = arma::sqrt(x);
  arma::mat us = (tU/s)/2.0;
  double    un = arma::norm(us, "fro");
  if (un < arma::datum::eps){
    y = x;
  } else {
    z = (std::cos(un)*s + (std::sin(un)/un)*us);
    y = z%z;
  }
  arma::mat output = y/arma::accu(y);
  return(output);
}
double multinomial_dist(arma::mat x, arma::mat y){
  arma::mat logxy = multinomial_log(x, y);
  double output   = std::sqrt(multinomial_metric(x, logxy, logxy));
  return(output);
}
arma::vec multinomial_equiv(arma::mat x, int m, int n){ // s = 2*sqrt(p)
  arma::vec out = arma::sqrt(arma::vectorise(x,0))*2.0;
  return(out);
}

arma::mat multinomial_invequiv(arma::vec x, int m, int n){ // p = (s^2)/4
  arma::mat out  = arma::reshape(x,m,n);
  arma::mat out2 = out%out;
  arma::mat out3 = out2/arma::accu(arma::abs(out2));
  return(out3);
}
double multinomial_distext(arma::mat x, arma::mat y){
  int mm = x.n_rows;
  int nn = x.n_cols;
  arma::vec xx = multinomial_equiv(x, mm, nn);
  arma::vec yy = multinomial_equiv(y, mm, nn);
  return(arma::norm(xx-yy,2));
}

// 08. SPD-K : Fixed-Rank K ====================================================
// https://www.manopt.org/reference/manopt/manifolds/symfixedrank/symfixedrankYYfactory.html
arma::mat spdk_initialize(arma::field<arma::mat> data, arma::vec weight){
  int N       = data.n_elem;
  double wsum = arma::accu(weight);
  
  arma::mat outmat(data(0).n_rows, data(0).n_cols, fill::zeros);
  for (int n=0; n<N; n++){
    outmat += (weight(n)/wsum)*data(n);
  }
  return(outmat);
}
double spdk_metric(arma::mat X, arma::mat eta, arma::mat zeta){
  return(arma::dot(arma::vectorise(eta), arma::vectorise(zeta)));
}
arma::mat spdk_log(arma::mat Y, arma::mat Z){
  arma::mat YtZ = Y.t()*Z;
  arma::mat U;  arma::vec s;  arma::mat V;
  arma::svd(U,s,V,YtZ);
  arma::mat Qt = V*U.t();
  arma::mat eta = Z*Qt - Y;
  return(eta);
}
arma::mat spdk_exp(arma::mat Y, arma::mat eta, double t){
  arma::mat Ynew = t*eta;
  return(Ynew);
}
double spdk_dist(arma::mat X, arma::mat Y){
  arma::mat logXY = spdk_log(X, Y);
  double distval  = std::sqrt(spdk_metric(X, logXY, logXY));
  return(distval);
}

// 10. EUCLIDEAN ===============================================================
arma::mat euclidean_initialize(arma::field<arma::mat> data, arma::vec weight){
  int N       = data.n_elem;
  double wsum = arma::accu(weight);
  
  arma::mat outmat(data(0).n_rows, data(0).n_cols, fill::zeros);
  for (int n=0; n<N; n++){
    outmat += (weight(n)/wsum)*data(n);
  }
  return(outmat);
}
arma::mat euclidean_exp(arma::mat x, arma::mat d, double t){
  arma::mat y = x + t*d;
  return(y);
}
arma::mat euclidean_log(arma::mat x, arma::mat y){
  return(y-x);
}
double euclidean_metric(arma::mat x, arma::mat d1, arma::mat d2){
  return(arma::dot(arma::vectorise(d1), arma::vectorise(d2)));
}
double euclidean_dist(arma::mat x, arma::mat y){
  return(arma::norm(x-y,"fro"));
}
arma::mat euclidean_proj(arma::mat x, arma::mat u){
  return(u);
}
double euclidean_distext(arma::mat x, arma::mat y){
  return(arma::norm(x-y,"fro"));
}
arma::vec euclidean_equiv(arma::mat x, int m, int n){
  arma::vec out = arma::vectorise(x,0);
  return(out);
}
arma::mat euclidean_invequiv(arma::vec x, int m, int n){
  arma::mat out = arma::reshape(x,m,n);
  return(out);
}

// 12. LANDMARK ================================================================
arma::mat landmark_aux_nearest(arma::mat x){ // auxiliary for landmark : nearest
  int n = x.n_rows;
  int p = x.n_cols;
  
  arma::rowvec xvec = arma::mean(x, 0);
  arma::mat y(n,p,fill::zeros);
  for (int i=0; i<n; i++){
    y.row(i) = x.row(i)-xvec;
  }
  return(y/arma::norm(y,"fro"));
}
arma::mat landmark_aux_matching(arma::mat x, arma::mat y){ // aux : matching
  arma::mat xy = arma::trans(x)*y;                         // OPS problem 
  
  arma::mat U; 
  arma::vec s;
  arma::mat V;
  arma::svd(U,s,V,xy);
  
  arma::mat output = y*V*arma::trans(U);
  return(output);
}
arma::mat landmark_initialize(arma::field<arma::mat> data, arma::vec weight){
  int N       = data.n_elem;
  double wsum = arma::accu(weight);
  
  arma::mat outmat(data(0).n_rows, data(0).n_cols, fill::zeros);
  for (int n=0; n<N; n++){
    outmat += (weight(n)/wsum)*data(n);
  }
  arma::mat output = landmark_aux_nearest(outmat);
  return(output);
}
arma::mat landmark_exp(arma::mat x, arma::mat d, double t){
  int n = x.n_rows;
  int p = x.n_cols;
  
  arma::mat xx = arma::reshape(x, n*p, 1);
  arma::mat dd = arma::reshape(d, n*p, 1);
  arma::mat zz = sphere_exp(xx, dd, t);
  arma::mat ZZ = arma::reshape(zz, n, p);
  arma::mat output = landmark_aux_nearest(ZZ);
  return(output);
}
arma::mat landmark_log(arma::mat X, arma::mat Y){
  int n = X.n_rows;
  int p = X.n_cols;
  arma::mat YY = landmark_aux_matching(X,Y);
  arma::mat xx = arma::reshape(X, n*p, 1);  // resize for sphere
  arma::mat yy = arma::reshape(YY, n*p, 1);
  
  arma::mat zz = sphere_log(xx, yy);
  arma::mat ZZ = arma::reshape(zz, n, p);
  return(ZZ);
}
double landmark_metric(arma::mat x, arma::mat d1, arma::mat d2){
  return(arma::dot(arma::vectorise(d1),arma::vectorise(d2)));  
}
double landmark_dist(arma::mat x, arma::mat y){
  arma::mat logxy = landmark_log(x, y);
  double output = std::sqrt(landmark_metric(x, logxy, logxy));
  return(output);
}
double landmark_distext(arma::mat x, arma::mat y){
  arma::mat yy  = landmark_aux_matching(x,y);
  double output = arma::norm(x-yy, "fro");
  return(output);
}
arma::vec landmark_equiv(arma::mat x, int m, int n){
  arma::vec out = arma::vectorise(x,0);
  return(out);
}
arma::mat landmark_invequiv(arma::vec x, int m, int n){
  arma::mat out = landmark_aux_nearest(arma::reshape(x,m,n));
  return(out);
}

#endif

