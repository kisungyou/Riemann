#ifndef RIEMANN_MANIFOLDS_H
#define RIEMANN_MANIFOLDS_H


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
arma::mat sphere_unit(arma::mat x){
  double nrm = arma::norm(x, "fro");
  if (!x.is_finite() || !std::isfinite(nrm) || nrm == 0.0 || std::abs(nrm-1.0) > 1e-8)
    Rcpp::stop("Sphere operation requires a finite unit-norm point.");
  return x/nrm;
}
arma::mat sphere_initialize(arma::field<arma::mat> data, arma::vec weight){
  arma::mat out(data(0).n_rows, data(0).n_cols, fill::zeros);
  for (arma::uword n=0; n<data.n_elem; ++n) if (weight(n)>0) out += weight(n)*data(n);
  double nrm = arma::norm(out, "fro");
  // This is an intrinsic optimizer initializer, not a unique extrinsic mean.
  if (nrm <= 64.0*arma::datum::eps) {
    for (arma::uword n=0; n<data.n_elem; ++n) if (weight(n)>0) return sphere_unit(data(n));
    Rcpp::stop("Sphere initialization needs positive weight.");
  }
  return sphere_unit(out/nrm);
}
arma::mat sphere_exp(arma::mat x, arma::mat d, double t){
  x = sphere_unit(x);
  if (!d.is_finite() || !std::isfinite(t)) Rcpp::stop("Nonfinite sphere tangent or step.");
  double radial = arma::accu(x%d);
  if (std::abs(radial)>1e-8*std::max(1.0,arma::norm(d,"fro")))
    Rcpp::stop("Sphere exponential requires a tangent vector.");
  d -= radial*x;
  double nrm = arma::norm(t*d, "fro");
  if (!std::isfinite(nrm)) Rcpp::stop("Sphere step exceeds numerical range.");
  if (nrm == 0.0) return x;
  arma::mat out = std::cos(nrm)*x + (std::sin(nrm)/nrm)*t*d;
  return out/arma::norm(out,"fro");
}
arma::mat sphere_proj(arma::mat x, arma::mat u){
  return u-x*arma::accu(x%u);
}
double sphere_dist(arma::mat x, arma::mat y){
  x = sphere_unit(x); y = sphere_unit(y);
  return 2.0*std::atan2(arma::norm(x-y,"fro"),arma::norm(x+y,"fro"));
}
arma::mat sphere_log(arma::mat x, arma::mat y){
  x = sphere_unit(x); y = sphere_unit(y);
  double difference = arma::norm(x-y,"fro");
  if (difference == 0.0) return arma::zeros<arma::mat>(x.n_rows,x.n_cols);
  double sum_norm = arma::norm(x+y,"fro");
  if (sum_norm <= 32.0*arma::datum::eps)
    Rcpp::stop("Sphere logarithm is nonunique at antipodal points.");
  arma::mat v = sphere_proj(x,y-x);
  double nv = arma::norm(v,"fro");
  if (nv == 0.0) Rcpp::stop("Sphere logarithm direction cannot be resolved numerically.");
  return v*(2.0*std::atan2(difference,sum_norm)/nv);
}
double sphere_metric(arma::mat x, arma::mat d1, arma::mat d2){
  return arma::accu(d1%d2);
}
double sphere_distext(arma::mat x, arma::mat y){
  return arma::norm(x-y,"fro");
}
arma::vec sphere_equiv(arma::mat x, int m, int n){
  return arma::vectorise(x);
}
arma::mat sphere_invequiv(arma::vec x, int m, int n){
  arma::mat out = arma::reshape(x,m,n);
  double nrm = arma::norm(out,"fro");
  if (!std::isfinite(nrm) || nrm <= 64.0*arma::datum::eps)
    Rcpp::stop("The projected sphere summary is nonunique for a zero resultant.");
  return out/nrm;
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
  return x/2.0+x.t()/2.0;
}
void spdaux_eigen(const arma::mat& x, arma::vec& values, arma::mat& vectors){
  if (!x.is_finite() || x.n_rows == 0 || x.n_rows != x.n_cols)
    Rcpp::stop("SPD operation requires a finite square matrix.");
  double scale = arma::abs(x).max();
  if (scale == 0.0 || arma::abs(x-x.t()).max() > 128.0*arma::datum::eps*scale)
    Rcpp::stop("SPD operation requires a symmetric positive-definite matrix.");
  if (!arma::eig_sym(values,vectors,spdaux_symm(x/scale)) || values.min() <= 0.0)
    Rcpp::stop("SPD matrix is not numerically positive definite.");
  values *= scale;
}
arma::mat spdaux_power(const arma::mat& x, double power){
  arma::vec values; arma::mat vectors;
  spdaux_eigen(x,values,vectors);
  arma::mat out = vectors*arma::diagmat(arma::pow(values,power))*vectors.t();
  if (!out.is_finite()) Rcpp::stop("SPD matrix power exceeds numerical range.");
  return spdaux_symm(out);
}
arma::mat spdaux_log(const arma::mat& x){
  arma::vec values; arma::mat vectors;
  spdaux_eigen(x,values,vectors);
  return spdaux_symm(vectors*arma::diagmat(arma::log(values))*vectors.t());
}
arma::mat spdaux_exp(const arma::mat& x){
  arma::vec values; arma::mat vectors;
  if (!x.is_finite() || !arma::eig_sym(values,vectors,spdaux_symm(x)))
    Rcpp::stop("Symmetric matrix exponential failed.");
  arma::vec ev = arma::exp(values);
  if (!ev.is_finite() || ev.min() <= 0.0)
    Rcpp::stop("SPD exponential exceeds representable positive-definite range.");
  return spdaux_symm(vectors*arma::diagmat(ev)*vectors.t());
}
arma::mat spd_exp(arma::mat x, arma::mat eta, double t){
  arma::mat half = spdaux_power(x,0.5);
  arma::mat invhalf = spdaux_power(x,-0.5);
  arma::mat out = half*spdaux_exp(invhalf*(t*eta)*invhalf)*half;
  if (!out.is_finite()) Rcpp::stop("SPD exponential produced nonfinite values.");
  return spdaux_symm(out);
}
double spd_dist(arma::mat x, arma::mat y){
  arma::mat invhalf = spdaux_power(x,-0.5);
  return arma::norm(spdaux_log(spdaux_symm(invhalf*y*invhalf)),"fro");
}
arma::mat spd_proj(arma::mat x, arma::mat u){
  return spdaux_symm(u);
}
arma::mat spd_log(arma::mat x, arma::mat y){
  arma::mat half = spdaux_power(x,0.5);
  arma::mat invhalf = spdaux_power(x,-0.5);
  return spdaux_symm(half*spdaux_log(spdaux_symm(invhalf*y*invhalf))*half);
}
double spd_metric(arma::mat x, arma::mat u, arma::mat v){
  arma::mat invhalf = spdaux_power(x,-0.5);
  arma::mat a = invhalf*u*invhalf;
  arma::mat b = invhalf*v*invhalf;
  return arma::accu(a%b);
}
arma::vec spd_equiv(arma::mat x, int m, int n){
  return arma::vectorise(spdaux_log(x));
}
arma::mat spd_invequiv(arma::vec x, int m, int n){
  return spdaux_exp(arma::reshape(x,m,n));
}
double spd_distext(arma::mat x, arma::mat y){
  return arma::norm(spdaux_log(x)-spdaux_log(y),"fro");
}

// Shared validation for legacy matrix geometries. These checks deliberately
// reject ambiguous inverse projections instead of manufacturing an ONB.
void legacy_orthonormal(const arma::mat& x, const char* name, bool square=false){
  if (!x.is_finite() || x.n_cols == 0 || x.n_rows < x.n_cols ||
      (square && x.n_rows != x.n_cols) ||
      arma::norm(x.t()*x-arma::eye<arma::mat>(x.n_cols,x.n_cols),"fro") >
        1e-8*std::sqrt(static_cast<double>(x.n_cols)))
    Rcpp::stop("%s requires finite orthonormal matrices of compatible dimensions.",name);
}
void legacy_pair(const arma::mat& x,const arma::mat& y,const char* name){
  if (x.n_rows != y.n_rows || x.n_cols != y.n_cols)
    Rcpp::stop("%s observation dimensions disagree.",name);
}
void legacy_factor(const arma::mat& x){
  if (!x.is_finite() || x.n_cols == 0 || x.n_rows <= x.n_cols)
    Rcpp::stop("Fixed-rank geometry requires finite n-by-k factors with 1 <= k < n.");
  arma::vec sv=arma::svd(x);
  if (sv.min() <= 64.0*arma::datum::eps*sv.max())
    Rcpp::stop("Fixed-rank geometry requires full-column-rank factors.");
}
void legacy_horizontal(const arma::mat& x,const arma::mat& u,const char* name){
  legacy_pair(x,u,name);
  if (!u.is_finite()) Rcpp::stop("%s tangent is not finite.",name);
  arma::mat cross=x.t()*u;
  if (arma::norm(cross-cross.t(),"fro") >
      512.0*arma::datum::eps*arma::norm(x,"fro")*arma::norm(u,"fro"))
    Rcpp::stop("%s requires a horizontal tangent representation.",name);
}
arma::mat legacy_horizontal_projection(const arma::mat& x,const arma::mat& u){
  arma::mat cross=x.t()*u;
  arma::vec values; arma::mat vectors;
  if (!arma::eig_sym(values,vectors,x.t()*x) || values.min()<=0.0)
    Rcpp::stop("Horizontal projection requires a full-rank representative.");
  arma::mat skew=vectors.t()*(cross-cross.t())*vectors;
  for (arma::uword i=0;i<skew.n_rows;++i)
    for (arma::uword j=0;j<skew.n_cols;++j) skew(i,j)/=(values(i)+values(j));
  return u-x*vectors*skew*vectors.t();
}
void legacy_probability(const arma::mat& x){
  if (!x.is_finite() || x.n_cols != 1 || x.n_rows < 2 ||
      arma::any(arma::vectorise(x) <= 0.0) || std::abs(arma::accu(x)-1.0)>1e-8)
    Rcpp::stop("Multinomial geometry requires strictly positive probability vectors.");
}
void legacy_rotation(const arma::mat& x){
  legacy_orthonormal(x,"Rotation",true);
  if (arma::det(x)<=0.0) Rcpp::stop("Rotation geometry requires determinant +1.");
}
void legacy_skew(const arma::mat& x,const arma::mat& u){
  legacy_pair(x,u,"Rotation");
  if (!u.is_finite() || arma::norm(u+u.t(),"fro") >
      512.0*arma::datum::eps*arma::norm(u,"fro"))
    Rcpp::stop("Rotation tangents are represented by body-coordinate skew-symmetric matrices.");
}
void legacy_shape(const arma::mat& x){
  if (!x.is_finite() || x.n_cols<2 || x.n_rows<=x.n_cols ||
      std::abs(arma::norm(x,"fro")-1.0)>1e-8 ||
      arma::norm(arma::mean(x,0),2)>1e-8)
    Rcpp::stop("Landmark operations require centered, unit-norm regular preshapes.");
  arma::vec sv=arma::svd(x);
  if (sv.min()<=64.0*arma::datum::eps*sv.max())
    Rcpp::stop("Singular landmark configurations are outside the regular shape domain.");
}

// 03. CORRELATION =============================================================
// NOTE : find geodesic minimizing D for d(A, DBD) with some default parameters
arma::mat corr_airm_findD(arma::mat C1, arma::mat C2){
Rcpp::stop("Correlation quotient geometry is unavailable pending a compatible, accuracy-controlled implementation.");
}
arma::mat correlation_initialize(arma::field<arma::mat> data, arma::vec weight){
Rcpp::stop("Correlation quotient geometry is unavailable pending a compatible, accuracy-controlled implementation.");
}
double correlation_dist(arma::mat A, arma::mat B){
Rcpp::stop("Correlation quotient geometry is unavailable pending a compatible, accuracy-controlled implementation.");
}
arma::mat correlation_exp(arma::mat X, arma::mat eta, double t){
Rcpp::stop("Correlation quotient geometry is unavailable pending a compatible, accuracy-controlled implementation.");
}
arma::mat correlation_log(arma::mat X, arma::mat Y){
Rcpp::stop("Correlation quotient geometry is unavailable pending a compatible, accuracy-controlled implementation.");
}
double correlation_metric(arma::mat X, arma::mat eta1, arma::mat eta2){
Rcpp::stop("Correlation quotient geometry is unavailable pending a compatible, accuracy-controlled implementation.");
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
if (!x.is_finite() || x.n_cols==0 || x.n_rows<x.n_cols)
    Rcpp::stop("Stiefel projection requires a finite n-by-k matrix with k <= n.");
  arma::mat U,V; arma::vec sv;
  if (!arma::svd_econ(U,sv,V,x) || sv.min()<=64.0*arma::datum::eps*sv.max())
    Rcpp::stop("Stiefel inverse projection is nonunique for a rank-deficient embedding.");
  return U*V.t();
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
Rcpp::stop("Intrinsic Stiefel geometry is unavailable: the legacy logarithm and metric use incompatible conventions.");
}
double stiefel_dist(arma::mat x, arma::mat y){
Rcpp::stop("Intrinsic Stiefel geometry is unavailable: the legacy logarithm and metric use incompatible conventions.");
}
arma::vec stiefel_equiv(arma::mat x, int m, int n){
legacy_orthonormal(x,"Stiefel");
  return arma::vectorise(x);
}
arma::mat stiefel_invequiv(arma::vec x, int m, int n){
if (x.n_elem!=static_cast<arma::uword>(m)*static_cast<arma::uword>(n))
    Rcpp::stop("Stiefel embedding dimensions disagree.");
  return stiefel_nearest(arma::reshape(x,m,n));
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
legacy_orthonormal(X,"Grassmann"); legacy_orthonormal(Y,"Grassmann");
  legacy_pair(X,Y,"Grassmann");
  arma::mat U,V; arma::vec cosine;
  arma::svd(U,cosine,V,X.t()*Y);
  arma::mat normal=Y*V-X*U*arma::diagmat(cosine);
  arma::vec angle(cosine.n_elem);
  for (arma::uword i=0;i<cosine.n_elem;++i)
    angle(i)=std::atan2(arma::norm(normal.col(i),2),std::min(1.0,std::max(0.0,cosine(i))));
  return arma::norm(angle,2);
}
arma::mat grassmann_proj(arma::mat x, arma::mat u){
  return(u-x*((x.t()*u)));
}
arma::mat grassmann_nearest(arma::mat x){
  return(stiefel_nearest(x));
}
arma::mat grassmann_exp(arma::mat x, arma::mat d, double t){
  legacy_orthonormal(x,"Grassmann"); legacy_pair(x,d,"Grassmann");
  if (!d.is_finite() || !std::isfinite(t) ||
      arma::norm(x.t()*d,"fro")>512.0*arma::datum::eps*arma::norm(d,"fro"))
    Rcpp::stop("Grassmann exponential requires a finite horizontal tangent and step.");
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
legacy_orthonormal(x,"Grassmann"); legacy_orthonormal(y,"Grassmann");
  legacy_pair(x,y,"Grassmann");
  arma::mat U,V; arma::vec cosine;
  arma::svd(U,cosine,V,x.t()*y);
  if (cosine.min()<=128.0*arma::datum::eps)
    Rcpp::stop("Grassmann logarithm is nonunique or ill-conditioned at a principal angle pi/2.");
  arma::mat normal=y*V-x*U*arma::diagmat(cosine);
  for (arma::uword i=0;i<cosine.n_elem;++i) {
    const double sine=arma::norm(normal.col(i),2);
    if (sine>0.0) normal.col(i)*=std::atan2(sine,cosine(i))/sine;
  }
  arma::mat out=normal*U.t();
  return out-x*(x.t()*out);
}
arma::vec grassmann_equiv(arma::mat x, int n, int p){
legacy_orthonormal(x,"Grassmann");
  return arma::vectorise(x*x.t());
}
arma::mat grassmann_invequiv(arma::vec x, int n, int p){
if (n<1 || p<1 || p>n || x.n_elem!=static_cast<arma::uword>(n)*static_cast<arma::uword>(n) || !x.is_finite())
    Rcpp::stop("Grassmann embedding dimensions are invalid.");
  arma::mat raw=arma::reshape(x,n,n);
  arma::mat sym=raw/2.0+raw.t()/2.0;
  arma::vec values; arma::mat vectors;
  if (!arma::eig_sym(values,vectors,sym)) Rcpp::stop("Grassmann inverse projection failed.");
  if (p<n && values(n-p)-values(n-p-1)<=128.0*arma::datum::eps*arma::abs(values).max())
    Rcpp::stop("Grassmann inverse projection is nonunique at a tied cutoff eigenvalue.");
  return arma::fliplr(vectors.tail_cols(p));
}
double grassmann_distext(arma::mat x, arma::mat y){
  int m = x.n_rows;
  int n = x.n_cols;
  
  arma::vec xext = grassmann_equiv(x, m, n);
  arma::vec yext = grassmann_equiv(y, m, n);
  
  return(arma::norm(xext-yext,2));
}
arma::mat grassmann_initialize(arma::field<arma::mat> data, arma::vec weight){
if (data.n_elem==0) Rcpp::stop("Grassmann initialization requires observations.");
  const int n=data(0).n_rows; const int p=data(0).n_cols;
  arma::vec embedded(n*n,fill::zeros);
  for (arma::uword i=0;i<data.n_elem;++i)
    if (weight(i)>0.0) embedded+=(weight(i)/arma::accu(weight))*grassmann_equiv(data(i),n,p);
  return grassmann_invequiv(embedded,n,p);
}


arma::mat rotation_nearest(arma::mat x);

// 06. ROTATION ================================================================
// https://www.manopt.org/reference/manopt/manifolds/rotations/rotationsfactory.html
arma::mat rotation_initialize(arma::field<arma::mat> data, arma::vec weight){
arma::mat average(data(0).n_rows,data(0).n_cols,fill::zeros);
  for (arma::uword i=0;i<data.n_elem;++i)
    if (weight(i)>0.0) average+=(weight(i)/arma::accu(weight))*data(i);
  return rotation_nearest(average);
}
double rotation_metric(arma::mat X, arma::mat eta1, arma::mat eta2){
legacy_rotation(X); legacy_skew(X,eta1); legacy_skew(X,eta2);
  return arma::accu(eta1%eta2);
}
arma::mat rotation_log(arma::mat X, arma::mat Y){
legacy_rotation(X); legacy_rotation(Y); legacy_pair(X,Y,"Rotation");
  arma::mat relative=X.t()*Y;
  arma::vec gap=arma::svd(relative+arma::eye<arma::mat>(X.n_rows,X.n_rows));
  if (gap.min()<=128.0*arma::datum::eps)
    Rcpp::stop("Rotation logarithm is nonunique at an eigenvalue -1 (angle pi).");
  arma::cx_mat value=arma::logmat(relative);
  if (!value.is_finite() || arma::norm(arma::imag(value),"fro")>
      1e-10*std::max(1.0,arma::norm(value,"fro")))
    Rcpp::stop("Rotation principal logarithm did not produce a stable real tangent.");
  arma::mat real=arma::real(value);
  return real/2.0-real.t()/2.0;
}
arma::mat rotation_exp(arma::mat X, arma::mat eta, double t){
legacy_rotation(X); legacy_skew(X,eta);
  if (!std::isfinite(t)) Rcpp::stop("Rotation exponential requires a finite step.");
  arma::mat out=X*arma::expmat(t*eta);
  legacy_rotation(out);
  return out;
}
double rotation_dist(arma::mat X, arma::mat Y){
legacy_rotation(X); legacy_rotation(Y); legacy_pair(X,Y,"Rotation");
  arma::cx_vec values;
  if (!arma::eig_gen(values,X.t()*Y)) Rcpp::stop("Rotation distance eigendecomposition failed.");
  arma::vec angle(values.n_elem);
  for (arma::uword i=0;i<values.n_elem;++i)
    angle(i)=std::atan2(values(i).imag(),values(i).real());
  return arma::norm(angle,2);
}
arma::vec rotation_equiv(arma::mat x, int m, int n){
legacy_rotation(x);
  return arma::vectorise(x);
}
arma::mat rotation_nearest(arma::mat x){
if (!x.is_finite() || x.n_rows!=x.n_cols || x.n_rows==0)
    Rcpp::stop("Rotation inverse projection requires a finite square embedding.");
  if (x.n_rows==1) return arma::ones<arma::mat>(1,1);
  arma::mat U,V; arma::vec sv;
  if (!arma::svd(U,sv,V,x)) Rcpp::stop("Rotation inverse projection failed.");
  const double sign=arma::det(U*V.t())<0.0 ? -1.0 : 1.0;
  if (sv(sv.n_elem-2)+sign*sv(sv.n_elem-1)<=128.0*arma::datum::eps*sv.max())
    Rcpp::stop("Rotation inverse projection is nonunique at the smallest singular values.");
  arma::mat correction=arma::eye<arma::mat>(x.n_rows,x.n_rows);
  correction(x.n_rows-1,x.n_rows-1)=sign;
  return U*correction*V.t();
}
arma::mat rotation_invequiv(arma::vec x, int m, int n){
if (m!=n || x.n_elem!=static_cast<arma::uword>(m)*static_cast<arma::uword>(n))
    Rcpp::stop("Rotation embedding dimensions disagree.");
  return rotation_nearest(arma::reshape(x,m,n));
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
legacy_probability(x); legacy_pair(x,d1,"Multinomial"); legacy_pair(x,d2,"Multinomial");
  if (!d1.is_finite() || !d2.is_finite() ||
      std::abs(arma::accu(d1))>512.0*arma::datum::eps*arma::norm(d1,"fro") ||
      std::abs(arma::accu(d2))>512.0*arma::datum::eps*arma::norm(d2,"fro"))
    Rcpp::stop("Multinomial tangents must be finite and sum to zero.");
  return arma::accu((d1%d2)/x);
}
arma::mat multinomial_log(arma::mat x, arma::mat y){
legacy_probability(x); legacy_probability(y); legacy_pair(x,y,"Multinomial");
  arma::mat rootx=arma::sqrt(x);
  arma::mat out=2.0*rootx%sphere_log(rootx,arma::sqrt(y));
  // Remove accumulated roundoff in the probability tangent constraint.
  return out-arma::accu(out)*x;
}
arma::mat multinomial_exp(arma::mat x, arma::mat u, double t){
legacy_probability(x); legacy_pair(x,u,"Multinomial");
  if (!u.is_finite() || !std::isfinite(t) ||
      std::abs(arma::accu(u))>512.0*arma::datum::eps*arma::norm(u,"fro"))
    Rcpp::stop("Multinomial exponential requires a finite zero-sum tangent and step.");
  arma::mat rootx=arma::sqrt(x);
  arma::mat velocity=(t*u)/(2.0*rootx);
  const double angle=arma::norm(velocity,"fro");
  if (angle==0.0) return x;
  if (!std::isfinite(angle)) Rcpp::stop("Multinomial step exceeds numerical range.");
  arma::mat direction=velocity/angle;
  for (arma::uword i=0;i<x.n_elem;++i) {
    const double boundary=std::atan2(rootx(i),-direction(i));
    if (angle>=boundary)
      Rcpp::stop("Multinomial exponential reaches the boundary of the positive simplex.");
  }
  arma::mat outroot=std::cos(angle)*rootx+std::sin(angle)*direction;
  arma::mat out=outroot%outroot;
  out/=arma::accu(out);
  legacy_probability(out);
  return out;
}
double multinomial_dist(arma::mat x, arma::mat y){
legacy_probability(x); legacy_probability(y); legacy_pair(x,y,"Multinomial");
  return 2.0*sphere_dist(arma::sqrt(x),arma::sqrt(y));
}
arma::vec multinomial_equiv(arma::mat x, int m, int n){
legacy_probability(x);
  return 2.0*arma::sqrt(arma::vectorise(x));
}

arma::mat multinomial_invequiv(arma::vec x, int m, int n){
if (n!=1 || m<2 || x.n_elem!=static_cast<arma::uword>(m) || !x.is_finite() || arma::any(x<=0.0))
    Rcpp::stop("Multinomial inverse embedding requires strictly positive compatible coordinates.");
  arma::vec scaled=x/x.max();
  arma::mat out=arma::reshape(scaled%scaled,m,n);
  out/=arma::accu(out);
  legacy_probability(out);
  return out;
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
const arma::mat reference=data(0);
  legacy_factor(reference);
  arma::mat out(reference.n_rows,reference.n_cols,fill::zeros);
  for (arma::uword i=0;i<data.n_elem;++i) {
    if (weight(i)==0.0) continue;
    legacy_factor(data(i)); legacy_pair(reference,data(i),"Fixed-rank");
    arma::mat U,V; arma::vec sv;
    arma::svd(U,sv,V,reference.t()*data(i));
    out+=(weight(i)/arma::accu(weight))*data(i)*V*U.t();
  }
  legacy_factor(out);
  return out;
}
double spdk_metric(arma::mat X, arma::mat eta, arma::mat zeta){
legacy_factor(X); legacy_horizontal(X,eta,"Fixed-rank"); legacy_horizontal(X,zeta,"Fixed-rank");
  return arma::accu(eta%zeta);
}
arma::mat spdk_log(arma::mat Y, arma::mat Z){
legacy_factor(Y); legacy_factor(Z); legacy_pair(Y,Z,"Fixed-rank");
  arma::mat U,V; arma::vec sv;
  arma::svd(U,sv,V,Y.t()*Z);
  if (sv.min()<=128.0*arma::datum::eps*arma::norm(Y,2)*arma::norm(Z,2))
    Rcpp::stop("Fixed-rank logarithm is nonunique at a singular cross-Gram matrix.");
  return legacy_horizontal_projection(Y,Z*V*U.t()-Y);
}
arma::mat spdk_exp(arma::mat Y, arma::mat eta, double t){
legacy_factor(Y); legacy_horizontal(Y,eta,"Fixed-rank");
  if (!std::isfinite(t)) Rcpp::stop("Fixed-rank exponential requires a finite step.");
  arma::mat out=Y+t*eta;
  legacy_factor(out);
  arma::mat cross=Y.t()*out;
  arma::mat factor;
  if (!arma::chol(factor,cross/2.0+cross.t()/2.0))
    Rcpp::stop("Fixed-rank exponential left the supported local horizontal domain.");
  return out;
}
double spdk_dist(arma::mat X, arma::mat Y){
legacy_factor(X); legacy_factor(Y); legacy_pair(X,Y,"Fixed-rank");
  arma::mat U,V; arma::vec sv;
  arma::svd(U,sv,V,X.t()*Y);
  return arma::norm(Y*V*U.t()-X,"fro");
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
arma::mat landmark_aux_nearest(arma::mat x){
if (!x.is_finite() || x.n_cols<2 || x.n_rows<=x.n_cols)
    Rcpp::stop("Landmark configuration dimensions are invalid.");
  const double scale=arma::abs(x).max();
  if (scale==0.0) Rcpp::stop("Constant landmark configurations have no shape.");
  arma::mat centered=x/scale;
  centered.each_row()-=arma::mean(centered,0);
  const double size=arma::norm(centered,"fro");
  if (size==0.0) Rcpp::stop("Constant landmark configurations have no shape.");
  arma::mat out=centered/size;
  legacy_shape(out);
  return out;
}
arma::mat landmark_aux_matching(arma::mat x, arma::mat y){
legacy_shape(x); legacy_shape(y); legacy_pair(x,y,"Landmark");
  arma::mat U,V; arma::vec sv;
  arma::svd(U,sv,V,x.t()*y);
  return y*V*U.t();
}
arma::mat landmark_initialize(arma::field<arma::mat> data, arma::vec weight){
const arma::mat reference=data(0);
  arma::mat out(reference.n_rows,reference.n_cols,fill::zeros);
  for (arma::uword i=0;i<data.n_elem;++i)
    if (weight(i)>0.0) out+=(weight(i)/arma::accu(weight))*landmark_aux_matching(reference,data(i));
  return landmark_aux_nearest(out);
}
arma::mat landmark_exp(arma::mat x, arma::mat d, double t){
legacy_shape(x); legacy_horizontal(x,d,"Landmark");
  const double norm=arma::norm(d,"fro");
  if (!std::isfinite(t) || arma::norm(arma::mean(d,0),2)>512.0*arma::datum::eps*norm ||
      std::abs(arma::accu(x%d))>512.0*arma::datum::eps*norm)
    Rcpp::stop("Landmark exponential requires a centered horizontal sphere tangent.");
  if (std::abs(t)*norm>=arma::datum::pi/2.0)
    Rcpp::stop("Landmark exponential is restricted to the regular local domain of radius pi/2.");
  arma::mat out=sphere_exp(x,d,t);
  legacy_shape(out);
  arma::mat cross=x.t()*out;
  arma::mat factor;
  if (!arma::chol(factor,cross/2.0+cross.t()/2.0))
    Rcpp::stop("Landmark exponential left the supported regular horizontal domain.");
  return out;
}
arma::mat landmark_log(arma::mat X, arma::mat Y){
legacy_shape(X); legacy_shape(Y); legacy_pair(X,Y,"Landmark");
  arma::mat U,V; arma::vec sv;
  arma::svd(U,sv,V,X.t()*Y);
  if (sv.min()<=128.0*arma::datum::eps*arma::norm(X,2)*arma::norm(Y,2))
    Rcpp::stop("Landmark logarithm is nonunique at a singular cross-Gram matrix.");
  arma::mat aligned=Y*V*U.t();
  arma::mat out=legacy_horizontal_projection(X,sphere_log(X,aligned));
  out.each_row()-=arma::mean(out,0);
  return out-arma::accu(X%out)*X;
}
double landmark_metric(arma::mat x, arma::mat d1, arma::mat d2){
legacy_shape(x); legacy_horizontal(x,d1,"Landmark"); legacy_horizontal(x,d2,"Landmark");
  return arma::accu(d1%d2);
}
double landmark_dist(arma::mat x, arma::mat y){
arma::mat aligned=landmark_aux_matching(x,y);
  return sphere_dist(x,aligned);
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
