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
namespace {
// Both the loss and its derivatives use the same normalized distances. A
// projected collision is excluded while optimizing because the loss is not
// differentiable there, even when its limiting value is finite.
double sammon_loss(const arma::mat& distances, const arma::mat& embedding,
                   double total, arma::mat& embedded_distances) {
  double loss = 0.0;
  const arma::uword n = distances.n_rows;
  for (arma::uword i = 0; i + 1 < n; ++i) {
    for (arma::uword j = i + 1; j < n; ++j) {
      const double d = arma::norm(embedding.row(i) - embedding.row(j), 2);
      if (!std::isfinite(d) || d <= 0.0) return arma::datum::inf;
      embedded_distances(i, j) = embedded_distances(j, i) = d;
      const double residual = distances(i, j) - d;
      loss += residual * (residual / distances(i, j));
    }
  }
  return loss / total;
}
}

// [[Rcpp::export]]
Rcpp::List visualize_sammon(std::string mfd, std::string geo, Rcpp::List& data,
                           int ndim, int maxiter, double abstol) {
  const int N = data.size();
  if (N < 2 || ndim < 1 || ndim >= N || maxiter < 1 ||
      !std::isfinite(abstol) || abstol < 0.0) {
    Rcpp::stop("Sammon mapping requires at least two observations, 1 <= ndim < N, maxiter >= 1, and finite eps >= 0.");
  }
  const arma::mat first = Rcpp::as<arma::mat>(data[0]);
  arma::cube mydata(first.n_rows, first.n_cols, N, arma::fill::zeros);
  for (int n = 0; n < N; ++n) mydata.slice(n) = Rcpp::as<arma::mat>(data[n]);

  arma::mat distances(N, N, arma::fill::zeros);
  for (int i = 0; i < N - 1; ++i) {
    for (int j = i + 1; j < N; ++j) {
      const double d = geo == "intrinsic" ?
        riem_dist(mfd, mydata.slice(i), mydata.slice(j)) :
        riem_distext(mfd, mydata.slice(i), mydata.slice(j));
      if (!std::isfinite(d) || d <= 0.0) {
        Rcpp::stop("Sammon mapping requires strictly positive finite distances between distinct observations; remove coincident points.");
      }
      distances(i, j) = distances(j, i) = d;
    }
  }
  const double distance_scale = distances.max();
  distances /= distance_scale;
  double total = 0.0;
  for (int i = 0; i < N - 1; ++i) {
    for (int j = i + 1; j < N; ++j) {
      if (distances(i, j) <= 0.0) {
        Rcpp::stop("The relative distance scale is outside the positive numerical range.");
      }
      total += distances(i, j);
    }
  }

  // Classical scaling on normalized distances. Negative eigenvalues in a
  // non-Euclidean dissimilarity matrix are not real coordinate directions.
  const arma::mat centering = arma::eye<arma::mat>(N, N) -
    arma::ones<arma::mat>(N, N) / static_cast<double>(N);
  arma::mat gram = -0.5 * centering * arma::square(distances) * centering;
  gram = 0.5 * (gram + gram.t());
  arma::vec values;
  arma::mat vectors;
  if (!arma::eig_sym(values, vectors, gram)) {
    Rcpp::stop("Sammon initialization failed to diagonalize the distance matrix.");
  }
  arma::mat embedding = vectors.tail_cols(ndim) *
    arma::diagmat(arma::sqrt(arma::clamp(values.tail(ndim), 0.0, arma::datum::inf)));
  arma::mat embedded_distances(N, N, arma::fill::zeros);
  double loss = sammon_loss(distances, embedding, total, embedded_distances);
  if (!std::isfinite(loss)) {
    // A lower-dimensional projection can coincide even for distinct inputs.
    // Separate those coordinates deterministically, without consuming R's RNG.
    const arma::mat initial = embedding;
    double perturbation = std::sqrt(arma::datum::eps);
    for (int trial = 0; trial < 10 && !std::isfinite(loss); ++trial) {
      embedding = initial;
      for (int i = 0; i < N; ++i) {
        embedding(i, 0) += perturbation * static_cast<double>(i) / N;
      }
      embedding.each_row() -= arma::mean(embedding, 0);
      loss = sammon_loss(distances, embedding, total, embedded_distances);
      perturbation *= 10.0;
    }
  }
  if (!std::isfinite(loss)) {
    Rcpp::stop("Sammon initialization is outside the finite numerical range.");
  }

  for (int iteration = 0; iteration < maxiter; ++iteration) {
    Rcpp::checkUserInterrupt();
    arma::mat gradient(N, ndim, arma::fill::zeros);
    arma::mat hessian(N, ndim, arma::fill::zeros);
    for (int i = 0; i < N - 1; ++i) {
      for (int j = i + 1; j < N; ++j) {
        const double d = embedded_distances(i, j);
        const double relative_error = d / distances(i, j) - 1.0;
        for (int q = 0; q < ndim; ++q) {
          const double direction = (embedding(i, q) - embedding(j, q)) / d;
          const double g = (2.0 / total) * relative_error * direction;
          // Algebraically 2/c * [1/D - 1/d + delta^2/d^3]. The legacy
          // expression omitted the pair-specific denominator D*d entirely.
          const double h = (2.0 / total) *
            (relative_error + direction * direction) / d;
          gradient(i, q) += g;
          gradient(j, q) -= g;
          hessian(i, q) += h;
          hessian(j, q) += h;
        }
      }
    }
    if (!gradient.is_finite() || !hessian.is_finite()) {
      Rcpp::stop("Sammon derivatives are outside the finite numerical range.");
    }
    if (arma::abs(gradient).max() == 0.0) break;
    hessian = arma::abs(hessian);
    const double hessian_scale = hessian.max();
    const double floor = hessian_scale > 0.0 ?
      std::max(1e-12 * hessian_scale, arma::datum::eps) : 1.0;
    arma::mat direction = -gradient / arma::clamp(hessian, floor, arma::datum::inf);
    const double direction_scale = arma::abs(direction).max();
    if (!direction.is_finite()) {
      Rcpp::stop("The Sammon update is outside the finite numerical range.");
    }
    if (direction_scale > 1.0) direction /= direction_scale;
    const double slope = arma::accu(gradient % direction);
    if (!std::isfinite(slope) || slope >= 0.0) break;

    double step = 0.3;
    bool accepted = false;
    arma::mat candidate, candidate_distances(N, N, arma::fill::zeros);
    double candidate_loss = loss;
    for (int backtrack = 0; backtrack < 60; ++backtrack) {
      candidate = embedding + step * direction;
      candidate.each_row() -= arma::mean(candidate, 0);
      candidate_loss = sammon_loss(distances, candidate, total, candidate_distances);
      if (std::isfinite(candidate_loss) && candidate_loss < loss &&
          candidate_loss <= loss + 1e-4 * step * slope) {
        accepted = true;
        break;
      }
      step *= 0.5;
    }
    if (!accepted) break;
    const double increment = arma::norm(candidate - embedding, "fro") /
      std::sqrt(static_cast<double>(N) * ndim);
    embedding = candidate;
    embedded_distances = candidate_distances;
    loss = candidate_loss;
    if (increment <= abstol) break;
  }

  const double stress = engine_stress(distances, embedded_distances);
  embedding *= distance_scale;
  if (!embedding.is_finite() || !std::isfinite(stress)) {
    Rcpp::stop("The Sammon result is outside the finite numerical range.");
  }
  return Rcpp::List::create(Rcpp::Named("embed") = embedding,
                            Rcpp::Named("stress") = stress);
}
