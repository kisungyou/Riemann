// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// basic_pdist
arma::mat basic_pdist(std::string mfdname, Rcpp::List& data, std::string dtype);
RcppExport SEXP _Riemann_basic_pdist(SEXP mfdnameSEXP, SEXP dataSEXP, SEXP dtypeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type mfdname(mfdnameSEXP);
    Rcpp::traits::input_parameter< Rcpp::List& >::type data(dataSEXP);
    Rcpp::traits::input_parameter< std::string >::type dtype(dtypeSEXP);
    rcpp_result_gen = Rcpp::wrap(basic_pdist(mfdname, data, dtype));
    return rcpp_result_gen;
END_RCPP
}
// basic_pdist2
arma::mat basic_pdist2(std::string mfdname, Rcpp::List& data1, Rcpp::List& data2, std::string dtype);
RcppExport SEXP _Riemann_basic_pdist2(SEXP mfdnameSEXP, SEXP data1SEXP, SEXP data2SEXP, SEXP dtypeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type mfdname(mfdnameSEXP);
    Rcpp::traits::input_parameter< Rcpp::List& >::type data1(data1SEXP);
    Rcpp::traits::input_parameter< Rcpp::List& >::type data2(data2SEXP);
    Rcpp::traits::input_parameter< std::string >::type dtype(dtypeSEXP);
    rcpp_result_gen = Rcpp::wrap(basic_pdist2(mfdname, data1, data2, dtype));
    return rcpp_result_gen;
END_RCPP
}
// basic_interpolate
arma::cube basic_interpolate(std::string mfdname, std::string dtype, arma::mat mat1, arma::mat mat2, arma::vec vect);
RcppExport SEXP _Riemann_basic_interpolate(SEXP mfdnameSEXP, SEXP dtypeSEXP, SEXP mat1SEXP, SEXP mat2SEXP, SEXP vectSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type mfdname(mfdnameSEXP);
    Rcpp::traits::input_parameter< std::string >::type dtype(dtypeSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type mat1(mat1SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type mat2(mat2SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type vect(vectSEXP);
    rcpp_result_gen = Rcpp::wrap(basic_interpolate(mfdname, dtype, mat1, mat2, vect));
    return rcpp_result_gen;
END_RCPP
}
// inference_mean_intrinsic
Rcpp::List inference_mean_intrinsic(std::string mfdname, Rcpp::List& data, arma::vec myweight, int myiter, double myeps);
RcppExport SEXP _Riemann_inference_mean_intrinsic(SEXP mfdnameSEXP, SEXP dataSEXP, SEXP myweightSEXP, SEXP myiterSEXP, SEXP myepsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type mfdname(mfdnameSEXP);
    Rcpp::traits::input_parameter< Rcpp::List& >::type data(dataSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type myweight(myweightSEXP);
    Rcpp::traits::input_parameter< int >::type myiter(myiterSEXP);
    Rcpp::traits::input_parameter< double >::type myeps(myepsSEXP);
    rcpp_result_gen = Rcpp::wrap(inference_mean_intrinsic(mfdname, data, myweight, myiter, myeps));
    return rcpp_result_gen;
END_RCPP
}
// inference_mean_extrinsic
Rcpp::List inference_mean_extrinsic(std::string mfdname, Rcpp::List& data, arma::vec myweight, int myiter, double myeps);
RcppExport SEXP _Riemann_inference_mean_extrinsic(SEXP mfdnameSEXP, SEXP dataSEXP, SEXP myweightSEXP, SEXP myiterSEXP, SEXP myepsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type mfdname(mfdnameSEXP);
    Rcpp::traits::input_parameter< Rcpp::List& >::type data(dataSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type myweight(myweightSEXP);
    Rcpp::traits::input_parameter< int >::type myiter(myiterSEXP);
    Rcpp::traits::input_parameter< double >::type myeps(myepsSEXP);
    rcpp_result_gen = Rcpp::wrap(inference_mean_extrinsic(mfdname, data, myweight, myiter, myeps));
    return rcpp_result_gen;
END_RCPP
}
// inference_median_intrinsic
Rcpp::List inference_median_intrinsic(std::string mfdname, Rcpp::List& data, arma::vec myweight, int myiter, double myeps);
RcppExport SEXP _Riemann_inference_median_intrinsic(SEXP mfdnameSEXP, SEXP dataSEXP, SEXP myweightSEXP, SEXP myiterSEXP, SEXP myepsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type mfdname(mfdnameSEXP);
    Rcpp::traits::input_parameter< Rcpp::List& >::type data(dataSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type myweight(myweightSEXP);
    Rcpp::traits::input_parameter< int >::type myiter(myiterSEXP);
    Rcpp::traits::input_parameter< double >::type myeps(myepsSEXP);
    rcpp_result_gen = Rcpp::wrap(inference_median_intrinsic(mfdname, data, myweight, myiter, myeps));
    return rcpp_result_gen;
END_RCPP
}
// inference_median_extrinsic
Rcpp::List inference_median_extrinsic(std::string mfdname, Rcpp::List& data, arma::vec myweight, int myiter, double myeps);
RcppExport SEXP _Riemann_inference_median_extrinsic(SEXP mfdnameSEXP, SEXP dataSEXP, SEXP myweightSEXP, SEXP myiterSEXP, SEXP myepsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type mfdname(mfdnameSEXP);
    Rcpp::traits::input_parameter< Rcpp::List& >::type data(dataSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type myweight(myweightSEXP);
    Rcpp::traits::input_parameter< int >::type myiter(myiterSEXP);
    Rcpp::traits::input_parameter< double >::type myeps(myepsSEXP);
    rcpp_result_gen = Rcpp::wrap(inference_median_extrinsic(mfdname, data, myweight, myiter, myeps));
    return rcpp_result_gen;
END_RCPP
}
// clustering_nmshift
Rcpp::List clustering_nmshift(std::string mfdname, Rcpp::List& data, double h, int iter, double eps);
RcppExport SEXP _Riemann_clustering_nmshift(SEXP mfdnameSEXP, SEXP dataSEXP, SEXP hSEXP, SEXP iterSEXP, SEXP epsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type mfdname(mfdnameSEXP);
    Rcpp::traits::input_parameter< Rcpp::List& >::type data(dataSEXP);
    Rcpp::traits::input_parameter< double >::type h(hSEXP);
    Rcpp::traits::input_parameter< int >::type iter(iterSEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    rcpp_result_gen = Rcpp::wrap(clustering_nmshift(mfdname, data, h, iter, eps));
    return rcpp_result_gen;
END_RCPP
}
// clustering_kmeans_lloyd
Rcpp::List clustering_kmeans_lloyd(std::string mfdname, std::string geotype, Rcpp::List& data, int iter, double eps, arma::uvec initlabel);
RcppExport SEXP _Riemann_clustering_kmeans_lloyd(SEXP mfdnameSEXP, SEXP geotypeSEXP, SEXP dataSEXP, SEXP iterSEXP, SEXP epsSEXP, SEXP initlabelSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type mfdname(mfdnameSEXP);
    Rcpp::traits::input_parameter< std::string >::type geotype(geotypeSEXP);
    Rcpp::traits::input_parameter< Rcpp::List& >::type data(dataSEXP);
    Rcpp::traits::input_parameter< int >::type iter(iterSEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type initlabel(initlabelSEXP);
    rcpp_result_gen = Rcpp::wrap(clustering_kmeans_lloyd(mfdname, geotype, data, iter, eps, initlabel));
    return rcpp_result_gen;
END_RCPP
}
// clustering_kmeans_macqueen
Rcpp::List clustering_kmeans_macqueen(std::string mfdname, std::string geotype, Rcpp::List& data, int iter, double eps, arma::uvec initlabel);
RcppExport SEXP _Riemann_clustering_kmeans_macqueen(SEXP mfdnameSEXP, SEXP geotypeSEXP, SEXP dataSEXP, SEXP iterSEXP, SEXP epsSEXP, SEXP initlabelSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type mfdname(mfdnameSEXP);
    Rcpp::traits::input_parameter< std::string >::type geotype(geotypeSEXP);
    Rcpp::traits::input_parameter< Rcpp::List& >::type data(dataSEXP);
    Rcpp::traits::input_parameter< int >::type iter(iterSEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type initlabel(initlabelSEXP);
    rcpp_result_gen = Rcpp::wrap(clustering_kmeans_macqueen(mfdname, geotype, data, iter, eps, initlabel));
    return rcpp_result_gen;
END_RCPP
}
// clustering_clrq
Rcpp::List clustering_clrq(std::string mfdname, Rcpp::List& data, arma::uvec init_label, double par_a, double par_b);
RcppExport SEXP _Riemann_clustering_clrq(SEXP mfdnameSEXP, SEXP dataSEXP, SEXP init_labelSEXP, SEXP par_aSEXP, SEXP par_bSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type mfdname(mfdnameSEXP);
    Rcpp::traits::input_parameter< Rcpp::List& >::type data(dataSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type init_label(init_labelSEXP);
    Rcpp::traits::input_parameter< double >::type par_a(par_aSEXP);
    Rcpp::traits::input_parameter< double >::type par_b(par_bSEXP);
    rcpp_result_gen = Rcpp::wrap(clustering_clrq(mfdname, data, init_label, par_a, par_b));
    return rcpp_result_gen;
END_RCPP
}
// clustering_sup_intrinsic
Rcpp::List clustering_sup_intrinsic(std::string mfdname, Rcpp::List& data, arma::vec weight, double multiplier, int maxiter, double eps);
RcppExport SEXP _Riemann_clustering_sup_intrinsic(SEXP mfdnameSEXP, SEXP dataSEXP, SEXP weightSEXP, SEXP multiplierSEXP, SEXP maxiterSEXP, SEXP epsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type mfdname(mfdnameSEXP);
    Rcpp::traits::input_parameter< Rcpp::List& >::type data(dataSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type weight(weightSEXP);
    Rcpp::traits::input_parameter< double >::type multiplier(multiplierSEXP);
    Rcpp::traits::input_parameter< int >::type maxiter(maxiterSEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    rcpp_result_gen = Rcpp::wrap(clustering_sup_intrinsic(mfdname, data, weight, multiplier, maxiter, eps));
    return rcpp_result_gen;
END_RCPP
}
// clustering_kmeans18B
Rcpp::List clustering_kmeans18B(std::string mfdname, std::string geotype, Rcpp::List& data, int K, int M, int maxiter);
RcppExport SEXP _Riemann_clustering_kmeans18B(SEXP mfdnameSEXP, SEXP geotypeSEXP, SEXP dataSEXP, SEXP KSEXP, SEXP MSEXP, SEXP maxiterSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type mfdname(mfdnameSEXP);
    Rcpp::traits::input_parameter< std::string >::type geotype(geotypeSEXP);
    Rcpp::traits::input_parameter< Rcpp::List& >::type data(dataSEXP);
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    Rcpp::traits::input_parameter< int >::type M(MSEXP);
    Rcpp::traits::input_parameter< int >::type maxiter(maxiterSEXP);
    rcpp_result_gen = Rcpp::wrap(clustering_kmeans18B(mfdname, geotype, data, K, M, maxiter));
    return rcpp_result_gen;
END_RCPP
}
// cpp_scNJW
Rcpp::List cpp_scNJW(arma::mat& D, int K, double sigma, bool usekmeans, int maxiter);
RcppExport SEXP _Riemann_cpp_scNJW(SEXP DSEXP, SEXP KSEXP, SEXP sigmaSEXP, SEXP usekmeansSEXP, SEXP maxiterSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type D(DSEXP);
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    Rcpp::traits::input_parameter< double >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< bool >::type usekmeans(usekmeansSEXP);
    Rcpp::traits::input_parameter< int >::type maxiter(maxiterSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_scNJW(D, K, sigma, usekmeans, maxiter));
    return rcpp_result_gen;
END_RCPP
}
// cpp_scUL
Rcpp::List cpp_scUL(arma::mat& D, int K, double sigma, bool usekmeans, int maxiter);
RcppExport SEXP _Riemann_cpp_scUL(SEXP DSEXP, SEXP KSEXP, SEXP sigmaSEXP, SEXP usekmeansSEXP, SEXP maxiterSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type D(DSEXP);
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    Rcpp::traits::input_parameter< double >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< bool >::type usekmeans(usekmeansSEXP);
    Rcpp::traits::input_parameter< int >::type maxiter(maxiterSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_scUL(D, K, sigma, usekmeans, maxiter));
    return rcpp_result_gen;
END_RCPP
}
// cpp_scSM
Rcpp::List cpp_scSM(arma::mat& D, int K, double sigma, bool usekmeans, int maxiter);
RcppExport SEXP _Riemann_cpp_scSM(SEXP DSEXP, SEXP KSEXP, SEXP sigmaSEXP, SEXP usekmeansSEXP, SEXP maxiterSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type D(DSEXP);
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    Rcpp::traits::input_parameter< double >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< bool >::type usekmeans(usekmeansSEXP);
    Rcpp::traits::input_parameter< int >::type maxiter(maxiterSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_scSM(D, K, sigma, usekmeans, maxiter));
    return rcpp_result_gen;
END_RCPP
}
// cpp_sc05Z
Rcpp::List cpp_sc05Z(arma::mat& D, int K, int nnbd, bool usekmeans, int maxiter);
RcppExport SEXP _Riemann_cpp_sc05Z(SEXP DSEXP, SEXP KSEXP, SEXP nnbdSEXP, SEXP usekmeansSEXP, SEXP maxiterSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type D(DSEXP);
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    Rcpp::traits::input_parameter< int >::type nnbd(nnbdSEXP);
    Rcpp::traits::input_parameter< bool >::type usekmeans(usekmeansSEXP);
    Rcpp::traits::input_parameter< int >::type maxiter(maxiterSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_sc05Z(D, K, nnbd, usekmeans, maxiter));
    return rcpp_result_gen;
END_RCPP
}
// visualize_pga
Rcpp::List visualize_pga(std::string mfdname, Rcpp::List& data);
RcppExport SEXP _Riemann_visualize_pga(SEXP mfdnameSEXP, SEXP dataSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type mfdname(mfdnameSEXP);
    Rcpp::traits::input_parameter< Rcpp::List& >::type data(dataSEXP);
    rcpp_result_gen = Rcpp::wrap(visualize_pga(mfdname, data));
    return rcpp_result_gen;
END_RCPP
}
// visualize_kpca
Rcpp::List visualize_kpca(std::string mfdname, Rcpp::List& data, double sigma, int ndim);
RcppExport SEXP _Riemann_visualize_kpca(SEXP mfdnameSEXP, SEXP dataSEXP, SEXP sigmaSEXP, SEXP ndimSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type mfdname(mfdnameSEXP);
    Rcpp::traits::input_parameter< Rcpp::List& >::type data(dataSEXP);
    Rcpp::traits::input_parameter< double >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< int >::type ndim(ndimSEXP);
    rcpp_result_gen = Rcpp::wrap(visualize_kpca(mfdname, data, sigma, ndim));
    return rcpp_result_gen;
END_RCPP
}
// visualize_isomap
arma::mat visualize_isomap(std::string mfdname, Rcpp::List& data, std::string geometry, int nnbd);
RcppExport SEXP _Riemann_visualize_isomap(SEXP mfdnameSEXP, SEXP dataSEXP, SEXP geometrySEXP, SEXP nnbdSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type mfdname(mfdnameSEXP);
    Rcpp::traits::input_parameter< Rcpp::List& >::type data(dataSEXP);
    Rcpp::traits::input_parameter< std::string >::type geometry(geometrySEXP);
    Rcpp::traits::input_parameter< int >::type nnbd(nnbdSEXP);
    rcpp_result_gen = Rcpp::wrap(visualize_isomap(mfdname, data, geometry, nnbd));
    return rcpp_result_gen;
END_RCPP
}
// learning_seb
Rcpp::List learning_seb(std::string mfdname, Rcpp::List& data, int myiter, double myeps, std::string method);
RcppExport SEXP _Riemann_learning_seb(SEXP mfdnameSEXP, SEXP dataSEXP, SEXP myiterSEXP, SEXP myepsSEXP, SEXP methodSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type mfdname(mfdnameSEXP);
    Rcpp::traits::input_parameter< Rcpp::List& >::type data(dataSEXP);
    Rcpp::traits::input_parameter< int >::type myiter(myiterSEXP);
    Rcpp::traits::input_parameter< double >::type myeps(myepsSEXP);
    Rcpp::traits::input_parameter< std::string >::type method(methodSEXP);
    rcpp_result_gen = Rcpp::wrap(learning_seb(mfdname, data, myiter, myeps, method));
    return rcpp_result_gen;
END_RCPP
}
// learning_rmml
arma::mat learning_rmml(std::string mfdname, Rcpp::List& data, double lambda, arma::uvec label);
RcppExport SEXP _Riemann_learning_rmml(SEXP mfdnameSEXP, SEXP dataSEXP, SEXP lambdaSEXP, SEXP labelSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type mfdname(mfdnameSEXP);
    Rcpp::traits::input_parameter< Rcpp::List& >::type data(dataSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type label(labelSEXP);
    rcpp_result_gen = Rcpp::wrap(learning_rmml(mfdname, data, lambda, label));
    return rcpp_result_gen;
END_RCPP
}
// learning_coreset18B
Rcpp::List learning_coreset18B(std::string mfdname, std::string geoname, Rcpp::List& data, int M, int myiter, double myeps);
RcppExport SEXP _Riemann_learning_coreset18B(SEXP mfdnameSEXP, SEXP geonameSEXP, SEXP dataSEXP, SEXP MSEXP, SEXP myiterSEXP, SEXP myepsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type mfdname(mfdnameSEXP);
    Rcpp::traits::input_parameter< std::string >::type geoname(geonameSEXP);
    Rcpp::traits::input_parameter< Rcpp::List& >::type data(dataSEXP);
    Rcpp::traits::input_parameter< int >::type M(MSEXP);
    Rcpp::traits::input_parameter< int >::type myiter(myiterSEXP);
    Rcpp::traits::input_parameter< double >::type myeps(myepsSEXP);
    rcpp_result_gen = Rcpp::wrap(learning_coreset18B(mfdname, geoname, data, M, myiter, myeps));
    return rcpp_result_gen;
END_RCPP
}
// runif_sphere
arma::mat runif_sphere(int n, int p);
RcppExport SEXP _Riemann_runif_sphere(SEXP nSEXP, SEXP pSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    rcpp_result_gen = Rcpp::wrap(runif_sphere(n, p));
    return rcpp_result_gen;
END_RCPP
}
// runif_stiefel
arma::cube runif_stiefel(int p, int k, int N);
RcppExport SEXP _Riemann_runif_stiefel(SEXP pSEXP, SEXP kSEXP, SEXP NSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    rcpp_result_gen = Rcpp::wrap(runif_stiefel(p, k, N));
    return rcpp_result_gen;
END_RCPP
}
// mat_rank
arma::uword mat_rank(arma::mat A);
RcppExport SEXP _Riemann_mat_rank(SEXP ASEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type A(ASEXP);
    rcpp_result_gen = Rcpp::wrap(mat_rank(A));
    return rcpp_result_gen;
END_RCPP
}
// mat_symm
arma::mat mat_symm(arma::mat A, bool diag);
RcppExport SEXP _Riemann_mat_symm(SEXP ASEXP, SEXP diagSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type A(ASEXP);
    Rcpp::traits::input_parameter< bool >::type diag(diagSEXP);
    rcpp_result_gen = Rcpp::wrap(mat_symm(A, diag));
    return rcpp_result_gen;
END_RCPP
}
// mat_diaghalf
arma::mat mat_diaghalf(arma::mat D);
RcppExport SEXP _Riemann_mat_diaghalf(SEXP DSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type D(DSEXP);
    rcpp_result_gen = Rcpp::wrap(mat_diaghalf(D));
    return rcpp_result_gen;
END_RCPP
}
// mat_diaginvhalf
arma::mat mat_diaginvhalf(arma::mat D);
RcppExport SEXP _Riemann_mat_diaginvhalf(SEXP DSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type D(DSEXP);
    rcpp_result_gen = Rcpp::wrap(mat_diaginvhalf(D));
    return rcpp_result_gen;
END_RCPP
}
// mat_cov2cor
arma::mat mat_cov2cor(arma::mat A);
RcppExport SEXP _Riemann_mat_cov2cor(SEXP ASEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type A(ASEXP);
    rcpp_result_gen = Rcpp::wrap(mat_cov2cor(A));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_Riemann_basic_pdist", (DL_FUNC) &_Riemann_basic_pdist, 3},
    {"_Riemann_basic_pdist2", (DL_FUNC) &_Riemann_basic_pdist2, 4},
    {"_Riemann_basic_interpolate", (DL_FUNC) &_Riemann_basic_interpolate, 5},
    {"_Riemann_inference_mean_intrinsic", (DL_FUNC) &_Riemann_inference_mean_intrinsic, 5},
    {"_Riemann_inference_mean_extrinsic", (DL_FUNC) &_Riemann_inference_mean_extrinsic, 5},
    {"_Riemann_inference_median_intrinsic", (DL_FUNC) &_Riemann_inference_median_intrinsic, 5},
    {"_Riemann_inference_median_extrinsic", (DL_FUNC) &_Riemann_inference_median_extrinsic, 5},
    {"_Riemann_clustering_nmshift", (DL_FUNC) &_Riemann_clustering_nmshift, 5},
    {"_Riemann_clustering_kmeans_lloyd", (DL_FUNC) &_Riemann_clustering_kmeans_lloyd, 6},
    {"_Riemann_clustering_kmeans_macqueen", (DL_FUNC) &_Riemann_clustering_kmeans_macqueen, 6},
    {"_Riemann_clustering_clrq", (DL_FUNC) &_Riemann_clustering_clrq, 5},
    {"_Riemann_clustering_sup_intrinsic", (DL_FUNC) &_Riemann_clustering_sup_intrinsic, 6},
    {"_Riemann_clustering_kmeans18B", (DL_FUNC) &_Riemann_clustering_kmeans18B, 6},
    {"_Riemann_cpp_scNJW", (DL_FUNC) &_Riemann_cpp_scNJW, 5},
    {"_Riemann_cpp_scUL", (DL_FUNC) &_Riemann_cpp_scUL, 5},
    {"_Riemann_cpp_scSM", (DL_FUNC) &_Riemann_cpp_scSM, 5},
    {"_Riemann_cpp_sc05Z", (DL_FUNC) &_Riemann_cpp_sc05Z, 5},
    {"_Riemann_visualize_pga", (DL_FUNC) &_Riemann_visualize_pga, 2},
    {"_Riemann_visualize_kpca", (DL_FUNC) &_Riemann_visualize_kpca, 4},
    {"_Riemann_visualize_isomap", (DL_FUNC) &_Riemann_visualize_isomap, 4},
    {"_Riemann_learning_seb", (DL_FUNC) &_Riemann_learning_seb, 5},
    {"_Riemann_learning_rmml", (DL_FUNC) &_Riemann_learning_rmml, 4},
    {"_Riemann_learning_coreset18B", (DL_FUNC) &_Riemann_learning_coreset18B, 6},
    {"_Riemann_runif_sphere", (DL_FUNC) &_Riemann_runif_sphere, 2},
    {"_Riemann_runif_stiefel", (DL_FUNC) &_Riemann_runif_stiefel, 3},
    {"_Riemann_mat_rank", (DL_FUNC) &_Riemann_mat_rank, 1},
    {"_Riemann_mat_symm", (DL_FUNC) &_Riemann_mat_symm, 2},
    {"_Riemann_mat_diaghalf", (DL_FUNC) &_Riemann_mat_diaghalf, 1},
    {"_Riemann_mat_diaginvhalf", (DL_FUNC) &_Riemann_mat_diaginvhalf, 1},
    {"_Riemann_mat_cov2cor", (DL_FUNC) &_Riemann_mat_cov2cor, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_Riemann(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
