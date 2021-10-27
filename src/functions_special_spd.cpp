#include <RcppArmadillo.h>
#include "riemann_src.h"
#include <map>

using namespace Rcpp;
using namespace arma;
using namespace std;





// SPECIAL FUNCTIONS ON SPD MANIFOLD
// (01) src_spd_dist    : compute distance of two SPD matrices : ADD ON HERE!
//      src_spd_pdist   : pairwise distances




// (01) spd_dist  : compute distance of two SPD matrices -----------------------
double src_spd_dist(arma::mat X, arma::mat Y, std::string geometry){
  double output = 0.0;
  if (geometry=="airm"){                                              // 1. AIRM
    output = riem_dist("spd",X,Y);
  } else if (geometry=="lerm"){                                       // 2. LERM
    output = riem_distext("spd",X,Y);
  } else if (geometry=="jeffrey"){                                    // 3. Jeffrey
    double term1 = arma::trace(arma::solve(X,Y))/2.0;
    double term2 = arma::trace(arma::solve(Y,X))/2.0;
    double term3 = static_cast<double>(X.n_rows);
    output = term1 + term2 - term3;
  } else if (geometry=="stein"){                                      // 4. Stein
    output = std::sqrt(std::log(arma::det((X+Y)/2.0)) - 0.5*std::log(arma::det(X*Y)));
  } else if (geometry=="wasserstein"){                                // 5. Wasserstein
    arma::mat Xsqrt = arma::sqrtmat_sympd(X);
    output = std::sqrt(arma::trace(X+Y-2.0*arma::sqrtmat_sympd((Xsqrt*Y*Xsqrt))));
  }
  return(output);
}
// [[Rcpp::export]]
arma::mat src_spd_pdist(arma::cube &data, std::string geometry){
  // PRELIMINARY
  int p = data.n_rows;
  int N = data.n_slices;
  
  //arma::mat exmat  = Rcpp::as<arma::mat>(data[0]);
  //int p = exmat.n_rows;
  //int N = data.size();
  //double pp = static_cast<double>(p);
  
  // COMPUTE
  arma::mat distance(N,N,fill::zeros);
  for (int i=0; i<(N-1); i++){
    for (int j=(i+1); j<N; j++){
      distance(i,j) = src_spd_dist(data.slice(i), data.slice(j), geometry);
      distance(j,i) = distance(i,j);
    }
  }
  
  // RETURN
  return(distance);
}