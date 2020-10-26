#include <RcppArmadillo.h>
#include "riemann_src.h"
#include <map>

using namespace Rcpp;
using namespace arma;
using namespace std;





// SPECIAL FUNCTIONS ON SPD MANIFOLD
// (01) spd_dist    : compute distance of two SPD matrices
//      spd_pdist   : pairwise distances




// (01) spd_dist  : compute distance of two SPD matrices -----------------------
// [[Rcpp::export]]
double spd_dist(arma::mat X, arma::mat Y, std::string geometry){
  double output = 0.0;
  if (geometry=="airm"){                 // 1. AIRM
    output = riem_dist("spd",X,Y);
  } else if (geometry=="lerm"){          // 2. LERM
    output = riem_distext("spd",X,Y);
  } else if (geometry=="jeffrey"){       // 3. Jeffrey
    double term1 = arma::trace(arma::solve(X,Y))/2.0;
    double term2 = arma::trace(arma::solve(Y,X))/2.0;
    double term3 = static_cast<double>(X.n_rows);
    output = term1 + term2 - term3;
  } else if (geometry=="stein"){        // 4. Stein
    output = std::sqrt(std::log(arma::det((X+Y)/2.0)) - 0.5*std::log(arma::det(X*Y)));
  }
  return(output);
}
// [[Rcpp::export]]
arma::mat spd_pdist(Rcpp::List& data, std::string geometry){
  // PRELIMINARY
  arma::mat exmat  = Rcpp::as<arma::mat>(data[0]);
  int p = exmat.n_rows;
  int N = data.size();
  double pp = static_cast<double>(p);
  
  arma::cube mydata(p,p,N,fill::zeros);
  for (int n=0; n<N; n++){
    mydata.slice(n) = Rcpp::as<arma::mat>(data[n]);
  }
  
  // COMPUTE
  arma::mat distance(N,N,fill::zeros);
  if (geometry=="jeffrey"){                         // special case : 3. Jeffrey
    arma::cube invdata(p,p,N,fill::zeros);
    for (int n=0; n<N; n++){
      invdata.slice(n) = arma::inv_sympd(mydata.slice(n));
    }
    for (int i=0; i<(N-1); i++){
      for (int j=(i+1); j<N; j++){
        distance(i,j) = (arma::trace(invdata.slice(i)*mydata.slice(j))*0.5) + (arma::trace(invdata.slice(j)*mydata.slice(i))*0.5) - pp;
        distance(j,i) = distance(i,j);
      }
    }
  } else if (geometry=="stein"){                      // special case : 4. Stein
    arma::vec logdetX(N,fill::zeros);
    for (int n=0; n<N; n++){
      logdetX(n) = std::log(arma::det(mydata.slice(n)));
    }
    for (int i=0; i<(N-1); i++){
      for (int j=(i+1); j<N; j++){
        distance(i,j) = std::sqrt(std::log(arma::det((mydata.slice(i)+mydata.slice(j))/2.0)) - 0.5*(logdetX(i) + logdetX(j)));
        distance(j,i) = distance(i,j);
      }
    }
    
  } else {                                                 // all the other case
    for (int i=0; i<(N-1); i++){
      for (int j=(i+1); j<N; j++){
        distance(i,j) = spd_dist(mydata.slice(i), mydata.slice(j), geometry);
        distance(j,i) = distance(i,j);
      }
    } 
  }
  
  // RETURN
  return(distance);
}