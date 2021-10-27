#' Supported Geometries on SPD Manifold
#' 
#' SPD manifold is a well-studied space in that there have been many geometries 
#' proposed on the space. For special functions on under SPD category, this 
#' function finds whether there exists a matching name that is currently 
#' supported in \pkg{Riemann}. If there is none, it will return an error message.
#' 
#' @param geometry name of supported geometries, including
#' \describe{
#' \item{AIRM}{Affine-Invariant Riemannian Metric.}
#' \item{LERM}{Log-Euclidean Riemannian Metric.}
#' \item{Jeffrey}{Jeffrey's divergence.}
#' \item{Stein}{Stein's metric.}
#' \item{Wasserstein}{2-Wasserstein geometry.}
#' }
#' 
#' @return a matching name in lower-case.
#' 
#' @examples 
#' # it just returns a small-letter string.
#' mygeom = spd.geometry("stein")
#' 
#' @concept spd
#' @export
spd.geometry <- function(geometry){
  name.all = tolower(c("AIRM","LERM","Jeffrey","Stein","Wasserstein"))
  name.tgt = tolower(geometry)
  
  if (name.tgt %in% name.all){
    return(match.arg(name.tgt, name.all))
  } else {
    stop(paste0("* spd.geometry : ",geometry," is not currently supported."))
  }
}


# auxiliary : spd.wrap3d for convenient use -------------------------------
#' @keywords internal
#' @noRd
spd.wrap3d <- function(riemdata){
  N = length(riemdata)
  p = base::nrow(riemdata[[1]])
  
  output = array(0,c(p,p,N))
  for (n in 1:N){
    output[,,n] = riemdata[[n]]
  }
  return(output)
}
