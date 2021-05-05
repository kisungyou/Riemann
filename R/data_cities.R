#' Data : Populated Cities in the U.S.
#' 
#' As of January 2006, there are 60 cities in the contiguous U.S. with population size 
#' larger than \eqn{300000}. We extracted information of the cities from the data 
#' delivered by \pkg{maps} package. Variables \code{coord} and \code{cartesian} are 
#' two identical representations of locations, which can be mutually converted by 
#' \code{\link{sphere.convert}}.
#' 
#' @usage data(cities)
#' 
#' @examples
#' \donttest{
#' ## LOAD THE DATA AND WRAP AS RIEMOBJ
#' data(cities)
#' myriem = wrap.sphere(cities$cartesian)
#' 
#' ## COMPUTE INTRINSIC/EXTRINSIC MEANS
#' intmean = as.vector(riem.mean(myriem, geometry="intrinsic")$mean)
#' extmean = as.vector(riem.mean(myriem, geometry="extrinsic")$mean)
#' 
#' ## CONVERT TO GEOGRAPHIC COORDINATES (Lat/Lon)
#' geo.int = sphere.xyz2geo(intmean)
#' geo.ext = sphere.xyz2geo(extmean)
#' }
#' 
#' @format a named list containing\describe{
#' \item{names}{a length-\eqn{60} vector of city names.}
#' \item{coord}{a \eqn{(60\times 2)} matrix of latitude and longitude.}
#' \item{cartesian}{a \eqn{(60\times 3)} matrix of cartesian coordinates on the unit sphere.}
#' \item{population}{a length-\eqn{60} vector of cities' populations.}
#' }
#' 
#' @seealso \code{\link{wrap.sphere}}
#' @concept data
"cities"
# 
# # https://stackoverflow.com/questions/1185408/converting-from-longitude-latitude-to-cartesian-coordinates
# library(maps)
# data("us.cities")
# idbig    = which(us.cities$pop >= 300000)
# idremove = which(us.cities$name=="Honolulu HI")
# idbig    = setdiff(idbig, idremove)
# 
# dat.long = us.cities$long[idbig]
# dat.lati = us.cities$lat[idbig]
# 
# globe3d  = array(0,c(length(idbig),3))
# for (i in 1:length(idbig)){
#   globe3d[i,] = sphere.geo2xyz(dat.lati[i], dat.long[i])
# }
# colnames(globe3d) = c("x","y","z")
# 
# coords = cbind(dat.lati, dat.long); colnames(coords)=c("latitude","longitude")
# names  = us.cities$name[idbig]
# papa   = us.cities$pop[idbig]
# 
# cities = list()
# cities$names = names
# cities$coord = coords
# cities$cartesian = globe3d
# cities$population = papa
# library(ggplot2)
# library(usmap)