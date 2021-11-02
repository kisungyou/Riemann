#' Data : EEG Covariances for Event-Related Potentials
#' 
#' This dataset delivers 216 covariance matrices from EEG ERPs with 4 different 
#' known classes by types of sources. Among 60 channels, only 32 channels are 
#' taken and sample covariance matrix is computed for each participant. The 
#' data is taken from a Python library \href{https://mne.tools/stable/generated/mne.datasets.sample.data_path.html#mne.datasets.sample.data_path}{mne}'s 
#' sample data.
#' 
#' @usage data(ERP)
#' 
#' @examples
#' \donttest{
#' ## LOAD THE DATA AND WRAP AS RIEMOBJ
#' data(ERP)
#' myriem = wrap.spd(ERP$covariance)
#' }
#' 
#' @format a named list containing\describe{
#' \item{covariance}{an \eqn{(32\times 32\times 216)} array of covariance matrices.}
#' \item{label}{a length-\eqn{216} factor of 4 different classes.}
#' }
#' 
#' @seealso \code{\link{wrap.spd}}
#' @concept data
"ERP"

# data(ERP)
# wrapERP = wrap.spd(ERP$covariance)
# embed2d = riem.mds(wrapERP)