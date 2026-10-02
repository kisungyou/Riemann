#' Method Contracts and Validation Status
#'
#' Riemann distinguishes validated input representations, available geometry
#' primitives, algorithm contracts, and statistical calibration. A wrapped point
#' or available distance does not establish every algorithm's assumptions.
#'
#' @details The installed file \file{method-contracts.csv} contains one disposition
#'   for every exported function and registered S3 method, including its domain,
#'   assumptions, restrictions, and evidence scope. Read it with the example below.
#'   The four possible dispositions are:
#'   \describe{
#'     \item{core}{The consolidated input, geometry, summary, clustering, tangent
#'       PCA, or scalar regression contract, within its documented capabilities.}
#'     \item{validated_restricted}{A specific legacy contract checked against
#'       independent formulas or reference fixtures; the stated restrictions apply.}
#'     \item{experimental}{A retained legacy interface whose full numerical,
#'       mathematical, or statistical contract has not been independently verified.}
#'     \item{disabled}{An unavailable interface with a documented incompatible
#'       implementation. Current geometry restrictions can disable combinations
#'       without disabling an entire exported function.}
#'   }
#'
#'   A test-file reference is an evidence pointer, not branch coverage or a proof.
#'   No disposition guarantees a global optimum, an identifiable model, finite
#'   sample calibration, or validity for every manifold. Experimental methods are
#'   retained for compatibility and exploration; their source references and
#'   examples should not be read as completed validation.
#'
#'   Common distinctions are particularly relevant: nonnegative graph affinities
#'   or regression weights need not form a positive-semidefinite kernel; a metric
#'   dissimilarity need not be Euclidean; a local tangent representation need not
#'   solve exact nonlinear principal geodesic analysis; and a permutation test
#'   requires exchangeability under its specified null. Rayleigh/Bingham and
#'   Frechet mean/variance tests need not detect all distributional alternatives.
#'
#'   \code{riem.capabilities()} records primitive availability and geometry status.
#'   The table here separately records method status. The package maintenance
#'   inventory records source and dependency ownership; the release/check logs
#'   supply actual execution evidence for each environment.
#'
#' @examples
#' contracts <- read.csv(system.file("method-contracts.csv", package = "Riemann"))
#' contracts[contracts$name %in% c("riem.kpca", "riem.fanova"),
#'           c("name", "disposition", "limitations")]
#'
#' @seealso \code{\link{riem.capabilities}}, \code{\link{riem.geometry}}
#' @name riem-method-contracts
#' @concept validation
NULL
