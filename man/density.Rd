% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mixture_methods.R
\name{density}
\alias{density}
\title{Evaluate density of given data for a fitted model}
\usage{
density(object, newdata)
}
\arguments{
\item{object}{a fitted mixture model of \code{riemmix} class.}

\item{newdata}{data of \eqn{n} objects (vectors, matrices) that can be wrapped by one of \code{wrap.*} functions in the \pkg{Riemann} package.}
}
\value{
a length-\eqn{n} vector of class labels.
}
\description{
Compute density for a fitted mixture model.
}