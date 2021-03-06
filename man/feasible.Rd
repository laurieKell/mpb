% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/biodyn-feasible.R
\docType{methods}
\name{feasible}
\alias{feasible}
\alias{feasible,FLPar,FLQuant-method}
\alias{feasible,biodyn,missing-method}
\alias{feasible,data.frame,FLQuant-method}
\title{feasible}
\usage{
\S4method{feasible}{FLPar,FLQuant}(object, catch, min = 0.01)
}
\arguments{
\item{object}{an \code{FLPar} or \code{data.frame} object with parameters for the production function, r, k, p and b0.}

\item{catch}{an \code{FLQuant} with a time series of catch, iters must be equal to 1}

\item{min}{the minimum permissable population level, used to check that the catch can be taken.}

\item{...}{any other parameters}
}
\value{
a \code{FLPar} a subset of params with parameter values that can explain the catch
}
\description{
Checks that the parameters can be taken given the catch time series
}
\examples{
\dontrun{
params=feasibleFn(catch,params)
}
}
