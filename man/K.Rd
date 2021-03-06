% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/biodyn-msy.R
\name{K}
\alias{K}
\title{Carrying capacity}
\usage{
K(msy, r = 0.6, p = 1, params = FLPar(r = r, p = p))
}
\arguments{
\item{msy}{a guess for MSY}

\item{r}{a guess for $r$ the population growth rate}

\item{p}{a guess for $p$ the shape parameter}

\item{params}{provide $r$ and $p$ as \code{FLPar}}
}
\value{
an \code{FLPar} object with an estimate for $k$
}
\description{
Calculates $k$ given msy, r and K for a Pella-Tomlinson biomass dynamic model
}
\examples{
\dontrun{
 K(5000,r=.6,p=1.0)
 }
}
\seealso{
\code{\link{msy}} and \code{\link{bmsy}}
}
