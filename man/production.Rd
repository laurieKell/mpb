% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FLBRP-production.R
\name{productionFn}
\alias{productionFn}
\alias{production,biodyn-method}
\alias{production,FLStock,missing}
\alias{production,biodyn,FLQuant}
\title{production}
\usage{
productionFn(
  object,
  biomass = "missing",
  slots = c("landings.sel", "discards.sel", "stock.wt", "landings.wt", "discards.wt",
    "m", "mat", "harvest.spwn", "m.spwn")
)
}
\arguments{
\item{object}{either \emph{biodyn} or \emph{FLBRP} class}

\item{biomass}{an \emph{FLQuant}}

\item{...}{any other parameter}
}
\description{
Estimates production for a given biomass
}
