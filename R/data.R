#' An example Operating Model, based on North Atlantic Albacore
#'
#' A dataset containing the historic stock estimates, derived from the ICCAT multifan-CL assessment 
#'
#' @format An \code{FLStock} object
#' \describe{
#'   \item{stock.n}{numbers-at.age, \code{FLQuant}}
#'   \item{stock.wt}{mass-at-age, \code{FLQuant}}
#'   ...
#' }
#' @source \url{http://http://iccat.int/Documents/CVSP/CV070_2014/n_4/CV070041365.pdf/}
"om"

#' Reference points and expected dynamics, for Operating Model based on North Atlantic Albacore
#'
#' A dataset containing the average estimates of biological parameters and selection patterns, 
#' derived from the ICCAT multifan-CL assessment 
#'
#' @format An \code{FLBRP} object
#' \describe{
#'   \item{stock.wt}{mass-at-age, \code{FLQuant}}
#'   ...
#' }
#' @source \url{http://http://iccat.int/Documents/CVSP/CV070_2014/n_4/CV070041365.pdf/}
"eql"