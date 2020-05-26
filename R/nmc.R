setGeneric('nmc',    function(object,...) standardGeneric('nmc'))

#' @title nmc
#'
#' @description Calculates the number of iterations needed in Monte Carlo Simulation to achieve a given precision.
#' 
#' @param x \code{FLQuant} 
#' @param z \code{numeric}  Value of Z for a given confidence level for a normally distributed random variable, default is 1.96 for a 95% CI 
#' @param E \code{numeric}  The required percentage error of the mean
#' @param ... other arguments
#' 
#' @export
#' @rdname nmc
#' 
#' @return \code{FLQuant} with number of iters required by number of iters 
#' 
#' @examples
# \dontrun{
#'  data(ple4) 
#'  ssb=rlnorm(2000,log(ssb(ple4)),0.5)
#'  nmc(ssb) 
#'  }
#'  
nmcFn<-function(x,s,E=5,z=1.96){
  ((100*z*s)/(E*x))^2}

setMethod('nmc',   signature(object='FLQuant'),
          function(object,E=5,z=1.96){

            res=aaply(object,1:5,function(x) 
                   maply(10:length(x), function(i) nmcFn(x=mean(x[seq(i)]),s=var(x[seq(i)])^0.5,E=E,z=z)))
            
            names(dimnames(res))[length(dim(res))]="iter"
            dimnames(res)[[2]]=seq(10:dim(object)[6])
            dmn=dim(object)
            dmn[6]=dmn[6]-9
            res2=as.FLQuant(array(res,dmn))

    res2})
