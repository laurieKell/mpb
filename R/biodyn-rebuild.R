#' Calculate Rebuild Trajectories
#'
#' @description Projects stock rebuilding trajectories from different initial depletion levels
#'
#' @param object A biodyn object
#' @param target Target biomass (default = BMSY)
#' @param nInitial Number of initial depletion levels (default = 100)
#' @param growthRate Growth rate for depletion sequence (default = 0.3)
#' @param minVal Minimum depletion value (default = 1e-6)
#' @param maxVal Maximum depletion value (default = 1)
#' @param nx Number of interpolation points (default = 101)
#' @return A data frame with columns:
#'   \item{year}{Projection year}
#'   \item{initial}{Initial depletion level relative to BMSY}
#' @export
#' @examples
#' bd=biodyn(FLPar(r=0.5, k=1000, p=1))
#' rebuild_data=rebuild(bd)
setGeneric("rebuild", function(object, targetF=NULL, targetSSB=NULL,
                               nInitial=100, growthRate=0.25, minVal=1e-6, maxVal=1,
                               burnin=20, truncate=TRUE,...) {
  standardGeneric("rebuild")
})

setMethod("rebuild", signature(object="biodyn"), 
          function(object, target=refpts(object)["bmsy"], nInitial=100, 
                   growthRate=0.3, minVal=1e-6, maxVal=1, nx=101) {
            if (!is(object, "biodyn"))
              stop("object must be a biodyn object")
            
            if (!is.numeric(target) || length(target) != 1)
              stop("target must be a single numeric value")
            
            if (!all(sapply(list(nInitial, nx), function(x) is.numeric(x) && x > 0)))
              stop("nInitial, and nx must be positive integers")
            
            if (!all(sapply(list(growthRate, minVal, maxVal), is.numeric)))
              stop("growthRate, minVal, and maxVal must be numeric")
            
            if (minVal >= maxVal)
              stop("minVal must be less than maxVal")
            
            bmsy=c(refpts(object)["bmsy"])
            
            rtn=propagate(object, nInitial)
            
            # Create stock projection
            target_seq=c(target)*seq(minVal^growthRate, maxVal^growthRate, length.out=nInitial)^(1/growthRate)
            rtn@stock=FLQuant(rep(target_seq, each=dim(rtn)[2]), dimnames=dimnames(stock(rtn)))
            rtn=fwd(rtn, catch=catch(rtn)[,-1]%=%0.0)

            # Transform data
            dat=as.data.frame(stock(rtn), drop=TRUE)
            dat$initial=c(stock(rtn)[,1])[an(dat$iter)]
            dat=dat[,-2]
            
            # Interpolate results
            dat=as.data.frame(with(dat, interp(initial, data, year, yo=bmsy, duplicate="mean", nx=nx, jitter=1e-6)))[,c(3,1)]
            names(dat)=c("year", "initial")
            
            transform(dat,initial=initial/bmsy)})
