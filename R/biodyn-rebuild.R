#' Create biodyn Object from FLPar Parameters
#'
#' @description Creates a biodyn object initialized with biomass at BMSY and catch at MSY
#'
#' @param params An FLPar object containing parameters r, k, and p
#' @param nyrs Number of years for projection (default = 50)
#' @return A biodyn object initialized for the specified number of years
#' @export
#' @examples
#' params=FLPar(r=0.5, k=1000, p=1)
#' bd=biodyn(params)
setMethod("biodyn", signature(object="FLPar","missing"), 
          function(object, params=50) {
            if (!is(object, "FLPar"))
              stop("params must be an FLPar object")
            
            if (!all(c("r", "k", "p") %in% dimnames(object)$params))
              stop("params must contain r, k, and p")
            
            nyrs=params
            if (!is.numeric(nyrs) || nyrs <= 0)
              stop("nyrs must be a positive integer")
            
            rtn=new("biodyn")
            params(rtn)[c("r","k","p")]=object[c("r","k","p")]

            # Set initial states
            rtn@stock=window(FLQuant(), end=nyrs)
            rtn@stock[]=refpts(rtn)["bmsy"]
            rtn@catch=window(FLQuant(), end=nyrs)
            rtn@catch[]=refpts(rtn)["msy"]
            
            range(rtn)=unlist(dims(stock(rtn))[c("minyear","maxyear")])
            
            rtn=fwd(rtn,catch=catch[,-1])
            
            return(rtn)})

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
            dat=as.data.frame(with(dat, interp(initial, data, year, yo=bmsy, duplicate="mean", nx=nx)))[,c(3,1)]
            names(dat)=c("year", "initial")
            
            transform(dat,initial=initial/bmsy)})
