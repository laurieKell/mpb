utils::globalVariables(c("aaply"))

#' rand
#'
#' Simulates a \code{biodyn} object for a catch series, given the parameter estimates  in the \code{param} 
#' slot and variance covariance matrix 
#' http://young.physics.ucsc.edu/jackboot.pdf
#' 
#' @param n  \code{numeric} with number of simulations 
#' @param object \code{biodyn} 
#' @param ... other arguments
#' 
#' @return \code{biodyn} with estimates of stock based on catch time series
#' 
#' @aliases rand-method rand,numeric,FLPar,FLPar-method rand,numeric,FLQuant,FLQuant-method 
#' 
#' @export
#' @rdname rand
#'
#' @aliases rand-method rand,numeric,biodyn-method rand, numeric,FLPar,FLPar-method rand, numeric,FLQuant,FLQuant-method rand, numeric,biodyn,biodyn-method
#'
#' @examples
#' \dontrun{
#' #simulate an object with known properties
#' bd=sim()
#' bd=window(bd,end=49)
#' 
#' #simulate a proxy for stock abundance
#' cpue=(stock(bd)[,-dims(bd)$year]+stock(bd)[,-1])/2
#' cpue=rlnorm(1,log(cpue),.2)
#' 
#' #set parameters
#' setParams(bd) =cpue
#' setControl(bd)=params(bd)
#' control(bd)[3:4,"phase"]=-1
#' 
#' #fit
#' bd=fit(bd,cpue)
#' 
#' bdJK=fit(bd,jackknife(cpue))
#' 
#' 
#' bdRnd  =rand(100,bd,bdJK)
#' 
#' plot(rand(100,stock(bd)[,40:45],stock(bdJK)[,40:45]))
#' }
setGeneric('rand',   function(n,object,sim,...)    standardGeneric('rand'))

randFn<-function(n,object,sim){
  res=jackSummary(params(object),params(sim))
  
  cov=res$cov[modelParams(model(object)),modelParams(model(object))]
  object@params=propagate(object@params,n)
  object@catch =propagate(object@catch, n)
  object@catch[]=object@catch[,,,,,1]
  object@stock =propagate(object@stock, n)
  
  nms=dimnames(cov)[[1]]
  nms=nms[aaply(cov,1,function(x) !all(is.na(x)))]
    
  #FLParBug in drop so need c()
  params(object)[nms]=
    t(maply(seq(n),function(x) 
         mvrnorm(1,c(params(object)[nms,x,1]),cov[nms,nms,drop=T])))
  
  rtn=fwd(object,catch=object@catch)
  
  return(rtn)}

setMethod('rand', signature(n="numeric",object='biodyn',sim="biodyn"),  
           function(n,object,sim) randFn(n,object,sim))

setMethod('rand', signature(n="numeric",object='FLQuant',sim="FLQuant"),  
          function(n,object,sim){ 
              res=jackSummary(object,sim)
               
              if (is.null(res$cov))
                rnorm(n,res$mean,res$se)
              else
                mvrnorm(n,res$hat,res$cov)
              })
          
setMethod('rand', signature(n="numeric",object='FLPar',sim="FLPar"),  
          function(n,object,sim){ 
            res=jackSummary(object,sim)
            
            if (sum(dim(res$cov))>0)
              mvrnorm(n,res$hat,res$cov)
            else  
              rlnorm(n,res$hat,res$se)
          })
