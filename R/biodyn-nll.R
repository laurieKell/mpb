#' nll
#' 
#' @description 
#' Checks that the parameters can be taken given the catch time series
#'       
#' @param object an \code{FLQuant} with a time series of catch, iters must be equal to 1
#' @param params an \code{FLPar} object with parameters for the production function, r, k, p and b0.
#' @param min, the minimum permissable population level, used to check that the catch can be taken.
#' @param ... any other parameters
#' @return a \code{FLPar} a subset of params with parameter values that can explain the catch
#' 
#' @aliases
#' nll,FLQuant,FLQuant,FLPar-method
#' nll,FLQuant,FLQuants,FLPar-method
#' nll,biodyn,FLQuant,missing-method
#' nll,biodyn,FLQuants,missing-method
#' 
#' @export
#' @docType methods
#' @rdname nll
#' 
#' @examples
#' \dontrun{
#' params=nllFn(catch,params)
#' }
setMethod('nll',  signature(object='FLQuant',index="FLQuant",params="FLPar"), 
   function(object,index,params,min=0.01)
     nllFn(object,index,params,min))

setMethod('nll',  signature(object='biodyn',index="FLQuant",params="missing"), 
    function(object,index,params,min=0.01)
      nllFn(catch(object),index,params(object),min))

setMethod('nll',  signature(object='FLQuant',index="FLQuants",params="FLPar"), 
    function(object,index,params,min=0.01){
            
    rtn=laply(index,function(x) nllFn(object,x,params,min))
    dmns=list(params="ll",index=names(index),iter=seq(dim(params)[2]))
    FLPar(array(c(rtn),laply(dmns,length),dmns))})

setMethod('nll',  signature(object='biodyn',index="FLQuants",params="missing"), 
   function(object,index,params,min=0.01){
            
   rtn=laply(index,function(x) nllFn(catch(object),x,params(object),min))
   dmns=list(params="ll",index=names(index),iter=seq(dim(params(object))[2]))
   FLPar(array(c(rtn),laply(dmns,length),dmns))})

nllFn<-function(ctc,idx,par,min=0.01){       
  
  full=FLQuant(NA,dimnames=dimnames(ctc))
  full[,dimnames(idx)$year]=idx
  full[full<=0]<-NA
  
  res=nllCpp(c(ctc),c(full),t(par@.Data))
  rtn=FLPar(array(res,c(1,dim(par)[2]),
                  dimnames=list(params="ll",iter=seq(dim(par)[2]))))
  rtn}