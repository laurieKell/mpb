utils::globalVariables('feasibleCpp')

#' feasible
#' 
#' Checks that parameters are feasible given data.
#' 
#' @description Checks that the parameters can be taken given the catch time series
#'       
#' @param object an \code{FLPar} or \code{data.frame} object with parameters for the production function, r, k, p and b0.
#' @param catch an \code{FLQuant} with a time series of catch, iters must be equal to 1
#' @param min the minimum permissable population level, used to check that the catch can be taken.
#' @param ... any other parameters
#' 
#' @return a \code{FLPar} a subset of params with parameter values that can explain the catch
#' 
#' @aliases 
#' feasible,FLPar,FLQuant-method
#' feasible,biodyn,missing-method
#' feasible,data.frame,FLQuant-method
#' 
#' @export
#' @docType methods
#' @rdname feasible
#' 
#' @examples
#' \dontrun{
#' params=feasibleFn(catch,params)
#' }
setMethod('feasible',  signature(object="FLPar",catch="FLQuant"), 
   function(object,catch,min=0.01){
     feasibleFn(params=object,catch=catch,min=min)})

setMethod('feasible',  signature(object="data.frame",catch="FLQuant"), 
          function(object,catch,min=0.01){
            feasibleFn(x=as(object,"FLPar"),y=catch,min=min)})


setMethod('feasible',  signature(object='biodyn',catch="missing"), 
          function(object,min=0.01)
            feasibleFn(params=params(object),catch=catch,min))

feasibleFn<-function(params,catch,min=0.01){          
  
  flg=feasibleCpp(catch,t(params@.Data))>0

  params[,flg]}

setMethod('grid',  signature(object='FLQuant'), 
          function(object, 
                   r=seq(0.01, 1.0, length.out=101),
                   k=seq(1,    100, length.out= 11)*signif(mean(catch),1),
                   p=0.001,b0=0.95,
                   maxk=max(k),mult=1.0,step=signif(max(k)/100,1))
            gridFn(object,
                   r=r,k=k,p=p,b0=b0,
                   maxk=maxk,mult=mult,step=step)) 
              
setMethod('grid',  signature(object='biodyn'), 
          function(object,
                   r=seq(0.01, 1.0, length.out=101),
                   k=seq(1,    100, length.out= 11)*signif(mean(catch(object)),1),
                   p=0.001,b0=0.95,
                   maxk=max(k),mult=1.0,step=signif(max(k)/100,1))
            gridFn(catch(object),
                   r=r,k=k,p=p,b0=b0,
                   maxk=maxk,mult=mult,step=step)) 

gridFn<-function(catch,
                 r=seq(0.01, 1.0, length.out=101),
                 k=seq(1,    100, length.out= 11)*signif(mean(catch),1),
                 p=0.001,b0=0.95,
                 maxk=max(k),mult=1.0,step=signif(max(k)/100,1)) {
  
  par=as(data.frame(expand.grid(r=r,k=k),p=p,b0=b0),"FLPar")
  par=feasible(par,catch)
  fsb=with(model.frame(par), coefficients(lm(1/k~r)))
  names(fsb)=c("intercept","slope")
  
  grd=data.frame(r=r,p=p,b0=b0,maxk=maxk,mult=mult,step=step)
  grd=mdply(grd, function(r,p,b0,maxk,mult,step)
    data.frame(k=seq(1/(r*fsb["slope"]*mult+fsb["intercept"]),
               maxk,
               step)))
  as(grd[,c("r","k","p","b0")],"FLPar")}
