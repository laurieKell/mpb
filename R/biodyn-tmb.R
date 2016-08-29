#' feasible, biodyn.
#'
#' @description 
#' Checks that parameters are feasible by ensuring catch can be taken given parameters. 
#'
#' @param object; a \code{biodyn} object 
#' @seealso \code{\link{biodyn}}
#' 
#' @aliases feasible-method feasible,biodyn-method
#' 
#' @export
#' 
#'@importFrom Rcpp evalCpp
#'@useDynLib mpb
#'
#' @rdname feasible
#' @examples
#' \dontrun{
#'     dat(asap)
#'     asp=boot(asp)}
setGeneric('feasible',      function(object,catch,...)  standardGeneric('feasible'))
setGeneric('setFeasible',   function(object,catch,...)  standardGeneric('setFeasible'))
setGeneric('setFeasible<-', function(object,value,...)  standardGeneric('setFeasible<-'))

setMethod('feasible', signature(object='biodyn',catch='missing'), function(object) feasibleFn(params(object),catch(object)))
setMethod('feasible', signature(object='FLPar', catch='FLQuant'), function(object,catch) feasibleFn(object,catch))

#setMethod('setFeasible<-', signature(object='FLPar', value='FLQuant'), function(object,value) feasibleFn2(object,catch))
#setMethod('setFeasible<-', signature(object='biodyn',value='missing'), function(object,value) feasibleFn3(params(object),catch(object))

setGeneric('fitmb',      function(object,index,...)  standardGeneric('fitmb'))


feasibleFn2<-function(params,catch){
  
  res=testGrid(catch,t(params@.Data))
  apply(res,6,function(x) all(x>=0.001))
  }
  
feasibleFn2<-function(params,catch){
  
  res=testGrid(catch,t(params@.Data))
  pararms[,apply(res,6,function(x) all(x>=0.001))]
  }

feasibleFn3<-function(object,catch){
  
  res=testGrid(catch,t(object@params@.Data))
  flag=apply(res,6,function(x) all(x>=0.001))

  object@params=object@pararms[,flag]
  }

#' nll, biodyn.
#'
#' @description 
#' Checks that parameters are feasible by ensuring catch can be taken given parameters. 
#'
#' @param object; a \code{biodyn} object 
#' @seealso \code{\link{biodyn}}
#' 
#' @aliases feasible-method feasible,biodyn-method
#' 
#' @export
#' 
#' @rdname feasible
#' @examples
#' \dontrun{
#'     data(asp)
#'     asp=boot(asp)}
setGeneric('nll', function(object,index,...)  standardGeneric('nll'))

setMethod('nll', signature(object='biodyn',index='FLQuant'),  
          function(object,index) nllFn(params(object),catch(object),index))
setMethod('nll', signature(object='biodyn',index='FLQuants'),  
          function(object,index) nllFn(params(object),catch(object),index))
setMethod('nll', signature(object='FLPar',index='FLQuant'),  
          function(object,index) nllFn(params(object),catch=catch,index))
setMethod('nll', signature(object='FLPar',index='FLQuants'),  
          function(object,index) nllFn(params(object),catch=catch,index))
  
#fitFn<-function(object,index){
#   
#   cat(run," ")
#   
#   ctrl=biodyn()@control[drop=T]
#   ctrl[1:4,1]=c(1,1,-1,-1)
#   ctrl[1:4,3]=unlist(c(runs[run,c("r","k","p","b0")]))
#   #ctrl[1,2]=0.1
#   #ctrl[2,2]=2.5e9
#   
#   ctrl[1:4,2]=ctrl[1:4,3]*0.01
#   ctrl[1:4,4]=ctrl[1:4,3]*10.0
#   
#   if (any(ctrl[,1]>0))
#     ctrl[ctrl[,1]>0,3]=logit(ctrl[ctrl[,1]>0,2],ctrl[ctrl[,1]>0,3],ctrl[ctrl[,1]>0,4])
#   
#   ctc=c(iter(window(catch(o),  start=runs[run,"yr0"]),as.numeric(runs[run,"scenario"])))
#   idx=c(iter(window(u[["bio"]],start=runs[run,"yr0"]),as.numeric(runs[run,"scenario"])))
#   idx=matrix(idx,c(length(idx),1))
#   
#   f=MakeADFun(data      =list(ctc=ctc,idx=idx,ctl=ctrl),
#               parameters=list(par=unlist(c(ctrl[,3]))),
#               DLL = "pellaTmb")
#   
#   hat=optimx(f$par,f$fn,control=list(trace=0,fnscale=-1),method="BFGS")
#   res=hat[c(1:5)]
#   names(res)=c("r","k","p","b0","ll")
#   
#   ctrl[,3][]=unlist(c(res[1:4]))
#   
#   if (any(ctrl[,1]>0))
#     res[-5][ctrl[,1]>0]=invLogit(ctrl[ctrl[,1]>0,2],ctrl[ctrl[,1]>0,3],ctrl[ctrl[,1]>0,4])
#   
#   res
# })