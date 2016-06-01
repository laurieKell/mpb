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
#' @rdname feasible
#' @examples
#' \dontrun{
#'     data(asp)
#'     asp=boot(asp)}
setGeneric('feasible',    function(object,catch,...)  standardGeneric('feasible'))
setGeneric('setFeasible', function(object,catch,...)  standardGeneric('setFeasible'))

setMethod('feasible', signature(object='biodyn',catch='missing'), function(object) feasibleFn(params(object),catch(object))
setMethod('feasible', signature(object='FLPar', catch='FLQuant'), function(object) feasibleFn(object,catch))

setMethod('setFeasible<-', signature(object='FLPar', value='FLQuant'), function(object) feasibleFn2(object,catch))
setMethod('setFeasible<-', signature(object='biodyn',value='missing'), function(object) feasibleFn3(params(object),catch(object))
          
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
setGeneric('nll', function(object,catch,...)  standardGeneric('nll'))

setMethod('nll', signature(object='biodyn',index='FLQuant'),  
          function(object) nllFn(params(object),catch(object),index)
setMethod('nll', signature(object='biodyn',index='FLQuants'),  
          function(object) nllFn(params(object),catch(object),index)
setMethod('nll', signature(object='FLPar',index='FLQuant'),  
          function(object) nllFn(params(object),catch=catch,index)
setMethod('nll', signature(object='FLPar',index='FLQuants'),  
          function(object) nllFn(params(object),catch=catch,index)
                    
nllFn<-function(object,catch,index){
  
       }          