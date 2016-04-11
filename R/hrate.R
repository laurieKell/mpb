#' hrate
#'
#' @description harvest rate
#' 
#' @param object either \emph{biodyn} or \emph{aspic} classes
#' 
#' @rdname hrate
#' @export
#' 
#' @aliases hrate hrate,biodyn-method  hrate,aspic-method 
#' 
setGeneric('hrate',   function(object)  standardGeneric('hrate'))
setMethod('hrate', signature(object='biodyn'),
  function(object){
    yrs=dimnames(stock(object))$year
    yrs=yrs[yrs%in%dimnames(catch(object))$year]
    
    res=catch(object)[,yrs]/stock(object)[,yrs]
    units(res)="per y"
    
    res})