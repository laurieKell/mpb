#' biodyns creates a list of biodyn objects
#' 
#' @slot object
#' 
setMethod('biodyns', signature(object='biodyn'), 
          function(object) {
  lst <- c(object, list(...))
  biodyns(lst)})

setMethod('biodyns', signature(object='missing'),
          function(object,...) {
            # empty
            if(missing(...)){
              new('biodyns')
              # or not
            } else {
              args <- list(...)
              object <- args[!names(args)%in%c('names', 'desc', 'lock')]
              args <- args[!names(args)%in%names(object)]
              do.call('biodyns',  c(list(object=object), args))
            }
          }
)

#' Class biodyns
#' 
#' \code{biodyns} is a class that extends \code{list} but implements a set of
#' features that give a little more structure to list objects. First the
#' elements of \code{biodyns} must all be of the same class \code{biodyn}. 
#' Second it implements lock mechanism that, when turned on, does not allow 
#' the user to increase
#' or decrease the object length.
#' 
#' @name biodyns
#' @aliases biodyns-class biodyns biodyns-methods biodyns,ANY-method
#' biodyns,missing-method biodyns,list-method
#' @docType class
#' @section Slots: \describe{
#'   \item{.Data}{The data. \code{list}.}
#'   \item{names}{Names of the list elements. \code{character}.}
#'   \item{desc}{Description of the object. \code{character}.}
#'   \item{lock}{Lock mechanism, if turned on the length of the list can not be
#'     modified by adding or removing elements. \code{logical}.} }
#' @template biodyns-constructors
#' @author The FLR Team
#' @seealso \link[base]{[}, \link[base]{[<-}, \link[base]{[[<-},
#' \link[base]{$<-}, \link[methods]{coerce}, \link[base]{lapply},
#' \link[stats]{window}, \link[base]{list}
#' @keywords classes
#' @examples
#' \dontrun{
#' bds <- biodyns("1"=sim(),"2"=biodyn())}
#' 
setMethod('biodyns', signature(object='list'),    
          function(object, ...) {
            
            args <- list(...)
            
            # names in args, ... 
            if('names' %in% names(args)) {
              names <- args[['names']]
            } else {
              # ... or in object,
              if(!is.null(names(object))) {
                names <- names(object)
                # ... or in elements, ...
              } else {
                names <- unlist(lapply(object, name))
                # ... or 1:n
                idx <- names == 'NA' | names == ''
                if(any(idx))
                  names[idx] <- as.character(length(names))[idx]
              }
            }
            
            # desc & lock
            args <- c(list(Class='biodyns', .Data=object, names=names),
                      args[!names(args)%in%'names'])
            
            return(
              do.call('new', args)
            )
            
          }) 

#setMethod('biodyns', signature(object='aspics'), 
as2bs<-function(object, ...) {
  bd=llply(object,function(x){
    b=as(x,"biodyn")
    #b=aspic2biodyn(x)
    
    b@control[substr(dimnames(b@control)[[1]],1,1)=="q",3]=
      x@params[substr(dimnames(x@params)[[1]],1,1)=="q"]
    b@control[substr(dimnames(b@control)[[1]],1,1)=="q",2]=
      b@control[substr(dimnames(b@control)[[1]],1,1)=="q",3]*0.1
    b@control[substr(dimnames(b@control)[[1]],1,1)=="q",4]=
      b@control[substr(dimnames(b@control)[[1]],1,1)=="q",3]*1.0
    
    b@control["p","phase"]=-10
    
    #b@control["r",3:4]=b@control["r",3:4]*2.0 
    b@control["k",3:4]=b@control["k",3:4]*2.0
    b@control[c("r","k","q1","sigma1"),"phase"]=c(4,3,2,1)

    u=index(x,FALSE)
    
    b=fit(b,u)
    b})
  
  biodyns(bd)
  }
#)
   

biodyns <- setClass('biodyns', contains='FLComps',
                    validity=function(object) {
                      # All items are biodyn
                      if(!all(unlist(lapply(object, is, 'biodyn'))))
                        return('Components must be biodyn')  
                      
                      return(TRUE)})