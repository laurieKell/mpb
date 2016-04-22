#' biodyns
#' 
#' @description Create a list of biodyn objects
#' 
#' #'  biodyns-method biodyns,biodyn-method  biodyns,list-method  biodyns,missing-method biodyns,aspics-method
#' @name biodyns
#' @param object can be \code{biodyn} object or a \code{list} of \code{biodyn} objects
#' @param ... additional \code{biodyn} objects
#' 
#' @return \code{biodyns} object
#' 
#' @export
#' @rdname biodyns
#' 
#' 
#' @examples 
#' \dontrun{
#'    biodyns(biodyn())
#'    }

setGeneric('biodyns',   function(object,...)  standardGeneric('biodyns'))

# constructor
setMethod('biodyns', signature(object='biodyn'), function(object) {
  lst <- c(object, list(...))
  biodyns(lst)})

setMethod('biodyns', signature(object='missing'),function(object,...) {
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

setMethod('biodyns', signature(object='list'),    function(object, ...) {
            
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
   

setGeneric('biodyns', function(object, ...) standardGeneric('biodyns'))
biodyns <- setClass('biodyns', contains='FLComps',
                    validity=function(object) {
                      # All items are biodyn
                      if(!all(unlist(lapply(object, is, 'biodyn'))))
                        return('Components must be biodyn')  
                      
                      return(TRUE)})