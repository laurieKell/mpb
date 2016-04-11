# FLQuantJKs {{{
# validity
vFLQs <- function(object){
  # Make sure the list contains all items of the same class
  for(i in 1:length(object)){
    if(!is(object[[i]], "FLQuant")) stop("Components must be FLQuant")	
  }
  # Everything is fine
  return(TRUE)
}

#' Class FLQuantJKs
#' 
#' \code{FLQuantJKs} is a \code{list} of \code{FLQuant} objects.
#' It is very similar to the standard \code{list} class.
#' It implements a lock mechanism that, when turned on, does
#' not allow the user to increase or decrease the object length.
#' The elements of \code{FLQuantJKs} must all be of class  \code{FLQuant}. 
#' 
#' @name FLQuantJKs
#' @aliases FLQuantJKs-class FLQuantJKs FLQuantJKs-methods FLQuantJKs,ANY-method
#' FLQuantJKs,missing-method FLQuantJKs,list-method FLQuantJKs,FLQuantJKs-method
#' @docType class
#' @section Slots: \describe{
#'     \item{.Data}{The data. \code{list}.}
#'     \item{names}{Names of the list elements. \code{character}.}
#'     \item{desc}{Description of the object. \code{character}.}
#'     \item{lock}{Lock mechanism, if turned on the length of the list can not be modified by adding or removing elements. \code{logical}.}
#' }
#' @author The FLR Team
#' @seealso \link[base]{*}, \link[methods]{Arith}, \link[base]{as.data.frame},
#' \link{bubbles}, \link{catch<-}, \link{iter}, \link[stats]{model.frame},
#' \link[methods]{show}, \link[base]{summary}, \link[lattice]{xyplot},
#' \link{FLlst}, \link[base]{list}
#' @keywords classes
# class
setClass("FLQuantJKs", contains="FLlst",
         validity=vFLQs
)

# constructor
setGeneric("FLQuantJKs", function(object, ...){
  standardGeneric("FLQuantJKs")
}
)

setMethod("FLQuantJKs", signature(object="ANY"), function(object, ...){
  lst1 <- list(...)
  nlst <- length(lst1)
  lst <- list()
  length(lst) <- nlst + 1
  lst[[1]] <- object
  lst[-1] <- lst1
  new("FLQuantJKs", lst)
})


setMethod("FLQuantJKs", signature(object="FLComp"),
          function(object, ...) {
            
            args <- list(...)
            
            # SPLIT into list if a character vector 
            if(length(args) == 1 & length(args[[1]]) > 1)
              args <- as.list(args[[1]])
            
            # CHECK args are char or function
            chr <- unlist(lapply(args, function(x) is(x, 'character')))
            fun <- unlist(lapply(args, function(x) is(x, 'function')))
            
            if(sum(chr + fun) != length(args))
              stop("Arguments in ... must be of class 'character' or 'function'")
            
            # CHECK function elements have names
            if(any(names(args[fun]) == ""))
              stop("Function calls must be named, e.g. catch=catch")
            
            # GET names
            nms <- names(args)
            nms[chr] <- unlist(args[chr])
            
            # DO.CALL list elements
            res <- lapply(args, function(x) do.call(x, args=list(object)))
            
            # ASSIGN names
            names(res) <- nms
            
            return(new("FLQuantJKs", res))
          })

setMethod("FLQuantJKs", "missing", function(...){
  if(missing(...)){
    new("FLQuantJKs")
  } else { 
    lst <- list(...)
    new("FLQuantJKs", lst)
  }
})

setMethod("FLQuantJKs", "list", function(object){
  res=new("FLQuantJKs", object)
  names(res)=names(object)
  res
})

setMethod("FLQuantJKs", "FLQuantJKs", function(object){
  return(object)
}) # }}}
