#' aspics a class that contains a list of ASPIC biomass dynamic stock assessment model classes
#' 
#' @slot object
#' 
setClass("aspics",
         representation(
           "FLlst"))

setMethod("aspics", signature(object="missing"),
  function(object,...) {
    # empty
    if(missing(...)){
	  	new("aspics")
    # or not
  	} else {
      args <- list(...)
      object <- args[!names(args)%in%c('names', 'desc', 'lock')]
      args <- args[!names(args)%in%names(object)]
      do.call('aspics',  c(list(object=object), args))
	  }
  }
)

#' Class aspics
#' 
#' \code{aspics} is a class that extends \code{list} but implements a set of
#' features that give a little more structure to list objects. First the
#' elements of \code{aspics} must all be of the same class \code{biodyn}. 
#' Second it implements lock mechanism that, when turned on, does not allow 
#' the user to increase
#' or decrease the object length.
#' 
#' @name aspics
#' @aliases aspics-class aspics aspics-methods aspics,ANY-method
#' aspics,missing-method aspics,list-method
#' @docType class
#' @section Slots: \describe{
#'   \item{.Data}{The data. \code{list}.}
#'   \item{names}{Names of the list elements. \code{character}.}
#'   \item{desc}{Description of the object. \code{character}.}
#'   \item{lock}{Lock mechanism, if turned on the length of the list can not be
#'     modified by adding or removing elements. \code{logical}.} }
#' @template aspics-constructors
#' @author The FLR Team
#' @seealso \link[base]{[}, \link[base]{[<-}, \link[base]{[[<-},
#' \link[base]{$<-}, \link[methods]{coerce}, \link[base]{lapply},
#' \link[stats]{window}, \link[base]{list}
#' @keywords classes
#' @examples
#' \dontrun{
#' asp <- aspics("1"=aspic(),"2"=aspic())}
#' 
setMethod("aspics", signature(object="list"),
  function(object, ...) {
    
    args <- list(...)
    
    # names in args, ... 
    if("names" %in% names(args)) {
      names <- args[['names']]
    } else {
    # ... or in object,
      if(!is.null(names(object))) {
        names <- names(object)
    # ... or in elements, ...
      } else {
        names <- unlist(lapply(object, name))
        # ... or 1:n
        idx <- names == "NA" | names == ""
        if(any(idx))
          names[idx] <- as.character(length(names))[idx]
      }
    }

    # desc & lock
    args <- c(list(Class="aspics", .Data=object, names=names),
      args[!names(args) %in% 'names'])

    return(
      do.call('new', args)
      )}) 

setMethod("aspics", signature(object="character"),
          function(object) {
           
          aspics(mlply(object,aspic))})
