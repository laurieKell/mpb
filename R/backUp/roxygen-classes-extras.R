#' @title biodyns
#' @description Create a list of biodyn objects
#' @name biodyns
#' @param object can be \code{biodyn} object or a \code{biodyn} of \code{biodyn} objects
#' @param ... additional \code{biodyn} objects
#' 
#' @return \code{biodyns} object
#' @export
#' @rdname biodynsConstructors2
#' 
#' @aliases 
#' biodyns-method 
#' biodyns,biodyn-method 
#' biodyns,missing-method 
#' biodyns,list-method
#' 
#' @examples 
#' \dontrun{
#' biodys(biodyns())
#' }


#' @title biodyns
#' 
#' @description 
#' The \code{biodyns} is a class that extends \code{list} but implements a set of
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
#' @rdname biodynsConstructors3
#' @author The FLR Team
#' @seealso \link[base]{[}, \link[base]{[<-}, \link[base]{[[<-},
#' \link[base]{$<-}, \link[methods]{coerce}, \link[base]{lapply},
#' \link[stats]{window}, \link[base]{list}
#' @keywords classes
#' @examples
#' \dontrun{
#' bds <- biodyns("1"=sim(),"2"=biodyn())}
#' 
