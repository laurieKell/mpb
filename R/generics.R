#' biodyn constructor
#' 
#' @name biodyn
#' 
#' @description Creates an object of the biodyn class representing a biomass dynamic stock assessment model.
#' @param object can be \code{character} object or a \code{FLPar}, \code{FLStock}, \code{aspic} objects
#' @param ...  named parameter being passed to slots
#' 
#' @return biodyn object
#' 
#' @aliases biodyn-method biodyn,ANY-method  biodyn,ANY,ANY-method biodyn,FLBRP,FLStock-method
#' 
#' @rdname biodyn-constructors
#' @export
#' 
#' @examples 
#' \dontrun{
#' bd=biodyn(params=FLPar(r=0.6,k=50000,p=1,b0=1))
#' }
setGeneric('biodyn',      function(object,params,...)            standardGeneric('biodyn'))


#' aspic constructor
#' 
#' @name aspic
#' 
#' @description Create an object of the \strong{aspic} class representing a biomass dynamic stock assessment model.
#' @param object can be \code{character} or \code{biodyn} objects
#' @param ...  named parameter being passed to slots
#' 
#' @return aspic object
#' 
#' @aliases aspic-method aspic,ANY-method aspic,character-method aspic,data.frame-method aspic,missing-method aspic,ANY,ANY-method aspic,FLBRP,FLStock-method
#' 
#' @export
#' @rdname aspic-constructors
#' 
#' @examples 
#' \dontrun{
#' asp=aspic()
#' }
setGeneric('aspic',       function(object,value,...)   standardGeneric('aspic'))

#' aspics
#' @description Create a list of aspic objects
#' @name aspics
#' @param object can be \code{aspic} object or a \code{list} of \code{aspic} objects
#' @param ... additional \code{aspic} objects
#' 
#' @return \code{aspics} object
#' @export
#' @rdname aspics-constructors
#' 
#' @aliases aspics-method aspics,missing-method aspics,list-method
#' 
#' @examples 
#' \dontrun{
#' aspics(aspic())
#' }
setGeneric('aspics',      function(object,...)            standardGeneric('aspics'))

#' fit
#'
#' @description 
#' A generic method for fitting catch and index of relative abundance for both \emph{biodyn} and \emph{aspic}.
#' 
#' @param object either \emph{biodyn} or \emph{aspic} class
#' @param index with relative abundance, \emph{FLQuant} or \emph{FLQuants}, if \strong{object} is of type \emph{bodyn}
#' @param ... any other parameter
#' 
#' @rdname fit
#' @export
#' 
#' @aliases fit,aspic-method fit,aspics-method fit,biodyn,FLQuant-method fit,biodyn,FLQuants-method  
#' @seealso  \code{\link{aspic}}, \code{\link{biodyn}}, \code{\link{jk}}, \code{\link{boot}}
setGeneric('fit',       function(object,index,...)  standardGeneric('fit'), package='mpb')

#' jk
#'
#' @description 
#' A generic method for jack knifing a fit to catch and index of relative abundance for both \emph{biodyn} and \emph{aspic}.
#' 
#' @param object either \emph{biodyn} or \emph{aspic} class
#' @param index with relative abundance of type \emph{FLQuantJK} or \emph{FLQuantJKs}, if object is of type \emph{biodyn}
#' @param ... any other parameter
#' 
#' @rdname jk
#' @export
#' 
#' @aliases jk,aspic-method jk,biodyn,FLQuant-method jk,biodyn,FLQuants-method  
#' @seealso  \code{\link{aspic}}, \code{\link{biodyn}}, \code{\link{fit}}, \code{\link{boot}}
setGeneric('jk',        function(object,index,...)        standardGeneric('jk'))

setGeneric('boot',      function(object,...)        standardGeneric('boot'))

#' readAspic
#'
#' @description 
#' Read ASPIC text files, either \code{inp} for inputs, or output files, produced by the executable version of ASPIC
#'
#' @param object an \strong{file path} 
#' @param ... any other parameter
#'
#' @examples
#' \dontrun{
#'     readAspic("aspic.inp")}
#' 
#' @rdname readAspic
#' @export
#' 
#' @aliases readAspic,character-method
#'  
#' @seealso  \code{\link{readFLStock}}, \code{\link{writeAspic}}
setGeneric('readAspic', function(object,...)        standardGeneric('readAspic'))

#' writeAspic
#'
#' @description 
#' Writes the ASPIC text input file \code{inp} to a file or connection.
#' The executable version of ASPIC uses an input file, this method generates that file
#'
#' @param object an \strong{apic} object
#' @param ... any other parameter
#'
#' @examples
#' \dontrun{
#'     writeAspic(albn,"aspic.inp")}
#' 
#' @rdname writeAspic
#' @export
#' 
#' @aliases writeAspic,aspic-method
#'  
#' @seealso  \code{\link{writeFLStock}}, \code{\link{readAspic}}
setGeneric('writeAspic',function(object,...)        standardGeneric('writeAspic'))

#setMethod('writeAspic',    signature(object="aspic"),       
#          function(object,index=object@index,what="FIT",niter=1,fl="aspic.inp",...)        .writeAspicInp(object,index,what,niter,fl=fl,...))

setGeneric('msy',      function(object,params,...)  standardGeneric('msy'))
setGeneric('fmsy',     function(object,params,...)  standardGeneric('fmsy'))
setGeneric('bmsy',     function(object,params,...)  standardGeneric('bmsy'))
setGeneric('refpts',   function(object,params,...)  standardGeneric('refpts'))
setGeneric('refptSE',  function(object,params,...)  standardGeneric('refptSE'))
setGeneric('kobe',     function(object,method,...)  standardGeneric('kobe'))

setGeneric('fwd',      function(object, ctrl, ...)    standardGeneric('fwd'))
setGeneric('hcr',      function(object,refs,...)      standardGeneric('hcr'))
setGeneric('tac',      function(object, harvest, ...) standardGeneric('tac'))
setGeneric('hcr<-',    function(object,value)         standardGeneric('hcr<-'))

#setGeneric("hcr",      function(object, ctrl, ...) standardGeneric("hcr"))
#setGeneric("tac",      function(object, ctrl, ...) standardGeneric("tac"))

setGeneric('diags',    function(object,method,...)    standardGeneric('diags'))
setGeneric('diags<-',  function(object,value)         standardGeneric('diags<-'))

setGeneric('survey',      function(object,...)       standardGeneric('survey'))
setGeneric('index',       function(object,...)       standardGeneric('index'))

#' control
#'
#' @description 
#' Slot for fitting options, i.e. for starting guesses, to fix parameters or set bounds
#' 
#' @param object either \emph{biodyn} or \emph{aspic} class
#' @param value 
#' @param ... any other parameter
#' 
#' @rdname control
#' @export
#' 
#' @aliases control-method control,biodyn,missing-method 
#' 
setGeneric('setControl<-',function(object,...,value)  standardGeneric('setControl<-'))

#' params
#'
#' @description 
#' Slot for estimates
#' 
#' @param object either \emph{biodyn} or \emph{aspic} class
#' @param value 
#' @param ... any other parameter
#' 
#' @rdname control
#' @export
#' 
#' @aliases setParams<-,biodyn,data.frame-method setParams<-,biodyn,FLPar-method setParams<-,biodyn,FLQuant-method setParams<-,biodyn,FLQuants-method setParams<-,aspic,data.frame-method
#' setParams<-,biodyn,FLPar-method  setParams<-,biodyn,FLQuant-method setParams<-,biodyn,FLQuants-method  setParams<-,biodyn,data.frame-method setParams<-  setParams<-,biodyn,FLBRP-method 
#' 
setGeneric('setParams<-', function(object,value,...)     standardGeneric('setParams<-'))
           
setGeneric('computePrd',function(object,biomass,...) standardGeneric('computePrd'))

setGeneric('plot',       function(x,y)               standardGeneric('plot'))
setGeneric('plotIndex',  function(data,...)          standardGeneric('plotIndex'))
setGeneric('plotDiags',  function(data,...)          standardGeneric('plotDiags'))
setGeneric('plotHcr',    function(data, ...)         standardGeneric('plotHcr'))
setGeneric('plotPrd',    function(data,biomass,...)  standardGeneric('plotPrd'))
setGeneric('plotEql',    function(data,biomass,...)  standardGeneric('plotEql'))
setGeneric('plotMSE',    function(x,y,z,...)         standardGeneric('plotMSE'))

setGeneric("profile", useAsDefault = profile)

