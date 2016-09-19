#grep -r "setGeneric">>../generic.txt

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
#' @aliases aspic-method 
#' aspic,ANY-method 
#' aspic,character-method 
#' aspic,data.frame-method 
#' aspic,missing-method 
#' aspic,ANY,ANY-method 
#' aspic,FLBRP,FLStock-method
#' 
#' @export
#' @rdname aspic-constructors
#' 
#' @examples 
#' \dontrun{
#' asp=aspic()
#' }
setGeneric('aspic',       function(object,value,...)   standardGeneric('aspic'))

#' biodyns
#' @description Create a list of biodyn objects
#' @name biodyn
#' @param object can be \code{biodyn} object or a \code{biodyn} of \code{biodyn} objects
#' @param ... additional \code{biodyn} objects
#' 
#' @return \code{biodyns} object
#' @export
#' @rdname biodyns-constructors
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
setGeneric('biodyns', function(object, ...) standardGeneric('biodyns'))

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
#' @aliases aspics-method 
#' aspics,missing-method 
#' aspics,aspic-method
#' aspics,list-method
#' 
#' @examples 
#' \dontrun{
#' aspics(aspic())
#' }
setGeneric('aspics',      function(object,...)            standardGeneric('aspics'))

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

#' msy reference points
#'
#' Benchmarks for yield, biomass and harvest rate
#' 
#' Calculates maximum sustainable yield (MSY) reference points given the model parameters, 
#'
#' @param object an object of class \code{biodyn} or \code{FLPar} with Pella-Tomlinson production
#' function parameters
#' @param ... any other parameters
#' 
#' @return an \code{FLPar} object with benchmark
#' 
#' @seealso \code{\link{refpts}} 
#'
#' @aliases 
#' msy,biodyn-method  
#' fmsy,biodyn-method  
#' bmsy,biodyn-method  
#' msy,aspic-method  
#' fmsy,aspic-method  
#' bmsy,aspic-method  
#' msy,FLPar-method  
#' fmsy,FLPar-method  
#' bmsy,FLPar-method  
#'   
#' @export
#' @rdname msy
#'
#' @examples \dontrun{ msy('logistic',FLPar(msy=100,k=500))}
#'   
setGeneric('msy',      function(object,...) standardGeneric('msy'))
setGeneric('fmsy',     function(object,...) standardGeneric('fmsy'))
setGeneric('bmsy',     function(object,...) standardGeneric('bmsy'))
setGeneric('refpts',   function(object,params,...) standardGeneric('refpts'))
setGeneric('refptSE',  function(object,params,...) standardGeneric('refptSE'))

#' setParams<-
#'
#' @description Sets catchability \code{q} and CV \code{sigma} in the params \code{FLPar} slot 
#' for CPUE provided either an \code{FLQuant} or \code{FLQuants}
#'
#' @param object \code{biodyn}
#' @param value  CPUE as \code{FLQuant} or \code{FLQuants} 
#' @param ... any other parameter
#'
#' 
#' @export
#' @rdname setParams
#'
#' @aliases 
#' setControl<-,biodyn,FLQuant-method  
#' setControl<-,biodyn,FLQuants-method
#' 
#' @examples
#' \dontrun{
#' setParams(bd)=cpue
#' }
#'  
setGeneric('setParams<-', function(object,value,...)  standardGeneric('setParams<-'))

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
setGeneric('control',     function(object,...)        standardGeneric('control'))
setGeneric('control<-',   function(object,value)      standardGeneric('control<-'))
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
#' @rdname params
#' @export
#' 
#' @aliases setParams<-,biodyn,data.frame-method setParams<-,biodyn,FLPar-method setParams<-,biodyn,FLQuant-method setParams<-,biodyn,FLQuants-method setParams<-,aspic,data.frame-method
#' setParams<-,biodyn,FLPar-method  setParams<-,biodyn,FLQuant-method setParams<-,biodyn,FLQuants-method  setParams<-,biodyn,data.frame-method setParams<-  setParams<-,biodyn,FLBRP-method 
#' 
setGeneric('index',       function(object,...)    standardGeneric('index'))
setGeneric("setIndex<-",  function(object,value)  standardGeneric('setIndex<-'))


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

setGeneric('diags',    function(object,method,...)    standardGeneric('diags'))
setGeneric('diags<-',  function(object,value)         standardGeneric('diags<-'))

setGeneric("profile", useAsDefault = profile)

setGeneric('fwd', function(object,ctrl,...) standardGeneric('fwd'))
setGeneric('fwd',      function(object, ctrl, ...)    standardGeneric('fwd'))
setGeneric('hcr',      function(object,refs,...)      standardGeneric('hcr'))
setGeneric('tac',      function(object, harvest, ...) standardGeneric('tac'))
setGeneric('hcr<-',    function(object,value)         standardGeneric('hcr<-'))

#' production
#'
#' @description 
#' Estimates production for a given biomass
#' 
#' @param object either \emph{biodyn} or \emph{FLBRP} class
#' @param biomass 
#' @param ... any other parameter
#' 
#' @rdname production
#' @export
#' 
#' @aliases production,biodyn-method 
#'  
setGeneric('production',function(object,biomass,...) standardGeneric('production'))

setGeneric('pellat',    function(object,...)         standardGeneric('pellat'))

setGeneric('plot',           function(x,y)               standardGeneric('plot'))
setGeneric('plotIndex',      function(data,...)          standardGeneric('plotIndex'))
setGeneric('plotDiags',      function(data,...)          standardGeneric('plotDiags'))
setGeneric('plotHcr',        function(data, ...)         standardGeneric('plotHcr'))
setGeneric('plotProduction', function(data,biomass,...)  standardGeneric('plotProduction'))
setGeneric('plotEql',        function(data,biomass,...)  standardGeneric('plotEql'))
setGeneric('plotMSE',        function(x,y,z,...)         standardGeneric('plotMSE'))
setGeneric('plotHcr',        function(object,...)        standardGeneric('plotHcr'))
setGeneric('plotCc',         function(data,...)          standardGeneric('plotCc'))

setGeneric("rate",         function(object,...)   standardGeneric('rate'))
setGeneric('hrate',        function(object)       standardGeneric('hrate'))
setGeneric("instantaneous",function(object,...)   standardGeneric('instantaneous'))
setGeneric('f',            function(object,...) standardGeneric('f'))
setGeneric('computeStock', function(object,...) standardGeneric('computeStock'))


setGeneric('nll',           function(object,index,params,...) standardGeneric('nll'))
setGeneric('grid',          function(object,...)              standardGeneric('grid'))
setGeneric('feasible',      function(object,catch,...)        standardGeneric('feasible'))
setGeneric('feasible',      function(object,catch,...)        standardGeneric('feasible'))
setGeneric('setFeasible',   function(object,catch,...)        standardGeneric('setFeasible'))
setGeneric('setFeasible<-', function(object,value,...)        standardGeneric('setFeasible<-'))

setGeneric('sim',   function(stock,brp,...)     standardGeneric('sim'))
setGeneric('xval',  function(object,index,...) standardGeneric('xval'))

setGeneric('oem',    function(object,...) standardGeneric('oem'))
setGeneric('survey', function(object,...) standardGeneric('survey'))

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
setGeneric('boot',        function(object,...)        standardGeneric('boot'))
setGeneric('randJack',    function(n,object,sim,...)  standardGeneric('randJack'))
setGeneric("jackSummary", function(object, ...)       standardGeneric("jackSummary"))
setGeneric("FLQuantJKs",  function(object, ...)     standardGeneric("FLQuantJKs"))

setGeneric('rand',        function(n,object,sim,...)  standardGeneric('rand'))
setGeneric('resample',    function(x,...)             standardGeneric('resample'))

setGeneric('timeSeries', function(object,params,...) standardGeneric('timeSeries'))
setGeneric('mng',        function(object,params,...) standardGeneric('mng'))
setGeneric('kobe',       function(object,method,...) standardGeneric('kobe'))

setGeneric('fitmb',      function(object,index,...)  standardGeneric('fitmb'))

# #' Sum of vector elements.
# #' 
# #' \code{sum} returns the sum of all the values present in its arguments.
# #' 
# #' This is a generic function: methods can be defined for it directly
# #' or via the \code{\link{Summary}} group generic. For this to work properly,
# #' the arguments \code{...} should be unnamed, and dispatch is on the
# #' first argument.
# #' 
# #' @param
# #' 
# #' @return
# #' 
# #' @export
# #' 
# #' @rdname
# #' 
# #' @examples
# #' \dontrun{
# #' }
