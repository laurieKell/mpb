#grep -r "setGeneric">>../generic.txt

utils::globalVariables(c("d_ply","dlply","mdply","filter","weighted.mean",
                         "stockCPP","ts","yrs","fCPP","filter","stockCPP",".",
"adply","aes","as.FLQuant","av","cast","catch<-","catch.wt","catch.wt<-","computeCatch",
"computeStock","ddply","desc","dims","discards.wt","discards.wt<-","dlply","dmns",
"element_line","element_rect","element_text","end","f0.1","FLIndex","FLIndices",
"FLQuant","FLQuants","geom_abline","geom_point","ggplot","grid.layout","grid.newpage",
"harvest","harvest<-","harvest.spwn","harvest.spwn<-","hat","laf_open_csv","ldply",
"llply","lm","m","m<-","maply","mat","mat<-","mcmc","mdply","melt","mlply","model.frame",
"modifyList","m_ply","m.spwn","m.spwn<-","n","opts","predict","propagate","pushViewport",
"qqnorm","quantile","read.csv","read.table","readVPA2Box","rec","rel","%+replace%","residuals",
"rstandard","rstudent","scale_x_continuous","scale_y_continuous","sd","series","setPlusGroup",
"ssb","ssb0.1","start","stat_smooth","stock","stock<-","stock.n","stock.n<-","stock.wt","stock.wt<-",
"str_trim","theme","theme_blank","theme_grey","theme_line","theme_rect","theme_segment","theme_text",
"trim","type","unit","V2","viewport","window","write.table","year","yrw"))

#' @title biodyn constructor
#' 
#' @name biodyn
#' 
#' @description Creates an object of the biodyn class representing a biomass dynamic stock assessment model.
#' @param object can be \code{character} object or a \code{FLPar}, \code{FLStock}, \code{aspic} objects
#' @param ...  named parameter being passed to slots
#' 
#' @return biodyn object
#' 
# @aliases biodyn-method biodyn,ANY-method  biodyn,ANY,ANY-method biodyn,FLBRP,FLStock-method
#' 
#' @rdname biodynConstructors
#' @export
#' 
#' @import ggplotFL
#' 
#' @examples 
#' \dontrun{
#' bd=biodyn(params=FLPar(r=0.6,k=50000,p=1,b0=1))
#' }
setGeneric('biodyn',      function(object,params,...)            standardGeneric('biodyn'))

#' @title aspic constructor
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
#  aspic,FLBRP,FLStock-method
#' 
#' @export
#' @rdname aspicConstructors
#' 
#' @examples 
#' \dontrun{
#' asp=aspic()
#' }
setGeneric('aspic',       function(object,value,...)   standardGeneric('aspic'))

setGeneric('biodyns', function(object, ...) standardGeneric('biodyns'))

#' @title aspics
#' @description Create a list of aspic objects
#' @name aspics
#' @param object can be \code{aspic} object or a \code{list} of \code{aspic} objects
#' @param ... additional \code{aspic} objects
#' 
#' @return \code{aspics} object
#' @export
#' @rdname aspicsConstructors
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

#' @title readAspic
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

#' @title writeAspic
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

#' @title msy reference points
#'
#' @description 
#' Calculates maximum sustainable yield (MSY) reference points given the model parameters, for yield, 
#' biomass and harvest rate 
#'
#' @param object an object of class \code{biodyn} or \code{FLPar} with Pella-Tomlinson production
#' function parameters
#' @param ... any other parameters
#' 
#' @return an \code{FLPar} object with benchmark
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
#' refptSD,biodyn,missing-method
#' refptSD,character,FLPar-method
#' refptSD,factor,FLPar-method
#' refpts,FLPar,missing-method
#' refpts,aspic,missing-method
#' refpts,biodyn,missing-method
#' refpts,character,FLPar-method
#' refpts,factor,FLPar-method
#'   
#' @export msy fmsy bmsy refpts refptSD
#' @rdname msy
#'
#' @examples \dontrun{ msy('logistic',FLPar(msy=100,k=500))}
#'   
setGeneric('refptSD',  function(object,params,...) standardGeneric('refptSD'))

#' setParams<-
#' 
#' @title setParams<-
#'
#' @description Sets catchability \code{q} and CV \code{sigma} in the params \code{FLPar} slot 
#' for CPUE provided either an \code{FLQuant} or \code{FLQuants}
#'
#' @param object \code{biodyn}
#' @param value  CPUE as \code{FLQuant} or \code{FLQuants} 
#'
#' @aliases   setParams<-, biodyn,FLPar-method
#'            setParams<-, biodyn,biodyn,FLQuant-method
#'            setParams<-, biodyn,biodyn,FLQuants-method
#'            setParams<-, biodyn,biodyn,data.frame-method
#'                        
#' @export
#' @rdname setParams
#'
#' @examples
#' \dontrun{
#' setParams(bd)=cpue
#' }
#'  
#'  
setGeneric('setParams<-', function(object,value)  standardGeneric('setParams<-'))
#' @title control
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

#' @title params
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
#  setParams<-,biodyn,FLPar-method  setParams<-,biodyn,FLQuant-method setParams<-,biodyn,FLQuants-method  setParams<-,biodyn,data.frame-method setParams<-  setParams<-,biodyn,FLBRP-method 

setGeneric('diags',    function(object)               standardGeneric('diags'))
setGeneric('diags<-',  function(object,value)         standardGeneric('diags<-'))

setGeneric("profile", useAsDefault = profile)

setGeneric('hcr',      function(object,refs,...)      standardGeneric('hcr'))
setGeneric('tac',      function(object, harvest, ...) standardGeneric('tac'))
setGeneric('hcr<-',    function(object,value)         standardGeneric('hcr<-'))

setGeneric('pellat',    function(object,...)         standardGeneric('pellat'))

setGeneric("rate",         function(object,...)   standardGeneric('rate'))
setGeneric('hrate',        function(object)       standardGeneric('hrate'))
setGeneric("instantaneous",function(object,...)   standardGeneric('instantaneous'))
setGeneric('f',            function(object,...) standardGeneric('f'))
setGeneric('computeStock', function(object,...) standardGeneric('computeStock'))


setGeneric('nll',           function(object,index,params,...) standardGeneric('nll'))
setGeneric('grid',          function(object,...)              standardGeneric('grid'))
setGeneric('feasible',      function(object,catch,...)        standardGeneric('feasible'))
setGeneric('setFeasible',   function(object,catch,...)        standardGeneric('setFeasible'))
setGeneric('setFeasible<-', function(object,value)        standardGeneric('setFeasible<-'))

setGeneric('sim',   function(stock,brp,...)     standardGeneric('sim'))
setGeneric('xval',  function(object,index,...) standardGeneric('xval'))

setGeneric('oem',    function(object,...)       standardGeneric('oem'))
setGeneric('oem<-',  function(object,value,...) standardGeneric('oem<-'))
setGeneric('survey', function(object,...)       standardGeneric('survey'))

#' @title jk
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
setGeneric("FLQuantJKs",  function(object, ...)     standardGeneric("FLQuantJKs"))

setGeneric('rand',        function(n,object,sim,...)  standardGeneric('rand'))
setGeneric('resample',    function(x,...)             standardGeneric('resample'))

setGeneric('timeSeries', function(object,params,...) standardGeneric('timeSeries'))
setGeneric('mng',        function(object,params,...) standardGeneric('mng'))
setGeneric('kobe',       function(object,method,...) standardGeneric('kobe'))

setGeneric('fit',        function(object,index,...)  standardGeneric('fit'))
setGeneric('fitmb',      function(object,index,...)  standardGeneric('fitmb'))

setGeneric('production',function(object,biomass,...) standardGeneric('production'))

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
