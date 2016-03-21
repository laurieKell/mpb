udms=function(x) unlist(dims(x))

utils::globalVariables('validParams')

# setParams<-function(model='pellat',its=1)
#   return(FLPar(as.numeric(NA),dimnames=list(params=c(validParams(model),'b0','q','sigma'),iter=its)))

   
   defaultParams<-function(object) {
     params(object)<-FLPar(as.numeric(NA),dimnames=list(params=c(validParams(model(object)),'b0','q','sigma'),iter=1:dims(object)$iter))
     
     unt<-NULL
     if ('r'     %in% dimnames(params(object))$params){
       params(object)['r',    ]<-0.5
       unt<-c(unt,'')}
     if ('k'     %in% dimnames(params(object))$params){
       params(object)['k',    ]<-mean(catch(object))*10
       unt<-c(unt,units(catch(object)))}
     if ('p'     %in% dimnames(params(object))$params){
       params(object)['p',    ]<-2
       unt<-c(unt,'')}
     if ('msy'   %in% dimnames(params(object))$params){
       params(object)['msy',  ]<-mean(catch(object))
       unt<-c(unt,units(catch(object)))}
     if ('b0'    %in% dimnames(params(object))$params){
       params(object)['b0',   ]<-1
       unt<-c(unt,'')}
     if ('m'     %in% dimnames(params(object))$params){
       params(object)['m',    ]<-0.5
       unt<-c(unt,'')}
     if ('q'     %in% dimnames(params(object))$params){
       params(object)['q',    ]<-1.0
       unt<-c(unt,'')}
     if ('sigma' %in% dimnames(params(object))$params){
       params(object)['sigma',]<-0.3
       unt<-c(unt,'')}
     
     units(params(object))<-unt
     
     invisible(params(object))}
   
   # setParams<-function(model='pellat',its=1)
   #   return(FLPar(as.numeric(NA),dimnames=list(params=c(validParams(model),'b0','q','sigma'),iter=its)))
   
   getParams<-function(params,nm){
     if (nm %in% dimnames(params)$params)
       return(c(params[nm,]))
     else
       return(rep(as.numeric(NA),length=dims(params)$iter))}
   

validity<-function(object) {
  return(TRUE)
  ## Catch must be continous
  yrs<-dimnames(catch(object))$year
  
  if (!all(yrs == ac(dims(catch(object))$minyear:dims(catch(object))$maxyear)))
    return('years in catch not continous')
  
  # range
  dims <-dims(object)
  range<-as.list(object@range)
  
  return(TRUE)}

#' control
#' @description sets initial guess and lower and upper bounds
#'
#' @return FLPar
#' @export
#' @examples
#' \dontrun{control(biodyn())}
setGeneric('control',     function(object,...)        standardGeneric('control'))
setMethod( 'control',   signature(object='biodyn'),function(object, ...)   object@control)
   
#' control<-
#'
#' @description sets in \code{biodyn} initial guess and lower and upper bounds
#' 
#' @return \code{biodyn} with new control slot
#' @export
#' 
#' @rdname biodyn-accessors
#' @docType methods
#' 
#' @examples
#'  
#' \dontrun{control(biodyn())}
setGeneric('control<-',   function(object,value)      standardGeneric('control<-'))
setMethod('control<-', signature(object='biodyn', value='FLPar'),
             function(object, value){
               updateFLComp(object, 'control', value)
               return(object)})   
   

setMethod('catch',   signature(object='biodyn'),function(object, ...)   object@catch)
setMethod('catch<-', signature(object='biodyn', value='FLQuant'),
          function(object, value){
            updateFLComp(object, 'catch', value)
            return(object)})

createFLAccesors <- function(class, exclude=character(1), include=missing) {
  
  object <- class

  if(!missing(include))
  	slots <- getSlots(class)[include]
  else
  	slots <- getSlots(class)[!names(getSlots(class))%in%exclude]

	defined <- list()

	for (x in names(slots)) {
		# check method is defined already and signatures match
		eval(
		substitute(if(isGeneric(x) && names(formals(x)) != 'object') {warning(paste('Accesor
			method for', x, 'conflicts with a differently defined generic. Type', x,
			'for more information')); break}, list(x=x))
			)
		# create new generic and accesor method
		eval(
		substitute(if(!isGeneric(x)) setGeneric(x, function(object, ...) standardGeneric(x)),
		list(x=x))
		)
		eval(
		substitute(setMethod(x, signature(y), function(object) return(slot(object, x))),
      list(x=x, y=class))
		)
		# create replacement method
		xr <- paste(x, '<-', sep='')
		eval(
		substitute(if(!isGeneric(x)) setGeneric(x,
			function(object, ..., value) standardGeneric(x)), list(x=xr))
		)
		eval(
		substitute(setMethod(x, signature(object=y, value=v), function(object, value)
			{slot(object, s) <- value; if(validObject(object)) object else stop('')}),
      list(x=xr, y=class, s=x, v=unname(slots[x])))
		)
    if(any(unname(slots[x]) %in% c('FLArray', 'FLQuant', 'FLCohort', 'FLPar')))
    eval(
		substitute(setMethod(x, signature(object=y, value='numeric'), function(object, value)
			{slot(object, s)[] <- value; object}), list(x=xr, y=object, s=x))
		)
		defined[[x]] <- c(x, xr, paste('alias{',x,',', class,'-method}', sep=''),
			paste('\alias{',xr,',', class,',',unname(slots[x]), '-method}', sep=''),
			paste('\alias{',x,'-methods}', sep=''),
			paste('\alias{',xr, '-methods}', sep='')
		)
	}
	return(defined)
}	# }}}

invisible(createFLAccesors('biodyn', 
            exclude=c('desc','range','priors','objFn','mng','diags','stock',
                      'refpts','hessian','hst','ref',
                      'profile','mngVcov','ll'))) #,'priors','diags','objFn','control','mng')))
