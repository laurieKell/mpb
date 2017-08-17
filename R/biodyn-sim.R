## Observation Error Model
globalVariables("ctrl")
globalVariables("prrs")
globalVariables("cpue")
globalVariables("phaseQ")
globalVariables("bounds")
globalVariables("uCV")

#' @title sim
#'
#' @description Creates a biodyn object with known properties
#' 
#' @param model character corresponding to model
#' @param params surplus production parameters
#' @param harvest \code{FLQuant} with harvest rate
#' 
#' @param bounds on \code{control}
#' @param ... other arguments
#' 
#' @export
#' @rdname sim
#' 
#' @return biodyn object with simulated time series
#' 
#' @aliases sim-method sim,FLStock,ANY-method sim,missing,missing-method 
#' 
#' @export
#' 
#' @examples
#' \dontrun{
#'  bd=sim() 
#'  }
setMethod( 'sim',   signature(stock='missing',brp='missing'),
           function(params=FLPar(r=0.5, k=1000, p=1, b0=1.0),
                    harvest=FLQuant(c(seq(0,1.5,length.out=30), 
                                              rev(seq(0.5,1.5,length.out=15))[-1],
                                              rep(0.5,5))*(params['r']*(1/(1+params['p'])))),
                    bounds =c(0.1,10),
                    p=NULL,b0=NULL,...) {

  args <- list(...)

  if (!is.null(p))  params["p"] =p
  if (!is.null(b0)) params["b0"]=b0
  
  nyr=dim(harvest)[2]
  stock =FLQuant(rep(params['k'], nyr), dimnames=dimnames(harvest))

  nyr <- dims(harvest)$year
  object = biodyn('pellat',params=params,
                  stock=stock)
  object@control['r',     'val']=params['r']
  object@control['k',     'val']=params['k']
  object@control['p',     'val']=params['p']
  object@control['b0',    'val']=params['b0']
  
  object@control[,'min']=object@control[,'val']*bounds[1]
  object@control[,'max']=object@control[,'val']*bounds[2]
  
  object@control['p', 'phase']=-1
  object@control['b0','phase']=-1
  object@priors[,1]=-1
  
  # Load given slots
  for(i in names(args))
    slot(object, i) <- args[[i]]

  object <- fwd(object, harvest=harvest)
  
  return(object)}) 

setMethod('sim', signature(stock='FLStock',brp='ANY'),function(stock,brp) {

  bd=biodyn(stock)
  
  params(bd)[dimnames(ctrl)$param]=ctrl[dimnames(ctrl)$param,'val']
  
  bd@priors=prrs
  setParams( bd)=cpue
  setControl(bd)=params(bd)
  bd@control[dimnames(ctrl)$params,'phase'][]=ctrl[dimnames(ctrl)$params,'phase']
  bd@control['q1','phase']=phaseQ
  bd@control['q1','val']  =1
  
  nyr <- dims(harvest)$year
  object = biodyn(model ='pellat',
                  stock =FLQuant(rep(params['k'], nyr), dimnames=dimnames(harvest)),
                  params=params)
  
  object@control['r',     'val']=params['r']
  object@control['k',     'val']=params['k']
  object@control['p',     'val']=params['p']
  object@control['b0',    'val']=params['b0']
  
  object@control[,'min']=object@control[,'val']*bounds[1]
  object@control[,'max']=object@control[,'val']*bounds[2]
  
  object@control['p', 'phase']=-1
  object@control['b0','phase']=-1
  object@priors[,1]=-1
  
  # Load given slots
  for(i in names(args))
    slot(object, i) <- args[[i]]
  
  object <- fwd(object, harvest=harvest)
  
  return(object)})