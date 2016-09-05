utils::globalVariables('validParams')
models=factor(c('fox',      'schaefer',
                'pellat',   'gulland',
                'fletcher', 'shepherd',
                'logistic', 'genfit'))

modelParams=function(mdl) {
  if (is.factor(mdl)) mdl=as.character(mdl)
  list(fox       =c('r','k'),
       schaefer  =c('r','k'),
       pellat    =c('r','k','p'),
       shepherd  =c('r','k','m'),
       gulland   =c('r','k'),
       fletcher  =c('k','msy','p'),
       logistic  =c('k','msy'),
       genfit    =c('r','k','p'))[[mdl]]}

validity<-function(object) {
  return(TRUE)
  ## Catch must be continous
  yrs<-dimnames(catch(object))$year
  
  if (!all(yrs == ac(dims(catch(object))$minyear:dims(catch(object))$maxyear)))
    return("years in catch not continous")
  
  # range
  dims <-dims(object)
  range<-as.list(object@range)
  
  return(TRUE)}

#' biodyn class
#'
#' @description A class that implement a biomass dynamic stock assessment model. 
#' 
#' @details 
#' The Class is intended to be used as part of an MSE and includes methods for diagnostics, calculating  reference points and other quantities used when providing management advice.
#' 
#' @slot name    {A \code{character} with the name of the stock}
#' @slot desc    {A \code{character} providing a fuller description of the object}       
#' @slot range   {A \code{numeric} vector containing the quant and year ranges}
#' @slot model   {A \code{factor} giving name of production function, for now this is only `pellat`}
#' @slot obj     { \code{factor} that determines the objective function type -LL or LAV}  
#' @slot catch   {An \code{FLQuant} with total catch by year}        
#' @slot stock   {An \code{FLQuant} which will hold the estimated stock by year}       
#' @slot control {An \code{FLPar} which sets initial guess (val) and bounds (min and max) for each parameter. The phase allows a parameter to be fixed if less <0 and for paramters to be estimated sequentially}       
#' @slot hcr     {A \code{data.frame} with harvest control rule options}       
#' @slot priors  {An \code{array} which sets penalties for parameters}         
#' @slot params  {An \code{FLPar} with parameter estmates}
#' @slot vcov    {An \code{FLPar} with the covariance matrix of the parameter estimates} 
#' @slot hessian {An \code{FLPar} with the hessian of the estimated parameters} 
#' @slot ref     {A \code{numeric} with parameters for estimating mng quantities}
#' @slot mng     {\code{FLPar} with derived quatities of management interest}
#' @slot mngVcov {An \code{FLPar} with the variance matrix of management quanties}
#' @slot diags   {A \code{data.frame} with residuals and covariates from fit of CPUE to stock }     
#' @slot objFn   {\code{FLPar} with objective function} 
#' @slot ll      {\code{FLPar} with negative log likelihood by data component}
#' @slot profile {\code{data.frame} not yet implemented} 
#' @slot hst     {\code{data.frame} not yet implemented} 
#' 
#' All slots in the class have accessor and replacement methods that provide validation and protection of their data.
#' 
#' @export
#' @importFrom plyr ddply ldply laply mdply maply alply mlply llply daply d_ply dlply
#' @importFrom stringr str_trim
#' @importFrom reshape melt
#' @import FLCore FLBRP
#' @import methods
#'
#' @aliases biodyn-class
#' 
#' @rdname biodyn-class
#' @docType class
#'  
#' @examples
#' \dontrun{biodyn()}
 setClass('biodyn', representation(
    'FLComp',
    model         ='factor',   
    obj           ='factor',   
    catch         ='FLQuant',
    stock         ='FLQuant',
    diags         ='data.frame',
    params        ='FLPar',
    control       ='FLPar',
    priors        ='array',
    hcr           ='data.frame',
    vcov          ='FLPar',
    hessian       ='FLPar',
    objFn         ='FLPar',
    ll            ='FLPar',
    ref           ='numeric',
    mng           ='FLPar',
    mngVcov       ='FLPar',
    profile       ='data.frame',
    hst           ='data.frame'
    ),
  prototype(
    range       =unlist(list(minyear=as.numeric(NA), maxyear=as.numeric(NA))),
    ref         =c(nyr=3,nreg=5,yr=NA,then=3),
    catch       =FLQuant(),
    stock       =FLQuant(),
    model       =models[3],
    params      =FLPar(c(.5,as.numeric(NA),1,1),dimnames=list(params=c('r','k','p','b0'),iter=1)),
    control     =FLPar(array(c(1,1,-1,-1,rep(NA,each=6),1,rep(NA,5)), dim=c(4,4,1), dimnames=list(params=c('r','k','p','b0'),option=c('phase','min','val','max'),iter=1))),
    priors      =array(rep(c(0,0,0.3,1),       each=7), dim=c(7,4),   dimnames=list(params=c('r','k','p','b0','msy','bmsy','fmsy'),c('weight','a','b','type'))),
    vcov        =FLPar(array(as.numeric(NA), dim=c(4,4,1), dimnames=list(params=c('r','k','p','b0'),params=c('r','k','p','b0'),iter=1))),
    hessian     =FLPar(array(as.numeric(NA), dim=c(4,4,1), dimnames=list(params=c('r','k','p','b0'),params=c('r','k','p','b0'),iter=1))),
    objFn       =FLPar(array(as.numeric(NA),dim=c(2,1),dimnames=list('value'=c('ll','rss'),iter=1))),
    ll          =FLPar(array(rep(as.numeric(NA),4),dim=c(4,1,1),
                         dimnames=list('params'=c('ll','rss',"sigma","n"),
                                       'index' =1,
                                       'iter'  =1)))), 
	validity=validity) 
