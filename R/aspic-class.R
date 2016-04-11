utils::globalVariables(c("year","fwd","setParams<-","setControl<-","K","FLBRP2biodyn","fwd"))


#### Convergence ###############################################################
# 1 100000     					## 0=no MC search, 1=search, 2=repeated srch; N trials #
# 1.00000e-08 					## Convergence crit. for simplex                       #  
# 3.00000e-08 6					## Convergence crit. for restarts, N restarts          #
# 1.00000e-04 0 				## Convergence crit. for estimating effort; N steps/yr #
# 8.00000 							## Maximum F allowed in estimating effort              #
################################################################################

model=factor(c("LOGISTIC", #Schaefer
               "GENGRID",  #generalized model at grid of values or at one specified value
               "FOX",      #Fox
               "GENFIT"))  #Fit the generalized model and estimate its exponent directly.

conditioning=c("YLD", #Condition fitting on yield (recommended for most analyses).
               "EFT") #Condition fitting on fishing-effort rate

objFn=c("SSE",     #Sum of squared errors (recommended default).
        "WTDSSE",  #SSE with annual data weighting
        "LAV")     #Least absolute values (robust objective function).

indexCode=c("CE", "Fishing effort rate, catch (weight)",                  "Effort rate: annual average,Catch: annual total",
            "CC", "Index (weight-based), catch (weight)",                 "Index: annual average,Catch: annual total",
            "B0", "Estimate of biomass Effort rate: annual average",      "Start of year",
            "B1", "Estimate of biomass: annual total",                    "Annual average",
            "B2", "Estimate of biomass: annual average",                  "End of year",
            "I0", "Index of biomass: annual total",                       "Start of year",
            "I1", "Index of biomass Start of year",                       "Annual average", 
            "I2", "Index of biomass Annual average",                      "End of year")

indexCode=t(array(indexCode, dim=c(3,8),dimnames=list(c("code","desc","timing"),NULL)))
dimnames(indexCode)[[1]]=indexCode[,1]
indexCode=transform(indexCode[,-1],startf=c(0,0,0,0,1,0,0,1),
                                 endf  =c(1,1,0,1,1,0,1,1),
                                 ncol  =c(3,3,2,2,2,2,2,2),
                                 col2  =c("effort","index","biomass","biomass","biomass","index","index","index"),
                                 col3  =c("catch", "catch","",       "",       "",       "",     "",     ""))


validAspic <- function(object) {
  ## Catch must be continous
  yrs<-dimnames(catch(object))$year
  
  if (!all(yrs == ac(dims(catch(object))$minyear:dims(catch(object))$maxyear)))
      return("years in catch not continous")
  
  #model         ="factor"
  if (!("factor" %in% is(object@model)) || !(model %in% model))
    stop()
  
  #obj           ="factor",
  #conditioning  ="factor",
  #options       ="numeric",

  # range
  dims <-dims(object)
  range<-as.list(object@range)

  return(TRUE)}


#' ASPIC Biomass Dynamic Model Class
#' 
#' @description A class that represents the ASPIC biomass dynamic stock assessment model.
#' @return aspics object
#' @export
#' @examples
#' 
#' \dontrun{aspic()}
#' 
#' @details 
#' The Class is intended to be used as part of an MSE and includes methods for diagnostics, calculating  reference points and other quantities used when providing management advice.
#' 
#' @slot name         {\code{character} with the name of the stock}      
#' @slot conditioning { \code{factor}}
#' @slot options      { \code{numeric}}     
#' @slot index        { \code{data.frame}}   
#' @slot stopmess     { \code{character}}  
#' @slot rnd          { \code{numeric}}       
#' @slot model        { \code{factor}}        
#' @slot catch        { \code{FLQuant}}      
#' @slot stock        { \code{FLQuant}}       
#' @slot diags        { \code{data.frame}}   
#' @slot params       { \code{FLPar}}    
#' @slot control      { \code{FLPar}}      
#' @slot priors       { \code{array}}        
#' @slot vcov         { \code{FLPar}}         
#' @slot hessian      { \code{FLPar}}       
#' @slot objFn        { \code{FLPar}}        
#' @slot mng          { \code{FLPars}}          
#' @slot name         { \code{character}}     
#' @slot desc         { \code{character}}      
#' @slot range        { \code{numeric}}      
#' 
#' @aliases harvest harvest,aspic-method catch,aspic-method catch<-,aspic,FLQuant-method
#' @rdname aspic-class
#' @docType class
#'  
setClass('aspic', representation(
    "biodyn",
    conditioning  ="factor",
    options       ="numeric",     
    
    index          ="data.frame",
    
    stopmess      ="character",
    rnd           ="numeric"),
  prototype(
    range         =unlist(list(minyear=as.numeric(NA),   maxyear=as.numeric(NA))),
    model         =factor("LOGISTIC",levels=model,       labels=model),
    obj           =factor("SSE",     levels=objFn,       labels=objFn),
    conditioning  =factor("YLD",     levels=conditioning,labels=conditioning),
    options       =c(search=1,trials=100000,simplex=1e-8,restarts=3e-8,nrestarts=6,effort=1e-4,nsteps=0,maxf=8.0),
   
    params        =FLPar(as.numeric(NA),dimnames=list(params=c("b0","msy","k"),iter=1)),
    control       =FLPar(as.numeric(NA),c(length(c(c("b0","msy","k"),paste("q",seq(1),sep=""))),5),
                         dimnames=list(params=c(c("b0","msy","k"),paste("q",seq(1),sep="")),
                                              c("fit","min","val","max","lambda"),iter=1)),
    objFn         =FLPar(array(as.numeric(NA),dim=c(2,1),dimnames=list("value"=c("rss","rsq"),iter=1))),
    vcov          =FLPar(as.numeric(NA),dimnames=list(params=c("b0","msy","k"),param=c("b0","msy","k"),iter=1)),
    hessian       =FLPar(as.numeric(NA),dimnames=list(params=c("b0","msy","k"),param=c("b0","msy","k"),iter=1)),
    stopmess      ="not ran"),
  validity=validAspic)

# printSlot=function(x){
#    res=data.frame(getSlots(x))
#    c("#'@section Slots:", "#'  \\describe{",
#      paste("#' \\item{\\code{", dimnames(res)[[1]], "}: \\code{", as.character(res[,1]),"}}",sep=""),
#      "#'  }")}




