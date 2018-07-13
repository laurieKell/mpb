utils::globalVariables(c('refJacobian'))
utils::globalVariables(c('jacobian', 'swon'))
utils::globalVariables(c('data.x','data.y','X..x'))
globalVariables("paramas")
utils::globalVariables('qnorm')
utils::globalVariables('fnJ')
utils::globalVariables('data')

#' @title msy
#'
#' @description 
#' Calculates msy based reference points
#'
#' @param \code{object} a \code{biodyn} object  
#' 
#' @return an \code{FLPar} object with an estimate for $msy$
#' 
#' @aliases msy msy-method msy,biodyn-method bmsy bmsy-method bmsy,biodyn-method fmsy fmsy-method fmsy,biodyn-method msy msy-method msy,aspic-method msy,aspic,missing-method bmsy bmsy-method bmsy,aspic-method fmsy fmsy-method fmsy,aspic-method 
#' 
#' @export
#' @rdname msy
#'
#' @examples
#' \dontrun{
#'  msy(bd)
#'  }
setMethod('msy', signature(x='biodyn'), function(x,...){
      params=params(x) 
      res=params['r']*params['k']*(1/(1+params['p']))^(1/params['p']+1)
      dimnames(res)$params="msy"
      res})

# setMethod('fmsy', signature(object='biodyn'), function(object,...)                             
#   fmsyPellaT(params(object)))
# setMethod('bmsy', signature(object='biodyn'), function(object,...)                             
#   bmsyPellaT(params(object)))

setMethod('refpts', signature(object='character'),   
          function(object=factor(object), params=params)          refptsFn(object, params))
setMethod('refpts', signature(object='factor'),   
          function(object=       object)          refptsFn(object, params))
setMethod('refpts', signature(object='biodyn'), 
          function(object)  {  
            model=model( object)
            par  =params(object)
            
            sim=refptsFn(model,par)
            
            if (!("FLParJK"%in%is(par)))
              return(sim)
            else{
              org=refptsFn(model,par@orig)
              return(FLParJK(sim,orig=org))}
            })
setMethod('refpts', signature(object='FLPar'),    
          function(object=object)  refptsFn(factor("pellaT"),params=object))

setMethod('refptSD', signature(object='character', params='FLPar'),   
          function(object=factor(object), params=params)         refptsFn(object, params))
setMethod('refptSD', signature(object='factor',    params='FLPar'),   
          function(object=       object,  params=params)         refptsFn(object, params))
setMethod('refptSD', signature(object='biodyn',    params='missing'), 
          function(object)  {  model=model(object)  
                               par  =params(object)
                               refptsFn(model,par)})

# Fox
msyFox  <- function(params){
  res=params['r']*params['k']*exp(-1)
  dimnames(res)$params="msy"
  res}

# Schaefer
msySchaefer <- function(params){
  res=params['r']*params['k']/4
  dimnames(res)$params="msy"
  res}

msyLogistic<-function(params){
  r=4*params['msy']/params['k']
  res=r*params['k']/2
  dimnames(res)$params="msy"
  res}

# PellaT
msyPellaT<-function(params){
  res=params['r']*params['k']*(1/(1+params['p']))^(1/params['p']+1)
  dimnames(res)$params="msy"
  res}

msyGenfit <- function(params){
  res=params['r']*params['k']*(1/(1+params['p']))^(1/params['p']+1)
  dimnames(res)$params="msy"
  res}

# Shepherd
msyShepherd<-function(params) {
  aPrime<-params['r']/params['m'] - 1
  Bmax  <-params['k']*aPrime
  .bmsy <- 0 #bmsy('shepherd',param)
  
  res=aPrime*params['m']*.bmsy*(1-.bmsy/Bmax)/(1+aPrime)^.5
  dimnames(res)$params="msy"
  res}

# Gulland
msyGulland  <- function(params){
  res=(params['r']*params['k']^2)/4
  dimnames(res)$params="msy"
  res}

msyFletcher <- function(params){
  res=params['msy']
  dimnames(res)$params="msy"
  res}

# Fox
bmsyFox  <- function(params){
  res=params['k']*exp(-1)
  dimnames(res)$params="bmsy"
  res}

# Schaefer
bmsySchaefer <- function(params){
  res=params['k']/2
  dimnames(res)$params="bmsy"
  res}

bmsyLogistic <- function(params){
  res=params['k']/2
  dimnames(res)$params="bmsy"
  res}

# PellaT
bmsyPellaT <- function(params){
  res=params['k']*(1/(1+params['p']))^(1/params['p'])
  dimnames(res)$params="bmsy"
  res}

bmsyGenfit <- function(params){
  res=params['k']*(1/(1+params['p']))^(1/params['p'])
  dimnames(res)$params="bmsy"
  res}

  # Shepherd
bmsyShepherd <- function(params) {
  aPrime <- params['r']/params['m'] - 1
  Bmax  <- params['k']*aPrime
  
  res=Bmax*((1+aPrime)^.5-1)/aPrime
  
  dimnames(res)$params="bmsy"
  res}

# Gulland
bmsyGulland <-function(params){
  res=params['k']/2
  dimnames(res)$params="bmsy"
  res}

# Fletcher
bmsyFletcher <- function(params){
  res=params['k']*(1/(params['p']+1)^(1/(params['p'])))
  dimnames(res)$params="bmsy"
  res}

fmsyPellaT  <-function(params) {
  res=params['r']*(1/(1+params['p']))
  dimnames(res)$params="fmsy"
  res}

fmsyFox     <-function(params) {
  res=msyFox(params)/bmsyFox(params)
  dimnames(res)$params="fmsy"
  res}

fmsySchaefer<-function(params) {
  res=params['r']/2
  dimnames(res)$params="fmsy"
  res}

fmsyShepherd<-function(params) {
  res=msyShepherd(params)/bmsyShepherd(params)
  dimnames(res)$params="fmsy"
  res}

fmsyGulland <-function(params) {
  res=params['r']*params['k']/2
  dimnames(res)$params="fmsy"
  res}

fmsyFletcher<-function(params) {
  res=msyFletcher(params)/bmsyFletcher(params)
  dimnames(res)$params="fmsy"
  res}

fmsyLogistic<-function(params) {r=4*params['msy']/params['k']; r/2}

fmsyFn=function(object,params,probs=0.5){
  
  if (probs!=0.5){
    res=qnorm(probs, object@mng['fmsy','hat'],object@mng['fmsy','sd'])
    return(res)}
  
  object=tolower(object)

  res<-switch(as.character(object),
              
              fox     =fmsyFox(     params),
              schaefer=fmsySchaefer(params),
              gulland =fmsyGulland( params),
              fletcher=fmsyFletcher(params),
              pellat  =fmsyPellaT(  params),
              shepherd=fmsyShepherd(params),
              logistic=fmsyLogistic(params),
              genfit  =fmsyPellaT(params))
  
  dimnames(res)$params='fmsy'
  
  return(res)}

msyFn=function(object,params,probs=0.5) {
   
  if (probs!=0.5){
    res=qnorm(probs, object@mng['msy','hat'],object@mng['msy','sd'])
    return(res)}
  
  
  object=tolower(object)
  res<-switch(object,
              fox      = msyFox(params),
              schaefer = msySchaefer(params),
              gulland  = msyGulland(params),
              fletcher = msyFletcher(params),
              pellat   = msyPellaT(params),
              shepherd = msyShepherd(params),
              logistic = msyLogistic(params),
              genfit   = msyPellaT(params))
  
  dimnames(res)$params='msy'
  
  return(res)}

bmsyFn=function(object,params,probs=0.5) {
  
  if (probs!=0.5){    
      res=qnorm(probs, object@mng['bmsy','hat'],object@mng['bmsy','sd'])
      return(res)}
  
  object=tolower(object)
  res<-switch(object,
              fox     =bmsyFox(params) ,
              schaefer=bmsySchaefer(params),
              gulland =bmsyGulland(params),
              fletcher=bmsyFletcher(params),
              pellat  =bmsyPellaT(params),
              shepherd=bmsyShepherd(params),
              logistic=bmsySchaefer(params),
              genfit  =bmsyPellaT(params))
  
  dimnames(res)$params='bmsy'
  
  return(res)}

#' @title Carrying capacity
#'
#' @description 
#' Calculates $k$ given msy, r and K for a Pella-Tomlinson biomass dynamic model
#'
#' @param msy a guess for MSY
#' @param r a guess for $r$ the population growth rate
#' @param p a guess for $p$ the shape parameter
#' @param params provide $r$ and $p$ as \code{FLPar}
#'
#' @return an \code{FLPar} object with an estimate for $k$
#' 
#' @seealso \code{\link{msy}} and \code{\link{bmsy}} 
#' 
#' @export
#' @rdname K
#'
#' @examples
#' \dontrun{
#'  K(5000,r=.6,p=1.0)
#'  }
K <- function(msy,r=.6,p=1,params=FLPar(r=r,p=p)){
  res=msy/(params['r']*(1/(1+params['p']))^(1/params['p']+1))
  
  dimnames(res)$params='k'
  
  return(res)} 

refptsFn=function(model,params){
  model=tolower(model)
  dmns<-dimnames(params)
  names(dmns)[1]<-'refpts'
  dmns[[1]]<-c('msy','fmsy','bmsy')
  res<-FLPar(as.numeric(NA),dimnames=dmns)
  
  obj=as.character(model)
  
  res[ 'msy']<- msyFn(model, params)
  res['fmsy']<-fmsyFn(model, params)
  res['bmsy']<-bmsyFn(model, params)
  
  return(res)}

funcB=function(x,bd,year=range(bd)['maxyear']){
  if (is.numeric(year)) year=ac(year)
  params(bd)[c('r','k')]=x
  stock(bd[,year])/refpts(bd)['bmsy']}

funcF=function(x,bd,year=range(bd)['maxyear']){
  if (is.numeric(year)) year=ac(year)
  params(bd)[c('r','k')]=x
  harvest(bd[,'2008'])/refpts(bd)['fmsy']}

funcY=function(x,bd,year=range(bd)['maxyear']){
  if (is.numeric(year)) year=ac(year)
  params(bd)[c('r','k')]=x
  catch(bd[,'2008'])/refpts(bd)['msy']}

funcJ=function(x,bd,year=range(bd)['maxyear']){
  if (is.numeric(year)) year=ac(year)
  bd@params[activeParams(bd)]=x
  c(stock  =stock(  bd[,'2008'])/refpts(bd)['bmsy'],
    harvest=harvest(bd[,'2008'])/refpts(bd)['fmsy'],
    catch  =catch(  bd[,'2008'])/refpts(bd)[ 'msy'])}  


refJacobian=function(bd,year=range(bd)['maxyear']){

  x0=bd@params[activeParams(bd)]
  
  res <- jacobian(func=fnJ, x=x0, bd=swon)
  
  nms=dimnames(bd@params[activeParams(bd)])[[1]]
  
  res=FLPar(array(t(res),dim=c(length(nms),3,1),dimnames=list('params'=nms,' '=c('stock','harvest','catch'),iter=1)))
  
  return(res)}

relVar=function(bd,year=range(bd)['maxyear']){
  
  nms=dimnames(bd@params[activeParams(bd)])[[1]]
  res=refJacobian(bd,year=year)
  cov=bd@vcov[nms,nms]
  
  
  res.=as.data.frame(res)
  
  FLCore::iter(cov,1)[lower.tri(iter(cov,1)[drop=T])] =NA
  cov.=as.data.frame(cov)
  cov.=cov.[!is.na(cov.$data),]
  
  vars=merge(cov., merge(res.,res.,by='params',all=T),by='params',all=T)
  
  vars=ddply(transform(vars,var=data.x*data.y*data)[,c('X..x','iter','var')],.(X..x,iter), function(x) data.frame(var=sum(x$var)))
  
  return(res)}  



