utils::globalVariables('end')
utils::globalVariables('fwdCpp')
utils::globalVariables('optim')
utils::globalVariables('par')

# fwd.R - 
# biodyn/R/fwd.R

# Copyontright 2003-2007 FLR Team. Distributed under the GPL 2 or later
# Maintainer: Iago Mosqueira, JRC
# $Id:  $

utils::globalVariables('finite')
globalVariables("control")

setMethod('production', signature(object='biodyn',biomass='FLQuant'),
          function(object,biomass) prdFn("pellat",params(object),biomass))
          
prdFn=function(mdl,params,biomass=0) {
  if (!is.FLQuant(biomass)) biomass=FLQuant(biomass)  
  
  if (dims(params)$iter==1 & dims(biomass)$iter >1) params =propagate(params,dims(biomass)$iter)
  if (dims(params)$iter> 1 & dims(biomass)$iter==1) biomass=propagate(biomass,dims(params)$iter)
  
  mdl=tolower(as.character(mdl))
  
  fox <-function(biomass, params){
    res=exp(1)*(params['msy']%/%params['k'])%*%biomass%*%log(params['k'])%*%(1-log(biomass)%/%log(params['k']))
    res}

  schaefer <- function(biomass, params)
    params['r']*biomass*(1-biomass/params['k'])
  
  pellat <- function(biomass, params){
    a=sweep(biomass,6,params['r']/params['p'],'*')
    b=sweep(biomass,6,params['k'],'/')
    c=sweep(b,      6,params['p'],'^')
    a*(1-c)}
  
  shepherd <- function(biomass,params)
    params['r']*biomass/(1+biomass/params['k'])-params['m']*biomass
  
  gulland <- function(biomass,params)
    params['r']*biomass*(params['k']-biomass)
  
  fletcher <- function(biomass,params) {
    params['p']=params['p']+1
    lambda <- (params['p']^(params['p']/(params['p']-1)))/(params['p']-1)
    lambda*params['msy']*(biomass/params['k'])-lambda*params['msy']*(biomass/params['k'])^params['p']
  }
  
  logistic <- function(biomass, params){
    r=4*params['msy']/params['k']
    r*biomass%*%(1-biomass%/%params['k'])}
  
  genfit <- function(biomass, params)
    params['r']*biomass*(1-biomass/params['k'])
  
  res <- switch(substr(mdl,1,2),
                fo=fox(     biomass,params),
                sc=schaefer(biomass,params),
                gu=gulland( biomass,params),
                fl=fletcher(biomass,params),
                pe=pellat(  biomass,params=params),
                sh=shepherd(biomass,params),
                ge=pellat(  biomass,params),
                lo=logistic(biomass,params))
  
  return(res)}

# calcs a relative value
relFn=function(val,ratio,lag=0){

  if (lag<=0) return(ratio)
  
  yrs=ac(as.numeric(dimnames(ratio)$year)-lag)
  yrs=yrs[yrs %in% dimnames(val)$year]
  yr2=ac(as.numeric(yrs)+lag)
  res=ratio[,yr2]*val[,yrs]
  
  return(res)}

## sets a minimum F
minFn=function(val,bnd){
  val[is.na(val)]=0
  
  res=val
  res[val<bnd]=bnd[val<bnd]
  
  return(res)}

## sets a maximum F
maxFn=function(val,bnd){
  val[is.na(val)]=0
  
  res=val
  res[val>bnd]=bnd[val>bnd]
 
  return(res)}

## calculates inter-annual bounds
iavFn=function(val,bnd,lag=1){
  if (all(!is.finite(bnd)) | lag<1) stop('cant have lags < 1')
  
  ref=apply(val[,dims(val)[2]-lag],2:6,mean)
  
  bnd=val[,dims(val)[2]]*(1-bnd[1])    
  if (finite(bnd[1])) val[,dims(val)[2]][ref<bnd[1]]=ref[ref<bnd[1]]
  
  bnd=val[,dims(val)[2]]*(1+bnd[2])    
  if (finite(bnd[2])) val[,dims(val)[2]][ref>bnd[2]]=ref[ref>bnd[2]]
 
  return(val[,dims(val)[2]])}


#val=FLQuant(rlnorm(200),dimnames=list(age=1:5,year=1:10,iter=1:4))
#bnd=FLQuant(1,          dimnames=list(age=1:5,year=1:10,iter=1:4))


#' @title fwd
#'
#' @description Projects a \code{biodyn} object (i.e. a biomass dynamic model)
#' for a given future catch, harvest rate, or stock biomass. Only
#' one of these has to be supplied as an argument.
#'
#' @param object an object of class \code{biodyn} or  \code{biodyns}
#' @param control missing
# @param catch   an \code{FLQuant} or \code{FLQuants} containing future catches
# @param harvest an \code{FLQuant} or \code{FLQuants} containing future harvest
# @param stock   an \code{FLQuant} or \code{FLQuants} containing future stock    
# @param minF    minimum harvest, \code{FLQuant} or single numeric value    
# @param maxF    maximum harvest, \code{FLQuant} or single numeric value
# @param pe      process error term, an \code{FLQuant}
# @param peMult  logical, i.e. by default equals true so multiplicative, otherwise additive
#' @param ... any other parameters
#'
#' @aliases 
#' fwd-method 
#' fwd,biodyn,FLQuants-method  
#' fwd,biodyn,missing-method 
#' fwd,FLPar,missing-method
#' 
#' @export
#' @rdname fwd
#'
#' @examples
#' \dontrun{
#' bd=simBidyn()
#' harvest=rlnorm(100,log(harvest(bd))[,-dims(bd)$year],.1)
#' 
#' bdHat =fwd(bdHat,harvest=harvest)
#' 
#' plot(bdHat,worm=c(2,8))+
#' theme(legend.position="bottom")
#'  }
setMethod( 'fwd', signature(object='biodyn',fishery='missing',control='missing'),
#setMethod( 'fwd', signature(object='biodyn',control='missing'),
   function(object, fishery, control,
                    catch  =NULL, 
                    harvest=NULL, 
                    f      =NULL, 
                    stock  =NULL, 
                    hcr    =NULL, 
                    pe     =NULL, peMult=TRUE,
                    minF   =0,    maxF  =2,
                    # these next arguments are not currently used, as 
                    # the code havnt been enabled
                    bounds =list(catch=c(Inf,Inf)),
                    lag    =0,
                    end    =NULL,
                    starvationRations=0.75,
                    ...) {
     
     ## target arg is an FLQuant
     if (is.null(harvest)&('FLQuant' %in% class(f)|'FLQuant' %in% class(f))) harvest=f
     
     if ('FLQuant' %in% class(stock)   |
         'FLQuant' %in% class(harvest) |
         'FLQuant' %in% class(catch))
        res=fwdFn(object,control=control,
              catch,harvest,stock,pe,peMult,minF,maxF,bounds,lag,end,
              starvationRations=starvationRations,...)
     else if ('FLQuants' %in% class(stock))  
        res=mpb:::biodyns(llply(stock, function(x) fwdFn(object,control=control,
                 catch,harvest,x,pe,peMult,minF,maxF,bounds,lag,end,...)))
      else if ('FLQuants' %in% class(harvest))
        res=mpb:::biodyns(llply(harvest, function(x) fwdFn(object,control=control,
                 catch,x,stock,pe,peMult,minF,maxF,bounds,lag,end,...))) 
      else if ('FLQuants' %in% class(catch))  
        res=mpb:::biodyns(llply(catch, function(x) fwdFn(object,control=control,
                x,harvest,stock,pe,peMult,minF,maxF,bounds,lag,end,...)))
                        
     res})
     
fwdFn=function(object, 
               catch, 
               harvest, 
               stock, 
               pe, peMult,
               minF,    maxF,
               starvationRations,
               ...){
        
      lag=0
      object@stock=FLQuant(object@stock,quant=names(object@catch)[1])   
    
      ## catch, harvest or stock?
      ctcTrgt=FALSE
      hvtTrgt=FALSE
      stkTrgt=FALSE
      hcrTrgt=FALSE

      if (!is.null(catch))   ctcTrgt=TRUE 
      if (!is.null(harvest)) hvtTrgt=TRUE
      if (!is.null(stock))   stkTrgt=TRUE
      ## hcr option not implemented, see hcr method
      
      if (hvtTrgt) f=harvest
      
      if(!ctcTrgt & hvtTrgt & stkTrgt & hcrTrgt)
          stop('must supply catch, harvest or stock as a target')

      if (ctcTrgt) if (dims(object@stock)$maxyear   < dims(  catch)$maxyear) object=window(object,end=dims(  catch)$maxyear)
      if (hvtTrgt) if (dims(object@stock)$maxyear-1 < dims(harvest)$maxyear) object=window(object,end=dims(harvest)$maxyear+1)
      if (stkTrgt) if (dims(object@stock)$maxyear   < dims(  stock)$maxyear) object=window(object,end=dims(  stock)$maxyear)

      if (stkTrgt) catch=stock*0
     
      ## check year range
      if (ctcTrgt | stkTrgt) {
        if (!(all(dimnames(catch)$year %in% dimnames(object@catch)$year)))
             object = window(object,end=dims(catch)$maxyear)
          if (dims(object@catch)$iter==1 & dims(catch)$iter>1) object@catch=propagate(object@catch, dims(catch)$iter)
          object@catch[,dimnames(catch)$year] <- catch
          yrs <- dimnames(catch)$year
      } else if (hvtTrgt) {
          if (!(all(dimnames(harvest)$year %in% dimnames(object@stock)$year))){
            
            stop('years in harvest & stock dont match')}
          yrs <- dimnames(harvest)$year
      } else if (hcrTrgt) {
         yrs=ac(dims(object@stock)$maxyear:end) 
         object=window(object,end=end)
      }  
     
      if (stkTrgt) yrs = yrs[-length(yrs)]
      ## B0 in year 1?
      if (as.numeric(yrs[1]) == range(object,'minyear')){
         if (!('year' %in% names(dimnames(params(object)['k']))))  
           object@stock[,ac(range(object,'minyear'))] = params(object)['k'] * params(object)['b0'] else
           object@stock[,ac(range(object,'minyear'))] = params(object)['k',ac(range(object,'minyear'))] * params(object)['b0',ac(range(object,'minyear'))]          
         }
      
      ## maxyear
      if (max(as.numeric(yrs)) == range(object,'maxyear'))
         object@stock <- window(object@stock,end=range(object,'maxyear')+1)

      ## niters
      ow=options("warn"); options(warn=-1)
      nits=dims(object)$iter
      options(ow)
      
      if (hvtTrgt) nits=max(nits,dims(harvest)$iter)
      if (ctcTrgt) nits=max(nits,dims(  catch)$iter)
      if (stkTrgt) nits=max(nits,dims(  stock)$iter)
      if (!is.null(pe)) nits=max(nits,dims(pe)$iter)

      if (hvtTrgt) nits=max(nits,dims(harvest)$iter) else
      if (ctcTrgt) nits=max(nits,dims(catch  )$iter) else
      if (stkTrgt) nits=max(nits,dims(stock  )$iter) 
      if (nits>1){ 
         if(dim(object@catch)[6]==1) object@catch =propagate(object@catch,nits)
         if(dim(object@stock)[6]==1) object@stock =propagate(object@stock,nits)
         if (hvtTrgt) if(dim(harvest)[6]==1) harvest=propagate(harvest,nits)
       
         if (dims(params(object))$iter==1) params(object)=propagate(params(object),nits)
         #if (!is.null(pe))                 pe            =propagate(pe            ,nits)
         } 
      
      ## projections
      if (!hvtTrgt) harvest=harvest(object)
      for(y in as.numeric(yrs)) {
    
         ## sp & process error
         if (!is.null(pe)) {    
            if (peMult) sp.=mpb::production(object,object@stock[, ac(y)])%*%pe[, ac(y)] 
            else        sp.=mpb::production(object,object@stock[, ac(y)])%+%pe[, ac(y)]
         } else sp.=mpb::production(object,object@stock[, ac(y)])
         #} else sp.=production(object,object@stock[, ac(y)]*(1-ptYr)+object@stock[, ac(y+1)]*(ptYr))
         
          ## targets, if lag<0 then the targets are relative 
          if (hcrTrgt)
              harvest[,ac(y)]=hcr(window(object,end=y))
          if (hvtTrgt | hcrTrgt){
              object@catch[,ac(y)]=object@stock[,ac(y)]*relFn(harvest(object),harvest[,ac(y)],lag)
          } else if (ctcTrgt){           
              object@catch[,ac(y)]=relFn(object@catch,catch[,ac(y)],lag)
          } else {
              object@catch[,ac(y)]=relFn(object@stock,stock[,ac(y+1)],lag) - object@stock[,ac(y)] - sp.}
         
         object@catch[,ac(y)]=qmin(object@stock[,ac(y)]*starvationRations,
                                   object@catch[,ac(y)])
                                   
         object@stock[,ac(y+1)] = object@stock[,ac(y)] - 
                                  object@catch[,ac(y)] + sp.
         
    ## Not implemented     
    ### bounds
    #      if ('catch' %in% bounds){ 
    #           object@catch[,y]=iavFn(window(object@catch,end=y),bnd,lag=lag)
    #           object@stock[,ac(y+1)] = object@stock[,ac(y)] - object@catch[,ac(y)] + sp.
    #           }
        
          ### min/max
          harvest[,ac(y)]       =maxFn(minFn(harvest(object)[,ac(y)],minF),maxF)
          object@catch[,ac(y)]  =object@stock[,ac(y)]*harvest[,ac(y)]
          object@stock[,ac(y+1)]=object@stock[,ac(y)] - object@catch[,ac(y)] + sp.
          }

        object@stock[stock(object) < 0] = 0.001
        object@catch[catch(object) < 0] = 0.0000001
    
        return(object)}

setMethod( 'fwd', signature(object='biodyn',fishery='FLQuants',control='missing'),
#setMethod('fwd', signature(object='biodyn',control='FLQuants'),
  function(object,fishery,control, pe=NULL, peMult=TRUE,minF=0,maxF=2,lag=0,
           bounds=list(catch=c(Inf,Inf)),...) {
  control=fishery
  res=mlply(seq(length(names(control))),
      function(x,object,control,pe,peMult,minF,maxF,lag,bounds){
        if (names(control)[x]=='catch')  return(fwd(object,catch  =control[[x]],pe=pe,peMult=peMult,minF=minF,maxF=maxF,lag=lag,bounds=bounds))
        if (names(control)[x]=='harvest')return(fwd(object,harvest=control[[x]],pe=pe,peMult=peMult,minF=minF,maxF=maxF,lag=lag,bounds=bounds))
        if (names(control)[x]=='stock')  return(fwd(object,stock  =control[[x]],pe=pe,peMult=peMult,minF=minF,maxF=maxF,lag=lag,bounds=bounds))},
       object=object,control=control,pe=pe,peMult=peMult,minF=minF,maxF=maxF,lag=lag,bounds=bounds)
 
  names(res)=names(control)
        
  return(mpb:::biodyns(res))})

## this is wrapper for a HCR projection, but now being done by the HCR method
# setMethod('fwd', signature(object='biodyn',control='list'),
#   function(object, control, pe=NULL, peMult=TRUE,minF=0,maxF=2,lag=0,
#            bounds=list(catch=c(Inf,Inf)),end=range(object,'maxyear')+15,...) {
#     
#   res=mlply(seq(length(names(control))),
#       function(x,object,control,pe,peMult,minF,maxF,lag,bounds,end){
#         return(fwd(object,hcr=control[[x]],pe=pe,peMult=peMult,minF=minF,maxF=maxF,lag=lag,bounds=bounds,end=end))},
#        object=object,control=control,pe=pe,peMult=peMult,minF=minF,maxF=maxF,lag=lag,bounds=bounds,end=end)
#   
#   names(res)=names(control)
#         
#   return(biodyns(res))})
      
setMethod( 'fwd', signature(object='FLPar',fishery='missing',control='missing'),
#setMethod( 'fwd', signature(object='FLPar',control='missing'),
  function(object,fishery,control,...){
  
   args=list(...)
   if ("catch"%in%names(args)){
     catch=args[["catch"]]
     
     if (any(dim(catch)[c(1,3:5)]!=1)) stop("Only for year and iter >1 in catch")
     if (!all(c("r","k","p","b0")%in%dimnames(object)$params)) stop(cat(c("r","k","p","b0"), "not all in params"))
     if (length(dim(object))!=2) stop("only params and iter allowed in params")

     nits=c(dims(catch)$iter,dims(object)$iter)
     if (length(unique(nits))==2&min(nits)!=1) stop("iter have to be 1 or n")

     par=matrix(object[c("r","k","p","b0"),],nrow=4,            ncol=max(nits))
     ctc=matrix(c(catch),                    nrow=dim(catch)[2],ncol=max(nits))
     
     res=fwdCpp(ctc,par)
     res=FLQuant(unlist(c(res)),dimnames=list(year=dimnames(catch)$year,iter=seq(max(nits))))
     }
          
   res})
           

# names(parSim1)[15:16]=c("ref85","current")
# parSim1=transform(parSim1,b85=ref85/current,
#                   bref=ref85/bmsy)

