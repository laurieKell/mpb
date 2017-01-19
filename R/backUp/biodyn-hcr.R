utils::globalVariables('laply')
utils::globalVariables('ages')

#' @title tac , 
#' 
#' @description 
#' Calculates the Total Allowable Catch for a \code{biodyn} object and target harvest rate
#' by projecting the last year.
#'
#' @param  object an object of class \code{biodyn} or
#' @param  harvest an \code{FLQuant} object with harvest rate
#' @param ... other arguments
#' 
#' @return FLQuant object with TAC value(s)
#' 
#' @export
#' @rdname tac
#' @aliases tac tac-method tac,biodyn-method
#' 
#' @examples
#' \dontrun{
#' tac(bd,FLQuant(0.1,dimnames=list(year=dims(bd)$maxyear)))
#' }
setMethod( 'tac', signature(object='biodyn'),
           function(object,harvest,...){

             yr  =dimnames(harvest)$year  
             #maxY =max(as.numeric(yr))
          
             #stock(object)=window(stock(object),end=maxY)
             #stock(object)[,ac(maxY)]=stock(object)[,ac(maxY-1)]-catch(object)[,ac(maxY-1)]+production(object,stock(object)[,ac(maxY-1)])
             
             #catch(object)=propagate(catch(object),dims(object)$iter)  
             #harvest      =window(harvest,start=dims(object)$year-1)
             #harvest[,ac(dims(object)$year-1)]=harvest(object)[,ac(dims(object)$year-1)]
             
             #object=fwd(object, harvest=harvest(object)[,ac(dimnames(object)$year-1)])
             
             object=window(object, end=max(as.numeric(yr)))
             object=fwd(object,harvest=harvest(object)[,ac(dims(object)$maxyear-1)])
             object=fwd(object, harvest=harvest)
             
             return(catch(object)[,yr])})

#' @title hcrParam
#' 
#' @description 
#' Combines reference points into the HCR breakpts
#'
#' @param ftar an object of class \code{FLPar}
#' @param btrig an object of class \code{FLPar}
#' @param fmin an object of class \code{FLPar}
#' @param blim an object of class \code{FLPar}
#' 
#' @export
#' @rdname hcrParam
#'
#' @examples
#' \dontrun{
#' tac('logistic',FLPar(msy=100,k=500))
#' }
hcrParam=function(ftar,btrig,fmin,blim){
  
  setNms=function(x,nm,nits){
    
    names(dimnames(x))[1]='params'
    dimnames(x)[[1]]     =nm
    if (nits!=dims(x)$iter)
      x=propagate(x,nits)
    
    return(x)}
  
  nits=max(laply(list(ftar,btrig,fmin,blim), function(x) dims(x)$iter))
  
  ftar =setNms(ftar, nm='ftar', nits)
  btrig=setNms(btrig,nm='btrig',nits)
  fmin =setNms(fmin, nm='fmin', nits)
  blim =setNms(blim, nm='blim', nits)
  
  if (nits==1) res=FLPar(  array(c(ftar,btrig,fmin,blim),c(4,nits),dimnames=list(params=c('ftar','btrig','fmin','blim'),iter=seq(nits)))) else
               res=FLPar(t(array(c(ftar,btrig,fmin,blim),c(nits,4),dimnames=list(iter=seq(nits),params=c('ftar','btrig','fmin','blim')))))
  
  #units(res)='harvest'
  return(res)}
  #return(as(res,'FLQuant'))}
  
#' @title hcr
#' 
#' @description 
#' Harvest Control Rule, calculates F, or Total Allowable Catch (TAC) based on a hockey stock harvest control rule.
#'
#' @param object an object of class \code{biodyn} or
#' @param ... other parameters, i.e.
#' params \code{FLPar} object with hockey stick HCR parameters, see hcrParam
#' yr numeric vector with years used to values in HCR
#' byr numeric vector with years used for bounds
#' hyr numeric vector with years to use in projection
#' tac \code{logical} should return value be TAC rather than F?
#' bndF \code{vector} with bounds (i.e.min and max values) on iter-annual variability on  F
#' bndTac \code{vector} with bounds (i.e. min and max values) on iter-annual variability on TAC
#'  
#' @aliases hcr,biodyn-method
#' 
#' @return \code{FLPar} object with value(s) for F or TAC if tac==TRUE
#' 
#' @export
#' @rdname hcr
#'
#' @examples
#' \dontrun{
#' bd   =sim()
#' 
#' bd=window(bd,end=29)
#' for (i in seq(29,49,1))
#' bd=fwd(bd,harvest=hcr(bd,yr=i,yr=i+1)$hvt)
#' }
#setMethod('hcr', signature(object='biodyn',refs='missing'),
hcr=function(object, 
           params=hcrParam(ftar =0.70*refpts(object)['fmsy'],
                           btrig=0.80*refpts(object)['bmsy'],
                           fmin =0.01*refpts(object)['fmsy'],
                           blim =0.40*refpts(object)['bmsy']),
    
           yr =max(as.numeric(dimnames(catch(object))$year)),
           byr=yr-1,
           hyr=yr+1:3,
           
           tac   =FALSE,
           tacMn =TRUE,
          
           bndF  =NULL, #c(1,Inf),
           bndTac=NULL, #c(1,Inf),
           iaF   =TRUE, 
           iaTac =TRUE, 
           maxF  =1,
           ...) {
  
  ## HCR
  dimnames(params)$params=tolower(dimnames(params)$params)
  params=as(params,'FLQuant')  
  #if (blim>=btrig) stop('btrig must be greater than blim')
  a=(params['ftar']-params['fmin'])/(params['btrig']-params['blim'])
  b=params['ftar']-a*params['btrig']

  ## Calc F
  # bug
  #val=(SSB%*%a) %+% b
  stk=FLCore::apply(stock(object)[,ac(yr)],6,mean)
  
  rtn=(stk%*%a)  
  rtn=FLCore::sweep(rtn,2:6,b,'+')

  fmin=as(params['fmin'],'FLQuant')
  ftar=as(params['ftar'],'FLQuant')
  for (i in seq(dims(object)$iter)){
    FLCore::iter(rtn,i)[]=max(FLCore::iter(rtn,i),FLCore::iter(fmin,i))
    FLCore::iter(rtn,i)[]=min(FLCore::iter(rtn,i),FLCore::iter(ftar,i))} 
  
  rtn=window(rtn,end=max(hyr))
  #dimnames(rtn)$year=min(hyr)  
  #if (length(hyr)>1){
  rtn=window(rtn,end=max(hyr))
  rtn[,ac(hyr)]=rtn[,dimnames(rtn)$year[1]]
  #}
  
  ### Bounds ##################################################################################
  ## F
  if (!is.null(bndF)){  
      ref=FLCore::apply(harvest(object)[,ac(byr)],6,mean)

      rtn[,ac(min(hyr))]=qmax(rtn[,ac(min(hyr))],ref*bndF[1])
      rtn[,ac(min(hyr))]=qmin(rtn[,ac(min(hyr))],ref*bndF[2])
    
      if (length(hyr)>1)        
        for (i in hyr[-1]){
          if (iaF){
            rtn[,ac(i)]=qmax(rtn[,ac(i)],rtn[,ac(i-1)]*bndF[1])
            rtn[,ac(i)]=qmin(rtn[,ac(i)],rtn[,ac(i-1)]*bndF[2])
          }else{
            rtn[,ac(i)]=rtn[,ac(i-1)]}
  
      if (!is.null(maxF)) rtn=qmin(rtn,maxF)}}
 
   hvt=rtn
   
   ## TAC
   if (tac){     
      ref=FLCore::apply(catch(object)[,ac(byr)],6,mean)
      
      object=window(object, end=max(as.numeric(hyr)))
      object=fwd(object,harvest=harvest(object)[,ac(min(as.numeric(hyr)-1))])
     
      rtn   =catch(mpb::fwd(object, harvest=rtn))[,ac(hyr)]

      plot(rtn)
      if (!is.null(bndTac)){  
        rtn[,ac(min(hyr))]=qmax(rtn[,ac(min(hyr))],ref*bndTac[1])
        rtn[,ac(min(hyr))]=qmin(rtn[,ac(min(hyr))],ref*bndTac[2])

        if (length(hyr)>1)        
          for (i in hyr[-1]){
            if (iaTac){
              rtn[,ac(i)]=qmax(rtn[,ac(i)],rtn[,ac(i-1)]*bndTac[1])
              rtn[,ac(i)]=qmin(rtn[,ac(i)],rtn[,ac(i-1)]*bndTac[2])
            }else{
              rtn[,ac(i)]=rtn[,ac(i-1)]}}
      }
      if (tacMn) rtn[]=c(apply(rtn,3:6,mean))}
     
  #if (tac) rtn=list(hvt=hvt,tac=rtn,stock=stk) else rtn=list(hvt=hvt,stock=stk)
  if (tac) {
    rtn=window(rtn,start=hyr[1]-1)
    rtn[,ac(hyr[1]-1)]=catch(object)[,ac(hyr[1]-1)]
    return(rtn)
  }else{
    hvt=hvt[,ac(c(hyr[1]-1,hyr))]
    hvt[,ac(hyr[1]-1)]=harvest(object)[,ac(hyr[1]-1)]

    return(hvt)}
  
  return(hvt)}
#)
