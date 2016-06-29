#' hcr
#' 
#' Harvest Control Rule, calculates F, or Total Allowable Catch (TAC) based on a hockey stock harvest control rule.
#'
#' @param object an object of class \code{FLStock} or
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
#' @seealso \code{\link{bmsy}}, \code{\link{fmsy}}, \code{\link{fwd}} and \code{\link{hcrParam}}
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
#' bd=fwd(bd,harvest=hcr(bd,byr=i,yrs=i+1)$hvt)
#' }
setGeneric('hcr', function(object,refs,...) standardGeneric('hcr'))
setMethod('hcr', signature(object='FLStock',refs='FLBRP'),
 function(object,refs, 
           hpar=hcrParam(ftar =0.70*refpts(refs)["msy",'harvest'],
                         btrig=0.80*refpts(refs)["msy",'ssb'],
                         fmin =0.01*refpts(refs)["msy",'harvest'],
                         blim =0.40*refpts(refs)["msy",'ssb']),
           yr =max(as.numeric(dimnames(catch(object))$year))-1,
           byr=yr,
           hyr=yr+2:4,
           sr.residuals=FLQuant(1,dimnames=list(year=hyr)),
           tac   =FALSE,
           tacMn =TRUE,
           bndF  =NULL, #c(1,Inf), #not needed as fmin and maxF sort this
           bndTac=NULL, #c(1,Inf), #absolute
           iaF   =TRUE,            #relative
           iaTac =TRUE,            #relative
           maxF  =2,
           bEr   =NULL,
           tEr   =NULL,
           params=NULL,
           ref   =stock.n(refs)[,1],
           mdd   =FALSE,
           matdd =FALSE,  
           ...) {
  
  ## HCR use hockey stick in yr to set F
  dimnames(hpar)$params=tolower(dimnames(hpar)$params)
  hpar=as(hpar,'FLQuant')  
  #if (blim>=btrig) stop('btrig must be greater than blim')
  a=(hpar['ftar']-hpar['fmin'])/(hpar['btrig']-hpar['blim'])
  b=hpar['ftar']-a*hpar['btrig']

  ## Calc F
  # bug
  #val=(SSB%*%a) %+% b
  bNow=FLCore::apply(ssb(object)[,ac(yr)],6,mean)

  if (!is.null(bEr))
    bNow=bNow*rlnorm(1,0,bEr)
  
  rtn=(bNow%*%a)  
  rtn=FLCore::sweep(rtn,2:6,b,'+')

  # set all iter
  fmin=as(hpar['fmin'],'FLQuant')
  ftar=as(hpar['ftar'],'FLQuant')
  for (i in seq(dims(object)$iter)){
    FLCore::iter(rtn,i)[]=max(FLCore::iter(rtn,i),FLCore::iter(fmin,i))
    FLCore::iter(rtn,i)[]=min(FLCore::iter(rtn,i),FLCore::iter(ftar,i))} 
  
  dimnames(rtn)$year=hyr[1]
  #dimnames(rtn)$year=min(hyr)  
  if (length(hyr)>1){
    rtn=window(rtn,end=max(hyr))
    rtn[,ac(hyr)]=rtn[,ac(min(hyr))]}
  
  ### Bounds ##################################################################################
  ## F
  if (!is.null(bndF)){  

      ref=FLCore::apply(harvest(object)[,ac(byr-1)],6,mean)
    
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
   
  #hvt=rtn
   
   ## TAC
   if (tac){
      #object=FLash:::fwd(object,harvest=harvest(object)[,ac(min(as.numeric(hyr)-1))])
      #object=FLash:::fwd(object,f=fbar(object)[,ac(min(as.numeric(hyr)-1))],sr=eq)
      rtn   =catch(FLash:::fwd(object, f=rtn, sr=refs))[,ac(hyr)]                                  
      
      if (!is.null(bndTac)){  
        ref=FLCore::apply(catch(object)[,ac(byr)],6,mean)
        
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

      if (!is.null(tEr))
        rtn=rtn*rlnorm(1,0,tEr)
      
      if (tac)
        res=FLash:::fwd(object,catch=rtn,sr=refs,sr.residuals=sr.residuals)
      else 
        res=FLash:::fwd(object,    f=rtn,sr=refs,sr.residuals=sr.residuals)
  
      if (mdd|matdd){

        for (i in dimnames(rtn)$year){ 
          scl=(stock.n(res)[,i]%-%scale)%/%scale
          if(mdd){
             m(res)[,i]=mdd(exp(log(stock.wt(res)[,i]%/%params["a"])%/%log(params["b"])),params[c("m1","m2","ddm")],scale=scl)
             }
          if(matdd)
             mat(res)[,i]=matdd(ages(m(res)[,i]),params[c("a50","ato95","asym","ddmat")],scale=scl) 

          if (tac) 
             res=FLash:::fwd(res,catch=rtn[i],sr=refs,sr.residuals=sr.residuals)
           else
             res=FLash:::fwd(res,    f=rtn,sr=refs,sr.residuals=sr.residuals)
          }
      }  
  
     res})

  #matdd(age,params,scale)
  
  #mdd(wt,params,scale) 
  