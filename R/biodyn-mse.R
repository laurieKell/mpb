utils::globalVariables(c('eql','srDev'))
utils::globalVariables('priors')
utils::globalVariables('mou')
utils::globalVariables('rcvPeriod')
utils::globalVariables('priors')
utils::globalVariables('cmdOps')
utils::globalVariables('ftar')
utils::globalVariables('btrig')
utils::globalVariables('fmin')
utils::globalVariables('blim')
utils::globalVariables('FLStock2biodyn')

#' @title mseBiodyn
#' 
#' @description Runs a full MSE using an \code{FLStock} object as the Operating Model and \code{biodyn} as the Mangement Procedure
#'           
#' @aliases mse
#' 
#' @param om an \code{FLStock} object 
#  @param eql an \code{FLBRP} object that holds the biological parameters for use in the projections
#' @param mp an \code{biodyn} object that holds the options for the biomass dynamic assessment model
#' @param range a \code{vector} the starting and end years for the projections, and the interval for running the MP
#' @param srDev  a \code{FLQuant} with recruitment deviates
#' @param uDev an \code{FLQuant} or \code{FLQuants} with CPUE residuals
#' @param ftar a \code{numeric} with target F in HCR
#' @param fmin a \code{numeric} with minimum F in HCR
#' @param blim a \code{numeric} with biomass limit for HCR
#' @param btrig a \code{numeric} with biomass trigger (i.e. when to reduce F) in HCR 
#' @param what a \code{character} that specifies what is to be used for the reference point in the HCR, recycled as required
#' @param mult a \code{logical} that specifies whether quantity in HCR options is a multiplier or probability, recycled as required
#'
#' @return  a list of \code{data.frame}s with performance measures from OM and summaries from MP, if \code{con!=NULL} will
#' also write to a MYSQL database
#'  
#' @export
#' @rdname runMSE
#' 
#' @seealso \code{\link{biodyn}}
#' 
#' @examples
#' \dontrun{
#' library(mpb)
#' library(FLash)
#' library(FLBRP)
#' 
#' load(om)
#' load(eql)
#' 
#' om=mpb::fwdWindow(om,eql,end=2030)
#' om=propagate(om,100)
#' 
#' srDev=FLQuant(0,dimnames=list(year=2000:2030))
#' srDev=rlnorm(100,srDev,0.3)
#' 
#' om=mpb::fwd(om,catch=catch(om)[,ac(2000:2011)],sr=eql,sr.residuals=srDev)
#' 
#' library(popbio)
#' 
#' mp=mpb::FLBRP2biodyn(  eql,"biodyn")
#' mp=mpb::FLStock2biodyn(om, "biodyn")
#' }
mseBiodyn<-function(om,eql,
                    srDev,uDev,
                    mp,
                    start=range(om)["maxyear"],end=start+30,interval=3,
                    oem   =oem,
                    hcrPar=function(mp,ftar=0.70,btrig=0.60,fmin=0-01,blim=0.001){
                                hcrParam(ftar =ftar *fmsy(mp),
                                         btrig=btrig*bmsy(mp),
                                         fmin =fmin *fmsy(mp), 
                                         blim =blim *bmsy(mp))},
                    bndF=NULL,bndTac=NULL,maxF=1.0,     
                    omega =1,refB  =1,
                    qTrend=0){
 
  ## Get number of iterations in OM
  nits=c(om=dims(om)$iter, eql=dims(params(eql))$iter, rsdl=dims(srDev)$iter)
  if (length(unique(nits))>=2 & !(1 %in% nits)) ("Stop, iters not '1 or n' in OM")
  if (nits['om']==1) stock(om)=propagate(stock(om),max(nits))

  ## Cut in capacity
  maxF=mean(apply(fbar(window(om,end=start)),6,max)*maxF)
  
  ## open loop feed forward
  mou=om
   
  #### Observation Error (OEM) setup #######################
  ## Random variation for CPUE
  cpue=oem(window(om,end=start),cv=uDev,fishDepend=TRUE)
  if (!("FLQuant"%in%is(qTrend)))
    qTrend=FLQuant(cumprod(rep(1+qTrend,(end+interval-as.numeric(dims(cpue)["minyear"])+1))),
                   dimnames=list(year=range(om)["minyear"]:(end+interval)))
  cpue=cpue*qTrend[,dimnames(cpue)$year]
  
  ## Loop round years
  mp =NULL
  hcr=NULL
  for (iYr in seq(start,range(om,'maxyear')-interval,interval)){
    #iYr = seq(start+rcvPeriod,range(om,'maxyear')-interval,interval)[1]
    cat('\n===================', iYr, '===================\n')
    
    ## use data from last year
    cpue=window(cpue,end=iYr-1)
    cpue[,ac(iYr-(interval:1))]=oem(om[,ac(iYr-(interval:1))],uDev,fishDepend=TRUE)
    cpue[,ac(iYr-(interval:1))]=cpue[,ac(iYr-(interval:1))]*qTrend[,ac(iYr-(interval:1))]
    
    #### Management Procedure
    ## Set up assessment parameter options
    bd=FLStock2biodyn(window(om,end=iYr-1))
    pnms=dimnames(control)$param[dimnames(control)$param%in%dimnames(params(bd))$params]
    params(bd)[pnms]=control[pnms,'val']
    
    bd@priors=priors
    setParams( bd)=cpue 
    setControl(bd)=params(bd)
    bd@control[dimnames(control)$params,'phase'][]=control[dimnames(control)$params,'phase']
    bd@control['q1','phase']=phaseQ
    bd@control['q1','val']  =1
    
    ## fit
    bd =fit(bd,cpue,cmdOps=cmdOps)
    bd =mpb::fwd(bd,catch=catch(om)[,ac(iYr)])
        
    ## HCR
    hcrPar=hcrParam(ftar =ftar *fmsy(bd),
                     btrig=btrig*bmsy(bd),
                     fmin =fmin *fmsy(bd), 
                     blim =blim *bmsy(bd))
    hcrOutcome=hcr(bd,hcrPar,
                   hcrYrs=iYr+seq(interval),
                   bndF  =bndF,
                   bndTac=bndTac,
                   tac =TRUE)
            
    ## TACs for next year (iYtr+1) for n=interval years
    TAC  =hcrOutcome$tac
    TAC[]=rep(apply(TAC,6,mean)[drop=T],each=interval)
    
    #### Operating Model Projectionfor TAC
    om =mpb::fwd(om,catch=TAC,maxF=maxF,sr=eql,sr.residuals=srDev)  

    #### Summary Statistics
    ## HCR actions, i.e. is biomass<Btrig?, what is F?, ..
    hcr =rbind(hcr,data.frame(yearHcr=min(as.numeric(dimnames(hcrOutcome$hvt)$year)),
                              #yearAss=rep(range(bd)[2],dims(bd)$iter),
                              model.frame(           hcrPar,drop=T)[,-5],
                              tac    =as.data.frame(apply(hcrOutcome$tac,6,mean),drop=T)[,'data'],
                              harvest=as.data.frame(apply(hcrOutcome$hvt,6,mean),drop=T)[,'data'],
                              stock  =as.data.frame(hcrOutcome$stock,drop=T)[,2]))
    
    ## Assessment parameters and reference points
    mp =rbind(mp,cbind(cbind(year=iYr,model.frame(params(bd))),
                       model.frame(refpts(bd))[,-4],
                       hcr))
    }
  
  ## save OM, projection without feedback, last assessment and MP summary
  return(list(om=om,mou=mou,bd=bd,mp=mp,oem=mcf(FLQuants(cpue=cpue,catch=catch(om)))))}


hcrFn2=function(om,btrig,blim,ftar,fmin,start,end,interval,lag=seq(interval)){    
  
  a=(ftar-fmin)/(btrig-blim)
  b=ftar-a*btrig
  
  for (iYr in seq(start+rcvPeriod,range(om,'maxyear')-interval,interval)){
    stk=FLCore::apply(stock(om)[,ac(iYr-1)],6,mean)
    
    trgt=(stk%*%a)+b  
    
    for (i in seq(dims(om)$iter)){
      FLCore::iter(trgt,i)[]=max(FLCore::iter(trgt,i),FLCore::iter(fmin,i))
      FLCore::iter(trgt,i)[]=min(FLCore::iter(trgt,i),FLCore::iter(ftar,i))} 
    
    dmns     =dimnames(trgt)
    dmns$year=as.character(iYr+lag)
    
    trgt=FLQuant(rep(c(trgt),each=length(lag)),dimnames=dmns)
    
    om=mpb::fwd(om,f=trgt,sr=eql,sr.residuals=srDev)
  }
  
  return(om)}

  
#with FLStock OM
demoBiodyn<-function(om,mp,
                     eDev  =0.3,
                     uDev  =0.2,
                     u     =oem,
                     ftar  =0.70,    btrig =0.60,
                     fmin  =0.01,    blim  =0.01,
                     bndF  =NULL,
                     start   =range(mp)["maxyear"],          
                     end     =start+30,         
                     interval=3,
                     maxF    =2.0,     
                     cmdOps=paste('-maxfn 500 -iprint 0 -est')){
  
  ## Get number of iterations in OM
  nits=c(om=dims(om)$iter, eql=dims(params(eql))$iter, rsdl=dims(srDev)$iter)
  if (length(unique(nits))>=2 & !(1 %in% nits)) ("Stop, iters not '1 or n' in OM")
  if (nits['om']==1) stock(om)=propagate(stock(om),max(nits))
  
  ## Cut in capacity
  maxF=mean(apply(fbar(window(om,end=start)),6,max)*maxF)
  
  
  #### Observation Error (OEM) setup #######################
  ## Random variation for CPUE  
  cpue=oem(window(om,end=start),uDev,fishDepend=TRUE)
  
  ## Loop round years
  mp =NULL
  hcr=NULL
  for (iYr in seq(start,range(om,'maxyear')-interval,interval)){
    #iYr = seq(start+rcvPeriod,range(om,'maxyear')-interval,interval)[1]
    cat('\n===================', iYr, '===================\n')
    
    ## use data from last year
    cpue=window(cpue,end=iYr-1)
    cpue[,ac(iYr-(interval:1))]=oem(om[,ac(iYr-(interval:1))],uDev,fishDepend=TRUE)
    
    #### Management Procedure
    ## Set up assessment parameter options
    bd=FLStock2biodyn(window(om,end=iYr-1))
    pnms=dimnames(control)$param[dimnames(control)$param%in%dimnames(params(bd))$params]
    params(bd)[pnms]=control[pnms,'val']
    
    bd@priors=priors
    setParams( bd)=cpue 
    setControl(bd)=params(bd)
    bd@control[dimnames(control)$params,'phase'][]=control[dimnames(control)$params,'phase']
    bd@control['q1','phase']=phaseQ
    bd@control['q1','val']  =1
    
    ## fit
    bd =fit(bd,cpue,cmdOps=cmdOps)
    bd =mpb::fwd(bd,catch=catch(om)[,ac(iYr)])
    
    ## HCR
    hcrPar=hcrParam(ftar =ftar *fmsy(bd),
                    btrig=btrig*bmsy(bd),
                    fmin =fmin *fmsy(bd), 
                    blim =blim *bmsy(bd))
    hcrOutcome=hcr(bd,hcrPar,
                           hcrYrs=iYr+seq(interval),
                           bndF=bndF,
                           tac =TRUE)
    
    ## TACs for next year (iYtr+1) for n=interval years
    TAC  =hcrOutcome$tac
    TAC[]=rep(apply(TAC,6,mean)[drop=T],each=interval)
    
    #### Operating Model Projectionfor TAC
    om =mpb::fwd(om,catch=TAC,maxF=maxF,sr=eql,sr.residuals=srDev)  
    
    #### Summary Statistics
    ## HCR actions, i.e. is biomass<Btrig?, what is F?, ..
    hcr =rbind(hcr,data.frame(yearHcr=min(as.numeric(dimnames(hcrOutcome$hvt)$year)),
                              #yearAss=rep(range(bd)[2],dims(bd)$iter),
                              model.frame(           hcrPar,drop=T)[,-5],
                              tac    =as.data.frame(apply(hcrOutcome$tac,6,mean),drop=T)[,'data'],
                              harvest=as.data.frame(apply(hcrOutcome$hvt,6,mean),drop=T)[,'data'],
                              stock  =as.data.frame(hcrOutcome$stock,drop=T)[,2]))
    
    ## Assessment parameters and reference points
    mp =rbind(mp,cbind(cbind(year=iYr,model.frame(params(bd))),
                       model.frame(refpts(bd))[,-4],
                       hcr))
  }
  
  ## save OM, projection without feedback, last assessment and MP summary
  return(list(om=om,mou=mou,bd=bd,mp=mp,oem=mcf(FLQuants(cpue=cpue,catch=catch(om)))))}



#with biodyn OM
demo<-function(om,mp,pe,
               uDev  =0.3,
               oem   =function(x,uDev){
                            dmns=list(year=dimnames(stock(x))$year[-dim(stock(x))[2]])
                            res =rlnorm(dim(stock(x))[6],
                                   FLQuant(0,dimnames=dmns),uDev)
  
                            res*(stock(x)[,-1]+stock(x)[,-dim(stock(x))[2]])/2},
               ftar    =0.70,    
               btrig   =0.60,
               fmin    =0.01,   
               blim    =0.01,
               bndF    =NULL,
               start   =range(mp)["maxyear"],          
               end     =start+30,
               interval=3,
               maxF    =.75, 
               qKnown  =FALSE,
               cmdOps=paste('-maxfn 500 -iprint 0'),
               silent=FALSE){

  om=window(om,end=start+1)

  ## Get number of iterations in OM
  nits=dims(om)$iter
  if (dims(mp)$iter==1&nits>1) mp=propagate(mp,nits)
  
  ## Cut in capacity
  #maxF=apply(harvest(window(om,end=start)),6,max,na.rm=T)*maxF
  
  #### Observation Error (OEM) setup #######################
  ## Random variation for CPUE  
  cpue=oem(window(om,end=start),uDev)
  setParams( mp)=cpue  
  setControl(mp)=params(mp)  
  
  if (qKnown){
    control(mp)["q1",c("phase","val")]=c(-1,1)
    print(control(mp))
    }   
  
  ## Loop round years
  for (iYr in seq(start,end-interval,interval)){
    #iYr = seq(start+rcvPeriod,range(om,'maxyear')-interval,interval)[1]
    if (!silent)
    cat('\n===================', iYr, '===================\n')
   
    ## use data from last year
    cpue=window(cpue,end=max(iYr-(interval:0)))
    cpue[,ac(iYr-(interval:1))] =oem(om[,ac(iYr-(interval:0))],uDev)
    mp=window(mp,end=iYr-1)
    catch(mp)=catch(om)
    #### Management Procedure
    
    ## fit
    mp =fit(window(mp,end=iYr-1),cpue,
            cmdOps=ifelse(nits>1,cmdOps,'-maxfn 500 -iprint 0'))
    mp =mpb::fwd(mp,catch=catch(om)[,ac(iYr)],maxF=maxF)

    ## HCR
    hcrPar=mpb::hcrParam(ftar =ftar *fmsy(mp),
                    btrig=btrig*bmsy(mp),
                    fmin =fmin *fmsy(mp), 
                    blim =blim *bmsy(mp))
    TAC  =mpb::hcr(mp,params=hcrPar,
                       hyr =iYr+seq(interval),
                       bndF=bndF,
                       tac =TRUE,
                       maxF=maxF)
    
    #TAC[]=rep(apply(TAC,6,mean)[drop=T],each=interval)
    TAC[is.na(TAC)]=10
    TAC=qmax(TAC,10)
    TAC=qmin(TAC,bmsy(mp)*fmsy(mp))
    
    #### Operating Model Projectionfor TAC
    om =mpb::fwd(om,catch=TAC,pe=pe,maxF=maxF)
    
    #print((1:100)[is.na(stock(om)[,ac(iYr)])])
    }
  
  if (!silent)
    print("bug biodyns needs a list")  
  
  biodyns(list("om"=window(om,end=end-interval),"mp"=window(mp,end=end-interval)))}

