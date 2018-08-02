hcrConstantCatch<-function(yrs,catch,...){
  res=FLQuant(c(apply(catch,6,mean)), 
              dimnames=list(year=yrs,iter=dimnames(catch)$iter)) 
  res}

hcrSBT1<-function(yrs,control=c(k1=2.0,k2=3.0,gamma=1),catch,cpue,...){
  
  lambda=as.FLQuant(ddply(as.data.frame(cpue%/%apply(cpue,6,mean)), 
                          .(iter), with,  data.frame(data=coefficients(lm(data~year))[2])))
  
  flag  =lambda<0
  lambda=abs(lambda)
  res= 1+ifelse(flag,-control["k1"],control["k2"])*exp(log(lambda)*ifelse(flag,control["gamma"],1))
  res=res%*%apply(catch,6,mean)
  
  dmns=dimnames(catch)
  dmns$year=yrs
  dmns$iter=dimnames(cpue)$iter
  
  res=FLQuant(rep(c(res),each=length(yrs)),dimnames=dmns)
  
  return(res)}

hcrSBT2=function(yrs,index,ref,
                 ei=0.25,ed=0.25,
                 tac,target){
  
  flag    =index<ref
  cBit    =target*(index/ref)*(1+ifelse(flag,-ed,ei))
  res     =tac+bit
  
  #   cat('TAC:',        as.integer(mean(tac)),
  #       '\t ratio:',   as.integer((mean(index/ref))),
  #       '\t delta:',   as.integer((mean(bit))),
  #       '\t New TAC:', as.integer(mean(res)),
  #       '\t bit:',    mean(bit),'\n')
  
  dmns=dimnames(tac)
  dmns$year=yrs
  dmns$iter=dimnames(index)$iter
  
  res=FLQuant(rep(c(res),each=interval),dimnames=dmns)
  
  return(res)}


hcrSBT2Orig=function(adult,juve,
                 yrAdult,yrJuve,
                 refJuve=-(1:5),
                 tac,tarCatch,eb=0.25,er=0.75,lag=1,interval=3){
  
  adultIdx=adult[,ac(dims(adult)$maxyear)]
  adultRef=aaply(adult[,ac(yrAdult)],3:6,mean)
  flag    =adultIdx<adultRef
  cBit    =tarCatch*(adultIdx/adultRef)*(1+ifelse(flag,-eb,eb))
  
  juveIdx =aaply(juve[,ac(dims(juve)$maxyear+refJuve)],3:6,mean)
  juveRef =aaply(juve[,ac(yrJuve) ],3:6,mean)
  flag    =juveIdx<juveRef
  rBit    =(juveIdx/juveRef)*(1+ifelse(flag,er,-er))
  
  # cat('ref Juve:',   as.integer(mean(refJuve)),
  #     '\t Juve:',    as.integer(mean(juve)),
  #     '\t ratio:',   mean(juve/refJuve),
  #     '\t rBit:',    mean(rBit),'\n')
  
  res =0.5*(tac+cBit*rBit)
  
  #   cat('TAC:',        as.integer(mean(tac)),
  #       '\t ratio:',   as.integer((mean(adult/refAdult))),
  #       '\t delta:',   as.integer((mean(cBit))),
  #       '\t New TAC:', as.integer(mean(res)),
  #       '\t rBit:',    mean(rBit),'\n')
  
  dmns=dimnames(tac)
  dmns$year=as.character(as.integer(dmns$year)+lag+seq(interval)-1)
  dmns$iter=dimnames(adult)$iter
  
  res=FLQuant(rep(c(res),each=interval),dimnames=dmns)
  
  return(res)}

mseEMP<-function(
  #OM as FLStock and FLBRP
  om,eq,
  
  #MP
  mp,
  control="missing",
  
  srDev,
  uDev,
  
  #years over which to run MSE, doesnt work if interval==1, this is a bug
  interval=1,start=range(om)["maxyear"]-30,end=range(om)["maxyear"]-interval,
  
  #Capacity, i.e. F in OM can not be greater than this
  maxF=1.5){
  
  ##So you dont run to the end then crash
  end=min(end,range(om)["maxyear"]-interval)
  
  ## Make sure number of iterations are consistent
  nits=c(om=dims(om)$iter, eq=dims(params(eq))$iter, rsdl=dims(srDev)$iter)
  if (length(unique(nits))>=2 & !(1 %in% nits)) ("Stop, iters not '1 or n' in om")
  if (nits['om']==1) stock(om)=propagate(stock(om),max(nits))
  
  ## Limit on capacity, add to fwd(om) if you want
  maxF=median(FLQuant(1,dimnames=dimnames(srDev))%*%apply(fbar(window(om,end=start)),6,max)*maxF)
  
  ## Observation Error (OEM) setup
  cpue=window(stock(om),end=start)
  cpue=cpue%*%uDev[,dimnames(cpue)$year]
  
  ## Loop round years
  cat('\n==')
  for (iYr in seq(start,end,interval)){
    cat(iYr,", ",sep="")
    
    ## Observation Error, using data from last year back to the last assessment
    ## CPUE
    cpue=window(cpue,end=iYr-1)
    cpue[,ac(iYr-(interval:1))]=stock(om)[,ac(iYr-(interval:1))]%*%uDev[,ac(iYr-(interval:1))]
    
    #### Management Procedure
    ##Constant catch
    #tac=hcrConstantCatch(iYr+1,catch=catch(om)[,ac(iYr-(2:1))]) 
    tac=hcrSBT1(iYr+1,control=control,
                catch(om)[,ac(iYr-(2:1))],
                cpue[,ac(ac(iYr-(3:1)))])
      
    
    #### Operating Model update
    om =fwd(om,catch=tac,sr=eq,residual=srDev,effort_max=mean(maxF))
    
    print(plot(window(om,end=iYr+interval)))
    }
  cat('==\n')
  
  return(om)}
