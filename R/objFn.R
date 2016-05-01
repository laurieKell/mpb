utils::globalVariables(c("aaply"))

logLikeFn<-function(residual){
  
  SS  =sum(residual^2,   na.rm=TRUE)
  n   =sum(!is.na(residual))
  se  =(SS/n)^.5
  res =(log(1/(2*pi))-n*log(se)-SS/(2*se^2))/2
  
  return(res)}

objLog=function(residual){
  
  dmn=list(params=c("ss","lav","ll","n"),
           iter  =dimnames(residual)$iter)
  rtn=FLPar(array(as.numeric(NA),dim=unlist(laply(dmn, length)),dimnames=dmn))
  
  rtn["ss", ]=aaply(residual%*%residual,6,sum,na.rm=TRUE)
  rtn["lav",]=aaply(abs(residual),      6,sum,na.rm=TRUE)
  rtn["ll",] =aaply(residual,6,function(x) logLikeFn(x))
  rtn[ "n",] =aaply(residual,6,function(x) sum(!is.na(x)))
  
  rtn}

calcQFLQuant<-function(stock,index,when=0.5,na.rm=T){
  
  when =pmin(pmax(when,0),1)
  stock=stock[,-dim(stock)[2]]*(1-when)+stock[,-1]*when

  yrs=dimnames(index)$year
  yrs=yrs[yrs%in%dimnames(stock)$year]
  index=index[,yrs]
  stock=stock[,yrs]
  
  q    =exp(apply(log(index)%-%log(stock),6,sum, na.rm=na.rm)/dim(index)[2])
  #sigma=calcSigma(log(index),log(stock%*%q))
              
  return(q)}

calcObjFn=function(object,index="missing",when=0.5,calcq=TRUE,na=10e6){
  if (missingArg(index))
    index=FLQuants(index(object,FALSE))
  else if (is.FLQuant(index))
    index=FLQuants(index)
  
  if (!calcq)
      q=params(object)[grep("q",dimnames(params(object))$params)]
  else
      q=FLPar(q=laply(index,function(x) calcQFLQuant(stock(object),x)))
  
  res=FLPars(mlply(data.frame(component=seq(length(index))), function(component){
    obs=index[[component]]%/%q[component]
    hat=stock(object,when=when)
    yrs=dimnames(obs)$year
    yrs=yrs[yrs%in%dimnames(hat)$year]
    residuals=log(obs[,yrs]/hat[,yrs])
    residuals[is.na(residuals)]=na              
    objLog(residuals)}))
    
  attributes(res)=NULL
  
  res}

