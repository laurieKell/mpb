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

calcObjFn=function(object,index="missing",when=0.5){
  if (missingArg(index))
    index=FLQuants(index(object,FALSE))
  else if (is.FLQuant(index))
    index=FLQuants(index)
  
  q=params(object)[grep("q",dimnames(params(object))$params)]
  
  res=FLPars(mlply(data.frame(component=seq(length(index))), function(component){
    obs=index[[component]]%/%q[component]
    hat=stock(object,when=when)
    objLog(log(obs/hat))}))
    
  attributes(res)=NULL
  
  res}
