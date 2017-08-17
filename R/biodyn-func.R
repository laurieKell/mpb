jackknifeKey=function(object){
  ## Points that were jackknifed
  orig=ldply(object@indices,function(x) {res=as.data.frame(orig(x),drop=TRUE)
  res[!is.na(res$data),]})
  sim =ldply(object@indices,function(x) {res=as.data.frame(x,drop=TRUE)
  res[is.na(res$data),]})
  key=merge(sim,orig,by=c(".id","year"))[,-4]
  names(key)[c(1,4)]=c("name","obs")
  
  key=transform(key,iter=seq(length(unique(iter))))
  key}

press=function(object){
  key=jackknifeKey(object)
  
  prsdl=mdply(key,function(name,year,iter,obs) {
    hat=stock(object,0.5)[,ac(year),,,,iter]%*%
      params(object)[paste("q",(as.numeric(factor(name))),sep=""),iter]
    
    c(log(obs/hat))})
  
  names(prsdl)[5]="residual.p"
  prsdl}

cooksdFn=function(x){
  dgs=subset(diags(x),!is.na(obs))[,c("year","obs")]
  dat=merge(dgs,as.data.frame(stock(x,0.5),drop=TRUE))
  
  dat=data.frame(year    =dgs$year,
                 cooksd  =cooks.distance(lm(obs~data-1,dat)),
                 hat     =hatvalues(lm(obs~data-1,dat)),
                 residual=residuals(lm(obs~data-1,dat)))
  dat}

dif<-function(object){
  n=seq(length(dim(object))-1)
  sweep(as(object,is(object)[2]),n,orig(object),"-")}

relDif<-function(object){
  n=seq(length(dim(object))-1)
  sweep(sweep(as(object,is(object)[2]),n,orig(object),"-"),n,orig(object),"/")}

dfbetaFn=function(object)
  relDif(refs(object))

refs=function(object){
  
  ##sims
  rf=FLPar(refpts(object))
  names(rf)[1]="params"
  pr=FLPar(params(object))
  sim=rbind(pr,rf)
  
  b =FLQuant(stock(  object)[,ac(dims(catch(object))$maxyear)])
  f =FLQuant(harvest(object)[,ac(dims(catch(object))$maxyear)])
  
  bbmsy=b%/%refpts(object)["bmsy"]
  ffmsy=f%/%refpts(object)["fmsy"]
  ymsy =catch(object)[,ac(dims(catch(object))$maxyear)]%/%refpts(object)["msy"]
  
  #bug
  #res=as.data.frame(FLQuants(b=b,f=f,bbmsy=bbmsy,ffmsy=ffmsy,ymsy=ymsy),drop=T)
  #names(res)[3]="params"
  #res=as(res,"FLPar")
  
  res=model.frame(FLQuants(b=b,f=f,bbmsy=bbmsy,ffmsy=ffmsy,ymsy=ymsy),drop=T)
  res=FLife:::mf2FLPar(res)
  #names(res)[3]="params"
  #res=as(res,"FLPar")
  sim=rbind(sim,res)
  
  ##orig
  rf=orig(refpts(object))
  names(rf)[1]="params"
  pr=orig(params(object))
  orig=rbind(pr,rf)
  
  b =orig(stock(  object))[,ac(dims(catch(object))$maxyear)]
  f =orig(harvest(object))[,ac(dims(catch(object))$maxyear)]
  
  bbmsy=b%/%orig(refpts(object))["bmsy"]
  ffmsy=f%/%orig(refpts(object))["fmsy"]
  ymsy =catch(object)[,ac(dims(catch(object))$maxyear)]%/%orig(refpts(object))["msy"]
  
  res=as.data.frame(FLQuants(b=b,f=f,bbmsy=bbmsy,ffmsy=ffmsy,ymsy=ymsy))[,c("iter","data","qname")]
  names(res)[3]="params"
  res=as(res,"FLPar")
  orig=rbind(orig,res)
  
  FLParJK(sim,orig=orig)}


jackSmry<-function(object) {
  
  orig=object@orig
  sim =FLPar(object)
  
  nms =names(dimnames(orig))
  idx =seq(length(nms))[nms != 'iter']
  n   =dims(sim)$iter 
  
  mn  =orig
  u   =sim
  mnU =apply(u, idx, mean)   
  
  SS =apply(sweep(FLPar(sim), idx, mnU,"-")^2, idx, sum)
  
  bias = (n-1)*(mn-mnU)
  se   = sqrt(((n-1)/n)*SS)
  
  influence=2*se/sqrt(n)
  
  cov  =FLPar(cov(model.frame(u)[,dimnames(u)[[1]]])*(n-1)*(n-1)/n)
  cor  =FLPar(cor(cov[drop=T]))
  
  res=FLPars("hat"=mn, "mean"=mnU, "se"=se, "cv"=se%/%mnU,
             "bias"=bias, "biasRel"=bias/mn,"influence"=influence,
             "cov"=cov,"cor"=cor)
  return(res)}
