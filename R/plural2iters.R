utils::globalVariables('FLBRPs')
utils::globalVariables('FLBRP')
utils::globalVariables('refpts<-')

FLQuants2FLQuant<-function(object){
  mat=ldply(object,as.data.frame)[,-7]
  names(mat)[1]="iter"
  mat$iter=as.numeric(as.factor((mat$iter)))
  
  as.FLQuant(mat)}

flqs2flq<-function(object,flq){
  nit=laply(object,function(x) dims(x)$iter)
  nit=c(0,cumsum(nit))[seq(length(nit))]
  
   mat=mdply(seq(length(nit)),function(i) {
    transform(as.data.frame(slot(object[[i]],flq)),
              iter=as.numeric(as.character(iter))+nit[i])})[,-1]
  
  as.FLQuant(mat)}

FLStocks2FLStock<-function(object){
  dms=ldply(object,function(x) unlist(dims(x)))
  
  dms$min=as.numeric(dms$min)
  dms$max=as.numeric(dms$max)
  dms$minyear=as.numeric(dms$minyear)
  dms$maxyear=as.numeric(dms$maxyear)
  
  dmns=list(age =min(dms$min):max(dms$max),
            year=min(dms$minyear):max(dms$maxyear))
  
  nts=length(object)
  res=propagate(FLStock(catch.n=FLQuant(NA,dimnames=dmns)),nts)
  
  minyear=unique(laply(object,function(x) dims(x)$minyear))
  if (length(minyear)>1)
    object=FLStocks(llply(object,window,start=min(minyear)))
  
  slt=getSlots("FLStock")  
  for (s in names(slt[slt%in%"FLQuant"]))
    slot(res,s)=flqs2flq(object,s)
  
  units(harvest(res))=units(harvest(object[[1]]))
  
  res}


FLBRPs2FLBRP<-function(object){
  
  if (dim(unique(laply(object,function(x) ac(model(x)))))[1]>1)
    stop("model has to be the same in all objects")
  
  dms=ldply(object,function(x) unlist(dims(x)))
  
  dms$min=as.numeric(dms$min)
  dms$max=as.numeric(dms$max)
  dms$minyear=as.numeric(dms$minyear)
  dms$maxyear=as.numeric(dms$maxyear)
  
  dmns=list(age =min(dms$min):max(dms$max),
            year=min(dms$minyear):max(dms$maxyear))
  
  obs=c("fbar.obs","landings.obs","discards.obs","rec.obs","ssb.obs","stock.obs","profit.obs")   
  
  minyear=unique(c(laply(object, function(x) maply(obs,function(y) dims(slot(x,y))$minyear))))
  if (length(minyear)>1)
    object=FLBRPs(llply(object,function(x) {
        for (i in obs)
           slot(x,i)=window(slot(x,i),start=min(minyear))
      x}))
  
  res=FLBRP()
  slt=getSlots("FLBRP")
  for (s in names(slt[slt%in%"FLQuant"]))
    slot(res,s)=flqs2flq(object,s)
  
  model(res)=model(object[[1]])
  
  par=aperm(laply(object,params),c(2:1))
  names(dimnames(par))[2]="iter"
  dimnames(par)[[2]]=seq(dim(par)[2])
  params(res)=FLPar(as.array(par))
  
  refs=aperm(laply(object,refpts),c(2:3,1))
  names(dimnames(refs))[3]="iter"
  dimnames(refs)[[3]]=seq(dim(refs)[3])
  refpts(res)=FLPar(as.array(refs))
  
  name( res)=""
  desc( res)="created from FLBRPs"
  rng=ldply(object,range)
  r=c(min=unlist(min(rng$min)),
      max=unlist(max(rng$max)),
      minfbar=unlist(max(rng$minfbar)),
      maxfbar=unlist(min(rng$maxfbar)),
      plusgroup=NA)
  
  if (length(unique(rng$plusgroup))==1)
    r["plusgroup"]=min(rng$plusgroup)
  
  range(res)=r
  
  units(landings.sel(res))=units(landings.sel(object[[1]]))
  units(discards.sel(res))=units(discards.sel(object[[1]]))
  
  res}

#fls=FLStocks2FLStock(om)
#flb=FLBRPs2FLBRP(eql)
