#om.=fwd(om,f=fbar(om)[,-1,1]%=%refpts(eq)["msy","harvest"],sr=eq)

getPath <- function(file) {
  if (!grepl(.Platform$file.sep,file))
    res <- getwd()
  else
    res <- substr(file,1,max(gregexpr(.Platform$file.sep,file)[[1]])-1)
  return(res)}

getFile <- function(file) {
  res <- substr(file,max(gregexpr(.Platform$file.sep,file)[[1]])+1,
                nchar(file))
  return(res)}

getBiodyn<-function(x,what="biodyn"){
  y1=ls()
  load(x)
  y2=ls()
  y2=y2[!(y2%in%c(y1,"y1"))]
  
  rtn=y2[maply(y2,function(x,y=what) y%in%is(get(x)))]
  rtn=get(rtn[1])
  switch(what,
        FLStock=as(rtn,"biodyn"),
        list  =jabba2biodyn(rtn),
        rtn)}

getOM<-function(x,what="FLStock"){
  y1=ls()
  load(x)
  y2=ls()
  y2=y2[!(y2%in%c(y1,"y1"))]
  
  rtn=y2[maply(y2,function(x,y=what) y%in%is(get(x)))]
  get(rtn[1])}

bdTimeSeries<-function(om,eq){
  FLQuants(om,"Stock." =function(x) apply(stock(x),2,mean)/refpts(eq)["msy","biomass"],
              "Blim."  =function(x) apply(stock(x),2,mean)/refpts(eq)["virgin","biomass"],
              "Stock"  =function(x) ssb(x)[,,"F"]/refpts(eq)["msy","ssb"],
              "Blim"   =function(x) ssb(x)[,,"F"]/refpts(eq)["virgin","ssb"],
              "Harvest"=function(x) apply(catch(x)/stock(x),2,mean)/(refpts(eq)["msy","yield"]/refpts(eq)["msy","biomass"]),
              "Catch"  =function(x) 0.5*apply(catch(x),2,sum)/refpts(eq)["msy","yield"])}

refTimeSeries<-function(om,eq,mp,historical=NULL){
  om=model.frame(bdTimeSeries(om,eq),drop=TRUE)

  mp=model.frame(mcf(FLQuants(mp,"Stock"  =function(x) stock(  x)/refpts(x)["bmsy"],
                                 "Blim"   =function(x) stock(  x)/params(x)["k"],
                                 "Harvest"=function(x) harvest(x)/refpts(x)["fmsy"],
                                 "Catch"  =function(x) catch(  x)/refpts(x)["msy"])),drop=TRUE)

  rtn=rbind.fill(cbind(what="OM",om),
                 cbind(what="MP",mp))

  if ("FLStock"%in%is(historical)) historical=FLStocks(list(Historical=historical))  
  if (!(is.null(historical))){
    rtn2=ldply(historical,function(x) model.frame(bdTimeSeries(x,eq),drop=TRUE))
    rtn=rbind.fill(rtn,transform(rtn2,what=.id)[,-1])}
  
  rtn}

smryBiodyns<-function(jb,om,eq){
  mp=jabba2biodyn(jb)
  rf=c(refpts(mp),params(mp)[c("k","r")])
  names(rf)=c("msy","fmsy","bmsy","virgin","r")
  rf=t(as.data.frame(rf))
  
  rf2=c(refpts(eq)["msy",c("yield","ssb","biomass","harvest")],refpts(eq)["virgin",c("ssb","biomass")])
  names(rf2)=c("msy","ssbmsy","bmsy","fmsy","virginssb","virgin")
  rf2=t(as.data.frame(rf2))
  
  rbind.fill(
    cbind(cbind("what"="OM", model.frame(bdTimeSeries(om, eq),drop=TRUE)),rf2),
    cbind(mutate(model.frame(FLQuants(mp,
                                      "Blim"   =function(x) stock(  x)/params(x)["k"],
                                      "Stock"  =function(x) stock(  x)/refpts(x)["bmsy"],
                                      "Harvest"=function(x) harvest(x)/refpts(x)["fmsy"],
                                      "Catch"  =function(x) catch(x)/refpts(x)["msy"]),drop=T)),rf,what="JABBA"))}

if (FALSE){
  dirMy="/home/laurence-kell/Desktop/papers/albio/inputs/ss"
  load("/home/laurence-kell/Desktop/submission/albio-1/results/baseMainCorners.RData")
  
  library(plyr)
  library(dplyr)
  library(reshape)
  library(ggplot2)
  library(JABBA)
  library(mpb)

  biodyns=mpb:::biodyns
  
  run="161-M0202_sigmaR0.4_steepness0.8_cpuecv0.3_ess50_llq1.0000_llselLog_number161"
  x  =c("mp_1_1999_.RData","mp_1_2002_.RData","mp_1_2005_.RData","mp_1_2008_.RData","mp_1_2011_.RData")
  
  bds=biodyns(mlply(x,function(x) getBiodyn(file.path(dirMy,run,"Schaeffer",1,x))))
  
  mse=getOM(file.path(dirMy,run,"mse.Schaeffer.1.1.RData"))
  om =getOM(file.path(dirMy,run,"om.RData"))
  eq =getOM(file.path(dirMy,run,"om.RData"),"FLBRP")
  dat=rbind.fill(
        cbind("what"="OM", as.data.frame(bdTimeSeries(om, eq),drop=TRUE)),
        cbind("what"="MSE",as.data.frame(bdTimeSeries(mse,eq),drop=TRUE)),
        mutate(ldply(bds,function(x){ as.data.frame(FLQuants(x,
                                          "Stock"  =function(x) stock(  x)/refpts(x)["bmsy"],
                                          "Harvest"=function(x) harvest(x)/refpts(x)["fmsy"],
                                          "Catch"  =function(x) catch(  x)/refpts(x)["msy"]),drop=T)}),
                            what=paste("MP",seq(1999,2011,3))[as.numeric(.id)]))[,-5]
     
  ggplot(dat)+
    geom_hline(aes(yintercept=1),col="red")+
    geom_line(aes(year,data,col=what))+
    facet_grid(qname~.,scale="free")
}

if (FALSE){
  library(plyr)
  library(dplyr)
  library(reshape)
  library(ggplot2)
  library(JABBA)
  library(mpb)
  library(mydas)
  
  biodyns=mpb:::biodyns
  oemFn=mydas:::oemFn
  jabba2biodyn=mydas:::jabba2biodyn

  dat=mdply(subset(scen, what%in%c("base","main"))$run, function(run) {
  #dat=mdply(scen$run, function(run) {
    dir="/home/laurence-kell/Desktop/papers/albio/inputs/ss"
    om =getOM(file.path(dir,run,"om.RData"))
    eq =getOM(file.path(dir,run,"om.RData"),"FLBRP")
    
    #om =fwd(om,f=fbar(om)[,ac(1999:2014),1]*2.5,sr=eq)
    dat     =oemFn(om,2014)
    sa=list(assessment="mp",
            scenario  ="Schaeffer",
            model.type="Pella",
            r.prior   =c(0.5,0.2),
            BmsyK     =0.37,
            K.prior   =c(300000,1),
            psi.prior =c(1,0.05),
            sigma.proc=0.1, 
            
            proc.dev.all=!TRUE,
            sigma.est   =!TRUE, 
            fixed.obsE  =0.02)  
    sa$catch=transmute(dat,Yr=year,Total=catch)
    sa$cpue =transmute(dat,Yr=year,Inxdex=index) 
    sa$se   =sa$cpue
    sa$se[,-1]=0.1
    
    jb=do.call("build_jabba",sa)
    jb=fit_jabba(jb,
                 init.values=TRUE,
                 init.K=sa$K.prior[1],
                 init.r=sa$r.prior[1],
                 init.q=rep(mean(dat$index/dat$ssb),dim(sa$cpue)[2]-1)/4,
                 ni =5500,
                 nt =1,
                 nb =500,
                 nc =2)
    mp=jabba2biodyn(jb)
    
    rf=c(refpts(mp),params(mp)[c("k","r")])
    names(rf)=c("msy","fmsy","bmsy","virgin","r")
    rf=t(as.data.frame(rf))
    
    rf2=c(refpts(eq)["msy",c("yield","ssb","biomass","harvest")],refpts(eq)["virgin",c("ssb","biomass")])
    names(rf2)=c("msy","ssbmsy","bmsy","fmsy","virginssb","virgin")
    rf2=t(as.data.frame(rf2))
    
    dat=rbind.fill(
      cbind(cbind("what"="OM", model.frame(bdTimeSeries(om, eq),drop=TRUE)),rf2),
      cbind(mutate(model.frame(FLQuants(mp,
                                "Stock"  =function(x) stock(  x)/refpts(x)["bmsy"],
                                "Harvest"=function(x) harvest(x)/refpts(x)["fmsy"],
                                "Catch"  =function(x) catch(x)/refpts(x)["msy"]),drop=T)),rf,
            what="JABBA"))
    dat})

  ggplot(melt(dat[,1:6],id=c("X1","what","year")))+
    geom_hline(aes(yintercept=1),col="red")+
    geom_line(aes(year,value,col=what,group=paste(what,X1)))+
    facet_grid(variable~.,scale="free")
  }


if (FALSE){
  
  dir ="/home/laurence-kell/Desktop/papers/albio/inputs/ss/172-M0303_sigmaR0.4_steepness0.9_cpuecv0.3_ess50_llq1.0000_llselLog_number172"
  runs=paste("mse",adply(key,1,paste,collapse=".")[,4],".RData",sep="")
  x=file.path(dir,runs)
  
  dat=mdply(x, function(x)  as.data.frame(bdTimeSeries(getOM(x), getOM(file.path(getPath(x),"om.RData"),"FLBRP")),drop=TRUE))
  dat=cbind(key[dat$X1,],dat)
  
  ggplot(subset(dat,year>1990&dev==0))+
    geom_line(aes(year,data,group=X1,col=ac(ftar)))+
    facet_grid(qname~model)
}
  
  
  