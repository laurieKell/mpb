fn<-function(id,ts=NULL,...){
  print(id)
  args=list(...)
  
  inp=subset(ts,assessid==id)
  inp=inp[order(inp$year),]
  inp=with(inp,list(obsC=catch,
                    timeC=year,
                    obsI =biomass,
                    timeI=year))
  
  parNm=c("log"=c("logn","logm","logK","logq"),
          "par"=c("n",   "m",   "K",   "q"))
  
  #n: Parameter determining the shape of the production curve as in the generalised form of Pella & Tomlinson (1969).
  #m: Log of maximum sustainable yield.
  #K: Log of carrying capacity.
  #q: Log of catchability vector.
  #sdb: Log of standard deviation of biomass process error.
  #sdf: Log of standard deviation of fishing mortality process error.
  #sdi: Log of standard deviation of index observation error.
  #sdc: Log of standard deviation of catch observation error.
  print(names(args))
  for (i in names(args)){
    inp$priors[[i]]=c(args[[i]],1e-3)   
    inp$phase[[i]] =-1}
  
  res=try(fit.spict(inp))
  
  return(res)}

fnFLQuant<-function(x,val="logB"){
  dat=data.frame(year  =c(rep(seq(x$inp$timerange[1],x$inp$timerange[2]),each=16)),
                 season=c(rep(seq(16),x$inp$timerange[2]-x$inp$timerange[1]+1)),
                 data  =exp(x$report[[val]])[-length(x$report[[val]])])
  as.FLQuant(dat)}

fnMSY<-function(x){
  with(x$report,FLPar(msy=MSY,bmsy=Bmsy,fmsy=Fmsy))}

fnKobe<-function(x){
  FLQuants(stock  =fnFLQuant(x,"logBBmsy"),
           harvest=fnFLQuant(x,"logFFmsy"))}

fnPf<-function(x){
  FLPar(r=x$value["r"],K=x$value["K"],p=x$value["p"])}

