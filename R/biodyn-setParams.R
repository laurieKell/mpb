utils::globalVariables(c('ldply','melt','variable'))
globalVariables("lambda")
utils::globalVariables(c("eql"))
utils::globalVariables('optimise')

setMethod('setParams<-', signature(object='biodyn',value='data.frame'), function(object,value) {
  nms=c(modelParams(as.character(object@model)),'b0')
  object@params=object@params[nms]
  
  object@params =setQ(object,value)
  
  return(object)})

setMethod('setParams<-', signature(object='biodyn',value='FLPar'), function(object,value) {
  object@params=value 
  
  return(object)})

setMethod('setParams<-', signature(object='biodyn',value='FLQuant'), function(object,value) {

  nms=c(modelParams(tolower(as.character(object@model))),'b0')

  #value=FLCore::apply(value,2,mean)
  object@params =setQ(object,value)
  
  return(object)})

setMethod('setParams<-', signature(object='biodyn',value='FLQuants'), function(object,value) {
  nms=c(modelParams(as.character(object@model)),'b0')
  object@params=object@params[nms]
  
  object@params =setQ(FLCore::iter(object,1),
                             FLQuants(lapply(value, function(x) FLCore::apply(x,2,mean,na.rm=T))))
  
  nms=dimnames(params(object))$param
  n  = as.numeric(summary(substr(nms,1,1)=='q')['TRUE'])
  nms[substr(nms,1,1)=='q']    =paste('q',    seq(n),sep='')
  nms[substr(nms,1,5)=='sigma']=paste('sigma',seq(n),sep='')
  
  dimnames(params(object))$params=nms
  
  return(object)})


# setMethod('setParams<-', signature(object='biodyn',value='FLQuants'), function(object,value,msy=TRUE) {
#   nms=c(modelParams(as.character(object@model)),'b0')
#   object@params=object@params[nms]
#   
#   #object@params =setQ(FLCore::iter(object,1),FLQuants(lapply(value, function(x) FLCore::apply(x,2,mean,na.rm=T))))
#   object@params =setQ(object,value)
#   
#   nms=dimnames(params(object))$param
#   n  = as.numeric(summary(substr(nms,1,1)=='q')['TRUE'])
#   nms[substr(nms,1,1)=='q']    =paste('q',    seq(n),sep='')
#   nms[substr(nms,1,5)=='sigma']=paste('sigma',seq(n),sep='')
#   
#   dimnames(params(object))$params=nms
#   
#   return(object)})

getP<-function(bmsy,k,p=c(0.001,5)){
  optimise(function(p,bmsy,k) 
    (bmsy-k*(1/(1+p))^(1/p))^2, p, bmsy=bmsy,k=k)$minimum}

getRK<-function(msy,bmsy,k,p=c(0.0001,5)){
  p=getP(bmsy,k,p)
  k=bmsy/((1/(1+p))^(1/p))
  r=msy/(k*(1/(1+p))^(1/p+1))
  
  FLPar(c(r=r,k=k,p=p))}

#getRK(100,1000,2000)

# setMethod('setParams<-', signature(object='biodyn',value='FLBRP'), function(object,value,msy=TRUE) {
#   
#   #if (!(is(value)%in%'FLBRP')) return(object)
#   
#   if (msy){
#     params(object)[c("r","k","p")]=getRK(value@refpts["msy",   "yield"],
#                                          value@refpts["msy",   "biomass"],
#                                          value@refpts["virgin","biomass"])
#     object@priors[c("r","k","p"),"a"]=params(object)[c("r","k","p")]
#   }else{
#     require(popbio)
#     params(object)["k"]=value@refpts["virgin","biomass"]
#     params(object)["p"]=getP(c(value@refpts["msy","biomass"]),
#                              c(params(object)["k"]),c(0.001,5))
#     params(object)["r"]=log(lambda(leslie(eql,fbar=c(value@refpts["crash","harvest"]))))
#   }
#   
#   object=fwd(object,catch=catch(object)[,-1])
#   
#   return(object)})

