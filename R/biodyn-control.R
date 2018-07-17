utils::globalVariables(c('ldply','melt','variable'))
globalVariables("lambda")
utils::globalVariables(c("eql"))
utils::globalVariables('optimise')

setGeneric('control',     function(object,...)        standardGeneric('control'))
setGeneric('control<-',   function(object,value)      standardGeneric('control<-'))
setGeneric('setControl<-',function(object,...,value)  standardGeneric('setControl<-'))

#' @title control
#' 
#' @description sets initial guess and lower and upper bounds
#'
#' @param object \code{biodyn} object
#' 
#' @return FLPar
#' @export
#' @examples
#' \dontrun{control(biodyn())}
setMethod( 'control',   signature(object='biodyn'),function(object, ...)   object@control)


#' @title setControl<-
#'
#' @description Sets the control slot in a biodyn object given the parameters in the \code{params}
#' slot. The starting values \code{val} are set to those in \code{params} and the \code{min} and
#' \code{max} bounds to .1 and 10 times these.
#'
#' @param object an object of class \code{biodyn}
#' @param value  \code{params} object
#' @param ... any other parameter
#'
#' @seealso \code{\link{controlFn}}
#' 
#' @export
#' @rdname setControl
#'
#' @aliases 
#' setControl<-,biodyn,FLPar-method  
#' setControl<-,biodyn,FLQuant-method  
#' setControl<-,biodyn,FLQuants-method
#' setControl<-,aspic,FLPar-method 
#' 
#' @examples
#' \dontrun{
#' bd=sim()
#' setControl(bd)=params(bd)
#' }
#'  

setMethod( 'setControl<-', signature(object='biodyn',value='FLPar'), function(object,value,min=0.1,max=10.0) {
  
  if (dims(value)$iter>1 & dims(object@control)$iter==1)
    object@control=propagate(control(object),dims(value)$iter)
  
  ctr=object@control
  nms=dimnames(object@params)$params
  
  object@control=FLPar(array(rep(c(1,NA,NA,NA),each=length(nms)), dim=c(length(nms),4,dims(value)$iter), dimnames=list(params=nms,option=c('phase','min','val','max'),iter=seq(dims(value)$iter))))
  nms.=nms[nms %in% dimnames(ctr)$params]
  
  nits=seq(dims(value)$iter)
  object@control[nms.,'phase',nits]=ctr[nms.,"phase",nits]
  nms.=nms[nms %in% dimnames(object@params)$params]
  
  object@control[,"val",nits]=object@params[nms.,nits]
  object@control[,"min",nits]=object@params[nms.,nits]*min
  object@control[,"max",nits]=object@params[nms.,nits]*max
 
  #if (!is.na(any(value[nms]<0)) & any(value[nms]<0))
  #  object@control[nms[value[nms]<0],c('min','max')]=object@control[nms[value[nms]<0],c('max','min')]
  #object@control[nms[value[nms]<0],c('min','max')]=object@control[nms[value[nms]<0],c('max','min')]

  max=object@control["p","min",object@control["p","val"]<0]
  min=object@control["p","max",object@control["p","val"]<0]
  
  object@control["p","max",object@control["p","val"]<0]=max
  object@control["p","min",object@control["p","val"]<0]=min
  
  object@control[-(1:4),"phase"]=-2
  
  prr=object@priors
  nms=c(nms[!substr(dimnames(object@params)$params,1,1) %in% c('s','q')],
        'msy','bmsy','fmsy',
        nms[substr(dimnames(object@params)$params,1,1) %in% c('s','q')])
  object@priors=array(rep(c(0,0,0.3,1),each=length(nms)), dim=c(length(nms),4),   dimnames=list(params=nms,c('weight','a','b','type')))
  nms=dimnames(prr)$params[dimnames(prr)$params %in%  dimnames(object@priors)$params]
  
  object@priors[nms,]=prr[nms,]
  
  return(object)})

.calcSigma <- function(obs,hat=rep(0,length(obs)),na.rm=T){
  n  =length(obs[!is.na(obs+hat)])
  SS =sum((obs-hat)^2,na.rm=na.rm)
  
  return((SS/n)^.5)}

calcB0<-function(index,q,k,nyrB0=3,error='log'){
  if (is.null(nyrB0)) return(params['b0'])
  
  if (error=='log'){
    t.<-sweep(log(index[,1:nyrB0,,,,,drop=FALSE]),c(1,6),q,'/')
    return(qmax(qmin(exp(apply(t.,c(1,6),mean))/k,1),0))}
  if (error=='normal'){
    t.<-sweep(index[,1:nyrB0,,,,,drop=FALSE],c(1,6),q,'/')
    return(qmax(qmin(apply(t.,c(1,6),mean)/k,1),0))}       
}

calcQ<-function(stock,index,error='log',na.rm=T){
  
  stock<-(stock[-length(stock)]+stock[-1])/2
  n    <-length(stock)
  index<-index[seq(n)]
  if (na.rm)
    n=length(seq(n)[!is.na(index+stock)])
  
  res=switch(error,
             normal={q    =sum(stock*index, na.rm=T)/sum(stock*stock, na.rm=na.rm)
             sigma=.calcSigma(index/(q*stock))
             data.frame(q=q,sigma=sigma)
             },
             log   ={q    =exp(sum(log(index)-log(stock), na.rm=na.rm)/n)
             sigma=.calcSigma(log(index),log(q*stock))
             data.frame(q=q,sigma=sigma)},
             cv   ={res   <-sum(index/stock)
             sigma2<-.calcSigma(res,na.rm=na.rm)
             q     <-(-res+(res^2+4*length(index)*sigma2*sum((index/stock)^2)))/(2*length(index)*sigma2)
             data.frame(q=q,sigma=sigma)})
  
  return(res)}

setQ=function(object,value,error='log'){
  if (!is.FLQuant(value))
    names(value)=seq(length(value))
  
  fn=function(value,stock){
    if (dims(value)$iter==1)
      dimnames(value)$iter=1
    if (dims(stock)$iter==1)
      dimnames(stock)$iter=1
    
    if (dims(stock)$iter==1 & dims(value)$iter>1)
      stock=propagate(stock,dims(value)$iter)
    
    if (dims(stock)$iter>1 & dims(value)$iter==1)
      value=propagate(value,dims(stock)$iter)
    
    model.frame(mcf(FLQuants(stock=stock,value=value)))}
  
  res=switch(is(value)[1],
             FLQuant   ={res=fn(value,stock(object));data.frame(name=1,res)},
             #FLQuants  =ldply(value, fn, model.frame(mcf(FLQuants(stock=stock,value=x))),stock=stock(object)),
             FLQuants  =ldply(value, function(x,stock) fn(x,stock), stock=stock(object)),
             data.frame=merge(model.frame(FLQuants('stock'=stock(object))),value,by='year',all=T))
  
  res=res[!is.na(res$iter),]
  
  if (!('name' %in% names(res))) 
    names(res)[1]='name'
  
  res=res[!is.na(res$name),]
  res=ddply(res, .(name,iter), function(x,log) data.frame(calcQ(x$stock,x$value)),log='log')
  its=max(as.numeric(ac(res$iter)))
  
  res.=transform(melt(res,id=c('name','iter')),params=paste(variable,name,sep=''))[,c('params','value','iter')]
  names(res.)[2]='data'
  
  #bug
  res=as(res.,'FLPar')
  #res=FLCore::iter(res,seq(its))
  units(res)='NA'
  res=res.[with(res.,order(iter,params)),]
  
  if (dims(object@params)$iter==1)
    object@params=propagate(object@params,its)
  t.=try(rbind(object@params,FLPar(as.FLQuant(res)[,1,drop=T])),silent=TRUE)
  
  if (is(t., "try-error"))
    t.=as(rbind(as.data.frame(object@params),
                as.data.frame(FLPar(as.FLQuant(res)[,1,drop=T]))),"FLPar")
  
  dmns=dimnames(t.)
  names(dmns)=c('params','iter')
  t.=FLPar(array(t.,dim=unlist(lapply(dmns,length)),dimnames=dmns))
  units(t.)='NA'
  
  object@params=t.
  #object@params=FLPar(rbind(FLPar(object@params),FLPar(res)))
  
  object@params}

setMethod('setControl<-', signature(object='biodyn',value='FLQuant'), 
          function(object,value,min=0.1,max=10.0) {
            
            setParams(object)<-value
            setControl(object,min=min,max=max)<-params(object)
            
            return(object)})


setMethod('setControl<-', signature(object='biodyn',value='FLQuants'), 
          function(object,value,min=0.1,max=10.0) {
            
            setParams(object)<-value
            setControl(object,min=min,max=max)<-params(object)
            
            return(object)})
# 
# setMethod('control<-', signature(object='biodyn',value='FLPar'), function(object,value) {
#   object@control=value
#     
#   return(object)})
# 
# setMethod('control', signature(object='biodyn'), function(object,...) {
#   object@control})

#' @title controlFn
#' 
#' @description A utility function to help set up the \code{control} slot in \code{biodyn} 
#' 
#' @param r a \code{numeric} value with best guess
#' @param k a \code{numeric} value with best guess      
#' @param p a \code{numeric} value with best guess, default=1      
#' @param b0 a \code{numeric} value with best guess default=1
#' @param phaseR a \code{numeric} value for phase, default=1
#' @param phaseK a \code{numeric}  value for phase, default=1,
#' @param phaseP a \code{numeric}  value for phase, default=-1,
#' @param phaseB0 a \code{numeric}  value for phase, default=-1,
#' @param min a \code{numeric} a multipler for the best guess  
#' @param max \code{numeric} a multipler for the best guess
#'
#' @return a \code{control} object
#'   
#' @export
#' @rdname controlFn
#' 
#' @seealso \code{\link{biodyn}}  \code{\link{control}} 
#'   
#' @examples
#' \dontrun{
#' sim()
#'    }
## utility function for setting control object
controlFn=function(r,       k,       p=1,      b0=1,
                   phaseR=1,phaseK=1,phaseP=-1,phaseB0=-1,
                   min=.5,  max=2){ 
  
  dmns=list(params=c('r','k','p','b0'),
            option=c('phase','min','val','max'),  
            iter  =1)
  
  res=FLPar(array(0,unlist(laply(dmns,length)),dmns))
  res[,'val'][]=c(r,k,p,b0)
  res[,'min']=res[,'val']*min
  res[,'max']=res[,'val']*max
  
  res[,'phase']=c(phaseR,phaseK,phaseP,phaseB0)
  
  res}  


#' @title control<-
#'
#' @description sets in \code{biodyn} initial guess and lower and upper bounds
#' 
#' @return \code{biodyn} with new control slot
#' @export
#' 
#' @rdname biodynAccessors
#' @docType methods
#' 
#' @examples
#'  
#' \dontrun{control(biodyn())}
setMethod('control<-', signature(object='biodyn', value='FLPar'),
          function(object, value){
            updateFLComp(object, 'control', value)
            return(object)})   

