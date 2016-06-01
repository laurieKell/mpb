#' resample
#' 
#' @description 
#' Resamples from a distribution to generate a frequency table
#'       
#' @param object an \code{FLQuant}
#' @param dim dimensions of \code{FLQuant} to sample within, i.e. if equal to \code{2} then sampling is within \code{year}
#' @param size sum of frequency disribution
#' @param replace sample with replacement?, defaults to \code{FALSE} only useful if values in \code{object} are integers
#' @params ...
#' @return a \code{FLQuant} with simulated frequency distribution
#' @export
#' @docType methods
#' @rdname lk-methods
#' 
#' @examples
#' \dontrun{
#'   library(FLCore)
#'   data(ple4sex)
#'   
#'   apply(resample(catch.n(ple4sex)[,1:4],2,  100),2,  sum)
#'   apply(resample(catch.n(ple4sex)[,1:4],2:3,100),2:3,sum)
#'   apply(resample(catch.n(ple4sex)[,1:4],1:2,100),1:2,sum)
#'   }
setGeneric('resample',   function(x,...) 
  standardGeneric('resample'))

setMethod('resample',  signature(x='FLQuant'), 
    function(x,dim=2:6,size,replace=T,...) {
      
      resampleFn=function(x,size,replace=T,...){
        
        sampleFn= function(x, ...) x[sample.int(length(x),...)]
        
        y  =c(x)
        nms=seq(length(y))
        
        res=sampleFn(nms,size,replace,prob=y,...)
        res=factor(res,levels=nms)
        
        summary(res)}
      
      res=apply(x,dim,resampleFn,size=size,replace=replace,...)
        
        ## rearranges array to get dims in FLQuant correct
        d=seq(length(dim(x)))
        i=c(d[!(d %in% dim)],dim)
        names(i)=d
        o=as.numeric(names(i[order(i)]))
        
        dmns=dimnames(x)[i]
        res=apply(x,dim,resampleFn,size=100)
        
        res=array(c(res),unlist(lapply(dmns,length)),dimnames=dmns)
        
        as.FLQuant(aperm(res,o))})
  

#dim(apply(m(ple4)[,1:4],1:6, function(x) rbeta(4,x,x)))

setMethod('resample',  signature(x='FLPar'), 
          function(x,dim=2:length(dim(x)),size,replace=T,...) {
            
            resampleFn=function(x,size,replace=T,...){
              
              sampleFn= function(x, ...) x[sample.int(length(x),...)]

              y  =c(x)
              nms=seq(length(y))
              
              res=sampleFn(nms,size,replace,prob=y,...)
              res=factor(res,levels=nms)
              
              summary(res)}
            
            res=apply(x,dim,resampleFn,size=size,replace=replace,...)
            
            ## rearranges array to get dims in FLQuant correct
            d=seq(length(dim(x)))
            i=c(d[!(d %in% dim)],dim)
            names(i)=d
            o=as.numeric(names(i[order(i)]))
            
            dmns=dimnames(x)[i]
            res=apply(x,dim,resampleFn,size=100)
            
            res=array(c(res),unlist(lapply(dmns,length)),dimnames=dmns)
            
            FLPar(aperm(res,o))})

#' mixnorm
#' 
#' @description 
#' Simulates a normal density function from an \code{FLQuant}
#'       
#' @param mn an \code{FLQuant} with means
#' @param sd an \code{FLQuant} with standard deviation
#' @return an \code{FLPar} with expected probability for lengths-at-age
#' @export
#' @docType methods
#' @rdname lk-methods
#' 
#' @examples
#' \dontrun{
#' mn=FLQuant( 1:5,     dimnames=list(age=1:5,year=1991:2000,unit=1:2))
#' sd=FLQuant((1:5)/5,dimnames=list(age=1:5,year=1991:2000))
#' md=mixnorm(mn,sd,(1:20)/4)
#' ggplot(as.data.frame(md))+
#'    geom_line(aes(len,data,group=age))+
#'    facet_wrap(year~unit)
#' ggplot(as.data.frame(apply(md,1,sum)))+
#'    geom_line(aes(as.numeric(as.character(params)),data))     
#'  }
mixnorm=function(n,mn,sd,bin,left=T){
  require(plyr)
  
  res=FLPar(maply(data.frame(len=bin), function(len) n*pnorm(len,mn,sd)))
  
  rtn=res
  
  if (!left){
    rtn[-1]          =res[-1]-res[-dim(res)[1]]
    attributes(rtn)[["boundary"]]="right"
  }else{
    rtn[-dim(res)[1]]=res[-1]-res[-dim(res)[1]]
    rtn[ dim(res)[1]]=1-res[ dim(res)[1]]
    attributes(rtn)[["boundary"]]="left"
  }
  
  if (length(unique(diff(bin)))==1) attributes(rtn)[["bin"]]=unique(diff(bin))
  
  return(rtn)}
  
rLen<-function(len,sd,n,smp=n,samplesize=1000,a=1,b=1){
  
  smp=smp/sum(smp)*samplesize
  res=mdply(data.frame(l=len,s=sd,m=smp),
            function(l,s,c,m)
              rnorm(1,l, s/(m^0.5)))
  sum(n*res$V1)/sum(n)}

