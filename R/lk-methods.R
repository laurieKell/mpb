#' resample
#' 
#' @description 
#' Resamples from a distribution to generate a frequency table
#'       
#' @param object an \code{FLQuant}
#' @param dim dimensions of \code{FLQuant} to sample within, i.e. if equal to \code{2} then sampling is within \code{year}
#' @param size sum of frequency disribution
#' @param replace sample with replacement?, defaults to \code{FALSE} only useful if values in \code{object} are integers
#' @param ... any other parameters
#' @return a \code{FLQuant} with simulated frequency distribution
#' @export
#' @docType methods
#' @rdname lk-methods-resample
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
  