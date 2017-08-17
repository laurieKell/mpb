scale<-function(x,y,...){
  args=list(...)
  
  if (length(args)==0) group=rep(1,length(x)) else group=args[[1]]  
  
  gm=gam(y~lo(x)+group,data=data.frame(x=x,y=y,group=group))
  
  res=data.frame(hat =predict(gm),
                 y     =gm$y,
                 x     =x,
                 group =group,
                 scl   =c(0,coefficients(gm)[-(1:2)])[as.numeric(as.factor(group))]
  )
  res$y  =res$y  -res$scl
  res$hat=res$hat-res$scl
  
  if (length(args)==1) names(res)[4]=names(args)[1]  
  
  names(res)[1:4]=c("hat","data","year","name")
  res[,c(4,3,2,1)]}

plotIndexFn<-function(object,...){
  
  dat=FLQuants(llply(object,stdz)) 
  dat=subset(as.data.frame(dat),!is.na(data))
  rng=range(dat$year)
  dat=with(dat,scale(year,data,name=qname))
  
  ggplot(dat)+
    geom_line( aes(year,hat),data=dat[,c("year","hat")],col="grey60")+
    geom_line( aes(year,hat))+
    geom_line( aes(year,data))+
    geom_point(aes(year,data),shape=21,fill="grey")+
    facet_wrap(~name)
  }


plotIndexResidualFn<-function(object,...){
  
  dat=FLQuants(llply(object,stdz)) 
  dat=subset(as.data.frame(dat),!is.na(data))
  rng=range(dat$year)
  dat=with(dat,scale(year,data,name=qname))
  dat2=merge(dat,expand.grid(year=rng[1]:rng[2],name=unique(dat$name)),by=c("year","name"),all.y=TRUE)
  
  ggplot(transform(dat,residual=data-hat))+
    geom_hline(aes(yintercept=0))+
    geom_point(aes(year,residual),position=position_dodge(width = 1),shape=21,fill="grey")+
    geom_line(aes(year,residual),alpha=0.5)+
    geom_linerange(aes(year,ymin=0,ymax=residual),position=position_dodge(width = 1))+
    facet_wrap(~name)
  
}

#' @title plotIndex
#' 
#' @description 
#'
#' @aliases plotIndex,biodyn,missing-method plotIndex,aspic,missing-method  plotIndex,FLQuants,missing-method
#' 
#' @param data an object of class \code{FLQuants} \code{biodyn} \code{aspic} 
#' 
#' @return an \code{ggplot2} object
#' 
#' @export
#' @rdname plotIndex
#' 
#' @examples
#' \dontrun{
#' 
#' } 
setMethod('plotIndex', signature(data='FLQuants'),
          function(data,...) plotIndexFn(data,...))

setMethod('plotIndex', signature(data='aspic'),
          function(data,...){
            
            u=FLQuants(dlply(data,.(name), with,
                             as.FLQuant(data.frame(year=year,data=index))))
            plotIndexFn(u,...)})

setMethod('plotIndex', signature(data='biodyn'),
          function(data,...){
            
            plotIndexResidualFn(data@index,...)})

setMethod('plotIndex', signature(data='aspics'),
          function(data,facet=facet_wrap(~qname,ncol=1,scales="free_y"),...){
            
            u=ldply(data,index)
            u=u[!duplicated(u[,c("name","year")]),]
            u=FLQuants(dlply(u,.(name), with,
                             as.FLQuant(data.frame(year=year,data=index))))
            plotIndexFn(u,...)})

#' @title plotIndexResidual
#' 
#' @description 
#'
#' @aliases plotIndexResidual,biodyn,missing-method plotIndexResidual,aspic,missing-method  plotIndexResidual,FLQuants,missing-method
#' 
#' @param data an object of class \code{FLQuants} \code{biodyn} \code{aspic} 
#' 
#' @return an \code{ggplot2} object
#' 
#' @export
#' @rdname plotIndexResidual
#' 
#' @examples
#' \dontrun{
#' 
#' } 
setMethod('plotIndexResidual', signature(data='FLQuants'),
          function(data,...) plotIndexFn(data,...))

setMethod('plotIndexResidual', signature(data='aspic'),
          function(data,...){
            
            u=FLQuants(dlply(data,.(name), with,
                             as.FLQuant(data.frame(year=year,data=index))))
            plotIndexResidualFn(u,...)})

setMethod('plotIndexResidual', signature(data='biodyn'),
          function(data,...){
            
            plotIndexResidualFn(data@index,...)})

setMethod('plotIndexResidual', signature(data='aspics'),
          function(data,facet=facet_wrap(~qname,ncol=1,scales="free_y"),...){
            
            u=ldply(data,index)
            u=u[!duplicated(u[,c("name","year")]),]
            u=FLQuants(dlply(u,.(name), with,
                             as.FLQuant(data.frame(year=year,data=index))))
            plotIndexResidualFn(u,...)})

if (FALSE){
object=FLQuants(dlply(uN,.(name), with, as.FLQuant(data.frame(year=year,data=index))))
plotIndex(object)
plotIndexResidual(object)
}
