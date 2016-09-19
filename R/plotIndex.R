utils::globalVariables(c("gam","geom_smooth","theme_bw"))
utils::globalVariables('hat')

setMethod('plotIndex', signature(data='FLQuants'),
          function(data,
                   facet=facet_wrap(~qname,ncol=1,scale="free_y"),...){

  require(gam)
  
  u=subset(as.data.frame(data,drop=TRUE),!is.na(data))
  u=cbind(u,hat=predict(gam(log(data)~lo(year)+qname,data=u)))
  
  ggplot(u)                                  +   
    geom_point(aes(year,data,col=qname))     +
    geom_line( aes(year,exp(hat)),col="red") +
    theme(legend.position="none")            +
    geom_smooth(aes(year,data),se=FALSE)     +           
    theme_bw()                               +
    theme(legend.position="none")            +
    facet
})

setMethod('plotIndex', signature(data='aspic'),
          function(data,facet=facet_wrap(~qname,ncol=1,scale="free_y"),...){
            
            u=FLQuants(dlply(data,.(name), with,
                             as.FLQuant(data.frame(year=year,data=index))))
            plotIndex(u,facet,...)})

setMethod('plotIndex', signature(data='aspics'),
          function(data,facet=facet_wrap(~qname,ncol=1,scale="free_y"),...){
    
    u=ldply(data,index)
    u=u[!duplicated(u[,c("name","year")]),]
    u=FLQuants(dlply(u,.(name), with,
                             as.FLQuant(data.frame(year=year,data=index))))
    plotIndex(u,facet,...)})