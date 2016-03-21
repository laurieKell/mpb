setGeneric('plotCc',  function(data,...)     standardGeneric('plotCc'))
setMethod('plotCc', signature(data='FLQuants'),
          function(data,...) plotCcFn(data,...))

plotCcFn<-function(data,...){            
  cc=mdply(expand.grid(a=names(cpue),b=names(cpue)),
           function(a,b){
             res=model.frame(mcf(FLQuants(cpue[c(a,b)])))
             res=subset(res,!is.na(res[,7])&!is.na(res[,8]))
             
             if (dim(res)[1]>10){
               res=data.frame(lag=-10:10,data=ccf(res[,7],res[,8],plot=F,
                                                  lag.max=10)$acf)
               return(res)}else{return(NULL)}})
  
  # cc=transform(subset(cc,a%in%c("us","jll1","jll2")&b%in%c("us","jll1","jll2")),
  #                     a=factor(a,levels=rev(names(cpue))))
  
  ggplot(cc)+
    geom_linerange(aes(x=lag,ymin=0,ymax=data))+
    facet_grid(a~b)+
    geom_vline(aes(xintercept=0))                                   
  }                                            
