utils::globalVariables(c('ggplot','geom_line','aes','yield',
                         'geom_point','cast','xlab','ylab',
                         'jack.mean','facet_wrap','geom_ribbon',
                         'jack.se','expand_limits','jack.bias',
                         'theme','element_blank','y','se'))

utils::globalVariables(c('laply','ldply','geom_path','.id',
                         'scale_size_manual','scale_linetype_manual'))
utils::globalVariables(c("scale_alpha_manual"))

utils::globalVariables(c('geom_path','scale_size_manual','scale_linetype_manual'))

utils::globalVariables(c('qname','What'))

globalVariables(c("ymax","ymin","CI","50%","ymax","ymin","CI"))

utils::globalVariables(c("gam","geom_smooth","theme_bw"))
utils::globalVariables('hat')

utils::globalVariables(c("x"))
utils::globalVariables('data')
utils::globalVariables(c("geom_linerange","facet_grid","geom_vline"))
utils::globalVariables('ccf')
utils::globalVariables('lag')


setGeneric('plotIndex',      function(data,...)          standardGeneric('plotIndex'))
setGeneric('plotDiags',      function(data,...)          standardGeneric('plotDiags'))
setGeneric('plotCcf',        function(data,...)          standardGeneric('plotCcf'))
setGeneric('plotEql',        function(data,biomass,...)  standardGeneric('plotEql'))
setGeneric('plotProduction', function(data,biomass,...)  standardGeneric('plotProduction'))
setGeneric('plotMSE',        function(x,y,z,...)         standardGeneric('plotMSE'))
setGeneric('plotHcr',        function(object,...)        standardGeneric('plotHcr'))
setGeneric('plotJack',       function(object,...)        standardGeneric('plotJack'))

#' @title plot
#' 
#' @description 
#' Plots time series of biomass, harvest rate and catch for a \code{biodyn} object, using \code{ggplot2}.
#'
#' @aliases plot,biodyns,missing-method
#' @param x an object of class \code{biodyn} 
#' @param y second argument
#' @param worm iter(s) to plot as lines 
#' @param probs numeric vector of probabilities with values in [0,1].  
#' @param na.rm	a logical value indicating whether NA values should be stripped before computation.
#' @param type an integer between 1 and 9 selecting one of the quantile algorithms to be used.
#' @param fn  functions
#' @param facet facet for panelling
#'
#' @importFrom reshape cast
#' 
#' @return an \code{ggplot2} object
#' 
#' @aliases plot,biodyn,missing-method plot,biodyns,missing-method plot,aspics,missing
#' 
#' @export
#' @rdname plot
#' 
#' @examples
#' \dontrun{
#' 
#' bd =sim()
#' plot(bd)
#' } 
setMethod('plot', signature(x='biodyn', y='missing'),
  function(x, y, probs=c(0.95,0.75,0.50,0.25,0.05), na.rm=FALSE, type = 7, 
    worm =NULL,
    fn   =list('Stock'  =function(x) stock(x), 
              'Harvest'=function(x) harvest(x),
              'Yield'  =function(x) catch(x)),
    facet=facet_wrap(~qname,scales='free',ncol=1),...){
    
      ## contains iters
      if (dims(x)$iter>=length(probs)){  
        res=whooow(x,fn,probs)
        res=transform(res,iter=factor(iter,labels= paste(as.integer(probs*100),"%",sep="")))
        if ("quant"%in%names(res))
          res=cast(res,quant+year+unit+season+area+qname~iter,value="data")
        else
          res=cast(res,age+year+unit+season+area+qname~iter,value="data")
        
        nprob=length(probs)
        if (any(0.5 %in% probs)){
          val=grep("%",names(res))

          mdn=res[,c(1:6,val[ (val==grep("50%",names(res)))])]
          
          res=res[,c(1:6,val[!(val==grep("50%",names(res)))])]}
        
        rbn=mdply(data.frame(ymax=dim(res)[2]-seq(nprob/2)+1,
                             ymin=6+seq(nprob/2)),
                    function(ymax,ymin){
                       data.frame(res[,c(1:6)],ymin=res[,ymin],ymax=res[,ymax],
                                  CI=paste(names(res)[ymin],names(res)[ymax],sep="-"))
                    })
          
          p=ggplot(rbn)+
            geom_ribbon(aes(year,ymax=ymax,ymin=ymin,alpha=CI),fill="red")
        
         if (any(0.5 %in% probs))
           p=p+geom_line(  aes(year,`50%`),data=mdn)

      }else{
        res=whooow(x,fn,0.5)
        p  =ggplot(res) + geom_path(aes(x=year,y=data))    }
      
      p = p + facet + 
            expand_limits(y = 0) +
            xlab('Year') + ylab('')
      
      if (length(worm) > 0)
        if (length(worm)<=dims(x)$iter)
          if (!(length(worm)==1 & is.na(worm[1])))  
            p=p+geom_path(aes(year,data,group=iter,colour=iter),
                            data=transform(subset(as.data.frame(FLQuants(lapply(fn,function(f,x) f(x), x=x))),iter %in% worm),iter=factor(iter)))
      
      p}) 
    
setMethod('plot', signature(x='biodyns', y='missing'),
  function(x, y, probs=c(0.95,0.75,0.50,0.25,0.05),type=7,na.rm=FALSE,
    facet=facet_wrap(~qname,scales='free',ncol=1),
           fn=list('Stock'  =function(x) stock(x), 
                   'Harvest'=function(x) harvest(x),
                   'Yield'  =function(x) catch(x)),...)
    {  
    res=ldply(x,function(x) plot(x,fn=fn,probs=probs,type=type,na.rm=na.rm)$data)
    mdn=ldply(x,function(x) whooow(x,fn,probs=.5))
    
    if (names(res)[1]=="X1") names(res)[1]=".id"
    if (names(mdn)[1]=="X1") names(mdn)[1]=".id"
    
    if (length(unique(res$CI))>=1){
      p=ggplot(res)+
        geom_ribbon(aes(year,ymax=ymax,ymin=ymin,alpha=CI,group=paste(.id,CI),fill=.id))+
        geom_line(  aes(year,data,col=.id,group=.id),data=mdn)+
        scale_alpha_manual(values=seq(.1,.75,length.out=length(unique(res$CI))))
      
    }else{
      p=ggplot(res)+
        geom_line(aes(year,data,col=.id,group=.id),data=mdn)
      
    }
    
    p = p + facet + 
      expand_limits(y = 0) +
      xlab('Year') + ylab('')
  
    p})

# @param  \code{fn}, a list of functions that estimate the quantities for plotting
# @param  \code{probs}, a vector specifying the percentiles for plotting, these are c(0.95,0.50,0.05) by default.
# @param  \code{size}, thinkness of percentile lines
# @param  \code{lty}, line type for percentiles
# @param \code{facet}, a layer that determines the facetting of the plot

#' @title plotProduction
#' 
#' @description 
#' Creates a \code{ggplot2} object that plots equilibrium values of biomass, harvest rate and catch against each other.
#' The basic object can then be modified by adding ggpot2 layers.
#'
#' @param  object an object of class \code{biodyn} 
#' @param  biomass optional argument, an FLQuant with biomass at beginning of year 
#' @param ... other arguments
#'
#' @return an \code{ggplot2} object
#' 
#' @export
#' @rdname plotProduction
#'
#' @aliases 
#  plotProduction,biodyn,FLBRP-method 
#' plotProduction,biodyn,FLQuant-method 
#' plotProduction,biodyn,missing-method
#' plotProduction,biodyns,missing-method 
#' 
#' @examples
#' \dontrun{
#'  refpts('logistic',FLPar(msy=100,k=500))
#' }
setMethod('plotProduction',signature(data='biodyn',biomass='missing'),
          function(data,biomass=FLQuant(seq(0,max(params(data)['k']),length.out=101)),...)
            plotProductionfn(data,biomass,...))
setMethod('plotProduction',signature(data='biodyn',biomass='FLQuant'),
          function(data,biomass,...) plotProductionfn(data,biomass,...))
# setMethod('plotProduction',signature(x='biodyn',biomass='FLBRP'),  
#           function(x,biomass,II=FALSE,...) plotProdfn(bd=x,brp=biomass,II=II,...))

setMethod('plotProduction',signature(data='biodyns',biomass='missing'),
          function(data,biomass,...) {
    res=ldply(data,function(x) {
           biomass=FLQuant(seq(0,max(params(x)['k']),length.out=101))
           model.frame(FLQuants(stock=biomass, yield=FLQuant(production(x,biomass))))})

    msy=ldply(data,function(x)
         cast(as.data.frame(refpts(x)),iter~refpts,value='data'))

    ggplot(res)+
        geom_line( aes(stock, yield, group=.id, col=.id))+
        geom_point(aes(bmsy,  msy,   col=.id),size=2,data=msy) +
        xlab('Stock') + ylab('Surplus Production')})

plotProductionfn=function(data,biomass=FLQuant(seq(0,max(params(data)['k']),length.out=101)),...) {
  warn=options()$warn
  options(warn=-1)
  
  if ((dims(data)$iter>1 | dims(params(data))$iter>1) & dims(biomass)$iter==1) 
    biomass=propagate(biomass,max(dims(data)$iter,dims(params(data))$iter))
  
  p <-  ggplot(model.frame(FLQuants(stock=biomass, yield=FLQuant(production(data,biomass))))) +
    geom_line(aes(stock, yield, group=iter, col=iter)) +
    geom_point(aes(bmsy,msy,col=iter),size=2,data=cast(as.data.frame(refpts(data)),iter~refpts,value='data')) +
    xlab('Stock') + ylab('Surplus Production')
  
  options(warn=warn)
  
  p} 

plotProdfn=function(bd,brp,II=FALSE){
  res=cbind(What='Age I',  model.frame(FLQuants('catch'=catch(brp),'stock'=stock(brp)),drop=T))
  
  if (II){
    landings.wt(brp)=stock.wt(brp)*mat(brp)
    res=rbind(res, cbind(What='Age II',   model.frame(FLQuants('catch'=catch(brp),'stock'=stock(brp)),drop=T)))}
  
  if (!is.null(bd)){
    bm =FLQuant(seq(0, max(params(bd)['k']), length.out = 101))
    res=rbind(res,cbind(What='Biomass', 
                        model.frame(mcf(FLQuants(catch=production(bd,bm),stock=bm)),drop=T)))}
  
  ggplot(res)}

plotProdfn=function(bd,brp,II=FALSE){
  res=cbind(What='Age I',  model.frame(FLQuants('catch'=catch(brp),'stock'=stock(brp)),drop=T))
  
  if (II){
    landings.wt(brp)=stock.wt(brp)*mat(brp)
    res=rbind(res, cbind(What='Age II',   model.frame(FLQuants('catch'=catch(brp),'stock'=stock(brp)),drop=T)))}
  
  if (!is.null(bd)){
    bm =FLQuant(seq(0, max(params(bd)['k']), length.out = 101))
    res=rbind(res,cbind(What='Biomass', 
                        model.frame(mcf(FLQuants(catch=production(bd,bm),stock=bm)),drop=T)))}
  
  ggplot(res)}

#' @title plotEql
#' 
#' @description 
#' Creates a \code{ggplot2} object that plots time series of biomass, harvest rate and catch. The basic object can then be modified by adding ggpot2 layers.
#'
#' @param  x an object of class \code{biodyn} 
#' @param  biomass an object of holding biomass at beginning of year 
#' @param ... other arguments
#'
#' @return an \code{ggplot2} object
#' 
#' @export
#' @rdname plotEql
#'
# @aliases plotEql,biodyn,FLBRP-method plotEql,biodyn,FLQuant-method plotEql,biodyn,missing-method
#' 
#' @examples
#' \dontrun{
#'  refpts('logistic',FLPar(msy=100,k=500))
#' }
setMethod('plotEql',signature(data='biodyn',biomass='missing'),  
          function(data,biomass,...) plotEqlfn(data,biomass=FLQuant(seq(0,max(params(data)['k']),length.out=101)),...))
setMethod('plotEql',signature(data='biodyn',biomass='FLQuant'),  
          function(data,biomass,...) plotEqlfn(data,biomass,...))
# setMethod('plotEql',signature(x='biodyn',biomass='FLBRP'),  
#           function(x,biomass,II=FALSE,...) plotProdfn(bd=x,brp=biomass,II=II,...))

plotEqlfn=function(data,biomass=FLQuant(seq(0,max(params(data)['k']),length.out=101)),...) {
  if ((dims(data)$iter>1 | dims(params(data))$iter>1) & dims(biomass)$iter==1) 
    biomass=propagate(biomass,max(dims(data)$iter,dims(params(data))$iter))
  
  res=model.frame(FLQuants(stock=biomass, 
                           yield=FLQuant(production(data,biomass))))
  res=transform(res,harvest=yield/stock)
  
  res=rbind(cbind(pnl="Equilibrium Stock v Harvest Rate",
                  transform(res,x=harvest,y=stock)[,-(7:9)]),
            cbind(pnl="Equilibrium Production v Harvest Rate",
                  transform(res,x=harvest,y=yield)[,-(7:9)]),
            cbind(pnl="Equilibrium Production v Stock",
                  transform(res,x=stock,  y=yield)[,-(7:9)]))
  
    p=ggplot(res)+geom_line(aes(x,y))+
      facet_wrap(~pnl,scales="free",ncol=1)+    
      xlab('') + ylab('')
  
  p} 


#' @title plotMSE
#'
#' @description 
#' Creates a \code{ggplot2} object that plots absolute and relative to MSY benchmarks time series of 
#' ssb, biomass, harvest rate and catch for  \code{FLStock} and  \code{biodyn}  objects
#' The basic object can then be modified by adding ggpot2 layers.
#'
#' @param  x, an object of class \code{biodyn} 
#' @param  y, an object of class \code{FLStock} 
#' @param  z, an object of class \code{FLBRP} 
#' @param ... other arguments
#'
#' @return an \code{ggplot2} object
#' 
#' @export
#' @rdname plotMSE
#'
#' @aliases plotMSE,biodyn,FLStock,FLBRP-method
#' 

# setMethod('plotMSE',signature(x='biodyn',y='FLStock',z='FLBRP'),  
#           function(x,y,z,...) plotMSEfn(x,y,z,...))


## compares age and biomass based time series      
plotMSEfn=function(mp,om,brp){
  ### OM ######################################
  ## absolute
  omAbs=cbind(Type='Absolute',
              rbind(as.data.frame(FLQuants(om,'stock','ssb','catch'),          drop=T),
                    as.data.frame(FLQuants(harvest=catch(om)/stock(om)),drop=T)))
  
  ## relative
  omRel=cbind(Type='Relative',as.data.frame(mcf(
    FLQuants('stock'  =stock(om)%/%refpts(brp)['msy','biomass'],
             'ssb'    =ssb(  om)%/%refpts(brp)['msy','biomass'],
             'harvest'=fbar( om)%/%refpts(brp)['msy','harvest'],
             #'harvest'=catch(om)/stock(om)%/%(refpts(brp)['msy','yield']/refpts(brp)['msy','biomass']),
             'catch'  =catch(om)%/%refpts(brp)['msy','yield'])),drop=T))
  
  ### MP ######################################
  ## absolute
  mpAbs=cbind(Type='Absolute',as.data.frame(FLQuants(mp,'stock','harvest','catch'),drop=T))
  
  ## relative
  mpRel=cbind(Type='Relative',as.data.frame(FLQuants('stock'  =stock(  mp)%/%bmsy(mp),
                                                     'harvest'=harvest(mp)%/%fmsy(mp),
                                                     'catch'  =catch(  mp)%/%msy( mp)),drop=T))
  
  ggplot(subset(rbind(cbind(What='OM',rbind(omAbs,omRel)),
                      cbind(What='MP',rbind(mpAbs,mpRel))),qname!='ssb'))+
    geom_line(aes(year,data,group=What,col=What))+
    facet_wrap(qname~Type,scales='free',ncol=2)  }

#' @title plotJack
#' 
#' @description 
#' Create a \code{ggplot2} plot based on a jack knifed biodyn and plots 
#' time series of biomass and harvest rate. 
#' The basic object can then be modified by adding ggpot2 layers.
#'
#' @param  x an object of class \code{biodyn} that has been jack knifed, i.e. by 
#' providing a jack knifed CPUE series to fit
#' @param y the original \code{biodyn} object
#' @param ncol number of colums in plot panel
#'
#' @return an \code{ggplot2} object
#' 
#' @export
#' @rdname plotJack
#'
#' @examples
#' \dontrun{
#' #simulate an object with known properties
#' bd=sim()
#' bd=window(bd,end=49)
#' 
#' #simulate a proxy for stock abundance
#' cpue=(stock(bd)[,-dims(bd)$year]+stock(bd)[,-1])/2
#' cpue=rlnorm(1,log(cpue),.2)
#' 
#' #set parameters
#' setParams(bd) =cpue
#' setControl(bd)=params(bd)
#' control(bd)[3:4,"phase"]=-1
#' 
#' #fit
#' bd=fit(bd,cpue)
#' 
#' bdJK=fit(bd,jackknife(cpue))
#' plotJack(bdJK)
#' bd  =randJack(100,bd)
#' }
plotJack=function(x,y,ncol=1){
   
  df=rbind(cbind(qname='stock',  
                 model.frame(jackSummary(stock(  x),stock(  y)),drop=TRUE)),
           cbind(qname='harvest',
                 model.frame(mcf(jackSummary(harvest(x),harvest(y))),drop=TRUE)))
  
  # basic plot data vs. year
  p=ggplot(data=df, aes(x=year, y=mean))+
    facet_wrap(~qname,ncol=ncol,scales='free_y')+
    geom_ribbon(aes(x=year, 
                  ymin=mean-2*se, 
                  ymax=mean+2*se),
                fill='blue', alpha = .20)+
    # line + xlab + ylab + limits to include 0 +
    geom_line(colour='red') + xlab('Year') + ylab('') + expand_limits(y=0) +
    #geom_line(aes(year,jack.mean+jack.bias),colour='black')   +
    # no legend
    theme(legend.title = element_blank())
  
  return(p)}

#' @title plotHcr
#'
#' @description 
#' Plots a hockey stick HCR with break pointts
#'
#' @param object an object of class \code{biodyn} or
#' @param params \code{FLPar} object with hockey stock HCR parameters
#' @param maxB  =1
#' @param rel   =TRUE
#' 
#' @return a \code{FLPar} object with value(s) for HCR
#' 
#' @rdname plotHcr
#' @aliases plotHcr-method  plotHcr,biodyn-method
#'
#' @export
#' @examples
#' \dontrun{
#' simBiodyn()
#' }
#' 
setMethod('plotHcr', signature(object='biodyn'),
          function(object,params=FLPar(ftar=0.7, btrig=0.7, fmin=0.01, blim=0.20),maxB=1,rel=TRUE){
            
            pts=rbind(cbind(refpt='Target',model.frame(rbind(bmsy(object)*c(params['btrig']),
                                                             fmsy(object)*c(params['ftar'])))),
                      cbind(refpt='Limit', model.frame(rbind(bmsy(object)*c(params['blim']),
                                                             fmsy(object)*c(params['fmin'])))))
            pts.=pts
            pts.[1,'bmsy']=params(object)['k']*maxB
            pts.[2,'bmsy']=0
            pts.[,1]=c('')
            
            
            pts=try(rbind(pts.[1,],pts[1:2,],pts.[2,]),silent=TRUE)
            if (is(pts)=="try-error")
              t.=as(rbind(as.data.frame(pts.[1,]),
                          as.data.frame(pts[1:2,]),
                          as.data.frame(pts.[2,])),"FLPar")
            
            names(pts)[2:3]=c('stock','harvest')
            
            if (rel){
              pts[,'stock']=pts[,'stock']/mpb::bmsy(object)
              pts[,'harvest']=pts[,'harvest']/mpb::fmsy(object)}
            
            pts})

# plotcf=function(x,y,qname="stock"){
#   ggplot(rbind.fill(cbind(What="MP",plot(FLQuants(llply(FLQuants(x,qname))))$data),
#                     cbind(What="OM",plot(FLQuants(llply(FLQuants(y,qname)))$data)))+
#     geom_ribbon(aes(year,min=`25%`,max=`75%`,fill=What),alpha=.5)+
#     geom_line(aes(year,`50%`,col=What))+
#     facet_wrap(~qname)
#   }

#' @title plotCcf
#' 
#' @description 
#'
#' @aliases plotCcf,biodyns,missing-method
#' 
#' @param x an object of class \code{biodyn} 
#' 
#' @return an \code{ggplot2} object
#' 
#' @export
#' @rdname plot
#' 
#' @examples
#' \dontrun{
#' 
#' } 
setMethod('plotCcf', signature(data='FLQuants'),
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

utils::globalVariables(c('geom_path','scale_size_manual','scale_linetype_manual'))
utils::globalVariables('data')
utils::globalVariables('data')

# plotComp.R - 
# ggplotFL/R/plotComp.R

# Copyright 2003-2007 FLR Team. Distributed under the GPL 2 or later
# Maintainer: Iago Mosqueira, JRC, Laurie Kell, ICCAT
# $Id:  $

# whooow {{{
whooow  =function(x,fn,probs)
  as.data.frame(FLQuants(lapply(fn,
                                function(fn,x)
                                  quantile(fn(x), probs=probs, na.rm=T), x=x))) # }}}

# plotComp {{{
plotComp = function(x, fn=NULL, probs=c(0.75,0.50,0.25), size=c(0.5,1.0,0.5),
                    lty=c(2,1,2), facet=facet_wrap(~qname, scales='free'),worm=NA) {
  
  if (dims(x)$iter>=length(probs)){  
    res = whooow(x,fn,probs)
    p1  = ggplot(res) + geom_path(aes(x=year,y=data,group=iter,size=iter,lty=iter)) +
      scale_size_manual(    values=size, name='Quantile') +
      scale_linetype_manual(values=lty , name='Quantile')
  }else{
    res = whooow(x,fn, 0.5)
    p1  = ggplot(res) + geom_path(aes(x=year,y=data))    }
  
  p1 = p1 +  expand_limits(y = 0) +
    xlab('Year') + ylab('') +
    facet
  
  if (length(worm) > 0)
    if (length(worm)<=dims(x)$iter)
      if (!(length(worm)==1 & is.na(worm[1])))  
        p1=p1+geom_path(aes(year,data,group=iter,colour=iter),
                        data=transform(subset(as.data.frame(FLQuants(lapply(fn,function(f,x) f(x), x=x))),iter %in% worm),iter=factor(iter)))
  
  p1
} # }}}

# plotComps {{{
plotComps = function(x, fn=NULL, probs=c(0.75,0.50,0.25), size=c(0.5,1.0,0.5),
                     lty=c(2,1,2), facet=facet_wrap(~qname,scales='free')) {
  
  x=x[seq(length(x))]
  
  if (max(laply(x,function(x) dims(x)$iter))>=length(probs))
    res = ldply(x, whooow, fn=fn, probs=probs)
  else
    res = ldply(x, whooow, fn=fn, probs=0.5)
  
  if ('X1' %in% names(res)) 
    names(res)[names(res)=='X1']='.id'
  
  if ('.id' %in% names(res)) 
    res$.id  = factor(res$.id)
  
  res$iter = factor(res$iter)
  
  if (length(unique(res$iter))>=length(probs)){
    p1  = ggplot(res) + geom_path(aes(x=year,y=data,group=.id:iter,size=iter,col=.id)) +
      scale_size_manual(    values=size, name='Quantile') +
      scale_linetype_manual(values=lty , name='Quantile')
  }else{
    p1  = ggplot(res) + geom_path(aes(x=year,y=data,group=.id,col=.id)) 
  }  
  
  p1= p1 + expand_limits(y = 0) +
    xlab('Year') + ylab('') +
    facet
  
  p1} 
# }}}

#' @title plotIndex
#' 
#' @description 
#'
#' @aliases plotIndex,biodyns,missing-method
#' 
#' @param x an object of class \code{biodyn} 
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
          function(data,
                   facet=facet_wrap(~qname,ncol=1,scales="free_y"),...){
            
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
          function(data,facet=facet_wrap(~qname,ncol=1,scales="free_y"),...){
            
            u=FLQuants(dlply(data,.(name), with,
                             as.FLQuant(data.frame(year=year,data=index))))
            plotIndex(u,facet,...)})

setMethod('plotIndex', signature(data='aspics'),
          function(data,facet=facet_wrap(~qname,ncol=1,scales="free_y"),...){
            
            u=ldply(data,index)
            u=u[!duplicated(u[,c("name","year")]),]
            u=FLQuants(dlply(u,.(name), with,
                             as.FLQuant(data.frame(year=year,data=index))))
            plotIndex(u,facet,...)})
