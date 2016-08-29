utils::globalVariables(c('admbCor','admbProfile',
                         'daply','residual', 'count'))
utils::globalVariables(c("%dopar%","foreach","i"))

utils::globalVariables(c('laply','llply','maply','mlply'))
utils::globalVariables('alply')

setMethod('fit',signature(object='biodyn',index='FLQuant'),
          function(object,index=index,exeNm='pella',package='mpb', 
                   dir=tempdir(),
                   cmdOps=paste('-maxfn 500 -iprint 0'),
                   lav=FALSE){

            #sink(file = "/home/laurie/Desktop/temp/output.txt")          
            
            res=fitPella(object,index=index,exeNm=exeNm,package=package, 
                         dir=dir,cmdOps=cmdOps,lav=lav)
            #sink()
            
            res})

setMethod('fit',signature(object='biodyn',index='FLQuants'),
          function(object,index=index,exeNm='pella',package='mpb', 
                   dir=tempdir(),
                   cmdOps=paste('-maxfn 500 -iprint 0'),lav=FALSE)
            fitPella(object,index,exeNm,package, 
                     dir=dir,
                     cmdOps=cmdOps,lav=lav))

setMethod('fit',signature(object='biodyn',index='FLQuantJKs'),
          function(object,index=index, 
                   dir=tempdir(),
                   cmdOps=paste('-maxfn 500 -iprint 0'),
                   lav=FALSE){
           
            nits  =max(laply(index,function(x) as.numeric(dims(x)$iter)))
            control(object)=propagate(control(object),nits)
            object@params=propagate(object@params,nits)
            object@stock   =propagate(stock(  object),nits)
            index =FLQuants(llply(index,as.FLQuant))
            object=fitPella(object,index=index, 
                            dir=dir,cmdOps=cmdOps,lav=lav)
            attributes(object)['jk']=TRUE
            object})

setMethod('fit',signature(object='biodyn',index='FLQuantJK'),
          function(object,index=index,exeNm='pella',package='mpb', 
                   dir=tempdir(),
                   cmdOps=paste('-maxfn 500 -iprint 0'),
                   lav=FALSE){
            
            object=propagate(object,dims(index)$iter)
            
            index =as.FLQuant(index)
            
            object=fitPella(object,index=index,exeNm=exeNm,package=package, 
                            dir=dir,cmdOps=cmdOps,lav=lav)
            
            attributes(object)['jk']=TRUE
            
            object})

#' fit
#'
#' Estimates parameters \code{biodyn} class by fitting catch and CPUE indices
#' 
#' @param   object an object of class \code{biodyn}
#' @param   index an \code{FLQuant}, \code{FLQuants} or  \code{data.frame} object with CPUE indices
#' @param   ... other arguments
#'
#' @export
#' @rdname fit-biodyn
#' 
#' @examples
#' \dontrun{
#'    #simulate an object with known properties
#'    bd=sim()
#'    bd=window(bd,end=49)
#'    
#'    #simulate a proxy for stock abundance
#'    cpue=(stock(bd)[,-dims(bd)$year]+stock(bd)[,-1])/2
#'    cpue=rlnorm(1,log(cpue),.2)
#'    
#'    #set parameters
#'    setParams(bd) =cpue
#'    setControl(bd)=params(bd)
#'    control(bd)[3:4,"phase"]=-1
#'    
#'    #fit
#'    bd=fit(bd,cpue)
#' }

diagsFn=function(res){
  ow=options("warn");options(warn=-1)
  
  res$residualLag <- c(res$residual[-1],NA)
  
  
  qqLine <- function(x,y){ 
    qtlx <- quantile(x, prob=c(0.25,0.75), na.rm=T)
    qtly <- quantile(y, prob=c(0.25,0.75), na.rm=T)
    
    a <- (qtly[1]- qtly[2]) / (qtlx[1] - qtlx[2])
    b <- qtly[1] - qtlx[1] * a
    
    res <- c(a,b)
    
    names(res) <- NULL
    names(res) <- c("a","b")
    
    return(res)}
  
  try({qq.     <- qqnorm(res$residual,plot.it=FALSE,na.rm=T)
       res$qqx <- qq.$x
       res$qqy <- qq.$y
       
       qqpar <- qqLine(qq.$x,qq.$y)[c("a","b")]
       
       res$qqHat=qqpar["a"]*res$qqx+qqpar["b"]})
  
  options(ow)
  
  res}

#### ADMB ###################################################################################
## 1) setExe copies exe from bin and creates a temo dir
## 2) Write data files
## 3) Run exe
## 4) Read output files
#############################################################################################

#tmp=model.frame(FLPar(array(NA,dim=c(3,2,1),dimnames=list(param=c('a','b','c'),var=c('x','y'),iter=1))))
#tmp=as.data.frame(FLPar(array(NA,dim=c(3,2,1),dimnames=list(param=c('a','b','c'),var=c('x','y'),iter=1))))

## copies exe into temp dir ready to run
setExe=function(exeNm,package,dir=tempdir()){
  ##### set up temp dir with exe for data files
  
  # Linux
  if (R.version$os=='linux-gnu') {
    exe = paste(system.file('bin', 'linux', package=package, mustWork=TRUE),exeNm, sep='/')
    if (length(grep("-rwxrwxr-x",system(paste("ls -l",exe),intern=TRUE)))==0)
      warning("Executable privilege not set for \n",exe,call.=FALSE)
    
    file.copy(exe, dir)
    dir = paste(dir, '/', sep='')
     
    # Windows
  } else if (.Platform$OS.type=='windows') {
    exe = paste(system.file('bin', 'windows', package=package, mustWork=TRUE), paste(exeNm, '.exe', sep=''), sep='/')
    file.copy(exe, dir)

    print(exe)
    print(dir)
    
    dir = paste(dir, '\\', sep='')
    
    # Mac OSX
  }else 
    stop()
  
  oldwd = getwd()
  
  # change wd to avoid exe case bug
  setwd(dir)
  
  oldwd}

setPella=function(obj, exeNm='pella', dir=tempdir(), lav=FALSE) {
  # create input files ################################
  dgts=options()$digits
  options(digits=22)

  # cpue
  if (is.FLQuant(obj[[2]])){  
     idx=model.frame(FLQuants(index=obj[[2]]), drop=TRUE)
     idx=data.frame(idx,name=1)
  }else if ('FLQuants' %in% class(obj[[2]])){  
     idx=as.data.frame(obj[[2]], drop=TRUE)
     names(idx)[2:3]=c('index','name')
     idx=transform(idx,name=as.numeric(name))
  }
  idx=idx[!is.na(idx$index),]

  bd.=obj[[1]]
 
  nms=c(modelParams('pellat'),'b0')
  if (length(unique(idx$name))>0)
    nmIdx=paste(c('q','sigma'), rep(unique(idx$name),each=2),sep='')
  else  
    nmIdx=c('q','sigma')

  ctl = bd.@control[nms,]
  ctl[,2:4] = ctl[,c(2,4,3)]
  ctl=alply(array(c(ctl),dim=dim(ctl),dimnames=dimnames(ctl))[,,1,drop=T],1)
#  ctl = alply(ctl,1)
  names(ctl) = nms
#   
#   # ctl file
#   ctl        =bd.@control[nms,,1]
#   ctl[,2:4,1]=ctl[,c(2,4,3),1]
# print("set 1")
#   ctl        =alply(ctl,1)
# print("set 2")
#   names(ctl) =nms

  writeADMB(ctl, paste(dir, '/', exeNm, '.ctl', sep=''),append=FALSE)
  writeADMB(ifelse(lav,0,1), paste(dir, '/', exeNm, '.obj', sep=''),append=FALSE)

  cat('# q ####################\n', file=paste(dir, '/', exeNm, '.ctl', sep=''),append=TRUE)

  ctl           = bd.@control[nmIdx[grep('q',nmIdx)],,1]
  ctl[,2:4,1]   = ctl[,c(2,4,3),1]
  ctl           = alply(t(matrix(ctl,dim(ctl))),1)
  names(ctl)    = c('phase','lower','upper','guess')

  writeADMB(ctl, paste(dir, '/', exeNm, '.ctl', sep=''),append=TRUE)
  
  cat('# sigma ################\n', file=paste(dir, '/', exeNm, '.ctl', sep=''),append=TRUE)
  ctl           = bd.@control[nmIdx[grep('s',nmIdx)],,1]
  ctl[,2:4,1]   = ctl[,c(2,4,3),1]
  ctl           = alply(t(matrix(ctl,dim(ctl))),1)
  names(ctl)    = c('phase','lower','upper','guess')

  writeADMB(ctl, paste(dir, '/', exeNm, '.ctl', sep=''),append=TRUE)

    # prr file
  prr = bd.@priors[c(nms,c('msy','bmsy','fmsy'),nmIdx),] 
  prr = alply(prr,1)
  names(prr) = dimnames(bd.@priors)$params[1:9]
  writeADMB(prr, paste(dir, '/', exeNm, '.prr', sep=''),append=FALSE)
 
  # ref file
  if (is.na(bd.@ref["yr"])) 
       bd.@ref["yr"]=as.integer(range(bd.)["maxyear"]-(range(bd.)["minyear"]+range(bd.)["maxyear"])/2)
  
  writeADMB(bd.@ref, paste(dir, '/', exeNm, '.ref', sep=''),append=FALSE)

  # write data
  ctc = as.list(model.frame(FLQuants("catch"=bd.@catch), drop=TRUE))
  ctc = c(nYrs=length(ctc[[1]]), ctc)
  res = c(ctc, c(nIdxYrs=dim(idx)[1], nIdx=length(unique(idx$name)), idx))

  writeADMB(res, paste(dir, '/', exeNm, '.dat', sep=''),append=FALSE)

#   # propagate as required
#   its = dims(bd)$iter
#   
#   # params
#   bd@params = propagate(params(bd)[,1], its)
#   # stock
#   stock(bd)  = FLQuant(dimnames=dimnames(stock(bd))[1:5], iter=its)
#   
  
  # vcov
  vcov(bd.)=FLPar(array(NA, dim     =c(dim(params(bd.))[1],dim(params(bd.))[1],1), 
                           dimnames=list(params=dimnames(params(bd.))[[1]],
                                         params=dimnames(params(bd.))[[1]],iter=1)))
  
  options(digits=dgts)
  return(bd.)}

getPella=function(obj, exeNm='pella') {
  ow=options("warn");options(warn=-1)
  
  t1 = read.table(paste(exeNm,'.rep',sep=''),skip =18,header=T) 

  # params
  t2 = unlist(c(read.table(paste(exeNm,'.rep',sep=''),nrows=4)))
  q. = unlist(c(read.table(paste(exeNm,'.rep',sep=''),nrows=1,skip=8)))
  s. = unlist(c(read.table(paste(exeNm,'.rep',sep=''),nrows=1,skip=10)))
  
  nms=c('r','k','b0','p')
  obj@params[nms,] = t2

  obj@params[grep('q',dimnames(obj@params)$params),]=q. 
  obj@params[grep('s',dimnames(obj@params)$params),]=s. 
    
  err=try(t3<-t(array(read.table('lls.txt',sep="\n"))))

  if (!(any(is(err)=="try-error"))){
    
    t3=t3[,length(t3)]
    t3=as.numeric(unlist(strsplit(str_trim(t3)," ")))
t3<<-t3
    obj@ll =FLPar(array(t(array(t3,c(length(s.),5))),
                dim=c(5,length(s.),1),
                dimnames=list(params=c("ll","rss","sigma","n","q"),
                              index =seq(length(s.)),
                              iter=1)))[-5,]
    }

  # stock biomass
  obj@stock[,1:dim(t1)[1]] = unlist(c(t1['stock'])) 

  options(ow)
  
  return(obj)} 

#FLParBug sim@control['r','val',1]=c(.5,.6)
#exe(object, exeNm='biodyn', dir=tempdir(), set=set, get=biodyn::set, cmdOps=paste('-maxfn 500'))
  
activeParams=function(obj) dimnames(obj@control)$params[c(obj@control[,'phase']>-1)]

fitPella=function(object,index=index,exeNm='pella',package='mpb', 
                  dir=tempdir(),cmdOps=paste('-maxfn 500 -iprint 0'),lav=FALSE,maxF=2.5)          
  {
  ow=options("warn");options(warn=-1)
  
  first=TRUE   
  catch=NULL
  if ('FLQuant'%in%is(index))
    index=FLQuants(index)
  
  its=max(laply(index,function(x) dims(x)$iter),dims(catch(object))$iter)
              
  catch(object)=propagate(catch(object),its)   

  max=min(dims(catch(object))$maxyear,max(laply(index,function(x) dims(x)$maxyear)))
  if (!is.na(range(object)['maxyear'])) max=min(max,range(object)['maxyear']) 
  min=min(dims(catch(object))$minyear,max(laply(index,function(x) dims(x)$minyear)))
  if (!is.na(range(object)['minyear'])) min=max(min,range(object)['minyear'])

  object=window(object,start=min,end=max)
  
  index=FLQuants(llply(index, window,start=min,end=max))
  
  its=max(its,dims(object@control)$iter)
  
  slts=getSlots('biodyn')
  slts=slts[slts %in% c('FLPar','FLQuant')]
 
  oldwd =setExe(exeNm,package,dir)
  oldwd=getwd()
  setwd(dir)
  #exe()
  
  object=list(object,index)
  bd =object[[1]]
  its=max(maply(names(slts), function(x) { 
          dims(slot(bd,x))$iter 
          }))
  
  nms=dimnames(params(bd))$params
  bd@vcov   =FLPar(array(as.numeric(NA), dim=c(length(nms),length(nms),its), dimnames=list(params=nms,params=nms,iter=seq(its))))
  bd@hessian=bd@vcov

  bd@ll=FLPar(NA,dimnames=list(params=c("ll","rss","sigma","n"),
                               index =paste('u',seq(length(dimnames(params(bd))$params[grep('q',dimnames(params(bd))$params)])),sep=''),
                               iter  =seq(1)))

  
  if (its>1){
   
      ## these are all results, so doesnt loose anything
      bd@stock  =FLCore::iter(bd@stock,  1)
      bd@params =FLCore::iter(bd@params, 1)
      bd@objFn  =FLCore::iter(bd@objFn,  1)
      bd@vcov   =FLCore::iter(bd@vcov,   1)
      bd@ll     =FLCore::iter(bd@ll,     1)
      bd@hessian=FLCore::iter(bd@hessian,1)
      bd@mng    =FLPar(a=1)
      bd@mngVcov=FLPar(a=1,a=1)
      
      if (dim(bd@stock)[6]==1) bd@stock=propagate(bd@stock, iter=its, fill.iter=TRUE)      
      if (dim(bd@catch)[6]==1) bd@stock=propagate(bd@catch, iter=its, fill.iter=TRUE)     

      for(i in c("params","control","vcov","hessian","objFn")){
        slot(bd, i)=FLCore::iter(slot(bd, i),1)
        slot(bd, i)<-propagate(slot(bd, i), its)}
        }

  slot(bd, "ll")<-propagate(slot(bd, "ll"), its)

  cpue=object[[2]]
  bd2 =object[[1]]
  
  for (i in seq(its)){       
     object[[2]] = FLCore::iter(cpue,i) 
     
     for (s in names(slts)[-(7:8)]){      
        slot(object[[1]],s) = FLCore::iter(slot(bd2,s),i) 
        }  

     object[[1]]=setPella(object,exeNm,dir,lav=lav)

     exe = paste(system.file('bin', 'linux', package="mpb", mustWork=TRUE),exeNm, sep='/')      

    #bug in windows
    try(
      if (length(grep("-rwxrwxr-x",system(paste("ls -l",exe),intern=TRUE)))==0)
        warning("Executable privilege not set for \n",exe,call.=FALSE) )

     # run
     system(paste('./', exeNm, ' ', cmdOps, sep=''),ignore.stdout=TRUE)
     #system(paste(exeNm, ' ', cmdOps, sep=''),ignore.stdout=TRUE)

     # gets results
     object[[1]]=getPella(object[[1]], exeNm)     
  
     s=names(slts)[slts%in%c('FLQuant','FLPar')]
    
     for (s in s[!("hessian"%in%s)]){

       try(FLCore::iter(slot(bd,s),i) <- slot(object[[1]],s)) 
       }

     if (its<=1 & file.exists(paste(dir,'admodel.hes',sep='/'))){
       ##hessian
       x<-file(paste(dir,'admodel.hes',sep='/'),'rb')
       nopar<-readBin(x,'integer',1)
       H<-matrix(readBin(x,'numeric',nopar*nopar),nopar)
       
       try(bd@hessian@.Data[activeParams(object[[1]]),activeParams(object[[1]]),i] <- H, silent=TRUE)
       close(x)
     
       ## vcov
       #print(file.exists(paste(dir,'admodel.cov',sep='/')))
       
       if (file.exists(paste(dir,'admodel.cov',sep='/')))
         try(bd@vcov@.Data[activeParams(object[[1]]),activeParams(object[[1]]),i] <- cv(paste(dir,'admodel.hes',sep='/')), silent=TRUE) 
        
       #if (file.exists(paste(dir,'admodel.cov',sep='/'))){
       #   x<-file(paste(dir,'admodel.cov',sep='/'),'rb')
       #   nopar<-readBin(x,'integer',1)
       #   H<-matrix(readBin(x,'numeric',nopar*nopar),nopar)
       #   try(bd@vcov@.Data[activeParams(object[[1]]),activeParams(object[[1]]),i] <- H, silent=TRUE)
       #close(x)}
       
       if (file.exists(paste(dir,'pella.hst',sep='/')))
          bd@hst=admbProfile(paste(dir,'pella.hst',sep='/'))$profile
       if (file.exists(paste(dir,'lpr.plt',sep='/')))
          bd@profile=mdply(data.frame(var=c("r", "k","bnow","fnow",
                                            "bnow","fnow","bnowthen","fnowthen",
                                            "msy","bmsy","fmsy","cmsy",
                                            "bmsy","ffmsy","bk","fr",
                                            "bratio","fratio","slopeb","slopef")),
                                function(var){
                                    #print(var)
                                    fl=paste("lp",var,".plt",sep="")
                                    if (file.exists(fl))
                                      admbPlt(fl)})
       }
 
     bd@params@.Data[  ,i] = object[[1]]@params
     bd@control@.Data[,,i] = object[[1]]@control

     bd@objFn@.Data[   ,i] = object[[1]]@objFn
     
     bd@ll@.Data[,,i] = object[[1]]@ll
     
     if (file.exists('pella.std')){
       err1=try(mng.<-read.table('pella.std',header=T)[,-1])
    
       #err2=try(mngVcov.<-fitFn(paste(dir,'pella',sep='/'))$vcov)
       
     ## FLPar hack
     if (first) {

       if (any(is(err1)!='try-error')) 
         bd@mng=FLPar(array(unlist(c(mng.[   ,-1])), dim     =c(dim(mng.)[1],2,its),
                                                     dimnames=list(param=mng.[,1],var=c('hat','sd'),iter=seq(its))))
        
       #if (any(is(err2)!='try-error')) 
       #   bd@mngVcov<-FLPar(array(unlist(c(mngVcov.)),dim     =c(dim(mng.)[1],dim(mng.)[1],its),
       #                                               dimnames=list(param=dimnames(mngVcov.)[[1]],
       #                                                             param=dimnames(mngVcov.)[[1]],iter=seq(its))))
 
       first=!first  
    }else{    
      
       try(if (all(is(err1)!='try-error')) bd@mng@.Data[    ,,i][]=unlist(c(mng.[,-1])))
       #try(if (all(is(err2)!='try-error')) bd@mngVcov@.Data[,,i][]=unlist(c(mngVcov.)))
       }}
  }

  units(bd@mng)='NA'  

  bd=mpb:::fwd(bd,catch=catch(bd)[,rev(dimnames(catch(bd))$year)[1]],starvationRations=2) 

  if (length(grep('-mcmc',cmdOps))>0 & length(grep('-mcsave',cmdOps))>0){
    #'-mcmc 100000 -mcsave 100'
     setMCMC=function(obj,dir){
       ps=read.psv(paste(dir,'pella.psv',sep='/'))

       dmns=list(params=activeParams(obj),iter=seq(dim(ps)[1]))

       ps=array(t(ps),dim=unlist(llply(dmns,length)),dimnames=dmns)
       ps=FLPar(ps)
       
       units(ps)='NA'
       ps}
     
    par=setMCMC(bd,dir)  
    
    cmd=strsplit(cmdOps,',')
    grp=unlist(gregexpr('-mcmc',cmd[[1]])) 
    mcmc =sub(' +', '', cmd[[1]][grp>0]) 
    mcmc =as.numeric(substr(mcmc,6,nchar(mcmc)))
    grp=unlist(gregexpr('-mcsave',cmd[[1]])) 
    mcsave=sub(' +', '', cmd[[1]][grp>0])
    mcsave=sub(' +', '', mcsave)
    mcsave=as.numeric(substr(mcsave,8,nchar(mcsave)))
    
    bd@params=propagate(bd@params[,1],dims(par)$iter)
    bd@objFn =propagate(bd@objFn[ ,1],dims(par)$iter)

    bd@params[dims(par)$params,]=par
    bd@stock=propagate(bd@stock,dim(params(bd))[2])
    bd=fwd(bd,catch=catch(bd),starvationRations=2)  

    attributes(bd@params)[['mcmc']]  =mcmc
    attributes(bd@params)[['mcsave']]=mcsave 
    } 

  if ("FLQuant"%in%class(index)) index=FLQuants("1"=index)

  if (its<=1){
      bd@diags=mdply(seq(length(index)),function(i,index){
        stockHat=(stock(bd)[,-dims(stock(bd))$year]+stock(bd)[,-1])/2
        hat     =stockHat*params(bd)[paste("q",i,sep="")]
        res=model.frame(mcf(FLQuants(
          obs   =index[[i]],
          #stock   =stock(bd),
          #stockHat=stockHat,
          hat     =hat,
          residual=log(index[[i]]/hat))),drop=T)
        
        diagsFn(res)},index=index)
  
      names(bd@diags)[1]="name"
  }else 
    bd@diags=dgsFn(bd,index)
    
  
  setwd(oldwd) 

#  if (!is.null(catch)) catch(object)=catch

  bd@stock=stock(mpb:::fwd(bd,catch=catch(bd)))
  
  options(ow)
  return(bd)}
            
#library(matrixcalc)
#is.positive.definite(H)

#vc=vcov(bd)[activeParams(bd),activeParams(bd),drop=T]
#vc[lower.triangle(vc)==0]=t(vc)[lower.triangle(vc)==0]
#is.positive.definite(vc,tol=1e-8)
#(vc-t(vc))/vc


admbCor=function(fl='pella.cor'){
  if (!file.exists(fl)) return
  
  ow=options("warn");  options(warn=-1)
  
  res=scan(fl,sep='\n',what=as.character(),quiet=TRUE)[-1]
  res=maply(res, str_trim)
  res=strsplit(res, ' +')
     
  nms=laply(res,function(x) x[2])[-1]
  hat=laply(res,function(x) as.numeric(x[3]))[-1]
  sd =laply(res,function(x) as.numeric(x[4]))[-1]
  
  npar=length(res)-1
  
  val =as.numeric(unlist(llply(res[-1], function(x) strsplit(x,' +')[-(1:4)])))
  cor=upper.triangle(array(1,dim=c(npar,npar),dimnames=list(nms,nms)))
  cor[upper.triangle(cor)==1]=val
  cor=cor+t(cor)

  diag(cor)=1
  vcov=cor*sd%o%sd
  
  options(ow)
  
  hat =FLPar(hat, units='NA')
  dimnames(hat)$params=dimnames(vcov)[[1]]
  vcov=FLPar(vcov,units='NA')
  names(vcov)[1:2]=c('params','params')
  
  
  return(FLPars(hat =hat,
                vcov=vcov))}

getDiags=function(fl='pella.rep'){
  skip=grep('Model',scan(fl,sep='\n',what=as.character(),quiet=TRUE))
      
  res=read.table(fl,skip=skip,header=TRUE)
  res=transform(res,residual=log(index/hat))
  res=diagsFn(res)
  
  res$harvest=res$catch/res$stock
  res$stock  =res$stock.
   
  res[,-7]}

calcElasticity=function(bd,mn=3,rg=5){
  
  elasFn=function(x,dmns,bd,mn,rg) {
    
    params(bd)[dmns]=exp(x)
    bd=fwd(bd,catch=catch(bd),starvationRations=2) 
    
    maxYr =ac(range(bd)['maxyear'])
    max1  =ac(range(bd)['maxyear']-(seq(mn)-1))
    max2  =ac(range(bd)['maxyear']-(seq(mn)-1+mn))
    maxR  =ac(range(bd)['maxyear']-(seq(rg)-1))
     
    smy=c(  stock(bd)[,maxYr],
          harvest(bd)[,maxYr],
              msy(bd),
             fmsy(bd),
             bmsy(bd),
            catch(bd)[,maxYr]/msy( bd),
            stock(bd)[,maxYr]/bmsy(bd),
            harvest(bd)[,maxYr]/fmsy(bd),
            mean(stock(bd)[,max1])/mean(stock(bd)[,max2]),
            mean(harvest(bd)[,max1])/mean(harvest(bd)[,max2]),
            abs(coefficients(lm(data~year,as.data.frame(stock(  bd)[,maxR],drop=T)))['year']),
            abs(coefficients(lm(data~year,as.data.frame(harvest(bd)[,maxR],drop=T)))['year']))     
      
    return(log(smy))}
  
  parNms=c(modelParams(model(bd)),'b0')

  jbn=jacobian(elasFn,log(c(bd@params[parNms])),dmns=parNms,bd=bd,mn=mn,rg=rg)
  
  
  dimnms=list(params=parNms,
              stat  =c('stock',    'harvest',
                       'msy',      'bmsy',     'fmsy',
                       'catchMsy', 'stockBmsy','harvestFmsy',
                       'stockMn',  'harvestMn',
                       'stockRg',  'harvestRg'))
  
  jbn=FLPar(array(t(jbn),dim=dim(t(jbn)),dimnames=dimnms))
  units(jbn)='NA'
  
  return(jbn)}

exe=function(package='mpb'){

  sep =  function() if (R.version$os=='linux-gnu') ':' else if (.Platform$OS=='windows') ';' else ','
  
  # Linux
  if (R.version$os=='linux-gnu') {
    path = paste(system.file('bin', 'linux',  package=package, mustWork=TRUE), sep='/')
    # Windows
  } else if (.Platform$OS.type == 'windows') {
    
    path = paste(system.file('bin', 'windows', package=package, mustWork=TRUE), sep='/')
    # Mac OSX
  }else 
    stop()
  
  #path='C:\\R\\R-2.15.1\\library\\aspic\\bin\\windows'
  
  path <- paste(path, sep(), Sys.getenv('PATH'),sep='')
   
  Sys.setenv(PATH=path)

  return(path)}

#asp=aspic('http://gbyp-sam.googlecode.com//svn//trunk//tests//aspic//swon//2009//high//aspic.inp')
#asp=fit(asp)

calcSS=function(x) daply(x@diags, .(name), 
                         with, sum(residual^2,na.rm=T)/sum(count(!is.na(residual))))

#FLParBug
#tst=FLPar(as.numeric(NA),dimnames=list(params=c('a','b'),col=c('x','y'),iter=1:2))
#tst['a','x',2]=3

fitFn=function(file){

  res=admbFit(file)

  #est        
  hat =FLPar(array(c(res$est),dim     =c(length(res$names),1),
                              dimnames=list(params=res$names,iter=1)))
  units(hat)='NA'  
    
  #std        
  std =FLPar(array(c(res$std),dim     =c(length(res$names),1),
                              dimnames=list(params=res$names,iter=1)))
  units(std)='NA'  
    
  #cor 
  cor =FLPar(array(c(res$cor),dim     =c(length(res$names),length(res$names),1),
                              dimnames=list(params=res$names,params=res$names,iter=1)))
  units(cor)='NA'  
    
  #cov
  vcov=FLPar(array(c(res$cov),dim     =c(length(res$names),length(res$names),1),
                              dimnames=list(params=res$names,params=res$names,iter=1)))
    
  units(vcov)='NA'  
    
  return(list(std       =std,hat=hat,cor=cor,vcov=vcov,  
              nlogl     =res$nlogl,      
              maxgrad   =res$maxgrad,    
              npar      =res$npar,       
              logDetHess=res$logDetHess))}

#file='/tmp/Rtmp6dBbkc/pella'
#fitFn(file)

#as.data.frame(t(array(read.csv("/tmp/RtmptQ93hN/lls.txt",sep=""))))

dgsFn=function(bd,index){
  
  stockHat=(stock(bd)[,-dims(stock(bd))$year]+stock(bd)[,-1])/2
  
  itB=dims(stockHat)$iter
  diags=mdply(seq(length(index)),function(i,index){
    
    itI=dims(index[[i]])$iter
    mdply(data.frame(iter=seq(max(itB,itI))), function(iter) {
      
      obs=iter(index[[i]],min(iter,itB))
      hat=iter(stockHat,min(iter,itB))*iter(params(bd),min(iter,itB))[paste("q",i,sep="")]
      yrs=dimnames(obs)$year[dimnames(obs)$year%in%dimnames(hat)$year]
      
    res=as.data.frame(log(obs[,yrs]/hat[,yrs]),drop=T)})},index=index)
  
  res=transform(diags[!is.na(diags$data),],name=names(index)[X1])[,c("name","year","iter","data")]
  names(res)[4]="residual"
  res}


logit<-function(min,x,max){
  x=(x-min)/(max-min)
  log(x/(1-x))}

invLogit<-function(min,x,max){
  x=1/(1+exp(-x))
  return(x*(max-min)+min)}

setControlFn<-function(params,min=0.01,max=100){
  dmns=dimnames(params)[1]
  dmns[["option"]]=c("phase","min", "val","max")
  dmns[["iter"]]  =seq(dim(params)[2])
  
  control=FLPar(array(c(1,1,-1,-1,rep(NA,12)),dim=c(4,4,dim(params)[2]),dimnames=dmns))
  control[1:4,"val"]=c(params)
  control[,"min"]=control[,"val"]*min
  control[,"max"]=control[,"val"]*max 
  
  control}

setLogitFn<-function(control){
  
  control=aaply(control,3,function(x){
    flg=x[,"phase"]>0
    if (any(flg))
      x[flg,"val"]=logit(x[flg,"min"],x[flg,"val"],x[flg,"max"])
    x})
  
  FLPar(aperm(control,c(2:3,1)))}

#control=setControlFn(params[,1],min=0.0001,max=10000)
#control["r","phase"]=-1

setMethod('fit',signature(object='FLPar',index='FLQuant'),
  function(object,index=index,...){
            
  args=list(...)
  catch=args[["catch"]]           
  nits=max(dims(object)$iter,dims(catch)$iter,dims(index)$iter)
  
  #DATA_VECTOR(ctc);      nyr
  #DATA_MATRIX(idx);      nyr, nidx
  
  nyr =dim(catch)[2]
  if (nyr!=dim(index)[2]) stop("catch and index dont match")
  
  #DATA_MATRIX(ctl);      4,   3
  if (dim(object)[1]!=4) stop("control needs 4 parameters")
  if (dim(object)[2]!=4) stop("control needs 4 options")
  
  warn=options()$warn
  options(warn=-1)
  sink("/dev/null")
              
  res=mdply(seq(nits),function(i) {
      print(i)
        
      params=FLPar(iter(object,i)[,"val",drop=TRUE])
    
      ctc=c(window(iter(catch,i)))
      idx=matrix(iter(index,i),dim(index)[2],1)
               
      #cat(i,file="/home/laurie/Desktop/tmp/chk2.txt")
                
      flg=(iter(object[,"phase"],i)>0)[drop=T]
        
      if (any(flg))
        val=c(logit(iter(object[flg,"min"],i),iter(object[flg,"val"],i),iter(object[flg,"max"],i)))
                
          
      f=MakeADFun(data      =list(ctc=ctc,idx=idx,ctl=iter(object,i)[drop=TRUE]),
                  parameters=list(par=val),
                  DLL       ="pellaTmb")
      hat=try(optimx(f$par,f$fn,f$gr,control=list(trace=0,fnscale=-1),method="BFGS"))
                
      if ("character" %in% mode(hat))
        hat=try(optimx(f$par,f$fn,control=list(trace=0,fnscale=-1),method="BFGS"))
      if ("character" %in% mode(hat))
        return(NULL)
                
      if (sum(flg)>0){
        params[flg]=invLogit(iter(object[flg,2],i),
                             unlist(c(hat[seq(sum(flg))])),
                             iter(object[flg,4],i))
          }
                
      nms=c("value","fevals","gevals","niter","convcode","kkt1","kkt2","xtimes")
      par=unlist(c(model.frame(params)[-5],hat[nms]))
      names(par)[length(flg)+1]="ll"
      par})
            
  sink();sink(); 
  options(warn=warn)
  
  
  res})
