utils::globalVariables(c('admbCor','admbProfile',
                         'daply','residual', 'count'))
utils::globalVariables(c("%dopar%","foreach","i"))

utils::globalVariables(c('laply','llply','maply','mlply'))
utils::globalVariables('alply')
utils::globalVariables('coefficients')
utils::globalVariables('lm')
utils::globalVariables('X1')
utils::globalVariables('qqnorm')
utils::globalVariables('read.table')
utils::globalVariables('hat')
utils::globalVariables('read.table')
utils::globalVariables('MakeADFun')
utils::globalVariables('optimx')

#' @title fit
#'
#' @description 
#' A generic method for fitting catch and index of relative abundance for both \emph{biodyn} and \emph{aspic}.
#' 
#' @param object either \emph{biodyn} or \emph{aspic} class
#' @param index with relative abundance, \emph{FLQuant} or \emph{FLQuants}, if \strong{object} is of type \emph{bodyn}
#' @param ... any other parameter
#' 
#' @rdname fit
#' @export
#' 
#' @aliases fit,aspic-method fit,aspics-method fit,biodyn,FLQuant-method fit,biodyn,FLQuants-method  
#' @seealso  \code{\link{aspic}}, \code{\link{biodyn}}, \code{\link{jk}}, \code{\link{boot}}
setGeneric('fit',       function(object,index,...)  standardGeneric('fit'))
  setMethod('fit',signature(object='biodyn',index='FLQuant'),
          function(object,index=index, 
                   dir=tempfile("tmpdir"),
                   cmdOps=paste('-maxfn 500 -iprint 0'),
                   lav=FALSE,maxF=2.5,silent=TRUE){

            #sink(file = "/home/laurie/Desktop/temp/output.txt")          

            object@indices=FLQuants("1"=index)
            res=fitPella(object, 
                         dir=dir,cmdOps=cmdOps,lav=lav,maxF=maxF,silent=silent)
            #sink()
            
            res})

setMethod('fit',signature(object='biodyn',index='FLQuants'),
          function(object,index=index, 
                   dir=tempfile("tmpdir"),
                   cmdOps=paste('-maxfn 500 -iprint 0'),lav=FALSE,maxF=2.5,silent=!TRUE){
          
            object@indices=index
            fitPella(object, 
                     dir=dir,
                     cmdOps=cmdOps,lav=lav,maxF=maxF,silent=silent)})


setMethod('fit',signature(object='biodyn',index="missing"),
          function(object,index=index, 
                   dir=tempfile("tmpdir"),
                   cmdOps=paste('-maxfn 500 -iprint 0'),lav=FALSE,maxF=2.5,silent=!TRUE){
          
            fitPella(object, 
                     dir=dir,
                     cmdOps=cmdOps,lav=lav,maxF=maxF,silent=silent)})

setMethod('jackknife',signature(object='biodyn'),
         function(object,index="missing", 
                   dir=tempfile("tmpdir"),
                   cmdOps=paste('-maxfn 500 -iprint 0'),
                   lav=FALSE,silent=!TRUE){

            orig=object
            if (index!="missing"){
              if      (is(index)=="FLQuantJK") object@indices=FLQuants("1"=index)
              else if (is(index)=="FLQuant")   object@indices=FLQuants("1"=jackknife(index))
              else if (is(index)=="FLQuants") {
                if (any(laply(index,function(x) (is(x)=="FLQuantJK")))) object@indices=index
                else object@indices=jackknife(index)}
              else stop("index either has to be missing or of FLQuant, FLQuants, or FLQuantJK")
              }
            
            nits  =max(laply(object@indices,function(x) as.numeric(dims(x)$iter)))
            if (dim(control(object))[3]==1)
               control(object)=propagate(control(object),nits)
            if (dim(object@params)[2]==1)
              object@params=propagate(object@params,nits)
            if (dims(object@stock)["iter"]==1)
              object@stock   =propagate(stock(  object),nits)
            object=fitPella(object,dir=dir,cmdOps=cmdOps,lav=lav,maxF=maxF,silent=silent)
            
            object@stock=as(stock(object),"FLQuantJK")
            object@stock@orig=stock(orig)
            object@catch     =orig@catch
            
            for (iSlot in names(getSlots("biodyn")[getSlots("biodyn")=="FLPar"])){
              slot(object,iSlot)      =as(slot(object,iSlot),"FLParJK")
              slot(object,iSlot)@orig=slot(orig,iSlot)
              }

            attributes(object)['jk']=TRUE
            
            object})


#' @title fit
#' 
#' @description 
#' Estimates parameters \code{biodyn} class by fitting catch and CPUE indices
#' 
#' @param   object an object of class \code{biodyn}
#' @param   index an \code{FLQuant}, \code{FLQuants} or  \code{data.frame} object with CPUE indices
#' @param   ... other arguments
#'
#' @export
#' @rdname biodyn-fit
#' 
#' @aliases fit,FLPar,FLQuant-method fit,aspics,missing-method fit,biodyn,FLQuantJK-method
#'          fit,biodyn,FLQuantJKs-method
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
setExe=function(dir=tempfile("tmpdir")){
  ##### set up temp dir with exe for data files
  
  # Linux
  if (R.version$os=='linux-gnu') {
    exe = paste(system.file('bin', 'linux', package="mpb", mustWork=TRUE),"pella", sep='/')
    if (length(grep("-rwxrwxr-x",system(paste("ls -l","pella"),intern=TRUE)))==0)
      warning("Executable privilege not set for \n","pella",call.=FALSE)
    
    file.copy(exe, dir)
    dir = paste(dir, '/', sep='')
     
    # Windows
  } else if (.Platform$OS.type=='windows') {
    exe = paste(system.file('bin', 'windows', package="mpb", mustWork=TRUE), "pella.exe", sep='/')
    
    file.copy(exe, file.path(dir,"pella.exe"))


    #dir = paste(dir, '\\', sep='')
    # Mac OSX
  }else 
    stop()
  oldwd = getwd()
  
  # change wd to avoid exe case bug
  try(setwd(dir))
  
  oldwd}

setPella=function(x,dir=tempfile("tmpdir"),lav=FALSE) {
  
  if (file.exists(dir)) system(paste("rm", dir))
  dir.create(dir)
  
  if ((.Platform$OS=='windows'))
     setExe(dir)
  else if (.Platform$OS.type=="unix")
    setExePath()
  else if (.Platform$OS.type=="Mac OSX")
    stop("Mac not set up yet")
  else stop("what OS?")
  
  try(setwd(dir))
  
  # create input files ################################
  
  dgts=options()$digits
  options(digits=22)

  # cpue
  idx=as.data.frame(x@indices, drop=TRUE)
  names(idx)[-1:0+dim(idx)[2]]=c('index','name')
  idx=transform(idx,name=as.numeric(name))
  idx=idx[!is.na(idx$index),]

  nms=c(mpb:::modelParams('pellat'),'b0')
  if (length(unique(idx$name))>0)
    nmIdx=paste(c('q','sigma'), rep(unique(idx$name),each=2),sep='')
  else nmIdx=c('q','sigma')

  ctl = x@control[nms,]
  ctl[,2:4] = ctl[,c(2,4,3)]
  ctl=alply(array(c(ctl),dim=dim(ctl),dimnames=dimnames(ctl))[,,1,drop=T],1)
  names(ctl) = nms

  mpb:::writeADMB(ctl, paste(dir, '/', 'pella', '.ctl', sep=''),append=FALSE)
  mpb:::writeADMB(ifelse(lav,0,1), paste(dir, '/', 'pella', '.obj', sep=''),append=FALSE)

  cat('# q ####################\n', file=paste(dir, '/', 'pella', '.ctl', sep=''),append=TRUE)

  ctl           = x@control[nmIdx[grep('q',nmIdx)],,1]
  ctl[,2:4,1]   = ctl[,c(2,4,3),1]
  ctl           = alply(t(matrix(ctl,dim(ctl))),1)
  names(ctl)    = c('phase','lower','upper','guess')

  mpb:::writeADMB(ctl, paste(dir, '/', 'pella', '.ctl', sep=''),append=TRUE)

  cat('# sigma ################\n', file=paste(dir, '/', 'pella', '.ctl', sep=''),append=TRUE)
  ctl           = x@control[nmIdx[grep('s',nmIdx)],,1]
  ctl[,2:4,1]   = ctl[,c(2,4,3),1]
  ctl           = alply(t(matrix(ctl,dim(ctl))),1)
  names(ctl)    = c('phase','lower','upper','guess')

  mpb:::writeADMB(ctl, paste(dir, '/', 'pella', '.ctl', sep=''),append=TRUE)

  # prr file
  prr = x@priors[c(nms,c('msy','bmsy','fmsy'),nmIdx),] 
  prr = alply(prr,1)
  names(prr) = dimnames(x@priors)$params[1:9]
  mpb:::writeADMB(prr, paste(dir, '/', 'pella', '.prr', sep=''),append=FALSE)
 
  # ref file
  if (is.na(x@ref["yr"])) 
       x@ref["yr"]=as.integer(range(x)["maxyear"]-(range(x)["minyear"]+range(x)["maxyear"])/2)
  
  mpb:::writeADMB(x@ref, paste(dir, '/', 'pella', '.ref', sep=''),append=FALSE)

  # write data
  ctc = as.list(model.frame(FLQuants("catch"=x@catch), drop=TRUE))
  ctc = c(nYrs=length(ctc[[1]]), ctc)
  res = c(ctc, c(nIdxYrs=dim(idx)[1], nIdx=length(unique(idx$name)), idx))

  mpb:::writeADMB(res, paste(dir, '/', 'pella', '.dat', sep=''),append=FALSE)

  # vcov
  vcov(x)=FLPar(array(NA, dim     =c(dim(params(x))[1],dim(params(x))[1],1), 
                           dimnames=list(params=dimnames(params(x))[[1]],
                                         params=dimnames(params(x))[[1]],iter=1)))
  
  options(digits=dgts)
  return(x)}

getPella=function(obj) {
  ow=options("warn");options(warn=-1)
  
  t1 = read.table(paste('pella.rep',sep=''),skip =18,header=T) 

  # params
  t2 = unlist(c(read.table(paste('pella.rep',sep=''),nrows=4)))
  q. = unlist(c(read.table(paste('pella.rep',sep=''),nrows=1,skip=8)))
  s. = unlist(c(read.table(paste('pella.rep',sep=''),nrows=1,skip=10)))
  
  nms=c('r','k','b0','p')
  obj@params[nms,] = t2

  obj@params[grep('q',dimnames(obj@params)$params),]=q. 
  obj@params[grep('s',dimnames(obj@params)$params),]=s. 

  err=try(t3<-t(array(read.table('lls.txt',sep="\n"))),silent=TRUE)

  if (!(any(is(err)=="try-error"))){
    
    t3=t3[,length(t3)]
    t3=as.numeric(unlist(strsplit(str_trim(t3)," ")))

    obj@ll =FLPar(array(t(array(t3,c(length(s.),5))),
                dim=c(5,length(s.),1),
                dimnames=list(params=c("ll","rss","sigma","n","q"),
                              index =seq(length(s.)),
                              iter=1)))[-5,]
    }

  # stock biomass
  obj@stock[,1:dim(t1)[1]] = unlist(c(t1['stock'])) 

  #unlink(FLCore:::getFile(getwd()), recursive=TRUE, force=TRUE)
   
  system("rm pella.*")
  system("rm admodel.dep admodel.hes debug.txt fmin.log lls.txt mcmc_bio.csv mcmc_par.csv trace.txt")
  whereNow=getwd()
  setwd(FLCore:::getDir(whereNow))
  system(paste("rm -r",FLCore:::getFile(whereNow)))
  
  options(ow)
  
  return(obj)} 

#FLParBug sim@control['r','val',1]=c(.5,.6)
#exe(object, exeNm='biodyn', dir=tempfile("tmpdir"), set=set, get=biodyn::set, cmdOps=paste('-maxfn 500'))
  
activeParams=function(obj) dimnames(obj@control)$params[c(obj@control[,'phase']>-1)]
            
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

  jbn=jacobian(elasFn,log(c(object@params[parNms])),dmns=parNms,bd=bd,mn=mn,rg=rg)
  
  
  dimnms=list(params=parNms,
              stat  =c('stock',    'harvest',
                       'msy',      'bmsy',     'fmsy',
                       'catchMsy', 'stockBmsy','harvestFmsy',
                       'stockMn',  'harvestMn',
                       'stockRg',  'harvestRg'))
  
  jbn=FLPar(array(t(jbn),dim=dim(t(jbn)),dimnames=dimnms))
  units(jbn)='NA'
  
  return(jbn)}

setExePath=function(){
  
  sep =  function() if (R.version$os=='linux-gnu') ':' else if (.Platform$OS=='windows') ';' else ','
  
  # Linux
  if (R.version$os=='linux-gnu') {
    path = paste(system.file('bin', 'linux', package="mpb", mustWork=TRUE), sep='/')
    # Windows
  } else if (.Platform$OS.type == 'windows') {
    
    path = paste(system.file('bin', 'windows', package="mpb", mustWork=TRUE), sep='/')
    # Mac OSX
  }else 
    stop()
  
  #path='C:\\R\\R-2.15.1\\library\\aspic\\bin\\windows'
  
  if (length(grep("/mpb/bin/",Sys.getenv('PATH')))==0){

    Sys.setenv(PATH=paste(path, sep(), Sys.getenv('PATH'),sep=''))
    } 
  
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

runExe<-function(bd,indices=bd@indices,wkdir=tempfile(),cmdOps=paste('-maxfn 500 -iprint 0')){
    if (os.type("linux")) {
              
        wkdir=tempfile("tmpdir")
        system(paste("mkdir",wkdir))
        setwd(wkdir)
        
        mpb:::setPella(bd,wkdir)
        system(paste("cp",system.file(package="mpb","bin/linux/pella"),file.path(wkdir,"pella")))
        system(paste("./pella"))}
}


fitPella=function(object,
                  dir=tempdir(),
                  cmdOps=paste('-maxfn 500 -iprint 0'),lav=FALSE,maxF=2.5,silent=!TRUE){
  
  first=TRUE
  
  if (silent){
    ow=options("warn");options(warn=-1)
    
    tmp <- tempfile()
    sink(tmp)
    on.exit(sink())
    on.exit(file.remove(tmp),add=TRUE)
    }
  
  oldwd=getwd()
  
  its=max(laply(object@indices,function(x) dims(x)$iter),dims(catch(object))$iter)
  its=max(its,dims(object@control)$iter)

  catch(object)=propagate(catch(object),its)   
  
  max=min(dims(catch(object))$maxyear,max(laply(object@indices,function(x) dims(x)$maxyear)))
  if (!is.na(range(object)['maxyear'])) max=min(max,range(object)['maxyear']) 
  min=min(dims(catch(object))$minyear,max(laply(object@indices,function(x) dims(x)$minyear)))
  if (!is.na(range(object)['minyear'])) min=max(min,range(object)['minyear'])
  
  object=window(object,start=min,end=max)
  
  #index=FLQuants(llply(object@indices, window,start=min,end=max))
  
  slts=getSlots('biodyn')
  slts=slts[slts %in% c('FLPar','FLQuant')]
  
  its=max(maply(names(slts), function(x) {  dims(slot(object,x))$iter}))
  
  nms=dimnames(params(object))$params
  object@vcov   =FLPar(array(as.numeric(NA), dim=c(length(nms),length(nms),its), dimnames=list(params=nms,params=nms,iter=seq(its))))
  object@hessian=object@vcov
  
  object@ll=FLPar(NA,dimnames=list(params=c("ll","rss","sigma","n"),
                                   index =paste('u',seq(length(dimnames(params(object))$params[grep('q',dimnames(params(object))$params)])),sep=''),
                                   iter  =seq(1)))
  
  if (its>1){
    
    ## these are all results, so doesnt loose anything
    object@stock  =FLCore::iter(object@stock,   1)
    object@params =FLCore::iter(object@params,  1)
    object@objFn  =FLCore::iter(object@objFn,   1)
    object@vcov   =FLCore::iter(object@vcov,    1)
    object@ll     =FLCore::iter(object@ll,      1)
    object@hessian=FLCore::iter(object@hessian, 1)
    object@mng    =FLPar(a=1)
    object@mngVcov=FLPar(a=1,a=1)
  
    if (dim(object@stock)[6]==1) object@stock=propagate(object@stock, iter=its, fill.iter=TRUE)      
    if (dim(object@catch)[6]==1) object@stock=propagate(object@catch, iter=its, fill.iter=TRUE)     
    
    for(i in c("params","vcov","hessian","objFn")){
      slot(object, i)=FLCore::iter(slot(object, i),1)
      slot(object, i)<-propagate(slot(object, i), its)}
  }
  
  slot(object, "ll")<-propagate(slot(object, "ll"), its)
  
  
  for (i in seq(its)){       
    
    tmpObj         = FLCore::iter(object,i) 
    tmpObj@indices = FLQuants(llply(object@indices, FLCore::iter,i)) 
    
    for (s in names(slts)[-(7:8)])
      slot(tmpObj,s) = FLCore::iter(slot(object,s),i) 
    
    ##set data
    tmpObj=setPella(tmpObj,dir=dir,lav=lav)
    
    #bug in windows
    # try(if (length(grep("-rwxrwxr-x",system(paste("ls -l","pella"),intern=TRUE)))==0)
    #        warning("Executable privilege not set for \n","pella",call.=FALSE),silent=FALSE)
    
    #run
    if ((.Platform$OS=='windows'))
      system(paste(file.path(dir,'pella.exe'), ' ', cmdOps, sep='')) #,ignore.stdout=TRUE, ignore.stderr=TRUE)
    else if (.Platform$OS.type=="unix")
      system2('pella', args=cmdOps, stdout=NULL, stderr=NULL) #,ignore.stdout=TRUE, ignore.stderr=TRUE)
    else if (.Platform$OS.type=="Mac OSX")
      stop("Mac not set up yet")
    else stop("what OS?")
    
    #gets results
    tmpObj=getPella(tmpObj)  
    
    s=names(slts)[slts%in%c('FLQuant','FLPar')]
    
    for (s in s[!("hessian"%in%s)]){
      try(FLCore::iter(slot(object,s),i) <- slot(tmpObj,s)) 
    }
    
    if (its<=1 & file.exists(paste(dir,'admodel.hes',sep='/'))){
      ##hessian
      x<-file(paste(dir,'admodel.hes',sep='/'),'rb')
      nopar<-readBin(x,'integer',1)
      H<-matrix(readBin(x,'numeric',nopar*nopar),nopar)
      
      try(object@hessian@.Data[activeParams(tmpObj),activeParams(tmpObj),i] <- H, silent=TRUE)
      close(x)
      
      ## vcov
      #print(file.exists(paste(dir,'admodel.cov',sep='/')))
      
      if (file.exists(paste(dir,'admodel.cov',sep='/')))
        try(object@vcov@.Data[activeParams(tmpObj),activeParams(tmpObj),i] <- cv(paste(dir,'admodel.hes',sep='/')), silent=TRUE) 
      
      #if (file.exists(paste(dir,'admodel.cov',sep='/'))){
      #   x<-file(paste(dir,'admodel.cov',sep='/'),'rb')
      #   nopar<-readBin(x,'integer',1)
      #   H<-matrix(readBin(x,'numeric',nopar*nopar),nopar)
      #   try(object@vcov@.Data[activeParams(tmpObj),activeParams(tmpObj),i] <- H, silent=TRUE)
      #close(x)}
      
      if (file.exists(paste(dir,'pella.hst',sep='/')))
        object@hst=admbProfile(paste(dir,'pella.hst',sep='/'))$profile
      if (file.exists(paste(dir,'lpr.plt',sep='/')))
        object@profile=mdply(data.frame(var=c("r", "k","bnow","fnow",
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
    
    object@params@.Data[  ,i] = tmpObj@params
    #object@control@.Data[,,i] = tmpObj@control
    
    wrn=options()$warn
    options(warn=-1)
    err=try(object@objFn@.Data[,i]<-tmpObj@objFn,silent=TRUE)
    if (!(any(is(err)=="try-error"))) {
      #print(dim(object@objFn@.Data))
      #print(dim(tmpObj@objFn))
      #print("warning bug when objFn has iters")
    }
    options(warn=wrn)
    
    object@ll@.Data[,,i] = tmpObj@ll
    
    if (file.exists('pella.std')){
      err1=try(mng.<-read.table('pella.std',header=T)[,-1])
      
      #err2=try(mngVcov.<-fitFn(paste(dir,'pella',sep='/'))$vcov)
      
      ## FLPar hack
      if (first) {
        
        if (any(is(err1)!='try-error')) 
          object@mng=FLPar(array(unlist(c(mng.[   ,-1])), dim     =c(dim(mng.)[1],2,its),
                                 dimnames=list(param=mng.[,1],var=c('hat','sd'),iter=seq(its))))
        
        #if (any(is(err2)!='try-error')) 
        #   object@mngVcov<-FLPar(array(unlist(c(mngVcov.)),dim     =c(dim(mng.)[1],dim(mng.)[1],its),
        #                                               dimnames=list(param=dimnames(mngVcov.)[[1]],
        #                                                             param=dimnames(mngVcov.)[[1]],iter=seq(its))))
        
        first=!first  
      }else{    
        
        try(if (all(is(err1)!='try-error')) object@mng@.Data[    ,,i][]=unlist(c(mng.[,-1])))
        #try(if (all(is(err2)!='try-error')) object@mngVcov@.Data[,,i][]=unlist(c(mngVcov.)))
      }}
  }
  
  units(object@mng)='NA'  
  
  object=fwd(object,catch=catch(object)[,rev(dimnames(catch(object))$year)[1]],starvationRations=2) 
  
  if (length(grep('-mcmc',cmdOps))>0 & length(grep('-mcsave',cmdOps))>0){
    #'-mcmc 100000 -mcsave 100'
    setMCMC=function(obj,dir){
      ps=read.psv(paste(dir,'pella.psv',sep='/'))
      
      dmns=list(params=activeParams(obj),iter=seq(dim(ps)[1]))
      
      ps=array(t(ps),dim=unlist(llply(dmns,length)),dimnames=dmns)
      ps=FLPar(ps)
      
      units(ps)='NA'
      ps}
    
    par=setMCMC(object,dir)  
    
    cmd=strsplit(cmdOps,',')
    grp=unlist(gregexpr('-mcmc',cmd[[1]])) 
    mcmc =sub(' +', '', cmd[[1]][grp>0]) 
    mcmc =as.numeric(substr(mcmc,6,nchar(mcmc)))
    grp=unlist(gregexpr('-mcsave',cmd[[1]])) 
    mcsave=sub(' +', '', cmd[[1]][grp>0])
    mcsave=sub(' +', '', mcsave)
    mcsave=as.numeric(substr(mcsave,8,nchar(mcsave)))
    
    object@params=propagate(object@params[,1],dims(par)$iter)
    object@objFn =propagate(object@objFn[ ,1],dims(par)$iter)
    
    object@params[dims(par)$params,]=par
    object@stock=propagate(object@stock,dim(params(object))[2])
    object=fwd(object,catch=catch(object),starvationRations=2)  
    
    attributes(object@params)[['mcmc']]  =mcmc
    attributes(object@params)[['mcsave']]=mcsave 
  } 
  
  if ("FLQuant"%in%class(index)) index=FLQuants("1"=index)
  
  object=fwd(object,catch=catch(object)) 
  #stock(object)=fwd(params(object),catch=catch(object)) 
  
  # system("rm pella.*")
  # system("rm admodel.dep admodel.hes debug.txt fmin.log lls.txt mcmc_bio.csv mcmc_par.csv trace.txt")
  # whereNow=getwd()
  # setwd(FLCore:::getDir(whereNow))
  # system(paste("rm -r",FLCore:::getFile(whereNow)))
  # setwd(oldwd) 
  # 
  if (its<=1){
     object@diags=mdply(seq(length(object@indices)),function(i,index){
       stockHat=(stock(object)[,-dims(stock(object))$year]+stock(object)[,-1])/2
       hat     =stockHat*params(object)[paste("q",i,sep="")]
       
        
       res=model.frame(mcf(FLQuants(
         obs     =index[[i]],
         hat     =hat,
         residual=log(index[[i]]%/%hat[,dimnames(index[[i]])$year]))),drop=T)
       
       diagsFn(res)},index=object@indices)
     
     names(object@diags)[1]="name"
  }else 
     object@diags=dgsFn(object,object@indices)
  
  #  if (!is.null(catch)) catch(object)=catch
  
  object@stock=stock(mpb::fwd(object,catch=catch(object)))
  
  if (silent) options(ow)
  
  return(object)}

  
