utils::globalVariables(c("year","swon","year","B","obs"))
utils::globalVariables(c("m_ply","b"))
utils::globalVariables(c("%dopar%","foreach","i"))

setGeneric('fit',       function(object,index,...)  standardGeneric('fit'), package='mpb')

#' @title fit
#'
#' @description
#' Fits the aspic model to catch and catch per unit effort data
#'
#' @rdname aspic-fit
#' @export
#' @docType methods
#'
#' @param object; an \code{aspic} object
#' @param dir; an optional \code{dir} where aspic text files used for fitting can be found
#' 
#' @return An aspic object with fitted values and parameter estimates
#' @seealso \code{\link{aspic},\link{biodyn},\link{boot},\link{jk}}
#'
#' @aliases fit,aspic,missing-method

#' @title boot, 
#' 
#' @description 
#' Bootstraps the ASPIC biomass dynamic model.
#'
#' @name boot
#' @export
#' @docType methods
#'
#' @description
#' Bootstraps the aspic model
#'
#' @param object; an \code{aspic} object or a character string giving an aspic "inp" file
#' @param dir; an optional \code{dir} where aspic text files used for fitting can be found
#' @return An aspic object with fitted values and parameter estimates
#' @seealso \code{\link{biodyn},\link{boot},\link{jk}}
#'
#' @examples
#' \dontrun{
#'     data(asp)
#'     asp=boot(asp)}


#' @title jk
#'
#' @description
#' Fits the aspic model to catch and catch per unit effort data removing 1 cpue observation at a time
#'
#' @name jk
#' @rdname jk
#' @export
#' @docType methods
#'
#' @param object; an \code{aspic} object or
#' @param object; a character string giving an aspic "inp" file
#' @param dir; an optional \code{dir} where aspic text files used for fitting can be found
#' @return An aspic object with fitted values and parameter estimates
#' @seealso \code{\link{biodyn},\link{boot},\link{fit}}
#'
#' @examples
#' \dontrun{
#'     data(asp)
#'     asp=jk(asp)}


chkIters=function(object){
  N=max(dims(object)$iter,dims(object@control)$iter,dims(object@params)$iter)
  
  #params(object)=FLPar(object@control[,"val",drop=T])
  
  if (N>1){
    object@stock=propagate(stock( object),N)
    if (dims(object@params)[1]!=N)
      object@params=propagate(object@params,N)
    
    object@ll   =propagate(object@ll,   N)
    object@objFn=propagate(object@objFn,N)
  }
  
  object@stock=propagate(stock(object),dims(object)$iter)
  object@stock=window(stock(object),end=max(as.numeric(dimnames(catch(object))$year))+1)
  
  return(object)}

jkIdx=function(x) dimnames(x)[[1]][ !is.na(x$index)]


runExe=function(object,package="mpb",exeNm="aspic",dir=tempdir(),jk=FALSE){
  ow=options("warn");options(warn=-1)
  
  object@index=object@index[object@index$year %in% range(object)["minyear"]:range(object)["maxyear"],]
  
  if (any(is.na(object@catch))){
    tmp=ddply(object@index, .(year), with, sum(catch,na.rm=TRUE))
    object@catch=as.FLQuant(tmp[,"V1"], dimnames=list(year=tmp[,"year"]))
    dmns=dimnames(object@catch)
    dmns$year=c(dmns$year,as.numeric(max(dmns$year))+1)
    object@stock=FLQuant(NA,dimnames=dmns)
  }
  
  oldwd=getwd()
  setwd(dir)
  #path=exe(package)
  
  if (.Platform$OS.type == "windows")
    file.copy(paste(paste(system.file("bin", "windows", package="mpb", mustWork=TRUE),exeNm, sep="/"),"exe",sep="."), dir)
  else  
    file.copy(     paste(system.file("bin", "linux", package="mpb", mustWork=TRUE),exeNm, sep="/"),                 dir)
  
  ## Jack knife if wished
  j=1
  if (jk){
    object=propagate(object,length(jkIdx(object@index)))
    object@params=propagate(object@params,length(jkIdx(object@index)))
    j   = jkIdx(object@index)
    index=object@index}
  
  us=paste("u",seq(length(dimnames(params(object))$params[grep("q",dimnames(params(object))$params)])),sep="")
  dmns=list(params=us,quantity=c("ll","ss","n"),iter=seq(1))
  object@ll=FLPar(array(NA,dim=unlist(lapply(dmns,length)),dimnames=dmns))
  
  object=chkIters(object)
  for (i in seq(dims(object)$iter)){
    m_ply(c("prn","rdat","bio","inp","fit","sum","rdatb","det","sum","bot"), function(x)
      if (file.exists(paste(exeNm,".",x,sep=""))) system(paste("rm ",exeNm,".",x,sep="")))
    
    if (jk){
      object@index=index
      object@index[j[i],"index"]=NA
    }
  
    # create exe input files
    .writeAspicInp(FLCore::iter(object,i),what="FIT",niter=1,
                   fl=paste(exeNm,".inp",sep=""))
     
    # run
    #exe=file.path(system.file("bin", "linux", package="mpb", mustWork=TRUE),"aspic")
    #inp=file.path(getwd(),"aspic.inp")
    #system(paste(exeNm, paste(" ",exeNm,".inp",sep=""),sep=""))
    
    system(paste("aspic.exe aspic.inp"),ignore.stdout=TRUE)
    
    rdat=dget(file.path(getwd(),paste(exeNm,"rdat",sep=".")))
    
    #rdat$estimates
    object@params[c("b0","msy","k"),i]=rdat$estimates[c("B1.K","MSY","K")]
    object@params[4:dim(object@params)[1],i]=rdat$estimates[substr(names(rdat$estimates),1,2)=="q."]
    
    names(rdat$t.series)=tolower(names(rdat$t.series))
    
    FLCore::iter(object@stock,i)=as.FLQuant(transform(rdat$t.series[,c("year","b")],data=b)[c("year","data")])[,dimnames(object@stock)$year]
  
    if (.Platform$OS!="windows"){
      try(object@objFn[1,i]<-rdat$diagnostics$obj.fn.value)
      
      #try(object@objFn[2,i]<-rdat$diagnostics$rsquare)
      rtn=aspicPrn(paste(exeNm,"prn",sep="."))
      
      pos=seq(dim(params(object))[1])[substr(dimnames(params(object))[[1]],1,1)=="q"]
      object@diags<-rtn
      #           object@diags=transform(object@diags,
      #                 stock.  =hat/c(object@params[pos,i])[name],
      #                 stockHat=obs/c(object@params[pos,i])[name])
      
      object@diags$name=factor(unique(object@index$name)[as.integer(object@diags$name)])
      object@diags=merge(object@diags,model.frame(mcf(FLQuants(stock  =object@stock)),drop=TRUE),all=T)
      object@diags$stock=object@diags$stock.
      
      object@diags=object@diags[,-10]
      object@diags=object@diags[!is.na(object@diags$name),]
      
    } else {
      rtn=try(readAspic(paste(exeNm,"prn",sep=".")))
      
      if (is.data.frame(rtn)) object@diags=rtn[!is.na(rtn$residual),]
      
      #           object@diags=transform(object@diags,stock.  =hat/c(object@params[grep("q",dimnames(params(object))$params),i])[name],
      #                                               stockHat=obs/c(object@params[grep("q",dimnames(params(object))$params),i])[name])
      object@diags$name=ac(factor(object@diags$name, labels=unique(object@index$name)))
      
      object@diags=merge(object@diags,model.frame(mcf(FLQuants(stock=object@stock,harvest=harvest(object))),drop=TRUE),all=T)
      #object@diags$stock=object@diags$stock.
      object@diags=object@diags[,-10]
      try(object@objFn[1,i]<-sum(diags(object)$residual^2,na.rm=T))
      
      object@diags=object@diags[!is.na(object@diags$name),]
    }
    
    #         dgs=subset(object@diags,!is.na(object@diags$residual))
    #         try(object@ll@.Data[,"ll",i]<-daply(dgs, .(name), with,
    #                                             calcLogLik(as.FLQuant(data.frame(data=residual,year=year)),
    #                                                        as.FLQuant(data.frame(data=1,year=year)),type=3)))
    #
    #         try(object@ll@.Data[,"ss",i]<-daply(dgs, .(name), with, sum(residual^2)))
    #         try(object@ll@.Data[,"n", i]<-daply(dgs, .(name), function(x) dim(x)[1]))
    #print(i)
    #print(daply(object@diags,.(name),with,sum(residual^2,na.rm=TRUE)))
    object@ll[,"ss",i]=daply(object@diags,.(name),with,sum(residual^2,na.rm=TRUE))
    
    #print(object@ll[,"ss",i])
  }
  
  #if (dims(object)$iter!=1) object@diags=data.frame(NULL)
  setwd(oldwd)
  
  options(ow)
  
  return(object)}

runBoot=function(object, package="mpb", exeNm="aspic", dir=tempdir(),boot=500){
  
  print(system.file('bin', 'linux', package=package, mustWork=TRUE))
  
  object@index=object@index[object@index$year %in% range(object)["minyear"]:range(object)["maxyear"],]
  
  if (boot<3) {
    boot=3
    warning("Requires a minimum of 3 bootstraps so boot option changed")
  }
  
  ## add catch baed on index catches if missing
  if (any(is.na(object@catch))){
    tmp=ddply(object@index, .(year), with, mean(catch,na.rm=TRUE))
    object@catch=as.FLQuant(tmp[,"V1"], dimnames=list(year=tmp[,"year"]))
  }
 
  ## add catch baed on index catches if missing
  if (dim(object@control)[3] >1)    stop("control can only have iter dim of 1")
  if (dim(object@params)[2]==1)    object@params=propagate(object@params,boot)
  if (dim(object@params)[2]!=boot) stop("params iters either have to be 1 or same as number of boot")
  if (dims(object@catch)$iter>1)   stop("catch must only have 1 iter")
   
  dmns=dimnames(object@catch)
  dmns$year=c(dmns$year,as.numeric(max(dmns$year))+1)
  object@stock=propagate(FLQuant(NA,dimnames=dmns),boot)
  
  oldwd=setExe(exeNm,package,dir)
  m_ply(c("prn","rdat","bio","inp","fit","sum","rdatb","det","sum","bot"), function(x)
    if (file.exists(paste(exeNm,".",x,sep=""))) system(paste("rm ",exeNm,".",x,sep="")))
  
  # create exe input files
  .writeAspicInp(object,what="BOT",niter=boot,fl=paste(exeNm,".inp",sep=""))

  # run
  system(paste("./", exeNm, paste(" ",exeNm,".inp",sep=""),sep=""))
  
  det=aspicDet(paste(exeNm,"det",sep="."))
  bio=aspicBio(paste(exeNm,"bio",sep="."))
  
  coerceDP=function(x)  FLPar(unlist(c(t(x))),params=names(x),iter=dim(x)[1])
  
  parNms=c("b0",modelParams(tolower(model(object))))
  
  object@params[parNms,]=coerceDP(det[,parNms])
  qNms=names(det)[!(names(det) %in% c("iter","stock","harvest","r","trial","loss","msy","bmsy","brel","frel","b1.k",parNms))]
  
  object@params[-(seq(length(parNms)))][]=unlist(c(det[,qNms]))
  object@stock=bio[["stock"]]%*%bio[["bmsy"]]
  
  setwd(oldwd)
  
  return(object)}

setMethod('fit',signature(object='aspic',index="missing"),
          function(object,dir=tempdir(), package="mpb", exeNm="aspic",jk=FALSE)
            runExe(object=object, dir=dir, package="mpb", exeNm="aspic",jk=jk))

setMethod('fit',signature(object='aspics',index="missing"),
          function(object, dir=tempdir(), package="mpb", exeNm="aspic",jk=FALSE,
                   .combine=list,
                   .multicombine=TRUE,.maxcombine=length(object),
                   .packages=c("mpb","plyr","reshape")){

            if (is.null(.combine)) ..combine=list else ..combine=.combine

            res=foreach(i=seq(length(object)), .combine=..combine,
                        .multicombine=.multicombine,
                        .maxcombine  =.maxcombine,
                        .packages    =.packages) %dopar% {

                        fit(object[[i]],dir=tempdir())
                        }

            #             if (is.null(.combine)) {
#               res=aspics(list(unlist(res)))
#               #names(res)=names(object)
#             }

            names(res)=names(object)

      return(aspics(res))})

setMethod('boot', signature(object='aspic'),
          function(object, dir=tempdir(), package="mpb", exeNm="aspic",boot=500)
            runBoot(object=object, dir=dir,package="mpb", exeNm="aspic",boot=boot))
setMethod('boot',  signature(object='aspics'),
          function(object, dir=tempdir(), package="mpb", exeNm="aspic",boot=500,
                   .combine=NULL,
                   .multicombine=T,.maxcombine=10,.packages=c("mpb","plyr")){

            if (is.null(.combine)) ..combine=list else ..combine=.combine
            res=foreach(i=names(object), .combine=..combine,
                        .multicombine=.multicombine,
                        .maxcombine  =.maxcombine,
                        .packages    =.packages) %dopar% {
                          wkdir=tempfile('file', dir)
                          dir.create(wkdir, showWarnings = FALSE)
                          boot(object[[i]],dir=wkdir,boot=boot)}

            if (is.null(.combine)) {
              res=aspics(res)
              names(res)=names(object)}

            res})

setMethod('jk',  signature(object='aspic'),
          function(object, dir=tempdir(), package="mpb", exeNm="aspic")
            runExe(object=object, dir=dir, package="mpb", exeNm="aspic",jk=TRUE))


