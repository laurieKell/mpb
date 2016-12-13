utils::globalVariables('data')

biodynCpue2aspic=function(object,value){
  res=aspic()
  
  ### New index #################################
  ## aspic
  # a) expunge old index
  # b) count number of indices
  # c) expand relevant slots
  # d) calculate default values
  
  idxs=names(value)
  
  ## a) expunge old index
  newIdx=ldply(value,as.data.frame,drop=TRUE)
  newIdx=transform(newIdx,name=.id,index=data)[c("name","year","index")]
  catch=window(object@catch,start=range(object)["minyear"],end=range(object)["maxyear"])
  catch=data.frame(as.data.frame(catch,drop=TRUE),name=unique(newIdx$name)[1])[,c("year","data","name")]
  names(catch)[2]="catch"
  dummy=data.frame(year=range(object)["minyear"]:range(object)["maxyear"])
  newIdx=merge(newIdx,dummy,by="year",all.y=TRUE)
  newIdx=merge(newIdx,catch,by=c("year","name"),all.x=TRUE)
  newIdx=newIdx[do.call("order",(newIdx[,c("name","year")])),]
  newIdx$code="CC"
  
  res@index=newIdx
  
  # stopmess  
  res@stopmess="not ran"
  
  # catch
  tmp=ddply(res@index, .(year), with, sum(catch,na.rm=TRUE))
  res@catch=as.FLQuant(tmp[,"V1"], dimnames=list(year=tmp[,"year"]))
  dmns=dimnames(res@catch)
  dmns$year=c(dmns$year,as.numeric(max(dmns$year))+1)
  
  # stock        
  res@stock=FLQuant(NA,dimnames=dmns)
  
  range(res)=range(as.numeric(dimnames(res@catch)$year))
  
  # diags       
  res@diags=data.frame(NULL)
  
  # params      
  res@params       =FLPar(as.numeric(NA),dimnames=list(params=c("b0","msy","k"),iter=1))
  res@params["msy"]=msyFn("pellat",object@params)
    
  # control 
  res@control=FLPar(as.numeric(NA),c(length(c(c("b0","msy","k"),paste("q",seq(length(idxs)),sep=""))),5),
                    dimnames=list(params=c(c("b0","msy","k"),paste("q",seq(length(idxs)),sep="")),
                                  c("fit","min","val","max","lambda"),iter=1))
  res@control[,    "fit"]=1
  res@control["b0","fit"]=0

  # vcov      
  res@vcov=FLPar(as.numeric(NA),dimnames=list(params=c("b0","msy","k"),param=c("b0","msy","k"),iter=1))
  
  # hessian       
  res@hessian=FLPar(as.numeric(NA),dimnames=list(params=c("b0","msy","k"),param=c("b0","msy","k"),iter=1))
  
  # objFn
  res@objFn  =FLPar(array(as.numeric(NA),dim=c(2,1),dimnames=list("value"=c("rss","rsq"),iter=1)))
  
  # mng 
  res@rnd=9785
  
  # mng 
  res@mng=FLPar()
  
  # mngVcov      
  res@mngVcov=FLPar()
  
  # profile         
  res@profile=data.frame(NULL)
  
  # desc        
  res@desc=paste(res@desc,"new index")
  
  nms=dimnames(res@params)$params[dimnames(res@params)$params %in% dimnames(object@params)$params]
  
  res@params[ nms]=object@params[ nms]
  #res@control[nms,c("min","val","max")]=object@control[nms,c("min","val","max")]
  #res@control[c("k","msy","b0"),"fit"]=c(1,1,0)
  
  res@stock=object@stock
  
  ## add q??s
  #setParams( res)=value
  
  nms=c("b0","msy","k")
  
  res@params=res@params[nms]
  idx<-index(res)
  stk<-as.data.frame(stock(res),drop=T)
  
  q=ddply(idx,.(name), function(idx){
    with(merge(idx,stk),
         mean(index/data,na.rm=T))})
  
  q.=FLPar(q[,2])
  dimnames(q.)$params=paste("q",seq(dim(q)[1]),sep="")
  res@params=rbind(params(res),q.)
  
  setControl(res)=params(res)
  
  res}

setMethod('aspic', signature(object="missing",value="missing"),
          function(object,value,...)
            new("aspic"))
setMethod('aspic', signature(object="biodyn",value="FLQuants"),
          function(object,value) biodynCpue2aspic(object,value))

setMethod('aspic', signature(object="biodyn",value="FLQuant"),
          function(object,value=FLQuants(value)) biodynCpue2aspic(object,value))

setMethod('aspic', signature(object="missing",value="missing"),
    function(object,value,...)
          {
            #args <- list(...)
            
            # if no FLQuant argument given, then use empty FLQuant
            #slots <- lapply(args, class)
            #slots <- names(slots)[slots == 'FLQuant']
            
          return(new("aspic"))
          })

setMethod('aspic', signature(object="data.frame",value="missing"),
    function(object,value,r=0.25,k=as.numeric(NA),msy=as.numeric(NA),...){
         
            args <- list(...)
            
            nms=names(object)
         
            res=new("aspic")
            
            if (all(c("year","catch") %in% nms)){
              o=ddply(object, .(year), with, sum(catch,na.rm=T))
              res@catch <- FLQuant(o$V1,dimnames=list(year=o$year))
              }
            
            ## CC
            if (all(c("year","catch") %in% nms) & !("index" %in% nms))
               object=transform(object, index=catch/effort)
            range(res)=unlist(list(minyear=min(object$year), maxyear=max(object$year)))
            
            res@index=object
            # Load given slots
            for(i in names(args))
              slot(res, i) <- args[[i]]
            
            nms=dimnames(res@params)
            nms$params=c(nms$params,paste("q",seq(length(unique(object$name))),sep=""))
            
            res@params=FLPar(as.numeric(NA),dimnames=nms)
            res@params["b0"] =1.0
            if (is.na(msy)) res@params["msy"]=mean(res@catch,na.rm=T)       else res@params["msy"]=msy
            if (is.na(k))   res@params["k"]  =mean(res@params["msy"])*4.0/r else res@params["k"]  =k
            
            res=fwd(res,catch=res@catch)
            
            setParams( res)=res@index
            setControl(res)=res@params
            
            res@control["b0","fit"]   =0
            res@control[1:3, "lambda"]=0
            
            res@rnd=2062012
            
            return(res)})

setMethod('readAspic', signature(object="character"),
          function(object,...)
            return(.readAspic(object)))
setMethod('aspic', signature(object="character",value="missing"),
    function(object,value,...)
          return(.readAspic(object)))


setMethod('aspic', signature(object="FLStock",value="missing"),
          function(object,value,min=0.5,max=1.5,...){
            
            res      =new("aspic")
            res@catch=catch(object)
            res@stock=window(catch(object),end=dims(catch(object))$maxyear+1)
            res@index =data.frame(model.frame(FLQuants(catch=catch(object),
                                                       index=catch(object)/fbar(object)/mean(catch(object)/fbar(object))),drop=T),type="CC",name="1")
            
            dmns=dimnames(res@params)
            dmns$params=c(dmns$params,"q1")
            
            res@params=FLPar(as.numeric(NA),dimnames=dmns)
            res@params[]=c(1,mean(res@catch),4*mean(res@catch), mean(res@index$index/res@index$catch)*.2)
            
            res@control[,"val"]=res@params
            res@control[,"min"]=res@control[,"val"]*min
            res@control[,"max"]=res@control[,"val"]*max
            res@control[,"lambda"]=1.0
            res@control[,"fit"]=c(0,1,1,1)
            
            range(res)[]=range(res@index$year)
            
            res@rnd=99999
            res})

#bdModel=attributes(model(new("FLBioDym")))$levels
# 
# setMethod('aspic', signature(object="data.frame"),
# asAspic=function(object,...){
#   
#             args <- list(...)
#             
#             res=new("aspic")
#             
#             ## The same
#             slot(res,"desc")    =slot(object,"desc")
#             slot(res,"name")    =slot(object,"name")
#             slot(res,"range")   =slot(object,"range")
#             
#             slot(res,"catch")   =slot(object,"catch")
#             slot(res,"stock")   =slot(object,"stock")
#             
#             slot(res,"stopmess")=slot(object,"stopmess")
#       
#             ## model
#             slot(res,"model")=switch(model(object),
#                                      schaefer=factor("LOGISTIC"),
#                                      fox     =factor("FOX"),
#                                      pellat  =factor("GENFIT"))
#           
#             slot(res,"params")  =paramFn(object)
#             slot(res,"control")  =controlFn(res)
#  
#             # Load given slots
#             for(i in names(args))
#               slot(res, i) <- args[[i]]
#                         
#             return(res)}


#' is.aspic
#'
#' @description Checks class type and returns TRUE if object is of type biodyn
#' @param x biodyn class
#' 
#' @return TRUE or FALSE
#' 
#' @export
#' @examples
#' \dontrun{
#'  is.aspic(aspic()) 
#'  }
is.aspic = function(x)
  return(inherits(x, 'aspic'))

