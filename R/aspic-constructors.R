setMethod('aspic', signature(object="missing"),
    function(object,...)
          {
            #args <- list(...)
            
            # if no FLQuant argument given, then use empty FLQuant
            #slots <- lapply(args, class)
            #slots <- names(slots)[slots == 'FLQuant']
            
          return(new("aspic"))
          })

setMethod('aspic', signature(object="data.frame"),
    function(object,r=0.25,k=as.numeric(NA),msy=as.numeric(NA),...){
         
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
setMethod('aspic', signature(object="character"),
    function(object,...)
          return(.readAspic(object)))
