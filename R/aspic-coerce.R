utils::globalVariables(c("bdModel","biomass"))

setGeneric('as.biodyn', function(object, ...) standardGeneric('as.biodyn'))
setGeneric('as.aspic',  function(object, ...) standardGeneric('as.aspic'))

paramFn=function(object){
  
#   fox       =c("r","K")
#   schaefer  =c("r","K")
#   pellat    =c("r","K","p")
  
  b2aParams=function(model,params) {
    
    #model=model(bd)
    
    if(!(model %in% bdModel)) stop("has to be one of", bdModel)
    
    schaeferFn = function(biomass, params) { #logistic
      params["r"]=params["r"]*params["K"]/4
      dimnames(params)$params[1]="msy"
      
      params}
    
    foxFn =function(biomass, params) params
    
    pellatFn = function(biomass, params) params
    
    res = switch(model,
                 "fox"     =foxFn(     biomass,params),
                 "schaefer"=schaeferFn(biomass,params),
                 "pellat"  =pellatFn(  biomass,params),
                 stop("has to be either 'fox', 'schaefer' or 'pellat'"))
    
    res@rnd=9999
    return(res)}
   
   b2aParams(model(object),params(object))}
     
controlFn=function(object){
  
  schaefer2Logistic=function(x){FLPar()}
  fox2fox          =function(x){FLPar()}
  pellat2Genfit    =function(x){FLPar()}
  
  switch(model(object),
         pellat  =schaefer2Logistic(params(object)),
         fox     =fox2fox(          params(object)),
         schaefer=pellat2Genfit(    params(object)),
         gulland =stop("Gulland not available in ASPIC"),
         fletcher=stop("Gulland not available in ASPIC"),
         shepherd=stop("Gulland not available in ASPIC"))
  }


aspic2biodyn<-function(from,
                       phase=c("b0"=-1,"r"=4,"k"=3,"p"=-2,"q"=2,"sigma"=1),
                       min=0.5,max=1.5){
  sA=getSlots("aspic")
  sB=getSlots("biodyn")
  
  sA=sA[!(names(sA) %in% c("model","params"))]
  sB=sB[!(names(sB) %in% c("model","params"))]
  
  nits=dim(from@params)[2]
  par=FLPar(array(NA,dim=c(4,nits),dimnames=list(params=c("r","k","b0","p"),iter=seq(nits))))
  par["p"] =1
  par["b0"]=from@params["b0"]
  par["k"]=from@params["k"]
  par["r"]=c((from@params["msy"]%/%par["k"]))
  par["r"]=par["r"]%/%((1/(1+par["p"]))^(1/par["p"]+1))
  
  r<-function(msy,p,k)
    msy/k/(1/(1+p))^(1/p+1)
  
  res=biodyn()
  res@params=par
  
  nms=names(sA[(names(sA) %in% names(sB))])
  nms=nms[nms!="control"]
  for (i in nms)
    slot(res,i)=slot(from,i)
  
  cpue=index(from,F)
  if ("FLQuant"%in%is(cpue)) cpue=FLQuants("1"=cpue) 
  
  #bug
  setParams( res)        =cpue      
  setControl(res,min=min,max=max)=res@params
  
  for (i in names(phase[names(phase)%in%dimnames(control(res))[[1]]]))
    control(res)[i,"phase"]=phase[i]
  
  if ("q"%in%names(phase))
    control(res)[grep("q",dimnames(control(res))[[1]]),"phase"]=phase["q"]
  
  if ("sigma"%in%names(phase))
    control(res)[grep("s",dimnames(control(res))[[1]]),"phase"]=phase["sigma"]
  
  dimnames(res@objFn)$value=c("rss","ll")
  
  return(res)}

setMethod("as.biodyn", signature(object="aspic"),
          function(object, ...) {
            
            sA=getSlots("aspic")
            sB=getSlots("biodyn")
            
            sA=sA[!(names(sA) %in% c("model","params"))]
            sB=sB[!(names(sB) %in% c("model","params"))]
            
            par=FLPar("r"=.6,"k"=c(object@params["k"]),"b0"=c(object@params["b0"]),"p"=1)
            par["r"]=c(object@params["msy"]/(par["k"]*(1/(1+par["p"]))^(1/par["p"]+1)))
            
            nms=c(mpb:::modelParams("pellat"),'b0')
            par=par[nms]
            units(par)="NA"
            
            res=biodyn(factor("pellat"),par)
            
            res@control[c("p","b0"),1,1]=-1
            
            nms=names(sA)
            nms=nms[!(nms%in%c("control","priors"))]
            for (i in nms[nms %in% names(sB)])
              slot(res,i)=slot(object,i)
            
            # cpue=index(from,F)        
            # setParams( res)<-cpue      
            # setControl(res)            =params(res)
            
            dimnames(res@objFn)$value=c("rss","ll")
            
            return(res)})

setAs(from='aspic', to='biodyn',  def=function(from) as.biodyn(from))

biodyn2aspic=function(from){
        sA=getSlots("aspic")
        sB=getSlots("biodyn")
        
        sA=sA[!(names(sA) %in% c("model","params","hcr"))]
        sB=sB[!(names(sB) %in% c("model","params","hcr"))]
        
        res=aspic()
        
        model(res) =factor("LOGISTIC")
        params(res)=FLPar("msy"=msyFn("pellat",from@params),"k"=c(from@params["k"]),"b0"=c(from@params["b0"]))
        
        res@control[c("b0","k"),c("min","val","max")]=from@control[c("b0","k"),c("min","val","max")]
        res@control["msy",c("min","val","max")]=from@control["r",-1]/c(from@control["r",3])*c(params(res)["msy"])
        
        for (i in names(sA[(names(sA) %in% names(sB))]))
          slot(res,i)=slot(from,i)
        
        dimnames(res@objFn)$value=c("rss","ll")
        
        return(res)}

setMethod("as.aspic", signature(object="biodyn"),
    function(object, ...) {
            
        sA=getSlots("aspic")
        sB=getSlots("biodyn")
        
        sA=sA[!(names(sA) %in% c("model","params","hcr"))]
        sB=sB[!(names(sB) %in% c("model","params","hcr"))]
        
        res=aspic()
        
        model(res) =factor("LOGISTIC")
        params(res)=FLPar("msy"=msyFn("pellat",object@params),"k"=c(object@params["k"]),"b0"=c(object@params["b0"]))
        
        res@control[c("b0","k"),c("min","val","max")]=object@control[c("b0","k"),c("min","val","max")]
        res@control["msy",c("min","val","max")]=object@control["r",-1]/c(object@control["r",3])*c(params(res)["msy"])
        
        for (i in names(sA[(names(sA) %in% names(sB))]))
          slot(res,i)=slot(object,i)
        
        dimnames(res@objFn)$value=c("rss","ll")
        
        return(res)})

setAs(from='biodyn', to='aspic', def=function(from) as.aspic(from))

aspics2biodyns<-function(from,
                         phase=c("b0"=-1,"r"=4,"k"=3,"p"=-1,"q"=2,"sigma"=1),
                         min=0.5,max=1.5)
  biodyns(llply(from,aspic2biodyn,phase=phase,min=min,max=max))
  
setAs(from='aspics', to='biodyns', def=function(from) aspics2biodyns(from))