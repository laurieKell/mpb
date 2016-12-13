utils::globalVariables(c("bdModel","biomass"))

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

setAs('biodyn','aspic',
      function(from){
        sA=getSlots("aspic")
        sB=getSlots("biodyn")
        
        sA=sA[!(names(sA) %in% c("model","params","hcr"))]
        sB=sB[!(names(sB) %in% c("model","params","hcr"))]
        
        res=aspic()
        
        model(res) =factor("LOGISTIC")
        params(res)=FLPar("msy"=msyFn("pellat",params(from)),"k"=c(params(from)["k"]),"b0"=c(params(from)["b0"]))
        
        res@control[c("b0","k"),c("min","val","max")]=from@control[c("b0","k"),c("min","val","max")]
        res@control["msy",c("min","val","max")]=from@control["r",-1]/c(from@control["r",3])*c(params(res)["msy"])
        
        for (i in names(sA[(names(sA) %in% names(sB))]))
          slot(res,i)=slot(from,i)
        
        dimnames(res@objFn)$value=c("rss","ll")
        
        return(res)})

biodyn2aspic=function(from){
        sA=getSlots("aspic")
        sB=getSlots("biodyn")
        
        sA=sA[!(names(sA) %in% c("model","params","hcr"))]
        sB=sB[!(names(sB) %in% c("model","params","hcr"))]
        
        res=aspic()
        
        model(res) =factor("LOGISTIC")
        params(res)=FLPar("msy"=msyFn("pellat",params(from)),"k"=c(params(from)["k"]),"b0"=c(params(from)["b0"]))
        
        res@control[c("b0","k"),c("min","val","max")]=from@control[c("b0","k"),c("min","val","max")]
        res@control["msy",c("min","val","max")]=from@control["r",-1]/c(from@control["r",3])*c(params(res)["msy"])
        
        for (i in names(sA[(names(sA) %in% names(sB))]))
          slot(res,i)=slot(from,i)
        
        dimnames(res@objFn)$value=c("rss","ll")
        
        return(res)}

aspic2biodyn<-function(from,
                       phase=c("b0"=-1,"r"=4,"k"=3,"p"=-2,"q"=2,"sigma"=1),
                       min=0.5,max=1.5){
  sA=getSlots("aspic")
  sB=getSlots("biodyn")
  
  sA=sA[!(names(sA) %in% c("model","params"))]
  sB=sB[!(names(sB) %in% c("model","params"))]
  
  nits=dim(params(from))[2]
  par=FLPar(array(NA,dim=c(4,nits),dimnames=list(params=c("r","k","b0","p"),iter=seq(nits))))
  par["p"] =1
  par["b0"]=params(from)["b0"]
  par["k"]=params(from)["k"]
  par["r"]=c((params(from)["msy"]%/%par["k"]))
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
  biodyn:::setParams( res)        =cpue      
  setControl(res,min=min,max=max)=res@params
  
  for (i in names(phase[names(phase)%in%dimnames(control(res))[[1]]]))
      control(res)[i,"phase"]=phase[i]
  
  if ("q"%in%names(phase))
     control(res)[grep("q",dimnames(control(res))[[1]]),"phase"]=phase["q"]
  
  if ("sigma"%in%names(phase))
    control(res)[grep("s",dimnames(control(res))[[1]]),"phase"]=phase["sigma"]
  
  dimnames(res@objFn)$value=c("rss","ll")
  
  return(res)}

setAs('aspic', 'biodyn', function(from) 
  aspic2biodyn(from,
               phase=c("b0"=-1,"r"=4,"k"=3,"p"=-2,"q"=2,"sigma"=1),
               min=0.5,max=1.5))

#   dat <- edat(harvest = albsa$catch_tonnes, 
#               index   = cbind(fleet1 = albsa$fleet1_cpue,
#                             fleet2 = albsa$fleet2_cpue, fleet3 = albsa$fleet3_cpue, fleet4 = albsa$fleet4_cpue,
#                             fleet8 = albsa$fleet8_cpue), time = rownames(albsa))
#   time
#   n
#   sigmao
#   sigmap
#   renormalise=TRUE
  
setAs('aspic', 'biodyn',
      function(from){
        
        sA=getSlots("aspic")
        sB=getSlots("biodyn")
        
        sA=sA[!(names(sA) %in% c("model","params"))]
        sB=sB[!(names(sB) %in% c("model","params"))]
        
        par=FLPar("r"=.6,"k"=c(params(from)["k"]),"b0"=c(params(from)["b0"]),"p"=1)
        par["r"]=c(params(from)["msy"]/(par["k"]*(1/(1+par["p"]))^(1/par["p"]+1)))
        res=biodyn(factor("pellat"),par)
        
        res@control[c("p","b0"),1,1]=-1
        
        for (i in names(sA[(names(sA) %in% names(sB))]))
          slot(res,i)=slot(from,i)
        
        cpue=index(from,F)        
        setParams( res)            =cpue      
        setControl(res)            =params(res)
        
        dimnames(res@objFn)$value=c("rss","ll")
        
        return(res)})

aspics2biodyns<-function(from,
                         phase=c("b0"=-1,"r"=4,"k"=3,"p"=-1,"q"=2,"sigma"=1),
                         min=0.5,max=1.5)
  biodyns(llply(from,aspic2biodyn,phase=phase,min=min,max=max))
  

setAs('aspic', 'biodyn', function(from) aspics2biodyns(from))