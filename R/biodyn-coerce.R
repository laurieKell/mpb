globalVariables("ctrl")

#' as("biodyn", "FLStock")
#'
#' @description
#' Coerces an \code{FLStock} into a \code{biodyn} object.
#'
#' @name as
#' @family biodyn
#'
#' @rdname as


FLStock2biodyn=function(from,
                        control=biodyn()@control,
                        mult   =FLPar(rep(c(0.5,1,1/0.5),each=4),
                                      dimnames=dimnames(biodyn()@control[,-1]))){
  
  stock(from)=computeStock(from)
  bd         =biodyn()
  range(bd)[]=range(from)[c('minyear','maxyear')]
  bd@catch   =catch(from)
  
  eql     =FLBRP(from,ny=dims(from)$year)
  sr      =as.FLSR(from,model="bevholtSV")
  
  ow      =options("warn");options(warn=-1)
  sr      =fmle(sr,fixed=list(spr0=spr0(eql),s=0.95),control=list(silent=TRUE))
  options(ow)
  
  eql     =brp(FLBRP(from,ny=dims(from)$year,sr=ab(sr)))
  
  params(bd)["k" ] =refpts(eql)["virgin","biomass"]
  params(bd)["b0"] =mean(stock(from)[,1:5])/params(bd)["k" ]
  
  if (control["p","val"]==1){
    params(bd)["p"]=1
    params(bd)["r"]=4*refpts(eql)["msy","yield"]%/%params(bd)["k"]       
  }else{
    params(bd)["r" ] =log(lambda(leslie(eql,c(refpts(eql)["crash","harvest"]))))
    
    fn  =function(p,bmsy,k) (bmsy-k*(1/(1+p))^(1/p))^2
    
    params(bd)["p"]=optimise(fn, c(0.001,5),    
                             bmsy=c(refpts(eql)["msy","biomass"]),
                             k   =c(refpts(eql)["virgin","biomass"]))$minimum
  }
  
  nms =dimnames(params( bd))$params
  nms.=dimnames(mult)$params
  nms =nms.[nms.%in%nms]
  
  params(bd)[nms]=params(bd)[nms]*mult[nms,"val",drop=T]
  bd=mpb::fwd(bd,catch=catch(bd))
  setControl(bd)=params(bd)
  
  params(bd)[nms][!is.na(control(bd)[nms,"val"])]=
    control(bd)[nms,"val"][!is.na(control(bd)[nms,"val"])]
  
  #setParams(bd) =index
  setControl(bd)=params(bd)
  
  control(bd)[dimnames(control)$params,dimnames(control)[[2]], ][!is.na(control)]=control[!is.na(control)]
  
  nms=dimnames(mult)$params
  control(bd)[nms,"min"]=mult[nms,"min"]*control(bd)[nms,"val"]
  control(bd)[nms,"max"]=mult[nms,"max"]*control(bd)[nms,"val"]
  control(bd)[nms,"val"]=mult[nms,"val"]*control(bd)[nms,"val"]
  
  #bd=fit(bd,index)
  
  bd}

FLStock2biodynSimple=function(from){
  
  res      =biodyn()
  res@catch=catch(from)
  res@stock=window(stock(from),end=dims(from)$maxyear+1)
  
  dmns=dimnames(res@params)
  
  res@params=FLPar(as.numeric(NA),dimnames=dmns)
  res@params[]=c(.6,4*mean(res@catch),1,1)
  
  res@control[,'val']=res@params
  res@control[,'min']=res@control[,'val']*.1
  res@control[,'max']=res@control[,'val']*10
  res@control[,'phase']=c(1,1,-1,-1)
  
  range(res)[]=range(from)[c('minyear','maxyear')]
  
  res}

setAs('FLStock','biodyn',
      function(from) FLStock2biodyn(from))

FLBRP2biodyn=function(from,what=c("ssb","biomass","exploitable")[1],fix=c("bmsy","msy")[1]){
  
  warn=options()$warn
  options(warn=-1)
  
  r=lambda(leslie(from,f=c(refpts(from)["crash","harvest"]))[drop=T])-1
  msy =c(from@refpts["msy",   "yield"])
  
  fbar(from)[,1]=0
  k=switch(tolower(substr(what[1],1,1)),
      s=from@refpts["virgin","ssb"],
      b=from@refpts["virgin","biomass"],
      e=apply(stock.n(from)[,1]%*%stock.wt(from)[,1]%*%catch.sel(from)%/%
        apply(catch.sel(from),c(2,6),max),6,sum))

  b0=switch(tolower(substr(what[1],1,1)),
        s=ssb.obs(from)[,1]%/%k,
        b=from@stock.obs[,1]%/%k,
        e=from@stock.obs[,1]%/%from@refpts["virgin","biomass"])
  
  fbar(from)[,1]=from@refpts["msy","harvest"]  
  bmsy=switch(tolower(substr(what[1],1,1)),
           s=from@refpts["msy","ssb"],
           b=from@refpts["msy","biomass"],
           e=apply(stock.n(from)[,1]%*%stock.wt(from)[,1]%*%catch.sel(from)%/%
             apply(catch.sel(from),c(2,6),max),6,sum))
  
  msy=from@refpts["msy","yield"]
  
  if (tolower(substr(fix,1,1))=="m")
    p=optimise(function(p,msy,r,k) {
              res=(msy/(r*k)-(1/(1+p))^(1/p+1))^2
              res}, 
               c(1e-12,1e12),    
               msy=msy,
               r  =r,
               k  =k)$minimum
  else  
    p=mpb:::getP(bmsy,k,p=c(0.001,5))
  
  bd=biodyn(catch=catch.obs(from))
  
  bd@params=FLPar(c(r=r,k=k,p=p,b0=b0))
  
  bd@control["r", c("min","val","max")][]=c(params(bd)["r"])
  bd@control["k", c("min","val","max")]=params(bd)["k"]
  bd@control["p", c("min","val","max")]=params(bd)["p"]
  bd@control["b0",c("min","val","max")]=params(bd)["b0"]
  
  bd@control[,"min"]=bd@control[,"min"]*.1
  bd@control[,"max"]=bd@control[,"max"]*10
  
  bd@priors["r",   "a"]=params(bd)["r"]
  bd@priors["k",   "a"]=params(bd)["k"]
  bd@priors["p",   "a"]=params(bd)["p"]
  bd@priors["b0",  "a"]=params(bd)["b0"]
  bd@priors[ "msy","a"]=mpb:::refpts(bd)["msy"]
  bd@priors["bmsy","a"]=mpb:::refpts(bd)["bmsy"]
  bd@priors["fmsy","a"]=mpb:::refpts(bd)["fmsy"]
  
  bd=mpb:::fwd(bd,catch=catch(bd))
  
  range(bd)["minyear"]=dims(bd@catch)$minyear
  range(bd)["maxyear"]=dims(bd@catch)$maxyear
  
  options(warn=warn)
  
  bd}

setAs('FLBRP','biodyn', 
      function(from) FLBRP2biodyn(from,to))

setMethod('biodyn', signature(object='FLStock',params='missing'),
    function(object,params){
                  
      res      =new('biodyn')
      res@catch=catch(object)
      res@stock=window(stock(object),end=dims(object)$maxyear+1)
      
      dmns=dimnames(res@params)
      
      res@params=FLPar(as.numeric(NA),dimnames=dmns)
      res@params[]=c(.6,4*mean(res@catch),1,1)
      
      res@control[,'val']=res@params
      res@control[,'min']=res@control[,'val']*.1
      res@control[,'max']=res@control[,'val']*10
      res@control[,'phase']=c(1,1,-1,-1)
      
      range(res)[]=range(object)[c('minyear','maxyear')]
      
      res})

# setMethod('as.biodyn',signature(x='FLStock'),
#           as.biodyn<-function(x){FLStock2biodyn(x)})

p<-function(from){
  
  ## age based regference points
  msy =c(from@refpts["msy",   "yield"])
  bmsy=c(from@refpts["msy",   "biomass"])
  k   =c(from@refpts["virgin","biomass"])
  
  # biomass based reference points
  optimise(function(p,bmsy,k) 
    (bmsy-k*(1/(1+p))^(1/p))^2, c(0.001,5), bmsy=bmsy,k=k)$minimum}
  


