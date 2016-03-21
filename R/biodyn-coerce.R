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
FLStock2biodyn=function(from){
  
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

FLBRP2biodyn=function(from){
  
  ## age based regference points
  msy =c(from@refpts["msy",   "yield"])
  bmsy=c(from@refpts["msy",   "biomass"])
  k   =c(from@refpts["virgin","biomass"])
  
  # biomass based reference points
  p   =optimise(function(p,bmsy,k) 
    (bmsy-k*(1/(1+p))^(1/p))^2, c(0.001,5), bmsy=bmsy,k=k)$minimum
  k=bmsy/((1/(1+p))^(1/p))
  r=msy/(k*(1/(1+p))^(1/p+1))
  b0=mean(stock(from)[,1:5]/k)
  
  bd=biodyn()
  
  bd@params=FLPar(c(r=r,k=k,p=p,b0=b0))
  control(bd)["r", c("min","val","max")]=params(bd)["r"]
  control(bd)["k", c("min","val","max")]=params(bd)["k"]
  control(bd)["p", c("min","val","max")]=params(bd)["p"]
  control(bd)["b0",c("min","val","max")]=params(bd)["b0"]
  
  control(bd)[,"min"]=control(bd)[,"min"]*.1
  control(bd)[,"max"]=control(bd)[,"max"]*10
  
  bd@priors["r",   "a"]=params(bd)["r"]
  bd@priors["k",   "a"]=params(bd)["k"]
  bd@priors["p",   "a"]=params(bd)["p"]
  bd@priors["b0",  "a"]=params(bd)["b0"]
  bd@priors[ "msy","a"]=refpts(bd)["msy"]
  bd@priors["bmsy","a"]=refpts(bd)["bmsy"]
  bd@priors["fmsy","a"]=refpts(bd)["fmsy"]
  
  bd@catch=catch.obs(from)
  bd@stock=from@stock.obs
  
  range(bd)["minyear"]=dims(bd@catch)$minyear
  range(bd)["maxyear"]=dims(bd@catch)$maxyear
  
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
  
       
