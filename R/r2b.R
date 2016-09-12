#' 
#' @description Biomass dynamic and age based model are equivalent if they provide similar advice on 
#' the response of a stock to exploitation and management, i.e. are 
#'   i) limit and target reference points  and 
#'   ii) the properties of the predicted time series
#' comparable.   
#'   
#' The main processes influencing the dynamics of exploited populations are gains due to growth and 
#' recruitment and losses due to fishing and natural mortality. In a biomass dynamic model recruitment, 
#' growth and natural mortality are simplified into a single production function, i.e. that of Pella 
#' and Tomlinson, which is parameterised by population growth rate, carrying capacity and its shape. The 
#' later is determined by density dependence.
#' 
#'
r2b=function(from){
  
  warn=options()$warn
  options(warn=-1)
  
  params=params(biodyn())
  refpts=refpts(biodyn())
  
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
