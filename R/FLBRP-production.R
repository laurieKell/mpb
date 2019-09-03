utils::globalVariables('biomass.obs')
utils::globalVariables('catch.obs')
utils::globalVariables('refpts<-')
utils::globalVariables('computeRefpts')
utils::globalVariables('ssb.obs')

#' @title production
#'
#' @description 
#' Estimates production for a given biomass
#' 
#' @param object either \emph{biodyn} or \emph{FLBRP} class
#' @param biomass an \emph{FLQuant} 
#' @param ... any other parameter
#' 
#' @rdname production
#' @export production
#' 
#' @aliases production,biodyn-method 
#'          production,FLStock,missing
#'          production,biodyn,FLQuant
#'  

productionFn=function(object,biomass="missing",
                      slots=c("landings.sel","discards.sel",
                              "stock.wt","landings.wt","discards.wt",
                              "m","mat",
                              "harvest.spwn","m.spwn")){
  
  nyr=dim(biomass.obs(object))[2]
  
  prd=biomass.obs(object)[,-nyr]+
    catch.obs(  object)[,-nyr]-
    biomass.obs(object)[,-1]
  
  #ref               =FLBRP::refpts(object)[1,]
  dimnames(ref)[[1]]="biomass"
  ref               =propagate(ref,nyr-1)
  ref[]             =NA
  ref[,"biomass"]   =biomass.obs(object)[,-nyr]
  refpts(object)    =ref
  
  if (!missing(biomass)){
    for (i in slots[!(slots%in%c("landings.sel","discards.sel"))]){
      slot(object,i)=propagate(slot(object,i),dims(ref)$iter)
      slot(object,i)[]=slot(biomass,i)[,-nyr]
    }
    
    if (any(c("landings.sel","discards.sel")%in%slots)){
      sel=harvest(biomass)[,-nyr]
      sel=sel%/%apply(sel[ac(range(biomass)["minfbar"]:range(biomass)["maxfbar"])],2,mean)
      
      if ("landings.sel"%in%slots){
        landings.sel(object)  =propagate(landings.sel(object),dims(ref)$iter)
        landings.sel(object)[]=sel*landings.n(biomass)[,-nyr]/catch.n(biomass)[,-nyr]}
      
      if ("discards.sel"%in%slots){
        discards.sel(object)  =propagate(discards.sel(object),dims(ref)$iter)
        discards.sel(object)[]=sel*discards.n(biomass)[,-nyr]/catch.n(biomass)[,-nyr]}
    }
  }
  
  refpts(object)=computeRefpts(object)
  
  res=cbind(as.data.frame(biomass.obs(object)[,-nyr],drop=T)[,c("year","data")],
            as.data.frame(ssb.obs(    object)[,-nyr],drop=T)[,c("data")],
            as.data.frame(catch.obs(  object)[,-nyr],drop=T)[,"data"],
            as.data.frame(rec.obs(    object)[,-nyr],drop=T)[,"data"],
            as.data.frame(biomass.obs(object)[,-nyr]-
                          ssb.obs(    object)[,-nyr],drop=T)[,"data"],
            as.data.frame(prd,drop=T)[,"data"]
            #model.frame(FLBRP::refpts(object)[,"yield"])[,"biomass"]
            )
  
  names(res)=c("year","biomass","ssb","catch","rec","juve","obs","hat")
  
  res[,c("year","biomass","ssb","juve","rec","catch","obs","hat")]}

# 
# setMethod('production', signature(object='FLBRP',biomass='missing'),
#           function(object,biomass) productionFn(object,biomass))
# 
# setMethod('production', signature(object='FLBRP',biomass='FLStock'),
#          function(object,biomass,slots=c("landings.sel","discards.sel",
#                                          "stock.wt","landings.wt","discards.wt",
#                                          "m","mat",
#                                          "harvest.spwn","m.spwn"))
#          productionFn(object,biomass,slots))
          
productionFn=function(object,biomass="missing",
                      slots=c("landings.sel","discards.sel",
                              "stock.wt","landings.wt","discards.wt",
                              "m","mat",
                              "harvest.spwn","m.spwn")){
  
  nyr=dim(biomass.obs(object))[2]
  
  prd=biomass.obs(object)[,-nyr]+
    catch.obs(  object)[,-nyr]-
    biomass.obs(object)[,-1]
  
  ref               =refpts(object)[1,]
  dimnames(ref)[[1]]="biomass"
  ref               =propagate(ref,nyr-1)
  ref[]             =NA
  ref[,"biomass"]   =biomass.obs(object)[,-nyr]
  refpts(object)    =ref
  
  if (!missing(biomass)){
    for (i in slots[!(slots%in%c("landings.sel","discards.sel"))]){
      slot(object,i)=propagate(slot(object,i),dims(ref)$iter)
      slot(object,i)[]=slot(biomass,i)[,-nyr]
    }
    
    if (any(c("landings.sel","discards.sel")%in%slots)){
      sel=harvest(biomass)[,-nyr]
      sel=sel%/%apply(sel[ac(range(biomass)["minfbar"]:range(biomass)["maxfbar"])],2,mean)
      
      if ("landings.sel"%in%slots){
        landings.sel(object)  =propagate(landings.sel(object),dims(ref)$iter)
        landings.sel(object)[]=sel*landings.n(biomass)[,-nyr]/catch.n(biomass)[,-nyr]}
      
      if ("discards.sel"%in%slots){
        discards.sel(object)  =propagate(discards.sel(object),dims(ref)$iter)
        discards.sel(object)[]=sel*discards.n(biomass)[,-nyr]/catch.n(biomass)[,-nyr]}
    }
  }
  
  refpts(object)=computeRefpts(object)
  
  res=cbind(as.data.frame(biomass.obs(object)[,-nyr],drop=T)[,c("year","data")],
            as.data.frame(ssb.obs(    object)[,-nyr],drop=T)[,c("data")],
            as.data.frame(catch.obs(  object)[,-nyr],drop=T)[,"data"],
            as.data.frame(rec.obs(    object)[,-nyr],drop=T)[,"data"],
            as.data.frame(biomass.obs(object)[,-nyr]-
                            ssb.obs(    object)[,-nyr],drop=T)[,"data"],
            as.data.frame(prd,drop=T)[,"data"],
            model.frame(refpts(object)[,"yield"])[,"biomass"])
  
  names(res)=c("year","biomass","ssb","catch","rec","juve","obs","hat")
  
  res[,c("year","biomass","ssb","juve","rec","catch","obs","hat")]}

# setMethod('production', signature(object='FLBRP',biomass='missing'),
#           function(object,biomass) productionFn(object,biomass))
# 
# setMethod('production', signature(object='FLBRP',biomass='FLStock'),
#           function(object,biomass,slots=c("landings.sel","discards.sel",
#                                       "stock.wt","landings.wt","discards.wt",
#                                       "m","mat",
#                                       "harvest.spwn","m.spwn"))
#             productionFn(object,biomass,slots))

setMethod('production', signature(object='FLStock',biomass='missing'),
          function(object,biomass)  {
            
            object=window(object,start=range(object)["minyear"],
                          end  =range(object)["maxyear"])
            
            res=stock(object)[,-dims(object)$year]-
              stock(object)[,-1]+catch(object)[,-dims(object)$year]
            
            units(res)=units(stock(object))
            res})
