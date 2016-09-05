productionFn=function(object,stk="missing",
                      slots=c("landings.sel","discards.sel",
                              "stock.wt","landings.wt","discards.wt",
                              "m","mat",
                              "harvest.spwn","m.spwn")){
  
  nyr=dim(biomass.obs(object))[2]
  
  prd=biomass.obs(object)[,-nyr]+
    catch.obs(  object)[,-nyr]-
    biomass.obs(object)[,-1]
  
  ref               =FLBRP::refpts(object)[1,]
  dimnames(ref)[[1]]="biomass"
  ref               =propagate(ref,nyr-1)
  ref[]             =NA
  ref[,"biomass"]   =biomass.obs(object)[,-nyr]
  refpts(object)    =ref
  
  if (!missing(stk)){
    for (i in slots[!(slots%in%c("landings.sel","discards.sel"))]){
      slot(object,i)=propagate(slot(object,i),dims(ref)$iter)
      slot(object,i)[]=slot(stk,i)[,-nyr]
    }
    
    if (any(c("landings.sel","discards.sel")%in%slots)){
      sel=harvest(stk)[,-nyr]
      sel=sel%/%apply(sel[ac(range(stk)["minfbar"]:range(stk)["maxfbar"])],2,mean)
      
      if ("landings.sel"%in%slots){
        landings.sel(object)  =propagate(landings.sel(object),dims(ref)$iter)
        landings.sel(object)[]=sel*landings.n(stk)[,-nyr]/catch.n(stk)[,-nyr]}
      
      if ("discards.sel"%in%slots){
        discards.sel(object)  =propagate(discards.sel(object),dims(ref)$iter)
        discards.sel(object)[]=sel*discards.n(stk)[,-nyr]/catch.n(stk)[,-nyr]}
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
            model.frame(FLBRP::refpts(object)[,"yield"])[,"biomass"])
  
  names(res)=c("year","biomass","ssb","catch","rec","juve","obs","hat")
  
  res[,c("year","biomass","ssb","juve","rec","catch","obs","hat")]}


setGeneric('production',   function(object,stk,...) standardGeneric('production'))

setMethod('production', signature(object='FLBRP',stk='missing'),
          function(object,stk) productionFn(object,stk))

setMethod('production', signature(object='FLBRP',stk='FLStock'),
         function(object,stk,slots=c("landings.sel","discards.sel",
                                     "stock.wt","landings.wt","discards.wt",
                                     "m","mat",
                                     "harvest.spwn","m.spwn"))
         productionFn(object,stk,slots))
          
productionFn=function(object,stk="missing",
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
  
  if (!missing(stk)){
    for (i in slots[!(slots%in%c("landings.sel","discards.sel"))]){
      slot(object,i)=propagate(slot(object,i),dims(ref)$iter)
      slot(object,i)[]=slot(stk,i)[,-nyr]
    }
    
    if (any(c("landings.sel","discards.sel")%in%slots)){
      sel=harvest(stk)[,-nyr]
      sel=sel%/%apply(sel[ac(range(stk)["minfbar"]:range(stk)["maxfbar"])],2,mean)
      
      if ("landings.sel"%in%slots){
        landings.sel(object)  =propagate(landings.sel(object),dims(ref)$iter)
        landings.sel(object)[]=sel*landings.n(stk)[,-nyr]/catch.n(stk)[,-nyr]}
      
      if ("discards.sel"%in%slots){
        discards.sel(object)  =propagate(discards.sel(object),dims(ref)$iter)
        discards.sel(object)[]=sel*discards.n(stk)[,-nyr]/catch.n(stk)[,-nyr]}
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

#' production
#'
#' @description Estimates estimates surplus production based on equilibrium analysis 
#' 
#' @param   object an object of class \code{FLBRP}
#' @param   stk an \code{FLStock}
#' @param   slots if provides then allows for non-stationarity by using the observed 
#' parameters each year, rather than the averages in  \code{FLBRP}
#'
#' @export
#' @rdname production
#'
#' @details Returns a data.frame with catch, biomass, observed (obs) and
#' expected (hat) production by year 
#' @examples
#' \dontrun{
#'    prd=production(x)
#' }

setGeneric('production',   function(object,stk,...) standardGeneric('production'))

setMethod('production', signature(object='FLBRP',stk='missing'),
          function(object,stk) productionFn(object,stk))

setMethod('production', signature(object='FLBRP',stk='FLStock'),
          function(object,stk,slots=c("landings.sel","discards.sel",
                                      "stock.wt","landings.wt","discards.wt",
                                      "m","mat",
                                      "harvest.spwn","m.spwn"))
            productionFn(object,stk,slots))

setMethod('production', signature(object='FLStock',stk='missing'),
          function(object,stk)  {
            
            object=window(object,start=range(object)["minyear"],
                          end  =range(object)["maxyear"])
            
            res=stock(object)[,-dims(object)$year]-
              stock(object)[,-1]+catch(object)[,-dims(object)$year]
            
            units(res)=units(stock(object))
            res})
