inline double pella(double r, double k, double p=1, double biomass=0) {
    return(biomass*r/p*(1-pow(biomass/k,p)));}

FLQuant fwd(FLQuant F, FLQuant C, FLQuant B, FLQuant params) {
     
     if ('FLQuant' %in% class(stock)   |
         'FLQuant' %in% class(harvest) |
         'FLQuant' %in% class(catch))
        res=fwdFn(object,ctrl=ctrl,
              catch,harvest,stock,pe,peMult,minF,maxF,bounds,lag,end,
              starvationRations=starvationRations,...) 
     else if ('FLQuants' %in% class(stock))  
        res=biodyns(llply(stock, function(x) fwdFn(object,ctrl=ctrl,
                 catch,harvest,x,pe,peMult,minF,maxF,bounds,lag,end,...)))
      else if ('FLQuants' %in% class(harvest))  
        res=biodyns(llply(harvest, function(x) fwdFn(object,ctrl=ctrl,
                 catch,x,stock,pe,peMult,minF,maxF,bounds,lag,end,...))) 
      else if ('FLQuants' %in% class(catch))  
        res=biodyns(llply(catch, function(x) fwdFn(object,ctrl=ctrl,
                x,harvest,stock,pe,peMult,minF,maxF,bounds,lag,end,...)))
                        
     res})
     
fwdFn(object, 
               catch, 
               harvest, 
               stock, 
               pe, peMult,
               minF,    maxF,
               starvationRations,
               ...){
        
      lag=0
      object@stock=FLQuant(object@stock,quant=names(object@catch)[1])   
    
      // catch, harvest or stock?
      ctcTrgt=FALSE
      hvtTrgt=FALSE
      stkTrgt=FALSE
      hcrTrgt=FALSE
      
      if (!is.null(catch))   ctcTrgt=TRUE 
      if (!is.null(harvest)) hvtTrgt=TRUE
      if (!is.null(stock))   stkTrgt=TRUE
      // hcr option not implemented, see hcr method
      
      if (hvtTrgt) f=harvest
      
      if(!ctcTrgt & hvtTrgt & stkTrgt & hcrTrgt)
          stop('must supply catch, harvest or stock as a target')
    
      if (ctcTrgt) if (dims(object@stock)$maxyear   < dims(  catch)$maxyear) object=window(object,end=dims(  catch)$maxyear)
      if (hvtTrgt) if (dims(object@stock)$maxyear-1 < dims(harvest)$maxyear) object=window(object,end=dims(harvest)$maxyear+1)
      if (stkTrgt) if (dims(object@stock)$maxyear   < dims(  stock)$maxyear) object=window(object,end=dims(  stock)$maxyear)
       
      if (stkTrgt) catch=stock*0
     
      // check year range
      if (ctcTrgt | stkTrgt) {
        if (!(all(dimnames(catch)$year %in% dimnames(object@catch)$year)))
             object = window(object,end=dims(catch)$maxyear)
          if (dims(object@catch)$iter==1 & dims(catch)$iter>1) object@catch=propagate(object@catch, dims(catch)$iter)
          object@catch[,dimnames(catch)$year] <- catch
          yrs <- dimnames(catch)$year
      } else if (hvtTrgt) {
          if (!(all(dimnames(harvest)$year %in% dimnames(object@stock)$year))){
            
            stop('years in harvest & stock dont match')}
          yrs <- dimnames(harvest)$year
      } else if (hcrTrgt) {
         yrs=ac(dims(object@stock)$maxyear:end) 
         object=window(object,end=end)
      }  
     
      if (stkTrgt) yrs = yrs[-length(yrs)]
      // B0 in year 1?
      if (as.numeric(yrs[1]) == range(object,'minyear')){
         if (!('year' %in% names(dimnames(params(object)['k']))))  
           object@stock[,ac(range(object,'minyear'))] = params(object)['k'] * params(object)['b0'] else
           object@stock[,ac(range(object,'minyear'))] = params(object)['k',ac(range(object,'minyear'))] * params(object)['b0',ac(range(object,'minyear'))]          
         }
      
      // maxyear
      if (max(as.numeric(yrs)) == range(object,'maxyear'))
         object@stock <- window(object@stock,end=range(object,'maxyear')+1)
    
      // niters
      nits=dims(object)$iter
      if (hvtTrgt) nits=max(nits,dims(harvest)$iter)
      if (ctcTrgt) nits=max(nits,dims(  catch)$iter)
      if (stkTrgt) nits=max(nits,dims(  stock)$iter)
      if (!is.null(pe)) nits=max(nits,dims(pe)$iter)
      
      if (hvtTrgt) nits=max(nits,dims(harvest)$iter) else
      if (ctcTrgt) nits=max(nits,dims(catch  )$iter) else
      if (stkTrgt) nits=max(nits,dims(stock  )$iter) 
      if (nits>1){ 
         if(dim(object@catch)[6]==1) object@catch =propagate(object@catch,nits)
         if(dim(object@stock)[6]==1) object@stock =propagate(object@stock,nits)
         if (hvtTrgt) if(dim(harvest)[6]==1) harvest=propagate(harvest,nits)
       
         if (dims(params(object))$iter==1) params(object)=propagate(params(object),nits)
         //if (!is.null(pe))                 pe            =propagate(pe            ,nits)
         } 
      
      // projections
      if (!hvtTrgt) harvest=harvest(object)
      for(y in as.numeric(yrs)) {
    
         // sp & process error
         if (!is.null(pe)) {    
            if (peMult) sp.=computePrd(object,object@stock[, ac(y)])%*%pe[, ac(y)] 
            else        sp.=computePrd(object,object@stock[, ac(y)])%+%pe[, ac(y)]
         } else sp.=computePrd(object,object@stock[, ac(y)])
         //} else sp.=computePrd(object,object@stock[, ac(y)]*(1-ptYr)+object@stock[, ac(y+1)]*(ptYr))
         
          // targets, if lag<0 then the targets are relative 
          if (hcrTrgt)
              harvest[,ac(y)]=hcr(window(object,end=y))
          if (hvtTrgt | hcrTrgt){
              object@catch[,ac(y)]=object@stock[,ac(y)]*relFn(harvest(object),harvest[,ac(y)],lag)
          } else if (ctcTrgt){           
              object@catch[,ac(y)]=relFn(object@catch,catch[,ac(y)],lag)
          } else {
              object@catch[,ac(y)]=relFn(object@stock,stock[,ac(y+1)],lag) - object@stock[,ac(y)] - sp.}
         
         object@catch[,ac(y)]=qmin(object@stock[,ac(y)]*starvationRations,
                                   object@catch[,ac(y)])
                                   
         object@stock[,ac(y+1)] = object@stock[,ac(y)] - 
                                  object@catch[,ac(y)] + sp.
    
    // Not implemented     
    // bounds
    //      if ('catch' %in% bounds){ 
    //           object@catch[,y]=iavFn(window(object@catch,end=y),bnd,lag=lag)
    //           object@stock[,ac(y+1)] = object@stock[,ac(y)] - object@catch[,ac(y)] + sp.
    //           }
        
          // min/max
          harvest[,ac(y)]       =maxFn(minFn(harvest(object)[,ac(y)],minF),maxF)
          object@catch[,ac(y)]  =object@stock[,ac(y)]*harvest[,ac(y)]
          object@stock[,ac(y+1)]=object@stock[,ac(y)] - object@catch[,ac(y)] + sp.
          }
      
        object@stock[stock(object) < 0] = 0.001
        object@catch[catch(object) < 0] = 0.0000001
    
        return(object)}


setMethod('fwd', signature(object='biodyn',ctrl='FLQuants'),
  function(object, ctrl, pe=NULL, peMult=TRUE,minF=0,maxF=2,lag=0,
           bounds=list(catch=c(Inf,Inf)),...) {
    
  res=mlply(seq(length(names(ctrl))),
      function(x,object,ctrl,pe,peMult,minF,maxF,lag,bounds){
        if (names(ctrl)[x]=='catch')  return(fwd(object,catch  =ctrl[[x]],pe=pe,peMult=peMult,minF=minF,maxF=maxF,lag=lag,bounds=bounds))
        if (names(ctrl)[x]=='harvest')return(fwd(object,harvest=ctrl[[x]],pe=pe,peMult=peMult,minF=minF,maxF=maxF,lag=lag,bounds=bounds))
        if (names(ctrl)[x]=='stock')  return(fwd(object,stock  =ctrl[[x]],pe=pe,peMult=peMult,minF=minF,maxF=maxF,lag=lag,bounds=bounds))},
       object=object,ctrl=ctrl,pe=pe,peMult=peMult,minF=minF,maxF=maxF,lag=lag,bounds=bounds)
  
  names(res)=names(ctrl)
        
  return(biodyns(res))})
