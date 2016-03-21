setMethod('refpts', signature(object='aspic', params='missing'),  
  function(object,params=NULL){

    if (object@model=="LOGISTIC"){
      msy =params(object)["msy"]
      fmsy=msy/params(object)["k"]*4
      dimnames(fmsy)[[1]]="fmsy"
      bmsy=msy/fmsy
      dimnames(bmsy)[[1]]="bmsy"
      
      res=rbind(msy,bmsy,fmsy)}
    else  if (object@model=="FOX"){
      msy =params(object)["msy"]
      bmsy=params(object)["k"]*exp(-1)
      fmsy=msy/bmsy
      
      dimnames( msy)$params="msy"
      dimnames(bmsy)$params="bmsy"
      dimnames(fmsy)$params="fmsy"
      res=rbind(msy,bmsy,fmsy)}
     
  res})

setGeneric('msy', function(object,params,...) standardGeneric('msy')) 
setMethod('msy', signature(object='aspic', params='missing'),  
          function(object,params=NULL){
            
            if (object@model=="LOGISTIC"){
              msy =params(object)["msy"]}
            else  if (object@model=="FOX"){
              msy =params(object)["msy"]}
            
            msy})

setGeneric('fmsy', function(object,params,...) standardGeneric('fmsy')) 
setMethod('fmsy', signature(object='aspic', params='missing'),  
          function(object,params=NULL){
            
            if (object@model=="LOGISTIC"){
              fmsy=params(object)["msy"]/params(object)["k"]*4
              }
            else  if (object@model=="FOX"){
              msy =params(object)["msy"]
              bmsy=params(object)["k"]*exp(-1)
              fmsy=msy/bmsy}
            
            dimnames(fmsy)$params="fmsy"
            
            fmsy})

setGeneric('bmsy', function(object,params,...) standardGeneric('bmsy')) 
setMethod('bmsy', signature(object='aspic', params='missing'),  
          function(object,params=NULL){
            
            if (object@model=="LOGISTIC"){
              msy =params(object)["msy"]
              fmsy=msy/params(object)["k"]*4
              bmsy=msy/fmsy
            }else  if (object@model=="FOX"){
              bmsy=params(object)["k"]*exp(-1)
            }
            
            dimnames(bmsy)$params="bmsy"
            
            bmsy})
