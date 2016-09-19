setMethod('refpts', signature(object='aspic', params='missing'),  
  function(object,params=NULL){

    if (object@model=="LOGISTIC"){
      msy =params(object)["msy"]
      fmsy=msy/params(object)["k"]*2
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

setMethod('msy', signature(object='aspic'),  
          function(object,params=NULL){
            
            if (object@model=="LOGISTIC"){
              msy =params(object)["msy"]}
            else  if (object@model=="FOX"){
              msy =params(object)["msy"]}
            
            msy})

setMethod('fmsy', signature(object='aspic'),  
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

setMethod('bmsy', signature(object='aspic'),  
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
