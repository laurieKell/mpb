setMethod('refpts', signature(object='aspic'),
  function(object){

    if (x@model=="LOGISTIC"){
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

setMethod('msy', signature(x='aspic'),
          function(x){

            if (x@model=="LOGISTIC"){
              msy =params(x)["msy"]}
            else  if (x@model=="FOX"){
              msy =params(x)["msy"]}

            msy})

setMethod('fmsy', signature(x='aspic'),
          function(x){

            if (x@model=="LOGISTIC"){
              fmsy=params(x)["msy"]/params(x)["k"]*4
              }
            else  if (x@model=="FOX"){
              msy =params(x)["msy"]
              bmsy=params(x)["k"]*exp(-1)
              fmsy=msy/bmsy}

            dimnames(fmsy)$params="fmsy"

            fmsy})

setMethod('bmsy', signature(x='aspic'),
          function(x){

            if (x@model=="LOGISTIC"){
              msy =params(x)["msy"]
              fmsy=msy/params(x)["k"]*4
              bmsy=msy/fmsy
            }else  if (x@model=="FOX"){
              bmsy=params(x)["k"]*exp(-1)
            }

            dimnames(bmsy)$params="bmsy"

            bmsy})
