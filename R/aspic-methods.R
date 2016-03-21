setMethod('harvest', signature(object='aspic'),function(object,when=.5,...) {
             
  yrs1=  dimnames(stock(object))$year
  yrs2=c(dimnames(stock(object))$year[-1],as.numeric(max(dimnames(stock(object))$year))+1)
             
  #res <- catch(object)/(stock(object)[,yrs1]*(1-when)+
  #                        stock(object)[,yrs2]*when)
             
  yrs=dimnames(catch(object))$year[dimnames(catch(object))$year %in% dimnames(catch(object))$year]
  res <- catch(object)[,yrs]/stock(object)[,yrs]
  units(res) <- 'hr'
  return(res)
  })

