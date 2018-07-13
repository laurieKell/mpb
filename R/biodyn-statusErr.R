statusErrFn<-function(x,y,rfs){
  
  lastYear=ac(dims(harvest(y))$maxyear)
  mp =(stock(y)%/%refpts(y)["bmsy"])[,lastYear]
  om =(ssb(x)%/%rfs["msy","ssb"])[,lastYear]
  bmsy=model.frame(FLQuants(om=om,mp=mp),drop=TRUE)
  
  mp =(harvest(y)%/%refpts(y)["fmsy"])[,lastYear]
  om =(fbar(x)%/%rfs["msy","harvest"])[,lastYear]
  fmsy=model.frame(FLQuants(om=om,mp=mp),drop=TRUE)
  
  res=rbind(cbind(quantity="harvest",fmsy),cbind(quantity="stock",bmsy))
  
  res}