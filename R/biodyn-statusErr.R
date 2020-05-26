statusErrFn<-function(x,y,rfs){
  
  lastYear=ac(dims(harvest(y))$maxyear)
  mp =(stock(y)%/%refpts(y)["bmsy"])[,lastYear]
  om =(ssb(x)%/%rfs["msy","ssb"])[,lastYear]
  bmsy=model.frame(FLQuants(om=om,mp=mp),drop=TRUE)
  
  mp =(harvest(y)%/%refpts(y)["fmsy"])[,lastYear]
  om =(fbar(x)%/%rfs["msy","harvest"])[,lastYear]
  fmsy=model.frame(FLQuants(om=om,mp=mp),drop=TRUE)
  
  res=cbind("relative"=TRUE,rbind(cbind(quantity="harvest",fmsy),cbind(quantity="stock",bmsy)))
  
  mp   =stock(y)[,lastYear]
  om   =ssb(x)[,lastYear]
  stock=model.frame(FLQuants(om=om,mp=mp),drop=TRUE)
  
  mp     =harvest(y)[,lastYear]
  om     =fbar(x)[,lastYear]
  harvest=model.frame(FLQuants(om=om,mp=mp),drop=TRUE)
  
  res2=cbind("relative"=FALSE,rbind(cbind(quantity="harvest",harvest),cbind(quantity="stock",stock)))
  
  rbind(res,res2)}