setMP<-function(om,r,k,p,b0){
  
  mp=as(om,"biodyn")
  
  params(mp)["r"] =r
  params(mp)["k"] =k 
  params(mp)["p"] =p
  params(mp)["b0"]=b0
  params(mp)      =iter(params(mp),1)
  
  u               =stock(mp,0.5)
  setParams( mp)  =u
  setControl(mp)  =params(mp)
  
  mp@control["q1",c("phase","val")]=c(-1,1)

  hat             =fit(mp,u)
  
  print(plot(as(list("MP"=hat,"OM"=mp),"biodyns"),probs=c(0.25,0.75)))
  
  params(mp)=iter(params(hat),1)
  params(mp)=apply(params(mp),1,median)
  setControl(mp)=params(mp)
  
  mp}
