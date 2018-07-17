if (FALSE){
    r/(m-1)*B*(1-(B/K)^(m-1))
    
    k=1000
    r=0.5
    p=-0.75
    B=seq(1,1000,10)
    r/p*B*(1-(B/k)^p)
    
    plot(r/p*B*(1-(B/k)^p))
    }

p<-function(shape){
  
  calcP<-function(shape){
    
    fn<-function(x,y)
      (y-(1/(1+x))^(1/x))^2
  
    optimise(fn,c(-0.9999,10),y=shape)$minimum}
  
  res=aaply(shape,seq(length(dim(shape)))[-1], calcP)
  
  dmns=dimnames(shape)
  dmns[[1]]="p"
  
  FLPar(array(res,dim=unlist(laply(dmns,length)),dimnames=dmns))}

if (FALSE){
  params=FLPar(t(array(c(seq(125,750,length.out=10)/rep(1000,10)),
                       dim=c(10,1),
                       dimnames=list(iter=seq(10),params=c("shape")))))
  
  p(params)}
