bounder<-function(x,min,max,scale){
  y=x/scale;
  if (y<20.0){
    z=exp(y)/(1.0+exp(y))
  }else{
    z=1.0/(1+exp(-y))}
  
  return(min+(max-min)*z)}

x=seq(-1000,1000,length.out=201)

ggplot(data.frame(x=x,y=bounder(x,0,100,10)))+
  geom_line(aes(x,y))+geom_vline(xintercept=c(0,500))

boundPar<-function(t., fmin, fmax){
  x=asin((2*t.-fmin)/(fmax-fmin)-1)/1.570795
  
  return(x)}


boundPar<-function(t., fmin, fmax){
  x=asin((2*t.-fmin)/(fmax-fmin)-1)/1.570795
  
  return(x)}

unboundPar<-function(x, fmin, fmax){
  t=(fmin+(fmax-fmin)*(sin(x*1.570795)+1))/2
  
  return(t)}

x=seq(-10000,20000,length.out=201)

ggplot(data.frame(x=x,
                  y=unboundPar(boundPar(x,0,10000),0,10000)))+
  geom_line(aes(x,y))+
  geom_vline(xintercept=c(0,500))
