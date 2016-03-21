##PELLA , J. J., 1967 A study of methods to estimate the Schaefer model parameters with special reference to the yellowfin tuna fishery in the eastern tropical Pacific ocean. University of Washington, Seattle.
##PELLA , J. J., and P. K. T OMLINSON , 1969 A generalized stock production model. Bulletin of the Inter-American Tropical Tuna Commission 13: 419-496.
##PRAGER , M. H., 1994 A suite of extensions to a nonequilibrium surplus-production model. U. S. Fishery Bulletin 92: 374-389.

alpha=function(r,F) r-F
beta =function(r,K) r/K

yield<-function(F,B,r,K)
    F/(r/K)*log(1-(r/K)*B*(1-exp((r-F)))/(r-F))

fnY   =function(F,C,B,r,K)
  C-yield(F,B,r,K)

minY   =function(F,C,B,r,K)
  fnY(F,C,B,r,K)^2

gradY =function(F,C,B,r,K){
  .expr1 <- r/K
  .expr2 <- F/.expr1
  .expr3 <- .expr1 * B
  .expr4 <- r - F
  .expr5 <- exp(.expr4)
  .expr7 <- .expr3 * (1 - .expr5)
  .expr9 <- 1 - .expr7/.expr4
  .expr10 <- log(.expr9)
  .value <- C - .expr2 * .expr10
  -(1/.expr1 * .expr10 - .expr2 * ((.expr3 *.expr5/.expr4 + .expr7/.expr4^2)/.expr9))
  }

gradMinY=function(F,C,B,r,K)
  2*minY(F,C,B,r,K)*gradY(F,C,B,r,K)

fn<-function(F,C,B,r,K)
  C*(r/K)/(log(r/K*B*(exp(r-F)-1)/(r-F)+1))

fnF   =function(F,C,B,r,K)
  F-fn(F,C,B,r,K)

minF   =function(F,C,B,r,K)
  fnF(F,C,B,r,K)^2

gradF=function(F,C,B,r,K){
  .expr1 <- r/K
  .expr2 <- C * .expr1
  .expr3 <- .expr1 * B
  .expr4 <- r - F
  .expr5 <- exp(.expr4)
  .expr7 <- .expr3 * (.expr5 - 1)
  .expr9 <- .expr7/.expr4 + 1
  .expr10 <- log(.expr9)
  1 - .expr2 * ((.expr3 * .expr5/.expr4 - .expr7/.expr4^2)/.expr9)/.expr10^2
  }

gradMinF=function(F,C,B,r,K)
  2*minF(F,C,B,r,K)*gradF(F,C,B,r,K)

NewRhap=function(x, func, grad)
  x - func / grad

prj<-function(F,C,B,r,K)
  alpha(r,F)*B*exp(alpha(r,F))/(alpha(r,F)+beta(r,K)*B*(exp(alpha(r,F))-1))

nr<-function(C,F,B,r,K,tolVal=1e-5,niter=100,yield=TRUE){
  
  for (t in rev(rev(dimnames(B)$year)[-1])){
       iters=0;tol=1
       while (tol>tolVal&iters<niter){
         iters  =iters+1
         if (yield){
           func   = fnY(F[,t],C[,t],B[,t],r,K) 
           grad   = gradY(F[,t],C[,t],B[,t],r,K)
         }else{
           func   = fnF(F[,t],C[,t],B[,t],r,K) 
           grad   = gradF(F[,t],C[,t],B[,t],r,K)
           }
         
         val=NewRhap(F[,t],func,grad)
        
         if (is.na(val)) 
           tol=tolVal
         else{
         tol=abs(F[,t]-val)
         F[,t]=val}
         }
       
   B[,ac(as.numeric(t)+1)]=prj(F[,t],C[,t],B[,t],r,K)
   }
  
  return(list(F=F,B=B))}

slv<-function(C,F,B,r,K){
  for (t in rev(rev(dimnames(B)$year)[-1])){
    res=ucminf(c(F[,t]), fn = minF, #gr = gradMinF, 
               control = list(trace = -1), 
               C=C[,t],B=B[,t],r=r,K=K)
    
    F[,t]=res$par
    
    B[,ac(as.numeric(t)+1)]=prj(F[,t],C[,t],B[,t],r,K)
    }
  
  return(FLQuants(F=F,B=B))}

if (FALSE){
  grad(yield, .01,         B=B[,t],r=r,K=K)
  c(  gradY(  .01, C=C[,t],B=B[,t],r=r,K=K))
  
  grad(fnF,.01, C=C[,t],B=B[,t],r=r,K=K)
  c( gradF(.01, C=C[,t],B=B[,t],r=r,K=K))
  
  ggplot(mdply(data.frame(F=seq(0,.03,length.out=101)), 
               function(F) c(fnY(F,C[,t],B[,t],r,K)/dfdY(F,C[,t],B[,t],r,K))))+
    geom_line(aes(F,V1))
  
  ggplot(mdply(data.frame(F=seq(0,.7,length.out=101)), 
               function(F) c(fnY(F,C[,t],B[,t],r,K))))+
    geom_line(aes(F,V1))
  
  ggplot(mdply(data.frame(F=seq(0,.7,length.out=101)), 
               function(F) c(gradY(F,C[,t],B[,t],r,K))))+
    geom_line(aes(F,V1))
  }


