globalVariables(c("slv","ucminf"))

#PELLA , J. J., 1967 A study of methods to estimate the Schaefer model parameters with special reference to the yellowfin tuna fishery in the eastern tropical Pacific ocean. University of Washington, Seattle.
#PELLA , J. J., and P. K. TOMLINSON , 1969 A generalized stock production model. Bulletin of the Inter-American Tropical Tuna Commission 13: 419-496.
#PRAGER , M. H., 1994 A suite of extensions to a nonequilibrium surplus-production model. U. S. Fishery Bulletin 92: 374-389.

alpha=function(r,F) r-F
beta =function(r,k) r/k

yield<-function(F,B,r,k)
    (F/(r/k))*log(1-(r/k)*B*(1-exp((r-F)))/(r-F))

fnY   =function(F,C,B,r,k)
  C-yield(F,B,r,k)

minY   =function(F,C,B,r,k)
  fnY(F,C,B,r,k)^2

gradY =function(F,C,B,r,k){
  .expr1  <- r/k
  .expr2  <- F/.expr1
  .expr3  <- .expr1*B
  .expr4  <- r-F
  .expr5  <- exp(.expr4)
  .expr7  <- .expr3 * (1 - .expr5)
  .expr9  <- 1 - .expr7/.expr4
  .expr10 <- log(.expr9)
  .value  <- C - .expr2 * .expr10
  -(1/.expr1 * .expr10 - .expr2 * ((.expr3 *.expr5/.expr4 + .expr7/(.expr4^2))/.expr9))
  }

gradMinY=function(F,C,B,r,k)
  2*minY(F,C,B,r,k)*gradY(F,C,B,r,k)

fn<-function(F,C,B,r,k)
  C*(r/k)/(log(r/k*B*(exp(r-F)-1)/(r-F)+1))

fnF   =function(F,C,B,r,k)
  F-fn(F,C,B,r,k)

minF   =function(F,C,B,r,k)
  fnF(F,C,B,r,k)^2

gradF=function(F,C,B,r,k){
  .expr1 <- r/k
  .expr2 <- C*.expr1
  .expr3 <- .expr1*B
  .expr4 <- r-F
  .expr5 <- exp(.expr4)
  .expr7 <- .expr3*(.expr5 - 1)
  .expr9 <- .expr7/.expr4 + 1
  .expr10 <- log(.expr9)
  1 - .expr2 * ((.expr3*.expr5/.expr4 - .expr7/(.expr4^2))/.expr9)/.expr10^2
  }

gradMinF=function(F,C,B,r,k)
  2*minF(F,C,B,r,k)*gradF(F,C,B,r,k)

NewRhap=function(x, func, grad)
  x - func / grad

prj<-function(F,C,B,r,k)
  alpha(r,F)*B*exp(alpha(r,F))/(alpha(r,F)+beta(r,k)*B*(exp(alpha(r,F))-1))

nr<-function(C,B,r,k,b0,tolVal=1e-10,niter=200,yieldFlag=!TRUE){
  
  F=-log(1-C/B[,dimnames(C)$year])*.5
  for (t in rev(rev(dimnames(C)$year)[-1])){
    iters=0;tol=1;val=1
    
    while (abs(val)>tolVal&iters<niter){
      iters  =iters+1
      if (yieldFlag){
        func   = fnY(F[,t],C[,t],B[,t],r,k) 
        grad   = gradY(F[,t],C[,t],B[,t],r,k)
      }else{
        func   = fnF(F[,t],C[,t],B[,t],r,k) 
        grad   = gradF(F[,t],C[,t],B[,t],r,k)
      }
      
      val=func/grad
      
      F[,t]=F[,t]-val
      
      #print(c(func,val,c(F[,t])))
    }
    
    #print(c(B[,ac(as.numeric(t)+1)],prj(F[,t],C[,t],B[,t],r,k)))     
    B[,ac(as.numeric(t)+1)]=prj(F[,t],C[,t],B[,t],r,k)
  }
  
  return(list(F=F,B=B))}

slv<-function(C,r,k,b0){
  B=window(FLQuant(c(k*b0),dimnames=dimnames(C)),end=max(as.numeric(dimnames(C)$year))+1)
  F=-log(1-C/B[,-dim(B)[2]])
  for (t in rev(rev(dimnames(B)$year)[-1])){
    res=ucminf(c(F[,t]), fn = minF, #gr = gradMinF, 
               control = list(trace = -1), 
               C=C[,t],B=B[,t],r=r,k=k)
    
    F[,t]=res$par
    
    B[,ac(as.numeric(t)+1)]=prj(F[,t],C[,t],B[,t],r,k)
    }
  
  return(FLQuants(F=F,B=B))}

setGeneric('f', function(object,...) standardGeneric('f'))
setMethod( 'f',signature(object='biodyn'),
      function(object,tolVal=1e-10,niter=200){
          fCPP(catch(object),params(object)[c("r","k","b"),tolVal,niter])
        })

setGeneric('computeStock', function(object,...) standardGeneric('computeStock'))
setMethod( 'computeStock',signature(object='biodyn'),
           function(object,tolVal=1e-10,niter=200){
             stockCPP(harvest(object),catch(object),stock(object),
                      params(object)[c("r","k","b0"),tolVal,niter])
             })


if (FALSE){
  grad(yield, .01,         B=B[,t],r=r,k=k)
  c(  gradY(  .01, C=C[,t],B=B[,t],r=r,k=k))
  
  grad(fnF,.01, C=C[,t],B=B[,t],r=r,k=k)
  c( gradF(.01, C=C[,t],B=B[,t],r=r,k=k))
  
  ggplot(mdply(data.frame(F=seq(0,.03,length.out=101)), 
               function(F) c(fnY(F,C[,t],B[,t],r,k)/gradY(F,C[,t],B[,t],r,k))))+
    geom_line(aes(F,V1))
  
  ggplot(mdply(data.frame(F=seq(0,.7,length.out=101)), 
               function(F) c(fnY(F,C[,t],B[,t],r,k))))+
    geom_line(aes(F,V1))
  
  ggplot(mdply(data.frame(F=seq(0,.7,length.out=101)), 
               function(F) c(gradY(F,C[,t],B[,t],r,k))))+
    geom_line(aes(F,V1))
  }
