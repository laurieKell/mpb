setGeneric("jackSummary", function(object,sim,...) standardGeneric("jackSummary"))

setMethod("jackSummary", signature(object="FLQuant",sim="FLQuant"),
  function(object,sim,...) {
            
    n   <- dims(sim)$iter
            
    mn  <- object
    u   <- sim
    mnU <- apply(u, 1:5, mean)   
            
    SS  <- apply((u%-%mnU)^2, 1:5, sum)
    #SS  <- apply(sweep(u, 1:5, mnU,"-")^2, 1:5, sum)
    
    bias<- (n-1)*(mn-mnU)
    se  <- sqrt(((n-1)/n)*SS)
            
    return(FLQuants(mean=mn, se=se, bias=bias))
    })

setMethod("jackSummary", signature(object="FLQuant",sim="FLQuant"),
  function(object=as(object,"FLQuantJK"),sim,...) {
    jackSummary(object,sim)})

setMethod("jackSummary", signature(object="FLQuant"),
  function(object,sim,...) {
            
    n   <- dims(sim)$iter
            
    mn  <- object
    u   <- sim
    mnU <- apply(u, 1:5, mean)   
            
    SS  <- apply(sweep(u, 1:5, mnU,"-")^2, 1:5, sum)
            
    bias<- (n-1)*(mn-mnU)
    se  <- sqrt(((n-1)/n)*SS)
            
   return(FLQuants(mean=mn, se=se, bias=bias))
   })

setMethod("jackSummary", signature(object="FLPar",sim="FLPar"),
          function(object,sim,...) {
            
            nms <-names(dimnames(object))
            idx <-seq(length(nms))[nms != 'iter']
            n   <-dims(sim)$iter 
            
            mn  <-object
            u   <-sim
            mnU <-apply(u, idx, mean)   
            
            SS <-apply(sweep(u, idx, mnU,"-")^2, idx, sum)
            
            bias <- (n-1)*(mn-mnU)
            se   <- sqrt(((n-1)/n)*SS)
            
            cov  <-FLPar(cov(model.frame(u)[,dimnames(u)[[1]]])*(n-1)*(n-1)/n)
            cor  =FLPar(cor(cov[drop=T]))
            
            return(FLPars("hat"=mn, "mean"=mnU, "se"=se, "cv"=se%/%mnU,
                          "bias"=bias, "biasRel"=bias%/%mnU,
                          "cov"=cov, "cor"=cor))})

setMethod("jackSummary", signature(object="FLPar"),
  function(object,sim,...) {
            
   nms <-names(dimnames(object))
   idx <-seq(length(nms))[nms != 'iter']
   n   <-dims(sim)$iter 
            
   mn  <-object
   u   <-sim
   mnU <-apply(u, idx, mean)   
            
   SS <-apply(sweep(u, idx, mnU,"-")^2, idx, sum)
            
   bias <- (n-1)*(mn-mnU)
   se   <- sqrt(((n-1)/n)*SS)
   
   cov  <-FLPar(cov(model.frame(u)[,dimnames(u)[[1]]])*(n-1)*(n-1)/n)
   
   return(FLPars(hat=mn, mean=mnU, se=se, bias=bias, cov=cov))})

jackSmryFn<-
  function(x, theta, ...)
  {
    call <- match.call()
    n <- length(x)
    u <- rep(0, n)
    for(i in 1:n) {
      u[i] <- theta(x[ - i], ...)
    }
    thetahat <- theta(x, ...)
    jack.bias <- (n - 1) * (mean(u) - thetahat)
    jack.se <- sqrt(((n - 1)/n) * sum((u - mean(u))^2))
    return(list(jack.se=jack.se, 
                jack.bias=jack.bias, 
                jack.values = u, 
                call=call))
  }

setGeneric("jackSummary", function(object,sim,...) standardGeneric("jackSummary"))


if (FALSE){
  # jackknife values for the sample mean 
  # (this is for illustration;  # since "mean" is  a 
  #  built in function,  jackknife(x,mean) would be simpler!)
  x       <- rnorm(20)               
  theta   <- function(x){mean(x)}
  results <- jackknife(x,theta)        
  
  # To jackknife functions of more  complex data structures, 
  # write theta so that its argument x
  #  is the set of observation numbers  
  #  and simply  pass as data to jackknife the vector 1,2,..n. 
  # For example, to jackknife
  # the correlation coefficient from a set of 15 data pairs:      
  
  xdata   <- matrix(rnorm(30),ncol=2)
  n       <- 15
  theta   <- function(x,xdata){ cor(xdata[x,1],xdata[x,2]) }
  results <- jackknife(1:n,theta,xdata)
}


covJK=function(object,...) {
            
  n  =dims(object)$iter 

  FLPar(cov(model.frame(object)[,dimnames(object)$params])*(n-1)*(n-1)/n)}

corJK=function(object,...)   {
  res=covJK(object)
  FLPar(cor(res[drop=T]))}
