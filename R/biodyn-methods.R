#' harvest
#'
#' @description harvest rate
#' 
#' @param object either \emph{biodyn} or \emph{aspic} classes
#' 
#' @rdname hrate
#'
#' @export
#' 
#' @aliases harvest harvest,biodyn-method harvest,aspic-method 
#' 
setMethod('harvest', signature(object='biodyn',catch="missing"),function(object) {
             
  yrs1=  dimnames(stock(object))$year
  yrs2=c(dimnames(stock(object))$year[-1],as.numeric(max(dimnames(stock(object))$year))+1)
             
  #res <- catch(object)/(stock(object)[,yrs1]*(1-when)+
  #                        stock(object)[,yrs2]*when)
             
  yrs=dimnames(catch(object))$year[dimnames(catch(object))$year %in% dimnames(catch(object))$year]
  res <- catch(object)[,yrs]/stock(object)[,yrs]
  units(res) <- 'hr'
  return(res)
  })

setMethod('stock', signature(object='biodyn'),
          function(object,when=0) {
            
            when=max(min(when,1),0)
            if (when<=0) return(object@stock)
            
            yrs =  dimnames(stock(object))$year
            yrs1=  rev(rev(yrs)[-1])
            yrs2=  yrs[-1]
             
            (1-when)*stock(object)[,yrs1]+when*stock(object)[,yrs2]})

#' computePrd
#'
#' Calculates the surplus production for a biomass dynamic model given a level of stock biomass
#' 
#' @param object an object of class \code{biodyn} 
#' @param biomass stock biomaas, may be a \code{numerix},  \code{FLQuant} or missing. In the latte case the stock slot will be used.
#' @param ... other arguments
#'
#' @return an \code{FLPar} object
#' 
#' @seealso \code{\link{plotPrd}}
#' 
#' @export
#' @rdname sp
#'
#' @aliases computePrd,biodyn,FLQuant-method  computePrd,biodyn,missing-method  computePrd,biodyn,numeric-method
#' @examples
#' \dontrun{ computePrd(bd,seq(0,params(bd)['k'])) }
#'  
setGeneric('computePrd',   function(object,biomass,...)        standardGeneric('computePrd'))
setMethod( 'computePrd', signature(object='biodyn',   biomass='missing'),     function(object,biomass=stock(object))  prdFn(model(object),params(object),biomass))
setMethod( 'computePrd', signature(object='biodyn',   biomass='numeric'),     function(object,biomass)                prdFn(model(object),params(object),biomass))
setMethod( 'computePrd', signature(object='biodyn',   biomass='FLQuant'),     function(object,biomass)                prdFn(model(object),params(object),biomass))


# calcLogLik

calcSigma <- function(obs,hat=rep(1,length(obs)),error='log'){
  yrs=dimnames(obs)$year
  yrs=yrs[yrs %in% dimnames(hat)$year]
  hat=hat[,yrs]
  obs=obs[,yrs]
  
  if (error=='log'){
    hat=log(hat)
    obs=log(obs)}
  
  SS =sum((obs-hat)^2,na.rm=T)
  
  return((SS/length(hat))^.5)}

loglFn<-function(obs,se,hat=rep(0,length(obs))){
  flag=!is.na(obs) & !is.na(hat)
  obs =obs[flag]
  hat =hat[flag]
  
  SS<-sum((obs-hat)^2)
  
  n   <-length(obs)
  res <-(log(1/(2*pi))-n*log(se)-SS/(2*se^2))/2
  
  return(res)}

calcLogLik<-function(obs,hat=rep(0,length(obs)),error='log',type=1){
  
  yrs=dimnames(obs)$year
  yrs=yrs[yrs %in% dimnames(hat)$year]
  hat=hat[,yrs]
  obs=obs[,yrs]
  
  hat=hat[!(is.na(obs)&is.na(hat))]
  obs=obs[,dimnames(hat)$year]
  
  if (error=='log'){
    hat=log(hat)
    obs=log(obs)}
  
  se<-calcSigma(obs,hat)
  
  if (type==1) return(loglFn(se,obs,hat)) else
    if (type==2) return(-sum(dnorm(obs, hat, se, log=(error=='log')))) else
      if (type==3) return(sum((obs-hat)^2))}


utils::globalVariables('calcSigmaFLQ')

# calcSigma
calcSigma=function(obs,hat=rep(0,length(obs)),error='log'){
  SS   <-sum((obs-hat)^2,na.rm=T)
  
  return((SS/length(hat))^.5)}

llSigma=function(obs,hat=obs*0,dims=c(1,3:6)){
  
  hat=hat[,dimnames(obs)$year]
  
  SS =apply(obs-hat, dims, function(x) sum(x^2,na.rm=T))
  n  =apply(obs-hat, dims, function(x) sum(!is.na(x)))
  
  return((SS/n)^.5)}

llQ=function(obs,hat,dims=c(1,3:6),error='log'){
  
  yrs=dimnames(obs)$year[dimnames(obs)$year %in% dimnames(hat)$year]
  obs=obs[,yrs]
  
  res=switch(error,
             normal={q    =  apply(hat*obs,dims, function(x) sum(x))
                     q    =q/apply(hat,    dims, function(x) sum(hat*hat))
                     sigma=calcSigmaFLQ(obs/(q%*%hat))
                     
                     FLQuants(q=q,sigma=sigma)},
             log   ={q    =apply(log(obs)-log(hat), dims, 
                                 function(x) exp(sum(x,na.rm=T)/sum(!is.na(x))))
                     sigma=llSigma(log(obs),log(q%*%hat))
                     
                     FLQuants(q=q,sigma=sigma)},
             cv   ={res   =apply(obs/hat, dims, sum) #bug!
                    sigma2=llSigma(res)
                    q     =(-res+(res^2+4*apply(obs,dims,function(x) sum(!is.na()))*sigma2*
                                    apply(obs/hat, dims, function(x) sum((x)^2)))/
                              (2*apply(obs,dims,function(x) sum(!is.na()))*sigma2))          
                    
                    FLQuants(q=q,sigma=sigma)})
  
  return(res)}

# calcLogLik
calcLl<-function(obs,hat=obs*0,error='log',type=1){
  
  hat=hat[,dimnames(obs)$year]
  
  logl<-function(se,obs,hat=obs*0,dims=c(1,3:6)){
    
    SS  =apply(obs-hat, dims, function(x) sum((x)^2,na.rm=T))    
    n   =apply(obs, dims, function(x) sum(!is.na(x)))
    
    res =(log(1/(2*pi))-n*log(se)-SS/(2*se^2))/2
    
    return(res)}
  
  se=llSigma(obs,hat)
  
  if (type==1) return(logl(se,obs,hat))
  if (type==2) return(apply(obs-hat, dims, 
                            function(x) {
                              se=llSigma(x)  
                              -sum(dnorm(x, 0, se, log=(error=='log')))}))
  if (type==3) return(apply(obs-hat, dims, function(x) sum(x^2)))
  
}

