#' @title boot, Bootstraps biodyn.
#'
#' @description 
#' Bootstraps the aspic model
#'
#' @param object; a \code{biodyn} object 
#' 
#' @seealso \code{\link{biodyn},\link{boot},\link{jk}}
#' 
#' @aliases boot-method boot,biodyn-method boot,aspic-method boot,aspics-method
#' 
#' @export
#' 
#' @rdname boot
#' @examples
#' \dontrun{
#'     data(asp)
#'     asp=boot(asp)}
setMethod('boot', signature(object='biodyn'),
          function(object,run=TRUE)
            bootFn(object=object, n=500, run=run))

bootFn<-function(object=object, n=500, run=TRUE){
   dgs=object@diags
   
   fn<-function(x){   
     hat=as.FLQuant(transform(dgs,data=hat)[,     c("data","year")])
     rsd=as.FLQuant(transform(dgs,data=residual)[,c("data","year")])
   
     smp=FLQuant(sample(c(rsd),n*dim(rsd)[2],replace=TRUE),
             dimnames=list(year=dimnames(rsd)$year,iter=seq(n)))
   
     hat*exp(smp)}
   
   FLQuants(dlply(dgs,.(name),fn))}


bootPar<-function(object){
  
}

bootFn2<-function(n,object,index){
  
  q=object@params[grep("q",dimnames(object@params)[["params"]]),]
  s=object@params[grep("s",dimnames(object@params)[["params"]]),]
  
  res=FLQuants(mlply(seq(length(index)), function(x) {
    res=q[x]%*%rlnorm(n,FLQuant(0,dimnames=dimnames(index[[x]])),c(s[x]))
    res=res%*%stock(object)[,dimnames(index[[x]])$year]
  }))
  
  names(res)=names(index)               
  
  res}

# data(bd)
# uBoot=bootFn(100,bd,cpue[4])

   

