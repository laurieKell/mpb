xvalFnAspic=function(object,minyear=range(object)["maxyear"]-4,maxyear=range(object)["maxyear"]){

  index=index(object,TRUE)

  res=aspics(mlply(data.frame(year=(minyear:maxyear)),function(year){
    if (maxyear==year) object else{
    bd=window(object,end=year)
    bd=fit(bd)
    fwd(bd,catch=catch(object)[,ac((year+1):maxyear)])}}))
  
  names(res)=minyear:maxyear
  
  attributes(res)$split_labels=NULL
  
  # if (mf){
  #   q=ldply(res,function(x) 
  #           params(x)[substr(dimnames(params(x))[[1]],1,1)=="q"])
  #   hat=ldply(res,function(x) (stock(x)[,-dims(stock(x)[2]]+stock(x)[,-1)/2)
  #   obs=ldply(res,index)
  #   }


  res}

#' @title xval
#'
#' @description Performs a cross-validation uisng a hindcast
#' 
#' @param   object an object of class \code{biodyn}
#' @param   index an \code{FLQuant} or \code{FLQuants} with index of relative stock abundance
#' @param   minyear last year to fit data to
#' @param   maxyear last year to project to, by default is the last year in the catch
#'
#' @aliases xval-method xval,biodyn-method
#'
#' @export
#' @rdname xval
#'
#' @details Returns a data.frame with index, year, obs and hat
#' @examples
#' \dontrun{
#'  data(bd)  
#' }
setMethod('xval', signature(object='biodyn',index="missing"),
          function(object,minyear=range(object)["maxyear"]-4,maxyear=range(object)["maxyear"],mf=TRUE) 
            xvalFnAspic(object,minyear,maxyear))


