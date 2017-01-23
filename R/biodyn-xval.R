xvalFnBiodyn=function(object,index,minyear=range(object)["maxyear"]-4,maxyear=range(object)["maxyear"]){

  u =window(index,   end=minyear)
  bd=window(object,  end=minyear)
  bd=fit(bd,u)
  bd=fwd(bd,catch=catch(object)[,ac((minyear+1):maxyear)])
  
  nU     =(dim(params(object))[1]-4)/2
  biomass=(stock(bd)[,-dim(stock(bd))[2]]+stock(bd)[,-1])/2
  
  if (!("FLQuants"%in%is(index))) index=FLQuants(index)
  
  res=mdply(seq(nU), function(i)
    model.frame(mcf(FLQuants(
      hat=biomass%*%params(bd)[4+i],
      obs=index[[i]])),drop=T))
  
  res=subset(res,year>=minyear)
  names(res)[1]="index"
  
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
#' @rdname biodyn-xval
#'
#' @details Returns a data.frame with index, year, obs and hat
#' @examples
#' \dontrun{
#'  data(bd)  
#' }
setMethod('xval', signature(object='biodyn'),
          function(object,index,minyear=range(object)["maxyear"]-4,maxyear=range(object)["maxyear"]) 
            xvalFnBiodyn(object,index,minyear,maxyear))


