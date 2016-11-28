utils::globalVariables('FLBRP')
utils::globalVariables('computeRefpts')
utils::globalVariables('rbind.fill')
utils::globalVariables('data')
utils::globalVariables('quantity')
utils::globalVariables('FLBRP')
utils::globalVariables('computeRefpts')
utils::globalVariables('optimise')
utils::globalVariables('catch.obs')
utils::globalVariables('ages')

ma<-function(x,n=5,sides=2){filter(x,rep(1/n,n), sides=sides)}

#' @title ts
#'
#' @description calculates time series of quantities useful for management
#' 
#' @param object either  \emph{FLStock} \emph{biodyn} or \emph{aspic} classes
#' 
#' @rdname timeSeries2
#' 
#' @aliases 
#' ts,FLStock,FLBRP-method 
#' ts,FLStock,FLPar-method 
#' ts,FLStocks,FLBRP-method 
#' ts,FLStocks,FLBRPs-method 
#' ts,biodyn,FLPar-method
#'
#' @export
#' 
#' 
setMethod( 'timeSeries', signature(object='FLStock',params="FLPar"), 
  function(object,params,df=TRUE) {

  object=window(object,start=range(object)["minyear"],end=range(object)["maxyear"])
  maxyr=unlist(dims(object)["year"])

	abs=FLQuants(
		ssb     =ssb( object),
 		f       =fbar(object),
 		h       =catch(object)%/%stock(object),
 		prod    =catch(object)[,-maxyr]%-%(-stock(object)[,-maxyr]+stock(object)[,-1]),
 		catch   =catch(object),
 		rec     =rec(object),
 		biomass =stock(object),
		juvenile=(stock(object)%-%ssb( object)))


  rel=FLQuants(
		ssb     =ssb( object)%/%params[1,"ssb"],
 		f       =fbar(object)%/%params[1,"harvest"],
 		h       =(catch(object)%/%stock(object))%/%(params[1,"yield"]%/%params[1,"biomass"]),
 		prod    =(catch(object)[,-maxyr]%-%(-stock(object)[,-maxyr]+stock(object)[,-1]))%/%params[1,"yield"],
 		catch   =catch(object)%/%params[1,"yield"],
 		rec     =rec(object)%/%params[1,"rec"],
 		biomass =stock(object)%/%params[1,"biomass"],
		juvenile=(stock(object)%-%ssb( object))%/%(params[1,"biomass"]%-%params[1,"ssb"]))

  dimnames(abs[["rec"]])[[1]]="all"
  dimnames(rel[["rec"]])[[1]]="all"
  
  if (df){
    
  res=rbind(data.frame(quantity="relative",model.frame(mcf(rel),drop=TRUE)),
            data.frame(quantity="absolute",model.frame(mcf(abs),drop=TRUE)))
  
  return(res)}
  else return(list(relative=rel,absolute=abs))
  })

setMethod( 'timeSeries', signature(object='FLStock',params="FLBRP"), 
      function(object,params,df=TRUE,ref="msy") {
          ts(object,FLBRP::refpts(params)[ref])})
  
setMethod( 'timeSeries', signature(object='biodyn',params="missing"), 
  function(object,params,df=TRUE,relative=TRUE) {

  object=window(object,start=range(object)["minyear"],end=range(object)["maxyear"])

  pFn<-function(x) {
    sYrs=dimnames(stock(x))$year
    cYrs=dimnames(catch(x))$year
    cYrs=cYrs[cYrs%in%yrs[-length(cYrs)]]
     
    stock(x)[,sYrs[-length(sYrs)]]-stock(x)[,-1]+catch(x)[,cYrs]}
  
  if (relative)
    res=model.frame(mcf(FLQuants(
  		h       =harvest(object)%/%refpts(object)["fmsy"],
  		prod    =pFn(object)%/%production(object,stock(object)[,-dim(stock(object))[2]]),
  		catch   =catch(object)%/%refpts(object)["msy"],
  		biomass=stock(object)%/%refpts(object)["bmsy"])))
	else	
  	res=model.frame(mcf(FLQuants(
  		h       =harvest(object),
  		prod    =pFn(object),
  		catch   =catch(object),
  		biomass=stock(object))))
             
  if (df)
    return(data.frame(res))
  else  
    return(res)
  })

setMethod( 'timeSeries', signature(object='aspic',params="missing"), 
    function(object,params,df=TRUE,relative=TRUE) {
             
      object=window(object,start=range(object)["minyear"],end=range(object)["maxyear"])
             
      pFn<-function(x) {
        stock(x)[,-dim(stock(x))[2]]-stock(x)[,-1]+catch(x)[,-dim(stock(x))[2]]}
      
      if (relative)
       res=model.frame(mcf(FLQuants(
           h       =harvest(object)%/%refpts(object)["fmsy"],
           prod    =pFn(object)%/%production(object,stock(object)[,-dim(stock(object))[2]]),
           catch   =catch(object)%/%refpts(object)["msy"],
           biomass=stock(object)%/%refpts(object)["bmsy"])))
      else	
       res=model.frame(mcf(FLQuants(
           h       =harvest(object),
           prod    =production(object,stock(object)),
           catch   =catch(object),
           biomass=stock(object))))
             
      if (df)
        return(data.frame(res))
      else  
        return(res)
      })

setMethod( 'timeSeries', signature(object='FLStocks',params="FLBRP"), 
           function(object,params,df=TRUE,ref="msy") {
             mdply(data.frame(.id=names(object)),function(.id) 
               ts(object[[.id]],FLBRP::refpts(params)[ref],df=TRUE))})

setMethod( 'timeSeries', signature(object='FLStocks',params="FLBRPs"), 
           function(object,params,ref="msy") {
             mdply(data.frame(.id=names(object)),function(.id) 
               ts(object[[.id]],FLBRP::refpts(params[[.id]])[ref]))})

#' @title mng
#'
#' @description calculates time series of quantities useful for management
#' 
#' @param object either  \emph{FLStock} \emph{biodyn} or \emph{aspic} classes
#' @aliases  timeSeries timeSeries-method timeSeries,FLStock,FLBRP-method timeSeries,FLStocks,FLBRP-method timeSeries,FLStocks,FLBRPs-method timeSeries,aspic,missing-method timeSeries,biodyn,missing-method
#' 
#' @rdname timeSeries1
#'
#' @export

# mngFn<-function(object,lyr=NULL,ryr=-(10:12),syr=-(0:4)) {
#   
#   if(is.null(lyr))
#     lyr=dim(object)$maxyear
#   
#   if (ryr<0) 
#     ryr=ac(as.numeric(lyr)+ryr)
#   
#   if (syr<0) 
#     syr=ac(as.numeric(lyr)+syr)
#   
#   xy=0.0; 
#   x =0.0;
#   y=0.0; 
#   xx=0.0; 
#   for (int i=nc; i>nc-ref[1]; i--){
#     x +=i
#     xx+=i*i  
#     y +=B[i]  
#     xy+=i*B[i]  
#     }
#   
#   slopeb = -(ref[1]*xy - x*y)/(ref[1]*xx - x*x);
# 
#   
#   #"bnow","bnowthen","slopeb"
#              
#     }
# 
# setMethod( 'mng', signature(object='FLStock',params="FLPar"), 
#     function(object,params,lyr=NULL,ryr=-(10:12),syr=-(0:4)) {
# 
#   res=ts(object,params)
# 	
#   lyr=NULL
# 	ryr=-(10:12)
# 	syr=-( 0:4)
# 	
# 	nms=c("r", "k",
# 		  "bnow","fnow","bnowthen","fnowthen",
# 		  "msy","bmsy","fmsy","cmsy",
# 		  "bmsy","ffmsy","bk","fr",
# 		  "slopeb","slopef")
# 
# 	if (ryr<0) 
#     ryr=range(params)["maxyear"]-ryr
# 
# 	#r k
# 	#bnow fnow
# 	#bnowthen fnowthen
# 	#msy bmsy fmsy cmsy
# 	#bmsy ffmsy 
# 	#bk fr
# 	#slopeb slopef
# 
# 	})

fapexAge<-function(object){
  tmp=harvest(object)
  tmp[]=fapex(object)
  tmp=FLQuant(ages(tmp)[harvest(object)==tmp],
              dimnames=dimnames(fapex(object)))
  tmp}



