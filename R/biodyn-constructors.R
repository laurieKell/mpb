setMethod('biodyn', signature(object='FLQuant',params='FLPar'),
    function(object,
             params,
             min=.1,max=10,...){
    
  model=tolower("pellat")

  args = list(...)

  res=new("biodyn")

  # Load given slots
  for(i in names(args))
    slot(res, i) = args[[i]]
     
  res@catch =object
  res@params=params
 
  if (any(is.na(res@params["k"]))){
    flag=is.na(res@params["k"])
    k=  guessK(mean(c(res@params["r",,flag],na.rm=TRUE)),
               mean(c(res@catch),na.rm=TRUE),
               mean(c(res@params["p",,flag]),na.rm=TRUE))
    
    res@params["k",,is.na(res@params["k"])]=k
    
    
    }
      
  res@control=propagate(res@control,dims(res@params)$iter)
  nms=dimnames(res@control)$param[dimnames(res@control)$param %in% dimnames(res@params)$param]
  res@control[nms,  'val']=res@params[nms,]
  res@control[nms,  'min']=res@params[nms,]*min
  res@control[nms,  'max']=res@params[nms,]*max
  res@control[c("b0","p"),'phase']=-1

  range(res)=unlist(dims(catch(res))[c("minyear","maxyear")])

  if (!("stock"%in%names(args)))
    res@stock=FLQuant(NA,dimnames=list(year= range(res)["minyear"]:(range(res)["maxyear"]+1)))
        
  res=fwd(res,catch=res@catch)

  return(res)})

guessK=function(r,catch,p=1)
    mean(catch)/(r*(1/(1+p))^(1/p+1))

setMethod('biodyn', signature(object='missing',params='FLPar'),
  function(object,params,model="pellat",min=.1,max=10,msy=NULL,r=NULL,...){
  args=list(...)
          
  res=biodyn(object=factor(model),
             params=params,
             #catch=ifelse("catch"%in%names(args),args[["catch"]],NA),
             min=min,max=max,...)
                               
  res})

setMethod('biodyn', signature(object='FLQuant',params='missing'),
  function(object,params=FLPar(r=0.5,k=NA,p=1,b0=1),model="pellat",min=.1,max=10,msy=NULL,r=NULL,...){
    args=list(...)
            
    res=biodyn(object=object,
               params=params,
               min=min,max=max,...)
            
    res})

# setMethod('biodyn', signature(object='FLBRP',params='FLStock'),
#           function(object,params,model="pellat",min=.1,max=10,msy=NULL,r=NULL,...){
#         
#   #if (!(is(y)%in%"FLBRP")) return(NULL)
#           
#   res=FLBRP2biodyn(object)
#   res=fwd(res,catch=catch(params)[,-1])
#         
#   res})
            

##########################################
setMethod('biodyn', signature(object='factor',params='FLPar'),
          function(object,params,min=0.1,max=10,catch=NULL,stock=NULL,msy=NULL,...){

  model=object

  if (is.null(msy) & !is.null(catch)) 
    msy=mean(catch,na.rm=TRUE)
            
  args = list(...)
            
  dimnames(params)$params=tolower(dimnames(params)$params)
  
  if (!('b0' %in%  dimnames(params)$params)) 
    params=rbind(params,propagate(FLPar('b0'=1),dims(params)$iter))

  if (model=='pellat' & !('p' %in%  dimnames(params)$params)) 
    params=rbind(params,propagate(FLPar('p'=1),dims(params)$iter))
            
  if (!('k' %in%  dimnames(params)$params))
    if (!is.null(msy) & model=='pellat') 
      params=rbind(params,'k'=FLPar(K(msy,params)))
    else  
      params=rbind(params,'k'=FLPar(k=as.numeric(NA)))

  if (model=='pellat')
      params=params[c('r','k','p','b0'),]
  
  res        =new("biodyn")
            
  if (!('factor' %in% is(model)))
      model=factor(model)
  
  res@model  =model
  res@params =params 
  if (!is.null(stock))
    res@stock[]=params(res)['k']*params(res)['b0']
            
  if (!is.null(catch)){
    res@catch=catch

    if (dim(params(res))[2]>1&dim(params(res))[2]!=dim(stock(res))[6]){
      res@stock=propagate(stock(res),dim(params(res))[2])}

    res=fwd(res,catch=catch)
  }else  if (!is.null(stock)) {
    res@stock=stock
              
    res@catch=window(stock,end=dims(stock)$maxyear-1)
    res@catch[]=NA}
            
    res@control=propagate(res@control,dims(params)$iter)
    nms=dimnames(res@control)$param[dimnames(res@control)$param %in% dimnames(res@params)$param]
    res@control[nms,  'val']=res@params[nms,]
    res@control[nms,  'min']=res@params[nms,]*min
    res@control[nms,  'max']=res@params[nms,]*max
            
    if (!('b0' %in% nms))
      res@control['b0',c('min','max','val')]=c(0.75,1,1)
            
    # Load given slots
    for(i in names(args))
      slot(res, i) = args[[i]]
            
    return(res)})

setMethod('biodyn', signature(object='character',params='FLPar'),
           function(object,params,min=0.1,max=10,catch=NULL,stock=NULL,...) {
  biodyn(object=factor(object,levels=models),params=params,min=min,max=max,catch=catch,stock=stock,...)
  })
 
setMethod('biodyn', signature(object='factor',params='missing'),
          function(object,params,min=min,max=max,catch=NULL,stock=NULL,...){

  args     =list(...)
  res      =new("biodyn")
  res@model=object
            
  nms=c(modelParams(object),'b0')
  par=rep(NA,length(nms))
  names(par)=nms
            
  res@params =FLPar(par) 
  res@params['b0']=1
            
  if (!is.null(stock))
    res@stock[]=params(res)['k']*params(res)['b0']
            
  if (!is.null(catch))
    res@catch=catch
  else  if (!is.null(stock)) {
    res@catch=window(res@stock,end=dims(res@catch)$maxyear-1)
    res@catch[]=NA}
            
  # Load given slots
  for(i in names(args))
    slot(res, i) = args[[i]]
          
  return(res)})

setMethod('biodyn', signature(object='character',params='missing'),
          function(object=object,min=0.1,max=10.0,catch=NULL,index=NULL,stock=NULL,...) 
  biodyn(object=factor(object,levels=models),min=min,max=max,catch=catch,stock=stock,...))

setMethod('biodyn', signature(object='missing',params='missing'),
          function(object,params,min=0.1,max=10.0,msy=NULL,...) {
  args = list(...)

  res=new('biodyn')
            
  # Load given slots
  for(i in names(args))
    slot(res, i) = args[[i]]
            
  range(res)=unlist(dims(catch(res))[c("minyear","maxyear")])
            
  return(res)})

setMethod('biodyn', signature(object='FLStock',params='FLPar'),
          function(object,params,
                   model="pellat",min=.1,max=10,...){
            args=list(...)

            res=biodyn(factor(model),
                       params=params,
                       #catch=ifelse("catch"%in%names(args),args[["catch"]],NA),
                       min=min,max=max,...)

            res@catch=catch(object)
            range(res)[c("minyear","maxyear")]=unlist(dims(catch(object)))[c("minyear","maxyear")]
            res@catch=catch(object)
            res@stock=object@stock
            res@stock=window(stock(res),end=range(res)["maxyear"]+1)
            
            res})

#' is.biodyn
#'
#' @description Checks class type and returns TRUE if object is of type biodyn
#' @param x biodyn class
#' 
#' @return TRUE or FALSE
#' 
#' @export
#' @examples
#' \dontrun{
#'  is.biodyn(biodyn()) 
#'  }
is.biodyn = function(x)
  return(inherits(x, 'biodyn'))

setMethod('biodyn', signature(object='missing',params='FLPar'),
     function(params,model="pellat",...){

     biodyn(factor("pellat"),params,...)
     })
