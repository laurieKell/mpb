utils::globalVariables(c('profileGrid'))
utils::globalVariables('maply')

setMethod('profile', signature(fitted='biodyn'),
      function(fitted, 
               which,
               index="missing",
               range=seq(0.5,1.5,length.out=21),
               fn   =function(x) mdply(data.frame(index=dimnames(x@ll)$index),
                                   function(index)      
                                     cbind(model.frame(params(x)),
                                           model.frame(x@ll[,index]),
                                           model.frame(refpts(x))[,-4],
                                           stock  =c(stock(  x)[,ac(range(x)['maxyear'])]%/%refpts(x)["bmsy"]),
                                           harvest=c(harvest(x)[,ac(range(x)['maxyear'])]%/%refpts(x)["fmsy"]))),
                   run  =TRUE,
                   comp =FALSE,...){
  
  if (is.FLQuant(index)) index=FLQuants(index)
  
  for (i in seq(length(index)))
    if (dims(index[[i]])$maxyear>=dims(stock(fitted))$maxyear) stop('index years greater in length than stock')
        
  if (dims(catch(fitted))$iter>1) catch(fitted)=iter(catch(fitted),1) #stop('can only be done for a single iter')
        
  if (dim(fitted@control)[3]==1){
     fitted@control=propagate(fitted@control,length(range)^length(which))
          
     sq=list(range)
     sq=do.call('expand.grid',sq[rep(1,length(which))])
       
     for (i in seq(length(which))){
       fitted@control[which[i],'val']=     params(fitted)[which[i]]*sq[,i]
       fitted@control[which[i],'min']=min(fitted@control[which[i],'val'])*range[1]
       fitted@control[which[i],'max']=max(fitted@control[which[i],'val'])*range[2]*2}

     fitted@control[which,'phase']=-1
  }else{
  
     fitted@catch      =iter(catch(fitted),1)
     fitted@stock      =iter(stock(fitted),1)
     setControl(fitted)=params(fitted)
           
     params(fitted)    =iter(params(fitted),1)

     fitted@control    =profileGrid(fitted@control,which,range)
           
     vcov(fitted)      =propagate(iter(vcov(fitted),  1),dims(fitted@control)$iter)
     fitted@hessian    =propagate(iter(fitted@hessian,1),dims(fitted@control)$iter)
     fitted@mng        =propagate(iter(fitted@mng,    1),dims(fitted@control)$iter)
     fitted@mngVcov    =propagate(iter(fitted@mngVcov,1),dims(fitted@control)$iter)
  }
        
  if (!run) return(fitted)
        
  f=fitted
  f@catch=propagate(f@catch,dim(f@control)[3])
  res=fit(f,index)
  res@catch=FLCore::iter(res@catch,1)
  rtn=fn(res)
      
  if (comp){
    rsd=mdply(data.frame(iter=seq(dims(f)$iter)), 
            function(iter){
              ft =fit(iter(fitted,iter),index)
              dgs=ft@diags[,c(".id","year","residual")]
              names(dgs)=c("name","year","residual")
                    
              dgs})
          
    return(list(comp=rtn,residuals=rsd))}
        
  return(rtn)})

# # CIs
# cis <- max(surface) - qchisq(ci, 2)
# 
# do.call('contour', list(x=sort(profiled[[1]]), y=sort(profiled[[2]]), z=surface,
#                             levels=cis, add=TRUE, labcex=0.8, labels=ci))
# 
  

profileGrid=function(object,which,range=seq(0.95,1.05,length.out=11)){
    
    res=maply(range,function(x,which,ctl) {
      ctl[which,'val']  =ctl[which,'val']*x
      ctl[which,'phase']=-1
      ctl}, which=which, ctl=object)
    
    res=aperm(res,c(2:4,1))
    
    dmns=dimnames(object)
    dmns$iter=seq(dim(object)[3]*length(range))
    
    object=FLPar(array(c(res),unlist(lapply(dmns,length)),dmns))
    
    return(object)}
  
  