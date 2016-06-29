utils::globalVariables(c('profileGrid'))
utils::globalVariables('maply')

profileFn<-function(r=0.5,k=1000,p=0.001,b0=0.75,min=0.1,max=10){
  
  nit=max(length(r),length(k),length(p),length(b0))
  
  if (length(r) <nit) r_=rep(r, nit) else r_=r
  if (length(k) <nit) k_=rep(k, nit) else k_=k
  if (length(p) <nit) p_=rep(p, nit) else p_=p
  if (length(b0)<nit) b_=rep(b0,nit) else b_=b0
  
  val=cbind(r=r_,
            k=k_,
            p=p_,
            b0=b_)

  val=FLPar(array(c(t(val)),dim     =c(4,nit),
                         dimnames=list(params=c("r","k","p","b0"),
                                          iter  =seq(nit))))
  
  ctl=propagate(control(biodyn()),dim(val)[2])
  ctl[,c("min","val","max")]=val
  
  if (length(r)>1){
    ctl["r","phase"]=-1
    ctl["r",c("min","max")]=ctl["r","val"]}
  if (length(k)>1){
    ctl["k","phase"]=-1
    ctl["k",c("min","max")]=ctl["k","val"]}
  if (length(p)>1){
    ctl["p","phase"]=-1
    ctl["p",c("min","max")]=ctl["p","val"]}
  if (length(b0)>1){
    ctl["b0","phase"]=-1
    ctl["b0",c("min","max")]=ctl["b0","val"]}
  
  if (any(ctl["r","phase"]>0)){
    ctl["r","min",]=ctl["r","val",]*min
    ctl["r","max",]=ctl["r","val",]*max}
    
  if (any(ctl["k","phase"]>0)){
    ctl["k","min",]=ctl["k","val",]*min
    ctl["k","max",]=ctl["k","val",]*max}
  
  if (any(ctl["p","phase"]>0)){
    ctl["p","min",]=ctl["p","val",]*min
    ctl["p","max",]=ctl["p","val",]*max}

  if (any(ctl["b0","phase"]>0)){
    ctl["b0","min",]=ctl["b0","val",]*min
    ctl["b0","max",]=ctl["b0","val",]*max}
  
  ctl}

profileFn2<-function(catch,index,control){

    full=FLQuant(NA,dimnames=dimnames(ctc))
    full[,dimnames(idx)$year]=idx
    
    res=nllCpp(c(ctc),c(full),t(params@.Data))
    rtn=FLPar(array(res,c(1,dim(par)[2]),
                    dimnames=list(params="ll",iter=seq(dim(par)[2]))))
    rtn}

tmp<-function(fitted, 
           which,
           index="missing",
           range=seq(0.5,1.5,length.out=21),
           fn   =function(x) cbind(model.frame(params(x)),
                                   ll     =model.frame(x@ll),
                                   model.frame(refpts(x))[,-4],
                                   stock  =c(stock(  x)[,ac(range(x)['maxyear'])]%/%bmsy(x)),
                                   harvest=c(harvest(x)[,ac(range(x)['maxyear'])]%/%fmsy(x))),
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
    rsd=mdply(data.frame(iter=seq(dims(f)$iter)), function(iter){
        ft =fit(iter(fitted,iter),index)
        dgs=ft@diags[,c(".id","year","residual")]
        names(dgs)=c("name","year","residual")
            
        dgs})
          
    return(list(comp=rtn,residuals=rsd))}
        
  return(rtn)}

# # CIs
# cis <- max(surface) - qchisq(ci, 2)
# 
# do.call('contour', list(x=sort(profiled[[1]]), y=sort(profiled[[2]]), z=surface,
#                             levels=cis, add=TRUE, labcex=0.8, labels=ci))
# 
  
