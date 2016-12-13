series=function(x){
  #pop numbers
  sn=FLQuants(r=stock.n(x)[1],
              j=apply(stock.n(x)%*%(1-mat(x)),c(2,6),sum),
              m=apply(stock.n(x)%*%mat(x),c(2,6),sum))
  
  #pop biomass
  sb=FLQuants(r=stock.n(x)[1]%*%stock.wt(x)[1],
              j=apply(stock.n(x)%*%(1-mat(x))%*%stock.wt(x),c(2,6),sum),
              m=apply(stock.n(x)%*%mat(x)%*%stock.wt(x),c(2,6),sum))
  
  #catch numbers
  cn=FLQuants(r=catch.n(x)[1],
              j=apply(catch.n(x)%*%(1-mat(x)),c(2,6),sum),
              m=apply(catch.n(x)%*%mat(x),c(2,6),sum))
  
  #catch biomass
  cb=FLQuants(r=catch.n(x)[1]%*%stock.wt(x)[1],
              j=apply(catch.n(x)%*%(1-mat(x))%*%catch.wt(x),c(2,6),sum),
              m=apply(catch.n(x)%*%mat(x)%*%catch.wt(x),c(2,6),sum))
  
  #pop wt
  sw=FLQuants(r=stock.wt(x)[1],
              j=apply(stock.n(x)%*%(1-mat(x))%*%stock.wt(x),c(2,6),sum)%/%
                apply(stock.n(x)%*%(1-mat(x)),c(2,6),sum),
              m=apply(stock.n(x)%*%mat(x)%*%stock.wt(x),c(2,6),sum)%/%
                apply(stock.n(x)%*%mat(x),c(2,6),sum))
  
  #catch wt
  cw=FLQuants(r=stock.wt(x)[1],
              j=apply(catch.n(x)%*%(1-mat(x))%*%catch.wt(x),c(2,6),sum)%/%
                apply(catch.n(x)%*%(1-mat(x)),c(2,6),sum),
              m=apply(catch.n(x)%*%mat(x)%*%catch.wt(x),c(2,6),sum)%/%
                apply(catch.n(x)%*%mat(x),c(2,6),sum))
  
  rbind(
    cbind(quantity="sn",model.frame(sn,drop=T)),
    cbind(quantity="sb",model.frame(sb,drop=T)),
    cbind(quantity="sw",model.frame(sw,drop=T)),
    cbind(quantity="cn",model.frame(cn,drop=T)),
    cbind(quantity="cb",model.frame(cb,drop=T)),
    cbind(quantity="cw",model.frame(cw,drop=T)))}
