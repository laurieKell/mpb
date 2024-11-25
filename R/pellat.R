
setMethod('pellat',  signature(object="FLPar"),  
          function(object){
            
            fn=function(x,object){
              res=optim(x,msy2pellatFn,refs=object,control=list(trace=0),
                        method="L-BFGS-B",
                        lower=c(0,0,1e-12),
                        upper=c(10,Inf,1))$par
              
              names(res)=c("r","k","p")
              res}
            
            res=aaply(object,seq(length(dim(object)))[-1],function(x) fn(c(0.5,object["bmsy"]*2,1),x))
            dmns=dimnames(par)
            dmns[[1]]=c("r","k","p")
            names(dmns)[1]="params"
            
            if (class(res)=="numeric")
              return(FLPar(res))
            
            if ((dims(par)$iter>1))
              prm=c(length(dim(par)),seq(length(dim(par)))[-length(dim(par))])
            else 
              prm=2:1
            
            res=as(array(aperm(res,prm),
                         dim=unlist(laply(dmns,length)),dmns),"FLPar")
            
            res})
