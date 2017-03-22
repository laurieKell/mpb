## ----knitr_init, echo=FALSE, results="hide"------------------------------
library(knitr)
## Global options
opts_chunk$set(echo    =!TRUE,
               eval    =TRUE,
               cache   =FALSE,
               cache.path="cache/",
               prompt  =FALSE,
               comment =NA,
               message =FALSE,
               tidy    =FALSE,
               warnings=FALSE,
               fig.height=4.5,
               fig.width =4.5,
               fig.path  ="tex/")

## ---- pkgs, echo=FALSE, message=FALSE------------------------------------
library(ggplot2)

theme_set(theme_bw())
options(digits=3)

## ----install,eval=FALSE--------------------------------------------------
## install.packages("mpb", repos = "http://flr-project.org/R")

## ----lib,echo=TRUE-------------------------------------------------------
library(ggplot2)
library(FLCore)
library(ggplotFL)
library(mpb)
library(FLife)
library(plyr)

## ----data,echo=TRUE------------------------------------------------------
data(ple4)

## ----plot,echo=TRUE------------------------------------------------------
plot(FLQuants(ple4,"Stock"=stock,"Index"=function(x) apply(oem(x),2,sum)))

## ----example1,echo=TRUE,eval=TRUE----------------------------------------
apply(oem(ple4),2,sum)

## ----plot-2,echo=TRUE----------------------------------------------------
ggplot(model.frame(mcf(FLQuants(ple4,"Stock"=stock,"Index"=function(x) apply(oem(x),2,sum)))))+
  geom_point( aes(Stock,Index))+
  geom_smooth(aes(Stock,Index),method="lm")+
  facet_null()

## ----example6,echo=TRUE--------------------------------------------------
timing=0.5
fish.dependent=FALSE
effort    =c("f","h")
mass      =TRUE

## ----example2,echo=TRUE--------------------------------------------------
cv=rlnorm(100,log(stock(ple4)),0.3)
ggplot(cv)+
  geom_boxplot(aes(factor(year),data))

## ----example4,echo=TRUE--------------------------------------------------
sel  =apply(harvest(ple4),1,mean)
ggplot(sel)+
  geom_line(aes(age,data))

## ----example3,echo=TRUE--------------------------------------------------
q    =FLQuant(cumprod(1+rep(.02,dim(fbar(ple4))[2])),dimnames=dimnames(fbar(ple4)))
plot(q)

## ----example5,echo=TRUE--------------------------------------------------
cpue   =stock(ple4)/mean(stock(ple4))
stable =cpue^0.1
deplete=cpue^2

plot(FLQuants("CPUE"          =cpue,
              `Hyper_Stable`  =stable,
              `Hyper_Depleted`=deplete))

## ----example7,echo=TRUE--------------------------------------------------
trend=FLQuant(seq(1,2,length.out=dim(stock(ple4))[2]),dimnames=dimnames(stock(ple4)))
var  =trend*abs(cv)*sign(cv)

ggplot(FLQuants("Log Normal" =cv,
                "Trend in CV"=var))+
  geom_boxplot(aes(factor(year),data))+
  facet_grid(qname~.)

## ----example8,echo=TRUE--------------------------------------------------
bias=FLPar(omega=1,ref=mean(stock(ple4)),q=0)

hyperstability<-function(object,omega=1,ref=apply(object,c(1,3:6),mean)) 
  ref%*%((object%/%ref)^omega)

bias<-function(object,bias=0.02) 
     FLQuant(cumprod(1+rep(bias,dim(object)[2])),dimnames=dimnames(object))


## ----example9,echo=TRUE--------------------------------------------------
set.seed(1234)
u     =FLQuants("Unbiased"      =rlnorm(100,log(apply(oem(ple4),2:6,sum)),.3),
                "Hyperstability"=rlnorm(100,log(apply(oem(ple4),2:6,sum)%*%
                                                  hyperstability(stock(ple4),0.52)),.3),
                "Trend"         =rlnorm(100,log(apply(oem(ple4),2:6,sum)%*%bias(stock(ple4),0.02)),.3),
                "AR"            =apply(oem(ple4),2:6,sum)%*%
                                   exp(rnoise(100,apply(oem(ple4),2:6,sum)*0,.3,b=.7)),
                "Variable"      =var,
                "Juvenile"      =rlnorm(100,log(apply(oem(ple4,sel=mat(ple4)),2:6,sum)),.3),
                "Mature"        =rlnorm(100,log(apply(oem(ple4,sel=1-mat(ple4)),2:6,sum)),.3),
                "Numbers"       =rlnorm(100,log(apply(oem(ple4,mass=FALSE),2:6,sum),.3)))

u=FLQuants(llply(u,function(x) x/mean(x)))
u=ldply(u,as.data.frame)

u.=ddply(u,.(year,.id), with, quantile(data))
ggplot()+
  geom_line(aes(year,data,col=factor(iter)),
            data=subset(u,iter%in%c(2,11)))+
  geom_ribbon(aes(year,ymin=`25%`,ymax=`75%`),data=u.,col="grey",alpha=.5)+
  facet_wrap(~.id,ncol=2)+
  theme_bw()+theme(legend.position="none")

## ---- devtools, echo=TRUE, eval=FALSE------------------------------------
## 	library(devtools)
## 	install_github('flr/FLPKG')

