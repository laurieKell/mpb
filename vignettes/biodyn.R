## ----initKnitr, echo=FALSE, results="asis"-------------------------------
library(knitr)

#output: rmarkdown::tufte_handout

# output:
#   rmdformats::html_clean:
#     fig_width: 6
#     fig_height: 6
#     highlight: pygments

## Global options
options(max.print="75")
opts_chunk$set(fig.path="out/",
               echo =TRUE,
               cache=!TRUE,
               cache.path="cache/",
               prompt=FALSE,
               tidy=FALSE,
               comment=NA,
               message=FALSE,
               warning=FALSE)

opts_knit$set(width=75)

## ----initPurl,eval=FALSE,echo=FALSE--------------------------------------
## library(knitr)
## purl(    "/home/laurie/Desktop/flr/git/biodyn/stuff/vignettes/biodyn.Rmd",
##      out="/home/laurie/Desktop/flr/git/biodyn/vignettes/biodyn.R")
## 

## ----initLib-------------------------------------------------------------
library(biodyn)
library(ggplot2)
library(plyr)
library(diags)
library(kobe)
library(ggplotFL)

## ----classCreate---------------------------------------------------------
bd =biodyn()

## ----classCreate2--------------------------------------------------------
bd=biodyn(catch=FLQuant(100,dimnames=list(year=1990:2010)))

## ----classCoerce,eval=FALSE----------------------------------------------
## data(ple4)
## bd =as(ple4,"biodyn")

## ----classCoerce2,eval=FALSE---------------------------------------------
## library(aspic)
## asp=aspic("/home/laurie/Desktop/rfmos/iccat/kobe/Inputs/swon/2013/aspic/run2/aspic.inp")
## bd =as(asp,"biodyn")

## ----classSim,eval=FALSE-------------------------------------------------
## bd=sim()

## ----plot, fig.margin=TRUE, fig.width=4, fig.height=8, fig.cap="Simulated time series"----
bd=sim()
bd=window(bd,end=49)
plot(bd)

## ----plotEql, fig.margin=TRUE, eval=FALSE, fig.height=6, fig.width=4, echo=FALSE, fig.cap="Simulated CPUE series"----
## plotEql(bd)

## ----plotPrd, fig.margin=TRUE, fig.cap="Simulated CPUE series"-----------
library(reshape)
x=sim()
plotPrd(x)+
  geom_path( aes(stock,catch),
             model.frame(FLQuants(x,"stock","catch")))+
  geom_point(aes(stock,catch),
             model.frame(FLQuants(x,"stock","catch")))

## ----plotMSE, fig.margin=TRUE, eval=FALSE, fig.cap="Simulated CPUE series"----
## plotMSE()

## ----plotMC, fig.margin=TRUE, fig.height=6, fig.width=4,fig.cap="Monte Carlo simuation time series with confidence intervals and a worm (single simulation)."----
harvest=rlnorm(200,log(harvest(bd)[,-1]),.2)

bd=fwd(bd,harvest=harvest)

plot(bd,worm=3)+
  theme(legend.position="bottom")

## ----fit, echo=TRUE, fig.margin=TRUE, fig.height=4, fig.cap="Simulated stock"----
bd=sim()

## ----fitU, fig.margin=TRUE, fig.cap="Simulated CPUE series"--------------
cpue=(stock(bd)[,-dims(bd)$year]+
      stock(bd)[,-1])/2
cpue=rlnorm(1,log(cpue),.2)

ggplot(as.data.frame(cpue))+
  geom_point(aes(year,data))+
  geom_line(aes(year,data),data=as.data.frame(stock(bd)),col="salmon")

## ----fitGuess,eval=FALSE-------------------------------------------------
## bd=biodyn(catch=catch(bd))

## ----fitParams,fig.margin=TRUE,fig.width=4,fig.height=6------------------
setParams( bd)=cpue
params(bd)

## ----fitControl,fig.margin=TRUE,fig.width=4,fig.height=6-----------------
setControl(bd)=params(bd)
  
bd@control

## ----fitRun, eval=TRUE, fig.margin=TRUE, fig.height=6,fig.cap="A comparison of the true and fitted time series"----
bd@control[3:4,"phase"]=-1
bdHat=fit(bd,cpue)

# plot(biodyns("True"=bd,"Hat"=bdHat))+
#   theme(legend.position="bottom")

## ----fitCheck,fig.margin=TRUE,fig.width=4,fig.height=6-------------------
params(bdHat)
params(bdHat)/params(bd)

## ----diag,echo=TRUE------------------------------------------------------
head(bdHat@diags)

## ----diagQQ, fig.width=4,fig.height=4,fig.margin=TRUE, fig.cap="Quantile-quantile plot to compare residual distribution with the normal distribution."----
rsdl=bdHat@diags
ggplot(rsdl)                                           +
  geom_point( aes(qqx,qqy))                            +
  stat_smooth(aes(qqx,qqHat),method="lm",se=T,fill="blue", alpha=0.1)         +
  theme_ms(14,legend.position="bottom")               

## ----diagHat, fig.margin=TRUE, fig.height=4, figwidth=4, fig.cap="Observed CPUE verses fitted, blue line is a linear resgression fitted to points, black the y=x line."----
ggplot(with(rsdl, data.frame(obs=stdz(index),hat=stdz(hat)))) +
    geom_abline(aes(0,1))                                     +
    geom_point( aes(obs,hat))                                 +
    stat_smooth(aes(obs,hat),method="lm", se=F)               +
    theme_ms(14,legend.position="bottom")                     +
    xlab("Fitted") + ylab("Observed")

## ----diagYr,fig.height=3, fig.margin=TRUE, fig.cap="Residuals by year, with lowess smoother"----
dat=transform(subset(rsdl,!is.na(residual), 
                     residual=stdz(residual,na.rm=T)))

ggplot(aes(year,residual),data=dat)  +
  geom_hline(aes(yintercept=0))      +
  geom_point()                       +
  stat_smooth(method="loess",se=F)   +
  theme_ms(14,legend.position="bottom") 

## ----diagVar,fig.height=3, fig.margin=TRUE, fig.cap="Plot of residuals against fitted value, to check variance relationship."----
ggplot(aes(hat, residual),
       data=subset(rsdl,!is.na(hat) & !is.na(residual)))   +
  geom_hline(aes(yintercept=0))         +
  geom_point()                          +
  stat_smooth(method="loess",se=F)      +
  theme_ms(14,legend.position="bottom")

## ----diagAR, fig.width=4,fig.width=4,fig.margin=TRUE, fig.cap="Plot of autocorrelation, i.e. $residual_{t+1}$ verses $residual_{t}$."----
ggplot(rsdl)                                              +
  geom_point( aes(residual,residualLag))                  +
  stat_smooth(aes(residual,residualLag),method="lm",se=F) +
  geom_hline(aes(yintercept=0))     +
  xlab(expression(Residual[t]))     + 
  ylab(expression(Residual[t+1]))   +
  theme_ms(14,legend.position="bottom")  

## ----prfl, fig.margin=TRUE, fig.cap="Likelihood profile for r"-----------
bdHat=fit(bdHat,cpue)
setControl(bdHat)=params(bdHat)
res=profile(bdHat,which='r',fixed=c('b0','p'),
            cpue,range=seq(0.95,1.03,.002))
ggplot(subset(res,ll.u1<0))+geom_line(aes(r,ll.u1))

## ----prfl2, eval=FALSE, fig.margin=TRUE, fig.cap="Likelihood profile for r"----
## res=profile(bdHat,which=c('r','k'),fixed=c('b0','p'),
##             cpue,range=seq(0.97,1.03,.02))
## ggplot(res, aes(r, k, z=ll.u1))+
##   stat_contour(aes(colour = ..level..), size = 1)

## ----prflLike, fig.margin=TRUE, fig.height=6, fig.width=4, fig.cap="Likelihood profile by data conmponent, i.e. CPUE series"----

bd=sim()

Us  =FLQuants("Unbiased"     =
                rlnorm(1,log((stock(bd)[,-dims(bd)$year]+
                              stock(bd)[,-1])/2),0.2),
              "Increase in q"=
                rlnorm(1,log((stock(bd)[,-dims(bd)$year]+
                              stock(bd)[,-1])/2),0.2))

bds=bd

setParams( bds)=Us
setControl(bds)=params(bds)

bds@control[3:4,"phase"]=-1
bds=fit(bds,index=Us)
bds@control[,c("min")]=bds@params*0.1
bds@control[,c("val")]=bds@params
bds@control[,c("max")]=bds@params*10

fn=function(x) cbind(model.frame(params(x)["r"]),
                     ll=model.frame(x@ll)[,-3],
                     ll=apply(x@ll,2,sum))

prfl=profile(bds,which='r',index=Us,
             range=seq(0.70,1.05,.001),fn=fn)

ggplot(subset(melt(prfl[,c("r","ll.u1","ll.u2","1")],id="r"),
              value<10))+
  geom_line(aes(r,value,group=variable,col=variable))+
  facet_wrap(~variable,scale="free",ncol=1)          +
  theme(legend.position="bottom")

## ----prflADMB,echo=FALSE,eval=FALSE--------------------------------------
## bd2=fit(bdHat,cpue,cmdOps="-lprof")
## prf=subset(bd@profile, param %in% c("bbmsy","ffmsy"))
## prf=data.frame(What="Profile",t(daply(prf, .(param), with, sample(value,500,prob=p,replace=T))))
## names(prf)[2:3]=c("Stock","Harvest")

## ----uncertainty, fig.margin=TRUE----------------------------------------
bd   =window(sim(),end=39)
cpue=(stock(bd)[,-dims(bd)$year]+
      stock(bd)[,-1])/2
cpue=rlnorm(1,log(cpue),.2)
bdHat=bd

setParams( bdHat)=cpue
setControl(bdHat)=params(bdHat)
bdHat@control[3:4,"phase"]=-1
bdHat=fit(bdHat,cpue)

sims=biodyns("True"=bd,"Best Fit"=bdHat)

## ----uncertaintyCov,fig.height=6, fig.margin=TRUE------------------------
vcov(bdHat)

## ----uncertaintyBoot, fig.height=4, fig.margin=TRUE, fig.cap="Bootstrapped CPUE series",eval=FALSE----
## cpueSim=bdHat@diags[,c("year","hat")]
## names(cpueSim)[2]="data"
## cpueSim=as.FLQuant(cpueSim)
## 
## cv=sd(sims[["Best Fit"]]@diags[,"residual"])
## 
## cpueSim=rlnorm(100,log(cpueSim),cv)
## 
## cpueSim[cpueSim==0]=NA
## 
## plot(cpueSim,na.rm=TRUE)
## 
## sims[["CPUE"]]=fit(propagate(bdHat,100),cpueSim)

## ----uncertaintyJackknife,fig.height=4,fig.margin=TRUE, fig.cap="Plot predicted stock trend by index"----

bdJK =fit(bdHat,jackknife(cpue))

sims[["Jack Knife"]]=bdJK

#plotJack(bdJK)

## ----uncertaintyMCMC-----------------------------------------------------
sims[["MCMC"]]=fit(bdHat,cpue,cmdOps=c("-mcmc 1000000, -mcsave 5000"))

## ----uncertaintyMCMC2,fig.height=4,fig.margin=TRUE-----------------------
acf(c(params(sims[["MCMC"]])["r"]))

## ----uncertaintySave,eval=FALSE,echo=FALSE-------------------------------
## save(sims,cpue,bdJK,bdHat,file="/home/laurie/Desktop/flr/git/biodyn/stuff/data/sims.RData")

## ----ref,eval=FALSE------------------------------------------------------
## bdHat@mng

## ----ref2,eval=FALSE-----------------------------------------------------
## bdHat@mngVcov

## ----ref3,eval=FALSE,fig.margin=TRUE,fig.width=4,fig.height=6,fig.cap=""----
## currentState   =bdHat@mng[c("bbmsy","ffmsy"),"hat",drop=T]
## currentStateVar=bdHat@mngVcov[c("bbmsy","ffmsy"),
##                               c("bbmsy","ffmsy"),drop=T]
## 
## refs=mvrnorm(100,currentState,currentStateVar)
## 
## ggplot(data=as.data.frame(refs))+
##    geom_histogram(aes(x=bbmsy))

## ----refBmsy,fig.fullwidth=TRUE, fig.width=4, fig.height=6,fig.cap="Densities of Stock from different methods for estimating uncertainty.",eval=FALSE----
## #load("/home/laurie/Desktop/flr/git/biodyn/stuff/data/sims.RData")
## sims[["Jack Knife"]]=bdJK
## 
## c("Best Fit","CPUE","Jack Knife","MCMC")
## 
## boot=stock(sims[["CPUE"]])[,39]
## # ggplot(as.data.frame(boot))+
## #   geom_density(aes(x=data, y=..count..), position = "stack",fill="red")+
## #   scale_x_continuous(limits=c(0,700))
## 
## mcmc=stock(sims[["MCMC"]])[,39]
## # ggplot(as.data.frame(mcmc))+
## #   geom_density(aes(x=data, y=..count..), position = "stack",fill="red")+
## #   scale_x_continuous(limits=c(0,700))
## 
## vcov=rnorm(500,sims[["Best Fit"]]@mng["bnow","hat"],
##                sims[["Best Fit"]]@mng["bnow", "sd"])
## # ggplot(as.data.frame(vcov))+
## #   geom_density(aes(x=data, y=..count..), position = "stack",fill="red")+
## #   scale_x_continuous(limits=c(0,1000))
## 
## jack=randJack(500,stock(sims[[  "Best Fit"]])[,39],
##                   stock(sims[["Jack Knife"]])[,39])
## # ggplot(as.data.frame(jack))+
## #   geom_density(aes(x=data, y=..count..), position = "stack",fill="red")+
## #   scale_x_continuous(limits=c(0,700))
## 
## bnow=rbind(data.frame(Method="boot",stock=c(boot)),
##            data.frame(Method="mcmc",stock=c(mcmc)),
##            data.frame(Method="vcov",stock=c(vcov)),
##            data.frame(Method="jack",stock=c(jack)))
## 
## ggplot(bnow)+
##   geom_density(aes(x=stock, y=..count..), position = "stack",fill="red")+
##   facet_wrap(~Method,scale="free_y",ncol=1)+
##   geom_vline(aes(xintercept=c(stock(sims[["Best Fit"]])[,"39"])))

## ----refFMSY,eval=FALSE,fig.fullwidth=TRUE, fig.width=6, fig.height=8,fig.cap="Densities of Stock/BMSY from different methods for estimating uncertainty."----
## boot=boot%/%bmsy(sims[["CPUE"]])
## mcmc=mcmc%/%bmsy(sims[["MCMC"]])
## vcov=rnorm(500,sims[["Best Fit"]]@mng["bbmsy","hat"],
##                sims[["Best Fit"]]@mng["bbmsy", "sd"])
## jack=randJack(500,stock(sims[[  "Best Fit"]])[,39]%/%bmsy(sims[[  "Best Fit"]]),
##                   stock(sims[["Jack Knife"]])[,39]%/%bmsy(sims[["Jack Knife"]]))
## 
## bbmsy=rbind(data.frame(Method="boot",stock=c(boot)),
##             data.frame(Method="mcmc",stock=c(mcmc)),
##             data.frame(Method="vcov",stock=c(vcov)),
##             data.frame(Method="jack",stock=c(jack)))
## 
## ggplot(bnow)+
##   geom_density(aes(x=stock, y=..count..), position = "stack",fill="red")+
##   facet_wrap(~Method,scale="free_y",ncol=1)+
##   geom_vline(aes(xintercept=c(stock(sims[["Best Fit"]])[,"39"])))+
##   scale_x_continuous(limits=c(0,2.0))

## ----kobe,fig.margin=TRUE,fig.width=4,fig.height=4,fig.caption="Kobe Phase Plots",eval=FALSE----
## library(kobe)
## 
## kb=rbind(data.frame(Method="Boot",kobe(sims[["CPUE"]],      what="pts")),
##          data.frame(Method="MCMC",kobe(sims[["MCMC"]],      what="pts")),
##          data.frame(Method="Vcov",kobe(sims[["Best Fit"]],  what="pts")),
##          data.frame(Method="Jack",kobe(sims[["Jack Knife"]],what="pts")))
## 
## ggplot(kb)+
##   geom_point(aes(stock,harvest),data=subset(df,year==39))+
##   facet_wrap(~Method,scale="free_y",ncol=1)

## ----uncertaintyJackknife2,,eval=FALSE,fig.height=4,fig.margin=TRUE,eval=FALSE----
## source('~/Desktop/flr/git/biodyn/R/biodyn-jackRand.R')
## source('~/Desktop/flr/git/biodyn/R/biodyn-jackSummary.R')
## 
## sims[["Jack Knife"]]=randJack(500,bdHat,bdJK)

## ----uncertaintyCov2,fig.height=6, fig.margin=TRUE,eval=FALSE------------
## sims[["Vcov"]]=bdHat
## sims[["Vcov"]]=mvn(bdHat,500,nms=c("r","k"),fwd=TRUE)

## ----fdwd, fig.margin=TRUE,fig.width=4, fig.height=6,fig.cap="Projection"----
harvest=rlnorm(100,log(harvest(bdHat))[,-dims(bdHat)$year],.1)
bdHat =fwd(bdHat,harvest=harvest)

plot(bdHat,worm=c(2,8))+    
  theme(legend.position="bottom")

## ----hcr1----------------------------------------------------------------
bd   =sim()

bd=window(bd,end=29)
for (i in seq(29,49,1))
  bd=fwd(bd,harvest=hcr(bd,refYrs=i,yrs=i+1)$hvt)
simHCR=biodyns("1"=bd)

## ----hcr3----------------------------------------------------------------
bd=window(bd,end=29)
for (i in seq(29,49,3))
  bd=fwd(bd,harvest=hcr(bd,refYrs=i,yrs=i+1:3)$hvt)
simHCR[["3"]]=bd

## ----hcrF----------------------------------------------------------------
bd=window(bd,end=29)
for (i in seq(29,49,1))
  bd=fwd(bd,harvest=hcr(bd,refYrs=i,yrs=i+1,bndF=c(0.9,1.1))$hvt)
simHCR[["bound F"]]=bd

## ----hcrY----------------------------------------------------------------
bd=window(bd,end=30)
for (i in seq(29,49,1))
  bd=fwd(bd,catch=hcr(bd,refYrs=i,yrs=i+1,tac=T,bndTac=c(0.9,1.1))$tac)
simHCR[["bound TAC"]]=bd

## ----hcrPlot,fig.fullwidth=TRUE,fig.width=6,fig.height=6,fig.cap="Plots of projections"----
plot(simHCR)+
  theme(legend.position="bottom")


## ----MC,fig.margin=TRUE,fig.width=6,fig.height=6-------------------------
pe=rlnorm(500,FLQuant(0,dimnames=list(year=1:50)),0.5)

bd=window(sim(),end=30)
bd.=bd
bd@stock =propagate(bd@stock, 500)
bd=fwd(bd,harvest=harvest(bd)[,2:30],pe=pe)

for (i in seq(30,48,1))
  bd=fwd(bd,
         catch=hcr(bd,refYrs=i,yrs=i+1,tac=T,bndTac=c(0.9,1.1))$tac,
         pe   =pe)

plot(bd)

## ----MCkobe,eval=FALSE,fig.margin=TRUE,fig.width=6,fig.height=6----------
## trks=biodyn::kobe(bd,what="trks")
## trks=mdply(data.frame(Year=seq(33,49,3)),
##            function(Year) subset(trks,year<=Year))
## 
## pts =transform(biodyn::kobe(bd,what="pts",year=seq(33,49,3)),
##                  Year=year)[,c("stock","harvest","Year")]
## 
## kobePhase()+
##     geom_line(aes(stock,harvest),data=hcrPlot(bd.),
##               col="brown",size=1.5)                             +
##     geom_path( aes(stock,harvest),data=subset(trks,pctl=="50%"))+
##     geom_point(aes(stock,harvest),data=subset(pts,Year>=33))    +
##     facet_wrap(~Year)

## ----mse,eval=FALSE------------------------------------------------------
## mseBiodyn

## ----,eval=FALSE, results='asis'-----------------------------------------
## library(xtable)
## options(xtable.comment = FALSE)
## options(xtable.booktabs = TRUE)
## xtable(head(mtcars[,1:6]), caption = "First rows of mtcars")

