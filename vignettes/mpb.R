## ----knitr_init, echo=FALSE, results="hide"------------------------------
library(knitr)
## Global options
opts_chunk$set(echo    =TRUE,
               eval    =TRUE,
               cache   =FALSE,
               prompt  =FALSE,
               comment =NA,
               message =FALSE,
               tidy    =FALSE,
               warnings=FALSE,
               fig.height=4.5,
               fig.width =4.5,
               fig.path  ="tex/")

## ----echo=FALSE----------------------------------------------------------
sink(NULL)
warn=options()$warn
options(warn=-1)
library(ggplotFL)
library(plyr)
library(reshape)
library(mpb)

theme_set(theme_bw())
options(digits=3)
options(warn=warn)

sink()

## ----classCreate---------------------------------------------------------
bd=biodyn()

## ----classCoerce2,eval=FALSE,echo=TRUE-----------------------------------
#  asp=aspic("aspic.inp")
#  bd =as(asp,"biodyn")

## ----classSim,eval=FALSE-------------------------------------------------
#  bd=sim()

## ----plot, fig.margin=TRUE, fig.width=4, fig.height=5, fig.cap="Production function with simulated time series"----
bd=window(sim(),end=49)
plot(bd)+
  theme_bw()

## ----plotPrd, fig.margin=TRUE, fig.height=2.5, fig.cap="Simulated CPUE series"----
plotProduction(bd)+
  geom_path( aes(stock,catch),
             model.frame(FLQuants(bd,"stock","catch")))+
  geom_point(aes(stock,catch),
             model.frame(FLQuants(bd,"stock","catch")))+
  theme_bw()+theme(legend.position="none")

## ----fitU, fig.margin=TRUE, fig.cap="Simulated CPUE series"--------------
cpue=(stock(bd)[,-dims(bd)$year]+
      stock(bd)[,-1])/2
set.seed(7890)
cpue=rlnorm(1,log(cpue),.2)

ggplot(as.data.frame(cpue))+
  geom_point(aes(year,data))+
  geom_line(aes(year,data),col="salmon",
            data=as.data.frame(stock(bd)))+
  theme_bw()

## ----fitGuess2,size="tiny"-----------------------------------------------
params(bd)["k"]=mpb:::guessK(params(bd)["r"],mean(catch(bd),na.rm=T),params(bd)["p"])

## ----fitParams,fig.margin=TRUE,fig.width=4,fig.height=6------------------
setParams(bd)=cpue
params(bd)

## ----fitControl,fig.margin=TRUE,fig.width=4,fig.height=6-----------------
setControl(bd)=params(bd)
control(bd)

## ----fitrun,fig.margin=TRUE,fig.width=4,fig.height=6---------------------
control(bd)[c("p","b0"),"phase"]=-1
control(bd)["r","max"]=1.0
bdHat=fit(bd,cpue)

## ----fitcheck, eval=TRUE, fig.margin=TRUE, fig.height=6,fig.cap="A comparison of the true and fitted time series"----
params(bdHat)
params(bdHat)/params(bd)

plot(mpb:::biodyns(list("True"=bd,"Hat"=bdHat)))+
  theme(legend.position="bottom")+
  theme_bw()

## ----diag,echo=TRUE------------------------------------------------------
print(head(bdHat@diags),digits=3)

## ----diagQQ, fig.width=4,fig.height=4,fig.margin=TRUE, fig.cap="Quantile-quantile plot to compare residual distribution with the normal distribution."----
rsdl=bdHat@diags
ggplot(rsdl)                                           +
  geom_point( aes(qqx,qqy))                            +
  stat_smooth(aes(qqx,qqHat),method="lm",se=T,fill="blue", alpha=0.1) +
  theme_bw()+theme(legend.position="bottom")               

## ----diagHat, fig.margin=TRUE, fig.height=4, figwidth=4, fig.cap="Observed CPUE verses fitted, blue line is a linear resgression fitted to points, black the y=x line."----
library(diags)

ggplot(with(rsdl, data.frame(obs=stdz(obs),hat=stdz(hat))))   +
    geom_abline(aes(slope=1,intercept=0))                     +
    geom_point( aes(obs,hat))                                 +
    stat_smooth(aes(obs,hat),method="lm", se=F)               +
    theme_bw()+theme(legend.position="bottom")                +
    xlab("Fitted") + ylab("Observed")

## ----diagYr,fig.height=3, fig.margin=TRUE, fig.cap="Residuals by year, with lowess smoother"----
dat=transform(subset(rsdl,!is.na(residual), 
                     residual=stdz(residual,na.rm=T)))

ggplot(aes(year,residual),data=dat)  +
  geom_hline(aes(yintercept=0))      +
  geom_point()                       +
  stat_smooth(method="loess",se=F)   +
  theme_bw()+theme(legend.position="bottom")                

## ----diagVar,fig.height=3, fig.margin=TRUE, fig.cap="Plot of residuals against fitted value, to check variance relationship."----
ggplot(aes(hat, residual),
       data=subset(rsdl,!is.na(hat) & !is.na(residual)))   +
  geom_hline(aes(yintercept=0))         +
  geom_point()                          +
  stat_smooth(method="loess",se=F)      +
  theme_bw()+theme(legend.position="bottom")               

## ----diagAR, fig.width=4,fig.width=4,fig.margin=TRUE, fig.cap="Plot of autocorrelation, i.e. $residual_{t+1}$ verses $residual_{t}$."----
sum(rsdl$residual^2)

ggplot(rsdl)                                              +
  geom_point( aes(residual,residualLag))                  +
  stat_smooth(aes(residual,residualLag),method="lm",se=F) +
  geom_hline(aes(yintercept=0))     +
  xlab(expression(Residual[t]))     + 
  ylab(expression(Residual[t+1]))   +
  theme_bw()+theme(legend.position="bottom")                 

## ----prfl, fig.margin=TRUE, fig.height=3, fig.cap="Likelihood profile for r"----
bdHat=fit(bdHat,cpue)
setControl(bdHat)=params(bdHat)
res=profile(bdHat,which='r',fixed=c('b0','p'),
            cpue,range=seq(0.95,1.03,.002))
ggplot(subset(res,ll<0))+
  geom_line(aes(r,ll))  +
  theme_bw()

## ----prfl2, eval=FALSE, fig.margin=TRUE, fig.cap="Likelihood profile for r"----
#  res=profile(bdHat,which=c('r','k'),fixed=c('b0','p'),
#              cpue,range=seq(0.97,1.03,.02))
#  ggplot(res, aes(r, k, z=ll))+
#    stat_contour(aes(colour = ..level..), size = 1)+
#    theme_bw()

## ----prflLike, fig.margin=TRUE, fig.height=4, fig.width=4, fig.cap="Likelihood profile by data conmponent, i.e. CPUE series"----

bd=sim()

set.seed(7890)
Us  =FLQuants("Unbiased"     =
                rlnorm(1,log((stock(bd)[,-dims(bd)$year]+
                              stock(bd)[,-1])/2),0.2),
              "Increase in q"=
                rlnorm(1,log((stock(bd)[,-dims(bd)$year]+
                              stock(bd)[,-1])/2),0.2))

setParams( bd)=Us
setControl(bd)=params(bd)

bd@control[3:4,"phase"]=-1
bd=fit(bd,index=Us)
bd@control[,c("min")]=bd@params*0.1
bd@control[,c("val")]=bd@params
bd@control[,c("max")]=bd@params*10

prfl=profile(bd,which='r',index=Us,
             range=seq(0.975,1.05,.001))

ggplot(prfl)+
  geom_path(aes(r,ll,group=index,col=index))+
  facet_wrap(~index,scale="free",ncol=1)          +
  theme(legend.position="bottom")+
  theme_bw()

## ----prflADMB,echo=FALSE,eval=FALSE--------------------------------------
#  bd2=fit(bdHat,cpue,cmdOps="-lprof")
#  prf=subset(bd@profile, param %in% c("bbmsy","ffmsy"))
#  prf=data.frame(What="Profile",t(daply(prf, .(param), with, sample(value,500,prob=p,replace=T))))
#  names(prf)[2:3]=c("Stock","Harvest")

## ----uncertainty, fig.margin=TRUE----------------------------------------
bd   =window(sim(),end=39)
cpue=(stock(bd)[,-dims(bd)$year]+
      stock(bd)[,-1])/2
set.seed(7890)
cpue=rlnorm(1,log(cpue),.2)
bdHat=bd

setParams( bdHat)=cpue
setControl(bdHat)=params(bdHat)
bdHat@control[3:4,"phase"]=-1
bdHat=fit(bdHat,cpue)

sims=mpb:::biodyns(list("True"=bd,"Best Fit"=bdHat))

## ----uncertaintyCov,fig.height=6, fig.margin=TRUE------------------------
v=vcov(  bdHat)[c("r","k"),c("r","k"),1]
p=params(bdHat)[c("r","k")]
#refs=mvn(500,p,v)

## ----uncertaintyBoot, fig.height=4, fig.margin=TRUE, fig.cap="Bootstrapped CPUE series"----
set.seed(7890)
cpueBoot      =boot(bdHat)

sims[["Bootstrap"]]=fit(bdHat,cpueBoot)

## ----uncertaintyJackknife,fig.height=4,fig.margin=TRUE, fig.cap="Plot predicted stock trend by index"----
bdJK =fit(bdHat,FLQuant(jackknife(cpue)))

sims[["Jack Knife"]]=bdJK

## ----uncertaintyMCMC, eval=FALSE-----------------------------------------
#  sims[["MCMC"]]=fit(bdHat,cpue,cmdOps=c("-mcmc 1000000, -mcsave 5000"))

## ----uncertaintyMCMC2, eval=FALSE,fig.height=4,fig.margin=TRUE-----------
#  acf(c(params(sims[["MCMC"]])["r"]))

## ----fig.height=4,fig.margin=TRUE, fig.cap=""----------------------------
plot(sims[-4])+
  theme_bw()

## ----ref,eval=FALSE------------------------------------------------------
#  bdHat@mng

## ----ref2,eval=FALSE-----------------------------------------------------
#  bdHat@mngVcov

## ----ref3,eval=FALSE,fig.margin=TRUE,fig.width=4,fig.height=6,fig.cap=""----
#  currentState   =bdHat@mng[c("bbmsy","ffmsy"),"hat",drop=T]
#  currentStateVar=bdHat@mngVcov[c("bbmsy","ffmsy"),
#                                c("bbmsy","ffmsy"),drop=T]
#  
#  refs=mvrnorm(100,currentState,currentStateVar)
#  
#  ggplot(data=as.data.frame(refs))+
#    geom_histogram(aes(x=bbmsy))+
#    theme_bw()

## ----refBmsy,fig.margin=TRUE, fig.width=4, fig.height=3,fig.cap="Densities of Stock from different methods for estimating uncertainty."----

boot=stock(sims[["Bootstrap"]])[,39]

set.seed(7890) 
jack=mpb:::randJack(500,stock(sims[[  "Best Fit"]])[,39],
                  stock(sims[["Jack Knife"]])[,39])

bnow=rbind(data.frame(Method="boot",stock=c(boot)),
           data.frame(Method="jack",stock=c(jack)))
  
ggplot(bnow)+ 
  geom_density(aes(x=stock, y=..count..), position = "stack",fill="red")+
  facet_wrap(~Method,scale="free_y",ncol=1)+
  geom_vline(aes(xintercept=c(stock(sims[["Best Fit"]])[,"39"])))+
  theme_bw()

## ----kobe,eval=FALSE,fig.margin=TRUE,fig.width=4,fig.height=5,fig.caption="Kobe Phase Plots"----
#  library(kobe)
#  
#  kb=rbind(data.frame(Method="Boot",kobe(sims[["Bootstrap"]], what="pts")),
#           data.frame(Method="Jack",kobe(sims[["Jack Knife"]],what="pts")))
#  
#  ggplot(kb)+
#    geom_point(aes(stock,harvest))+
#    facet_wrap(~Method,scale="free_y",ncol=1)+
#    theme_bw()

## ----fdwd, fig.margin=TRUE,fig.width=4, fig.height=6,fig.cap="Projection"----
set.seed(7890)
harvest=rlnorm(100,log(harvest(bdHat))[,-dims(bdHat)$year],.1)
bdHat =fwd(bdHat,harvest=harvest)

plot(bdHat,worm=c(2,8))+    
  theme(legend.position="bottom")+
  theme_bw()

## ----hcr1----------------------------------------------------------------
bd=window(sim(),end=29)
for (i in seq(29,49,1))
  bd=fwd(bd,harvest=mpb:::hcr(bd,yr=i-1,hyr=i+1:2))
simHCR=mpb:::biodyns(list("Annual"=bd))

## ----hcr3----------------------------------------------------------------
bd=window(bd,end=29)
for (i in seq(29,49,3))
  bd=fwd(bd,harvest=mpb:::hcr(bd,yr=i,hyr=i+1:3))
simHCR[["Triennial"]]=bd

## ----hcrF----------------------------------------------------------------
bd=window(bd,end=29)
for (i in seq(29,49,3))
  bd=fwd(bd,harvest=mpb:::hcr(bd,yr=i,byr=i,hyr=i+1:3,bndF=c(0.9,1.1)))
simHCR[["bound F"]]=bd

## ----hcrY----------------------------------------------------------------
bd=window(bd,end=29)
for (i in seq(29,49,3))
  bd=fwd(bd,catch=mpb:::hcr(bd,yr=i,hyr=i+1:3,tac=TRUE,bndTac=c(0.9,1.1)))
simHCR[["bound TAC"]]=bd

## ----hcrPlot, fig.margin=TRUE,fig.width=4, fig.height=4,fig.cap="Plots of projections"----
plot(simHCR)+
  theme_bw()+
  theme(legend.position="bottom")

## ----MC,fig.margin=TRUE,fig.width=6,fig.height=6-------------------------
set.seed(7890)
pe=rlnorm(500,FLQuant(0,dimnames=list(year=1:50)),0.5)

bd=window(sim(),end=30)
bd.=bd
bd@stock =propagate(bd@stock, 500)
bd=fwd(bd,harvest=harvest(bd)[,2:30],pe=pe)

for (i in seq(30,48,1))
  bd=fwd(bd,
         catch=hcr(bd,yr=i,hyr=i+1,tac=TRUE,bndTac=c(0.9,1.1)),
         pe   =pe)

plot(bd)+
  theme_bw()

## ----MCkobe,fig.margin=TRUE,fig.width=6,fig.height=6---------------------
library(plyr)
library(mpb)
library(reshape)
library(kobe)
bd=sim()

#source('~/Desktop/flr/mpb/R/biodyn-msy.R')
#source('~/Desktop/flr/mpb/R/biodyn-hcr.R')

trks=kobe(bd,what="trks")
trks=mdply(data.frame(Year=seq(33,49,3)), 
           function(Year) subset(trks,year<=Year))

pts =mdply(data.frame(Year=seq(33,49,3)),function(Year)
                 kobe(bd,year=Year,what="pts"))

kobePhase()+    
    geom_line(aes(stock,harvest),data=plotHcr(bd),
              col="brown",size=1.5)                             +    
    geom_path( aes(stock,harvest),data=subset(trks,pctl=="50%"),col="blue")+
    geom_point(aes(stock,harvest),data=subset(pts,year>=33),size=.3,col="cyan")    +
    facet_wrap(~Year)

## ----mse,eval=FALSE------------------------------------------------------
#  mseBiodyn

