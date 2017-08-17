library(FLCore)
library(plyr)
library(reshape)
library(ggplotFL)

dirMed="/home/laurie/ownCloud/GVYP_VPA/BFTAss/Analysis/vpa/east/reported/med"
file  =file.path(dirMed,"retro/bfte2014.d1")
cfile =file.path(dirMed,"retro/bfte2014.c1")

m  =c(0.49,0.24,0.24,0.24,0.24,0.20,0.175,0.15,0.125,0.10) 
vpa=readVPA2Box(cfile,m=m)
names(vpa)=seq(2014,2004,-1)

plot(vpa)
#lets pretend


