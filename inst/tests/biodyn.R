library(mpb)
library(plyr)

dirTest="/home/laurie/Desktop/flr/mpb/inst/aspic"

runs=rbind(expand.grid(stock="albn",run=paste("run",1:7,sep="")),
           expand.grid(stock="bet", run=paste("run",c(3,5,6),sep="")))
asp=aspics(alply(runs,1,with, aspic(ac(file.path(dirTest,stock,run,"aspic.inp")))))
#fits=aspics(alply(runs,1,with, aspic(ac(file.path(dirTest,stock,run,"aspic.rdat")))))

attributes(runs)$split_labels=NULL

library(doParallel)
library(foreach)

cl=makeCluster(4)
registerDoParallel(cl)

asp=fit(asp)
plot(asp)

us=llply(asp,index,FALSE)
bd=llply(asp,as,"biodyn")

bd=llply(asp,function(x) as(x,"biodyn"))