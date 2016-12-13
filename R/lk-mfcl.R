utils::globalVariables(c('par','polygon','lines','segments','box','text','ploc','axis','dev.copy'))
utils::globalVariables('X2')
utils::globalVariables('X1')
utils::globalVariables('as.FLPar')
utils::globalVariables('par')
utils::globalVariables('polygon')
utils::globalVariables('lines')
utils::globalVariables('segments')
utils::globalVariables('box')
utils::globalVariables('text')
utils::globalVariables('ploc')
utils::globalVariables('axis')
utils::globalVariables('dev.copy')
utils::globalVariables('getfrqfile')
utils::globalVariables('getoutputparfile')
utils::globalVariables('getfrqsampledat')

getfrqheader <-function (frqfile = getfrqfile("plot.rep")) 
  {
    dat <- scan(frqfile, comment.char="#", quiet = TRUE)
    x <- numeric(12)
    names(x) <- c("nreg", "nflt", "year1", "nevents", "nlengths", 
                  "startlens", "lenbinsize", "lenfactor",
                  "nweights", "startwts", "wtbinsize", "wtfactor")
    x["nreg"] <- dat[1]
    x["nflt"] <- dat[2]
    x["year1"] <- dat[5]
    version <- dat[10]
    tnma <- 10 + x["nreg"] + x["nflt"] + x["nreg"]*(x["nreg"] - 1)/2 +1
    tnm <- tnma + 5*x["nflt"]
    fdatflags <- matrix(dat[tnma:(tnm-1)], nrow=5, byrow=TRUE)
    if (version == 4) {
      tnm <- tnm + dat[tnm] +1
    }
    else if (version > 4) {
      stop("DIE_YUPPIE_SCUM")
    }
    x["nevents"]    <- dat[tnm]
    x["nlengths"]   <- dat[tnm + 1]
    x["startlens"]  <- dat[tnm + 2]
    x["lenbinsize"] <- dat[tnm + 3]
    x["lenfactor"]  <- dat[tnm + 4]
    if(version >= 4) {
      x["nweights"]   <- dat[tnm + 5]
      x["startwts"]   <- dat[tnm + 6]
      x["wtbinsize"]  <- dat[tnm + 7]
      x["wtfactor"]   <- dat[tnm + 8]
    } else {
      x["nweights"] <- x["startwts"] <- x["wtbinsize"] <-
        x["wtbinsize"] <- x["wtfactor"] <- NA
    }
    x <- as.list(x)
    x$regionsize <- dat[11:(10 + dat[1])]
    x$fltregion <- dat[(11 + dat[1]):(10 + dat[1] + dat[2])]
    x$fdatflags <- fdatflags
    x
  }


freqplot=function (length.fit = "length.fit", path, flt = 0, xlab = "", 
                   ylab = "", sampids, fltcodes, lwd = 2, hcol = "yellow", lcol = "red", 
                   titleloc = "in", nplotcols=1, nplotrows=30, plot=FALSE, ...) 
{
  if (!missing(path)) 
    length.fit <- paste(path, "/", length.fit, sep = "")
  if (is.character(flt) | flt <= 0) 
    flt <- "\"totals\""
  ## 2005/06/27 Y.T.
  qq <- if(.Platform$OS.type!="unix") '\" '
  else "' "
  
  
  awk <- paste("awk ",qq,"(NR==2){print;", "getline;print;", "getline;print;", 
               "getline;print}", "(/^# fishery/ && $3==", flt, "){read=1;next}", 
               "(/^# fishery/ && read){exit}", "(read){print}",qq, length.fit)
  dat <- scan(pipe(awk), quiet = TRUE)
  nbins <- dat[1]
  bin1 <- dat[2]
  binwid <- dat[3]
  nflt1 <- dat[4]
  if (missing(fltcodes)) 
    fltcodes <- 1:nflt1
  if (is.character(flt)) 
    flt <- nflt1
  nsamp <- dat[4 + flt]
  if (nsamp <= 0) 
    stop(paste("No size samples for fleet ", flt), call. = F)
  nages <- dat[5 + nflt1]
  xp <- seq(bin1, by = binwid, len = nbins + 1)
  xl <- seq(bin1 + binwid/2, by = binwid, len = nbins)
  dpt <- 4 + nages + (2 + nages) * nbins
  id <- character(nsamp)
  lenatage <- matrix(nrow = nages, ncol = nsamp)
  ofrq <- matrix(nrow = nbins, ncol = nsamp)
  tfrq <- matrix(nrow = nbins, ncol = nsamp)
  frq <- array(dim = c(nages, nbins, nsamp))
  for (i in 1:nsamp) {
    offset <- 5 + nflt1 + (i - 1) * dpt
    id[i] <- paste(dat[offset + (1:3)], collapse = "|")
    lenatage[, i] <- dat[offset + 3 + (1:nages)] * binwid + 
      bin1
    ofrq[, i] <- dat[offset + 4 + nages + (1:nbins)]
    tfrq[, i] <- dat[offset + 4 + nages + nbins + (1:nbins)]
    frq[, , i] <- matrix(dat[offset + 4 + nages + nbins * 
                               2 + (1:(nages * nbins))], byrow = T, ncol = nbins)
  }
  if (!missing(sampids)) {
    k <- match(sampids, id)
  }
  else k <- 1:nsamp
  
  nplt <- length(k)
  if (missing(nplotcols)) nplotcols <- min(5,floor(sqrt(nplt)))
  nplotrows <- min(ceiling(nplt/nplotcols),nplotrows)
  
  res=list(hat=tfrq,      #fitted cas   
           obs=ofrq,      #observed cas
           laa=lenatage,  #laa
           caa=frq,       #caa
           id =id)    
  
  if (!plot) return(res)
  
  
  if(nplt>1) {
    opar <- par(mfcol = c(nplotrows+1, nplotcols),
                mar = c(3/(nplotrows+1), 2, 0, 1) + 0.1, xpd = T)
    on.exit(par(opar))
  } else {
    nplotrows <- nplotcols <- 1
  }
  
  for (jj in seq(along=k)) {
    j <- k[jj]
    ymax <- max(ofrq[, j], tfrq[, j], frq[, , j])
    if (ymax > 0) {
      plot(range(xp), range(0, ymax), type = "n", axes = F, ann=F) 
      polygon(rep(xp, rep(2, nbins + 1)), c(0, rep(ofrq[,j], rep(2, nbins)), 0), col = hcol, border=hcol,lwd = 0.1)
      lines(xl, tfrq[, j], type = "l", col = "red", lwd = lwd)
      if (flt < nflt1) {
        for (i in 1:nages) lines(xl, frq[i, , j], type = "l", col = lcol)
        segments(lenatage[, j], rep(0, nages), lenatage[,j], rep(ymax, nages), lty = 3)
      }
      box(bty = "l", lwd = 0.1)
      if (titleloc == "in") {
        title <- paste(ifelse(flt < nflt1, paste(fltcodes[flt], 
                                                 ", ", id[j], sep = ""), paste(fltcodes[j])))
        usr <- par("usr")
        text(ploc(0.98, 0.85), title, adj = 1, cex = 1, 
             ...)
      }
      else {
        title(if (flt < nflt1) 
          paste(fltcodes[flt], ", ", id[j], sep = "")
              else paste(fltcodes[j]), ...)
      }
      if ( jj%%nplotrows == 0 | jj == nplt ) {
        axis(1, lwd = 0.1)
        if(nplt>1) plot(c(0,1),c(0,1),type='n',axes=F,ann=F)
      } else {
        axis(1, labels = F, lwd = 0.1)
      }
      if(jj%%(nplotrows*nplotcols)==0 & j<nplt) {
        if(readline("Save plot? (y/n): ") == "y") {
          dev.copy(eval(parse(text=options("device"))))
        }
      }
      #cat(paste("\"", id[j], "\" ", sep = ""))
    } else plot(c(0,1),c(0,1),type='n',axes=F,ann=F)
  }
  
  #cat("\n")
  
  return(res)}


mfclLfd=function(ffile){
  ln =with(getfrqheader(ffile),startlens+cumsum(c(0,rep(lenbinsize,nlengths-1))))
  frq=getfrqsampledat(ffile)
  frq=melt(frq,id=c("yr","mo","wk","flt"))
  frq=transform(frq,params=ln[variable],data=value)[,c(7,1,2,4,8)]
  frq=frq[frq$data>0,]
  names(frq)[c(2:4)]=c("year","season","unit")
  
  #frq=as(frq,"FLQuant")
  
  return(frq)}

#' Returns length age data from Multifan-CL as a list with four data frames
#' 
#' @param lfile the \code{length.fit} file
#' @param ffile the \code{freq} file
#' @param df return data as a data.frame \code{logical} default is \code{TRUE}
#' @param i fleet to read in
#' 
#' @return a list with elements \code{caa} numbers by age, len, season and year, \code{laa} length by age season and year, 
#' \code{obs} observed length samples by season and year and \code{hat} fitted lengths by season and year,
#' 
alkMFCL=function(lfile,ffile,i=1,df=TRUE){

  #' lfile="/home/laurie/Desktop/gcode/mse4mfcl/ALB/papers/SCRS/SCRS2013-ALK/Inputs/length09.fit"
  #' rfile="/home/laurie/Desktop/gcode/mse4mfcl/ALB/papers/CSRS/SCRS2013-ALK/Inputs/plot-09.par.rep"
  #' ffile="/home/laurie/Desktop/Dropbox/collaboration/Shelton/ALBN/4B/2009/albN.frq"
  
  wrn=options()$warn
  options(warn=-1)
  
  flt=freqplot(lfile,flt=i)
  len=with(getfrqheader(ffile),startlens+cumsum(c(0,rep(lenbinsize,nlengths-1))))
  
  ## catch-at-age
  age   =rep(1:15,prod(dim(flt$caa)[-1]))
  year  =rep(laply(strsplit(flt$id,"\\|"),function(x) x[1]),each=prod(dim(flt$caa)[-3]))
  season=rep(laply(strsplit(flt$id,"\\|"),function(x) x[2]),each=prod(dim(flt$caa)[-3]))
  
  caa=transform(melt(flt$caa),len   =len[X2],
                              age   =seq(dim(flt$caa)[1])[X1],
                              year  =laply(strsplit(flt$id,"\\|"),function(x) x[1]),
                              season=laply(strsplit(flt$id,"\\|"),function(x) x[2]))[,-(1:3)]
  names(caa)[1]="data"                           
  #ggplot(caa)+geom_histogram(aes(len,weight=freq),binwidth=2)+facet_wrap(~year)
  #ggplot(caa)+geom_histogram(aes(age,weight=freq),binwidth=1)+facet_wrap(~year)
  
  if (!df) {
    names(caa)[4]="params"
    caa=as(caa,"FLPar")
    }
  
  
  ## length-at-age
  age   =rep(1:15,prod(dim(flt$laa)[-1]))
  year  =rep(laply(strsplit(flt$id,"\\|"),function(x) x[1]),each=prod(dim(flt$laa)[-2]))
  season=rep(laply(strsplit(flt$id,"\\|"),function(x) x[2]),each=prod(dim(flt$laa)[-2]))
  laa   =data.frame(age=age,year=year,season=season,data=c(flt$laa),iter=1)
  if (!df) {
    #laa=as.FLQuant(laa)
  }
  
  ## observed size sample
  len   =with(getfrqheader(ffile),startlens+cumsum(c(0,rep(lenbinsize,nlengths-1))))
  year  =rep(laply(strsplit(flt$id,"\\|"),function(x) x[1]),prod(dim(flt$obs)[-2]))
  season=rep(laply(strsplit(flt$id,"\\|"),function(x) x[2]),prod(dim(flt$obs)[-2]))
  obs   =data.frame(len=len,year=year,season=season,data=c(flt$obs),iter=1)
  if (!df) {
    names(obs)[1]="params"
    obs=as(obs,"FLPar")
  }
  
  ## fitted size sample
  hat   =data.frame(len=len,year=year,season=season,data=c(flt$hat),iter=1)
  if (!df) {
    names(hat)[1]="params"
    hat=as.FLPar(hat)
  }
  
  options(warn=wrn)
  
  if (df){
    caa=caa[caa$data>0 & !is.na(caa$data),c("year","age","season","len","data")]
    laa=laa[laa$data>0 & !is.na(laa$data),c("year",      "season",      "data")]
    obs=obs[obs$data>0 & !is.na(obs$data),c("year",      "season","len","data")]
    hat=hat[hat$data>0 & !is.na(hat$data),c("year",      "season","len","data")]
    }
  
  return(list(caa=caa,laa=laa,obs=obs,hat=hat))}

#################################################################################
# MFCL growth curve Get function
#################################################################################

#' Returns von Bertalannfy growth parameters from Multifan-CL
#' 
#' @param parfile the par file output from Multifan-cl
#' 
#' @return a numeric vector
#' 
#' @export
#' 
mfclGrw=function(parfile = getoutputparfile("plot.rep")){
  #==============================================================
  # Get avg size in 1st and last age classes plus the Brody rho
  #==============================================================
  tt <- getplotdat4("# The von Bertalanffy parameters",parfile)[c(1,4,7)]
  tt <- list(size1=tt[1],size2=tt[2],rho=tt[3])
  dt <- getplotdat1("# The number of age classes",parfile)-1
  tt$Linf <- with(tt,(size2-size1*exp(-rho*dt))/(1-exp(-rho*dt)))
  tt$t0 <- with(tt,1-log(Linf/(Linf-size1))/rho)
  tt=unlist(tt)

  names(tt)[3:5]=c("k","linf","t0")
  tt[c(4,3,5,1,2)]}


scanText<-function(string, what = character(0), ...){
  ## Like scan() but reading from a vector of character strings
  tc <- textConnection(string)
  result <- scan(tc, what = what, quiet = TRUE, ...)
  close(tc)
  return(result)}

getplotdat1 <- function(h="",plotrepfile,skip=1) {
  ##=================================================
  ## List single line after header h. 
  ##=================================================
  dat <- readLines(plotrepfile)
  recnum <- grep(h, dat)
  scanText(dat[recnum + skip],what=0)
}

getplotdat4 <- function(h="",plotrepfile) {
  ##=================================================
  ## Start listing after header h.  Quit if encounter
  ##  "^#"
  ##=================================================
  dat <- readLines(plotrepfile)
  rec1 <- grep(h, dat)
  if(length(rec1) <= 0)
    stop(paste('"',h,'"',"not found in",plotrepfile," Die yuppie scum!"))
  recnum <- rec1+1
  tt <- numeric(0)
  for(i in recnum:length(dat)) {
    if (regexpr("^#", dat[i]) != -1) break
    tt <- c(tt, scanText(dat[i], what = 0))
  }
  tt
}


