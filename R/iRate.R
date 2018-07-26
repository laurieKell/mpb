#file:///home/laurence/Downloads/IOTC-2016-WPM07-08.pdf
#file:///home/laurence/Downloads/IOTC-2015-WPM06-09_-_SKJ_OM.pdf

#  Biomass threshold
#  Biomass limit
#  Fishing mortality target
#  Spawning biomass at MSY estimate
#  Lags in data collection, between years of application and last year of CPUE data, and management implementation,
# between proposed TAC and actual implementation to the shery.
#  Implementation error

 # bthreshold = 1
 # blim = 0.4 * SBMSY
 # ftarget = FMSY
 # DLAG=1, MLAG=1, SFREQ=1
 # errcpue=~0, effcpue=~0, errimp=~0

# Parameter Symbol Description
# Values
# evaluated
# Responsiveness Degree of smoothing in biomass index 0.5
# Target harvest
# rate muliplier
# Target harvest rate relative to historic levels
# i.e 0.9 = 90% of historic average
# 0.8, 0.9,
# 1.0, 1.1
# Threshold
# biomass index
# Biomass index at which the harvest rate is
# reduced relative to historic levels i.e. 0.7 =
#   reduce harvest rate when the biomas index is
# at 70% of historic levels
# 0.5, 0.6,
# 0.7, 0.8
# Limit biomass
# index
# Biomass index at which harvest rate is zero
# relative to historic levels i.e. 0.2 = close the
# fishery when the biomas index is at 20% of
# historic levels
# 0.05, 0.1,
# 0.2
# Maximum
# change
# Maximum allowable percenatge change in
# effort
# 0.4
# Maximum TAC
# Maximum total allowable catch (thousand
#                                tonnes)
# 300, 400,
# 500, 600

iRate<-function(yrs,cpue,catch,
                control=c(smooth=NA,responsiveness=0.5,
                          dlag=1,sfreq=2,mlag=1,
                          threshold=0.5,limit=0.2,hrMult= 1.1,
                          maxTac=1.5)){
  
  tc=catch[,ac(yrs - control["dlag"])]
  
  # Get control["smooth"]ed CPUE using last available timestep
  if (is.na(control["smooth"])){
    control["smooth"] <- cpue[,ac(yrs - control["dlag"])]
  } else {
    control["smooth"] <- control["responsiveness"] * cpue[,ac(yrs - control["dlag"])] +
      (1 - control["responsiveness"]) * control["smooth"]
  }
  
  # Calc recommended catch scalar
  rate   <- control["smooth"]
  rate[] <- (control["smooth"] - control["limit"]) * control["hrMult"] /
    (control["threshold"] - control["limit"])
  rate[control["smooth"] < control["limit"]] <- 0
  rate[control["smooth"] > control["threshold"]] <- control["hrMult"][control["smooth"] > control["threshold"]]
  
  # Set catch for next control["sfreq"] years, starting in yrs + control["mlag"]
  tac <- FLCore::expand(tc, year=seq(yrs + control["mlag"], yrs + control["mlag"] + control["sfreq"] - 1))
  tac[] <- rep(pmin(c(tc * rate), c(control["maxTac"])), each=control["sfreq"])
  
  tac}