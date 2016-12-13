iRate<-function(cpue,catch,y,
                cpue_smooth=NULL,responsiveness=0.5,
                DLAG=1,SFREQ=2,MLAG=1,
                biomass_threshold=0.5,biomass_limit=0.2,hr_multiplier= 1.1,
                maxTAC=1.5){
  
  tc=catch[,ac(y - DLAG)]
  
  # Get smoothed CPUE using last available timestep
  if (is.null(cpue_smooth)){
    cpue_smooth <- cpue[,ac(y - DLAG)]
  } else {
    cpue_smooth <- responsiveness * cpue[,ac(y - DLAG)] +
      (1 - responsiveness) * cpue_smooth
  }
  
  # Calc recommended catch scalar
  rate <- cpue_smooth
  rate[] <- (cpue_smooth - biomass_limit) * hr_multiplier /
    (biomass_threshold - biomass_limit)
  rate[cpue_smooth < biomass_limit] <- 0
  rate[cpue_smooth > biomass_threshold] <- hr_multiplier[cpue_smooth > biomass_threshold]
  
  # Set catch for next SFREQ years, starting in y + MLAG
  tac <- FLCore::expand(tc, year=seq(y + MLAG, y + MLAG + SFREQ - 1))
  tac[] <- rep(pmin(c(tc * rate), c(maxTAC)), each=SFREQ)
  
  tac}