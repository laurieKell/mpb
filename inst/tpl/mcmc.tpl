  int mcmc_iteration = 0;

  ofstream mcmc_par("mcmc_par.csv");
  ofstream mcmc_bio("mcmc_bio.csv");

FUNCTION write_mcmc
  mcmc_iteration++;
  // Parameters
  if(mcmc_iteration == 1){
    mcmc_par<<"neglogL,r,k,b0,p,q,s"<<endl;
    mcmc_par<<neglogL<<","
          <<r      <<","
          <<k      <<","
          <<a      <<","
          <<p      <<","
          <<q      <<","
          <<s  <<endl;}
  // Biomass
  if(mcmc_iteration == 1){
    mcmc_bio<<Cyear[1];
    for(int t=2; t<=nc; t++)
      mcmc_bio<<","<<Cyear[t];
     mcmc_bio<<endl;
     }

  mcmc_bio<<B[1];
  for(int t=2; t<=nc; t++)
    mcmc_bio<<","<<B[t];
  mcmc_bio<<endl;
