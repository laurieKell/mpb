
  ofstream priors("prr.txt");

  // Switch to prior file
  !! change_datafile_name((adstring)run_name.c_str() + ".prr");

  if (false){
  // weighted likelihood priors
  dvariable _msy  =r*k*pow(1/(1+p),1/p+1);
  dvariable _bmsy =(k*pow((1/(1+p)),(1/p)));
  dvariable _fmsy =msy/bmsy;
 
  if (r_prior[1]>0)    neglogL += r_prior[1]*dnorm(r, r_prior[2], r_prior[3]); 
  if (k_prior[1]>0)    neglogL += k_prior[1]*dnorm(k, k_prior[2], k_prior[3]);
  if (p_prior[1]>0)    neglogL += p_prior[1]*dnorm(p, p_prior[2], p_prior[3]);
  if (a_prior[1]>0)    neglogL += a_prior[1]*dnorm(a, a_prior[2], a_prior[3]);

  if ( msy_prior[1]>0)  
    neglogL +=  msy_prior[1]*dnorm(_msy,   msy_prior[2],  msy_prior[3]); 
  if (bmsy_prior[1]>0)  
    neglogL += bmsy_prior[1]*dnorm(_bmsy, bmsy_prior[2], bmsy_prior[3]); 
  if (fmsy_prior[1]>0)  
    neglogL += fmsy_prior[1]*dnorm(_fmsy, fmsy_prior[2], fmsy_prior[3]); 

  //Reference year in "a" of priors slot if type==10
  if ( msy_prior[1]==10) {
      double a=value(C[(int)msy_prior[2]]);
      neglogL +=  msy_prior[1]*dnorm(_msy,a,a*msy_prior[2]); 
      }
  if (fmsy_prior[1]==10) {
      double a=value(C[(int)msy_prior[2]]/B[(int)  fmsy_prior[2]]);
      neglogL +=  fmsy_prior[1]*dnorm(_fmsy,a,a*fmsy_prior[2]); 
      }
  if (bmsy_prior[1]==10) {
      double a=value(B[(int)  bmsy_prior[2]]);
      neglogL += bmsy_prior[1]*dnorm(_bmsy,a,a*bmsy_prior[2]); 
      }
  }

  //penalty on catch, i.e. if catch cant be taken 
  neglogL += pen;
    
