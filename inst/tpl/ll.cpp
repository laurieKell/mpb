
FUNCTION get_neglogL
  neglogL = halfnlog2pi;
  for (int j=1; j<=ni; j++)
     neglogL += log(_s[Idx[j]])   
             +  pow(log(I[j])-log(Ifit[j]),2.0)/(2*_s[Idx[j]]*_s[Idx[j]]);

  
  // weighted likelihood priors
  dvariable _msy  =r*k*pow(1/(1+p),1/p+1);
  dvariable _bmsy =(k*pow((1/(1+p)),(1/p)));
  dvariable _fmsy =msy/bmsy;
 
  if (r_prior[1]>0)    neglogL += r_prior[1]*dnorm(r, r_prior[2], r_prior[3]); 
  if (k_prior[1]>0)    neglogL += k_prior[1]*dnorm(k, k_prior[2], k_prior[3]);
  if (p_prior[1]>0)    neglogL += p_prior[1]*dnorm(p, p_prior[2], p_prior[3]);
  if (a_prior[1]>0)    neglogL += a_prior[1]*dnorm(a, a_prior[2], a_prior[3]);

  if ( msy_prior[1]>0)  neglogL +=  msy_prior[1]*dnorm(_msy,   msy_prior[2],  msy_prior[3]); 
  if (bmsy_prior[1]>0)  neglogL += bmsy_prior[1]*dnorm(_bmsy, bmsy_prior[2], bmsy_prior[3]); 
  if (fmsy_prior[1]>0)  neglogL += fmsy_prior[1]*dnorm(_fmsy, fmsy_prior[2], fmsy_prior[3]); 

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

  //bound on maximum F
  if(_p_plui[1]<-10){
     for (int t=1; t<=nc; t++){
       dvariable alpha=sfabs(r-F[t]);
     
       neglogL += (2*.1*.1)*(log(C[t])-
       log((F[t]/(r/k)*log(1-(r/k)*B[t]*(1-exp((alpha)))/
       (alpha)))))*(log(C[t])-log((F[t]/(r/k)*
       log(1-(r/k)*B[t]*(1-exp((alpha)))/(alpha)))));}}
  