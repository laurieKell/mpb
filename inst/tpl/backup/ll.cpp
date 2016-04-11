
FUNCTION get_neglogL
  neglogL = halfnlog2pi;
  for (int j=1; j<=ni; j++){  
     //change from log
     //s_       = mfexp(logs[Idx[j]]);
     s_       = _s[Idx[j]];
     neglogL += log(s_)   
             +  pow(log(I[j])-log(Ifit[j]),2.0)/(2*s_*s_);
     }

  //neglogL = halfnlog2pi + ni*log(s(1)) + RSS[1]/(2*s(1)*s(1));
  
  // weighted likelihood priors
  dvariable _msy  =r*k*pow(1/(1+p),1/p+1);
  dvariable _bmsy =(k*pow((1/(1+p)),(1/p)));
  dvariable _fmsy =msy/bmsy;
 
  if (r_prior[1]>0)    neglogL += r_prior[1]*dnorm(r, r_prior[2], r_prior[3]); // /dnorm(r_prior[2], r_prior[2], r_prior[3]);
  if (k_prior[1]>0)    neglogL += k_prior[1]*dnorm(k, k_prior[2], k_prior[3]); // /dnorm(k_prior[2], k_prior[2], k_prior[3]);
  if (p_prior[1]>0)    neglogL += p_prior[1]*dnorm(p, p_prior[2], p_prior[3]); // /dnorm(p_prior[2], p_prior[2], p_prior[3]);
  if (a_prior[1]>0)    neglogL += a_prior[1]*dnorm(a, a_prior[2], a_prior[3]); // /dnorm(a_prior[2], a_prior[2], a_prior[3]);

  if ( msy_prior[1]>0)  neglogL +=  msy_prior[1]*dnorm(_msy,   msy_prior[2],  msy_prior[3]); // /dnorm( msy_prior[2],  msy_prior[2],  msy_prior[3]);
  if (bmsy_prior[1]>0)  neglogL += bmsy_prior[1]*dnorm(_bmsy, bmsy_prior[2], bmsy_prior[3]); // /dnorm(bmsy_prior[2], bmsy_prior[2], bmsy_prior[3]);
  if (fmsy_prior[1]>0)  neglogL += fmsy_prior[1]*dnorm(_fmsy, fmsy_prior[2], fmsy_prior[3]); // /dnorm(fmsy_prior[2], fmsy_prior[2], fmsy_prior[3]);

  //reference year in "a" of priors slot
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

  if(_p_plui[1]<-10){
     for (int t=1; t<=nc; t++){
       dvariable alpha=sfabs(r-F[t]);//-sfabs(F[t]-r);
     
       neglogL += (2*.1*.1)*(log(C[t])-log((F[t]/(r/k)*log(1-(r/k)*B[t]*(1-exp((alpha)))/(alpha)))))*(log(C[t])-log((F[t]/(r/k)*log(1-(r/k)*B[t]*(1-exp((alpha)))/(alpha)))));}}
   
  //if (bmsy_prior[1]==10) neglogL += bmsy_prior[1]*dnorm(_bmsy, C[(int)bmsy_prior[2]]/B[(int)bmsy_prior[2]], 
  //                                                             C[(int)bmsy_prior[2]]/B[(int)bmsy_prior[2]]*bmsy_prior[3]); 
  //if (fmsy_prior[1]==10) neglogL += fmsy_prior[1]*dnorm(_fmsy, C[(int)fmsy_prior[2]]/B[(int)fmsy_prior[2]], 
  //                                                             C[(int)fmsy_prior[2]]/B[(int)fmsy_prior[2]]*fmsy_prior[3]); 


 //for (i=1; i<=nIdx; i++){
  //  if (q_prior[i]>0) neglogL += q_prior[1]*dnorm(q, q_prior[2], q_prior[3]); // /dnorm(q_prior[2], q_prior[2], q_prior[3]);
  //  if (s_prior[i]>0) neglogL += s_prior[1]*dnorm(s, s_prior[2], s_prior[3]); // /dnorm(s_prior[2], s_prior[2], s_prior[3]);}