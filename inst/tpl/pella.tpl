//=========================================================================================================================
// File:        pella.tpl
// Model:       Pella-Tomlinson model, with Binit=k*a
// Parameters:  r, k, a, p, q, s
// Fitted data: Abundance index
// Likelihood:  Log-transformed normal
// References:  Polacheck et al. (1993)
// Notes:       q and s are free parameters, to allow full uncertainty in MCMC analysis
// History:      9 Mar 2010 Arni Magnusson created, to benchmark against R optimizers
//               7 Oct 2010 Arni Magnusson improved string handling and comments
//               22 may 2016 Checked fo Alb WG
//=========================================================================================================================
// Implementation notes
//   Abundance index may not exist for all years
//   Vectors that include all years: B, C
//   Vectors that include abundance index years: I, Ifit, X
//   X links long and short vectors
//=========================================================================================================================

GLOBALS_SECTION
  #include "admodel.h"
  #include <string>
  #include <dnorm.cpp> //include functions from custom library
  #include <nr.cpp>    //include functions from custom library

  using std::string;

  const double pi = 3.141592654;
  int mcmc_iteration = 0;
  int phz;    // phase
  double lb;  // lower bound
  double ub;  // upper bound

  ofstream mcmc_par("mcmc_par.csv");
  ofstream mcmc_bio("mcmc_bio.csv");
  ofstream priors("prr.txt");
  ofstream bounds("bnd.txt");
  ofstream    lls("lls.txt");
  ofstream  trace("trace.txt");
  ofstream  debug("debug.txt");
  
DATA_SECTION
  // Read data [ile
  init_int nc
  init_matrix Cdata(1,2,1,nc)  // Year | C
  init_int ni
  init_int nIdx
  init_matrix Idata(1,3,1,ni)  // Year | I
  
  // Vectors
  ivector Cyear(1,nc)
  vector      C(1,nc)

  ivector Iyear(1,ni)
  ivector   Idx(1,ni)
  vector      I(1,ni)		
  ivector     X(1,ni)  // years with abundance index: 1995 | 1998 | ...
  vector   logI(1,ni)
  
  // Constants
  number halfnlog2pi
  number stepN
  number stepSz

  // Switch to control file
  !! string run_name = string(adprogram_name);
  !! if(option_match(argc,argv,"-ind") > -1){
  !!   run_name = argv[option_match(argc,argv,"-ind") + 1];
  !!   run_name = run_name.substr(0, run_name.rfind("."));}
  
  // Read control file (phase, lower, upper, init)
  !! change_datafile_name((adstring)run_name.c_str() + ".ctl");

  init_vector _r_plui(1,4)
  init_vector _k_plui(1,4)
  init_vector _p_plui(1,4)
  init_vector _a_plui(1,4)

  init_ivector qPh(1,nIdx)
  init_vector  qLo(1,nIdx)
  init_vector  qHi(1,nIdx)
  init_vector  qPr(1,nIdx)
  
  init_ivector sPh(1,nIdx)
  init_vector  sLo(1,nIdx)
  init_vector  sHi(1,nIdx)
  init_vector  sPr(1,nIdx)
 
  // Switch to prior file
  !! change_datafile_name((adstring)run_name.c_str() + ".prr");

  // Read prior file (wt, mean, sd)
  init_vector r_prior(1,4)
  init_vector k_prior(1,4)
  init_vector p_prior(1,4)
  init_vector a_prior(1,4)
  init_vector msy_prior(1,4)
  init_vector bmsy_prior(1,4)
  init_vector fmsy_prior(1,4)
  init_vector q_prior(1,4)
  init_vector s_prior(1,4)

  // Switch to ref file for reference year etc
  !! change_datafile_name((adstring)run_name.c_str() + ".ref");
  
  // c(nyr=3,nreg=5,refyr=NA)
  init_vector ref(1,4)

  // Read flag for objective function
  !! change_datafile_name((adstring)run_name.c_str() + ".obj");
  init_number lav

PARAMETER_SECTION
  // Estimated
  !! phz = (int) _r_plui[1];
  !! lb  =       _r_plui[2];
  !! ub  =       _r_plui[3];
  init_bounded_number _r(lb,ub,phz)
  !! phz = (int) _k_plui[1];
  !! lb  =       _k_plui[2];
  !! ub  =       _k_plui[3];
  init_bounded_number _k(lb,ub,phz)
  !! phz = (int) _a_plui[1];
  !! lb  =       _a_plui[2];
  !! ub  =       _a_plui[3];
  init_bounded_number _a(lb,ub,phz)
  !! phz = (int) _p_plui[1];
  !! lb  =       _p_plui[2];
  !! ub  =       _p_plui[3];
  init_bounded_number _p(lb,ub,phz)
  
  //init_bounded_number_vector F(1,nc,0,1)

  //change from log
  //init_bounded_number_vector logq(1,nIdx,qLo,qHi,qPh)  
  //init_bounded_number_vector logs(1,nIdx,sLo,sHi,sPh)  
  init_bounded_number_vector _q(1,nIdx,qLo,qHi,qPh)  
  init_bounded_number_vector _s(1,nIdx,sLo,sHi,sPh)  
  
  // Derived
  sdreport_number r
  sdreport_number k
  sdreport_number a
  sdreport_number p
  sdreport_vector q(1,nIdx)
  sdreport_vector s(1,nIdx)
  sdreport_number cnow
  sdreport_number bnow
  sdreport_number fnow
  sdreport_number bthen
  sdreport_number fthen
  sdreport_number bnowthen
  sdreport_number fnowthen
  sdreport_number msy
  sdreport_number bmsy
  sdreport_number fmsy
  sdreport_number cmsy
  sdreport_number bbmsy
  sdreport_number ffmsy
  sdreport_number bk
  sdreport_number fr
  sdreport_number bratio
  sdreport_number fratio
  sdreport_number slopeb;
  sdreport_number slopef;
  sdreport_vector ll(1,nIdx)
  
  number pen
  number SS
  number s_
  vector nii(1,2)
  vector tmpn(1,nIdx)
  vector tmpq(1,nIdx)
  vector tmps(1,nIdx)
  vector nll(1,nIdx)
  vector rss(1,nIdx)
  vector n(1,nIdx)
  
  number _bratio
  number _fratio

  number  xy  
  number  x   
  number  y 
  number  xx

  // Updated
  vector B(1,nc+1)
  vector F(1,nc)
  vector Bfit(1,ni)
  vector Ifit(1,ni)
  vector RSS(1,nIdx)
    
  // Report
  matrix summary(1,nc,1,8)  // Year | B | C | I | Ifit | stockHat | stock | F.

  // likelihood profile
  likeprof_number lpbnow;
  likeprof_number lpfnow;
  likeprof_number lpbthen;
  likeprof_number lpfthen;
  likeprof_number lpbnowthen;
  likeprof_number lpfnowthen;
  likeprof_number lpmsy;
  likeprof_number lpbmsy;
  likeprof_number lpfmsy;
  likeprof_number lpcmsy;
  likeprof_number lpbbmsy;
  likeprof_number lpffmsy;
  likeprof_number lpbratio;
  likeprof_number lpfratio;
  likeprof_number lpslopeb;
  likeprof_number lpslopef;
  likeprof_number lpr; 
  likeprof_number lpk;
  likeprof_number lpfr;
  likeprof_number lpbk;
  
  // Objfun
  objective_function_value neglogL

PRELIMINARY_CALCS_SECTION
  
   pen=0.0;
   
   // Parameters
  _r    = _r_plui[4];
  _k    = _k_plui[4];
  _a    = _a_plui[4];
  _p    = _p_plui[4];
 
  halfnlog2pi = 0.5*ni*log(2*pi);
  stepN =50;
  stepSz=0.05;

  // Data
  Cyear = (ivector) row(Cdata,1);
  C     =           row(Cdata,2);
  
  Iyear = (ivector) row(Idata,1);
  I     =           row(Idata,2);
  Idx   = (ivector) row(Idata,3);
  logI  = log(I);

  X     = Iyear - Cyear[1] + 1;


  for(int i=1; i<=nc; i++) F[i]=_r*0.2;
 
  for (int j=1; j<=nIdx; j++){
    // change from log
    //logq(j) = qPr[j];
    //logs(j) = sPr[j];
    _q(j) = qPr[j];
    _s(j) = sPr[j];
    }

  // likelihood profile 
  lpbnow  =bnow;
  lpfnow  =fnow;
  lpbnow  =bthen;
  lpfnow  =fthen;
  lpbnow  =bnowthen;
  lpfnow  =fnowthen;

  if (1!=2){
	lpmsy   =msy;
	lpbmsy  =bmsy;
	lpfmsy  =fmsy;
	lpcmsy  =cmsy;
	lpbbmsy =bbmsy;
	lpffmsy =ffmsy;
	lpbratio=bratio;
	lpfratio=fratio;
	lpslopeb=slopeb;
	lpslopef=slopef;
	lpr     =_r;
	lpk     =_k;
	lpfr    =fr;
	lpbk    =bk;
	}

  lpbnow.set_stepnumber(stepN);
  lpbnow.set_stepsize(stepSz);

  lpfnow.set_stepnumber(stepN);
  lpfnow.set_stepsize(stepSz);

  lpbthen.set_stepnumber(stepN);
  lpbthen.set_stepsize(stepSz);

  lpfnowthen.set_stepnumber(stepN);
  lpfnowthen.set_stepsize(stepSz);

  lpbnowthen.set_stepnumber(stepN);
  lpbnowthen.set_stepsize(stepSz);

  lpfnow.set_stepnumber(stepN);
  lpfnow.set_stepsize(stepSz);

  if (1!=2){

  lpmsy.set_stepnumber(stepN);
  lpmsy.set_stepsize(stepSz);

  lpfmsy.set_stepnumber(stepN);
  lpfmsy.set_stepsize(stepSz);

  lpbmsy.set_stepnumber(stepN);
  lpbmsy.set_stepsize(stepSz);

  lpcmsy.set_stepnumber(stepN);
  lpcmsy.set_stepsize(stepSz);

  lpbbmsy.set_stepnumber(stepN);
  lpbbmsy.set_stepsize(stepSz);

  lpffmsy.set_stepnumber(stepN);
  lpffmsy.set_stepsize(stepSz);

  lpbratio.set_stepnumber(stepN);
  lpbratio.set_stepsize(stepSz);

  lpfratio.set_stepnumber(stepN);
  lpfratio.set_stepsize(stepSz);

  lpslopeb.set_stepnumber(stepN);
  lpslopeb.set_stepsize(stepSz);

  lpslopef.set_stepnumber(stepN);
  lpslopef.set_stepsize(stepSz);

  lpr.set_stepnumber(stepN);
  lpr.set_stepsize(stepSz);

  lpk.set_stepnumber(stepN);
  lpk.set_stepsize(stepSz);

  lpfr.set_stepnumber(stepN);
  lpfr.set_stepsize(stepSz);

  lpbk.set_stepnumber(stepN);
  lpbk.set_stepsize(stepSz);
  }
  
  trace<<"ll"<<","<<"ss"<<","<<"k"<<","<<"r"<<","<<"p"<<endl;
 
PROCEDURE_SECTION
  get_fit();
  get_neglogL();

  if(mceval_phase())
    write_mcmc();

  // likelihood profile 
  lpbnow  =bnow;
  lpfnow  =fnow;
  lpbthen =bthen;
  lpfthen =fthen;
  lpbnowthen  =bnowthen;
  lpfnowthen  =fnowthen;

  if (1!=2){
  lpmsy   =msy;
  lpbmsy  =bmsy;
  lpfmsy  =fmsy;
  lpcmsy  =cmsy;
  lpbbmsy =bbmsy;
  lpffmsy =ffmsy;
  lpbratio=bratio;
  lpfratio=fratio;
  lpslopeb=slopeb;
  lpslopef=slopef;
  lpr     =_r;
  lpk     =_k;
  lpfr    =fr;
  lpbk    =bk;
  }

REPORT_SECTION
  write_bounds();
  write_priors();
  
  summary.initialize();
  get_summary();
  get_fit();
   
  report<<setprecision(12)
        <<"# r"      <<endl<<r      <<endl
        <<"# k"      <<endl<<k      <<endl
        <<"# b0"     <<endl<<a      <<endl
        <<"# p"      <<endl<<p      <<endl
        <<"# q"      <<endl<<q      <<endl
        <<"# s"      <<endl<<_s      <<endl
        <<"#-ll"     <<endl<<nll    <<endl
        <<"# neglogL"<<endl<<neglogL<<endl<<endl;
  report<<setprecision(12)
        <<"# Model summary"<<endl
        <<" year stock catch index hat stockHat stock F."<<endl
        <<summary<<endl;

  write_ll();
 
 
FUNCTION get_fit
  //r = mfexp(logr);
  
  r=_r;
  k=_k;
  a=_a;
  p=_p;
  
  for (int j=1; j<=nIdx; j++){
    //change from log
    //q[j]   = mfexp(logq[j]);
    //s[j]   = mfexp(logs[j]);
    q[j]   = _q[j];
    s[j]   = _s[j];
    }
  
  pen=0.0;
  B[1] = k*a;
  for(int t=1; t<=nc; t++){
    if (_p_plui[1]<(-1)){

       F[t]=nr(-log(1-C[t]/B[t])*.5, C[t], B[t], r, k);

       //(r-F[t])
       //(sfabs(r-F[t])-sfabs(F[t]-r))
       dvariable alpha=sfabs(r-F[t]);//-sfabs(F[t]-r));

       B[t+1]=((r-F[t]))*B[t]*exp((alpha))/(alpha+(r/k)*B[t]*(exp(alpha)-1));
       B[t+1]=sfabs(B[t+1]);
     }else{
       dvariable now=posfun(B[t]-C[t],.001,pen);
  
       B[t+1]=now+r/p*B[t]*(1-pow(B[t]/k,p));
       }
    }
     
  // constricted likelihood
  bool flag=false;
  
  //initialise
  for (int i=1; i<=nIdx; i++){
    tmpn(i)=0;
    tmpq(i)=0;
    tmps(i)=0;
    if ((qPh(i)< -1) || (sPh(i)< -1))
       flag=true;
    }

  // calculate
  if (flag){
    for (int j=1; j<=ni; j++){
      tmpn(Idx[j])=tmpn(Idx[j])+1;
      tmpq(Idx[j])=tmpq(Idx[j])+log(I[Idx[j]])-log(0.5*(B(X[j])+B(X[j]+1)));
      }

    for (int i=1; i<=nIdx; i++)
      tmpq(i)=exp(tmpq(i)/tmpn[i]);
    
    for (int j=1; j<=ni; j++)
      tmps[Idx[j]]=+log(I[Idx[j]])-log(0.5*(B(X[j])+B(X[j]+1))*tmpq(Idx[j]));

    for (int i=1; i<=nIdx; i++){
      if (qPh[i]<-1) q[i]=tmpq(i);
      if (sPh[i]<-1) s[i]=pow((tmps(i)/tmpn[i]),.5);
      }
    }

  for (int j=1; j<=ni; j++)
     //Ifit[j] = B(X[j])*q(Idx[j]);
     Ifit[j] = 0.5*(B(X[j])+B(X[j]+1))*q(Idx[j]);
  
  write_trace();
  
  cnow =C[nc];
  fnow =C[nc]/B[nc];
  bnow =B[nc];
  fthen =C[nc-ref[4]+1]/B[nc-ref[4]+1];
  bthen =B[nc-ref[4]+1];
  fnowthen  =fnow/fthen;
  bnowthen =bnow/bthen;
  
  msy  =r*k*pow(1/(1+p),1/p+1);
  bmsy =k*pow((1/(1+p)),(1/p));
  fmsy =msy/bmsy;
  cmsy =C[nc]	/msy;
  bbmsy=bnow/bmsy;
  ffmsy=fnow/fmsy;
  
  bk=bnow/k;
  fr=fnow/r;
  
   bratio=0.0;
   fratio=0.0;
  _bratio=0.0;
  _fratio=0.0;
  
  for (int i=nc; i>nc-ref[2]; i--){
    bratio+=B[i];
    fratio+=C[i]/B[i];
    
    _bratio+=B[i-ref[3]];
    _fratio+=C[i-ref[3]]/B[i-ref[3]];
    }
 bratio=bratio/_bratio;
 fratio=fratio/_fratio;

  xy=0.0; 
  x =0.0;
   y=0.0; 
  xx=0.0; 
  for (int i=nc; i>nc-ref[1]; i--){
    x +=i;
    xx+=i*i;  
    y +=B[i];  
    xy+=i*B[i];  
    }
  slopeb = -(ref[1]*xy - x*y)/(ref[1]*xx - x*x);
   
  xy=0.0; 
  x =0.0;
   y=0.0; 
  xx=0.0; 
  for (int i=nc; i>nc-ref[1]; i--){
    x +=i;
    xx+=i*i;  
    y +=C[i]/B[i];  
    xy+=i*C[i]/B[i];  
    }
  slopef = -(ref[1]*xy - x*y)/(ref[1]*xx - x*x);
   
  //neglogL = 0; //halfnlog2pi;
  //for (int j=1; j<=ni; j++)
  //  // neglogL += log(_s[Idx[j]])+  
  //  //            pow(log(I[j])-log(Ifit[j]),2.0)/(2*_s[Idx[j]]*_s[Idx[j]]);
  //  neglogL += pow(log(I[j])-log(Ifit[j]),2.0);
  
FUNCTION get_neglogL

  neglogL = 0.5*ni*log(2*pi);
  for (int j=1; j<=ni; j++){
     if (lav==1)
	 	  neglogL += pow(log(I[j])-log(Ifit[j]),2.0);

     if (lav!=1)
	 	  neglogL += sfabs(log(I[j])-log(Ifit[j]));
	 }
     

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
    
FUNCTION get_summary
  summary.colfill(1,(dvector)Cyear);
  summary.colfill(2,B);
  summary.colfill(3,C);
  summary.colfill(8,F);

  for(int i=1; i<=ni; i++)  // allow missing years in abundance index
    {
    summary(X[i],4) = I[i];
    summary(X[i],5) = Ifit[i];
    summary(X[i],6) = Ifit[i]/q(Idx(i));
    summary(X[i],7) = (B[X[i]]+B[X[i]+1])/2.0;
    }

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

FUNCTION write_priors
  priors<<r_prior<<endl;
	
FUNCTION write_bounds
  bounds<<q	
        <<s<<endl;	

FUNCTION write_ll
  
  for (int i=1; i<=nIdx; i++){
      n[i]  =0;
      rss[i]=0;
      nll[i]=0;}
       
  for (int i=1; i<=ni; i++){
	n[  Idx[i]]+=1;
    rss[Idx[i]]+=pow(log(I[i])-log(Ifit[i]),2.0);}
    
  for (int i=1; i<=nIdx; i++)
    nll[i] = halfnlog2pi*n[i]/2+
              rss[i]/(2*_s[i]*_s[i])*n[i]/2+
              (2*_s[i]*_s[i])*n[i]/2;
  
  //  ll(i)= n/2*log(3.14159265359*2)
  //        +n/2*(log(se*se))
  //        +ss/(2*se*se);
          

   lls<<nll<<n<<endl;
 
FUNCTION write_debug
  debug<<_p_plui[1]<<" "<<(_p_plui[1]<(-1))<<endl;

FUNCTION write_trace
  get_neglogL();

  dvariable ss=0.0;
  for (int j=1; j<=ni; j++)
      ss += pow(log(I[j])-log(Ifit[j]),2.0);
  
  neglogL=ss;
       
  trace<<neglogL<<","<<ss<<","<<k<<","<<r<<","<<p<<endl;


TOP_OF_MAIN_SECTION
  arrmblsize = 40000000L;
  gradient_structure::set_GRADSTACK_BUFFER_SIZE(3000000);
  gradient_structure::set_CMPDIF_BUFFER_SIZE(200000);
  gradient_structure::set_MAX_NVAR_OFFSET(10000);
  gradient_structure::set_MAX_DLINKS(100000);

