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
  
#include <admodel.h>
#include <contrib.h>

  extern "C"  {
    void ad_boundf(int i);
  }
#include <pella.htp>

model_data::model_data(int argc,char * argv[]) : ad_comm(argc,argv)
{
  nc.allocate("nc");
  Cdata.allocate(1,2,1,nc,"Cdata");
  ni.allocate("ni");
  nIdx.allocate("nIdx");
  Idata.allocate(1,3,1,ni,"Idata");
  Cyear.allocate(1,nc);
  C.allocate(1,nc);
  Iyear.allocate(1,ni);
  Idx.allocate(1,ni);
  I.allocate(1,ni);
  X.allocate(1,ni);
  logI.allocate(1,ni);
 string run_name = string(adprogram_name);
 if(option_match(argc,argv,"-ind") > -1){
   run_name = argv[option_match(argc,argv,"-ind") + 1];
   run_name = run_name.substr(0, run_name.rfind("."));}
 change_datafile_name((adstring)run_name.c_str() + ".ctl");
  _r_plui.allocate(1,4,"_r_plui");
  _k_plui.allocate(1,4,"_k_plui");
  _p_plui.allocate(1,4,"_p_plui");
  _a_plui.allocate(1,4,"_a_plui");
  qPh.allocate(1,nIdx,"qPh");
  qLo.allocate(1,nIdx,"qLo");
  qHi.allocate(1,nIdx,"qHi");
  qPr.allocate(1,nIdx,"qPr");
  sPh.allocate(1,nIdx,"sPh");
  sLo.allocate(1,nIdx,"sLo");
  sHi.allocate(1,nIdx,"sHi");
  sPr.allocate(1,nIdx,"sPr");
 change_datafile_name((adstring)run_name.c_str() + ".prr");
  r_prior.allocate(1,4,"r_prior");
  k_prior.allocate(1,4,"k_prior");
  p_prior.allocate(1,4,"p_prior");
  a_prior.allocate(1,4,"a_prior");
  msy_prior.allocate(1,4,"msy_prior");
  bmsy_prior.allocate(1,4,"bmsy_prior");
  fmsy_prior.allocate(1,4,"fmsy_prior");
  q_prior.allocate(1,4,"q_prior");
  s_prior.allocate(1,4,"s_prior");
 change_datafile_name((adstring)run_name.c_str() + ".ref");
  ref.allocate(1,4,"ref");
 change_datafile_name((adstring)run_name.c_str() + ".obj");
  lav.allocate("lav");
}

model_parameters::model_parameters(int sz,int argc,char * argv[]) : 
 model_data(argc,argv) , function_minimizer(sz)
{
  initializationfunction();
 phz = (int) _r_plui[1];
 lb  =       _r_plui[2];
 ub  =       _r_plui[3];
  _r.allocate(lb,ub,phz,"_r");
 phz = (int) _k_plui[1];
 lb  =       _k_plui[2];
 ub  =       _k_plui[3];
  _k.allocate(lb,ub,phz,"_k");
 phz = (int) _a_plui[1];
 lb  =       _a_plui[2];
 ub  =       _a_plui[3];
  _a.allocate(lb,ub,phz,"_a");
 phz = (int) _p_plui[1];
 lb  =       _p_plui[2];
 ub  =       _p_plui[3];
  _p.allocate(lb,ub,phz,"_p");
  _q.allocate(1,nIdx,qLo,qHi,qPh,"_q");
  _s.allocate(1,nIdx,sLo,sHi,sPh,"_s");
  r.allocate("r");
  k.allocate("k");
  a.allocate("a");
  p.allocate("p");
  q.allocate(1,nIdx,"q");
  s.allocate(1,nIdx,"s");
  cnow.allocate("cnow");
  bnow.allocate("bnow");
  fnow.allocate("fnow");
  bthen.allocate("bthen");
  fthen.allocate("fthen");
  bnowthen.allocate("bnowthen");
  fnowthen.allocate("fnowthen");
  msy.allocate("msy");
  bmsy.allocate("bmsy");
  fmsy.allocate("fmsy");
  cmsy.allocate("cmsy");
  bbmsy.allocate("bbmsy");
  ffmsy.allocate("ffmsy");
  bk.allocate("bk");
  fr.allocate("fr");
  bratio.allocate("bratio");
  fratio.allocate("fratio");
  slopeb.allocate("slopeb");
  slopef.allocate("slopef");
  ll.allocate(1,nIdx,"ll");
  pen.allocate("pen");
  #ifndef NO_AD_INITIALIZE
  pen.initialize();
  #endif
  SS.allocate("SS");
  #ifndef NO_AD_INITIALIZE
  SS.initialize();
  #endif
  s_.allocate("s_");
  #ifndef NO_AD_INITIALIZE
  s_.initialize();
  #endif
  nii.allocate(1,2,"nii");
  #ifndef NO_AD_INITIALIZE
    nii.initialize();
  #endif
  tmpn.allocate(1,nIdx,"tmpn");
  #ifndef NO_AD_INITIALIZE
    tmpn.initialize();
  #endif
  tmpq.allocate(1,nIdx,"tmpq");
  #ifndef NO_AD_INITIALIZE
    tmpq.initialize();
  #endif
  tmps.allocate(1,nIdx,"tmps");
  #ifndef NO_AD_INITIALIZE
    tmps.initialize();
  #endif
  nll.allocate(1,nIdx,"nll");
  #ifndef NO_AD_INITIALIZE
    nll.initialize();
  #endif
  rss.allocate(1,nIdx,"rss");
  #ifndef NO_AD_INITIALIZE
    rss.initialize();
  #endif
  n.allocate(1,nIdx,"n");
  #ifndef NO_AD_INITIALIZE
    n.initialize();
  #endif
  _bratio.allocate("_bratio");
  #ifndef NO_AD_INITIALIZE
  _bratio.initialize();
  #endif
  _fratio.allocate("_fratio");
  #ifndef NO_AD_INITIALIZE
  _fratio.initialize();
  #endif
  xy.allocate("xy");
  #ifndef NO_AD_INITIALIZE
  xy.initialize();
  #endif
  x.allocate("x");
  #ifndef NO_AD_INITIALIZE
  x.initialize();
  #endif
  y.allocate("y");
  #ifndef NO_AD_INITIALIZE
  y.initialize();
  #endif
  xx.allocate("xx");
  #ifndef NO_AD_INITIALIZE
  xx.initialize();
  #endif
  B.allocate(1,nc+1,"B");
  #ifndef NO_AD_INITIALIZE
    B.initialize();
  #endif
  F.allocate(1,nc,"F");
  #ifndef NO_AD_INITIALIZE
    F.initialize();
  #endif
  Bfit.allocate(1,ni,"Bfit");
  #ifndef NO_AD_INITIALIZE
    Bfit.initialize();
  #endif
  Ifit.allocate(1,ni,"Ifit");
  #ifndef NO_AD_INITIALIZE
    Ifit.initialize();
  #endif
  RSS.allocate(1,nIdx,"RSS");
  #ifndef NO_AD_INITIALIZE
    RSS.initialize();
  #endif
  summary.allocate(1,nc,1,8,"summary");
  #ifndef NO_AD_INITIALIZE
    summary.initialize();
  #endif
  lpbnow.allocate("lpbnow");
  lpfnow.allocate("lpfnow");
  lpbthen.allocate("lpbthen");
  lpfthen.allocate("lpfthen");
  lpbnowthen.allocate("lpbnowthen");
  lpfnowthen.allocate("lpfnowthen");
  lpmsy.allocate("lpmsy");
  lpbmsy.allocate("lpbmsy");
  lpfmsy.allocate("lpfmsy");
  lpcmsy.allocate("lpcmsy");
  lpbbmsy.allocate("lpbbmsy");
  lpffmsy.allocate("lpffmsy");
  lpbratio.allocate("lpbratio");
  lpfratio.allocate("lpfratio");
  lpslopeb.allocate("lpslopeb");
  lpslopef.allocate("lpslopef");
  lpr.allocate("lpr");
  lpk.allocate("lpk");
  lpfr.allocate("lpfr");
  lpbk.allocate("lpbk");
  neglogL.allocate("neglogL");
  prior_function_value.allocate("prior_function_value");
  likelihood_function_value.allocate("likelihood_function_value");
}

void model_parameters::preliminary_calculations(void)
{

#if defined(USE_ADPVM)

  admaster_slave_variable_interface(*this);

#endif
  
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
 
}

void model_parameters::userfunction(void)
{
  neglogL =0.0;
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
}

void model_parameters::report(const dvector& gradients)
{
 adstring ad_tmp=initial_params::get_reportfile_name();
  ofstream report((char*)(adprogram_name + ad_tmp));
  if (!report)
  {
    cerr << "error trying to open report file"  << adprogram_name << ".rep";
    return;
  }
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
        <<"# s"      <<endl<<s      <<endl
        <<"#-ll"     <<endl<<nll    <<endl
        <<"# neglogL"<<endl<<neglogL<<endl<<endl;
  report<<setprecision(12)
        <<"# Model summary"<<endl
        <<" year stock catch index hat stockHat stock F."<<endl
        <<summary<<endl;
  write_ll();
}

void model_parameters::get_fit(void)
{
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
}

void model_parameters::get_neglogL(void)
{
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
}

void model_parameters::get_summary(void)
{
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
}

void model_parameters::write_mcmc(void)
{
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
}

void model_parameters::write_priors(void)
{
  priors<<r_prior<<endl;
}

void model_parameters::write_bounds(void)
{
  bounds<<q	
        <<s<<endl;	
}

void model_parameters::write_ll(void)
{
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
}

void model_parameters::write_debug(void)
{
  debug<<_p_plui[1]<<" "<<(_p_plui[1]<(-1))<<endl;
}

void model_parameters::write_trace(void)
{
  get_neglogL();
  dvariable ss=0.0;
  for (int j=1; j<=ni; j++)
      ss += pow(log(I[j])-log(Ifit[j]),2.0);
  neglogL=ss;
  trace<<neglogL<<","<<ss<<","<<k<<","<<r<<","<<p<<endl;
}

model_data::~model_data()
{}

model_parameters::~model_parameters()
{}

void model_parameters::final_calcs(void){}

void model_parameters::set_runtime(void){}

#ifdef _BORLANDC_
  extern unsigned _stklen=10000U;
#endif


#ifdef __ZTC__
  extern unsigned int _stack=10000U;
#endif

  long int arrmblsize=0;

int main(int argc,char * argv[])
{
    ad_set_new_handler();
  ad_exit=&ad_boundf;
  arrmblsize = 40000000L;
  gradient_structure::set_GRADSTACK_BUFFER_SIZE(3000000);
  gradient_structure::set_CMPDIF_BUFFER_SIZE(200000);
  gradient_structure::set_MAX_NVAR_OFFSET(10000);
  gradient_structure::set_MAX_DLINKS(100000);
    gradient_structure::set_NO_DERIVATIVES();
    gradient_structure::set_YES_SAVE_VARIABLES_VALUES();
    if (!arrmblsize) arrmblsize=15000000;
    model_parameters mp(arrmblsize,argc,argv);
    mp.iprint=10;
    mp.preliminary_calculations();
    mp.computations(argc,argv);
    return 0;
}

extern "C"  {
  void ad_boundf(int i)
  {
    /* so we can stop here */
    exit(i);
  }
}
