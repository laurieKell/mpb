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
  nU.allocate("nU");
  Idata.allocate(1,3,1,ni,"Idata");
  Cyear.allocate(1,nc);
  C.allocate(1,nc);
  Iyear.allocate(1,ni);
  u.allocate(1,ni);
  uNm.allocate(1,ni);
  uYr.allocate(1,ni);
 string run_name = string(adprogram_name);
 if(option_match(argc,argv,"-ind") > -1){
   run_name = argv[option_match(argc,argv,"-ind") + 1];
   run_name = run_name.substr(0, run_name.rfind("."));}
 change_datafile_name((adstring)run_name.c_str() + ".ctl");
  r_plui.allocate(1,4,"r_plui");
  k_plui.allocate(1,4,"k_plui");
  p_plui.allocate(1,4,"p_plui");
  a_plui.allocate(1,4,"a_plui");
  qPh.allocate(1,nU,"qPh");
  qLo.allocate(1,nU,"qLo");
  qHi.allocate(1,nU,"qHi");
  qPr.allocate(1,nU,"qPr");
  sPh.allocate(1,nU,"sPh");
  sLo.allocate(1,nU,"sLo");
  sHi.allocate(1,nU,"sHi");
  sPr.allocate(1,nU,"sPr");
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
 phz = (int) r_plui[1];
 lb  =       r_plui[2];
 ub  =       r_plui[3];
  r.allocate(lb,ub,phz,"r");
 phz = (int) k_plui[1];
 lb  =       k_plui[2];
 ub  =       k_plui[3];
  k.allocate(lb,ub,phz,"k");
 phz = (int) a_plui[1];
 lb  =       a_plui[2];
 ub  =       a_plui[3];
  a.allocate(lb,ub,phz,"a");
 phz = (int) p_plui[1];
 lb  =       p_plui[2];
 ub  =       p_plui[3];
  p.allocate(lb,ub,phz,"p");
  q.allocate(1,nU,qLo,qHi,qPh,"q");
  s.allocate(1,nU,sLo,sHi,sPh,"s");
  pen.allocate("pen");
  #ifndef NO_AD_INITIALIZE
  pen.initialize();
  #endif
  nll.allocate(1,nU,"nll");
  #ifndef NO_AD_INITIALIZE
    nll.initialize();
  #endif
  ss.allocate(1,nU,"ss");
  #ifndef NO_AD_INITIALIZE
    ss.initialize();
  #endif
  n.allocate(1,nU,"n");
  #ifndef NO_AD_INITIALIZE
    n.initialize();
  #endif
  sm.allocate(1,nU,"sm");
  #ifndef NO_AD_INITIALIZE
    sm.initialize();
  #endif
  _q.allocate(1,nU,"_q");
  #ifndef NO_AD_INITIALIZE
    _q.initialize();
  #endif
  B.allocate(1,nc+1,"B");
  #ifndef NO_AD_INITIALIZE
    B.initialize();
  #endif
  F.allocate(1,nc,"F");
  #ifndef NO_AD_INITIALIZE
    F.initialize();
  #endif
  I.allocate(1,nc,"I");
  #ifndef NO_AD_INITIALIZE
    I.initialize();
  #endif
  summary.allocate(1,nc,1,5,"summary");
  #ifndef NO_AD_INITIALIZE
    summary.initialize();
  #endif
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
  r = r_plui[4];
  k = k_plui[4];
  p = p_plui[4];
  a = a_plui[4];
  
  stepN =50;
  stepSz=0.05;
  // Data
  Cyear = (ivector) row(Cdata,1);
  C     =           row(Cdata,2);
  
  Iyear = (ivector) row(Idata,1);  
  u     = (dvector) row(Idata,2); //index
  uNm   = (ivector) row(Idata,3); //name
  uYr   = Iyear - Cyear[1] + 1;
  
  for(int i=1; i<=nc; i++) F[i]=r*0.2;
 
  for (int j=1; j<=nU; j++){
    q(j) = qPr[j];
    s(j) = sPr[j];
    }
  trace<<"neglogL"<<"nll"<<" "<<"ss"<<" "<<"k"<<" "<<"r"<<endl;
 
}

void model_parameters::userfunction(void)
{
  neglogL =0.0;
  fwd();
  qs();
  ll();
  trace<<neglogL<<" "<<nll<<" "<<ss<<" "<<r<<" "<<k<<endl;
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
  setSummary();
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
        <<" year stock catch F index"<<endl
        <<summary<<endl;
  lls<<setprecision(12)
     <<nll<<ss<<s<<n<<endl;
}

void model_parameters::fwd(void)
{
  pen  =0.0;
  B[1] =k*a;
  for(int t=1; t<=nc; t++){
    //instantaneous
    if (p_plui[1]<(-1)){
       F[t]=nr(-log(1-C[t]/B[t])*.5, C[t], B[t], r, k);
       dvariable alpha=sfabs(r-F[t]);
       B[t+1]=((r-F[t]))*B[t]*exp((alpha))/(alpha+(r/k)*B[t]*(exp(alpha)-1));
       B[t+1]=sfabs(B[t+1]);
    }else{
    // harvest rate
       dvariable now=posfun(B[t]-C[t],.001,pen);
       B[t+1]=now+r/p*B[t]*(1-pow(B[t]/k,p));
       }
    //index   
    I(t)=0.5*(B[t]+B[t+1]);   
    }
}

void model_parameters::qs(void)
{
  //initialise
  for (int i=1; i<=nU; i++){
    nll(i)=0.0;
     ss(i)=0.0;
      n(i)=0.0;
    sm(i)=0.0;
     _q(i)=0.0;
    }
  for (int i=1; i<=ni; i++){
       n(uNm[i])+=1;
      _q(uNm[i])+=log(u(i)/I(uYr[i]));
     }
 for (int i=1; i<=nU; i++)
     _q[i]= exp(_q(i)/n(i));
  for (int i=1; i<=ni; i++)
    ss(uNm[i]) +=pow(log(I(uYr[i])*q(uNm(i)))-log(u(i)),2);
 for (int i=1; i<=nU; i++){
   if (qPh[i]<-1) q[i]=_q(i);
   if (sPh[i]<-1) s[i]=pow((ss(i)/n[i]),.5);
   }
}

void model_parameters::ll(void)
{
  double halflog2pi = 0.5*log(2*pi);
  if (lav!=1)
	  for (int i=1; i<=ni; i++)
   	  neglogL += sfabs((log(I(uYr(i))*q(uNm(i)))-log(u(i))));
	else{
    for (int i=1; i<=nU; i++){
	    nll(i)=0;
	    sm(i)=0;
	    ss(i) =0;}
	  for (int i=1; i<=ni; i++){
	    sm(uNm(i))+=log(u(i));
      ss[uNm(i)] +=pow(log(I(uYr[i])*q(uNm(i)))-log(u(i)),2.0);
      }
    neglogL=pen;    
    for (int i=1; i<=nU; i++){
	    nll[i] = sm(i)+
	             n[i]*log(2.0*pi)/2.0+
	             n(i)*log(s(i))+
	             ss(i)/(2.0*s(i));
      neglogL+=nll[i];
      //neglogL+=ss[i];
      }               
	 }
}

void model_parameters::setSummary(void)
{
  summary.colfill(1,(dvector)Cyear);
  summary.colfill(2,B);
  summary.colfill(3,C);
  summary.colfill(4,F);
  summary.colfill(5,I);
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
