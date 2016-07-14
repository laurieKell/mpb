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
//   Vectors that include abundance index years: u, X
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
  ofstream    lls("lls.txt");
  ofstream  trace("trace.txt");
  ofstream  debug("debug.txt");
  
DATA_SECTION

  // Read data file
  //catch time series which has to be complete
  init_int nc
  init_matrix Cdata(1,2,1,nc)  // Year | C

  //indices of abundance, multiple with missing years
  init_int ni
  init_int nU
  init_matrix Idata(1,3,1,ni)  // Year | I
  
  // Vectors to hold catch and year after read in above
  ivector Cyear(1,nc)
   vector     C(1,nc)

  // Vectors to hold indices and year
  ivector Iyear(1,ni)
   vector     u(1,ni)
  ivector   uNm(1,ni)
  ivector   uYr(1,ni)  // years with abundance index: 1995 | 1998 | ...
  
  // Constants for profile
  number stepN
  number stepSz

  // Read in data
  // Switch to control file
  !! string run_name = string(adprogram_name);
  !! if(option_match(argc,argv,"-ind") > -1){
  !!   run_name = argv[option_match(argc,argv,"-ind") + 1];
  !!   run_name = run_name.substr(0, run_name.rfind("."));}
  
  // Read control file (phase, lower, upper, init)
  !! change_datafile_name((adstring)run_name.c_str() + ".ctl");

  init_vector r_plui(1,4)
  init_vector k_plui(1,4)
  init_vector p_plui(1,4)
  init_vector a_plui(1,4)

  init_ivector qPh(1,nU)
  init_vector  qLo(1,nU)
  init_vector  qHi(1,nU)
  init_vector  qPr(1,nU)
  
  init_ivector sPh(1,nU)
  init_vector  sLo(1,nU)
  init_vector  sHi(1,nU)
  init_vector  sPr(1,nU)
 
  // Switch to prior file
  !! change_datafile_name((adstring)run_name.c_str() + ".prr");
  // Read wt, mean, sd
  init_vector r_prior(1,4)
  init_vector k_prior(1,4)
  init_vector p_prior(1,4)
  init_vector a_prior(1,4)
  init_vector msy_prior(1,4)
  init_vector bmsy_prior(1,4)
  init_vector fmsy_prior(1,4)
  init_vector q_prior(1,4)
  init_vector s_prior(1,4)

  // Switch to ref file for mng reference year etc
  !! change_datafile_name((adstring)run_name.c_str() + ".ref");
  
  // c(nyr=3,nreg=5,refyr=NA)
  init_vector ref(1,4)

  // Read flag for objective function, i.e. LAV
  !! change_datafile_name((adstring)run_name.c_str() + ".obj");
  init_number lav

PARAMETER_SECTION
  
  // Estimated
  !! phz = (int) r_plui[1];
  !! lb  =       r_plui[2];
  !! ub  =       r_plui[3];
  init_bounded_number r(lb,ub,phz)
  !! phz = (int) k_plui[1];
  !! lb  =       k_plui[2];
  !! ub  =       k_plui[3];
  init_bounded_number k(lb,ub,phz)
  !! phz = (int) a_plui[1];
  !! lb  =       a_plui[2];
  !! ub  =       a_plui[3];
  init_bounded_number a(lb,ub,phz)
  !! phz = (int) p_plui[1];
  !! lb  =       p_plui[2];
  !! ub  =       p_plui[3];
  init_bounded_number p(lb,ub,phz)
  
  //init_bounded_number_vector F(1,nc,0,1)

  init_bounded_number_vector q(1,nU,qLo,qHi,qPh)  
  init_bounded_number_vector s(1,nU,sLo,sHi,sPh)  

  //likelihood stuff  
  number pen

  vector nll(1,nU)
  vector  ss(1,nU)
  vector   n(1,nU)
  vector sm(1,nU)
  vector  _q(1,nU)
  
  // Updated pop time series
  vector B(1,nc+1)
  vector F(1,nc)
  vector I(1,nc)
    
  // Report matrix
  matrix summary(1,nc,1,5)  // Year | B | C | F | I 

  // Objfun
  objective_function_value neglogL

PRELIMINARY_CALCS_SECTION
  
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
 
PROCEDURE_SECTION
  fwd();
  qs();
  ll();
  
  trace<<neglogL<<" "<<nll<<" "<<ss<<" "<<r<<" "<<k<<endl;

REPORT_SECTION

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
 
FUNCTION fwd
  
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

FUNCTION qs
    
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

FUNCTION ll

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

FUNCTION setSummary
  summary.colfill(1,(dvector)Cyear);
  summary.colfill(2,B);
  summary.colfill(3,C);
  summary.colfill(4,F);
  summary.colfill(5,I);

TOP_OF_MAIN_SECTION
  arrmblsize = 40000000L;
  gradient_structure::set_GRADSTACK_BUFFER_SIZE(3000000);
  gradient_structure::set_CMPDIF_BUFFER_SIZE(200000);
  gradient_structure::set_MAX_NVAR_OFFSET(10000);
  gradient_structure::set_MAX_DLINKS(100000);

