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
  init_int nIdx
  init_matrix Idata(1,3,1,ni)  // Year | I
  
  // Vectors to hold catch and year after read in above
  ivector Cyear(1,nc)
  vector      C(1,nc)

  // Vectors to hold indices and year
  ivector Iyear(1,ni)
  ivector   Idx(1,ni)
  vector      I(1,ni)		
  ivector     X(1,ni)  // years with abundance index: 1995 | 1998 | ...
  vector   logI(1,ni)
  vector   lhat(1,ni)
  
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

  init_bounded_number_vector q(1,nIdx,qLo,qHi,qPh)  
  init_bounded_number_vector s(1,nIdx,sLo,sHi,sPh)  
  
  number pen
  
  vector nll(1,nIdx)
  vector rss(1,nIdx)
  vector n(  1,nIdx)
  vector inv(1,nIdx)
  
  // Updated pop time series
  vector B(   1,nc+1)
  vector F(   1,nc)
  vector Ihat(1,nc)
    
  // Report matrix
  matrix summary(1,nc,1,4)  // Year | B | C | F | I 

  // Objfun
  objective_function_value neglogL

PRELIMINARY_CALCS_SECTION
  
  pen=0.0;
   
  // Parameters
  r    = r_plui[4];
  k    = k_plui[4];
  p    = p_plui[4];
  a    = a_plui[4];
  
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
    q(j) = qPr[j];
    s(j) = sPr[j];
    }

  trace<<"nll"<<" "<<"ss"<<" "<<"k"<<" "<<"r"<<endl;
 
PROCEDURE_SECTION
  fwd();
  qs();
  ll();
  
  trace<<neglogL<<" "<<rss<<" "<<r<<" "<<k<<endl;

REPORT_SECTION

  summary();
  
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

  lls<<endl;
 
 
FUNCTION fwd
  
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

FUNCTION qs
     
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
  
FUNCTION ll

  neglogL = 0.5*ni*log(2*pi);
  for (int j=1; j<=ni; j++){
     if (lav==1)
	 	  neglogL += pow(log(I[j])-log(Ifit[j]),2.0);

     if (lav!=1)
	 	  neglogL += sfabs(log(I[j])-log(Ifit[j]));
	 }


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
     
  //halflog2pi = 0.5*log(2*pi);
  
    
FUNCTION summary

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


TOP_OF_MAIN_SECTION
  arrmblsize = 40000000L;
  gradient_structure::set_GRADSTACK_BUFFER_SIZE(3000000);
  gradient_structure::set_CMPDIF_BUFFER_SIZE(200000);
  gradient_structure::set_MAX_NVAR_OFFSET(10000);
  gradient_structure::set_MAX_DLINKS(100000);

