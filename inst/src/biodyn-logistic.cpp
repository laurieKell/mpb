// Implements the logistic function
//PELLA , J. J., 1967 A study of methods to estimate the Schaefer model parameters with special reference to the yellowfin tuna fishery in the eastern tropical Pacific ocean. University of Washington, Seattle.
//PELLA , J. J., and P. K. T OMLINSON , 1969 A generalized stock production model. Bulletin of the Inter-American Tropical Tuna Commission 13: 419-496.
//PRAGER , M. H., 1994 A suite of extensions to a nonequilibrium surplus-production model. U. S. Fishery Bulletin 92: 374-389.

#include <biodyn-logistic.h>
  
logistic::logistic(FLQuant F, FLQuant C, FLQuant B, FLQuant params){
  FLQuant ctc=C;
  FLQuant hvt=F;
  FLQuant stk=B;
  
  FLQuant par=params;
  }
     
void logistic::nr(double tolVal, int niter){
  
  for (int i=1; i<=stk.get_niter(); i++){
     int iters=0;
     int tol  =1;
     for (int yr=1; yr<=ctc.get_nyear(); yr++){
        while (tol>tolVal&iters<niter){
         iters++;
         double func = flq(ctc,1,yr,1,1,1,i)-yield(flq(hvt,1,yr,1,1,1,i),flq(stk,1,yr,1,1,1,i),flq(par1,yr,1,1,1,i),flq(par2,yr,1,1,1,i));
         double grad = gradY(flq(hvt,1,yr,1,1,1,i),flq(ctc,1,yr,1,1,1,i),flq(stk,1,yr,1,1,1,i),flq(par1,yr,1,1,1,i),flq(par2,yr,1,1,1,i));
         double val  =NewRhap(flq(hvt,1,yr,1,1,1,i),func,grad);
        
         tol=abs(flq(hvt,1,yr,1,1,1,i)-val);
         hvt(1,yr,1,1,1,i)=val;
         }
       
   stk(1,yr+1,1,1,1,i)=stock(flq(hvt,1,yr,1,1,1,i),flq(ctc,1,yr,1,1,1,i),
                             flq(stk,1,yr,1,1,1,i),flq(par1,yr,1,1,1,i),
                             flq(par2,yr,1,1,1,i));
   }}}
     
logistic::~logistic(void){
;  
}     

double logistic::gradY(double F, double C, double  B, double r, double K){

  double expr1  = r/K;
  double expr2  = F/expr1;
  double expr3  = expr1 * B;
  double expr4  = r - F;
  double expr5  = exp(expr4);
  double expr7  = expr3 * (1 - expr5);
  double expr9  = 1 - expr7/expr4;
  double expr10 = log(expr9);

  return(-(1/expr1*expr10-expr2*((expr3*expr5/expr4 + expr7/(expr4*expr4))/expr9)));
  }

 