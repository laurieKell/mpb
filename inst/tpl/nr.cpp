dvariable alpha(dvariable, dvariable);
dvariable beta( dvariable, dvariable);
dvariable yield(dvariable, dvariable, dvariable, dvariable); 
dvariable stock(dvariable, dvariable, dvariable, dvariable, dvariable );
dvariable gradY(dvariable, dvariable, dvariable, dvariable, dvariable);
dvariable nr(dvariable, dvariable, dvariable, dvariable, dvariable, double, int);

double alpha(double, double);
double beta( double, double);
double yield(double, double, double, double); 
double stock(double, double, double, double, double );
double gradY(double, double, double, double, double);
double _nr(double, double, double, double, double, double, int);

dvariable alpha(prevariable& r, prevariable& F);
dvariable beta( prevariable& r, prevariable& K);
dvariable yield(  prevariable& F,prevariable& B,prevariable& r,prevariable& K); 
dvariable stock(prevariable& F, prevariable&  C, prevariable&  B, prevariable& r, prevariable& K);
dvariable gradY(prevariable& F, prevariable&  C, prevariable& B, prevariable& r, prevariable& K);
dvariable _nr(prevariable& F, prevariable&  C, prevariable& B, prevariable& r, prevariable& K,double tolVal,int);
   
inline dvariable alpha(dvariable r, dvariable F) {return(r-F);}
inline dvariable beta( dvariable r, dvariable K) {return(r/K);}

inline dvariable yield(dvariable F,dvariable B,dvariable r,dvariable K){ 
     return((F/(r/K))*log(1-(r/K)*B*(1-exp((r-F)))/(r-F)));}

inline dvariable stock(dvariable F, dvariable C, dvariable  B, dvariable r, dvariable K){
    return(alpha(r,F)*B*exp(alpha(r,F))/(alpha(r,F)+beta(r,K)*B*(exp(alpha(r,F))-1)));}
 
dvariable gradY(dvariable F, dvariable C, dvariable B, dvariable r, dvariable K){

  dvariable expr1  = r/K;
  dvariable expr2  = F/expr1;
  dvariable expr3  = expr1*B;
  dvariable expr4  = (r-F);//-sfabs(F-r); //r-F; 
  dvariable expr5  = exp(expr4);
  dvariable expr7  = expr3*(1-expr5);
  dvariable expr9  = 1-expr7/expr4;
  dvariable expr10 = log(expr9);

  return(-(1/expr1*expr10-expr2*((expr3*expr5/expr4+expr7/(expr4*expr4))/expr9)));}

dvariable nr(dvariable F, dvariable C, dvariable B, 
             dvariable r, dvariable K,
             double tolVal=1e-10,int niter=20){

  dvariable val=1;
  int iters=0;
  while ((val*val)>1e-10&&iters<niter)
     {
     iters++;   
  
     val = (C-yield(F,B,r,K))/(gradY(F,C,B,r,K));
        
     F-=val;
     }
  
  return(F);}


inline double alpha(double r, double F) {return(r-F);}
inline double beta( double r, double K) {return(r/K);}

inline double yield(double F,double B,double r,double K){ 
     return(F/(r/K)*log(1-(r/K)*B*(1-exp((r-F)))/(r-F)));}
  
inline double stock(double F, double C, double  B, double r, double K){
    return(alpha(r,F)*B*exp(alpha(r,F))/(alpha(r,F)+beta(r,K)*B*(exp(alpha(r,F))-1)));}

double gradY(double F, double C, double B, double r, double K){

  double expr1  = r/K;
  double expr2  = F/expr1;
  double expr3  = expr1*B;
  double expr4  = (r-F);
  double expr5  = exp(expr4);
  double expr7  = expr3*(1-expr5);
  double expr9  = 1-expr7/expr4;
  double expr10 = log(expr9);

  return(-(1/expr1*expr10-expr2*((expr3*expr5/expr4+expr7/(expr4*expr4))/expr9)));}


double _nr(double F, double C, double B, 
          double r, double K,
          double tolVal=1e-10,int niter=20){
               
  double val=1;
  int iters=0;
  while ((val*val)>1e-10&&iters<niter)
     {
     iters++;   
  
     val = (C-yield(F,B,r,K))/(gradY(F,C,B,r,K));
        
     F-=val;
     }
  
  return(F);}
 

inline dvariable alpha(prevariable& r, prevariable& F) {return(r-F);}
inline dvariable beta( prevariable& r, prevariable& K) {return(r/K);}

inline dvariable yield(  prevariable& F,prevariable& B,prevariable& r,prevariable& K){ 
     return(F/(r/K)*log(1-(r/K)*B*(1-exp((r-F)))/(r-F)));}
inline dvariable stock(prevariable& F, double C, prevariable&  B, prevariable& r, prevariable& K){
    return(alpha(r,F)*B*exp(alpha(r,F))/(alpha(r,F)+beta(r,K)*B*(exp(alpha(r,F))-1)));}


dvariable gradY(prevariable& F, prevariable& C, prevariable& B, prevariable& r, prevariable& K){

  dvariable expr1  = r/K;
  dvariable expr2  = F/expr1;
  dvariable expr3  = expr1*B;
  dvariable expr4  = (r-F);//-sfabs(F-r); //r-F; 
  dvariable expr5  = exp(expr4);
  dvariable expr7  = expr3*(1-expr5);
  dvariable expr9  = 1-expr7/expr4;
  dvariable expr10 = log(expr9);

  return(-(1/expr1*expr10-expr2*((expr3*expr5/expr4+expr7/(expr4*expr4))/expr9)));}
  
  
dvariable _nr(prevariable& F, prevariable& C, prevariable& B, 
             prevariable& r, prevariable& K,
             double tolVal=1e-10,int niter=20){
                
  dvariable val=1;
  int iters=0;
  while ((val*val)>1e-10&&iters<niter)
     {
     iters++;   
     val = (C-yield(F,B,r,K))/(gradY(F,C,B,r,K));
        
     F-=val;
     }
  
  return(F);}
   
  
