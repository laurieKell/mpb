#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() () {
  // data:
  DATA_VECTOR(ctc);
  DATA_VECTOR(idx);
  
  // parameters:
  PARAMETER(log_r); 
  PARAMETER(log_k); 
  PARAMETER(p); 
  PARAMETER(b0); 
  
  vector<Type> b(ctc.size()+1);
  
  // procedures: (transformed parameters)
  Type r=exp(log_r);
  Type k=exp(log_k);
  
  int n=ctc.size(); 
  
  Type nll=0.0; 
  Type pen=0.0;
  Type eps=0.001;
  
  b(0)=k*b0;
  for(int i=0; i<(n-1); i++){
    b(i+1) =b(i)-ctc(i)+r/p*b(i)*(1-pow(b(i)/k,p));
    if(b(i+1)<=eps){
      pen+=.01*(b(i+1)-eps)*(b(i+1)-eps);
      b(i+1) = eps/(2-b(i+1)/eps);}
   }
  
  //likelihood
  Type ss=0.0, q=0.0, se=0.0;
  for (int j=0; j<ctc.size(); j++)
    q+=log(idx(j))-log((b(j)+b(j+1))/2.0);
   q=exp(q/ctc.size());
   
  for (int j=0; j<ctc.size(); j++)
    ss  +=pow(log(idx(j))-log(b(j)*q),2.0);
   
  se =pow(ss/ctc.size(),0.5);
  nll=ctc.size()/2*log(3.14159265359*2)
     +ctc.size()/2*(log(se*se))
     +ss/(2*se*se);
     
  //for (int k=0; k<ctc.size(); k++)
  //   nll+=log((idx(k)-q*(b(k)+b(k+1))/2.0)*(idx(k)-q*(b(k)+b(k+1))/2.0));
  
  return nll+pen;}

