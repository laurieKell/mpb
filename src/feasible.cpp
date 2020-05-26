#include <Rcpp.h>



using namespace Rcpp;

// [[Rcpp::export]]
RObject feasibleCpp(const Rcpp::NumericVector& ctc, const Rcpp::NumericMatrix& par,
                    const Rcpp::NumericVector& eps=0.001){
  
  NumericVector rtn(par.nrow());
  NumericMatrix b(ctc.size(),par.nrow());
  
  for (int i=0; i<par.nrow(); i++){
    rtn(i)=1;
    b(0,i)=par(i,1)*par(i,3);
    for (int j=0; j<ctc.size()-1; j++){
      b(j+1,i)=b(j,i)-ctc(j)+par(i,0)/par(i,2)*b(j,i)*(1-pow(b(j,i)/par(i,1),par(i,2)));
      if(b(j+1,i)<=eps[0]){
        b(j+1,i) = eps[0]/(2-b(j+1,i)/eps[0]);
        rtn(i)=0;}
    }}
  
  return rtn;}
