#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export]]
RObject fwdCpp(const Rcpp::NumericMatrix& ctc, 
               const Rcpp::NumericMatrix& par,
               const Rcpp::NumericVector& eps=0.001){
  
  int nit=fmax(ctc.ncol(),par.ncol());
  NumericMatrix rtn;
  NumericMatrix b(ctc.nrow(),nit);
  
  for (int i=0; i<nit; i++){
    b(0,i)=par(1,i)*par(3,i);
    for (int j=0; j<ctc.nrow()-1; j++){
      b(j+1,i)=b(j,i)-ctc(j,i)+par(0,i)/par(2,i)*b(j,i)*(1-pow(b(j,i)/par(1,i),par(2,i)));
      if(b(j+1,i)<=eps[0]){
        b(j+1,i) = eps[0]/(2-b(j+1,i)/eps[0]);}
      }}
  
  rtn=b;
  return rtn;}

// [[Rcpp::export]]
RObject fwdTimeVaryCpp(const Rcpp::NumericVector& ctc, const Rcpp::NumericMatrix& par,
                       const Rcpp::NumericVector& eps=0.001){
  
  NumericVector rtn;
  NumericVector b(ctc.size());
  
  b(0)=par(0,1)*par(0,3);
  for (int i=0; i<ctc.size()-1; i++){
      b(i+1)=b(i)-ctc(i)+par(i,0)/par(i,2)*b(i)*(1-pow(b(i)/par(i,1),par(i,2)));
      if(b(i)<=eps[0]){
         b(i)=eps[0]/(2-b(i)/eps[0]);}
      }
  
  rtn=b;
  return rtn;}

// [[Rcpp::export]]
RObject dBdp(const Rcpp::NumericVector& ctc,
             const Rcpp::NumericVector& par,
             const Rcpp::NumericVector& eps=0.001){
  
  NumericVector rtn(1);
  NumericVector b(ctc.size());
  
  b(0)=par(1)*par(3);
  for (int i=0; i<ctc.size()-1; i++){
    b(i+1)=b(i)-ctc(i)+par(1)/par(2)*b(i)*(1-pow(b(i)/par(1),par(2)));
    if(b(i+1)<=eps[0]){
      b(i+1) = eps[0]/(2-b(i+1)/eps[0]);}
    }
  
  rtn(0)=b[ctc.size()]/exp(log((par(1)*(1/(1+par(2))))*(1/par(2))));
  return rtn;}
