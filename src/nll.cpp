#include <Rcpp.h>

using namespace Rcpp;

double posfun2(const double &x, const double eps, const double& _pen){
  double& pen=(double&)_pen;
  
  if (x>=eps) {
    return x;}
  else{
    double y=eps-x;
    double tmp=y/eps;
    double tmp2=tmp*tmp;
    double tmp3=tmp2*tmp;
    pen+=.01*tmp3*tmp3*tmp3;
    return eps/(1.0+tmp+tmp2+tmp3);}
  }

//https://en.wikipedia.org/wiki/Maximum_likelihood

// [[Rcpp::export]]
RObject nll(const Rcpp::NumericVector& ctc, const Rcpp::NumericVector& idx, 
            const Rcpp::NumericMatrix& par, const Rcpp::NumericVector& eps=0.001){
  
  NumericVector rtn;
  NumericVector ll(par.nrow());
  NumericMatrix b(ctc.size(),par.nrow());
  
  for (int i=0; i<par.nrow(); i++){
    double pen=0.0;
    b(0,i)=par(i,1)*par(i,3);
    for (int j=0; j<ctc.size()-1; j++){
       b(j+1,i)=b(j,i)-ctc(j)+
                 par(i,0)/par(i,2)*b(j,i)*(1-pow(b(j,i)/par(i,1),par(i,2)));

       if(b(j+1,i)<=eps[0]){
          pen+=.01*(b(j+1,i)-eps[0])*(b(j+1,i)-eps[0]);
          b(j+1,i) = eps[0]/(2-b(j+1,i)/eps[0]);}
         
     double ss=0.0, q=0.0, se=0.0;
     for (int j=0; j<ctc.size()-1; j++)
       q +=log(idx(j))-log(b(j,i));
     q =exp(q/ctc.size());
     
    for (int j=0; j<ctc.size(); j++)
       ss  +=pow(log(idx(j))-log(b(j,i)*q),2.0);
     
    se =pow(ss/ctc.size(),0.5);
    //ll(i)=-(log(1/(2*3.14159265359))-ctc.size()*log(se)-ss/(2*pow(se,2)))/2;
    ll(i)=ctc.size()/2*log(3.14159265359*2)
          +ctc.size()/2*(log(se*se))
          +ss/(2*se*se)
          +pen;
     }}
  
  rtn=ll;
  return rtn;}
