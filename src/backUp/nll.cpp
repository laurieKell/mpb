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
RObject nllCpp(const Rcpp::NumericVector& ctc, const Rcpp::NumericVector& idx, 
               const Rcpp::NumericMatrix& par, const Rcpp::NumericVector& eps=0.001){
  
  NumericVector rtn;
  NumericVector ll(par.nrow());
  NumericMatrix b(ctc.size()+1,par.nrow());
  NumericVector I(ctc.size());
  
  Rcout << par <<"\n";
  
  double ss=0.0, q=0.0, se=0.0, n=0.0;
  
  // loop over parameters
  for (int i=0; i<par.nrow(); i++){
    double pen=0.0;
    
    // Set population
    b(0,i)=par(i,1)*par(i,3);
    for (int j=0; j<ctc.size(); j++){
       b(j+1,i)=b(j,i)-ctc(j)+
                par(i,0)/par(i,2)*b(j,i)*(1-exp(log(b(j,i)/par(i,1))*par(i,2)));
       I(j)=0.5*(b(j,i)+b(j+1,i));
       }
    
       //if(b(j+1,i)<=eps(0)){
       //    pen+=.01*(b(j+1,i)-eps(0))*(b(j+1,i)-eps(0));
       //    b(j+1,i) = eps(0)/(2-b(j+1,i)/eps(0));}

    ss=0.0; q=0.0; se=0.0; n=0.0;
    for (int j=0; j<idx.size(); j++){
      if (!NumericVector::is_na(idx(j))){
         q+=log(idx(j))-log(I(j));
         n+=1.0;
         }
      //Rcout<<j<<"\t"<<ctc(j)<<"\t"<<idx(j)<<"\t"<<b(j,i)<<"\t"<<I(j)<<"\t"<<q<<"\n";
      }
    
    q =exp(q/n);

    for (int j=0; j<ctc.size(); j++)
      if (!NumericVector::is_na(idx(j)))
        ss+=pow(log(idx(j))-log(q*I(j)),2.0);

    se =pow(ss/n,0.5);

    //ll(i)=-(log(1/(2*3.14159265359))-ctc.size()*log(se)-ss/(2*pow(se,2)))/2;
    ll(i)= n/2*log(3.14159265359*2)
          +n/2*(log(se*se))
          +ss/(2*se*se);
  }

  //Rcout << "\n" << n <<"\t"<< q <<"\t"<< se <<"\t"<< ss <<"\n";
  
  rtn=ll;
  return rtn;}
