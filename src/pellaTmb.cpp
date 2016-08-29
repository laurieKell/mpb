#include <TMB.hpp>

// logit<-function(x,min=0,max=1){
//   x=(x-min)/(max-min)
//   x=log(x/(1-x))
//   
//   return(x)}
// 
// invLogit<-function(x,min=0,max=1){
//   
//   x=1/(1+exp(-x))
//   
//   return(x*(max-min)+min)}


template<class Type>
Type invLogit(Type x, Type min, Type max){
  x=1/(1+exp(-x));
  
  return(x*(max-min)+min);}
  
template<class Type>
Type logit(Type x, Type min, Type max){
  x=(x-min)/(max-min);
  x=log(x/(1-x));
  
  return x;}

template<class Type>
Type posfun(Type x, Type eps, Type &pen){
  pen += CppAD::CondExpLt(x,eps,Type(0.01)*pow(x-eps,2),Type(0));
  return CppAD::CondExpGe(x,eps,x,eps/(Type(2)-x/eps));}

template<class Type>
Type objective_function<Type>::operator() () {
  // data:
  DATA_VECTOR(ctc);
  DATA_MATRIX(idx);
  
  DATA_MATRIX(ctl);
  
  PARAMETER_VECTOR(par);
  
  vector<Type> b(ctc.size()+1);
  vector<Type> I(ctc.size());
  
  //transformed parameters
  Type r,k,p,b0;
  int i=0;
  if(ctl(0,0)<0) r =ctl(0,2); else if(ctl(0,0)==0) r  =par(i++); else r =invLogit(par(i++), ctl(0,1), ctl(0,3));
  if(ctl(1,0)<0) k =ctl(1,2); else if(ctl(1,0)==0) k  =par(i++); else k =invLogit(par(i++), ctl(1,1), ctl(1,3));
  if(ctl(2,0)<0) p =ctl(2,2); else if(ctl(2,0)==0) p  =par(i++); else p =invLogit(par(i++) ,ctl(2,1), ctl(2,3));
  if(ctl(3,0)<0) b0=ctl(3,2); else if(ctl(3,0)==0) b0 =par(i++); else b0=invLogit(par(i++), ctl(3,1), ctl(3,3));
        
  Type nll=0.0; 
  Type pen=0.0;
  Type eps=0.00001;
  
  vector<Type> n( idx.cols()+1);
  vector<Type> q( idx.cols()+1);
  vector<Type> ss(idx.cols()+1);
  vector<Type> se(idx.cols()+1);
  vector<Type> ll(idx.cols()+1);
  Rcout << r <<"\t"<< k <<"\t"<< p <<"\t"<< b0 <<"\n";
  
  //Population dynamics
  b(0)=k*b0;
  for(int i=0; i<ctc.size(); i++){
    Type now=b(i)-ctc(i);//+r/p*b(i)*(1-exp(log(b(i)/k)*p));
    now=posfun(now, eps, pen);
    b(i+1)=now+r/p*b(i)*(1-exp(log(b(i)/k)*p));
    I(i)  =(b(i)+b(i+1))/2.0;
    }
  
  //likelihoods by Index
  int ncol=idx.cols();
  if (ncol>ctc.size()) ncol=ctc.size();
  for (int j=0; j<ncol; j++){
	  ss(j)=0.0; q(j)=0.0, se(j)=0.0, n(j)=0.0;
	
	  //cathability
	  for (int i=0; i<ctc.size(); i++){
	    if (idx(i,j)>0.0){
	      q(j)+=log(idx(i,j))-log(I(i));
	      n(j)+=1.0;
	      }
	    //Rcout<<i<<"\t"<<ctc(i)<<"\t"<<idx(i,j)<<"\t"<<b(i)<<"\t"<<I(i)<<"\t"<<q(j)<<"\n";
	    }
	  q(j)=exp(q(j)/n(j));
	  
      //sum ogf squares
	  for (int i=0; i<ctc.size(); i++)
	    if (idx(i,j)>0.0){
		    ss(j)+=(log(idx(i,j))-log(q(j)*I(i)))*
			       (log(idx(i,j))-log(q(j)*I(i)));}
	    
	  // SEs
	  se(j) =exp(log(ss(j)/n(j))*0.5);
	  
	  //Negative likelihood
	  nll+=n(j)/2*log(3.14159265359*2)
		    +n(j)/2*(log(se(j)*se(j)))
		    +ss(j)/(2*se(j)*se(j));
		
		Rcout <<"\n"<<  n(j) <<"\t"<< q(j) <<"\t"<< se(j) <<"\t"<< ss(j) <<"\t" << nll <<"\t" << pen <<"\n\n";
    }
  
  return -nll-pen;}
