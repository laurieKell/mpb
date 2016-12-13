
  cnow =C[nc];
  fnow =C[nc]/B[nc];
  bnow =B[nc];
  fthen =C[nc-ref[4]+1]/B[nc-ref[4]+1];
  bthen =B[nc-ref[4]+1];
  fnowthen  =fnow/fthen;
  bnowthen =bnow/bthen;
  
  msy  =r*k*pow(1/(1+p),1/p+1);
  bmsy =k*pow((1/(1+p)),(1/p));
  fmsy =msy/bmsy;
  cmsy =C[nc]	/msy;
  bbmsy=bnow/bmsy;
  ffmsy=fnow/fmsy;
  
  bk=bnow/k;
  fr=fnow/r;
  
   bratio=0.0;
   fratio=0.0;
  _bratio=0.0;
  _fratio=0.0;
  
  for (int i=nc; i>nc-ref[2]; i--){
    bratio+=B[i];
    fratio+=C[i]/B[i];
    
    _bratio+=B[i-ref[3]];
    _fratio+=C[i-ref[3]]/B[i-ref[3]];
    }
 bratio=bratio/_bratio;
 fratio=fratio/_fratio;

  xy=0.0; 
  x =0.0;
   y=0.0; 
  xx=0.0; 
  for (int i=nc; i>nc-ref[1]; i--){
    x +=i;
    xx+=i*i;  
    y +=B[i];  
    xy+=i*B[i];  
    }
  slopeb = -(ref[1]*xy - x*y)/(ref[1]*xx - x*x);
   
  xy=0.0; 
  x =0.0;
   y=0.0; 
  xx=0.0; 
  for (int i=nc; i>nc-ref[1]; i--){
    x +=i;
    xx+=i*i;  
    y +=C[i]/B[i];  
    xy+=i*C[i]/B[i];  
    }
  slopef = -(ref[1]*xy - x*y)/(ref[1]*xx - x*x);
   
