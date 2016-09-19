globalVariables('fuel')

ee<-function(i,e,B,price,irate,fcost,vcost,age){
  #fleet exit/entering responding to biomass y-1
  x1 = exp((-81.0017961) + 
             (  2.737885 *log(fuel)) + 
             ( -2.194898 *log(age)) + 
             (  5.7208019*log(B[i-1])) + 
             ( -2.292542 *log(price)) + 
             ( -0.4412905*log(irate)))
  
  x3 = exp((-0.3439559)  + 
             (-1.430336 * log(age)) + 
             ( 0.8241269* log(B[i-1])) + 
             ( 2.807732 * log(fuel)) + 
             (-3.481020 * log(price)) + 
             (-0.5133259* log(irate)))
  
  denom = x1 + x3
  entry = e[i-1]*(x1/1+denom)
  exit  = e[i-1]*(x3/1+denom)
  
  e[i] = round(e[i-1] + entry - exit)
  
  e}
