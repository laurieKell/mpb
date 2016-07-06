
  // constricted likelihood
  bool flag=false;
  
  //initialise
  for (int i=1; i<=nIdx; i++){
    tmpn(i)=0;
    tmpq(i)=0;
    tmps(i)=0;
    if ((qPh(i)< -1) || (sPh(i)< -1))
       flag=true;
    }

  // calculate
  if (flag){
    for (int j=1; j<=ni; j++){
      tmpn(Idx[j])=tmpn(Idx[j])+1;
      tmpq(Idx[j])=tmpq(Idx[j])+log(I[Idx[j]])-log(0.5*(B(X[j])+B(X[j]+1)));
      }

    for (int i=1; i<=nIdx; i++)
      tmpq(i)=exp(tmpq(i)/tmpn[i]);
    
    for (int j=1; j<=ni; j++)
      tmps[Idx[j]]=+log(I[Idx[j]])-log(0.5*(B(X[j])+B(X[j]+1))*tmpq(Idx[j]));

    for (int i=1; i<=nIdx; i++){
      if (qPh[i]<-1) q[i]=tmpq(i);
      if (sPh[i]<-1) s[i]=pow((tmps(i)/tmpn[i]),.5);
      }
    }


