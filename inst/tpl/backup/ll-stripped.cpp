FUNCTION get_neglogL

  neglogL = halfnlog2pi;
  for (int j=1; j<=ni; j++)
     neglogL += log(s_)+  
                pow(log(I[j])-log(Ifit[j]),2.0)/(2*_s[Idx[j]]*s_);
     


