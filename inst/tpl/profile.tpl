
  // likelihood profile
  likeprof_number lpbnow;
  likeprof_number lpfnow;
  likeprof_number lpbthen;
  likeprof_number lpfthen;
  likeprof_number lpbnowthen;
  likeprof_number lpfnowthen;
  likeprof_number lpmsy;
  likeprof_number lpbmsy;
  likeprof_number lpfmsy;
  likeprof_number lpcmsy;
  likeprof_number lpbbmsy;
  likeprof_number lpffmsy;
  likeprof_number lpbratio;
  likeprof_number lpfratio;
  likeprof_number lpslopeb;
  likeprof_number lpslopef;
  likeprof_number lpr; 
  likeprof_number lpk;
  likeprof_number lpfr;
  likeprof_number lpbk;

  // likelihood profile 
  lpbnow  =bnow;
  lpfnow  =fnow;
  lpbnow  =bthen;
  lpfnow  =fthen;
  lpbnow  =bnowthen;
  lpfnow  =fnowthen;

  if (1!=2){
	lpmsy   =msy;
	lpbmsy  =bmsy;
	lpfmsy  =fmsy;
	lpcmsy  =cmsy;
	lpbbmsy =bbmsy;
	lpffmsy =ffmsy;
	lpbratio=bratio;
	lpfratio=fratio;
	lpslopeb=slopeb;
	lpslopef=slopef;
	lpr     =_r;
	lpk     =_k;
	lpfr    =fr;
	lpbk    =bk;
	}

  lpbnow.set_stepnumber(stepN);
  lpbnow.set_stepsize(stepSz);

  lpfnow.set_stepnumber(stepN);
  lpfnow.set_stepsize(stepSz);

  lpbthen.set_stepnumber(stepN);
  lpbthen.set_stepsize(stepSz);

  lpfnowthen.set_stepnumber(stepN);
  lpfnowthen.set_stepsize(stepSz);

  lpbnowthen.set_stepnumber(stepN);
  lpbnowthen.set_stepsize(stepSz);

  lpfnow.set_stepnumber(stepN);
  lpfnow.set_stepsize(stepSz);

  if (1!=2){

  lpmsy.set_stepnumber(stepN);
  lpmsy.set_stepsize(stepSz);

  lpfmsy.set_stepnumber(stepN);
  lpfmsy.set_stepsize(stepSz);

  lpbmsy.set_stepnumber(stepN);
  lpbmsy.set_stepsize(stepSz);

  lpcmsy.set_stepnumber(stepN);
  lpcmsy.set_stepsize(stepSz);

  lpbbmsy.set_stepnumber(stepN);
  lpbbmsy.set_stepsize(stepSz);

  lpffmsy.set_stepnumber(stepN);
  lpffmsy.set_stepsize(stepSz);

  lpbratio.set_stepnumber(stepN);
  lpbratio.set_stepsize(stepSz);

  lpfratio.set_stepnumber(stepN);
  lpfratio.set_stepsize(stepSz);

  lpslopeb.set_stepnumber(stepN);
  lpslopeb.set_stepsize(stepSz);

  lpslopef.set_stepnumber(stepN);
  lpslopef.set_stepsize(stepSz);

  lpr.set_stepnumber(stepN);
  lpr.set_stepsize(stepSz);

  lpk.set_stepnumber(stepN);
  lpk.set_stepsize(stepSz);

  lpfr.set_stepnumber(stepN);
  lpfr.set_stepsize(stepSz);

  lpbk.set_stepnumber(stepN);
  lpbk.set_stepsize(stepSz);
  }
  
