struct KMF
{
  int iplane;
  int idwi;
  double KMFpar0[5][22],KMFpar1[5][22],SMTpar[5][22];
  double Fplane[5][5];
  double reso[22],mulsth[22],muls[22],meas[22],measz[22],zcoor[22],path[22];
  double eloss[22];
  double kmfchi2,kmfchi2new[22],kmfresi[22],kmferr[22];
  double smtchi2,smtchi2new[22],smtresi[22],smterr[22];
};


  
