#ifndef MMMINUITFIT_H
#define MMMINUITFIT_H 1
#include "track.h"

typedef void ( * MinuitFUtil )();
typedef void ( * MinuitFCN )( int*,
			      double*, double*, double*, int*, MinuitFUtil );

//TpcParam* tpcPar4minuit;

int minuitFit(double* initPara, float* finalChi2, double* para,
              Track* aTrack);



#endif
