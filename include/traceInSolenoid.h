#ifndef __TRACEINSOLENOID_H__
#define __TRACEINSOLENOID_H__


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <unistd.h>

#include "readTosca.h"
#include "arraySize.h"
#include "track.h"

/*#define MAX_STEP 2000*/
#define MAX_STEP 100

int traceInSolenoid(const int n, double* parm, double* res,Track* aTrack,  int sw, TpcParam* tpcPar);
int EquMot(double s, double* x, double* xp, double* f);

int trackCrossLayer(double R0, double R1, double R2);

#endif /* __TRACEINSOLENOID_H__ */

