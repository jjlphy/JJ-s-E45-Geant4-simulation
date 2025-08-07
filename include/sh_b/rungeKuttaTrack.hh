#ifndef RUNGEKUTTATRACK_H
#define RUNGEKUTTATRACK_H 1

//#include <stdio.h>
//#include <stdlib.h>
//#include <math.h>
//#include <limits.h>
//#include <unistd.h>
#include "track.h"
//#include <switch.h>

int setInitPara(Track* aTrack, double* initPara);
int  rungeKuttaTrack(Track* aTrack);

typedef struct {
  double origin[3];
  double e1[3];
  double e2[3];
} Plane;


/* example
typedef struct TrackGroup{
  int numTracks;
  int numGene;
  int idTrack[MAX_TRACK_IN_GROUP];
  int TrackParent[MAX_TRACK_IN_GROUP];
  int id;
  int quality;
  int nVtx;
  float Vtx[3];
  float Mom[4];
} TrackGroup;
*/


 /*
 a plane on which a starting point (P) moves .
  vec{P} = vec{origin}+ s*vec{e1} + t*vec{e2}
 */


/*NTPC*/
//Hit virtualFirstHit; /*entity of hit[0](virtual hit)*/
#endif
