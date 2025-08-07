#ifndef __TRACK_H__
#define __TRACK_H__
//#include <arraySize.h>
//#include <cluster.h>
//#include <hit.h>
//#include "minuitFit.h"
//#include <switch.h>

#define MAX_HIT_IN_TRACK 500
#define MAX_ITERATION 1000
#define NUM_PARA_RK 5
#define MAX_TRACK_IN_GROUP 30
#define MAX_TRACK 30

typedef struct Track{
  int igroup;
  int trkQual;
  int numHits;
  //  int numSectors; /*NTPC*/
  int numLayers;  /*total number NTPC*/
  //  int numLayersinSec[NUM_SECTOR]; /*NTPC  in each sector*/
  int sector[MAX_HIT_IN_TRACK+1]; /*NTPC*/
  int lay[MAX_HIT_IN_TRACK+1];
  /////  Hit* hit[MAX_HIT_IN_TRACK+1];
  int charge;
  int numStep;  /* # of steps in RK trace */
  double totalLength;
  double totalLengthTOF; /*NTPC TOF*/
  int CrossOuter;       /*NTPC TOF crossing outer scinti*/
  double meanAdc;
  double chi2Pad;
  double chi2Z;
  double chi2Prob;
  double chi2;
  double x[MAX_HIT_IN_TRACK+1][3];
  double uv[MAX_HIT_IN_TRACK+1][2];
  double xlocal[MAX_HIT_IN_TRACK+1][3]; /*NTPC*/
  double xEvaluate[MAX_HIT_IN_TRACK+1][3]; /* a point on the tarck to calcurate residual */
  double err[MAX_HIT_IN_TRACK+1][3]; /* err[0][]=errX, err[1][] = errY, err[1][] = errZ */  
  double radius;
  double center[2];
  double arcLen[MAX_HIT_IN_TRACK+1];
  double resPad[MAX_HIT_IN_TRACK+1];
  double resPady[MAX_HIT_IN_TRACK+1]; /*NTPC*/
  double resZ[MAX_HIT_IN_TRACK+1];
  double res[MAX_HIT_IN_TRACK+1];
  int idwi[MAX_HIT_IN_TRACK+1];
  double rKresXYZ[MAX_HIT_IN_TRACK+1][4];  /* 0:X 0:Y 0:Z  !!! in LOCAL sector !!*/
  /*1st  double rKresRPhiZ[MAX_HIT_IN_TRACK+1][4]; */ /* [0] : Res-phi, [1] : Res-Z  */
  double rKNumIter;
  double rKChi2[MAX_ITERATION];
  double rKInitPara[NUM_PARA_RK];
  double rKFinalPara[NUM_PARA_RK];
  double parRZFit[2];
  double mom[4];
  double xOnTrack[3];
  double resVirtual[4];  /* resPhi, resZ on the virtual plane    */
  double initRes[MAX_HIT_IN_TRACK+1][3];  /* [0] : Res-phi, [1] : Res-Z  */
  /*test 3*/
  double RKPFinal[3];

  /**/
  double ydist[MAX_HIT_IN_TRACK+1];
  int whilet[MAX_HIT_IN_TRACK+1];
  double x_org[MAX_HIT_IN_TRACK+1][3];
  double xlocal_org[MAX_HIT_IN_TRACK+1][3];
  double phi_local[MAX_HIT_IN_TRACK+1];
  int nout;
  int ngood;
  int ibad[MAX_HIT_IN_TRACK+1];
  int zbad[MAX_HIT_IN_TRACK+1];
} Track;


typedef struct TrackGroup{
  int numTracks;
  int numGene;
  int idTrack[MAX_TRACK_IN_GROUP];
  int TrackParent[MAX_TRACK_IN_GROUP];
  int id;
  int quality;
  int nVtx;
  double Vtx[3];
  double Mom[4];
} TrackGroup;



/*NTPC hit[sector]*/
/*
int findTrackGroup(Hit hits[][NUM_LAY][MAX_SELECTED_HITS], Track* tracks,
		   TrackGroup* trackGroups, int MODE, TpcParam* tpcPar,
		   Switch* sw, Wire* pwire);

int getInitialMom(Track* tracks,Switch *sw);

int setSpatialErr(Track* aTrack);

int minuitFit( double* initPara, double* chi2, double* finalPara,
               Track* aTrack, int sector, int lay,Switch *sw);


int readpar_trk_(); 
*/


  /*
  ** x = r*cos(phi)+center[0]
  ** y = r*sin(phi)+center[1]
  ** z = lambda*phi + x[2]
  */
typedef struct Helix {
  int charge;
  int signBz;

  double r;
  double center[2];
  double x[3];      /* (x,y,z) global coor. on the helix where phi=0 */
  double xlocal[3]; /* (x,y,z) local coor. on the helix where phi=0  NTPC*/
  double tanLambda;
} Helix;


#endif /* __TRACK_H__ */
