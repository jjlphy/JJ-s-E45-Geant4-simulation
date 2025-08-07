#include <math.h>
#include <float.h>
 
#include "traceInSolenoid.h"
#include "arraySize.h"
#include "tpcGeometry.h"
#include "libTpc.h"
#include "readTosca.h"
#include "minuitFit.h"
#include "rungeKuttaTrack.h"
#include "driftTable.h"
#include <error.h>
#include <ntTpc.h>

#define SQ(a) ((a)*(a))

#define DIM 3
#define SDIM (6*(DIM))

double RKCharge;          /* Charge for the tracked particle   */
double RKTransMomentum;  /* transverse Momentum for the tracked particle */
double RKMomentum;  /* Momentum for the tracked particle */

extern double RKResXYZ[MAX_HIT_IN_TRACK][4];
/* [0] : res-X, [1] : res-Y,[2] : res-Z */
/*extern double RKResRPhiZ[MAX_HIT_IN_TRACK][4];*/
/* [0] : res-phi, [1] : res-phi,[2] : res-Z */

extern Plane RKVirtualPlane;
extern double RKInitRes[MAX_HIT_IN_TRACK][3];

extern int RKNumIter;
extern int eveNum;

extern float B_FIELD;

int storeFinalRes(double res[MAX_HIT_IN_TRACK][4],int numHits);
int storeInitialRes(double res[MAX_HIT_IN_TRACK][4],int numHits);

double calcResidualLine(double *res,Track *aTrack,int ii
			,double *curXlocal,double *curXPlocal,double step,int sw);

typedef void (* SubDRKNYS)( double* x, double* y, double* yp, double* f);

void drknys_( int* n, double* h, double* x, double* y, double* yp,
              SubDRKNYS sub, double* w );
extern void equf2_( double* x, double* y, double* yp, double* f );
extern void equmot_(double* A1, double* A2,double* A3, double* A4) {
  EquMot( *A1, A2, A3, A4 );
  return;
}
/* equmot is not declared in other files. Then, next declaration is probably ok.
void equmot_(double* A1, double* A2,double* A3, double* A4) {
  EquMot( *A1, A2, A3, A4 );
  return;
}
*/

extern int chi2prb_(float* chi2, int* ndf, float* prb);

double calcStep(double *Xlocal,double *XPlocal,Track *aTrack,int iHit);

extern int setSpatialErrRK(Track* aTrack,int iHit);
double chiz[MAX_HIT_IN_TRACK][5];
double dyCor[9][10][5];

/*********************************************
 *********************************************/
int traceInSolenoid(const int nPara, double* param, double* chiSq, Track* aTrack, int sw, TpcParam* tpcPar)
{
  int i,ii;

  int numHits = aTrack->numHits;
  aTrack->nout =0;
  const int dim = DIM;

  double work[SDIM];
  double curX[DIM],curXP[DIM];
  double preX[DIM],preXP[DIM];

  double curXlocal[DIM],preXlocal[DIM];
  double curXPlocal[DIM],preXPlocal[DIM];

  double curPathLength = 0.0;
  double prevPathLength = 0.0;

  double start[DIM];

  double step = 5.;/* default*/ /* mm */
  /*double step = 20.;*/ /* mm */

  int numStep = 0;

  FILE *fp2 = NULL;

  double residual = 0.0;
  double chi[MAX_HIT_IN_TRACK];
  double tmpRes[MAX_HIT_IN_TRACK][4]={{0.}};
  double res[MAX_HIT_IN_TRACK][4]   ={{0.}};

  /*NTPC SelectHits >>*/
  int sele[MAX_HIT_IN_TRACK];
  int whilet=0;
  double RKX[MAX_STEP+1][3];
  double RKXlocal[MAX_STEP+1][3];
  double RKXP[MAX_STEP+1][3];
  int wbreak=0;
  /*<< */

  double ydist;
  
  /* 1st
     static double matrix[MAX_HIT_IN_TRACK][DIM][DIM] = {{{0.}}};
  */
  int foundEstimatePoint[MAX_HIT_IN_TRACK];

  int dumpFile = 0;

  int exitFirstIter=0;

  double t0,t1;
  double bVec[DIM];
  Helix helix;
  int flag;

  /*NTPC test>>*/
  int stepIDcross[MAX_HIT_IN_TRACK];
  for( i=0 ; i<MAX_HIT_IN_TRACK ; i++ ){
    stepIDcross[i]=-10;
    sele[i]=0;
  }
  for(i=0;i<MAX_STEP+1;i++){
    RKX[i][0]=0.;
    RKX[i][1]=0.;
    RKX[i][2]=0.;
    RKXP[i][0]=0.;
    RKXP[i][1]=0.;
    RKXP[i][2]=0.;
  }
  /*<< */


  if(exitFirstIter == 1)
    sw=1;

  curPathLength = 0.0;
  prevPathLength = 0.0;

  for( i=0 ; i < MAX_HIT_IN_TRACK; i++){
    tmpRes[i][0]= tmpRes[i][1]= tmpRes[i][2]= tmpRes[i][3]=
      res[i][0]= res[i][1]= res[i][2]= res[i][3]=
      RKResXYZ[i][0] = RKResXYZ[i][1] = RKResXYZ[i][2] = RKResXYZ[i][3] =
      aTrack->xEvaluate[i][0] = aTrack->xEvaluate[i][1] = aTrack->xEvaluate[i][2] = 
      10000.;

    /*chi[i] = DBL_MAX-1;*/
    chi[i] = 1000000000000000000000.;/*2010Aug06*/
    foundEstimatePoint[i] = -1;

  }

  if( numHits >= MAX_HIT_IN_TRACK){
    reportErr(0,"traceInSolenoid RKNumHit %d >=MAX_HIT_IN_TRACK\n",
           numHits);
    exit(2);
  }
  

  if( dumpFile == 1){
    if(sw == 1){
      if((fp2 = fopen("dat2","w")) == NULL){
        fprintf(stderr,"Can't open dat2\n");
        exit(1);
      }
    }
  }
  
  if(fabs(param[2])>100) param[2]= 100;

  isNanInf(param, NUM_PARA_RK, "traceInSolenoid -1", "param");

  
  param[2] = fabs(param[2]);

  RKMomentum = param[2];
  RKCharge   = (double)aTrack->charge;
  RKTransMomentum = fabs(param[2]*sin(param[3]));

  
  for( i=0; i<3; i++){
    start[i] = curX[i] = RKVirtualPlane.origin[i]
      +param[0]*RKVirtualPlane.e1[i]+param[1]*RKVirtualPlane.e2[i];

    preX[i] = curX[i];
  }
  
  curXP[0] = sin(param[3])*cos(param[4]);
  curXP[1] = sin(param[3])*sin(param[4]);
  curXP[2] = cos(param[3]);


  for(i=0;i<3;i++){
    preXP[i] = curXP[i];
  }

  if( isNanInf(curX, 3, "traceInSolenoid 0", "curX")== 1
      || isNanInf(curXP, 3, "traceInSolenoid 0", "curXP")== 1
      || isInBMap(curX, "traceInSolenoid 0",sw) == 0
      || isInTpcRKtrk(curX,&(tpcPar4minuit->align),No)==No/*2010Aug05*/
      || RKTransMomentum > 100./*2010Aug06 100GeV/c*/
      ){
    double chiSqorg;
    fflush(stdout);
    /*
     *chiSq = 0;
     for(i=0; i < numHits; i++)
     *chiSq += chi[i];
     */
    chiSqorg = *chiSq;
    *chiSq = 1000000000000000000000.;/*2010Aug05*/
    
    if(sw == 2 ){
      storeInitialRes(res, aTrack->numHits); 
    }

    if( sw == 1){

      storeFinalRes(res, aTrack->numHits);   

      reportErr(0,"error at RK start point x:%lf y:%lf z:%lf Pt:%lf radius:%lf\n"
		,curX[0],curX[1],curX[2],RKTransMomentum,aTrack->radius
		);
      
    }
    
    return numStep;
  }

  if(dumpFile == 1 ){
    if(sw==1){
      fprintf(fp2,"%f %f %f %f\n",curX[0],curX[1],curX[2],
              sqrt(curX[0]*curX[0]+curX[1]*curX[1]));

      fflush(fp2);
    }
  }

  /*
  while ( isInTpcRKtrk(curX,&(tpcPar4minuit->align),No)==Yes && numStep < MAX_STEP){
    whilet++;
  */
  for(ii=0;ii<numHits;ii++){
    
    for( i=0; i<3; i++){
      preX[i] = curX[i];
      preXP[i] = curXP[i];
      RKX[ii][i]=curX[i];
      RKXP[ii][i]=curXP[i];
    }
    prevPathLength = curPathLength;

    SectorRotation(aTrack->sector[ii],curXlocal,curX,3,G2L);
    SectorRotation(aTrack->sector[ii],curXPlocal,curXP,3,G2L);

    /**/
    if(curXlocal[1]>aTrack->xlocal[ii][1]){
      SectorRotation(aTrack->sector[ii],curXPlocal,curXP,3,G2L);
      chi[ii] = calcResidualLine(res[ii],aTrack,ii
				 ,curXlocal,curXPlocal,step,sw);
      continue;
    }
    /**/

    step = calcStep(curXlocal,curXPlocal,aTrack,ii);
    
    flag = 0;
    whilet = 0;
    ydist = 10000.;
    while(flag!=1){
      double norm;
      
      if( isNanInf(curX, 3, "traceInSolenoid 2.0", "curX")== 1
	  || isInBMap(curX, "traceInSolenoid 2.0",sw) == 0
	  || ( numStep > MAX_STEP) ){
	wbreak=1;
	break;
      }

      drknys_( ( int* )( & dim ), & step, &curPathLength,
	       curX, curXP, equf2_, work );

      if(sw==1){
	norm = sqrt(SQ(curXP[0])+SQ(curXP[1])+SQ(curXP[2]));
	if( fabs(norm-1.)>0.01 ){
	  reportErr(1,"curXP far from 1. %lf , step=%lf x:%lf y:%lf z:%lf sw:%d\n",norm,step,curX[0],curX[1],curX[2],sw);
	}
      }
      
      SectorRotation(aTrack->sector[ii],curXlocal,curX,3,G2L);
      SectorRotation(aTrack->sector[ii],curXPlocal,curXP,3,G2L);

      if( fabs(curXlocal[1]-aTrack->xlocal[ii][1])<0.5 ){
	/*if( fabs(curXlocal[1]-aTrack->xlocal[ii][1])<0.1 ){*/
	ydist = curXlocal[1]-aTrack->xlocal[ii][1];
	flag = 1;
      }else if( curXlocal[1]<aTrack->xlocal[ii][1] ){

	for( i=0; i<3; i++){
	  preX[i] = curX[i];
	  preXP[i] = curXP[i];
	}
	prevPathLength = curPathLength;

	step = calcStep(curXlocal,curXPlocal,aTrack,ii);
	
      }else if( curXlocal[1]>aTrack->xlocal[ii][1] ){
	step *= 0.5;
	for(i=0;i<3;i++){
	  curX[i]  = preX[i];
	  curXP[i] = preXP[i];
	}
      }    

      whilet++;
      
    } /*while*/


    if(wbreak == 1) {
      if(sw==1){
	reportErr(0,"wbreak numStep=%d, x=%lf, y=%lf, z=%lf step=%lf:hitx=%lf y=%lf z=%lf   hitxp=%lf yp=%lf zp=%lf\n"
		  ,numStep,curX[0],curX[1],curX[2],step
		  ,aTrack->x[ii][0],aTrack->x[ii][1],aTrack->x[ii][2]
		  ,curXP[0],curXP[1],curXP[2]
		  );

      }
      break;
    }

    chi[ii] = calcResidualLine(res[ii],aTrack,ii
			       ,curXlocal,curXPlocal,step,sw);

    if(sw==1){
      aTrack->ydist[ii]=ydist;
      aTrack->whilet[ii]=whilet;
    }
    
    numStep++;
  }/*for*/

  /*end point of RK*/
  for(i=0;i<3;i++){
    RKX[ii][i]=curX[i];
    RKXP[ii][i]=curXP[i];
  }
  
  *chiSq = 0;
  float maxchi=0;
  int   bid   =-1;
  for(i=1; i < numHits; i++){
    if(aTrack->ibad[i-1]>0){
      if(chi[i]>maxchi){
	maxchi = chi[i];
	bid    = i;
      }
    }
  }

  for(i=1; i < numHits; i++){
    if(sw!=1){
      if(aTrack->ibad[i-1]>0){
	*chiSq += chi[i];
      }
    }else if(sw==1){
      if(maxchi<25){
	if(aTrack->ibad[i-1]>0){
	  *chiSq += chi[i];
	}
      }else if(i!=bid&&aTrack->ibad[i-1]>0){
	*chiSq += chi[i];
      }else if(i==bid){
	aTrack->nout++;
	aTrack->ibad[i-1] =-10;
      }
    }
  }
      /*
      if(chi[i]<25){
	*chiSq += chi[i];
      }else{
	aTrack->nout++;
	aTrack->ibad[i-1] =-10;
      }
      */
  
  if(RKNumIter>50){
    int  ndf = (2*(aTrack->numHits-aTrack->nout)-5);
    float prob,chi2;
    chi2 = *chiSq;
    if(ndf>0){
      chi2prb_( &chi2, &ndf, &prob);
      if(prob<0.02){
	float zmax=0;
	int   zid;
	for(i=1; i < numHits; i++){
	  if(fabs(chiz[i][0])>zmax){
	    zmax = fabs(chiz[i][0]);
	    zid  = i;
	  }
	}
	
	if(aTrack->zbad[zid-1]>0&&chiz[zid][0]<-2&&param[3]<0.9){
	  int pad   = aTrack->hit[zid]->peakTdc;
	  if(aTrack->hit[zid]->wirepulse[2]!=NULL&&aTrack->hit[zid]->wirepulse[3]!=NULL){
	    int wpad3 = aTrack->hit[zid]->wirepulse[2]->peakTdc;
	    int wpad4 = aTrack->hit[zid]->wirepulse[3]->peakTdc;
	    //if(fabs(wpad4-pad-1)<3.5&&wpad3-pad<1){
	    if(fabs(wpad4-pad-1)<3.5){
	      for(i=1; i < numHits; i++){
		float tdrift = ((aTrack->hit[i]->peakTdc - 5.5 - tpcPar->fadcT0Global)*FADC_CLOCK/1000.0);
		float newz   = tdrift * tpcPar->vDrift + tpcPar->align.zShield ;
		aTrack->hit[i]->x[2]      = tdrift * tpcPar->vDrift + tpcPar->align.zShield;
		aTrack->hit[i]->xlocal[2] = tdrift * tpcPar->vDrift + tpcPar->align.zShield;
		aTrack->x[i][2]           = tdrift * tpcPar->vDrift + tpcPar->align.zShield;
		aTrack->xlocal[i][2]      = tdrift * tpcPar->vDrift + tpcPar->align.zShield;
		aTrack->zbad[i-1] = -10;
	      }
	      
	      /*
	      float tdrift = ((aTrack->hit[zid]->wirepulse[2]->cfTdc - tpcPar->fadcT0Global)*FADC_CLOCK/1000.0);
	      float newz   = tdrift * tpcPar->vDrift + tpcPar->align.zShield ;
	      float newsz  = (newz - chiz[zid][1])/chiz[zid][2];
	      if(fabs(newsz)<fabs(chiz[zid][0])){
		aTrack->hit[zid]->x[2]      = tdrift * tpcPar->vDrift + tpcPar->align.zShield;
		aTrack->hit[zid]->xlocal[2] = tdrift * tpcPar->vDrift + tpcPar->align.zShield;
		aTrack->x[zid][2]           = tdrift * tpcPar->vDrift + tpcPar->align.zShield;
		aTrack->xlocal[zid][2]      = tdrift * tpcPar->vDrift + tpcPar->align.zShield;
		aTrack->zbad[zid-1] = -10;
	      }
	      */
	    }//padcut
	  }
	}//hitcut
      }//hit loop
    }
  }
  
    /*
  float ndf = (float)(2*(aTrack->numHits-aTrack->nout)-5);
  if(ndf>0){
    *chiSq = *chiSq/ndf;
  }else{
    *chiSq = 1000000000000000000000.;
  }
  */

  if(wbreak==1 && sw==1){
    *chiSq = 1000000000000000000000.;
    aTrack->trkQual = 0;
  }

  
  if(dumpFile == 1 ){
    if(sw==1){
      fclose(fp2);
    }
  }

  if(sw == 2 ){
    extern int trkNum; 
    storeInitialRes(res, aTrack->numHits);
  }


  if( sw == 1){
    extern int trkNum; 

    storeFinalRes(res, aTrack->numHits);
    
  }

  if(exitFirstIter == 1){
    reportErr(0,"chi2 %f\n",*chiSq);
    exit(0);
  }

  
  /*test 3*/
  if(sw==1){
    double theta;
    double phi;
    double xpt;
    
    theta = acos(curXP[2]);
    xpt   = sqrt(SQ(curXP[0])+SQ(curXP[1]));
    phi   = atan2(curXP[1],curXP[0]);
    
    aTrack->RKPFinal[0] = RKMomentum*sin(theta)*cos(phi);
    aTrack->RKPFinal[1] = RKMomentum*sin(theta)*sin(phi);
    aTrack->RKPFinal[2] = RKMomentum*cos(theta);

  }
  /*******/

  
  return numStep;
}

/*
 * EquMot
 * Func.     : Equation of Motion in TPC Solenoid magnet
 *             df2/dX2 value set for DRKSTP :: One Step
 *             s: path length, *x:x,y,z, *xp:dx/ds,dy/ds,dz/ds
 *
 * return    : 0 : OK
 *             1 : globalMom = 0 :: exit();
 *             2 : Cannot Get MagField
 *
 */

int EquMot(double s, double* x, double* xp, double* f) {
  /*  Hep3VectorD posVec, magVec; */
  double magVec[3];
  double Const;
  extern int eveNum;

  /*test*/ int sw=-1;/*sw is dummy.*/

  f[0] = f[1] = f[2] = 0.;

  if( RKMomentum == 0 ) {
    reportErr(0,"###### in EquMot Div0 ######\n");
    exit(1);
  }

  if(isNanInf(x, 3, "EquMot", "x") == 1
     /*     || isInBMap(x, "EquMot",sw) == 0)*/
     || isInBMap(x, "EquMot",sw) == 0)
    return 2;

  Const = 0.299792458*1.0E-3/RKMomentum/RKCharge;
  /* ---------- "Get Magnetic Field" ----- */

  if(getBField( x, magVec) == -1){
    reportErr(0,"ev:%d EquMot: Magnetic field out of range %f %f %f \n",eveNum,x[0],x[1],x[2]);
    return 2;
  }

  /* ---------- "Lorentz Force " -------- */

  f[0] = Const * (xp[1]*magVec[2] - xp[2]*magVec[1]);  /* d2x/ds2 */
  f[1] = Const * (xp[2]*magVec[0] - xp[0]*magVec[2]);  /* d2y/ds2 */
  f[2] = Const * (xp[0]*magVec[1] - xp[1]*magVec[0]);  /* d2z/ds2 */

  return 0;
}

/*******************************************
 *******************************************/
int storeFinalRes(double res[MAX_HIT_IN_TRACK][4], int numHits){

  int iHit;

  for(iHit=0; iHit < numHits; iHit++){
    RKResXYZ[iHit][0] = res[iHit][0];
    RKResXYZ[iHit][1] = res[iHit][1];
    RKResXYZ[iHit][2] = res[iHit][2];
    RKResXYZ[iHit][3] = res[iHit][3];
  }

  return 0;

}
/*******************************************
 *******************************************/


/*******************************************
 *******************************************/
int storeInitialRes(double res[MAX_HIT_IN_TRACK][4], int numHits){

  int iHit;

  for(iHit=0; iHit < numHits; iHit++){
    RKInitRes[iHit][0] = res[iHit][0];
    RKInitRes[iHit][1] = res[iHit][1];
    RKInitRes[iHit][2] = res[iHit][2];
  }

  return 0;

}
/*******************************************
 *******************************************/


/*************************************
 *************************************
static  double f(double t, double dX, double dY,
                 double r,double tanLam, double eX2, double eY2,double eZ2){

  return dX*sin(t)/eX2-dY*cos(t)/eY2+(1./eX2-1./eY2)*r*sin(t)*cos(t)-r*tanLam*tanLam*t/eZ2;

}
*************************************
 *************************************/


/******************************
 ******************************/
double calcResidualLine(double *res, Track *aTrack, int iHit,
			double *Xlocal,double *XPlocal,double step,int sw){
  int i;
  double chi;
  double p1,q1,r1;
  double x0,y0,z0;
  double y1;
  double xEvlocal[3];
  double xcor[3];
  
  x0=Xlocal[0];
  y0=Xlocal[1];
  z0=Xlocal[2];
  y1=aTrack->xlocal[iHit][1];
  p1 = XPlocal[0]/XPlocal[1];
  q1 = XPlocal[2]/XPlocal[1];
  r1 = XPlocal[2]/sqrt(SQ((XPlocal[0]))+SQ((XPlocal[1])));
  
  xEvlocal[0] = p1*(y1-y0)+x0;
  xEvlocal[1] = y1;
  xEvlocal[2] = q1*(y1-y0)+z0;

  SectorRotation(aTrack->sector[iHit]
		 ,xEvlocal,aTrack->xEvaluate[iHit],3,L2G);

  aTrack->hit[iHit]->tanphiRK = p1;
  aTrack->hit[iHit]->tandipyzRK = q1;
  aTrack->hit[iHit]->tandipRK = r1;
  if(sw==2){
    setSpatialErrRK(aTrack,iHit);
  }
  dyHitCor(aTrack,iHit,xcor);
  
  chi = 0.;
  for(i=0;i<3;i++){
    res[i] = xcor[i] - xEvlocal[i];
    chi += SQ(res[i])/SQ(aTrack->err[iHit][i]);
    if(i==2) {
      chiz[iHit][0] = res[i]/aTrack->err[iHit][i];
      chiz[iHit][1] = xEvlocal[i];
      chiz[iHit][2] = aTrack->err[iHit][i];
      chiz[iHit][3] = - aTrack->xlocal[iHit][2]  + xcor[i];
    }
    /*
    if(i<2) {
      chi += SQ(res[i])/SQ(0.2);
    }else{
      chi += SQ(res[i])/SQ(0.4);
    }
    */
  }

  return chi;
}
/************************
 ************************/


/************************
 ************************/
double calcStep(double *Xlocal,double *XPlocal,Track *aTrack,int iHit){
  double step;
  int i;
  double dy,yp;
  double step_max;
  
  step_max=40.;
  
  dy= fabs(Xlocal[1]-aTrack->xlocal[iHit][1]);
  yp = fabs(XPlocal[1]);
  if(yp<DBL_EPSILON){
    step = dy;
  }else{
    step = dy/yp;
  }

  
  if(step>step_max){
    step = step_max;
  }

  return step;
}
/************************
 ************************/


/************************
 ************************/
int readdycorpar(){
  FILE *fp;
  int i;
  int ilay,idip;
  double dlay,ddip;
  double p[5];
  char line[500];
  
  fp = fopen("dydep.par","r");

  while( fgets(line,sizeof(line),fp)!=NULL ){
    sscanf(line,"%lf %lf %lf %lf %lf %lf %lf"
	   ,&dlay,&ddip,&p[0],&p[1],&p[2],&p[3],&p[4]);
    ilay = (int)dlay;
    idip = (int)ddip;

    for(i=0;i<5;i++){
      dyCor[ilay][idip][i] = p[i];
    }

  }
  
  fclose(fp);

  return 0;
}
/************************
 ************************/


/************************
 ************************/
int dyHitCor(Track *aTrack,int iHit,double *xcor){
  int i;
  double dip,phi;
  int ilay,idip;
  int ndip;
  double mindip,maxdip;
  
  ndip = 10;
  mindip = -35.;
  maxdip = 65.;
  
  xcor[0]=aTrack->xlocal[iHit][0];
  xcor[1]=aTrack->xlocal[iHit][1];
  xcor[2]=aTrack->xlocal[iHit][2];

  ilay = aTrack->lay[iHit];
  if(ilay<0) return 0; /*virtual hit*/
  
  dip = atan(aTrack->hit[iHit]->tandipRK)*180./M_PI;
  /*phi = atan(aTrack->hit[iHit]->tanphiRK)*180./M_PI;*/
  phi  = aTrack->phi_local[iHit];

  if(dip<mindip){
    idip = 0;
  }else if(dip>maxdip){
    idip = ndip-1;
  }else{
    idip = (int)((dip-mindip)/10.);
    if(idip>=ndip){
      idip = ndip-1;
    }
  }
    
  for(i=0;i<5;i++){
    xcor[2] -= dyCor[ilay][idip][i]*pow(aTrack->hit[iHit]->dx[1],(double)i);
  }

  return 0;
}
