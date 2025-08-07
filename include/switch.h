#ifndef __SWITCH_H__
#define __SWITCH_H__

enum {DATA, MC};

typedef struct {
  float BField;
  int sol_current; /*NTPC*/
  int uniformBForRK;
  int rkTracking;
  int ExBCorrection;
  int makeDriftTable;
  int uniformEForExB;
  int uniformBForExB;
  int anodeYcorrection; /*NTPC anode*/
  int normPad;   /*NTPC normalization of pad gain*/
  int normAnode; /*NTPC anode*/
  int use5wires; /*NTPC anode*/
  int dipCor;    /*dip angle correction for hit position*/
  int wireselectdip;/*dip angle-dependent selection of wire hit*/
  int pulsepeakfit;
  int anaTpc;
  int anaTag;
  int anaNDC;
  int anaTpcTrig;
  int ntTpc;
  int ntHitRaw;
  int ntWireRaw;
  int ntTag;
  int ntNDC;
  int ntTpcTrigRaw;
  int ntTpcTrig;
  int wireOnTimePulse;
  int resolution;
  int ntDebug;
} Switch;

int readSwitch(char* fname, Switch* sw, char* dataFile);

#endif
