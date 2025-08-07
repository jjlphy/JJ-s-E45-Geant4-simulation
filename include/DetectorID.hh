// -*- C++ -*-

#ifndef DETECTOR_ID_HH
#define DETECTOR_ID_HH

#include <iostream>

#include <Rtypes.h>
#include <TString.h>

// Counters ___________________________________________________________________
const Int_t DetIdBH1     =  1;
const Int_t DetIdBH2     =  2;
const Int_t DetIdBAC     =  3;
const Int_t DetIdHTOF    =  4;
const Int_t DetIdKVC     =  5;
const Int_t DetIdFTOF    =  6;
const Int_t DetIdBVH_U     =  11; //BVH upstream 추가
const Int_t DetIdBVH_D     =  12; //BVH upstream 추가
const Int_t DetIdT0     = 13; // T0 추가

const Int_t NumOfSegBH1  = 11;
const Int_t NumOfSegBH2  = 12;
const Int_t NumOfSegBAC  =  1;
const Int_t NumOfSegHTOF = 32;
const Int_t NumOfSegKVC  = 8;
const Int_t NumOfSegFTOF = 24;
const Int_t NumOfPlaneHTOF = 8;
const Int_t NumOfSegHTOFOnePlane = NumOfSegHTOF / NumOfPlaneHTOF;
const Int_t NumOfSegBVH_U  = 10; // BVH up
const Int_t NumOfSegBVH_D  = 60; // BVH downstream
const Int_t NumOfSegT0  = 5; // T0

const Int_t DetIdPVAC    =  5;
const Int_t DetIdFAC     =  6;
const Int_t DetIdSCH     =  7;
const Int_t DetIdLAC     =  9;
const Int_t DetIdWC      = 10;
const Int_t NumOfSegFBH  = 16;
const Int_t NumOfSegCFBH = 31;
const Int_t NumOfSegPVAC =  1;
const Int_t NumOfSegFAC  =  1;
const Int_t NumOfSegSCH  = 64;
const Int_t NumOfSegLAC  = 30;
const Int_t NumOfSegWC   = 20;

// Misc _______________________________________________________________________
const Int_t DetIdTrig      =  21;
const Int_t DetIdScaler    =  22;
const Int_t DetIdEMC       =  23;
const Int_t DetIdSSDT      =  24;
const Int_t DetIdMsT       =  25;
const Int_t DetIdHulTOF    =  26;
const Int_t DetIdHulSCH    =  27;
const Int_t DetIdHulFBH    =  28;
const Int_t DetIdVmeRm     =  81;
const Int_t DetIdMsTRM     =  82;
const Int_t DetIdVP        =  83;
const Int_t NumOfSegTrig   =  32;
const Int_t NumOfSegScaler =  32;
const Int_t NumOfSegEMC    =   1;
const Int_t NumOfSegSSDT   =  16;
const Int_t SpillEndFlag   =  20;// 0-based
const Int_t NumOfPlaneVmeRm = 11;
const Int_t NumOfSegVP     =   5;

enum eTriggerFlag
  {
    kUB       =  0,
    kUBUB     =  1,
    kKPI      =  2,
    kKK       =  3,
    kUBPS     =  4,
    kUBUBPS   =  5,
    kKPIPS    =  6,
    kKKPS     =  7,
    kMtx1     =  8,
    kMtx2     =  9,
    kMtx3     = 10,
    kMtx4     = 11,
    kKIn      = 12,
    kPiIn     = 13,
    kKOut     = 14,
    kPiOut    = 15,
    kKBeam    = 16,
    kPiBeam   = 17,
    kPBeam    = 18,
    kClock    = 19,
    kSpillEnd = 20
  };

static TString sTriggerFlag[NumOfSegTrig] =
  {
    "UB",
    "UBUB",
    "KPI",
    "KK",
    "UBPS",
    "UBUBPS",
    "KPIPS",
    "KKPS",
    "Mtx1",
    "Mtx2",
    "Mtx3",
    "Mtx4",
    "KIn",
    "PiIn",
    "KOut",
    "PiOut",
    "KBeam",
    "PiBeam",
    "PBeam",
    "Clock",
    "SpillEnd"
  };

const Int_t DetIdVmeCalib      = 999;
const Int_t NumOfPlaneVmeCalib =   5;
const Int_t NumOfSegVmeCalib   =  32;

// Trackers ___________________________________________________________________
const Int_t DetIdBC1  = 101;
const Int_t DetIdBC2  = 102;
const Int_t DetIdBC3  = 103;
const Int_t DetIdBC4  = 104;
const Int_t DetIdSDC1 = 105;
const Int_t DetIdSDC3 = 106;
const Int_t DetIdSDC4 = 107;
const Int_t DetIdBFT  = 110;
const Int_t DetIdSFT  = 111;
const Int_t DetIdSSD1 = 121;
const Int_t DetIdSSD2 = 122;

const Int_t PlMinBcIn    =   1;
const Int_t PlMaxBcIn    =  12;
const Int_t PlMinBcOut   =  13;
const Int_t PlMaxBcOut   =  24;
const Int_t PlMinSdcIn   =   1;
const Int_t PlMaxSdcIn   =  14;
const Int_t PlMinSdcOut  =  31;
const Int_t PlMaxSdcOut  =  38;
const Int_t PlMinSsdIn   =   1;
const Int_t PlMaxSsdIn   =   4;
const Int_t PlMinSsdOut  =   5;
const Int_t PlMaxSsdOut  =   8;
const Int_t PlMinSsdOutSdcIn = 1;
const Int_t PlOffsBc     = 100;
const Int_t PlOffsSdcIn  =   0;
const Int_t PlOffsSdcOut =  30;
const Int_t PlOffsSsd    =   6;

const Int_t NumOfLayersBc   = 6;
const Int_t NumOfLayersSSD1 = 4;
const Int_t NumOfLayersSSD2 = 4;
const Int_t NumOfLayersSDC1 = 6;
const Int_t NumOfLayersSFT  = 3;
const Int_t NumOfLayersSDC3 = 4;
const Int_t NumOfLayersSDC4 = 4;
const Int_t NumOfLayersBcIn   = PlMaxBcIn   - PlMinBcIn   + 1;
const Int_t NumOfLayersBcOut  = PlMaxBcOut  - PlMinBcOut  + 1;
const Int_t NumOfLayersSsdIn  = PlMaxSsdIn  - PlMinSsdIn  + 1;
const Int_t NumOfLayersSsdOut = PlMaxSsdOut - PlMinSsdOut + 1;
const Int_t NumOfLayersSdcIn  = PlMaxSdcIn  - PlMinSdcIn  + 1;
const Int_t NumOfLayersSdcOut = PlMaxSdcOut - PlMinSdcOut + 1;

const Int_t MaxWireBC1      = 256;
const Int_t MaxWireBC2      = 256;
const Int_t MaxWireBC3      =  64;
const Int_t MaxWireBC4      =  64;

const Int_t MaxWireSDC1     =  64;
const Int_t MaxWireSDC3     = 128;
const Int_t MaxWireSDC4X    =  96;
const Int_t MaxWireSDC4Y    =  64;

const Int_t NumOfSampleSSD  =    8;
const Double_t SamplingIntervalSSD = 25.; //[ns]
const Int_t NumOfSegSSD0    = 1536;
const Int_t NumOfSegSSD1    = 1536;
const Int_t NumOfSegSSD2    = 1536;

const Int_t NumOfPlaneBFT   =   2;
const Int_t NumOfSegBFT     = 160;
// SFT X layer has U D plane.
// SFT UV layers have only U plnane.
enum SFT_PLANE{ SFT_X1, SFT_X2, SFT_V, SFT_U };
const Int_t NumOfPlaneSFT   =   4;
const Int_t NumOfSegSFT_X   = 256;
const Int_t NumOfSegSFT_UV  = 320;

// HypTPC _____________________________________________________________________
const Int_t DetIdTPC  = 201;
const Int_t NumOfPadTPCIn  = 10;
const Int_t NumOfPadTPCOut = 22;
const Int_t NumOfPadTPC = NumOfPadTPCIn + NumOfPadTPCOut;

#endif
