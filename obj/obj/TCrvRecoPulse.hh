#ifndef STNTUPLE_TCrvRecoPulse
#define STNTUPLE_TCrvRecoPulse

#include "TClonesArray.h"

namespace mu2e {
  class CrvRecoPulse;
};

class TCrvRecoPulse: public TObject {
public:
  
  int       fSbid;
  int       fSipm;
  int       fRoc;
  int       fFeb;
  int       fFebCh;               // channel within the feb

  float     fPes;
  float     fPesPh;
  float     fTime;
  float     fPh;
  float     fBeta;
  float     fChi2;
  float     fLeTime;
  float     fPed;

  mu2e::CrvRecoPulse* fOfflineCrvp; // !
//-----------------------------------------------------------------------------
//  functions
//-----------------------------------------------------------------------------
public:
					// ****** constructors and destructor
  TCrvRecoPulse();
  virtual ~TCrvRecoPulse();
//-----------------------------------------------------------------------------
// accessors
//-----------------------------------------------------------------------------
  int   Sbid     () { return fSbid; }
  int   Sipm     () { return fSipm; }
  int   Roc      () { return fRoc;  }
  int   Feb      () { return fFeb;  }
  int   FebCh    () { return fFebCh;}
  
  int   OfflineChID() { return fSbid*4+fSipm; }

  float Time  ()   { return fTime; }
  float Pes   ()   { return fPes; }
  float PesPh ()   { return fPesPh; }
  float Ph    ()   { return fPh; }
  float Beta  ()   { return fBeta; }
  float Chi2  ()   { return fChi2; }
  float LeTime()   { return fLeTime; }

//-----------------------------------------------------------------------------
// modifiers
//-----------------------------------------------------------------------------
  void Set(int Sbid, int Sipm, int Roc, int Feb, int FebCh,
           float Pes, float PesPh, float Time, float Ph, float Beta,
           float Chi2, float LeTime, float Pedestal);
//-----------------------------------------------------------------------------
// overloaded methods of TObject
//-----------------------------------------------------------------------------
  virtual void Clear(Option_t* opt="") override;
  virtual void Print(Option_t* opt="") const override;

  ClassDefOverride(TCrvRecoPulse,1)	         // CRV reco pulse
};


#endif
