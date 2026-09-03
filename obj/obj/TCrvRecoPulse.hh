#ifndef STNTUPLE_TCrvRecoPulse
#define STNTUPLE_TCrvRecoPulse

#include "TClonesArray.h"


class TCrvRecoPulse: public TObject {
public:
  int       fIndex;                   // index in the list of _all_ pulses
  int       fNPe;
  int       fNPeHeight;
  int       fNDigis;
  int       fBar;
  int       fSipm;

  float     fTime;
  float     fHeight;
  float     fWidth;
  float     fChi2;
  float     fLeTime;
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
  int   Index    () { return fIndex; }
  int   NPe      () { return fNPe   ; }
  int   NPeHeight() { return fNPeHeight ; }
  int   NDigis   () { return fNDigis; }
  int   Bar      () { return fBar; }
  int   Sipm     () { return fSipm; }

  float Time  ()   { return fTime; }
  float Height()   { return fHeight; }
  float Width ()   { return fWidth; }
  float Chi2  ()   { return fChi2; }
  float LeTime()   { return fLeTime; }

//-----------------------------------------------------------------------------
// modifiers
//-----------------------------------------------------------------------------
  void Set(int I, int NPes, int NPesHeight, int NInd, int Bar, int Sipm, 
	   float Time, float Height, float Width, float Chi2, float LeTime);
//-----------------------------------------------------------------------------
// overloaded methods of TObject
//-----------------------------------------------------------------------------
  void Clear(Option_t* opt="");
  void Print(Option_t* opt="") const;

  ClassDefOverride(TCrvRecoPulse,1)	         // CRV reco pulse
};


#endif
