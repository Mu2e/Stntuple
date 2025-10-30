///////////////////////////////////////////////////////////////////////////////
// tracekr envelope
///////////////////////////////////////////////////////////////////////////////
#ifndef Stntuple__geom_TStnTracker_hh
#define Stntuple__geom_TStnTracker_hh

#include "TString.h"
#include "TFolder.h"
#include "TFile.h"

#include "TGeoTube.h"
#include "TGeoVolume.h"

class TStnTracker : public TObject {
public:
  struct GeomData_t {
    float fRMin;
    float fRMax;
    float fX0;
    float fY0;
    float fZ0;
    float fHalfLength;
  };

  TGeoTube*   fEnvelope;
  TGeoVolume* fMotherVolume;
  
public:
  int      fInitialized;
  
//-----------------------------------------------------------------------------
// methods
//-----------------------------------------------------------------------------
  TStnTracker();
  TStnTracker(GeomData_t* Geom);

  ~TStnTracker();
//-----------------------------------------------------------------------------
// accessors
//-----------------------------------------------------------------------------
  int     Initialized() { return fInitialized; }
//-----------------------------------------------------------------------------
// other methods
//-----------------------------------------------------------------------------
  int Init(GeomData_t* Geom);
  int InitEvent();
//-----------------------------------------------------------------------------
// overloaded methods of TObject
//-----------------------------------------------------------------------------
  void    Clear(Option_t* Opt = "") override ;
  void    Print(Option_t* Opt = "") const override;

  ClassDef(TStnTracker,0)

};

//}  // end namespace

#endif
