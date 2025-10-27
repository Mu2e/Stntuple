///////////////////////////////////////////////////////////////////////////////
// vis node displays one wedge
///////////////////////////////////////////////////////////////////////////////
#ifndef TEvdTracker_hh
#define TEvdTracker_hh

#include "Gtypes.h"
#include "TClonesArray.h"
#include "TH1.h"
#include "TPad.h"
#include "TArc.h"

namespace mu2e {
  class Tracker;
}
  
namespace stntuple {

class TEvdStation;

class TEvdTracker: public TNamed {
public:
  
protected:

  int        fNStations;
  TObjArray* fListOfStations;

  const mu2e::Tracker*  fTracker;

public:
//-----------------------------------------------------------------------------
// constructors and destructor
//-----------------------------------------------------------------------------
  TEvdTracker(const mu2e::Tracker* Tracker = NULL);
  //  TEvdTracker(const char* Name); 

  virtual ~TEvdTracker();
//-----------------------------------------------------------------------------
// accessors
//-----------------------------------------------------------------------------
  int          NStations     () { return fNStations;      }
  TObjArray*   ListOfStations() { return fListOfStations; }

  const mu2e::Tracker* GetMu2eTracker() { return fTracker; }

  TEvdStation* Station  (int I) { 
    return (TEvdStation*) fListOfStations->UncheckedAt(I); 
  }
  
  static void  ConvertPanelGeoIndices(int Station, int ZFace, int IPFace, int& Plane, int& Panel);
//-----------------------------------------------------------------------------
// modifiers
//-----------------------------------------------------------------------------
  virtual void  Paint   (Option_t* option = "");

  void  PaintXY   (Option_t* option = "");
  void  PaintRZ   (Option_t* option = "");
  void  PaintVST  (Option_t* option = "");

  //  virtual void  ExecuteEvent(Int_t event, Int_t px, Int_t py);

  virtual Int_t DistancetoPrimitive   (Int_t px, Int_t py);
  virtual Int_t DistancetoPrimitiveXY (Int_t px, Int_t py);
  virtual Int_t DistancetoPrimitiveRZ (Int_t px, Int_t py);
  virtual Int_t DistancetoPrimitiveVST(Int_t px, Int_t py);

  //  virtual void   Print(const char* Opt = "") const ; // **MENU**

  ClassDef(stntuple::TEvdTracker,0)
};
}

#endif
