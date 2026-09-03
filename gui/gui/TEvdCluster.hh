///////////////////////////////////////////////////////////////////////////////
// vis node displays one wedge
///////////////////////////////////////////////////////////////////////////////
#ifndef TEvdCluster_hh
#define TEvdCluster_hh

#include "Gtypes.h"
#include "TClonesArray.h"
#include "TH1.h"
#include "TPad.h"
#include "TArc.h"

#include "Stntuple/base/TVisNode.hh"
#include "Stntuple/gui/TEvdCrystal.hh"

#ifndef __CINT__
#include "Offline/RecoDataProducts/inc/CaloCluster.hh"
#include "Offline/CalorimeterGeom/inc/Disk.hh"
#else
namespace mu2e {
  class CaloCluster;
  class Disk;
};
#endif

class TEvdCluster: public TObject {
public:
  
protected:

  TEllipse*                fTrkEllipse;
  TEllipse*                fCalEllipse;

  const mu2e::CaloCluster* fCluster;

  TObjArray*               fListOfCrystals;
  const mu2e::Disk*        fDisk;
public:
//-----------------------------------------------------------------------------
// constructors and destructor
//-----------------------------------------------------------------------------
  TEvdCluster() {}
  TEvdCluster(const mu2e::CaloCluster* fCluster); 

  virtual ~TEvdCluster();
//-----------------------------------------------------------------------------
// accessors
//-----------------------------------------------------------------------------
  const mu2e::CaloCluster*  Cluster() { return fCluster; }

  TEvdCrystal* Crystal  (int I) { return (TEvdCrystal*) fListOfCrystals->At(I); }
  int          NCrystals()      { return fListOfCrystals->GetEntriesFast(); }
//-----------------------------------------------------------------------------
// modifiers
//-----------------------------------------------------------------------------
  void    AddCrystal(TEvdCrystal* Crystal) { fListOfCrystals->Add(Crystal); }

  virtual void  Paint   (Option_t* Option = "");
  virtual void  PaintXY (Option_t* Option = "");
  virtual void  PaintCal(Option_t* Option = "");

  //  virtual void  ExecuteEvent(Int_t event, Int_t px, Int_t py);

  virtual Int_t DistancetoPrimitive  (Int_t px, Int_t py);
  virtual Int_t DistancetoPrimitiveXY(Int_t px, Int_t py);
  virtual Int_t DistancetoPrimitiveRZ(Int_t px, Int_t py);

  virtual void  Clear(const char* Opt = "");
  virtual void  Print(const char* Opt = "") const ; // **MENU**

  ClassDefOverride(TEvdCluster,0)
};


#endif
