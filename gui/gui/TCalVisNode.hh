// vis node displays one wedge

#ifndef TCalVisNode_hh
#define TCalVisNode_hh

#include "Gtypes.h"
#include "TClonesArray.h"
#include "TH1.h"
#include "TPad.h"
#include "TArc.h"

#include "Stntuple/gui/TStnVisNode.hh"
#include "Stntuple/gui/TEvdCluster.hh"

#ifndef __CINT__

#include "Offline/RecoDataProducts/inc/CaloCluster.hh"
#include "Offline/RecoDataProducts/inc/CaloHit.hh"
#include "Offline/CalorimeterGeom/inc/Disk.hh"
#include "Offline/CalorimeterGeom/inc/DiskCalorimeter.hh"

#else

namespace mu2e {
  class CaloCrystalHitCollection;
  class CaloClusterCollection;
  class CaloSection;
  class DiskCalorimeter;
  class Disk;
};

#endif

class TEvdCalSection;
class TEvdCrystal;


class TCalVisNode: public TStnVisNode {
public:
  enum {
    kPickClusters = 0,
    kPickTracks   = 1,
    kPickStrips   = 2
  };
protected:
  const mu2e::CaloClusterCollection**     fListOfClusters;
  const mu2e::CaloHitCollection**         fListOfCrystalHits;

  TObjArray**        fListOfTracks;
  int                fSectionID;

  TEvdCalSection*    fEvdCalSection; 	// disk ?

  TObjArray*         fListOfEvdClusters;
  int                fNClusters;
  TObjArray*         fListOfEvdCrystals;

  Int_t              fDisplayHits;
  Int_t              fPickMode;
  double             fMinClusterEnergy;
  double             fMinCrystalEnergy;

  int                fNCrystals;
  int                fFirst;

public:
					// ****** constructors and destructor

  TCalVisNode(): TStnVisNode() {}
  TCalVisNode(const char* Name, const mu2e::Disk* Disk, int SectionID); 

  virtual ~TCalVisNode();
//-----------------------------------------------------------------------------
// accessors
//-----------------------------------------------------------------------------
  const mu2e::CaloClusterCollection* GetListOfClusters() { return *fListOfClusters; }

  TObjArray* GetListOfEvdClusters() { return fListOfEvdClusters; }

  int NClusters() { return fListOfEvdClusters->GetEntries(); }

  int NCrystals() { return fListOfEvdCrystals->GetEntries(); }

  TEvdCluster*  EvdCluster(int I) { return (TEvdCluster*) fListOfEvdClusters->At(I); }

  TEvdCrystal*  EvdCrystal(int I) { return (TEvdCrystal*) fListOfEvdCrystals->At(I); }

  int  SectionID() { return fSectionID; }

  int  LocalCrystalID (int CrystalID);

//-----------------------------------------------------------------------------
// modifiers
//-----------------------------------------------------------------------------
  void  SetListOfClusters (const mu2e::CaloClusterCollection** List) { 
    fListOfClusters  = List; 
  }

  void  SetListOfCrystalHits(const mu2e::CaloHitCollection** List) { 
    fListOfCrystalHits  = List; 
  }

  void SetMinClusterEnergy(float MinEnergy) { fMinClusterEnergy = MinEnergy; }
  void SetMinCrystalEnergy(float MinEnergy) { fMinCrystalEnergy = MinEnergy; }

  TEvdCluster*   NewEvdCluster(const mu2e::CaloCluster* Cluster) {
    TEvdCluster* cl = new ((*fListOfEvdClusters)[fNClusters]) TEvdCluster(Cluster);
    fNClusters++;
    return cl;
  }

  void  SetSectionID(int SectionID) { fSectionID = SectionID; }

  void  SetListOfTracks   (TObjArray** List) { fListOfTracks    = List; }

  void  SetPickMode   (Int_t Mode) { fPickMode    = Mode; }
  void  SetDisplayHits(Int_t Mode) { fDisplayHits = Mode; }

  //  void  Set(Int_t Side, Int_t Wedge) ; // **MENU**


  virtual int   InitEvent();

  virtual void  PaintXY  (Option_t* option = "");
  virtual void  PaintRZ  (Option_t* option = "");
  virtual void  PaintTZ  (Option_t* option = "");
  virtual void  PaintPhiZ(Option_t* option = "");
  virtual void  PaintCal (Option_t* option = "");
  virtual void  PaintCrv (Option_t* option = "");
  virtual void  PaintVST (Option_t* option = "");
  virtual void  PaintVRZ (Option_t* option = "");

  //  virtual void  ExecuteEvent(Int_t event, Int_t px, Int_t py);

  //  virtual Int_t DistancetoPrimitive  (Int_t px, Int_t py);
  virtual Int_t DistancetoPrimitiveXY (Int_t px, Int_t py);
  virtual Int_t DistancetoPrimitiveRZ (Int_t px, Int_t py);
  virtual Int_t DistancetoPrimitiveCal(Int_t px, Int_t py);
//-----------------------------------------------------------------------------
// overloaded functions of TObject
//-----------------------------------------------------------------------------
  virtual void   Clear(Option_t* Opt = "");
  virtual void   Print(Option_t* Opt = "") const ; // **MENU**

  //  ClassDef(TCalVisNode,0)
};


#endif
