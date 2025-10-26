#ifndef Stntuple_gui_TEvdTimeClusterVisNode_hh
#define Stntuple_gui_TEvdTimeClusterVisNode_hh

#include "Gtypes.h"
#include "TClonesArray.h"

#include "canvas/Persistency/Common/Ptr.h"

#include "Stntuple/gui/TStnVisNode.hh"
#include "Offline/RecoDataProducts/inc/TimeCluster.hh"
#include "Offline/RecoDataProducts/inc/ComboHit.hh"
#include "Offline/RecoDataProducts/inc/StrawHitFlag.hh"

class TStnTrackBlock;
class TSimpBlock;
class TStnHelixBlock;
class TStnTimeClusterBlock;

namespace stntuple {
  class TEvdTimeCluster;
}

class TEvdTimeClusterVisNode: public TStnVisNode {
public:
  
protected:

  const mu2e::TimeClusterCollection*  fTcColl;
  const mu2e::TimeClusterCollection*  fPcColl;
  const mu2e::ComboHitCollection*     fChColl;
  const mu2e::ComboHitCollection*     fSschColl;

  TStnTimeClusterBlock*              fTcBlock;

  art::InputTag                      fTcCollTag;   // time cluster coll tag
  art::InputTag                      fPcCollTag;   // phi cluster coll tag (also time clusters)
  art::InputTag                      fChCollTag;   // combo hit collection
  art::InputTag                      fSdmcCollTag; // straw digi collection
  art::InputTag                      fSschCollTag; // single-straw combo hit collection

  Color_t                            fTcColor;
  Color_t                            fPcColor;

  TObjArray*                         fListOfTimeClusters; // list of EvdTimeCLusters (Time)
  TObjArray*                         fListOfPhiClusters;  // list of EvdTimeCLusters (Phi)

public:
//-----------------------------------------------------------------------------
// constructors and destructor
//-----------------------------------------------------------------------------
  TEvdTimeClusterVisNode();
  TEvdTimeClusterVisNode(const char* Name, TStnTimeClusterBlock* TcBlock = nullptr);

  virtual ~TEvdTimeClusterVisNode();
//-----------------------------------------------------------------------------
// accessors
//-----------------------------------------------------------------------------
  TObjArray* GetListOfTimeClusters() { return fListOfTimeClusters; }
  TObjArray* GetListOfPhiClusters () { return fListOfPhiClusters;  }

  Color_t    GetTcColor           () { return fTcColor;   }
  Color_t    GetPcColor           () { return fPcColor;   }

  int        NTimeClusters        ()      { return fListOfTimeClusters->GetEntriesFast(); }
  int        NPhiClusters         ()      { return fListOfPhiClusters->GetEntriesFast() ; }

  stntuple::TEvdTimeCluster*  EvdTimeCluster(int I) { 
    return (stntuple::TEvdTimeCluster*)    fListOfTimeClusters->At(I); 
  }

  stntuple::TEvdTimeCluster*  EvdPhiCluster(int I) { 
    return (stntuple::TEvdTimeCluster*)    fListOfPhiClusters->At(I); 
  }

  art::InputTag& TcCollTag()   { return fTcCollTag   ; }
  art::InputTag& PcCollTag()   { return fPcCollTag   ; }
  art::InputTag& ChCollTag()   { return fChCollTag   ; }
  art::InputTag& SschCollTag() { return fSschCollTag ; }
  art::InputTag& SdmcCollTag() { return fSdmcCollTag ; }

  const mu2e::TimeClusterCollection*  TcColl () { return fTcColl; }
  const mu2e::TimeClusterCollection*  PcColl () { return fPcColl; }
  const mu2e::ComboHitCollection*     ChColl () { return fChColl; }
  //  const mu2e::StrawHitFlagCollection* ChfColl() { return fChfColl; }
//-----------------------------------------------------------------------------
// modifiers
//-----------------------------------------------------------------------------
  void SetChCollTag  (art::InputTag& CollTag) { fChCollTag   = CollTag; }
  void SetSschCollTag(art::InputTag& CollTag) { fSschCollTag = CollTag; }
  void SetSdmcCollTag(art::InputTag& CollTag) { fSdmcCollTag = CollTag; }
  void SetTcCollTag  (art::InputTag& CollTag) { fTcCollTag   = CollTag; }
  void SetPcCollTag  (art::InputTag& CollTag) { fPcCollTag   = CollTag; }

  void SetTcColor    (Color_t       Color   ) { fTcColor    = Color;   }
  void SetPcColor    (Color_t       Color   ) { fPcColor    = Color;   }
//-----------------------------------------------------------------------------
// overloaded methods of TVisNode
//-----------------------------------------------------------------------------
  virtual int  InitEvent();

  virtual void  PaintXY  (Option_t* option = "");
  virtual void  PaintRZ  (Option_t* option = "");
  virtual void  PaintTZ  (Option_t* option = "");
  virtual void  PaintPhiZ(Option_t* option = "");
  virtual void  PaintCrv (Option_t* option = "");
  virtual void  PaintCal (Option_t* option = "");
  virtual void  PaintVST (Option_t* option = "");
  virtual void  PaintVRZ (Option_t* option = "");
//-----------------------------------------------------------------------------
// overloaded methods of TObject
//-----------------------------------------------------------------------------
  virtual int  DistancetoPrimitiveTZ(Int_t px, Int_t py);

  virtual void Clear(const char* Opt = "")       ;

  virtual void NodePrint(const void* Object, const char* ClassName);

  virtual void Print(const char* Opt = "") const ; // **MENU**

  ClassDef(TEvdTimeClusterVisNode,0)
};


#endif
