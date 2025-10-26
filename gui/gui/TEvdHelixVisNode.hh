#ifndef Stntuple_gui_TEvdHelixVisNode_hh
#define Stntuple_gui_TEvdHelixVisNode_hh

#include "Gtypes.h"
#include "TClonesArray.h"
#include "TArc.h"
					// needed for rootcint... ComboHit
#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/Ptr.h"

#include "Offline/RecoDataProducts/inc/HelixSeed.hh"

#include "Stntuple/gui/TStnVisNode.hh"

class TStnTrackBlock;
class TSimpBlock;
class TStnHelixBlock;

namespace stntuple {
  class TEvdHelixSeed;
  class TEvdSimParticle;
  class TEvdStrawTracker;
}

class TEvdHelixVisNode: public TStnVisNode {
public:
  
protected:

  const mu2e::HelixSeedCollection*  fHsColl;

  art::InputTag                     fHsCollTag;
  art::InputTag                     fShCollTag;   // straw/combo hit collection
  art::InputTag                     fSdmcCollTag; // straw digi collection

  TStnHelixBlock*                   fHelixBlock;
  Color_t                           fHelixColor;

  TArc*                             fArc;
  const mu2e::TimeCluster*          fTimeCluster;

  Int_t                             fTimeWindow;
  TObjArray*                        fListOfHelixSeeds;

  stntuple::TEvdStrawTracker*       fTracker;
public:
//-----------------------------------------------------------------------------
// constructors and destructor
//-----------------------------------------------------------------------------
  TEvdHelixVisNode();
  TEvdHelixVisNode(const char* Name, TStnHelixBlock* fHelixBlock);
  virtual ~TEvdHelixVisNode();

//-----------------------------------------------------------------------------
// accessors
//-----------------------------------------------------------------------------
  TObjArray* GetListOfHelices() { return fListOfHelixSeeds; }
  Color_t    GetHelixColor   () { return fHelixColor;   }

  int        GetNHelices()      { return fListOfHelixSeeds->GetEntriesFast(); }

  stntuple::TEvdHelixSeed*  GetEvdHelixSeed   (int I) { 
    return (stntuple::TEvdHelixSeed*)    fListOfHelixSeeds->At(I); 
  }

  art::InputTag& HsCollTag  () { return fHsCollTag   ; }
  art::InputTag& ShCollTag  () { return fShCollTag   ; }
  art::InputTag& SdmcCollTag() { return fSdmcCollTag ; }
//-----------------------------------------------------------------------------
// modifiers
//-----------------------------------------------------------------------------
  void SetHelixSeedCollTag (art::InputTag& Tag) { fHsCollTag   = Tag; }
  void SetSdmcCollTag      (art::InputTag& Tag) { fSdmcCollTag = Tag; }
  void SetShCollTag        (art::InputTag& Tag) { fShCollTag   = Tag; }

//  void  SetListOfTracks(TObjArray*  List ) { fListOfTracks = List ; }
  void  SetHelixColor  (Color_t     Color) { fHelixColor   = Color; }

//-----------------------------------------------------------------------------
// overloaded methods of TVisNode
//-----------------------------------------------------------------------------
  virtual int  InitEvent();

  virtual void NodePrint(const void* Object, const char* ClassName) ;

  virtual void  PaintCal (Option_t* option = "");
  virtual void  PaintCrv (Option_t* option = "");
  virtual void  PaintPhiZ(Option_t* option = "");
  virtual void  PaintRZ  (Option_t* option = "");
  virtual void  PaintTZ  (Option_t* option = "");
  virtual void  PaintVST (Option_t* option = "");
  virtual void  PaintVRZ (Option_t* option = "");
  virtual void  PaintXY  (Option_t* option = "");
//-----------------------------------------------------------------------------
// overloaded methods of TObject
//-----------------------------------------------------------------------------
  virtual int  DistancetoPrimitiveXY(Int_t px, Int_t py);
  virtual int  DistancetoPrimitiveRZ(Int_t px, Int_t py);
  virtual int  DistancetoPrimitiveTZ(Int_t px, Int_t py);

  virtual void Clear(const char* Opt = "")       ; // **MENU**
  virtual void Print(const char* Opt = "") const ; // **MENU**

  ClassDef(TEvdHelixVisNode,0)
};


#endif
