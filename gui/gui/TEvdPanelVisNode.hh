//-----------------------------------------------------------------------------
// P.Murat: node for a side view of a panel
//-----------------------------------------------------------------------------
#ifndef Stntuple_gui_TEvdPanelVisNode_hh
#define Stntuple_gui_TEvdPanelVisNode_hh

#include "Gtypes.h"
#include "TClonesArray.h"
#include "TArc.h"
#include "TH1F.h"
//------------------------------------------------------------------------------
#include "canvas/Utilities/InputTag.h"
#include "art/Framework/Principal/Event.h"

#include "Stntuple/gui/TStnVisNode.hh"
#include "Offline/RecoDataProducts/inc/ComboHit.hh"

class TStrawHitBlock;

namespace stntuple {

  class TEvdPanel;

class TEvdPanelVisNode: public TStnVisNode {
public:
  
  stntuple::TEvdPanel*      fPanel;	      // not owned

  TObjArray*                fListOfStrawHits; // to be displayed

  int                       fWfIndex;	      // index of the waveform to be displayed
  TH1F*                     fWfHist;	      // display waveformas a histogram

public:
//-----------------------------------------------------------------------------
// constructors and destructor
//-----------------------------------------------------------------------------
  TEvdPanelVisNode(const char* Name = "", TEvdPanel* Panel = nullptr);
  virtual ~TEvdPanelVisNode();
//-----------------------------------------------------------------------------
// accessors
//-----------------------------------------------------------------------------
//  TObjArray* GetListOfTracks() { return fListOfTracks; }
  // Color_t    GetTrackColor  () { return fTrackColor;   }

  // int	      GetNStrawHits  () { return fListOfStrawHits->GetEntriesFast(); }
  //  TVstStrawHit* VstStrawHit(int I) { return (TVstStrawHit*) fListOfStrawHits->UncheckedAt(I); }
//-----------------------------------------------------------------------------
// modifiers
//-----------------------------------------------------------------------------
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
  virtual int  DistancetoPrimitive   (Int_t px, Int_t py);
  virtual int  DistancetoPrimitiveXY (Int_t px, Int_t py);
  virtual int  DistancetoPrimitiveRZ (Int_t px, Int_t py);
  virtual int  DistancetoPrimitiveVST(Int_t px, Int_t py);
  virtual int  DistancetoPrimitiveVRZ(Int_t px, Int_t py);
  //  virtual int  DistancetoPrimitiveWF (Int_t px, Int_t py);

  // virtual void Clear(const char* Opt = "")       ; // **MENU**
  virtual void Print(const char* Opt = "") const ; // **MENU**

  ClassDef(TEvdPanelVisNode,0)
};
}

#endif
