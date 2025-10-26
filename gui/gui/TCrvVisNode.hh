///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
#ifndef Stntuple_gui_TCrvVisNode_hh
#define Stntuple_gui_TCrvVisNode_hh

#include "Gtypes.h"
#include "TClonesArray.h"
#include "TH1.h"	
#include "TPad.h"
#include "TBox.h"

#include "Stntuple/gui/TStnVisNode.hh"
#include "Stntuple/gui/TEvdCrvBar.hh"

#ifndef __CINT__

#include "Offline/CosmicRayShieldGeom/inc/CosmicRayShield.hh"
#include "Offline/CosmicRayShieldGeom/inc/CRSScintillatorBar.hh"
#include "Offline/CosmicRayShieldGeom/inc/CRSScintillatorLayer.hh"
#include "Offline/CosmicRayShieldGeom/inc/CRSScintillatorModule.hh"
#include "Offline/CosmicRayShieldGeom/inc/CRSScintillatorShield.hh"
#include "Offline/DataProducts/inc/CRSScintillatorBarIndex.hh"
#include "Offline/RecoDataProducts/inc/CrvRecoPulse.hh"

#else

namespace mu2e {
  class CosmicRayShield;
  class CRSScintillatorBar;
  class CRSScintillatorLayer;
  class CRSScintillatorModule;
  class CRSScintillatorShield;
  class CrvRecoPulseCollection;
  class CRSScintillatorBarIndex;
  class CrvRecoPulse;
  struct CrvRecoPulse::CrvSingleRecoPulse;
};

#endif

class TCrvVisNode : public TStnVisNode {
protected:
  mu2e::CrvRecoPulseCollection** fCrvRecoPulsesCollection;

  int		fSectionID;
  TObjArray*	fListOfBars;

  float		fMinPulsePEs;
  float		ftimeLow;
  float		ftimeHigh;
  
  Int_t         colorPalette[1000];
  
public:
//-----------------------------------------------------------------------------
// Constructors and Destructor
 //-----------------------------------------------------------------------------
  TCrvVisNode() {}
  TCrvVisNode(const char* Name, /*const mu2e::CosmicRayShield* CRV,*/ int SectionID);

  virtual ~TCrvVisNode();
//-----------------------------------------------------------------------------
// Accessors
//-----------------------------------------------------------------------------
  int		SectionID() { return fSectionID; }
  
  int           NBars()    const { return fListOfBars->GetEntriesFast(); }
  TEvdCrvBar*   EvdBar(int I) { return (TEvdCrvBar*) fListOfBars->UncheckedAt(I); }
  TEvdCrvBar*   EvdBarWithIndex(int BarIndex);
//-----------------------------------------------------------------------------
// Modifiers
//-----------------------------------------------------------------------------
  void	SetRecoPulsesCollection(mu2e::CrvRecoPulseCollection** List) { fCrvRecoPulsesCollection = List; }

  void	SetMinPulsePEs(float minPulsePEs){ fMinPulsePEs = minPulsePEs; }
  void	SetTimeWindow(float timeLow, float timeHigh);
//-----------------------------------------------------------------------------
// Overloaded methods of TStnVisNode
//-----------------------------------------------------------------------------
  virtual void  PaintCrv(Option_t* option = "");
  
  virtual void  PaintXY  (Option_t* option = "");
  virtual void  PaintRZ  (Option_t* option = "");
  virtual void  PaintTZ  (Option_t* option = "");
  virtual void  PaintPhiZ(Option_t* option = "");
  virtual void  PaintCal (Option_t* option = "");
  virtual void  PaintVST (Option_t* option = "");
  virtual void  PaintVRZ (Option_t* option = "");
//-----------------------------------------------------------------------------
// Overloaded methods of TVisNode
//-----------------------------------------------------------------------------
  virtual int   InitEvent();
  void	        UpdateEvent();
//-----------------------------------------------------------------------------
// Overloaded methods of TObject
//-----------------------------------------------------------------------------
  virtual void	Clear(Option_t* Opt = "");
  virtual void	Print(Option_t* Opt = "") const; // **MENU**

  virtual Int_t	DistancetoPrimitive(Int_t px, Int_t py);
  virtual Int_t	DistancetoPrimitiveXY(Int_t px, Int_t py);
  virtual Int_t	DistancetoPrimitiveRZ(Int_t px, Int_t py);

  int getCRVSection(int shieldNumber);

  ClassDef(TCrvVisNode, 0)
};


#endif
