#ifndef Stntuple_gui_TTrkVstView_hh
#define Stntuple_gui_TTrkVstView_hh


#include "TNamed.h"
#include "TPad.h"
#include "Stntuple/base/TStnView.hh"


class TTrkVstView: public TStnView {
protected:
  //  Int_t               fSectionToDisplay;	// a disk or vane number
  //  TVirtualPad*        fPad;
public:

  TTrkVstView();
  virtual ~TTrkVstView();

  //  TVirtualPad*  GetPad() { return fPad; }

  // int    SectionToDisplay() { return fSectionToDisplay; }
  // void   SetSectionToDisplay(int Section) { fSectionToDisplay = Section; }

  //  void   SetPad (TVirtualPad* Pad) { fPad = Pad; }
//-----------------------------------------------------------------------------
// menu
//-----------------------------------------------------------------------------
  // void   SetMinClusterEnergy(float MinEnergy);  // *MENU*
  // void   SetMinCrystalEnergy(float MinEnergy);  // *MENU*
  // void   PrintClosestCrystal();                 // *MENU*
//-----------------------------------------------------------------------------
// overloaded virtual functions of TObject
//-----------------------------------------------------------------------------
  // virtual void  Paint              (Option_t* Option = "");
  // virtual void  ExecuteEvent       (Int_t event, Int_t px, Int_t py);
  //  virtual Int_t DistancetoPrimitive(Int_t px, Int_t py);

  ClassDef(TTrkVstView,0)
};

#endif
