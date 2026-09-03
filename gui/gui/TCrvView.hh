#ifndef Stntuple_gui_TCrvView_hh
#define Stntuple_gui_TCrvView_hh


#include "TNamed.h"
#include "TPad.h"

#include "Stntuple/base/TStnView.hh"

//-----------------------------------------------------------------------------
class TCrvView : public TStnView {
protected:
  Int_t		fSectionToDisplay;
  TVirtualPad*	fPad;		 

public:
  TCrvView(int Section);
  virtual ~TCrvView();
  
  TVirtualPad* GetPad() { return fPad; }

  int	       SectionToDisplay() { return fSectionToDisplay; }

  void	       SetSectionToDisplay(int Section) { fSectionToDisplay = Section; }

  void	       SetPad(TVirtualPad* Pad)			{ fPad = Pad; }
//-----------------------------------------------------------------------------
// menu
//-----------------------------------------------------------------------------
  //void	SetMinPulseHeight(float MinHeight);	// *MENU*
  void	        SetMinPulsePEs(float MinPEs);               // *MENU*
  void	        PrintClosestBar();			    // *MENU*
//-----------------------------------------------------------------------------
// overloaded virtual functions of TObject
//-----------------------------------------------------------------------------
  virtual void Paint              (Option_t* Option = "");
  virtual void ExecuteEvent       (Int_t event, Int_t px, Int_t py);
  virtual int  DistancetoPrimitive(Int_t px, Int_t py);
  
  ClassDefOverride(TCrvView, 0)
};

#endif
