#ifndef TEvdFrame_hh
#define TEvdFrame_hh

#include "TGFrame.h"
#include "TCanvasImp.h"
#include "TRootEmbeddedCanvas.h"

class TControlBar;
class TGMenuBar;
class TGPopupMenu;
class TGLayoutHints;
class TGStatusBar;
class TGGroupFrame;
class TGCheckButton;
class TGTextButton;
class TGRadioButton;

class TVisManager;

//_____________________________________________________________________________
class TEvdFrame: public TGMainFrame, public TCanvasImp {
public:

  enum { fGroupFrameWidth = 150 };

protected:

  TGMenuBar           *fMenuBar;	   // !
  TGPopupMenu         *fMenuFile;	   // !
  TGPopupMenu         *fMenuEdit;	   // !
  TGPopupMenu         *fMenuOption;	   // !
  TGPopupMenu         *fMenuOpen;	   // !
  TGPopupMenu         *fMenuPrint;	   // !
  TGPopupMenu         *fMenuHelp;	   // !

  TGLayoutHints       *fMenuBarLayout;	   // !
  TGLayoutHints       *fMenuBarItemLayout; // !
  TGLayoutHints       *fMenuBarHelpLayout; // !

  TGHorizontalFrame   *fHorizontalFrame;   // !
  TGGroupFrame        *fGroupFrame;        // !
  TGTextButton        *fNextEvent;         // ! "Next Event" command

  TControlBar         *fEditorBar;         // ! Editor control bar

  TRootEmbeddedCanvas *fEmbeddedCanvas;	   // ! canvas widget
  TGLayoutHints       *fCanvasLayout;	   // ! layout for canvas widget

  TGStatusBar         *fStatusBar;         // !
  TGLayoutHints       *fStatusBarLayout;   // ! layout for the status bar

  //  TGRadioButton       *fRb[2];             // ! [0]: SH [1]

  TGCheckButton*      fDisplayHelices;     // ! 

  TVisManager*        fVisManager;	   // !
  int                 fView;		   // !

public:
  TEvdFrame(const char* name, const char* title, 
	    TVisManager* VisManager, Int_t View, 
	    UInt_t w, UInt_t h,
	    UInt_t options = kMainFrame | kVerticalFrame);

  virtual ~TEvdFrame();

					// ****** accessors

  TVisManager*    GetVisManager() { return fVisManager; }
  int             GetView      () { return fView; }

					// ****** setters

					// ****** slots
  void DoOK();
  void DoCancel();
  void DoTab(Int_t id);

  void HandleButtons(Int_t id = -1);
  void EditorBar();
  virtual Bool_t ProcessMessage(Long_t msg, Long_t parm1, Long_t parm2);

  void  ShowStatusBar(Bool_t show);
  void  SetStatusText(const char* txt = 0, Int_t partidx = 0);

  TCanvas*  GetCanvas() { return fEmbeddedCanvas->GetCanvas(); }

  virtual void CloseWindow();

  void OpenLegendPopup();
  void OpenLegendCalView();
  void OpenLegendCrvView();

  void PrintToPNG(TCanvas* running_canvas);

  ClassDef(TEvdFrame,1)
};


#endif
