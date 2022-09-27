//-----------------------------------------------------------------------------
//  Mu2e Event display - stolen from TGeant visualization
//-----------------------------------------------------------------------------
#include <time.h>

#include "TROOT.h"
#include "TApplication.h"
#include "TVirtualX.h"

#include "TPad.h"
#include "TText.h"
#include "TGMenu.h"
#include "TGMsgBox.h"
#include "TGFrame.h"
#include "TGButton.h"
#include "TGFileDialog.h"
#include "TControlBar.h"
#include "TInterpreter.h"
#include "TGStatusBar.h"
#include "TGShutter.h"

#include "TRootEmbeddedCanvas.h"
#include "TCanvas.h"

#include "Stntuple/gui/TStnFrame.hh"
#include "Stntuple/gui/TStnVisManager.hh"

ClassImp(TStnFrame)

//-----------------------------------------------------------------------------
enum TGeantCommandIdentifiers {
  M_FILE_OPEN,
  M_FILE_SAVE,
  M_FILE_SAVEAS,
  M_FILE_EXIT,

  M_TEST_DLG,
  M_TEST_MSGBOX,
  M_TEST_SLIDER,
  M_TEST_SHUTTER,
  M_TEST_PROGRESS,

  M_EDIT_EDITOR,
  M_EDIT_UNDO,
  M_EDIT_CLEARPAD,
  M_EDIT_CLEARCANVAS,

  M_OPTION_EVENT_STATUS,
  M_OPTION_AUTO_EXEC,
  M_OPTION_AUTO_RESIZE,
  M_OPTION_RESIZE_CANVAS,
  M_OPTION_MOVE_OPAQUE,
  M_OPTION_RESIZE_OPAQUE,
  M_OPTION_REFRESH,
  M_OPTION_STATISTICS,
  M_OPTION_HIST_TITLE,
  M_OPTION_FIT_PARAMS,
  M_OPTION_CAN_EDIT,

  M_HELP_CONTENTS,
  M_HELP_SEARCH,
  M_HELP_ABOUT,

  M_OPEN_XY,
  M_OPEN_RZ,
  M_OPEN_TZ,
  M_OPEN_CAL,
  M_OPEN_CRV,
  M_OPEN_VST,

  M_PRINT_STRAW_H,
  M_PRINT_COMBO_H

};

//-----------------------------------------------------------------------------
Int_t mb_button_id[9] = { kMBYes, kMBNo, kMBOk, kMBApply,
                          kMBRetry, kMBIgnore, kMBCancel,
                          kMBClose, kMBDismiss };

EMsgBoxIcon mb_icon[4] = { kMBIconStop, kMBIconQuestion,
                           kMBIconExclamation, kMBIconAsterisk };

const char *filetypes[] = { "All files",     "*",
                            "ROOT files",    "*.root",
                            "ROOT macros",   "*.C",
                            0,               0 };

//-----------------------------------------------------------------------------
TStnFrame::TStnFrame(const char* Name,
		     const char* Title, 
		     Int_t       View,
		     UInt_t      w,
		     UInt_t      h, 
		     UInt_t      Options):

  TGMainFrame(gClient->GetRoot(),w, h, Options),
  fView(View)
{

  TStnVisManager* vm = TStnVisManager::Instance();
//-----------------------------------------------------------------------------
//  create menu bar
//-----------------------------------------------------------------------------
  fMenuBarLayout     = new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX,
					 0, 0, 1, 1);
  fMenuBarItemLayout = new TGLayoutHints(kLHintsTop | kLHintsLeft, 0, 4, 0, 0);
  fMenuBarHelpLayout = new TGLayoutHints(kLHintsTop | kLHintsRight);

  fMenuFile = new TGPopupMenu(gClient->GetRoot());
  fMenuFile->AddEntry("&Open...", M_FILE_OPEN);
  fMenuFile->AddEntry("&Save", M_FILE_SAVE);
  fMenuFile->AddEntry("S&ave as...", M_FILE_SAVEAS);
  fMenuFile->AddEntry("&Close", -1);
  fMenuFile->AddSeparator();
  fMenuFile->AddEntry("&Print", -1);
  fMenuFile->AddEntry("P&rint setup...", -1);
  fMenuFile->AddSeparator();
  fMenuFile->AddEntry("E&xit", M_FILE_EXIT);

  fMenuFile->DisableEntry(M_FILE_SAVEAS);

  fMenuHelp = new TGPopupMenu(gClient->GetRoot());
  fMenuHelp->AddEntry("&Contents", M_HELP_CONTENTS);
  fMenuHelp->AddEntry("&Search...", M_HELP_SEARCH);
  fMenuHelp->AddSeparator();
  fMenuHelp->AddEntry("&About", M_HELP_ABOUT);
//-----------------------------------------------------------------------------
// EDIT menu item on top 
//-----------------------------------------------------------------------------
  fMenuEdit = new TGPopupMenu(gClient->GetRoot());
  fMenuEdit->AddEntry("&Editor",             M_EDIT_EDITOR);
  fMenuEdit->AddEntry("&Undo",               M_EDIT_UNDO);
  fMenuEdit->AddEntry("Clear &Pad",          M_EDIT_CLEARPAD);
  fMenuEdit->AddEntry("&Clear Canvas",       M_EDIT_CLEARCANVAS);
//-----------------------------------------------------------------------------
// OPTION menu item on top 
//-----------------------------------------------------------------------------
  fMenuOption = new TGPopupMenu(gClient->GetRoot());
  fMenuOption->AddEntry("&Event Status",         M_OPTION_EVENT_STATUS);
  fMenuOption->AddEntry("&Pad Auto Exec",        M_OPTION_AUTO_EXEC);
  fMenuOption->AddSeparator();
  fMenuOption->AddEntry("&Auto Resize Canvas",   M_OPTION_AUTO_RESIZE);
  fMenuOption->AddEntry("&Resize Canvas",        M_OPTION_RESIZE_CANVAS);
  fMenuOption->AddEntry("&Move Opaque",          M_OPTION_MOVE_OPAQUE);
  fMenuOption->AddEntry("Resize &Opaque",        M_OPTION_RESIZE_OPAQUE);
  fMenuOption->AddEntry("R&efresh",              M_OPTION_REFRESH);
  fMenuOption->AddSeparator();
  fMenuOption->AddEntry("Show &Statistics",      M_OPTION_STATISTICS);
  fMenuOption->AddEntry("Show &Histogram Title", M_OPTION_HIST_TITLE);
  fMenuOption->AddEntry("Show &Fit Parameters",  M_OPTION_FIT_PARAMS);
  fMenuOption->AddEntry("Can Edit Histograms",   M_OPTION_CAN_EDIT);
//-----------------------------------------------------------------------------
// OPEN menu item on top 
//-----------------------------------------------------------------------------
  fMenuOpen = new TGPopupMenu(gClient->GetRoot());
  fMenuOpen->AddEntry("&XY View",            M_OPEN_XY);
  fMenuOpen->AddEntry("&RZ View",            M_OPEN_RZ);
  fMenuOpen->AddEntry("&TZ View",            M_OPEN_TZ);
  fMenuOpen->AddEntry("&Cal View",           M_OPEN_CAL);
  fMenuOpen->AddEntry("&CRV View",           M_OPEN_CRV);
  fMenuOpen->AddEntry("&VST View",           M_OPEN_VST);
//-----------------------------------------------------------------------------
// PRINT menu item on top 
//-----------------------------------------------------------------------------
  fMenuPrint = new TGPopupMenu(gClient->GetRoot());
  fMenuPrint->AddEntry("Print &Straw Hits",  M_PRINT_STRAW_H);
  fMenuPrint->AddEntry("Print &Combo Hits",  M_PRINT_COMBO_H);
//-----------------------------------------------------------------------------
// define menu handlers
//-----------------------------------------------------------------------------
  fMenuFile  ->Associate(this);
  fMenuEdit  ->Associate(this);
  fMenuOption->Associate(this);
  fMenuHelp  ->Associate(this);
  fMenuOpen  ->Associate(this);

  fMenuBar = new TGMenuBar(this, 1, 1, kHorizontalFrame);

  fMenuBar->AddPopup("&File"  , fMenuFile  , fMenuBarItemLayout);
  fMenuBar->AddPopup("&Edit"  , fMenuEdit  , fMenuBarItemLayout);
  fMenuBar->AddPopup("&Option", fMenuOption, fMenuBarItemLayout);
  fMenuBar->AddPopup("O&pen"  , fMenuOpen  , fMenuBarItemLayout);
  fMenuBar->AddPopup("&Print" , fMenuPrint , fMenuBarItemLayout);
  fMenuBar->AddPopup("&Help"  , fMenuHelp  , fMenuBarHelpLayout);
  
  AddFrame(fMenuBar, fMenuBarLayout);
//-----------------------------------------------------------------------------
// left Group frame for commands, options and such
//-----------------------------------------------------------------------------
//  SetLayoutManager(new TGHorizontalLayout(this));

   // horizontal frame
  fHorizontalFrame = new TGHorizontalFrame(this, GetWidth()+4,GetHeight()+4,kHorizontalFrame);
  fHorizontalFrame->SetName("HorizontalFrame");

  fGroupFrame = new TGGroupFrame(fHorizontalFrame,"GroupFrame");
  fGroupFrame->SetLayoutBroken(kTRUE);
//-----------------------------------------------------------------------------
// "Next Event" button
//-----------------------------------------------------------------------------
  TGTextButton *tb(nullptr);
  tb = new TGTextButton(fGroupFrame,"Next Event",-1,
			TGTextButton::GetDefaultGC()(),
			TGTextButton::GetDefaultFontStruct(),
			kRaisedFrame);
  tb->SetTextJustify(36);
  tb->SetMargins(0,0,0,0);
  tb->SetWrapLength(-1);
  fGroupFrame->AddFrame(tb, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
  tb->MoveResize(10,20,120,25);
  tb->Connect("Pressed()", "TStnVisManager", vm, "NextEvent()");
//-----------------------------------------------------------------------------
// "Next Event" button
//-----------------------------------------------------------------------------
  tb = new TGTextButton(fGroupFrame,"Quit",-1,
			TGTextButton::GetDefaultGC()(),
			TGTextButton::GetDefaultFontStruct(),
			kRaisedFrame);
  tb->SetTextJustify(36);
  tb->SetMargins(0,0,0,0);
  tb->SetWrapLength(-1);
  tb->MoveResize(10,50,120,25);
  fGroupFrame->AddFrame(tb, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
  tb->Connect("Pressed()", "TStnVisManager", vm, "Quit()");
//-----------------------------------------------------------------------------
// shutter - below the "Next Event" button
//-----------------------------------------------------------------------------
  TGShutter *shutter = new TGShutter(fGroupFrame,kSunkenFrame);
//-----------------------------------------------------------------------------
// configure shutter items
// "Display" : choose collections to display
//-----------------------------------------------------------------------------
  TGShutterItem*    sh_item = new TGShutterItem(shutter, new TGHotString("Display"),1000,kVerticalFrame);
  TGCompositeFrame* frame   = (TGCompositeFrame *) sh_item->GetContainer();

  // TGPictureButton *fPictureButton805  = new TGPictureButton(frame,gClient->GetPicture("profile_s.xpm"),-1,TGPictureButton::GetDefaultGC()(),kRaisedFrame);
  // fPictureButton805->SetToolTipText("TProfile");
  // frame->AddFrame(fPictureButton805, new TGLayoutHints(kLHintsNormal));

  // TGPictureButton *fPictureButton802 = new TGPictureButton(frame,gClient->GetPicture("h3_s.xpm"),-1,TGPictureButton::GetDefaultGC()(),kRaisedFrame);
  // fPictureButton802->SetToolTipText("TH3");
  // frame->AddFrame(fPictureButton802, new TGLayoutHints(kLHintsNormal));

  // TGPictureButton *fPictureButton799 = new TGPictureButton(frame,gClient->GetPicture("h2_s.xpm"),-1,TGPictureButton::GetDefaultGC()(),kRaisedFrame);
  // fPictureButton799->SetToolTipText("TH2");
  // frame->AddFrame(fPictureButton799, new TGLayoutHints(kLHintsNormal));

  // TGPictureButton *fPictureButton796 = new TGPictureButton(frame,gClient->GetPicture("h1_s.xpm"),-1,TGPictureButton::GetDefaultGC()(),kRaisedFrame);
  // fPictureButton796->SetToolTipText("TH1");
  // frame->AddFrame(fPictureButton796, new TGLayoutHints(kLHintsNormal));

  // tb = new TGTextButton(frame,"test_01",-1,TGTextButton::GetDefaultGC()(),TGTextButton::GetDefaultFontStruct(),kRaisedFrame);
  // frame->AddFrame(tb, new TGLayoutHints(kLHintsNormal));
  
  fRb[0] = new TGRadioButton(frame,"display SH", 21);
  // rb->SetTextJustify(36);
  // rb->SetMargins    (0,0,0,0);
  // rb->SetWrapLength (-1);
  // rb->MoveResize    (8,64,100,40);
  fRb[0]->Connect("Clicked()", "TStnFrame", this, "DoRadio()");
  frame->AddFrame(fRb[0], new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));

  fRb[1] = new TGRadioButton(frame,"display CH", 22);
  // rb->SetTextJustify(36);
  // rb->SetMargins    (0,0,0,0);
  // rb->SetWrapLength (-1);
  // rb->MoveResize    (8,64,100,40);
  fRb[1]->Connect("Clicked()", "TStnFrame", this, "DoRadio()");
  frame->AddFrame(fRb[1], new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));

  TGCheckButton* cbtn;

  cbtn = new TGCheckButton(frame, "Helices", kDisplayHelices);
  cbtn->Connect("Clicked()", "TStnFrame", this, "DoCheckButtons()");
  frame->AddFrame(cbtn, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));

  cbtn = new TGCheckButton(frame, "Tracks", kDisplayTracks);
  cbtn->Connect("Clicked()", "TStnFrame", this, "DoCheckButtons()");
  frame->AddFrame(cbtn, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));

  cbtn = new TGCheckButton(frame, "SimParticles", kDisplaySimParticles);
  cbtn->Connect("Clicked()", "TStnFrame", this, "DoCheckButtons()");
  frame->AddFrame(cbtn, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));

  shutter->AddItem(sh_item);
//-----------------------------------------------------------------------------
// "print" : print different collections
//-----------------------------------------------------------------------------
  TGShutterItem *sh02 = new TGShutterItem(shutter, new TGHotString("print"),1001,kVerticalFrame);
  frame = (TGCompositeFrame *)sh02->GetContainer();

//-----------------------------------------------------------------------------
// "print KalSeedColls" button
//-----------------------------------------------------------------------------
  tb = new TGTextButton(frame,"p_KalSeedColls",-1,
			TGTextButton::GetDefaultGC()(),
			TGTextButton::GetDefaultFontStruct(),
			kRaisedFrame);
  tb->SetTextJustify(36);
  tb->SetMargins(0,0,0,0);
  tb->SetWrapLength(-1);
  tb->MoveResize(5,50,120,30);
  frame->AddFrame(tb, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
  tb->Connect("Pressed()", "TStnVisManager", vm, "PrintColls(=\"kalseed_colls\")");
//-----------------------------------------------------------------------------
// "print StrawDigiMCs" button
//-----------------------------------------------------------------------------
  tb = new TGTextButton(frame,"StrawDigiMCs",-1,
			TGTextButton::GetDefaultGC()(),
			TGTextButton::GetDefaultFontStruct(),
			kRaisedFrame);
  tb->SetTextJustify(36);
  tb->SetMargins(0,0,0,0);
  tb->SetWrapLength(-1);
  tb->MoveResize(5,850,120,30);
  frame->AddFrame(tb, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));

  tb->Connect("Pressed()", "TStnVisManager", vm, "PrintColls(=\"sdmc_colls\")");

  // TGPictureButton *fPictureButton821 = new TGPictureButton(fVerticalFrame820,gClient->GetPicture("f1_s.xpm"),-1,TGPictureButton::GetDefaultGC()(),kRaisedFrame);
  // fPictureButton821->SetToolTipText("print_straw_hits");
  // fVerticalFrame820->AddFrame(fPictureButton821, new TGLayoutHints(kLHintsCenterX | kLHintsTop,5,5,5,0));

  // TGPictureButton *fPictureButton824 = new TGPictureButton(fVerticalFrame820,gClient->GetPicture("f2_s.xpm"),-1,TGPictureButton::GetDefaultGC()(),kRaisedFrame);
  // fPictureButton824->SetToolTipText("TF2");
  // fVerticalFrame820->AddFrame(fPictureButton824, new TGLayoutHints(kLHintsCenterX | kLHintsTop,5,5,5,0));

  shutter->AddItem(sh02);
//-----------------------------------------------------------------------------
// "Trees" : no use yet
//-----------------------------------------------------------------------------
  TGShutterItem* shi03 = new TGShutterItem(shutter, new TGHotString("Trees"),1002,kVerticalFrame);
  frame = (TGCompositeFrame *) shi03->GetContainer();

  // TGPictureButton *fPictureButton840 = new TGPictureButton(fVerticalFrame839,gClient->GetPicture("ntuple_s.xpm"),-1,TGPictureButton::GetDefaultGC()(),kRaisedFrame);
  // fPictureButton840->SetToolTipText("TNtuple");
  // fVerticalFrame839->AddFrame(fPictureButton840, new TGLayoutHints(kLHintsCenterX | kLHintsTop,5,5,5,0));

  // TGPictureButton *fPictureButton843 = new TGPictureButton(fVerticalFrame839,gClient->GetPicture("tree_s.xpm"),-1,TGPictureButton::GetDefaultGC()(),kRaisedFrame);
  // fPictureButton843->SetToolTipText("TTree");
  // fVerticalFrame839->AddFrame(fPictureButton843, new TGLayoutHints(kLHintsCenterX | kLHintsTop,5,5,5,0));

  // TGPictureButton *fPictureButton846 = new TGPictureButton(fVerticalFrame839,gClient->GetPicture("chain_s.xpm"),-1,TGPictureButton::GetDefaultGC()(),kRaisedFrame);
  // fPictureButton846->SetToolTipText("TChain");
  // fVerticalFrame839->AddFrame(fPictureButton846, new TGLayoutHints(kLHintsCenterX | kLHintsTop,5,5,5,0));

  shutter->AddItem(shi03);

  //  sh->SetSelectedItem(fShutter02);

  shutter->Resize(100,300);
  shutter->MoveResize(10,140,120,300);

  fGroupFrame->AddFrame(shutter, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
//-----------------------------------------------------------------------------
// 'group frame' has vertical layout
//-----------------------------------------------------------------------------
  fGroupFrame->SetLayoutManager(new TGVerticalLayout(fGroupFrame));
  fGroupFrame->Resize(150,770);
  fGroupFrame->Resize(fGroupFrameWidth,GetHeight()+4);
//-----------------------------------------------------------------------------
// Create canvas and canvas container that will host the ROOT graphics
// fHorizontalFrame holds them
//-----------------------------------------------------------------------------
  fEmbeddedCanvas = new TRootEmbeddedCanvas(Name, 
					    fHorizontalFrame, 
					    GetWidth ()+4, 
					    GetHeight()+4,
					    kSunkenFrame | kDoubleBorder);

  fCanvasLayout = new TGLayoutHints(kLHintsExpandX | kLHintsExpandY);
  TCanvas* c = fEmbeddedCanvas->GetCanvas();

					// pad for the event display
  c->cd();
  TString name1(Name);
  name1 += "_1";
  TPad *p1 = new TPad(name1, "p1",0.0,0.0,1,0.96);
  p1->Draw();
  c->cd();

  TString name2(Name);
  name2 += "_2";
  TPad *p2 = new TPad(name2, "p2",0.7,0.96,1,1);
  p2->Draw();
  p2->cd();

  time_t t2 = time(0);
  tm* t22 = localtime(&t2);
  TText *text = new TText(0.05,0.3,asctime(t22));
  text->SetTextFont(22);
  text->SetTextSize(0.4);
  text->Draw();
  p2->Modified();

  c->cd();

  TString name_title(Name);
  name2 += "_title";
  TPad *title = new TPad(name_title, "title",0.,0.96,0.7,1);
  title->Draw();

  title->Modified();
 
  fHorizontalFrame->AddFrame(fGroupFrame, new TGLayoutHints(kLHintsLeft | kLHintsTop |  kLHintsExpandY,2,2,2,2));
  fHorizontalFrame->AddFrame(fEmbeddedCanvas, fCanvasLayout);

  AddFrame(fHorizontalFrame, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY,1,1,1,1));
//-----------------------------------------------------------------------------
// no editor bar by default
//-----------------------------------------------------------------------------
  fEditorBar     = NULL;
//-----------------------------------------------------------------------------
// status bar (hidden by default)
//-----------------------------------------------------------------------------
   int parts[] = { 33, 10, 10, 47 };
   fStatusBar  = new TGStatusBar(this, 10, 10);
   fStatusBar->SetParts(parts, 4);
   fStatusBarLayout = new TGLayoutHints(kLHintsBottom | kLHintsLeft | 
					kLHintsExpandX, 2, 2, 1, 1);
   AddFrame(fStatusBar, fStatusBarLayout);
   HideFrame(fStatusBar);
//-----------------------------------------------------------------------------
// final actions
//-----------------------------------------------------------------------------
  MapSubwindows();
  Resize(GetDefaultSize());

  SetWindowName(Title);

  MapWindow();
}

//-----------------------------------------------------------------------------
TStnFrame::~TStnFrame() {
  // Delete window, as TStnFrame's do not exist by themselves, but they are
  // always managed by TVisManager, we need to erase this frame from the
  // list of frames

  TVisManager::Instance()->GetListOfCanvases()->Remove(GetCanvas());

  delete fStatusBar;
  delete fStatusBarLayout;

  delete fEmbeddedCanvas;
  delete fCanvasLayout;

  delete fMenuFile;
  delete fMenuEdit;
  delete fMenuOption;
  delete fMenuOpen;
  delete fMenuPrint;
  delete fMenuHelp;

  delete fMenuBarItemLayout;
  delete fMenuBarHelpLayout;
  delete fMenuBarLayout;

  delete fMenuBar;

  if (fEditorBar) delete fEditorBar;

}

//-----------------------------------------------------------------------------
void TStnFrame::EditorBar() {
  // Create the Editor Controlbar

   TControlBar *ed = new TControlBar("vertical", "Editor");
   ed->AddButton("Arc",       "gROOT->SetEditorMode(\"Arc\")",       "Create an arc of circle");
   ed->AddButton("Line",      "gROOT->SetEditorMode(\"Line\")",      "Create a line segment");
   ed->AddButton("Arrow",     "gROOT->SetEditorMode(\"Arrow\")",     "Create an Arrow");
   ed->AddButton("Button",    "gROOT->SetEditorMode(\"Button\")",    "Create a user interface Button");
   ed->AddButton("Diamond",   "gROOT->SetEditorMode(\"Diamond\")",   "Create a diamond");
   ed->AddButton("Ellipse",   "gROOT->SetEditorMode(\"Ellipse\")",   "Create an Ellipse");
   ed->AddButton("Pad",       "gROOT->SetEditorMode(\"Pad\")",       "Create a pad");
   ed->AddButton("Pave",      "gROOT->SetEditorMode(\"Pave\")",      "Create a Pave");
   ed->AddButton("PaveLabel", "gROOT->SetEditorMode(\"PaveLabel\")", "Create a PaveLabel (prompt for label)");
   ed->AddButton("PaveText",  "gROOT->SetEditorMode(\"PaveText\")",  "Create a PaveText");
   ed->AddButton("PavesText", "gROOT->SetEditorMode(\"PavesText\")", "Create a PavesText");
   ed->AddButton("PolyLine",  "gROOT->SetEditorMode(\"PolyLine\")",  "Create a PolyLine (TGraph)");
   ed->AddButton("CurlyLine", "gROOT->SetEditorMode(\"CurlyLine\")", "Create a Curly/WavyLine");
   ed->AddButton("CurlyArc",  "gROOT->SetEditorMode(\"CurlyArc\")",  "Create a Curly/WavyArc");
   ed->AddButton("Text/Latex","gROOT->SetEditorMode(\"Text\")",      "Create a Text/Latex string");
   ed->AddButton("Marker",    "gROOT->SetEditorMode(\"Marker\")",    "Create a marker");
   ed->Show();
   fEditorBar = ed;
}

//-----------------------------------------------------------------------------
Bool_t TStnFrame::ProcessMessage(Long_t msg, Long_t parm1, Long_t parm2) {
  // Handle menu items.

  TCanvas* c;
  int message = GET_MSG(msg);
  //   double     x,y;
  //   int        px, py;
  
  TVisManager* vm = TVisManager::Instance();

  c = GetCanvas();

  switch (message) {
  case kC_COMMAND:
    int submessage = GET_SUBMSG(msg);

    printf(" *** TStnFrame::ProcessMessage SUBMESSAGE: %i \n",submessage);

    switch (submessage) {
    case kCM_MENU:
      switch (parm1) {
 //-----------------------------------------------------------------------------
//  FILE menu
//-----------------------------------------------------------------------------
      case M_FILE_OPEN: {
	TGFileInfo fi;
	fi.fFileTypes = (const char **)filetypes;
	new TGFileDialog(gClient->GetRoot(), this, kFDOpen,&fi);
      }
      break;

      case M_FILE_SAVE:
	printf("M_FILE_SAVE\n");
	break;
	
      case M_FILE_EXIT:
	CloseWindow();   // terminate theApp no need to use SendCloseMessage()
	break;
//-----------------------------------------------------------------------------
//  EDIT menu
//-----------------------------------------------------------------------------
      case M_EDIT_EDITOR:
	EditorBar();
	break;
//-----------------------------------------------------------------------------
//  OPTIONS menu
//-----------------------------------------------------------------------------
      case M_OPTION_EVENT_STATUS:

	printf(" *** TStnFrame::ProcessMessage M_OPTION_EVENT_STATUS: msg = %li parm1 = %li parm2 = %li\n", 
	       msg,parm1,parm2);
	c->ToggleEventStatus();
	if (c->GetShowEventStatus()) {
	  ShowFrame(fStatusBar);
	  fMenuOption->CheckEntry(M_OPTION_EVENT_STATUS);
	} 
	else {
	  HideFrame(fStatusBar);
	  fMenuOption->UnCheckEntry(M_OPTION_EVENT_STATUS);
	}
	break;

      case M_OPTION_CAN_EDIT:
	if (gROOT->GetEditHistograms()) {
	  gROOT->SetEditHistograms(kFALSE);
	  fMenuOption->UnCheckEntry(M_OPTION_CAN_EDIT);
	} 
	else {
	  gROOT->SetEditHistograms(kTRUE);
	  fMenuOption->CheckEntry(M_OPTION_CAN_EDIT);
	}
	break;
//-----------------------------------------------------------------------------
//  OPEN menu
//-----------------------------------------------------------------------------
      case M_OPEN_XY:
	if (vm->DebugLevel() > 0) {
	  printf(" *** TStnFrame::ProcessMessage M_OPEN_XY: msg = %li parm1 = %li parm2 = %li\n", 
		 msg,parm1,parm2);
	}
	TStnVisManager::Instance()->OpenTrkXYView();
	break;
      case M_OPEN_RZ:
	if (vm->DebugLevel() > 0) {
	  printf(" *** TStnFrame::ProcessMessage M_OPEN_RZ: msg = %li parm1 = %li parm2 = %li\n", 
		 msg,parm1,parm2);
	}
	TStnVisManager::Instance()->OpenTrkRZView();
	break;
      case M_OPEN_TZ:
	
	printf(" *** TStnFrame::ProcessMessage M_OPEN_TZ: msg = %li parm1 = %li parm2 = %li\n", 
	       msg,parm1,parm2);
	TStnVisManager::Instance()->OpenTrkTZView();
	break;
      case M_OPEN_CAL:
	
	if (vm->DebugLevel() > 0) {
	  printf(" *** TStnFrame::ProcessMessage M_OPEN_CAL: msg = %li parm1 = %li parm2 = %li\n", 
		 msg,parm1,parm2);
	}
	TStnVisManager::Instance()->OpenCalView();
	break;
      case M_OPEN_CRV:
	
	if (vm->DebugLevel() > 0) {
	  printf(" *** TStnFrame::ProcessMessage M_OPEN_CRV: msg = %li parm1 = %li parm2 = %li\n", 
		 msg,parm1,parm2);
	}
	TStnVisManager::Instance()->OpenCrvView();
	break;
      case M_OPEN_VST:
	
	if (vm->DebugLevel() > 0) {
	  printf(" *** TStnFrame::ProcessMessage M_OPEN_VST: msg = %li parm1 = %li parm2 = %li\n", 
		 msg,parm1,parm2);
	}
	TStnVisManager::Instance()->OpenVSTView();
	break;
//-----------------------------------------------------------------------------
//  default
//-----------------------------------------------------------------------------
      default:
	if (vm->DebugLevel() > 0) {
	  printf(" *** TStnFrame::ProcessMessage msg = %li parm1 = %li parm2 = %li\n", 
		 msg,parm1,parm2);
	}
	break;
      }
    default:
      if (vm->DebugLevel() > 0) {
	printf(" *** TStnFrame::ProcessMessage msg = %li parm1 = %li parm2 = %li\n", 
	       msg,parm1,parm2);
      }
      break;
    }
  }

//   px = gPad->GetEventX();
//   py = gPad->GetEventY();
  
//   x = ((TPad*) gPad)->AbsPixeltoX(px);
//   y = ((TPad*) gPad)->AbsPixeltoY(py);
  
//   fStatusBar->SetText(Form("Z = %8.3f",x),0);
//   fStatusBar->SetText(Form("R = %8.3f",y),1);
//   fStatusBar->SetText("ccc",2);
//   fStatusBar->SetText("ddd",3);
  
  return true;
}

//-----------------------------------------------------------------------------
void TStnFrame::CloseWindow() {
  // Called when window is closed via the window manager.
  
  TVirtualPad *savepad = gPad;
  gPad = 0;        // hide gPad from CINT
  gInterpreter->DeleteGlobal(fEmbeddedCanvas->GetCanvas());
  gPad = savepad;  // restore gPad for ROOT
  
  delete this;
}

//-----------------------------------------------------------------------------
void TStnFrame::ShowStatusBar(Bool_t show) {
  // Show or hide statusbar.

  if (show) {
    ShowFrame(fStatusBar);
    fMenuOption->CheckEntry(M_OPTION_EVENT_STATUS);
  } else {
    HideFrame(fStatusBar);
    fMenuOption->UnCheckEntry(M_OPTION_EVENT_STATUS);
  }
}

//-----------------------------------------------------------------------------
void TStnFrame::SetStatusText(const char *txt, Int_t partidx) {
  // Set text in status bar.

  fStatusBar->SetText(txt, partidx);
}

//-----------------------------------------------------------------------------
void TStnFrame::DoOK() {
  printf("\nTerminating dialog: OK pressed\n");

  // Send a close message to the main frame. This will trigger the
  // emission of a CloseWindow() signal, which will then call
  // TStnFrame::CloseWindow(). Calling directly CloseWindow() will cause
  // a segv since the OK button is still accessed after the DoOK() method.
  // This works since the close message is handled synchronous (via
  // message going to/from X server).
  
  SendCloseMessage();
  
  // The same effect can be obtained by using a singleshot timer:
  //TTimer::SingleShot(50, "TStnFrame", this, "CloseWindow()");
}

//_____________________________________________________________________________
void TStnFrame::DoCancel() {
  printf("\nTerminating dialog: Cancel pressed\n");
  SendCloseMessage();
}

//_____________________________________________________________________________
void TStnFrame::HandleButtons(Int_t id) {
  // Handle different buttons.
}

//-----------------------------------------------------------------------------
void TStnFrame::DoTab(Int_t id) {
   printf("*** TStnFrame::DoTab : Tab item %d activated\n", id);
}

//-----------------------------------------------------------------------------
void TStnFrame::DoRadio() {
   // Handle radio buttons.

  TGButton *btn = (TGButton *) gTQSender;
  Int_t id = btn->WidgetId();

  printf(" TStnFrame::DoRadio radio button ID: %i\n",id);
  
  TStnVisManager* vm = TStnVisManager::Instance();
  if (id >= 21 && id <= 22) {

    if (id == 21) vm->SetDisplayStrawHitsXY(1);
    else          vm->SetDisplayStrawHitsXY(0);

    for (int i = 0; i < 2; i++) {
      if (fRb[i]->WidgetId() != id) fRb[i]->SetState(kButtonUp);
    }
  }
}

//-----------------------------------------------------------------------------
// TGCheckButton changes its state on its own
//-----------------------------------------------------------------------------
void TStnFrame::DoCheckButtons() {
   // Handle radio buttons.

  TGButton *btn = (TGButton *) gTQSender;
  Int_t id = btn->WidgetId();

  printf(" TStnFrame::DoCheckButtons check button ID: %i state: %i\n",id,btn->IsOn());
  
  TStnVisManager* vm = TStnVisManager::Instance();

  int doit = (int) btn->IsOn(); 

  if      (id == kDisplayHelices     ) vm->SetDisplayHelices     (doit);
  else if (id == kDisplayTracks      ) vm->SetDisplayTracks      (doit);
  else if (id == kDisplaySimParticles) vm->SetDisplaySimParticles(doit);

}
