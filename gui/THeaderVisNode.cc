//_____________________________________________________________________________
// Feb 10 2001 P.Murat
//_____________________________________________________________________________
#include "TVirtualX.h"
#include "TPad.h"

#include "Stntuple/gui/THeaderVisNode.hh"
#include "Stntuple/gui/TStnVisManager.hh"
#include "Stntuple/obj/TStnHeaderBlock.hh"

ClassImp(THeaderVisNode)

//_____________________________________________________________________________
THeaderVisNode::THeaderVisNode(const char* name, TStnHeaderBlock* h): 
  TStnVisNode(name) 
{
  fHeader = h;
  fText = new TText(0.1,0.1,"elki-palki");
  fText->SetTextFont(52);
}

//_____________________________________________________________________________
THeaderVisNode::~THeaderVisNode() {
  delete fText;
}


//_____________________________________________________________________________
int THeaderVisNode::InitEvent() {
  //
  return 0;
}

//-----------------------------------------------------------------------------
// override Paint
//-----------------------------------------------------------------------------
void THeaderVisNode::Paint(Option_t* option) {
  PaintXY(option);
}

//_____________________________________________________________________________
void THeaderVisNode::PaintXY(Option_t* option) {
  // draw event/run

  TStnVisManager* vm = TStnVisManager::Instance();

  const art::Event* ev = vm->Event();
  char  text[100];

  if (ev != nullptr) sprintf(text,"Event: %7i:%06i:%7i",ev->run(),ev->subRun(),ev->event());
  else               sprintf(text,"Event: %7i:%06i:%7i",0,0,0);

  fText->SetTextSize(0.4);
  fText->SetText(0.02,0.3,text);

  fText->Paint(option);
  gPad->Modified();
}


//_____________________________________________________________________________
void THeaderVisNode::PaintRZ(Option_t* option) {
}
//-----------------------------------------------------------------------------
void THeaderVisNode::PaintTZ(Option_t* Option) {
}
//-----------------------------------------------------------------------------
void THeaderVisNode::PaintPhiZ(Option_t* Option) {
}
//-----------------------------------------------------------------------------
void THeaderVisNode::PaintCal(Option_t* Option) {
}
//-----------------------------------------------------------------------------
void THeaderVisNode::PaintCrv(Option_t* Option) {
}
//-----------------------------------------------------------------------------
void THeaderVisNode::PaintVST(Option_t* Option) {
}
//-----------------------------------------------------------------------------
void THeaderVisNode::PaintVRZ(Option_t* Option) {
}
//_____________________________________________________________________________
Int_t THeaderVisNode::DistancetoPrimitive(Int_t px, Int_t py) {
  return 9999;
}

//_____________________________________________________________________________
Int_t THeaderVisNode::DistancetoPrimitiveXY(Int_t px, Int_t py) {
  return 9999;
}

//_____________________________________________________________________________
Int_t THeaderVisNode::DistancetoPrimitiveRZ(Int_t px, Int_t py) {
  return 9999;
}

