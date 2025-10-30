//-----------------------------------------------------------------------------
//  Jan 04 2001 P.Murat: cut and paste to create a base class for 
//  TGEANT visualization manager
//
//-----------------------------------------------------------------------------
#include "TROOT.h"
#include "TGraph.h"
#include "TColor.h"
#include "TLine.h"
#include "TPolyLine.h"
#include "TPolyMarker.h"
#include "TPaveLabel.h"
#include "TText.h"
#include "TEnv.h"
#include <string.h> 
#include <ctype.h> 
#include "TCanvas.h"

//#include "TG3Volume.hh"
// #include "TGeometryManager.hh"
#include "Stntuple/base/TVisManager.hh"
// #include "TDetectorElement.hh"


TVisManager*  TVisManager::fgInstance = 0;

ClassImp(TVisManager)
//_____________________________________________________________________________
TVisManager::TVisManager(const char* name, const char* title): 
  TNamed(name,title)
{
  fListOfCanvases = new TList();
  fListOfNodes    = new TObjArray();
  fListOfViews    = new TObjArray();
  fgInstance      = this;
  fDebugLevel     = gEnv->GetValue("TVisManager.DebugLevel",0);
  fTitleNode      = nullptr;
}


//_____________________________________________________________________________
TVisManager::~TVisManager() {
  fListOfNodes->Delete();
  delete fListOfNodes;

  fListOfViews->Delete();
  delete fListOfViews;

  fListOfCanvases->Delete();
  delete fListOfCanvases;

  if (fTitleNode) delete fTitleNode;

  fgInstance  = 0;
}

//-----------------------------------------------------------------------------
// AddNode checks for the presence of the node being present
//-----------------------------------------------------------------------------
void TVisManager::AddView(TStnView* View) { 
  if (fListOfViews->FindObject(View) == 0) {
    fListOfViews->Add(View); 
    
    int n = View->GetNNodes();
    for (int i=0; i<n; i++) {
      TVisNode* node = View->GetNode(i);
      AddNode(node);
    }
  }
}


//-----------------------------------------------------------------------------
// dist - in points
//-----------------------------------------------------------------------------
void TVisManager::SetClosestDetElement(TDetectorElement* Det, Int_t Dist) {
  fClosestDetElement = Det;
  fMinDistDetElement = Dist;
  if (Dist < fMinDist) {
					// new closest object
    SetClosestObject((TObject*) Det,Dist);
  }
}


//-----------------------------------------------------------------------------
void TVisManager::DeclareCanvas(TCanvas* c) { 
  // add new canvas to the list of canvases and populate its list of primitives

  fListOfCanvases->Add(c);
}

//-----------------------------------------------------------------------------
TStnView* TVisManager::FindView(int Type, int Index) { 
  TStnView* v(nullptr);
  int nv = GetNViews();
  for (int i=0; i<nv; i++) {
    TStnView* vi = (TStnView*) GetView(i);
    if ((vi->Type() == Type) and ((Index < 0) or (vi->Index() == Index))) {
      v = vi;
      break;
    }
  }
  return v;
}

//-----------------------------------------------------------------------------
void TVisManager::MarkModified(TPad* Pad) {
  Pad->Modified();
  
  TObject     *obj;
  TObjOptLink *lnk = 0;

  TList       *pList = Pad->GetListOfPrimitives();
  if (pList) lnk = (TObjOptLink*)pList->FirstLink();
  
  while (lnk) {
    obj = lnk->GetObject();
    if (obj->InheritsFrom(TPad::Class())) {
      MarkModified((TPad*) obj);
    }
    lnk = (TObjOptLink*)lnk->Next();
  }
}

//-----------------------------------------------------------------------------
// to be overloaded
//-----------------------------------------------------------------------------
int TVisManager::GetViewID(const char* Name) {
  return -1;
}

//-----------------------------------------------------------------------------
void TVisManager::OpenView(TStnView* Mother, int Px1, int Py1, int Px2, int Py2) {
}

//-----------------------------------------------------------------------------
void TVisManager::SetStations(int IMin, int IMax) {
}

//-----------------------------------------------------------------------------
void TVisManager::SetTimeWindow(float TMin, float TMax) {
}

//-----------------------------------------------------------------------------
void TVisManager::SetMinEDep(float E) {
}
//-----------------------------------------------------------------------------
void TVisManager::SetMaxEDep(float E) {
}

//-----------------------------------------------------------------------------
// to be overloaded, if needed, and there already is a need for that
//-----------------------------------------------------------------------------
void TVisManager::InitEvent() { 
}

//-----------------------------------------------------------------------------
// redraw all open windows
//-----------------------------------------------------------------------------
void TVisManager::DisplayEvent() { 
//-----------------------------------------------------------------------------
// reinitialize graphics objects in the event
// vis manager may need some specific initialization, add a virtual function for that
//-----------------------------------------------------------------------------
  InitEvent();

  if (fTitleNode) fTitleNode->InitEvent();

  int nn = fListOfNodes->GetEntries();
  for (int i=0; i<nn; i++) {
    TVisNode* vn = (TVisNode*) fListOfNodes->At(i);
    vn->InitEvent();
  }

//-----------------------------------------------------------------------------
// mark all views as modified
//-----------------------------------------------------------------------------
  TIter it(fListOfCanvases);
  while (TCanvas* c = (TCanvas*) it.Next()) {
//     printf("TVisManager::DisplayEvent CANVAS: class=%-30s name=%-30s\n",
// 	   c->ClassName(),c->GetName());
    TIter it1(c->GetListOfPrimitives());
    while(TObject* o = it1.Next()) {
//       printf("TVisManager::DisplayEvent object: class=%-30s name=%-30s\n",
// 	     o->ClassName(),o->GetName());
      if (o->InheritsFrom("TPad")) {
	TPad* pad = (TPad*) o;
	
	MarkModified(pad);
      }
    }
    c->Modified();
    c->Update();
  }
}

//-----------------------------------------------------------------------------
// returns 1 if end-of-run and don't need to redraw
//-----------------------------------------------------------------------------
int TVisManager::EndRun() { 
  return 0;
}

//_____________________________________________________________________________
void TVisManager::Gdhead(Int_t isel, const char *name, Float_t chrsiz) { 
}


//_____________________________________________________________________________
void TVisManager::Gdman(Float_t u, Float_t v, const char *type) { 
}
 
//_____________________________________________________________________________
void TVisManager::Gdtree(const char *name,Int_t levmax, Int_t isel) { 
} 

//_____________________________________________________________________________
void TVisManager::GdtreeParent(const char *name,Int_t levmax, Int_t isel)
{ 
  //
  //  NAME   Volume name
  //  LEVMAX Depth level
  //  ISELT  Options
  //
  //  This function draws the logical tree of the parent of name.
  //  

  Clear();
//   TG3Shape* mother = fGeometryManager->GetMother(name);
//   if (mother) {
//     Gdtree(mother->GetName(),levmax,isel);
//   }
} 


//____________________________________________________________________________ 
void TVisManager::DefaultRange()
{ 
}


//_____________________________________________________________________________
void TVisManager::Gsatt(const char *name, const char *att, Int_t val)
{ 
}

//_____________________________________________________________________________
void TVisManager::Gdopen(Int_t iview) { 
}
 
//_____________________________________________________________________________
void TVisManager::Gdclose() { 
}
 
//_____________________________________________________________________________
void TVisManager::Gdshow(Int_t iview) { 
} 

//_____________________________________________________________________________
void TVisManager::Gdopt(const char *name,const char *value) { 
} 
 
//_____________________________________________________________________________
void TVisManager::Gdraw(const char *name,
		  Float_t theta, Float_t phi, Float_t psi,
		  Float_t u0,Float_t v0,Float_t ul,Float_t vl)
{ 
}

 
//_____________________________________________________________________________
void TVisManager::Gdrawc(const char *name,Int_t axis, Float_t cut,Float_t u0,
			 Float_t v0,Float_t ul,Float_t vl) 
{ 
} 
 
//_____________________________________________________________________________
void TVisManager::Gdrawx(const char *name,Float_t cutthe, Float_t cutphi,
			 Float_t cutval, Float_t theta, Float_t phi, Float_t u0,
			 Float_t v0,Float_t ul,Float_t vl)
{ 
}
 

//_____________________________________________________________________________
void TVisManager::Gdelete(Int_t iview) { 
}
 
//_____________________________________________________________________________
void  TVisManager::Gdxyz(Int_t it) {
}

//_____________________________________________________________________________
void  TVisManager::Gdcxyz() {
}

//_____________________________________________________________________________
void TVisManager::Gdspec(const char *name) { 
} 

//_____________________________________________________________________________
void TVisManager::DoCheckButton() { 
} 

//_____________________________________________________________________________
void TVisManager::DoRadioButton() { 
} 

//_____________________________________________________________________________
void TVisManager::OpenView(const char* View) { 
  printf("TVisManager::OpenView is meant to be overloaded, view=\"%s\"\n",View);
}

//_____________________________________________________________________________
void TVisManager::Quit() { 
  printf("TVisManager::Quit is meant to be overloaded\n");
}

//_____________________________________________________________________________
void TVisManager::DrawOneSpec(const char *name) { 
} 


//______________________________________________________________________________
void TVisManager::SetBomb(Float_t boom) {
}
 
//_____________________________________________________________________________
void  TVisManager::SetClipBox(const char *name, 
			      Float_t xmin, Float_t xmax,
			      Float_t ymin, Float_t ymax,
			      Float_t zmin, Float_t zmax) 
{
} 
