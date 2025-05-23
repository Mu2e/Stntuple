//-----------------------------------------------------------------------------
//  Jan 04 2001 P.Murat: base class for Geant visualization manager
//  THis class also provides low-level graphics routines to resolve
//  references from FORTRAN part of GEANT visualization package
//-----------------------------------------------------------------------------
#ifndef Stntuple_base_TVisManager_hh
#define Stntuple_base_TVisManager_hh

#include "TList.h"
#include "TCanvas.h"
#include "TObjArray.h"

#include "Stntuple/base/TVisNode.hh"
#include "Stntuple/base/TStnView.hh"

class TGeometryManager;
class TDetectorElement;

class TVisManager: public TNamed {
//-----------------------------------------------------------------------------
//  static class members
//-----------------------------------------------------------------------------
protected:
  static TVisManager*  fgInstance;
//-----------------------------------------------------------------------------
//  data members
//-----------------------------------------------------------------------------
protected:
  TGeometryManager* fGeometryManager;   // ! geometry (do we need to know? )
  TList*            fListOfCanvases;	// ! list of canvases
  TCanvas*          fActiveCanvas;	// ! pointer to the active canvas

  TObjArray*        fListOfViews;       // ! list of views
  TObjArray*        fListOfNodes;       // ! list of vis nodes
  TObject*          fClosestObject;	// ! object closest to the cursor
  Int_t             fMinDist;		// ! distance to the closest object
  TDetectorElement* fClosestDetElement;	// !
  Int_t             fMinDistDetElement;	// !
  TVisNode*         fTitleNode;         // !
  TStnView*         fCurrentView;       // ! 
  int               fDebugLevel; 	// !
//-----------------------------------------------------------------------------
// methods
//-----------------------------------------------------------------------------
public:
					// ****** constructors and destructor

  TVisManager(const char* name  = "TVisManager", 
	      const char* title = "TVisManager"); 

  virtual ~TVisManager();

  static TVisManager* Instance() {
    return (fgInstance) ? fgInstance : fgInstance = new TVisManager();
  }
					// ****** registry
  void  DeclareCanvas(TCanvas* c);
					// ****** accessors

  TGeometryManager*  GetGeometryManager() { return fGeometryManager; }

  TCanvas*           GetCanvas() { return fActiveCanvas; }

  TList*             GetListOfCanvases   () { return fListOfCanvases; }
  TObjArray*         GetListOfNodes      () { return fListOfNodes; }
  TObject*           GetClosestObject    () { return fClosestObject; }
  Int_t              GetMinDist          () { return fMinDist; }
  Int_t              GetMinDistDetElement() { return fMinDistDetElement; }
  TDetectorElement*  GetClosestDetElement() { return fClosestDetElement; }
  TStnView*          GetCurrentView      () { return fCurrentView;       }

  Int_t              GetNNodes           () { 
    return fListOfNodes->GetEntriesFast(); 
  }

  void  AddNode(TVisNode* node) { 
				        // make sure we're not adding the same node twice
    if (fListOfNodes->FindObject(node) == 0) {
      fListOfNodes->Add(node); 
    }
  }

  TVisNode*          GetNode      (Int_t i) { return (TVisNode*) fListOfNodes->UncheckedAt(i); }

  TVisNode*          FindNode     (const char* Name) { return (TVisNode*) fListOfNodes->FindObject(Name); }

  int                GetNViews()    { return fListOfViews->GetEntriesFast(); }
  TStnView*          GetView(int I) { return (TStnView*) fListOfViews->UncheckedAt(I); } 

  void               AddView(TStnView* View);

  TStnView*          FindView(int Type, int Index);

  int                DebugLevel() { return fDebugLevel; }
//-----------------------------------------------------------------------------
// modifiers
//-----------------------------------------------------------------------------
  void  SetGeometryManager(TGeometryManager* gm) { fGeometryManager = gm; }

  void  SetCurrentView(TStnView* View) { fCurrentView = View; }

  void  SetClosestObject    (TObject* o, Int_t Dist) { 
    fClosestObject = o; fMinDist = Dist;
  }

  void  SetClosestDetElement(TDetectorElement* det, Int_t Dist) ;

  void  SetTitleNode (TVisNode* node) { fTitleNode = node; }
  void  SetDebugLevel(int      Level) { fDebugLevel = Level; }


  virtual int  GetViewID     (const char* View); // to be overloaded
  virtual void SetStations   (int IMin, int IMax);

  virtual void SetTimeWindow (float TMin, float TMax);
  virtual void SetMinEDep    (float E);
  virtual void SetMaxEDep    (float E);
//-----------------------------------------------------------------------------
// GUI comands / drawing functions
//-----------------------------------------------------------------------------
  virtual void InitEvent   ();		// self-initialization, called from DisplayEvent
  virtual void DisplayEvent();

  virtual void DefaultRange ();
  virtual void DoCheckButton();
  virtual void DoRadioButton();

  virtual void DrawOneSpec(const char* name);

					// returns 1 if end run
  virtual int  EndRun();

  virtual TCanvas*       NewCanvas(const char* Name, 
				   const char* Title, 
				   Int_t       SizeX, 
				   Int_t       SizeY) { 
    return NULL;
  }

  virtual void OpenView(const char* View);
  virtual void OpenView(TStnView* Mother, int Px1, int Py1, int Px2, int Py2);

					// update all the open windows

  virtual void MarkModified(TPad* Pad);
  virtual void Quit        ();   //
//-----------------------------------------------------------------------------
// Geant3 legacy 
//-----------------------------------------------------------------------------
  virtual void Gdclose();
  virtual void Gdelete(Int_t view);
  virtual void Gdhead(Int_t isel, const char* name, Float_t csize=0.6);   
  virtual void Gdman(Float_t u0, Float_t v0, const char* type="MAN");
  virtual void Gdopen(Int_t view);
  virtual void Gdraw(const char* name,
		     Float_t theta = 30  , Float_t phi = 30, Float_t psi = 0, 
		     Float_t u0    = 10  , Float_t v0  = 10,
		     Float_t ul    = 0.01, Float_t vl = 0.01); //*MENU*

  virtual void Gdrawc(const char *name,
		      Int_t axis=1, Float_t cut=0,Float_t u0=10,
		      Float_t v0=10,Float_t ul=0.01,Float_t vl=0.01); // *MENU*

  void Gdrawx(const char *name,Float_t cutthe, Float_t cutphi,
	      Float_t cutval, Float_t theta, Float_t phi, Float_t u0,
	      Float_t v0,Float_t ul,Float_t vl); // *MENU*

  virtual void Gdopt(const char* name,const char* opt);	// *MENU*
  virtual void Gdshow(Int_t view);
  virtual void Gdspec(const char *name);                // *MENU*
  virtual void Gdtree(const char *name,Int_t levmax=15,Int_t ispec=0); // *MENU*
  virtual void GdtreeParent(const char *name,Int_t levmax=15,Int_t ispec=0);
  virtual void Gdxyz(Int_t it); 
  virtual void Gdcxyz(); 

  virtual void Gsatt(const char* name, const char* att, Int_t val); //*MENU*

  virtual void SetBomb(Float_t boom = 1.);

  virtual void SetClipBox(const char* name, 
			  Float_t xmin=-9999, Float_t xmax=0, 
			  Float_t ymin=-9999, Float_t ymax=0,
			  Float_t zmin=-9999, Float_t zmax=0);

  virtual void SetColors() {}
					// ****** dummies

  virtual void Divide(Int_t nx=1, Int_t ny=1, 
		      Float_t xmargin=0.01, Float_t ymargin=0.01, 
		      Int_t color=0) {};


  virtual void SetGrid (Int_t valuex = 1, Int_t valuey = 1) {}
  virtual void SetGridx(Int_t value = 1) {}
  virtual void SetGridy(Int_t value = 1) {}
  virtual void SetLogx (Int_t value = 1) {}
  virtual void SetLogy (Int_t value = 1) {}
  virtual void SetLogz (Int_t value = 1) {}
  virtual void SetTickx(Int_t value = 1) {}
  virtual void SetTicky(Int_t value = 1) {}

  ClassDef(TVisManager,1)		// Abs visualization manager for TGeant
};

#endif
