///////////////////////////////////////////////////////////////////////////////
// 2D view 
///////////////////////////////////////////////////////////////////////////////
#ifndef Stntuple_base_TStnView_hh
#define Stntuple_base_TStnView_hh

#include "TMarker.h"
#include "TNamed.h"
#include "TObjArray.h"
#include "TGaxis.h"
#include "TGeoMatrix.h"
#include "TVector3.h"

#include "Stntuple/base/TVisNode.hh"

class TStnView: public TNamed {
protected:
  int                 fType;  // view type
  int                 fIndex; // for calorimeter - 2 views, for example

  int                 fDebugLevel;
  
  TGeoCombiTrans      fCombiTrans;   // rotate, then translate
  TVector3            fUDir;         // temporary
  TVector3            fVDir;
  TVector3            fWDir;
  
  Int_t               fPx1;
  Int_t               fPy1;
  Int_t               fPx2;
  Int_t               fPy2;
					// displayed range in world coordinates
  double              fXMin;
  double              fXMax;
  double              fYMin;
  double              fYMax;

  TGaxis*             fXAxis;
  TGaxis*             fYAxis;

  TMarker*            fCenter;

  int                 fCursorX;
  int                 fCursorY;

  TObjArray*          fListOfNodes;	// list of nodes to visualize

protected:
  float     fTMin;
  float     fTMax;

public:
  TStnView(int Type = -1, int Index = -1); 

  TStnView(int Type, int Index, const char* Name, const char* Title);

  virtual ~TStnView();
//-----------------------------------------------------------------------------
// accessors
//-----------------------------------------------------------------------------
  int           Type () { return fType;  }
  int           Index() { return fIndex; }

  TGeoCombiTrans*  GetCombiTrans() { return &fCombiTrans; }

  int           GetNNodes()      { return fListOfNodes->GetEntriesFast(); }
  TVisNode*     GetNode  (int I) { return (TVisNode*) fListOfNodes->UncheckedAt(I);   }
  TObjArray*    GetListOfNodes() { return fListOfNodes; }

  TVector3*     UDir() { return & fUDir; }
  TVector3*     VDir() { return & fVDir; }
  TVector3*     WDir() { return & fWDir; }

  void          AddNode(TVisNode* Node) { fListOfNodes->Add(Node); }
//-----------------------------------------------------------------------------
// setters
//-----------------------------------------------------------------------------
  void          SetType       (int Type ) { fType  = Type;  } 
  void          SetIndex      (int Index) { fIndex = Index; } 

  void          SetTimeWindow (float TMin, float TMax);  // *MENU* 
  void          SetMinEDep    (float E);                 // *MENU* 
  void          SetMaxEDep    (float E);                 // *MENU* 
  void          SetStations   (int I1, int I2);          // *MENU* 
  void          SetDebugLevel (int Level);               // *MENU*
//-----------------------------------------------------------------------------
// overloaded functions of TObject
//-----------------------------------------------------------------------------
  virtual void  Paint               (Option_t* option = "");
  virtual Int_t DistancetoPrimitive (Int_t px, Int_t py);
  virtual void  ExecuteEvent        (Int_t event, Int_t px, Int_t py);
  virtual char* GetObjectInfo       (Int_t Px, Int_t Py) const ;

  virtual void  Print               (Option_t* option = "") const ;  // *MENU* 

  ClassDef(TStnView,0)
};

#endif
