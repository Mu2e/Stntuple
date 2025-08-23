#ifndef __Stntuple_gui_TStnVisNode_hh__
#define __Stntuple_gui_TStnVisNode_hh__

#include "TObject.h"
#include "TString.h"

#include "Stntuple/base/TVisNode.hh"

class TStnVisNode: public TVisNode {
protected:
  int fID;                             // just to test...
public:
					// ****** constructors and destructor
  TStnVisNode(const char* name = "TStnVisNode");
  virtual ~TStnVisNode();
//-----------------------------------------------------------------------------
// accessors
//-----------------------------------------------------------------------------
  int GetID() { return fID; }
//-----------------------------------------------------------------------------
// these are the views. Each view has its own view manager handling multiple, potentially,
// windows with this view
//-----------------------------------------------------------------------------
  virtual void  Paint    (Option_t* option = "");
  virtual void  PaintXY  (Option_t* option = "");
  virtual void  PaintRZ  (Option_t* option = "");
  virtual void  PaintTZ  (Option_t* option = "");
  virtual void  PaintPhiZ(Option_t* option = "");
  virtual void  PaintCal (Option_t* option = "");
  virtual void  PaintCrv (Option_t* option = "");
  virtual void  PaintVST (Option_t* option = "");
  virtual void  PaintVRZ (Option_t* option = "");
  virtual void  PaintVYZ (Option_t* option = "");

  virtual int   DistancetoPrimitive    (Int_t px, Int_t py);
  virtual int   DistancetoPrimitiveXY  (Int_t px, Int_t py);
  virtual int   DistancetoPrimitiveRZ  (Int_t px, Int_t py);
  virtual int   DistancetoPrimitiveTZ  (Int_t px, Int_t py);
  virtual int   DistancetoPrimitivePhiZ(Int_t px, Int_t py);
  virtual int   DistancetoPrimitiveCal (Int_t px, Int_t py);
  virtual int   DistancetoPrimitiveCrv (Int_t px, Int_t py);
  virtual int   DistancetoPrimitiveVST (Int_t px, Int_t py);
  virtual int   DistancetoPrimitiveVRZ (Int_t px, Int_t py);
  virtual int   DistancetoPrimitiveVYZ (Int_t px, Int_t py);
};

#endif
