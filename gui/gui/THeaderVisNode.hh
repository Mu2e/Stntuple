#ifndef THeaderVisNode_hh
#define THeaderVisNode_hh

#include "Gtypes.h"
#include "TText.h"

#include "Stntuple/gui/TStnVisNode.hh"

class TStnHeaderBlock;

class THeaderVisNode: public TStnVisNode {

protected:
  TStnHeaderBlock*    fHeader;
  TText*              fText;

public:
					// ****** constructors and destructor

  THeaderVisNode(const char* name = "",
		 TStnHeaderBlock* h = 0);

  virtual ~THeaderVisNode();
					// ****** accessors

					// ****** modifiers

  int   InitEvent();

  virtual void Paint  (Option_t* option = "");

  virtual void  PaintXY  (Option_t* option = "");
  virtual void  PaintRZ  (Option_t* option = "");
  virtual void  PaintTZ  (Option_t* option = "");
  virtual void  PaintPhiZ(Option_t* option = "");
  virtual void  PaintCal (Option_t* option = "");
  virtual void  PaintCrv (Option_t* option = "");
  virtual void  PaintVST (Option_t* option = "");
  virtual void  PaintVRZ (Option_t* option = "");
  
  virtual Int_t DistancetoPrimitive  (Int_t px, Int_t py);
  virtual Int_t DistancetoPrimitiveXY(Int_t px, Int_t py);
  virtual Int_t DistancetoPrimitiveRZ(Int_t px, Int_t py);

  ClassDef(THeaderVisNode,0)
};


#endif
