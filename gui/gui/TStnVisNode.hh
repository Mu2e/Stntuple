#ifndef Stntuple_gui_TStnVisVode_hh
#define Stntuple_gui_TStnVisVode_hh

#include "TObject.h"
#include "TString.h"

#include "Stntuple/base/TVisNode.hh"

class TStnVisNode: public TVisNode {
protected:
  // int        fSectionToDisplay;
  // int        fDiskID;                   // calorimenter disk ID 
public:
					// ****** constructors and destructor
  TStnVisNode(const char* name = "TStnVisNode");
  virtual ~TStnVisNode();
//-----------------------------------------------------------------------------
// accessors
//-----------------------------------------------------------------------------
//  int   SectionToDisplay() { return fSectionToDisplay; }

//  int   DebugLevel      () { return fDebugLevel;       }

					// called by TEvdManager::DisplayEvent
  virtual int   InitEvent();
//-----------------------------------------------------------------------------
// these are the views. Each view has its own view manager handling multiple, potentially,
// windows with this view
//-----------------------------------------------------------------------------
  virtual void  Paint   (Option_t* option = "");
  virtual void  PaintXY (Option_t* option = "");
  virtual void  PaintRZ (Option_t* option = "");
  virtual void  PaintTZ (Option_t* option = "");
  virtual void  PaintCal(Option_t* option = "");
  virtual void  PaintCrv(Option_t* option = "");

  virtual int   DistancetoPrimitive   (Int_t px, Int_t py);
  virtual int   DistancetoPrimitiveXY (Int_t px, Int_t py);
  virtual int   DistancetoPrimitiveRZ (Int_t px, Int_t py);
  virtual int   DistancetoPrimitiveTZ (Int_t px, Int_t py);
  virtual int   DistancetoPrimitiveCal(Int_t px, Int_t py);
  virtual int   DistancetoPrimitiveCrv(Int_t px, Int_t py);

  // void SetSectionToDisplay(int Section) { fSectionToDisplay= Section; }

  void SetDebugLevel(int Level) { fDebugLevel = Level; }

  ClassDef(TStnVisNode,0)
};

#endif
