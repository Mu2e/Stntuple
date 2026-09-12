///////////////////////////////////////////////////////////////////////////////
// vis node displays one wedge
///////////////////////////////////////////////////////////////////////////////
#ifndef TStnCrystal_hh
#define TStnCrystal_hh

#include "Gtypes.h"
#include "TObjArray.h"
#include "TH1.h"
#include "TPad.h"
#include "TArc.h"
#include "TVector3.h"

#include "Stntuple/geom/TDiskIndex.hh"
#include "Stntuple/base/TStnShape.hh"
#include "Stntuple/obj/TCalHitData.hh"

class TDisk;

class TStnCrystal: public TObject {
public:
  
protected:
  int                   fIndex;
  TDiskIndex            fDiskIndex;

  TVector3              fCenter;
  double                fSize;

  TObjArray*            fListOfHits;    // keeps [not-owned] pointers to TCalDataBlock
					// display in XY view
  TStnShape*            fShape;

  int                   fNEdges;
  int                   fFillStyle;
  int                   fFillColor;
  int                   fLineColor;
  int                   fNHits;		// number of (resolved) crystal hits

  float                 fEnergy;  	// total deposited energy
  TDisk*                fDisk;          // backward pointer to the disk

public:
//-----------------------------------------------------------------------------
// constructors and destructor
//-----------------------------------------------------------------------------
  TStnCrystal() {}
  TStnCrystal(TDiskIndex* Index, double X0, double Y0, double Z0, int NEdges, double Size); 

  virtual ~TStnCrystal();
//-----------------------------------------------------------------------------
// accessors
//-----------------------------------------------------------------------------
  TObjArray*     ListOfHits   () { return fListOfHits; }
  TCalHitData*   CalHitData(int I) { return (TCalHitData*) fListOfHits->UncheckedAt(I); }
  TVector3*      Center       () { return &fCenter   ; }
  int            Index        () { return fIndex     ; }
  int            NHits        () { return fNHits; }
  float          Energy       () { return fEnergy;     }
  double         Radius       () { return fShape->Radius(); }
  double         Size         () { return fSize; }

  //  const mu2e::Crystal* Crystal() { return fCrystal;    }
  TDisk*         Disk         () { return fDisk; }

  TDiskIndex*    DiskIndex     () { return &fDiskIndex; }
//-----------------------------------------------------------------------------
// modifiers
//-----------------------------------------------------------------------------
  void  SetFillStyle(int Style) { fShape->fFillStyle = Style; }
  void  SetFillColor(int Color) { fShape->fFillColor = Color; }
  void  SetLineColor(int Color) { fShape->fLineColor = Color; }
  void  SetDisk     (TDisk* Disk) { fDisk = Disk; }

  void   AddHit(TCalHitData* CalHit) { 
    fListOfHits->Add(CalHit); 
    fEnergy += CalHit->Energy();
    fNHits++;
  }
					// use for rad damage estimates

  void   AddEnergy(float E) { fEnergy += E; }

  //  virtual void  Draw    (Option_t* option = "");

  virtual void  Paint   (Option_t* option = "") override;
  virtual void  PaintXY (Option_t* option = "");
  virtual void  PaintCal(Option_t* option = "");

  //  virtual void  ExecuteEvent(Int_t event, Int_t px, Int_t py);

  virtual Int_t DistancetoPrimitive  (Int_t px, Int_t py) override;
  virtual Int_t DistancetoPrimitiveXY(Int_t px, Int_t py);
  virtual Int_t DistancetoPrimitiveRZ(Int_t px, Int_t py);
//-----------------------------------------------------------------------------
// overloaded methods of TObject
//-----------------------------------------------------------------------------
  virtual void   Clear(const char* Opt = "") override;
  virtual void   Print(const char* Opt = "") const override; // **MENU**

  ClassDefOverride(TStnCrystal,0)
};


#endif
