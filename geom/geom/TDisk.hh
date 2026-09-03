//
// Read the tracks added to the event by the track finding and fitting code.
//
// $Id: TDisk.hh,v 1.1 2014/06/13 06:14:48 murat Exp $
// $Author: murat $
// $Date: 2014/06/13 06:14:48 $
//
// Contact person Pavel Murat
//
#ifndef Stntuple_geom_TDisk_hh
#define Stntuple_geom_TDisk_hh

// storable objects (data products)
// 

// C++ includes.
#include <iostream>
#include <math.h>

#include "TMath.h"
#include "TString.h"
#include "TFolder.h"
#include "TFile.h"
#include "TObjArray.h"
#include "TVector2.h"

#include "Stntuple/geom/TDiskIndex.hh"
#include "Stntuple/geom/TDiskCrystalMap.hh"

class TStnCrystal;

namespace mu2e {
  class Disk;
}

class TDisk : public TObject {
public:
  int                fSectionID;
					
  TDiskCrystalMap*   fCrystalMap;
  TObjArray*         fListOfCrystals;	// a sparse array of crystals 

  int                fNEdges;		// 6:hex or 4:square crystals

  int                fFirstChanOffset;	// offset of the first channel

  double             fRMin;
  double             fRMax;
  double             fZ0;
  double             fSize;             // hex size for hex crystals, side - for square
  double             fDeadSpace;	// wrapper+shell
  double             fMinFraction;	// min fraction of area to add a crystal

  double             fEnergy;           // total energy deposited in the disk
//-----------------------------------------------------------------------------
// methods
//-----------------------------------------------------------------------------
  TDisk();
  TDisk(int      SectionID  , 
	double   RMin       , 
	double   RMax       , 
	double   Z0         ,
	int      NEdges     ,
	double   Size       , 
	double   DeadSpace  , 
	double   MinFraction);

  ~TDisk();
//-----------------------------------------------------------------------------
// accessors
//-----------------------------------------------------------------------------
  int            SectionID() { return fSectionID; }
  TObjArray*     ListOfCrystals() { return fListOfCrystals; }
  int            First    (int Ir) { return fCrystalMap->GetFirst(Ir); }
  TStnCrystal*   Crystal  (int I ) { return (TStnCrystal*) fListOfCrystals->UncheckedAt(I); }

  void           GetHexIndex(int I, TDiskIndex* HexIndex);
					// for crystal number I return its ring number
  int            GetRing (int I);
  int            GetRing (TDiskIndex* Index);

  double         GetRMin() const { return fRMin; }
  double         GetRMax() const { return fRMax; }

  int            GetNRings   () const { return fCrystalMap->GetNRings();    }
  int            NCrystals() const { return fCrystalMap->GetNInsideTotal(); }

  void           GetPosition(int I           , TVector2* Pos);
  void           GetPosition(TDiskIndex* Index, TVector2* Pos);

  double         Energy() { return fEnergy; }

  int            GetNCrystalsPerRing(int Ring) const { return fCrystalMap->GetNCrystalsPerRing(Ring); }
  //  int            GetNInside         (int Ring) const { return fCrystalMap->fNInside         [Ring]; }

  double         GetActiveArea() { return GetCrystalArea()*NCrystals(); }
  double         GetTotalArea () { return (fRMax*fRMax-fRMin*fRMin)*TMath::Pi() ; }

  double         GetRadius(int I);
  double         GetRadius(TDiskIndex* Index);

  double         GetCrystalArea() { return  fSize*fSize*sqrt(3.)/2.; }

// two crystal cross sections require different initialization

  int            InitHexCrystals();
  int            InitSqrCrystals();
  int            InitSsqCrystals();

  int            FirstChanOffset() { return fFirstChanOffset; }
//-----------------------------------------------------------------------------
// setters
//-----------------------------------------------------------------------------
  void           SetFirstChanOffset(int Offset) { fFirstChanOffset = Offset; }

  void           SetEnergy(double E) { fEnergy = E; }
//-----------------------------------------------------------------------------
// overloaded functions of TObject
//-----------------------------------------------------------------------------
  virtual void   Paint(Option_t* Opt = "");
  virtual void   Clear(Option_t* Opt = "") ;
  virtual void   Print(Option_t* Opt = "") const ;

  ClassDefOverride(TDisk,0)
    
};

#endif
