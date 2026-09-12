//
// Read the tracks added to the event by the track finding and fitting code.
//
// $Id: TDiskCalorimeter.hh,v 1.1 2014/06/13 06:14:48 murat Exp $
// $Author: murat $
// $Date: 2014/06/13 06:14:48 $
//
// Contact person Pavel Murat
//
#ifndef Stntuple_inc_TDiskCalorimeter_hh
#define Stntuple_inc_TDiskCalorimeter_hh

// storable objects (data products)
// 

// C++ includes.
#include <iostream>

#include "TString.h"
#include "TFolder.h"
#include "TFile.h"

//namespace murat {

class TDisk;
class TCalDataBlock;

class TDiskCalorimeter : public TObject {
public:

  struct GeomData_t {
    int    fNDisks;
    int    fNEdges;
    int    fNCrystals[4];
    double fRMin     [4];
    double fRMax     [4];
    double fZ0       [4];
    double fHexSize;
    double fMinFraction;
    double fWrapperThickness;
    double fShellThickness;
  };

  int      fInitialized;
  int      fNDisks;
  TDisk*   fDisk[2];

  double   fEnergy; // total energy in the calorimeter
//-----------------------------------------------------------------------------
// methods
//-----------------------------------------------------------------------------
  TDiskCalorimeter();
  TDiskCalorimeter(GeomData_t* Geom);

  ~TDiskCalorimeter();
//-----------------------------------------------------------------------------
// accessors
//-----------------------------------------------------------------------------
  int     Initialized() { return fInitialized; }
  int     NDisks()    { return fNDisks; }
  TDisk*  Disk(int I) { return fDisk[I]; }

  int     DiskNumber   (int HitID);
  double  CrystalRadius(int HitID);
  double  Energy()     { return fEnergy; }
//-----------------------------------------------------------------------------
// other methods
//-----------------------------------------------------------------------------
  int Init(GeomData_t* Geom);

  int InitEvent(TCalDataBlock* CalDataBlock);
  
//-----------------------------------------------------------------------------
// overloaded methods of TObject
//-----------------------------------------------------------------------------
  virtual void  Clear(Option_t* Opt = "") override;
  virtual void  Print(Option_t* Opt = "") const override;

  ClassDefOverride(TDiskCalorimeter,0)

};

//}  // end namespace

#endif
