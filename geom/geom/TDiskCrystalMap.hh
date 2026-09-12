#ifndef Stntuple_base_TDiskCrystalMap_hh
#define Stntuple_base_TDiskCrystalMap_hh

#include "Stntuple/base/TStnArrayI.hh"
#include "Stntuple/geom/TDiskIndex.hh"

#include "TClonesArray.h"

class TVector2;

class TDiskCrystalMap {
protected:

  double fSize;
  double fRMin;
  double fRMax;
  int    fNRings;
  int    fNCrystals;

  double fMinFraction;

  TClonesArray*  fIndex;                // a list of disk indices, NULL means that the crystal 
					// in not in the ring

  TStnArrayI*    fIndexMap;		// fIndexMap[i] gives 1-D index of the crystal in [sparse] fIndex
					// fIndexMap is non-sparse, runs from 0 to N-1, where N is the total 
					// number of crystals _inside_ the disk
  int*           fNInside ;

public:

  TDiskCrystalMap();
  TDiskCrystalMap(double Size, double RMin, double RMax);

  virtual ~TDiskCrystalMap();

  virtual int GetRing(int I) = 0;
  virtual int GetRing(TDiskIndex* Index) = 0;

  virtual int GetNRings() const { return fNRings; }

					// total number of crystals per ring (including ones outside the disk)

  virtual int GetNCrystalsPerRing(int I) const = 0;

  virtual int GetNTotal() const = 0;
					// number of crystals in a Ring included into a disk

  int GetNInside(int Ring) const { return fNInside[Ring]; }

  int GetNInsideTotal() const { return fIndexMap->NDataWords(); }

  int GetIndexMap(int I) const { return fIndexMap->At(I); }
					// index of the first crystal in a sparse array
  virtual int  GetFirst(int Ir) const = 0;

  virtual TDiskIndex* GetDiskIndex(int I) const { return (TDiskIndex*) fIndex->At(I) ; }

  virtual void GetPosition(int I, TVector2* Pos) = 0;
  virtual void GetPosition(TDiskIndex* Index, TVector2* Pos) = 0;

  virtual double GetRadius(TDiskIndex* Index) = 0;
  virtual double GetRadius(int I) = 0;

  virtual int    InsideCode(TDiskIndex* Index, double* Fraction) = 0;

  virtual void   Print(Option_t* Opt = "") const ;
};

#endif
