#ifndef Stntuple_base_TSsqCrystalMap_hh
#define Stntuple_base_TSsqCrystalMap_hh

#include "Stntuple/geom/TDiskIndex.hh"
#include "Stntuple/geom/TDiskCrystalMap.hh"

class TSsqCrystalMap: public TDiskCrystalMap {
protected:
  static int fgStep[12];

public:

  TSsqCrystalMap();
  TSsqCrystalMap(double Size, double RMin, double RMax); 

  virtual ~TSsqCrystalMap() override;

//-----------------------------------------------------------------------------
// virtual functions of TDiskCrystalMap
//-----------------------------------------------------------------------------
  virtual int  GetFirst(int Ir) const override { 
    if (Ir == 0) return 0;
    else         return 3*Ir*(Ir-1)+1;
  }
					// total number of crystals per ring (including ones outside the disk)

  virtual int GetNCrystalsPerRing(int I) const override { 
    if (I == 0) return 1;
    else        return 6*I; 
  }

  virtual int    GetNTotal() const override { return 3*fNRings*(fNRings-1)+1 ; }

  virtual int    GetRing(int I) override;

  virtual int    GetRing(TDiskIndex* Index) override;

  virtual void   GetPosition(int I, TVector2* Pos) override;
  virtual void   GetPosition(TDiskIndex* Index, TVector2* Pos) override;

  virtual double GetRadius(TDiskIndex* Index) override;
  virtual double GetRadius(int I) override;

  virtual int    InsideCode(TDiskIndex* Index, double* Fraction) override;
};

#endif
