#ifndef Stntuple_base_TSqrCrystalMap_hh
#define Stntuple_base_TSqrCrystalMap_hh

#include "Stntuple/geom/TDiskIndex.hh"
#include "Stntuple/geom/TDiskCrystalMap.hh"

class TSqrCrystalMap: public TDiskCrystalMap {
protected:

public:

  TSqrCrystalMap();
  TSqrCrystalMap(double Size, double RMin, double RMax); 

  virtual ~TSqrCrystalMap() override;

//-----------------------------------------------------------------------------
// virtual functions of TDiskCrystalMap
//-----------------------------------------------------------------------------
  virtual int  GetFirst(int Ir) const override { 
    if (Ir == 0) return 0;
    else         {
      printf(" TSqrCrystalMap::GetFirst ERROR: not implemented yet\n");
      return -1;
    }
  }
					// total number of crystals per ring (including ones outside the disk)

  virtual int GetNCrystalsPerRing(int I) const override { 
    if (I == 0) return 1;
    else        return 4*(I+1); 
  }

  virtual int    GetNTotal() const override { return 1+4*(fNRings*(fNRings+1)/2-1) ; }

  virtual int    GetRing(int I) override;
  virtual int    GetRing(TDiskIndex* Index) override;

  virtual void   GetPosition(int I, TVector2* Pos) override;
  virtual void   GetPosition(TDiskIndex* Index, TVector2* Pos) override;

  virtual double GetRadius(TDiskIndex* Index) override;
  virtual double GetRadius(int I) override;

  virtual int    InsideCode(TDiskIndex* Index, double* Fraction) override;
};

#endif
