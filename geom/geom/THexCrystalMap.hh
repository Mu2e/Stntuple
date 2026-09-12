#ifndef Stntuple_base_THexCrystalMap_hh
#define Stntuple_base_THexCrystalMap_hh

#include "Stntuple/geom/TDiskIndex.hh"
#include "Stntuple/geom/TDiskCrystalMap.hh"

class THexCrystalMap : public TDiskCrystalMap {
protected:
					// hexagon vertices
  static TDiskIndex fgPos[6];

public:

  THexCrystalMap(double Size, double RMin, double RMax);
  THexCrystalMap(); 

  virtual ~THexCrystalMap() override;

  virtual int  GetFirst(int Ir) const override{ 
    if (Ir == 0) return 0;
    else         {
      printf(" THexCrystalMap::GetFirst ERROR: not implemented yet\n");
      return -1;
    }
  }
					// total number of crystals per ring (including ones outside the disk)

  virtual int GetNCrystalsPerRing(int I) const override { 
    if (I == 0) return 1;
    else        return 6*I; 
  }

  virtual int    GetNTotal() const override { return 3*fNRings*(fNRings-1)+1 ; }

  virtual int GetRing(int I) override;
  virtual int GetRing(TDiskIndex* Index) override;

  virtual void GetPosition(int I, TVector2* Pos) override;
  virtual void GetPosition(TDiskIndex* Index, TVector2* Pos) override;

  virtual double GetRadius(TDiskIndex* Index) override;
  virtual double GetRadius(int I) override;

  virtual TDiskIndex&    Pos(int I) { return fgPos[I]; }

  virtual int InsideCode(TDiskIndex* Index, double* Fraction) override;

};

#endif
