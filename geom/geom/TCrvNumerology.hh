///////////////////////////////////////////////////////////////////////////////
// CRV numerology
///////////////////////////////////////////////////////////////////////////////
#ifndef Stntuple_inc_TCrvNumerology_hh
#define Stntuple_inc_TCrvNumerology_hh

#include "TObject.h"
#include "TString.h"
#include "TVector3.h"

class TCrvNumerology : public TObject {
public:

  struct SectorData_t {
    int     fNumber;
    TString fName;
    int     fNModules;
    int     fNLayers;
    int     fNBarsPerLayer;
    int     fFirstIndex;
    float   fBarLength;     //bar length in meters
  };

  struct BarData_t {
    int     fBarIndex;
    int     fSector;
    int     fSectorType;
    TVector3 fBarPos;
    int     fBarU;
    int     fBarV;
    int     fBarK;
    // float   fBarX;
    // float   fBarY;
    // float   fBarZ;
    // float   fWidthDir;
    float   localX() {
      float dir(1);
      if((fSector >= 0) && (fSector <= 4))
	dir = -1.;
      return dir*fBarPos[fBarU];
    }
    float   localY() { return fBarPos[fBarV];}
  };

  enum { kNSectors = 22, kNBars = 5504 } ;

  SectorData_t fSector[kNSectors];
  BarData_t    fBar[kNBars];

  static  TCrvNumerology* fgInstance;

  class  Cleaner {
  public:
    Cleaner ();
    ~Cleaner();
  };
  friend class Cleaner;
//-----------------------------------------------------------------------------
// methods
//-----------------------------------------------------------------------------
  TCrvNumerology();
  ~TCrvNumerology();
//-----------------------------------------------------------------------------
// accessors
//-----------------------------------------------------------------------------
  int         NSectors          () { return kNSectors;                 }
  const char* Name         (int I) { return fSector[I].fName.Data();   }
  int         NModules     (int I) { return fSector[I].fNModules;      }
  int         NBarsPerLayer(int I) { return fSector[I].fNBarsPerLayer; }
  int         FirstIndex   (int I) { return fSector[I].fFirstIndex;    }
  float       BarLength    (int I) { return fSector[I].fBarLength;     }
  float       LocalBarX    (int I) { return fBar[I].localX(); }
  float       LocalBarY    (int I) { return fBar[I].localY(); }
//-----------------------------------------------------------------------------
// other methods
//-----------------------------------------------------------------------------
  int        GetBarInfo(int BarIndex, int& Sector, int& Module, int& Layer, int& Bar);

  static TCrvNumerology*  Instance();
//-----------------------------------------------------------------------------
// overloaded methods of TObject
//-----------------------------------------------------------------------------
  void    Clear(Option_t* Opt = "") ;
  void    Print(Option_t* Opt = "") const ;

  ClassDefOverride(TCrvNumerology,0)

};

//}  // end namespace

#endif
