#ifndef STNTUPLE_TCrvCoincidenceCluster
#define STNTUPLE_TCrvCoincidenceCluster

#include "TClonesArray.h"
#include "TArrayI.h"
#include "TVector3.h"

class TCrvCoincidenceCluster: public TObject {
public:
  enum {
    kNFreeInts = 10,
    kNFreeFloats = 9
  };

  int          fIndex;                   // index in the list
  int          fSectorType;
  int          fNPulses;
  int          fNPe;
  int          fSimID; //most likely SIM particle
  int          fMCNPulses;
  int          fFreeInts[kNFreeInts];
  float        fStartTime;
  float        fEndTime;
  float        fMCEnergyDep;
  float        fMCAvgHitTime;
  float        fSlope;
  float        fFreeFloats[kNFreeFloats];
  TVector3     fPosition;
  TVector3     fMCAvgPosition;

//-----------------------------------------------------------------------------
//  functions
//-----------------------------------------------------------------------------
public:
					// ****** constructors and destructor
  TCrvCoincidenceCluster();
  virtual ~TCrvCoincidenceCluster();
					// ****** accessors

  int     Index           () const { return fIndex;          }
  int     SectorType      () const { return fSectorType;      }
  int     NPe             () const { return fNPe;             }
  int     NPulses         () const { return fNPulses;         }

  float   StartTime       () const { return fStartTime;       }
  float   EndTime         () const { return fEndTime;         }
  float   Slope           () const { return fSlope;           }

  const TVector3* Position() const { return &fPosition;       }
//-----------------------------------------------------------------------------
// modifiers
//-----------------------------------------------------------------------------
  void Set(int Index, int SectorType, int Np, int NPe, 
	   float X, float Y, float Z, float T1, float T2, float Slope);
  void SetMC(int SimID, int Np, float EnergyDep, float AvgTime,
	   float X, float Y, float Z);
//-----------------------------------------------------------------------------
// overloaded methods of TObject
//-----------------------------------------------------------------------------
  void Clear(Option_t* Opt = "");
  void Print(Option_t* Opt = "") const;

  ClassDef(TCrvCoincidenceCluster,3)	         // STNTUPLE representation of CrvCoincidenceCluster
};


#endif
