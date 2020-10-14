#ifndef __stntuple_ana_TrackPar_t__
#define __stntuple_ana_TrackPar_t__

#include "ParBase_t.hh"

namespace stntuple {
class TrackPar_t : public ParBase_t {
public:
  int     fNHPl;
  int     fNEPl;
  int     fNDPl;
  int     fDiskID;			// 0: first disk, 1:second disk
  int     fIDWord[20];

  int     fIDWord_RMC;                  // ID word for mu- --> E+ analysis (RMC rejection)

  double  fMVAOut[20];			// outputs of different MVA classifiers
  double  fProb;
  
  double  fP   ; 			// momentum, corrected to set DPF at 0
  double  fPStOut;                      // MC true mom on exit from ST
  double  fPFront;                      // MC true mom at the tracker front
  float   fDpF ;			// tracker-only resolution (corrected)
  float   fXDpF;                        // delta(Pfront)/sigmaP
  float   fDp0 ;
  float   fDp2 ;
  float   fDpFSt;
  double  fDioWt;			// DIO LO weight
  double  fDioWtRC;			// DIO LO weight with rad corrections
  double  fLumWt;			// luminosity weight
  double  fTotWt;			// total weight
  double  fTotWtRC;			// total weight with rad corrections

  double  fDtZ0;			// delta(T) at z=0
  double  fDtBack;			// delta(T) at z=Z(TT_Back)
  double  fXtZ0;

  double  fEcl;
  double  fEp;
  double  fDrDzCal;			// dr/dz at the calorimeter intersection point
  double  fDtClZ0;
  double  fEDiskP;
  double  fECalP;
  double  fDx;
  double  fDy;
  double  fDz;
  double  fDt;
  double  fDu;		         	// rotated residuals
  double  fDv;
  double  fChi2Tcm;
  double  fChi2XY;
  double  fChi2T;
  double  fPath;
  
  double  fSinTC;			// angle between the cluster and the track
  double  fDrTC;
  double  fSInt;

  double  fLogLHDedm;			// downstram electron (DE) vs downstream muon (DM)

  double  fTMean;			// mean time over the track hits

};
}
#endif
