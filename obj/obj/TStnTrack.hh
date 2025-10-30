//
// Read the tracks added to the event by the track finding and fitting code.
//
// $Id: TStnTrack.hh,v 1.3 2014/09/13 01:16:28 murat Exp $
// $Author: murat $
// $Date: 2014/09/13 01:16:28 $
//
// Contact person Pavel Murat
//
// v9 : use of of the ints to store NHits
//-----------------------------------------------------------------------------
#ifndef murat_inc_TStnTrack_hh
#define murat_inc_TStnTrack_hh

// Mu2e includes.

// C++ includes.
#include <iostream>

#include "TString.h"
#include "TFolder.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "Stntuple/base/TBitset.hh"
// #include "TStnCluster.hh"

					// 'KalRep' is a BaBar class
class TStnCluster;

namespace mu2e {
  class KalSeed ;
  class CaloCluster;
  class TrkToCaloExtrapol;
}

// namespace murat {

class TStnTrack : public TObject {

  enum {
    kNFreeIntsV1   =  9,
    kNFreeFloatsV1 = 10,

    kNFreeIntsV2   =  9,
    kNFreeFloatsV2 =  6,

    kNFreeIntsV3   =  9,
    kNFreeFloatsV3 =  6,

    kNFreeIntsV4   =  7,		// hopefully, correct
    kNFreeFloatsV4 =  6,		// hopefully, correct

    kNFreeIntsV5   =  7,
    kNFreeFloatsV5 =  6,

    kNFreeIntsV6   =  7,
    kNFreeFloatsV6 =  6,

    kNFreeIntsV7   =  7,
    kNFreeFloatsV7 =  5,

    kNFreeIntsV8   =  7,
    kNFreeFloatsV8 =  4,

    kNFreeIntsV9   =  5,
    kNFreeFloatsV9 =  4,

    kNFreeIntsV10  =  5,     //         in V10 2 words per intersection are added
    kNFreeFloatsV10=  4,

    kNFreeIntsV11  =  5,     //         V11: add Dave's TrkQual, I/O doesn't change
    kNFreeFloatsV11=  3,

    kNFreeIntsV12  =  3,     //         V12: add helix and seed indices, I/O doesn't change
    kNFreeFloatsV12=  3,

                             //         V13: add interData for housing the TrkCaloHit info

    kNFreeIntsV15  =  2,     //         V15: add additional track momentum sampling points and MC trajectory
    kNFreeFloatsV15=  3,
    kNFreeFloats2V15= 20,

    kNFreeInts     =  1,     //         V16: add ST/IPA/OPA intersection counts
    kNFreeFloats   =  2,
    kNFreeFloats2  = 20,
  };

  //  enum { kMaxNLayers = 88 }; // 22x2*2
  enum { kMaxNLayers = 144 }; // 18(stations)x2(planes/station)*2(faces/plane)*2(layers/face)

					// added 3 chi*2's
public:
  enum { kNVanesV4 = 4,
	 kNVanesV5 = 4,
	 kNDisks   = 2			// there are only two disks, no vane support
  };

  struct InterData_t {
    int          fID;			// = -1 if no intersection
    int          fClusterIndex;         // cluster index in the list of clusters
    float        fTime;			// extrapolated track time, not corrected by _dtOffset
    float        fEnergy;		// closest cluster energy
    float        fXTrk;
    float        fYTrk;
    float        fZTrk;
    float        fNxTrk;		// track direction cosines in the intersection point
    float        fNyTrk;
    float        fNzTrk;
    float        fXCl;			// cluster coordinates
    float        fYCl;
    float        fZCl;
    float        fDx;			// TRK-CL
    float        fDy;			// TRK-CL
    float        fDz;
    float        fDt;			// TRK-CL , _corrected_ by _dTOffset (!)
    float        fDu;			// ** added in V6
    float        fDv;			// ** added in V6
    float        fChi2Match;		// track-cluster match chi&^2 (coord)
    float        fChi2Time;		// track-cluster match chi&^2 (time)
    float        fPath;			// track path in the disk
    float        fIntDepth;             // ** added in V6 :assumed interaction depth
    float        fDr;                   // ** added in V10: DR(cluster-track), signed
    float        fSInt;                 // ** added in V10: interaction length, calculated
    const mu2e::CaloCluster*       fCluster;
    const mu2e::TrkToCaloExtrapol* fExtrk;
  };
    
  TLorentzVector            fMomentum;        // this assumes DEM fit hypothesis
  
  TBitset                   fHitMask;	      // bit #i: 1 if there is an active hit 
  TBitset                   fExpectedHitMask; // bit #i: 1 if expect to have a hit at this Z

  int                       fNumber;          // track index in the list of reconstructed tracks
  int                       fNHyp;            // number of hyp's with successfull fits
  int                       fBestHyp[2];      // hypothesis with the best chi2/ndof
  int                       fIDWord;	      // now - for selection "C"

  int                       fNActive;	      // NWrong << 16 + NActive
  //  int                       fVaneID;	      // (Nmat ) + (nmatactive << 16) 
  int                       fNMatSites;	      // (Nmat ) + (nmatactive << 16) 
  int                       fDiskID;	      // 
  int                       fPdgCode;         // PDG code of the particle produced most hits
  int                       fNGoodMcHits;     // (NDOF << 16) + (Nhits on the track produced by the associated MC particle)
  int                       fPartID;          // MC particle ID (number in the list)
  int                       fNMcStrawHits;    // Nhits by associated particle in the straw tracker
  int                       fAlgorithmID;     // bit-packed : (alg_mask << 16 ) | best
  int                       fNHits;           // undefined before V9: total number of hits associated with the track | (nbend << 16)
  int                       fNDoublets;       // undefined before V9: nd_os | (nd_ss << 8) | (nhitsambig0 << 16) | (nda << 24)
  int                       fHelixIndex;      // added in V12
  int                       fSeedIndex;
  int                       fMcDirection;     // -1 for upstream, +1 for downstream
  int                       fInterCounts;     // ST/IPA/OPA/TSdA intersection (packed)
  int                       fInt[kNFreeInts]; // provision for future I/O expansion
  
  float                     fChi2;
  float                     fTBack;           // T at Z(TT_Back)
  float                     fFitCons;
  float                     fT0;
  float                     fT0Err;
  float                     fFitMomErr;
  float                     fTanDip;	      // at Z=Z0
  float                     fP;		      // total momentum in the first point (Z=Z1)
  float                     fCharge;
  float                     fPt;	      // transverse momentum in the first point (Z=Z1)
  float                     fD0;	      // at Z=Z0
  float                     fZ0;

  float                     fPStOut;    // MC momentum in the VD ST_Out 
  float                     fPFront;    // MC momentum in the VD front of the tracker

  float                     fClusterE;	// energy of the associated cluster
  float                     fDt;
  float                     fEp;
  float                     fDx;	// about 0 for vanes
  float                     fDy;	// 
  float                     fDz;	// about 0 for disks
  
  float                     fEleLogLHCal;         // log likelihood of the electron hypothesis
  float                     fMuoLogLHCal;	  // log lilelihood of the muon     hypothesis

  float                     fRSlope;	          // timing residual slope dres(T)/dZ
  float                     fRSlopeErr;

  float                     fLogLHRXs;            // XSlope-only-based likelihood

  float                     fEleLogLHDeDx;        // dE/dX LH calculated by Vadim based 
  float                     fMuoLogLHDeDx;        // 
  float                     fX1;	          // X1,Y1,Z1: track coordinates at Z of the first hit
  float                     fY1;
  float                     fZ1;
  float                     fP0;                  // momentum defined at Z=Z(TT_FrontPA) - tracker front
  float                     fP2;                  // momentum defined at Z=Z(TT_BACK)    - tracker back
  float                     fC0;	          // curvature at Z0
  float                     fPhi0;	          // phi0 at Z0 **specify in V8 , no I/O changes***
  float                     fTrkQual;             // ** added in V11
  float                     fTFront;              // time at tracker front, added in V17
  float                     fFloat[kNFreeFloats]; // provision for future I/O expansion

  InterData_t               fDisk [kNDisks];      // track intersections with disks

  // new floats added in v15
  float fPSTFront; // momentum at the front of the stopping target
  float fPSTBack; // at back
  float fPTrackerEntrance; // momentum at tracker entrance (i.e. entering to stopping target)
  float fPTrackerMiddle; // at middle
  float fPTrackerExit; // at exit (i.e. exiting towards calorimeter)
  float fFloat2[kNFreeFloats2]; // for future I/O expansion for versions > 15

  // floats added 
  
  InterData_t               fTrkCaloHit;          // TrkCaloHit info
//-----------------------------------------------------------------------------
//  transient data members, all persistent ones should go above
//-----------------------------------------------------------------------------
  const mu2e::TrkToCaloExtrapol* fExtrk;              //!
  const mu2e::CaloCluster*       fClosestCaloCluster; //!

  InterData_t*                   fVMinS;	      //! intersection with min S
  InterData_t*                   fVMaxEp;	      //! intersection with max E/P
  InterData_t*                   fVTCH;               //! TrackCaloHit info

  mu2e::KalSeed*                 fKalRep[4];          //! different fits, sequence: e-, e+, mu-, mu+

  int                            fITmp[ 5];           //!
  Float_t                        fTmp [10];           //! for temporary analysis needs
  int                            fNHPerStation[50];   //! currently - 44 x 3
  TStnCluster*                   fCluster;            //! 
//-----------------------------------------------------------------------------
// 
//-----------------------------------------------------------------------------
  TStnTrack(int i = -1);
  ~TStnTrack();

  TLorentzVector* Momentum() { return &fMomentum; }

  TBitset*        HitMask        () { return &fHitMask; }
  TBitset*        ExpectedHitMask() { return &fExpectedHitMask; }
  
  int    Number   () const { return fNumber; }
  int    NHits    () const { return (fNHits        ) & 0xffff; }
  int    NBend    () const { return (fNHits   >> 16) & 0xffff; }
  int    NActive  () const { return (fNActive      ) & 0xffff; }
  int    NWrong   () const { return (fNActive >> 16) & 0xffff; }

  int    NOSDoublets  () const { return (fNDoublets      ) & 0xff;   }
  int    NSSDoublets  () const { return (fNDoublets >> 8 ) & 0xff;   }
  int    NDoublets    () const { return NOSDoublets()+NSSDoublets(); }
  int    NHitsAmbZero () const { return (fNDoublets >> 16) & 0xff;   }
  int    NDoubletsAct () const { return (fNDoublets >> 24) & 0xff;   }

  int    NMat         () const { return (fNMatSites      ) & 0xffff; }
  int    NMatActive   () const { return (fNMatSites >> 16) & 0xffff; }

  int    NSTIntersections ()  const { return (fInterCounts      ) & 0xff;}
  int    NIPAIntersections()  const { return (fInterCounts >>  8) & 0xff;}
  int    NOPAIntersections()  const { return (fInterCounts >> 16) & 0xff;}
  int    NTSdAIntersections() const { return (fInterCounts >> 24) & 0xff;}

  int    NClusters();
  int    NMcStrawHits() const { return fNMcStrawHits; }

                                        // P.M. pack the number of degrees of freedom into the same word with
                                        // the number of good MC hits, not sure if the latter is still used

  int    NGoodMcHits () const { return (fNGoodMcHits      ) & 0xffff; }
  int    NDoF        () const { return (fNGoodMcHits >> 16) & 0xffff; }

  int    AlgorithmID() const { return fAlgorithmID; }
  int    BestAlg    () const { return (fAlgorithmID      ) & 0xffff; }
  int    AlgMask    () const { return (fAlgorithmID >> 16) & 0xffff; }

  int    HelixIndex     () const { return fHelixIndex; }
  int    TrackSeedIndex () const { return fSeedIndex; }

  int    McDirection () const { return fMcDirection;  }

  int    PDGCode() const { return fPdgCode; }
                                        // ID of the corresponding TSimParticle
  int    SimID  () const { return fPartID;  }

  float  T0       () const { return fT0       ; }
  float  T0Err    () const { return fT0Err    ; }
  float  FitCons  () const { return fFitCons  ; }
  float  FitMomErr() const { return fFitMomErr; }
  float  TanDip   () const { return fTanDip   ; }
  float  Charge   () const { return fCharge   ; }

  float  D0       () const { return fD0; }
  float  Z0       () const { return fZ0; }
  float  Dt       () const { return fDt; }
  float  Ep       () const { return fEp; }
  float  Dy       () const { return fDy; }
  float  Dz       () const { return fDz; }
                                                          // for KK fit, NDoF is a tricky thing, rely on Dave
  float  Chi2Dof  () const { return fChi2/NDoF(); }
  float  ClusterE () const { return fClusterE;    }

  int    TCH_diskId () const { if (fVTCH) { return fTrkCaloHit.fID;} else {return -1e6;} }
  float  TCH_Ep     () const { if (fVTCH) { return fTrkCaloHit.fEnergy/fP2;} else {return -1.e6;} }
  float  TCH_Dt     () const { if (fVTCH) { return fTrkCaloHit.fDt;} else {return -1.e6;} }
  float  TCH_Doca   () const { if (fVTCH) { return fTrkCaloHit.fDr;} else {return -1.e6;} }
  float  TCH_HitLen () const { if (fVTCH) { return fTrkCaloHit.fPath;} else {return -1.e6;} }

  float  TBack    () const { return fTBack; }
  float  TFront   () const { return fTFront; }

  float  EleLogLHCal() const { return fEleLogLHCal; }
  float  MuoLogLHCal() const { return fMuoLogLHCal; }
  float  LogLHRCal  () const { return fEleLogLHCal-fMuoLogLHCal; }

  float  RSlope   () const { return fRSlope;    }
  float  RSlopeErr() const { return fRSlopeErr; }
  float  XSlope   () const { return fRSlope/fRSlopeErr; }

					// track-only likelihood based on timing residuals

  float  EleLogLHDeDx() const { return fEleLogLHDeDx; }
  float  MuoLogLHDeDx() const { return fMuoLogLHDeDx; }
  float  LogLHRDeDx  () const { return fEleLogLHDeDx-fMuoLogLHDeDx; }
  float  LogLHRXs    () const { return fLogLHRXs; }

  float  LogLHRTrk   () const { return LogLHRXs()+LogLHRDeDx(); }

  TStnCluster*   Cluster() const { return fCluster; }

  float  P () const { return fP; }
  float  Pt() const { return fPt; }
  Int_t  GetMomentum  (TLorentzVector* Momentum) ;
  
  float  PSTFront         () const { return fPSTFront; }
  float  PSTBack          () const { return fPSTBack; }
  float  PTrackerEntrance () const { return fPTrackerEntrance; }
  float  PTrackerMiddle   () const { return fPTrackerMiddle; }
  float  PTrackerExit     () const { return fPTrackerExit; }
  
  mu2e::KalSeed*   GetKalRep() { return fKalRep[0]; }
  
  float  P0         () const { return fP0;      }
  float  P2         () const { return fP2;      }
  float  C0         () const { return fC0;      }
  float  Phi0       () const { return fPhi0;    }
  float  DaveTrkQual() const { return fTrkQual; }  // kept for not hurting previous Ana codes
  float  TrkQual    () const ;

  float  RMax       () const { return fabs(fD0+2/fC0); }
  float  Tmp   (int I) const { return fTmp[I]; }
//-----------------------------------------------------------------------------
// setters
//-----------------------------------------------------------------------------
  void   SetNumber        (int I ) { fNumber      = I ; }
  void   SetAlgorithmID   (int ID) { fAlgorithmID = ID; }
  void   SetHelixIndex    (int I ) { fHelixIndex  = I;  }
  void   SetTrackSeedIndex(int I ) { fSeedIndex   = I;  }

  void   SetITmp          (int I, int   Val) { fITmp[I] = Val; }
  void   SetTmp           (int I, float Val) { fTmp [I] = Val; }
//-----------------------------------------------------------------------------
// overloaded methods of TObject
//-----------------------------------------------------------------------------
  virtual void Clear (Option_t* Opt = "") ;
  virtual void Print (Option_t* Opt = "") const ;

//-----------------------------------------------------------------------------
// schema evolution
//-----------------------------------------------------------------------------
  void ReadV4 (TBuffer& R__b);
  void ReadV5 (TBuffer& R__b);
  void ReadV6 (TBuffer& R__b);
  void ReadV7 (TBuffer& R__b);
  void ReadV8 (TBuffer& R__b);
  void ReadV9 (TBuffer& R__b);
  void ReadV10(TBuffer& R__b);

  ClassDef(TStnTrack,17)

};

#endif
