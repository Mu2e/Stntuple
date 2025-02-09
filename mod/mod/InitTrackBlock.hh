///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
#ifndef __Stntuple_mod_InitTrackBlock__
#define __Stntuple_mod_InitTrackBlock__

#include <string.h>

#include "canvas/Utilities/InputTag.h"

#include "Stntuple/obj/TStnInitDataBlock.hh"
#include "Stntuple/obj/TStnTrackBlock.hh"

#ifndef __CINT__

#include "Offline/TrkReco/inc/DoubletAmbigResolver.hh"
#include "Offline/TrackerGeom/inc/Tracker.hh"
#include "Offline/RecoDataProducts/inc/AlgorithmID.hh"
#include "Offline/RecoDataProducts/inc/ComboHit.hh"
#include "Offline/RecoDataProducts/inc/KalSeed.hh"
#include "Offline/RecoDataProducts/inc/TrkCaloIntersect.hh"
#include "Offline/RecoDataProducts/inc/MVAResult.hh"
#include "Offline/RecoDataProducts/inc/PIDProduct.hh"
#include "Offline/MCDataProducts/inc/StrawDigiMC.hh"
#include "Offline/MCDataProducts/inc/StepPointMC.hh"

#else
namespace mu2e {
  class AlgorithmIDCollection;
  class DoubletAmbigResolver;
  class Tracker;
  class ComboHitCollection;
  class KalSeedCollection;
  class PIDProductCollection;
  class StrawDigiMCCollection;
  class StepPointMCCollection;
  class MVAResultCollection;
  class TrkCaloIntersectCollection;
};
#endif

class StntupleInitTrackBlock : public TStnInitDataBlock {
public:

  struct ZMap_t {
    int    fMap[44][6][2];		// 44 means "by plane"
    double fZ  [88];
  };

  art::InputTag   fAlgorithmIDCollTag;
  art::InputTag   fCaloClusterCollTag;
  art::InputTag   fSsChCollTag;           // single-straw combo hit collection tag
  art::InputTag   fKFFCollTag;            // everything (all tracks) is a KalSeed now
  art::InputTag   fPIDProductCollTag;
  art::InputTag   fStrawDigiMCCollTag;
  art::InputTag   fVdhCollTag;                     // VDH = virtual detector hit
  art::InputTag   fTciCollTag;                     // TCI = track calo    intersection
  art::InputTag   fTcmCollTag;                     // TCM = track cluster match
  art::InputTag   fTrkQualCollTag;

  TString         fTrackTsBlockName;
  art::InputTag   fTrackTsCollTag;

  TString         fTrackHsBlockName;
  int             fVerbose;

  mu2e::AlgorithmIDCollection*             list_of_algs               ;
  const mu2e::KalSeedCollection*           list_of_kffs               ;
  const mu2e::MVAResultCollection*         list_of_trk_qual           ;
  const mu2e::StrawDigiMCCollection*       list_of_mc_straw_hits      ;
  const mu2e::ComboHitCollection*          fSschColl                  ;
  const mu2e::TrkCaloIntersectCollection*  list_of_extrapolated_tracks;
  const mu2e::PIDProductCollection*        list_of_pidp               ;

  const mu2e::Tracker*                     tracker;
  ZMap_t                                   zmap;

  mu2e::DoubletAmbigResolver*              _dar;
//-----------------------------------------------------------------------------
// functions
//-----------------------------------------------------------------------------
public:

  void   SetAlgorithmIDCollTag      (std::string& Tag) { fAlgorithmIDCollTag   = art::InputTag(Tag); }
  void   SetCaloClusterCollTag      (std::string& Tag) { fCaloClusterCollTag   = art::InputTag(Tag); }
  void   SetSsChCollTag             (art::InputTag& Tag) { fSsChCollTag        = Tag; }
  void   SetKFFCollTag              (std::string& Tag) { fKFFCollTag           = art::InputTag(Tag); }
  void   SetPIDProductCollTag       (std::string& Tag) { fPIDProductCollTag    = art::InputTag(Tag); }
  void   SetVdhCollTag              (art::InputTag& Tag) { fVdhCollTag         = Tag; }
  void   SetStrawDigiMCCollTag      (art::InputTag& Tag) { fStrawDigiMCCollTag = Tag; }
  void   SetTciCollTag              (std::string& Tag) { fTciCollTag           = art::InputTag(Tag); }
  void   SetTcmCollTag              (std::string& Tag) { fTcmCollTag           = art::InputTag(Tag); }
  void   SetTrkQualCollTag          (std::string& Tag) { fTrkQualCollTag       = art::InputTag(Tag); }
  void   SetTrackTsCollTag          (std::string& Tag) { fTrackTsCollTag       = art::InputTag(Tag); }

  void   SetDoubletAmbigResolver    (mu2e::DoubletAmbigResolver* Dar) { _dar   = Dar; }

  void   SetTrackTsBlockName        (const char* Name) { fTrackTsBlockName     = Name; }
  void   SetTrackHsBlockName        (const char* Name) { fTrackHsBlockName     = Name; }

  void   InitTrackerZMap(const mu2e::Tracker* Tracker, ZMap_t* Map);
  void   get_station    (const mu2e::Tracker* Tracker, ZMap_t* Map, double Z, int* Plane, int* Offset);
  double s_at_given_z   (const mu2e::KalSeed* KSeed, double Z);
  
  
  virtual int InitDataBlock(TStnDataBlock* Block, AbsEvent* Evt, int Mode);
  virtual int ResolveLinks (TStnDataBlock* Block, AbsEvent* Evt, int Mode);

};

#endif
