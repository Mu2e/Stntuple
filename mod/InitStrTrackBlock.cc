//-----------------------------------------------------------------------------
//  2026-09-09 PM 
//-----------------------------------------------------------------------------
#include <cstdio>
#include <format>

#include "TROOT.h"
#include "TFolder.h"
#include "TLorentzVector.h"
#include "TVector2.h"

#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Event.h"

#include "Stntuple/obj/TStnDataBlock.hh"
#include "Stntuple/obj/TStnEvent.hh"

#include "Offline/RecoDataProducts/inc/KalSeed.hh"
#include "Offline/RecoDataProducts/inc/KalSegment.hh"
#include "Stntuple/mod/InitStrTrackBlock.hh"
//-----------------------------------------------------------------------------
// assume that the collection name is set, so we could grab it from the event
// ComboHitCollection and StrawHitCollection are of the same time, the first one 
// contains real combo hits (one combo hit could be made out of more than one straw hit),
// the other one has a combo hit per straw digi
//-----------------------------------------------------------------------------
namespace stntuple {

//-----------------------------------------------------------------------------
int  InitStrTrackBlock::InitDataBlock(TStnDataBlock* Block, AbsEvent* Evt, int Mode) {
  const char* oname = {"StntupleInitStrTrackBlock::InitDataBlock"};

  TStrTrackBlock* tb = (TStrTrackBlock*) Block;
  
  tb->Clear();
  
  const mu2e::KalSeedCollection*  ksc(nullptr);
  int                             ntracks(0);
  
  if (! fTrackCollTag.empty()) {
    art::Handle<mu2e::KalSeedCollection> ksch;
    bool ok = Evt->getByLabel(fTrackCollTag,ksch);
    if (ok) {
      ksc     = ksch.product();
      ntracks = ksc->size();
    }
    else {
      // no cal digi collection: print diagnostics but do nothing else, just leave the data block empty
      mf::LogWarning(oname) << std::format("ERROR: no KalSeedCollection tag={} found. BAIL OUT",
                                           fTrackCollTag.encode().data());
      return 0;
    }
  }

  if (ntracks == 0) return 0;
  
  const mu2e::KalSeed* ks0 =  &ksc->at(0);
  for (int i=0; i<ntracks; i++) {
    const mu2e::KalSeed* ks = &ksc->at(i);
    // use the second one - seems to be the right one...
    const mu2e::KalSegment* seg = &ks->segments()[1];
    int ind         = ks-ks0;
    TStrTrack* tst  = tb->NewTrack(ind);
    tst->fNHits     = ks->nHits();
    tst->fChi2      = ks->chisquared();
    tst->fT0        = ks->t0().t0();
    tst->fNDof      = ks->nDOF();
    tst->fX0        = seg->position3().x();
    tst->fY0        = seg->position3().y();
    tst->fZ0        = seg->position3().z();
    double tmom     = seg->mom();
    tst->fNx        = seg->momentum3().x()/tmom;
    tst->fNy        = seg->momentum3().y()/tmom;
    tst->fNz        = seg->momentum3().z()/tmom;
    tst->fOfflineKs = ks;
  }

  tb->f_RunNumber    = Evt->run();
  tb->f_EventNumber  = Evt->event();
  tb->f_SubrunNumber = Evt->subRun();

  return 0;
}

//_____________________________________________________________________________
Int_t InitStrTrackBlock::ResolveLinks(TStnDataBlock* Block, AbsEvent* AnEvent, int Mode) {
  // Mu2e version, do nothing

//   Int_t  ev_number, rn_number;

//   ev_number = AnEvent->event();
//   rn_number = AnEvent->run();

//   if (! Block->Initialized(ev_number,rn_number)) return -1;

// 					// do not do initialize links 2nd time

//   if (Block->LinksInitialized()) return 0;

//   TStnEvent*                 ev;
//   TStnTimeClusterBlock*      hb;

//   TStnHelixBlock*            tsb;
//   TStnHelix*                 helixseed;

//   const mu2e::TimeCluster*   ktcluster, *fktcluster;
//   const mu2e::HelixSeed*     kseed;

//   char                       short_tcluster_block_name[100];

//   ev     = Block->GetEvent();
//   hb     = (TStnTimeClusterBlock*) Block;
  
//   hb->GetModuleLabel("mu2e::TimeClusterCollection"  , short_tcluster_block_name);

//   tsb    = (TStnHelixBlock*) ev->GetDataBlock(short_tcluster_block_name);
  
//   int    ntc(0);
//   if (hb!=nullptr){
//     ntc = hb ->NTimeClusters();
//   }
//   int    nhelixseed(0);
//   if (tsb !=nullptr){
//     nhelixseed = tsb->NHelices();
//   }

//   for (int i=0; i<ntc; ++i){
//     TStnTimeCluster* tc = hb->TimeCluster(i);
//     ktcluster = tc->fTimeCluster;
//     int      helixseedIndex(-1);
//     for (int j=0; j<nhelixseed; ++j){
//       helixseed   = tsb->Helix(j);
//       kseed       = helixseed->fHelix;
//       fktcluster  = kseed->timeCluster().get();
//       if (fktcluster == ktcluster) {
// 	helixseedIndex = j;
// 	break;
//       }
//     }
    
//     if (helixseedIndex < 0) {
//       printf(">>> ERROR: TimeClusterFinder timeCluster %i -> no HelixSeed associated\n", i);//FIXME!
// 	  continue;
//     }
    
//     tc->SetHelixSeedIndex(helixseedIndex);
//   }
// //-----------------------------------------------------------------------------
// // mark links as initialized
// //-----------------------------------------------------------------------------
//   hb->fLinksInitialized = 1;

  return 0;
}

}
