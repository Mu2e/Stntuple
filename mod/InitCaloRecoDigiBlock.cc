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

#include "Offline/RecoDataProducts/inc/CaloRecoDigi.hh"
#include "Stntuple/mod/InitCaloRecoDigiBlock.hh"
//-----------------------------------------------------------------------------
// assume that the collection name is set, so we could grab it from the event
// ComboHitCollection and StrawHitCollection are of the same time, the first one 
// contains real combo hits (one combo hit could be made out of more than one straw hit),
// the other one has a combo hit per straw digi
//-----------------------------------------------------------------------------
namespace stntuple {

//-----------------------------------------------------------------------------
int  InitCaloRecoDigiBlock::InitDataBlock(TStnDataBlock* Block, AbsEvent* Evt, int Mode) {
  const char* oname = {"StntupleInitCaloRecoDigiBlock::InitDataBlock"};

  int evn = Evt->event();
  int run = Evt->run();
  int srn = Evt->subRun();
  
  if (Block->Initialized(evn,run)) return 0;

  TCaloRecoDigiBlock* tcrdb = (TCaloRecoDigiBlock*) Block;
  tcrdb->Clear();
  
  tcrdb->f_RunNumber    = run;
  tcrdb->f_EventNumber  = evn;
  tcrdb->f_SubrunNumber = srn;
  
  const mu2e::CaloRecoDigiCollection*  crdc(nullptr);

  int ndigis(0); // tcrdb->fNDigis is incremented in the constructor
  if (! fCaloRecoDigiCollTag.empty()) {
    art::Handle<mu2e::CaloRecoDigiCollection> crdch;
    bool ok = Evt->getByLabel(fCaloRecoDigiCollTag,crdch);
    if (ok) {
      crdc   = crdch.product();
      ndigis = crdc->size();
    }
    else {
      // no cal digi collection: print diagnostics but do nothing else, just leave the data block empty
      mf::LogWarning(oname) << std::format("ERROR: no CaloDigiCollection tag={} found. BAIL OUT",
                                           fCaloRecoDigiCollTag.encode().data());
      return 0;
    }
  }

  if (ndigis == 0)  return 0;
  
  const mu2e::CaloRecoDigi* crd0 =  &crdc->at(0);
  for (int i=0; i<ndigis; i++) {
    const mu2e::CaloRecoDigi* crd = &crdc->at(i);
    // index in the original list of CaloRecoDigis
    int ind = crd-crd0;
    TCaloRecoDigi* tcrd  = tcrdb->NewCaloRecoDigi(ind);
    tcrd->fSipmID   = crd->SiPMID();
    tcrd->fNdf      = crd->ndf();
    tcrd->fPileup   = crd->pileUp() ? 1 : 0;
    tcrd->fTime     = crd->time();
    tcrd->fSigT     = crd->timeErr();
    tcrd->fEDep     = crd->energyDep();
    tcrd->fSigE     = crd->energyDepErr();
    tcrd->fChi2     = crd->chi2();
  }
  // at this point tcrdb->fNDigis is defined
  return 0;
}

//_____________________________________________________________________________
Int_t InitCaloRecoDigiBlock::ResolveLinks(TStnDataBlock* Block, AbsEvent* AnEvent, int Mode) {
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
