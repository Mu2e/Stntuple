//-----------------------------------------------------------------------------
//  Apr 2016 G. Pezzullo: initialization of the MU2E STNTUPLE TimePeak block
//  2020-10-18 P.M. rewrite 
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

#include "Offline/RecoDataProducts/inc/CrvDigi.hh"
#include "Stntuple/mod/InitCrvDigiBlock.hh"
//-----------------------------------------------------------------------------
// assume that the collection name is set, so we could grab it from the event
// ComboHitCollection and StrawHitCollection are of the same time, the first one 
// contains real combo hits (one combo hit could be made out of more than one straw hit),
// the other one has a combo hit per straw digi
//-----------------------------------------------------------------------------
int  StntupleInitCrvDigiBlock::InitDataBlock(TStnDataBlock* Block, AbsEvent* Evt, int Mode) {
  const char* oname = {"StntupleInitCrvDigiBlock::InitDataBlock"};

  TCrvDigiBlock* block = (TCrvDigiBlock*) Block;
  
  block->Clear();
  
  const mu2e::CrvDigiCollection*        cdc(nullptr);
  
  art::Handle<mu2e::CrvDigiCollection>  cdch;
  int                                   ndigis(0);
  
  if (! fCrvDigiCollTag.empty()) {
    bool ok = Evt->getByLabel(fCrvDigiCollTag,cdch);
    if (ok) {
      cdc    = cdch.product();
      ndigis = cdc->size();
    }
    else {
      // no crv digi collection: print diagnostics but do nothing else, just leave the data block empty
      mf::LogWarning(oname) << std::format("ERROR: no CrvDigiCollection tag={} found. BAIL OUT",
                                           fCrvDigiCollTag.encode().data());
      return 0;
    }
  }
  
  for (int i=0; i<ndigis; i++) {
    const mu2e::CrvDigi* crvd = &cdc->at(i);
    int ns = crvd->GetADCs().size();

    TCrvDigi* tcrvd = (TCrvDigi*) block->NewDigi();
    tcrvd->Init(ns);
    
    tcrvd->fSbid        = crvd->GetScintillatorBarIndex().asInt();
    tcrvd->fTdc         = crvd->GetStartTDC();
    tcrvd->fNzs         = crvd->IsNZS();
    tcrvd->fOddTs       = crvd->HasOddTimestamp();
    tcrvd->fSipm        = crvd->GetSiPMNumber();
    tcrvd->fRoc         = crvd->GetROC();
    tcrvd->fFeb         = crvd->GetFEB();
    tcrvd->fCh          = crvd->GetFEBchannel();
//-----------------------------------------------------------------------------
// store the waveform, at this point the ADC vector of the TCrvDigi is already resized
//-----------------------------------------------------------------------------
    for (int is=0; is<ns; is++) {
      tcrvd->fAdc[is] = crvd->GetADCs()[is];
    }
  }

  return 0;
}

//_____________________________________________________________________________
Int_t StntupleInitCrvDigiBlock::ResolveLinks(TStnDataBlock* Block, AbsEvent* AnEvent, int Mode) {
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

