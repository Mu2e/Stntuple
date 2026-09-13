//-----------------------------------------------------------------------------
//  Apr 2016 G. Pezzullo: initialization of the MU2E STNTUPLE TimePeak block
//  2020-10-18 P.M. rewrite 
//-----------------------------------------------------------------------------
#include <cstdio>
#include "TROOT.h"
#include "TFolder.h"
#include "TLorentzVector.h"
#include "TVector2.h"

#include "messagefacility/MessageLogger/MessageLogger.h"

#include "Stntuple/obj/TStnDataBlock.hh"
#include "Stntuple/obj/TStnEvent.hh"

#include "Stntuple/obj/TStnTimeCluster.hh"
#include "Stntuple/obj/TStnTimeClusterBlock.hh"

#include "Stntuple/obj/TStnHelix.hh"
#include "Stntuple/obj/TStnHelixBlock.hh"

#include "Stntuple/mod/InitTimeClusterBlock.hh"

#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Event.h"

#include "Offline/GeometryService/inc/GeometryService.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"

#include "Offline/CalorimeterGeom/inc/Calorimeter.hh"

#include "Offline/RecoDataProducts/inc/TimeCluster.hh"
#include "Offline/RecoDataProducts/inc/HelixSeed.hh"
#include "Offline/RecoDataProducts/inc/StrawHitIndex.hh"

#include "Offline/MCDataProducts/inc/StrawDigiMC.hh"
#include "Offline/MCDataProducts/inc/StrawGasStep.hh"
#include "Offline/MCDataProducts/inc/SimParticle.hh"

#include "Offline/RecoDataProducts/inc/CaloCluster.hh"

namespace stntuple {
//-----------------------------------------------------------------------------
// assume that the collection name is set, so we could grab it from the event
// ComboHitCollection and StrawHitCollection are of the same time, the first one 
// contains real combo hits (one combo hit could be made out of more than one straw hit),
// the other one has a combo hit per straw digi
//-----------------------------------------------------------------------------
int  InitTimeClusterBlock::InitDataBlock(TStnDataBlock* Block, AbsEvent* Event, int Mode) {
  const char* oname = {"StntupleInitTimeClusterBlock::InitDataBlock"};

  int ev_number, rn_number; // , mc_flag(0); /*,n_combo_hits(0), n_straw_hits(0)*/

  ev_number = Event->event();
  rn_number = Event->run();
  //  if (rn_number < 100000) mc_flag = 1;

  TStnTimeClusterBlock* tcb = (TStnTimeClusterBlock*) Block;
  tcb->Clear();

  const mu2e::TimeClusterCollection* tcc(nullptr);
  int                                ntc(0);

  if (! fTcCollTag.empty()) {
    art::Handle<mu2e::TimeClusterCollection> tccH;
    bool ok = Event->getByLabel(fTcCollTag,tccH);
    if (ok) {
      tcc = tccH.product();
      ntc               = tcc->size();
    }
  }

  art::Handle<mu2e::ComboHitCollection>    chcH;
  const mu2e::ComboHitCollection*          chc(nullptr);
//-----------------------------------------------------------------------------
// combohits in ntuples are mostly needed for debugging and MC-specific purpose
//-----------------------------------------------------------------------------
  if (! fChCollTag.empty()) {
    bool ok = Event->getByLabel(fChCollTag,chcH);
    if (ok) {
      chc          = chcH.product();
    }
    else {
      mf::LogWarning(oname) << " ERROR: no ComboHitCollection tag=" 
			    << fChCollTag.encode().data() <<  " found. BAIL OUT";
      // return -1;
    }
  }
//-----------------------------------------------------------------------------
// single straw hit collection (also ComboHit's
//-----------------------------------------------------------------------------
  // if (! fShCollTag.empty()) {
  //   bool ok = Event->getByLabel(fShCollTag,sschcH);
  //   if (ok) {
  //     sschc          = sschcH.product();
  //   }
  // }

  art::Handle<mu2e::StrawDigiMCCollection> sdmccH;
  const mu2e::StrawDigiMCCollection*       mcdigis(nullptr);

  if (! fSdmcCollTag.empty()) {
    bool ok = Event->getByLabel(fSdmcCollTag,sdmccH);
    if (ok) {
      mcdigis = sdmccH.product();
    }
  }

  const mu2e::CaloCluster     *cluster(0);
  
  for (int i=0; i<ntc; i++) {
    TStnTimeCluster* tc = tcb->NewTimeCluster();
    const mu2e::TimeCluster* otc = &tcc->at(i);
    cluster             = otc->caloCluster().get();
    if (cluster != 0) {
      mu2e::GeomHandle<mu2e::Calorimeter> ch;
      const mu2e::Calorimeter* _calorimeter = ch.get();      
      
      tc->fClusterTime    = cluster->time();
      tc->fClusterEnergy  = cluster->energyDep();
      CLHEP::Hep3Vector         gpos = _calorimeter->diskToMu2e(cluster->diskID(),cluster->cog3Vector());
      CLHEP::Hep3Vector         tpos = _calorimeter->mu2eToTracker(gpos);
      tc->fClusterX       = tpos.x();
      tc->fClusterY       = tpos.y();
      tc->fClusterZ       = tpos.z();
    }

    tc->fOfflineTc    = otc;
    tc->fNComboHits   = otc->hits().size();
    tc->fNHits        = otc->nStrawHits();
    tc->fT0           = otc->t0()._t0;
    tc->fT0Err        = otc->t0()._t0err;     
    tc->fPosX         = otc->position().x();     
    tc->fPosY         = otc->position().y();     
    tc->fPosZ         = otc->position().z();
//-----------------------------------------------------------------------------
// loop over combo hits to determine the matching MC particle
//-----------------------------------------------------------------------------
    const mu2e::StrawGasStep* step (nullptr);

    int const max_np(200);
    int     np(0), sim_nh[max_np], sim_id[max_np], pdg_code[max_np]; 

    const mu2e::SimParticle*  simm[max_np];
    simm[0] = nullptr;

    if (mcdigis and chc) {
      for (int ih=0; ih<tc->fNComboHits; ih++) {
        StrawHitIndex hit_index   = otc->hits().at(ih);
        const mu2e::ComboHit* hit = &chc->at(hit_index);
//-----------------------------------------------------------------------------
// loop over straw hits of one combo hit
//-----------------------------------------------------------------------------
        int nsh = hit->nStrawHits();
        for (int ish=0; ish<nsh; ish++) {
          int ind = hit->index(ish);
          const mu2e::StrawDigiMC* mcdigi = &mcdigis->at(ind);
          
          step = mcdigi->earlyStrawGasStep().get();

          int id(-1);
          const mu2e::SimParticle* sim(nullptr);
          if (step) {
            sim = step->simParticle().get(); 
            id  = sim->id().asInt();
          }
//-----------------------------------------------------------------------------
// accumulate list of sim particle ID's for this time cluster
//-----------------------------------------------------------------------------
          int found = 0;
          for (int ip=0; ip<np; ip++) {
            if (id == sim_id[ip]) {
              found       = 1;
              sim_nh[ip] += 1;
              break;
            }
          }
	
          if (sim && (found == 0)) {
            sim_id  [np] = id;
            simm    [np] = sim;
            pdg_code[np] = sim->pdgId();
            sim_nh  [np] = 1;
            np          += 1;
          }
        }
      }
//-----------------------------------------------------------------------------
// identify time cluster with the particle which produced most hits
//-----------------------------------------------------------------------------
      int ipart  = -1;
      int max_nh = sim_nh[0];
    
      const mu2e::SimParticle* best_sim(nullptr);

      for (int ip=1; ip<np; ip++) {
        if (sim_nh[ip] > max_nh) {
          max_nh = sim_nh[ip];
          ipart  = ip;
        }
      }
    
      if (ipart >= 0) {
        best_sim        = simm    [ipart];
        tc->fSimID      = sim_id  [ipart];
        tc->fPdgID      = pdg_code[ipart];
        tc->fNHitsSimID = max_nh;
        tc->fMcMom      = best_sim->startMomentum().mag() ;
      }
    }
  }

  tcb->f_RunNumber   = rn_number;
  tcb->f_EventNumber = ev_number;

  return 0;
}

//-----------------------------------------------------------------------------
Int_t InitTimeClusterBlock::ResolveLinks(TStnDataBlock* Block, AbsEvent* AnEvent, int Mode) {
  // Mu2e version, do nothing

  Int_t  ev_number, rn_number;

  ev_number = AnEvent->event();
  rn_number = AnEvent->run();

  if (! Block->Initialized(ev_number,rn_number)) return -1;

					// do not do initialize links 2nd time

  if (Block->LinksInitialized()) return 0;
//-----------------------------------------------------------------------------
// determine the helix corresponding to this time cluster
// what if there is more than one ? and what if several helix finders used
// the same list of time clusters as input ?
// - this makes having only one helix index a very questionable proposition
// - disable that
// instead build a list of links to the combohits
//-----------------------------------------------------------------------------
//  TStnHelix*                 helixseed;

  // const mu2e::TimeCluster*   ktcluster, *fktcluster;
  // const mu2e::HelixSeed*     kseed;

  // char                       short_tc_block_name[100];

  //  TStnEvent* ev = Block->GetEvent();
  auto tc_block = (TStnTimeClusterBlock*) Block;

  int ntc = tc_block->NTimeClusters();

  for (int i=0; i<ntc; i++) {
    TStnTimeCluster* tc = tc_block->TimeCluster(i);
    const mu2e::TimeCluster* otc = tc->OfflineTc();
    int nch = otc->hits().size();
    for (int j=0; j<nch; j++) {
      uint16_t ind = otc->hits().at(j);
      tc_block->ListOfChLinks()->Add(i,ind);
    }
  }
  
//-----------------------------------------------------------------------------
// mark links as initialized
//-----------------------------------------------------------------------------
  tc_block->fLinksInitialized = 1;

  return 0;
}

}
