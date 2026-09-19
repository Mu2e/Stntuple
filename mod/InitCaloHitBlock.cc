///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "Stntuple/mod/InitCaloHitBlock.hh"

#include "Offline/RecoDataProducts/inc/CaloHit.hh"

#include <iostream>
#include <format>
#include <vector>

using std::vector ;
//-----------------------------------------------------------------------------
// in this case AbsEvent is just not used
//-----------------------------------------------------------------------------
namespace stntuple {

//-----------------------------------------------------------------------------
InitCaloHitBlock::InitCaloHitBlock() : TStnInitDataBlock() {
  fLastRun     = -1;
}

//-----------------------------------------------------------------------------
int InitCaloHitBlock::InitDataBlock(TStnDataBlock* Block, AbsEvent* Event, int Mode) {
  std::string oname("InitCaloHitBlock::InitDataBlock");

  int ev_number, rn_number; // , mc_flag(0); /*,n_combo_hits(0), n_straw_hits(0)*/

  ev_number = Event->event();
  rn_number = Event->run();

  if (Block->Initialized(ev_number,rn_number)) return 0;

  TCaloHitBlock* chb = (TCaloHitBlock*) Block;
  chb->Clear();
  
  chb->f_RunNumber   = rn_number;
  chb->f_EventNumber = ev_number;

  // if (rn_number < 100000) mc_flag = 1; 

  if (fLastRun != rn_number) {
    fLastRun  = rn_number;
  }
//-----------------------------------------------------------------------------
// assume that straw hits and combo hits are created by the same module - why ?
//-----------------------------------------------------------------------------
  int nhits(0);
  
  if (! fCaloHitCollTag.empty() != 0) {
    art::Handle<mu2e::CaloHitCollection> chch;
    bool ok = Event->getByLabel(fCaloHitCollTag,chch);
    if (ok) { 
      fCaloHitColl = chch.product();
      nhits        = fCaloHitColl->size();
    }
    else {
      mf::LogWarning(oname) << " ERROR:" << __LINE__ 
                            << " : mu2e::CaloHitCollection " 
                            << fCaloHitCollTag.encode().data() 
                            << " not found.";
      return -1;
    }
  }

  if (nhits == 0) return 0;

  //  int ncrdtot(0);
  fCaloRecoDigiColl = nullptr;
  if (! fCaloRecoDigiCollTag.empty() != 0) {
    art::Handle<mu2e::CaloRecoDigiCollection> crdch;
    bool ok = Event->getByLabel(fCaloRecoDigiCollTag,crdch);
    if (ok) { 
      fCaloRecoDigiColl = crdch.product();
      // ncrdtot           = fCaloRecoDigiColl->size();
    }
    else {
      mf::LogWarning(oname) << " WARNING:" << __LINE__ 
                            << " : mu2e::CaloHitCollection " 
                            << fCaloRecoDigiCollTag.encode().data() 
                            << " not found.";
    }
  }

  if (nhits == 0) return 0;
  
  const mu2e::CaloHit* ch0 = &fCaloHitColl->at(0);
  for (int i=0; i<nhits; i++) {
    const mu2e::CaloHit* ch = &fCaloHitColl->at(i);
    int ind = ch-ch0;
    // float corrected_time = ch->correctedTime();
    // if (fMinDt > 0) {
    //   if (ch->energyDep() < fMinEDep) {
    //     bool close_enough = CloseEnough(corrected_time);
    //     if (not close_enough) continue;
    //   }
    // }

    //    TCaloHit* nt_ch = chb->NewHit(i);

    // int pln = ch->strawId().plane();
    // int pnl = ch->strawId().panel();
    // const mu2e::TrkPanelMap::Row* tpmd = fTpm->panel_map_by_offline_ind(pln,pnl);
    
    // ind is a hit index in the original Mu2e hit collection
    TCaloHit* tch   = chb->NewCaloHit(ind);
    tch->fCid       = ch->crystalID();
    int ncrd        = ch->recoCaloDigis().size();
    tch->fNSipms    = (ch->nSiPMs() & 0xff) | ((ncrd & 0xff) << 8) ;
    if (fCaloRecoDigiColl != nullptr) {
      auto crd0 = &fCaloRecoDigiColl->at(0);
      for (int icrd=0; icrd<ncrd; ++icrd) {
        const mu2e::CaloRecoDigi* crd = ch->recoCaloDigis().at(icrd).get();
        int index = crd-crd0;
        if (icrd < 2) {
          tch->fCrdIndex[icrd] = index;
        }
        else {
          mf::LogWarning(oname) << " ERROR:" << __LINE__ 
                                << " : calo hit number " << i 
                                << " is made from " << ncrd << "CaloRecoDidis.";
        }
      }
    }
    tch->fTime      = ch->time();
    tch->fSigT      = ch->timeErr();
    tch->fEDep      = ch->energyDep();
    tch->fSigE      = ch->energyDepErr();

    int disk = 0;
    if (ch->crystalID() >= 674) disk = 1;
    chb->fEDep[disk] += tch->fEDep;
  }
  
  return 0;
}

//-----------------------------------------------------------------------------
int InitCaloHitBlock::ResolveLinks(TStnDataBlock* Block, AbsEvent* Event, int Mode) {
  std::string oname("InitCaloHitBlock::ResolveLinks");

  const int evn = Event->event();
  const int rn  = Event->run();
  const int srn = Event->subRun();

  if (! Block->Initialized(evn,rn,srn)) return -1;
  if (  Block->LinksInitialized()     ) return  0;
  
  TCaloHitBlock* chb = (TCaloHitBlock*) Block;

  // int nch = ch_block->NHits();
  // for (int i=0; i<nch; i++) {
  //   auto ch = ch_block->Hit(i);
  //   // this one has a list of pointers to straw hits
  //   // rely on that the Mu2e collectiosn being vectors of things
  //   const mu2e::CaloHit* o_ch = ch->OfflineCh();
  //   int nsh = o_ch->nStrawHits();
  //   for (int j=0; j<nsh; j++) {
  //     uint16_t ind = o_ch->index(j);
  //     // assuming all straw hits are stored
  //     // to make that work in case not all straw hits are stored, need to store
  //     // offline index of each straw hit - why not ?
  //     // TObject should have space
  //     // same shoudl hold for all objects stored in Stntuple
  //     ch_block->ListOfShIndices()->Add(i,ind);
  //   }
  // }
//-----------------------------------------------------------------------------
// mark links as initialized
//-----------------------------------------------------------------------------
  chb->fLinksInitialized = 1;

  return 0;
}
}
