///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "Stntuple/mod/InitComboHitBlock.hh"

#include "Offline/RecoDataProducts/inc/ComboHit.hh"

#include "Offline/MCDataProducts/inc/StrawDigiMC.hh"
#include "Offline/MCDataProducts/inc/StrawGasStep.hh"

#include <vector>

using std::vector ;
//-----------------------------------------------------------------------------
// in this case AbsEvent is just not used
//-----------------------------------------------------------------------------
namespace stntuple {

  //-----------------------------------------------------------------------------
InitComboHitBlock::InitComboHitBlock() : TStnInitDataBlock() {
  fSaveShLinks = 0;
  fMinDt       = -1;     // if < 0, save all hits
  fMinEDep     = 0.5e-3; // 0.5 keV
  fLastRun     = -1;
}
  
//-----------------------------------------------------------------------------
// compressing ntuple for specialized studies
// check whether the hist time is close enough to one of the time clusters
//-----------------------------------------------------------------------------
bool InitComboHitBlock::CloseEnough(float Time) {
  bool close_enough = false;
  for (int itc=0; itc<fNTc; itc++) {
    const mu2e::TimeCluster* tc = &fTcc->at(itc);
    float dt = Time - tc->t0().t0();
    if (fabs(dt) < fMinDt) {
      close_enough = true;
      break;
    }
  }
  return close_enough;
}

//-----------------------------------------------------------------------------
int InitComboHitBlock::InitDataBlock(TStnDataBlock* Block, AbsEvent* Event, int Mode) {
  std::string oname("InitComboHitBlock::InitDataBlock");

  int ev_number, rn_number, mc_flag(0); /*,n_combo_hits(0), n_straw_hits(0)*/

  ev_number = Event->event();
  rn_number = Event->run();
  if (rn_number < 100000) mc_flag = 1; 

  if (Block->Initialized(ev_number,rn_number)) return 0;

  TComboHitBlock* data = (TComboHitBlock*) Block;
  data->Clear();
  
  if (fLastRun != rn_number) {
    fTpm      = &fTpmh.get(Event->id());
    fLastRun  = rn_number;
  }
    
//-----------------------------------------------------------------------------
// straw hit information
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
// assume that straw hits and combo hits are created by the same module - why ?
//-----------------------------------------------------------------------------
  if (! fShCollTag.empty() != 0) {
    art::Handle<mu2e::StrawHitCollection> shch;
    bool ok = Event->getByLabel(fShCollTag,shch);
    if (ok) { 
      fShc = shch.product();
      fNSh = fShc->size();
    }
    else {
      fNSh = 0;
      mf::LogWarning(oname) << " ERROR:" << __LINE__ 
                            << " : mu2e::StrawHitCollection " 
                            << fShCollTag.encode().data() 
                            << " not found. No CH->SH links will be stored ";
    }
  }

  if (! fChCollTag.empty() != 0) {
    art::Handle<mu2e::ComboHitCollection> chch;
    bool ok = Event->getByLabel(fChCollTag,chch);
    if (ok) { 
      fChc = chch.product();
      fNCh = fChc->size();
    }
    else {
      fNCh = 0;
      mf::LogWarning(oname) << " ERROR:" << __LINE__ 
                            << " : mu2e::ComboHitCollection " 
                            << fChCollTag.encode().data() 
                            << " not found, BAIL OUT. rc = -1";
      return -1;
    }
  }

  if (! fTcCollTag.empty() != 0) {
    art::Handle<mu2e::TimeClusterCollection> tcch;
    bool ok = Event->getByLabel(fTcCollTag,tcch);
    if (ok) { 
      fTcc = tcch.product();
      fNTc = fTcc->size();
    }
    else {
      fNTc = 0;
      mf::LogWarning(oname) << " ERROR:" << __LINE__ 
                            << " : mu2e::ComboHitCollection " 
                            << fChCollTag.encode().data() 
                            << " not found, BAIL OUT. rc = -1";
    }
  }

  const mu2e::ComboHit* ch0 = &fChc->at(0);
  int ich = 0;
  for (int i=0; i<fNCh; i++) {
    const mu2e::ComboHit* ch = &fChc->at(i);
    int ind = ch-ch0;
    float corrected_time = ch->correctedTime();
    if (fMinDt > 0) {
      if (ch->energyDep() < fMinEDep) {
        bool close_enough = CloseEnough(corrected_time);
        if (not close_enough) continue;
      }
    }

    //    TComboHit* nt_ch = data->NewHit(i);

    int pln = ch->strawId().plane();
    int pnl = ch->strawId().panel();
    const mu2e::TrkPanelMap::Row* tpmd = fTpm->panel_map_by_offline_ind(pln,pnl);
    // ind is a hit index in the original Mu2e hit collection
    TComboHit* nt_ch   = data->NewHit(ind);
    ich++;
    nt_ch->fStrawID    = ch->strawId().asUint16();
    nt_ch->fNsh        = ch->nStrawHits();
    nt_ch->fZface      = tpmd->zface();
    nt_ch->fMnid       = tpmd->mnid();
    nt_ch->fTime       = corrected_time; // ch->correctedTime();
    nt_ch->fDrTime     = ch->driftTime();
    nt_ch->fX          = ch->pos().x();
    nt_ch->fY          = ch->pos().y();
    nt_ch->fZ          = ch->pos().z();
    nt_ch->fUx         = ch->uDir().x();
    nt_ch->fUy         = ch->uDir().y();
    nt_ch->fUres       = ch->uRes();
    nt_ch->fVres       = ch->vRes();
    nt_ch->fEDep       = ch->energyDep();

    // if (_debugMode  > 0) {
    //   if (_debugBit[1] != 0) {
    //     // printf("%8i %5i %3i %3i %3i %12.4f %12.4f %6.1f %6.1f %7.4f\n",
    //     //        _event->evn,
    //     //        (int) nt_sh->sid, pln, pnl, nt_sh->mnid,
    //     //        nt_sh->time, nt_sh->dt,
    //     //        nt_sh->tot0, nt_sh->tot1,
    //     //        nt_sh->edep);
    //   }
    // }
  }

  data->f_RunNumber   = rn_number;
  data->f_EventNumber = ev_number;
  
  return 0;
}

//-----------------------------------------------------------------------------
int InitComboHitBlock::ResolveLinks(TStnDataBlock* Block, AbsEvent* Event, int Mode) {
  std::string oname("InitComboHitBlock::ResolveLinks");

  if (fSaveShLinks == 0) return 0;

  const int evn = Event->event();
  const int rn  = Event->run();
  const int srn = Event->subRun();

  if (! Block->Initialized(evn,rn,srn)) return -1;
  if (  Block->LinksInitialized()     ) return  0;
  
  TComboHitBlock* ch_block = (TComboHitBlock*) Block;

  int nch = ch_block->NHits();
  for (int i=0; i<nch; i++) {
    auto ch = ch_block->Hit(i);
    // this one has a list of pointers to straw hits
    // rely on that the Mu2e collectiosn being vectors of things
    const mu2e::ComboHit* o_ch = ch->OfflineCh();
    int nsh = o_ch->nStrawHits();
    for (int j=0; j<nsh; j++) {
      uint16_t ind = o_ch->index(j);
      // assuming all straw hits are stored
      // to make that work in case not all straw hits are stored, need to store
      // offline index of each straw hit - why not ?
      // TObject should have space
      // same shoudl hold for all objects stored in Stntuple
      ch_block->ListOfShIndices()->Add(i,ind);
    }
  }
//-----------------------------------------------------------------------------
// mark links as initialized
//-----------------------------------------------------------------------------
  ch_block->fLinksInitialized = 1;

  return 0;
}
}
