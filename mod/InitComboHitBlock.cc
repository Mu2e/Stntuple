///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "Stntuple/mod/InitComboHitBlock.hh"
#include "Stntuple/mod/StntupleUtilities.hh"

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
  fSaveShLinks =  1;
  fMinDt       = -1;     // if < 0, save all hits
  fMinEDep     = 0.5e-3; // 0.5 keV
  fLastRun     = -1;
}
  
//-----------------------------------------------------------------------------
int InitComboHitBlock::InitDataBlock(TStnDataBlock* Block, AbsEvent* Event, int Mode) {
  std::string oname("InitComboHitBlock::InitDataBlock");

  int evn, run, srn; // , mc_flag(0); /*,n_combo_hits(0), n_straw_hits(0)*/

  evn = Event->event();
  run = Event->run();
  srn = Event->subRun();

  //  if (rn_number < 100000) mc_flag = 1; 

  if (Block->Initialized(evn,run,srn)) return 0;

  TComboHitBlock* data = (TComboHitBlock*) Block;
  data->Clear();
  
  if (fLastRun != run) {
    fTpm      = &fTpmh.get(Event->id());
    fLastRun  = run;
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

  if (fNCh > 0) {
    const mu2e::ComboHit* ch0 = &fChc->at(0);
    int ich = 0;
    for (int i=0; i<fNCh; i++) {
      const mu2e::ComboHit* ch = &fChc->at(i);
      int ind = ch-ch0;
      float corrected_time = ch->correctedTime();
      if (fStoreCloseHitsOnly) {
        if (ch->energyDep() < fMinEDep) {
          //        bool close_enough = stntuple::CloseEnough(corrected_time,fTcc,fMinDt);
          bool close_enough = stntuple::CloseEnough(ch->correctedTime(),fTcc,fMinDt);
          if (not close_enough) {
            // accidentally, the time clustering window may be large enough
            // so hits included into the time cluster are outside the checked window
            // make sure the hit is not in the list explicitly
            int ntc = fTcc->size();
            for (int itc=0; itc<ntc; itc++) {
              const mu2e::TimeCluster* tc = &fTcc->at(itc);
              int nch = tc->hits().size();
              for (int j=0; j<nch; j++) {
                // store the CH index in the initial collection
                uint16_t tc_ind = tc->hits().at(j);
                if (ind == tc_ind) {
                  close_enough = true;
                  break;
                }
              }
              if (close_enough) break;
            }
            if (not close_enough) continue;
          }
        }
      }
//-----------------------------------------------------------------------------
// store combo hit
//-----------------------------------------------------------------------------
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

      nt_ch->fOfflineCh  = ch;

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
  }
  data->f_RunNumber    = run;
  data->f_EventNumber  = evn;
  data->f_SubrunNumber = srn;
  
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
    // rely on that the Mu2e collections being vectors of things
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
