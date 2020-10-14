///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////

#include "Stntuple/mod/InitStrawDataBlock.hh"

#include "RecoDataProducts/inc/ComboHit.hh"
#include "RecoDataProducts/inc/StrawHit.hh"

#include "MCDataProducts/inc/StrawDigiMC.hh"
#include "MCDataProducts/inc/StrawGasStep.hh"

#include "vector"

using std::vector ;
//-----------------------------------------------------------------------------
// in this case AbsEvent is just not used
//-----------------------------------------------------------------------------
int StntupleInitStrawDataBlock::InitDataBlock(TStnDataBlock* Block, AbsEvent* Event, int Mode) {

  int ev_number, rn_number, nhits(0);

  ev_number = Event->event();
  rn_number = Event->run();

  if (Block->Initialized(ev_number,rn_number)) return 0;

  TStrawDataBlock* data = (TStrawDataBlock*) Block;
  data->Clear();
//-----------------------------------------------------------------------------
//  straw hit information
//-----------------------------------------------------------------------------
  art::Handle<mu2e::ComboHitCollection> chch;
  const mu2e::ComboHitCollection*       chc (nullptr);

  art::Handle<mu2e::StrawHitCollection> shch;
  const mu2e::StrawHitCollection*       shc (nullptr);

  art::Handle<mu2e::StrawDigiMCCollection> sdmcch;
  const mu2e::StrawDigiMCCollection*       sdmcc (nullptr);

  //  int   nhits(0);

  if (! fStrawHitCollTag.empty() != 0) {
    bool ok = Event->getByLabel(fStrawHitCollTag,chch);
    if (ok) chc = chch.product();
    
    ok = Event->getByLabel(fStrawHitCollTag,shch);
    if (ok) { 
      shc   = shch.product();
      nhits = shc->size();
    }
  }
  
  if (! fStrawDigiMCCollTag.empty() != 0) {
    bool ok = Event->getByLabel(fStrawDigiMCCollTag,sdmcch);
    if (ok) sdmcc = sdmcch.product();
  }
  
  if (shc == nullptr) {
    printf(" >>> ERROR in StntupleInitMu2eStrawDataBlock: no StrawHitCollection with tag: %s. BAIL OUT\n",fStrawHitCollTag.encode().data());
    return -1;
  }
  else if (sdmcc == nullptr) {
    printf(" >>> ERROR in StntupleInitMu2eStrawDataBlock: no StrawDigiMCCollection with tag: %s. BAIL OUT\n",fStrawDigiMCCollTag.encode().data());
    return -1;
  }
//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------
  const mu2e::StrawGasStep* step (nullptr);
  const mu2e::SimParticle* sim;

  TStrawHitData*           hit; 

  int   pdg_id, mother_pdg_id, sim_id, gen_index;
  float mc_mom;

  const mu2e::ComboHit* ch0 = &chc->at(0);
 
  for (int i=0; i<nhits; i++) {
    const mu2e::StrawHit* sh = &shc->at(i);
    const mu2e::ComboHit* ch = &chc->at(i);

    size_t ih  = ch-ch0;
    vector<StrawDigiIndex> shids;
    chc->fillStrawDigiIndices((const art::Event&)*Event,ih,shids);

    const mu2e::StrawDigiMC* mcdigi = &sdmcc->at(shids[0]);

    if (mcdigi->wireEndTime(mu2e::StrawEnd::cal) < mcdigi->wireEndTime(mu2e::StrawEnd::hv)) {
      step = mcdigi->strawGasStep(mu2e::StrawEnd::cal).get();
    }
    else {
      step = mcdigi->strawGasStep(mu2e::StrawEnd::hv ).get();
    }

    hit = data->NewHit();

    if (step) {
      art::Ptr<mu2e::SimParticle> const& simptr = step->simParticle(); 
      art::Ptr<mu2e::SimParticle> mother = simptr;

      while (mother->hasParent()) mother = mother->parent();

      sim = mother.get();

      pdg_id        = simptr->pdgId();
      mother_pdg_id = sim->pdgId();

      if (simptr->fromGenerator()) gen_index = simptr->genParticle()->generatorId().id();
      else                         gen_index = -1;
      
      sim_id        = simptr->id().asInt();
      mc_mom        = step->momvec().mag();
    }
    else {
      pdg_id        = -1;
      mother_pdg_id = -1;
      gen_index     = -1;
      sim_id        = -1;
      mc_mom        = -1.;
    }

    hit->Set(sh->strawId().asUint16(), sh->time(), sh->dt(), sh->energyDep(),
	     pdg_id, mother_pdg_id, gen_index, sim_id, mc_mom);
    
  }

  data->f_RunNumber   = rn_number;
  data->f_EventNumber = ev_number;
  
  return 0;

}


