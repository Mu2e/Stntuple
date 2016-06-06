///////////////////////////////////////////////////////////////////////////////
// 2014-01-26 P.Murat
///////////////////////////////////////////////////////////////////////////////

#include "Stntuple/mod/InitStntupleDataBlocks.hh"
#include "Stntuple/obj/TTrackStrawHitBlock.hh"

#include "RecoDataProducts/inc/StrawHitCollection.hh"

#include "RecoDataProducts/inc/KalRepPtrCollection.hh"
#include "BTrkData/inc/TrkStrawHit.hh"

#include "MCDataProducts/inc/SimParticle.hh"
#include "MCDataProducts/inc/StepPointMC.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "MCDataProducts/inc/PtrStepPointMCVectorCollection.hh"

//-----------------------------------------------------------------------------
Int_t StntupleInitMu2eTrackStrawHitBlock(TStnDataBlock* Block, AbsEvent* AnEvent, int Mode) 
{
  // initialize COT data block with the `event' data
  // return -1, if bank doesn't exist, 0, if everything is OK

  int ev_number, rn_number, nhits;

  mu2e::KalRepPtrCollection*               list_of_kreps(0);

  static char   strh_module_label[100], strh_description[100];
  static char   krep_module_label[100], krep_description[100];

  ev_number = AnEvent->event();
  rn_number = AnEvent->run();

  if (Block->Initialized(ev_number,rn_number)) return 0;

  TTrackStrawHitBlock* data = (TTrackStrawHitBlock*) Block;
  data->Clear();
//-----------------------------------------------------------------------------
//  straw hit information
//-----------------------------------------------------------------------------
  data->GetModuleLabel("mu2e::StrawHitCollection",strh_module_label);
  data->GetDescription("mu2e::StrawHitCollection",strh_description );

  data->GetModuleLabel("mu2e::KalRepCollection"  ,krep_module_label);
  data->GetDescription("mu2e::KalRepCollection"  ,krep_description );

  art::Handle<mu2e::StrawHitCollection>       strh_handle;
  const mu2e::StrawHitCollection*             list_of_hits(0);

  art::Handle<mu2e::PtrStepPointMCVectorCollection> mcptr_handle;
  const mu2e::PtrStepPointMCVectorCollection* list_of_steps(0);

  art::Handle<mu2e::KalRepPtrCollection> krepsHandle;
  if (krep_module_label[0] != 0) {
    if (krep_description[0] == 0) AnEvent->getByLabel(krep_module_label,krepsHandle);
    else                          AnEvent->getByLabel(krep_module_label,krep_description, krepsHandle);
    if (krepsHandle.isValid())    list_of_kreps = (mu2e::KalRepPtrCollection*) krepsHandle.product();
  }

  if (strh_module_label[0] != 0) {
    if (strh_description[0] == 0) AnEvent->getByLabel(strh_module_label,strh_handle);
    else                          AnEvent->getByLabel(strh_module_label,strh_description,strh_handle);

    if (strh_handle.isValid()) list_of_hits = strh_handle.product();

    AnEvent->getByLabel(strh_module_label,"StrawHitMCPtr",mcptr_handle);
    if (mcptr_handle.isValid()) list_of_steps = mcptr_handle.product();
  }

  if (list_of_hits == NULL) {
    printf(" >>> ERROR in StntupleInitMu2eCalDataBlock: no list_of_hits. BAIL OUT\n");
    return -1;
  }
//-----------------------------------------------------------------------------
// 
//-----------------------------------------------------------------------------
  const mu2e::StrawHit*    sh; 
  const mu2e::StepPointMC* step;
  const mu2e::SimParticle* sim;

  TStrawHitData*           hit; 
  data->fNTracks = list_of_kreps->size();

  const KalRep*            krep;
  const TrkHitVector*       krep_hits; 

  int   pdg_id, mother_pdg_id, sim_id, gen_index;
  float mc_mom;

  data->fNTrackHits->Set(data->fNTracks);
  data->fFirst->Set(data->fNTracks);

  int nhtot = 0;
  for (int i=0; i<data->fNTracks; i++) {
    krep      = list_of_kreps->at(i).get();
    krep_hits = &krep->hitVector();
    nhits     = krep_hits->size();
    
    (*data->fNTrackHits)[i] = nhits;
    (*data->fFirst)     [i] = nhtot;

    nhtot += nhits;

    for (int ih=0; ih<nhits; ih++) {
      sh  = &list_of_hits->at(ih);
      mu2e::PtrStepPointMCVector const& mcptr (list_of_steps->at(ih));
      step = mcptr[0].operator ->();
      hit  = data->NewHit();

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
	mc_mom        = step->momentum().mag();
      }
      else {
	pdg_id        = -1;
	mother_pdg_id = -1;
	gen_index     = -1;
	sim_id        = -1;
	mc_mom        = -1.;
      }
      
      hit->Set(sh->strawIndex().asInt(), sh->time(), sh->dt(), sh->energyDep(),
	       pdg_id, mother_pdg_id, gen_index, sim_id, mc_mom);
    }
  }

  data->fNStrawHits   = nhtot;
  
  data->f_RunNumber   = rn_number;
  data->f_EventNumber = ev_number;
  
  return 0;
}

