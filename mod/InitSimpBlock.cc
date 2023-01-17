///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
#include "TLorentzVector.h"
#include "TDatabasePDG.h"
#include "TParticlePDG.h"
#include <vector>

#include "messagefacility/MessageLogger/MessageLogger.h"

#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Selector.h"
#include "canvas/Utilities/InputTag.h"

#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/VirtualDetector.hh"

#include "DataProducts/inc/VirtualDetectorId.hh"

#include "MCDataProducts/inc/GenParticle.hh"
#include "MCDataProducts/inc/SimParticle.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "MCDataProducts/inc/StrawDigiMCCollection.hh"

#include "RecoDataProducts/inc/StrawHitCollection.hh"

#include "Stntuple/obj/TSimpBlock.hh"

#include "Stntuple/mod/InitStntupleDataBlocks.hh"


#include "Stntuple/mod/THistModule.hh"
#include "Stntuple/base/TNamedHandle.hh"

#include "Stntuple/mod/InitSimpBlock.hh"

#include "messagefacility/MessageLogger/MessageLogger.h"

//-----------------------------------------------------------------------------
// fill SimParticle's data block
//-----------------------------------------------------------------------------
int StntupleInitSimpBlock::InitDataBlock(TStnDataBlock* Block, AbsEvent* AnEvent, int mode) 
{
  const char* oname = {"StntupleInitMu2eSimpBlock"};

  //  static int    initialized(0);
  //  static char   g4_module_label  [100], g4_description  [100];

  // static double _min_energy;

  std::vector<art::Handle<mu2e::SimParticleCollection>> list_of_sp;

  const mu2e::SimParticleCollection*       simp_coll(0);
  const mu2e::SimParticle*                 sim(0);
  const mu2e::StrawHitCollection*          list_of_straw_hits(0);
  const mu2e::StrawDigiMCCollection*       mcdigis(nullptr);

  double        px, py, pz, energy;
  int           id, parent_id, generator_id, n_straw_hits(0), nhits;
  int           pdg_code, start_vol_id, end_vol_id, creation_code, termination_code;
  TSimParticle* simp;

  TSimpBlock* simp_block = (TSimpBlock*) Block;
  simp_block->Clear();

  art::Handle<mu2e::StrawHitCollection> shHandle;
  if (! fStrawHitCollTag.empty()) {
    bool ok = AnEvent->getByLabel(fStrawHitCollTag,shHandle);
    if (ok) {
      list_of_straw_hits = shHandle.product();
      n_straw_hits      = list_of_straw_hits->size();
    }
  }

  art::Handle<mu2e::StrawDigiMCCollection> mcdH;

  if (! fStrawDigiMCCollTag.empty()) {
    bool ok = AnEvent->getByLabel(fStrawDigiMCCollTag,mcdH);
    if (ok) {
      mcdigis = mcdH.product();
    }
  }

  mu2e::GeomHandle<mu2e::VirtualDetector> vdg;

  art::Handle<mu2e::SimParticleCollection> simp_handle;
  if (! fSimpCollTag.empty()) {
    AnEvent->getByLabel(fSimpCollTag,simp_handle);
  }
//-----------------------------------------------------------------------------
// load simulation time offsets for this event: timeOffsets may not be defined yet (stages 1, 2, 3)
//-----------------------------------------------------------------------------
  if (fTimeOffsets) fTimeOffsets->updateMap(*AnEvent);
//-----------------------------------------------------------------------------
// figure out how many straw hits each particle has produced
//-----------------------------------------------------------------------------
  std::vector<int> vid, vin;
  int np_with_straw_hits(0);   // number of particles with straw hits

  if (n_straw_hits > 0) {

    vid.reserve(n_straw_hits);
    vin.reserve(n_straw_hits);

    for (int i=0; i<n_straw_hits; i++) vin[i] = 0;

    const mu2e::StrawGasStep* step(nullptr);

    for (int i=0; i<n_straw_hits; i++) {
      const mu2e::StrawDigiMC* mcdigi = &mcdigis->at(i);

      if (mcdigi->wireEndTime(mu2e::StrawEnd::cal) < mcdigi->wireEndTime(mu2e::StrawEnd::hv)) {
      	step = mcdigi->strawGasStep(mu2e::StrawEnd::cal).get();
      }
      else {
      	step = mcdigi->strawGasStep(mu2e::StrawEnd::hv ).get();
      }
    
      if (step) {
//------------------------------------------------------------------------------
// looking at the ppbar annihilation events - 
// searching for the mother works for delta-electrons, but not otherwise, do we really need 
// that? - it worked before, though... comment outp for teh time being
//-----------------------------------------------------------------------------
	art::Ptr<mu2e::SimParticle> const& simptr = step->simParticle(); 
	art::Ptr<mu2e::SimParticle> mother        = simptr;
	//	while(mother->hasParent())  mother        = mother->parent();
	const mu2e::SimParticle*    sim           = mother.get();
	  
	int sim_id = sim->id().asInt();

	int found  = 0;
	for (int ip=0; ip<np_with_straw_hits; ip++) {
	  if (sim_id == vid[ip]) {
	    vin[ip] += 1;
	    found    = 1;
	    break;
	  }
	}

	if (found == 0) {
	  vid[np_with_straw_hits] = sim_id;
	  vin[np_with_straw_hits] = 1;
	  np_with_straw_hits      = np_with_straw_hits+1;
	}
      }
    }
  }

  if (simp_handle.isValid()) {
    simp_coll = simp_handle.product();

    for (mu2e::SimParticleCollection::const_iterator ip = simp_coll->begin(); ip != simp_coll->end(); ip++) {
      sim      = &ip->second;
      const mu2e::GenParticle* genp;

      id        = sim->id().asInt();
      parent_id = -1;
//------------------------------------------------------------------------------
// count number of straw hits produced by the particle
//-----------------------------------------------------------------------------
      nhits = 0;
      
      for (int i=0; i<np_with_straw_hits; i++) {
	if (vid[i] == id) {
	  nhits = vin[i];
	  break;
	}
      }
//-----------------------------------------------------------------------------
// a semi-kludge: store e+ and e- from an external photon conversion
//-----------------------------------------------------------------------------
      if (sim->parent()) { 
	parent_id = sim->parent()->id().asInt();
	genp      = sim->parent()->genParticle().get();
      }
      else {
	genp      = sim->genParticle().get();
      }

      pdg_code         = (int) sim->pdgId();

      if (genp) generator_id = genp->generatorId().id();   // ID of the MC generator
      else      generator_id = -1;

      if ((fGenProcessID > 0) && (generator_id != fGenProcessID)) continue;
      if ((fPdgID       != 0) && (pdg_code     != fPdgID       )) continue;
      
      creation_code    = sim->creationCode();
      termination_code = sim->stoppingCode();

      start_vol_id     = sim->startVolumeIndex();
      end_vol_id       = sim->endVolumeIndex();
      
      px     = sim->startMomentum().x();
      py     = sim->startMomentum().y();
      pz     = sim->startMomentum().z();
      energy = sim->startMomentum().e();
//-----------------------------------------------------------------------------
// by default, do not store low energy SimParticles not making hits in the tracker
//-----------------------------------------------------------------------------
      const CLHEP::Hep3Vector* sp = &sim->startPosition();

      if (fMinSimpEnergy >= 0) {
	if ((nhits == 0) and (energy < fMinSimpEnergy))       continue;
      }

      double time_offsets = 0;
      if (fTimeOffsets) time_offsets = fTimeOffsets->totalTimeOffset(sim);

      simp   = simp_block->NewParticle(id, parent_id, pdg_code, 
				       creation_code, termination_code,
				       start_vol_id, end_vol_id,
				       generator_id);
      simp->SetStartMom(px, py, pz, energy);
//-----------------------------------------------------------------------------
// for SimParticles, do not wrap the time around the microbunch, do that only for hits
//-----------------------------------------------------------------------------
      simp->SetStartPos(sp->x(),sp->y(),sp->z(),sim->startGlobalTime()+time_offsets);
      simp->SetStartProperTime(sim->startProperTime());
      simp->SetEndMom  (sim->endMomentum().x(),
			sim->endMomentum().y(),
			sim->endMomentum().z(),
			sim->endMomentum().e());
      const CLHEP::Hep3Vector* ep = &sim->endPosition();
      simp->SetEndPos(ep->x(),ep->y(),ep->z(),sim->endGlobalTime()+time_offsets);
      simp->SetEndProperTime(sim->endProperTime());

      simp->SetNStrawHits(nhits);
//-----------------------------------------------------------------------------
// particle parameters at virtual detectors
//-----------------------------------------------------------------------------
      if (vdg->nDet() > 0) {
	art::Handle<mu2e::StepPointMCCollection> vdhits;
	AnEvent->getByLabel(fVDHitsCollTag,vdhits);
	if (!vdhits.isValid()) {
	  char warning[100];
	  sprintf(warning,"StepPointMCCollection %s NOT FOUND\n",fVDHitsCollTag.encode().data());
	  mf::LogWarning(oname) << warning;
	}
	else {
	  int nvdhits = vdhits->size();
	  for (int i=0; i<nvdhits; i++) {
	    const mu2e::StepPointMC* hit = &(*vdhits)[i];
	    
	    mu2e::VirtualDetectorId vdid(hit->volumeId());
	    
	    if (vdid.id() == mu2e::VirtualDetectorId::ST_Out) {
	      
	      const mu2e::SimParticle* sim = hit->simParticle().get();
	      
	      if (sim == NULL) {
		printf(">>> ERROR: %s sim == NULL\n",oname);
	      }
	      int sim_id = sim->id().asInt();
	      if (sim_id == id) {
		simp->SetMomTargetEnd(hit->momentum().mag());
	      }
	    }
	    else if (vdid.isTrackerFront()) {
	      art::Ptr<mu2e::SimParticle> const& simptr = hit->simParticle();
	      const mu2e::SimParticle* sim = simptr.get();
	      
	      if (sim == NULL) {
		printf("[%s] ERROR: sim == NULL. CONTINUE.\n",oname);
	      }
	      else {
		int sim_id = sim->id().asInt();
		if (sim_id == id) {
		  simp->SetMomTrackerFront(hit->momentum().mag());
		}
	      }
	    }
	    
	  }
	}
      }
    }
  }
  else {
//-----------------------------------------------------------------------------
// no SIMP collection
//-----------------------------------------------------------------------------
    mf::LogWarning(oname) << " WARNING in " << oname << ": SimParticleCollection " 
			  << fSimpCollTag.encode() << " not found, BAIL OUT\n";
    return -1;
  }

  return 0;
}
