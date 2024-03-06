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

#include "Offline/GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "Offline/GlobalConstantsService/inc/ParticleDataList.hh"

#include "Offline/GeometryService/inc/GeometryService.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/GeometryService/inc/VirtualDetector.hh"

#include "Offline/DataProducts/inc/VirtualDetectorId.hh"

#include "Offline/MCDataProducts/inc/GenParticle.hh"
#include "Offline/MCDataProducts/inc/SimParticle.hh"
#include "Offline/MCDataProducts/inc/StepPointMC.hh"
#include "Offline/MCDataProducts/inc/StrawDigiMC.hh"
#include "Offline/MCDataProducts/inc/PrimaryParticle.hh"

#include "Offline/RecoDataProducts/inc/StrawHit.hh"

#include "Stntuple/obj/TSimpBlock.hh"

#include "Stntuple/mod/StntupleUtilities.hh"
#include "Stntuple/mod/InitStntupleDataBlocks.hh"


// #include "Stntuple/mod/THistModule.hh"
#include "Stntuple/base/TNamedHandle.hh"

#include "Stntuple/mod/InitSimpBlock.hh"

#include "messagefacility/MessageLogger/MessageLogger.h"

//-----------------------------------------------------------------------------
// fill SimParticle's data block
//-----------------------------------------------------------------------------
int StntupleInitSimpBlock::InitDataBlock(TStnDataBlock* Block, AbsEvent* AnEvent, int Mode) {
  const char* oname = {"Stntuple::InitSimpBlock"};

  std::vector<art::Handle<mu2e::SimParticleCollection>> list_of_sp;

  const mu2e::SimParticleCollection*       simp_coll(nullptr);
  const mu2e::SimParticle*                 sim      (nullptr);
  const mu2e::StrawHitCollection*          shColl   (nullptr);
  const mu2e::StrawDigiMCCollection*       mcdigis  (nullptr);

  double        px, py, pz, energy;
  int           id, parent_id, process_id, nsh(0), nhits;
  int           pdg_code, start_vol_id, end_vol_id, creation_code, termination_code;
  TSimParticle* simp;

  TSimpBlock* simp_block = (TSimpBlock*) Block;
  simp_block->Clear();

  auto const& ptable = mu2e::GlobalConstantsHandle<mu2e::ParticleDataList>();

  art::Handle<mu2e::StrawHitCollection> shcH;
  if (! fShCollTag.empty()) {
    bool ok = AnEvent->getByLabel(fShCollTag,shcH);
    if (ok) {
      shColl = shcH.product();
      nsh    = shColl->size();
    }
    else {
      mf::LogWarning(oname) << " WARNING line " << __LINE__ << ": no StrawHitCollection tag=" 
			    << fShCollTag.encode().data() <<  " found";
    }
  }

  art::Handle<mu2e::StrawDigiMCCollection> mcdH;

  if (! fStrawDigiMCCollTag.empty()) {
    bool ok = AnEvent->getByLabel(fStrawDigiMCCollTag,mcdH);
    if (ok) mcdigis = mcdH.product();
    else {
      mf::LogWarning(oname) << " WARNING line " << __LINE__ << ": no StrawDigiMCCollection tag=" 
			    << fShCollTag.encode().data() <<  " found";
    }
  }

  art::Handle<mu2e::SimParticleCollection> simp_handle;
  if (! fSimpCollTag.empty()) {
    bool ok = AnEvent->getByLabel(fSimpCollTag,simp_handle);
    if (! ok) {
      mf::LogWarning(oname) << " WARNING line " << __LINE__ 
			    << ": SimpCollection:" 
			    << fSimpCollTag.encode().data() << " NOT FOUND";
      return -1;
    }
  }

  mu2e::GeomHandle<mu2e::VirtualDetector> vdg;

  int nvdhits(0);
  art::Handle<mu2e::StepPointMCCollection> vdhcH;
  const mu2e::StepPointMCCollection*       vdhc(nullptr);

  if (! fVDHitsCollTag.empty()) {
    AnEvent->getByLabel(fVDHitsCollTag,vdhcH);
    if (!vdhcH.isValid()) {
      mf::LogWarning(oname) << " WARNING in " << oname << ":" << __LINE__ 
                            << ": StepPointMCCollection:" 
                            << fVDHitsCollTag.encode().data() << " NOT FOUND";
    }
    else {
      vdhc    = vdhcH.product();
      nvdhits = vdhc->size();
    }
  }
//-----------------------------------------------------------------------------
// figure out how many straw hits each particle has produced
// StrawHitCollection, StrawDigiCollection, and StrawDigiMCCollection are all parallel,
// indexed with the same index
//-----------------------------------------------------------------------------
  std::vector<int>  vid, vin;
  int               np_with_straw_hits(0);  // number of particles with straw hits
                                            // straw hit ID's, per particle
  std::vector< std::vector<int>* > vshid(nsh+1);

  vid.reserve(nsh);
  vin.reserve(nsh);

  for (int i=0; i<nsh; i++) {
    vin  [i] = 0;
    vshid[i] = nullptr;
  }

  if ((nsh > 0) and (mcdigis != nullptr)) {

    const mu2e::StrawGasStep* step(nullptr);

    for (int i=0; i<nsh; i++) {
      const mu2e::StrawDigiMC* mcdigi = &mcdigis->at(i);

      step = mcdigi->earlyStrawGasStep().get();
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
	    vshid[ip]->push_back(i);
	    found    = 1;
	    break;
	  }
	}

	if (found == 0) {
	  vshid[np_with_straw_hits]      = new std::vector<int>;
	  vshid[np_with_straw_hits]->push_back(i);
	  vid[np_with_straw_hits]        = sim_id;
	  vin[np_with_straw_hits]        = 1;
	  np_with_straw_hits             = np_with_straw_hits+1;
	}
      }
    }
  }
//-----------------------------------------------------------------------------
// figure out the primary particle
//-----------------------------------------------------------------------------
  art::Handle<mu2e::PrimaryParticle> pp_handle;
  const mu2e::PrimaryParticle*       pp(nullptr);
  const mu2e::SimParticle*           primary(nullptr);

  if (! fPrimaryParticleTag.empty()) {
    AnEvent->getByLabel(fPrimaryParticleTag,pp_handle);

    if (pp_handle.isValid()) {
      pp            = pp_handle.product();
      if (pp->primarySimParticles().size() > 0) {
	primary       = pp->primarySimParticles().front().get();
	fGenProcessID = primary->creationCode();
	fPdgID        = primary->pdgId();
      }
    }
  }

  if (simp_handle.isValid()) {
    simp_coll = simp_handle.product();

    for (mu2e::SimParticleCollection::const_iterator ip = simp_coll->begin(); ip != simp_coll->end(); ip++) {
      sim      = &ip->second;

      id        = sim->id().asInt();
      parent_id = -1;
//------------------------------------------------------------------------------
// count number of straw hits produced by the particle
//-----------------------------------------------------------------------------
      nhits  = 0;
      std::vector<int>* v_shid(nullptr);

      for (int i=0; i<np_with_straw_hits; i++) {
	if (vid[i] == id) {
	  nhits    = vin  [i];
	  v_shid   = vshid[i];
	  vshid[i] = nullptr;
	  break;
	}
      }
//-----------------------------------------------------------------------------
// need to store e+ and e- from an external photon conversion
// it is possible that will need to lop over the primary particles
//-----------------------------------------------------------------------------
      pdg_code         = (int) sim->pdgId();
      process_id       = sim->creationCode();

      if (pp != nullptr) {
	int found = 0;
	for (auto pr : pp->primarySimParticles()) {
	  if (pr.get() == sim) {
	    found = 1;
	    break;
	  }
	}
	if (found == 0)                                     continue;
      }
//-----------------------------------------------------------------------------
// if primary particle is not defined, or defined incorrectly (!) 
// store all particles
//-----------------------------------------------------------------------------
      creation_code    = sim->creationCode();
      termination_code = sim->stoppingCode();

      start_vol_id     = sim->startVolumeIndex();
      end_vol_id       = sim->endVolumeIndex();
      
      px               = sim->startMomentum().x();
      py               = sim->startMomentum().y();
      pz               = sim->startMomentum().z();
      double ptot      = sim->startMomentum().vect().mag();
      energy           = sim->startMomentum().e();
//-----------------------------------------------------------------------------
// by default, do not store low energy SimParticles not making hits in the tracker
//-----------------------------------------------------------------------------
      const CLHEP::Hep3Vector sp = sim->startPosition();

      if (fMinSimpMomentum >= 0) {
	if ((nhits == 0) or  (ptot < fMinSimpMomentum))     continue;
      }

      simp   = simp_block->NewParticle(id, parent_id, pdg_code        , 
				       creation_code, termination_code,
				       start_vol_id , end_vol_id      ,
				       process_id);
      simp->SetStartMom(px, py, pz, energy);
      simp->SetStartPos(sp.x(),sp.y(),sp.z(),sim->startGlobalTime());
//-----------------------------------------------------------------------------
// proper time: if StepPointMC collection is defined, make an attempt to correct 
// things for multistage simulation, where the proper time is calculated 
// for each stage separately
//-----------------------------------------------------------------------------
      // double tau               = gc->getParticleLifetime(sim->pdgId());
      double tau = ptable->particle(sim->pdgId()).lifetime()*1.e9;  // convert to ns
      if (tau == 0) tau = 1.e20;

      double start_proper_time = sim->startProperTime();
      double dpt               = stntuple::get_proper_time(sim);

      start_proper_time = start_proper_time+dpt;

      simp->SetStartProperTime(start_proper_time/tau);

      simp->SetEndMom  (sim->endMomentum().x(),
			sim->endMomentum().y(),
			sim->endMomentum().z(),
			sim->endMomentum().e());
      const CLHEP::Hep3Vector ep = sim->endPosition();
      simp->SetEndPos(ep.x(),ep.y(),ep.z(),sim->endGlobalTime());

      double end_proper_time = sim->endProperTime();

      end_proper_time = end_proper_time+dpt;

      simp->SetEndProperTime(end_proper_time/tau);

      simp->SetNStrawHits(nhits);
      simp->SetSimStage(sim->simStage());
      simp->SetSimParticle(sim);
      simp->SetShid(v_shid);
//-----------------------------------------------------------------------------
// particle parameters at virtual detectors -stored only for those which have 
// VD hits stored
//-----------------------------------------------------------------------------
      for (int i=0; i<nvdhits; i++) {
        const mu2e::StepPointMC* hit = &(*vdhc)[i];
	    
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
//-----------------------------------------------------------------------------
// memory clean up , reassigned ones will be deleted together with particles ? 
//-----------------------------------------------------------------------------
    for (int i=0; i<np_with_straw_hits; i++) {
      if (vshid[i]) delete vshid[i];
    }
  }
  else {
//-----------------------------------------------------------------------------
// no SIMP collection
//-----------------------------------------------------------------------------
    mf::LogWarning(oname) << " WARNING in " << oname << ":" << __LINE__ 
			  << " : SimParticleCollection " 
			  << fSimpCollTag.encode() << " not found";
    return -1;
  }

  return 0;
}
