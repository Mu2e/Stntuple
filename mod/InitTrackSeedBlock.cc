//-----------------------------------------------------------------------------
//  Apr 2016 G. Pezzullo: initialization of the MU2E STNTUPLE TrackSeed block
//
//-----------------------------------------------------------------------------
#include <cstdio>
#include "TROOT.h"
#include "TFolder.h"
#include "TLorentzVector.h"
#include "TDatabasePDG.h"
#include "TParticlePDG.h"


#include "Stntuple/obj/TStnDataBlock.hh"
#include "Stntuple/obj/TStnEvent.hh"

#include "Stntuple/obj/TStnHelix.hh"
#include "Stntuple/obj/TStnHelixBlock.hh"

#include "Stntuple/obj/TStnTrackSeed.hh"
#include "Stntuple/obj/TStnTrackSeedBlock.hh"

#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Event.h"

#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"

#include "RecoDataProducts/inc/KalSeed.hh"
#include "RecoDataProducts/inc/CaloCluster.hh"
#include "RecoDataProducts/inc/StrawHitCollection.hh"
#include "RecoDataProducts/inc/StrawHit.hh"

#include "MCDataProducts/inc/SimParticle.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "MCDataProducts/inc/StrawGasStep.hh"
#include "MCDataProducts/inc/StrawDigiMCCollection.hh"
//-----------------------------------------------------------------------------
// assume that the collection name is set, so we could grab it from the event
//-----------------------------------------------------------------------------
int  StntupleInitMu2eTrackSeedBlock(TStnDataBlock* Block, AbsEvent* Evt, int Mode) {

  const mu2e::KalSeedCollection*  list_of_trackSeeds;
  const mu2e::ComboHitCollection* list_of_combo_hits(0);

  char                      trkSeed_module_label[100], trkSeed_description[100]; 
  char                      cmbh_module_label   [100], cmbh_description[100];

  int  ev_number, rn_number;

  ev_number = Evt->event();
  rn_number = Evt->run();

  if (Block->Initialized(ev_number,rn_number)) return 0;

  TStnTrackSeed*              trackSeed(0);
  TStnTrackSeedBlock*         data = (TStnTrackSeedBlock*) Block;

  data->Clear();

  data->GetModuleLabel("mu2e::KalSeedCollection", trkSeed_module_label);
  data->GetDescription("mu2e::KalSeedCollection", trkSeed_description );

  art::Handle<mu2e::KalSeedCollection>               trackSeed_handle;
  if (trkSeed_description[0] == 0) Evt->getByLabel(trkSeed_module_label, trackSeed_handle);
  else                          Evt->getByLabel(trkSeed_module_label, trkSeed_description, trackSeed_handle);
  list_of_trackSeeds = (mu2e::KalSeedCollection*) &(*trackSeed_handle);

  data->GetModuleLabel("mu2e::ComboHitCollection",cmbh_module_label);
  data->GetDescription("mu2e::ComboHitCollection",cmbh_description );

  art::Handle<mu2e::ComboHitCollection> shHandle;
  if (cmbh_module_label[0] != 0) {
    if (cmbh_description[0] == 0) Evt->getByLabel(cmbh_module_label,shHandle);
    else                          Evt->getByLabel(cmbh_module_label,cmbh_description,shHandle);
    if (shHandle.isValid()) list_of_combo_hits = shHandle.product();
  }
  
  
  const mu2e::KalSeed         *trkSeed(0);
  int                           ntrkseeds(0);
  
  const mu2e::CaloCluster       *cluster(0);
 
  ntrkseeds = list_of_trackSeeds->size();    

  char                 makeSD_module_label[100];
  data->GetModuleLabel("mu2e::StrawDigiMCCollection",makeSD_module_label);

  const mu2e::StrawDigiMCCollection* mcdigis(0);
  art::Handle<mu2e::StrawDigiMCCollection> mcdH;
  Evt->getByLabel(makeSD_module_label, mcdH);
  mcdigis = mcdH.product();
  
  TParticlePDG*        part(0);
  TDatabasePDG*        pdg_db = TDatabasePDG::Instance();
 
  for (int i=0; i<ntrkseeds; i++) {
    std::vector<int>     hits_simp_id, hits_simp_index, hits_simp_z;
    
    trackSeed                  = data->NewTrackSeed();
    trkSeed                    = &list_of_trackSeeds->at(i);
    cluster                    = trkSeed->caloCluster().get();

    if (cluster != 0){
      trackSeed->fClusterTime    = cluster->time();
      trackSeed->fClusterEnergy  = cluster->energyDep();
      trackSeed->fClusterX       = cluster->cog3Vector().x();
      trackSeed->fClusterY       = cluster->cog3Vector().y();
      trackSeed->fClusterZ       = cluster->cog3Vector().z();
    }else {
      trackSeed->fClusterTime    = 0; 
      trackSeed->fClusterEnergy  = -1; 
      trackSeed->fClusterX       = -1.e6; 
      trackSeed->fClusterY       = -1.e6; 
      trackSeed->fClusterZ       = -1.e6; 
    }
    
    mu2e::KalSegment kalSeg  = trkSeed->segments().at(0);//take the KalSegment closer to the entrance of the tracker
    trackSeed->fTrackSeed    = trkSeed;
    trackSeed->fNHits        = trkSeed->hits().size();
    trackSeed->fT0           = trkSeed->t0()._t0;
    trackSeed->fT0Err        = trkSeed->t0()._t0err;     
    trackSeed->fD0           = kalSeg.helix().d0();
    trackSeed->fFlt0         = trkSeed->flt0();     
    trackSeed->fTanDip       = kalSeg.helix().tanDip();
    trackSeed->fP            = kalSeg.mom();
    trackSeed->fPt           = trackSeed->fP*std::cos( std::atan(trackSeed->fTanDip));
    trackSeed->fZ0           = kalSeg.helix().z0();

    trackSeed->fChi2         = trkSeed->chisquared();
    trackSeed->fFitCons      = trkSeed->fitConsistency();

    //now loop over the hits to search the particle that generated the track
    int                      nsh = trkSeed->hits().size();
    const mu2e::ComboHit*    hit(0), *hit_0(0);
    hit_0     = &list_of_combo_hits->at(0);

    int                   loc(-1);

    float    first_hit_z(0), last_hit_z(0);
    float    lambda = kalSeg.helix().tanDip()/kalSeg.helix().omega();

    for (int j=0; j<nsh; ++j){
      int  hitIndex  = int(trkSeed->hits().at(j).index());
      hit            = &list_of_combo_hits->at(hitIndex);
      loc            = hit - hit_0;
      
      if(j==0) first_hit_z = hit->pos().z();
      else if(j==nsh-1) last_hit_z = hit->pos().z();
   
      const mu2e::StrawGasStep* step(nullptr);
      if (mcdigis) {
	const mu2e::StrawDigiMC* sdmc = &mcdigis->at(loc);
	if (sdmc->wireEndTime(mu2e::StrawEnd::cal) < sdmc->wireEndTime(mu2e::StrawEnd::hv)) {
	  step = sdmc->strawGasStep(mu2e::StrawEnd::cal).get();
	}
	else {
	  step = sdmc->strawGasStep(mu2e::StrawEnd::hv ).get();
	}
      }

      if (step) {
	art::Ptr<mu2e::SimParticle> const& simptr = step->simParticle(); 
	int sim_id        = simptr->id().asInt();

	hits_simp_id.push_back   (sim_id); 
	hits_simp_index.push_back(loc);
	hits_simp_z.push_back(step->position().z());
     }
    }

//-----------------------------------------------------------------------------
// calculate the number of loops made
//-----------------------------------------------------------------------------
    trackSeed->fNLoops = (last_hit_z - first_hit_z)/(fabs(lambda)*2.*M_PI);
 
//-----------------------------------------------------------------------------
// find the simparticle that created the majority of the hits
//-----------------------------------------------------------------------------
    int     max(0), mostvalueindex(-1), mostvalue= hits_simp_id[0];
    float   dz_most(1e4);
    for (int k=0; k<(int)hits_simp_id.size(); ++k){
      int co = (int)std::count(hits_simp_id.begin(), hits_simp_id.end(), hits_simp_id[k]);
      if ( (co>0) &&  (co>max)) {
	float  dz      = std::fabs(hits_simp_z[k]);
	if (dz < dz_most){
	  dz_most        = dz;
	  max            = co;
	  mostvalue      = hits_simp_id[k];
	  mostvalueindex = hits_simp_index[k];
	}
      }
    }

    trackSeed->fSimpId1Hits = max;
    trackSeed->fSimpId2Hits = -1;

    const mu2e::StrawGasStep* step(nullptr);
    const mu2e::StrawDigiMC* sdmc = &mcdigis->at(mostvalueindex);
    if (mcdigis) {
      if (sdmc->wireEndTime(mu2e::StrawEnd::cal) < sdmc->wireEndTime(mu2e::StrawEnd::hv)) {
	step = sdmc->strawGasStep(mu2e::StrawEnd::cal).get();
      }
      else {
	step = sdmc->strawGasStep(mu2e::StrawEnd::hv ).get();
      }
    }
    const mu2e::SimParticle * sim (0);

    if (step) {
      art::Ptr<mu2e::SimParticle> const& simptr = step->simParticle(); 
      trackSeed->fSimpPDG1    = simptr->pdgId();
      art::Ptr<mu2e::SimParticle> mother = simptr;
      part   = pdg_db->GetParticle(trackSeed->fSimpPDG1);

      while(mother->hasParent()) mother = mother->parent();
      sim = mother.operator ->();

      trackSeed->fSimpPDGM1   = sim->pdgId();
      
      double   px = simptr->startMomentum().x();
      double   py = simptr->startMomentum().y();
      double   pz = simptr->startMomentum().z();
      double   mass  (-1.);
      double   energy(-1.);
      if (part) {
	mass   = part->Mass();
	energy = sqrt(px*px+py*py+pz*pz+mass*mass);
      }
      trackSeed->fMom1.SetPxPyPzE(px,py,pz,energy);

      const CLHEP::Hep3Vector* sp = &simptr->startPosition();
      trackSeed->fOrigin1.SetXYZT(sp->x(),sp->y(),sp->z(),simptr->startGlobalTime());
    }
    
    //look for the second most frequent hit
    if (max != nsh && max > 0){
      int   secondmostvalueindex(-1);
      max     = 0;//reset max
      dz_most = 1e4;

      for (int k=0; k<(int)hits_simp_id.size(); ++k){
	int value = hits_simp_id[k];
	int co = (int)std::count(hits_simp_id.begin(), hits_simp_id.end(), value);
	if ( (co>0) && (co>max) && (value != mostvalue)) {
	  float  dz      = std::fabs(hits_simp_z[k]);
	  if (dz < dz_most){
	    max                  = co;
	    dz_most              = dz;
	    secondmostvalueindex = hits_simp_index[k];
	  }
	}
      }
      trackSeed->fSimpId2Hits = max;

      if (secondmostvalueindex >=0){
	const mu2e::StrawGasStep* step(nullptr);
	const mu2e::StrawDigiMC* sdmc = &mcdigis->at(secondmostvalueindex);
	if (sdmc->wireEndTime(mu2e::StrawEnd::cal) < sdmc->wireEndTime(mu2e::StrawEnd::hv)) {
	  step = sdmc->strawGasStep(mu2e::StrawEnd::cal).get();
	}
	else {
	  step = sdmc->strawGasStep(mu2e::StrawEnd::hv ).get();
	}
		
	const mu2e::SimParticle * sim (nullptr);

	if (step) {
	  art::Ptr<mu2e::SimParticle> const& simptr = step->simParticle(); 
	  trackSeed->fSimpPDG2    = simptr->pdgId();
	  art::Ptr<mu2e::SimParticle> mother = simptr;
	  part   = pdg_db->GetParticle(trackSeed->fSimpPDG2);

	  while(mother->hasParent()) mother = mother->parent();
	  sim = mother.operator ->();

	  trackSeed->fSimpPDGM2   = sim->pdgId();
      
	  double   px = simptr->startMomentum().x();
	  double   py = simptr->startMomentum().y();
	  double   pz = simptr->startMomentum().z();
	  double   mass  (-1.);
	  double   energy(-1.);
	  if (part) {
	    mass   = part->Mass();
	    energy = sqrt(px*px+py*py+pz*pz+mass*mass);
	  }
	  trackSeed->fMom2.SetPxPyPzE(px,py,pz,energy);

	  const CLHEP::Hep3Vector* sp = &simptr->startPosition();
	  trackSeed->fOrigin2.SetXYZT(sp->x(),sp->y(),sp->z(),simptr->startGlobalTime());
	}      
      }
    }
  }
					// on return set event and run numbers
					// to mark block as initialized
  data->f_RunNumber   = rn_number;
  data->f_EventNumber = ev_number;

  return 0;
}

//_____________________________________________________________________________
Int_t StntupleInitMu2eTrackSeedBlockLinks(TStnDataBlock* Block, AbsEvent* AnEvent, int Mode) {
  // Mu2e version, do nothing

  // Int_t  ev_number, rn_number;

  // ev_number = AnEvent->event();
  // rn_number = AnEvent->run();
					// do not do initialize links 2nd time

  if (Block->LinksInitialized()) return 0;

  TStnTrackSeedBlock* tsb;
  TStnTrackSeed*      trkseed;
  TStnHelixBlock*     hb;
  TStnHelix          *helix;
  TStnEvent*          ev;

  const mu2e::HelixSeed* khelix, *fkhelix;
  const mu2e::KalSeed*   kseed; 

  tsb    = (TStnTrackSeedBlock*) Block;

  char                kseed_module_label[100], helix_block_name[100];
  tsb->GetModuleLabel("mu2e::KalSeedCollection"  , kseed_module_label);
  tsb->GetModuleLabel("HelixBlockName"           , helix_block_name);

  ev     = Block->GetEvent();
  hb     = (TStnHelixBlock*)     ev->GetDataBlock(helix_block_name);
  
  int    ntrkseed = tsb ->NTrackSeeds();
  int    nhelix(0);

  if (hb) nhelix = hb->NHelices();
  
  for (int i=0; i<ntrkseed; i++) {
    trkseed = tsb->TrackSeed(i);
    kseed   = trkseed->fTrackSeed;

    fkhelix  = kseed->helix().get();
    int  helixIndex(-1);
    for (int j=0; j<nhelix; ++j){
      helix  = hb->Helix(j);
      khelix = helix->fHelix;
      if ( khelix == fkhelix){
	helixIndex = j;
	break;
      }    
    }
    
    if (helixIndex < 0) {
      printf(">>> InitTrackSeedBlock ERROR: %s trackseed %i -> no HelixSeed associated\n", kseed_module_label, i);
      continue;
    }
    trkseed->SetHelixIndex(helixIndex);
  }
//-----------------------------------------------------------------------------
// mark links as initialized
//-----------------------------------------------------------------------------
  Block->SetLinksInitialized();

  return 0;
}

