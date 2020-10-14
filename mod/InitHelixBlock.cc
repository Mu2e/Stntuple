//-----------------------------------------------------------------------------
//  Apr 2016 G. Pezzullo: initialization of the MU2E STNTUPLE Helix block
//
//-----------------------------------------------------------------------------
#include <cstdio>
#include "TROOT.h"
#include "TFolder.h"
#include "TLorentzVector.h"
#include "TVector2.h"
#include "TDatabasePDG.h"
#include "TParticlePDG.h"


#include "Stntuple/obj/TStnDataBlock.hh"
#include "Stntuple/obj/TStnNode.hh"
#include "Stntuple/obj/TStnEvent.hh"

#include "Stntuple/obj/TStnTimeCluster.hh"
#include "Stntuple/obj/TStnTimeClusterBlock.hh"

#include "Stntuple/obj/TStnHelix.hh"
#include "Stntuple/obj/TStnHelixBlock.hh"

#include "Stntuple/obj/TStnTrackSeed.hh"
#include "Stntuple/obj/TStnTrackSeedBlock.hh"

#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Event.h"

#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"

#include "CalorimeterGeom/inc/Calorimeter.hh"

#include "DataProducts/inc/XYZVec.hh"
#include "RecoDataProducts/inc/TimeCluster.hh"
#include "RecoDataProducts/inc/HelixSeed.hh"
#include "RecoDataProducts/inc/HelixHit.hh"
#include "RecoDataProducts/inc/KalSeed.hh"

// #include "RecoDataProducts/inc/AlgorithmIDCollection.hh"
#include "RecoDataProducts/inc/CaloCluster.hh"
#include "RecoDataProducts/inc/StrawHitPositionCollection.hh"

// #include "MCDataProducts/inc/PtrStepPointMCVectorCollection.hh"
#include "MCDataProducts/inc/StrawDigiMC.hh"
#include "MCDataProducts/inc/SimParticle.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "MCDataProducts/inc/StepPointMC.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "MCDataProducts/inc/StrawDigiMCCollection.hh"
#include "TrkDiag/inc/TrkMCTools.hh"

#include "Mu2eUtilities/inc/LsqSums4.hh"
#include "Math/VectorUtil.h"

using namespace ROOT::Math::VectorUtil;

//-----------------------------------------------------------------------------
// assume that the collection name is set, so we could grab it from the event
//-----------------------------------------------------------------------------
int  StntupleInitMu2eHelixBlock(TStnDataBlock* Block, AbsEvent* Evt, int Mode) {

  const mu2e::HelixSeedCollection* list_of_helices(0);
  //  mu2e::AlgorithmIDCollection*     aid_coll    (0);

  char                 helix_module_label[100], helix_description[100]; 
  char                 makeSD_module_label[100];

  int ev_number = Evt->event();
  int rn_number = Evt->run();
  int sr_number = Evt->subRun();

  if (Block->Initialized(ev_number,rn_number,sr_number)) return 0;

  TStnHelixBlock*      cb = (TStnHelixBlock*) Block;
  TStnHelix*           helix(0);

  cb->Clear();
  cb->GetModuleLabel("mu2e::StrawDigiMCCollection",makeSD_module_label);
//-----------------------------------------------------------------------------
// algorithm ID's are determined by the same modules as the helices themselves
//-----------------------------------------------------------------------------
  cb->GetModuleLabel("mu2e::HelixSeedCollection", helix_module_label);
  cb->GetDescription("mu2e::HelixSeedCollection", helix_description );

  art::Handle<mu2e::HelixSeedCollection>   hch;

  if (helix_module_label[0] != 0) {
    if (helix_description[0] == 0) {
      Evt->getByLabel(helix_module_label,hch);
    }
    else {
      Evt->getByLabel(helix_module_label,helix_description,hch);
    }

    if (hch.isValid()) list_of_helices = hch.product();
  }

  const mu2e::StrawDigiMCCollection* mcdigis(0);
  art::Handle<mu2e::StrawDigiMCCollection> mcdH;
  Evt->getByLabel(makeSD_module_label, mcdH);
  mcdigis = mcdH.product();
  
  const mu2e::HelixSeed     *tmpHel(0);
  int                        nhelices(0);
  const mu2e::RobustHelix   *robustHel(0);
  const mu2e::CaloCluster   *cluster(0);
 
  if (list_of_helices) nhelices = list_of_helices->size();
  
  TParticlePDG* part(nullptr);
  TDatabasePDG* pdg_db = TDatabasePDG::Instance();

  static XYZVec zaxis(0.0,0.0,1.0); // unit in z direction

  for (int i=0; i<nhelices; i++) {
    std::vector<int>     hits_simp_id, hits_simp_index, hits_simp_z;
    
    helix                  = cb->NewHelix();
    tmpHel                 = &list_of_helices->at(i);
    cluster                = tmpHel->caloCluster().get();
    robustHel              = &tmpHel->helix();
    if (cluster != 0){
      mu2e::GeomHandle<mu2e::Calorimeter> ch;
      const mu2e::Calorimeter* _calorimeter = ch.get();      
      
      helix->fClusterTime    = cluster->time();
      helix->fClusterEnergy  = cluster->energyDep();
      CLHEP::Hep3Vector gpos = _calorimeter->geomUtil().diskToMu2e(cluster->diskId(),cluster->cog3Vector());
      CLHEP::Hep3Vector tpos = _calorimeter->geomUtil().mu2eToTracker(gpos);
      helix->fClusterX       = tpos.x();
      helix->fClusterY       = tpos.y();
      helix->fClusterZ       = tpos.z();
    } else {
      helix->fClusterTime    = 0; 
      helix->fClusterEnergy  = 0; 
      helix->fClusterX       = 0; 
      helix->fClusterY       = 0; 
      helix->fClusterZ       = 0; 
    }
    
    helix->fHelix        = tmpHel;
    helix->fHelicity     = robustHel->helicity()._value;
    helix->fT0           = tmpHel->t0()._t0;
    helix->fT0Err        = tmpHel->t0()._t0err;     
    helix->fRCent        = robustHel->rcent  ();
    helix->fFCent        = robustHel->fcent  ();     
    helix->fRadius       = robustHel->radius ();
    helix->fLambda       = robustHel->lambda ();     
    helix->fFZ0          = robustHel->fz0    ();
    helix->fD0           = robustHel->rcent  () - robustHel->radius ();
    
    const mu2e::HelixHitCollection* hits = &tmpHel->hits();
    int nhits = hits->size();

    const mu2e::ComboHit*           hit(0);
//-----------------------------------------------------------------------------
// figure out which algorithm was used to reconstruct the helix - helix collections 
// could be merged
// 1: RobustHelixFinder 2: CalHelixFinder
//-----------------------------------------------------------------------------
    if (tmpHel->status().hasAnyProperty(mu2e::TrkFitFlagDetail::CPRHelix)) {
      helix->fAlgorithmID = 0x20002;
    }
    else {
      helix->fAlgorithmID = 0x10001;
    }

    int      nStrawHits(0);
    float    first_hit_z(0), last_hit_z(0);
    
    for (int j=0; j<nhits; ++j) {      //this is a loop is over ComboHits
      hit       = &hits->at(j);

      if(j==0) first_hit_z = hit->pos().z();
      else if(j==nhits-1) last_hit_z = hit->pos().z();
      
      // get the MC truth info
      if (hit->_flag.hasAnyProperty(mu2e::StrawHitFlag::outlier))         continue;

      std::vector<StrawDigiIndex> shids;
      hits->fillStrawDigiIndices(*(Evt),j,shids);
      
      for (size_t k=0; k<shids.size(); ++k) {
      	const mu2e::StrawDigiMC* sdmc = &mcdigis->at(shids[k]);

	const mu2e::StrawGasStep* step (nullptr);
	if (sdmc->wireEndTime(mu2e::StrawEnd::cal) < sdmc->wireEndTime(mu2e::StrawEnd::hv)) {
	  step = sdmc->strawGasStep(mu2e::StrawEnd::cal).get();
	}
	else {
	  step = sdmc->strawGasStep(mu2e::StrawEnd::hv ).get();
	}

      	if (step) {
      	  art::Ptr<mu2e::SimParticle> const& simptr = step->simParticle(); 
      	  int sim_id        = simptr->id().asInt();
      	  float   dz        = step->position().z();// - trackerZ0;
      	  hits_simp_id.push_back   (sim_id); 
      	  hits_simp_index.push_back(shids[k]);
      	  hits_simp_z.push_back(dz);
      	  break;
      	}
      }
      nStrawHits += hit->nStrawHits();
    } 

//-----------------------------------------------------------------------------
// calculate the number of loops made
//-----------------------------------------------------------------------------
    helix->fNLoops = (last_hit_z - first_hit_z)/(fabs(helix->fLambda)*2.*M_PI);

//-----------------------------------------------------------------------------
// find the SimParticle that created the majority of the hits
//-----------------------------------------------------------------------------
    int     max(0), mostvalueindex(-1), mostvalue(-1), id_max(0);
    float   dz_most(1e4);
    if (hits_simp_id.size()>0) mostvalue = hits_simp_id[0];
    
    for (int  k=0; k<(int)hits_simp_id.size(); ++k){
      int co = (int)std::count(hits_simp_id.begin(), hits_simp_id.end(), hits_simp_id[k]);
      if ( (co>0) && (co>=max)) {
	float  dz      = std::fabs(hits_simp_z[k]);
	if (dz < dz_most){
	  dz_most        = dz;
	  max            = co;
	  id_max         = k;
	  mostvalue      = hits_simp_id[k];
	  mostvalueindex = hits_simp_index[k];
	}
      }
    }

    helix->fSimpId1Hits = max;
    helix->fSimpId2Hits = -1;
    helix->fSimpPDGM1   = -1;
    helix->fSimpPDGM2   = -1;
    
    if (hits_simp_id.size()>0) {
      if ( (mostvalueindex != hits_simp_index[id_max]) ){
	printf(">>> ERROR: event %i helix %i no MC found. MostValueindex = %i hits_simp_index[id_max] =%i \n", 
	       Evt->event(), i, mostvalueindex, hits_simp_index[id_max]);
      }
    }
    
    if ( (mostvalueindex<0) || (mostvalueindex >= (int)mcdigis->size()))        {
      printf(">>> ERROR: event %i helix %i no MC found. MostValueindex = %i hits_simp_index[id_max] = %i mcdigis_size =%li \n", 
	     Evt->event(), i, mostvalueindex, hits_simp_index[id_max], mcdigis->size());
      continue;
    }
    
    const mu2e::StrawDigiMC* sdmc = &mcdigis->at(mostvalueindex);

    const mu2e::StrawGasStep* step(nullptr);
    if (sdmc->wireEndTime(mu2e::StrawEnd::cal) < sdmc->wireEndTime(mu2e::StrawEnd::hv)) {
      step = sdmc->strawGasStep(mu2e::StrawEnd::cal).get();
    }
    else {
      step = sdmc->strawGasStep(mu2e::StrawEnd::hv ).get();
    }

    const mu2e::SimParticle* sim(nullptr);

    if (step) {
      art::Ptr<mu2e::SimParticle> const& simptr = step->simParticle(); 
      helix->fSimpPDG1    = simptr->pdgId();
      art::Ptr<mu2e::SimParticle> mother = simptr;
      part   = pdg_db->GetParticle(helix->fSimpPDG1);

      while(mother->hasParent()) mother = mother->parent();
      sim = mother.operator ->();

      helix->fSimpPDGM1   = sim->pdgId();
    
      double   px = step->momvec().x();
      double   py = step->momvec().y();
      double   pz = step->momvec().z();
      double   mass  (-1.);
      double   energy(-1.);
      if (part) {
    	mass   = part->Mass();
    	energy = sqrt(px*px+py*py+pz*pz+mass*mass);
      }
      helix->fMom1.SetPxPyPzE(px,py,pz,energy);

      const CLHEP::Hep3Vector* sp = &simptr->startPosition();
      helix->fOrigin1.SetXYZT(sp->x(),sp->y(),sp->z(),simptr->startGlobalTime());
    }
    
    //look for the second most frequent hit
    if (max != int(hits_simp_id.size()) ){  //nhits){
      int   secondmostvalueindex(-1);
      max     = 0;//reset max
      dz_most = 1e4;

      for (int k=0; k<(int)hits_simp_id.size(); ++k){
	int value = hits_simp_id[k];
	int co = (int)std::count(hits_simp_id.begin(), hits_simp_id.end(), value);
	if ( (co>0) && (co>=max) && (value != mostvalue)) {
	  float  dz      = std::fabs(hits_simp_z[k]);
	  if (dz < dz_most){
	    max                  = co;
	    dz_most              = dz;
	    secondmostvalueindex = hits_simp_index[k];
	  }
	}
      }
      helix->fSimpId2Hits = max;

      if (secondmostvalueindex >= 0) {

       	const mu2e::StrawGasStep* step(nullptr);

      	const mu2e::StrawDigiMC*  sdmc = &mcdigis->at(secondmostvalueindex);
	if (sdmc->wireEndTime(mu2e::StrawEnd::cal) < sdmc->wireEndTime(mu2e::StrawEnd::hv)) {
	  step = sdmc->strawGasStep(mu2e::StrawEnd::cal).get();
	}
	else {
	  step = sdmc->strawGasStep(mu2e::StrawEnd::hv ).get();
	}

	const mu2e::SimParticle* sim (nullptr);

       	if (step) {
	  art::Ptr<mu2e::SimParticle> const& simptr = step->simParticle(); 
	  helix->fSimpPDG2    = simptr->pdgId();
	  art::Ptr<mu2e::SimParticle> mother = simptr;
	  part   = pdg_db->GetParticle(helix->fSimpPDG2);

	  while(mother->hasParent()) mother = mother->parent();
	  sim = mother.operator ->();

	  helix->fSimpPDGM2   = sim->pdgId();
      
	  double   px = simptr->startMomentum().x();
	  double   py = simptr->startMomentum().y();
	  double   pz = simptr->startMomentum().z();

	  double   mass  (-1.), energy(-1.);
	  if (part) {
	    mass   = part->Mass();
	    energy = sqrt(px*px+py*py+pz*pz+mass*mass);
	  }
	  helix->fMom2.SetPxPyPzE(px,py,pz,energy);

	  const CLHEP::Hep3Vector* sp = &simptr->startPosition();
	  helix->fOrigin2.SetXYZT(sp->x(),sp->y(),sp->z(),simptr->startGlobalTime());
	}      
      }
    }
    
    helix->fNComboHits   = tmpHel->hits().size();
    helix->fNHits        = nStrawHits;
    helix->fChi2XYNDof   = robustHel->chi2dXY();
    helix->fChi2PhiZNDof = robustHel->chi2dZPhi();
  }
//-----------------------------------------------------------------------------
// on return set event and run numbers to mark block as initialized
//-----------------------------------------------------------------------------
  cb->f_RunNumber    = rn_number;
  cb->f_EventNumber  = ev_number;
  cb->f_SubrunNumber = sr_number;

  return 0;
}

//-----------------------------------------------------------------------------
Int_t StntupleInitMu2eHelixBlockLinks(TStnDataBlock* Block, AbsEvent* AnEvent, int Mode) 
{

  Int_t  evt, run, srn;

  evt = AnEvent->event();
  run = AnEvent->run();
  srn = AnEvent->subRun();

  if (! Block->Initialized(evt,run,srn)) return -1;

					// do not do initialize links 2nd time

  if (Block->LinksInitialized()) return 0;

  TStnEvent*               ev;
  TStnHelixBlock*          hb;
  TStnHelix*               helix;
  TStnTrackSeedBlock*      tsb;
  TStnTrackSeed*           trkseed;
  TStnTimeClusterBlock*    tcb;
  TStnTimeCluster*         tp;

  const mu2e::HelixSeed*   khelix, *fkhelix;
  const mu2e::KalSeed*     kseed;
  const mu2e::TimeCluster* ktimepeak, *fktimepeak;

  char                     ts_block_name[100], tc_block_name[100];

  ev     = Block->GetEvent();
  hb     = (TStnHelixBlock*) Block;
  int nhelices = hb->NHelices();
//-----------------------------------------------------------------------------
// this is a hack, to be fixed soon
//-----------------------------------------------------------------------------
  hb->GetModuleLabel("TrackSeedBlockName"  ,ts_block_name);
  hb->GetModuleLabel("TimeClusterBlockName",tc_block_name);

  int ntseeds(0), ntpeaks(0);

  tsb    = (TStnTrackSeedBlock*  ) ev->GetDataBlock(ts_block_name);
  tcb    = (TStnTimeClusterBlock*) ev->GetDataBlock(tc_block_name);

  if (tsb) ntseeds  = tsb->NTrackSeeds();
  if (tcb) ntpeaks = tcb->NTimeClusters();

  for (int i=0; i<nhelices; ++i){
    helix  = hb   ->Helix(i);
    khelix = helix->fHelix;
    int      trackseedIndex(-1);
    for (int j=0; j<ntseeds; ++j){
      trkseed = tsb->TrackSeed(j);
      kseed   = trkseed->fTrackSeed;
      fkhelix = kseed->helix().get();
      if (fkhelix == khelix) {
	trackseedIndex = j;
	break;
      }
    }
//-----------------------------------------------------------------------------
// don't print diagnostics if the block name is empty
// non-merged helices and not chozen helices do not have associated seeds
//-----------------------------------------------------------------------------
    // if ((trackseedIndex < 0) && (ts_block_name[0] != 0)) {
    //   printf(">>> ERROR: block %s HelixFinder helix %i -> no TrackSeed associated\n", hb->GetNode()->GetName(),i);
    //   continue;
    // }
    helix->SetTrackSeedIndex(trackseedIndex);

    ktimepeak = khelix->timeCluster().get();
    int      timepeakIndex(-1);
    for (int j=0; j<ntpeaks;++j){
      tp         = tcb->TimeCluster(j);
      fktimepeak = tp->fTimeCluster;
      if (fktimepeak == ktimepeak){
	timepeakIndex = j;
	break;
      }
    }
//-----------------------------------------------------------------------------
// don't print diagnostics if the time cluster block name is empty
// for TimeClusterBlockTpr clusters are not defined
//-----------------------------------------------------------------------------
    // if ((timepeakIndex < 0) && (tc_block_name[0] != 0)) {
    //   printf(">>> ERROR: block %s  %s CalHelixFinder helix %i -> no TimeCluster associated\n", hb->GetNode()->GetName(), tc_block_name, i);//FIXME!
    //   continue;
    // }

    helix->SetTimeClusterIndex(timepeakIndex);
  }
//-----------------------------------------------------------------------------
// mark links as initialized
//-----------------------------------------------------------------------------
  hb->fLinksInitialized = 1;

  return 0;
}

