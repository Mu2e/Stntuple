////////////////////////////////////////////////////////////////////////////////

#include "Stntuple/print/TAnaDump.hh"
#include "Stntuple/print/Stntuple_print_functions.hh"

#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/Selector.h"

#include "Offline/RecoDataProducts/inc/StrawHit.hh"
#include "Offline/RecoDataProducts/inc/KalSeed.hh"
#include "Offline/RecoDataProducts/inc/HelixSeed.hh"
#include "Offline/RecoDataProducts/inc/ComboHit.hh"
#include "Offline/RecoDataProducts/inc/TrkStrawHitSeed.hh"
#include "Offline/RecoDataProducts/inc/TimeCluster.hh"
#include "Offline/RecoDataProducts/inc/TrkFitFlag.hh"

using namespace std;

//-----------------------------------------------------------------------------
// KKLine is unique
//-----------------------------------------------------------------------------
void TAnaDump::printKalSeed_Line(const mu2e::KalSeed* KalSeed      ,
                                 const char*          Opt          , 
                                 const char*          ShCollTag    ,
                                 const char*          SdmcCollTag) {
  TString opt = Opt;

  if (not KalSeed->status().hasAnyProperty(mu2e::TrkFitFlag::KKLine)) return;
  
  if ((opt == "") || (opt == "banner")) {
    printf("--------------------------------------------------------------------------------");
    printf("----------------------------------------------------------------------------\n");
    printf("  TrkID       Address    N  chi2/ndof    T0      T0err       X           Y         Z       ");
    printf("    Nx         Ny        Nz     tmin     tmax    fmin      fmax\n");
    printf("--------------------------------------------------------------------------------");
    printf("----------------------------------------------------------------------------\n");
  }

  int    nhits   = KalSeed->hits().size();
  if ((opt == "") || (opt.Index("data") >= 0)) {
    
    double t0     = KalSeed->t0()._t0;
    double t0err  = KalSeed->t0()._t0err;
    double ndof   = KalSeed->nDOF();
    // auto status   = KalSeed->status();
    float chi2    = KalSeed->chisquared()/(ndof+1.e-12);
    
    for (const mu2e::KalSegment& kalSeg : KalSeed->segments()) {

      KinKal::VEC3 pos   = kalSeg.kinematicLine().pos0();
      KinKal::VEC3 dir   = kalSeg.kinematicLine().direction();
      double tmin        = kalSeg.tmin();
      double tmax        = kalSeg.tmax();
      double fmin        = kalSeg.fmin();
      double fmax        = kalSeg.fmax();

      // const mu2e::CaloCluster*cluster = KalSeed->caloCluster().get();
      // double clusterEnergy(-1);
      // if (cluster != 0) clusterEnergy = cluster->energyDep();
      printf("%5i %16p %3i %8.3f %8.3f %8.5f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %8.2f %8.2f %8.2f %8.2f\n",
	     -1,
	     static_cast<const void*>(KalSeed),
	     nhits,chi2,
             t0, t0err, pos.x(),pos.y(),pos.z(), dir.x(),dir.y(),dir.z(),
             tmin,tmax,fmin, fmax);
    }
  }

  if ((opt == "") || (opt.Index("hits") >= 0)) {
  
    const mu2e::ComboHit* hit(0), *hit_0(0);
    const mu2e::StrawGasStep* step(0);

    art::Handle<mu2e::ComboHitCollection> shcolH;
    const mu2e::ComboHitCollection*       shcol(0);
    fEvent->getByLabel<mu2e::ComboHitCollection>(art::InputTag(ShCollTag),shcolH);
    if (shcolH.isValid()) shcol = shcolH.product();
    else {
      printf("ERROR in TAnaDump::%s: no ComboHitCollection with tag=%s, BAIL OUT\n",__func__,ShCollTag);
      return;
    }

    art::InputTag tag; // sdmccT(fSdmcCollTag.Data());
    art::Handle<mu2e::StrawDigiMCCollection> sdmccH;
    const mu2e::StrawDigiMCCollection* sdmcc(nullptr);

    if      (SdmcCollTag != nullptr) tag = SdmcCollTag;
    else if (fSdmcCollTag       != ""     ) tag = fSdmcCollTag;
    else {
      printf("ERROR in TAnaDump::%s: no SdmcCollTag specified, BAIL OUT\n",__func__);
      return;
    }
    fEvent->getByLabel(tag,sdmccH);
    if (sdmccH.isValid()) sdmcc = sdmccH.product();
    else {
      printf("ERROR in TAnaDump::%s: no StrawDigiMCCollection with tag=%s,",__func__,tag.encode().data());
      printf(" available collections are:\n");
 
      vector<art::Handle<mu2e::StrawDigiMCCollection>> list;
      const  art::Handle<mu2e::StrawDigiMCCollection>*  handle;
      const art::Provenance*                       prov;
      
      list  = fEvent->getMany<mu2e::StrawDigiMCCollection>();

      for (auto it = list.begin(); it != list.end(); it++) {
	handle = it.operator -> ();
	if (handle->isValid()) {
	  prov = handle->provenance();
	
	  printf("moduleLabel: %-20s, productInstanceName: %-20s, processName: %-30s nHelices: %3li\n" ,
		 prov->moduleLabel().data(),
		 prov->productInstanceName().data(),
		 prov->processName().data(),
		 handle->product()->size()
		 );
	}
      }
      return;
    }
//-----------------------------------------------------------------------------
// StrawDigiMCCollection is present
//-----------------------------------------------------------------------------
    hit_0    = &shcol->at(0);
    int      loc(-1);
    int banner_printed(0);
    for (int i=0; i<nhits; ++i){
      const mu2e::TrkStrawHitSeed* hit_seed = &KalSeed->hits().at(i);
      int  hitIndex  = int(hit_seed->index());
      hit            = &shcol->at(hitIndex);
      loc            = hit - hit_0;
//-----------------------------------------------------------------------------
// fake hit flag - FIXME
//-----------------------------------------------------------------------------
      int straw_hit_flag = hit_seed->flag().hasAllProperties(mu2e::StrawHitFlagDetail::active);

      if (sdmcc) {
	const mu2e::StrawDigiMC* sdmc = &sdmcc->at(loc);
	step = nullptr;
	if (sdmc->wireEndTime(mu2e::StrawEnd::cal) < sdmc->wireEndTime(mu2e::StrawEnd::hv)) {
	  step = sdmc->strawGasStep(mu2e::StrawEnd::cal).get();
	}
	else {
	  step = sdmc->strawGasStep(mu2e::StrawEnd::hv ).get();
	}
      }

      if (banner_printed == 0){
	printComboHit(hit, step, "banner", -1, straw_hit_flag);
	banner_printed = 1;
      } 
      else {
	printComboHit(hit, step, "data"  , -1, straw_hit_flag);
      }
    }
  }
}

//-----------------------------------------------------------------------------
// KKLine is unique
//-----------------------------------------------------------------------------
void TAnaDump::printKalSeed(const mu2e::KalSeed* KalSeed      ,
			    const char*          Opt          , 
			    const char*          ShCollTag    ,
			    const char*          SdmcCollTag) {
  TString opt = Opt;

  if (KalSeed->status().hasAnyProperty(mu2e::TrkFitFlag::KKLine)) {
    printKalSeed_Line(KalSeed,Opt,ShCollTag,SdmcCollTag);
    return;
  }
  
  if ((opt == "") || (opt == "banner")) {
    printf("------------------------------------------------------------------------------------");
    printf("----------------------------------------------------------------------------\n");
    printf("  TrkID       Address    N  Q       P      pT        T0      T0err      fmin      fmax");
    printf("       D0        Z0   Phi0   TanDip    radius      Ecl      chi2   FitCon  \n");
    printf("------------------------------------------------------------------------------------");
    printf("----------------------------------------------------------------------------\n");
  }

  if ((opt == "") || (opt.Index("data") >= 0)) {
    int    nhits   = KalSeed->hits().size();
    
    double t0     = KalSeed->t0()._t0;
    double t0err  = KalSeed->t0()._t0err;

    for (const mu2e::KalSegment& kalSeg : KalSeed->segments()) {

      double d0     = kalSeg.helix().d0();
      double z0     = kalSeg.helix().z0();
      double phi0   = kalSeg.helix().phi0();
      double tandip = kalSeg.helix().tanDip();
      double mm2MeV = 3/10.;
      double mom    = kalSeg.mom();
      double pt     = mom*cos(atan(tandip));
      double radius = pt/mm2MeV;
      double q      = kalSeg.centralHelix().charge();


      const mu2e::CaloCluster*cluster = KalSeed->caloCluster().get();
      double clusterEnergy(-1);
      if (cluster != 0) clusterEnergy = cluster->energyDep();
      printf("%5i %16p %3i %2.0f %8.3f %8.5f %7.3f %6.3f %9.3f %9.3f",
	     -1,
	     static_cast<const void*>(KalSeed),
	     nhits,q,
	     mom, pt, t0, t0err, kalSeg.fmin(), kalSeg.fmax() );

      float chi2    = KalSeed->chisquared()/double(nhits - 5.);
      float fitCons = KalSeed->fitConsistency();
      
      printf(" %8.3f %9.3f %6.3f %8.4f %10.4f %8.3f %8.3f %10.3e\n",
	     d0,z0,phi0,tandip,radius,clusterEnergy,chi2,fitCons);
    }
  }

  if ((opt == "") || (opt.Index("hits") >= 0)) {
    int nsh = KalSeed->hits().size();

    const mu2e::ComboHit* hit(0), *hit_0(0);
    const mu2e::StrawGasStep* step(0);

    art::Handle<mu2e::ComboHitCollection> shcolH;
    const mu2e::ComboHitCollection*       shcol(0);
    fEvent->getByLabel<mu2e::ComboHitCollection>(art::InputTag(ShCollTag),shcolH);
    if (shcolH.isValid()) shcol = shcolH.product();
    else {
      printf("ERROR in TAnaDump::printTrackSeed: no ComboHitCollection with tag=%s, BAIL OUT\n",ShCollTag);
      return;
    }

    art::InputTag tag; // sdmccT(fSdmcCollTag.Data());
    art::Handle<mu2e::StrawDigiMCCollection> sdmccH;
    const mu2e::StrawDigiMCCollection* sdmcc(nullptr);

    if      (SdmcCollTag != nullptr) tag = SdmcCollTag;
    else if (fSdmcCollTag       != ""     ) tag = fSdmcCollTag;
    else {
      printf("ERROR in TAnaDump::printTrackSeed: no SdmcCollTag specified, BAIL OUT\n");
      return;
    }
    fEvent->getByLabel(tag,sdmccH);
    if (sdmccH.isValid()) sdmcc = sdmccH.product();
    else {
      printf("ERROR in TAnaDump::printTrackSeed: no StrawDigiMCCollection with tag=%s,",tag.encode().data());
      printf(" available collections are:\n");
 
      vector<art::Handle<mu2e::StrawDigiMCCollection>> list;
      const  art::Handle<mu2e::StrawDigiMCCollection>*  handle;
      const art::Provenance*                       prov;
      
      list  = fEvent->getMany<mu2e::StrawDigiMCCollection>();

      for (auto it = list.begin(); it != list.end(); it++) {
	handle = it.operator -> ();
	if (handle->isValid()) {
	  prov = handle->provenance();
	
	  printf("moduleLabel: %-20s, productInstanceName: %-20s, processName: %-30s nHelices: %3li\n" ,
		 prov->moduleLabel().data(),
		 prov->productInstanceName().data(),
		 prov->processName().data(),
		 handle->product()->size()
		 );
	}
      }
      return;
    }
//-----------------------------------------------------------------------------
// StrawDigiMCCollection is present
//-----------------------------------------------------------------------------
    hit_0    = &shcol->at(0);
    int      loc(-1);
    int banner_printed(0);
    for (int i=0; i<nsh; ++i){
      const mu2e::TrkStrawHitSeed* hit_seed = &KalSeed->hits().at(i);
      int  hitIndex  = int(hit_seed->index());
      hit            = &shcol->at(hitIndex);
      loc            = hit - hit_0;
//-----------------------------------------------------------------------------
// fake hit flag - FIXME
//-----------------------------------------------------------------------------
      int straw_hit_flag = hit_seed->flag().hasAllProperties(mu2e::StrawHitFlagDetail::active);

      if (sdmcc) {
	const mu2e::StrawDigiMC* sdmc = &sdmcc->at(loc);
	step = nullptr;
	if (sdmc->wireEndTime(mu2e::StrawEnd::cal) < sdmc->wireEndTime(mu2e::StrawEnd::hv)) {
	  step = sdmc->strawGasStep(mu2e::StrawEnd::cal).get();
	}
	else {
	  step = sdmc->strawGasStep(mu2e::StrawEnd::hv ).get();
	}
      }

      if (banner_printed == 0){
	printComboHit(hit, step, "banner", -1, straw_hit_flag);
	banner_printed = 1;
      } 
      else {
	printComboHit(hit, step, "data"  , -1, straw_hit_flag);
      }
    }
  }
}

//-----------------------------------------------------------------------------
// coll tag format: ModuleLabel[:CollName[:ProcessName]]
//-----------------------------------------------------------------------------
void TAnaDump::printKalSeedCollection(const char* KalSeedCollTag    ,
				      int         hitOpt            ,
				      const char* ShCollTag   ,
				      const char* SdmcCollTag) {

  art::Handle<mu2e::KalSeedCollection> kseedHandle;
  
  bool ok = fEvent->getByLabel(art::InputTag(KalSeedCollTag),kseedHandle);
//-----------------------------------------------------------------------------
// make sure collection exists
//-----------------------------------------------------------------------------
  if (! ok) {
    printf("ERROR in TAnaDump::printKalSeedCollection: no KalSeedCollection tag=%s\n",KalSeedCollTag);
    print_kalseed_colls();
    return;
  }

  art::Handle<mu2e::StrawDigiMCCollection> sdmccH;
  fEvent->getByLabel(art::InputTag(SdmcCollTag),sdmccH);

  if (sdmccH.isValid()) {
    _mcdigis     = sdmccH.product();
    fSdmcCollTag = SdmcCollTag;
  }
  else {
    printf("ERROR in TAnaDump::printKalSeedCollection: no StrawDigiMCCollection tag=%s, available are\n",SdmcCollTag);
    print_sdmc_colls();
    _mcdigis = nullptr;
  }

  mu2e::KalSeedCollection*  list_of_kseeds;
  list_of_kseeds = (mu2e::KalSeedCollection*) &(*kseedHandle);

  int nks = list_of_kseeds->size();

  const mu2e::KalSeed *ks;

  int banner_printed = 0;
  for (int i=0; i<nks; i++) {
    ks = &list_of_kseeds->at(i);
    if (banner_printed == 0) {
      printKalSeed(ks,"banner");
      banner_printed = 1;
    }
    printKalSeed(ks,"data");
    if (hitOpt > 0) printKalSeed(ks,"hits",ShCollTag,SdmcCollTag);

  }
}

