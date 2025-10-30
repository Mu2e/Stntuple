////////////////////////////////////////////////////////////////////////////////

#include "Stntuple/print/TAnaDump.hh"
#include "Stntuple/print/Stntuple_print_functions.hh"

#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/Selector.h"

#include "Offline/RecoDataProducts/inc/StrawHit.hh"
#include "Offline/RecoDataProducts/inc/CosmicTrackSeed.hh"
#include "Offline/RecoDataProducts/inc/HelixSeed.hh"
#include "Offline/RecoDataProducts/inc/ComboHit.hh"
#include "Offline/RecoDataProducts/inc/TrkStrawHitSeed.hh"
#include "Offline/RecoDataProducts/inc/TimeCluster.hh"
#include "print/TAnaDump.hh"

using namespace std;
//-----------------------------------------------------------------------------
void TAnaDump::printCosmicTrackSeed(const mu2e::CosmicTrackSeed* Cts        ,
                                    const char*                  Opt        , 
                                    const char*                  SHCollTag  ,
                                    const char*                  SdmcCollTag) {
  TString opt = Opt;
  
  if ((opt == "") || (opt == "banner")) {
    printf("-----------------------------------------------------------------------------------------------------------\n");
    printf("  TrkID  Address   Npt    Status      T0        T0Err     A0        A1       B0        B1    Chi2X    Chi2Y\n");
    printf("-----------------------------------------------------------------------------------------------------------\n");
  }

  int    nhits   = Cts->hits().size();
  
  if ((opt == "") || (opt.Index("data") >= 0)) {
    double t0     = Cts->t0()._t0;
    double t0err  = Cts->t0()._t0err;
    mu2e::TrkFitFlag status = Cts->status();

    double a0     = Cts->_track.FitParams.A0;
    double a1     = Cts->_track.FitParams.A1;
    double b0     = Cts->_track.FitParams.B0;
    double b1     = Cts->_track.FitParams.B1;
    double chi2x  = Cts->_track.Diag.FinalChiX;
    double chi2z  = Cts->_track.Diag.FinalChiY;
    
    printf(" %5i 0x%8p %4i 0x%08x %8.3f %8.3f %8.3f %8.4f %8.3f %8.4f %8.3f %8.3f\n",
           -1, (void*) Cts, nhits, *((int*) &status),
           t0,t0err,a0,a1,b0,b1,chi2x,chi2z);
  }

  if ((opt == "") || (opt.Index("hits") >= 0)) {

    const mu2e::ComboHit* hit(0), *hit_0(0);
    const mu2e::StrawGasStep* step(0);

    art::Handle<mu2e::ComboHitCollection> shcolH;
    const mu2e::ComboHitCollection*       shcol(0);
    fEvent->getByLabel<mu2e::ComboHitCollection>(art::InputTag(SHCollTag),shcolH);
    if (shcolH.isValid()) shcol = shcolH.product();
    else {
      printf("ERROR in TAnaDump::%s: no ComboHitCollection with tag=%s, BAIL OUT\n",__func__,SHCollTag);
      return;
    }

    art::InputTag tag; // sdmccT(fStrawDigiMCCollTag.Data());
    art::Handle<mu2e::StrawDigiMCCollection> sdmccH;
    const mu2e::StrawDigiMCCollection* sdmcc(nullptr);

    if      (SdmcCollTag != nullptr     ) tag = SdmcCollTag;
    else if (strcmp(SdmcCollTag,"") != 0) tag = SdmcCollTag;
    else {
      printf("ERROR in TAnaDump::%s: no SdmcCollTag specified, BAIL OUT\n",__func__);
      return;
    }
    fEvent->getByLabel(tag,sdmccH);
    if (sdmccH.isValid()) sdmcc = sdmccH.product();
    else {
      printf("ERROR in TAnaDump::%s: no StrawDigiMCCollection with tag=%s,",__func__,tag.encode().data());
      print_sdmc_colls();
      return;
    }
//-----------------------------------------------------------------------------
// StrawDigiMCCollection is present
//-----------------------------------------------------------------------------
    hit_0    = &shcol->at(0);
    int      loc(-1);
    int banner_printed(0);
    for (int i=0; i<nhits; ++i){
      const mu2e::ComboHit* hit_seed = &Cts->hits().at(i);
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
void TAnaDump::printCosmicTrackSeedCollection(const char* CTSCollTag        ,
                                              int         hitOpt            ,
                                              const char* StrawHitCollTag   ,
                                              const char* StrawDigiMCCollTag) {

  art::Handle<mu2e::CosmicTrackSeedCollection> ctsch;
  
  bool ok = fEvent->getByLabel(art::InputTag(CTSCollTag),ctsch);
//-----------------------------------------------------------------------------
// make sure collection exists
//-----------------------------------------------------------------------------
  if (! ok) {
    printf("ERROR in TAnaDump::%s: no CosmicTrackSeedCollection tag=%s\n",__func__,CTSCollTag);
    print_kalseed_colls();
    return;
  }

  art::Handle<mu2e::StrawDigiMCCollection> sdmccH;
  fEvent->getByLabel(art::InputTag(StrawDigiMCCollTag),sdmccH);

  _mcdigis = nullptr;
  if (StrawDigiMCCollTag != nullptr) {
    if (sdmccH.isValid()) {
      _mcdigis     = sdmccH.product();
      fSdmcCollTag = StrawDigiMCCollTag;
    }
    else {
      printf("ERROR in TAnaDump::%s: no StrawDigiMCCollection tag=%s, available are\n",__func__,StrawDigiMCCollTag);
      print_sdmc_colls();
    }
  }
    
  mu2e::CosmicTrackSeedCollection*  list_of_cts;
  list_of_cts = (mu2e::CosmicTrackSeedCollection*) &(*ctsch);

  int ncts = list_of_cts->size();

  //  const mu2e::CosmicTrackSeed *cts;

  int banner_printed = 0;
  for (int i=0; i<ncts; i++) {
    const mu2e::CosmicTrackSeed* cts = &list_of_cts->at(i);
    if (banner_printed == 0) {
      printCosmicTrackSeed(cts,"banner");
      banner_printed = 1;
    }
    printCosmicTrackSeed(cts,"data");
    if (hitOpt > 0) printCosmicTrackSeed(cts,"hits",StrawHitCollTag,StrawDigiMCCollTag);
  }
}

