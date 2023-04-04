////////////////////////////////////////////////////////////////////////////////
#include "Stntuple/print/TAnaDump.hh"
#include "Stntuple/print/Stntuple_print_functions.hh"

#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/Selector.h"

#include "Offline/RecoDataProducts/inc/StrawHit.hh"
#include "Offline/RecoDataProducts/inc/ComboHit.hh"
#include "Offline/RecoDataProducts/inc/TrkStrawHitSeed.hh"

using namespace std;

//-----------------------------------------------------------------------------
void TAnaDump::printComboHit(const mu2e::ComboHit* Hit, const mu2e::StrawGasStep* Step, const char* Opt, int IHit, int Flags) {
  TString opt = Opt;
  opt.ToLower();

  if ((opt == "") || (opt.Index("banner") >= 0)) {
    printf("#-------------------------------------------------------------------------------------------");
    printf("--------------------------------------------------------------------------------------------\n");
    printf("#   I nsh   SID   Flags  Pln:Pnl:Lay:Str      X        Y        Z       Time      eDep   End");
    printf("  DrTime  PrTime  TRes    WDist     WRes        PDG     PDG(M) GenID simID       p        pz\n");
    printf("#-------------------------------------------------------------------------------------------");
    printf("--------------------------------------------------------------------------------------------\n");
  }

  if (opt == "banner") return;

  const mu2e::SimParticle * sim (0);
    
  int      pdg_id(-1), mother_pdg_id(-1), generator_id(-1), sim_id(-1);
  double   mc_mom(-1.);
  double   mc_mom_z(-1.);

  mu2e::GenId gen_id;

  if (Step) {
    art::Ptr<mu2e::SimParticle> const& simptr = Step->simParticle(); 
    art::Ptr<mu2e::SimParticle> mother        = simptr;

    while(mother->hasParent()) mother = mother->parent();

    sim           = mother.operator ->();

    pdg_id        = simptr->pdgId();
    mother_pdg_id = sim->pdgId();

    if (simptr->fromGenerator()) generator_id = simptr->genParticle()->generatorId().id();
    else                         generator_id = -1;

    sim_id        = simptr->id().asInt();
    mc_mom        = Step->momvec().mag();
    mc_mom_z      = Step->momvec().z();
  }
    
  if ((opt == "") || (opt.Index("data") >= 0)) {
    if (IHit  >= 0) printf("%5i " ,IHit);
    else            printf("      ");

    printf("%3i ",Hit->nStrawHits());

    printf("%5u",Hit->strawId().asUint16());

    if (Flags >= 0) printf(" %08x",Flags);
    else            printf("        ");

    printf(" %3i %3i %3i %3i  %8.3f %8.3f %9.3f %8.3f %8.5f  %3i %7.2f %7.2f %5.2f %8.3f %8.3f %10i %10i %5i %5i %8.3f %8.3f\n",
	   Hit->strawId().plane(),
	   Hit->strawId().panel(),
	   Hit->strawId().layer(),
	   Hit->strawId().straw(),
	   Hit->pos().x(),Hit->pos().y(),Hit->pos().z(),
	   Hit->time(),
	   Hit->energyDep(),

	   (int) Hit->endTime(),
	   Hit->driftTime(),
	   Hit->propTime(),
	   Hit->transRes(),
	   Hit->wireDist(),
	   Hit->wireRes(),

	   pdg_id,
	   mother_pdg_id,
	   generator_id,
	   sim_id,
	   mc_mom,
	   mc_mom_z);
  }
}

//-----------------------------------------------------------------------------
// 
//-----------------------------------------------------------------------------
void TAnaDump::printComboHitCollection(const char* StrawHitCollTag   , 
				       const char* StrawDigiMCCollTag,
				       double TMin, double TMax) {

  //  const char* oname = "TAnaDump::printComboHitCollection";
//-----------------------------------------------------------------------------
// get straw hits
//-----------------------------------------------------------------------------
  art::Handle<mu2e::ComboHitCollection> shcH;
  const mu2e::ComboHitCollection* shc(nullptr);
  fEvent->getByLabel(StrawHitCollTag,shcH);

  if (shcH.isValid()) shc = shcH.product();
  else {
    printf("ERROR: cant find StrawHitCollection tag=%s, EXIT\n",StrawHitCollTag);
    print_sh_colls();
    return;
  }

  TString sdmc_tag = StrawDigiMCCollTag;
  if (sdmc_tag == "") sdmc_tag = fStrawDigiMCCollTag;

  art::Handle<mu2e::StrawDigiMCCollection> mcdH;
  fEvent->getByLabel<mu2e::StrawDigiMCCollection>(sdmc_tag.Data(),mcdH);
  const mu2e::StrawDigiMCCollection*  mcdigis(nullptr);
  if (mcdH.isValid())   mcdigis = mcdH.product();
  else {
    printf("ERROR: cant find StrawDigiMCCollection tag=%s, EXIT\n",sdmc_tag.Data());
    print_sdmc_colls();
    return;
  }

  int nhits = shc->size();

  const mu2e::ComboHit* hit;

  int flags;
  
  const mu2e::ComboHit* hit0 = &shc->at(0);
 
  int banner_printed = 0;
  for (int i=0; i<nhits; i++) {
    hit         = &shc->at(i);
    size_t ish  = hit-hit0;
    vector<StrawDigiIndex> shids;
    shc->fillStrawDigiIndices(*(fEvent),ish,shids);

    const mu2e::StrawDigiMC* mcdigi = &mcdigis->at(shids[0]);
    const mu2e::StrawGasStep* step (nullptr);

    if (mcdigi->wireEndTime(mu2e::StrawEnd::cal) < mcdigi->wireEndTime(mu2e::StrawEnd::hv)) {
      step = mcdigi->strawGasStep(mu2e::StrawEnd::cal).get();
    }
    else {
      step = mcdigi->strawGasStep(mu2e::StrawEnd::hv ).get();
    }

    flags = *((int*) &hit->flag());
    if (banner_printed == 0) {
      printComboHit(hit, step, "banner");
      banner_printed = 1;
    }
    if ((hit->time() >= TMin) && (hit->time() <= TMax)) {
      printComboHit(hit, step, "data", i, flags);
    }
  }
 
}
