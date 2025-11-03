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

#include "Offline/MCDataProducts/inc/GenParticle.hh"
#include "Offline/MCDataProducts/inc/EventWeight.hh"

#include "Stntuple/obj/TGenpBlock.hh"
#include "Stntuple/mod/InitGenpBlock.hh"

// #include "Stntuple/mod/THistModule.hh"


#include "messagefacility/MessageLogger/MessageLogger.h"

//-----------------------------------------------------------------------------
// fill GenParticle's data block
//-----------------------------------------------------------------------------
int StntupleInitGenpBlock::InitDataBlock(TStnDataBlock* Block, AbsEvent* AnEvent, int mode) 
{
  //  const char* oname = {"StntupleInitGenpBlock"};

  std::vector<art::Handle<mu2e::GenParticleCollection>> list_of_gp;

  const mu2e::GenParticleCollection*     coll(0);
  const mu2e::GenParticle*               gp(0);
  const art::Provenance*                 prov;
  const art::Handle<mu2e::GenParticleCollection>* handle;

  const int verbose(0);

  art::Selector  selector(art::ProductInstanceNameSelector(""));

  double px, py, pz, mass, energy;
  int    pdg_code, gen_id;

  TGenpBlock* genp_block = (TGenpBlock*) Block;
  genp_block->Clear();
//-----------------------------------------------------------------------------
// MC event weight (by "generate") is not always present
//-----------------------------------------------------------------------------
  art::InputTag ewT("generate");
  art::Handle<mu2e::EventWeight> ewH;

  AnEvent->getByLabel(ewT,ewH);

  if (ewH.isValid()) genp_block->fWeight = ewH->weight();
  else               genp_block->fWeight = 1.;

//-----------------------------------------------------------------------------
// MC event energy (by "generate") is not always present
//-----------------------------------------------------------------------------
  art::InputTag geeT("generate:photon");
  art::Handle<mu2e::GenParticle> geeH;

  AnEvent->getByLabel(geeT,geeH);

  if (geeH.isValid())   genp_block->fGenEnergy = geeH->momentum().e();
  else {
    art::Handle<mu2e::GenParticleCollection> geecH; //newer versions input GenParticleCollection

    AnEvent->getByLabel(geeT,geecH);
    if(geecH.isValid()) genp_block->fGenEnergy = geecH->at(0).momentum().e();
    else                genp_block->fGenEnergy = -1.;
  }

  genp_block->SetGenProcessID(fGenProcessID);
//-----------------------------------------------------------------------------
// loop over existing GENP collections, there could be many of them
//-----------------------------------------------------------------------------
  list_of_gp = AnEvent->getMany<mu2e::GenParticleCollection>(selector);

  TDatabasePDG* pdg_db = TDatabasePDG::Instance();
  TParticlePDG* part;

  for (std::vector<art::Handle<mu2e::GenParticleCollection>> ::const_iterator it = list_of_gp.begin();
       it != list_of_gp.end(); it++) {
    handle = it.operator -> ();

    if (handle->isValid()) {
      coll = handle->product();
      prov = handle->provenance();

//       printf("moduleLabel = %-20s, producedClassname = %-30s, productInstanceName = %-20s\n",
// 	     prov->moduleLabel().data(),
// 	     prov->producedClassName().data(),
// 	     prov->productInstanceName().data());

      if (fGenpCollTag.label().data()[0] != 0) {
//-----------------------------------------------------------------------------
// module label is defined, assume there is only one generator module to save
//-----------------------------------------------------------------------------
	if (strcmp(fGenpCollTag.label().data(),prov->moduleLabel().data()) != 0) continue;
      }

      for (std::vector<mu2e::GenParticle>::const_iterator ip = coll->begin();
	   ip != coll->end(); ip++) {
	gp       = ip.operator -> ();
	gen_id   = (int) gp->generatorId().id();
	pdg_code = (int) gp->pdgId();
        if(verbose > 1) printf("InitGenpBlock::%s: Gen particle: PDG = %5i, gen_id = %3i\n", __func__, pdg_code, gen_id);
//-----------------------------------------------------------------------------
// conditionally. store only particles corresponding to the requested process
//-----------------------------------------------------------------------------
	if ((fGenProcessID > 0) && (gen_id   != fGenProcessID)) continue;
	if ((fPdgID       != 0) && (pdg_code != fPdgID       )) continue;
        if ((fPdgID == 0) && (pdg_code == 2212))                continue; // ignore POT unless requested
        if(verbose > 0) printf("InitGenpBlock::%s: Accepted gen particle: PDG = %5i, gen_id = %3i\n", __func__, pdg_code, gen_id);

	part     = pdg_db->GetParticle(pdg_code);

	px     = gp->momentum().x();
	py     = gp->momentum().y();
	pz     = gp->momentum().z();
	if (part != nullptr) mass   = part->Mass();
	else                 mass   = 0.;// FIXME!
	energy = sqrt(px*px+py*py+pz*pz+mass*mass);

	genp_block->NewParticle(pdg_code,
				gen_id,
				-1, -1, -1, -1,
				px, py, pz, energy,
				gp->position().x(),
				gp->position().y(),
				gp->position().z(),
				gp->time(),
				gp->properTime());
      }
    }
    else {
      printf(">>> ERROR in StntupleInitMu2eGenpBlock: failed to locate collection %s. BAIL OUT.\n",
	     fGenpCollTag.encode().data());
      return -1;
    }
  }

  return 0;
}


