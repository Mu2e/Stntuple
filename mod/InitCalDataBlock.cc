///////////////////////////////////////////////////////////////////////////////
// 2014-01-26 P.Murat Mu2e version of TCalDataBlock initialization
///////////////////////////////////////////////////////////////////////////////
#include <cstdio>
#include "TROOT.h"
#include "TFolder.h"
#include "TLorentzVector.h"

#include "art/Framework/Principal/Handle.h"

#include "Offline/GeometryService/inc/GeometryService.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"

#include "Offline/CalorimeterGeom/inc/Calorimeter.hh"
#include "Offline/CalorimeterGeom/inc/DiskCalorimeter.hh"

#include "Offline/RecoDataProducts/inc/CaloHit.hh"

#include "Stntuple/obj/TCalDataBlock.hh"
#include "Stntuple/mod/InitCalDataBlock.hh"

//-----------------------------------------------------------------------------
namespace stntuple {

//-----------------------------------------------------------------------------
InitCalDataBlock::InitCalDataBlock() {
}

//-----------------------------------------------------------------------------
int InitCalDataBlock::InitDataBlock(TStnDataBlock* Block, AbsEvent* AnEvent, int Mode) {
  // initialize CAL data block with the `event' data

  int ev_number = AnEvent->event();
  int rn_number = AnEvent->run();

  if (Block->Initialized(ev_number,rn_number)) return 0;

  TCalDataBlock* data = (TCalDataBlock*) Block;
  data->Clear();

       // Get handles to calorimeter crystal hits

  const mu2e::CaloHitCollection* list_of_hits(nullptr);

  if (not fCaloHitCollTag.empty()) {
    art::Handle<mu2e::CaloHitCollection> chch;
    bool ok = AnEvent->getByLabel(fCaloHitCollTag,chch);
    if (ok) {
      list_of_hits = chch.product();
    }
  }

  if (list_of_hits == nullptr) {
    printf(" >>> ERROR in stntuple::InitCalDataBlock: no list_of_hits. BAIL OUT\n");
    return -1;
  }

  int nhits = list_of_hits->size();

  TCalHitData*   hit;

  // reminder: data->fNHits is set to 0 by TCalDataBlock::Clear(), should be this way

  for (int i=0; i<nhits; i++) {
    const mu2e::CaloHit* calo_hit = &list_of_hits->at(i);
    hit      = data->NewCalHitData();

    hit->Set(calo_hit->crystalID(),
	     calo_hit->nSiPMs(),
	     calo_hit->time(),
	     calo_hit->energyDep());
  }
//-----------------------------------------------------------------------------
// store geometry data (this is, obviously, a kludge)
//-----------------------------------------------------------------------------
  art::ServiceHandle<mu2e::GeometryService> geom;
  mu2e::GeomHandle<mu2e::DiskCalorimeter>   dc;

  const mu2e::DiskCalorimeter*              cal;
  const mu2e::Disk*                         disk;

  cal = dc.get();

  data->fNDisks = cal->nDisks();
  for (int i=0; i<data->fNDisks; i++) {
    disk = &cal->disk(i);
    data->fNCrystals[i] = disk->nCrystals();
    data->fRMin     [i] = disk->diskInfo().innerEnvelopeR();
    data->fRMax     [i] = disk->diskInfo().outerEnvelopeR();
    data->fZ0       [i] = disk->diskInfo().origin().z();
  }

  data->fCrystalSize = cal->G4Info().get<double>("crystalXYLength")/2.; // crystalHalfTrans();

				        // also a dummy line
  data->fMinFraction      = 1.0;
  data->fWrapperThickness = cal->G4Info().get<double>("wrapperThickness");      // wrapperThickness();
  data->fShellThickness   = -1.; // 2026-08-08: currently undefined // cal->G4Info().get<double>("crystalFrameThickness"); // caseThickness  ();

					// on return set event and run numbers
					// to mark block as initialized
  data->f_RunNumber   = rn_number;
  data->f_EventNumber = ev_number;

  return 0;
}

//-----------------------------------------------------------------------------
int InitCalDataBlock::ResolveLinks(TStnDataBlock* Block, AbsEvent* AnEvent, int Mode)  {
  return 0;
}

}
