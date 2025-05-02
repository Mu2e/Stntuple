///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
#include "TLorentzVector.h"
#include "TDatabasePDG.h"
#include "TParticlePDG.h"
#include <vector>

#include "messagefacility/MessageLogger/MessageLogger.h"

#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Provenance.h"
#include "canvas/Utilities/InputTag.h"

#include "Stntuple/obj/TStnHeaderBlock.hh"
#include "Stntuple/mod/InitHeaderBlock.hh"
// #include <Stntuple/mod/StntupleUtilities.hh>

#include "Offline/RecoDataProducts/inc/StrawHit.hh"
#include "Offline/RecoDataProducts/inc/ComboHit.hh"

// #include "Stntuple/mod/THistModule.hh"
#include "Offline/MCDataProducts/inc/ProtonBunchIntensity.hh"
#include "Offline/MCDataProducts/inc/EventWeight.hh"

void stntuple_get_version(char* ver, char* test);

//-----------------------------------------------------------------------------
// fill header  block
//-----------------------------------------------------------------------------
namespace stntuple {

int InitHeaderBlock::InitDataBlock(TStnDataBlock* Block, AbsEvent* AnEvent, int Mode) {
  const char oname []  = {"stntuple::InitHeaderBlock::InitDataBlock"};
  const int verbose(0);

  TStnHeaderBlock* data = (TStnHeaderBlock*) Block;

  data->fEventNumber   = AnEvent->event();
  data->fRunNumber     = AnEvent->run();
  data->fSectionNumber = AnEvent->subRun();
//-----------------------------------------------------------------------------
// MC is supposed to have run numbers < 100000
//-----------------------------------------------------------------------------
  if (data->fRunNumber < 100000) data->fMcFlag = 1;
  else                           data->fMcFlag = 0;

  data->fNTracks       = -1;
//-----------------------------------------------------------------------------
// instantaneous luminosity
//-----------------------------------------------------------------------------
  data->fInstLum       = -1.;
  data->fMeanLum       = -1.;

  art::Handle<mu2e::ProtonBunchIntensity> pbiHandle;
  AnEvent->getByLabel("PBISim", pbiHandle);

  if (!pbiHandle.isValid()) {
    mf::LogWarning(oname) << " WARNING: no ProtonBunchIntensity objects found\n";
  }
  else {
//-----------------------------------------------------------------------------
// Loop over ProtonBunchIntensity objects
//-----------------------------------------------------------------------------
    const mu2e::ProtonBunchIntensity* pbi= pbiHandle.product();
    data->fInstLum = pbi->intensity();
    data->fMeanLum = -1. ;
  }
//-----------------------------------------------------------------------------
//  STNTUPLE version
//-----------------------------------------------------------------------------
  char  ver[200], text[200];
  stntuple_get_version(ver,text);
  data->fStnVersion = ver;
//-----------------------------------------------------------------------------
// number of straw hits
//-----------------------------------------------------------------------------
  art::Handle<mu2e::StrawHitCollection> shcH;
  const mu2e::StrawHitCollection*       shColl(nullptr);

  if (! fShCollTag.empty()) {
    AnEvent->getByLabel(fShCollTag,shcH);

    if (shcH.isValid()) {
      shColl = shcH.product();
      data->fNStrawHits = shColl->size();
    }
  }
//-----------------------------------------------------------------------------
// number of combo hits
//-----------------------------------------------------------------------------
  art::Handle<mu2e::ComboHitCollection> chcH;
  const mu2e::ComboHitCollection*       chColl(nullptr);

  if (! fChCollTag.empty()) {
    AnEvent->getByLabel(fChCollTag,chcH);

    if (chcH.isValid()) {
      chColl = chcH.product();
      data->fNComboHits = chColl->size();
    }
  }
//-----------------------------------------------------------------------------
// event weight
//-----------------------------------------------------------------------------
  data->fEventWeight      = 1.f;

  auto handles = AnEvent->getMany<mu2e::EventWeight>();
  for(auto handle : handles) {
    if(!handle || !handle.isValid()) continue;
    if(verbose) std::cout << " Weight = " << handle->weight() << std::endl;
    data->fEventWeight *= handle->weight();
  }
  if(verbose) std::cout << " Total Weight = " << data->fEventWeight << std::endl;

  return 0;
}

}
