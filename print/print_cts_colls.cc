///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
#include "Offline/RecoDataProducts/inc/CosmicTrackSeed.hh"

#include "Stntuple/print/TAnaDump.hh"
#include "Stntuple/print/Stntuple_print_functions.hh"
//-----------------------------------------------------------------------------
// print all CosmicTrackSeed collections in the event
//-----------------------------------------------------------------------------
void print_cts_colls() {
  printf("Available CosmicTrackSeedCollections: ");

  const art::Event* event = TAnaDump::Instance()->Event();

  art::Selector  selector(art::ProductInstanceNameSelector(""));

  auto vcoll = event->getMany<mu2e::CosmicTrackSeedCollection>(selector);

  int ncoll = vcoll.size();
  if (ncoll == 0) {
    printf("none\n");
    return;
  }
  printf("\n");
  
  for (auto handle = vcoll.begin(); handle != vcoll.end(); handle++) {
    if (handle->isValid()) {
      const art::Provenance* prov = handle->provenance();
      
      printf("moduleLabel: %-20s, productInstanceName: %-20s, processName:= %-30s\n" ,
	     prov->moduleLabel().data(),
	     prov->productInstanceName().data(),
	     prov->processName().data()
	     );
    }
  }
  return;
}
