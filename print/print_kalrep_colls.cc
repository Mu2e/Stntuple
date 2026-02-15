///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
// #include "Offline/RecoDataProducts/inc/KalRepPtrCollection.hh"

#include "Stntuple/print/TAnaDump.hh"
#include "Stntuple/print/Stntuple_print_functions.hh"
//--------------------------------------------------------------------------------------------------------------
void print_kalrep_colls() {
  printf("Available KalRepPtrCollections: \n");

  // const art::Event* event = TAnaDump::Instance()->Event();

  // art::Selector  selector(art::ProductInstanceNameSelector(""));

  // // std::vector<art::Handle<mu2e::KalRepPtrCollection>> vcoll;
  // auto vcoll = event->getMany<mu2e::KalRepPtrCollection>(selector);

  // for (auto handle = vcoll.begin(); handle != vcoll.end(); handle++) {
  //   if (handle->isValid()) {
  //     const art::Provenance* prov = handle->provenance();
      
  //     printf("moduleLabel: %-20s, productInstanceName: %-20s, processName:= %-30s\n" ,
  //            prov->moduleLabel().data(),
  //            prov->productInstanceName().data(),
  //            prov->processName().data()
  //            );
  //   }
  // }
  // return;
}
