#include <iostream>
#include <iomanip>

#include "obj/TStnTimeClusterBlock.hh"
#include "obj/TStnTimeCluster.hh"

ClassImp(TStnTimeClusterBlock)

//-----------------------------------------------------------------------------
void TStnTimeClusterBlock::ReadV1(TBuffer &R__b) {
  R__b >> fNTimeClusters;
  fListOfTimeClusters->Streamer(R__b);
  // and keep fLiskOfChLinks untouched
}

//_____________________________________________________________________________
void TStnTimeClusterBlock::Streamer(TBuffer &R__b) {

  if (R__b.IsReading()) {
    Version_t R__v = R__b.ReadVersion(); if (R__v) { }
    if (R__v == 1) {
      ReadV1(R__b);
    }
    else {
      R__b >> fNTimeClusters;
      if (fNTimeClusters > 0) {
        fListOfTimeClusters->Streamer(R__b);
        fListOfChLinks->Streamer(R__b);
      }
    }
  } 
  else {
    R__b.WriteVersion(TStnTimeClusterBlock::IsA());
    R__b << fNTimeClusters;
    if (fNTimeClusters > 0) {
      fListOfTimeClusters->Streamer(R__b);
      fListOfChLinks->Streamer(R__b);
    }
  }
}


//-----------------------------------------------------------------------------
TStnTimeClusterBlock::TStnTimeClusterBlock() {
  fNTimeClusters      = 0;
  fListOfTimeClusters = new TClonesArray("TStnTimeCluster",100);
  //  fListOfTimeClusters->BypassStreamer(kFALSE);
  fListOfTimeClusters->BypassStreamer(kTRUE);
  fCollName  = "default";
  fListOfChLinks = new TStnLinkBlock();
}


//_____________________________________________________________________________
TStnTimeClusterBlock::~TStnTimeClusterBlock() {
  fListOfTimeClusters->Delete();
  delete fListOfTimeClusters;
  delete fListOfChLinks;
}


//_____________________________________________________________________________
void TStnTimeClusterBlock::Clear(Option_t* opt) {
  fNTimeClusters = 0;
  fListOfTimeClusters->Clear(opt);
  fListOfChLinks->Clear(opt);

  f_EventNumber       = -1;
  f_RunNumber         = -1;
  f_SubrunNumber      = -1;
  fLinksInitialized   =  0;
}

//------------------------------------------------------------------------------
void TStnTimeClusterBlock::Print(Option_t* opt) const {

  int banner_printed = 0;
  for (int i=0; i<fNTimeClusters; i++) {
    TStnTimeCluster* t = ((TStnTimeClusterBlock*) this)->TimeCluster(i);
    if (! banner_printed) {
      t->Print("banner");
      banner_printed = 1;
    }
    t->Print("data");
    // afer that try to print combo hit indices
    int nch = fListOfChLinks->NLinks(i);
    int ip = 0;
    for (int l=0; l<nch; l++) {
      std::cout << std::format("{:6}",fListOfChLinks->Index(i,l));
      ip++;
      if (ip == 20) {
        std::cout << std::endl;
        ip = 0;
      }
    }
    if (ip > 0) {
      std::cout << std::endl;
    }
  }
}
