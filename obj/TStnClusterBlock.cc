#include <iostream>
#include <iomanip>

#include "obj/TStnClusterBlock.hh"
#include "obj/TStnCluster.hh"

ClassImp(TStnClusterBlock)
//_____________________________________________________________________________
void TStnClusterBlock::ReadV1(TBuffer &R__b) {

  R__b >> fNClusters;
  fListOfClusters->Streamer(R__b);

  for (int i=0; i<fNClusters; i++) {
    Cluster(i)->SetNumber(i);
  }
  
  fListOfHitLinks->Clear();
}

//______________________________________________________________________________
void TStnClusterBlock::Streamer(TBuffer &R__b) {
  // Stream an object of class TStnClusterBlock.

  if (R__b.IsReading()) {
    Version_t R__v = R__b.ReadVersion();
    if (R__v == 1) {
      ReadV1(R__b);
    }
    else {
      // current version : V2
      R__b >> fNClusters;
      if (fNClusters > 0) {
        fListOfClusters->Streamer(R__b);
        for (int i=0; i<fNClusters; i++) {
          Cluster(i)->SetNumber(i);
        }
        fListOfHitLinks->Streamer(R__b);    // added in V2
      }
    }
  }
  else {
    R__b.WriteVersion(TStnClusterBlock::IsA());
    R__b << fNClusters;
    if (fNClusters > 0) {
      fListOfClusters->Streamer(R__b);
      fListOfHitLinks->Streamer(R__b);  // added in V2
    }
  }
}

//_____________________________________________________________________________
TStnClusterBlock::TStnClusterBlock() {
  fNClusters   = 0;
  fListOfClusters = new TClonesArray("TStnCluster",100);
  fListOfClusters->BypassStreamer(kFALSE);
  fListOfHitLinks = new TStnLinkBlock();
  fCollName  = "default";
}


//_____________________________________________________________________________
TStnClusterBlock::~TStnClusterBlock() {
  fListOfClusters->Delete();
  delete fListOfClusters;
  delete fListOfHitLinks;
}


//_____________________________________________________________________________
void TStnClusterBlock::Clear(Option_t* opt) {
  fNClusters = 0;
  fListOfClusters->Clear(opt);
  fListOfHitLinks->Clear(opt);

  f_EventNumber       = -1;
  f_RunNumber         = -1;
  f_SubrunNumber      = -1;
  fLinksInitialized   =  0;
}

//------------------------------------------------------------------------------
void TStnClusterBlock::Print(Option_t* opt) const {

  int banner_printed = 0;
  for (int i=0; i<fNClusters; i++) {
    TStnCluster* t = ((TStnClusterBlock*) this)->Cluster(i);
    if (! banner_printed) {
      t->Print("banner");
      banner_printed = 1;
    }
    t->Print("data");
  }
}
