//
// 2026-09-07 P.Murat
//
#include <iostream>
#include <format>

#include "Stntuple/obj/TComboHitBlock.hh"

ClassImp(TComboHitBlock)

//-----------------------------------------------------------------------------
void TComboHitBlock::ReadV1(TBuffer& R__b) {
  R__b >> fNHits;
  fListOfHits->Streamer(R__b);
  fListOfShIndices->Clear();
}

//-----------------------------------------------------------------------------
// R_v is so far unused
//-----------------------------------------------------------------------------
  void TComboHitBlock::Streamer(TBuffer &R__b) {
  if(R__b.IsReading()) {
    Version_t R__v = R__b.ReadVersion();
    if (R__v == 1) ReadV1(R__b);
    else {
      R__b >> fNHits;
      fListOfHits->Streamer(R__b);
      fListOfShIndices->Streamer(R__b);
    }
  }
  else {
//-----------------------------------------------------------------------------
// current version = 1
//-----------------------------------------------------------------------------
    R__b.WriteVersion(TComboHitBlock::IsA());
    R__b << fNHits;
    fListOfHits->Streamer(R__b);
    fListOfShIndices->Streamer(R__b);
  }
}

//______________________________________________________________________________
TComboHitBlock::TComboHitBlock() {
  fListOfHits = new TClonesArray("TComboHit",10000);
  fListOfHits->BypassStreamer(kFALSE);
  fListOfShIndices = new TStnLinkBlock();
  Clear();
}

//______________________________________________________________________________
TComboHitBlock::~TComboHitBlock() {
  fListOfHits->Delete();
  delete fListOfHits;
  delete fListOfShIndices;
}

//-----------------------------------------------------------------------------
int TComboHitBlock::FindHitIndex(TComboHit* Ch) {
  int found(-1);
  for (int i=0; i<fNHits; i++) {
    TComboHit* hit = Hit(i);
    if (hit == Ch) {
      found = i;
      break;
    }
  }
  return found;
}

//-----------------------------------------------------------------------------
TComboHit* TComboHitBlock::GetHitByID(int ID) {
  TComboHit* found(nullptr);
  for (int i=0; i<fNHits; i++) {
    TComboHit* hit = Hit(i);
    if (hit->GetUniqueID() == (uint32_t) ID) {
      found = hit;
      break;
    }
  }
  return found;
}

//______________________________________________________________________________
void TComboHitBlock::Clear(Option_t* opt) {
  fListOfHits->Clear();
  fNHits      = 0;

  f_EventNumber       = -1;
  f_RunNumber         = -1;
  f_SubrunNumber      = -1;
  fLinksInitialized   =  0;
}

//-----------------------------------------------------------------------------
//  print all hits in the straw tracker
//-----------------------------------------------------------------------------
void TComboHitBlock::Print(Option_t* opt) const {
  TComboHitBlock* blk = (TComboHitBlock*) this;
  int banner_printed = 0;
  std::cout << std::format(" *** N(combo hits): {}\n",fNHits);
  for (int i=0; i<fNHits; i++) {
    TComboHit* hit = blk->Hit(i);
    if (banner_printed == 0) {
      hit->Print("banner");
      banner_printed = 1;
    }
    hit->Print("data");
  }
}
