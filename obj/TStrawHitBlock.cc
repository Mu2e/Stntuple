//
// 2014-01-27 P.Murat
//

#include "Stntuple/obj/TStrawHitBlock.hh"

ClassImp(TStrawHitBlock)


//-----------------------------------------------------------------------------
void TStrawHitBlock::ReadV1(TBuffer& R__b) {
  R__b >> fNHits;
  fListOfHits->Streamer(R__b);
  fNWaveforms      = 0;
  // fListOfWaveforms = nullptr;
}

//-----------------------------------------------------------------------------
// R_v is so far unused
//-----------------------------------------------------------------------------
  void TStrawHitBlock::Streamer(TBuffer &R__b) {
  if(R__b.IsReading()) {
    Version_t R__v = R__b.ReadVersion();
    if (R__v == 1) ReadV1(R__b);
    else {
      R__b >> fNHits;
      R__b >> fNWaveforms;
      fListOfHits->Streamer(R__b);
      fListOfWaveforms->Streamer(R__b);
    }
  }
  else {
//-----------------------------------------------------------------------------
// current version = 2
//-----------------------------------------------------------------------------
    R__b.WriteVersion(TStrawHitBlock::IsA());
    R__b << fNHits;
    R__b << fNWaveforms;
    fListOfHits->Streamer(R__b);
    fListOfWaveforms->Streamer(R__b);
  }
}

//______________________________________________________________________________
TStrawHitBlock::TStrawHitBlock() {
  fListOfHits = new TClonesArray("TStrawHit",30240);
  fListOfHits->BypassStreamer(kFALSE);
  fListOfWaveforms = new TClonesArray("TStrWaveform",30240);
  fListOfWaveforms->BypassStreamer(kFALSE);
  Clear();
}

//______________________________________________________________________________
TStrawHitBlock::~TStrawHitBlock() {
  fListOfHits->Delete();
  delete fListOfHits;

  if (fListOfWaveforms) {
    fListOfWaveforms->Delete();
    delete fListOfWaveforms;
  }
}

//______________________________________________________________________________
void TStrawHitBlock::Clear(Option_t* opt) {
  fListOfHits->Clear();
  if (fListOfWaveforms) fListOfWaveforms->Clear();
  fNHits      = 0;
  fNWaveforms = 0;

  f_EventNumber       = -1;
  f_RunNumber         = -1;
  f_SubrunNumber      = -1;
  fLinksInitialized   =  0;
}

//-----------------------------------------------------------------------------
//  print all hits in the straw tracker
//-----------------------------------------------------------------------------
void TStrawHitBlock::Print(Option_t* opt) const {
  TStrawHitBlock* blk = (TStrawHitBlock*) this;
  int banner_printed = 0;
  printf(" *** N(straw hits): %d, N(waveforms): %d\n",fNHits,fNWaveforms);
  for (int i=0; i<fNHits; i++) {
    TStrawHit* hit = blk->Hit(i);
    if (banner_printed == 0) {
      hit->Print("banner");
      banner_printed = 1;
    }
    hit->Print("data");
  }

  for (int i=0; i<fNWaveforms; i++) {
    TStrWaveform* wf = blk->Waveform(i);
    if (banner_printed == 0) {
      wf->Print("banner");
      banner_printed = 1;
    }
    wf->Print("data");
  }
}
