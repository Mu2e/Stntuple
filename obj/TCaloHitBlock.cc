///////////////////////////////////////////////////////////////////////////////
//  Dec 07 2001 P.Murat: start putting in some comments
//  ---------------------------------------------------
// TCaloHitBlock: ROOT-parseable description of TCalData to be stored in 
//                STNTUPLE
///////////////////////////////////////////////////////////////////////////////
#include "TVector2.h"

#include "Stntuple/obj/TCaloHit.hh"
#include "Stntuple/obj/TCaloHitBlock.hh"

ClassImp(TCaloHitBlock)

//______________________________________________________________________________
void TCaloHitBlock::Streamer(TBuffer &R__b) {

  if (R__b.IsReading()) {
    // Version_t R__v = R__b.ReadVersion(); 
    R__b.ReadVersion(); 
    // if      (R__v == 1) ReadV1(R__b);
    // else if (R__v == 2) ReadV2(R__b);
    R__b >> fNHits;
    R__b >> fEDep[0];
    R__b >> fEDep[1];
    if (fNHits > 0) {
      fListOfCaloHits->Streamer(R__b);
    }
  }
  else {
    R__b.WriteVersion(TCaloHitBlock::IsA());

    R__b << fNHits;
    R__b << fEDep[0];
    R__b << fEDep[1];
    
    if (fNHits > 0) {
      fListOfCaloHits->Streamer(R__b);
    }
  }
}

//_____________________________________________________________________________
TCaloHitBlock::TCaloHitBlock() {
				// commented out pieces are available starting 
				// from version 2.24/05

  fListOfCaloHits    = new TClonesArray("TCaloHit",100);
  fListOfCaloHits->BypassStreamer(kFALSE);
  //  fAdcThreshold = 0;
  Clear();
}

//_____________________________________________________________________________
TCaloHitBlock::~TCaloHitBlock() {
  fListOfCaloHits->Delete();
  delete fListOfCaloHits;
}


//_____________________________________________________________________________
void TCaloHitBlock::Clear(Option_t* opt) {
  fListOfCaloHits->Clear();
  fNHits      = 0;
  fEDep[0]    = 0;
  fEDep[1]    = 0;

  f_EventNumber       = -1;
  f_RunNumber         = -1;
  f_SubrunNumber      = -1;
  fLinksInitialized   =  0;
}

//-----------------------------------------------------------------------------
void TCaloHitBlock::Print(Option_t* opt) const {
  // print all the towers in the list
  if (fNHits > 0) {
    fListOfCaloHits->At(0)->Print("banner");
    for (int i=0; i<fNHits; i++) {
      fListOfCaloHits->At(i)->Print("data");
    }
  }
}
