///////////////////////////////////////////////////////////////////////////////
//  Dec 07 2001 P.Murat: start putting in some comments
//  ---------------------------------------------------
// TStrTrackBlock: ROOT-parseable description of TCalData to be stored in 
//                STNTUPLE
///////////////////////////////////////////////////////////////////////////////
#include <format>
#include "TVector2.h"

#include "Stntuple/obj/TStrTrackBlock.hh"

ClassImp(TStrTrackBlock)

// //_____________________________________________________________________________
// void TStrTrackBlock::ReadV1(TBuffer &R__b) {

//   struct TStrTrackBlockV1_t {
//     int            fNHits;		// number of hit crystals
//     int            fNDisks;             // 
//     int            fNCrystals  [4];	// 
//     float          fRMin       [4];	// as a temporary measure, store 
//     float          fRMax       [4];
//     float          fZ0         [4];     // 
//     float          fCrystalSize;
//     float          fMinFraction;        // min fr of the included crystal area

//     TClonesArray*  fListOfStrTracks;	// list of crystal hit data 
//   };

//   TStrTrackBlockV1_t data; 

//   int nwi = ((int*  ) data.fRMin       ) - &data.fNHits;
//   int nwf = ((float*) &data.fListOfStrTracks) - data.fRMin;

//   R__b.ReadFastArray(&fNHits,nwi);
//   R__b.ReadFastArray(fRMin  ,nwf);

//   if (fNHits > 0) {
//     fListOfStrTracks->Streamer(R__b);
//   }
// 				// initialize V2 variables 
//   fWrapperThickness = 0.065;    // 65 microns
//   fShellThickness   = 0;
// }



//______________________________________________________________________________
void TStrTrackBlock::Streamer(TBuffer &R__b) {
  // Stream an object of class TStrTrackBlock.

  // int nwi = ((int*  ) &fListOfStrTracks) - &fNDigis;
  // int nwf = 0; // ((float*) &fListOfStrTracks) - fRMin;

  if (R__b.IsReading()) {
    Version_t R__v = R__b.ReadVersion(); 
    if      (R__v == 1) {         // ReadV1(R__b);
                                        // else if (R__v == 2) ReadV2(R__b);
      R__b >> fNTracks;
      if (fNTracks > 0) {
        fListOfTracks->Streamer(R__b);
      }
    }
    else {
//-----------------------------------------------------------------------------
// read version > 1 ???
//-----------------------------------------------------------------------------
      std::cout << std::format(">>> ERROR: TStrTrackBlock::Streamer read version:{}\n",R__v);
    } 
  }
  else {
    R__b.WriteVersion(TStrTrackBlock::IsA());
    R__b << fNTracks;

    if (fNTracks > 0) {
      fListOfTracks->Streamer(R__b);
    }
  }
}

//_____________________________________________________________________________
TStrTrackBlock::TStrTrackBlock() {

  fListOfTracks = new TClonesArray("TStrTrack",100);
  fListOfTracks->BypassStreamer(kFALSE);
  Clear();
}

//_____________________________________________________________________________
TStrTrackBlock::~TStrTrackBlock() {
  fListOfTracks->Delete();
  delete fListOfTracks;
}

//_____________________________________________________________________________
void TStrTrackBlock::Clear(Option_t* opt) {
  fListOfTracks->Clear();
  fNTracks            = 0;

  f_EventNumber       = -1;
  f_RunNumber         = -1;
  f_SubrunNumber      = -1;
  fLinksInitialized   =  0;
}

//_____________________________________________________________________________
void TStrTrackBlock::Print(Option_t* opt) const {
  // print all the towers in the list
  if (fNTracks > 0) {
    fListOfTracks->At(0)->Print("banner");
    for (int i=0; i<fNTracks; i++) {
      fListOfTracks->At(i)->Print();
    }
  }
}
