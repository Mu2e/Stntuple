///////////////////////////////////////////////////////////////////////////////
//  Dec 07 2001 P.Murat: start putting in some comments
//  ---------------------------------------------------
// TCaloDigiBlock: ROOT-parseable description of TCalData to be stored in 
//                STNTUPLE
///////////////////////////////////////////////////////////////////////////////
#include <format>
#include "TVector2.h"

#include "Stntuple/obj/TCaloDigiBlock.hh"

ClassImp(TCaloDigiBlock)

// //_____________________________________________________________________________
// void TCaloDigiBlock::ReadV1(TBuffer &R__b) {

//   struct TCaloDigiBlockV1_t {
//     int            fNHits;		// number of hit crystals
//     int            fNDisks;             // 
//     int            fNCrystals  [4];	// 
//     float          fRMin       [4];	// as a temporary measure, store 
//     float          fRMax       [4];
//     float          fZ0         [4];     // 
//     float          fCrystalSize;
//     float          fMinFraction;        // min fr of the included crystal area

//     TClonesArray*  fListOfCaloDigis;	// list of crystal hit data 
//   };

//   TCaloDigiBlockV1_t data; 

//   int nwi = ((int*  ) data.fRMin       ) - &data.fNHits;
//   int nwf = ((float*) &data.fListOfCaloDigis) - data.fRMin;

//   R__b.ReadFastArray(&fNHits,nwi);
//   R__b.ReadFastArray(fRMin  ,nwf);

//   if (fNHits > 0) {
//     fListOfCaloDigis->Streamer(R__b);
//   }
// 				// initialize V2 variables 
//   fWrapperThickness = 0.065;    // 65 microns
//   fShellThickness   = 0;
// }



//______________________________________________________________________________
void TCaloDigiBlock::Streamer(TBuffer &R__b) {
  // Stream an object of class TCaloDigiBlock.

  // int nwi = ((int*  ) &fListOfCaloDigis) - &fNDigis;
  // int nwf = 0; // ((float*) &fListOfCaloDigis) - fRMin;

  if (R__b.IsReading()) {
    Version_t R__v = R__b.ReadVersion(); 
    if      (R__v == 1) {         // ReadV1(R__b);
                                        // else if (R__v == 2) ReadV2(R__b);
      R__b >> fNDigis;
      if (fNDigis > 0) {
        fListOfCaloDigis->Streamer(R__b);
      }
    }
    else {
//-----------------------------------------------------------------------------
// read version > 1 ???
//-----------------------------------------------------------------------------
      std::cout << std::format(">>> ERROR: TCaloDigiBlock::Streamer read version:{}\n",R__v);
    } 
  }
  else {
    R__b.WriteVersion(TCaloDigiBlock::IsA());
    R__b << fNDigis;

    if (fNDigis > 0) {
      fListOfCaloDigis->Streamer(R__b);
    }
  }
}

//_____________________________________________________________________________
TCaloDigiBlock::TCaloDigiBlock() {

  fListOfCaloDigis    = new TClonesArray("TCaloDigi",100);
  fListOfCaloDigis->BypassStreamer(kFALSE);
  Clear();
}

//_____________________________________________________________________________
TCaloDigiBlock::~TCaloDigiBlock() {
  fListOfCaloDigis->Delete();
  delete fListOfCaloDigis;
}

//_____________________________________________________________________________
void TCaloDigiBlock::Clear(Option_t* opt) {
  fListOfCaloDigis->Clear();
  fNDigis             = 0;

  f_EventNumber       = -1;
  f_RunNumber         = -1;
  f_SubrunNumber      = -1;
  fLinksInitialized   =  0;
}

//_____________________________________________________________________________
void TCaloDigiBlock::Print(Option_t* opt) const {
  // print all the towers in the list
  if (fNDigis > 0) {
    fListOfCaloDigis->At(0)->Print("banner");
    for (int i=0; i<fNDigis; i++) {
      fListOfCaloDigis->At(i)->Print();
    }
  }
}
