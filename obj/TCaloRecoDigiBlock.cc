///////////////////////////////////////////////////////////////////////////////
//  Dec 07 2001 P.Murat: start putting in some comments
//  ---------------------------------------------------
// TCaloRecoDigiBlock: ROOT-parseable description of TCalData to be stored in 
//                STNTUPLE
///////////////////////////////////////////////////////////////////////////////
#include <format>
#include "TVector2.h"

#include "Stntuple/obj/TCaloRecoDigiBlock.hh"

ClassImp(TCaloRecoDigiBlock)

// //_____________________________________________________________________________
// void TCaloRecoDigiBlock::ReadV1(TBuffer &R__b) {

//   struct TCaloRecoDigiBlockV1_t {
//     int            fNHits;		// number of hit crystals
//     int            fNDisks;             // 
//     int            fNCrystals  [4];	// 
//     float          fRMin       [4];	// as a temporary measure, store 
//     float          fRMax       [4];
//     float          fZ0         [4];     // 
//     float          fCrystalSize;
//     float          fMinFraction;        // min fr of the included crystal area

//     TClonesArray*  fListOfCaloRecoDigis;	// list of crystal hit data 
//   };

//   TCaloRecoDigiBlockV1_t data; 

//   int nwi = ((int*  ) data.fRMin       ) - &data.fNHits;
//   int nwf = ((float*) &data.fListOfCaloRecoDigis) - data.fRMin;

//   R__b.ReadFastArray(&fNHits,nwi);
//   R__b.ReadFastArray(fRMin  ,nwf);

//   if (fNHits > 0) {
//     fListOfCaloRecoDigis->Streamer(R__b);
//   }
// 				// initialize V2 variables 
//   fWrapperThickness = 0.065;    // 65 microns
//   fShellThickness   = 0;
// }



//______________________________________________________________________________
void TCaloRecoDigiBlock::Streamer(TBuffer &R__b) {
  // Stream an object of class TCaloRecoDigiBlock.

  // int nwi = ((int*  ) &fListOfCaloRecoDigis) - &fNDigis;
  // int nwf = 0; // ((float*) &fListOfCaloRecoDigis) - fRMin;

  if (R__b.IsReading()) {
    Version_t R__v = R__b.ReadVersion(); 
    if      (R__v == 1) {         // ReadV1(R__b);
                                        // else if (R__v == 2) ReadV2(R__b);
      R__b >> fNDigis;
      if (fNDigis > 0) {
        fListOfCaloRecoDigis->Streamer(R__b);
      }
    }
    else {
//-----------------------------------------------------------------------------
// read version > 1 ???
//-----------------------------------------------------------------------------
      std::cout << std::format(">>> ERROR: TCaloRecoDigiBlock::Streamer read version:{}\n",R__v);
    } 
  }
  else {
    R__b.WriteVersion(TCaloRecoDigiBlock::IsA());
    R__b << fNDigis;

    if (fNDigis > 0) {
      fListOfCaloRecoDigis->Streamer(R__b);
    }
  }
}

//_____________________________________________________________________________
TCaloRecoDigiBlock::TCaloRecoDigiBlock() {

  fListOfCaloRecoDigis = new TClonesArray("TCaloRecoDigi",100);
  fListOfCaloRecoDigis->BypassStreamer(kFALSE);
  Clear();
}

//_____________________________________________________________________________
TCaloRecoDigiBlock::~TCaloRecoDigiBlock() {
  fListOfCaloRecoDigis->Delete();
  delete fListOfCaloRecoDigis;
}

//_____________________________________________________________________________
void TCaloRecoDigiBlock::Clear(Option_t* opt) {
  fListOfCaloRecoDigis->Clear();
  fNDigis             = 0;

  f_EventNumber       = -1;
  f_RunNumber         = -1;
  f_SubrunNumber      = -1;
  fLinksInitialized   =  0;
}

//_____________________________________________________________________________
void TCaloRecoDigiBlock::Print(Option_t* opt) const {
  // print all the towers in the list
  if (fNDigis > 0) {
    fListOfCaloRecoDigis->At(0)->Print("banner");
    for (int i=0; i<fNDigis; i++) {
      fListOfCaloRecoDigis->At(i)->Print();
    }
  }
}
