///////////////////////////////////////////////////////////////////////////////
//  2014-01-26 P.Murat TCaloHit
///////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <format>

#include "Stntuple/obj/TCaloHit.hh"

ClassImp(TCaloHit)

//-----------------------------------------------------------------------------
// frst version of the I/O with the schema evolution
//-----------------------------------------------------------------------------
// void TCaloHit::ReadV1(TBuffer &R__b) {

//   struct TCaloHit_t_V1 {
//     int        fID;         // hit ID, cods disk,  x1, x2
//     int        fNChannels;  // number of readout channels, 1 or 2 (kludge)
//     float      fTime; 
//     float      fEnergy;
//     void*      fDummy;      // !
//   } data;
  
//   int nwi = ((int*  ) &data.fTime ) - &data.fID;
//   int nwf = ((float*) &data.fDummy) - &data.fTime  ;
  
//   R__b.ReadFastArray(&data.fID  ,nwi);
//   R__b.ReadFastArray(&data.fTime,nwf);

//   fCid    = data.fID;
//   fNSipms = data.fNChannels;
//   fTime   = data.fTime;
//   fEDep   = data.fEnergy;
//   fSigT   = -1;
//   fSigE   = -1;
// }

//_____________________________________________________________________________
void TCaloHit::Streamer(TBuffer &R__b) {
  int nwi = ((int*) &fTime)             - &fCid;
  int nwf = ((float*) &fOfflineCaloHit) - &fTime;
  
  if(R__b.IsReading()) {
    Version_t R__v = R__b.ReadVersion();
    if (R__v == 1) { // ReadV1(R__b);
//-----------------------------------------------------------------------------
// current version 1
//-----------------------------------------------------------------------------
      TObject::Streamer(R__b);
      R__b.ReadFastArray(&fCid ,nwi);
      R__b.ReadFastArray(&fTime,nwf);
    }
  }
  else {
    R__b.WriteVersion(TCaloHit::IsA());
    TObject::Streamer(R__b);
    R__b.WriteFastArray(&fCid ,nwi);
    R__b.WriteFastArray(&fTime,nwf);
  } 
}

//_____________________________________________________________________________
TCaloHit::TCaloHit(): TObject() {
  Clear();
}

//_____________________________________________________________________________
TCaloHit::TCaloHit(int ID): TObject() {
  SetUniqueID(ID);
  Clear();
}

//_____________________________________________________________________________
TCaloHit::~TCaloHit() {
}

//_____________________________________________________________________________
void TCaloHit::Clear(Option_t* opt) {
  fCid       = -1;
  fNSipms    = -1;
  fTime      = -1.;
  fEDep      = -1.;
  fSigT      = -1.;
  fSigE      = -1.;
}

//_____________________________________________________________________________
void TCaloHit::Print(Option_t* opt) const {
  std::cout << std::format("WARNING: TCaloHit::Print not implemented yet\n");
}
