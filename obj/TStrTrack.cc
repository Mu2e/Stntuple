//
#include <iostream>
#include <format>

#include "TBuffer.h"
#include "Stntuple/obj/TStrTrack.hh"

ClassImp(TStrTrack)

//-----------------------------------------------------------------------------
void TStrTrack::Streamer(TBuffer& R__b) {
  
  int nwi = ((int*  ) &fT0       ) - &fNHits;
  int nwf = ((float*) &fOfflineKs) - &fT0;

  if (R__b.IsReading()) {
    Version_t R__v = R__b.ReadVersion();  // not used but want to look at 
    
    if      (R__v == 1) { // ReadV1(R__b);
    // else if (R__v == 2) ReadV2(R__b);
    // else {
                                        // current version: V1
      TObject::Streamer(R__b);
      R__b.ReadFastArray(&fNHits,nwi);
      R__b.ReadFastArray(&fT0   ,nwf);
    }
  }
  else {
    R__b.WriteVersion(TStrTrack::IsA());
    TObject::Streamer(R__b);
    R__b.WriteFastArray(&fNHits,nwi);
    R__b.WriteFastArray(&fT0   ,nwf);
  }
}

//-----------------------------------------------------------------------------
TStrTrack::TStrTrack() : TObject() {
  Clear();
}

//-----------------------------------------------------------------------------
TStrTrack::TStrTrack(int ID) : TObject() {
  SetUniqueID(ID);
  Clear();
}

//-----------------------------------------------------------------------------
TStrTrack::~TStrTrack() {
}


//-----------------------------------------------------------------------------
void TStrTrack::Clear(const char* Opt) {
  fNHits  = -1;
  fNDof   = -1;
  fTcIndex = -1;
  
  fT0     = -1;
  fChi2   = -1;
  fX0     = -1;
  fY0     = -1;
  fZ0     = -1;
  fNx     = -1;
  fNy     = -1;
  fNz     = -1;
}

//-----------------------------------------------------------------------------
void TStrTrack::Print(const char* Opt) const {
  std::cout << std::format("TStrTrack::Print not implemented yet\n");
}
