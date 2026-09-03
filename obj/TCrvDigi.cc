//

#include "TBuffer.h"
#include "Stntuple/obj/TCrvDigi.hh"

ClassImp(TCrvDigi)

//-----------------------------------------------------------------------------
void TCrvDigi::Streamer(TBuffer& R__b) {
  
  int nwi = ((int*) &fAdc   ) - &fNs;

  if (R__b.IsReading()) {
    Version_t R__v = R__b.ReadVersion();  // not used but want to look at 
    
    // if      (R__v == 1) ReadV1(R__b);
    // else if (R__v == 2) ReadV2(R__b);
    // else {
                                        // current version: V3
    R__b.ReadFastArray(&fNs, nwi);
    if (fNs != fAdc.size()) {
      fAdc.resize(fNs);
    }
    if (fNs != 0) {
      R__b.ReadFastArray(fAdc.data(),fNs);
    }
  }
  else {
    R__b.WriteVersion(TCrvDigi::IsA());

    R__b.WriteFastArray(&fNs,nwi);
    if (fNs != 0) {
      R__b.WriteFastArray(fAdc.data(),fNs);
    }
  }
}

//-----------------------------------------------------------------------------
TCrvDigi::TCrvDigi() : TObject() {
  fNs = 0;
}

//-----------------------------------------------------------------------------
TCrvDigi::~TCrvDigi() {
}

//-----------------------------------------------------------------------------
int TCrvDigi::Init(int Ns) {
  int rc(0);
  if (fNs != Ns) {
    fNs = Ns;
    fAdc.resize(fNs);
  }
  return 0;
}

//-----------------------------------------------------------------------------
void TCrvDigi::Clear(const char* Opt) {
}
