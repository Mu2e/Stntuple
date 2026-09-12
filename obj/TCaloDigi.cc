//

#include "TBuffer.h"
#include "Stntuple/obj/TCaloDigi.hh"

ClassImp(TCaloDigi)

//-----------------------------------------------------------------------------
void TCaloDigi::Streamer(TBuffer& R__b) {
  
  int nwi = ((int*) &fWf   ) - &fNs;

  if (R__b.IsReading()) {
    Version_t R__v = R__b.ReadVersion();  // not used but want to look at 
    
    if      (R__v == 1) { // ReadV1(R__b);
    // else if (R__v == 2) ReadV2(R__b);
    // else {
                                        // current version: V1
      R__b.ReadFastArray(&fNs, nwi);
      if (fNs != (int) fWf.size()) {
        fWf.resize(fNs);
      }
      if (fNs != 0) {
        R__b.ReadFastArray(fWf.data(),fNs);
      }
    }
  }
  else {
    R__b.WriteVersion(TCaloDigi::IsA());

    R__b.WriteFastArray(&fNs,nwi);
    if (fNs != 0) {
      R__b.WriteFastArray(fWf.data(),fNs);
    }
  }
}

//-----------------------------------------------------------------------------
TCaloDigi::TCaloDigi() : TObject() {
  fNs = 0;
}

//-----------------------------------------------------------------------------
TCaloDigi::~TCaloDigi() {
}

//-----------------------------------------------------------------------------
void TCaloDigi::Set(int SipmID, float T0, float PeakPos, const std::vector<int>* Wf) {
  
  fSipmID = SipmID;
  fT0     = T0;
  fPpos = PeakPos;
  fNs      = Wf->size();
  if (fNs != (int) fWf.size()) {
    fWf.resize(fNs);
  }
  for (int i=0; i<fNs; i++) {
    fWf[i] = (*Wf)[i];
  }
}

//-----------------------------------------------------------------------------
int TCaloDigi::Init(int Ns) {
  int rc(0);
  if (fNs != Ns) {
    fNs = Ns;
    fWf.resize(fNs);
  }
  return rc;
}

//-----------------------------------------------------------------------------
void TCaloDigi::Clear(const char* Opt) {
}
