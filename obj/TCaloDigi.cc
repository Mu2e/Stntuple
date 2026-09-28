//
#include <iostream>
#include <format>

#include "TBuffer.h"
#include "Stntuple/obj/TCaloDigi.hh"

ClassImp(TCaloDigi)

//-----------------------------------------------------------------------------
// first version of the I/O with the schema evolution
// doesn't read the TObject back
//-----------------------------------------------------------------------------
void TCaloDigi::ReadV1(TBuffer &R__b) {

  struct TCaloDigi_V1 {
    int                   fNs;
    int                   fSipmID;
    int                   fT0;
    int                   fPPos;          // peak position
    std::vector<uint16_t> fWf;
  } data;

  int nwi = ((int*) &fWf   ) - &fNs;
  R__b.ReadFastArray(&data.fNs,nwi);

  fNs     = data.fNs;
  fSipmID = data.fSipmID;
  fT0     = data.fT0;
  fPPos   = data.fPPos;
  
  if (fNs != (int) fWf.size()) {
    fWf.resize(fNs);
  }
  if (fNs != 0) {
    R__b.ReadFastArray(fWf.data(),fNs);
  }

}

//-----------------------------------------------------------------------------
void TCaloDigi::Streamer(TBuffer& R__b) {
  
  int nwi = ((int*) &fWf   ) - &fNs;

  if (R__b.IsReading()) {
    Version_t R__v = R__b.ReadVersion();  // not used but want to look at 
    
    if      (R__v == 1) {
      ReadV1(R__b);
    }
    else if (R__v == 2) {
                                        // current version: V2 - adds TObject part
      TObject::Streamer(R__b);
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
    TObject::Streamer(R__b);
    R__b.WriteFastArray(&fNs,nwi);
    if (fNs != 0) {
      R__b.WriteFastArray(fWf.data(),fNs);
    }
  }
}

//-----------------------------------------------------------------------------
TCaloDigi::TCaloDigi() : TObject() {
  fNs     =  0;
  fSipmID = -1;
  fT0     = -1;
  fPPos   = -1;
}

//-----------------------------------------------------------------------------
TCaloDigi::TCaloDigi(int ID) : TObject() {
  SetUniqueID(ID);
  fNs     = 0;
  fSipmID = -1;
  fT0     = -1;
  fPPos   = -1;
}

//-----------------------------------------------------------------------------
TCaloDigi::~TCaloDigi() {
}

//-----------------------------------------------------------------------------
void TCaloDigi::Set(int SipmID, float T0, float PeakPos, const std::vector<int>* Wf) {
  
  fSipmID = SipmID;
  fT0     = T0;
  fPPos = PeakPos;
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

//-----------------------------------------------------------------------------
void TCaloDigi::Print(const char* Opt) const {
  TString opt = Opt;
  opt.ToLower();
  
  if ((opt == "") || (opt.Index("banner") >= 0)) {
    printf("------------------------------------------------------------------\n");
    printf("   ID SipmID   Mask   T0     N s  PPos                  Wf              \n");
    printf("------------------------------------------------------------------\n");
  }
  
  if ((opt == "") || (opt.Index("data") >= 0)) {

    std::cout << std::format("{:5d} {:5d} 0x{:04x} {:6d} {:6d} {:6d} *** ",
                             GetUniqueID(),SipmID(),Mask(),fT0,fNs,fPPos);
    int ns = fWf.size();
    if (ns == fNs) {
      int pos = 0;
      for (int i=0; i<ns; i++) {
        printf("%5d",fWf[i]);
        pos++;
        if (pos == 20) {
          printf("\n");
          pos = 0;
          if (i == ns-1) break;
          printf("%43s","");
        }
      }
      if (pos != 0) {
        printf("\n");
      }
    }
    else {
      std::cout << std::format("ERROR: ns:{} not equal to fNs:{}. BAIL OUT\n",
                               ns,fNs);
    }
  }
}
