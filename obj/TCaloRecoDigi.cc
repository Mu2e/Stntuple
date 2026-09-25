//
#include <iostream>
#include <format>

#include "TBuffer.h"
#include "Stntuple/obj/TCaloRecoDigi.hh"

ClassImp(TCaloRecoDigi)

//-----------------------------------------------------------------------------
void TCaloRecoDigi::Streamer(TBuffer& R__b) {
  
  int nwi = ((int*  ) &fEDep      ) - &fSipmID;
  int nwf = ((float*) &fOfflineCrd) - &fEDep;

  if (R__b.IsReading()) {
    Version_t R__v = R__b.ReadVersion();  // not used but want to look at 
    
    if      (R__v == 1) { // ReadV1(R__b);
    // else if (R__v == 2) ReadV2(R__b);
    // else {
                                        // current version: V1
      TObject::Streamer(R__b);
      R__b.ReadFastArray(&fSipmID, nwi);
      R__b.ReadFastArray(&fEDep  , nwf);
    }
  }
  else {
    R__b.WriteVersion(TCaloRecoDigi::IsA());
    TObject::Streamer(R__b);
    R__b.WriteFastArray(&fSipmID,nwi);
    R__b.WriteFastArray(&fEDep  ,nwf);
  }
}

//-----------------------------------------------------------------------------
TCaloRecoDigi::TCaloRecoDigi() : TObject() {
  Clear();
}

//-----------------------------------------------------------------------------
TCaloRecoDigi::TCaloRecoDigi(int ID) : TObject() {
  SetUniqueID(ID);
  Clear();
}

//-----------------------------------------------------------------------------
TCaloRecoDigi::~TCaloRecoDigi() {
}


//-----------------------------------------------------------------------------
void TCaloRecoDigi::Clear(const char* Opt) {
  fSipmID = -1;
  fNdf    = -1;
  fPileup = -1;
  fEDep   = -1;
  fSigE   = -1;
  fTime   = -1;
  fSigT   = -1;
  fChi2   = -1;
}

//-----------------------------------------------------------------------------
void TCaloRecoDigi::Print(const char* Opt) const {
  TString opt = Opt;
  opt.ToLower();

  if ((opt == "") || (opt.Index("banner") >= 0)) {
    printf("------------------------------------------------------------------\n");
    printf("   ID SipmID    Time        EDep     SigT    SigE ndf  chi2 Pileup\n");
    printf("------------------------------------------------------------------\n");
  }

  if (opt.Index("data") < 0) return;

  std::cout << std::format("{:5d} {:5d} {:10.2f} {:10.3f} {:7.3f} {:7.3f} {:3} {:5.1f} {:4}\n",
                           GetUniqueID(),fSipmID, fTime,fEDep,fSigT,fSigE,fNdf,fChi2,fPileup );  
}
