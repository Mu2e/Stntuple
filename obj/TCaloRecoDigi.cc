//
#include <iostream>
#include <format>

#include "TBuffer.h"
#include "Stntuple/obj/TCaloRecoDigi.hh"

ClassImp(TCaloRecoDigi)

//-----------------------------------------------------------------------------
// frst version of the I/O with the schema evolution
//-----------------------------------------------------------------------------
void TCaloRecoDigi::ReadV1(TBuffer &R__b) {

  struct TCaloRecoDigi_V1 {
    int    fSipmID;
    int    fNdf;
    int    fPileup;
    float  fEDep;
    float  fSigE;
    float  fTime;
    float  fSigT;
    float  fChi2;
    const mu2e::CaloRecoDigi*  fOfflineCrd; //!
  } data;

  int nwi = ((int*  ) &data.fEDep      ) - &data.fSipmID;
  int nwf = ((float*) &data.fOfflineCrd) - &data.fEDep;
  
  TObject::Streamer(R__b);
  R__b.ReadFastArray(&data.fSipmID,nwi);
  R__b.ReadFastArray(&data.fEDep  ,nwf);

  fSipmID      = data.fSipmID;
  fNdf         = data.fNdf;
  fPileup      = data.fPileup;
  fCdIndex     = -1;                    // == added in V2 ==
  fEDep        = data.fEDep;
  fSigE        = data.fSigE;
  fTime        = data.fTime;
  fSigT        = data.fSigT;
  fChi2        = data.fChi2;
}

//-----------------------------------------------------------------------------
void TCaloRecoDigi::Streamer(TBuffer& R__b) {
  
  int nwi = ((int*  ) &fEDep) - &fSipmID;
  int nwf = (         &fNext) - &fEDep;

  if (R__b.IsReading()) {
    Version_t R__v = R__b.ReadVersion();  // not used but want to look at 
    
    if      (R__v == 1) {
      ReadV1(R__b);
    }
    else if (R__v == 2) {
                                        // current version: V2
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
  fSipmID  =  0;
  fNdf     = -1;
  fPileup  = -1;
  fCdIndex =  0;
  fEDep    = -1;
  fSigE    = -1;
  fTime    = -1;
  fSigT    = -1;
  fChi2    = -1;
}

//-----------------------------------------------------------------------------
void TCaloRecoDigi::Print(const char* Opt) const {
  TString opt = Opt;
  opt.ToLower();

  if ((opt == "") || (opt.Index("banner") >= 0)) {
    printf("---------------------------------------------------------000000---------\n");
    printf("   ID SipmID CdIndex   Time        EDep    SigT    SigE ndf  chi2 Pileup\n");
    printf("------------------------------------------------------------------------\n");
  }

  if (opt.Index("data") < 0) return;

  std::cout << std::format("{:5d} {:5d} {:5d} {:10.2f} {:10.3f} {:7.3f} {:7.3f} {:3} {:5.1f} {:4}\n",
                           GetUniqueID(),SipmID(), CdIndex(), fTime,fEDep,fSigT,fSigE,fNdf,fChi2,fPileup );  
}
