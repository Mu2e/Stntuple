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
void TCaloHit::ReadV1(TBuffer &R__b) {

  struct TCaloHit_t_V1 {
    int            fCid;                  // crystal ID
    int            fNSipms;               // (number of R/O channels used, 1 or 2) || (n)digis) << 8
    float          fTime;                 // 
    float          fEDep;                 //
    float          fSigT;                 // uncertainty on T.
    float          fSigE;                 // uncertainty on E
  
    mu2e::CaloHit* fOfflineCaloHit;       //! transient
  } data;

  int nwi = ((int*  ) &data.fTime          ) - &data.fCid;
  int nwf = ((float*) &data.fOfflineCaloHit) - &data.fTime  ;

  R__b.ReadFastArray(&data.fCid ,nwi);
  R__b.ReadFastArray(&data.fTime,nwf);

  fCid         = data.fCid;
  fNSipms      = data.fNSipms;
  fCrdIndex[0] = -1;                    // == added in V2 ==
  fCrdIndex[1] = -1;                    // == added in V2 ==
  fTime        = data.fTime;
  fEDep        = data.fEDep;
  fSigT        = data.fSigT;
  fSigE        = data.fSigE;
}

//_____________________________________________________________________________
void TCaloHit::Streamer(TBuffer &R__b) {
  int nwi = ((int*  ) &fTime          ) - &fCid;
  int nwf = ((float*) &fOfflineCaloHit) - &fTime;
  
  if(R__b.IsReading()) {
    Version_t R__v = R__b.ReadVersion();
    if (R__v == 1) {
      ReadV1(R__b);
    }
    else {
//-----------------------------------------------------------------------------
// current version 2
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

//-----------------------------------------------------------------------------
// Options: "banner", "data"
//-----------------------------------------------------------------------------
void TCaloHit::Print(Option_t* Option) const {
  TString opt = Option;
  opt.ToLower();

  if ((opt == "") || (opt.Index("banner") >= 0)) {
    printf("-----------------------------------------------------\n");
    printf("   ID  CID NSipms   Time        EDep     SigT    SigE\n");
    printf("-----------------------------------------------------\n");
  }

  if (opt.Index("data") < 0) return;

  std::cout << std::format("{:5d} {:5d} {:4d} {:10.2f} {:10.3f} {:7.3f} {:7.3f}\n",
                           GetUniqueID(),fCid,fNSipms,fTime,fEDep,fSigT,fSigE);  
}
