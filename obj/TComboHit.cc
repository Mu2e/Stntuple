///////////////////////////////////////////////////////////////////////////////
// 2020-09-07-26 P.Murat TComboHit
// - what code is setting the MCFlag for MC hits ?
///////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <format>
#include "TString.h"

#include "Stntuple/obj/TComboHit.hh"

ClassImp(TComboHit)

// //-----------------------------------------------------------------------------
// // V1: doesn't store the hit ID
// //-----------------------------------------------------------------------------
// void TComboHit::ReadV1(TBuffer &R__b) {
//   int nw_data = &fSimID-&fStrawID;

//   R__b.ReadFastArray(&fStrawID,nw_data);
//   R__b >> fEDep;
//   if (MCFlag()) {
//     int nwi_mc  = ((int*) &fEDep) - &fSimID;
//     R__b.ReadFastArray(&fSimID ,nwi_mc);
//     R__b >> fMcMom;
//   }
//   else {
//     fGenID       = -1;
//     fSimID       = -1;
//     fPdgID       = -1;
//     fMotherPdgID = -1;
//     fMcMom       = -1.;
//   }
// }

//_____________________________________________________________________________
void TComboHit::Streamer(TBuffer &R__b) {

  int nwi = ((int*  ) &fTime     )-&fStrawID;
  int nwf = ((float*) &fOfflineCh)-&fTime;
					// those are only in MC
  if (R__b.IsReading()) {
    // Version_t R__v = R__b.ReadVersion();
    R__b.ReadVersion();
    // if      (R__v == 1) ReadV1(R__b);
    // else if (R__v == 2) ReadV2(R__b);
    // else {
//-----------------------------------------------------------------------------
// curent version: V1
//-----------------------------------------------------------------------------
    TObject::Streamer(R__b);
    R__b.ReadFastArray(&fStrawID,nwi);
    R__b.ReadFastArray(&fTime   ,nwf);
  }
  else {
//-----------------------------------------------------------------------------
// write V3
//-----------------------------------------------------------------------------
    R__b.WriteVersion(TComboHit::IsA());
    TObject::Streamer(R__b);
    R__b.WriteFastArray(&fStrawID,nwi);
    R__b.WriteFastArray(&fTime   ,nwf);
  }
}

//-----------------------------------------------------------------------------
// comboHit ID is the hit index in the original Mu2e ComboHit collection
//-----------------------------------------------------------------------------
TComboHit::TComboHit(int ID): TObject() {
  SetUniqueID(ID);
  Clear();
}

//_____________________________________________________________________________
TComboHit::~TComboHit() {
}

//_____________________________________________________________________________
void TComboHit::Clear(Option_t* opt) {
}

//-----------------------------------------------------------------------------
// Options: "banner", "data"
//-----------------------------------------------------------------------------
void TComboHit::Print(Option_t* Option) const {
  // print Straw hit properties
  //  printf("Superlayer: %d, Wire: %d, Cell: %d,\n",fSuperLayer,fWire,fCell);
  
  TString opt = Option;
  opt.ToLower();

  if ((opt == "") || (opt.Index("banner") >= 0)) {
    printf("------------------------------------------------------------------------------------\n");
    printf("   ID   SID  MnID Nsh zface    Time     DrTime  EDep(keV)     X        Y         Z  \n");
    printf("------------------------------------------------------------------------------------\n");
  }

  if (opt.Index("data") < 0) return;

  std::cout << std::format("{:5} {:5} MN{:03d} {:2d} {:5d} {:10.2f} {:8.2f} {:8.2f} {:9.3f} {:9.3f} {:9.3f}\n",
                           GetUniqueID(),fStrawID,fMnid,fNsh,fZface, fTime,fDrTime,fEDep*1e3,fX,fY,fZ);  
}
