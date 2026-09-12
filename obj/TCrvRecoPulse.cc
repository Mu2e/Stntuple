///////////////////////////////////////////////////////////////////////////////
//  2014-01-26 P.Murat TCrvRecoPulse
///////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <format>

#include "TString.h"
#include "TBuffer.h"

#include "Stntuple/obj/TCrvRecoPulse.hh"

ClassImp(TCrvRecoPulse)

//_____________________________________________________________________________
void TCrvRecoPulse::Streamer(TBuffer &R__b) {

  int nwi = ((int*) &fPes) - &fSbid;
  int nwf = ((float*) (&fOfflineCrvp)) - &fPes;
  
  if (R__b.IsReading()) {
    //    Version_t R__v = R__b.ReadVersion();
    R__b.ReadVersion();
//-----------------------------------------------------------------------------
// curent version: V1
//-----------------------------------------------------------------------------
    R__b.ReadFastArray(&fSbid,nwi);
    R__b.ReadFastArray(&fPes ,nwf);
  }
  else {
    R__b.WriteVersion(TCrvRecoPulse::IsA());
    R__b.WriteFastArray(&fSbid,nwi);
    R__b.WriteFastArray(&fPes ,nwf);
  } 
}

//_____________________________________________________________________________
TCrvRecoPulse::TCrvRecoPulse(): TObject() {
  Clear();
}

//_____________________________________________________________________________
TCrvRecoPulse::~TCrvRecoPulse() {
}

//_____________________________________________________________________________
void TCrvRecoPulse::Set(int Sbid, int Sipm, int Roc, int Feb, int FebCh,
                        float Pes, float PesPh, float Time, float Ph, float Beta,
                        float Chi2, float LeTime, float Ped) {
  fSbid      = Sbid;
  fSipm      = Sipm;
  fRoc       = Roc;
  fFeb       = Feb;
  fFebCh     = FebCh;

  fPes       = Pes;
  fPesPh     = PesPh;
  fTime      = Time;
  fPh        = Ph;
  fBeta      = Beta;
  fChi2      = Chi2;
  fLeTime    = LeTime;
  fPed       = Ped;
}

//_____________________________________________________________________________
void TCrvRecoPulse::Clear(Option_t* opt) {
  fSbid      = -1;
  fSipm      = -1;
  fRoc       = -1;
  fFeb       = -1;
  fFebCh     = -1;
  
  fPes       = 0;
  fPesPh     = 0;
  fTime      = 0;
  fPh        = 0;
  fBeta      = 0;
  fChi2      = -1;
  fLeTime    = 0;
  fPed       = -1;
}

//_____________________________________________________________________________
void TCrvRecoPulse::Print(Option_t* Option) const {
  // print Straw hit properties
  //  printf("Superlayer: %d, Wire: %d, Cell: %d,\n",fSuperLayer,fWire,fCell);
  
  TString opt = Option;
  opt.ToLower();

  if ((opt == "") || (opt.Index("banner") >= 0)) {
    printf("---------------------------------------------------------------------------\n");
    printf(" Sbid Sipm Roc Feb FebCh  Time    Pes   PesPh    Beta    Chi2    LeTime    Ped  \n");
    printf("---------------------------------------------------------------------------\n");
  }
 
  if ((opt == "") || (opt.Index("data") >= 0)) {

    std::cout << std::format("{:5} {:5} {:4} {:4} {:4} {:8.3f} {:8.3f} {:8.3f} {:8.2f} {:8.2f} {:8.2f}\n",
                             fSbid, fSipm, fRoc, fFeb, fFebCh,
                             fTime, fPes,  fPesPh, fBeta, fChi2, fLeTime, fPed);
  }
}
