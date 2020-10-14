///////////////////////////////////////////////////////////////////////////////
//  2014-01-26 P.Murat TCrvRecoPulse
///////////////////////////////////////////////////////////////////////////////
#include "TString.h"
#include "TBuffer.h"

#include "Stntuple/obj/TCrvRecoPulse.hh"

ClassImp(TCrvRecoPulse)

//_____________________________________________________________________________
void TCrvRecoPulse::Streamer(TBuffer &R__b) {

  int nwi = ((int*) &fTime) - &fIndex;
  int nwf = &fLeTime - &fTime +1;
  
  if (R__b.IsReading()) {
    //    Version_t R__v = R__b.ReadVersion();
    R__b.ReadVersion();
//-----------------------------------------------------------------------------
// curent version: V1
//-----------------------------------------------------------------------------
    R__b.ReadFastArray(&fIndex,nwi);
    R__b.ReadFastArray(&fTime ,nwf);
  }
  else {
    R__b.WriteVersion(TCrvRecoPulse::IsA());
    R__b.WriteFastArray(&fIndex,nwi);
    R__b.WriteFastArray(&fTime ,nwf);
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
void TCrvRecoPulse::Set(int I, int NPe, int NPeHeight, int NDigis, int Bar, int Sipm, 
			float Time, float Height, float Width, float Chi2, float LeTime)
{
  fIndex    = I;
  fNPe       = NPe;
  fNPeHeight = NPeHeight;
  fNDigis    = NDigis;
  fBar       = Bar;
  fSipm      = Sipm;
  fTime      = Time;
  fHeight    = Height;
  fWidth     = Width;
  fChi2      = Chi2;
  fLeTime    = LeTime;
}

//_____________________________________________________________________________
void TCrvRecoPulse::Clear(Option_t* opt) {
  fIndex          = -1;
  fNPe            = -1;
  fNPeHeight      = -1;
  fNDigis         = -1;
  fBar            = -1;
  fSipm           = -1;
  fTime           = -1;
  fHeight         = -1;
  fWidth          = -1;
  fChi2           = -1;
  fLeTime         = -1;
}

//_____________________________________________________________________________
void TCrvRecoPulse::Print(Option_t* Option) const {
  // print Straw hit properties
  //  printf("Superlayer: %d, Wire: %d, Cell: %d,\n",fSuperLayer,fWire,fCell);
  
  TString opt = Option;
  opt.ToLower();

  if ((opt == "") || (opt.Index("banner") >= 0)) {
    printf("---------------------------------------------------------------------------\n");
    printf(" NPE   HPE  NDigis Bar   Sipm  Time     Height     Width      Chi2   LeTime   \n");
    printf("---------------------------------------------------------------------------\n");
  }
 
  if ((opt == "") || (opt.Index("data") >= 0)) {

    printf("%5i %5i %5i %5i %5i %8.3f %8.3f %8.3f %10.3f %8.3f",
	   fNPe,
	   fNPeHeight,
	   fNDigis,
	   fBar,
	   fSipm,
	   fTime,
	   fHeight,
	   fWidth,
	   fChi2,
	   fLeTime);
    printf("\n");
  }
}
