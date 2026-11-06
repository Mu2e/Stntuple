//

#include "RtypesCore.h"
#include "Stntuple/obj/TStrawDigi.hh"
#include "Offline/DataProducts/inc/StrawId.hh"

ClassImp(TStrawDigi)
//-----------------------------------------------------------------------------
// V1: can read directly into the object, no need to rearrange anything
//-----------------------------------------------------------------------------
void TStrawDigi::ReadV1(TBuffer &R__b) {

  struct StrawDigiV1_t {
    int    fNs;                   //     N(ADC samples)
    int    fSid;
    int    fMnid;                          // 
    int    fTdc0;
    int    fTdc1;
    int    fTot0;
    int    fTot1;
    int    fPmp;
    int    fFlag;
    int    fFs;
    float  fBl;
    float  fPh;
    std::vector<uint16_t> fAdc;
  } dat;

  int nwi = ((int*)   &fBl ) - &fNs;
  int nwf = ((float*) &fAdc) - &fBl;

  SetUniqueID(kMaxUInt);     // such index should not exist
  R__b.ReadFastArray(&fNs, nwi);
  R__b.ReadFastArray(&fBl, nwf);
  if (fNs != (int) fAdc.size()) {
    fAdc.resize(fNs);
  }
  if (fNs != 0) {
    R__b.ReadFastArray(fAdc.data(),fNs);
  }
}

//-----------------------------------------------------------------------------
void TStrawDigi::Streamer(TBuffer& R__b) {
  
  int nwi = ((int*) &fBl   ) - &fNs;
  int nwf = ((float*) &fAdc) - &fBl;

  if (R__b.IsReading()) {
    Version_t R__v = R__b.ReadVersion();  // not used but want to look at 
    
    if      (R__v == 1) {
      ReadV1(R__b);
    }
    else {
                                        // current version: V2
      TObject::Streamer(R__b);
      R__b.ReadFastArray(&fNs, nwi);
      R__b.ReadFastArray(&fBl, nwf);
      if (fNs != (int) fAdc.size()) {
        fAdc.resize(fNs);
      }
      if (fNs != 0) {
        R__b.ReadFastArray(fAdc.data(),fNs);
      }
    }
  }
  else {
    R__b.WriteVersion(TStrawDigi::IsA());
    TObject::Streamer(R__b);
    R__b.WriteFastArray(&fNs,nwi);
    R__b.WriteFastArray(&fBl,nwf);
    if (fNs != 0) {
      R__b.WriteFastArray(fAdc.data(),fNs);
    }
  }
}

//-----------------------------------------------------------------------------
// 
//-----------------------------------------------------------------------------
TStrawDigi::TStrawDigi() : TObject(), fNs(0), fSid(-1), fMnid(-1), fTdc0(-1), fTdc1(-1),
                           fTot0(-1), fTot1(-1), fPmp(-1), fFlag(-1), fBl(-1.), fPh(-1.)
{
  SetUniqueID(kMaxUInt);
}

//-----------------------------------------------------------------------------
TStrawDigi::TStrawDigi(UInt_t ID) : TObject(), fNs(0), fSid(-1), fMnid(-1), fTdc0(-1), fTdc1(-1),
                                    fTot0(-1), fTot1(-1), fPmp(-1), fFlag(-1), fBl(-1.), fPh(-1.)
{
  SetUniqueID(ID);
}

//-----------------------------------------------------------------------------
TStrawDigi::~TStrawDigi() {
}

//-----------------------------------------------------------------------------
int TStrawDigi::Init(int Ns) {
  int rc(0);
  if (fNs != Ns) {
    fNs = Ns;
    fAdc.resize(fNs);
  }
  return rc;
}

//-----------------------------------------------------------------------------
void TStrawDigi::Clear(const char* Opt) {
  SetUniqueID(kMaxUInt);
}

//-----------------------------------------------------------------------------
void TStrawDigi::Print(const char* Option) const {
  TString opt = Option;
  opt.ToLower();

  if ((opt == "") || (opt.Index("banner") >= 0)) {
    printf("-----------------------------------------------------------------------------");
    printf("-----------------------------------------------------------------------------\n");
    printf("  SID pln pnl str  MNID    TDC0     TDC1  TOT0 TOT1  PMP    Flag     Fs    Bl      Ph");
    printf("   Ns                                ADC\n");
    printf("------v----------------------------------------------------------------------");
    printf("----------------------------------------------------------------------------\n");
  }

  if (opt.Index("data") < 0) return;

  mu2e::StrawId sid(fSid);
  
  printf("%5d  %2d  %2d  %2d MN%03d %8d %8d %4d %4d %4d 0x%08x  %2d %7.2f %7.2f %2d",
	 GetUniqueID(),
	 sid.plane(),
	 sid.panel(),
	 sid.straw(),
         fMnid,
	 fTdc0, fTdc1, // ## format
	 fTot0, fTot1,
	 fPmp,
	 fFlag,
         fFs,
         fBl,
         fPh,
         fNs);

  if (fNs <= 0) {
    printf("\n");
  }
  else {
    int pos=0;
    for (int i=0; i<fNs; i++) {
      printf(" %4d",fAdc[i]);
      pos++;
      if (pos >= 20) {
        // TODO: handle the second line - need to start not from the first position
        // ok for now
        printf("\n");
        pos = 0;
      }
    }
    if (pos > 0) printf("\n");
  }
  
}
