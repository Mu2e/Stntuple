///////////////////////////////////////////////////////////////////////////////
//  2014-01-26 P.Murat TStrawHit
// - what code is setting the MCFlag for MC hits ?
///////////////////////////////////////////////////////////////////////////////
#include "TString.h"

#include "Stntuple/obj/TStrawHit.hh"

ClassImp(TStrawHit)

//-----------------------------------------------------------------------------
// V2: for an unknown reason, has times integer
//-----------------------------------------------------------------------------
void TStrawHit::ReadV2(TBuffer &R__b) {
  struct TStrawHit_v2 {
    int     fStrawID;			// sid | (mcflag << 16)
    int     fTime[2];			// 
    int     fTOT;                       // TOT[0] + TOT[1]<<16
					// -------------------- MC info, don't write out for the data
    int     fGenID;                     // generator ID + StrawDigiFlag << 16
    int     fSimID;			// sim particle ID
    int     fPdgID;			// sim particle PDG ID
    int     fMotherPdgID;	        // mother PDG ID
	  				// ----------------------------- floats
    float   fEDep;			// energy deposition
    float   fMcMom;			// MC particle momentum
  } sh;

  int nw_data = &sh.fSimID-&sh.fStrawID;
  int nwi_mc  = ((int*) &fEDep) - &fSimID;

  TObject::Streamer(R__b);
  R__b.ReadFastArray(&sh.fStrawID,nw_data);
  R__b >> sh.fEDep;
  if (MCFlag()) {
    R__b.ReadFastArray(&sh.fSimID ,nwi_mc);
    R__b >> sh.fMcMom;
  }
  else {
    sh.fGenID       = -1;
    sh.fSimID       = -1;
    sh.fPdgID       = -1;
    sh.fMotherPdgID = -1;
    sh.fMcMom       = -1.;
  }

  fStrawID     = sh.fStrawID;
  fGenID       = sh.fGenID;
  fSimID       = sh.fSimID;
  fPdgID       = sh.fPdgID;
  fMotherPdgID = sh.fMotherPdgID;

  fTime[0]     = sh.fTime[0];
  fTime[1]     = sh.fTime[1];

  fTOT[0]      = (sh.fTOT      ) & 0XFFFF;
  fTOT[1]      = (sh.fTOT >> 16) & 0XFFFF;
  
  fEDep        = sh.fEDep;
  fMcMom       = sh.fMcMom;
}

//-----------------------------------------------------------------------------
// V1: doesn't store the hit ID
//-----------------------------------------------------------------------------
void TStrawHit::ReadV1(TBuffer &R__b) {
  int nw_data = &fSimID-&fStrawID;

  R__b.ReadFastArray(&fStrawID,nw_data);
  R__b >> fEDep;
  if (MCFlag()) {
    int nwi_mc  = ((int*) &fEDep) - &fSimID;
    R__b.ReadFastArray(&fSimID ,nwi_mc);
    R__b >> fMcMom;
  }
  else {
    fGenID       = -1;
    fSimID       = -1;
    fPdgID       = -1;
    fMotherPdgID = -1;
    fMcMom       = -1.;
  }
}

//_____________________________________________________________________________
void TStrawHit::Streamer(TBuffer &R__b) {

  int nwi_mc = (int*)fTime-&fGenID;
  int nwf    = &fMcMom-fTime;
					// those are only in MC
  if (R__b.IsReading()) {
    Version_t R__v = R__b.ReadVersion();
    if      (R__v == 1) ReadV1(R__b);
    else if (R__v == 2) ReadV2(R__b);
    else {
//-----------------------------------------------------------------------------
// curent version: V3 - 
//-----------------------------------------------------------------------------
      TObject::Streamer(R__b);
      R__b >> fStrawID;

      fGenID       = -1;
      fSimID       = -1;
      fPdgID       = -1;
      fMotherPdgID = -1;
      fMcMom       = -1.;
      if (MCFlag()) R__b.ReadFastArray(&fGenID ,nwi_mc);
      R__b.ReadFastArray(fTime,nwf);
      if (MCFlag()) R__b >> fMcMom;
    }
  }
  else {
//-----------------------------------------------------------------------------
// write V3
//-----------------------------------------------------------------------------
    R__b.WriteVersion(TStrawHit::IsA());
    TObject::Streamer(R__b);
    R__b << fStrawID;
    if (MCFlag()) R__b.WriteFastArray(&fGenID ,nwi_mc);
    
    R__b.WriteFastArray(fTime,nwf);
    
    if (MCFlag()) R__b << fMcMom;
  }
}

//-----------------------------------------------------------------------------
TStrawHit::TStrawHit(int ID): TObject() {
  SetUniqueID(ID);
  Clear();
}

//_____________________________________________________________________________
TStrawHit::~TStrawHit() {
}

//_____________________________________________________________________________
// straw ID has a MC flag hidden in it
//-----------------------------------------------------------------------------
void TStrawHit::Set(int StrawID, float* Time, float* TOT, 
		    int GenID, int SimID, int PdgID, int MotherPdgID, 
		    float EDep, float McMom) 
{
  fStrawID     = StrawID; 
  fTime[0]     = Time[0]; 
  fTime[1]     = Time[1]; 
  fTOT[0]      = TOT[0];
  fTOT[1]      = TOT[1];
  fGenID       = GenID;
  fSimID       = SimID;
  fPdgID       = PdgID;
  fMotherPdgID = MotherPdgID;
  fEDep        = EDep;
  fMcMom       = McMom;
  
}

//_____________________________________________________________________________
void TStrawHit::Clear(Option_t* opt) {
  fStrawID     = 0xffffffff;
  fTime[0]     = -1.;
  fTime[1]     = -1.;
  fTOT [0]     = -1.;
  fTOT [1]     = -1.;
  fGenID       = -1;
  fSimID       = -1;
  fPdgID       = -1;
  fMotherPdgID = -1;
  fEDep        = -1.;
  fMcMom       = -1.;
}

//-----------------------------------------------------------------------------
// Options: "banner", "data"
//-----------------------------------------------------------------------------
void TStrawHit::Print(Option_t* Option) const {
  // print Straw hit properties
  //  printf("Superlayer: %d, Wire: %d, Cell: %d,\n",fSuperLayer,fWire,fCell);
  
  TString opt = Option;
  opt.ToLower();

  if ((opt == "") || (opt.Index("banner") >= 0)) {
    printf("-------------------------------------------------------------------------------------------------------------------\n");
    printf("  SID stn face pnl str   T_0      T_1     TOT_0 TOT_1    EDep  Flag MC GenID SimID      PdgID    PdgID_M      McMom\n");
    printf("-------------------------------------------------------------------------------------------------------------------\n");
  }

  if (opt.Index("data") < 0) return;
  
  printf("%5i %2i   %2i  %2i  %2i %9.3f %9.3f %5.0f %5.0f %8.5f 0x%02x %2i %5i %5i %10i %10i %10.3f\n",
	 StrawID(), 
	 Station(),
	 Face(),
	 Panel(),
	 Straw(),
	 fTime[0], fTime[1], // ## format
	 TOT(0)  , TOT(1)  ,
	 fEDep,
	 StrawDigiFlag(),
	 MCFlag(), 
	 GenID(),
	 fSimID,
	 fPdgID,
	 fMotherPdgID,
	 fMcMom);
  
}
