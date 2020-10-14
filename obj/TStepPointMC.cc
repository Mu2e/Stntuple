//-----------------------------------------------------------------------------
// TStepPointMC: STNTUPLE description of the Mu2e step point
// fMom and fPos are in the beginning of the step
//----------------------------------------------------------------------------- 
#include <limits.h>	// UINT_MAX
#include "TDatabasePDG.h"
#include "Stntuple/obj/TStepPointMC.hh"

ClassImp(TStepPointMC)

//-----------------------------------------------------------------------------
// void TStepPointMC::ReadV1(TBuffer &R__b) {

//   struct TStepPointMCV01_t {
//     int             fParentID;
//     int             fPdgCode;
//     int             fCreationCode;
//     int             fStartVolumeIndex;
//     int             fTerminationCode;
//     int             fEndVolumeIndex;
//     int             fNStrawHits;
    
//     float           fMomTargetEnd;
//     float           fMomTrackerFront;		// entrance to ST

//     TLorentzVector  fStartPos;
//     TLorentzVector  fStartMom;
//   };

//   TStepPointMCV01_t data;

//   int nwi = ((int*  ) &data.fMomTargetEnd) - &data.fParentID;
//   int nwf = ((float*) &data.fStartPos    ) - &data.fMomTargetEnd ;

//   TObject::Streamer(R__b);

//   R__b.ReadFastArray(&data.fParentID   ,nwi);
//   R__b.ReadFastArray(&data.fMomTargetEnd,nwf);

//   fParentID         = data.fParentID;
//   fPdgCode          = data.fPdgCode;
//   fCreationCode     = data.fCreationCode;
//   fStartVolumeIndex = data.fStartVolumeIndex;
//   fTerminationCode  = data.fTerminationCode;
//   fEndVolumeIndex   = data.fEndVolumeIndex;
//   fNStrawHits       = data.fNStrawHits;

//   fGenpID           = -1;         // ** added in V2 **

//   fStartPos.Streamer(R__b);
//   fStartMom.Streamer(R__b);
  
// }

//-----------------------------------------------------------------------------
// we don't really need to write out TObject part - so far it is not used
//-----------------------------------------------------------------------------
void TStepPointMC::Streamer(TBuffer& R__b) {
  int nwi, nwf;

  nwi = ((int*  ) &fEDepTot) - &fVolumeID;
  nwf = ((float*) &fPos    ) - &fEDepTot ;

  if (R__b.IsReading()) {
    //    Version_t R__v = R__b.ReadVersion(); 
    R__b.ReadVersion(); 
    // if (R__v == 1) ReadV1(R__b);
    // else {

    TObject::Streamer(R__b);
    R__b.ReadFastArray(&fVolumeID,nwi);
    R__b.ReadFastArray(&fEDepTot ,nwf);
    
    fPos.Streamer(R__b);
    fMom.Streamer(R__b);
    //    }
  }
  else {
    R__b.WriteVersion(TStepPointMC::IsA());
    TObject::Streamer(R__b);

    R__b.WriteFastArray(&fVolumeID,nwi);
    R__b.WriteFastArray(&fEDepTot ,nwf);

    fPos.Streamer(R__b);
    fMom.Streamer(R__b);
  }
}

//_____________________________________________________________________________
TStepPointMC::TStepPointMC() {
}

//_____________________________________________________________________________
TStepPointMC::TStepPointMC(int VolumeID    , int GenIndex      , 
			   int SimID       , int PDGCode       ,
			   int ParentSimID , int ParentPDGCode , 
			   int CreationCode, int EndProcessCode, 
			   float EDepTot   , float EDepNio     , 
			   float Time      , float StepLength  ,
			   float X, float Y, float Z, float Px, float Py, float Pz):
  fPos( X, Y, Z),
  fMom(Px,Py,Pz)
{
  fVolumeID       = VolumeID;
  fGenIndex       = GenIndex;
  fSimID          = SimID;
  fPDGCode        = PDGCode;
  fParentSimID    = ParentSimID;
  fParentPDGCode  = ParentPDGCode;
  
  fCreationCode   = CreationCode;
  fEndProcessCode = EndProcessCode;

  fEDepTot        = EDepTot;
  fEDepNio        = EDepNio;
  fTime           = Time;
  fStepLength     = StepLength;
}

//-----------------------------------------------------------------------------
TStepPointMC::~TStepPointMC() {
}

//_____________________________________________________________________________
void TStepPointMC::Print(Option_t* Opt) const {
  TString opt = Opt;
  if ((opt == "banner") || (opt == "")) {
    // print banner
    printf("---------------------------------------------------------------------------------------------------------------");
    printf("-------------------------------------------------------------------------\n");
    printf("   Vol   Gen     PDG Particle          Creation   SimID   PPdg  PSimID  StopProc     X          Y          Z   ");
    printf("  Edep(Tot) Edep(NIO)  Step    Time      Px        Py        Pz     Ptot\n");
    printf("---------------------------------------------------------------------------------------------------------------");
    printf("------------------------------------------------------------------------\n");
  }

  TDatabasePDG* db = TDatabasePDG::Instance();

  TParticlePDG* pdg = db->GetParticle(fPDGCode);

  if ((opt == "data") || (opt == "")) {
     printf("%5i",VolumeID());
     printf("%8i"  ,fGenIndex);

     printf("%7i"  ,fPDGCode);
     if (pdg) printf(" %-19s",pdg->GetName());
     else          printf(" %-19s","*** unknown ***");

     printf("%7i"   ,fCreationCode);
     printf("%8i"   ,fSimID);
     printf("%7i"   ,fParentPDGCode);
     printf("%8i"   ,fParentSimID);
     printf("%8i"   ,fEndProcessCode);
     printf("%11.3f",fPos.X());
     printf("%11.3f",fPos.Y());
     printf("%11.3f",fPos.Z());

     printf("%9.5f" ,fEDepTot);
     printf("%9.5f" ,fEDepNio);

     printf("%9.3f" ,fStepLength);
     printf("%9.3f" ,fTime);
     
     printf("%9.3f" ,fMom.Px());
     printf("%9.3f" ,fMom.Py());
     printf("%9.3f" ,fMom.Pz());
     printf("%9.3f" ,fMom.Mag());
     printf("\n");
  }
}
