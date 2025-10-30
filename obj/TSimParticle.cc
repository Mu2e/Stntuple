//-----------------------------------------------------------------------------
// TSimParticle: STNTUPLE description of the generator-level MC particle
// It inherits from ROOT's TParticle and adds more functions to it
// particle number in the original list is saved in TObject::fUniqueID 
// Nov 25 2000 P.Murat (CDF/FNAL)
//----------------------------------------------------------------------------- 
#include <limits.h>	// UINT_MAX
#include "TDatabasePDG.h"
#include "Stntuple/obj/TSimParticle.hh"

ClassImp(TSimParticle)

//-----------------------------------------------------------------------------
void TSimParticle::ReadV1(TBuffer &R__b) {

  struct TSimParticleV01_t {
    int             fParentID;
    int             fPdgCode;
    int             fCreationCode;
    int             fStartVolumeIndex;
    int             fTerminationCode;
    int             fEndVolumeIndex;
    int             fNStrawHits;
    
    float           fMomTargetEnd;
    float           fMomTrackerFront;		// entrance to ST

    TLorentzVector  fStartPos;
    TLorentzVector  fStartMom;
  };

  TSimParticleV01_t data;

  int nwi = ((int*  ) &data.fMomTargetEnd) - &data.fParentID;
  int nwf = ((float*) &data.fStartPos    ) - &data.fMomTargetEnd ;

  TObject::Streamer(R__b);

  R__b.ReadFastArray(&data.fParentID   ,nwi);
  R__b.ReadFastArray(&data.fMomTargetEnd,nwf);

  fParentID         = data.fParentID;
  fPdgCode          = data.fPdgCode;
  fCreationCode     = data.fCreationCode;
  fStartVolumeIndex = data.fStartVolumeIndex;
  fTerminationCode  = data.fTerminationCode;
  fEndVolumeIndex   = data.fEndVolumeIndex;
  fNStrawHits       = data.fNStrawHits;

  fGeneratorID         = -1;             // ** added in V2

  fStartPos.Streamer(R__b);
  fStartMom.Streamer(R__b);
  
  fEndPos.SetXYZT(0.,0.,0.,0.);	         // ** added in V3
  fEndMom.SetXYZT(0.,0.,0.,0.);	         // ** added in V3

  fStartProperTime = -1.;                // ** added in v4
  fEndProperTime   = -1.;                // ** added in v4
}

//-----------------------------------------------------------------------------
void TSimParticle::ReadV2(TBuffer &R__b) {

  struct TSimParticleV02_t {
    int             fParentID;
    int             fPdgCode;
    int             fCreationCode;
    int             fStartVolumeIndex;
    int             fTerminationCode;
    int             fEndVolumeIndex;
    int             fNStrawHits;
    int             fGeneratorID;               // ** MC generator ID, added in V2
    
    float           fMomTargetEnd;
    float           fMomTrackerFront;		// entrance to ST

    TLorentzVector  fStartPos;
    TLorentzVector  fStartMom;
  };

  TSimParticleV02_t data;

  int nwi = ((int*  ) &data.fMomTargetEnd) - &data.fParentID;
  int nwf = ((float*) &data.fStartPos    ) - &data.fMomTargetEnd ;

  TObject::Streamer(R__b);

  R__b.ReadFastArray(&data.fParentID   ,nwi);
  R__b.ReadFastArray(&data.fMomTargetEnd,nwf);

  fParentID         = data.fParentID;
  fPdgCode          = data.fPdgCode;
  fCreationCode     = data.fCreationCode;
  fStartVolumeIndex = data.fStartVolumeIndex;
  fTerminationCode  = data.fTerminationCode;
  fEndVolumeIndex   = data.fEndVolumeIndex;
  fNStrawHits       = data.fNStrawHits;
  fGeneratorID      = data.fGeneratorID; // ** added in V2

  fStartPos.Streamer(R__b);
  fStartMom.Streamer(R__b);

  fEndPos.SetXYZT(0.,0.,0.,0.);		 // ** added in V3
  fEndMom.SetXYZT(0.,0.,0.,0.);		 // ** added in V3

  fStartProperTime = -1.;                // ** added in v4
  fEndProperTime   = -1.;                // ** added in v4
}

//-----------------------------------------------------------------------------
void TSimParticle::ReadV3(TBuffer &R__b) {

  struct TSimParticleV03_t {
    int             fParentID;
    int             fPdgCode;
    int             fCreationCode;
    int             fStartVolumeIndex;
    int             fTerminationCode;
    int             fEndVolumeIndex;
    int             fNStrawHits;
    int             fGeneratorID;               // ** MC generator ID, added in V2
    
    float           fMomTargetEnd;
    float           fMomTrackerFront;		// entrance to ST

    TLorentzVector  fStartPos;
    TLorentzVector  fStartMom;
    TLorentzVector  fEndPos;                    // ** added in V3
    TLorentzVector  fEndMom;                    // ** added in V3
  } data;

  int nwi = ((int*  ) &data.fMomTargetEnd) - &data.fParentID;
  int nwf = ((float*) &data.fStartPos    ) - &data.fMomTargetEnd ;

  TObject::Streamer(R__b);

  R__b.ReadFastArray(&data.fParentID   ,nwi);
  R__b.ReadFastArray(&data.fMomTargetEnd,nwf);

  fParentID         = data.fParentID;
  fPdgCode          = data.fPdgCode;
  fCreationCode     = data.fCreationCode;
  fStartVolumeIndex = data.fStartVolumeIndex;
  fTerminationCode  = data.fTerminationCode;
  fEndVolumeIndex   = data.fEndVolumeIndex;
  fNStrawHits       = data.fNStrawHits;
  fGeneratorID      = data.fGeneratorID;         // ** added in V2 **

  fStartPos.Streamer(R__b);
  fStartMom.Streamer(R__b);
  fEndPos.Streamer(R__b);                        // ** added in V3
  fEndMom.Streamer(R__b);                        // ** added in V3

  fStartProperTime = -1.;                        // ** added in v4
  fEndProperTime   = -1.;                        // ** added in v4
}

//-----------------------------------------------------------------------------
// we don't really need to write out TObject part - so far it is not used
//-----------------------------------------------------------------------------
void TSimParticle::Streamer(TBuffer& R__b) {
  int nwi, nwf;

  nwi = ((int*  ) &fMomTargetEnd) - &fParentID;
  nwf = ((float*) &fStartPos    ) - &fMomTargetEnd ;

  if (R__b.IsReading()) {
    Version_t R__v = R__b.ReadVersion(); 
    if (R__v == 1) ReadV1(R__b);
    if (R__v == 2) ReadV2(R__b);
    if (R__v == 3) ReadV3(R__b);
    else {

      TObject::Streamer(R__b);
      R__b.ReadFastArray(&fParentID   ,nwi);
      R__b.ReadFastArray(&fMomTargetEnd,nwf);

      fStartPos.Streamer(R__b);
      fStartMom.Streamer(R__b);
      fEndPos.Streamer  (R__b);
      fEndMom.Streamer  (R__b);
    }
  }
  else {
    R__b.WriteVersion(TSimParticle::IsA());
    TObject::Streamer(R__b);

    R__b.WriteFastArray(&fParentID    ,nwi);
    R__b.WriteFastArray(&fMomTargetEnd,nwf);

    fStartPos.Streamer(R__b);
    fStartMom.Streamer(R__b);
    fEndPos.Streamer  (R__b);
    fEndMom.Streamer  (R__b);
  }
}

//_____________________________________________________________________________
TSimParticle::TSimParticle() {
  SetUniqueID(UINT_MAX);
  fGeneratorID   = -1;
  fShid          = nullptr;
  fSimParticle   = nullptr;
}

//_____________________________________________________________________________
TSimParticle::TSimParticle(Int_t ID, Int_t ParentID, Int_t PDGCode, 
			   int CreationCode, int TerminationCode,
			   int StartVolumeIndex, int EndVolumeIndex,
			   int GeneratorID): TObject()
  
{
  Init(ID,ParentID,PDGCode,CreationCode,TerminationCode,StartVolumeIndex,EndVolumeIndex,GeneratorID);
}

//-----------------------------------------------------------------------------
// fSimp is just cached
//-----------------------------------------------------------------------------
TSimParticle::~TSimParticle() {
  if (fShid) {
    delete fShid;
    fShid = nullptr;
  }
}

//_____________________________________________________________________________
int  TSimParticle::Init(Int_t ID, Int_t ParentID, Int_t PDGCode, 
			int CreationCode, int TerminationCode,
			int StartVolumeIndex, int EndVolumeIndex,
			int GeneratorID) 
{
  SetUniqueID(ID);
  fParentID         = ParentID;
  fPdgCode          = PDGCode;
  fCreationCode     = CreationCode;
  fTerminationCode  = TerminationCode;
  fStartVolumeIndex = StartVolumeIndex;
  fEndVolumeIndex   = EndVolumeIndex;
  fNStrawHits       = 0;
  fGeneratorID      = GeneratorID;
  fMomTargetEnd     = -1.;
  fMomTrackerFront  = -1.;
  fSimParticle      = nullptr;
  fShid             = nullptr;

  return 0;
}

//_____________________________________________________________________________
void TSimParticle::Clear(Option_t* Opt) {
  fSimParticle = nullptr;
  if (fShid) {
    delete fShid;
    fShid = nullptr;
  }
}


//_____________________________________________________________________________
void TSimParticle::Delete(Option_t* Opt) {
  fSimParticle = nullptr;
  if (fShid) {
    delete fShid;
    fShid = nullptr;
  }
}


//_____________________________________________________________________________
void TSimParticle::Print(Option_t* Opt) const {

  TString opt = Opt;
  if ((opt.Index("banner") >= 0) || (opt == "")) {
				// print banner
    printf("-----------------------------------------------------");
    printf("------------------------------------------------------------------------------------------");
    printf("-----------------------------------------------------------------------------------------------\n");
    printf("   i name                      PDG  ID/Stg GenID Prnt");
    printf("      p0x        p0y        p0z       p0     vol0     v0x        v0y      v0z      t0   t0p/tau");
    printf("      p1x        p1y        p1z       p1     vol1     v1x        v1y      v1z      t1   t1p/tau Nsh\n");
    printf("-----------------------------------------------------");
    printf("------------------------------------------------------------------------------------------");
    printf("-----------------------------------------------------------------------------------------------\n");
  }

  TDatabasePDG* db = TDatabasePDG::Instance();

  TParticlePDG* pdg = db->GetParticle(fPdgCode);

  if ((opt.Index("data") >= 0) || (opt == "")) {
    printf("%4i",Number());

    if (pdg) printf(" %-19s",pdg->GetName());
    else          printf(" %-19s","*** unknown ***");

    printf("%10i"  ,fPdgCode);
    printf("%4i/%2i",GetUniqueID(),SimStage());
    printf("%6i"   ,fGeneratorID);
    printf("%6i"   ,fParentID);
    printf("%11.3f",fStartMom.Px());
    printf("%11.3f",fStartMom.Py());
    printf("%11.3f",fStartMom.Pz());
    printf("%10.3f",fStartMom.P());
    printf("%6i"   ,fStartVolumeIndex);
    printf("%10.3f",fStartPos.X());
    printf("%10.3f",fStartPos.Y());
    printf("%10.3f",fStartPos.Z());
    printf("%10.3e",fStartPos.T());
    printf("%6.2f" ,fStartProperTime);
    printf("%11.3f",fEndMom.Px());
    printf("%11.3f",fEndMom.Py());
    printf("%11.3f",fEndMom.Pz());
    printf("%10.3f",fEndMom.P());
    printf("%6i"   ,fEndVolumeIndex);
    printf("%10.3f",fEndPos.X());
    printf("%10.3f",fEndPos.Y());
    printf("%10.3f",fEndPos.Z());
    printf("%10.3e",fEndPos.T());
    printf("%6.2f" ,fEndProperTime);
    printf("%4i"   ,NStrawHits());
    printf("\n");
  }
}

