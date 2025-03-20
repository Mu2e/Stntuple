///////////////////////////////////////////////////////////////////////////////
//  2019-06-10 P.Murat TCrvCoincidenceCluster
///////////////////////////////////////////////////////////////////////////////
#include "TString.h"
#include "TBuffer.h"

#include "Stntuple/obj/TCrvCoincidenceCluster.hh"

ClassImp(TCrvCoincidenceCluster)

//_____________________________________________________________________________
void TCrvCoincidenceCluster::Streamer(TBuffer &R__b) {

  int nwi = ((int*) &fStartTime) - &fIndex;
  int nwf = ((float*) &fPosition) - &fStartTime;
  
  if (R__b.IsReading()) {
    //    Version_t R__v = R__b.ReadVersion();
    R__b.ReadVersion();
//-----------------------------------------------------------------------------
// curent version: V2
//-----------------------------------------------------------------------------
    R__b.ReadFastArray(&fIndex    ,nwi);
    R__b.ReadFastArray(&fStartTime,nwf);
    fPosition.Streamer(R__b);
    fMCAvgPosition.Streamer(R__b);
  }
  else {
    R__b.WriteVersion(TCrvCoincidenceCluster::IsA());
    R__b.WriteFastArray(&fIndex    ,nwi);
    R__b.WriteFastArray(&fStartTime,nwf);
    fPosition.Streamer(R__b);
    fMCAvgPosition.Streamer(R__b);
  } 
}

//_____________________________________________________________________________
TCrvCoincidenceCluster::TCrvCoincidenceCluster(): TObject() {
  Clear();
}

//_____________________________________________________________________________
TCrvCoincidenceCluster::~TCrvCoincidenceCluster() {
}

//_____________________________________________________________________________
void TCrvCoincidenceCluster::Set(int Index, int SectorType, int NPulses, int NPe,
				 float X, float Y, float Z, float T1, float T2)
{
  fIndex      = Index;
  fSectorType = SectorType;
  fNPulses    = NPulses;
  fNPe        = NPe;
  fPosition.SetXYZ(X,Y,Z);
  fStartTime  = T1;
  fEndTime    = T2;
}

//_____________________________________________________________________________
void TCrvCoincidenceCluster::SetMC(int SimID, int NPulses, float EnergyDep, float AvgTime,
				 float X, float Y, float Z)
{
  fSimID        = SimID;
  fMCNPulses    = NPulses;
  fMCEnergyDep  = EnergyDep;
  fMCAvgHitTime = AvgTime;
  fMCAvgPosition.SetXYZ(X,Y,Z);
}

//_____________________________________________________________________________
void TCrvCoincidenceCluster::Clear(Option_t* opt) {
  fIndex = -1;
  fSectorType = -1;
  fNPulses = -1;
  fNPe = -1;
  fSimID = -1;
  fMCNPulses = -1;
  for(int i = 0; i < kNFreeInts; ++i) fFreeInts[i] = 0;
  fStartTime = 0.;
  fEndTime = 0.;
  fMCEnergyDep = -1.;
  fMCAvgHitTime = 0.;
  for(int i = 0; i < kNFreeFloats; ++i) fFreeFloats[i] = 0;
  fPosition.SetXYZ(0.,0.,0.);
  fMCAvgPosition.SetXYZ(0.,0.,0.);
}

//_____________________________________________________________________________
void TCrvCoincidenceCluster::Print(Option_t* Option) const {
  // print Straw hit properties
  //  printf("Superlayer: %d, Wire: %d, Cell: %d,\n",fSuperLayer,fWire,fCell);
  
  TString opt = Option;
  opt.ToLower();

  if ((opt == "") || (opt.Index("banner") >= 0)) {
    printf("-------------------------------------------------------------------------------------\n");
    printf("Index  SType   NPulses N(PE) StartTime  EndTime      X         Y         Z     SimID \n");
    printf("-------------------------------------------------------------------------------------\n");
  }

  if (opt == "banner") return;
  
  printf(" %3i ",fIndex);
  printf(" %5i ",fSectorType);
  printf(" %5i ",fNPulses);
  printf(" %5i ",fNPe);
  printf(" %8.2f ",fStartTime);
  printf(" %8.2f ",fEndTime);
  printf(" %8.2f ",fPosition.X());
  printf(" %8.2f ",fPosition.Y());
  printf(" %8.2f ",fPosition.Z());
  printf(" %5i ",fSimID);
  printf("\n");

}
