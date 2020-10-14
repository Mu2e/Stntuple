///////////////////////////////////////////////////////////////////////////////
// use of tmp:
//
// Tmp(0) : nax seg
// Tmp(1) : nst seg
// 
// use of debug bits: bits 0-2 are reserved
// 0  : all events
// 1  : passed events
// 2  : rejected events
// 
// 3  : events with SPMC T<100
// 4  : events with pbars ID=300000+I, X> 25cm
// 5  : events with pbars ID=400000+I
// 6  : events with pbars reaching the final stage
// 7  : events with pbars P > 100 MeV/c reaching the final stage
// 8  : events with pi (-211) in VD9 (before the target)
// 9  : events with pi (-211) at stage 4 (case of pbar simulation)
///////////////////////////////////////////////////////////////////////////////
#include "TF1.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TEnv.h"
#include "TSystem.h"

#include "Stntuple/loop/TStnAna.hh"
#include "Stntuple/obj/TStnHeaderBlock.hh"
#include "Stntuple/val/stntuple_val_functions.hh"
//------------------------------------------------------------------------------
// Mu2e offline includes
//-----------------------------------------------------------------------------
// #include "CalorimeterGeom/inc/HexMap.hh"

#include "ana/TStnSpmcAnaModule.hh"

#include "ana/StnInitVirtualDetectors.hh"

ClassImp(TStnSpmcAnaModule)
//-----------------------------------------------------------------------------
TStnSpmcAnaModule::TStnSpmcAnaModule(const char* name, const char* title):
  TStnModule(name,title)
{
  // fPdgCode       = 11;
  // fGeneratorCode = 28;
  fSpmcBlockName = "SpmcTrackerBlock";
  fVDetBlockName = "SpmcVDetBlock";

  fPdgDb = TDatabasePDG::Instance();

  for (int i=0; i<5000; i++) fParticleCache[i] = 0;

  SetParticleCache(   11,fPdgDb->GetParticle(   11)); // e
  SetParticleCache(  -11,fPdgDb->GetParticle(  -11)); // e
  SetParticleCache(   13,fPdgDb->GetParticle(   13)); // mu
  SetParticleCache(  -13,fPdgDb->GetParticle(  -13)); // mu
  SetParticleCache(   22,fPdgDb->GetParticle(   22)); // photon
  SetParticleCache(  211,fPdgDb->GetParticle(  211)); // pi^-
  SetParticleCache( -211,fPdgDb->GetParticle( -211)); // pi^-
  SetParticleCache( 2112,fPdgDb->GetParticle( 2112)); // neutron
  SetParticleCache( 2212,fPdgDb->GetParticle( 2212)); // proton
  SetParticleCache(-2212,fPdgDb->GetParticle(-2212)); // pbar

  fStnt = TStntuple::Instance();
}

//-----------------------------------------------------------------------------
TStnSpmcAnaModule::~TStnSpmcAnaModule() {
}


//-----------------------------------------------------------------------------
void TStnSpmcAnaModule::BookEventHistograms(EventHist_t* Hist, const char* Folder) {
  //  char name [200];
  //  char title[200];

  HBook1F(Hist->fEventNumber,"evtnum"  ,Form("%s: Event Number"        ,Folder), 1000, 0,  1.e4,Folder);
  HBook1F(Hist->fRunNumber  ,"runnum"  ,Form("%s: Run   Number"        ,Folder), 1000, 0,  1.e6,Folder);
  HBook1F(Hist->fNSimp      ,"nsimp"   ,Form("%s: N(sim particles)"    ,Folder),  200, 0,  200 ,Folder);
  HBook1F(Hist->fNSpmc      ,"nspmc"   ,Form("%s: N(step points)"      ,Folder),  500, 0,  500 ,Folder);
  HBook1F(Hist->fNSpmc150   ,"nspmc150",Form("%s: N(step points T>150)",Folder),  200, 0,  200 ,Folder);
}

//-----------------------------------------------------------------------------
void TStnSpmcAnaModule::BookStepPointMCHistograms(SpmcHist_t* Hist, const char* Folder) {
  //     char name [200];
  //     char title[200];

  HBook1F(Hist->fVolumeID,"vol_id"   ,Form("%s: VolumeID"       ,Folder), 200,    0, 10000,Folder);
  HBook1F(Hist->fGenIndex,"gen_index",Form("%s: GenIndex"       ,Folder), 200,    0, 10000,Folder);
  HBook1F(Hist->fSimID   ,"sim_id"   ,Form("%s: SimID"          ,Folder), 200,    0,  1000,Folder);

  HBook1F(Hist->fPDGCode[0] ,"pdg_0" ,Form("%s: PDG code"       ,Folder),  500,  -250,   250,Folder);
  HBook1F(Hist->fPDGCode[1] ,"pdg_1" ,Form("%s: PDG code"       ,Folder), 2000, -10000, 10000,Folder);

  HBook1F(Hist->fCreationCode   ,"cr_code",Form("%s: Creation code",Folder), 200,   0,   200,Folder);
  HBook1F(Hist->fParentSimID    ,"psim_id",Form("%s: Parent SimID",Folder), 200,   0,  1000,Folder);
  HBook1F(Hist->fParentPDGCode  ,"ppdg"   ,Form("%s: Parent PDG code",Folder), 200, -100,   100,Folder);
  HBook1F(Hist->fEndProcessCode ,"end_code",Form("%s: End process code",Folder), 200,   0,   200,Folder);
  HBook1F(Hist->fEDepTot        ,"edep_tot",Form("%s: EDEP tot"        ,Folder), 200,   0,   10 ,Folder);
  HBook1F(Hist->fEDepNio        ,"edep_nio",Form("%s: EDEP NIO"        ,Folder), 200,   0,   10 ,Folder);
  HBook1F(Hist->fTime           ,"time"    ,Form("%s: Time"            ,Folder), 200,   0,  2000,Folder);
  HBook1F(Hist->fStepLength     ,"step"    ,Form("%s: Ltep Length"     ,Folder), 200,   0,   100,Folder);
  HBook1F(Hist->fMom[0]         ,"mom"     ,Form("%s: Momentum[0]"     ,Folder), 500,   0,   500,Folder);
  HBook1F(Hist->fMom[1]         ,"mom_1"   ,Form("%s: Momentum[1]"     ,Folder), 500,   0,  5000,Folder);

  HBook2F(Hist->fCosThVsMom     ,"cth_vs_mom",Form("%s: Cos(Th) vs Mom",Folder), 250,   0,  5000,100,-1,1,Folder);

  HBook1F(Hist->fEKin           ,"ekin"    ,Form("%s: kinetic energy"  ,Folder), 400,   0,   100,Folder);

  HBook2F(Hist->fYVsX           ,"y_vs_x"     ,Form("%s: Y vs X"       ,Folder), 100, -250,  250, 100, -250, 250, Folder);
  HBook2F(Hist->fYVsZ           ,"y_vs_z"     ,Form("%s: Y vs Z"       ,Folder), 500, -250,  250, 500, -250, 250, Folder);

  HBook2F(Hist->fCosThVsMomPV   ,"cth_vs_mom_pv",Form("%s: Cos(Th):Mom PV" ,Folder), 250,   0,  5000,100,-1,1,Folder);
}

//-----------------------------------------------------------------------------
void TStnSpmcAnaModule::BookSimpHistograms(SimpHist_t* Hist, const char* Folder) {

  HBook1F(Hist->fVolumeID   ,"vol_id"   ,Form("%s: Volume ID"   ,Folder),1000,  2500, 3500,Folder);
  HBook1F(Hist->fStage      ,"stage"    ,Form("%s: Stage"       ,Folder),  10,     0,   10,Folder);
  HBook1F(Hist->fGeneratorID,"gen_id"   ,Form("%s: Generator ID",Folder), 200,   -10,  190,Folder);
  HBook1F(Hist->fTime       ,"time"     ,Form("%s: Stop Time"   ,Folder), 200,     0, 2000,Folder);
  HBook1F(Hist->fParentMom  ,"pmom"     ,Form("%s: Parent Mom"  ,Folder), 200,     0, 2000,Folder);
  HBook1F(Hist->fParentPDG  ,"ppdg"     ,Form("%s: Parent PDG"  ,Folder), 200, -1000, 1000,Folder);

  HBook1F(Hist->fStartMom[0],"mom"        ,Form("%s: start Mom[0]"  ,Folder), 500,     0,  500,Folder);
  HBook1F(Hist->fStartMom[1],"mom_1"      ,Form("%s: start Mom[1]"  ,Folder), 500,     0, 5000,Folder);
  HBook2F(Hist->fYVsX       ,"y_vs_x"     ,Form("%s: yend vs Xend " ,Folder), 250,  -250, 250, 250, -250, 250,Folder);
  HBook2F(Hist->fXEndVsZEnd ,"xe_vs_ze"   ,Form("%s: xend vs zend " ,Folder), 250,  -5000, 20000, 100, -5000, 5000,Folder);
  HBook2F(Hist->fYVsX_2480  ,"y_vs_x_2480",Form("%s: Y vs X [2480]" ,Folder), 250,  -250, 250, 250, -250, 250,Folder);
  HBook2F(Hist->fYVsX_2513  ,"y_vs_x_2513",Form("%s: Y vs X [2513]" ,Folder), 250,  -250, 250, 250, -250, 250,Folder);


  HBook2F(Hist->fCosThVsMom ,"cth_vs_mom",Form("%s: Cos(Th) vs Mom",Folder), 250,   0,  5000,100,-1,1,Folder);

  HBook2F(Hist->fCosThVsMomPV   ,"cth_vs_mom_pv",Form("%s: Cos(Th):Mom PV" ,Folder), 250,   0,  5000,100,-1,1,Folder);
}

//-----------------------------------------------------------------------------
void TStnSpmcAnaModule::BookVDetHistograms(VDetHist_t* Hist, const char* Folder) {

  HBook1F(Hist->fIndex   ,"index"   ,Form("%s: VD index"      ,Folder),1000, 0, 1000,Folder);
  HBook1F(Hist->fPDGCode ,"pdg_code",Form("%s: PDG code"      ,Folder),2000,-1000, 1000,Folder);
  HBook1F(Hist->fGenCode ,"gen_code",Form("%s: generator code",Folder), 100, -10, 90,Folder);
  HBook1F(Hist->fMom[0]  ,"mom"  ,Form("%s: Momentum[0]"   ,Folder), 500, 0, 500,Folder);
  HBook1F(Hist->fMom[1]  ,"mom_1",Form("%s: Momentum[1]"   ,Folder), 500, 0,5000,Folder);
  HBook1F(Hist->fTime    ,"time"    ,Form("%s: Hit Time  "    ,Folder), 200, 0,2000,Folder);
  HBook2F(Hist->fYVsX    ,"y_vs_x"  ,Form("%s: Y vs X (all)"  ,Folder), 250, -250, 250, 250, -250, 250,Folder);
  HBook2F(Hist->fYVsZ    ,"y_vs_z"  ,Form("%s: Y vs Z (all)"  ,Folder), 250, -250, 250, 250, -250, 250,Folder);
  HBook1F(Hist->fPt      ,"pt"      ,Form("%s: Pt"            ,Folder), 200, 0, 200,Folder);
  HBook1F(Hist->fPp      ,"pp"      ,Form("%s: P(parallel)"   ,Folder), 200, 0, 200,Folder);
  HBook1F(Hist->fTanTh   ,"tan_th"  ,Form("%s: tan(pitch ang)",Folder), 500, -1, 4,Folder);
  HBook1F(Hist->fEKin    ,"ekin"    ,Form("%s: kinetic energy",Folder), 400,  0, 100,Folder);

  HBook2F(Hist->fCosThVsMom,"cth_vs_mom",Form("%s: Cos(Th) vs Mom",Folder), 250,   0,  5000,100,-1,1,Folder);

  HBook2F(Hist->fCosThVsMomPV   ,"cth_vs_mom_pv",Form("%s: Cos(Th):Mom PV" ,Folder), 250,   0,  5000,100,-1,1,Folder);
}

//-----------------------------------------------------------------------------
void TStnSpmcAnaModule::BookHistograms() {

  //  char name [200];
  //  char title[200];

  TFolder* fol;
  TFolder* hist_folder;
  char     folder_name[200];

  DeleteHistograms();
  hist_folder = (TFolder*) GetFolder()->FindObject("Hist");
//-----------------------------------------------------------------------------
// book event histograms
//-----------------------------------------------------------------------------
  int book_event_histset[kNEventHistSets];
  for (int i=0; i<kNEventHistSets; i++) book_event_histset[i] = 0;

  book_event_histset[ 0] = 1;		// all events
  book_event_histset[ 1] = 1;		// events with at least one T>150 step point

  for (int i=0; i<kNEventHistSets; i++) {
    if (book_event_histset[i] != 0) {
      sprintf(folder_name,"evt_%i",i);
      fol = (TFolder*) hist_folder->FindObject(folder_name);
      if (! fol) fol = hist_folder->AddFolder(folder_name,folder_name);
      fHist.fEvent[i] = new EventHist_t;
      BookEventHistograms(fHist.fEvent[i],Form("Hist/%s",folder_name));
    }
  }
//-----------------------------------------------------------------------------
// book StepPointMC histograms
//-----------------------------------------------------------------------------
  int book_spmc_histset[kNStepPointMCHistSets];
  for (int i=0; i<kNStepPointMCHistSets; i++) book_spmc_histset[i] = 0;

  book_spmc_histset[0] = 1;		// all steps
  book_spmc_histset[1] = 1;		// electrons
  book_spmc_histset[2] = 1;		// positrons
  book_spmc_histset[3] = 1;		// mu-
  book_spmc_histset[4] = 1;		// mu+
  book_spmc_histset[5] = 1;		// photons
  book_spmc_histset[6] = 1;		// negative pions
  book_spmc_histset[7] = 1;		// positive pions
  book_spmc_histset[8] = 1;		// protons+antiprotons
  book_spmc_histset[9] = 1;		// everything else

  book_spmc_histset[ 900] = 1;		// all particles , VD=9
  book_spmc_histset[ 913] = 1;		// mu-   in VD=9
  book_spmc_histset[ 921] = 1;		// pbars in VD=9

  book_spmc_histset[9100] = 1;		// all particles , VD=91
  book_spmc_histset[9113] = 1;		// mu- in VD=91
  book_spmc_histset[9121] = 1;		// pbars in VD=91

  for (int i=0; i<kNStepPointMCHistSets; i++) {
    if (book_spmc_histset[i] != 0) {
      sprintf(folder_name,"spmc_%i",i);
      fol = (TFolder*) hist_folder->FindObject(folder_name);
      if (! fol) fol = hist_folder->AddFolder(folder_name,folder_name);
      fHist.fStepPointMC[i] = new SpmcHist_t;
      BookStepPointMCHistograms(fHist.fStepPointMC[i],Form("Hist/%s",folder_name));
    }
  }
//-----------------------------------------------------------------------------
// book SimParticle histograms ??  none so far
//-----------------------------------------------------------------------------
  int book_simp_histset[kNSimpHistSets];
  for (int i=0; i<kNSimpHistSets; i++) book_simp_histset[i] = 0;

  book_simp_histset[  0] = 1;		// all particles 0+ID
  book_simp_histset[  3] = 1;		// mu-
  book_simp_histset[  4] = 1;		// mu+
  book_simp_histset[  5] = 1;		// pi-
  book_simp_histset[ 21] = 1;		// pbars

  book_simp_histset[100] = 1;		// stage 0  (offset 0)
  book_simp_histset[103] = 1;		// mu-
  book_simp_histset[104] = 1;		// mu+
  book_simp_histset[105] = 1;		// pi-
  book_simp_histset[121] = 1;		// pbars

  book_simp_histset[200] = 1;		// stage 1  (offset 100000)
  book_simp_histset[203] = 1;		// mu-
  book_simp_histset[204] = 1;		// mu+
  book_simp_histset[205] = 1;		// pi-
  book_simp_histset[221] = 1;		// pbars

  book_simp_histset[300] = 1;		// stage 2
  book_simp_histset[303] = 1;		// mu-
  book_simp_histset[304] = 1;		// mu+
  book_simp_histset[305] = 1;		// pi-
  book_simp_histset[321] = 1;		// pbars

  book_simp_histset[400] = 1;		// stage 3
  book_simp_histset[403] = 1;		// mu-
  book_simp_histset[404] = 1;		// mu+
  book_simp_histset[405] = 1;		// pi-
  book_simp_histset[421] = 1;		// pbars

  book_simp_histset[500] = 1;		// stage 4
  book_simp_histset[503] = 1;		// mu-
  book_simp_histset[504] = 1;		// mu+
  book_simp_histset[505] = 1;		// pi-
  book_simp_histset[521] = 1;		// pbars

  for (int i=0; i<kNSimpHistSets; i++) {
    if (book_simp_histset[i] != 0) {
      sprintf(folder_name,"simp_%i",i);
      fol = (TFolder*) hist_folder->FindObject(folder_name);
      if (! fol) fol = hist_folder->AddFolder(folder_name,folder_name);
      fHist.fSimp[i] = new SimpHist_t;
      BookSimpHistograms(fHist.fSimp[i],Form("Hist/%s",folder_name));
    }
  }
//-----------------------------------------------------------------------------
// book VDet histograms - so far need only VDET=9 (in front of the stopping target)
// virtual detector ID's are defined in DataProducts/inc/VirtualDetectorId.hh
//-----------------------------------------------------------------------------
  int book_vdet_histset[kNVDetHistSets];
  for (int i=0; i<kNVDetHistSets; i++) book_vdet_histset[i] = 0;

  book_vdet_histset[  9] = 1;		// all particles, VDET= 9 , ST_In
  book_vdet_histset[ 13] = 1;		// all particles, VDET=13 , TT_FrontHollow
  book_vdet_histset[ 14] = 1;		// all particles, VDET=13 , TT_FrontHollow, r > 40

  for (int i=0; i<kNVDetHistSets; i++) {
    if (book_vdet_histset[i] != 0) {
      sprintf(folder_name,"vdet_%i",i);
      fol = (TFolder*) hist_folder->FindObject(folder_name);
      if (! fol) fol = hist_folder->AddFolder(folder_name,folder_name);
      fHist.fVDet[i] = new VDetHist_t;
      BookVDetHistograms(fHist.fVDet[i],Form("Hist/%s",folder_name));
    }
  }
}

//-----------------------------------------------------------------------------
void TStnSpmcAnaModule::FillEventHistograms(EventHist_t* Hist) {
//   double            cos_th, xv, yv, rv, zv, p;
//   TLorentzVector    mom;

  int event_number = GetHeaderBlock()->EventNumber();
  int run_number   = GetHeaderBlock()->RunNumber();

  Hist->fEventNumber->Fill(event_number);
  Hist->fRunNumber->Fill(run_number);
  Hist->fNSimp->Fill(fNSimp);
  Hist->fNSpmc->Fill(fNSpmc);
  Hist->fNSpmc150->Fill(fNSpmc150);
}

//-----------------------------------------------------------------------------
void TStnSpmcAnaModule::FillStepPointMCHistograms(SpmcHist_t* Hist, TStepPointMC* Step, SpmcData_t* SpmcData, double Weight) {

  int pdg_code = Step->PDGCode();

  int id = (Step->VolumeID());

  StnVDetData_t* vdd = fVDet+id;

  Hist->fVolumeID->Fill(id,Weight);
  Hist->fGenIndex->Fill(Step->GenIndex(),Weight);
  Hist->fSimID   ->Fill(Step->SimID(),Weight);
  Hist->fPDGCode[0] ->Fill(pdg_code,Weight);
  Hist->fPDGCode[1] ->Fill(pdg_code,Weight);
  Hist->fCreationCode->Fill(Step->CreationCode(),Weight);
  Hist->fParentSimID ->Fill(Step->ParentSimID(),Weight);
  Hist->fParentPDGCode->Fill(Step->ParentPDGCode(),Weight);
  Hist->fEndProcessCode->Fill(Step->EndProcessCode(),Weight);

  Hist->fEDepTot->Fill(Step->EDepTot(),Weight);
  Hist->fEDepNio->Fill(Step->EDepNio(),Weight);
  Hist->fTime   ->Fill(Step->Time(),Weight);
  Hist->fStepLength->Fill(Step->StepLength(),Weight);

  double p = Step->Mom()->Mag();
  Hist->fMom[0]->Fill(p,Weight);
  Hist->fMom[1]->Fill(p,Weight);
  
  double m(0);
  if (SpmcData->fParticle) {
    m = SpmcData->fParticle->Mass()*1e3;   // convert GeV -> MeV
  }

  double ekin = sqrt(p*p+m*m)-m;
  Hist->fEKin->Fill(ekin,Weight);

  float x = Step->Pos()->X();
  float y = Step->Pos()->Y();
  float z = Step->Pos()->Z();

  // hack for stage 2:

  if (fabs(z-2929.) < 0.1) x = x+3904;

  Hist->fYVsX->Fill(x,y,Weight);		// useful for stage 2
  Hist->fYVsZ->Fill(z,y,Weight);		// useful for stage 1

  float pp;
  if (vdd->fIZLocal == 3) pp = Step->Mom()->Pz();
  else                    pp = -Step->Mom()->Px();

  float cos_th = pp/p;
  Hist->fCosThVsMom->Fill(p,cos_th,Weight);

  Hist->fCosThVsMomPV->Fill(fPbarMomPV,fPbarCosThPV,Weight);
}

//-----------------------------------------------------------------------------
void TStnSpmcAnaModule::FillSimpHistograms(SimpHist_t* Hist, TSimParticle* Simp, SimpData_t* Sd, double Weight) {

  int stage  = Simp->GetUniqueID()/100000;
  
  Hist->fVolumeID->Fill(Simp->fEndVolumeIndex,Weight);
  Hist->fStage->Fill(stage,Weight);
  Hist->fGeneratorID->Fill(Simp->fGeneratorID,Weight);
  
  float xe = Simp->EndPos()->X()+3904.;
  float ye = Simp->EndPos()->Y();
  float ze = Simp->EndPos()->Z();
  float te = Simp->EndPos()->T();

  Hist->fTime->Fill(te,Weight);

  // Hist->fParentMom->Fill(fParent->StartMom()->P());
  // Hist->fParentPDG->Fill(fParent->PDGCode());

  float p = Simp->StartMom()->P();
  Hist->fStartMom[0]->Fill(p,Weight);
  Hist->fStartMom[1]->Fill(p,Weight);

  Hist->fYVsX->Fill(xe,ye,Weight);
  Hist->fXEndVsZEnd->Fill(ze,xe,Weight);
//-----------------------------------------------------------------------------
// looks like something to do with the stopping target - but this is still 34 foils..
//-----------------------------------------------------------------------------
  if (Simp->fEndVolumeIndex == 2480) Hist->fYVsX_2480->Fill(xe,ye,Weight);
  if (Simp->fEndVolumeIndex == 2513) Hist->fYVsX_2513->Fill(xe,ye,Weight);

  float cos_th = Simp->StartMom()->Pz()/p;
  Hist->fCosThVsMom->Fill(p,cos_th,Weight);

  Hist->fCosThVsMomPV->Fill(fPbarMomPV,fPbarCosThPV,Weight);
}

//-----------------------------------------------------------------------------
void TStnSpmcAnaModule::FillVDetHistograms(VDetHist_t* Hist, TStepPointMC* Step, double Weight) {

  int id = (Step->VolumeID());

  StnVDetData_t* vdd = fVDet+id;
  
  Hist->fIndex   ->Fill(id              ,Weight);
  Hist->fPDGCode ->Fill(Step->PDGCode() ,Weight);
  Hist->fGenCode ->Fill(Step->GenIndex(),Weight);

  float p = Step->Mom()->Mag();
  Hist->fMom[0]->Fill(p,Weight);
  Hist->fMom[1]->Fill(p,Weight);

  Hist->fTime    ->Fill(Step->Time(),Weight);

  // calculate local X and local Z

  double phi  = vdd->fPhiXZ*M_PI/180;

  double dx   = Step->Pos()->X()-vdd->fX;
  double dz   = Step->Pos()->Z()-vdd->fZ;

  double xloc =  dx*cos(phi) - dz*sin(phi);
  double zloc =  dx*sin(phi) + dz*cos(phi);

  Hist->fYVsX    ->Fill(xloc,Step->Pos()->Y(),Weight);
  Hist->fYVsZ    ->Fill(zloc,Step->Pos()->Y(),Weight);

  float py = Step->Mom()->Py();

  float px, pp, pt;
  if (vdd->fIZLocal == 3) {
    px = Step->Mom()->Px();
    pp = Step->Mom()->Pz();
  }
  else {
    px = Step->Mom()->Pz();
    pp = -Step->Mom()->Px();
  }

  pt = sqrt(px*px+py*py);

  Hist->fPt   ->Fill(pt   ,Weight);
  Hist->fPp   ->Fill(pp   ,Weight);
  Hist->fTanTh->Fill(pt/pp,Weight);

  float cos_th = pp/p;
  Hist->fCosThVsMom->Fill(p,cos_th,Weight);

  Hist->fCosThVsMomPV->Fill(fPbarMomPV,fPbarCosThPV,Weight);
}

//-----------------------------------------------------------------------------
// register data blocks and book histograms
//-----------------------------------------------------------------------------
int TStnSpmcAnaModule::BeginJob() {
//-----------------------------------------------------------------------------
// register data blocks 'SpmcBlock' or 'StepPointMCBlock' (old)
//-----------------------------------------------------------------------------
  RegisterDataBlock(fSpmcBlockName.Data(),"TStepPointMCBlock",&fSpmcBlock);
  RegisterDataBlock("SimpBlock"          ,"TSimpBlock"       ,&fSimpBlock       );
  RegisterDataBlock(fVDetBlockName.Data(),"TStepPointMCBlock",&fVDetBlock       );
//-----------------------------------------------------------------------------
// book histograms
//-----------------------------------------------------------------------------
  BookHistograms();

//-----------------------------------------------------------------------------
// initialize virtual detector offsets - a convenience for histogram filling
//-----------------------------------------------------------------------------
  StnInitVirtualDetectors(fVDet,&fNVDet);

  return 0;
}


//-----------------------------------------------------------------------------
void TStnSpmcAnaModule::FillHistograms() {

//-----------------------------------------------------------------------------
// event histograms
//
// EVT_0: all events
//-----------------------------------------------------------------------------
  FillEventHistograms(fHist.fEvent[0]);
//-----------------------------------------------------------------------------
// SIM_PARTICLE histograms
//-----------------------------------------------------------------------------
  SimpData_t* sd = nullptr;  // temp

  //  int pbar_stopped_in_st = 0;

  for (int i=0; i<fNSimp; i++) {
    TSimParticle* simp = fSimpBlock->Particle(i);
    int pdg_code = simp->PDGCode();
    int simp_id  = simp->GetUniqueID();
    //    int    vid1  = simp->EndVolumeIndex();
    //    double pend  = simp->EndMom()->P();
    //    double tend  = simp->EndPos()->T();

    FillSimpHistograms(fHist.fSimp[  0],simp,sd);
    if (pdg_code ==    13) FillSimpHistograms(fHist.fSimp[  3],simp,sd);
    if (pdg_code ==   -13) FillSimpHistograms(fHist.fSimp[  4],simp,sd);
    if (pdg_code ==  -211) FillSimpHistograms(fHist.fSimp[  5],simp,sd);
    if (pdg_code == -2212) FillSimpHistograms(fHist.fSimp[ 21],simp,sd);

    if      (simp_id < 100000) {
      FillSimpHistograms(fHist.fSimp[100],simp,sd);
      if (pdg_code ==    13) FillSimpHistograms(fHist.fSimp[103],simp,sd);
      if (pdg_code ==   -13) FillSimpHistograms(fHist.fSimp[104],simp,sd);
      if (pdg_code ==  -211) FillSimpHistograms(fHist.fSimp[105],simp,sd);
      if (pdg_code == -2212) FillSimpHistograms(fHist.fSimp[121],simp,sd);
    }
    else if (simp_id < 200000) {
      FillSimpHistograms(fHist.fSimp[200],simp,sd);
      if (pdg_code ==    13) FillSimpHistograms(fHist.fSimp[203],simp,sd);
      if (pdg_code ==   -13) FillSimpHistograms(fHist.fSimp[204],simp,sd);
      if (pdg_code ==  -211) FillSimpHistograms(fHist.fSimp[205],simp,sd);
      if (pdg_code == -2212) FillSimpHistograms(fHist.fSimp[221],simp,sd);
    }
    else if (simp_id < 300000) {
      FillSimpHistograms(fHist.fSimp[300],simp,sd);
      if (pdg_code ==    13) FillSimpHistograms(fHist.fSimp[303],simp,sd);
      if (pdg_code ==   -13) FillSimpHistograms(fHist.fSimp[304],simp,sd);
      if (pdg_code ==  -211) FillSimpHistograms(fHist.fSimp[305],simp,sd);
      if (pdg_code == -2212) FillSimpHistograms(fHist.fSimp[321],simp,sd);
    }
    else if (simp_id < 400000) {
      FillSimpHistograms(fHist.fSimp[400],simp,sd);
      if (pdg_code ==    13) FillSimpHistograms(fHist.fSimp[403],simp,sd);
      if (pdg_code ==   -13) FillSimpHistograms(fHist.fSimp[404],simp,sd);
      if (pdg_code ==  -211) FillSimpHistograms(fHist.fSimp[405],simp,sd);
      if (pdg_code == -2212) FillSimpHistograms(fHist.fSimp[421],simp,sd);
    }
    else if (simp_id < 500000) {
      FillSimpHistograms(fHist.fSimp[500],simp,sd);
      if (pdg_code ==    13) FillSimpHistograms(fHist.fSimp[503],simp,sd);
      if (pdg_code ==   -13) FillSimpHistograms(fHist.fSimp[504],simp,sd);
      if (pdg_code ==  -211) FillSimpHistograms(fHist.fSimp[505],simp,sd);
      if (pdg_code == -2212) FillSimpHistograms(fHist.fSimp[521],simp,sd);
    }
  }
//-----------------------------------------------------------------------------
// StepPointMC histograms
// for beamline studies, fStepPointMCBlock contains hits of particles for which 
// one of the stopping conditions has been satisfied, thus, this block contains 
// one StepPointMC per particle
//-----------------------------------------------------------------------------
  TStepPointMC* spmc;
  SpmcData_t    spmc_data;

  int nsteps = fSpmcBlock->NStepPoints();
  for (int i=0; i<nsteps; i++) {
    spmc             = fSpmcBlock->StepPointMC(i);
    //    float p          = spmc->Mom()->Mag();
    //    float t          = spmc->Time();
    int pdg_code     = spmc->PDGCode();
    int abs_pdg_code = abs(pdg_code);
//-----------------------------------------------------------------------------
// particles of interest are electrons, pions, muons and photons,
// speed up the mass extraction
//-----------------------------------------------------------------------------
    if (abs_pdg_code < 2500) spmc_data.fParticle = GetParticleCache(pdg_code);
    else                     spmc_data.fParticle = NULL;

    if (spmc_data.fParticle == NULL) {
      if (GetDebugBit(3) == 0) printf(">>> WARNING: no particle with PDF code=%i cached from ROOT particle DB\n",pdg_code);
    }

    FillStepPointMCHistograms(fHist.fStepPointMC[0],spmc,&spmc_data);

    if      (pdg_code ==   11) {
      FillStepPointMCHistograms(fHist.fStepPointMC[1],spmc,&spmc_data);
    }
    else if (pdg_code ==  -11) {
      FillStepPointMCHistograms(fHist.fStepPointMC[2],spmc,&spmc_data);
    }
    else if (pdg_code ==   13) {
      FillStepPointMCHistograms(fHist.fStepPointMC[3],spmc,&spmc_data);
    }
    else if (pdg_code ==   -13) FillStepPointMCHistograms(fHist.fStepPointMC[4],spmc,&spmc_data);
    else if (pdg_code ==    22) FillStepPointMCHistograms(fHist.fStepPointMC[5],spmc,&spmc_data);
    else if (pdg_code ==  -211) {
      FillStepPointMCHistograms(fHist.fStepPointMC[  6],spmc,&spmc_data);
    }
    else if (pdg_code     ==   211) FillStepPointMCHistograms(fHist.fStepPointMC[7],spmc,&spmc_data);
    else if (abs_pdg_code ==  2212) FillStepPointMCHistograms(fHist.fStepPointMC[8],spmc,&spmc_data);  // protons+pbars
    else                            FillStepPointMCHistograms(fHist.fStepPointMC[9],spmc,&spmc_data);  // everything else

    if (spmc->VolumeID() == 9) {
      FillStepPointMCHistograms(fHist.fStepPointMC[900],spmc,&spmc_data);
      if (pdg_code ==    13) FillStepPointMCHistograms(fHist.fStepPointMC[913],spmc,&spmc_data);
      if (pdg_code == -2212) FillStepPointMCHistograms(fHist.fStepPointMC[921],spmc,&spmc_data);
    }

    if (spmc->VolumeID() == 91) {
      FillStepPointMCHistograms(fHist.fStepPointMC[9100],spmc,&spmc_data);
      if (pdg_code ==    13) FillStepPointMCHistograms(fHist.fStepPointMC[9113],spmc,&spmc_data);
      if (pdg_code == -2212) FillStepPointMCHistograms(fHist.fStepPointMC[9121],spmc,&spmc_data);
    }
  }
  
//-----------------------------------------------------------------------------
// VDET histograms
//-----------------------------------------------------------------------------
  for (int i=0; i<fNVDetHits; i++) {
    TStepPointMC* step = fVDetBlock->StepPointMC(i);

    if (step->VolumeID() ==  9) FillVDetHistograms(fHist.fVDet[ 9],step);
    if (step->VolumeID() == 13) {
      FillVDetHistograms(fHist.fVDet[13],step);
      float x = step->Pos()->X()+3904.;
      float y = step->Pos()->Y();
      float r = sqrt(x*x+y*y);
      if ((r >= 400) && (r < 800)) { 
	FillVDetHistograms(fHist.fVDet[14],step);
      }
    }
  }
}



//_____________________________________________________________________________
int TStnSpmcAnaModule::BeginRun() {
  int rn = GetHeaderBlock()->RunNumber();
  TStntuple::Init(rn);
  return 0;
}


//_____________________________________________________________________________
int TStnSpmcAnaModule::Event(int ientry) {

  fSpmcBlock->GetEntry(ientry);
  fSimpBlock->GetEntry(ientry);
  fVDetBlock->GetEntry(ientry);
//-----------------------------------------------------------------------------
// assume electron in the first particle, otherwise the logic will need to 
// be changed
// if there are several hits, use the first one
//-----------------------------------------------------------------------------
  fNVDetHits = fVDetBlock->NStepPoints();
  fNSimp     = fSimpBlock->NParticles();

  fPbarCosThPV = -2.;
  fPbarMomPV   = -1.;
//-----------------------------------------------------------------------------
// number of step points with T>150
//-----------------------------------------------------------------------------
  fNSpmc     = fSpmcBlock->NStepPoints();
  fNSpmc150  = 0;

  for (int i=0; i<fNSpmc; i++) {
    TStepPointMC* spmc = fSpmcBlock->StepPointMC(i);
    if (spmc->Time() > 150) fNSpmc150++;
  }
//-----------------------------------------------------------------------------
// determine the cross section weight looking at the first particle with the PDG code of an antiproton
//-----------------------------------------------------------------------------
  fWeight = 1.;
  if (fNSimp > 0) {
//-----------------------------------------------------------------------------
// using the first antiproton in the list should work for old Bobs's dataset as well 
// as for the current ones
//-----------------------------------------------------------------------------
    for (int i=0; i<fNSimp; i++) {
      TSimParticle* sp0 = fSimpBlock->Particle(i);
      if (sp0->PDGCode() == -2212) {
//-----------------------------------------------------------------------------
// pbar production, assume Bob
// assuming parent particle exists, determine the production cross-section weight
// pbeam, nx, nz: the beam momentum and direction; beam ny = 0
//-----------------------------------------------------------------------------
	double pbeam(8.9), nx(-0.24192190), nz(-0.97029573);
	const TLorentzVector* sm = sp0->StartMom();
	double px    = sm->Px();
	double pz    = sm->Pz();
	double p     = sm->P();              // momentum in lab frame, MeV/c
	double costh = (px*nx+pz*nz)/p;
	double th    = TMath::ACos(costh);
//-----------------------------------------------------------------------------
// convert momentum to GeV/c
//-----------------------------------------------------------------------------
	double plab  = p/1000.;  

	fWeight      = fStnt->PBar_Striganov_d2N(pbeam,plab,th);
	fPbarMomPV   = p;
	fPbarCosThPV = costh;
	break;
      }
    }
  }
//-----------------------------------------------------------------------------
// determine simulation stage by looking at the last particle
//-----------------------------------------------------------------------------
  fStage = -1;
  if (fNSimp > 0) {
    TSimParticle* simp = fSimpBlock->Particle(fNSimp-1);
    int simp_id  = simp->GetUniqueID();
    fStage = simp_id / 100000;
  }

  for (int i=0; i<fNSimp; i++) {
    TSimParticle* simp = fSimpBlock->Particle(i);
    //    int pdg_code = simp->PDGCode();
    int simp_id  = simp->GetUniqueID();

    if (i < kMaxNSimp) {
      fSimData[i].fStage = simp_id/100000;
    }
    else {
      printf(" TStnSpmcAnaModule::Event ERROR: too many SimParticles\n");
    }
  }
//-----------------------------------------------------------------------------
// everything is precalculated, fill histograms
//-----------------------------------------------------------------------------
  FillHistograms();

  Debug();

  return 0;		       
}

//-----------------------------------------------------------------------------
// bit 4: pbars with X> 25cm at TS5 in
// bit 5: pbars ID=400000+I
// bit 6: pbar in the final state
// bit 7: pbar P > 100 MeV/c in the final state
//-----------------------------------------------------------------------------
void TStnSpmcAnaModule::Debug() {

  for (int i=0; i<fNSimp; i++) {
    TSimParticle* simp = fSimpBlock->Particle(i);
    int pdg_code = simp->PDGCode();
    int simp_id  = simp->GetUniqueID();

    if (GetDebugBit(4) == 1) {
      if ((pdg_code == -2212) && (simp_id > 300000) && (simp_id < 400000)) {
	float xe = simp->EndPos()->X()+3904.;
	float ye = simp->EndPos()->Y();
	float ze = simp->EndPos()->Z();
	float te = simp->EndPos()->T();
	if (xe > 250) GetHeaderBlock()->Print(Form("bit:4: xe, ye, ze, te = %10.3f %10.3f%10.3f %10.3f",xe,ye,ze,te));
      }
    }

    if (GetDebugBit(5) == 1) {
      if ((pdg_code == -2212) && (simp_id > 400000) && (simp_id < 500000)) {
	GetHeaderBlock()->Print(Form("bit:5: pbar ID=%10i",simp_id));
      }
    }
  }
//-----------------------------------------------------------------------------
// for all stages except S3, StepPointMC collection represents particles reaching 
// the "STOP" volume, one hit per particle
//-----------------------------------------------------------------------------
  if (GetDebugBit(6) == 1) {
    int nsteps = fSpmcBlock->NStepPoints();
    for (int i=0; i<nsteps; i++) {
      TStepPointMC* spmc = fSpmcBlock->StepPointMC(i);
      float p            = spmc->Mom()->Mag();
      float t            = spmc->Time();
      int pdg_code       = spmc->PDGCode();
      
      if (pdg_code          == -2212) {
	GetHeaderBlock()->Print(Form("bit:6: pbar in the end, p = %10.3f t= %10.3f",p,t));
      }
    }
  }    
  if (GetDebugBit(7) == 1) {
    int nsteps = fSpmcBlock->NStepPoints();
    for (int i=0; i<nsteps; i++) {
      TStepPointMC* spmc = fSpmcBlock->StepPointMC(i);
      float p            = spmc->Mom()->Mag();
      float t            = spmc->Time();
      int pdg_code       = spmc->PDGCode();
      
      if ((pdg_code          == -2212) && (p > 100)) {
	GetHeaderBlock()->Print(Form("bit:7: pbar P > 100 in the end, p = %10.3f t= %10.3f",p,t));
      }
    }
  }
//-----------------------------------------------------------------------------
// bit:8  events with pi- at VD9
//-----------------------------------------------------------------------------
  if (GetDebugBit(8) == 1) {
    for (int i=0; i<fNVDetHits; i++) {
      TStepPointMC* step = fVDetBlock->StepPointMC(i);
      if      ((step->PDGCode() == -211) && (step->VolumeID() ==  9)) {
	float p            = step->Mom()->Mag();
	float t            = step->Time();
	GetHeaderBlock()->Print(Form("bit:8: pi- at VD 09, p = %10.3f t= %10.3f",p,t));
      }
    }
  }
//-----------------------------------------------------------------------------
// bit:9  events with pi- at stage4
//-----------------------------------------------------------------------------
  if (GetDebugBit(9) == 1) {
    for (int i=0; i<fNSimp; i++) {
      TSimParticle* simp = fSimpBlock->Particle(i);
      int pdg_code = simp->PDGCode();
      int simp_id  = simp->GetUniqueID();

      if (simp_id > 400000) {
	if (pdg_code ==  -211) {
	  GetHeaderBlock()->Print(Form("bit:9: pi- at stage 4"));
	}
      }
    }
  }
}

//_____________________________________________________________________________
int TStnSpmcAnaModule::EndJob() {
  printf("----- end job: ---- %s\n",GetName());
  return 0;
}

//_____________________________________________________________________________
void TStnSpmcAnaModule::Test001() {
}

