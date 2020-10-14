///////////////////////////////////////////////////////////////////////////////
// use of tmp:
//
// Tmp(0) : 
// Tmp(1) : 
// 
// use of debug bits: bits 0-2 are reserved
// 0  : all events
// 1  : passed events
// 2  : rejected events
// 
// 3  : 
// 4  : 
// 5  : 
///////////////////////////////////////////////////////////////////////////////
#include "TF1.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TEnv.h"
#include "TSystem.h"

#include "Stntuple/loop/TStnAna.hh"
#include "Stntuple/obj/TStnHeaderBlock.hh"
#include "Stntuple/alg/TStntuple.hh"
#include "Stntuple/geom/TDisk.hh"
#include "Stntuple/geom/TStnCrystal.hh"
#include "Stntuple/val/stntuple_val_functions.hh"
//------------------------------------------------------------------------------
// Mu2e offline includes
//-----------------------------------------------------------------------------
#include "ana/TStnClusterAnaModule.hh"

ClassImp(TStnClusterAnaModule)
//-----------------------------------------------------------------------------
TStnClusterAnaModule::TStnClusterAnaModule(const char* name, const char* title):
  TStnModule(name,title)
{
  fDiskCalorimeter = new TDiskCalorimeter();
  fMinT0           = 0;                 // do not cut on time by default
}

//-----------------------------------------------------------------------------
TStnClusterAnaModule::~TStnClusterAnaModule() {
}


//-----------------------------------------------------------------------------
void TStnClusterAnaModule::BookClusterHistograms(ClusterHist_t* Hist, const char* Folder) {
//   char name [200];
//   char title[200];

  HBook1F(Hist->fDiskID ,"disk_id",Form("%s: Disk ID"       ,Folder), 10, 0,  10,Folder);
  HBook1F(Hist->fEnergy ,"energy" ,Form("%s: Cluster Energy",Folder),500, 0, 250,Folder);
  HBook1F(Hist->fT0     ,"t0"     ,Form("%s: cluster T0"    ,Folder),200, 0,2000,Folder);
  HBook1F(Hist->fRow    ,"row"    ,Form("%s: cluster Row"   ,Folder),200, 0, 200,Folder);
  HBook1F(Hist->fCol    ,"col"    ,Form("%s: cluster column",Folder),200, 0, 200,Folder);
  HBook1F(Hist->fX      ,"x"      ,Form("%s: cluster X"     ,Folder),200, -5000,5000,Folder);
  HBook1F(Hist->fY      ,"y"      ,Form("%s: cluster Y"     ,Folder),200,-1000,1000,Folder);
  HBook1F(Hist->fZ      ,"z"      ,Form("%s: cluster Z"     ,Folder),200, 11500,13500,Folder);
  HBook1F(Hist->fR      ,"r"      ,Form("%s: cluster Radius",Folder),100, 0,  1000,Folder);
  HBook1F(Hist->fYMean  ,"ymean"  ,Form("%s: cluster YMean" ,Folder),400,-200,200,Folder);
  HBook1F(Hist->fZMean  ,"zmean"  ,Form("%s: cluster ZMean" ,Folder),400,-200,200,Folder);
  HBook1F(Hist->fSigY   ,"sigy"   ,Form("%s: cluster SigY"  ,Folder),100, 0,100,Folder);
  HBook1F(Hist->fSigZ   ,"sigz"   ,Form("%s: cluster SigZ"  ,Folder),100, 0,100,Folder);
  HBook1F(Hist->fSigR   ,"sigr"   ,Form("%s: cluster SigR"  ,Folder),100, 0,100,Folder);
  HBook1F(Hist->fNCr0   ,"ncr0"   ,Form("%s: cluster NCR[0]",Folder),100, 0,100,Folder);
  HBook1F(Hist->fNCr1   ,"ncr1"   ,Form("%s: cluster NCR[1]",Folder),100, 0,100,Folder);
  HBook1F(Hist->fFrE1   ,"fre1"   ,Form("%s: E1/Etot"       ,Folder),200, 0,  1,Folder);
  HBook1F(Hist->fFrE2   ,"fre2"   ,Form("%s: (E1+E2)/Etot"  ,Folder),200, 0,  1,Folder);
  HBook1F(Hist->fSigE1  ,"sige1"   ,Form("%s: SigmaE/Etot"  ,Folder),200, 0, 10,Folder);
  HBook1F(Hist->fSigE2  ,"sige2"   ,Form("%s: SigmaE/Emean" ,Folder),200, 0, 10,Folder);
}

//-----------------------------------------------------------------------------
void TStnClusterAnaModule::BookEventHistograms(EventHist_t* Hist, const char* Folder) {
  //  char name [200];
  //  char title[200];

  HBook1F(Hist->fRv         ,"rv"      ,Form("%s: R(Vertex)"                       ,Folder), 100, 0, 200,Folder);
  HBook1F(Hist->fZv         ,"zv"      ,Form("%s: Z(Vertex)"                       ,Folder), 400, 5000,7000,Folder);
  HBook1F(Hist->fEleCosTh  ,"ce_costh" ,Form("%s: Conversion Electron Cos(Theta)"  ,Folder), 100,-1,    1,Folder);
  HBook1F(Hist->fNClusters ,"ncl"      ,Form("%s: Number of Reconstructed Clusters",Folder),1000, 0, 1000,Folder);
  HBook1F(Hist->fNCl20     ,"ncl20"    ,Form("%s: N(Clusters), E > 20"             ,Folder), 100, 0,  100,Folder);
  HBook1F(Hist->fNCl50     ,"ncl50"    ,Form("%s: N(Clusters), E > 50"             ,Folder), 100, 0,  100,Folder);
  HBook1F(Hist->fNCl70     ,"ncl70"    ,Form("%s: N(Clusters), E > 70"             ,Folder), 100, 0,  100,Folder);
  HBook1F(Hist->fEMax      ,"emax"     ,Form("%s: Max cluster energy"              ,Folder), 250, 0,  500,Folder);
  HBook1F(Hist->fNGenp     ,"ngenp"    ,Form("%s: N(Gen Particles)"                ,Folder), 500, 0,10000,Folder);

  for (int i=0; i<kNDisks; i++) {
    HBook1F(Hist->fNHits       [i],Form("nhits_%i"   ,i),Form("%s: N(Hits) [%i]"           ,Folder,i),  500,0,10000,Folder);
    HBook1F(Hist->fETot        [i],Form("etot_%i"    ,i),Form("%s: Etot[%i]"               ,Folder,i), 1000,0,10000,Folder);
    HBook1D(Hist->fECrVsR      [i],Form("ecr_vs_r_%i",i),Form("%s: E Cr Vs R [%i]"         ,Folder,i),  100,0, 1000,Folder);
    HBook1D(Hist->fNHitsVsR    [i],Form("nh_vs_r_%i" ,i),Form("%s: N hits Vs R [%i]"       ,Folder,i),  100,0, 1000,Folder);
    HBook1D(Hist->fNCrystalsVsR[i],Form("ncr_vs_r_%i",i),Form("%s: N Hit Crystals[%i] vs R",Folder,i),  100,0, 1000,Folder);

    HBook2F(Hist->fEHitVsR     [i],Form("eh_vs_r_%i",i),Form("%s: E Hit Vs R [%i]",Folder,i),100,0,1000,250,0,500,Folder);
  }

  HBook1F(Hist->fNHitsTot       ,"nhits_tot",Form("%s: Total number of hits",Folder), 400, 0,20000,Folder);
  HBook1F(Hist->fNHitCrystalsTot,"nhcr_tot" ,Form("%s: NHit Crystals Tot"   ,Folder), 800, 0,4000,Folder);

  HBook1F(Hist->fECaloTot       ,"ecal"     ,Form("%s: Total energy in the calorimeter",Folder), 500,0,250,Folder);
}

//_____________________________________________________________________________
void TStnClusterAnaModule::BookHistograms() {

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
  book_event_histset[ 1] = 1;		// EMin = 0.1
  book_event_histset[ 2] = 1;		// EMin = 0.5
  book_event_histset[ 3] = 1;		// EMin = 1.0
  book_event_histset[ 4] = 1;		// EMin = 1.5
  book_event_histset[ 5] = 1;		// EMin = 2.0
  book_event_histset[ 6] = 1;		// EMin = 1.0 , TMin = 400
  book_event_histset[ 7] = 1;		// EMin = 1.0 , TMin = 500
  book_event_histset[ 8] = 1;		// EMin = 1.0 , TMin = 600

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
// book cluster histograms
//-----------------------------------------------------------------------------
  int book_cluster_histset[kNClusterHistSets];
  for (int i=0; i<kNClusterHistSets; i++) book_cluster_histset[i] = 0;

  book_cluster_histset[0] = 1;		// all clusters
  book_cluster_histset[1] = 1;		// clusters disk=0
  book_cluster_histset[2] = 1;		// clusters disk=1
  book_cluster_histset[3] = 1;		// E > 20 MeV
  book_cluster_histset[4] = 1;		// E > 50 MeV
  book_cluster_histset[5] = 1;		// E > 70 MeV

  for (int i=0; i<kNClusterHistSets; i++) {
    if (book_cluster_histset[i] != 0) {
      sprintf(folder_name,"cls_%i",i);
      fol = (TFolder*) hist_folder->FindObject(folder_name);
      if (! fol) fol = hist_folder->AddFolder(folder_name,folder_name);
      fHist.fCluster[i] = new ClusterHist_t;
      BookClusterHistograms(fHist.fCluster[i],Form("Hist/%s",folder_name));
    }
  }
}

//-----------------------------------------------------------------------------
void TStnClusterAnaModule::FillClusterHistograms(ClusterHist_t* Hist, TStnCluster* Cluster) {
  int   row, col;
  float  x, y, z, r;

  row = Cluster->Ix1();
  col = Cluster->Ix2();

  x   = Cluster->fX; // +3904.;
  y   = Cluster->fY;
  z   = Cluster->fZ;
  r   = sqrt(x*x+y*y);

  if ((row < 0) || (row > 9999)) row = -9999;
  if ((col < 0) || (col > 9999)) col = -9999;

  Hist->fDiskID->Fill(Cluster->DiskID());
  Hist->fEnergy->Fill(Cluster->Energy());
  Hist->fT0->Fill(Cluster->Time());
  Hist->fRow->Fill(row);
  Hist->fCol->Fill(col);
  Hist->fX->Fill(x);
  Hist->fY->Fill(y);
  Hist->fZ->Fill(z);
  Hist->fR->Fill(r);

  Hist->fYMean->Fill(Cluster->fYMean);
  Hist->fZMean->Fill(Cluster->fZMean);
  Hist->fSigY ->Fill(Cluster->fSigY);
  Hist->fSigZ ->Fill(Cluster->fSigZ);
  Hist->fSigR ->Fill(Cluster->fSigR);
  Hist->fNCr0 ->Fill(Cluster->fNCrystals);
  Hist->fNCr1 ->Fill(Cluster->fNCr1);
  Hist->fFrE1 ->Fill(Cluster->fFrE1);
  Hist->fFrE2 ->Fill(Cluster->fFrE2);
  Hist->fSigE1->Fill(Cluster->fSigE1);
  Hist->fSigE2->Fill(Cluster->fSigE2);
}

//-----------------------------------------------------------------------------
// 
//-----------------------------------------------------------------------------
void TStnClusterAnaModule::FillEventHistograms(EventHist_t* Hist, double EMin, double TMin) {
  double            cos_th, xv, yv, rv, zv, p;
  TLorentzVector    mom;

  if (fElectron) {
    fElectron->Momentum(mom);

    p      = mom.P();
    cos_th = mom.Pz()/p;

    xv     = fElectron->Vx()+3904.;
    yv     = fElectron->Vy();
    rv     = sqrt(xv*xv+yv*yv);
    zv     = fElectron->Vz();

    Hist->fEleCosTh->Fill(cos_th);
    Hist->fRv->Fill(rv);
    Hist->fZv->Fill(zv);
  }

  Hist->fNClusters->Fill(fNClusters);
//   Hist->fNCl20->Fill(fNCl20);
//   Hist->fNCl50->Fill(fNCl50);
//   Hist->fNCl70->Fill(fNCl70);

  double emax = -1;

  TStnCluster* cluster(0);
  if (fNClusters > 0) cluster = fClusterBlock->Cluster(0);

  if (cluster) {
    emax   = cluster->Energy();
  }

  Hist->fEMax->Fill(emax);
  //  Hist->fNGenp->Fill(fNGenp);
//-----------------------------------------------------------------------------
// calorimeter
//-----------------------------------------------------------------------------
  int      ndisks, nc, nh, nhh, bin;

  int      n_hits_tot, n_hit_crystals_tot;
  int      n_hit_crystals[kNDisks], n_hits[kNDisks];
  int      n_hits_r[kNDisks][100], n_hit_crystals_r[kNDisks][100];

  float    etot[kNDisks];

  TCalHitData* hit;
  TStnCrystal* cr;

  float        ehit, r;

  ndisks = fDiskCalorimeter->NDisks();

  n_hits_tot         = 0;
  n_hit_crystals_tot = 0;

  for (int idisk=0; idisk<ndisks; idisk++) {

    TDisk* disk = fDiskCalorimeter->Disk(idisk);

    etot          [idisk] = 0;
    n_hits        [idisk] = 0;
    n_hit_crystals[idisk] = 0;

    for (int ib=0; ib<100; ib++) {
      n_hits_r        [idisk][ib] = 0;
      n_hit_crystals_r[idisk][ib] = 0;
    }

    nc = disk->NCrystals();
    for (int ic=0; ic<nc; ic++) {
      cr   = disk->Crystal(ic);
      r    = cr->Radius();
      nh   = cr->NHits ();
      ehit = cr->Energy();
      bin  = (int) (r/10.);
					// total energy deposited in the disk
      etot[idisk] += cr->Energy();

      nhh = 0;
      for (int ih=0; ih<nh; ih++) {
	hit = cr->CalHitData(ih);
	if ((hit->Energy() > EMin) && (hit->Time() > TMin)) {
					// number of hits above the threshold (EMin)

	  nhh                         += 1; // N(hits) per crystal
	  n_hits         [idisk]      += 1; // total N(hits) in the disk
	  n_hits_r       [idisk][bin] += 1; // total N(hits) in a given radial bin
	}
      }
    
      Hist->fECrVsR   [idisk]->Fill(r,ehit);
      Hist->fEHitVsR  [idisk]->Fill(r,ehit);

      if (nhh > 0) {
	n_hit_crystals  [idisk]      += 1;
	n_hit_crystals_r[idisk][bin] += 1;
      }
    }
    n_hit_crystals_tot += n_hit_crystals[idisk];
    n_hits_tot         += n_hits[idisk];
  }

  double ecalo(0);
  for (int idisk=0; idisk<ndisks; idisk++) {
    ecalo += etot[idisk];
//-----------------------------------------------------------------------------
// fill 'per-disk' histograms
//-----------------------------------------------------------------------------
    Hist->fETot [idisk]->Fill(etot  [idisk]);
    Hist->fNHits[idisk]->Fill(n_hits[idisk]);
//-----------------------------------------------------------------------------
// 100 is fixed by the number of bins in the radial distributions
//-----------------------------------------------------------------------------
    for (int ib=0; ib<100; ib++) {
      r = (ib+0.5)*10.;
      Hist->fNHitsVsR    [idisk]->Fill(r,n_hits_r        [idisk][ib]);
      Hist->fNCrystalsVsR[idisk]->Fill(r,n_hit_crystals_r[idisk][ib]);
    }
  }

  Hist->fECaloTot->Fill(ecalo);
  Hist->fNHitCrystalsTot->Fill(n_hit_crystals_tot);
  Hist->fNHitsTot->Fill(n_hits_tot);
}

//-----------------------------------------------------------------------------
// register data blocks and book histograms
//-----------------------------------------------------------------------------
int TStnClusterAnaModule::BeginJob() {
//-----------------------------------------------------------------------------
// register data blocks
//-----------------------------------------------------------------------------
//  RegisterDataBlock("TrackBlock"    ,"TStnTrackBlock"   ,&fTrackBlock  );
  RegisterDataBlock("ClusterBlock"  ,"TStnClusterBlock" ,&fClusterBlock);
  RegisterDataBlock("CalDataBlock"  ,"TCalDataBlock"    ,&fCalDataBlock);
  //  RegisterDataBlock("GenpBlock"     ,"TGenpBlock"       ,&fGenpBlock);
//   RegisterDataBlock("SimpBlock"     ,"TSimpBlock"       ,&fSimpBlock);
//-----------------------------------------------------------------------------
// book histograms
//-----------------------------------------------------------------------------
  BookHistograms();

  return 0;
}


//_____________________________________________________________________________
void TStnClusterAnaModule::FillHistograms() {

  //  static int first_entry(1);
//-----------------------------------------------------------------------------
// event histograms
//
// EVT_0: all events
//-----------------------------------------------------------------------------
  FillEventHistograms(fHist.fEvent[0], 0.0,300);
  FillEventHistograms(fHist.fEvent[1], 0.1,300);
  FillEventHistograms(fHist.fEvent[2], 0.5,300);
  FillEventHistograms(fHist.fEvent[3], 1.0,300);
  FillEventHistograms(fHist.fEvent[4], 1.5,300);
  FillEventHistograms(fHist.fEvent[5], 2.0,300);
//-----------------------------------------------------------------------------
// EVT_6: E> 1 MeV, T>400
// EVT_7: E> 1 MeV, T>500
// EVT_8: E> 1 MeV, T>600
//-----------------------------------------------------------------------------
  FillEventHistograms(fHist.fEvent[6], 1.0,400);
  FillEventHistograms(fHist.fEvent[7], 1.0,500);
  FillEventHistograms(fHist.fEvent[8], 1.0,600);
//-----------------------------------------------------------------------------
// cluster histograms
//-----------------------------------------------------------------------------
  TStnCluster* cl;

  for (int i=0; i<fNClusters; ++i ) {
    cl = fClusterBlock->Cluster(i);

    FillClusterHistograms(fHist.fCluster[0],cl);

    if (cl->DiskID() == 0)  FillClusterHistograms(fHist.fCluster[1],cl);
    if (cl->DiskID() == 1)  FillClusterHistograms(fHist.fCluster[2],cl);

    if (cl->Energy() > 20.) FillClusterHistograms(fHist.fCluster[3],cl);
    if (cl->Energy() > 50.) FillClusterHistograms(fHist.fCluster[4],cl);
    if (cl->Energy() > 70.) FillClusterHistograms(fHist.fCluster[5],cl);
  }
}



//_____________________________________________________________________________
int TStnClusterAnaModule::BeginRun() {
  int rn = GetHeaderBlock()->RunNumber();
  TStntuple::Init(rn);
  return 0;
}


//_____________________________________________________________________________
int TStnClusterAnaModule::Event(int ientry) {

  //  double                p;
  //  TEmuLogLH::CalData_t  dat;
  //  TStnTrack*            track;
  //  int                   id_word;
  TLorentzVector        mom;

  TDiskCalorimeter::GeomData_t disk_geom;

  fClusterBlock->GetEntry(ientry);
  fCalDataBlock->GetEntry(ientry);
  //  fGenpBlock->GetEntry(ientry);
  //  fSimpBlock->GetEntry(ientry);

  fElectron = NULL;
//-----------------------------------------------------------------------------
// assume electron in the first particle, otherwise the logic will need to 
// be changed
//-----------------------------------------------------------------------------
  if (fDiskCalorimeter->Initialized() == 0) {
    disk_geom.fNDisks = fCalDataBlock->NDisks();

    for (int i=0; i<disk_geom.fNDisks; i++) {
      //      disk_geom.fNCrystals[i] = fCalDataBlock->fNCrystals[i];
      disk_geom.fRMin[i]      = fCalDataBlock->fRMin[i];
      disk_geom.fRMax[i]      = fCalDataBlock->fRMax[i];
      disk_geom.fZ0  [i]      = fCalDataBlock->fZ0  [i];
    }

    disk_geom.fHexSize          = fCalDataBlock->CrystalSize()*2;
    // kludge , so far
    disk_geom.fMinFraction      = 1.; // fCalDataBlock->MinFraction();
    disk_geom.fWrapperThickness = fCalDataBlock->WrapperThickness();
    disk_geom.fShellThickness   = fCalDataBlock->ShellThickness();

    fDiskCalorimeter->Init(&disk_geom);
  }

  //  TStnCluster* cl;
  fNClusters  = fClusterBlock->NClusters();

  if (fNClusters == 0) fCluster = 0;
  else                 fCluster = fClusterBlock->Cluster(0);

//   fNCl20      = 0;
//   fNCl50      = 0;
//   fNCl70      = 0;
  for (int i=0; i<fNClusters; ++i ) {
    //    cl = fClusterBlock->Cluster(i);
//     if (cl->Energy() > 20.) fNCl20 += 1;
//     if (cl->Energy() > 50.) fNCl50 += 1;
//     if (cl->Energy() > 70.) fNCl70 += 1;
  }

  //  fNCalHits   = fCalDataBlock->NHits();

  fDiskCalorimeter->InitEvent(fCalDataBlock);

  FillHistograms();

  Debug();

  return 0;		       
}

//-----------------------------------------------------------------------------
void TStnClusterAnaModule::Debug() {

// //-----------------------------------------------------------------------------
// // bit 5: events with N(tracks) > 1
// //-----------------------------------------------------------------------------
//   if (GetDebugBit(5) == 1) {
//     int ntrk = fTrackBlock->NTracks();
//     if (ntrk > 1) {
//       GetHeaderBlock()->Print(Form("NTracks = %i5",ntrk));
//     }
//   }
}

//_____________________________________________________________________________
int TStnClusterAnaModule::EndJob() {
  printf("----- end job: ---- %s\n",GetName());
  return 0;
}

//_____________________________________________________________________________
void TStnClusterAnaModule::Test001() {

  // mu2e::HexMap* hmap      = new mu2e::HexMap();

  // mu2e::HexLK hex_index(0,0);

  // for (int i=0; i<40; i++) {
  //   hex_index = hmap->lk(i);
  //   printf(" i,l,k = %5i %5i %5i\n",i,hex_index._l,hex_index._k);
  // }
}

