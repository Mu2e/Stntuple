//-----------------------------------------------------------------------------
//  Dec 28 2000 P.Murat: base class for STNTUPLE analysis module
//-----------------------------------------------------------------------------
#include <iostream>
#include <format>

using namespace std;

#include "stdlib.h"
#include "time.h"

#include "TPad.h"
#include "TCanvas.h"
#include "TText.h"
#include "TROOT.h"
#include "TSystem.h"

#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TFitResult.h"

#include "Stntuple/base/TStnBookHist.hh"

#include "Stntuple/obj/TStnEvent.hh"
#include "Stntuple/obj/TStnGoodRunList.hh"

#include "Stntuple/loop/TStnModule.hh"
#include "Stntuple/loop/TStnAna.hh"
#include "TFolder.h"

ClassImp(TStnModule)
//_____________________________________________________________________________
TStnModule::TStnModule() 
{
  fEnabled           =  1;
  fInitialized       =  0;
  fLastRun           = -1;
  fMyronFlag         = -1;
  fPrintLevel        =  0;
  fFilteringMode     =  0;
				// by default all the events pass
  fPassed            =  1;
  fFolder            =  0;
  fListOfL3TrigNames =  0;
  for (int i=0; i<kNDebugBits; i++) fDebugBit[i] = 0;

  fBookHist          = nullptr;
}

//_____________________________________________________________________________
TStnModule::TStnModule(const char* name, const char* title):TNamed(name,title)
{
  // kludge around until ROOT v1.05+ (version after Feb 09 2003) comes to FNAL
  // struct Folder_t : public TNamed {
  //   TCollection* fFolders;
  //   Bool_t       fIsOwner;
  // };

  // Folder_t* x;

  fEnabled       = 1;
  fInitialized   = 0;
  fLastRun       = -1;
  fMyronFlag     = -1;
  fPrintLevel    =  0;
  fFilteringMode =  0;
				// by default all the events pass
  fPassed        =  1;

  fFolder = new TFolder(name,title);
				// here comes the hack (see above)
  // x = (Folder_t*) fFolder;
  // x->SetName(name);
  // x->SetTitle(title);
  // x->fFolders = new TList();

  fFolder->Add(this);

  // TFolder* hist = fFolder->AddFolder("Hist","ListOfHistograms");
  fFolder->AddFolder("Hist","ListOfHistograms");

  fBookHist          = new TStnBookHist(fFolder);

  fListOfL3TrigNames = new TObjArray();
  fListOfL3Triggers  = new TObjArray();

  for (int i=0; i<kNDebugBits; i++) fDebugBit[i] = 0;
}

//_____________________________________________________________________________
TStnModule::~TStnModule() {
  // destructor: module owns its histograms, but there could be other objects
  // added to it... derived classes should not delete objecs added to fFolder

  delete fBookHist;
  delete fFolder;
  
  fListOfL3TrigNames->Delete();
  delete fListOfL3TrigNames;
  fListOfL3Triggers->Delete();
  fListOfL3Triggers->Clear();
  delete fListOfL3Triggers;
}

//_____________________________________________________________________________
TStnHeaderBlock* TStnModule::GetHeaderBlock() {
  // return pointer to the header block

  return fAna->GetHeaderBlock();
}

//_____________________________________________________________________________
TStnEvent* TStnModule::GetEvent() {
  // return processed event

  return fAna->GetEvent();
}

//_____________________________________________________________________________
TStnGoodRunList* TStnModule::GetGoodRunList() {
  // return processed event

  return fAna->GetGoodRunList();
}

//_____________________________________________________________________________
void* TStnModule::RegisterDataBlock(const char* BranchName, 
				    const char* ClassName) 
{
  TStnModule* x = (TStnModule*) fAna->RegisterDataBlock(BranchName,ClassName);
  if (! x) {
    Error("RegisterDataBlock",
	  Form("No %-20s branch requested by %-20s in the input STNTUPLE",
	       BranchName, GetName()));
  }
  return x;
}

//_____________________________________________________________________________
void TStnModule::RegisterDataBlock(const char* BranchName, void* DataBlock) {
  fAna->RegisterDataBlock(BranchName,DataBlock);
  if (! *((void**)DataBlock)) {
    Error("RegisterDataBlock",
	  Form("No %-20s branch requested by %-20s in the input STNTUPLE",
	       BranchName, GetName()));
  }
}

//_____________________________________________________________________________
void TStnModule::RegisterDataBlock(const char* BranchName, 
				   const char* ClassName, 
				   void*       DataBlock) 
{
  fAna->RegisterDataBlock(BranchName,ClassName,DataBlock);
  if (! *((void**)DataBlock)) {
    Error("RegisterDataBlock",
	  Form("No %-20s branch requested by %-20s in the input STNTUPLE",
	       BranchName, GetName()));
  }
}

//_____________________________________________________________________________
Int_t TStnModule::AddOutputBlock(const char* BranchName, TStnDataBlock* Block)
{
  // add new branch to the output tree and associate it with the data block
  // 'Block'

  fAna->GetEvent()->AddOutputBlock(BranchName,Block);
  return 0;
}

//_____________________________________________________________________________
int TStnModule::BeginJob() {

  const char* env;
  env = gSystem->Getenv(Form("%s_PrintLevel",GetName()));
  if (env) {
    fPrintLevel = atoi(env);
    printf("%s: fPrintLevel = %i\n",GetName(),fPrintLevel);
  }
//-----------------------------------------------------------------------------
//  one more handle on filtering mode
//-----------------------------------------------------------------------------
  env = gSystem->Getenv(Form("%s_FilteringMode",GetName()));
  if (env) {
    fFilteringMode = atoi(env);
    printf("%s: fFilteringMode = %i\n",GetName(),fFilteringMode);
  }

  fFolder->SetName (GetName ());
  fFolder->SetTitle(GetTitle());

  return 0;
}

//_____________________________________________________________________________
int TStnModule::BeginRun() {
  return 0;
}

//_____________________________________________________________________________
int TStnModule::Event(Int_t i) {
  return 0;
}

//_____________________________________________________________________________
int TStnModule::EndRun() {
  return 0;
}

//_____________________________________________________________________________
int TStnModule::EndJob() {
  return 0;
}

//_____________________________________________________________________________
TCanvas* TStnModule::NewSlide(const char* name, 
			      const char* title, 
			      int nx, int ny) 
{
  // create new canvas with user-defined number of pads

  char canvas_name[120], p1_name[130], p2_name[130], p3_name[130];

  TCanvas* slide;

  if (name ) sprintf(canvas_name,"%s_%s",GetName(),name);
  else       sprintf(canvas_name,"%s_slide",GetName());

  sprintf(p1_name,"%s_p1",canvas_name);
  sprintf(p2_name,"%s_p2",canvas_name);
  sprintf(p3_name,"%s_p3",canvas_name);

  slide = new TCanvas(canvas_name,canvas_name,0,0,600,800);

  TPad *p1 = new TPad(p1_name,p1_name,0.01,0.01,0.99,0.96);

  p1->Divide(nx,ny);
  p1->Draw();
  p1->Range(0,0,1,1);
  slide->cd();

  TPad *p2 = new TPad(p2_name, p2_name,0.7,0.965,0.99,0.995);
  p2->Draw();
  p2->cd();

  time_t t2 = time(0);
  tm* t22 = localtime(&t2);
  TText *text = new TText(0.05,0.3,asctime(t22));
  text->SetTextSize(0.5);
  text->Draw();
  p2->Modified();

  slide->cd();
				// create title pad

  TPad *title_pad = new TPad(p3_name, "title",0.05,0.965,0.69,0.995);
  title_pad->Draw();
  title_pad->cd();

  if (title) {
    TText *text = new TText(0.05,0.3,title);
    text->SetTextSize(0.5);
    text->Draw();
  }
  title_pad->Modified();

  return slide;
}

//_____________________________________________________________________________
void     TStnModule::AddHistogram(TObject* hist, const char* FolderName) {
  if (hist == nullptr) return;
  
  TFolder* fol = (TFolder*) fFolder->FindObject(FolderName);
  if (fol) {
    fBookHist->AddHistogram(hist,fol);
  }
  else {
    printf("TStnModule::%s: Unable to find folder %s\n", __func__, FolderName);
    return;
  }
}

//-----------------------------------------------------------------------------
void TStnModule::HBook1F(TH1F*& Hist, const char* Name, const char* Title,
			 Int_t Nx, Double_t XMin, Double_t XMax,
			 const char* FolderName)
{
  // book 2D histogram, add it to the module's list of histograms and 
  // return pointer to it to the user

  Hist = new TH1F(Name,Title,Nx,XMin,XMax);
  AddHistogram(Hist,FolderName);
}

//_____________________________________________________________________________
void TStnModule::HBook1F(TH1F*& Hist, const char* Name, const char* Title,
			 Int_t Nx, const float* XMin,
			 const char* FolderName)
{
  // book 1D histogram with variable size bins, add it to the module's list of histograms and 
  // return pointer to it to the user

  Hist = new TH1F(Name,Title,Nx,XMin);
  AddHistogram(Hist,FolderName);
}

//_____________________________________________________________________________
void TStnModule::HBook1D(TH1D*& Hist, const char* Name, const char* Title,
			 Int_t Nx, Double_t XMin, Double_t XMax,
			 const char* FolderName)
{
  // this is introduced specifically for weighting the DIO Mu2e spectrum
  // book 1D histogram, add it to the module's list of histograms and 
  // return pointer to it to the user

  Hist = new TH1D(Name,Title,Nx,XMin,XMax);
  AddHistogram(Hist,FolderName);
}

//_____________________________________________________________________________
void TStnModule::HBook1D(TH1D*& Hist, const char* Name, const char* Title,
			 Int_t Nx, const double* XMin,
			 const char* FolderName)
{
  // book 1D histogram with variable size bins, add it to the module's list of histograms and 
  // return pointer to it to the user

  Hist = new TH1D(Name,Title,Nx,XMin);
  AddHistogram(Hist,FolderName);
}

//_____________________________________________________________________________
void TStnModule::HBook2F(TH2F*& Hist, const char* Name, const char* Title,
			 Int_t Nx, Double_t XMin, Double_t XMax,
			 Int_t Ny, Double_t YMin, Double_t YMax,
			 const char* FolderName)
{
  // book 2D histogram, add it to the module's list of histograms and 
  // return pointer to it to the user

  Hist = new TH2F(Name,Title,Nx,XMin,XMax,Ny,YMin,YMax);
  AddHistogram(Hist,FolderName);
}

//_____________________________________________________________________________
void TStnModule::HProf(TProfile*& Hist, const char* Name, const char* Title,
		       Int_t Nx, Double_t XMin, Double_t XMax,
		       Double_t YMin, Double_t YMax,
		       const char* FolderName)
{
  // book 2D histogram, add it to the module's list of histograms and 
  // return pointer to it to the user

  Hist = new TProfile(Name,Title,Nx,XMin,XMax,YMin,YMax);
  AddHistogram(Hist,FolderName);
}




//_____________________________________________________________________________
void TStnModule::DeleteHistograms(TFolder* Folder) {
  // internal method...
  fBookHist->DeleteHistograms(Folder);

  // if (((long int) Folder) == -1) Folder = fFolder;

  // TObject  *o1;

  // TIter    it1(Folder->GetListOfFolders());

  // while ((o1 = it1.Next())) {
  //   if (o1->InheritsFrom("TFolder")) {
  //     DeleteHistograms((TFolder*) o1);
  //   }
  //   else if (o1->InheritsFrom("TH1")) {
  //     Folder->Remove(o1);
  //     delete o1;
  //   }
  // }
}


//_____________________________________________________________________________
void TStnModule::Clear(const char* Opt) {
}

//_____________________________________________________________________________
void TStnModule::Delete(const char* Opt) {
  // whatever it means

  if (strcmp(Opt,"hist") == 0) {
					// delete all the histograms
    DeleteHistograms(fFolder);
  }
  else {
    fFolder->Clear();
  }
}

//_____________________________________________________________________________
void TStnModule::Print(const char* Opt) const {
}


//-----------------------------------------------------------------------------
// Ip1, Ip2 
//-----------------------------------------------------------------------------
int TStnModule::FitHistogram(TH1* Hist, fit_result_t* Fr, float XMin, float XMax, int MinSum) {

  //  fit_result_t* fr = &fFr[Ip2][Ip1];

  Fr->chi2dof = -1;
      
  for (int ip=0; ip<3; ip++) {
    Fr->p[ip] = 0;
    Fr->e[ip] = -1;
  }

  // TH1F* h = fHist->h_dt05[Ip2][Ip1];

  int nbins     = Hist->GetNbinsX();
  int integral  = Hist->Integral(1,nbins);
      
  if (integral < MinSum) {
    std::cout << std::format("ERROR: integral:{} < {}. BAIL OUT\n",integral,MinSum);
    return -1;
  }

  // find max bin

  int   imax = -1;
  float qmax = -1;
  for (int i=0; i<nbins; i++) {
    float q = Hist->GetBinContent(i+1);
    if (q > qmax) {
      imax = i+1;
      qmax = q;
    }
  }

  if (qmax < 3) {
    std::cout << std::format("ERROR: qmax:{} < 3. BAIL OUT\n",qmax);
    return -2;
  }

  // estimate integral of the expected gaussian

  int i=0;
  //  int imin(0), imax(0);

  float sum = qmax;
  
  while (1) {
    i += 1;
                                        // check to the right of the maximum
    int iplus  = imax+i;
    if (iplus <= nbins) {
      float qplus = Hist->GetBinContent(iplus);
      if (qplus/qmax > 0.2) {
        sum += qplus;
      }
      else {
        // done
        iplus = nbins+1;
      }
    }
                                        // check the left side
    int iminus = imax-i;
    if (iminus > 0) {
                                        // bins start from 1
      
      float qminus = Hist->GetBinContent(iminus);
      if (qminus/qmax > 0.2) {
        sum += qminus;
      }
      else {
        // done
        iminus = -1;
      }
    }
    if ((iplus > nbins) and (iminus < 0)) break;
  }

  if (sum < MinSum) {
    std::cout << std::format("ERROR: sum:{} < {}. BAIL OUT\n",sum,MinSum);
    return -3;
  }

  float t0 = Hist->GetBinCenter(imax);
  
  // for some reason, "sq" is required for tfr be defined
  float tmin{t0-50}, tmax{t0+50};
  if (XMax > XMin) {
    tmin = XMin;
    tmax = XMax;
  }
  
  TFitResultPtr tfr = Hist->Fit("gaus","sql","",tmin,tmax);
  
  if ((! tfr->IsValid()) or tfr->IsEmpty()) {
    // assume tha all indices are in the name/title
    std::cout << std::format("# FIT ERROR: Hist->name:{} Hist->title:{}\n",Hist->GetName(),Hist->GetTitle());
    return -4;
  }

  Fr->chi2dof = tfr->Chi2()/tfr->Ndf();
  double sf    = sqrt(Fr->chi2dof);
          
  for (int ip=0; ip<3; ip++) {
    Fr->p[ip] = tfr->Parameter(ip);
    Fr->e[ip] = tfr->Error(ip)*sf;
  }
  return 0;
}
