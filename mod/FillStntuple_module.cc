//--------------------------------------------------------------------------
// Description:
// -----------
// Class FillStntuple : fills Stntuple
//
//
//------------------------------------------------------------------------
#ifdef __GNUG__
#pragma implementation
#endif

#include "TTree.h"
#include "TBranchElement.h"
#include "TBranchObject.h"
#include "TLeafObject.h"
#include "TSystem.h"

#include <Stntuple/obj/TStnNode.hh>
#include <Stntuple/obj/TStnEvent.hh>
#include <Stntuple/obj/TStnDataBlock.hh>
#include <Stntuple/obj/TStnDBManager.hh>
#include <Stntuple/obj/TStnHeaderBlock.hh>

#include "Stntuple/mod/StntupleModule.hh"

namespace mu2e {

class FillStntuple : public StntupleModule {
//------------------------------------------------------------------------------
//  data members
//------------------------------------------------------------------------------
protected:
  Int_t             fLastRun;		// last run with events
//------------------------------------------------------------------------------
// function members
//------------------------------------------------------------------------------
public:
					// constructors and destructor

  FillStntuple(fhicl::ParameterSet const& PSet);

  ~FillStntuple();
//-----------------------------------------------------------------------------
// functions of the module
//-----------------------------------------------------------------------------
  int     ProcessNewRun      (int RunNumber);
//-----------------------------------------------------------------------------
// overloaded virtual functions of EDFilter
//-----------------------------------------------------------------------------
  virtual void beginRun(const art::Run&   r);
  virtual void endRun  (const art::Run&   r);
  virtual void analyze (const AbsEvent&   e);

  //  ClassDef(FillStntuple,0)
};


// ClassImp(FillStntuple)
//------------------------------------------------------------------------------
// constructors
//------------------------------------------------------------------------------
FillStntuple::FillStntuple(fhicl::ParameterSet const& PSet): 
  StntupleModule   (PSet.get<fhicl::ParameterSet>("THistModule"),"FillStntuple")
{
  fLastRun = -1;
  TTree::SetMaxTreeSize(8000000000LL);
}

//------------------------------------------------------------------------------
FillStntuple::~FillStntuple() {
}

//------------------------------------------------------------------------------
void FillStntuple::beginRun(const art::Run &  aRun) {

  THistModule::beforeBeginRun(aRun);

  TModule::fRun  = &aRun;

  int runnum = aRun.run();
  if (runnum != fLastRun) {
					// create new subdirectory and store 
					// calibration constants in there
    // [2025-05-12 PM bypass for the time being]  ProcessNewRun(runnum);
    fLastRun = runnum;
  }

  THistModule::afterBeginRun(aRun);
}

//------------------------------------------------------------------------------
void FillStntuple::endRun(const art::Run &  Rn) {

  THistModule::beforeEndRun(Rn);
  THistModule::afterEndRun (Rn);
}

//------------------------------------------------------------------------------
  Int_t FillStntuple::ProcessNewRun(int RunNumber)  {
  // create subdirectory with the name run_xxxxxxxx to store database-type
  // constants for this run

  char dir_name[200];

  TDirectory* olddir = gDirectory;
					// 

  sprintf(dir_name,"run_%i",RunNumber);

					// make sure that "db" directory exists
  if (fgFile->GetKey("db") == 0) {
    fgFile->mkdir("db");
  }
  fgFile->cd("db");
					// curent directory - "db"
					// check if subdirectory for this 
					// run has already been created 

  if (gDirectory->GetKey(dir_name) != 0) return 0;

  gDirectory->mkdir(dir_name);
  gDirectory->cd(dir_name);
//-----------------------------------------------------------------------------
// directory is created, time to store there whatever we need
// start from the beam position (note, that we cd'ed into the new directory)
//-----------------------------------------------------------------------------
  TStnDBManager* dbm = TStnDBManager::Instance();

  dbm->Write();

  gDirectory = olddir;
  gDirectory->cd();

  return 0;
}

//------------------------------------------------------------------------------
void FillStntuple::analyze(const AbsEvent& anEvent) {
  // it only fills the tree

  TTree* tree;
  //  TFile* file;
  TFile* old_file;

  char line[100];

  int rc, run_section;

  THistModule::beforeEvent(anEvent);

  static int old_rs = -1;
//-----------------------------------------------------------------------------
// at this point, all individual blocks are filled
// initialize information cross-linking the blocks
//-----------------------------------------------------------------------------
  TIter it(Event()->GetListOfNodes());

  unsigned long rtime = (unsigned long)(gSystem->Now());
  while(TStnNode* node = (TStnNode*) it.Next()) {
    TStnDataBlock* block = node->GetDataBlock();
    rc = block->ResolveLinks((AbsEvent*) &anEvent,0);
    if (rc != 0) {
					// pass all the messages/warnings
					// to the error logger
      TIter it(block->MessageList());
      // while (TObjString* mess = (TObjString*) it.Next()) {
      // 	ERRLOG(ELwarning,mess->GetString().Data()) << endmsg;
      // }
					// and don't forget to delete the 
					// messages after printing them
      block->MessageList()->Delete();
    }
  }
  // time to do resolve links in ms
  rtime = (unsigned long)(gSystem->Now()) - rtime;
  // add processing time to that from InitStntuple
  TStnHeaderBlock* fHeaderBlock = 
    (TStnHeaderBlock*) Event()->GetDataBlock("HeaderBlock");
  if(fHeaderBlock) {
    float t = fHeaderBlock->CpuTime(); // in s
    t += float(rtime)/1000.0;
    int it = int(t*10.0); // store time in 10*s
    if(it>(1<<24)) it=((1<<24)-1);
    int speed = (fHeaderBlock->fCpu & 0xFF);
    fHeaderBlock->fCpu = (it<<8 | speed);
  }

//-----------------------------------------------------------------------------
// if we need to close a file and to write a new one, make such a decision here
// fgOpenNextFile will be reset by FillStntuple after it writes the file
//-----------------------------------------------------------------------------
  run_section = anEvent.subRun();

  int mbytes_written = (int) (fgFile->GetBytesWritten()/1000000);
  if (mbytes_written >= fgMaxFileSize) {
    if (run_section != old_rs) {
      THistModule::fgOpenNextFile = 1;
    }
  }
//-----------------------------------------------------------------------------
// now fill the tree, however before doing it check if the output file is not
// too large and it is time to switch to writing a new file
//-----------------------------------------------------------------------------
  if (THistModule::fgOpenNextFile) {
					// open the new file
    THistModule::fgFileNumber++;
    sprintf(line,"_%i",THistModule::fgFileNumber);

    TString fn = THistModule::fgFileName;
    fn += line;
					// open new file 

    old_file = THistModule::fgFile;
    THistModule::fgFile = new TFile(fn.Data(),"RECREATE");

					// and create the tree in it

    tree = new TTree ("STNTUPLE","STNTUPLE");

					// loop over the branches of still
					// opened fgTree and duplicate them
					// in `tree'

    TObjArray* list = fgTree->GetListOfBranches();
    TIter it(list);

    TBranch*      input_branch;
    TBranch*      output_branch;
    TLeafObject*  leaf;

    const char* class_name;
    const char* input_branch_name;

    while ((input_branch = (TBranch*) it.Next())) {
      input_branch_name = input_branch->GetName();

      if (strcmp(input_branch->ClassName(),"TBranchElement") == 0) {
	class_name = ((TBranchElement*)input_branch)->GetClassName();
      }
      else {
	leaf = (TLeafObject*) ((TBranchObject*)input_branch)->GetLeaf(input_branch_name);
	class_name = leaf->GetTypeName();
      }
      TStnNode* node = fgEvent->FindNode(input_branch_name);

					// preserve split level

      output_branch = tree->Branch(input_branch_name,class_name,
				   node->GetDataBlockAddress(),
				   input_branch->GetBasketSize(),
				   input_branch->GetSplitLevel());

      output_branch->SetCompressionLevel(input_branch->GetCompressionLevel());
      output_branch->SetAutoDelete(kFALSE);

					// and also redefine branch in the node
      node->SetBranch(output_branch);
					// store calib consts for the last event
      ProcessNewRun(anEvent.run());
    }
					// close the old file
    old_file->Write();
    delete old_file;
					// redefine the static variables
    THistModule::fgTree = tree;

    THistModule::fgOpenNextFile = 0;
  }
					// and finally fill the tree
					// this is the first entry in the 
					// new file
  fgTree->Fill();
  old_rs = run_section;

  THistModule::afterEvent(anEvent);

}

} // end namespace mu2e

DEFINE_ART_MODULE(mu2e::FillStntuple)
