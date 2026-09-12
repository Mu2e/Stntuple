//------------------------------------------------------------------------------
// Implementation of an output module.
//
// Drop and Keep list strategy
// ---------------------------
// 1) The keep list is ignored as long as it is empty and everything
//    except for the objects mentioned in the drop list are kept.
// 2) If the keep list contains one or more elements the drop list is
//    ignored (even if it has something in it) and all blocks except
//    the one mentioned in the keep list ared dropped.
//------------------------------------------------------------------------------
#ifndef TStnOutputModule_hh
#define TStnOutputModule_hh

#include "TStnModule.hh"

class TStnEvent;
class TStnNode;
class TTree;
class TFile;

class TStnOutputModule: public TStnModule {
//-----------------------------------------------------------------------------
//  data members
//-----------------------------------------------------------------------------
protected:
  TString               fFileName;
  TFile*                fFile;		// output file
  TTree*                fTree;		// pointer to the output tree
  Int_t                 fFileNumber;	// sequential file number in a job
  Int_t                 fMaxFileSize;	// max file size in MBytes
  TObjArray*            fDropList;	// list of data blocks (names) to drop
  TObjArray*            fKeepList;	// list of data blocks (names) to keep
  int                   fLoRun;         // low and high run and subrun numbers
  int                   fLoSubrun;      // %08i is OK for an int
  int                   fHiRun;        
  int                   fHiSubrun;        
  TString               fDatasetID;     // dataset ID
  TString               fProject;       // project (description)
//-----------------------------------------------------------------------------
//  functions
//-----------------------------------------------------------------------------
public:
  TStnOutputModule(const char* FileName  = "stntuple_output.root");

  virtual ~TStnOutputModule();

  virtual int BeginJob      () override;
  virtual int BeginRun      () override;
  virtual int Event         (Int_t i) override;
  virtual int EndRun        () override;
  virtual int EndJob        () override;
//-----------------------------------------------------------------------------
// accessors
//-----------------------------------------------------------------------------
  TTree*         GetOutputTree () { return fTree; }
  Int_t          GetMaxFileSize() { return fMaxFileSize; }
  TObjArray*     GetDropList   () { return fDropList; }
  TObjArray*     GetKeepList   () { return fKeepList; }
  const TString& DatasetID     () const { return fDatasetID;  }
//-----------------------------------------------------------------------------
// modifiers
//-----------------------------------------------------------------------------
  void        SetMaxFileSize(int         size   ) { fMaxFileSize = size   ; }
  void        SetDatasetID  (const char* DsID   ) { fDatasetID   = DsID   ; }
  void        SetProject    (const char* Project) { fProject     = Project; }
//-----------------------------------------------------------------------------
// other methods
//-----------------------------------------------------------------------------
  void        DropDataBlock (const char* Name);
  void        KeepDataBlock (const char* Name);

  Int_t       OpenNewFile   (const char* Filename );

  ClassDefOverride(TStnOutputModule,0)
};
#endif
