///////////////////////////////////////////////////////////////////////////////
// Class StntupleMaker : fills Stntuple (P.Murat)
// ------------------------------------------
// order of the data blocks is essential - they are filled according to the
// order in which they are declared...
//
///////////////////////////////////////////////////////////////////////////////

#ifdef __GNUG__
#pragma implementation
#endif

#include <string>
#include <cstdio>

#include "fhiclcpp/ParameterSet.h"

#include <assert.h>
#include <iostream>
#include <iomanip>

#include "TNamed.h"
#include "TH1.h"
#include "TString.h"
#include "TProfile.h"
#include "TFolder.h"
#include "TSystem.h"

#include "Offline/GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "Offline/GlobalConstantsService/inc/PhysicsParams.hh"

#include "Stntuple/obj/TStnEvent.hh"
#include "Stntuple/obj/TStnNode.hh"
#include "Stntuple/obj/TStnErrorLogger.hh"
#include "Stntuple/obj/TStnTrackBlock.hh"
#include "Stntuple/obj/TStnHelixBlock.hh"
#include "Stntuple/obj/TStrawHitBlock.hh"
#include "Stntuple/obj/TCalDataBlock.hh"
#include "Stntuple/obj/TGenpBlock.hh"
#include "Stntuple/obj/TSimpBlock.hh"
#include "Stntuple/obj/TStepPointMCBlock.hh"
#include "Stntuple/obj/TStnHeaderBlock.hh"
//-----------------------------------------------------------------------------
// nothing new: link-wise, TModule depends on TAnaDump
//-----------------------------------------------------------------------------
#include "Stntuple/print/TAnaDump.hh"
#include "Stntuple/mod/StntupleModule.hh"
#include "Stntuple/mod/StntupleGlobals.hh"

#include "Stntuple/mod/InitCrvPulseBlock.hh"
#include "Stntuple/mod/InitCrvClusterBlock.hh"
#include "Stntuple/mod/InitGenpBlock.hh"
#include "Stntuple/mod/InitHeaderBlock.hh"
#include "Stntuple/mod/InitHelixBlock.hh"
#include "Stntuple/mod/InitSimpBlock.hh"
#include "Stntuple/mod/InitStrawHitBlock.hh"
#include "Stntuple/mod/InitStepPointMCBlock.hh"
#include "Stntuple/mod/InitTrackBlock.hh"
#include "Stntuple/mod/InitTrackBlock_KK.hh"
#include "Stntuple/mod/InitTrackSeedBlock.hh"
#include "Stntuple/mod/InitTrackStrawHitBlock.hh"
#include "Stntuple/mod/InitTriggerBlock.hh"
#include "Stntuple/mod/InitTimeClusterBlock.hh"

#include "Stntuple/mod/InitStntupleDataBlocks.hh"
#include "Stntuple/mod/StntupleUtilities.hh"

#include "Stntuple/base/TNamedHandle.hh"
#include "Stntuple/alg/TStntuple.hh"

#include "Offline/TrkReco/inc/DoubletAmbigResolver.hh"
#include "Offline/MCDataProducts/inc/GenId.hh"
#include "Offline/RecoDataProducts/inc/HelixSeed.hh"
// #include "TrkDiag/inc/KalDiag.hh"

using namespace std; 

// ClassImp(StntupleMaker)

static const char rcsid[] = "$Name:  $";
// stntuple_get_version is autogenerated, see Stntuple/scripts/create_print_header_routine.sh
void stntuple_get_version(char* ver, char* test);
namespace mu2e {
class StntupleMaker : public StntupleModule {
//------------------------------------------------------------------------------
//  data members
//------------------------------------------------------------------------------
protected:
					// process name, default - PROD
  std::string              fProcessName;
//-----------------------------------------------------------------------------
// switches for individual branches
//-----------------------------------------------------------------------------
  int                      fMakeCalData;
  int                      fMakeClusters;
  int                      fMakeGenp;
  int                      fMakePid;
  int                      fMakeSimp;         // 0:dont store, 1:all;
  int                      fMakeStepPointMC;
  int                      fMakeStrawHits;
  int                      fMakeStrawWaveforms;
  int                      fMakeTracks;
  int                      fMakeTrackStrawHits;
  int                      fMakeTimeClusters;
  int                      fMakeHelices;
  int                      fMakeTrackSeeds;
  int                      fMakeTrigger;
  int                      fMakeCrvPulses;
  int                      fMakeCrvClusters;
//-----------------------------------------------------------------------------
// module parameters
//-----------------------------------------------------------------------------
  art::InputTag            fGenpCollTag;
  art::InputTag            fSimpCollTag;

  art::InputTag            fChCollTag;
  art::InputTag            fShCollTag;
  string                   fStrawDigiCollTag;
  string                   fSdwfCollTag;
  art::InputTag            fSdmcCollTag;

  string                   fCrvRecoPulseCollTag;            //
  string                   fCrvCoincidenceCollTag;          //
  string                   fCrvCoincidenceClusterCollTag;   //

  art::InputTag            fVdhCollTag;                     // hits on virtual detectors (StepPointMCCollection)

  vector<string>           fTClBlockName;                   // time cluster block names
  vector<string>           fTClCollTag;                     // time cluster coll tags

  vector<string>           fHelixBlockName;
  vector<string>           fHelixSeedCollTag;
  vector<art::InputTag>    fHelixKsCollTag;

  art::InputTag            fPbiTag;
  art::InputTag            fPrimaryParticleTag;
  string                   fTriggerResultsTag;
  int                      fNTriggerBits;

  vector<string>           fKsfBlockName;
  vector<string>           fKsfCollTag;

  vector<string>           fTrackBlockName;
  vector<string>           fTrackCollTag;
  int                      fTrackFitType;

  vector<string>           fTrackTsBlockName;  // for each track block, the tag  of the corresponding TrackSeed coll
  vector<string>           fTrackTsCollTag;    // for each track block, the name of the corresponding TrackSeed block

  vector<string>           fTrackHsBlockName;  // KinKal fits: for each track block, 
                                               // the name of the corresponding Stntuple HelixSeed block

  vector<string>           fTciCollTag;        // collection produced by TrackCaloIntersection module
  vector<string>           fTcmCollTag;        // collection produced by TrackCaloMatching     module
  vector<string>           fTrkQualCollTag;    // collection produced by TrackQuality          module

  vector<string>           fPidBlockName;
  vector<string>           fPidCollTag;

  vector<string>           fTrackSHBlockName;

  vector<string>           fSpmcBlockName;
  vector<art::InputTag>    fSpmcCollTag;
  vector<string>           fStatusG4Tag;

  string                   fCaloCrystalHitMaker;
  string                   fCaloClusterMaker;
//-----------------------------------------------------------------------------
// initialization of various data blocks
//-----------------------------------------------------------------------------
  StntupleInitCrvPulseBlock*    fInitCrvPulseBlock;
  StntupleInitCrvClusterBlock*  fInitCrvClusterBlock;
  StntupleInitGenpBlock*        fInitGenpBlock;
  stntuple::InitHeaderBlock*    fInitHeaderBlock;
  StntupleInitSimpBlock*        fInitSimpBlock;
  stntuple::InitStrawHitBlock*  fInitStrawHitBlock;
  StntupleInitTriggerBlock*     fInitTriggerBlock;
  TObjArray*                    fInitTrackStrawHitBlock;
  TObjArray*                    fInitHelixBlock;
  TObjArray*                    fInitTrackBlock;
  TObjArray*                    fInitTrackSeedBlock;
  TObjArray*                    fInitStepPointMCBlock;
  TObjArray*                    fInitTimeClusterBlock;
//-----------------------------------------------------------------------------
// cut-type parameters
//-----------------------------------------------------------------------------
  GenId                    fGenId          ;     // generated process ID
  int                      fPdgId          ;     // PDG ID of the simparticle to be stored, 0 by default
  
  double                   fMinTActive     ;     // start of the active window
  double                   fMinECrystal    ;     // 
  double                   fMinSimpMomentum;     // min tot momentum of a particle to be stored in SIMP block
  double                   fSimpMaxZ       ;     // max Z of a particle to be stored in SIMP block
  int                      fMinNStrawHits  ;     // min number of straw hits produced by a SimParticle

  string                   fCutHelixSeedCollTag; // helix collection to cut on
  int                      fMinNHelices        ; // min number of helices (for cosmics)

  TNamed*                  fVersion;

  TNamedHandle*            fDarHandle;

  DoubletAmbigResolver*    fDar;

//------------------------------------------------------------------------------
// function members
//------------------------------------------------------------------------------
public:
					// constructors and destructor

  StntupleMaker(fhicl::ParameterSet const& pset);

  ~StntupleMaker();
//-----------------------------------------------------------------------------
// functions of the module
//-----------------------------------------------------------------------------
  void GetDefTrackCollName(char* Name);

					// ****** setters

//-----------------------------------------------------------------------------
// overloaded virtual functions of EDFilter
//-----------------------------------------------------------------------------
  virtual void beginRun(const art::Run& ARun);
  virtual void endRun  (const art::Run& ARun);
  virtual void beginJob();
  virtual void endJob  ();
  virtual void analyze (const AbsEvent& event);

  //  ClassDef(StntupleMaker,0)
};



//------------------------------------------------------------------------------
// constructors
//------------------------------------------------------------------------------
StntupleMaker::StntupleMaker(fhicl::ParameterSet const& PSet): 
  StntupleModule             (PSet,                   "StntupleMaker"         )
  , fProcessName             (PSet.get<string>        ("processName"         ))

  , fMakeCalData             (PSet.get<int>           ("makeCalData"         ))
  , fMakeClusters            (PSet.get<int>           ("makeClusters"        ))
  , fMakeGenp                (PSet.get<int>           ("makeGenp"            ))
  , fMakePid                 (PSet.get<int>           ("makePid"             ))
  , fMakeSimp                (PSet.get<int>           ("makeSimp"            ))
  , fMakeStepPointMC         (PSet.get<int>           ("makeStepPointMC"     ))
  , fMakeStrawHits           (PSet.get<int>           ("makeStrawHits"       ))
  , fMakeStrawWaveforms      (PSet.get<int>           ("makeStrawWaveforms"  ))
  , fMakeTracks              (PSet.get<int>           ("makeTracks"          ))
  , fMakeTrackStrawHits      (PSet.get<int>           ("makeTrackStrawHits"  ))
  , fMakeTimeClusters        (PSet.get<int>           ("makeTimeClusters"    ))
  , fMakeHelices             (PSet.get<int>           ("makeHelices"         ))
  , fMakeTrackSeeds          (PSet.get<int>           ("makeTrackSeeds"      ))
  , fMakeTrigger             (PSet.get<int>           ("makeTrigger"         ))
  , fMakeCrvPulses           (PSet.get<int>           ("makeCrvPulses"       ))
  , fMakeCrvClusters         (PSet.get<int>           ("makeCrvClusters"     ))
  
  , fGenpCollTag             (PSet.get<art::InputTag> ("genpCollTag"         ))
  , fSimpCollTag             (PSet.get<art::InputTag> ("simpCollTag"         ))

  , fChCollTag               (PSet.get<art::InputTag> ("chCollTag"     ))
  , fShCollTag               (PSet.get<art::InputTag> ("shCollTag"     ))
  , fStrawDigiCollTag        (PSet.get<string>        ("strawDigiCollTag"    ))
  , fSdwfCollTag             (PSet.get<string>        ("sdwfCollTag"         ))
  , fSdmcCollTag             (PSet.get<art::InputTag> ("strawDigiMCCollTag"  ))

  , fCrvRecoPulseCollTag         (PSet.get<string>    ("crvRecoPulseCollTag"         ))
  , fCrvCoincidenceCollTag       (PSet.get<string>    ("crvCoincidenceCollTag"       ))
  , fCrvCoincidenceClusterCollTag(PSet.get<string>    ("crvCoincidenceClusterCollTag"))

  , fVdhCollTag              (PSet.get<art::InputTag> ("vdHitsCollTag"       ))
  , fTClBlockName            (PSet.get<vector<string>>("timeClusterBlockName"))
  , fTClCollTag              (PSet.get<vector<string>>("timeClusterCollTag"  ))
  , fHelixBlockName          (PSet.get<vector<string>>("helixBlockName"      ))
  , fHelixSeedCollTag        (PSet.get<vector<string>>("helixCollTag"        ))
  , fHelixKsCollTag          (PSet.get<vector<art::InputTag>>("helixKsCollTag"))

  , fPbiTag                  (PSet.get<art::InputTag> ("pbiTag"    ))
  , fPrimaryParticleTag      (PSet.get<art::InputTag> ("primaryParticleTag"  ))

  , fTriggerResultsTag       (PSet.get<string>        ("triggerResultsTag"   ))
  , fNTriggerBits            (PSet.get<int>           ("nTriggerBits"        ))

  , fKsfBlockName            (PSet.get<vector<string>>("trackSeedBlockName"  ))
  , fKsfCollTag              (PSet.get<vector<string>>("trackSeedCollTag"    ))

  , fTrackBlockName          (PSet.get<vector<string>>("trackBlockName"      ))
  , fTrackCollTag            (PSet.get<vector<string>>("trackCollTag"        ))
  , fTrackFitType            (PSet.get<int>           ("trackFitType"        ))

  , fTrackTsBlockName        (PSet.get<vector<string>>("trackTsBlockName"    ))
  , fTrackTsCollTag          (PSet.get<vector<string>>("trackTsCollTag"      ))

  , fTrackHsBlockName        (PSet.get<vector<string>>("trackHsBlockName"    ))

  , fTciCollTag              (PSet.get<vector<string>>("tciCollTag"          ))
  , fTcmCollTag              (PSet.get<vector<string>>("tcmCollTag"          ))
  , fTrkQualCollTag          (PSet.get<vector<string>>("trkQualCollTag"      ))
  , fPidBlockName            (PSet.get<vector<string>>("pidBlockName"        ))
  , fPidCollTag              (PSet.get<vector<string>>("pidCollTag"          ))
  , fTrackSHBlockName        (PSet.get<vector<string>>("trackSHBlockName"    ))
  
  , fSpmcBlockName           (PSet.get<vector<string>>       ("spmcBlockName"       ))
  , fSpmcCollTag             (PSet.get<vector<art::InputTag>>("spmcCollTag"         ))
  , fStatusG4Tag             (PSet.get<vector<string>>       ("statusG4Tag"         ))

  , fCaloCrystalHitMaker     (PSet.get<string>        ("caloCrystalHitsMaker"))
  , fCaloClusterMaker        (PSet.get<string>        ("caloClusterMaker"    ))

  , fGenId(GenId::findByName (PSet.get<std::string>   ("genId"               ),"unknown"))
  , fPdgId                   (PSet.get<int>           ("pdgId"               ))

  , fMinTActive              (PSet.get<double>        ("minTActive"          ))
  , fMinECrystal             (PSet.get<double>        ("minECrystal"         ))
  , fMinSimpMomentum         (PSet.get<double>        ("minSimpMomentum"     ))
  , fSimpMaxZ                (PSet.get<double>        ("simpMaxZ"            ))
  , fMinNStrawHits           (PSet.get<int>           ("minNStrawHits"       ))
  , fCutHelixSeedCollTag     (PSet.get<string>        ("cutHelixSeedCollTag" ))
  , fMinNHelices             (PSet.get<int>           ("minNHelices"         ))
{

  char  ver[20], text[200];
  stntuple_get_version(ver,text);

  fVersion      = new TNamed(ver,text);
  TModule::fFolder->Add(fVersion);

  fInitCrvPulseBlock    = nullptr;
  fInitCrvClusterBlock  = nullptr;
  fInitGenpBlock        = nullptr;
  fInitHeaderBlock      = nullptr;
  fInitSimpBlock        = nullptr;
  fInitStrawHitBlock    = nullptr;
  fInitTriggerBlock     = nullptr;

  fInitStepPointMCBlock = new TObjArray();
  fInitStepPointMCBlock->SetOwner(kTRUE);

  fInitHelixBlock       = new TObjArray();
  fInitHelixBlock->SetOwner(kTRUE);

  fInitTrackBlock       = new TObjArray();
  fInitTrackBlock->SetOwner(kTRUE);

  fInitTrackSeedBlock   = new TObjArray();
  fInitTrackSeedBlock->SetOwner(kTRUE);

  fInitTrackStrawHitBlock = new TObjArray();
  fInitTrackStrawHitBlock->SetOwner(kTRUE);

  fInitTimeClusterBlock = new TObjArray();
  fInitTimeClusterBlock->SetOwner(kTRUE);
  fDar                  = new DoubletAmbigResolver (PSet.get<fhicl::ParameterSet>("DoubletAmbigResolver"),0.,0,0);
  fDarHandle            = new TNamedHandle("DarHandle",fDar);
  // fKalDiag              = new KalDiag     (PSet.get<fhicl::ParameterSet>("KalDiag",fhicl::ParameterSet()));
  // fKalDiagHandle        = new TNamedHandle("KalDiagHandle"      ,fKalDiag);

  fFolder->Add(fDarHandle);

}


//------------------------------------------------------------------------------
StntupleMaker::~StntupleMaker() {
  delete fDar;
  delete fDarHandle;
  delete fVersion;

  if (fInitCrvPulseBlock  ) delete fInitCrvPulseBlock;
  if (fInitCrvClusterBlock) delete fInitCrvClusterBlock;
  if (fInitSimpBlock      ) delete fInitSimpBlock;

  delete fInitStepPointMCBlock;
  delete fInitHelixBlock;
  delete fInitTrackBlock;
  delete fInitTrackSeedBlock;
  delete fInitHeaderBlock;
}

//------------------------------------------------------------------------------
void StntupleMaker::beginRun(const art::Run& aRun) {

  static int first_begin_run = 1;

  THistModule::beforeBeginRun(aRun);

  if (first_begin_run) {
//-----------------------------------------------------------------------------
// if we runnning stnmaker_prod.exe, save revision of the TCL file in STNTUPLE
//-----------------------------------------------------------------------------
    first_begin_run = 0;
    const char* c = gSystem->Getenv("STNMAKER_PROD_TCL");
    if (c) TModule::fFolder->Add(new TNamed("STNMAKER_PROD_TCL",c));
    else   TModule::fFolder->Add(new TNamed("STNMAKER_PROD_TCL","unknown"));
  }
//-----------------------------------------------------------------------------
// StepPointMC collections - set mbtime - have to do that at beginRun()
//-----------------------------------------------------------------------------
  float mbtime = GlobalConstantsHandle<PhysicsParams>()->getNominalDRPeriod();

  if (fMakeStepPointMC) {
    int nblocks = fSpmcBlockName.size();

    for (int i=0; i<nblocks; i++) {
      StntupleInitStepPointMCBlock* init_block = static_cast<StntupleInitStepPointMCBlock*> (fInitStepPointMCBlock->At(i));
      init_block->SetMbTime(mbtime);
    }
  }

  THistModule::afterBeginRun(aRun);
}

//------------------------------------------------------------------------------
void StntupleMaker::endRun(const art::Run& aRun ) {
  THistModule::beforeEndRun(aRun);
  THistModule::afterEndRun (aRun);
}


//------------------------------------------------------------------------------
void StntupleMaker::endJob() {

  THistModule::beforeEndJob();
  THistModule::afterEndJob ();

}

//------------------------------------------------------------------------------
void StntupleMaker::beginJob() {

  int split_mode, compression_level, buffer_size;

  THistModule::beforeBeginJob();

  // create data blocks and branches

  fgStntupleFolder->Add(new TNamed("ProcessName"     ,fProcessName));

					// for the moment do it by hands...
					// create default branches to go into 
					// STNTUPLE

  split_mode        = THistModule::SplitLevel();
  compression_level = THistModule::CompressionLevel();
  buffer_size       = THistModule::BufferSize();
//-----------------------------------------------------------------------------
// header block is always there
//-----------------------------------------------------------------------------
  fInitHeaderBlock = new stntuple::InitHeaderBlock();
  fInitHeaderBlock->SetPbiTag(fPbiTag);
  fInitHeaderBlock->SetChCollTag(fChCollTag);
  fInitHeaderBlock->SetShCollTag(fShCollTag);

  AddDataBlock("HeaderBlock","TStnHeaderBlock",
	       fInitHeaderBlock,
	       THistModule::BufferSize(),
	       0, // 99,                          // fSplitMode.value()
	       THistModule::CompressionLevel());

  //  SetResolveLinksMethod("HeaderBlock",StntupleInitMu2eHeaderBlockLinks);
//-----------------------------------------------------------------------------
// calorimeter hit data
// this is not RAW hit data yet...
//-----------------------------------------------------------------------------
  if (fMakeCalData) {
    TStnDataBlock* cal_data;

    cal_data = AddDataBlock("CalDataBlock","TCalDataBlock",
			    StntupleInitMu2eCalDataBlock,
			    buffer_size,
			    split_mode,
			    compression_level);
    if (cal_data) {
      cal_data->AddCollName("mu2e::CaloHitCollection",fCaloCrystalHitMaker.data());
    }
  }
//-----------------------------------------------------------------------------
// calorimeter clusters 
//-----------------------------------------------------------------------------
  if (fMakeClusters) {
    TStnDataBlock* db = AddDataBlock("ClusterBlock",
				     "TStnClusterBlock",
				     StntupleInitMu2eClusterBlock,
				     buffer_size,
				     split_mode,
				     compression_level);
    if (db) {
      db->AddCollName("mu2e::CaloClusterCollection",fCaloClusterMaker.data());
      SetResolveLinksMethod("ClusterBlock",StntupleInitMu2eClusterBlockLinks);
    }
  }
//-----------------------------------------------------------------------------
// CRV
//-----------------------------------------------------------------------------
  if (fMakeCrvPulses) {

    fInitCrvPulseBlock = new StntupleInitCrvPulseBlock();
    fInitCrvPulseBlock->SetCrvRecoPulseCollTag(fCrvRecoPulseCollTag);
    fInitCrvPulseBlock->SetCrvCoincidenceCollTag(fCrvCoincidenceCollTag);
    fInitCrvPulseBlock->SetCrvCoincidenceClusterCollTag(fCrvCoincidenceClusterCollTag);

    AddDataBlock("CrvPulseBlock","TCrvPulseBlock",
		 fInitCrvPulseBlock,
		 buffer_size,
		 split_mode,
		 compression_level);
  }

  if (fMakeCrvClusters) {

    fInitCrvClusterBlock = new StntupleInitCrvClusterBlock();
    fInitCrvClusterBlock->SetCrvRecoPulseCollTag(fCrvRecoPulseCollTag);
    fInitCrvClusterBlock->SetCrvCoincidenceClusterCollTag(fCrvCoincidenceClusterCollTag);

    AddDataBlock("CrvClusterBlock","TCrvClusterBlock",
		 fInitCrvClusterBlock,
		 buffer_size,
		 split_mode,
		 compression_level);
  }
//-----------------------------------------------------------------------------
// generator particles 
//-----------------------------------------------------------------------------
  if (fMakeGenp) {
    fInitGenpBlock = new StntupleInitGenpBlock();
    fInitGenpBlock->SetGenpCollTag(fGenpCollTag);
    fInitGenpBlock->SetGenProcessID (fGenId.id());
    fInitGenpBlock->SetPdgID        (fPdgId     );

    AddDataBlock("GenpBlock","TGenpBlock",fInitGenpBlock,buffer_size,split_mode,compression_level);
  }
//--------------------------------------------------------------------------------
// helix data
//--------------------------------------------------------------------------------
  if (fMakeHelices) {
    int nblocks    = fHelixBlockName.size();
    int nks_blocks = fKsfBlockName.size();

    for (int i=0; i<nblocks; i++) {
      const char* block_name = fHelixBlockName[i].data();

      StntupleInitHelixBlock* init_block(nullptr);
      init_block = new StntupleInitHelixBlock();

      fInitHelixBlock->Add(init_block);

      init_block->SetHSeedCollTag(fHelixSeedCollTag[i]);
      init_block->SetSdmcCollTag (fSdmcCollTag );
      init_block->SetTclBlockName(fTClBlockName   [i]);
      init_block->SetTrackFitType(fTrackFitType);

      if (i < nks_blocks) {
        init_block->SetKsCollTag   (fHelixKsCollTag[i]);
        init_block->SetKsfBlockName(fKsfBlockName  [i]);
      }

      AddDataBlock(block_name,"TStnHelixBlock",init_block,
                   buffer_size,split_mode,compression_level);
    }
  }
//-----------------------------------------------------------------------------
// PID - one PID block per track block
//-----------------------------------------------------------------------------
  if (fMakePid) {
    int nblocks = fTrackBlockName.size();

    for (int i=0; i<nblocks; i++) {
      const char* block_name  = fPidBlockName[i].data();
      if (block_name[0] == 0x0)                             continue;

      TStnDataBlock* db = AddDataBlock(block_name,
				       "TStnPidBlock",
				       StntupleInitMu2ePidBlock,
				       buffer_size,
				       split_mode,
				       compression_level);
      if (db) {
	db->AddCollName("mu2e::AvikPIDNewProductCollection",fPidCollTag[i].data());
	SetResolveLinksMethod(block_name,StntupleInitMu2ePidBlockLinks);
      }
    }
  }
//-----------------------------------------------------------------------------
// simulated particles 
//-----------------------------------------------------------------------------
  if (fMakeSimp) {
    fInitSimpBlock = new StntupleInitSimpBlock();

    fInitSimpBlock->SetSimpCollTag       (fSimpCollTag       );
    fInitSimpBlock->SetShCollTag         (fShCollTag         );
    fInitSimpBlock->SetSdmcCollTag       (fSdmcCollTag       );
    fInitSimpBlock->SetVDHitsCollTag     (fVdhCollTag        );
    fInitSimpBlock->SetPrimaryParticleTag(fPrimaryParticleTag);
    fInitSimpBlock->SetMinSimpMomentum   (fMinSimpMomentum   );
    fInitSimpBlock->SetMaxZ              (fSimpMaxZ          );
    fInitSimpBlock->SetGenProcessID      (fGenId.id()        );
    fInitSimpBlock->SetPdgID             (fPdgId             );
    fInitSimpBlock->SetMinNStrawHits     (fMinNStrawHits     );

    AddDataBlock("SimpBlock","TSimpBlock",fInitSimpBlock,buffer_size,split_mode,compression_level);
  }
//-----------------------------------------------------------------------------
// StepPointMC collections - could be several
//-----------------------------------------------------------------------------
  if (fMakeStepPointMC) {
    int nblocks = fSpmcBlockName.size();

    for (int i=0; i<nblocks; i++) {
      const char* block_name = fSpmcBlockName[i].data();

      StntupleInitStepPointMCBlock* init_block = new StntupleInitStepPointMCBlock();
      fInitStepPointMCBlock->Add(init_block);

      init_block->SetSpmcCollTag(fSpmcCollTag[i]);
      init_block->SetStatusG4Tag(fStatusG4Tag[i]);

      TStnDataBlock* db = AddDataBlock(block_name,
				       "TStepPointMCBlock",
				       init_block,
				       buffer_size,
				       split_mode,
				       compression_level);
      if (db) {
	((TStepPointMCBlock*) db)->SetGenProcessID(fGenId.id());
      }
    }
  }
//-----------------------------------------------------------------------------
// straw hit data
//-----------------------------------------------------------------------------
  if (fMakeStrawHits) {
    fInitStrawHitBlock = new stntuple::InitStrawHitBlock();

    fInitStrawHitBlock->SetShCollTag   (fShCollTag);
    fInitStrawHitBlock->SetStrawDigiCollTag  (fStrawDigiCollTag  );
    fInitStrawHitBlock->SetStrawDigiMCCollTag(fSdmcCollTag);

    if (fMakeStrawWaveforms) { 
      fInitStrawHitBlock->SetWriteSdwf(1);
      fInitStrawHitBlock->SetSdwfCollTag(fSdwfCollTag);
    }

    AddDataBlock("StrawHitBlock","TStrawHitBlock",fInitStrawHitBlock,buffer_size,split_mode,compression_level);
  }
//--------------------------------------------------------------------------------
// time clusters
//--------------------------------------------------------------------------------
  if (fMakeTimeClusters) {
    int nblocks = fTClBlockName.size();

    for (int i=0; i<nblocks; i++) {
      const char* block_name = fTClBlockName[i].data();
      if ((block_name[0] == 0) || (block_name[0] == ' ')) continue;

      StntupleInitTimeClusterBlock* init_block = new StntupleInitTimeClusterBlock();
      fInitTimeClusterBlock->Add(init_block);

      init_block->SetTimeClusterCollTag(fTClCollTag[i]);

      init_block->SetShCollTag   (fShCollTag);
      init_block->SetChCollTag   (fChCollTag);
      init_block->SetStrawDigiMCCollTag(fSdmcCollTag);

      AddDataBlock(block_name,"TStnTimeClusterBlock",init_block,buffer_size,split_mode,compression_level);
    }
  }
//--------------------------------------------------------------------------------
// trackSeed data
// pass to the track seed block the name of the single-straw combohit collection (fShCollTag)
//--------------------------------------------------------------------------------
  if (fMakeTrackSeeds) {
    int nb = fKsfBlockName.size();

    for (int i=0; i<nb; i++) {
      const char* block_name = fKsfBlockName[i].data();
      StntupleInitTrackSeedBlock* init_block = new StntupleInitTrackSeedBlock();
      fInitTrackSeedBlock->Add(init_block);

      init_block->SetHsBlockName(fHelixBlockName[i]);
      init_block->SetSschCollTag(fShCollTag       );
      init_block->SetKsfCollTag (fKsfCollTag[i]);
      init_block->SetSdmcCollTag(fSdmcCollTag);

      AddDataBlock(block_name,"TStnTrackSeedBlock",init_block,buffer_size,
                   split_mode,compression_level);
    }
  }
//-----------------------------------------------------------------------------
// track straw hits
//-----------------------------------------------------------------------------
  if (fMakeTrackStrawHits) {
    int nblocks = fTrackSHBlockName.size();

    for (int i=0; i<nblocks; i++) {
      const char* block_name = fTrackSHBlockName[i].data();

      stntuple::InitTrackStrawHitBlock* init_block = new stntuple::InitTrackStrawHitBlock();
      fInitTrackStrawHitBlock->Add(init_block);

      init_block->SetShCollTag    (fShCollTag   );
      init_block->SetKalSeedCollTag     (fTrackCollTag[i]   );  // tracks saved as lists of KalSeeds
      init_block->SetStrawDigiCollTag   (fStrawDigiCollTag  );
      init_block->SetStrawDigiMCCollTag (fSdmcCollTag);

      AddDataBlock(block_name,"TTrackStrawHitBlock",init_block,buffer_size,
		   split_mode,compression_level);
    }
  }
//-----------------------------------------------------------------------------
// track branches: for ROOT v3 to use streamers one has to specify split=-1
// KinKal fit links track directly to the helix, no more TrackSeeds
//-----------------------------------------------------------------------------
  if (fMakeTracks) {
    int nblocks = fTrackBlockName.size();

    StntupleInitTrackBlock* init_block(nullptr);
    for (int i=0; i<nblocks; i++) {
      const char* block_name = fTrackBlockName[i].data();

      if      (fTrackFitType == 1) init_block = new StntupleInitTrackBlock   ();
      else if (fTrackFitType == 2) init_block = new StntupleInitTrackBlock_KK();

      fInitTrackBlock->Add(init_block);

      // init_block->SetCaloClusterCollTag (fCaloClusterMaker);
      init_block->SetSsChCollTag        (fShCollTag       );
      init_block->SetKFFCollTag         (fTrackCollTag[i] );  // tracks saved as lists of KalSeeds
      init_block->SetPIDProductCollTag  (fPidCollTag[i]   );
      init_block->SetVdhCollTag         (fVdhCollTag      );  // 
      init_block->SetStrawDigiMCCollTag (fSdmcCollTag);
      init_block->SetTciCollTag         (fTciCollTag[i]);
      init_block->SetTcmCollTag         (fTcmCollTag[i]);
      init_block->SetTrkQualCollTag     (fTrkQualCollTag[i]);

      init_block->SetDoubletAmbigResolver(fDar);

      TStnDataBlock* db = AddDataBlock(block_name,"TStnTrackBlock",init_block,
				       buffer_size,split_mode,compression_level);
//-----------------------------------------------------------------------------
// each track points back to its seed
// if nshortblocks != 0, for each track we store an index to that tracks's seed
//-----------------------------------------------------------------------------
      if (db) {
        if (fTrackFitType == 1) {
//-----------------------------------------------------------------------------
// BTRK fits
//-----------------------------------------------------------------------------
          if (fTrackTsBlockName.size() > 0) {
            init_block->SetTrackTsBlockName(fTrackTsBlockName[i].data());
            init_block->SetTrackTsCollTag  (fTrackTsCollTag  [i]);
          }
        }
        else if (fTrackFitType == 2) {
//-----------------------------------------------------------------------------
// KinKal fits
//-----------------------------------------------------------------------------
          init_block->SetTrackHsBlockName(fTrackHsBlockName[i].data());
        }
      }
    }
  }
//-----------------------------------------------------------------------------
// CRV
//-----------------------------------------------------------------------------
  if (fMakeTrigger) {
    fInitTriggerBlock = new StntupleInitTriggerBlock();
    fInitTriggerBlock->SetTriggerResultsTag(fTriggerResultsTag);
    fInitTriggerBlock->SetNTriggerBits     (fNTriggerBits);

    AddDataBlock("TriggerBlock","TStnTriggerBlock",
		 fInitTriggerBlock,
		 buffer_size,
		 split_mode,
		 compression_level);
  }

  THistModule::afterEndJob();
}

//_____________________________________________________________________________
void StntupleMaker::analyze(const AbsEvent& AnEvent) {

  // when execution comes here al the registered data blocks are already
  // initialized with the event data. Left: variables in the data blocks
  // which depend on the variable defined in other blocks, like track number
  // for a muon or an electron - the idea is that these are defined during the
  // 2nd loop in FillStntupleModule, where ResolveLinks methods are called
  // for each data block

//-----------------------------------------------------------------------------
// connect to the error reporting facility
//-----------------------------------------------------------------------------
//  TStnErrorLogger* logger = Event()->GetErrorLogger();
//   logger->Connect("Report(Int_t, const char*)",
// 		  "StntupleModule",
// 		  this,
// 		  "LogError(const char*)");
//-----------------------------------------------------------------------------
// disconnect from the error reporting signal and return back to AC++
//-----------------------------------------------------------------------------
//   logger->Disconnect("Report(Int_t,const char*)",
// 		     this,"LogError(Int_t,const char*)");

  // bool passed(1);

  // if (fMinNHelices > 0) {
  //   auto hH = AnEvent.getValidHandle<mu2e::HelixSeedCollection>(fCutHelixSeedCollTag);
  //   int  nh = hH->size();
  //   passed  = (nh >= fMinNHelices);
  // }
}


//------------------------------------------------------------------------------
void StntupleMaker::GetDefTrackCollName(char* Name) {
  // put in a working kludge first

  strcpy(Name,"default");
}


// //_____________________________________________________________________________
// int StntupleMaker::InitCalDataBlock(TStnDataBlock* Block) {
//   int mode = 0;
//   AbsEvent* event = AbsEnv::instance()->theEvent();
//   return StntupleInitMu2eCalDataBlock(Block,event,mode);
// }

// //_____________________________________________________________________________
// int StntupleMaker::InitHeaderBlock(TStnDataBlock* Block) {
//   int mode = 0;
//   AbsEvent* event = AbsEnv::instance()->theEvent();
//   return StntupleInitMu2eHeaderBlock(Block,event,mode);
// }

// //_____________________________________________________________________________
// int StntupleMaker::InitTrackBlock(TStnDataBlock* Block) {
//   int mode = 0;
//   AbsEvent* event = AbsEnv::instance()->theEvent();
//   return StntupleInitMu2eTrackBlock(Block,event,mode);
// }

//_____________________________________________________________________________
// int StntupleMaker::InitTriggerBlock(TStnDataBlock* Block) {
//   int mode = 0;
//   AbsEvent* event = AbsEnv::instance()->theEvent();
//   return StntupleInitMu2eTriggerBlock(Block,event,mode);
// }

} // end namespace mu2e

using mu2e::StntupleMaker;

DEFINE_ART_MODULE(StntupleMaker)
