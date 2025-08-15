///////////////////////////////////////////////////////////////////////////////
// clang-format off 
// interactive 2D event display. 
//
// $Id: MuHitDisplay_module.cc,v 1.6 2014/09/20 17:54:06 murat Exp $
// $Author: murat $
// $Date: 2014/09/20 17:54:06 $
//
// Contact person:  Pavel Murat, Gianantonio Pezzulo
///////////////////////////////////////////////////////////////////////////////
#ifndef Stntuple_mod_MuHitDisplay
#define Stntuple_mod_MuHitDisplay

#include <iostream>
#include <string>

#include "Stntuple/mod/THistModule.hh"

#include "art/Framework/Principal/Event.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "Offline/TrackerConditions/inc/StrawResponse.hh"

#include "Offline/TrackerGeom/inc/Tracker.hh"
#include "Offline/CalorimeterGeom/inc/DiskCalorimeter.hh"
#include "Offline/CalorimeterGeom/inc/Calorimeter.hh"

#include "Offline/MCDataProducts/inc/GenParticle.hh"
#include "Offline/MCDataProducts/inc/SimParticle.hh"
#include "Offline/MCDataProducts/inc/StepPointMC.hh"
#include "Offline/MCDataProducts/inc/StrawDigiMC.hh"

#include "Offline/RecoDataProducts/inc/CaloHit.hh"
#include "Offline/RecoDataProducts/inc/CaloCluster.hh"
#include "Offline/RecoDataProducts/inc/ComboHit.hh"
#include "Offline/RecoDataProducts/inc/StrawHitFlag.hh"
#include "Offline/RecoDataProducts/inc/StrawDigi.hh"

#include "Offline/RecoDataProducts/inc/CrvRecoPulse.hh"
#include "Offline/RecoDataProducts/inc/TimeCluster.hh"

#include "Offline/BTrkData/inc/TrkStrawHit.hh"

class TApplication;

#include "Stntuple/gui/TStnVisManager.hh"
#include "Stntuple/gui/TStnGeoManager.hh"

#include "Stntuple/alg/TStnTrackID.hh"

#include "Stntuple/obj/TStnHeaderBlock.hh"
#include "Stntuple/obj/TSimpBlock.hh"

#ifndef __CLING__
#include "Offline/Mu2eUtilities/inc/McUtilsToolBase.hh"
#endif
// using namespace std;

class DoubletAmbigResolver;

namespace mu2e {

class TrackerPanelMap;

class MuHitDisplay : public THistModule {
public:
#ifndef __CLING__  
  struct VmConfig {
    using Name    = fhicl::Name;
    using Comment = fhicl::Comment;
    fhicl::Atom<int>         displayStrawDigiMC {Name("displayStrawDigiMC" )    , Comment("display sdmc") };
    fhicl::Atom<int>         displayStrawHitsXY {Name("displayStrawHitsXY" )    , Comment("display straw hits in XY") };
    fhicl::Atom<int>         debugLevel         {Name("debugLevel"         )    , Comment("debug level" ) };
    fhicl::Atom<float>       ewLength           {Name("ewLength"           )    , Comment("ewLength"    ) };
    fhicl::Atom<float>       tMin               {Name("tMin"               )    , Comment("tMin"        ) };
    fhicl::Atom<float>       tMax               {Name("tMax"               )    , Comment("tMax"        ) };
    fhicl::Atom<float>       minEDep            {Name("minEDep"            )    , Comment("minEDep"     ) };
    fhicl::Atom<float>       maxEDep            {Name("maxEDep"            )    , Comment("maxEDep"     ) };
    fhicl::Atom<std::string> defaultView        {Name("defaultView"        )    , Comment("defaultView" ) };
    fhicl::Sequence<int>     visibleStations    {Name("visibleStations"    )    , Comment("list of visible stations" ) };
    fhicl::Atom<int>         evdPanelDebugLevel {Name("evdPanelDebugLevel" )    , Comment("TEvdPanel::fgDebugLevel"  ) };
    fhicl::Atom<int>         stnViewDebugLevel  {Name("stnViewDebugLevel"  )    , Comment("TStnView::fgDebugLevel"   ) };
  };
                                        // Geometry Manager configuration
  struct GmConfig {
    using Name    = fhicl::Name;
    using Comment = fhicl::Comment;
    fhicl::Atom<float>       bField             {Name("bField"             )    , Comment("bField"      ) };
  };
    
  struct Config {
    using Name    = fhicl::Name;
    using Comment = fhicl::Comment;
    fhicl::DelegatedParameter    THistModule            {Name("THistModule"       )    , Comment("THistModule parameters") };
    fhicl::Atom<art::InputTag>   genpCollTag            {Name("genpCollTag"       )    , Comment("GenParticle collection tag") };
    fhicl::Atom<art::InputTag>   spmcCollTag            {Name("spmcCollTag"       )    , Comment("SPMC  collection tag") };
    fhicl::Atom<art::InputTag>   caloClusterCollTag     {Name("caloClusterCollTag" )   , Comment("calo cluster collection tag") };
    fhicl::Atom<art::InputTag>   crvRecoPulsesCollTag   {Name("crvRecoPulsesCollTag")  , Comment("CRV reco pulses collection tag") };
      
    fhicl::Atom<art::InputTag>   shCollTag              {Name("strawHitCollTag"    )   , Comment("StrawHit collection tag"   ) };
    fhicl::Atom<art::InputTag>   chCollTag              {Name("comboHitCollTag"    )   , Comment("ComboHit collection tag"   ) };

    fhicl::Atom<art::InputTag>   sdCollTag              {Name("sdCollTag"         )    , Comment("StrawDigi collection tag"   ) };
    fhicl::Atom<art::InputTag>   sdmcCollTag            {Name("sdmcCollTag"       )    , Comment("StrawDigiMC collection tag") };
    fhicl::Atom<art::InputTag>   swCollTag              {Name("swCollTag"         )    , Comment("Straw WF collection tag"   ) };

    fhicl::Atom<art::InputTag>   helixSeedCollTag       {Name("helixSeedCollTag"  )    , Comment("helix seed collection tag"   ) };
    fhicl::Atom<art::InputTag>   ksfCollTag             {Name("ksfCollTag"        )    , Comment("KSF   collection tag"   ) };
    fhicl::Atom<art::InputTag>   kffCollTag             {Name("kffCollTag"        )    , Comment("KFF   collection tag"   ) };
    fhicl::Atom<art::InputTag>   ctsCollTag             {Name("ctsCollTag"        )    , Comment("Cosmic Track Seed coll tag"   ) };
    fhicl::Atom<art::InputTag>   simpCollTag            {Name("simpCollTag"       )    , Comment("SIMP  collection tag"   ) };
    fhicl::Atom<art::InputTag>   timeClusterCollTag     {Name("timeClusterCollTag")    , Comment("time cluster collection tag") };
    fhicl::Atom<art::InputTag>   phiClusterCollTag      {Name("phiClusterCollTag" )    , Comment("phi  cluster collection tag") };
    fhicl::Atom<art::InputTag>   caloHitCollTag         {Name("caloHitCollTag"    )    , Comment("calo hit collection tag") };
    fhicl::Atom<art::InputTag>   primaryParticleTag     {Name("primaryParticleTag")    , Comment("primary particle tag") };
    fhicl::Atom<art::InputTag>   vdHitsCollTag          {Name("vdHitsCollTag"     )    , Comment("VD hits collection tag") };
    fhicl::Atom<mu2e::GenId>     generatorID            {Name("generatorID"       )    , Comment("generator ID") , mu2e::GenId::CeEndpoint};
    fhicl::Atom<float>           minSimpMomentum        {Name("minSimpMomentum"   )    , Comment("minSimpMomentum") };
    fhicl::Atom<float>           maxSimpMomentum        {Name("maxSimpMomentum"   )    , Comment("maxSimpMomentum") };
    fhicl::Atom<bool>            showCRVOnly            {Name("showCRVOnly"       )    , Comment("showCRVOnly"    ) };
    fhicl::Atom<bool>            showTracks             {Name("showTracks"        )    , Comment("showTracks"     ) };
    fhicl::Table<VmConfig>       visManager             {Name("visManager"        )    , Comment("vis manager config" ) };
    fhicl::Table<GmConfig>       geoManager             {Name("geoManager"        )    , Comment("geo manager config" ) };
  };
#endif
private:
//-----------------------------------------------------------------------------
// Input parameters: Module labels 
//-----------------------------------------------------------------------------
  art::InputTag _genpCollTag;
  art::InputTag _spmcCollTag;
  art::InputTag _caloClusterCollTag;
  art::InputTag _crvRecoPulseCollTag;

  art::InputTag _shCollTag;
  art::InputTag _comboHitCollTag;
  art::InputTag _sdCollTag;		// strawDigiCollTag
  art::InputTag _sdmcCollTag;
  art::InputTag _swCollTag;
  
  art::InputTag _helixSeedCollTag;	// helix seed collection tag
  art::InputTag _ksfCollTag;
  art::InputTag _kffCollTag;
  art::InputTag _ctsCollTag;            // cosmic track seed coll tag
  art::InputTag _trackCollTag;
  art::InputTag _simpCollTag;
  art::InputTag _timeClusterCollTag;
  art::InputTag _phiClusterCollTag;
  art::InputTag _caloHitCollTag;
  art::InputTag _ppTag;			// primary particle tag
  art::InputTag _vdHitsCollTag;
  std::string _defaultView; // view open in the first window
  std::vector<int> _listOfVisibleStations;
  
  GenId         _generatorID;

  // StrawHitFlag  fGoodHitMask;           // unused
  // StrawHitFlag  fBadHitMask;            // unused
  // size_t        _minHits;               // unused
  double        _minSimpMomentum;       // to initialize SimpBlock and TStnVisManager
  double        _maxSimpMomentum;       // to initialize SimpBlock and TStnVisManager
//-----------------------------------------------------------------------------
// Control for CRV-specific viewing
//-----------------------------------------------------------------------------
  bool				_showCRVOnly;
  bool				_showTracks;

  //  fhicl::ParameterSet           _vmConfig;
#ifndef __CLING__
  VmConfig           _vmConfig;
  GmConfig           _gmConfig;
#endif
//-----------------------------------------------------------------------------
// end of input parameters
// Options to control the display
// hit flag bits which should be ON and OFF
//-----------------------------------------------------------------------------
  TApplication*                               fApplication;

  const mu2e::GenParticleCollection*          _genpColl;        // 

  const mu2e::StrawDigiCollection*            _sdColl;   //
  const mu2e::StrawDigiMCCollection*          _sdmcColl; //
  const mu2e::StrawDigiADCWaveformCollection* _swColl;   //
  const mu2e::StrawHitFlagCollection*         _shfColl;  //
  
  const mu2e::CaloHitCollection*              _caloHitColl;      //
  const mu2e::CaloClusterCollection*          _caloClusterColl;  //

  const mu2e::StepPointMCCollection*          _spmcColl;         // on virtual detectors
  const mu2e::SimParticleCollection*          _simpColl;         //

  mu2e::CrvRecoPulseCollection*		      fCrvPulseColl_Right;
  mu2e::CrvRecoPulseCollection*		      fCrvPulseColl_Left;
  mu2e::CrvRecoPulseCollection*		      fCrvPulseColl_TopDS;
  mu2e::CrvRecoPulseCollection*		      fCrvPulseColl_TopTS;
  mu2e::CrvRecoPulseCollection*		      fCrvPulseColl_Dwnstrm;
  mu2e::CrvRecoPulseCollection*		      fCrvPulseColl_Upstrm;
  
#ifndef __CLING__
  std::unique_ptr<McUtilsToolBase>            _mcUtils;
#endif
//-----------------------------------------------------------------------------
// need to detect the first call to initialize fVisManager
//-----------------------------------------------------------------------------
  int                                         fFirstCall;    
  int                                         fLastRun;
  TStnVisManager*                             fVisManager;
//-----------------------------------------------------------------------------
// geometry manager may need to be reinitialized at  run boundary
//-----------------------------------------------------------------------------
  TStnGeoManager*                             fGeoManager;
//-----------------------------------------------------------------------------
// reuse STNTUPLE data blocks
//-----------------------------------------------------------------------------
  TStnHeaderBlock*       fHeaderBlock;
  TSimpBlock*            fSimpBlock;

  TStnTrackID*           fTrackID;

  const Tracker*         fTracker;    // straw tracker geometry

  // DoubletAmbigResolver*  fDar;
//-----------------------------------------------------------------------------
// run-dependent things
//-----------------------------------------------------------------------------
  const mu2e::TrackerPanelMap*  fTrkPanelMap;

public:
  // for some reason, this line is required by art to allow the command line help print
  MuHitDisplay();
  
#ifndef __CLING__
  typedef art::EDAnalyzer::Table<Config> Parameters;
  explicit MuHitDisplay(const art::EDAnalyzer::Table<Config>& config);
#endif
  // explicit MuHitDisplay(fhicl::ParameterSet const& pset);
  virtual ~MuHitDisplay();

  int      getData      (const art::Event* Evt);
  int      getCRVSection(int   shieldNumber   );

  const mu2e::GenParticleCollection*            GetGenppColl      () { return _genpColl       ; }
  const mu2e::SimParticleCollection*            GetSimpColl       () { return _simpColl       ; }
  const mu2e::StepPointMCCollection*            GetSpmcColl       () { return _spmcColl       ; }
  const mu2e::StrawDigiCollection*              GetSdColl         () { return _sdColl         ; }
  const mu2e::StrawDigiADCWaveformCollection*   GetSwColl         () { return _swColl         ; }
  const mu2e::StrawDigiMCCollection*            GetSdmcColl       () { return _sdmcColl       ; }

  TSimpBlock*                                   GetSimpBlock      () { return fSimpBlock      ; }

  void                                          InitVisManager();
  void                                          InitGeoManager();
//-----------------------------------------------------------------------------
// overloaded virtual methods of the base class - EDAnalyzer 
//-----------------------------------------------------------------------------
  virtual void     beginJob();
  virtual void     beginRun(const art::Run& aRun);
  virtual void     endRun  (const art::Run& aRun);
  virtual void     analyze (const art::Event& Evt);

  //  ClassDef(mu2e::MuHitDisplay,0)
};
}
#endif
