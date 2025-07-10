// ======================================================================
//
// ======================================================================
#ifndef __daqana_mod_StntupleRecoAna_hh__
#define __daqana_mod_StntupleRecoAna_hh__

// ROOT includes
#include "TH1F.h"
//#include "TFolder.h"

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileDirectory.h"
#include "art_root_io/TFileService.h"
#include "fhiclcpp/ParameterSet.h"

// #ifndef __CLING__ 
// #include "artdaq-core-mu2e/Overlays/FragmentType.hh"

// typedef artdaq::Fragment::type_t  type_t;

#include "artdaq-core-mu2e/Overlays/Decoders/TrackerDataDecoder.hh"
#include "artdaq-core/Data/Fragment.hh"
//  #else 
//  namespace mu2e {
//    class TrackerFragment;
//    class TrackerFragment::TrackerDataPacket;
//  }

// namespace artdaq {
//   class Fragment;
// }
// #endif

// Mu2e includes
#include "Offline/DataProducts/inc/StrawId.hh"
#include "Offline/DataProducts/inc/TrkTypes.hh"
#include "Offline/RecoDataProducts/inc/StrawHit.hh"
#include "Offline/RecoDataProducts/inc/ComboHit.hh"
#include "Offline/RecoDataProducts/inc/KalSeed.hh"
#include "Offline/RecoDataProducts/inc/HelixSeed.hh"

// #include "Stntuple/print/TAnaDump.hh"
#include "Stntuple/mod/mod/THistModule.hh"

#include <iostream>
#include <memory>
#include <string>

// ======================================================================
namespace mu2e {

  class StntupleRecoAna : public THistModule {
  public:

    struct Config {
      using Name    = fhicl::Name;
      using Comment = fhicl::Comment;
      fhicl::Atom<art::InputTag>   sschCollTag            {Name("sschCollTag"       )    , Comment("SS ComboHit collection tag") };
      fhicl::Atom<art::InputTag>   chCollTag              {Name("chCollTag"         )    , Comment("ComboHit collection tag"   ) };
      fhicl::Atom<art::InputTag>   tcCollTag              {Name("tcCollTag"         )    , Comment("TimeCluster collection tag") };
      fhicl::Atom<art::InputTag>   sdmcCollTag            {Name("sdmcCollTag"       )    , Comment("StrawDigiMC collection tag") };
      fhicl::Atom<int>             debugLevel             {Name("debugLevel"        )    , Comment("debug level"               ) };
    };

    enum {
      kMaxEvHistSets = 10,
      kMaxShHistSets = 10,
      kMaxChHistSets = 10,
      kMaxTcHistSets = 10
    };

    struct ShHist_t {
      TH1F*         time;
    };

    struct ChHist_t {
      TH1F*         time;
    };

    struct TcHist_t {
      TH1F*         nch;
      TH1F*         nsh;
      TH1F*         time;
    };
                                        // event histograms
    struct EvHist_t {
      TH1F*         nsh;
      TH1F*         nch;
      TH1F*         ntc;
    };

    struct EventData_t {
      int                          event_number;
      const ComboHitCollection*    sschColl;
      int                          nsh;
      const ComboHitCollection*    chColl;
      int                          nch;
      const TimeClusterCollection* tcColl;
      int                          ntc;
    } _data;
//-----------------------------------------------------------------------------
// data part
//
// talk-to parameters
//-----------------------------------------------------------------------------
    art::InputTag    _sschCollTag;
    art::InputTag    _chCollTag;
    art::InputTag    _tcCollTag;
    art::InputTag    _sdmcCollTag;

    int              _debugLevel;
//-----------------------------------------------------------------------------
// the rest
//-----------------------------------------------------------------------------
    const art::Event*  _event;
//-----------------------------------------------------------------------------
// forgetting, for now, about multiple DTC's
//-----------------------------------------------------------------------------
    struct Hist_t {
      EvHist_t ev[kMaxEvHistSets];
      ShHist_t sh[kMaxShHistSets];
      ChHist_t ch[kMaxChHistSets];
      TcHist_t tc[kMaxTcHistSets];
    } _Hist;

    explicit StntupleRecoAna(fhicl::ParameterSet const& pset);
    // explicit StntupleRecoAna(const art::EDAnalyzer::Table<Config>& config);
    virtual ~StntupleRecoAna() {}
    
    virtual void beginRun(const art::Run& cr_Run) override;

    virtual void beginJob() override;
    virtual void endJob  () override;
    
    virtual void analyze         (const art::Event& cr_Event) override;

    void         book_tc_histograms(art::TFileDirectory* Dir, int RunNumber, TcHist_t*   Hist);
    void         book_ch_histograms(art::TFileDirectory* Dir, int RunNumber, ChHist_t*   Hist);
    void         book_sh_histograms(art::TFileDirectory* Dir, int RunNumber, ShHist_t*   Hist);
    void         book_ev_histograms(art::TFileDirectory* Dir, int RunNumber, EvHist_t*   Hist);
    void         book_histograms             (int RunNumber);
  
    int          fill_tc_histograms(TcHist_t* Hist, const TimeCluster* Tc  );
    int          fill_ch_histograms(ChHist_t* Hist, const ComboHit*    Ch  );
    int          fill_sh_histograms(ShHist_t* Hist, const ComboHit*    Sh  );
    int          fill_ev_histograms(EvHist_t* Hist, EventData_t*       Data);

                                        // returns -1 if in trouble
    int          fill_histograms();
    bool         find_data(const art::Event& cr_Event);
  };
}
#endif
