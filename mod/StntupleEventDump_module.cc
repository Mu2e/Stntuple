///////////////////////////////////////////////////////////////////////////////
//
// $Author: murat $
// $Date: 2014/09/20 17:54:06 $
//
// Contact person:  PM
//
///////////////////////////////////////////////////////////////////////////////

// C++ includes.
#include <map>
#include <iostream>
#include <string>

// Framework includes.
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Principal/Selector.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// Mu2e includes.
#include "canvas/Persistency/Provenance/ProductID.h"
#include "art/Framework/Principal/Provenance.h"
#include "Offline/RecoDataProducts/inc/KalSeed.hh"

// ROOT includes
#include "TApplication.h"
#include "TDirectory.h"

#include "Stntuple/print/TAnaDump.hh"
#include "Stntuple/mod/TModule.hh"

using namespace std;
using CLHEP::Hep3Vector;

namespace mu2e {

  class StntupleEventDump : public TModule {
  struct Config {
    using Name = fhicl::Name;
    using Comment = fhicl::Comment;
    fhicl::Atom<int> printKalSeeds   {Name("printKalSeeds"   ), Comment("print KalSeeds")};
    fhicl::Atom<int> printKalSeedHits{Name("printKalSeedHits"), Comment("print KalSeed hits")};
    fhicl::DelegatedParameter TModule{Name("TModule"), Comment("TModule parameters")};
  };
  using Parameters = art::EDAnalyzer::Table<Config>;
  private:
//-----------------------------------------------------------------------------
// Input parameters: Module labels 
//-----------------------------------------------------------------------------
    int                _printKalSeeds;        // print KalSeed collections
    int                _printKalSeedHits;     // if >0, also print hits
//-----------------------------------------------------------------------------
// end of input parameters
// Options to control the display
// hit flag bits which should be ON and OFF
//-----------------------------------------------------------------------------
    TApplication*                               fApplication;

  public:

    explicit StntupleEventDump(const Parameters& config);
    virtual ~StntupleEventDump();

    int      getData(const art::Event* Evt);
    void     Init   (art::Event* Evt);
//-----------------------------------------------------------------------------
// overloaded virtual methods of the base class
//-----------------------------------------------------------------------------
    virtual void     beginJob();
    virtual void     beginRun(const art::Run& aRun);
    virtual void     analyze (const art::Event& Evt);
  };


//-----------------------------------------------------------------------------
  StntupleEventDump::StntupleEventDump(const Parameters& config) :
    TModule(config.get_PSet(), config().TModule.get<fhicl::ParameterSet>(), "StntupleEventDump"),
    _printKalSeeds(config().printKalSeeds()),
    _printKalSeedHits(config().printKalSeedHits())
  {
    fApplication      = 0;
  }

//-----------------------------------------------------------------------------
  StntupleEventDump::~StntupleEventDump() {
    if (fApplication) delete fApplication;
  }

//-----------------------------------------------------------------------------
  void StntupleEventDump::beginJob() {
    //     const char oname[] = "StntupleEventDump::beginJob";
    int    tmp_argc(0);
    char** tmp_argv(0);

    if (!gApplication) {
      fApplication = new TApplication("StntupleEventDump_module", &tmp_argc, tmp_argv);
    }
  }

//-----------------------------------------------------------------------------
  void StntupleEventDump::beginRun(const art::Run& Run) {
  }

//-----------------------------------------------------------------------------
// get data from the event record
//-----------------------------------------------------------------------------
  int StntupleEventDump::getData(const art::Event* Evt) {
    //    const char* oname = "StntupleEventDump::getData";
    return 0;
  }

//-----------------------------------------------------------------------------
  void StntupleEventDump::analyze(const art::Event& Evt) {
    //    const char* oname = "StntupleEventDump::filter";
//-----------------------------------------------------------------------------
// go into interactive mode, till '.q' is pressed
//-----------------------------------------------------------------------------
    TAnaDump* ad = TAnaDump::Instance();
    ad->SetEvent(&Evt);

//   TModule::analyze(Evt);
    if (_printKalSeeds) {
//-----------------------------------------------------------------------------
// print all KalSeed collections (any instance name) found in the event;
// TAnaDump::printKalSeedCollection expects tag = label:instance:process
//-----------------------------------------------------------------------------
      auto vcoll = Evt.getMany<mu2e::KalSeedCollection>();

      for (auto const& handle : vcoll) {
        if (! handle.isValid()) continue;
        const art::Provenance* prov = handle.provenance();

        printf("moduleLabel: %-20s, productInstanceName: %-20s, processName:= %-30s\n",
               prov->moduleLabel().data(),
               prov->productInstanceName().data(),
               prov->processName().data());

        std::string tag = prov->moduleLabel()+":"+prov->productInstanceName()+":"+prov->processName();
        ad->printKalSeedCollection(tag.data(),_printKalSeedHits);
      }
    }
  }
}

using mu2e::StntupleEventDump;
DEFINE_ART_MODULE(StntupleEventDump)
