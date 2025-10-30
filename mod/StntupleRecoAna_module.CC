///////////////////////////////////////////////////////////////////////////////
// this module makes key histograms which shoudl allow to diagnose 
// the reco failures, 
// so far, those are mostly due to misconfiguration of the reconstruction jobs
///////////////////////////////////////////////////////////////////////////////
#include "TRACE/tracemf.h"
#define TRACE_NAME "StntupleRecoAna_module"

#include "Stntuple/mod/StntupleRecoAna_module.hh"

namespace mu2e {
// ======================================================================

  StntupleRecoAna::StntupleRecoAna(fhicl::ParameterSet const& PSet) : 
    THistModule(PSet,"StntupleRecoAna"  ),

    // _sschCollTag           (config().sschCollTag()       ),
    // _chCollTag             (config().chCollTag()         ),
    // _sdmcCollTag           (config().sdmcCollTag()       ),
    // _writeStrawHitFlags    (config().writeStrawHitFlags()),

    _sschCollTag     (PSet.get<art::InputTag>   ("sschCollTag"      )),
    _chCollTag       (PSet.get<art::InputTag>   ("chCollTag"        )),
    _tcCollTag       (PSet.get<art::InputTag>   ("tcCollTag"        )),
    _sdmcCollTag     (PSet.get<art::InputTag>   ("sdmcCollTag"      )),
    _debugLevel      (PSet.get<int>             ("debugLevel"       ))
  {}


//-----------------------------------------------------------------------------
  void StntupleRecoAna::book_tc_histograms(art::TFileDirectory* Dir, int RunNumber, TcHist_t* Hist) {
    Hist->nch           = Dir->make<TH1F>("nch"  , Form("run %06i: N(ch)" ,RunNumber), 200, 0,  200);
    Hist->nsh           = Dir->make<TH1F>("nsh"  , Form("run %06i: N(sh)" ,RunNumber), 200, 0,  200);
    Hist->time          = Dir->make<TH1F>("time" , Form("run %06i: time"  ,RunNumber), 400, 0, 2000);
  }

//-----------------------------------------------------------------------------
  void StntupleRecoAna::book_ch_histograms(art::TFileDirectory* Dir, int RunNumber, ChHist_t* Hist) {
    Hist->time          = Dir->make<TH1F>("time" , Form("run %06i: CH time" ,RunNumber), 400, 0, 2000);
  }

//-----------------------------------------------------------------------------
  void StntupleRecoAna::book_sh_histograms(art::TFileDirectory* Dir, int RunNumber, ShHist_t* Hist) {
    Hist->time          = Dir->make<TH1F>("time" , Form("run %06i: SH time" ,RunNumber), 400, 0, 2000);
  }

//-----------------------------------------------------------------------------
  void StntupleRecoAna::book_ev_histograms(art::TFileDirectory* Dir, int RunNumber, EvHist_t* Hist) {
    Hist->nsh           = Dir->make<TH1F>("nsh" , Form("run %06i: nsh" ,RunNumber), 200, -0.5, 199.5);
    Hist->nch           = Dir->make<TH1F>("nch" , Form("run %06i: nch" ,RunNumber), 200, -0.5, 199.5);
    Hist->ntc           = Dir->make<TH1F>("ntc" , Form("run %06i: ntc" ,RunNumber), 200, -0.5, 199.5);
  }

//-----------------------------------------------------------------------------
  void StntupleRecoAna::book_histograms(int RunNumber) {
    art::ServiceHandle<art::TFileService> tfs;
    
    art::TFileDirectory ev_dir = tfs->mkdir("ev_00");
    book_ev_histograms(&ev_dir,RunNumber,&_Hist.ev[0]);

    art::TFileDirectory sh_dir = tfs->mkdir(Form("sh_00"));
    book_sh_histograms(&sh_dir,RunNumber,&_Hist.sh[0]);

    art::TFileDirectory ch_dir = tfs->mkdir(Form("ch_00"));
    book_ch_histograms(&ch_dir,RunNumber,&_Hist.ch[0]);

    art::TFileDirectory tc_dir = tfs->mkdir(Form("tc_00"));
    book_tc_histograms(&tc_dir,RunNumber,&_Hist.tc[0]);

    printf("[mu2e::StntupleRecoAna] pointer to the module: 0x%8p\n",(void*) this);
  }

//-----------------------------------------------------------------------------
// 
//-----------------------------------------------------------------------------
void StntupleRecoAna::beginRun(const art::Run& aRun) {
  int rn  = aRun.run();

  // if (_initialized != 0) return;
  // _initialized = 1;

  book_histograms(rn);
}

//--------------------------------------------------------------------------------
  void StntupleRecoAna::beginJob() {
  }

//-----------------------------------------------------------------------------
  void StntupleRecoAna::endJob() {
    printf("[mu2e::StntupleRecoAna] pointer to the module: 0x%8p\n",(void*) this);
  }


//-----------------------------------------------------------------------------
  int StntupleRecoAna::fill_ch_histograms(ChHist_t* Hist, const ComboHit* Ch) {
    Hist->time->Fill(Ch->correctedTime());
    return 0;
  }

//-----------------------------------------------------------------------------
  int StntupleRecoAna::fill_ev_histograms(EvHist_t* Hist,  EventData_t* Data) {
    Hist->nch->Fill(Data->nch);
    Hist->nsh->Fill(Data->nsh);
    Hist->ntc->Fill(Data->ntc);

    return 0;
  }

//-----------------------------------------------------------------------------
  int StntupleRecoAna::fill_sh_histograms(ShHist_t* Hist, const ComboHit* Ssch) {
    Hist->time->Fill(Ssch->correctedTime());
    return 0;
  }

//-----------------------------------------------------------------------------
  int StntupleRecoAna::fill_tc_histograms(TcHist_t* Hist, const TimeCluster* Tc) {
    Hist->time->Fill(Tc->t0().t0());
    Hist->nch->Fill(Tc->nhits());
    Hist->nsh->Fill(Tc->nStrawHits());
    return 0;
  }

//-----------------------------------------------------------------------------
// fill_roc_histograms also fills the channel histograms
// if in error, only histogram the error code
//-----------------------------------------------------------------------------
  int StntupleRecoAna::fill_histograms() {
//-----------------------------------------------------------------------------
// event histograms
//-----------------------------------------------------------------------------
    fill_ev_histograms(&_Hist.ev[0],&_data);
//-----------------------------------------------------------------------------
// straw hit histograms
//-----------------------------------------------------------------------------
    for (int i=0; i<_data.nsh; i++) {
      const ComboHit* ssch = &(*_data.sschColl )[i];
      fill_sh_histograms(&_Hist.sh[0],ssch);
    }
//----------------------------------------------------------------------------
// combo hit histograms
//-----------------------------------------------------------------------------
    for (int i=0; i<_data.nch; i++) {
      const ComboHit* ch = &(*_data.chColl )[i];
      fill_ch_histograms(&_Hist.ch[0],ch);
    }
//-----------------------------------------------------------------------------
// time cluster histograms
//-----------------------------------------------------------------------------
    for (int i=0; i<_data.ntc; i++) {
      const TimeCluster *tc = &(*_data.tcColl )[i];
      fill_tc_histograms(&_Hist.tc[0],tc);
    }

    return 0;
  }

//-----------------------------------------------------------------------------
  bool StntupleRecoAna::find_data(const art::Event& cr_Event) {

    auto chcH      = cr_Event.getValidHandle<mu2e::ComboHitCollection>(_chCollTag);
    _data.chColl   = chcH.product();
    
    auto sschcH    = cr_Event.getValidHandle<mu2e::ComboHitCollection>(_sschCollTag);
    _data.sschColl = sschcH.product();
    
    auto tccH      = cr_Event.getValidHandle<mu2e::TimeClusterCollection>(_tcCollTag);
    _data.tcColl   = tccH.product();
    
    return (_data.chColl != nullptr) and (_data.sschColl != nullptr);
  }

//--------------------------------------------------------------------------------
// assume that we only have tracker fragment(s)
//-----------------------------------------------------------------------------
void StntupleRecoAna::analyze(const art::Event& cr_Event) {

  _event = &cr_Event;

  TLOG(TLVL_DEBUG+10) << "event=" << _event->run() << ":" << _event->subRun() << ":" << _event->event() ;

  find_data(cr_Event);
//-----------------------------------------------------------------------------
// figure out event parameters
//-----------------------------------------------------------------------------
  _data.nsh = _data.sschColl->size();
  _data.nch = _data.chColl->size();
  _data.ntc = _data.tcColl->size();

  fill_histograms();
//-----------------------------------------------------------------------------
// finally, if requested, go into interactive mode, 
// fInteractiveMode = 0 : do not stop
// fInteractiveMode = 1 : stop after each event (event display mode)
// fInteractiveMode = 2 : stop only in the end of run, till '.q' is pressed
//-----------------------------------------------------------------------------
  TModule::analyze(cr_Event);
}

} // end namespace mu2e

DEFINE_ART_MODULE(mu2e::StntupleRecoAna)
