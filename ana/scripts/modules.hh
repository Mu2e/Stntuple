#ifndef Stntuple_ana_scripts_modules_hh
#define Stntuple_ana_scripts_modules_hh

// #include "Stntuple/ana/TDFCModule.hh"

class TStnDebugModule;
class TStnGeneratorModule;		// event generator base class 
class TStnLumiMonModule;

// #include "Stntuple/ana/TClusterAnaModule.hh"
// #include "Stntuple/ana/TGenAnaModule.hh"
// #include "Stntuple/ana/TSpmcAnaModule.hh"

namespace stntuple {

  class TCrvAnaModule;
  class TEventDisplayModule;
  class TGenAnaModule;
  class THelixAnaModule;
  class TPhotosAnaModule;
  class TTrackAnaModule;
  class TTrackSeedAnaModule;
  class TTriggerAnaModule;
  class TValidationModule;

  TClusterAnaModule*          m_cls   = nullptr;
  TCrvAnaModule*              m_crv   = nullptr;
  TEventDisplayModule*        m_evd   = nullptr;
  TGenAnaModule*              m_gen   = nullptr;
  THelixAnaModule*            m_hel   = nullptr;
  TPhotosAnaModule*           m_pho   = nullptr;
  TSpmcAnaModule*             m_spmc  = nullptr;
  TTrackAnaModule*            m_trk   = nullptr;
  TTrackSeedAnaModule*        m_trs   = nullptr;
  TTriggerAnaModule*          m_trig  = nullptr;
  TValidationModule*          m_val   = nullptr;

  TStnDebugModule*            m_dbg   = nullptr;
  TDFCModule*                 m_dfc   = nullptr;
  TStnLumiMonModule*          m_lum   = nullptr;
  TStnGeneratorModule*        m_stg   = nullptr;

};

#endif
