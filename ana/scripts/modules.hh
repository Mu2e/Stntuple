#ifndef stntuple_ana_scripts_modules_hh
#define stntuple_ana_scripts_modules_hh

// #include "Stntuple/ana/TDFCModule.hh"

class TStnCrvAnaModule;
class TStnDebugModule;
class TStnEventDisplayModule;
class TStnGeneratorModule;
class TStnHelixAnaModule;
class TStnLumiMonModule;
class TStnPhotosAnaModule;
class TStnTrackAnaModule;
class TStnPiplusenuAnaModule;
class TStnTrackSeedAnaModule;
class TStnTriggerAnaModule;
class TStnValidationModule;

namespace stntuple {

  TStnDebugModule*            m_dbg   = NULL;
  TDFCModule*                 m_dfc   = NULL;
  TStnLumiMonModule*          m_lum   = NULL;
  TStnEventDisplayModule*     m_evd   = NULL;
  TStnCrvAnaModule*           m_crv   = NULL;
  TStnHelixAnaModule*         m_hel   = NULL;
  TStnGeneratorModule*        m_stg   = NULL;
  TStnPhotosAnaModule*        m_pho   = NULL;
  TStnTrackAnaModule*         m_trk   = NULL;
  TStnPiplusenuAnaModule*     m_pen   = NULL;
  TStnTrackSeedAnaModule*     m_trs   = NULL;
  TStnTriggerAnaModule*       m_trig  = NULL;
  TStnValidationModule*       m_val   = NULL;
};

#endif
