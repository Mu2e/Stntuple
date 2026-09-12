///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
#include "Stntuple/scripts/global_vars.h"
#include "modules.hh"

def_name stn_crv_001("stn_crv_ana");
def_name stn_crv_002("stn_crv_time_ana");

//-----------------------------------------------------------------------------
void  stn_crv_ana(int UseAllPulses = 0, int DebugBit = -1) {
  stntuple::m_crv = (stntuple::TCrvAnaModule*) g.x->AddModule("stntuple::TCrvAnaModule",0);  

  stntuple::m_crv->SetUseAllPulses(UseAllPulses);

  if (DebugBit >= 0)stntuple::m_crv->SetDebugBit(DebugBit,1);
}

//-----------------------------------------------------------------------------
void  stn_crv_time_ana(int DebugBit = -1) {
  stntuple::m_crvt = (stntuple::TCrvTimeAnaModule*) g.x->AddModule("stntuple::TCrvTimeAnaModule",0);  

  if (DebugBit >= 0)stntuple::m_crvt->SetDebugBit(DebugBit,1);
}
