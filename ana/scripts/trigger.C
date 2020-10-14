///////////////////////////////////////////////////////////////////////////////
// 
///////////////////////////////////////////////////////////////////////////////
#include "Stntuple/scripts/global_vars.h"
#include "Stntuple/ana/scripts/modules.hh"

def_name stn_trigger_0010("stn_trigger_ana");
///////////////////////////////////////////////////////////////////////////////

//-----------------------------------------------------------------------------
void  stn_trigger_ana(int PdgCode = 11, int ProcessCode = 7, int DebugBit = -1) {
//-----------------------------------------------------------------------------
// configure analysis module
//-----------------------------------------------------------------------------
  stntuple::m_trig = (TStnTriggerAnaModule*) g.x->AddModule("TStnTriggerAnaModule",0);
  stntuple::m_trig->SetPdgCode(PdgCode);
  stntuple::m_trig->SetProcessCode(ProcessCode);
  if (DebugBit >= 0) {
    stntuple::m_trig->SetDebugBit(DebugBit,1);
  }
}

