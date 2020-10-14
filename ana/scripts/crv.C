///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
#include "Stntuple/scripts/global_vars.h"
#include "modules.hh"

def_name stn_crv_001("stn_crv_ana");

//-----------------------------------------------------------------------------
void  stn_crv_ana(int UseAllPulses = 0, int DebugBit = -1) {
  stntuple::m_crv = (TStnCrvAnaModule*) g.x->AddModule("TStnCrvAnaModule",0);  

  stntuple::m_crv->SetUseAllPulses(UseAllPulses);

  if (DebugBit >= 0)stntuple::m_crv->SetDebugBit(DebugBit,1);
}
