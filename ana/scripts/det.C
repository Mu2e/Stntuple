///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
#include "Stntuple/scripts/global_vars.h"
#include "modules.hh"

def_name stn_det_002("stn_det_time_ana");

//-----------------------------------------------------------------------------
void  stn_det_time_ana(int DebugBit = -1) {
  stntuple::m_det = (stntuple::TDetTimeAnaModule*) g.x->AddModule("stntuple::TDetTimeAnaModule",0);  

  if (DebugBit >= 0)stntuple::m_det->SetDebugBit(DebugBit,1);
}
