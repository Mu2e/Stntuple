///////////////////////////////////////////////////////////////////////////////
// 
///////////////////////////////////////////////////////////////////////////////
#include "modules.hh"

def_name stn_genp_001("stn_genp_ana");
def_name stn_spmc_002("stn_spmc_ana");
def_name stn_spmc_003("stn_g4s4_ana");
///////////////////////////////////////////////////////////////////////////////

//-----------------------------------------------------------------------------
void  stn_genp_ana(int DebugBit = -1) {
  stntuple::m_gen = (TStnGenAnaModule*) g.x->AddModule("TStnGenAnaModule",0);  
  if (DebugBit >= 0) stntuple::m_gen->SetDebugBit(DebugBit,1);
}

//-----------------------------------------------------------------------------
void  stn_spmc_ana(int DebugBit = -1) {
  stntuple::m_spmc = (TStnSpmcAnaModule*) g.x->AddModule("TStnSpmcAnaModule",0);  
  if (DebugBit >= 0) stntuple::m_spmc->SetDebugBit(DebugBit,1);
}

//-----------------------------------------------------------------------------
void  stn_g4s4_ana(int DebugBit = -1) {
  stntuple::m_spmc = (TStnSpmcAnaModule*) g.x->AddModule("TStnSpmcAnaModule",0);  
  stntuple::m_spmc->SetSpmcBlockName("SpmcTrackerBlock");
  stntuple::m_spmc->SetVDetBlockName("SpmcVDetBlock");
  if (DebugBit >= 0) stntuple::m_spmc->SetDebugBit(DebugBit,1);
}

