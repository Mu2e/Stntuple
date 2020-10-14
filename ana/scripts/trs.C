///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
#include "Stntuple/scripts/global_vars.h"
#include "Stntuple/ana/scripts/modules.hh"

def_name stn_trs_001("stn_trs_ana");

void stn_trs_ana(int DebugBit = -1, int PdgCode = 11, int GeneratorCode = 2) {
//-----------------------------------------------------------------------------
// configure validation module
//-----------------------------------------------------------------------------
  stntuple::m_trs = (TStnTrackSeedAnaModule*) g.x->AddModule("TStnTrackSeedAnaModule",0);  
  stntuple::m_trs->SetPdgCode      (PdgCode);
  stntuple::m_trs->SetGeneratorCode(GeneratorCode);

  if (DebugBit >= 0) stntuple::m_trs->SetDebugBit(DebugBit,1);
}
