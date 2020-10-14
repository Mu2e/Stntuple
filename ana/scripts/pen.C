///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
#include "Stntuple/scripts/global_vars.h"
#include "Stntuple/ana/scripts/modules.hh"

def_name piplusenu_001("piplusenu_ana");

void piplusenu_ana(int PdgCode = 11, int GeneratorCode = 2, int DebugBit = -1) {
//-----------------------------------------------------------------------------
// configure validation module
//-----------------------------------------------------------------------------
  stntuple::m_pen = (TStnPiplusenuAnaModule*) g.x->AddModule("TStnPiplusenuAnaModule",0);
  stntuple::m_pen->SetPdgCode      (PdgCode);
  stntuple::m_pen->SetGeneratorCode(GeneratorCode);

  if (DebugBit >= 0) {
    stntuple::m_pen->SetDebugBit(DebugBit,1);
  }
}
