///////////////////////////////////////////////////////////////////////////////
// 
///////////////////////////////////////////////////////////////////////////////
#include "Stntuple/scripts/global_vars.h"
#include "Stntuple/ana/scripts/modules.hh"

def_name stn_cal_001("stn_cluster_ana");
///////////////////////////////////////////////////////////////////////////////

//-----------------------------------------------------------------------------
// configure jobs
//-----------------------------------------------------------------------------
void  stn_cluster_ana(int DebugBit = -1) {
  stntuple::m_cls = (TStnClusterAnaModule*) g.x->AddModule("TStnClusterAnaModule",0);  
  if (DebugBit >= 0) stntuple::m_cls->SetDebugBit(DebugBit,1);
}





