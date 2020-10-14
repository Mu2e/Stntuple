///////////////////////////////////////////////////////////////////////////////
// classes in Stntuple/ana murat/ana have the same names, don't load 
// libmurat_ana.so and libStntuple_ana.so simultaneously
///////////////////////////////////////////////////////////////////////////////
#include "TInterpreter.h"
#include "modules.hh"

// TG3Pythia6*              py      = NULL;
//-----------------------------------------------------------------------------
// the first parameter is the script, the second - env.var telling whether 
// the script has to be loaded. If the corresponding env. variable is not defined,
// the script is not loaded. PWD is always defined
//-----------------------------------------------------------------------------
int load_stnana_scripts_Stntuple() {
  char        macro[200];

  const char* script[] = { 
    "catalog.C"    , "PWD",
    "crv.C"        , "PWD",
    "debug.C"      , "PWD",
    "helix.C"      , "PWD",
    "lumi.C"       , "PWD",
    "photos.C"     , "STNTUPLE_MC_GEN",
    "trk.C"        , "PWD",
    "pen.C"        , "PWD",
    "trs.C"        , "PWD",                  // sample track seed analysis module
    "stn_trigger.C", "PWD",
    "validation.C" , "PWD",
    0 
  };

  TString work_dir = gEnv->GetValue("Stnana.TestReleaseDir",gSystem->Getenv("PWD"));

  TInterpreter* cint = gROOT->GetInterpreter();
  
  for (int i=0; script[i] != 0; i+=2) {
    const char* dir = gSystem->Getenv(script[i+1]);
    if (dir) {
      sprintf(macro,"%s/Stntuple/ana/scripts/%s",dir,script[i]);
      if (! cint->IsLoaded(macro)) cint->LoadMacro(macro);
    }
  }
  
  return 0;
}
