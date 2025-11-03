///////////////////////////////////////////////////////////////////////////////
// hide STNTUPLE module classes in a namespace
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
  const char* pkg{"Stntuple"};
  char        macro[200];

  const char* script[] = { 
    "calorimeter.C"  , "PWD",		     // calorimeter, calorimeter clusters
    "catalog.C"      , "PWD",
    "crv.C"          , "PWD",		     // CRV
    "debug.C"        , "PWD",
    "genp.C"         , "PWD",
    "event_display.C", "PWD",
    "helix.C"        , "PWD",
    "lumi.C"         , "PWD",
    "photos.C"       , "STNTUPLE_MC_GEN",    // analyze output of the PHOTOS MC
    "spmc.C"         , "PWD",                // analyze StepPointMC's
    "trk.C"          , "PWD",		     // tracks
    "trs.C"          , "PWD",                // sample track seed analysis module
    "trigger.C"      , "PWD",
    "validation.C"   , "PWD",
    0 
  };

  TString code_dir = gEnv->GetValue("Stnana.TestReleaseDir",gSystem->Getenv("PWD"));

  TInterpreter* cint = gROOT->GetInterpreter();
  
  for (int i=0; script[i] != 0; i+=2) {
    const char* dir = gSystem->Getenv(script[i+1]);
    if (dir) {
      if (strcmp(script[i+1],"PWD") == 0) sprintf(macro,"%s/%s/ana/scripts/%s",code_dir.Data(),pkg,script[i]);
      else                                       sprintf(macro,"%s/%s/ana/scripts/%s",dir            ,pkg,script[i]);
      if (! cint->IsLoaded(macro)) cint->LoadMacro(macro);
    }
  }
  
  return 0;
}
