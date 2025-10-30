/////////////////////////////////////////////////////////////////////////////
#ifndef __Stntuple_gui_print_buttons_h
#define __Stntuple_gui_print_buttons_h
//-----------------------------------------------------------------------------
// these are not the collection names/types - just button names
// they have to correspond to the name list in TStnVisManager::PrintColls
// move them into an include file ?
//-----------------------------------------------------------------------------
struct print_button_t {
  const char* name;
  void  (*f)();
};

print_button_t print_buttons[] = {
  {"CosmicTrackSeedColls", print_cts_colls       },
  {"HelixSeedColls"      , print_helix_seed_colls},
  {"KalSeedColls"        , print_kalseed_colls   },
  {"StrawDigiColls"      , print_sd_colls        },
  {"StrawDigiMCColls"    , print_sdmc_colls      },
  {"TimeClusterColls"    , print_tc_colls        },
  {"ComboHits"           , nullptr},
  {"CosmicTrackSeeds"    , nullptr},
  {"Helices"             , nullptr},
  {"KalSeeds"            , nullptr},
  {"SimParticles"        , nullptr},
  {"StrawDigis"          , nullptr},
  {"StrawHits"           , nullptr},
  {"TimeClusters"        , nullptr},
  {nullptr               , nullptr}
};
#endif
