//////////////////////////////////////////////////////////////////////////////
#ifndef __Stntuple_print_functions_hh__
#define __Stntuple_print_functions_hh__

// #ifndef __CINT__

#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"

// #else

// namespace art {
//   class Event;
// }

// #endif

void     print_combo_hit_coll  (const char* ComboHitCollTag, const char* StrawHitCollTag = "makeSH");
void     print_combo_hit_coll_2(const char* ComboHitCollTag, const char* StrawHitCollTag = "makeSH");

void     print_calo_cluster_colls();   // print all CaloCluster collections
void     print_ch_colls          ();   // print all ComboHit     collections
void     print_cts_colls         ();   // 'cts' = CosmicTrack
void     print_genp_colls        ();   // print all GenParticle  collections
void     print_helix_seed_colls  ();   // print all HelixSeed    collections
void     print_kalrep_colls      ();   // print all KalRep       collections
void     print_kalseed_colls     ();   // print all KalSeed      collections
void     print_sd_colls          ();   // print all StrawDigi    collections
void     print_sdmc_colls        ();   // print all StrawDigiMC  collections
void     print_sh_colls          ();   // print all StrawHit     collections
void     print_shf_colls         ();   // print all StrawHitFlag collections
void     print_simp_colls        ();   // print all SimParticle  collections
void     print_spmc_colls        ();   // print all StepPointMC  collections
void     print_tc_colls          ();   // print all time cluster collections

#endif
