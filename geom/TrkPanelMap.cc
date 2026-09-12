//
#include <iostream>
#include <string>
#include <vector>
#include <format>

#include "cetlib/filepath_maker.h"
#include "fhiclcpp/ParameterSet.h"
// #include "fhiclcpp/make_ParameterSet.h"

#include "TEnv.h"

// #include "toml++/toml.hpp"
#include "Stntuple/geom/TrkPanelMap.hh"

//-----------------------------------------------------------------------------
TrkPanelMap::TrkPanelMap() {
}

//-----------------------------------------------------------------------------
TrkPanelMap::~TrkPanelMap() {
}

//-----------------------------------------------------------------------------
TrkPanelMap* TrkPanelMap::Instance() {
  static TrkPanelMap* instance(nullptr);

  if  (instance == nullptr) {
    instance = new TrkPanelMap();
  }
  return instance;
}
/*
//-----------------------------------------------------------------------------
int TrkPanelMap::Init(int RunNumber) {
  fRunNumber = RunNumber;
    // initialize
                                        // assume using spack
  std::string fn = std::format("{}/rundb/TrkPanelMap/TrkPanelMap.toml",getenv("SPACK_ENV"));
  toml::table tbl = toml::parse_file(fn);
    
  auto maps = tbl["TrkPanelMap"].as_array();
  //  int n_run_ranges = maps->size();

  for (auto&& node : *maps) {
    toml::table& entry_table = *node.as_table();

    auto* range = entry_table["run_range"].as_array();

    int min_run = range->at(0).as_integer()->get();
    int max_run = range->at(1).as_integer()->get();

    if ((RunNumber >= min_run) and (RunNumber <= max_run)) {
      // found run range - panel data - array of records

      auto panel_data_array = entry_table["panel_data"].as_array();
      for (auto&& item : *panel_data_array) {
        toml::table& p = *item.as_table();
        
        int mnid   = p["mnid" ].value_or(-1);
        Data_t* r = &_data[mnid];
        r->mnid = mnid;
        r->dtc_id = p["dtc_id"].value_or(-1);
        r->link   = p["link"  ].value_or(-1);
        r->plane  = p["plane" ].value_or(-1);
        r->ppid   = p["ppid"  ].value_or(-1);
        r->panel  = p["panel" ].value_or(-1);
        r->zface  = p["zface" ].value_or(-1);

        _panel_data_by_mnid[mnid] = r;
        _panel_data_by_online[r->dtc_id][r->link] = r;
        _panel_data_by_offline[r->plane][r->panel] = r;
      }
    }
  }
}
*/

//-----------------------------------------------------------------------------
int TrkPanelMap::Init(int RunNumber) {
  int rc(0);
  
  cet::filepath_lookup policy("FHICL_FILE_PATH");

  std::string rundb_dir = gEnv->GetValue("Stntuple.RunDb","");
  std::string fn        = std::format("{}/TrkPanelMap/TrkPanelMap.fcl",rundb_dir);
  auto const pset       = fhicl::ParameterSet::make(fn, policy);

  // TrkPanelMap is a sequence of tables
  auto trkPanelMap = pset.get<std::vector<fhicl::ParameterSet>>("TrkPanelMap");

  for (size_t i = 0; i < trkPanelMap.size(); ++i) {
    auto const& entry = trkPanelMap[i];
    auto run_range    = entry.get<std::vector<int>>("run_range");

    if ((RunNumber >= run_range[0]) and (RunNumber <= run_range[1])) {
      // found the needed run range, read the panel data

      auto panel_data = entry.get<std::vector<fhicl::ParameterSet>>("panel_data");

      // std::cout << "  panel_data size = " << panel_data.size() << "\n";

      for (size_t j = 0; j < panel_data.size(); ++j) {
        auto const& p = panel_data[j];
        int mnid   = p.get<int>("mnid");
        Data_t* r  = &_data[mnid];

        r->mnid    = mnid;
        r->dtc_id = p.get<int>("dtc_id");
        r->link   = p.get<int>("link");
        r->plane  = p.get<int>("plane");
        r->ppid   = p.get<int>("ppid");
        r->panel  = p.get<int>("panel");
        r->zface  = p.get<int>("zface");
        
        _panel_data_by_mnid   [mnid]                = r;
        _panel_data_by_online [r->dtc_id][r->link ] = r;
        _panel_data_by_offline[r->plane ][r->panel] = r;
      }
    }
  }

  return rc;
}
