//
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <format>
#include <deque>
#include <vector>
#include <algorithm>

#include "Stntuple/geom/TCrvChannelMap.hh"

TCrvChannelMap* TCrvChannelMap::fgInstance(nullptr);

namespace {
//-----------------------------------------------------------------------------
  std::string trim(std::string value) {
    auto notSpace = [](unsigned char c) { return !std::isspace(c); };
    
    value.erase(value.begin(), std::find_if(value.begin(), value.end(), notSpace));
    value.erase(std::find_if(value.rbegin(), value.rend(), notSpace).base(),value.end());
    
    return value;
  }

//-----------------------------------------------------------------------------
  std::vector<std::string> splitCommaSeparated(std::string const& line) {
    std::vector<std::string> fields;
    std::string field;
    
    for (char c : line) {
      if (c == ',') {
        fields.push_back(trim(field));
        field.clear();
      } else {
        field += c;
      }
    }

    fields.push_back(trim(field));
    return fields;
  }

//-----------------------------------------------------------------------------
  int toInt(std::string const& value, std::string const& fileName, unsigned lineNumber) {
    try {
      std::size_t used = 0;
      int result = std::stoi(value, &used);
      
      if (used != value.size()) {
        throw std::invalid_argument("extra characters");
      }
      
      return result;
    }
    catch (std::exception const& error) {
      throw std::runtime_error(
                               fileName + ":" + std::to_string(lineNumber) +
                               ": invalid integer '" + value + "': " + error.what());
    }
  }
}

//-----------------------------------------------------------------------------
TCrvChannelMap::TCrvChannelMap() {
}

//-----------------------------------------------------------------------------
int TCrvChannelMap::Init(int RunNumber) {
  int rc(0);

  fRunNumber = RunNumber;
  // initialize, assume using spack
  
  std::string fn = std::format("{}/Stntuple/geom/data/crv_channel_map_extracted.txt",getenv("MUSE_WORK_DIR"));
    
  std::ifstream input(fn);
  if (!input) {
    throw std::runtime_error("Unable to open CRV channel map: " + fn);
  }

  std::string line;
  unsigned    lineNumber = 0;

  while (std::getline(input, line)) {
    ++lineNumber;

    std::string stripped = trim(line);

    if (stripped.empty() || stripped.front() == '#') {
      continue;
    }

    auto fields = splitCommaSeparated(stripped);

    // chid,sbid,sipm,sectorName,module,side,roc,feb,fch
    if (fields.size() != 9) {
      throw std::runtime_error(
          fn + ":" + std::to_string(lineNumber) +
          ": expected 8 comma-separated fields, found " +
          std::to_string(fields.size()));
    }

    Data_t data;
    data.chid       = toInt(fields[0], fn, lineNumber);
    data.sbid       = toInt(fields[1], fn, lineNumber);
    data.sipm       = toInt(fields[2], fn, lineNumber);
    if      (fields[3] == "CRV_EX") data.sector = 0;
    else if (fields[3] == "CRV_T1") data.sector = 1;
    else if (fields[3] == "CRV_T2") data.sector = 2;
    else if (fields[3].find("CRV_M") == 0) data.sector = 3;
    else {
      throw std::runtime_error(
          fn + ":" + std::to_string(lineNumber) +
          ": unknown sector:" + std::to_string(fields.size()));
    }

    data.module     = toInt(fields[4], fn, lineNumber);
    if      (fields[5] == "-") data.side       = 0;
    else if (fields[5] == "+") data.side       = 1;
    else {
      throw std::runtime_error(
          fn + ":" + std::to_string(lineNumber) +
          ": unknown side:" + std::to_string(fields.size()));
    }
    data.roc        = toInt(fields[6], fn, lineNumber);
    data.feb        = toInt(fields[7], fn, lineNumber);
    data.fch        = toInt(fields[8], fn, lineNumber);

    _data.push_back(data);

    Data_t* p = &_data.back();
    _ch_data_by_offline[data.chid]                     = p;
    _ch_data_by_online [data.roc ][data.feb][data.fch] = p;
  }

  std::cout << "Read " << _data.size() << " channel mappings\n";
  return rc;
}

//-----------------------------------------------------------------------------
TCrvChannelMap::~TCrvChannelMap() {
}


//-----------------------------------------------------------------------------
TCrvChannelMap* TCrvChannelMap::Instance() {
  if  (fgInstance == nullptr) {
    fgInstance = new TCrvChannelMap();
  }
  return fgInstance;
}
