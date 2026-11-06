///////////////////////////////////////////////////////////////////////////////
// PM: this include is temporary and it will go away as soon
// as the DB-based approach is implemented
// in essence, it is a table prototype
///////////////////////////////////////////////////////////////////////////////
#ifndef __Stntuple_alg_TRunDb_hh__
#define __Stntuple_alg_TRunDb_hh__

#include "fhiclcpp/ParameterSet.h"

namespace stntuple {

class TRunDb {
public:
                                        // run types
  enum {
    e_BEAM            = 1,
    e_COSMICS         = 2,
    e_INTERNAL_PULSER = 3,
    e_PULSE_INJECTION = 4,
    e_MONTE_CARLO     = 5,
  };
     
  struct Data_t {
    int   run_number;                   //
    int   nevents;                      // total nevents on input (last event 
    int   run_type;                     // 1:beam 2:cosmics 3:internal_pulser 4:pulse_injection 5:noise 100:MC
    int   event_mode;                   // on/off spill
    int   cfo_rate;                     // event rate, Hz
    int   trigger_rate;                 // event rate to disk, Hz
    float ew_length;                    // event window length, sec
    int   max_ewtag;                    // max event window tag
    int   ref_channel;                  // for pulse injection runs
    int   pulsed_channel[12];           // assume one per 8
    int   plane_flag[36];               // 
  };

private:
  TRunDb();
  ~TRunDb();

public:
  static TRunDb*      Instance();
  
  int                 GetRunInfo(int RunNumber, TRunDb::Data_t* RunInfo);
  fhicl::ParameterSet GetTable(const char* TableName, int RunNumber);

};

}
#endif
