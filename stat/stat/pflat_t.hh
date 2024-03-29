#ifndef __Stntuple_stat_pflat_t_hh__
#define __Stntuple_stat_pflat_t_hh__

#include "Stntuple/stat/parameter_t.hh"

namespace stntuple {
  class  pflat_t : public parameter_t {
  public:
    double   fXMin;
    double   fXMax;
// -----------------------------------------------------------------------------
// functions
//-----------------------------------------------------------------------------

    // this is good for gaussian, log-normal, uniform
    pflat_t(const char* Name, double XMin, double XMax, int Debug = 0);
    
    virtual double XMin();
    virtual double XMax();
    
    virtual void InitValue();

    virtual void Print(const Option_t* Opt) const ;
    
  };
}
#endif
