#ifndef __su2020_sens_parameter_t_hh__
#define __su2020_sens_parameter_t_hh__

#include "TNamed.h"
#include "TH1.h"
#include "TH2.h"
#include "TRandom3.h"
#include "Math/GSLRandom.h"
#include "Math/Random.h"


namespace stntuple {
  class  parameter_t : public TNamed {
  public:
    double                    fValue;                   // current value
    double                    fMean;                    // mean (set parameter)
    ROOT::Math::RandomRanLux* fRng;
    TH1D*                     fHistPDF;			// PDF distribution
    int                       fDebug;
    int                       fFixed;                   // = 0 by default
// -----------------------------------------------------------------------------
// functions
//-----------------------------------------------------------------------------

    parameter_t (const char* Name = "undefined_parameter", int Debug = 0);
    ~parameter_t();

    double         InverseMean() { return 1./fMean; }
    double         Mean       () { return fMean   ; }

    virtual double  XMin();
    virtual double  XMax();

    virtual void    InitValue();
					// GetValue returns fluctuated mean

    double          GetValue    () { return fValue; };
    double          InverseValue() { return 1./fValue; };

    TH1D*           GetHistPDF() { return fHistPDF; };

    void            SetDebug(int    Value) { fDebug = Value; }
    void            SetFixed(int    Value) { fFixed = Value; }
    void            SetMean (double Mean ) { fMean  = Mean;  }
    
    virtual void    Print(const Option_t* Opt) const override;

    ClassDefOverride(stntuple::parameter_t,0)
    
  };
}

#endif
