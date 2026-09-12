#ifndef smooth_hh
#define smooth_hh

#include "TH1.h"
#include "TF1.h"
#include "TGraph.h"

class smooth {
public:

  int     fIndex;
  int     fN   ;
  double* fP0  ;
  double* fP1  ;
  double* fP2  ;
  double* fX   ;
  TF1*    fFunc;
//-----------------------------------------------------------------------------
// constructiors and destructor
//-----------------------------------------------------------------------------
  smooth();
  smooth(const TH1*    Hist , double XMin = -1., double XMax = 1.);
  smooth(const TGraph* Graph, double XMin = -1., double XMax = 1.);
  virtual ~smooth();

  double Eval(double* X);


  TF1*   GetFunc() { return fFunc; }

  static  double func(double* X, double* P);

  // ClassDef(smooth,0)
};

#endif
