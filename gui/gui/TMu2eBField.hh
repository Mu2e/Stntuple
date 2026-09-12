//-----------------------------------------------------------------------------
// apparenty , the field is in Tesla, not in Gauss
//-----------------------------------------------------------------------------
#ifndef Stntuple_ana_TMu2eBField
#define Stntuple_ana_TMu2eBField

#include "TObject.h"
#include "TEveTrackPropagator.h"
#include "TEveVector.h"

namespace mu2e {
  class BFieldManagerMaker;
  class BFieldConfig;
  class Beamline;
  class BFieldManager;
}

namespace stntuple {

class TMu2eBField: public TEveMagField {
protected:
  mu2e::Beamline*           fBeamline;
  mu2e::BFieldConfig*       fBfc;
  mu2e::BFieldManagerMaker* fBfmm;
  mu2e::BFieldManager*      fBmgr;
  
public:
  TMu2eBField();
  ~TMu2eBField();
  
  virtual Float_t    GetMaxFieldMag() const override { return 5. ; }
  virtual TEveVector GetField(float X, float Y, float Z) const override;

  ClassDefOverride(TMu2eBField, 0);
};
}
#endif
