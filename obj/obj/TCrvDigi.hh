#ifndef __daqana_obj_TCrvDigi_hh__
#define __daqana_obj_TCrvDigi_hh__

#include <vector>
#include "TClonesArray.h"
#include "TObject.h"

class TCrvDigi : public TObject {
public:
  int    fNs;
  int    fSbid;
  int    fTdc;
  int    fNzs;
  int    fOddTs;
  int    fSipm;
  int    fRoc;
  int    fFeb;
  int    fCh;

  std::vector<uint16_t> fAdc;
//-----------------------------------------------------------------------------
// functions
//-----------------------------------------------------------------------------
  TCrvDigi();
  TCrvDigi(int Ns);
  virtual ~TCrvDigi();

  int     Chid() { return 4*fSbid+fSipm; }
  int     Sbid() { return fSbid; }
  int     Sipm() { return fSipm; }
  int     Roc () { return fRoc;  }
  int     Ns  () { return fNs;   }
  int     Init(int Ns);

  std::vector<uint16_t>& Adc() { return fAdc; }

  virtual void Clear(const char* Opt) override ;

  ClassDefOverride(TCrvDigi,1);
};

#endif
