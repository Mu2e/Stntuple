#ifndef __daqana_obj_TCaloDigi_hh__
#define __daqana_obj_TCaloDigi_hh__

#include <vector>
#include "TClonesArray.h"
#include "TObject.h"

class TCaloDigi : public TObject {
public:
                                        // digi mask occupies high 16 bits of the fSipmID
  enum {
    kOverflowFlag = 0x0001,
  };
  
  int                   fNs;
  int                   fSipmID;
  int                   fT0;
  int                   fPPos;          // peak position
  std::vector<uint16_t> fWf;
//-----------------------------------------------------------------------------
// functions
//-----------------------------------------------------------------------------
  TCaloDigi();
  TCaloDigi(int ID);
  virtual ~TCaloDigi();

  int     Ns    () const { return fNs; }
  int     SipmID() const { return fSipmID & 0xffff ; }
  int     Mask  () const { return (fSipmID >> 16) & 0xffff; }
  int     T0    () const { return fT0; }
  int     PPos  () const { return fPPos; }

  int     Init  (int Ns);

  void    Set(int SipmID, float T0, float PeakPos, const std::vector<int>* Wf);

  void    SetMask(int AddedFlags) {
    int new_mask = Mask() | (AddedFlags & 0xffff);
    fSipmID      = (fSipmID & 0xffff) | (new_mask << 16);
  }

  std::vector<uint16_t>& Wf() { return fWf; }
//-----------------------------------------------------------------------------
// schema evolution
//-----------------------------------------------------------------------------
  void ReadV1(TBuffer &R__b);
//-----------------------------------------------------------------------------
// overloaded TObject functions
//-----------------------------------------------------------------------------
  virtual void Clear(const char* Opt = "") override ;
  virtual void Print(const char* Opt = "") const override ;

  ClassDefOverride(TCaloDigi,2);
};

#endif
