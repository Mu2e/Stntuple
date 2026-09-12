//-----------------------------------------------------------------------------
//  2026-09-07 P.Murat - TComboHit
//-----------------------------------------------------------------------------
#ifndef TComboHit_hh
#define TComboHit_hh

#include <math.h>
#include "TMath.h"
#include "TObject.h"
#include "TBuffer.h"

#include "Offline/DataProducts/inc/StrawId.hh"

#include "Stntuple/obj/TStnLinkBlock.hh"

namespace mu2e {
  class ComboHit;
};

class TComboHit : public TObject {
public:
  int     fStrawID;                     // includes MC flag (bit 16)
  int     fNsh;
  int     fZface;                        // z-ordered face
  int     fMnid;                         // Minnesota panel ID 
  float   fTime;                         // correctedTime
  float   fDrTime;                       // drift time
  float   fX;
  float   fY;
  float   fZ;
  float   fUx;                           // nz is always zero
  float   fUy;
  float   fUres;
  float   fVres;
  float   fEDep;

  mu2e::ComboHit* fOfflineCh;           //! backward pointer

public:
                                        // constructors and destructors
  TComboHit(int I = -1);
  virtual ~TComboHit();
//-----------------------------------------------------------------------------
// accessors
//-----------------------------------------------------------------------------
  mu2e::StrawId Sid() const { return mu2e::StrawId( uint16_t(fStrawID & 0xffff)); }
  int   MCFlag     () const { return (fStrawID & 0xffff0000 ) >> 16         ; }

  int   Nsh        () const { return fNsh;  }
  int   Mnid       () const { return fMnid; }
  float Time       () const { return fTime; }
  float DrTime     () const { return fDrTime; }
  float X          () const { return fX; }
  float Y          () const { return fY; }
  float Z          () const { return fZ; }
  float EDep       () const { return fEDep;        }

  mu2e::ComboHit* OfflineCh() { return fOfflineCh; }
//-----------------------------------------------------------------------------
// modifiers, assume TOT = tot[0] | (tot[1] << 16
//-----------------------------------------------------------------------------
  void    Set(int   StrawID, float Time , float* TOT  , 
              float EDep   , float  McMom);
//-----------------------------------------------------------------------------
// overloaded methods of TObject
//-----------------------------------------------------------------------------
  virtual void Clear(Option_t* opt = "") override;
  virtual void Print(Option_t* opt = "") const override;
//-----------------------------------------------------------------------------
// schema evolution
//-----------------------------------------------------------------------------
  // void ReadV1(TBuffer& R__b);
  // void ReadV2(TBuffer& R__b);

  ClassDefOverride (TComboHit,1)
};

#endif
