//-----------------------------------------------------------------------------
//  2022-09-24 P.Murat - TStrWaveform
//-----------------------------------------------------------------------------
#ifndef TStrWaveform_hh
#define TStrWaveform_hh

#include <math.h>
#include "TMath.h"
#include "TObject.h"
#include "TBuffer.h"

class TStrWaveform : public TObject {
public:
					// data
  int     fNWords;
  ushort* fData;                        // [fNWords] array of fNWords shorts

  float   fBaseline;	                //! baseline (based on first 5 samples)
  float   fQn;                          //! nsamples used to calculate the charge
  float   fVMin;			//! min and max of the waveform
  float   fVMax;

public:
                                        // constructors and destructors
  TStrWaveform(int ID = -1);
  virtual ~TStrWaveform();
//-----------------------------------------------------------------------------
// accessors
//-----------------------------------------------------------------------------
  int     NWords   () { return fNWords ; }
  ushort  Data(int I) { return fData[I]; }
  ushort* Data()      { return fData   ; }
//-----------------------------------------------------------------------------
// modifiers
//-----------------------------------------------------------------------------
  void    Set(int NWords, const ushort* Data);
//-----------------------------------------------------------------------------
// overloaded methods of TObject
//-----------------------------------------------------------------------------
  virtual void Clear(Option_t* opt = "")       override;
  virtual void Print(Option_t* opt = "") const override;
//-----------------------------------------------------------------------------
// schema evolution - no I/O changes from v1 to v2, only transient variables added
//-----------------------------------------------------------------------------
//  void ReadV1(TBuffer &R__b);

  ClassDefOverride (TStrWaveform,2)
};

#endif
