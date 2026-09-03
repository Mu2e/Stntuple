//-----------------------------------------------------------------------------
//  2014-01-26 P.Murat - TTrackStrawHit
//-----------------------------------------------------------------------------
#ifndef TTrackStrawHit_hh
#define TTrackStrawHit_hh

#include <math.h>
#include "TMath.h"
#include "TObject.h"
#include "TBuffer.h"

class TTrackStrawHit : public TObject {
public:
//-----------------------------------------------------------------------------
// constants imported from DataProducts/inc/StrawId.hh
// indexing of the calibration and HV ends 
//-----------------------------------------------------------------------------
  enum { eCal = 0, eHv = 1 } ;

  constexpr static ushort _layermsk    = 0x1;                 // mask for layer field
  constexpr static ushort _strawmsk    = 0x7F;                // mask for straw field
  constexpr static ushort _preampmsk   = 0x7E;                // mask for preamp
  constexpr static ushort _panelmsk    = 0x380;               // mask for panel field
  constexpr static ushort _preampsft   = 1;                   // shift for preamp field
  constexpr static ushort _panelsft    = 7;                   // shift for panel field
  constexpr static ushort _facemsk     = 0x80;                // mask for face field
  constexpr static ushort _facesft     = 7;                   // shift for face field
  constexpr static ushort _planemsk    = 0xFC00;              // mask for plane field
  constexpr static ushort _planesft    = 10;                  // shift for plane field
  constexpr static ushort _stationmsk  = 0xF800;              // mask for station field
  constexpr static ushort _stationsft  = 11;                  // shift for station field
  constexpr static ushort _invalid     = 0xFFFF;              // invalid identifier
  constexpr static ushort _nstraws     = 96;                  // number of straws per panel
  constexpr static ushort _nlayers     = 2;                   // number of layers per panel ; do we need it, see below
  constexpr static ushort _npanels     = 6;                   // number of panels per plane
  constexpr static ushort _nfaces      = 2;                   // number of faces in a plane
  constexpr static ushort _nplanes     = 36;                  // number of planes
  constexpr static ushort _nstations   = _nplanes/2;          // number of stations
  constexpr static ushort _nupanels    = _npanels * _nplanes; // number of unique panels
  constexpr static ushort _nustraws    = _nupanels* _nstraws; // number of unique straws
  constexpr static ushort _nustrawends = _nustraws*2;         // number of unique straw ends
  constexpr static ushort _ntotalfaces = _nfaces*_nplanes;    // tot N(faces)
  constexpr static ushort _maxval      = ((_nplanes -1) << _planesft) + ((_npanels -1) << _panelsft) + _nstraws; // maximum Id as uint16 value
public:
  int     fStrawID;
  int     fActive;
  int     fAmbig;
  int     fPdgCode;
  int     fMotherPdgCode;
  int     fGeneratorCode;
  int     fSimID;
  float   fTime;
  float   fDt;
  float   fEnergy;
  float   fDriftRadius;
  float   fMcDoca;
  float   fMcMomentum;			// MC particle momentum, assumed to be the last word
public:
                                        // constructors and destructors
  TTrackStrawHit();
  virtual ~TTrackStrawHit();
//-----------------------------------------------------------------------------
// accessors
//-----------------------------------------------------------------------------
  int     Index      () const { return fStrawID;  }
  int     StrawID    () const { return (fStrawID & 0x0000ffff )               ; }
  int     MCFlag     () const { return (fStrawID & 0xffff0000 ) >> 16         ; }
	  	     
  int     Station    () const { return (fStrawID & _stationmsk) >> _stationsft; }
  int     Panel      () const { return (fStrawID & _panelmsk  ) >> _panelsft  ; }
  int     Face       () const { return (fStrawID & _facemsk   ) >> _facesft   ; }
  int     Layer      () const { return (fStrawID & _layermsk  )               ; }
  int     Straw      () const { return (fStrawID & _strawmsk  )               ; }
  int     Preamp     () const { return (fStrawID & _preampmsk ) >> _preampsft ; }

  float   Time       () { return fTime;   }
  float   Dt         () { return fDt;     }
  float   Energy     () { return fEnergy; }
  float   DriftRadius() { return fDriftRadius; }
  int     Active     () { return fActive; }

  int     PdgCode      () const { return fPdgCode;       }
  int     MotherPdgCode() const { return fMotherPdgCode; }
  int     GeneratorCode() const { return fGeneratorCode & 0xFFFF; }
  int     StrawDigiFlag() const { return (fGeneratorCode >> 16) & 0xFF;  }
  int     SimID        () const { return fSimID;         }

  float   McDoca       () { return fMcDoca;        }
  float   McMomentum   () { return fMcMomentum;    }
//-----------------------------------------------------------------------------
// modifiers
//-----------------------------------------------------------------------------
  void    Set(int Index, float Time, float Dt, float EnergyDep,
	      int Active, int Ambig, float DriftRadius,
	      int PdgID, int MotherPdgID, int GenCode, int SimID, 
	      float McDoca, float McMomentum);
//-----------------------------------------------------------------------------
// overloaded methods of TObject
//-----------------------------------------------------------------------------
  void Clear(Option_t* opt = "");
  void Print(Option_t* opt = "") const;
//-----------------------------------------------------------------------------
// schema evolution
//-----------------------------------------------------------------------------
  void ReadV1(TBuffer &R__b);

  ClassDefOverride (TTrackStrawHit,1)
};

#endif
