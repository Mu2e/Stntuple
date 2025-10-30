///////////////////////////////////////////////////////////////////////////////
// vis node displays one wedge
///////////////////////////////////////////////////////////////////////////////
#ifndef TEvdTimeCluster_hh
#define TEvdTimeCluster_hh

#include "Gtypes.h"
#include "TClonesArray.h"
#include "TH1.h"
#include "TBox.h"
#include "TEllipse.h"

namespace mu2e {
  class TimeCluster;
  class ComboHit;
  class ComboHitCollection;
}

class TStnVisNode;

namespace stntuple {

class TEvdTrkStrawHit;

class TEvdTimeCluster: public TObject {
public:
  
protected:
  int                       fNumber;
  const mu2e::TimeCluster*         fTimeCluster;
  const mu2e::ComboHitCollection*  fChColl;
  TStnVisNode*              fVisNode;      // backward link to the note

  TObjArray*                fListOfHits;   // list of combo hits
  float                     fT0; 
  float                     fTMin;
  float                     fTMax;
  float                     fZMin;
  float                     fZMax;
  float                     fPhiMin;
  float                     fPhiMax;
  int                       fColor;
  TBox*                     fBox;
  TEllipse*                 fEllipse;
public:
//-----------------------------------------------------------------------------
// constructors and destructor
//-----------------------------------------------------------------------------
  TEvdTimeCluster();

  TEvdTimeCluster(int Number,
                  const mu2e::TimeCluster* aTimeCluster,
                  const mu2e::ComboHitCollection* ChColl,
                  float T0,
		  float TMin, float TMax, float ZMin, float ZMax,
                  float PhiMin, float PhiMax,
                  TStnVisNode* VisNode);

  virtual ~TEvdTimeCluster();
//-----------------------------------------------------------------------------
// accessors
//-----------------------------------------------------------------------------
  const mu2e::TimeCluster* TimeCluster() const { return fTimeCluster; }
  float                    T0         () const { return fT0    ; }
  float                    TMin       () const { return fTMin  ; }
  float                    TMax       () const { return fTMax  ; }
  float                    ZMin       () const { return fZMin  ; }
  float                    ZMax       () const { return fZMax  ; }
  float                    PhiMin     () const { return fPhiMin; }
  float                    PhiMax     () const { return fPhiMax; }
  TBox*                    Box        () const { return fBox   ; }
//-----------------------------------------------------------------------------
// modifiers
//-----------------------------------------------------------------------------
					// Type: "ch" or "sh"
  int                      HitOK(const mu2e::ComboHit* Hit, const char* Type);

  int TCHit(int Index);
//-----------------------------------------------------------------------------
// drawing functions
//-----------------------------------------------------------------------------
  virtual void  PaintXY  (Option_t* option = "");
  virtual void  PaintTZ  (Option_t* option = "");

  virtual void  Select();                        // **MENU**

  virtual Int_t DistancetoPrimitiveTZ (Int_t px, Int_t py);
//-----------------------------------------------------------------------------
// overloaded methods of TObject
//-----------------------------------------------------------------------------
  virtual void  Paint(Option_t* Opt = "");
  virtual void  Clear(Option_t* Opt = "");

  virtual void  PrintMe() const ; // **MENU**
  virtual void  Print  (Option_t* Opt = "") const ; // **MENU**

  virtual int   DistancetoPrimitive(Int_t px, Int_t py);
  virtual void  ExecuteEvent       (int Event, int Px, int Py);	


  ClassDef(stntuple::TEvdTimeCluster,0)
};

}
#endif
