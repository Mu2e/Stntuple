///////////////////////////////////////////////////////////////////////////////
// so far: a base class for a square and a hexagon
///////////////////////////////////////////////////////////////////////////////
#ifndef TStnShape_hh
#define TStnShape_hh

#include <math.h>

#include "TObject.h"
//-----------------------------------------------------------------------------
class TStnShape: public TObject {
public:
  double fSize;
  double fX0;
  double fY0;
  
  int    fLineColor;
  int    fFillColor;
  int    fFillStyle;
  int    fLineWidth;
//-----------------------------------------------------------------------------
// methods
//-----------------------------------------------------------------------------
  TStnShape();
  ~TStnShape();
  TStnShape(double Size, double X = 0, double Y = 0);

  double X0    () const { return fX0; }
  double Y0    () const { return fY0; }
  double Radius() const { return sqrt (fX0*fX0+fY0*fY0); }

  double Dist(double X, double Y) const { return sqrt((fX0-X)*(fX0-X)+(fY0-Y)*(fY0-Y)); }

  void   SetSize(double Size) { fSize = Size; }
  void   SetPos (double X, double Y) { fX0 = X; fY0 = Y; }

  void  SetLineColor(int Color) { fLineColor = Color; }
  void  SetLineWidth(int Width) { fLineWidth = Width; }
  void  SetFillColor(int Color) { fFillColor = Color; }
  void  SetFillStyle(int Style) { fFillStyle = Style; }

  virtual void  Paint(Option_t* Opt="") = 0;

  ClassDefOverride(TStnShape,0)

};


#endif
