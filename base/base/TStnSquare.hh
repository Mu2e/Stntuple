///////////////////////////////////////////////////////////////////////////////
// HexSize - distance between the 2 flat sides
///////////////////////////////////////////////////////////////////////////////
#ifndef TStnSquare_hh
#define TStnSquare_hh

#include "TStnShape.hh"
//-----------------------------------------------------------------------------
class TStnSquare: public TStnShape {
public:
//-----------------------------------------------------------------------------
// methods
//-----------------------------------------------------------------------------
  TStnSquare();
  ~TStnSquare();
  TStnSquare(double HexSize, double X = 0, double Y = 0);

  virtual void Paint(Option_t* Opt="") override;

  static int Test(double Size = 30., double RIn=360., double ROut=670.);

  ClassDefOverride(TStnSquare,0)

};


#endif
