///////////////////////////////////////////////////////////////////////////////
// HexSize - distance between the 2 flat sides
///////////////////////////////////////////////////////////////////////////////
#ifndef test_hexmap_hh
#define test_hexmap_hh

#include "TStnShape.hh"
//-----------------------------------------------------------------------------
class TStnHexagon: public TStnShape {
public:
//-----------------------------------------------------------------------------
// methods
//-----------------------------------------------------------------------------
  TStnHexagon();
  ~TStnHexagon();
  TStnHexagon(double HexSize, double X = 0, double Y = 0);

  void  Paint(Option_t* Opt="") ;

  static int Test(double HexSize = 30., double RIn=360., double ROut=670.);

  ClassDefOverride(TStnHexagon,0)

};


#endif
