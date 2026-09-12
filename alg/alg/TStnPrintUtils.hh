///////////////////////////////////////////////////////////////////////////////
// 2026-09-08: P.Murat 
///////////////////////////////////////////////////////////////////////////////
#ifndef __Stntuple_alg_PrintUtils_hh__
#define __Stntuple_alg_PrintUtils_hh__

#include "TObject.h"

class TStnTimeClusterBlock;
class TComboHitBlock;

class TStnPrintUtils: public TObject {
public:
//-----------------------------------------------------------------------------
// 
//-----------------------------------------------------------------------------
private:
  TStnPrintUtils();
  
public:
  static TStnPrintUtils* Instance();
  
  void PrintTimeClusterBlock(TStnTimeClusterBlock* TcB, TComboHitBlock* ChB);
};

#endif
