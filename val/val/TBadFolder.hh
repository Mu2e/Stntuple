#ifndef TBadFolder_hh
#define TBadFolder_hh

#include "TFolder.h"

class TBadFolder: public  TFolder {
public:
  TBadFolder() {};
  virtual ~TBadFolder(){}


  ClassDefOverride(TBadFolder,1)
};
#endif
