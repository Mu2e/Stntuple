#ifndef TStnErrorLogger_hh
#define TStnErrorLogger_hh

#include "TQObject.h"
#include "RQ_OBJECT.h"

#include "TNamed.h"
#include "TBuffer.h"

class TStnErrorLogger : public TNamed {
  RQ_OBJECT("TStnErrorLogger")
public:
//------------------------------------------------------------------------------
//  function members
//------------------------------------------------------------------------------
  TStnErrorLogger(const char* Name  = "StntupleErrorLogger",
		  const char* Title = "Error Logger");
  virtual ~TStnErrorLogger();
					// ****** overloaded methods of TObject

  void  Report(Int_t ErrorCode, const char* Message);  //*SIGNAL*

  virtual void  Clear(Option_t* opt="")       override;
  virtual void  Print(Option_t* opt="") const override;

  ClassDefOverride(TStnErrorLogger,0)

};

#endif




