///////////////////////////////////////////////////////////////////////////////
// description of the trigger stored in STNTUPLE
// Author:    P. Murat (CDF/FNAL)
// Date:      Nov 09 2001
//            Nov 09 2019 : adopt for Mu2e
///////////////////////////////////////////////////////////////////////////////
#ifndef Stntuple_obj_TStnTrigger_hh
#define Stntuple_obj_TStnTrigger_hh

#include <iostream>

#include "TFolder.h"
#include "TFile.h"
#include "TString.h"

class TStnTrigger: public TObject {

protected:
  Int_t     fBit;
  Int_t     fTag;
  Int_t     fId;
  TString   fName;
  //-----------------------------------------------------------------------------
  //  methods
  //-----------------------------------------------------------------------------
public:
  // ****** constructors and destructor
  TStnTrigger();
  TStnTrigger(Int_t Bit, Int_t Id, const char* Name, Int_t Tag);
  ~TStnTrigger();
  // ****** accessors

  Int_t           Bit  () const { return fBit;   }
  Int_t           Id   () const { return fId;    }
  Int_t           Tag  () const { return fTag;   }
  const TString&  Name () const { return *((const TString*) &fName);  }
  //TString         Name () const { return fName;  }

  // ****** modifiers

  void   SetBit  (Int_t       Bit  ) { fBit   = Bit;   }
  void   SetTag  (Int_t       Tag  ) { fTag   = Tag;   }
  void   SetId   (Int_t       Id   ) { fId    = Id;    }
  void   SetName (const char* Name ) { fName  = Name;  }

  // ****** overloaded methods of TObject

  virtual const char* GetName () const override { return (const char*) fName.Data(); }

  virtual void  Print  (Option_t* Option = "") const override;
  
  ClassDefOverride(TStnTrigger,1)

};

#endif
