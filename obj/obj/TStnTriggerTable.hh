#ifndef Stntuple_obj_TStnTriggerTable
#define Stntuple_obj_TStnTriggerTable
//-----------------------------------------------------------------------------
// description of the trigger table to be stored in STNTUPLE
// Author:    P. Murat (CDF/FNAL)
// Date:      Nov 2 2001
//            Nov 09 2019: adopt for Mu2e
//-----------------------------------------------------------------------------

#include "TObjArray.h"
#include "Stntuple/obj/TStnTrigger.hh"

class TStnTriggerTable: public TObject {
protected:
  TString    fObjName;                  //! name of this for list searching

protected:
  Int_t      fTag;			// trigger table tag
  TString    fName;			// trigger table name
  TObjArray* fListOfTriggers;	        // descriptions of the triggers
//-----------------------------------------------------------------------------
//  methods
//-----------------------------------------------------------------------------
public:
  TStnTriggerTable();
  ~TStnTriggerTable();
					// ****** accessors

  const char* GetTableName() const { return fName.Data(); }
  int         GetTableTag () const { return fTag;         }
  int         GetVersion  () const ; 

  Int_t NTriggers() const { return fListOfTriggers->GetEntriesFast(); }

  const TStnTrigger*  GetTrigger(int I) const {
    return (TStnTrigger*) fListOfTriggers->UncheckedAt(I);
  }
					// return list of triggers which names 
					// match given pattern

  void GetListOfTriggers(const char* Pattern, TObjArray* List);

					// ****** modifiers

  void SetName   (const char* Name) { fName = Name; }
  void SetTag    (Int_t       Tag ) { fTag  = Tag;  }

  void AddTrigger(const TStnTrigger* Trigger);

					// ****** overloaded methods of TObject

  const char* GetName() const { return fObjName.Data(); }
  void        Delete(Option_t* option = "");
  void        Print (Option_t* option = "") const; // *MENU*

  ClassDefOverride(TStnTriggerTable,1)
};

#endif
