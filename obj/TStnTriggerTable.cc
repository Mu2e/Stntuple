///////////////////////////////////////////////////////////////////////////////
//  STNTUPLE trigger table to be stored in the STNTUPLE DB subdirectory
//
//  Author:    P.Murat (CDF/FNAL)
//  Date:      Oct 15 2001
///////////////////////////////////////////////////////////////////////////////
#include <cstdlib>
#include "obj/TStnTriggerTable.hh"


ClassImp(TStnTriggerTable)

//_____________________________________________________________________________
TStnTriggerTable::TStnTriggerTable():fObjName("TStnTriggerTable") {
  fListOfTriggers = new TObjArray(128);
}

//_____________________________________________________________________________
TStnTriggerTable::~TStnTriggerTable() {
  fListOfTriggers->Delete();
  delete fListOfTriggers;
}

//_____________________________________________________________________________
void TStnTriggerTable::AddTrigger(const TStnTrigger* Trigger) {
  Int_t bit = Trigger->Bit();

  fListOfTriggers->AddAtAndExpand((TObject*) Trigger,bit);
}

//_____________________________________________________________________________
Int_t TStnTriggerTable::GetVersion() const {
  // parse name : AAAA_N_MM , N: high number of the version
  // MM - lower, return 100*N+MM

  int len = fName.Length();
  int ver = atoi(fName(len-4,1).Data())*100 + atoi(fName(len-2,2).Data());

  return ver;
}

//_____________________________________________________________________________
void TStnTriggerTable::GetListOfTriggers(const char* Pattern,
					 TObjArray*  List) 
{
  // do not assume that the trigger names are all uppercased

  List->Clear();
  Int_t nt = NTriggers();

  TString trigger_name;
  TString pat(Pattern);

  pat.ToUpper();
  for (int i=0; i<nt; i++) {
    const TStnTrigger* t = GetTrigger(i);

				// keep in mind that not all the triger bits
				// have corresponding triggers implemented, so
				// a "proof of existence" is necessary
    if (t) {
      trigger_name = t->Name();
      trigger_name.ToUpper();
      if (pat == "" || strstr(trigger_name,pat) != 0) {
	List->Add((TObject*)t);
      }
    }
  }
}

//_____________________________________________________________________________
void TStnTriggerTable::Delete(Option_t* Opt) {
  fListOfTriggers->Delete();
}

//_____________________________________________________________________________
void TStnTriggerTable::Print(Option_t* Opt) const {
  printf(" **** trigger table %s ",fName.Data());
  printf(" N(triggers) = %3i \n",fListOfTriggers->GetEntries());

  int banner_printed;

  int nt = fListOfTriggers->GetSize();
  banner_printed = 0;
  for (int i=0; i<nt; i++) {
    const TStnTrigger* trig = GetTrigger(i);

    // some of the trigger bits may not be used
    if (trig) {
      if (!banner_printed) {
	trig->Print("banner");
	banner_printed = 1;
      }
      trig->Print("data");
    }
  }
}

