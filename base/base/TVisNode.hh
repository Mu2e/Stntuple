#ifndef Stntuple_base_TVisVode_hh
#define Stntuple_base_TVisVode_hh

#include "TObject.h"
#include "TString.h"

class TVisNode: public TObject {
protected:
  TString    fName;
  TObject*   fClosestObject;
  int        fDist;

  static int fgDebugLevel;
public:
					// ****** constructors and destructor
  TVisNode(const char* name = "TVisNode");
  virtual ~TVisNode();
//-----------------------------------------------------------------------------
// accessors
//-----------------------------------------------------------------------------
  TObject*            GetClosestObject() { return fClosestObject; }
  
  virtual const char* GetName() const    { return fName.Data(); }

  static int          DebugLevel()       { return fgDebugLevel; }

					// generic print callback. to be overloaded, if needed

  virtual void        NodePrint(const void* Object, const char* ClassName) ;

				// called by TEvdManager::DisplayEvent. a must to overload

  virtual int         InitEvent() = 0;

  static  void        SetDebugLevel(int Level) { fgDebugLevel = Level; }

  void                SetClosestObject(TObject* Obj, int Dist) {
    fClosestObject = Obj;
    fDist          = Dist;
  }

  ClassDef(TVisNode,0)
};

#endif
