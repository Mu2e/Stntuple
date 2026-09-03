//
// - Serves as the link between ROOT "events" (e.g. mouse-clicks) and the ART
//   event display service by providing a receiver slot for signals generated
//   by the ROOT events.  A ROOT dictionary needs to be generated for this.
//

#ifndef __Stntuple_ana_TEventDisplayUtils_hh__
#define __Stntuple_ana_TEventDisplayUtils_hh__

#include <TObject.h>
#include <TApplication.h>
#include <TGTextBuffer.h>
#include <iostream>

namespace stntuple {

class TEventDisplayUtils: public TObject {
public:

  TGTextBuffer *fTbRun;
  TGTextBuffer *fTbEvt;

  TEventDisplayUtils();
  ~TEventDisplayUtils();
  
  void     PrevEvent         ();
  void     NextEvent         ();
  void     GotoEvent         ();

  ClassDefOverride(TEventDisplayUtils,0)
};
}
#endif 
