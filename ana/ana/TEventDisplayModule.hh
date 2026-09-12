//-----------------------------------------------------------------------------

#ifndef __Stntuple_ana_TEventDisplayModule_hh__
#define __Stntuple_ana_TEventDisplayModule_hh__

// ROOT includes
// ... libCore
#include <TApplication.h>
#include <TString.h>
#include <TSystem.h>
// ... libRIO
#include <TFile.h>
// ... libGui
#include <TGString.h>
#include <TGLabel.h>
#include <TGButton.h>
#include <TGButtonGroup.h>
#include <TGTextEntry.h>
#include <TGTextView.h>
#include <TGLayout.h>
#include <TGTab.h>
#include <TG3DLine.h>
// ... libGeom
#include <TGeoManager.h>
#include <TGeoTube.h>
#include <TGeoCompositeShape.h>
#include <TGeoBoolNode.h>
// ... libEG
#include <TParticle.h>
// ... libRGL
#include <TGLViewer.h>
// ... libEve
#include <TEveManager.h>
#include <TEveEventManager.h>
#include <TEveBrowser.h>
#include <TEveGeoNode.h>
#include "TEveGeoShape.h"
#include "TEveViewer.h"
#include <TEveScene.h>
#include <TEveProjectionManager.h>
#include <TEveProjectionAxes.h>
#include <TEvePointSet.h>
#include <TEveTrack.h>
#include <TEveTrackPropagator.h>

#include "Stntuple/gui/TEventDisplayUtils.hh"
#include "Stntuple/gui/TMu2eBField.hh"

#include "Stntuple/loop/TStnModule.hh"
#include "Stntuple/obj/TStnTrackBlock.hh"
#include "Stntuple/obj/TGenpBlock.hh"
#include "Stntuple/obj/TSimpBlock.hh"
#include "Stntuple/gui/TStnGeoManager.hh"

namespace stntuple {

class TEventDisplayModule : public TStnModule {
protected:

  // Set by parameter set variables.
  
  int                     fDisplayTracks;
  bool                    drawHits_;
  double                  hitMarkerSize_;
  double                  fTrkMaxR;
  double                  fTrkMaxZ;
  double                  fTrkMaxStepSize;
  double                  camRotateCenterH_;
  double                  camRotateCenterV_;
  double                  camDollyDelta_;
  
  TEventDisplayUtils*     fEvdUtils;
  TEveGeoShape*           fSimpleGeom;

  TStnGeoManager*         fGeoManager;
  TMu2eBField*            fBField;
  TEveTrackPropagator*    fTrackPropagator;
  TEveTrackPropagator*    fBeamPropagator;

  TEveViewer              *fXYView;
  TEveViewer              *fRZView;
  TEveProjectionManager   *fXYMgr;
  TEveProjectionManager   *fRZMgr;
  TEveScene               *fDetXYScene;
  TEveScene               *fDetRZScene;
  TEveScene               *fEvtXYScene;
  TEveScene               *fEvtRZScene;
  
  TGTextEntry             *fTeRun,*fTeEvt;
  TGLabel                 *fTlRun,*fTlEvt;
  
  TEveTrackList           *fBeamList;
  TEveTrackList           *fTrackList;
  TEveElementList         *fHitsList;

  TStnTrackBlock*         fTrackBlock;
  TGenpBlock*             fGenpBlock;
  TSimpBlock*             fSimpBlock;

  int                     fDisplayBeamMuMinus;
  int                     fDisplayBeamMuPlus;
//-----------------------------------------------------------------------------
// functions
//-----------------------------------------------------------------------------
  public:
  TEventDisplayModule(const char* Name = "Stntuple_EventDisplay", const char* Title = "Stntuple_EventDisplay");
  ~TEventDisplayModule();

  virtual int  BeginJob() override           ;
  virtual int  BeginRun() override           ;
  virtual int  EndJob  () override           ;
  virtual int  Event   (int IEntry) override ;
  
  void MakeNavPanel();
  
  void SetDisplayBeamMuMinus(int OnOff) { fDisplayBeamMuMinus = OnOff; }
  void SetDisplayBeamMuPlus (int OnOff) { fDisplayBeamMuPlus  = OnOff; }

  ClassDefOverride(TEventDisplayModule,0)
};
}
#endif
