///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
#ifndef __InitStrawDigiBlock__
#define __InitStrawDigiBlock__

#include <string.h>

#include "canvas/Utilities/InputTag.h"

#include "Offline/TrackerConditions/inc/TrackerPanelMap.hh"
#include "Offline/ProditionsService/inc/ProditionsHandle.hh"

#include "Stntuple/obj/TStnInitDataBlock.hh"

class StntupleInitStrawDigiBlock : public TStnInitDataBlock {
public:
  
  struct WfParam_t {
    int   fs;                           // first sample above _minPulseHeight
    float bl;                           // baseline
    float ph;                           // pulse height
    float q;                            // Q(positive)
    float qt;                           // Q(tail)
  };

  art::InputTag   fChCollTag;           // 
  art::InputTag   fSdCollTag;           // 
  art::InputTag   fShCollTag;           // straw hits are made out of straw digis,
  art::InputTag   fTcCollTag;           // time cluster collection

  int             fLastRun;
  int             fNAdcSamples;
  int             fNSamplesBL;
  float           fMinEDep;
  float           fMinDt;

  bool            fStoreCloseHitsOnly;

  float           fMinPulseHeight;
  float           fMinSDPHToSave;
  
  mu2e::ProditionsHandle<mu2e::TrackerPanelMap>   fTpm_h;
  const mu2e::TrackerPanelMap*                    fTrkPanelMap;
//-----------------------------------------------------------------------------
// functions
//-----------------------------------------------------------------------------
public:

  StntupleInitStrawDigiBlock();
  
  void        SetChCollTag(art::InputTag& Tag) { fChCollTag = Tag; }
  void        SetSdCollTag(art::InputTag& Tag) { fSdCollTag = Tag; }
  void        SetShCollTag(art::InputTag& Tag) { fShCollTag = Tag; }
  void        SetTcCollTag(art::InputTag& Tag) { fTcCollTag = Tag; }
  //   void   SetStrawDigiMCCollTag (art::InputTag& Tag) { fSdmcCollTag = Tag; }

  void        SetStoreCloseHitsOnly (bool  Flag) { fStoreCloseHitsOnly = Flag; }
  void        SetMinEDep            (float EDep) { fMinEDep            = EDep; }
  void        SetMinDt              (float Dt  ) { fMinDt              = Dt;   }

  int         ProcessAdcWaveform (float* Wf, WfParam_t* Wp);
  
  virtual int InitDataBlock      (TStnDataBlock* Block, AbsEvent* Evt, int Mode);
  virtual int ResolveLinks       (TStnDataBlock* Block, AbsEvent* Evt, int Mode);

};

#endif
