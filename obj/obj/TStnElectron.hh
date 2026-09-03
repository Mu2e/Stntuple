//
// Read the tracks added to the event by the track finding and fitting code.
//
// $Id: TStnElectron.hh,v 1.1 2014/06/13 06:14:48 murat Exp $
// $Author: murat $
// $Date: 2014/06/13 06:14:48 $
//
// Contact person Pavel Murat
//
#ifndef murat_inc_TStnElectron_hh
#define murat_inc_TStnElectron_hh

// C++ includes.
#include <iostream>

#include "TString.h"
#include "TFolder.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TBuffer.h"
					// 'KalRep' is a BaBar class
class TStnTrack;
class TStnCluster;

class TStnElectron : public TObject {
public:
  int                       fNumber;       // track index in the list of reconstructed tracks
  TStnTrack*                fTrack;
  TStnCluster*              fCluster;
  TLorentzVector            fMomentum;         // this assumes DELE fit hypothesis

  int                       fIDWord;	    // 
  int                       fITmp[5];
//-----------------------------------------------------------------------------
//  transient data members, all the persistent ones should go above
//-----------------------------------------------------------------------------
  Float_t                   fEp;
  Float_t                   fTmp[10];              // ! for temporary analysis needs
//-----------------------------------------------------------------------------
// methods
//-----------------------------------------------------------------------------
  TStnElectron(int i = -1);
  ~TStnElectron();

  TLorentzVector* Momentum() { return &fMomentum; }
  
  Int_t  GetMomentum  (TLorentzVector* Momentum);
  
  virtual void Print (Option_t* Opt = "") const ;

  ClassDefOverride(TStnElectron,0)
    
};

#endif
