///////////////////////////////////////////////////////////////////////////////
//  always written out in split mode, so Streamer should never be called
///////////////////////////////////////////////////////////////////////////////

#ifdef __GNUG__
#pragma implementation
#endif

#include <iostream>
#include <iomanip>

#include "Stntuple/obj/TStnEvent.hh"
#include "Stntuple/obj/TStnHeaderBlock.hh"

ClassImp(TStnHeaderBlock)

//______________________________________________________________________________
void TStnHeaderBlock::ReadV1(TBuffer &R__b) {



  struct TStnHeaderBlockDataV1_t {
    Int_t             fVersion;
    Int_t             fEventNumber;
    Int_t             fRunNumber;
    Int_t             fSectionNumber;	// section number within the run
    Int_t             fMcFlag;		// MC flag, 0 for real data
    Int_t             fGoodRun;		// run flag
    Int_t             fBrCode;		// 
    Int_t             fGoodTrig;	// 
    Int_t             fTrigWord;	// 
    Int_t             fNTracks;         //
    Int_t             fCpu;             // packed word with processing time

    float             fInstLum;		// instantaneous luminosity

    TString           fStnVersion;      // like dev_243_16

    Int_t             fLastNumber;	// !
    Int_t             fLastRunNumber;	// ! 
  };

  TStnHeaderBlockDataV1_t data;

  int nwi     = ((Int_t*  )&data.fInstLum   )-&data.fVersion;
  int nwf     = ((Float_t*)&data.fLastNumber)-&data.fInstLum;

  R__b.ReadFastArray(&data.fVersion,nwi);
  R__b.ReadFastArray(&data.fInstLum,nwf);

  fVersion       = data.fVersion    ;
  fEventNumber   = data.fEventNumber;
  fRunNumber     = data.fRunNumber  ;
  fSectionNumber = data.fSectionNumber;
  fMcFlag        = data.fMcFlag     ;
  fGoodRun       = data.fGoodRun    ;
  fBrCode        = data.fBrCode     ;
  fGoodTrig      = data.fGoodTrig   ;
  fTrigWord      = data.fTrigWord   ;
  fNTracks       = data.fNTracks    ;
  fNStrawHits    = -1               ;   // *** added in V2
  fNCaloHits     = -1               ;   // *** added in V2
  fNCRVHits      = -1               ;   // *** added in V2
  fCpu           = data.fCpu        ;
  fNComboHits    = -1               ;   // *** added in V3

  fInstLum       = data.fInstLum    ;
  fMeanLum       = -1.              ;   // *** added in V2
  fEventWeight   =  1.              ;   // *** added in V4

  fStnVersion    = "";
}

//______________________________________________________________________________
void TStnHeaderBlock::Streamer(TBuffer &R__b)
{
   // Stream an object of class TStnHeaderBlock: should never be called

  int nwi     = ((Int_t*  )&fInstLum   )-&fVersion;
  int nwf     = ((Float_t*)&fStnVersion)-&fInstLum;

  if (R__b.IsReading()) {
    Version_t R__v = R__b.ReadVersion(); 
    if (R__v == 1) { 
      ReadV1(R__b);
    }
    else {
//-----------------------------------------------------------------------------
//  current version : 4
//-----------------------------------------------------------------------------
      R__b.ReadFastArray(&fVersion,nwi);
      R__b.ReadFastArray(&fInstLum,nwf);
      if(R__v < 4) {
        if(R__v < 3) {
          fNComboHits = -1;
        }
        fEventWeight = 1.f;
      }
    }
  } 
  else {
    R__b.WriteVersion(TStnHeaderBlock::IsA());
    R__b.WriteFastArray(&fVersion,nwi);
    R__b.WriteFastArray(&fInstLum,nwf);
  }
}

//------------------------------------------------------------------------------
TStnHeaderBlock::TStnHeaderBlock() {
  fEventNumber   = -1;
  fRunNumber     = -1;
  fNTracks       = -1;
  fNStrawHits    = -1;
  fNCaloHits     = -1;
  fNCRVHits      = -1;
  fNComboHits    = -1;

  fInstLum       = -1;
  fMeanLum       = -1;
  fEventWeight   = 1.;
  fStnVersion    = "";

  for (int i=0; i<kNFreeInts  ; i++) fInt  [i] = -1;
  for (int i=0; i<kNFreeFloats; i++) fFloat[i] = -1;

  fLastNumber    = -1;
  fLastRunNumber = -1;
}


//_____________________________________________________________________________
TStnHeaderBlock::~TStnHeaderBlock() {
}

//_____________________________________________________________________________
void TStnHeaderBlock::Clear(Option_t* opt) {
  f_EventNumber       = -1;
  f_RunNumber         = -1;
  f_SubrunNumber      = -1;
  fLinksInitialized   =  0;
}

//_____________________________________________________________________________
void TStnHeaderBlock::Print(Option_t* opt) const {
  // don't print header for the same event 2 times in a raw
  static char f_last_opt[1000] = {0};

  TStnHeaderBlock* block = (TStnHeaderBlock*) this;

  if ((fLastNumber            == fEventNumber) && 
      (fLastRunNumber         == fRunNumber  ) &&
      (strcmp(opt,f_last_opt) == 0           )    ) return;

  block->fLastNumber    = fEventNumber;
  block->fLastRunNumber = fRunNumber;
  strncpy(f_last_opt,opt,1000-1);

  printf("*** Run,Subrun,Event: %6i,%7i,%-9i: %s\n", 
	 fRunNumber,fSectionNumber,fEventNumber,opt);
}
