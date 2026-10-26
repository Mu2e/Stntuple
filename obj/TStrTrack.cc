///////////////////////////////////////////////////////////////////////////////
// format : 1 (default) : track parameters only
//        : 
#include <iostream>
#include <format>

#include "Offline/RecoDataProducts/inc/TrkStrawHitSeed.hh"
#include "TBuffer.h"
#include "Stntuple/obj/TStrTrack.hh"

ClassImp(TStrTrack)

//-----------------------------------------------------------------------------
void TStrTrack::ReadV1(TBuffer& R__b) {
  struct TStrTrack_t_V1 {
    int     fNHits;
    int     fNDof ;
    int     fTcIndex;     // index of the corresponding time cluster
    int     fFormat;      // low byte: format
  
    float   fT0   ;
    float   fChi2 ;
    float   fX0   ;
    float   fY0   ;
    float   fZ0   ;
    float   fNx   ;
    float   fNy   ;
    float   fNz   ;
    const mu2e::KalSeed*  fOfflineKs; //!
  } data;
  
  int nwi = ((int*  ) &data.fT0       ) - &data.fNHits;
  int nwf = ((float*) &data.fOfflineKs) - &data.fT0;

  TObject::Streamer(R__b);
  R__b.ReadFastArray(&fNHits,nwi);
  R__b.ReadFastArray(&fT0   ,nwf);

//-----------------------------------------------------------------------------
// added in V2
//-----------------------------------------------------------------------------
    // fTrkStrawHitSeed; // written if fFormat & 0x2
    // fTrkStraw;        // written if fFormat & 0x4

}
//-----------------------------------------------------------------------------
void TStrTrack::Streamer(TBuffer& R__b) {
  
  int nwi = ((int*  ) &fT0       ) - &fNHits;
  int nwf = ((float*) &fOfflineKs) - &fT0;

  if (R__b.IsReading()) {
    Version_t R__v = R__b.ReadVersion();  // not used but want to look at 
    
    if      (R__v == 1) ReadV1(R__b);
    else if (R__v == 2) {
                                        // current version: V1
      TObject::Streamer(R__b);
      R__b.ReadFastArray(&fNHits,nwi);
      R__b.ReadFastArray(&fT0   ,nwf);

      if (fFormat & 0x2) {
        int n_hit_seeds;
        R__b >> n_hit_seeds;
        fTrkStrawHitSeed->resize(n_hit_seeds);
        for (int i=0; i<n_hit_seeds; i++) {
          mu2e::TrkStrawHitSeed* shs = &fTrkStrawHitSeed->at(i);
          uint16_t* w   = (uint16_t*) &shs->_sid;
          int*     flag = (int*     ) &shs->_flag;
          int* kkshflag = (int*     ) &shs->_kkshflag;
          R__b >> shs->_index >> *w;
          R__b >> *((char*) &shs->_eend._end);
          R__b >> *flag >> *kkshflag;
          R__b >> shs->_ambig >> shs->_algo >> shs->_frozen;

          int nwf1 =  &shs->_stime - &shs->_bkgqual + 1;
          R__b.ReadFastArray(&shs->_bkgqual,nwf1);
        }
      }
      if (fFormat & 0x4) {
// read back the straws
// straws need to be written out
        int n_straws;
        R__b >> n_straws;
        fTrkStraw->resize(n_straws);
        for (int i=0; i<n_straws; i++) {
          mu2e::TrkStraw* str = &fTrkStraw->at(i);
          uint16_t* sid = (uint16_t*) &str->_straw;
          int*     flag = (int*) &str->_flag;
          R__b >> *sid >> *flag;
          int nwf2 = &str->_phi-((float*) &str->_poca) + 1;
          R__b.ReadFastArray((float*) &str->_poca,nwf2);
        }
      }
    }
  }
  else {
    R__b.WriteVersion(TStrTrack::IsA());
    TObject::Streamer(R__b);
    R__b.WriteFastArray(&fNHits,nwi);
    R__b.WriteFastArray(&fT0   ,nwf);
    if (fFormat & 0x2) {
      // save hits 
      int n_hit_seeds = fTrkStrawHitSeed->size();
      R__b << n_hit_seeds;
      for (int i=0; i<n_hit_seeds; i++) {
        mu2e::TrkStrawHitSeed* shs = &fTrkStrawHitSeed->at(i);
        uint16_t* w        = (uint16_t*) &shs->_sid;
        int*      flag     = (int*) &shs->_flag;
        int*      kkshflag = (int*) &shs->_kkshflag;

        R__b << shs->_index << *w << shs->_eend._end << *flag << *kkshflag;
        R__b << shs->_ambig << shs->_algo << shs->_frozen;
        int nwf1 =  &shs->_stime-&shs->_bkgqual + 1;
        R__b.WriteFastArray(&shs->_bkgqual,nwf1);
      }
    }
    if (fFormat & 0x4) {
      // save straws
      int n_straws = fTrkStraw->size();
      R__b << n_straws;
      for (int i=0; i<n_straws; i++) {
        mu2e::TrkStraw* str = &fTrkStraw->at(i);
        uint16_t* sid = (uint16_t*) &str->_straw;
        int*     flag = (int*     ) &str->_flag;
        R__b << *sid << *flag;
        int nwf2 = &str->_phi-((float*) &str->_poca)+1;
        R__b.WriteFastArray((float*) &str->_poca,nwf2);
      }
    }
  }
}

//-----------------------------------------------------------------------------
TStrTrack::TStrTrack() :TObject() {
  fTrkStrawHitSeed = new std::vector<mu2e::TrkStrawHitSeed>;
  fTrkStraw        = new std::vector<mu2e::TrkStraw>;
  Clear();
}

//-----------------------------------------------------------------------------
TStrTrack::TStrTrack(int ID, int Format) : TObject() {
  SetUniqueID(ID);
  fTrkStrawHitSeed = new std::vector<mu2e::TrkStrawHitSeed>;
  fTrkStraw        = new std::vector<mu2e::TrkStraw>;
  Clear();
  // set format after clearing, Clear() sets it to 0
  fFormat = Format;
}

//-----------------------------------------------------------------------------
TStrTrack::~TStrTrack() {
}


//-----------------------------------------------------------------------------
void TStrTrack::Clear(const char* Opt) {
  fNHits  = -1;
  fNDof   = -1;
  fTcIndex = -1;
  fFormat  = 0;
  
  fT0     = -1;
  fChi2   = -1;
  fX0     = -1;
  fY0     = -1;
  fZ0     = -1;
  fNx     = -1;
  fNy     = -1;
  fNz     = -1;

  fTrkStrawHitSeed->clear();
  fTrkStraw->clear();
}

//-----------------------------------------------------------------------------
void TStrTrack::Print(const char* Opt) const {
  TString opt(Opt);

  opt.ToLower();
					// "non-const *this" for printing purposes
  // TStrTrack* t = (TStrTrack*) this;

  if ((opt == "") || (opt.Index("banner") >= 0)) {
//-----------------------------------------------------------------------------
// print banner
//-----------------------------------------------------------------------------
    printf("------------------------------------------------------------------------------------------------");
    printf("------------------------------------------------------------------\n");
    printf(" i   nh  ndof TcInd     chi2d     T0         X0         Y0        Z0          Nx         Ny       Nz\n");
    printf("------------------------------------------------------------------------------------------------");
    printf("------------------------------------------------------------------\n");
  }

  if ((opt == "") || (opt.Index("data") >= 0)) {
    printf("%2i  %3i  %4i  %4i  ",GetUniqueID(), fNHits, fNDof, fTcIndex);
    printf("%8.2f %9.2f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f",fChi2/fNDof,fT0,fX0,fY0,fZ0,fNx,fNy,fNz);
    printf("\n");
  }
}
