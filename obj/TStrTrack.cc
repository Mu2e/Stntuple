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
  
  int nwi = ((int*  ) &fT0      ) - &fNHits;
  int nwf = ((float*) &fHitSeeds) - &fT0;

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
        fHitSeeds->resize(n_hit_seeds);
        for (int i=0; i<n_hit_seeds; i++) {
          mu2e::TrkStrawHitSeed* shs = &fHitSeeds->at(i);
          uint16_t* w        = (uint16_t*) &shs->_sid;
          uint32_t* flag     = (uint32_t*) &shs->_flag;
          uint32_t* kkshflag = (uint32_t*) &shs->_kkshflag;

          R__b.ReadUShort(shs->_index);
          R__b.ReadUShort(*w);
          R__b.ReadUChar(*((u_char*)&shs->_eend._end));
          R__b.ReadUInt(*flag);
          R__b.ReadUInt(*kkshflag);
          R__b.ReadInt(shs->_ambig);
          R__b.ReadInt(shs->_algo);
          R__b.ReadBool(shs->_frozen);

          int nwf1 =  &shs->_stime - &shs->_bkgqual + 1;
          R__b.ReadFastArray(&shs->_bkgqual,nwf1);
        }
      }
      if (fFormat & 0x4) {
// read back the straws
// straws need to be written out
        int n_straws;
        R__b.ReadInt(n_straws);
        fStraws->resize(n_straws);
        for (int i=0; i<n_straws; i++) {
          mu2e::TrkStraw* str = &fStraws->at(i);
          uint16_t*   sid = (uint16_t*) &str->_straw;
          uint32_t*  flag = (uint32_t*) &str->_flag;
          R__b.ReadUShort(*sid);
          R__b.ReadUInt(*flag);
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
      int n_hit_seeds = fHitSeeds->size();
      R__b << n_hit_seeds;
      for (int i=0; i<n_hit_seeds; i++) {
        mu2e::TrkStrawHitSeed* shs = &fHitSeeds->at(i);
        uint16_t* w        = (uint16_t*) &shs->_sid;
        uint32_t* flag     = (uint32_t*) &shs->_flag;
        uint32_t* kkshflag = (uint32_t*) &shs->_kkshflag;

        R__b.WriteUShort(shs->_index);
        R__b.WriteUShort(*w);
        R__b.WriteUChar(shs->_eend._end);
        R__b.WriteUInt(*flag);
        R__b.WriteUInt(*kkshflag);
        R__b.WriteInt(shs->_ambig);
        R__b.WriteInt(shs->_algo);
        R__b.WriteBool(shs->_frozen);
        int nwf1 =  &shs->_stime-&shs->_bkgqual + 1;
        R__b.WriteFastArray(&shs->_bkgqual,nwf1);
      }
    }
    if (fFormat & 0x4) {
      // save straws
      int n_straws = fStraws->size();
      R__b.WriteInt(n_straws);
      for (int i=0; i<n_straws; i++) {
        mu2e::TrkStraw* str = &fStraws->at(i);
        uint16_t* sid  = (uint16_t*) &str->_straw;
        uint32_t* flag = (uint32_t*) &str->_flag;
        R__b.WriteUShort(*sid);
        R__b.WriteUInt(*flag);
        int nwf2 = &str->_phi-((float*) &str->_poca)+1;
        R__b.WriteFastArray((float*) &str->_poca,nwf2);
      }
    }
  }
}

//-----------------------------------------------------------------------------
TStrTrack::TStrTrack() :TObject() {
  fHitSeeds = new std::vector<mu2e::TrkStrawHitSeed>;
  fStraws   = new std::vector<mu2e::TrkStraw>;
  Clear();
}

//-----------------------------------------------------------------------------
TStrTrack::TStrTrack(int ID, int Format) : TObject() {
  SetUniqueID(ID);
  fHitSeeds = new std::vector<mu2e::TrkStrawHitSeed>;
  fStraws   = new std::vector<mu2e::TrkStraw>;
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

  fHitSeeds->clear();
  fStraws->clear();
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

    if (fFormat & 0x2) {
      // print hits
      int nh = fHitSeeds->size();
      for (int i=0; i<nh; i++) {
        mu2e::TrkStrawHitSeed* shs = &(*fHitSeeds)[i];
        std::cout << std::format("{:3d} {:4d} {:#04x} {} {:9.3f}\n",
                                 i,shs->index(),shs->_sid.asUint16(),shs->ambig(),shs->time());
      }
    }

    if (fFormat & 0x4) {
      // print hits
      std::cout << std::format(" I   SID   DOCA    X     Y    Z\n");
      int ns = fStraws->size();
      for (int i=0; i<ns; i++) {
        mu2e::TrkStraw* ts = &(*fStraws)[i];
        std::cout << std::format("{:3d} {:#04x} {:9.3f} {:9.3f} {:9.3f} {:9.3f}\n",
                                 i,ts->_straw.asUint16(),
                                 ts->_doca,
                                 ts->_poca.X(),ts->_poca.Y(),ts->_poca.Z());
      }
    }

    
  }
}
