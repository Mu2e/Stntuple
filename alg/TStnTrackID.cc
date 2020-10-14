///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
#include <cmath>

#include "TH1F.h"
#include "Stntuple/alg/TStnTrackID.hh"
#include "Stntuple/obj/TStnTrack.hh"

ClassImp(TStnTrackID)

//-----------------------------------------------------------------------------
// everything is commented out to avoid compilation warnings
//-----------------------------------------------------------------------------
void TStnTrackID::Streamer(TBuffer &R__b) {
// //  // Stream an object of class TStnTrackID.
// //
//  UInt_t R__s, R__c;
// //
// //  int nwi = ((Int_t*  )&fMinFitCons)-&fUseMask;
// //  int nwf = ((Float_t*)&fEOR       )-&fMinFitCons  ;
// //
  if (R__b.IsReading()) {
//       Version_t R__v = R__b.ReadVersion(&R__s, &R__c); 
// //    TNamed::Streamer(R__b);
// //    R__b.ReadFastArray(&fUseMask   ,nwi);
// //    R__b.ReadFastArray(&fMinFitCons,nwf);
// //    R__b.CheckByteCount(R__s, R__c, TStnTrackID::IsA());
  } 
  else {
//-----------------------------------------------------------------------------
//  write section
//-----------------------------------------------------------------------------
//    R__c = R__b.WriteVersion(TStnTrackID::IsA(), kTRUE);
// //    TNamed::Streamer(R__b);
// //    R__b.WriteFastArray(&fUseMask   ,nwi);
// //    R__b.WriteFastArray(&fMinFitCons,nwf);
// //    R__b.SetByteCount(R__c, kTRUE);
  }
}

//-----------------------------------------------------------------------------
// default: "Set C" cuts , as defined in 
// http://mu2e-docdb.fnal.gov:8080/cgi-bin/RetrieveFile?docid=2447;filename=Resolutions.pdf;version=2
//
// 2014-04-02: update according to 
//    http://mu2e-docdb.fnal.gov:8080/cgi-bin/RetrieveFile?docid=3996;filename=cutset_c.pdf;version=2
//-----------------------------------------------------------------------------
TStnTrackID::TStnTrackID(const char* Name): TNamed(Name,Name) {
  
  fUseMask         = 0xffffffff;
  fMinFitCons      = 2.e-3;
  fMinT0           = 700.;		// ns
  fMinNActive      = 25;
  fMaxNActive      = 200;
  fMaxT0Err        = 0.9;  		// in ns
  fMaxMomErr       = 0.25;  		// in MeV
  fMinTanDip       = tan(M_PI/6.);	// 0.5773
  fMaxTanDip       = 1.0;  
  fMinD0           = -80.;		// in mm
  fMaxD0           = 105.;

  fMinRMax         = 450;
  fMaxRMax         = 680;

  fMaxChi2Dof      = 1.e6;       // no chi2 cut by default
  fMaxDNa          = 1000;       // no DNa  cut by default
  fMinFNa          = -1.;        // no FNa  cut by default

  fMinTrkQual      = -100.;     // if undefined , = 100
				// initialize spare words to zero

  for (int i=0; i< kNFreeInts  ; i++) fInteger[i] = 0;
  for (int i=0; i< kNFreeFloats; i++) fFloat  [i] = 0;
}


//_____________________________________________________________________________
TStnTrackID::~TStnTrackID() {
}

//-----------------------------------------------------------------------------
//  need to implement the D0 cuts
//-----------------------------------------------------------------------------
int TStnTrackID::IDWord(TStnTrack* Track) {

  int              id_word(0), nactive, dna;
  double           fcons, fna, chi2dof, d0, t0, t0_err, mom_err, tan_dip, rmax, trk_qual;
  
  fcons          = Track->FitCons();
  chi2dof        = Track->Chi2Dof();
  t0             = Track->fT0;
  t0_err         = Track->fT0Err;
  nactive        = Track->NActive();
  dna            = Track->NHits()-nactive;
  fna            = double(nactive)/Track->NHits();
  tan_dip        = Track->fTanDip;
  mom_err        = Track->fFitMomErr;
  d0             = Track->fD0;              // signed impact parameter, convention: mu2e-781
  rmax           = fabs(d0+2./Track->C0()); // 2/c0 - signed diameter
  trk_qual       = Track->DaveTrkQual();

  if (fcons            <  fMinFitCons  ) id_word |= kFitConsBit ;
  if (chi2dof          >  fMaxChi2Dof  ) id_word |= kChi2DofBit ;
  if (t0               <  fMinT0       ) id_word |= kT0Bit      ;
  if (nactive          <  fMinNActive  ) id_word |= kNActiveBit ;
  if (nactive          >= fMaxNActive  ) id_word |= kNActiveBit ;
  if (dna              >= fMaxDNa      ) id_word |= kDNaBit     ;
  if (fna              <  fMinFNa      ) id_word |= kFNaBit     ;

  if (t0_err           >  fMaxT0Err    ) id_word |= kT0ErrBit   ;
  if (mom_err          >  fMaxMomErr   ) id_word |= kMomErrBit  ;

  if (fabs(tan_dip)    <  fMinTanDip   ) id_word |= kTanDipBit  ;
  if (fabs(tan_dip)    >  fMaxTanDip   ) id_word |= kTanDipBit  ;

  if (d0               <  fMinD0       ) id_word |= kD0Bit      ;
  if (d0               >  fMaxD0       ) id_word |= kD0Bit      ;

  if (rmax             <  fMinRMax     ) id_word |= kRMaxBit    ;
  if (rmax             >  fMaxRMax     ) id_word |= kRMaxBit    ;
  if (trk_qual         <  fMinTrkQual  ) id_word |= kTrkQualBit ;

  return  (id_word & fUseMask);
}

//_____________________________________________________________________________
int TStnTrackID::TightIDWord(TStnTrack* Track) {
  int  id_word;
  id_word = IDWord(Track);
  return id_word;
}


//-----------------------------------------------------------------------------
// so far, no use
//-----------------------------------------------------------------------------
int TStnTrackID::LooseIDWord(TStnTrack* Track) {
  return IDWord(Track);
}


//_____________________________________________________________________________
void TStnTrackID::FillHistograms(Hist_t* Hist, TStnTrack* Track, Int_t Mode) {
//-----------------------------------------------------------------------------
// distributions for ID variables (with n-1 ID cuts applied). 
// Plot 3 sets of histograms:
// --------------------------
// set 1: for all the muons
// set 2: for the parameter given all the rest cuts were successfull
// set 3&4: for the cuts in sequence they are applied
// Mode = 1: tight ID (default) 
//      = 2: loose ID
//-----------------------------------------------------------------------------
  int     id_word(0xffffffff);
 
  if      (Mode == 1) id_word = IDWord     (Track);
  else if (Mode == 2) id_word = LooseIDWord(Track);

  float nactive = Track->NActive();
  float dna     = Track->NHits()-Track->NActive();
  float fcons   = Track->FitCons();
  float chi2dof = Track->Chi2Dof();
  float t0      = Track->T0();
  float t0err   = Track->T0Err();
  //  float mom     = Track->P();
  float momerr  = Track->FitMomErr();
  float tandip  = Track->TanDip();
  float d0      = Track->D0();
  float rmax    = Track->RMax();
  float trkqual = Track->DaveTrkQual();
//-----------------------------------------------------------------------------
//  1. number of points
//-----------------------------------------------------------------------------
  Hist->fNActive[0]->Fill(nactive);
  if ((id_word & ~kNActiveBit) == 0) Hist->fNActive[1]->Fill(nactive);
  if (id_word == 0) Hist->fNActive[4]->Fill(nactive);

  Hist->fFitCons[0]->Fill(fcons);
  if ((id_word & ~kFitConsBit) == 0) Hist->fFitCons[1]->Fill(fcons);
  if (id_word == 0) Hist->fFitCons[4]->Fill(fcons);

  Hist->fChi2Dof[0]->Fill(chi2dof);
  if ((id_word & ~kChi2DofBit) == 0) Hist->fChi2Dof[1]->Fill(chi2dof);
  if (id_word == 0) Hist->fChi2Dof[4]->Fill(chi2dof);

  Hist->fT0[0]->Fill(t0);
  if ((id_word & ~kT0Bit) == 0) Hist->fT0[1]->Fill(t0);
  if (id_word == 0) Hist->fT0[4]->Fill(t0);

  Hist->fT0Err[0]->Fill(t0err);
  if ((id_word & ~kT0ErrBit) == 0) Hist->fT0Err[1]->Fill(t0err);
  if (id_word == 0) Hist->fT0Err[4]->Fill(t0err);

  Hist->fMomErr[0]->Fill(momerr);
  if ((id_word & ~kMomErrBit) == 0) Hist->fMomErr[1]->Fill(momerr);
  if (id_word == 0) Hist->fMomErr[4]->Fill(momerr);

  Hist->fDNa[0]->Fill(dna);
  if ((id_word & ~kDNaBit) == 0) Hist->fDNa[1]->Fill(dna);
  if (id_word == 0) Hist->fDNa[4]->Fill(dna);

  Hist->fTanDip[0]->Fill(tandip);
  if ((id_word & ~kTanDipBit) == 0) Hist->fTanDip[1]->Fill(tandip);
  if (id_word == 0) Hist->fTanDip[4]->Fill(tandip);

  Hist->fD0[0]->Fill(d0);
  if ((id_word & ~kD0Bit) == 0) Hist->fD0[1]->Fill(d0);
  if (id_word == 0) Hist->fD0[4]->Fill(d0);

  Hist->fRMax[0]->Fill(rmax);
  if ((id_word & ~kRMaxBit) == 0) Hist->fRMax[1]->Fill(rmax);
  if (id_word == 0) Hist->fRMax[4]->Fill(rmax);

  Hist->fTrkQual[0]->Fill(trkqual);
  if ((id_word & ~kRMaxBit) == 0) Hist->fTrkQual[1]->Fill(trkqual);
  if (id_word == 0) Hist->fTrkQual[4]->Fill(trkqual);

//-----------------------------------------------------------------------------
//  single histogram showing how often every particular cut failed
//-----------------------------------------------------------------------------
  for (int bit=0; bit<32; bit++) {
    if (((id_word >> bit) & 0x1) == 1) {
      Hist->fFailedBits->Fill(bit);
    }
  }
  Hist->fPassed->Fill(id_word == 0);
//-----------------------------------------------------------------------------
//  ***** now histogram effect of cuts in the order they are applied
//-----------------------------------------------------------------------------
  Hist->fNActive[2]->Fill(nactive);
  if ((id_word & kNActiveBit) != 0)                         goto END;
  Hist->fNActive[3]->Fill(nactive);
  
  Hist->fFitCons[2]->Fill(fcons);
  if ((id_word & kFitConsBit) != 0)                         goto END;
  Hist->fFitCons[3]->Fill(fcons);
  
  Hist->fChi2Dof[2]->Fill(chi2dof);
  if ((id_word & kChi2DofBit) != 0)                         goto END;
  Hist->fChi2Dof[3]->Fill(chi2dof);
  
  Hist->fT0[2]->Fill(t0);
  if ((id_word & kT0Bit) != 0)                              goto END;
  Hist->fT0[3]->Fill(t0);
  
  Hist->fT0Err[2]->Fill(t0err);
  if ((id_word & kT0ErrBit) != 0)                           goto END;
  Hist->fT0Err[3]->Fill(t0err);
  
  Hist->fMomErr[2]->Fill(momerr);
  if ((id_word & kMomErrBit) != 0)                          goto END;
  Hist->fMomErr[3]->Fill(momerr);
  
  Hist->fDNa[2]->Fill(dna);
  if ((id_word & kDNaBit) != 0)                             goto END;
  Hist->fDNa[3]->Fill(dna);
  
  Hist->fTanDip[2]->Fill(tandip);
  if ((id_word & kTanDipBit) != 0)                          goto END;
  Hist->fTanDip[3]->Fill(tandip);
  
  Hist->fD0[2]->Fill(d0);
  if ((id_word & kD0Bit) != 0)                             goto END;
  Hist->fD0[3]->Fill(d0);
  
  Hist->fRMax[2]->Fill(rmax);
  if ((id_word & kRMaxBit) != 0)                           goto END;
  Hist->fRMax[3]->Fill(rmax);
  
  Hist->fTrkQual[2]->Fill(trkqual);
  if ((id_word & kTrkQualBit) != 0)                        goto END;
  Hist->fTrkQual[3]->Fill(trkqual);
  
  
 END:;
}


//_____________________________________________________________________________
void TStnTrackID::Print(const char* Opt) const {
  printf("-----------------------------------------------------\n");
  printf("      track ID cuts                                  \n");
  printf("-----------------------------------------------------\n");
  printf(" bit  0: fMinFitCons   = %12.4f\n",fMinFitCons   );
  printf(" bit  1: fMinT0        = %12.4f\n",fMinT0        );
  printf(" bit  2: fNMinActive   = %12i  fNMaxNActive = %12i\n",fMinNActive,fMaxNActive);
  printf(" bit  3: fMaxT0Err     = %12.4f\n",fMaxT0Err     );
  printf(" bit  4: fMaxMomErr    = %12.4f\n",fMaxMomErr );
  printf(" bit  5: fTanDip       = %12.4f < tan(dip)   < %12.4f\n",fMinTanDip,fMaxTanDip);
  printf(" bit  6: fD0           = %12.4f < D0         < %12.4f\n",fMinD0    ,fMaxD0    );
  printf(" bit  7: fRMax         = %12.4f < D0+2/omega < %12.4f\n",fMinRMax  ,fMaxRMax  );
}


