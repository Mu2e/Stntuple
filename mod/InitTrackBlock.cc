//-----------------------------------------------------------------------------
//  Dec 26 2000 P.Murat: initialization of the STNTUPLE track block
//  2014-06-23: remove vane support
//-----------------------------------------------------------------------------
#include <cstdio>
#include <algorithm>
#include "TROOT.h"
#include "TFolder.h"
#include "TLorentzVector.h"
#include "TVector2.h"
					// Mu2e
#include "Stntuple/obj/TStnTrack.hh"
#include "Stntuple/obj/TStnTrackBlock.hh"
#include "Stntuple/obj/TStnEvent.hh"
#include "Stntuple/obj/TStnHelix.hh"
#include "Stntuple/obj/TStnHelixBlock.hh"

#include "Stntuple/obj/TStnTrackSeed.hh"
#include "Stntuple/obj/TStnTrackSeedBlock.hh"

#include "messagefacility/MessageLogger/MessageLogger.h"

#include "art/Framework/Principal/Selector.h"
#include "art/Framework/Principal/Handle.h"

#include "Offline/GeometryService/inc/GeometryService.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/GeometryService/inc/VirtualDetector.hh"
#include "Offline/GeometryService/inc/DetectorSystem.hh"

#include "Offline/TrackerGeom/inc/Tracker.hh"
#include "Offline/CalorimeterGeom/inc/Calorimeter.hh"
#include "Offline/CalorimeterGeom/inc/DiskCalorimeter.hh"

#include "Offline/RecoDataProducts/inc/HelixSeed.hh"
#include "Offline/RecoDataProducts/inc/KalSeed.hh"
#include "Offline/RecoDataProducts/inc/KalSeedAssns.hh"

#include "Offline/BTrkData/inc/TrkStrawHit.hh"
#include "Offline/BTrkData/inc/TrkCaloHit.hh"
#include "Offline/BTrkData/inc/Doublet.hh"

#include "Offline/RecoDataProducts/inc/TrkCaloIntersect.hh"
#include "Offline/RecoDataProducts/inc/TrackClusterMatch.hh"

#include "Offline/MCDataProducts/inc/GenParticle.hh"
#include "Offline/MCDataProducts/inc/SimParticle.hh"
#include "Offline/MCDataProducts/inc/StrawGasStep.hh"
#include "Offline/DataProducts/inc/VirtualDetectorId.hh"
#include "Offline/DataProducts/inc/SurfaceId.hh"

#include "Offline/RecoDataProducts/inc/StrawDigi.hh"
#include "Offline/RecoDataProducts/inc/StrawHit.hh"
#include "Offline/RecoDataProducts/inc/CaloHit.hh"
#include "Offline/RecoDataProducts/inc/CaloCluster.hh"
#include "Offline/RecoDataProducts/inc/AlgorithmID.hh"

					          // BaBar
#include "BTrk/ProbTools/ChisqConsistency.hh"
#include "BTrk/BbrGeom/BbrVectorErr.hh"
#include "BTrk/BbrGeom/TrkLineTraj.hh"
#include "BTrk/TrkBase/TrkPoca.hh"
#include "BTrk/KalmanTrack/KalHit.hh"

#include "Stntuple/mod/InitTrackBlock.hh"
// #include "Stntuple/mod/THistModule.hh"
// #include "Stntuple/base/TNamedHandle.hh"

//-----------------------------------------------------------------------------
// Map[iplane][ipanel][il]: index of the layer in Z-ordered sequence
//-----------------------------------------------------------------------------
void StntupleInitTrackBlock::InitTrackerZMap(const mu2e::Tracker* Tracker, ZMap_t* Map) {
  int      ix, loc;
  double   z0, z1;

  int nplanes = Tracker->nPlanes();

  for (int ipl=0; ipl<nplanes; ipl++) {
    for (int isec=0; isec<6; isec++) {
      ix  = isec % 2;
      loc = 2*ipl+ix;
      Map->fMap[ipl][isec][0] = loc;
      Map->fMap[ipl][isec][1] = loc;
    }
    // form the list of Z-coordinates
    const mu2e::Straw *s0, *s1;

    s0 = &Tracker->getPlane(ipl).getPanel(0).getStraw(0);
    s1 = &Tracker->getPlane(ipl).getPanel(0).getStraw(1);
    z0 = s0->getMidPoint().z();
    z1 = s1->getMidPoint().z();

    Map->fZ[2*ipl] = (z0+z1)/2.;

    s0 = &Tracker->getPlane(ipl).getPanel(1).getStraw(0);
    s1 = &Tracker->getPlane(ipl).getPanel(1).getStraw(1);
    z0 = s0->getMidPoint().z();
    z1 = s1->getMidPoint().z();

    Map->fZ[2*ipl+1] = (z0+z1)/2.;
  }
}

//-----------------------------------------------------------------------------
// for a given Z find closest Z-layer, returns 'Plane'
// 'Offset' is a 'face number' within the Plane
//-----------------------------------------------------------------------------
void StntupleInitTrackBlock::get_station(const mu2e::Tracker* Tracker, ZMap_t* Map, double Z, int* Plane, int* Offset) {

  double dz, dz_min(1.e10);
  int    iface(-1);
  // looks that Device == Plane
  int nplanes = Tracker->nPlanes();
  // a plane has 2 "faces", 2 layers in each
  int nfaces  = 2*nplanes;

  for (int i=0; i<nfaces; i++) {
    dz = Map->fZ[i]-Z;
    if (fabs(dz) < dz_min) {
      iface  = i;
      dz_min = fabs(dz);
    }
  }

  *Plane   = iface / 2;
  *Offset  = iface % 2;
}

//-----------------------------------------------------------------------------
// extrapolate track to a given Z
//-----------------------------------------------------------------------------
double StntupleInitTrackBlock::s_at_given_z(const mu2e::KalSeed* KSeed, double Z) {

  printf("ERROR in StntupleInitTrackBlock:::s_at_given_z : REIMPLEMENT! \n");
  return -1;

  // double  ds(10.), s0, s1, s2, z0, z1, z2, dzds, sz, sz1, z01;

  // // s1     = Kffs->firstHit()->kalHit()->hit()->fltLen();
  // // s2     = Kffs->lastHit ()->kalHit()->hit()->fltLen();

  // const TrkHitVector* hots = &Kffs->hitVector();
  // int nh = hots->size();

  // const TrkHit *first(nullptr), *last(nullptr);

  // for (int ih=0; ih<nh; ++ih) {
  //   const TrkHit* hit = hots->at(ih);
  //   if (hit  != nullptr) {
  //     if (first == nullptr) first = hit;
  //     last = hit;
  //   }
  // }

  // s1 = first->fltLen();
  // s2 = last ->fltLen();

  // z1     = Kffs->position(s1).z();
  // z2     = Kffs->position(s2).z();

  // dzds   = (z2-z1)/(s2-s1);
  // //-----------------------------------------------------------------------------
  // // iterate once, choose the closest point
  // //-----------------------------------------------------------------------------
  // if (fabs(Z-z1) > fabs(Z-z2)) {
  //   z0 = z2;
  //   s0 = s2;
  // }
  // else {
  //   z0 = z1;
  //   s0 = s1;
  // }

  // sz    = s0+(Z-z0)/dzds;

  // z0     = Kffs->position(sz).z();     // z0 has to be close to Z(TT_FrontPA)
  // z01    = Kffs->position(sz+ds).z();

  // dzds   = (z01-z0)/ds;
  // sz1    = sz+(Z-z0)/dzds;	          // should be good enough

  // return sz1;
}


//-----------------------------------------------------------------------------
int StntupleInitTrackBlock::InitDataBlock(TStnDataBlock* Block, AbsEvent* AnEvent, Int_t Mode) {
  const char*               oname = {"InitMu2eTrackBlock"};
//-----------------------------------------------------------------------------
// cached pointers, owned by the StntupleMaker_module
//-----------------------------------------------------------------------------
  static int                          initialized(0);

  const int verbose(fVerbose); //control output level for debugging

  int                       ntrk(0), nassns(0), ev_number, rn_number;
  TStnTrack*                track;
  TStnTrackBlock            *data(0);

  ev_number = AnEvent->event();
  rn_number = AnEvent->run();

  if (Block->Initialized(ev_number,rn_number)) {
    if(verbose > 0) printf("%s::%s: Block initialized already, exiting\n", typeid(*this).name(), __func__);
    return 0;
  }

  mu2e::GeomHandle<mu2e::Tracker> ttHandle;
  tracker = ttHandle.get();

  data = (TStnTrackBlock*) Block;
  data->Clear();

  if (initialized == 0) {
    initialized = 1;

    InitTrackerZMap(tracker,&zmap);
  }

  list_of_algs = 0;
  art::Handle<mu2e::AlgorithmIDCollection> algsHandle;
  AnEvent->getByLabel(fAlgorithmIDCollTag, algsHandle);
  if (algsHandle.isValid()) list_of_algs = (mu2e::AlgorithmIDCollection*) algsHandle.product();

  list_of_kffs = 0;
  art::Handle<mu2e::KalSeedCollection> kffcH;
  AnEvent->getByLabel(fKFFCollTag,kffcH);
  if (kffcH.isValid())    {
    list_of_kffs = kffcH.product();
    ntrk         = list_of_kffs->size();
    if(verbose > 0) printf("%s::%s: KalSeedCollection %15s has %2i tracks\n", typeid(*this).name(), __func__, fKFFCollTag.encode().c_str(), ntrk);
  } else {
    if(verbose > 0) printf("%s::%s: KalSeedCollection %s not found!\n", typeid(*this).name(), __func__, fKFFCollTag.encode().c_str());
  }

  const mu2e::KalHelixAssns* list_of_kff_assns = nullptr;
  art::Handle<mu2e::KalHelixAssns> kffAssnsH;
  AnEvent->getByLabel(fKFFCollTag,kffAssnsH);
  if (kffAssnsH.isValid()) {
    list_of_kff_assns = kffAssnsH.product();
    nassns = list_of_kff_assns->size();
    if(verbose > 0) printf("%s::%s: KalHelixAssns %15s has %2i associations\n", typeid(*this).name(), __func__, fKFFCollTag.encode().c_str(), nassns);
    if(nassns != ntrk) printf("%s::%s: KalHelixAssns %15s has a different number of tracks! %i assns %i trks\n",
                              typeid(*this).name(), __func__, fKFFCollTag.encode().c_str(), nassns, ntrk);
  } else {
    if(verbose > 0) printf("%s::%s: KalHelixAssns %s not found!\n", typeid(*this).name(), __func__, fKFFCollTag.encode().c_str());
  }

  list_of_trk_qual = 0;
  art::Handle<mu2e::MVAResultCollection> trkQualHandle;
  AnEvent->getByLabel(fTrkQualCollTag,trkQualHandle);
  if (trkQualHandle.isValid()) list_of_trk_qual = trkQualHandle.product();
  else if(fTrkQualCollTag != "") printf(" InitTrackBlock::%s: Track quality collection %s not found\n", __func__, fTrkQualCollTag.encode().c_str());

  fSschColl = 0;
  art::Handle<mu2e::ComboHitCollection> sschcH;
  AnEvent->getByLabel(fSsChCollTag,sschcH);
  if (sschcH.isValid()) fSschColl = sschcH.product();
  else printf(" WARNING InitTrackBlock::%s: ComboHitCollection %s not found\n", __func__, fSsChCollTag.encode().c_str());

  list_of_mc_straw_hits = 0;
  art::Handle<mu2e::StrawDigiMCCollection> sdmcHandle;
  AnEvent->getByLabel(fStrawDigiMCCollTag,sdmcHandle);
  if (sdmcHandle.isValid()) list_of_mc_straw_hits = sdmcHandle.product();
  else if(fStrawDigiMCCollTag != "") printf(" InitTrackBlock::%s: StrawDigiMC collection %s not found\n", __func__, fStrawDigiMCCollTag.encode().c_str());

  list_of_extrapolated_tracks = 0;
  art::Handle<mu2e::TrkCaloIntersectCollection>  texHandle;
  AnEvent->getByLabel(fTciCollTag,texHandle);
  if (texHandle.isValid()) list_of_extrapolated_tracks = texHandle.product();

  art::Handle<mu2e::TrackClusterMatchCollection>  tcmH;
  AnEvent->getByLabel(fTcmCollTag,tcmH);

  list_of_pidp = 0;
  art::Handle<mu2e::PIDProductCollection>  pidpHandle;
  AnEvent->getByLabel(fPIDProductCollTag,pidpHandle);
  if (pidpHandle.isValid()) list_of_pidp = pidpHandle.product();

  art::ServiceHandle<mu2e::GeometryService>   geom;
  mu2e::GeomHandle<mu2e::DetectorSystem>      ds;
  mu2e::GeomHandle<mu2e::VirtualDetector>     vdet;

  const mu2e::AlgorithmID*  alg_id;
  int                       mask;

  const mu2e::Calorimeter* bc(NULL);

  if (geom->hasElement<mu2e::DiskCalorimeter>() ) {
    mu2e::GeomHandle<mu2e::DiskCalorimeter> h;
    bc = (const mu2e::Calorimeter*) h.get();
  }

  for (int itrk=0; itrk<ntrk; itrk++) {
    if(verbose > 1) printf("%s::%s: KalSeedCollection %s track %2i:\n", typeid(*this).name(), __func__, fKFFCollTag.encode().c_str(), itrk);
    track          = data->NewTrack();
    const mu2e::KalSeed* kffs = &list_of_kffs->at(itrk);
    if(!kffs) {
      printf("%s::%s: KalSeed %s track %2i not defined!\n", typeid(*this).name(), __func__, fKFFCollTag.encode().c_str(), itrk);
      continue;
    }

    // Attempt to find the corresponding helix
    const mu2e::HelixSeed* helix = nullptr;
    if(list_of_kff_assns) {
      for(int iassn = 0; iassn < nassns; ++iassn) {
        const auto assn = list_of_kff_assns->at(iassn);
        const mu2e::KalSeed*   itrack = &(*(assn.first));
        const mu2e::HelixSeed* ihelix = &(*(assn.second));
        if(&(*(itrack)) == &(*kffs)) {
          if(verbose > 1) printf(" --> Associated helix found, association index %i\n", iassn);
          helix = ihelix;
          break;
        }
      }
      if(!helix) printf("%s::%s: KalSeedCollection %s track %2i: Associated helix not found! N(track) = %i N(Assns) = %i\n",
                        typeid(*this).name(), __func__, fKFFCollTag.encode().c_str(), itrk, ntrk, nassns);
    }

//-----------------------------------------------------------------------------
// track-only-based particle ID, initialization ahs already happened in the constructor
//-----------------------------------------------------------------------------
    if (list_of_pidp) {
      const mu2e::PIDProduct* pidp = &list_of_pidp->at(itrk);
      track->fEleLogLHDeDx = pidp->GetLogEProb();
      track->fMuoLogLHDeDx = pidp->GetLogMProb();
      track->fRSlope       = pidp->GetResidualsSlope();
      track->fRSlopeErr    = pidp->GetResidualsSlopeError();
    }

    track->fKalRep[0] = (mu2e::KalSeed*) kffs;
    mask = (0x0001 << 16) | 0x0000;

    if (list_of_algs) {
      alg_id = &list_of_algs->at(itrk);
      mask   = alg_id->BestID() | (alg_id->AlgMask() << 16);
    }

    track->fAlgorithmID = mask;
//-----------------------------------------------------------------------------
// in all cases define momentum at lowest Z - ideally, at the tracker entrance
// 'entlen' - trajectory length, corresponding to the first point in Z (?)
//-----------------------------------------------------------------------------
    // double  h1_fltlen(1.e6), hn_fltlen(1.e6), sent, sexit;
    // const mu2e::TrkStrawHitSeed *first(nullptr), *last(nullptr);

    const std::vector<mu2e::TrkStrawHitSeed>* hots = &kffs->hits();
    int n_kffs_hits = hots->size();
    if(verbose > 1) printf("  N(hits) = %2i:\n", n_kffs_hits);

    // const TrkHit* first = kffs->firstHit()->kalHit()->hit();
    // const TrkHit* last  = kffs->lastHit ()->kalHit()->hit();

    // for (int ih=0; ih<n_kffs_hits; ++ih) {
    //   const mu2e::TrkStrawHitSeed* hit =  &hots->at(ih);
    //   if ((hit != nullptr) and (hit->flag().hasAnyProperty(mu2e::StrawHitFlagDetail::active))) {
    //     if (first == nullptr) first = hit;
    //     last = hit;
    //   }
    // }

    // if (first) h1_fltlen = first->trkLen();
    // if (last ) hn_fltlen = last->trkLen();

    // sent        = std::min(h1_fltlen,hn_fltlen);
    // sexit       = std::max(h1_fltlen,hn_fltlen);
//-----------------------------------------------------------------------------
// find segments corresponding to entry and exit points in the tracker
//-----------------------------------------------------------------------------
    // Intersections with the tracker surfaces
    const mu2e::KalIntersection *kinter_front(nullptr), *kinter_mid(nullptr), *kinter_back(nullptr);
    // Intersections with the stopping target surfaces
    const mu2e::KalIntersection *kinter_st_front(nullptr), *kinter_st_back(nullptr);
    std::vector<const mu2e::KalIntersection *> kinter_st_foils;

    // Take the first intersection (in time) for each surface (for the ST foil, take the highest Z foil)
    for(size_t ikinter = 0; ikinter < kffs->intersections().size(); ++ikinter) {
      auto const& kinter = kffs->intersections()[ikinter];
      if (kinter.surfaceId() == mu2e::SurfaceIdDetail::ST_Front) {
        if(!kinter_st_front || kinter_st_front->time() < kinter.time()) kinter_st_front = &kinter;
      }
      if (kinter.surfaceId() == mu2e::SurfaceIdDetail::ST_Back) {
        if(!kinter_st_back || kinter_st_back->time() < kinter.time()) kinter_st_back = &kinter;
      }
      if (kinter.surfaceId() == mu2e::SurfaceIdDetail::ST_Foils) { // save all the foils, ask about them later
        kinter_st_foils.push_back(&kinter);
      }
      if (kinter.surfaceId() == mu2e::SurfaceIdDetail::TT_Front) {
        if(!kinter_front || kinter_front->time() < kinter.time()) kinter_front = &kinter;
      }
      if (kinter.surfaceId() == mu2e::SurfaceIdDetail::TT_Mid) {
        if(!kinter_mid || kinter_mid->time() < kinter.time()) kinter_mid = &kinter;
      }
      if (kinter.surfaceId() == mu2e::SurfaceIdDetail::TT_Back) {
        if(!kinter_back || kinter_back->time() < kinter.time()) kinter_back = &kinter;
      }
      if(verbose > 2) printf("  Surface %10s: p = %4.1f t = %6.1f:\n", kinter.surfaceId().name().c_str(), kinter.mom(), kinter.time());
    }

    // Store the momentum at each surface if found
    if(kinter_st_front) { track->fPSTFront         = kinter_st_front->mom(); }
    if(kinter_st_back ) { track->fPSTBack          = kinter_st_back ->mom(); }
    if(kinter_front   ) { track->fPTrackerEntrance = kinter_front   ->mom(); }
    if(kinter_mid     ) { track->fPTrackerMiddle   = kinter_mid     ->mom(); }
    if(kinter_back    ) { track->fPTrackerExit     = kinter_back    ->mom(); }

    // Decide which intersection to use for the defaults, using the front if available (only Mid available for Online tracks)
    const mu2e::KalIntersection* kinter((kinter_front) ? kinter_front : (kinter_mid) ? kinter_mid : nullptr);
    KinKal::VEC3 fitmom = (kinter) ? kinter->momentum3() : KinKal::VEC3();
    KinKal::VEC3 pos    = (kinter) ? kinter->position3() : KinKal::VEC3();

    track->fX1 = pos.x();
    track->fY1 = pos.y();
    track->fZ1 = pos.z();

    double px, py, pz;
    px = fitmom.x();
    py = fitmom.y();
    pz = fitmom.z();
//-----------------------------------------------------------------------------
// track parameters in the first point
//-----------------------------------------------------------------------------
    track->Momentum()->SetXYZM(px,py,pz,0.511);
    track->fP         = (kinter) ? kinter->mom() : 0.f;
    track->fPt        = track->Momentum()->Pt();
    track->fChi2      = kffs->chisquared();
    track->fFitCons   = kffs->fitConsistency();
    if(kinter_mid) {
      track->fT0 = kinter_mid->time();
      track->fT0Err = std::sqrt(kinter_mid->loopHelix().paramVar(KinKal::LoopHelix::t0_));
    }
//-----------------------------------------------------------------------------
// momentum error in the first point
//-----------------------------------------------------------------------------
//    ROOT::Math::XYZVector  momdir = fitmom.unit();

    track->fFitMomErr = (kinter) ? kinter->momerr() : 0.f;
//-----------------------------------------------------------------------------
// determine, approximately, 'sz0' - flight length corresponding to the
// virtual detector at the tracker front
//-----------------------------------------------------------------------------
    // Hep3Vector tfront = ds->toDetector(vdet->getGlobal(mu2e::VirtualDetectorId::TT_FrontPA));
    // double     zfront = tfront.z();
    // double     sz0    = s_at_given_z(kffs,zfront);
//-----------------------------------------------------------------------------
// fP0 : track momentum value at Z(TT_FrontPA) - the same as in the first point
// fP2 : track momentum at Z(TT_Back), just for fun, should not be used for anything
//-----------------------------------------------------------------------------
    track->fP0        = track->fP; // can reuse , if needed
    track->fP2        = (kinter_back) ? kinter_back->mom() : 0.f;
//-----------------------------------------------------------------------------
// helical parameters at Z(TT_FrontPA)
//-----------------------------------------------------------------------------
    if(kinter) {
      try {
        KinKal::CentralHelix helx = kinter->centralHelix();
        track->fC0        = helx.omega(); // old
        track->fD0        = helx.d0();
        track->fZ0        = helx.z0();
        track->fPhi0      = helx.phi0();
        track->fTanDip    = helx.tanDip(); // old
        track->fCharge    = helx.charge();
      }
      catch (...) {
        mf::LogWarning(oname) << " ERROR line " << __LINE__ << ": KinKal::CentralHelix trouble" ;
        continue;
      }
    }
    if(verbose > 1) printf("  p = %5.1f, pT = %5.1f, t0 = %6.1f, d0 = %6.1f, z0 = %7.1f, phi0 = %5.2f, tdip = %4.2f\n",
                           track->fP, track->fPt, track->fT0, track->fD0, track->fZ0, track->fPhi0, track->fTanDip);

//-----------------------------------------------------------------------------
// virtual detector at the tracker exit: Time at Z(TT_Back)
//-----------------------------------------------------------------------------
//    Hep3Vector vd_tt_back = ds->toDetector(vdet->getGlobal(mu2e::VirtualDetectorId::TT_Back));
    // double     zback      = vd_tt_back.z();
    // double     szb        = s_at_given_z(kffs,zback);

    const float tback = -1.f; //FIXME
    // rename later
    track->fTBack         = tback;
//-----------------------------------------------------------------------------
// the total number of planes is 36, use 40 for simplicity
//-----------------------------------------------------------------------------
    const mu2e::TrkStrawHitSeed  *hit; // , *closest_hit(NULL);

    //const TrkHitVector*       kffs_hits = &kffs->hitVector();

    for (int j=0; j<40; j++) {
      track->fNHPerStation[j] = 0;
    }

    int     loc, nss_ch, found, ntrkhits(0), nhitsambig0(0); // , pdg_code;
    int     ipart;
    const static int max_npart(100);
    int     id(-1),  npart(0), part_nh[max_npart], part_id[max_npart];
    int     part_pdg_code[max_npart];
    double  part_first_z[max_npart], part_first_z_p[max_npart], part_first_z_pz[max_npart];
    double  part_last_z [max_npart], part_last_z_p [max_npart];
    int     nwrong = 0;
    double  mcdoca;

    const mu2e::SimParticle   *sim   (nullptr);
    const mu2e::StrawGasStep  *stgs  (nullptr);

    nss_ch = (fSschColl) ? fSschColl->size() : -1;

    if (nss_ch <= 0) {
      printf(">>> ERROR in InitTrackBlock::%s ComboHitCollection by module %s is empty, NHITS = %i\n", __func__, fSsChCollTag.encode().c_str(), nss_ch);
    }
    else {
      for (int it=0; it<n_kffs_hits; it++) {
        if(verbose > 5) printf(" Checking hit %i\n", it);
        hit = &hots->at(it);
//-----------------------------------------------------------------------------
// skip calorimeter hit
//-----------------------------------------------------------------------------
        if (! hit) continue;
        mu2e::StrawId sid = hit->strawId();
        const mu2e::Straw* straw = &tracker->straw(sid);
        ++ntrkhits;
//-----------------------------------------------------------------------------
// the rest makes sense only for active hits
// all KalSeed hits are "active", figuring out non-active ones
// now requires comparing the outputs of the seed fit and the full fit
//-----------------------------------------------------------------------------
        if (1) { // hit->isActive()) { // all KalSeed hits are active
          loc   = hit->index();
          if ((loc >= 0) && (loc < nss_ch)) {
            if ((list_of_mc_straw_hits) && (list_of_mc_straw_hits->size() > 0)) {
              if(verbose > 5) printf(" --> MC digi found\n");
              const mu2e::StrawDigiMC* mcdigi = &list_of_mc_straw_hits->at(loc);

              stgs = mcdigi->earlyStrawGasStep().get();
//-----------------------------------------------------------------------------
// count number of active hits with R > 200 um and misassigned drift signs
//-----------------------------------------------------------------------------
              if (stgs) {
                if(verbose > 5) printf(" --> MC gas step found\n");
                if (hit->driftRadius() > 0.2) {
                  const CLHEP::Hep3Vector* v1 = &straw->getMidPoint();
                  HepPoint p1(v1->x(),v1->y(),v1->z());

                  CLHEP::Hep3Vector v2 = stgs->position();
                  HepPoint    p2(v2.x(),v2.y(),v2.z());

                  TrkLineTraj trstraw(p1,straw->getDirection()  ,0.,0.);

                  TrkLineTraj trstep (p2,stgs->momvec().unit(),0.,0.);

                  TrkPoca poca(trstep, 0., trstraw, 0.);

                  mcdoca = poca.doca();
//-----------------------------------------------------------------------------
// if mcdoca and hit->_iamb have different signs, the hit drift direction has wrong sign
//-----------------------------------------------------------------------------
                  if (hit->ambig()*mcdoca < 0) nwrong      += 1;
                  if (hit->ambig()       == 0) nhitsambig0 += 1;
                }

                sim = &(*stgs->simParticle());
              }

              // check if there's an associated SIM particle, store the ID
              if (sim != NULL) id = sim->id().asInt();
              else {
                printf(">>> ERROR in %s : sim is NULL, set PDG_CODE to -1\n",oname);
                id = -1;
              }

              // check if this SIM particle has already been seen
              found = 0;
              for (int ip=0; ip<npart; ip++) {
                if (id == part_id[ip]) { //increment N(hits) by this SIM if already seen
                  found        = 1;
                  part_nh[ip] += 1;
                  const double dz = stgs->position().z();
                  if(dz < part_first_z[ip]) {
                    part_first_z   [ip] = dz;
                    part_first_z_p [ip] = std::sqrt(stgs->momentum().mag2());
                    part_first_z_pz[ip] = stgs->momentum().z();
                  }
                  if(dz > part_last_z[ip]) {
                    part_last_z    [ip] = dz;
                    part_last_z_p  [ip] = std::sqrt(stgs->momentum().mag2());
                  }
                  break;
                }
              }

              // if this is the first occurrence, add it to the list
              if (found == 0) {
                const double dz = stgs->position().z();
                part_id        [npart] = id;
                part_pdg_code  [npart] = sim->pdgId();
                part_nh        [npart] = 1;
                part_first_z   [npart] = dz;
                part_first_z_p [npart] = std::sqrt(stgs->momentum().mag2());
                part_first_z_pz[npart] = stgs->momentum().z();
                part_last_z    [npart] = dz;
                part_last_z_p  [npart] = std::sqrt(stgs->momentum().mag2());
                npart                 += 1;
              }
            }
          }
          else {
            printf(">>> ERROR in StntupleInitMu2eTrackBlock: wrong hit collection used");
            printf(", loc = %10i, n_straw_hits = %10i\n", loc,nss_ch);
          }

          const mu2e::StrawId& straw_id = straw->id();

          int ist = straw_id.getStation();

          track->fNHPerStation[ist] += 1;

          int pan = straw_id.getPanel();
          int lay = straw_id.getLayer();
          int bit = zmap.fMap[ist][pan][lay];

          track->fHitMask.SetBit(bit,1);
        }
      }
    }
//-----------------------------------------------------------------------------
// Dave's variables calculated by KalDiag
//-----------------------------------------------------------------------------
    // printf("InitTrackBlock: ERROR: kalDiag is gone, FIXIT\n");
    // _kalDiag->kalDiag(kffs,false);
//-----------------------------------------------------------------------------
// total number of hits associated with the trackand the number of bend sites
//-----------------------------------------------------------------------------
    track->fNHits     = ntrkhits; // ntrkhits | (_kalDiag->_trkinfo._nbend << 16);
    track->fNMatSites = 0; // _kalDiag->_trkinfo._nmat | (_kalDiag->_trkinfo._nmatactive << 16);

    if (list_of_trk_qual) track->fTrkQual = list_of_trk_qual->at(itrk)._value;
    else                  track->fTrkQual = -1.e6;
//-----------------------------------------------------------------------------
// defined bit-packed fNActive word
//-----------------------------------------------------------------------------
    track->fNActive   = kffs->nHits() | (nwrong << 16);
    if(verbose > 1) printf("  N(hits) = %2i, N(active) = %2i, N(wrong) = %2i, trkqual = %5.2f\n",
                           track->fNHits, track->NActive(), track->NWrong(),
                           track->fTrkQual);

    mu2e::Doublet*                     d;
    mu2e::DoubletAmbigResolver::Data_t r;

    int   nd, nd_tot(0), nd_os(0), nd_ss(0), ns;
    vector<mu2e::Doublet> list_of_doublets;

    //    _dar->findDoublets(kffs,&list_of_doublets);

    nd = list_of_doublets.size();
//-----------------------------------------------------------------------------
// counting only 2+ hit doublets
//-----------------------------------------------------------------------------
    int nad(0);  // number of doublets with all hits active

    for (int i=0; i<nd; i++) {
      d  = &list_of_doublets.at(i);
      ns = d->fNStrawHits;

      if (ns > 1) {
	nd_tot += 1;
	if (d->isSameSign()) nd_ss += 1;
	else                 nd_os += 1;

	int active = 1;
	for (int is=0; is<ns; is++) {
	  if (!d->fHit[is]->isActive()) {
	    active = 0;
	    break;
	  }
	}

	if (active == 1) {
	  nad += 1;
	}
      }
    }

    track->fNDoublets = nd_os | (nd_ss << 8) | (nhitsambig0 << 16) | (nad << 24);
//-----------------------------------------------------------------------------
// given track parameters, build the expected hit mask
//-----------------------------------------------------------------------------
    double z, zw, dz, dz_min; //  s0, closest_z(-1.e6), s;
    int    iplane, offset;
    int    nz(88);

    for (int iz=0; iz<nz; iz++) {
      z = zmap.fZ[iz];
					// find the track hit closest to that Z
      dz_min = 1.e10;
      for (auto it : *hots) {

	hit = &it;
	if (! hit) continue;

	// s_hit = &hit->comboHit();
	loc   = hit->index(); // s_hit-s_hit0;
	mu2e::StrawId sid = hit->strawId();
	const mu2e::Straw* straw = &tracker->straw(sid);
	zw    = straw->getMidPoint().z();
	dz    = z-zw;

	if (fabs(dz) < dz_min) {
	  dz_min      = fabs(dz);
	  // closest_hit = hit;
	  // closest_z   = zw;
	}
      }
//-----------------------------------------------------------------------------
// found closest hit and the extrapolation length, then extrapolate track
//-----------------------------------------------------------------------------
      // s0  = closest_hit->trkLen();
      //      s   = (z-track->fZ0)/(closest_z-track->fZ0)*s0;

      HepPoint    pz(0.,0.,0.); // FIXME      = kffs->position(s);

      get_station(tracker,&zmap,z,&iplane,&offset);

      const mu2e::Panel*  panel0 = NULL;
      const mu2e::Panel*  panel;
      const mu2e::Plane*  plane;
      double              min_dphi(1.e10);
      double              dphi, nx, ny, wx, wy, wrho, rho, phi0;
      CLHEP::Hep3Vector   w0mid;
      int                 ipanel;

      plane = &tracker->getPlane(iplane);

      for (int i=0; i<3; i++) {
	ipanel = 2*i+offset;		// panel number
					// check if point pz(x0,y0) overlaps with the segment iseg
					// expected mask is set to zero
	panel = &plane->getPanel(ipanel);
	w0mid = panel->straw0MidPoint();
					// also calculate pho wrt the sector
	wx    = w0mid.x();
	wy    = w0mid.y();

	phi0  = w0mid.phi();            // make sure I understand the range

	wrho  = sqrt(wx*wx+wy*wy);

	nx    = wx/wrho;
	ny    = wy/wrho;

	rho = nx*pz.x()+ny*pz.y();

	dphi = TVector2::Phi_mpi_pi(phi0-pz.phi());
	if ((abs(dphi) < min_dphi) && (rho > wrho)) {
	  min_dphi = fabs(dphi);
	  panel0   = panel;
	}
      }
//-----------------------------------------------------------------------------
// OK, closest segment found, set the expected bit..
//-----------------------------------------------------------------------------
      if (panel0 != NULL) {
	track->fExpectedHitMask.SetBit(iz,1);
      }
    }
//-----------------------------------------------------------------------------
// identify track with the particle which produced most hits
//-----------------------------------------------------------------------------
    ipart = 0;
    int nh0 = part_nh[0];

    for (int ip=1; ip<npart; ip++) {
      if (part_nh[ip] > nh0) {
	nh0 = part_nh[ip];
	ipart = ip;
      }
    }

    track->fPdgCode     = part_pdg_code[ipart];
    track->fPartID      = part_id      [ipart];
    track->fNGoodMcHits = (kffs->nDOF() << 16) + nh0;

    if(verbose > 1) printf(" Track %i MC info: ID = %3i, PDG = %5i, N(hits) = %2i, z(first) = %6.1f, p(first) = %5.1f, z(last) = %6.1f, p(last) = %5.1f\n",
                           itrk, track->fPartID, track->fPdgCode, nh0,
                           part_first_z[ipart], part_first_z_p[ipart],
                           part_last_z [ipart], part_last_z_p [ipart]);
//-----------------------------------------------------------------------------
// particle parameters at virtual detectors
//-----------------------------------------------------------------------------
    mu2e::GeomHandle<mu2e::VirtualDetector> vdg;

    double t_front(1.e6), t_stout(1.e6);

    // initialize MC track front to an early sim hit in case no virtual detector hit is found
    track->fPFront = part_first_z_p[ipart];
    track->fPStOut = -1.;
    track->fMcDirection = (part_first_z_pz[ipart] >= 0.) ? 1 : -1;

    if(verbose > 1) printf(" N(virtual detectors) = %i\n", vdg->nDet());
    if (vdg->nDet() > 0) {
//-----------------------------------------------------------------------------
// no more step point MC's in the tracker - straw gas steps there
//-----------------------------------------------------------------------------
      art::Handle<mu2e::StepPointMCCollection> vdhits;
      AnEvent->getByLabel(fVdhCollTag,vdhits);
      if (!vdhits.isValid()) {
	char warning[500];
	snprintf(warning,500,"WARNING in InitTrackBlock::%s: StepPointMCCollection %s not found\n",
                 __func__,fVdhCollTag.encode().data());
	mf::LogWarning(oname) << warning;
      }
      else {
	int nvdhits = vdhits->size();
        if(verbose > 2) printf(" N(virtual detector hits) = %i with tag %s\n", nvdhits, fVdhCollTag.encode().c_str());
	for (int i=0; i<nvdhits; i++) {
	  const mu2e::StepPointMC* hit = &(*vdhits)[i];

	  //int vdid = hit.volumeId();
	  mu2e::VirtualDetectorId vdid(hit->volumeId());
          if(verbose > 3) printf("  Virtual detector hit %2i: %15s\n", i, vdid.name().c_str());

	  if (vdid.id() == mu2e::VirtualDetectorId::ST_Out) {

	    //	    TAnaDump::Instance()->printStepPointMC(hit,"");

	    art::Ptr<mu2e::SimParticle> const& simptr = hit->simParticle();
	    const mu2e::SimParticle* sim  = simptr.operator ->();

	    if (sim == NULL) {
	      char warning[100];
	      printf(">>> ERROR: %s sim == NULL\n",oname);
	      sprintf(warning,"WARNING: SimParticle for step %i = NULL\n",i);
	      mf::LogWarning(oname) << warning;
	    } else {
              int sim_id = sim->id().asInt();
              if(verbose > 4) printf("   simID = %5i\n", sim_id);
              if ((sim_id == track->fPartID)  && (hit->time() <  t_stout)) {
                track->fPStOut = hit->momentum().mag();
                t_stout        = hit->time();
                if(verbose > 3) printf(" Sim ST_Out hit found: p = %5.1f, t = %6.1f\n", track->fPStOut, t_stout);
              }
            }
	  }
	  else if (vdid.isTrackerFront()) {

	    //	    TAnaDump::Instance()->printStepPointMC(hit,"");

	    art::Ptr<mu2e::SimParticle> const& simptr = hit->simParticle();
	    const mu2e::SimParticle* sim  = simptr.operator ->();

	    if (sim == NULL) {
	      printf(">>> ERROR: %s sim == NULL\n",oname);
	    } else {
              int sim_id = sim->id().asInt();
              if(verbose > 4) printf("   simID = %5i\n", sim_id);
              if ((sim_id == track->fPartID) && (hit->time() < t_front)) {
                track->fPFront = hit->momentum().mag();
                t_front        = hit->time();
                track->fMcDirection = (hit->momentum().z() >= 0.) ? 1 : -1;
                if(verbose > 3) printf(" Sim TT_Front hit found: p = %5.1f, t = %6.1f MC trajectory = %2i\n", track->fPFront, t_front, track->fMcDirection);
              }
            }
	  }
	}
      }
    }
    if(verbose > 1) printf(" Track MC P(front) = %5.1f, MC P(ST Out) = %5.1f, MC trajectory = %2i\n", track->fPFront, track->fPStOut, track->fMcDirection);

//-----------------------------------------------------------------------------
// number of MC hits produced by the mother particle
//-----------------------------------------------------------------------------
    track->fNMcStrawHits = 0;

    if (list_of_mc_straw_hits->size() > 0) {
      for (int i=0; i<nss_ch; i++) {
	const mu2e::StrawDigiMC* mcdigi = &list_of_mc_straw_hits->at(i);

        const mu2e::StrawGasStep* step = mcdigi->earlyStrawGasStep().get();

	if (step) {
	  art::Ptr<mu2e::SimParticle> const& simptr = step->simParticle();
	  art::Ptr<mu2e::SimParticle> mother        = simptr;
	  // while(mother->hasParent())  mother        = mother->parent();
	  const mu2e::SimParticle*    sim           = mother.get();

	  int sim_id = sim->id().asInt();

	  if (sim_id == track->fPartID) {
	    track->fNMcStrawHits += 1;
	  }
	}
      }
    }
//-----------------------------------------------------------------------------
// consider half-ready case when can't use the extrapolator yet; turn it off softly
//-----------------------------------------------------------------------------
//    const mu2e::TrkCaloIntersect*  extrk;
    CLHEP::Hep3Vector                     x1, x2;
    HepPoint                       extr_point;
    CLHEP::Hep3Vector                     extr_mom;
    // int                            iv, next;
    // TStnTrack::InterData_t*        vint;

    // if (list_of_extrapolated_tracks != NULL) next = list_of_extrapolated_tracks->size();
    // else                                     next = 0;

//     for (int iext=0; iext<next; iext++) {
//       extrk = &list_of_extrapolated_tracks->at(iext);
//       kffs  = extrk->trk().get();
//       if (kffs == track->fKalRep[0]) {
// 	if (track->fExtrk == 0) {
// 	  track->fExtrk = (const mu2e::TrkToCaloExtrapol*) extrk;
// 	}
// 	if (bc) {
// //-----------------------------------------------------------------------------
// // store coordinates of the best intersection in a plane
// //-----------------------------------------------------------------------------
// 	  iv   = extrk->diskId();
// 	  vint = &(track->fDisk[iv]);

// 	  if (vint->fID == -1) {
// 	    vint->fID        = iv;
// 	    vint->fExtrk     = (const mu2e::TrkToCaloExtrapol*) extrk;
// 	    vint->fChi2Match = 1.e6;
// 	  }
// 	  else {
// 	    printf("Run:Event: %i:%i %s : ADDITIONAL EXTR POINT for track %i on vane = %i\n",
// 		   rn_number,ev_number,oname,itrk,iv);
// 	  }
// 	}
//       }
//    }
//-----------------------------------------------------------------------------
// now loop over track-cluster matches and find the right ones to associate
// with the track
//-----------------------------------------------------------------------------
//    unsigned int nm (0);

    // const mu2e::TrackClusterMatchCollection* tcmcoll;

    // if (tcmH.isValid()) {
    //   tcmcoll = tcmH.product();
    //   nm      = tcmcoll->size();
    // }

    // const mu2e::TrackClusterMatch* tcm;

    // double best_chi2_match(1.e6);

//    for (size_t im=0; im<nm; im++) {
//      tcm   = &tcmcoll->at(im);
//      extrk = tcm->textrapol();
//      kffs  = extrk->trk().get();
//      if (kffs == track->fKalRep[0]) {
//	const mu2e::CaloCluster* cl = tcm->caloCluster();
//	iv   = cl->diskID();
//	vint = &track->fDisk[iv];
//	if (bc == 0) {
//	  printf(">>> InitTrackBlock ERROR: %s calorimeter is not defined \n",oname);
//	  continue;
//	}
//
//	x1   = bc->geomUtil().mu2eToDisk(iv,cl->cog3Vector());
//
//	if ((track->fClosestCaloCluster == NULL) || (tcm->chi2() < best_chi2_match )) {
////-----------------------------------------------------------------------------
//// if closest cluster has not been defined or the energy of the new one is higher
//// depending on the calorimeter geometry choice either DX or DZ is meaningful
////-----------------------------------------------------------------------------
//	  track->fClosestCaloCluster = cl;
//	  track->fExtrk              = (const mu2e::TrkToCaloExtrapol*) extrk;
//	  best_chi2_match            = tcm->chi2();
//	}
//
//	vint->fXTrk  = tcm->xtrk();
//	vint->fYTrk  = tcm->ytrk();
//	vint->fZTrk  = tcm->ztrk();
//	vint->fTime  = tcm->ttrk();
//
//	vint->fNxTrk = tcm->nx();
//	vint->fNyTrk = tcm->ny();
//	vint->fNzTrk = tcm->nz();
//
//	if (vint->fCluster == 0) {
//	  vint->fCluster      = cl;
//	  vint->fClusterIndex = tcm->icl();
//	  vint->fEnergy       = cl->energyDep();
//	  vint->fXCl          = x1.x();
//	  vint->fYCl          = x1.y();
//	  vint->fZCl          = x1.z();
//	  vint->fDt           = tcm->dt();
//	  vint->fDx           = tcm->dx();
//	  vint->fDy           = tcm->dy();
//	  vint->fDz           = tcm->dz();
//	  vint->fDu           = tcm->du();
//	  vint->fDv           = tcm->dv();
//	  vint->fChi2Match    = tcm->chi2();
//	  vint->fChi2Time     = tcm->chi2_time();
//	  vint->fIntDepth     = tcm->int_depth();
//	  vint->fPath         = tcm->ds();
//	  vint->fDr           = tcm->dr();
//	  vint->fSInt         = tcm->sint();
//	}
//	else {
//	  printf("%s : ADDITIONAL MATCH for track %i on vane = %i\n", oname,itrk,iv);
//	}
//      }
//    }
//-----------------------------------------------------------------------------
// find intersections to use for electron ID,
// in this version both VMinS and VMaxEp are the same
//-----------------------------------------------------------------------------
    double                    min_chi2_match(1.e6);
    TStnTrack::InterData_t*   v;

    track->fVMinS  = 0;
    track->fVMaxEp = 0;

    int ndisks = bc->nDisks();

    for (int iv=0; iv<ndisks; iv++) {
      v = &track->fDisk[iv];
      if (v->fID >= 0) {
	if (v->fCluster) {
	  if (v->fChi2Match < min_chi2_match) {
	    track->fVMaxEp = v;
	    track->fVMinS  = v;
	    min_chi2_match = v->fChi2Match;
	  }
	}
      }
    }
//-----------------------------------------------------------------------------
// define E/P by the first intersection, if it exists, the second in the
// high-occupancy environment may be unreliable
//----------------------------------------------------
    track->fClusterE = -track->fP;
    track->fDt       = -1.e12;
    track->fDx       = -1.e12;
    track->fDy       = -1.e12;
    track->fDz       = -1.e12;

    const mu2e::CaloCluster* calo_cluster = kffs->caloHit().caloCluster().get();
    if (calo_cluster) track->fClusterE = calo_cluster->energyDep();
    track->fEp = track->fClusterE/track->fP2;

    const bool use_calo_hit(true); // use calo-hit info vs. intersections

    if(!use_calo_hit) {
      if (track->fVMinS != 0) {
        if (track->fVMinS->fCluster) {
          track->fClusterE = track->fVMinS->fCluster->energyDep();
          track->fDx       = track->fVMinS->fDx;
          track->fDy       = track->fVMinS->fDy;
          track->fDz       = track->fVMinS->fDz;
          track->fDt       = track->fVMinS->fDt;
        }
        else {
          //-----------------------------------------------------------------------------
          // intersection with minimal S doesn't have a cluster, check MaxP
          //-----------------------------------------------------------------------------
          if (track->fVMaxEp != 0) {
            if (track->fVMaxEp->fCluster) {
              track->fClusterE = track->fVMaxEp->fCluster->energyDep();
              track->fDx       = track->fVMaxEp->fDx;
              track->fDy       = track->fVMaxEp->fDy;
              track->fDz       = track->fVMaxEp->fDz;
              track->fDt       = track->fVMaxEp->fDt;
            }
          }
        }
      }
    }

//--------------------------------------------------------------------------------
// now set the parameters associated with the TrkCaloHit
//--------------------------------------------------------------------------------
    const mu2e::TrkCaloHitSeed*  tch;
    TStnTrack::InterData_t*      vtch;

    tch = &kffs->caloHit();
    // for (auto it : hots) {
    //   tch   = itdynamic_cast<const mu2e::TrkCaloHit*> (*it);
    //-----------------------------------------------------------------------------
    // skip TrkStrawHit hit
    //-----------------------------------------------------------------------------
    if (tch) {
      vtch = &(track->fTrkCaloHit);
      const mu2e::CaloCluster* cl = tch->caloCluster().get();

      if (cl) {
	CLHEP::Hep3Vector cpos = bc->geomUtil().mu2eToTracker(bc->geomUtil().diskFFToMu2e( cl->diskID(), cl->cog3Vector()));

	CLHEP::Hep3Vector pos;
	// tch->hitPosition(pos);

	vtch->fID           = cl->diskID();		//
	vtch->fClusterIndex = -1;         // cluster index in the list of clusters

      // the following includes the (Calibrated) light-propagation time delay.  It should eventually be put in the reconstruction FIXME!
      // This velocity should come from conditions FIXME!

	vtch->fTime         = tch->t0().t0();          // extrapolated track time, not corrected by _dtOffset
	vtch->fEnergy       = cl->energyDep();            // cluster energy
	vtch->fXTrk         = pos.x();
	vtch->fYTrk         = pos.y();
	vtch->fZTrk         = pos.z();
	// vtch->fNxTrk        = -9999.;		// track direction cosines in the intersection point
	// vtch->fNyTrk        = -9999.;
	// vtch->fNzTrk        = -9999.;
	vtch->fXCl          = cpos.x();			// cluster coordinates
	vtch->fYCl          = cpos.y();
	vtch->fZCl          = cpos.z();
	vtch->fDx           = vtch->fXTrk - vtch->fXCl;	// TRK-CL
	vtch->fDy           = vtch->fYTrk - vtch->fYCl;	// TRK-CL
	vtch->fDz           = vtch->fZTrk - vtch->fZCl;	// TRK-CL
	vtch->fDt           = tch->t0().t0() - tch->time();
	// vtch->fDu           = -9999.;			// ** added in V6
	// vtch->fDv           = -9999.;			// ** added in V6
	// vtch->fChi2Match    = -9999.;		// track-cluster match chi&^2 (coord)
	// vtch->fChi2Time     = -9999.;		// track-cluster match chi&^2 (time)
	vtch->fPath         = tch->hitLen();			// track path in the disk
	vtch->fIntDepth     = -9999.;                     // ** added in V6 :assumed interaction depth
	vtch->fDr           = tch->_udoca; // clusterAxisDOCA(); // tch->poca().doca();         // distance of closest approach
	// vtch->fSInt         = -9999.;                 // ** added in V10: interaction length, calculated
	vtch->fCluster      = cl;
	//    vtch->fExtrk        = NULL;

        if(use_calo_hit) {
          track->fDx       = vtch->fDx;
          track->fDy       = vtch->fDy;
          track->fDz       = vtch->fDz;
          track->fDt       = vtch->fDt;
        }
      }
    }
  }
					// on return set event and run numbers
					// to mark block as initialized
  data->f_RunNumber   = rn_number;
  data->f_EventNumber = ev_number;

  return 0;
}

//-----------------------------------------------------------------------------
// 2015-04-02: this routine is not finished yet
//-----------------------------------------------------------------------------
Int_t StntupleInitTrackBlock::ResolveLinks(TStnDataBlock* Block, AbsEvent* AnEvent, int Mode)
{
  int    ev_number, rn_number;

  ev_number = AnEvent->event();
  rn_number = AnEvent->run();

  if (! Block->Initialized(ev_number,rn_number)) return -1;
//-----------------------------------------------------------------------------
// do not do initialize links for 2nd time
//-----------------------------------------------------------------------------
  if (Block->LinksInitialized()) return 0;

  TStnTrackBlock* tb = (TStnTrackBlock*) Block;

  art::Handle<mu2e::CaloClusterCollection> cch;
  if (not fCaloClusterCollTag.empty()) {
    AnEvent->getByLabel(fCaloClusterCollTag,cch);
  }
//-----------------------------------------------------------------------------
// TStnTrackBlock is already initialized
//-----------------------------------------------------------------------------
  TStnEvent* ev = Block->GetEvent();

  if (not fTrackTsCollTag.empty()) {
//-----------------------------------------------------------------------------
// seeds are stored, fill links part: 'TrackTs' collection stores, for each track,
// its KalSeed
//-----------------------------------------------------------------------------
    art::Handle<mu2e::KalSeedCollection>  ksch;
    mu2e::KalSeedCollection*              list_of_kalSeeds;

    AnEvent->getByLabel(fTrackTsCollTag,ksch);
    list_of_kalSeeds = (mu2e::KalSeedCollection*) ksch.product();

    TStnTrackSeedBlock* tsb = (TStnTrackSeedBlock*) ev->GetDataBlock(fTrackTsBlockName.Data());

    int    ntrk = tb->NTracks();
    int    nts  = tsb->NTrackSeeds();

    for (int i=0; i<ntrk; i++) {
      TStnTrack* trk = tb->Track(i);
      const mu2e::KalSeed *ts = &list_of_kalSeeds->at(i); // seed corresponding to track # i
      int  loc(-1);
      for (int j=0; j<nts; ++j) {
	const mu2e::KalSeed* tss = tsb->TrackSeed(j)->fTrackSeed;
	if (ts == tss) {
	  loc = j;
	  break;
	}
      }

      if (loc < 0) {
	printf(">>> ERROR: %s track %i -> no TrackSeed associated\n",
               fKFFCollTag.encode().data(), i);
	continue;
      }

      trk->SetTrackSeedIndex(loc);
    }
  }
//-----------------------------------------------------------------------------
// mark links as initialized
//-----------------------------------------------------------------------------
  Block->SetLinksInitialized();

  return 0;
}
