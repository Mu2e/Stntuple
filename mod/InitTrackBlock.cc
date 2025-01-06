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
  
  int                       ntrk(0), ev_number, rn_number;
  TStnTrack*                track;
  TStnTrackBlock            *data(0);   

  ev_number = AnEvent->event();
  rn_number = AnEvent->run();

  if (Block->Initialized(ev_number,rn_number)) return 0;

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
  }

  list_of_trk_qual = 0;
  art::Handle<mu2e::TrkQualCollection> trkQualHandle;
  AnEvent->getByLabel(fTrkQualCollTag,trkQualHandle);
  if (trkQualHandle.isValid()) list_of_trk_qual = trkQualHandle.product();

  fSschColl = 0;
  art::Handle<mu2e::ComboHitCollection> sschcH;
  AnEvent->getByLabel(fSsChCollTag,sschcH);
  if (sschcH.isValid()) fSschColl = sschcH.product();

  list_of_mc_straw_hits = 0;
  art::Handle<mu2e::StrawDigiMCCollection> sdmcHandle;
  AnEvent->getByLabel(fStrawDigiMCCollTag,sdmcHandle);
  if (sdmcHandle.isValid()) list_of_mc_straw_hits = sdmcHandle.product();

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
    track          = data->NewTrack();
    const mu2e::KalSeed* kffs = &list_of_kffs->at(itrk);
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
    const mu2e::KalSegment *kseg(nullptr), *kseg_exit(nullptr);

    double zmin(1.e6), zmax(-1.e6);

    for(auto const& ks : kffs->segments() ) {
      double z = ks.position3().z();
      if (z < zmin) {
	kseg = &ks;
	zmin = z  ;  
      }

      if (z > zmax) {
	kseg_exit = &ks;
	zmax      = z;
      }
    }

    for(size_t ikinter = 0; ikinter < kffs->intersections().size(); ++ikinter) {
      auto const& kinter = kffs->intersections()[ikinter];
      if (kinter.surfaceId() == mu2e::SurfaceIdDetail::TT_Front) { track->fPTrackerEntrance = kinter.mom(); }
      if (kinter.surfaceId() == mu2e::SurfaceIdDetail::TT_Mid) { track->fPTrackerMiddle = kinter.mom(); }
      if (kinter.surfaceId() == mu2e::SurfaceIdDetail::TT_Back) { track->fPTrackerExit = kinter.mom();  }
    }

    KinKal::VEC3 fitmom = kseg->momentum3();
    KinKal::VEC3 pos    = kseg->position3();

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
    track->fP         = kseg->mom();
    track->fPt        = track->Momentum()->Pt();
    track->fChi2      = kffs->chisquared();
    track->fFitCons   = kffs->fitConsistency();
    track->fT0        = kffs->t0().t0();
    track->fT0Err     = kffs->t0().t0Err();
//-----------------------------------------------------------------------------
// momentum error in the first point
//-----------------------------------------------------------------------------
//    ROOT::Math::XYZVector  momdir = fitmom.unit();
    
    track->fFitMomErr = kseg->momerr();
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
    track->fP2        = kseg_exit->mom();
//-----------------------------------------------------------------------------
// helical parameters at Z(TT_FrontPA)
//-----------------------------------------------------------------------------
    try {
      KinKal::CentralHelix helx = kseg->centralHelix();
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

//-----------------------------------------------------------------------------
// virtual detector at the tracker exit: Time at Z(TT_Back)
//-----------------------------------------------------------------------------
//    Hep3Vector vd_tt_back = ds->toDetector(vdet->getGlobal(mu2e::VirtualDetectorId::TT_Back));
    // double     zback      = vd_tt_back.z();
    // double     szb        = s_at_given_z(kffs,zback);

    double     tback      = -1; // FIXME kffs->arrivalTime(szb);
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
    int     id(-1),  npart(0), part_nh[100], part_id[100];
    int     part_pdg_code[100]; 
    int     nwrong = 0;
    double  mcdoca;

    const mu2e::SimParticle   *sim   (nullptr); 
    const mu2e::StrawGasStep  *stgs  (nullptr);

    nss_ch = fSschColl->size();

    if (nss_ch <= 0) {
      printf(">>> ERROR in StntupleInitMu2eTrackBlock: ComboHitCollection by module XXXX is empty, NHITS = %i\n",
	     nss_ch);
    }
    else {
      for (int it=0; it<n_kffs_hits; it++) {
	hit = &hots->at(it);
	mu2e::StrawId sid = hit->strawId();
	const mu2e::Straw* straw = &tracker->straw(sid);
	++ntrkhits;
//-----------------------------------------------------------------------------
// skip calorimeter hit
//-----------------------------------------------------------------------------
	if (! hit) continue;
//-----------------------------------------------------------------------------
// the rest makes sense only for active hits
// all KalSeed hits are "active", figuring out non-active ones 
// now requires comparing the outputs of the seed fit and the full fit
//-----------------------------------------------------------------------------
	if (1) { // hit->isActive()) { // all KalSeed hits are active 
	  loc   = hit->index();
	  if ((loc >= 0) && (loc < nss_ch)) {
	    if ((list_of_mc_straw_hits != NULL) && (list_of_mc_straw_hits->size() > 0)) {

	      const mu2e::StrawDigiMC* mcdigi = &list_of_mc_straw_hits->at(loc);

	      stgs = mcdigi->earlyStrawGasStep().get();
//-----------------------------------------------------------------------------
// count number of active hits with R > 200 um and misassigned drift signs
//-----------------------------------------------------------------------------
	      if (stgs) {
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
	      if (sim != NULL) id = sim->id().asInt();
	      else {
		printf(">>> ERROR in %s : sim is NULL, set PDG_CODE to -1\n",oname);
		id = -1;
	      }

	      found = 0;
	      for (int ip=0; ip<npart; ip++) {
		if (id == part_id[ip]) {
		  found        = 1;
		  part_nh[ip] += 1;
		  break;
		}
	      }
	    
	      if (found == 0) {
		part_id      [npart] = id;
		part_pdg_code[npart] = sim->pdgId();
		part_nh      [npart] = 1;
		npart               += 1;
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

    if (list_of_trk_qual) track->fTrkQual = list_of_trk_qual->at(itrk).MVAOutput();
    else                  track->fTrkQual = -1.e6;
//-----------------------------------------------------------------------------
// defined bit-packed fNActive word
//-----------------------------------------------------------------------------
    track->fNActive   = kffs->hits().size() | (nwrong << 16);
    
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
//-----------------------------------------------------------------------------
// particle parameters at virtual detectors
//-----------------------------------------------------------------------------
    mu2e::GeomHandle<mu2e::VirtualDetector> vdg;

    double t_front(1.e6), t_stout(1.e6);

    track->fPFront = -1.;
    track->fPStOut = -1.;

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
	for (int i=0; i<nvdhits; i++) {
	  const mu2e::StepPointMC* hit = &(*vdhits)[i];
	  
	  //int vdid = hit.volumeId();
	  mu2e::VirtualDetectorId vdid(hit->volumeId());

	  if (vdid.id() == mu2e::VirtualDetectorId::ST_Out) {

	    //	    TAnaDump::Instance()->printStepPointMC(hit,"");

	    art::Ptr<mu2e::SimParticle> const& simptr = hit->simParticle();
	    const mu2e::SimParticle* sim  = simptr.operator ->();

	    if (sim == NULL) {
	      char warning[100];
	      printf(">>> ERROR: %s sim == NULL\n",oname);
	      sprintf(warning,"WARNING: SimParticle for step %i = NULL\n",i);
	      mf::LogWarning(oname) << warning;
	    }
	    int sim_id = sim->id().asInt();
	    if ((sim_id == track->fPartID)  && (hit->time() <  t_stout)) {
	      track->fPStOut = hit->momentum().mag();
	      t_stout        = hit->time();
	    }
	  }
	  else if (vdid.isTrackerFront()) {

	    //	    TAnaDump::Instance()->printStepPointMC(hit,"");

	    art::Ptr<mu2e::SimParticle> const& simptr = hit->simParticle();
	    const mu2e::SimParticle* sim  = simptr.operator ->();

	    if (sim == NULL) {
	      printf(">>> ERROR: %s sim == NULL\n",oname);
	    }
	    int sim_id = sim->id().asInt();
	    if ((sim_id == track->fPartID) && (hit->time() < t_front)) {
	      track->fPFront = hit->momentum().mag();
	      t_front        = hit->time();
	    }
	  }
	}
      }
    }
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
	  while(mother->hasParent())  mother        = mother->parent();
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

    int ndisks = bc->nDisk();

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

    track->fEp = track->fClusterE/track->fP2;
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

