///////////////////////////////////////////////////////////////////////////////
// May 04 2013 P.Murat
//
// in 'XY' mode draw calorimeter clusters as circles with different colors
// in 'Cal' mode draw every detail...
///////////////////////////////////////////////////////////////////////////////
#include "TVirtualX.h"
#include "TPad.h"
#include "TStyle.h"
#include "TVector3.h"
#include "TLine.h"
#include "TArc.h"
#include "TArrow.h"
#include "TBox.h"
#include "TDatabasePDG.h"
#include "TParticlePDG.h"

#include "messagefacility/MessageLogger/MessageLogger.h"

// #include "art/Framework/Principal/Event.h"
// #include "art/Framework/Principal/Handle.h"

#include "Stntuple/gui/TTrkVisNode.hh"

#include "Offline/GeometryService/inc/GeometryService.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/TrackerGeom/inc/Tracker.hh"

#include "Offline/ConditionsService/inc/ConditionsHandle.hh"
#include "Offline/TrackerConditions/inc/StrawResponse.hh"

#include "Offline/DataProducts/inc/StrawId.hh"

#include "Offline/RecoDataProducts/inc/StrawHit.hh"
#include "Offline/RecoDataProducts/inc/StrawHitFlag.hh"
#include "Offline/RecoDataProducts/inc/KalSegment.hh"
#include "Offline/RecoDataProducts/inc/TimeCluster.hh"

#include "Stntuple/gui/TEvdTimeCluster.hh"
#include "Stntuple/gui/TEvdComboHit.hh"
#include "Stntuple/gui/TEvdTrack.hh"
#include "Stntuple/gui/TEvdStraw.hh"
#include "Stntuple/gui/TEvdStrawHit.hh"
#include "Stntuple/gui/TEvdTrkStrawHit.hh"
#include "Stntuple/gui/TEvdStation.hh"
#include "Stntuple/gui/TEvdPanel.hh"
#include "Stntuple/gui/TEvdPlane.hh"
#include "Stntuple/gui/TEvdStrawTracker.hh"
#include "Stntuple/gui/TEvdSimParticle.hh"
#include "Stntuple/gui/TStnVisManager.hh"

#include "Stntuple/obj/TSimpBlock.hh"


// #include "DataProducts/inc/XYZVec.hh"

#include "Stntuple/print/TAnaDump.hh"

ClassImp(TTrkVisNode)

//-----------------------------------------------------------------------------
TTrkVisNode::TTrkVisNode() : TStnVisNode("") {
}

//_____________________________________________________________________________
TTrkVisNode::TTrkVisNode(const char* name, const mu2e::Tracker* Tracker, TStnTrackBlock* TrackBlock):
  TStnVisNode(name) {
  fTracker    = new stntuple::TEvdStrawTracker(Tracker);
  fTrackBlock = TrackBlock;

  fArc        = new TArc;
  fEventTime  = 0;
  fTimeWindow = 1.e6;

  fListOfStrawHits    = new TObjArray();
  fListOfComboHits    = new TObjArray();
  fUseStereoHits      = 0;

  fListOfTracks       = new TObjArray();
  fListOfSimParticles = new TObjArray();

  fChColl             = nullptr;
  fSchColl            = nullptr;
  fShColl             = nullptr;
  fSdmcColl           = nullptr;
  fKsColl             = nullptr;
  fSimpColl           = nullptr;
  fSpmcColl           = nullptr;
					// owned externally (by MuHitDisplay)
  fSimpBlock          = new TSimpBlock();
}

//-----------------------------------------------------------------------------
TTrkVisNode::~TTrkVisNode() {
  delete fArc;

  delete fListOfStrawHits;
  delete fListOfComboHits;

  fListOfTracks->Delete();
  delete fListOfTracks;

  fListOfSimParticles->Delete();
  delete fListOfSimParticles;
}

//-----------------------------------------------------------------------------
int TTrkVisNode::InitEvent() {

  //  const char* oname = "TTrkVisNode::InitEvent";

  mu2e::GeomHandle<mu2e::Tracker> ttHandle;
  const mu2e::Tracker* tracker = ttHandle.get();

  TStnVisManager* vm      = TStnVisManager::Instance();
  const art::Event* event = vm->Event();

  // Tracker calibration object.
  // mu2e::ConditionsHandle<mu2e::StrawResponse> srep = mu2e::ConditionsHandle<mu2e::StrawResponse>("ignored");

  const mu2e::ComboHit              *hit;
  stntuple::TEvdStrawHit            *evd_straw_hit;
  const CLHEP::Hep3Vector           *w;
  const mu2e::Straw                 *straw;

  int                               n_straw_hits, color, nl, ns; // , ipeak, ihit;
  bool                              isFromConversion, intime;
  double                            sigw(1000.), /*vnorm, v,*/ sigr;
  CLHEP::Hep3Vector                 vx0, vx1, vx2;
//-----------------------------------------------------------------------------
// first, clear the cached hit information from the previous event
//-----------------------------------------------------------------------------
  stntuple::TEvdStation*            station;
  stntuple::TEvdPlane*              plane;
  stntuple::TEvdPanel*              panel;

  int                               nst, nplanes, npanels/*, isec*/;

  nst = mu2e::StrawId::_nstations; // tracker->nStations();
  for (int ist=0; ist<nst; ist++) {
    station = fTracker->Station(ist);
    nplanes = station->NPlanes();
    for (int iplane=0; iplane<nplanes; iplane++) {
      plane = station->Plane(iplane);
      npanels = plane->NPanels();
      for (int ipanel=0; ipanel<npanels; ipanel++) {
	panel = plane->Panel(ipanel);
	nl    = panel->NLayers();
	for (int il=0; il<nl; il++) {
	  ns = panel->NStraws(il);
	  for (int is=0; is<ns; is++) {
	    panel->Straw(il,is)->Clear();
	  }
	}
      }
    }
  }
//-----------------------------------------------------------------------------
// offline collections
//-----------------------------------------------------------------------------
  art::Handle<mu2e::ComboHitCollection> chcH;
  event->getByLabel(art::InputTag(fChCollTag), chcH);
  if (chcH.isValid()) fChColl = chcH.product();
  else {
    mf::LogWarning("TTrkVisNode::InitEvent") << " WARNING:" << __LINE__
                                             << " : mu2e::ComboHitCollection "
                                             << fChCollTag << " not found";
    fChColl = nullptr;
  }

  // art::Handle<mu2e::StrawHitFlagCollection> chfcH;
  // event->getByLabel(art::InputTag(fChfCollTag), chfcH);
  // if (chfcH.isValid()) fChfColl = chfcH.product();
  // else {
  //   mf::LogWarning("TTrkVisNode::InitEvent") << " WARNING:" << __LINE__
  //                                            << " : mu2e::StrawHitFlagCollection "
  //                                            << fChfCollTag << " not found";
  //   fChfColl = nullptr;
  // }

  art::Handle<mu2e::ComboHitCollection> schcH;
  event->getByLabel(art::InputTag(fShCollTag), schcH);
  if (schcH.isValid()) fSchColl = schcH.product();
  else {
    mf::LogWarning("TTrkVisNode::InitEvent") << " WARNING:" << __LINE__
                                             << " : mu2e::ComboHitCollection "
                                             << fShCollTag << " not found";
    fSchColl = nullptr;
  }

  // art::Handle<mu2e::StrawHitFlagCollection> shfcH;
  // event->getByLabel(art::InputTag(fShfCollTag), shfcH);
  // if (shfcH.isValid()) fShfColl = shfcH.product();
  // else {
  //   mf::LogWarning("TTrkVisNode::InitEvent") << " WARNING:" << __LINE__
  //                                            << " : mu2e::StrawHitFlagCollection "
  //                                            << fShfCollTag << " not found";
  //   fShfColl = nullptr;
  // }

  art::Handle<mu2e::StrawHitCollection> shcH;
  event->getByLabel(art::InputTag(fShCollTag), shcH);
  if (shcH.isValid()) fShColl = shcH.product();
  else {
    mf::LogWarning("TTrkVisNode::InitEvent") << " WARNING:" << __LINE__
                                             << " : mu2e::StrawHitCollection "
                                             << fShCollTag << " not found";
    fShColl = nullptr;
  }

  art::Handle<mu2e::StrawDigiMCCollection> sdmccH;
  event->getByLabel(art::InputTag(fSdmcCollTag), sdmccH);
  if (sdmccH.isValid()) fSdmcColl = sdmccH.product();
  else {
    mf::LogWarning("TTrkVisNode::InitEvent") << " WARNING:" << __LINE__
                                             << " : mu2e::StrawDigiMCCollection "
                                             << fSdmcCollTag << " not found";
    fSdmcColl = nullptr;
  }

  art::Handle<mu2e::KalSeedCollection> kscH;
  event->getByLabel(art::InputTag(fKsCollTag), kscH);
  if (kscH.isValid()) fKsColl = kscH.product();
  else {
    mf::LogWarning("TTrkVisNode::InitEvent") << " WARNING:" << __LINE__
                                             << " : mu2e::KalSeedCollection "
                                             << fKsCollTag << " not found";
    fKsColl = nullptr;
  }
//-----------------------------------------------------------------------------
// display hits corresponding to a given time peak, or all hits,
// if the time peak is not found
//-----------------------------------------------------------------------------
  fListOfStrawHits->Delete();

  stntuple::TEvdStraw* evd_straw;

  n_straw_hits = 0;
  if (fSchColl != nullptr) n_straw_hits = fSchColl->size();

  for (int ihit=0; ihit<n_straw_hits; ihit++ ) {

    hit              = &fSchColl->at(ihit);
    straw            = &tracker->getStraw(hit->strawId());
    color            = kBlack;
    intime           = 0;
    isFromConversion = false;
//-----------------------------------------------------------------------------
// deal with MC information - later
//-----------------------------------------------------------------------------
    const mu2e::StrawDigiMC*             mcdigi(nullptr);
    if (fSdmcColl and (fSdmcColl->size() > 0)) {
      mcdigi = &fSdmcColl->at(ihit);
    }

    if (mcdigi) {

      const mu2e::StrawGasStep* step = mcdigi->earlyStrawGasStep().get();
      const mu2e::SimParticle*  sim  = &(*step->simParticle());

      if (sim->fromGenerator()) {
	mu2e::GenParticle* gen = (mu2e::GenParticle*) &(*sim->genParticle());
	//	    if ( gen->generatorId() == mu2e::GenId::conversionGun ){
	if ( gen->generatorId() == mu2e::GenId::StoppedParticleReactionGun ){
	  isFromConversion = true;
	}
      }
      int   pdg_id = sim->pdgId();
      float mc_mom = step->momvec().mag();

      intime = fabs(hit->time()-fEventTime) < fTimeWindow;

      if      (pdg_id == 11) {
	if    (mc_mom > 20  ) {
          // if (intime) color = kRed;
          // else        color = kAzure+1;
          color = kRed;
	}
	else                   { color = kRed+2;    }
      }
      else if (pdg_id ==  -11) { color = kBlue;     }
      else if (pdg_id ==   13) { color = kGreen+2;  }
      else if (pdg_id ==  -13) { color = kGreen-2;  }
      else if (pdg_id == 2212) { color = kBlue+2;   }
      else if (pdg_id ==  211) { color = kOrange-3; }
      else if (pdg_id == -211) { color = kOrange+7; }

      else                     { color = kBlack;    }
    }
//-----------------------------------------------------------------------------
// add a pointer to the hit to the straw
//-----------------------------------------------------------------------------
    sigw     = hit->wireRes();      // P.Murat
    sigr     = 2.5;                    // in mm
    int mask = 0;
    if (intime          ) mask |= stntuple::TEvdStrawHit::kInTimeBit;
    if (isFromConversion) mask |= stntuple::TEvdStrawHit::kConversionBit;

    int ist, ipl, ippl, /*ifc,*/ ipn, il, is;

    ipl  = straw->id().getPlane();      // plane number here runs from 0 to 2*NStations-1
    ist  = straw->id().getStation();
    ippl = ipl % 2 ;                    // plane number within the station
    ipn  = straw->id().getPanel();
    il   = straw->id().getLayer();
    is   = straw->id().getStraw();

    w             = &straw->getDirection();
    evd_straw     = fTracker->Station(ist)->Plane(ippl)->Panel(ipn)->Straw(il,is/2);
    evd_straw_hit = new stntuple::TEvdStrawHit(hit,
					       evd_straw,
					       mcdigi,
					       hit->pos().x(),
					       hit->pos().y(),
					       hit->pos().z(),
					       w->x(),w->y(),
					       sigw,sigr,
					       mask,color);
    evd_straw->AddHit(evd_straw_hit);
//-----------------------------------------------------------------------------
// so far, lists of straw hits and combo hits are the same, need to make them different
//-----------------------------------------------------------------------------
    fListOfStrawHits->Add(evd_straw_hit);
  }
//-----------------------------------------------------------------------------
// combo hits
//-----------------------------------------------------------------------------
  fListOfComboHits->Delete();
  int nch  = 0;
  if (fChColl != nullptr) nch = fChColl->size();
//-----------------------------------------------------------------------------
// the rest makes sense only if nhits > 0
//-----------------------------------------------------------------------------
  if (nch > 0) {
    //    const mu2e::ComboHit* hit0 = &fChColl->at(0);

    float                     mc_mom(-1.), mc_mom_z(-1.);
    int                       mother_pdg_id(0);
    const mu2e::SimParticle*  mother(nullptr);
    const mu2e::StrawGasStep* step  (nullptr);
    const mu2e::SimParticle*  sim   (nullptr);

    for (int ihit=0; ihit<nch; ihit++ ) {
      const mu2e::ComboHit* hit = &fChColl->at(ihit);
//-----------------------------------------------------------------------------
// handle MC truth, if that is present
//-----------------------------------------------------------------------------
      if (fSdmcColl != nullptr) {
        int ind = hit->indexArray().at(0);
	const mu2e::StrawDigiMC* mcdigi = &fSdmcColl->at(ind);

	step = mcdigi->earlyStrawGasStep().get();

	const art::Ptr<mu2e::SimParticle>& simptr = step->simParticle();
	sim = simptr.operator->();

	art::Ptr<mu2e::SimParticle>        momptr = simptr;

	while (momptr->hasParent()) momptr = momptr->parent();
	mother = momptr.operator->();

	mother_pdg_id = mother->pdgId();
	mc_mom        = step->momvec().mag();
	mc_mom_z      = step->momvec().z();

	// if (simptr->fromGenerator()) generator_id = simptr->genParticle()->generatorId().id();
	// else                         generator_id = -1;

      }
//-----------------------------------------------------------------------------
// store TEvdComboHit
//-----------------------------------------------------------------------------
      fListOfComboHits->Add(new stntuple::TEvdComboHit(hit,sim,step,mother_pdg_id,mc_mom,mc_mom_z));
    }
  }
//-----------------------------------------------------------------------------
// hit MC truth unformation from StepPointMC's
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// now initialize tracks
//-----------------------------------------------------------------------------
  stntuple::TEvdTrack      *trk;
  const mu2e::KalSeed      *kseed;
  // const mu2e::TrkStrawHit  *track_hit;

  fListOfTracks->Delete();
  int ntrk = 0;

  if (fKsColl) ntrk = fKsColl->size();

  for (int i=0; i<ntrk; i++) {
    kseed = &fKsColl->at(i);
    trk  = new stntuple::TEvdTrack(i,kseed);
//-----------------------------------------------------------------------------
// add hits, skip calorimeter clusters (TrkCaloHit's)
//-----------------------------------------------------------------------------
    //    const TrkHitVector* hits = &kseed->hitVector();
    const std::vector<mu2e::TrkStrawHitSeed>* hits = &kseed->hits();

    for (auto it=hits->begin(); it!=hits->end(); it++) {
      const mu2e::TrkStrawHitSeed* hit =  &(*it);
      // track_hit = dynamic_cast<mu2e::TrkStrawHit*> (*it);
      if (hit == nullptr) continue;
      // need to find this hit in the list of TEvdStrawHits (already existing) ... later

      const mu2e::Straw* straw = &tracker->straw(hit->strawId());
      stntuple::TEvdTrkStrawHit* evd_hit = new stntuple::TEvdTrkStrawHit(hit,straw);
      trk->AddHit(evd_hit);
    }

    fListOfTracks->Add(trk);
  }
//-----------------------------------------------------------------------------
// initialize SimParticles
//-----------------------------------------------------------------------------
  fListOfSimParticles->Delete();

  stntuple::TEvdSimParticle *esim;

  TDatabasePDG* pdb = TDatabasePDG::Instance();

  int np = fSimpBlock->NParticles();
  int ipp = 0;
  for (int i=0; i<np; i++) {
    TSimParticle* tsimp = fSimpBlock->Particle(i);

    TParticlePDG* p_pdg = pdb->GetParticle(tsimp->PDGCode());

    if (p_pdg == nullptr) {
      printf("WARNING: Unknown SimParticle with code %i, SKIP\n",tsimp->PDGCode());
      continue;
    }

    if (tsimp->NStrawHits() > 0) {
//-----------------------------------------------------------------------------
// taking the first and the last StrawGasStep's in hope that they are ordered in Z
//-----------------------------------------------------------------------------
      int first = tsimp->Shid()->front();
      int last  = tsimp->Shid()->back ();

      const mu2e::StrawGasStep* s1 = fSdmcColl->at(first).earlyStrawGasStep().get();
      const mu2e::StrawGasStep* s2 = fSdmcColl->at(last ).earlyStrawGasStep().get();

      esim = new stntuple::TEvdSimParticle(ipp,tsimp,s1,s2);
      fListOfSimParticles->Add(esim);
      ipp++;
    }
  }

  return 0;
}

//-----------------------------------------------------------------------------
// check if a hit with 'Index' in the fShColl collection belongs to the time cluster
// time cluster is made out of ComboHit's
// the code could be made a bit more efficient for combo hits at a cost of passing the
// the hit type explicitly
//-----------------------------------------------------------------------------
int TTrkVisNode::TCHit(const mu2e::TimeCluster* TimeCluster, int Index) {
  int ok(0);

  int nch = TimeCluster->nhits();
  for (int i=0; i<nch; i++) {
    int ind = TimeCluster->hits().at(i);
    const mu2e::ComboHit* ch = &fChColl->at(ind);
    int nsh = ch->nStrawHits();
    for (int ish=0; ish<nsh; ish++) {
      if (ch->index(ish) == Index) {
	ok = 1;
	break;
      }
    }
  }

  return ok;
}


//-----------------------------------------------------------------------------
// draw reconstructed tracks and STRAW hits, may want to display COMBO hits instead
//-----------------------------------------------------------------------------
void TTrkVisNode::PaintXY(Option_t* Option) {

  double                  time;
  int                     station, ntrk(0);

  const mu2e::Straw      *straw;

  //  int view_type = TVisManager::Instance()->GetCurrentView()->Type();

  mu2e::GeomHandle<mu2e::Tracker> ttHandle;
  const mu2e::Tracker* tracker = ttHandle.get();

  TStnVisManager* vm = TStnVisManager::Instance();

  double tmin = vm->TMin();
  double tmax = vm->TMax();

  stntuple::TEvdTimeCluster* etcl = vm->SelectedTimeCluster();
  if (etcl) {
    tmin = etcl->TMin(); // FIXME!
    tmax = etcl->TMax(); // FIXME!

    if (vm->DisplayStrawHitsXY()) {
//-----------------------------------------------------------------------------
// when displaying straw hits, tmin and tmax should be defined by the straw,
// not combo, hit times
//-----------------------------------------------------------------------------
      const mu2e::TimeCluster* tc = etcl->TimeCluster();
      int nch = tc->nhits();
      for (int i=0; i<nch; i++) {
	int ind = tc->hits().at(i);
	const mu2e::ComboHit* ch = &fChColl->at(ind);
	int nsh = ch->nStrawHits();
	for (int ish=0; ish<nsh; ish++) {
	  int loc = ch->index(ish);
	  const mu2e::ComboHit* sh = &fSchColl->at(loc);
	  if (sh->correctedTime() < tmin) tmin = sh->correctedTime();
	  if (sh->correctedTime() > tmax) tmax = sh->correctedTime();
	}
      }
    }
  }

  if (vm->DisplayStrawHitsXY()) {
//-----------------------------------------------------------------------------
// display straw hits
//-----------------------------------------------------------------------------
    int nhits = fListOfStrawHits->GetEntries();
    if (nhits > 0) {
      // const mu2e::ComboHit* sch0 = &fSchColl->at(0);
      for (int i=0; i<nhits; i++) {
	stntuple::TEvdStrawHit* evd_sh = GetEvdStrawHit(i);
	const mu2e::ComboHit*   sch    = evd_sh->StrawHit();
//-----------------------------------------------------------------------------
// see in flags need to be checked, use external flags
//-----------------------------------------------------------------------------
	// int loc  = sch-sch0;
	const mu2e::StrawHitFlag flag = sch->flag(); // &fShfColl->at(loc);
	if (vm->IgnoreComptonHits()) {
	  if (flag.hasAnyProperty(mu2e::StrawHitFlagDetail::bkg))      continue;
	}

	if (vm->IgnoreProtonHits()) {
	  if (! flag.hasAnyProperty(mu2e::StrawHitFlagDetail::energysel)) continue;
	}

	straw     = &tracker->getStraw(sch->strawId());
	station   = straw->id().getStation();
	time      = sch->correctedTime();

	if ((station >= vm->MinStation()) && (station <= vm->MaxStation())) {
	  if ((time >= tmin) && (time <= tmax)) {
//-----------------------------------------------------------------------------
// check if the hit belongs to the time cluster
//-----------------------------------------------------------------------------
	    int ok = 1;
	    if (etcl and vm->DisplayOnlyTCHits()) {
	      ok = TCHit(etcl->TimeCluster(),sch->index());
	    }
	    if (ok  ) evd_sh->PaintXY(Option);
	  }
	}
      }
    }
  }
  else {
//-----------------------------------------------------------------------------
// display combo hits
//-----------------------------------------------------------------------------
    int nch = fListOfComboHits->GetEntries();
    if (nch > 0) {
      // const mu2e::ComboHit* ch0 = &fChColl->at(0);
      for (int i=0; i<nch; i++) {
	stntuple::TEvdComboHit* evd_ch = (stntuple::TEvdComboHit*) fListOfComboHits->At(i);
	const mu2e::ComboHit*   ch     = evd_ch->ComboHit();

	// int loc  = ch-ch0;
	const mu2e::StrawHitFlag flag = ch->flag(); // &fChfColl->at(loc);
	if (vm->IgnoreComptonHits()) {
	  if (flag.hasAnyProperty(mu2e::StrawHitFlagDetail::bkg))      continue;
	}
	if (vm->IgnoreProtonHits()) {
	  if (! flag.hasAnyProperty(mu2e::StrawHitFlagDetail::energysel)) continue;
	}
//-----------------------------------------------------------------------------
// a combined hit doesn't have it's own measured time - that is calculated
// look at the first straw
//-----------------------------------------------------------------------------
        int index = ch->index(0);
        const mu2e::ComboHit* sh = &fSchColl->at(index);

        straw     = &tracker->getStraw(sh->strawId());           // first straw hit
        station   = straw->id().getStation();
        time      = ch->correctedTime();

        if ((station >= vm->MinStation()) && (station <= vm->MaxStation())) {
          if ((time >= tmin) && (time <= tmax)) {
            // check if the hit belongs to the time cluster
            int ok = 1;
            if (etcl and vm->DisplayOnlyTCHits()) {
              ok = TCHit(etcl->TimeCluster(),ch->index(0));
            }
            if (ok  ) evd_ch->PaintXY(Option);
          }
        }
      }
    }
  }
//-----------------------------------------------------------------------------
// now - tracks
//-----------------------------------------------------------------------------
  if (vm->DisplayTracks()) {
    stntuple::TEvdTrack* evd_trk;
  //  TAnaDump::Instance()->printKalRep(0,"banner");

    if ( (fListOfTracks) != 0 )  ntrk = fListOfTracks->GetEntriesFast();

    for (int i=0; i<ntrk; i++ ) {
      evd_trk = GetEvdTrack(i);
      evd_trk->Paint(Option);
    }
  }
//-----------------------------------------------------------------------------
// SimParticle's
//-----------------------------------------------------------------------------
  if (vm->DisplaySimParticles()) {

    stntuple::TEvdSimParticle* esp;

    int nsim = fListOfSimParticles->GetEntriesFast();

    for (int i=0; i<nsim; i++ ) {
      esp = GetEvdSimParticle(i);
      TSimParticle* simp = esp->SimParticle();

      if (vm->IgnoreProtons() != 0) {
        int pdg_code = simp->PDGCode();
        if (pdg_code > 2000)                                         continue;
      }

      float mom          = simp->StartMom()->P();
      if ((mom >= vm->MinSimpMomentum()) and (mom <= vm->MaxSimpMomentum())) {
        esp->PaintXY(Option);
      }
    }
  }
//-----------------------------------------------------------------------------
// seedfit, if requested - not implemented yet
//-----------------------------------------------------------------------------
//   TAnaDump::Instance()->printKalRep(0,"banner");

//   ntrk = fListOfTracks->GetEntriesFast();
//   for (int i=0; i<ntrk; i++ ) {
//     evd_trk = (TEvdTrack*) fListOfTracks->At(i);
//     evd_trk->Paint(Option);
//   }

  gPad->Modified();
}

//-----------------------------------------------------------------------------
// in RZ view can display only straw hits, they are painted by fTracker
//-----------------------------------------------------------------------------
void TTrkVisNode::PaintRZ(Option_t* Option) {
  int             ntrk(0), nhits;
  stntuple::TEvdTrack*      evd_trk;

  TStnVisManager* vm = TStnVisManager::Instance();

  fTracker->PaintRZ(Option);

  double tmin = vm->TMin();
  double tmax = vm->TMax();

  stntuple::TEvdTimeCluster* etcl = vm->SelectedTimeCluster();
  if (etcl) {
    tmin = etcl->TMin(); // FIXME!
    tmax = etcl->TMax(); // FIXME!

    if (vm->DisplayStrawHitsXY()) {
//-----------------------------------------------------------------------------
// for straw hit display, tmin and tmax should be defined by the straw,
// not combo, hit times
//-----------------------------------------------------------------------------
      const mu2e::TimeCluster* tc = etcl->TimeCluster();
      int nch = tc->nhits();
      for (int i=0; i<nch; i++) {
	int ind = tc->hits().at(i);
	const mu2e::ComboHit* ch = &fChColl->at(ind);
	int nsh = ch->nStrawHits();
	for (int ish=0; ish<nsh; ish++) {
	  int loc = ch->index(ish);
	  const mu2e::ComboHit* sh = &fSchColl->at(loc);
	  if (sh->correctedTime() < tmin) tmin = sh->correctedTime();
	  if (sh->correctedTime() > tmax) tmax = sh->correctedTime();
	}
      }
    }
  }
//-----------------------------------------------------------------------------
// do not draw all straw hits - just redraw straws in different color instead
//-----------------------------------------------------------------------------
//   int nhits = fListOfStrawHits->GetEntries();
//   for (int i=0; i<nhits; i++) {
//     hit       = (TEvdStrawHit*) fListOfStrawHits->At(i);

//     if ((station >= vm->MinStation()) && (station <= vm->MaxStation())) continue;
//     if ((time    <  vm->TMin()      ) || (time     > vm->TMax()      )) continue;

//     hit->Paint(Option);
//   }
//-----------------------------------------------------------------------------
// display tracks and track hits
//-----------------------------------------------------------------------------
  if (vm->DisplayTracks()) {
    if (fListOfTracks != 0)  ntrk = fListOfTracks->GetEntriesFast();

    for (int i=0; i<ntrk; i++ ) {
      evd_trk = (stntuple::TEvdTrack*) fListOfTracks->At(i);
      evd_trk->Paint(Option);

      nhits = evd_trk->NHits();
      for (int ih=0; ih<nhits; ih++) {
	stntuple::TEvdTrkStrawHit* hit = evd_trk->Hit(ih);
	double time = hit->TrkStrawHitSeed()->time();
	if ((time >= tmin) && (time <= tmax)) {
	  hit->PaintRZ(Option);
	}
      }
    }
  }
//-----------------------------------------------------------------------------
// SimParticle's : pT at  the ST is too large, need to use parameters at the tracker entrance ?
//-----------------------------------------------------------------------------
  if (vm->DisplaySimParticles()) {
    stntuple::TEvdSimParticle* esim;
    int nsim(0);

    if ( (fListOfSimParticles) != 0 )  nsim = fListOfSimParticles->GetEntriesFast();

    for (int i=0; i<nsim; i++ ) {
      esim = (stntuple::TEvdSimParticle*) fListOfSimParticles->At(i);
      TSimParticle* simp = esim->SimParticle();

      if (vm->IgnoreProtons() != 0) {
        int pdg_code = simp->PDGCode();
        if (pdg_code > 2000)                                           continue;
      }

      float mom          = simp->StartMom()->P();
      if ((mom >= vm->MinSimpMomentum()) and (mom <= vm->MaxSimpMomentum())) {
        esim->PaintRZ(Option);
      }
    }
  }

  gPad->Modified();
}


//-----------------------------------------------------------------------------
// TZ view is only for the pattern recognition / time cluster finding
// display reconstructed tracks and combo hits
//-----------------------------------------------------------------------------
void TTrkVisNode::PaintTZ(Option_t* Option) {

  TStnVisManager* vm = TStnVisManager::Instance();

  double tmin = vm->TMin();
  double tmax = vm->TMax();

  stntuple::TEvdTimeCluster* etcl = vm->SelectedTimeCluster();

  if (etcl) {
    tmin = etcl->TMin(); // FIXME!
    tmax = etcl->TMax(); // FIXME!
  }

  int nhits = fListOfComboHits->GetEntries();
  if (nhits > 0) {
    // const mu2e::ComboHit* ch0 = &fChColl->at(0);

    for (int i=0; i<nhits; i++) {
      stntuple::TEvdComboHit* ech = (stntuple::TEvdComboHit*) fListOfComboHits->At(i);
      const mu2e::ComboHit*   ch  = ech->ComboHit();

      // int loc = ch-ch0;
      const mu2e::StrawHitFlag flag = ch->flag(); // &fChfColl->at(loc);
      if (vm->IgnoreComptonHits()) {
        if (flag.hasAnyProperty(mu2e::StrawHitFlagDetail::bkg))      continue;
      }

      if (vm->IgnoreProtonHits()) {
        if (! flag.hasAnyProperty(mu2e::StrawHitFlagDetail::energysel)) continue;
      }

      float time  = ech->correctedTime();

      if ((time >= tmin) && (time <= tmax)) {
        // check if the hit belongs to the time cluster
        int ok = 1;
        if (etcl and vm->DisplayOnlyTCHits()) {
          ok = TCHit(etcl->TimeCluster(),ech->ComboHit()->index(0));
        }
        if (ok) ech->PaintTZ(Option);
      }
    }
  }
//-----------------------------------------------------------------------------
// SimParticle's
//-----------------------------------------------------------------------------
  if (vm->DisplaySimParticles()) {
    stntuple::TEvdSimParticle* esim;
    int nsim(0);

    if ( (fListOfSimParticles) != 0 )  nsim = fListOfSimParticles->GetEntriesFast();

    for (int i=0; i<nsim; i++ ) {
      esim = (stntuple::TEvdSimParticle*) fListOfSimParticles->At(i);

      if (vm->IgnoreProtons() != 0) {
        int pdg_code = esim->SimParticle()->PDGCode();
        if (pdg_code > 2000)                                         continue;
      }

      esim->PaintTZ(Option);
    }
  }

  gPad->Modified();
}

//-----------------------------------------------------------------------------
// PhiZ view is only for the pattern recognition / time cluster finding
// display reconstructed tracks and combo hits
//-----------------------------------------------------------------------------
void TTrkVisNode::PaintPhiZ(Option_t* Option) {

  TStnVisManager* vm = TStnVisManager::Instance();

  mu2e::GeomHandle<mu2e::Tracker> ttHandle;
  const mu2e::Tracker* tracker = ttHandle.get();

  double tmin   = vm->TMin();
  double tmax   = vm->TMax();

  double phimin = -M_PI; // vm->TMin();
  double phimax =  M_PI; // vm->TMax();

  stntuple::TEvdTimeCluster* etcl = vm->SelectedTimeCluster();

  if (etcl) {
    tmin = etcl->TMin(); // FIXME!
    tmax = etcl->TMax(); // FIXME!
  }

  int nhits = fListOfComboHits->GetEntries();
  if (nhits > 0) {
    // const mu2e::ComboHit* ch0 = &fChColl->at(0);

    for (int i=0; i<nhits; i++) {
      stntuple::TEvdComboHit* ech = (stntuple::TEvdComboHit*) fListOfComboHits->At(i);
      const mu2e::ComboHit*   ch  = ech->ComboHit();

      // int loc = ch-ch0;
      const mu2e::StrawHitFlag flag = ch->flag(); // &fChfColl->at(loc);
      if (vm->IgnoreComptonHits()) {
        if (flag.hasAnyProperty(mu2e::StrawHitFlagDetail::bkg))      continue;
      }

      if (vm->IgnoreProtonHits()) {
        if (! flag.hasAnyProperty(mu2e::StrawHitFlagDetail::energysel)) continue;
      }

      int index = ch->index(0);
      const mu2e::ComboHit* sh = &fSchColl->at(index);

      const mu2e::Straw* straw = &tracker->getStraw(sh->strawId()); // first straw hit
      int station = straw->id().getStation();
      double time = ch->correctedTime();

      if ((station >= vm->MinStation()) && (station <= vm->MaxStation())) {
        if ((time >= tmin) && (time <= tmax)) {
          float phi = ech->Pos()->Phi();

          if ((phi >= phimin) && (phi <= phimax)) {
                                        // check if the hit belongs to the time cluster
            int ok = 1;
            if (etcl and vm->DisplayOnlyTCHits()) {
              ok = TCHit(etcl->TimeCluster(),ech->ComboHit()->index(0));
            }
            if (ok) ech->PaintPhiZ(Option);
          }
        }
      }
    }
  }
//-----------------------------------------------------------------------------
// SimParticle's
//-----------------------------------------------------------------------------
  if (vm->DisplaySimParticles()) {
    stntuple::TEvdSimParticle* esim;
    int nsim(0);

    if ( (fListOfSimParticles) != 0 )  nsim = fListOfSimParticles->GetEntriesFast();

    for (int i=0; i<nsim; i++ ) {
      esim = (stntuple::TEvdSimParticle*) fListOfSimParticles->At(i);

      if (vm->IgnoreProtons() != 0) {
        int pdg_code = esim->SimParticle()->PDGCode();
        if (pdg_code > 2000)                                         continue;
      }

      esim->PaintPhiZ(Option);
    }
  }

  gPad->Modified();
}

//-----------------------------------------------------------------------------
// VST view : display all straws
//-----------------------------------------------------------------------------
void TTrkVisNode::PaintVST(Option_t* Option) {

  fTracker->PaintVST(Option);

  gPad->Modified();
}

//-----------------------------------------------------------------------------
int TTrkVisNode::DistancetoPrimitiveXY(Int_t px, Int_t py) {
  static TVector3 global;
  global.SetXYZ(gPad->AbsPixeltoX(px),gPad->AbsPixeltoY(py),0);

  TObject* closest(nullptr);

  int  x1, y1, dx1, dy1, min_dist(9999), dist;
  TStnVisManager* vm = TStnVisManager::Instance();

  if (vm->DisplayStrawHitsXY() == 1) {
    int nhits = fListOfStrawHits->GetEntries();
    for (int i=0; i<nhits; i++) {
      stntuple::TEvdStrawHit* hit = GetEvdStrawHit(i);

      x1  = gPad->XtoAbsPixel(hit->Pos()->X());
      y1  = gPad->YtoAbsPixel(hit->Pos()->Y());
      dx1 = px-x1;
      dy1 = py-y1;

      dist  = (int) sqrt(dx1*dx1+dy1*dy1);
      if (dist < min_dist) {
	min_dist = dist;
	closest  = hit;
      }
    }
  }
  else {
    int nhits = fListOfComboHits->GetEntries();
    for (int i=0; i<nhits; i++) {
      stntuple::TEvdComboHit* hit = GetEvdComboHit(i);

      x1  = gPad->XtoAbsPixel(hit->Pos()->X());
      y1  = gPad->YtoAbsPixel(hit->Pos()->Y());
      dx1 = px-x1;
      dy1 = py-y1;

      dist  = (int) sqrt(dx1*dx1+dy1*dy1);
      if (dist < min_dist) {
	min_dist = dist;
	closest  = hit;
      }
    }
  }

//-----------------------------------------------------------------------------
// tracks are represented by ellipses
//-----------------------------------------------------------------------------
  if (vm->DisplayTracks()) {
    int ntracks = fListOfTracks->GetEntries();
    for (int i=0; i<ntracks; i++) {
      stntuple::TEvdTrack* trk = GetEvdTrack(i);

      dist = trk->DistancetoPrimitiveXY(px,py);

      if (dist < min_dist) {
	min_dist = dist;
	closest  = trk;
      }
    }
  }
//-----------------------------------------------------------------------------
// simparticles are represented by ellipses
//-----------------------------------------------------------------------------
  if (vm->DisplaySimParticles()) {
    int nsim = fListOfSimParticles->GetEntries();
    for (int i=0; i<nsim; i++) {
      stntuple::TEvdSimParticle* sim = GetEvdSimParticle(i);

      if (vm->IgnoreProtons() != 0) {
        int pdg_code = sim->SimParticle()->PDGCode();
        if (pdg_code > 2000)                                         continue;
      }

      dist = sim->DistancetoPrimitiveXY(px,py);

      if (dist < min_dist) {
	min_dist = dist;
	closest  = sim;
      }
    }
  }
  SetClosestObject(closest,min_dist);

  return min_dist;
}

//-----------------------------------------------------------------------------
Int_t TTrkVisNode::DistancetoPrimitiveRZ(Int_t px, Int_t py) {
  return 9999;
}


//-----------------------------------------------------------------------------
Int_t TTrkVisNode::DistancetoPrimitiveTZ(Int_t px, Int_t py) {

  // static TVector3 global;
  // global.SetXYZ(gPad->AbsPixeltoX(px),gPad->AbsPixeltoY(py),0);

  TObject* closest(nullptr);

  int  x1, y1, dx1, dy1, min_dist(9999), dist;

  TStnVisManager* vm = TStnVisManager::Instance();

  int nhits = fListOfComboHits->GetEntries();
  for (int i=0; i<nhits; i++) {
    stntuple::TEvdComboHit* hit = (stntuple::TEvdComboHit*) fListOfComboHits->At(i);
    x1  = gPad->XtoAbsPixel(hit->Z());
    y1  = gPad->YtoAbsPixel(hit->correctedTime());
    dx1 = px-x1;
    dy1 = py-y1;

    dist  = (int) sqrt(dx1*dx1+dy1*dy1);
    if (dist < min_dist) {
      min_dist = dist;
      closest  = hit;
    }
  }
//-----------------------------------------------------------------------------
// simparticles are represented by lines
//-----------------------------------------------------------------------------
  if (vm->DisplaySimParticles()) {
    int nsim = fListOfSimParticles->GetEntries();
    for (int i=0; i<nsim; i++) {
      stntuple::TEvdSimParticle* esim = GetEvdSimParticle(i);

      if (vm->IgnoreProtons() != 0) {
        int pdg_code = (int) esim->SimParticle()->PDGCode();
        if (pdg_code > 2000)                                         continue;
      }

      dist = esim->DistancetoPrimitiveTZ(px,py);

      if (dist < min_dist) {
	min_dist = dist;
	closest  = esim;
      }
    }
  }

  SetClosestObject(closest,min_dist);

  return min_dist;
}

//-----------------------------------------------------------------------------
Int_t TTrkVisNode::DistancetoPrimitivePhiZ(Int_t px, Int_t py) {

  // static TVector3 global;
  // global.SetXYZ(gPad->AbsPixeltoX(px),gPad->AbsPixeltoY(py),0);

  TObject* closest(nullptr);

  int  x1, y1, dx1, dy1, min_dist(9999), dist;

  TStnVisManager* vm = TStnVisManager::Instance();

  int nhits = fListOfComboHits->GetEntries();
  for (int i=0; i<nhits; i++) {
    stntuple::TEvdComboHit* hit = (stntuple::TEvdComboHit*) fListOfComboHits->At(i);
    x1  = gPad->XtoAbsPixel(hit->Z());
    y1  = gPad->YtoAbsPixel(hit->Pos()->Phi());
    dx1 = px-x1;
    dy1 = py-y1;

    dist  = (int) sqrt(dx1*dx1+dy1*dy1);
    if (dist < min_dist) {
      min_dist = dist;
      closest  = hit;
    }
  }
//-----------------------------------------------------------------------------
// simparticles are represented by lines
//-----------------------------------------------------------------------------
  if (vm->DisplaySimParticles()) {
    int nsim = fListOfSimParticles->GetEntries();
    for (int i=0; i<nsim; i++) {
      stntuple::TEvdSimParticle* esim = GetEvdSimParticle(i);

      if (vm->IgnoreProtons() != 0) {
        int pdg_code = (int) esim->SimParticle()->PDGCode();
        if (pdg_code > 2000)                                         continue;
      }

      dist = esim->DistancetoPrimitivePhiZ(px,py);

      if (dist < min_dist) {
	min_dist = dist;
	closest  = esim;
      }
    }
  }

  SetClosestObject(closest,min_dist);

  return min_dist;
}


//-----------------------------------------------------------------------------
void TTrkVisNode::Clear(Option_t* Opt) {
  printf(">>> name: %s TTrkVisNode::Clear is not implemented yet\n",GetName());
}

//-----------------------------------------------------------------------------
void TTrkVisNode::Print(Option_t* Opt) const {

  TString opt = Opt;
  opt.ToLower();

  if (opt == "combo_hits") {
    TAnaDump* ad = TAnaDump::Instance();
    ad->printComboHitCollection(fChCollTag.data(),fSdmcCollTag.encode().data());
    return;
  }
//-----------------------------------------------------------------------------
// print SimParticles
//-----------------------------------------------------------------------------
  stntuple::TEvdSimParticle* sim;
  int nsim(0);

  if ( (fListOfSimParticles) != 0 )  nsim = fListOfSimParticles->GetEntriesFast();

  printf("n(sim particles) = %i\n",nsim);
  int banner_printed(0);
  for (int i=0; i<nsim; i++ ) {
    sim = (stntuple::TEvdSimParticle*) fListOfSimParticles->At(i);
    if (banner_printed == 0) {
      sim->Print("banner");
      banner_printed = 1;
    }
    sim->Print("data");
  }
//-----------------------------------------------------------------------------
// print ComboHits
//-----------------------------------------------------------------------------
  banner_printed = 0;
  int nch = fListOfComboHits->GetEntries();
  printf("n(combo hits) = %i\n",nch);

  for (int i=0; i<nch; i++) {
    stntuple::TEvdComboHit* hit = (stntuple::TEvdComboHit*) fListOfComboHits->At(i);
    if (banner_printed == 0) {
      hit->Print("banner");
      banner_printed = 1;
    }
    hit->Print("data");
  }
}


//-----------------------------------------------------------------------------
void TTrkVisNode::NodePrint(const void* Object, const char* ClassName) {
  TString class_name(ClassName);

  TAnaDump* ad = TAnaDump::Instance();

  if (class_name == "ComboHit") {
//-----------------------------------------------------------------------------
// print a ComboHit or a collection of those
//-----------------------------------------------------------------------------
    if (Object) ad->printComboHit          ((const mu2e::ComboHit*) Object,nullptr);
    else        ad->printComboHitCollection(fChCollTag.data(),fSdmcCollTag.encode().data());
  }
  else if (class_name == "KalSeed") {
//-----------------------------------------------------------------------------
// print a KalSeed or a KalSeed collection
//-----------------------------------------------------------------------------
    if (Object) ad->printKalSeed          ((const mu2e::KalSeed*) Object,"",fShCollTag.encode().data(),fSdmcCollTag.encode().data());
    else        ad->printKalSeedCollection(fKsCollTag.data(),1,fShCollTag.encode().data(),fSdmcCollTag.encode().data());
  }
  else if (class_name == "SimParticle") {
    ad->printSimParticleCollection(fSimpCollTag);
  }
  else if (class_name == "StrawHit") {
    if (Object) ad->printStrawHit          ((const mu2e::StrawHit*) Object,nullptr);
    else        ad->printStrawHitCollection(fShCollTag.encode().data(),fSdmcCollTag.encode().data());
  }
  else {
    printf("WARNING in TTrkVisNode::Print: print for %s not implemented yet\n",ClassName);
  }
}
