///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////

#include "Stntuple/mod/InitCrvClusterBlock.hh"
#include "Offline/RecoDataProducts/inc/CrvRecoPulse.hh"
#include "Offline/RecoDataProducts/inc/CrvCoincidence.hh"
#include "Offline/RecoDataProducts/inc/CrvCoincidenceCluster.hh"

//-----------------------------------------------------------------------------
// in this case AbsEvent is just not used
//-----------------------------------------------------------------------------
int StntupleInitCrvClusterBlock::InitDataBlock(TStnDataBlock* Block, AbsEvent* Event, int Mode) {

  const int ev = Event->event();
  const int rn = Event->run();
  const int sr = Event->subRun();

  if (Block->Initialized(ev,rn,sr)) return 0;

  TCrvClusterBlock* block = (TCrvClusterBlock*) Block;

  block->f_EventNumber  = ev;
  block->f_RunNumber    = rn;
  block->f_SubrunNumber = sr;
//-----------------------------------------------------------------------------
// initialize pointer to the pulse collection
//-----------------------------------------------------------------------------
  art::Handle<mu2e::CrvRecoPulseCollection> crpch;
  const mu2e::CrvRecoPulseCollection*       crpc(nullptr);
  int   ncrp(0);

  if(!fCrvRecoPulseCollTag.empty()) {
    if (Event->getByLabel(fCrvRecoPulseCollTag,crpch)) {
      crpc = crpch.product();
      ncrp = crpc->size();
    } else {
      printf("InitCrvClusterBlock::%s: No CRV pulse collection found\n", __func__);
    }
  }

  const mu2e::CrvRecoPulse* p0 = (ncrp > 0) ? &crpc->at(0) : nullptr;
//-----------------------------------------------------------------------------
// store CrvCoincidenceCluster's
//-----------------------------------------------------------------------------
  art::Handle<mu2e::CrvCoincidenceClusterCollection> cccch;
  const mu2e::CrvCoincidenceClusterCollection*       cccc(nullptr);
  int                                                nccc(0);

  if (!fCrvCoincidenceClusterCollTag.empty()) {
    if (Event->getByLabel(fCrvCoincidenceClusterCollTag,cccch)) {
      cccc = cccch.product();
      nccc = cccc->size();
    } else {
      printf("InitCrvClusterBlock::%s: No CRV coincidence cluster collection found\n", __func__);
    }
  }

  for (int iccc=0; iccc<nccc; iccc++) {
    const mu2e::CrvCoincidenceCluster* cluster = &cccc->at(iccc);

    TCrvCoincidenceCluster* ccc = block->NewCluster();

    const std::vector<art::Ptr<mu2e::CrvRecoPulse>>* list_of_pulses = &cluster->GetCrvRecoPulses();

    const int    sector = cluster->GetCrvSectorType();
    const int    np     = list_of_pulses->size();
    const int    npe    = cluster->GetPEs();

    const double x      = cluster->GetAvgHitPos().x();
    const double y      = cluster->GetAvgHitPos().y();
    const double z      = cluster->GetAvgHitPos().z();
    const float  t1     = cluster->GetStartTime();
    const float  t2     = cluster->GetEndTime();

    ccc->Set(iccc,sector,np,npe,x,y,z,t1,t2);
//-----------------------------------------------------------------------------
// now store pulses associated with the cluster
//-----------------------------------------------------------------------------
    for (int ip=0; ip<np; ip++) {
      const mu2e::CrvRecoPulse* pulse = list_of_pulses->at(ip).get();
      const int index = pulse-p0;
//-----------------------------------------------------------------------------
// check if the pulse is already stored
//-----------------------------------------------------------------------------
      int loc(-1);
      const int npulses = block->NPulses();

      for (int i=0; i<npulses; i++) {
        TCrvRecoPulse* p = block->Pulse(i);
        if (p->Index() == index) {
          loc   = i;
          break;
        }
      }
      if (loc == -1) {
//-----------------------------------------------------------------------------
// add pulse to the list
//-----------------------------------------------------------------------------
        loc                  = npulses;

        TCrvRecoPulse* new_pulse = block->NewPulse();

        const int   npes           = pulse->GetPEs();
        const int   npes_height    = pulse->GetPEsPulseHeight();
        const int   nind           = pulse->GetWaveformIndices().size();
        const int   bar            = pulse->GetScintillatorBarIndex().asInt();
        const int   sipm           = pulse->GetSiPMNumber();

        const float time           = pulse->GetPulseTime();
        const float height         = pulse->GetPulseHeight();
        const float width          = pulse->GetPulseBeta(); // was GetPulseWidth();
        const float chi2           = pulse->GetPulseFitChi2();
        const float le_time        = pulse->GetLEtime();

        new_pulse->Set(index,npes,npes_height,nind,bar,sipm,time,height,width,chi2,le_time);
      }
      block->fClusterPulseLinks->Add(iccc,loc);
    }
  }
  return 0;
}


//-----------------------------------------------------------------------------
// keep this function as an example, don't really need it
//-----------------------------------------------------------------------------
int StntupleInitCrvClusterBlock::ResolveLinks(TStnDataBlock* Block, AbsEvent* Event, int Mode) {

  const int ev = Event->event();
  const int rn = Event->run();
  const int sr = Event->subRun();

  if (! Block->Initialized(ev,rn,sr)) return -1;

  return 0;
}
