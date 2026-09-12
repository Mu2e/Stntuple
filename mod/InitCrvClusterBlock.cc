///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
#include <format>
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "Stntuple/mod/InitCrvClusterBlock.hh"
#include "Stntuple/obj/TCrvPulseBlock.hh"
#include "Stntuple/obj/TStnEvent.hh"
#include "Offline/RecoDataProducts/inc/CrvRecoPulse.hh"
#include "Offline/RecoDataProducts/inc/CrvCoincidenceCluster.hh"
#include "Offline/MCDataProducts/inc/CrvCoincidenceClusterMC.hh"
#include "Offline/MCDataProducts/inc/CrvCoincidenceClusterMCAssns.hh"

//-----------------------------------------------------------------------------
// in this case AbsEvent is just not used
//-----------------------------------------------------------------------------
int StntupleInitCrvClusterBlock::InitDataBlock(TStnDataBlock* Block, AbsEvent* Event, int Mode) {
  const char oname []  = {"stntuple::InitCrvClusterBlock::InitDataBlock"};

  const int ev = Event->event();
  const int rn = Event->run();
  const int sr = Event->subRun();

  if (Block->Initialized(ev,rn,sr)) return 0;

  const int verbose(0);

  TCrvClusterBlock* block = (TCrvClusterBlock*) Block;

  block->f_EventNumber  = ev;
  block->f_RunNumber    = rn;
  block->f_SubrunNumber = sr;
//-----------------------------------------------------------------------------
// initialize pointer to the pulse collection
//-----------------------------------------------------------------------------
  art::Handle<mu2e::CrvRecoPulseCollection> crpch;

  if (!fCrvRecoPulseCollTag.empty()) {
    if (not Event->getByLabel(fCrvRecoPulseCollTag,crpch)) {
      mf::LogWarning(oname) << std::format("WARNING: InitCrvClusterBlock::{} No CRV pulse collection (%s) found",
                                           __func__, fCrvRecoPulseCollTag.encode().c_str());
    }
  }
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
    }
    else {
      mf::LogWarning(oname) << std::format("WARNING: No CRV coincidence cluster collection (%s) found",
                                           fCrvCoincidenceClusterCollTag.encode().c_str());
    }
  }

  // Retrieve MC information if it's available
  art::Handle<mu2e::CrvCoincidenceClusterMCAssns> mc_ccc_assnsH;
  const mu2e::CrvCoincidenceClusterMCAssns*       mc_ccc_assns(nullptr);
  if (!fCrvCoincidenceClusterMCCollTag.empty()) {
    if (Event->getByLabel(fCrvCoincidenceClusterMCCollTag,mc_ccc_assnsH)) mc_ccc_assns = mc_ccc_assnsH.product();
    else {
      mf::LogWarning(oname) << std::format("WARNING: No MC <--> Reco CRV coincidence cluster associations (%s) found",
                                           fCrvCoincidenceClusterMCCollTag.encode().c_str());
    }
  }
  
  const int nmc_ccc_assns = (mc_ccc_assns) ? mc_ccc_assns->size() : 0;
  if(mc_ccc_assns && nmc_ccc_assns != nccc) {
    mf::LogWarning(oname) << std::format("WARNING: MC cluster associations (%i) and Reco clusters (%i) don't match",
                                         nmc_ccc_assns, nccc);
  }
//-----------------------------------------------------------------------------
// Loop over the CRV CC
//-----------------------------------------------------------------------------
  for (int iccc=0; iccc<nccc; iccc++) {
    if (verbose > 0) {
      printf("InitCrvClusterBlock::%s: Processing cluster %i\n", __func__, iccc);
    }
    const mu2e::CrvCoincidenceCluster* cluster = &cccc->at(iccc);
    const mu2e::CrvCoincidenceClusterMC* mc_cluster = nullptr;
    if (!cluster) {
      printf("InitCrvClusterBlock::%s: Cluster %i is not defined!\n", __func__, iccc);
      continue;
    }
    
    if (mc_ccc_assns) {
      for(int iassn = 0; iassn < nmc_ccc_assns; ++iassn) {
        const auto assn = mc_ccc_assns->at(iassn);
        if (!assn.second) {
          printf("InitCrvClusterBlock::%s: Cluster %2i: MC association %i is not valid! Reco = %o, MC = %o\n", __func__, iccc, iassn,
                 assn.first.isAvailable(), assn.second.isAvailable());
          continue;
        }
        const mu2e::CrvCoincidenceCluster*   ireco = (!assn.first) ? nullptr : &(*(assn.first));
        const mu2e::CrvCoincidenceClusterMC* imc   = &(*(assn.second));
        if (ireco && &(*(ireco)) == &(*cluster)) {
          if (verbose > 1) printf(" --> Associated MC cluster found, association index %i\n", iassn);
          mc_cluster = imc;
          break;
        }
      }
      
      if(!mc_cluster && nmc_ccc_assns == nccc) {
        // try by index if the numbers match
        const mu2e::CrvCoincidenceClusterMC* imc = &(*(mc_ccc_assns->at(iccc).second));
        mc_cluster = imc;
        if (verbose > 0) {
          printf("%s:%s: Cluster %2i: Associated MC cluster found by index %i\n",
                 typeid(*this).name(), __func__, iccc, iccc);
        }
      }
      if (!mc_cluster) {
        printf("%s::%s: Cluster %2i: Associated MC cluster not found! N(clusters) = %i N(Assns) = %i\n",
               typeid(*this).name(), __func__, iccc, nccc, nmc_ccc_assns);
      }
    }

    TCrvCoincidenceCluster* ccc = block->NewCluster();
    ccc->SetOfflineCrvc(cluster);

    const std::vector<art::Ptr<mu2e::CrvRecoPulse>>* list_of_pulses = &cluster->GetCrvRecoPulses();

    int    sector = cluster->GetCrvSectorType();
    int    np     = list_of_pulses->size();
    float  pes    = cluster->GetPEs();
    float  slope  = cluster->GetSlope();

    double x      = cluster->GetAvgHitPos().x();
    double y      = cluster->GetAvgHitPos().y();
    double z      = cluster->GetAvgHitPos().z();
    float  t1     = cluster->GetStartTime();
    float  t2     = cluster->GetEndTime();

    ccc->Set(iccc,sector,np,pes,x,y,z,t1,t2, slope);
    if (verbose > 1) {
      printf("  Cluster sector = %2i, N(pulses) = %2i, N(PE) = %4.1f, x = %7.1f, y = %7.1f, z = %8.1f, t1 = %6.1f, t2 = %6.1f, slope = %6.2f, mc_found = %o\n",
             sector, np, pes, x, y, z, t1, t2, slope, mc_cluster != nullptr);
    }
    if (mc_cluster) {
      auto sim = mc_cluster->GetMostLikelySimParticle();
      int   sim_id  = (sim) ? sim->id().asInt() : -1;
      int   mc_np   = mc_cluster->GetPulses().size();

      float mc_edep = mc_cluster->GetTotalEnergyDeposited();
      float mc_time = mc_cluster->GetAvgHitTime();
      float mc_x    = mc_cluster->GetAvgHitPos().x();
      float mc_y    = mc_cluster->GetAvgHitPos().y();
      float mc_z    = mc_cluster->GetAvgHitPos().z();

      ccc->SetMC(sim_id, mc_np, mc_edep, mc_time, mc_x, mc_y, mc_z);
      if (verbose > 1) {
        printf("  MC Info: SIM ID = %i, N(pulses) = %2i, E(dep) = %4.1f, x = %7.1f, y = %7.1f, z = %8.1f, tavg = %6.1f\n",
               sim_id, mc_np, mc_edep, mc_x, mc_y, mc_z, mc_time);
      }
    }
// //-----------------------------------------------------------------------------
// // now store pulses associated with the cluster
// //-----------------------------------------------------------------------------
//     if(fStorePulses) {
//       for (int ip=0; ip<np; ip++) {
//         const mu2e::CrvRecoPulse* pulse = list_of_pulses->at(ip).get();
//         const int index = pulse-p0;
// //-----------------------------------------------------------------------------
// // check if the pulse is already stored
// //-----------------------------------------------------------------------------
//         int loc(-1);
//         const int npulses = block->NPulses();

//         for (int i=0; i<npulses; i++) {
//           TCrvRecoPulse* p = block->Pulse(i);
//           if (p->Index() == index) {
//             loc   = i;
//             break;
//           }
//         }
//         if (loc == -1) {
// //-----------------------------------------------------------------------------
// // add pulse to the list
// //-----------------------------------------------------------------------------
//           loc                  = npulses;

//           TCrvRecoPulse* new_pulse = block->NewPulse();

//           const int   npes           = pulse->GetPEs();
//           const int   npes_height    = pulse->GetPEsPulseHeight();
//           const int   nind           = pulse->GetWaveformIndices().size();
//           const int   bar            = pulse->GetScintillatorBarIndex().asInt();
//           const int   sipm           = pulse->GetSiPMNumber();

//           const float time           = pulse->GetPulseTime();
//           const float height         = pulse->GetPulseHeight();
//           const float width          = pulse->GetPulseBeta(); // was GetPulseWidth();
//           const float chi2           = pulse->GetPulseFitChi2();
//           const float le_time        = pulse->GetLEtime();

//           new_pulse->Set(index,npes,npes_height,nind,bar,sipm,time,height,width,chi2,le_time);
//         }
//         block->fClusterPulseLinks->Add(iccc,loc);
//       }
//     }
  }
  if (verbose > 2) block->Print();
  return 0;
}


//-----------------------------------------------------------------------------
// keep this function as an example, don't really need it
//-----------------------------------------------------------------------------
int StntupleInitCrvClusterBlock::ResolveLinks(TStnDataBlock* Block, AbsEvent* Event, int Mode) {

  const int evn = Event->event();
  const int rn  = Event->run();
  const int srn = Event->subRun();

  if (! Block->Initialized(evn,rn,srn)) return -1;
  if (  Block->LinksInitialized()     ) return  0;
  
  TCrvClusterBlock* crvc_block = (TCrvClusterBlock*) Block;

  // TStnEvent* ev   = crvc_block->GetEvent();
  //  auto crvp_block = (TCrvPulseBlock*) ev->GetDataBlock("CrvpBlock");

  int ncrvc = crvc_block->NClusters();
  for (int i=0; i<ncrvc; i++) {
    auto crvc = crvc_block->Cluster(i);
    // this one has a list of pointers to crv reco pulses
    // rely on that the Mu2e o_crvc_c is a vector of things
    const mu2e::CrvCoincidenceCluster* o_crvc = crvc->OfflineCrvc();
    auto o_crvp_coll = &o_crvc->GetCrvRecoPulses();  // offline list of associated pulses
    auto crvp_0      = &o_crvp_coll->at(0);               // pointer to the first pulse
    int ncrvp        = o_crvp_coll->size();
    for (int j=0; j<ncrvp; j++) {
      auto crvp_j = &o_crvp_coll->at(j);
      int  pulse_index = crvp_j-crvp_0;
      crvc_block->ClusterToPulseLinks()->Add(i,pulse_index);
    }
  }
//-----------------------------------------------------------------------------
// mark links as initialized
//-----------------------------------------------------------------------------
  crvc_block->fLinksInitialized = 1;

  return 0;
}
