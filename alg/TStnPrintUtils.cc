///////////////////////////////////////////////////////////////////////////////
// 2026-09-08: P.Murat 
///////////////////////////////////////////////////////////////////////////////
#include "Stntuple/alg/TStnPrintUtils.hh"
#include "Stntuple/obj/TStnTimeClusterBlock.hh"
#include "Stntuple/obj/TComboHitBlock.hh"

//-----------------------------------------------------------------------------
TStnPrintUtils::TStnPrintUtils() {
}

//-----------------------------------------------------------------------------
// static 'instance' will be deleted in the end of the job
//-----------------------------------------------------------------------------
TStnPrintUtils* TStnPrintUtils::Instance() {
  static TStnPrintUtils* instance{nullptr};
  if (instance == nullptr) {
    instance = new TStnPrintUtils;
  }
  return instance;
}

//-----------------------------------------------------------------------------
void TStnPrintUtils::PrintTimeClusterBlock(TStnTimeClusterBlock* TcB,
                                           TComboHitBlock*       ChB) {
  
  int banner_printed = 0;

  int nchtot = ChB->NHits();

  // std::cout << std::format("nchtot:{} GetEntries:{}\n",nchtot,ChB->ListOfHits()->GetEntries());
  
  for (int i=0; i<TcB->fNTimeClusters; i++) {
    TStnTimeCluster* tc = TcB->TimeCluster(i);
    if (! banner_printed) {
      tc->Print("banner");
      banner_printed = 1;
    }
    tc->Print("data");
    
    // afer that try to print combo hits
    int nch = TcB->ListOfChLinks()->NLinks(i);
    bool ch_banner_printed(false);
    
    for (int l=0; l<nch; l++) {
      uint32_t ind = (uint32_t) TcB->ListOfChLinks()->Index(i,l);
      // search for the hit in ChB
      for (int k=0; k<nchtot; k++) {
        TComboHit* ch = ChB->Hit(k);
        if (ind == ch->GetUniqueID()) {
          // hit found, print it
          if (not ch_banner_printed) {
            ch->Print("banner");
            ch_banner_printed = true;
          }
          ch->Print("data");
        }
      }
      
    }
  }
  
};
