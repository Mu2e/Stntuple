//-----------------------------------------------------------------------------
// file: Stntuple/mod/scripts/display_001.C
//
// test macro : fill and display histogram of N(combo hits) in the event
// no includes needed
// hide variables - dos it have any effect ?
//-----------------------------------------------------------------------------
namespace {
  TH1F*               hist   (nullptr);	       // to be accessed interactively
  TCanvas*            c      (nullptr);
}

TTrkVisNode* tnode(nullptr);

//-----------------------------------------------------------------------------
// Mode = 0: begin job (run)
// Mode = 1: event
// Mode = 2: end job (run)
//-----------------------------------------------------------------------------
void display_001(int Mode, TModule* Module) {
  printf("display_001 called: Mode = %i, Module = %p\n",Mode,Module);

  // TStnVisManager* vm = TStnVisManager::Instance();
  
  //  m_disp = (mu2e::MuHitDisplay*) Module;

  if (Mode == 0) {  
//-----------------------------------------------------------------------------
// book histograms , perform other initializations
//-----------------------------------------------------------------------------
    c = new TCanvas;
    c->cd();
    hist    = new TH1F ("hist","hist",200,0,200);
    hist->Draw();
  }
  else if (Mode == 1) {
//-----------------------------------------------------------------------------
// fill histograms 
//-----------------------------------------------------------------------------
    if (tnode == nullptr) {
      tnode = (TTrkVisNode*) TStnVisManager::Instance()->FindNode("TrkVisNode");
    }

    if (tnode) {
                 
      const mu2e::ComboHitCollection* chc = tnode->GetCComboHitColl();  // const

      int nh = chc->size();

      hist->Fill(nh);

      printf("nh:  %i\n",nh);
    }
    else {
      std::cout << "ERROR: tnode == nullptr" << std::endl;
    }

    c->Modified();
    c->Update();
  }
  else if (Mode == 2) {
//-----------------------------------------------------------------------------
// end run
//-----------------------------------------------------------------------------
  }
}
