////////////////////////////////////////////////////////////////////////////////
// verbose = 0: default
//         = 1: fill histogram
//         = 2: print value
//         > 2: more detailed printout
////////////////////////////////////////////////////////////////////////////////
#include "Stntuple/stat/pgaus_t.hh"
namespace stntuple {

//-----------------------------------------------------------------------------
  pgaus_t::pgaus_t(const char* Name, double Mean, double Sigma, int Debug) : parameter_t(Name,Debug) {
    int nbins(1000);
    
    fMean  = Mean ;
    fSigma = Sigma;
    
    // double bin = (Mean+10*Sigma)/nbins;
    fHistPDF = new TH1D(Form("h_par_%s",Name),"parameter",nbins,0,Mean+10*Sigma);
  }

//-----------------------------------------------------------------------------
// initialize the parameter value just once per pseudoexperiment
//-----------------------------------------------------------------------------
  void pgaus_t::InitValue() {

    if (fFixed == 0) {
      do { fValue = fRng->Gaus(fMean,fSigma); } while (fValue < 0);
    }
    else             fValue = fMean;

    if (fDebug > 0) fHistPDF->Fill(fValue);
  }

  double pgaus_t::XMin() {
    // for histogramming
    return 0;
  }

  double pgaus_t::XMax() {
    // for histogramming
    return fMean + 10*fSigma;
  }
//-----------------------------------------------------------------------------
  void pgaus_t::Print(const Option_t* Opt) const {
  }

}
