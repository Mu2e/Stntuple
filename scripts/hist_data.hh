//-----------------------------------------------------------------------------
// a 'dataset' here is a histogram file plus data fields for plotting attributes
//-----------------------------------------------------------------------------
#ifndef __Stntuple_scripts_hist_data_hh__
#define __Stntuple_scripts_hist_data_hh__

#include "stn_book.hh"
#include "hist_file.hh"

class hist_data_t {
public:
  TH1*            fHist;                  // if non-null, fFile is not needed
  hist_file_t*    fFile;                  // it contains a pointer to the dataset
  const stn_book* fBook;

  TString       fLabel  ;		// label to appear on a plot
  float         fLabelXMin;
  float         fLabelYMin;
  float         fLabelFontSize;
  float         fLabelFont;

  TString       fModule;	        // module name
  TString       fName;	                // histogram name, as saved in a file
  TString       fNewName;	        // redefines histogram name at plotting time
  int           fLineColor;             // these are utility fields to be used as needed
  int           fLineWidth;
  int           fFillColor;             // these are utility fields to be used as needed
  int           fFillStyle;             // these are utility fields to be used as needed
  int           fMarkerStyle;
  float         fMarkerSize;
  int           fMarkerColor;
  float         fXMin;
  float         fXMax;
  float         fYMin;
  float         fYMax;
  int           fStats;
  int           fOptStat;
  float         fStatBoxXMin;
  float         fStatBoxYMin;
  float         fStatBoxXMax;
  float         fStatBoxYMax;
  int           fRebin;
  int           fYLogScale;
  TString       fDrawOpt;
  float         fScale;			// 1:normalize to NGenEvents, 2:normalize to integral
  TString       fPlotName;
  TString       fPlotLabel;
  TString       fXAxisTitle;
  TString       fYAxisTitle;
  float         fLegendXMin;
  float         fLegendYMin;
  float         fLegendXMax;
  float         fLegendYMax;
  TString       fCanvasName;
  int           fCanvasSizeX;
  int           fCanvasSizeY;
  TCanvas*      fCanvas;
  TString       fOutputFn;

  hist_data_t(const char* DsID = "", const char* JobName = "", const char* Module = "", const char* HistName = "") {
    fHist        = nullptr;
    if (DsID && (DsID[0] != 0)) fFile = get_hist_file(DsID,JobName);
    else                        fFile = nullptr;
    fBook        = nullptr;
    fName        = HistName;
    fModule      = Module;
    init();
  }
  
  hist_data_t(const stn_book* Book, const char* DsID = "", const char* JobName = "", const char* Module = "", const char* HistName = "") {
    fHist        = nullptr;
    fBook        = Book;
    fFile        = Book->FindHistFile(DsID,"",JobName);
    fName        = HistName;
    fModule      = Module;
    init();
  }
  
  hist_data_t(const stn_catalog* Catalog, const char* BookName, const char* DsID = "",
	      const char* JobName = "",   const char* Module = "", const char* HistName = "") {

    fHist        = nullptr;
    fBook        = Catalog->FindBook(BookName);
    fFile        = fBook->FindHistFile(DsID,"",JobName);
    fName        = HistName;
    fModule      = Module;
    init();
  }
  
  hist_data_t(TH1* Hist, const char* JobName = "", const char* Module = "") {
    fHist        = Hist;
    fFile        = nullptr;            // don't need it, fHist is already defined
    fBook        = nullptr;            // same reason
    fName        = Hist->GetName();
    fModule      = Module;
    init();
  }

  void  init() {
    fCanvasName  = "";
    fCanvasSizeX = 1200.;
    fCanvasSizeY =  800.;
    fNewName     = "";
    fLabel         = "";
    fLabelXMin     = 0.15;
    fLabelYMin     = 0.86;
    fLabelFontSize = 0.03;
    fLabelFont     = 52; // lower left corner
    fLineColor   = -1;
    fLineWidth   =  1;
    fMarkerSize  = -1;
    fMarkerStyle = -1;
    fMarkerColor = -1;
    fFillColor   = -1;
    fFillStyle   = -1;
    fXMin        =  0;
    fXMax        = -1;
    fYMin        =  0;
    fYMax        = -1;
    fStats       =  1;
    fOptStat     = -1;
    fStatBoxXMin =  0.65;
    fStatBoxYMin =  0.65;
    fStatBoxXMax =  0.90;
    fStatBoxYMax =  0.90;
    fRebin       = -1;
    fYLogScale   =  0;
    fDrawOpt     = "";
    fScale       = -1;
    fPlotName    = "";
    fPlotLabel   = "";
    fXAxisTitle  = "";
    fYAxisTitle  = "";
    fLegendXMin  = -1;
    fLegendYMin  = -1;
    fLegendXMax  = -1;
    fLegendYMax  = -1;
  }

  hist_file_t*  get_hist_file(const char* DsID, const char* JobName);
  
};

#endif
