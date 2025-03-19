//
#include "TEnv.h"
#include "Stntuple/geom/TCrvNumerology.hh"


TCrvNumerology*        TCrvNumerology::fgInstance       = 0;

int sector_type[22] = {
  1, 1, 1, 1, 1, 2, 2, 3, 3, 3,
  3, 4, 4, 5, 6, 6, 6, 6, 7, 8,
  9, 9
};
//-----------------------------------------------------------------------------
TCrvNumerology::TCrvNumerology() {

  const char* gdata{"Stntuple/geom/data"};

  const char* dir  = gEnv->GetValue("Stntuple.GeometryData",gdata);
//-----------------------------------------------------------------------------
// read file and construct sectors
//-----------------------------------------------------------------------------
  TString fn = Form("%s/crv_sectors.txt",dir);
  FILE* f  = fopen(fn.Data(),"r");
  if (f == 0) {
    Error("Init",Form("missing file %s\n",fn.Data()));
    // return -2;
  }

  int     sector, iwy, iwx, iwz ;
  float   dy, dx, dz;
  char    name[100];
  int done = 0;
  char c[10000];

  // fNChannels = 0;

  while ( ((c[0]=getc(f)) != EOF) && !done) {
                                        // check if it is a comment line
    if (c[0] != '#') {
      ungetc(c[0],f);
                                        // parse line
      fscanf(f,"%i" ,&sector       );
      fscanf(f,"%s" ,name          );

      SectorData_t* sec = &fSector[sector];
      sec->fNumber = sector;
      sec->fName   = name;

      fscanf(f,"%i" ,&sec->fNModules);
      fscanf(f,"%i" ,&sec->fNLayers );
      fscanf(f,"%i" ,&sec->fNBarsPerLayer);
      fscanf(f,"%i" ,&sec->fFirstIndex);
      fscanf(f,"%i" ,&iwy          );
      fscanf(f,"%i" ,&iwx          );
      fscanf(f,"%i" ,&iwz          );
      fscanf(f,"%f" ,&dy           );
      fscanf(f,"%f" ,&dx           );
      fscanf(f,"%f" ,&dz           ); // half-length
      sec->fBarLength = dz*2/100; // convert into meters....
//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------
    }
                                        // skip line
    fgets(c,1000,f);
  }

  fclose(f);
//-----------------------------------------------------------------------------
// now the same way read in the counter geom data
//-----------------------------------------------------------------------------
  TString fn1 = Form("%s/crv_counter_geom.txt",dir);
  FILE* f1  = fopen(fn1.Data(),"r");
  if (f1 == 0) {
    Error("Init",Form("missing file %s\n",fn1.Data()));
    // return -2;
  }

  int     index, module, layer, ibar;
  float   x0, y0, z0;

  //  fNChannels = 0;
  done = 0;
  while ( ((c[0]=getc(f1)) != EOF) && !done) {
                                        // check if it is a comment line
    if (c[0] != '#') {
      ungetc(c[0],f1);
                                        // parse line
      fscanf(f1,"%i" ,&index       );

      BarData_t* bar = &fBar[index];
      bar->fBarIndex = index;

      fscanf(f1,"%i" ,&sector);
      bar->fSector = sector;

      bar->fSectorType = sector_type[sector];

      fscanf(f1,"%i" ,&module );
      fscanf(f1,"%i" ,&layer);
      fscanf(f1,"%i" ,&ibar);
      fscanf(f1,"%i" ,&bar->fBarU );
      fscanf(f1,"%i" ,&bar->fBarV );
      fscanf(f1,"%i" ,&bar->fBarK );
      fscanf(f1,"%f" ,&x0           );
      fscanf(f1,"%f" ,&y0           );
      fscanf(f1,"%f" ,&z0           ); // half-length
      bar->fBarPos.SetXYZ(x0,y0,z0);

      fscanf(f1,"%f" ,&dx           );
      fscanf(f1,"%f" ,&dy           );
      fscanf(f1,"%f" ,&dz           ); // half-length
//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------
    }
                                        // skip line
    fgets(c,1000,f1);
  }

  fclose(f1);
}

//_____________________________________________________________________________
TCrvNumerology::~TCrvNumerology() {
}

//_____________________________________________________________________________
TCrvNumerology*  TCrvNumerology::Instance() {
  static Cleaner cleaner;
  return (fgInstance) ? fgInstance : (fgInstance = new TCrvNumerology());
}

//------------------------------------------------------------------------------
TCrvNumerology::Cleaner::Cleaner() {
}

//------------------------------------------------------------------------------
TCrvNumerology::Cleaner::~Cleaner() {
  if (TCrvNumerology::fgInstance) {
    delete TCrvNumerology::fgInstance;
    TCrvNumerology::fgInstance = 0;
  }
}

//------------------------------------------------------------------------------
int TCrvNumerology::GetBarInfo(int BarIndex, int& Sector, int& Module, int& Layer, int& Bar) {
  Sector = -1;
  Module = -1;
  Layer  = -1;
  Bar    = -1;
  for (int is=0; is<kNSectors; is++) {
    SectorData_t* s = fSector+is;
    int nbars = s->fNModules*s->fNLayers*s->fNBarsPerLayer;
    if ((BarIndex >= s->fFirstIndex) and (BarIndex < s->fFirstIndex+nbars)) {
      // found
      Sector  = is;
      int loc = BarIndex-s->fFirstIndex;
      int nbm = s->fNBarsPerLayer*s->fNLayers;   // n(bars per module)
      Module  = loc/nbm;
      int l2  = loc - nbm*Module;
      Layer   = l2/s->fNBarsPerLayer;
      Bar     = l2 - s->fNBarsPerLayer*Layer;
      break;
    }
  }
  return Sector;
}

//------------------------------------------------------------------------------
void TCrvNumerology::Clear(Option_t* Opt) {
}

//------------------------------------------------------------------------------
void TCrvNumerology::Print(Option_t* Opt) const {
  printf ("N(sectors): %i\n",kNSectors);
  for (int i=0; i<kNSectors; i++) {
    printf("%3i %-10s %5i\n",i,fSector[i].fName.Data(),fSector[i].fFirstIndex);
  }
}
