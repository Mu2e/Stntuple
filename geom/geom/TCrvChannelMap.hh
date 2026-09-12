///////////////////////////////////////////////////////////////////////////////
// PM: this include is temporary and it will go away as soon
// as the DB-based approach is implemented
// in essence, it is a table prototype
///////////////////////////////////////////////////////////////////////////////
#ifndef __Stntuple_geom_TCrvChannelMap_hh__
#define __Stntuple_geom_TCrvChannelMap_hh__

#include<deque>

class TCrvChannelMap {
public:
  enum {
    kNRocs           = 18,
    kNFebs           = 24,
    kNChannelsFeb    = 64,
    kNChannelsMax    = kNRocs*kNFebs*kNChannelsFeb,
  };
  
  struct Data_t {
    int  och;                           // offline channel
    int  roc;
    int  feb;
    int  fch;                           // FEB channel
  };

  static TCrvChannelMap* fgInstance;
                                        // there are 6000 smth bars < 25000 channels
  
  std::deque<Data_t> _data;
  
  int      fRunNumber;
  int      fNChannels;
  
  Data_t* _ch_data_by_offline   [kNChannelsMax];     // [mnid]
  Data_t* _ch_data_by_online    [kNRocs][kNFebs][kNChannelsFeb];  //

private:
  TCrvChannelMap(int RunNumber);
  ~TCrvChannelMap();

public:
  static TCrvChannelMap* Instance(int RunNumber);

  Data_t* ch_data_by_online (int roc, int feb, int fch) { return _ch_data_by_online [roc][feb][fch]; }
  Data_t* ch_data_by_offline(int och)                   { return _ch_data_by_offline[och]; }

};

#endif
