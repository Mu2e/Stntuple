///////////////////////////////////////////////////////////////////////////////
// PM: this include is temporary and it will go away as soon
// as the DB-based approach is implemented
// in essence, it is a table prototype
///////////////////////////////////////////////////////////////////////////////
#ifndef __daqana_mod_TCaloChannelMap_hh__
#define __daqana_mod_TCaloChannelMap_hh__

#include<array>

class TCaloChannelMap {
public:
  enum {
    kNDisks        = 2,
    kNCrates       = 10,
    kNBoards       = 20,  // per crate
    kNChPerBoard   = 20,  // guesswork
    kMaxNChannels  = kNDisks*kNCrates*kNBoards*kNChPerBoard,
  };
  
  struct Data_t {
    int  channel;
    int  disk;                           // offline channel
    int  crate;
    int  board;
    int  cid;                            // crystal ID
    int  sipm;
  };
                                        // there are 6000 smth bars < 25000 channels
  std::vector<Data_t> _data;
  
  int      fRunNumber;
  int      fNChannels;
  
  Data_t* _ch_data_by_offline   [kMaxNChannels];
  Data_t* _ch_data_by_online    [kNDisks][kNCrates][kNBoards][kNChPerBoard];  //

private:
  TCaloChannelMap();
  ~TCaloChannelMap();

public:
  static TCaloChannelMap* Instance();

  int Init(int RunNumber);

  Data_t* ch_data_by_online (int disk, int crate, int board, int ch) {
    return _ch_data_by_online [disk][crate][board][ch];
  }
  
  Data_t* ch_data_by_offline(int sipmid){ return _ch_data_by_offline[sipmid]; }

};

#endif
