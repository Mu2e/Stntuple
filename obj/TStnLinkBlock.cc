//_____________________________________________________________________________
//  base class for STNTUPLE link block
// 
//  Author:    Pasha Murat (CDF/FNAL)
//  Date:      Nov 10 2000
//
//  TStnLinkBlock contains one-directional links, stored as integers, between
//  the pieces of the data stored on different branches, for example, between
//  the electrons and the CES or CPR clusters (see TStnEleLinkBlock). 
//  The number of CPR clusters found in CES-based mode, which match to the 
//  electron number `i' can be found as follows:
//
//  TStnEleLinkBlock* link_block;  // suppose the block has been initialized
//
//  TStnLinkBlock* cpr_links = link_block->CprLinkCes();
//
//  int ncl = cpr_links->NLinks(i);
//
//  and the loop over the matching CPR clusters would look like this:
//
//  for (int j=0; j<ncl; j++) {
//                                      // `icl' is the index of `j'-th CPR
//                                      // cluster matching the `i'-th electron
//    int icl = cpr_links->Index(i,j);
//    ............ 
//  }
//_____________________________________________________________________________
#include "TClass.h"
#include <Stntuple/obj/TStnLinkBlock.hh>

ClassImp(TStnLinkBlock)
//_____________________________________________________________________________
void TStnLinkBlock::Streamer(TBuffer& R__b) {
   // Stream an object of class TStnLinkBlock as fast as possible

  int nw;

  if (R__b.IsReading()) {
    Version_t R__v = R__b.ReadVersion(); if (R__v) { }
    R__b >> nw;
    fLast = nw-2;
    if (fLast >= 0) {

      if (nw > fOffset.fN) {
	fOffset.Set(nw);
      }
      R__b.ReadStaticArray(fOffset.fArray);

      R__b >> fNLinksTotal;
      if (fNLinksTotal > fIndex.fN) {
	fIndex.Set(fNLinksTotal);
      }
      R__b.ReadStaticArray(fIndex.fArray);
      
    } 
    else {
//-----------------------------------------------------------------------------
//  as we're not reading anything, make sure link list looks empty
//-----------------------------------------------------------------------------
      memset(fOffset.fArray,0,fOffset.fN*sizeof(Int_t));
    }
  }
  else {
    R__b.WriteVersion(TStnLinkBlock::IsA());

    nw = fLast+2;
    R__b << nw;
    if (fLast >= 0) {
      R__b.WriteArray(fOffset.fArray,nw);
      R__b << fNLinksTotal;
      R__b.WriteArray(fIndex.fArray,fNLinksTotal);
    }
  }
}


//_____________________________________________________________________________
TStnLinkBlock::TStnLinkBlock() {
  fLast        = -1;
  fNLinksTotal = 0;

  fOffset.Set(10);
  fOffset.fArray[0] = 0;

  fIndex.Set(10);
}


//_____________________________________________________________________________
TStnLinkBlock::~TStnLinkBlock() {
}


//_____________________________________________________________________________
Int_t TStnLinkBlock::InitLinks(Int_t i) {
  // force initialization of links for object `i' (otherwise the last one 
  // may get lost)
  // assume that this routine is never called in non-sequential order

  int loc;

  if (i != fLast) {
    if (fOffset.fN <= i+1) {
      fOffset.Set(2*i+1);
    }
					// don't always expect i = fLast+1, but 
					// `i' is always greater than fLast
    loc = fOffset.fArray[fLast+1];

    for (int j=fLast+2; j<=i+1; j++) {
      fOffset.fArray[j] = loc;
    }
    fLast              = i;
  }
  return 0;
}


//_____________________________________________________________________________
Int_t TStnLinkBlock::Add(Int_t i, Int_t match) {
  // add one more link i --> match, where `i' is an index of an object stored
  // on one branch (electon, for example) and `match' is an index of an object 
  // stored on different branch (matching CES cluster, for example). 
  // It is assumed that link block is filled sequentially, meaning that
  // first all the electron#0 links are stored, then - all the electron#1 
  // links etc...

  int loc;

  if (i != fLast) {
    if (fOffset.fN <= i+1) {
      fOffset.Set(2*i+1);
    }
					// don't always expect i = fLast+1, but 
					// `i' is always greater than fLast
    loc = fOffset.fArray[fLast+1];

    for (int j=fLast+2; j<=i; j++) {
      fOffset.fArray[j] = loc;
    }

    fOffset.fArray[i+1] = fOffset.fArray[i]+1;
    fLast              = i;
  }
  else {
    fOffset.fArray[i+1]++;
  }

  if (fIndex.fN <= fNLinksTotal) {
    fIndex.Set(2*fIndex.fN+1);
  }
  fIndex.fArray[fNLinksTotal++] = match;

  return 0;
}


//_____________________________________________________________________________
void TStnLinkBlock::Print(Option_t* option) const {
  if (fLast < 0) {
    printf("%s->Print(): empty link block\n",GetName());
    return;
  }

  printf( "  i  n(links)   links \n");

  for (int i=0; i<fLast+1; i++) {
    int nl = NLinks(i);
    printf("%3i  %3i    ",i,nl);
    if (nl > 0) {
      printf("   ");
      for (int j=0; j<nl; j++) {
	printf("%5i ",Index(i,j));
      }
    }
    printf("\n");
  }
}

//_____________________________________________________________________________
void TStnLinkBlock::Clear(Option_t* option) {
  fLast             = -1;
  fNLinksTotal      = 0;
  fOffset.fArray[0] = 0;
}

//_____________________________________________________________________________
Int_t TStnLinkBlock::Init(Int_t N) {
  // initialize link block for N objects
  fLast             = -1;
  fNLinksTotal      = 0;

  if (fOffset.fN <= N+1) {
    fOffset.Set(2*N+1);
  }

  for (int i=0; i<=N; i++) {
    fOffset.fArray[i] = 0;
  }
  return 0;
}


