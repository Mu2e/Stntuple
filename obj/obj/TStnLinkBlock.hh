#ifndef STNTUPLE_TStnLinkBlock
#define STNTUPLE_TStnLinkBlock
//-----------------------------------------------------------------------------
//  definition of the STNTUPLE link block (holder for the data block)
//  Author:    Pasha Murat (CDF/FNAL)
//  Date:      Nov 10 2000
// 
//-----------------------------------------------------------------------------
#include "TArrayI.h"
#include "Stntuple/obj/TStnDataBlock.hh"

class TStnLinkBlock: public TStnDataBlock {
public:
//-----------------------------------------------------------------------------
//  data members
//-----------------------------------------------------------------------------
  Int_t            fInitialized;	// ! do we really need it? 
  Int_t            fLast;
  Int_t            fNLinksTotal;	// total number of links
  TArrayI          fOffset;		//
  TArrayI          fIndex;		//
//-----------------------------------------------------------------------------
//  functions
//-----------------------------------------------------------------------------
					// ****** constructors and destructor
  TStnLinkBlock();
  virtual ~TStnLinkBlock();
					// ****** init methods

					// initialize block for N elements
  Int_t   Init(Int_t N);
					// init links for object I
  Int_t   InitLinks(Int_t I);
					// ****** accessors

  Int_t   NLinksTotal() const { return fNLinksTotal; }
  Int_t   Last       () const { return fLast; }

					// number of links corresponding to the
					// element number `i'
  Int_t  NLinks (Int_t i) const {
    if(NLinksTotal() == 0) return 0;    // empty block, ignore the array
    return fOffset.fArray[i+1]-fOffset.fArray[i];
  }
					// index of the `j'-th match to `i'-th
					// element

  Int_t  Index (Int_t i, Int_t j) const { 
    

    return fIndex.fArray[fOffset.fArray[i]+j]; 
  }
					// ****** setters
  Int_t  Add(Int_t i, Int_t match);
					// ****** overloaded functions of 
					// TObject
  void Print(Option_t* option = "") const;
  void Clear(Option_t* option = "");
					// ****** schema evolution
					// nothing has changed I/O-wise between
					// versions 1 and 2, just a member 
					// function has been added and a data 
					// member has been renamed
  ClassDef(TStnLinkBlock,2)
};

#endif
