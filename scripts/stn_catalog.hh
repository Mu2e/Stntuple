//-----------------------------------------------------------------------------
// a 'dataset' here is a histogram file plus data fields for plotting attributes
//-----------------------------------------------------------------------------
#ifndef __Stntuple_scripts_stn_catalog_hh__
#define __Stntuple_scripts_stn_catalog_hh__

#include "Stntuple/scripts/stn_book.hh"
#include "Stntuple/scripts/stn_dataset.hh"
#include "Stntuple/scripts/hist_file.hh"
//-----------------------------------------------------------------------------
// book is a list of datasets, list owns its pointers
//-----------------------------------------------------------------------------
class stn_catalog {
public:
  TString              fName;
  vector<stn_book*>    fListOfBooks;

  stn_catalog(const char* Name = "") {
    fName    = Name;
  }
  
  ~stn_catalog() {
    int n = fListOfBooks.size();
    for (int i = 0; i<n; i++) {
      delete (fListOfBooks)[i];
    }
  }

  int AddBook(stn_book* Book) {
    int rc = fListOfBooks.size();
    fListOfBooks.push_back(Book);
    Book->fCatalog = this;
    return rc;
  }

  stn_book* Book(int I) { return fListOfBooks.at(I); }

  stn_book* FindBook(const char* Name) const {
    stn_book* book(nullptr);

    for(int i=0; i<fListOfBooks.size(); i++) {
      stn_book* b = fListOfBooks[i];
      if (b->fName == Name) {
	book = b;
	break;
      }
    }
    if (book == nullptr) {
      printf(" get_book ERROR: cant find book=%s, return NULL\n",Name);
    }
    return book;
  }
};

#endif
