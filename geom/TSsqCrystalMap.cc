//

#include "Stntuple/geom/TSsqCrystalMap.hh"
#include "TVector2.h"
#include "math.h"


int TSsqCrystalMap::fgStep[12] = {
  1, 1,
  0, 1,
  -1, 0,
  -1,-1,
  0,-1, 
  1, 0
};

//-----------------------------------------------------------------------------
TSsqCrystalMap::TSsqCrystalMap():TDiskCrystalMap() {
}


//-----------------------------------------------------------------------------
TSsqCrystalMap::~TSsqCrystalMap() {
}



//-----------------------------------------------------------------------------
TSsqCrystalMap::TSsqCrystalMap(double Size, double RMin, double RMax) : 
  TDiskCrystalMap(Size, RMin, RMax) {

  double fr;				// reserved for future

  TDiskIndex  di;
  
  int n_inside(1); // n(crystals) in a given ring inside the disk, seed the loop

  int nrings = 0;
  int ncrtot = 0;
  //  int icc    = 0;
//-----------------------------------------------------------------------------
// first loop : find number of rings
//-----------------------------------------------------------------------------
  for (int ir=0; n_inside>0; ir++) {
    n_inside = 0;

    //    printf(" ir  = %5i\n",ir);
    if (ir == 0) { 
//-----------------------------------------------------------------------------
// first ring
//-----------------------------------------------------------------------------
      di.fL = 0;
      di.fK = 0;

      int inside = InsideCode(&di,&fr);

      // printf("is = %2i ic=%2i  i1,i2 = %2i,%2i inside=%i\n",
      // 	     0, 0, di.fL,di.fK,inside);

      if (inside > -20) n_inside += 1; // end-loop counter

      if (inside == 0) {
//-----------------------------------------------------------------------------
// add crystal
//-----------------------------------------------------------------------------
	ncrtot++;
      }
    }
    else {
//-----------------------------------------------------------------------------
// rings with ir > 0 have 6*ir crystals each
//-----------------------------------------------------------------------------
      di.fL = 0;
      di.fK = -ir;

      for (int is=0; is<6; is++) {

	for (int ic=0; ic<ir; ic++) {

	  int inside = InsideCode(&di,&fr);
	  
	  // printf("is = %2i ic=%2i  i1,i2 = %2i,%2i inside=%i\n",
	  // 	 is, ic, di.fL,di.fK,inside);

	  if (inside > -20) n_inside += 1;
	  if (inside == 0) ncrtot++;

	  di.fL += fgStep[2*is+0];
	  di.fK += fgStep[2*is+1];
	}
      }
    }
    if (n_inside > 0) nrings ++;
  }

  fNRings   = nrings;
  fNInside  = new int[fNRings];

  int ntotal = GetNTotal();
  fIndex    = new TClonesArray("TDiskIndex",ntotal);
  fIndexMap = new TStnArrayI(ncrtot);
//-----------------------------------------------------------------------------
// at this point, NInsideTotal() is defined - it is the length of fIndexMap
// with nrings calculated, everything is much simpler 
//-----------------------------------------------------------------------------
  int nc = 0;				// running total number of crystals in all rings
  int ic = 0;
  TDiskIndex* pdi;

  for (int ir=0; ir<nrings; ir++) {
    fNInside[ir] = 0;
    pdi = new TDiskIndex(0,0);
    if (ir == 0) {
      int inside = InsideCode(pdi,&fr);
      if (inside == 0) {
	fNInside[ir] += 1;
	new ((*fIndex)[ic]) TDiskIndex(*pdi);
	  //	(*fIndex)[ic] = pdi;
	(*fIndexMap)[nc] = ic;
	nc += 1;
      }
      else {
	(*fIndex)[ic] = NULL;
      }
      ic += 1;
    }
    else {
//-----------------------------------------------------------------------------
// ring with ir>0
//-----------------------------------------------------------------------------
      int i1 = 0;
      int i2 = -ir;

      //      int ncr = 6*ir;			// number of crystals in this ring

      for (int is=0; is<6; is++) { 
					// new side, there are six sides, 'ir' crystals in each
	for (int icr=0; icr<ir; icr++) {
	  pdi = new TDiskIndex(i1,i2);
	  int inside = InsideCode(pdi,&fr);
	  if (inside == 0) {
	    fNInside[ir] += 1;
	    new ((*fIndex)[ic]) TDiskIndex(*pdi);
	      //	    (*fIndex)[ic] = pdi;
	    (*fIndexMap)[nc] = ic;
	    nc += 1;
	  }
	  // else {
	  //   (*fIndex)[ic] = NULL;
	  // }
	  ic += 1;
					// next crystal on the same side
	  i1 += fgStep[2*is+0];
	  i2 += fgStep[2*is+1];
	}
      }
    }
  }
}


//-----------------------------------------------------------------------------
int TSsqCrystalMap::GetRing(int I) {
  return -1;
}

//-----------------------------------------------------------------------------
int TSsqCrystalMap::GetRing(TDiskIndex* Index) {
  return -1;
}


//-----------------------------------------------------------------------------
// 'I' : sequential crystal number in (compressed) IndexMap
//-----------------------------------------------------------------------------
void TSsqCrystalMap::GetPosition(int I, TVector2* Pos) {
  int loc = GetIndexMap(I);
  TDiskIndex* di = GetDiskIndex(loc);
  GetPosition(di,Pos);
}

//-----------------------------------------------------------------------------
void TSsqCrystalMap::GetPosition(TDiskIndex* Index, TVector2* Pos) {
  double x = (Index->fL+Index->fK)/2.*fSize;
  double y = (Index->fL-Index->fK)*fSize;
  Pos->Set(x,y);
}

//-----------------------------------------------------------------------------
double TSsqCrystalMap::GetRadius(TDiskIndex* Index) {
  return -1.;
}

//-----------------------------------------------------------------------------
double TSsqCrystalMap::GetRadius(int I) {
  return -1.;
}

//-----------------------------------------------------------------------------
// square crystals. for partial overlaps, Fraction is undefined
//
// return codes: 
//   0 : crystal fully overlaps with the disk
// -12 : no overlap, crystal fully within the hole
// -22 : no overlap, crystal fully outside
// -11 : partial overlap with the inner edge
// -21 : partial overlap with the outer edge
//-----------------------------------------------------------------------------
int TSsqCrystalMap::InsideCode(TDiskIndex* Index, double* Fraction) {
  double xmin, xmax, ymin, ymax;
					// crystal center position
  double ix = (Index->fL+Index->fK)/2.;
  double iy = (Index->fL-Index->fK);

  double x0 = ix*fSize;
  double y0 = iy*fSize;

  if (x0 <= 0) {
    xmax = x0-fSize/2;
    xmin = x0+fSize/2;
  }
  else {
    xmax = x0+fSize/2;
    xmin = x0-fSize/2;
  }

  if (y0 <= 0) {
    ymax = y0-fSize/2;
    ymin = y0+fSize/2;
  }
  else {
    ymax = y0+fSize/2;
    ymin = y0-fSize/2;
  }

  double rmin = sqrt(xmin*xmin+ymin*ymin);
  double rmax = sqrt(xmax*xmax+ymax*ymax);

  *Fraction = 0;

  if (rmax < fRMin) {
					// crystal entirely in the hole
    return -12;
  }
  else if (rmin > fRMax) {
					// crystal entirely outside
    return -22;
  }
  else {
    if (rmin < fRMin) {
					// partial overlap with the inner edge
      return -11;
    }
    else {
      if (rmax > fRMax) {
					// partial overlap with the outer edge
	return -21;
      }
      else {
					// full overlap
	*Fraction = 1;
	return 0;
      }
    }
  }
}
