//

#include "Stntuple/geom/THexCrystalMap.hh"
#include "Stntuple/geom/TStnCrystal.hh"
#include "TMath.h"
#include "TVector2.h"

TDiskIndex THexCrystalMap::fgPos[6] = {
  TDiskIndex( 1, 0),  // up right
  TDiskIndex( 0, 1),  // down right
  TDiskIndex(-1, 1),  // down
  TDiskIndex(-1, 0),  // down left
  TDiskIndex( 0,-1),  // up left
  TDiskIndex( 1,-1)   // up 
};




//-----------------------------------------------------------------------------
THexCrystalMap::THexCrystalMap() {
}


//-----------------------------------------------------------------------------
THexCrystalMap::~THexCrystalMap() {
}


//-----------------------------------------------------------------------------
THexCrystalMap::THexCrystalMap(double Size, double RMin, double RMax): 
  TDiskCrystalMap(Size,RMin,RMax) 
{
  int         loc;
  TDiskIndex* index;
  TVector2    pos;
  double      fraction;

  for (int ir=0; ir<fNRings; ir++) {
    int ncr = GetNCrystalsPerRing(ir);
    for (int ic=0; ic<ncr; ic++) {
      loc       = GetFirst(ir)+ic;
      index = GetDiskIndex(loc);
      //      ring      = GetRing  (&index);
      int inside    = InsideCode(index,&fraction);
	
      if (inside == 0) {
	GetPosition(index,&pos);
	// *FIXME*	crystal = new TStnCrystal(&index,pos.X(),pos.Y(),Z0,fNEdges,fSize);
	// *FIXME*	crystal->SetDisk(this);
	// *FIXME*	fListOfCrystals->Add(crystal);
	fNInside[ir] += 1;
      }
    }
    
    fNCrystals += fNInside[ir];
  }

}


//-----------------------------------------------------------------------------
// void THexCrystalMap::GetDiskIndex(int I, TDiskIndex* Index) const {

//   int      ring, loc, l(0), k(0), n1, n2;

//   if (I > 0) {
// 				// find the ring number
//     for (int i=1; i>0; i++) {
//       n1 = 3*(i-1)*i+1;
//       n2 = 3*i*(i+1)+1;
//       if (I < n2) {
// 	ring = i;
// 	loc  = I-n1;
// 	break;
//       }
//     }

//     int seg1 = loc / ring; 
//     int pos  = loc % ring;

//     int seg2 = seg1+1;
//     if (seg2 == 6) seg2 = 0;

//     l = THexCrystalMap::fgPos[seg1].fL*ring + (THexCrystalMap::fgPos[seg2].fL-THexCrystalMap::fgPos[seg1].fL)*pos;
//     k = THexCrystalMap::fgPos[seg1].fK*ring + (THexCrystalMap::fgPos[seg2].fK-THexCrystalMap::fgPos[seg1].fK)*pos;
//   }
//   Index->Set(l,k);
// }


//-----------------------------------------------------------------------------
int THexCrystalMap::GetRing(TDiskIndex* Index) {

  int  ring;

  if (Index->fL*Index->fK > 0) {
    ring = std::abs(Index->fL+Index->fK);
  }
  else if ( std::abs(Index->fL) > std::abs(Index->fK) ) {
    ring = std::abs(Index->fL);
  }
  else {
    ring = std::abs(Index->fK);
  }

  return ring;
}

//-----------------------------------------------------------------------------
int THexCrystalMap::GetRing(int I) {
  TDiskIndex*  index;
  
  index = GetDiskIndex(I);

  return GetRing(index);
}

//-----------------------------------------------------------------------------
void THexCrystalMap::GetPosition(TDiskIndex* Index, TVector2* Pos) {
  double x, y, step;

  step = fSize;
  x    = step*(Index->fL+Index->fK)*sqrt(3.)/2.;
  y    = step*(Index->fL-Index->fK)/2.;
  Pos->Set(x,y);
}

//-----------------------------------------------------------------------------
void THexCrystalMap::GetPosition(int I, TVector2* Pos) {
  TDiskIndex*  index;
  
  index = GetDiskIndex(I);
  GetPosition(index,Pos);
}

//-----------------------------------------------------------------------------
double THexCrystalMap::GetRadius(TDiskIndex* Index) {
  double x, y, step, r;

  step = fSize;
  x    = step*(Index->fL+Index->fK)*sqrt(3.)/2.;
  y    = step*(Index->fL-Index->fK)/2.;
  r    = sqrt(x*x+y*y);

  return r;
}

//-----------------------------------------------------------------------------
double THexCrystalMap::GetRadius(int I) {
  TDiskIndex*  index;
  
  index = GetDiskIndex(I);
  return GetRadius(index);
}

//-----------------------------------------------------------------------------
int THexCrystalMap::InsideCode(TDiskIndex* Index, double* Fraction) {

  int       inside, nvin(0), nbelow(0), nabove(0);
  double    x0, y0, r, x, y, phi, s, s0, /*s1,*/ r0, dr, adr;
  double    fr, step;
					// crystal center position
  step = fSize;
  x0   = step*(Index->fL+Index->fK)*sqrt(3.)/2.;
  y0   = step*(Index->fL-Index->fK)/2.;

					// loop over 6 vertices and check if 
					// all of them are inside
  double rho = fSize/sqrt(3.);

  for (int i=0; i<6; i++) {
    phi = i*TMath::Pi()/3;
    x   = x0+rho*TMath::Cos(phi);
    y   = y0+rho*TMath::Sin(phi);

    r   = sqrt(x*x+y*y);

    if (r < fRMin) {
      nbelow += 1;
    }
    else if (r > fRMax) {
      nabove += 1;
    }
    else {
      nvin += 1;
    }
  }
//-----------------------------------------------------------------------------
// 'nvin' - number of vertices inside the disk annulus
//-----------------------------------------------------------------------------
  if (nvin == 0) {
    fr = 0;
  }
  else if (nvin == 6) {
    fr = 1.;
  }
  else { 
//-----------------------------------------------------------------------------
// at this point can calculate fraction of the crystal area inside the ring
//-----------------------------------------------------------------------------
    r0  = sqrt(x0*x0+y0*y0);
    dr  = r0-fRMin;

    if (nabove > 0) {
//-----------------------------------------------------------------------------
// crystal crosses the outer ring
//-----------------------------------------------------------------------------
      dr  = r0 -fRMax;
      adr = fabs(dr);
      if (adr > fSize/2) adr = fSize/2.;

      s   = fSize*fSize*sqrt(3)/2;
      //      s1  = (2*fSize-adr)*adr/sqrt(3.);
      s0  = (3*fSize-2*adr)*(fSize-2*adr)/4/sqrt(3);

      if (dr <= 0) {
	fr = 1-s0/s;
      }
      else {
        fr = s0/s;
      }
    }
    else {
//-----------------------------------------------------------------------------
// crystal crosses the inner ring
//-----------------------------------------------------------------------------
      adr = fabs(dr);
      if (adr > fSize/2) adr = fSize/2.;

      s   = fSize*fSize*sqrt(3)/2;
      //      s1  = (2*fSize-adr)*adr/sqrt(3.);
      s0  = (3*fSize-2*adr)*(fSize-2*adr)/4/sqrt(3);
      if (dr > 0) {
	fr = 1-s0/s;
      }
      else {
        fr = s0/s;
      }
    }
  }

  if (fr >= fMinFraction) inside = 1;
  else                    inside = 0;

  *Fraction = fr;

  return inside;
}

