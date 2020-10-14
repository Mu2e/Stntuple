#ifndef __Stntuple_ana_StnVDetData_t__
#define __Stntuple_ana_StnVDetData_t__
//------------------------------------------------------------------------------
//  IXLocal and IZLocal are either (1,3) or (3,1) 
//  vector normal to the VD is pointed along IZLocal
//------------------------------------------------------------------------------
  struct StnVDetData_t {
    int    fID;
    int    fIXLocal;			// orientation of the local axis to calculate Pt wrt to the field
    int    fIZLocal;
    double fX      ;		        // for X dist
    double fZ      ;
    double fBField ;                    // main component, orthogonal to the detector plane
    double fPhiXZ  ;
  };

#endif
