///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
#include "cmath"
#include "Stntuple/alg/smooth.hh"

// ClassImp(smooth)

//-----------------------------------------------------------------------------
double smooth::Eval(double* X) {

  int    i0, i1;
  double x, x0, x1, dx, dx0, dx1, f0(0), f1(0), f(0), w0(0), w1(0);

  if (fX == 0)   return 0;
//-----------------------------------------------------------------------------
// find the right bin, for the moment don't care about the algorithm speed
//-----------------------------------------------------------------------------
  int nx = fN;
  x  = X[0];
  i0 = -1;
  for (int ix=0; ix<nx; ix++) {
    if (x < fX[ix]) {
      i0 = ix;
      break;
    }
  }
  //  printf("smooth::Eval  x = %12.5e i0 = %i nx = %i\n",x,i0,nx);
//-----------------------------------------------------------------------------
// if X is outside the interpolation range, return zero
//-----------------------------------------------------------------------------
  if ((i0 == 0) || (i0 == -1))                              return 0;
//-----------------------------------------------------------------------------
//  here the regular part starts
//-----------------------------------------------------------------------------
  if (i0 == 1) {
    x0  = fX[i0];
    dx  = x-x0;
    f0  = fP0[i0] + dx*fP1[i0] + dx*dx*fP2[i0];
    w0  = 1;
  }
  else if (i0 == nx-1) {
    i1 = i0-1;
    x1 = fX[i1];
    dx = x-x1;
    f1 = fP0[i1] + dx*fP1[i1] + dx*dx*fP2[i1];
    w1 = 1;
    //   printf("smooth::Eval i1,x1,dx,f1,w1 = %i %12.5e %12.5e %12.5e %12.5e \n",i1,x1,dx,f1,w1);
  }
  else {
    x0  = fX[i0];
    dx0 = x-x0;
    f0  = fP0[i0] + dx0*fP1[i0] + dx0*dx0*fP2[i0];

    i1  = i0-1;
    x1  = fX[i1];
    dx1 = x-x1;
    f1  = fP0[i1] + dx1*fP1[i1] + dx1*dx1*fP2[i1];

    w0 = dx1/(x0-x1);
    w1 = fabs(dx0)/(x0-x1);
  }

  f   = f0*w0+f1*w1;

  //  printf("%10.3f %10.3f %3i %10.3f %3i %10.3f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f\n",
  //	 x,f,i0,x0,i1,x1,fP0[i0],fP1[i0],fP2[i0],dx0,f0,w0,dx1,f1,w1);

  double scale = fFunc->GetParameter(0);
  return f*scale;
}


//-----------------------------------------------------------------------------
// before calling f_smooth (or the corresponding function make sure that 
// smooth::gVar points to the data structure
//-----------------------------------------------------------------------------
double smooth::func(double* X, double* P) {

  double f(0);
					// a hack - P[1] carries the address...
  smooth* o = (smooth*) ((long int)P[1]);

  if (o) f = o->Eval(X);
  else   f = 0;

  return f;
}

//-----------------------------------------------------------------------------
smooth::smooth() {

  fP0   = NULL;
  fP1   = NULL;
  fP2   = NULL;
  fX    = NULL;
  fFunc = NULL;
  fN    = -1;
}
//-----------------------------------------------------------------------------
// before calling set smooth::gVar to the address of an actual structure
//-----------------------------------------------------------------------------
smooth::smooth(const TH1* Hist, double XMin, double XMax) {

  if (Hist == NULL) {
    fP0   = NULL;
    fP1   = NULL;
    fP2   = NULL;
    fX    = NULL;
    fFunc = NULL;
    printf("no hist to smooth\n");
  }
  else {

    int nx = Hist->GetNbinsX();
    fN = nx;
    double x1,x2,x3,y1,y2,y3;

    fP0   = new double[nx];
    fP1   = new double[nx];
    fP2   = new double[nx];
    fX    = new double[nx];

    if (XMin > XMax) {
      XMin = Hist->GetXaxis()->GetXmin();
      XMax = Hist->GetXaxis()->GetXmax();
    }

    //   printf("nx, XMin,XMax = %3i %12.5f %12.5f\n",nx,XMin,XMax);

    double  dx21, dx32, dx31, dy32, dy12, d;
    
    fX[0   ] = Hist->GetBinCenter(1);
    fX[nx-1] = Hist->GetBinCenter(nx);

    for (int i=2; i<=nx-1; i++) {
      y1 = Hist->GetBinContent(i-1);
      y2 = Hist->GetBinContent(i  );
      y3 = Hist->GetBinContent(i+1);

      dy32 = y3-y2;
      dy12 = y1-y2;

      x1 = Hist->GetBinCenter(i-1);
      x2 = Hist->GetBinCenter(i  );
      x3 = Hist->GetBinCenter(i+1);

      dx21 = x2-x1;
      dx31 = x3-x1;
      dx32 = x3-x2;
      d    = dx21*dx31*dx32;
      
      fP0[i] = y2;
      fP1[i] = (dy32*dx21*dx21-dy12*dx32*dx32)/d;
      fP2[i] = (dy32*dx21+dy12*dx32)/d;
      fX [i] = x2;
      
      //  printf(" i,Fx[i]: %2i %12.5e",i,fX[i]);
      // printf(" %12.5le %12.5le %12.5le ",y1,y2,y3);
      // printf(" %12.5le %12.5le %12.5le %12.5le ",dx21, dx31, dx32, d);
      // printf(" %12.5le %12.5le ",dy32, dy12);
      // printf(" %12.5le %12.5le %12.5le\n",fP0[i],fP1[i],fP2[i]);
    }
//-----------------------------------------------------------------------------
// parametrization is finished, how to use it? 
//-----------------------------------------------------------------------------
    fFunc = new TF1("hist_smooth_func",smooth::func,XMin,XMax,2);

    fFunc->SetParameter(0,1.); // scale factor
//------------------------------------------------------------------------------
// hack - passing the address to smooth::Eval()
//-----------------------------------------------------------------------------
//    printf(" fX = %p\n",fX);
    if (fX == 0) fFunc->SetParameter(1,0);
    else {
      double hidden_address = (double) ((long int) this);
      fFunc->FixParameter(1,hidden_address);
    }
					// fix the second parameter
    //    fFunc->SetParError(1,0.);
  }
}

//-----------------------------------------------------------------------------
// before calling set smooth::gVar to the address of an actual structure
//-----------------------------------------------------------------------------
smooth::smooth(const TGraph* Graph, double XMin, double XMax) {
  
  if (Graph == NULL) {
    fP0   = NULL;
    fP1   = NULL;
    fP2   = NULL;
    fX    = NULL;
    fFunc = NULL;
  }
  else {

    int nx = Graph->GetN();
    fN = nx;

    double* x  = Graph->GetX();
    double* y  = Graph->GetY();

    double y1, y2, y3, d;

    fX    = new double[nx];
    fP0   = new double[nx];
    fP1   = new double[nx];
    fP2   = new double[nx];

//   for (int i=1; i<=nx; i++) {
//     printf("%3i %10.5f %10.5f \n",i,Hist->GetBinCenter(i),Hist->GetBinContent(i));
//   }

    fX[0 ]   = x[0 ];
    fX[nx-1] = x[nx-1];
    
    for (int i=1; i<nx-1; i++) {
      fX [i] = x[i];

      y1 = y[i-1];
      y2 = y[i  ];
      y3 = y[i+1];

      double dx21 = x[i  ]-x[i-1];
      double dx31 = x[i+1]-x[i-1];
      double dx32 = x[i+1]-x[i  ];
      
      d    = dx21*dx31*dx32;

      double dy32 = y3-y2;
      double dy12 = y1-y2;

      fP0[i] = y2;
      fP1[i] = (dy32*dx21*dx21-dy12*dx32*dx32)/d;
      fP2[i] = (dy32*dx21+dy12*dx32)/d;

      // printf(" %2i %12.5e",i,fX[i]);
      // printf(" %12.5le %12.5le %12.5le ",y1,y2,y3);
      // printf(" %12.5le %12.5le %12.5le %12.5le ",dx21, dx31, dx32, d);
      // printf(" %12.5le %12.5le ",dy32, dy12);
      // printf(" %12.5le %12.5le %12.5le\n",fP0[i],fP1[i],fP2[i]);
    }

    if (XMin > XMax) {
      XMin = x[0];
      XMax = x[nx-1];
    }
//-----------------------------------------------------------------------------
// parametrization is finished, how to use it? 
//-----------------------------------------------------------------------------
    fFunc = new TF1("hist_smooth_func",smooth::func,XMin,XMax,2);

    fFunc->SetParameter(0,1.); // scale factor
//------------------------------------------------------------------------------
// hack - passing the address to smooth::Eval()
//-----------------------------------------------------------------------------
    if (fX == 0) {
      fFunc->SetParameter(1,0);
    }
    else {
      double hidden_address = (double) ((long int) this);
      fFunc->FixParameter(1,hidden_address);
    }
    // fix the second parameter
    //    fFunc->SetParError(1,0.);
  }
}

//-----------------------------------------------------------------------------
smooth::~smooth() {
  if (fX) {
    delete fP0;
    delete fP1;
    delete fP2;
    delete fX ;

    delete fFunc;
  }
}
