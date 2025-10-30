///////////////////////////////////////////////////////////////////////////////
// real base class for Stntuple VisNodes 
///////////////////////////////////////////////////////////////////////////////
#include "TPad.h"
#include "Stntuple/gui/TStnVisNode.hh"
#include "Stntuple/gui/TStnVisManager.hh"

//-----------------------------------------------------------------------------
TStnVisNode::TStnVisNode(const char* Name): TVisNode(Name) {
}


//-----------------------------------------------------------------------------
TStnVisNode::~TStnVisNode() {
}

//-----------------------------------------------------------------------------
void TStnVisNode::Paint   (Option_t* Option) {
				// parse option list

  int view_type = TStnVisManager::Instance()->GetCurrentView()->Type();

  if      (view_type == TStnVisManager::kCal ) PaintCal (Option);
  else if (view_type == TStnVisManager::kCrv ) PaintCrv (Option);
  else if (view_type == TStnVisManager::kPhiZ) PaintPhiZ(Option);
  else if (view_type == TStnVisManager::kRZ  ) PaintRZ  (Option);
  else if (view_type == TStnVisManager::kTZ  ) PaintTZ  (Option);
  else if (view_type == TStnVisManager::kVRZ ) PaintVRZ (Option);
  else if (view_type == TStnVisManager::kVST ) PaintVST (Option);
  else if (view_type == TStnVisManager::kVYZ ) PaintVRZ (Option);
  else if (view_type == TStnVisManager::kXY  ) PaintXY  (Option);
  else {
    Warning("Paint",Form("Unknown view_type %i",view_type));
  }
  
  gPad->Modified();
}

//-----------------------------------------------------------------------------
void TStnVisNode::PaintXY (Option_t* Option) {
  printf(Form(">>> ERROR in : TStnVisNode::%s: TXXXVisNode::%s is not implemented\n",__func__,__func__));
}

//-----------------------------------------------------------------------------
void TStnVisNode::PaintRZ (Option_t* Option) {
  printf(Form(">>> ERROR in : TStnVisNode::%s: TXXXVisNode::%s is not implemented\n",__func__,__func__));
}

//-----------------------------------------------------------------------------
void TStnVisNode::PaintTZ (Option_t* Option) {
  printf(Form(">>> ERROR in : TStnVisNode::%s: TXXXVisNode::%s is not implemented\n",__func__,__func__));
}

//-----------------------------------------------------------------------------
void TStnVisNode::PaintPhiZ(Option_t* Option) {
  printf(Form(">>> ERROR in : TStnVisNode::%s: TXXXVisNode::%s is not implemented\n",__func__,__func__));
}

//-----------------------------------------------------------------------------
void TStnVisNode::PaintCal(Option_t* Option) {
  printf(Form(">>> ERROR in : TStnVisNode::%s: TXXXVisNode::%s is not implemented\n",__func__,__func__));
}

//-----------------------------------------------------------------------------
void TStnVisNode::PaintCrv(Option_t* Option) {
  printf(Form(">>> ERROR in : TStnVisNode::%s: TXXXVisNode::%s is not implemented\n",__func__,__func__));
}

//-----------------------------------------------------------------------------
void TStnVisNode::PaintVST(Option_t* Option) {
  printf(Form(">>> ERROR in : TStnVisNode::%s: TXXXVisNode::%s is not implemented\n",__func__,__func__));
}

//-----------------------------------------------------------------------------
void TStnVisNode::PaintVRZ(Option_t* Option) {
  printf(Form(">>> ERROR in : TStnVisNode::%s: TXXXVisNode::%s is not implemented\n",__func__,__func__));
}

//-----------------------------------------------------------------------------
void TStnVisNode::PaintVYZ(Option_t* Option) {
  printf(Form(">>> ERROR in : TStnVisNode::%s: TXXXVisNode::%s is not implemented\n",__func__,__func__));
}

//-----------------------------------------------------------------------------
int  TStnVisNode::DistancetoPrimitive(Int_t px, Int_t py) {
  // by default, return a large number
  // decide how to deal with 3D views later

  int dist(9999);

  int view  = TStnVisManager::Instance()->GetCurrentView()->Type();
  //  int index = TStnVisManager::Instance()->GetCurrentView()->Index();
 
  if      (view == TStnVisManager::kCal ) dist = DistancetoPrimitiveCal (px,py);
  else if (view == TStnVisManager::kCrv ) dist = DistancetoPrimitiveCrv (px,py);
  else if (view == TStnVisManager::kPhiZ) dist = DistancetoPrimitivePhiZ(px,py);
  else if (view == TStnVisManager::kRZ  ) dist = DistancetoPrimitiveRZ  (px,py);
  else if (view == TStnVisManager::kTZ  ) dist = DistancetoPrimitiveTZ  (px,py);
  else if (view == TStnVisManager::kVST ) dist = DistancetoPrimitiveVST (px,py);
  else if (view == TStnVisManager::kVRZ ) dist = DistancetoPrimitiveVRZ (px,py);
  else if (view == TStnVisManager::kVYZ ) dist = DistancetoPrimitiveVYZ (px,py);
  else if (view == TStnVisManager::kXY  ) dist = DistancetoPrimitiveXY  (px,py);
  else {
    // what is the default?
    //    Warning("Paint",Form("Unknown option %s",option));
  }

  return dist;
}

//-----------------------------------------------------------------------------
int  TStnVisNode::DistancetoPrimitiveXY(Int_t px, Int_t py) {
  // by default, return a large number
  // decide how to deal with 3D views later

  return 10000000;
}

//-----------------------------------------------------------------------------
int  TStnVisNode::DistancetoPrimitiveRZ(Int_t px, Int_t py) {
  // by default, return a large number
  // decide how to deal with 3D views later

  return 10000000;
}

//-----------------------------------------------------------------------------
int  TStnVisNode::DistancetoPrimitiveTZ(Int_t px, Int_t py) {
  // by default, return a large number
  // decide how to deal with 3D views later

  return 10000;
}

//-----------------------------------------------------------------------------
int  TStnVisNode::DistancetoPrimitivePhiZ(Int_t px, Int_t py) {
  // by default, return a large number
  // decide how to deal with 3D views later

  return 10000;
}

//-----------------------------------------------------------------------------
int  TStnVisNode::DistancetoPrimitiveCal(Int_t px, Int_t py) {
  // by default, return a large number
  // decide how to deal with 3D views later

  return 10000000;
}

//-----------------------------------------------------------------------------
int  TStnVisNode::DistancetoPrimitiveCrv(Int_t px, Int_t py) {
  // by default, return a large number
  // decide how to deal with 3D views later

  return 10000000;
}

//-----------------------------------------------------------------------------
int  TStnVisNode::DistancetoPrimitiveVST(Int_t px, Int_t py) {
  // by default, return a large number
  // decide how to deal with 3D views later

  return 10000000;
}

//-----------------------------------------------------------------------------
int  TStnVisNode::DistancetoPrimitiveVRZ(Int_t px, Int_t py) {
  // by default, return a large number
  // decide how to deal with 3D views later

  return 10000000;
}

//-----------------------------------------------------------------------------
int  TStnVisNode::DistancetoPrimitiveVYZ(Int_t px, Int_t py) {
  // by default, return a large number
  // decide how to deal with 3D views later

  return 10000000;
}
