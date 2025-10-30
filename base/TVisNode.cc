///////////////////////////////////////////////////////////////////////////////
#include "Stntuple/base/TVisNode.hh"

// ClassImp(TVisNode)

//_____________________________________________________________________________
TVisNode::TVisNode(const char* name): fName(name) {
  fClosestObject = NULL;
  fDebugLevel    = 0;
}

//_____________________________________________________________________________
TVisNode::~TVisNode() {
}

//_____________________________________________________________________________
int TVisNode::InitEvent() {
  printf(">>> ERROR in TVisNode::%s: derived class TXXXNode::%s is not implemented\n",__func__,__func__);
  return 0;
}

//_____________________________________________________________________________
void TVisNode::NodePrint(const void* Object, const char* ClassName) {
  printf(">>> ERROR in TVisNode::%s: derived class TXXXNode::%s is not implemented\n",__func__,__func__);
}
