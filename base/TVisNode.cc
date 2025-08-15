///////////////////////////////////////////////////////////////////////////////
#include "Stntuple/base/TVisNode.hh"

ClassImp(TVisNode)


int TVisNode::fgDebugLevel(0);

//_____________________________________________________________________________
TVisNode::TVisNode(const char* name):
  fName(name)
{
  fClosestObject = NULL;
  fgDebugLevel    = 0;
}


//_____________________________________________________________________________
TVisNode::~TVisNode() {
}

//_____________________________________________________________________________
void TVisNode::NodePrint(const void* Object, const char* ClassName) {
}

