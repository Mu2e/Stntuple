///////////////////////////////////////////////////////////////////////////////
// draw different parts of Mu2e
//
// first example: 
//
// dmg = new DisplayMu2eGeometry("/projects/mu2e/geom/gdml/mu2e_geometry_v6_1_4.gdml")
// dmg->HideBuilding(1)
// dmg->gm->GetVolume("World")->Draw("ogl")
//
// for proper 2D ZX view, choose:
// - Camera/Orthographic ZnOX
// - Guides/Axes/Origin
// - Camera overlay/Show Mode; Grid Front or Axes
//
// comment: TGeoManager::Import chokes on filenames like "~/mu2e.gdml") 
///////////////////////////////////////////////////////////////////////////////
#include "Stntuple/gui/DisplayMu2eGeometry.hh"

//-----------------------------------------------------------------------------
DisplayMu2eGeometry::DisplayMu2eGeometry(const char* Fn, int KeepOriginalColors) {
  gm = new TGeoManager();
  gm->Import(Fn);

  fTop       = gGeoManager->GetTopNode();

  fDS1Vacuum = FindNodeByVolumeName(fTop,"DS1Vacuum");
  fDS2Vacuum = FindNodeByVolumeName(fTop,"DS2Vacuum");
  fDS3Vacuum = FindNodeByVolumeName(fTop,"DS3Vacuum");

  fTS1Vacuum = FindNodeByVolumeName(fTop,"TS1Vacuum");
  fTS2Vacuum = FindNodeByVolumeName(fTop,"TS2Vacuum");
  fTS3Vacuum = FindNodeByVolumeName(fTop,"TS3Vacuum");
  fTS4Vacuum = FindNodeByVolumeName(fTop,"TS4Vacuum");
  fTS5Vacuum = FindNodeByVolumeName(fTop,"TS5Vacuum");

  fSttMother = FindNodeByVolumeName(fDS2Vacuum,"StoppingTargetMother");

  fTrkMother = FindNodeByVolumeName(fDS3Vacuum,"TrackerMother");
  fCalMother = FindNodeByVolumeName(fDS3Vacuum,"CalorimeterMother");
  fMbsMother = FindNodeByVolumeName(fDS3Vacuum,"MBSMother");

  fTransp       = 40;
  fDefaultColor = kCyan-10;
  
  HideBuilding(KeepOriginalColors);
}

//-----------------------------------------------------------------------------
DisplayMu2eGeometry::~DisplayMu2eGeometry() {
  if (gm) delete gm;
}

//-----------------------------------------------------------------------------
TGeoNode* DisplayMu2eGeometry::FindNodeByName(TGeoNode* Top, const char* Name) {
  TGeoNode  *top, *found(0);

  if (Top) top = fTop;
  else     top = gm->GetTopNode();
  
  TObjArray* o =  top->GetNodes();

  int n = o->GetEntriesFast();
  
  for (int i=0; i<n; i++) {
    TGeoNode* node = (TGeoNode*) o->UncheckedAt(i);
    const char* name = node->GetName();
    if (strcmp(name,Name) == 0) {
      found = node;
      break;
    }
    else if (node->GetNodes() != NULL) {
      found = FindNodeByName(node,Name);
      if (found) break;
    }
  }
  return found;
}

//-----------------------------------------------------------------------------
TGeoNode* DisplayMu2eGeometry::FindNodeByVolumeName(TGeoNode* Top, const char* VolumeName) {
  TGeoNode  *top, *found(0);

  if (Top) top = Top;
  else     top = fTop;
  
  TObjArray* o =  top->GetNodes();

  int n = o->GetEntriesFast();
  
  for (int i=0; i<n; i++) {
    TGeoNode* node = (TGeoNode*) o->UncheckedAt(i);
    const char* name = node->GetVolume()->GetName();
    if (strcmp(name,VolumeName) == 0) {
      found = node;
      break;
    }
    else if (node->GetNodes() != NULL) {
      found = FindNodeByVolumeName(node,VolumeName);
      if (found) break;
    }
  }
  return found;
}

//-----------------------------------------------------------------------------
// assume that we're looking for one of the daughters
//-----------------------------------------------------------------------------
TGeoNode* DisplayMu2eGeometry::FindNodeByVolumeName(TGeoVolume* Top, const char* VolumeName) {
  TGeoVolume  *top;
  TGeoNode*    found(0);

  if (Top) top = Top;
  else     top = gm->GetTopVolume();

  TObjArray* o =  top->GetNodes();

  int n = o->GetEntriesFast();
  
  for (int i=0; i<n; i++) {
    TGeoNode* node = (TGeoNode*) o->UncheckedAt(i);
    const char* name = node->GetVolume()->GetName();
    if (strcmp(name,VolumeName) == 0) {
      found = node;
      break;
    }
    else if (node->GetNodes() != NULL) {
      found = FindNodeByVolumeName(node,VolumeName);
      if (found) break;
    }
  }
  return found;
}

//-----------------------------------------------------------------------------
void DisplayMu2eGeometry::SetRecursiveVisibilityByName(TGeoNode* Node, const char* Pattern, int OnOff) {

  TString name(Node->GetName());
  
  if (name.Index(Pattern) >= 0) {
    Node->SetVisibility(OnOff);
    //std::cout <<"hiding "<< name << std::endl;
  }
				        // Descend recursively into each daughter TGeoNode.
  int nd = Node->GetNdaughters();
  for (int i=0; i<nd; ++i) {
    TGeoNode * d = Node->GetDaughter(i);
    SetRecursiveVisibilityByName(d,Pattern,OnOff);
  }
}

//-----------------------------------------------------------------------------
void DisplayMu2eGeometry::SetRecursiveVisibility(TGeoNode* Node, int OnOff) {

  TString name(Node->GetName());
  
  Node->SetVisibility(OnOff);
    //std::cout <<"hiding "<< name << std::endl;
				        // Descend recursively into each daughter TGeoNode.
  int nd = Node->GetNdaughters();
  for (int i=0; i<nd; ++i) {
    TGeoNode * d = Node->GetDaughter(i);
    SetRecursiveVisibility(d,OnOff);
  }
}

//-----------------------------------------------------------------------------
void DisplayMu2eGeometry::SetRecursiveVisibilityByMaterial(TGeoNode* Node, const char* Material, int OnOff) {

  TString mat(Node->GetVolume()->GetMaterial()->GetName());
  
  if (mat.Index(Material) >= 0) Node->SetVisibility(OnOff);

				        // Descend recursively into each daughter TGeoNode.
  int ndau = Node->GetNdaughters();
  for ( int i=0; i<ndau; ++i ){
    TGeoNode * d = Node->GetDaughter(i);
    SetRecursiveVisibilityByMaterial(d,Material,OnOff);
  }
}

//-----------------------------------------------------------------------------
void DisplayMu2eGeometry::SetRecursiveColorTranspByName(TGeoNode* Node, const char* Name, Int_t Color, Int_t Transp) {

  
  TString node_name = Node->GetName();
  TGeoVolume*  vol = Node->GetVolume();

  if (node_name.Index(Name) >= 0) {
    vol->SetLineColor   (Color );
    vol->SetTransparency(Transp);
  }
     
  int nd = vol->GetNdaughters();
  for (int i=0; i<nd; i++) {
    TGeoNode* dn = vol->GetNode(i);
    SetRecursiveColorTranspByName(dn, Name, Color, Transp);
  }
}

//-----------------------------------------------------------------------------
void DisplayMu2eGeometry::SetRecursiveColorTranspByMaterial(TGeoNode* Node, const char* MatName, Int_t Color, Int_t Transp) {

  
  TGeoVolume*  vol = Node->GetVolume();
  TString mat_name = vol->GetMaterial()->GetName();

  if (mat_name.Index(MatName) >= 0) {
    vol->SetLineColor   (Color );
    vol->SetTransparency(Transp);
  }
     
  int nd = vol->GetNdaughters();
  for (int i=0; i<nd; i++) {
    TGeoNode* dn = vol->GetNode(i);
    SetRecursiveColorTranspByMaterial(dn, MatName, Color, Transp);
  }
}

//-----------------------------------------------------------------------------
void DisplayMu2eGeometry::SetRecursiveVisibilityColorTranspByNameAndMaterial(TGeoNode*   Top    ,
									  const char* Name   ,
									  const char* MatName,
									  int         Visibility,
									  int         Color  ,
									  int         Transp) {
  TGeoVolume*  vol = Top->GetVolume();
  TString name     = vol->GetName();
  TString mat_name = vol->GetMaterial()->GetName();

  if ((name.Index(Name) >= 0) && (mat_name.Index(MatName) >= 0)) {
    Top->SetVisibility  (Visibility);
    vol->SetLineColor   (Color );
    vol->SetTransparency(Transp);
  }
     
  int nd = vol->GetNdaughters();
  for (int i=0; i<nd; i++) {
    TGeoNode* node = vol->GetNode(i);
    SetRecursiveVisibilityColorTranspByNameAndMaterial(node,Name,MatName,Visibility,Color,Transp);
  }
}

//-----------------------------------------------------------------------------
void DisplayMu2eGeometry::SetRecursiveColorTransp(TGeoVolume *Vol, Int_t Color, Int_t Transp) {

  TString name = Vol->GetName();
  
  int col    = Color;
  int transp = Transp;

  if      (name.Index("TargetFoil") >= 0) { col = kBlue+4;  transp = 10; }
  else if (name.Index("CaloPipe"  ) >= 0) { col = kOrange+7; }
    
  if (col    >=0 ) Vol->SetLineColor   (col   );
  if (Transp >=0 ) Vol->SetTransparency(transp);
     
  int nd = Vol->GetNdaughters();
  for (int i=0; i<nd; i++) {
    TGeoVolume* vd = Vol->GetNode(i)->GetVolume();
    SetRecursiveColorTransp(vd, Color, transp);
  }
}

//-----------------------------------------------------------------------------
// everything is kCyan by default
// production tsrget is kRed+3
//-----------------------------------------------------------------------------
void DisplayMu2eGeometry::SetDefaultColorTransp() {
  SetRecursiveColorTransp(fTop->GetVolume(),fDefaultColor,fTransp);

//-----------------------------------------------------------------------------
// color production target
//-----------------------------------------------------------------------------
  TGeoVolume* ptarget = gm->GetVolume("ProductionTargetMother");
  int col             = kRed+2;
  int nd              = ptarget->GetNdaughters();
  const char*  name;
  
  for (int i=0; i<nd; i++) {
    TGeoVolume* vd = ptarget->GetNode(i)->GetVolume();
    name = vd->GetName();
    printf(" production target daughter: %s\n",name);
    if      (strcmp(name,"ProductionTargetSupportWheel") == 0) col = kGray;
    else if (strcmp(name,"ClampSupportWheel_R")          == 0) col = kGray+3;
    else                                                       col = kRed+2;
    SetRecursiveColorTransp(vd,col,fTransp);
  }
//-----------------------------------------------------------------------------
// color proton absorbers, start from TS1
//-----------------------------------------------------------------------------
  TGeoVolume* ts1_vacuum = gm->GetVolume("TS1Vacuum");
  col             = kRed+2;
  nd              = ts1_vacuum->GetNdaughters();
  
  int ts1_coll_transp = 30;
  for (int i=0; i<nd; i++) {
    TGeoVolume* vd = ts1_vacuum->GetNode(i)->GetVolume();
    name = vd->GetName();
    printf(" TS1Vacuum daughter: %s\n",name);
    if      (strcmp(name,"PbarAbsTS1Out") == 0) SetRecursiveColorTransp(vd,kRed+1       ,fTransp        );
    if      (strcmp(name,"Coll11"       ) == 0) SetRecursiveColorTransp(vd,kRed+3       ,ts1_coll_transp);
    else if (strcmp(name,"Coll12"       ) == 0) SetRecursiveColorTransp(vd,fDefaultColor,ts1_coll_transp);
    else if (strcmp(name,"Coll13"       ) == 0) SetRecursiveColorTransp(vd,kRed+3       ,ts1_coll_transp);
  }
//-----------------------------------------------------------------------------
// proceed with TS3
//-----------------------------------------------------------------------------
  TGeoVolume* ts3_vacuum = gm->GetVolume("TS3Vacuum");
  col             = kRed+2;
  nd              = ts3_vacuum->GetNdaughters();
  
  int ts3_coll_transp = 30;
  for (int i=0; i<nd; i++) {
    TGeoVolume* vd = ts3_vacuum->GetNode(i)->GetVolume();
    name = vd->GetName();
    printf(" TS3Vacuum daughter: %s\n",name);
    if      (strcmp(name,"PbarAbs"     ) == 0) SetRecursiveColorTransp(vd,kRed+2,fTransp);
    else if (strcmp(name,"PbarAbsWedge") == 0) SetRecursiveColorTransp(vd,kRed+2,fTransp);
    if      (strcmp(name,"Coll31"      ) == 0) SetRecursiveColorTransp(vd,kRed+2,ts3_coll_transp);
    else if (strcmp(name,"Coll32"      ) == 0) SetRecursiveColorTransp(vd,kRed+2,ts3_coll_transp);
  }
//-----------------------------------------------------------------------------
// color TS5 collimator
//-----------------------------------------------------------------------------
  TGeoVolume* ts5_vacuum = gm->GetVolume("TS5Vacuum");
  nd              = ts5_vacuum->GetNdaughters();
  
  int ts5_coll_transp = 30;
  for (int i=0; i<nd; i++) {
    TGeoVolume* vd = ts5_vacuum->GetNode(i)->GetVolume();
    name = vd->GetName();
    printf(" TS5Vacuum daughter: %s\n",name);
    if      (strcmp(name,"Coll51") == 0) SetRecursiveColorTransp(vd,kRed+2,ts5_coll_transp);
    else if (strcmp(name,"Coll52") == 0) SetRecursiveColorTransp(vd,kRed+2,ts5_coll_transp);
  }
//-----------------------------------------------------------------------------
// color things in the DS2
//-----------------------------------------------------------------------------
  // TGeoVolume* ds2_vacuum = gm->GetVolume("DS2Vacuum");
  // nd              = ds2_vacuum->GetNdaughters();
  
  // int ds2_coll_transp = 30;
  // for (int i=0; i<nd; i++) {
  //   TGeoVolume* vd = ds2_vacuum->GetNode(i)->GetVolume();
  //   name = vd->GetName();
  //   printf(" DS2Vacuum daughter: %s\n",name);
  //   if      (strcmp(name,"protonabs1") == 0) SetRecursiveColorTransp(vd,kRed+2 ,70);
  //   else if (strcmp(name,"protonabs3") == 0) SetRecursiveColorTransp(vd,kBlue+1,70);
  // }
//-----------------------------------------------------------------------------
// color things in the tracker and the calorimeter
//-----------------------------------------------------------------------------
  SetAbsorberColors();
  SetTrackerColors();
  SetCalorimeterColors();
}

//-----------------------------------------------------------------------------
void DisplayMu2eGeometry::HideBuilding(int KeepOriginalColors) {

  // Volumes will be made invisible if their name contains one
  // of these strings.
  //"CRS", "ExtShield""CRV"
  
  static TString name[] = {
    "Ceiling", "backfill", "dirt", "concrete",
    "VirtualDetector",
    "pipeType",
    "CRSAluminium",        // CRV
    "CRV","CRS","crv",     // CRV
    "ElectronicRackBox",   // electronics aside
    "ExtShield",
    "ExtMon",              // ExtMon
    "collimator1Channel",  // ExtMon
    "collimator2Channel",  // ExtMon
    "EMFPlane"          ,  // ExtMon
    "coll2Shielding",
    "pBendType22",         // who knows what it is ?
    "pBendType2",          // who knows what it is ?
    "dsAreaTrenchCover", 
    "ProtonBeam",
    "collimatorFOV",
    "FOVliner",
				// pieces inside DS
    "DSCoil",
    "DSSpacer",
    "DScenterRing",
    "DSleftSideRing",
    "DSrightSideRing",
    "BearingBlock",

    "stmMagnet",           // STM magnet and its support
    "stmDet",              // STM far behind
    "collimatorSS",	     // STM
    
    "PSEnclosureShell",    // part of PS
    "PSEnclosureWindow",   // part of PS
    
    "BearingBlock_DS2",
    "psAreaHatchLid",
    "remoteHandlingHatchLid",
    
    ""
  };

  for (int i=0; name[i] != ""; i++) {
    SetRecursiveVisibilityByName(fTop,name[i].Data(),0);
  }
//-----------------------------------------------------------------------------
// Volumes with these material names will be made invisible : 
//-----------------------------------------------------------------------------
  static TString material[] = {
    "MBOverburden",
    "CONCRETE",
    "BARITE",
    ""
  };

  for (int i=0; material[i] != ""; i++) {
    SetRecursiveVisibilityByMaterial(fTop,material[i].Data(),0);
  }
//-----------------------------------------------------------------------------
// inside DS3Vacuum: hide calorimeter electronics, MBS
//-----------------------------------------------------------------------------
  SetRecursiveVisibilityByName(fDS3Vacuum,"VPSP_"     ,0);
  SetRecursiveVisibilityByName(fDS3Vacuum,"IFB_"      ,0);
  SetRecursiveVisibilityByName(fDS3Vacuum,"protonabs4",0);

  SetRecursiveVisibilityByName(fDS3Vacuum,"CalorimeterFEB"    ,0);
  SetRecursiveVisibilityByName(fDS3Vacuum,"crateBoxLog"       ,0);
  SetRecursiveVisibilityByName(fDS3Vacuum,"crateSideLog"      ,0);
  SetRecursiveVisibilityByName(fDS3Vacuum,"crateTopLog"       ,0);
  SetRecursiveVisibilityByName(fDS3Vacuum,"crateBottomALog"   ,0);
  SetRecursiveVisibilityByName(fDS3Vacuum,"crateBottomBLog"   ,0);

      SetRecursiveVisibilityByName(fDS3Vacuum,"boardCrateLog"       ,0);
  SetRecursiveVisibilityByName(fDS3Vacuum,"radiatorBoardLog"    ,0);
  SetRecursiveVisibilityByName(fDS3Vacuum,"activeStripBoardLog" ,0);
  SetRecursiveVisibilityByName(fDS3Vacuum,"passiveStripBoardLog",0);
//-----------------------------------------------------------------------------
// hide things inside the tracker
//-----------------------------------------------------------------------------
					// tracker support
  SetRecursiveVisibilityByName(fTrkMother,"NorthRailDS"       ,0);
  SetRecursiveVisibilityByName(fTrkMother,"SouthRailDS"       ,0);
  SetRecursiveVisibilityByName(fTrkMother,"TrackerSupport"    ,0);
  SetRecursiveVisibilityByName(fTrkMother,"TTrackerEndRing"   ,0);
  SetRecursiveVisibilityByName(fTrkMother,"ThinSupportRing"   ,0);
//-----------------------------------------------------------------------------
// hide things inside the calorimeter
//-----------------------------------------------------------------------------
  SetRecursiveVisibilityByName(fCalMother,"DiskInnerRing"     ,0);
  SetRecursiveVisibilityByName(fCalMother,"DiskOuterRing"     ,0);
  //  SetRecursiveVisibilityByName(fCalMother,"DiskCase"          ,0);
  SetRecursiveVisibilityByName(fCalMother,"CrystalROLog"      ,0);
  SetRecursiveVisibilityByName(fCalMother,"ElectronicsROLog"  ,0);
  SetRecursiveVisibilityByName(fCalMother,"WrapLog"           ,0);
  // SetRecursiveVisibilityByName(fCalMother,"UnitLog"           ,0);
  
//-----------------------------------------------------------------------------
// dont show MBS
//-----------------------------------------------------------------------------
  SetRecursiveVisibility(fMbsMother,0);
//-----------------------------------------------------------------------------
// hide last saddle boxes
//-----------------------------------------------------------------------------
  char saddle_box_name[50];
  for (int i=79; i<120; i++) {
    sprintf(saddle_box_name,"SaddleBox_%i",i);
    SetRecursiveVisibilityByName(fTop,saddle_box_name,0);
  }
//-----------------------------------------------------------------------------
// colors
//-----------------------------------------------------------------------------
  if (KeepOriginalColors == 0) SetDefaultColorTransp();
}


//-----------------------------------------------------------------------------
void DisplayMu2eGeometry::SetAbsorberColors() {
  TGeoNode* abs_node = FindNodeByVolumeName(fDS2Vacuum,"protonabs1");
				      
  abs_node->GetVolume()->SetLineColor(804);
  //  gm->GetVolume("protonabs3")->SetLineColor(808);

  // gm->GetVolume("InternalNeutronAbsorber1")->SetLineColor(900); // default = 920
  // gm->GetVolume("InternalNeutronAbsorber2")->SetLineColor(850); // default = 920
}


//-----------------------------------------------------------------------------
// set crystal color to kOrange-2
//-----------------------------------------------------------------------------
void DisplayMu2eGeometry::SetCalorimeterColors() {
  SetRecursiveColorTranspByMaterial(fCalMother,"G4_CESIUM_IODIDE",kOrange-2  ,40);
}

//-----------------------------------------------------------------------------
void DisplayMu2eGeometry::DrawCalorimeterDisk() {
  gm->GetVolume("DiskCalorimeter_0")->Draw("ogl");
}

//-----------------------------------------------------------------------------
void DisplayMu2eGeometry::DrawCalorimeter() {

  HideBuilding(0);

  SetCalorimeterColors();
  gm->GetVolume("CalorimeterMother")->Draw("ogl");
}

//-----------------------------------------------------------------------------
void DisplayMu2eGeometry::DrawExtShielding() {

  HideBuilding(0);
  
  static TString name[] = {
    "ExtShield",
    ""
  };

  for (int i=0; name[i] != ""; i++) {
    SetRecursiveVisibilityColorTranspByNameAndMaterial(fTop,name[i].Data(),"BARITE"       ,1,kBlue+2   ,0);
    SetRecursiveVisibilityColorTranspByNameAndMaterial(fTop,name[i].Data(),"CONCRETE_CB4" ,1,kGray+2   ,0);
    SetRecursiveVisibilityColorTranspByNameAndMaterial(fTop,name[i].Data(),"CONCRETE_MARS",1,kMagenta+2,0);
  }

  // SetRecursiveColorTranspByMaterial(fTop,"BARITE"       ,kBlue+2   ,0);
  // SetRecursiveColorTranspByMaterial(fTop,"CONCRETE_CB4" ,kGray+2   ,0);
  // SetRecursiveColorTranspByMaterial(fTop,"CONCRETE_MARS",kMagenta-2,0);

  static TString crv_name[] = {
    "CRSAluminium","CRV","CRS","crv",
    ""
  };

  for (int i=0; crv_name[i] != ""; i++) {
    SetRecursiveVisibilityByName(fTop,crv_name[i].Data(),1);
  }
  SetRecursiveColorTranspByMaterial(fTop,"G4_POLYSTYRENE",kYellow-9  ,30);

  gm->GetVolume("HallAir")->Draw("ogl");
}

//-----------------------------------------------------------------------------
void DisplayMu2eGeometry::DrawCRV() {

  HideBuilding(0);
  
  static TString name[] = {
    //    "ExtShield",
    "CRSAluminium",        // CRV
    "CRV","CRS","crv",     // CRV
    ""
  };

  for (int i=0; name[i] != ""; i++) {
    SetRecursiveVisibilityByName(fTop,name[i].Data(),1);
  }

  SetRecursiveColorTranspByMaterial(fTop,"G4_POLYSTYRENE",kYellow-9  ,0);
  SetRecursiveColorTranspByMaterial(fTop,"G4_Al"         ,kGray      ,0);
  SetRecursiveColorTranspByMaterial(fTop,"ElectronicsFEB",kGray+2    ,0);

  SetRecursiveColorTranspByMaterial(fTop,"BARITE"       ,kBlue+2   ,0);
  SetRecursiveColorTranspByMaterial(fTop,"CONCRETE_CB4" ,kGray+2   ,0);
  SetRecursiveColorTranspByMaterial(fTop,"CONCRETE_MARS",kMagenta-2,0);

  gm->GetVolume("HallAir")->Draw("ogl");
}



//-----------------------------------------------------------------------------
void DisplayMu2eGeometry::SetTrackerColors() {
  static TString name[] = {
    "TTrackerSupport",
    "TTrackerEndRingUpstream",
    ""
  };
  
  SetRecursiveColorTranspByName(fTrkMother,"TTracker",kYellow   ,90);
  SetRecursiveColorTranspByName(fTrkMother,"Plane"   ,kYellow   ,99);
  SetRecursiveColorTranspByName(fTrkMother,"Panel"   ,kYellow   ,99);
  
  SetRecursiveColorTranspByName(fTrkMother,"TTrackerEndRingUpstream",kGray  ,0);
  SetRecursiveColorTranspByName(fTrkMother,"TTrackerSupport"        ,kGray  ,0);
  SetRecursiveColorTranspByName(fTrkMother,"TTrackerSupportBeam"    ,kGray+2,0);
  SetRecursiveColorTranspByName(fTrkMother,"TTrackerStrawGas"       ,kYellow,0);
}

//-----------------------------------------------------------------------------
void DisplayMu2eGeometry::DrawStrawTracker() {

  HideBuilding(0);

  SetTrackerColors();
  gm->GetVolume("TrackerMother")->Draw("ogl");
}


//-----------------------------------------------------------------------------
void DisplayMu2eGeometry::DrawProductionTarget() {
  gm->GetVolume("ProductionTargetMother")->Draw("ogl");
}

//-----------------------------------------------------------------------------
// the names are this is for v4_0_6
// clipping: xc = -387.93
//           zlen = 1429
//-----------------------------------------------------------------------------
void DisplayMu2eGeometry::DrawDetectorSolenoid() {

  TGeoVolume* hall = fTop->GetVolume();

  // TObjArray* list_of_nodes = hall->GetNodes();

  // int n_nodes = list_of_nodes->GetEntries();

  // TGeoVolume* node;

  // for (int i=0; i<n_nodes; i++) {
  //   node  = (TGeoVolume*) list_of_nodes->At(i);

  //   const char* name = node->GetName();

  //   //    printf(" -- node name: %s\n",name);

  //   if ((strstr(name,"DS2Vacuum") != 0) ||
  // 	(strstr(name,"DS3Vacuum") != 0)    ) {
  //     node->SetVisibility(1);
  //     node->SetVisDaughters(1);
  //     node->SetVisLeaves(1);
  //   }
  //   else {
  //     node->SetVisibility(0);
  //     node->SetVisDaughters(0);
  //     node->SetVisLeaves(0);
  //   }
  // }

  gm->GetVolume("DS1Vacuum")->SetVisibility(0);
  gm->GetVolume("DS2Vacuum")->SetVisibility(0);
  gm->GetVolume("DS3Vacuum")->SetVisibility(0);

  gm->GetVolume("protonabs1")->SetLineColor(807);


  SetRecursiveVisibilityByName(fTop,"MBSMother",0);
  
  // TGeoVolume* mbs = gm->GetVolume("MBSMother");
  // mbs->SetVisibility(0);
  // mbs->SetVisDaughters(0);
  // mbs->SetVisLeaves(0);

  gm->GetVolume("protonabs1")->SetLineColor(804);
  gm->GetVolume("protonabs3")->SetLineColor(808);

  // gm->GetVolume("InternalNeutronAbsorber1" )->SetLineColor(900); // default = 920
  // gm->GetVolume("InternalNeutronAbsorber2" )->SetLineColor(850); // default = 920
  // gm->GetVolume("InternalNeutronAbsorber3a")->SetLineColor(860); // default = 920

  //  gm->GetVolume("Foil")->SetLineColor(20);

  hall->Draw("ogl");
}

//-----------------------------------------------------------------------------
// the names are this is for v4_0_6
// clipping: xc = -387.93
//           zlen = 1429
//-----------------------------------------------------------------------------
void DisplayMu2eGeometry::DrawDetectorSolenoidDev2() {

  TGeoVolume* hall = fTop->GetVolume();

  TObjArray* list_of_nodes = hall->GetNodes();

  int n_nodes = list_of_nodes->GetEntries();

  TGeoNode    *node, *node2;
  const char  *name, *name2;

  for (int i=0; i<n_nodes; i++) {
    node  = hall->GetNode(i);

    name = node->GetName();

    //    printf(" -- node name: %s\n",name);

    if ((strstr(name,"DS2Vacuum") != 0) ||
	(strstr(name,"DS3Vacuum") != 0)    ) {

      node->SetVisibility(0);
      node->SetVisDaughters(1);
      node->SetVisLeaves(1);


      TObjArray* list_of_nodes_2 = node->GetNodes();
      int n_nodes_2 = list_of_nodes_2->GetEntries();

      for (int i2=0; i2<n_nodes_2; i2++) {
	node2  = (TGeoNode*) node->GetVolume()->GetNode(i2);
	name2 = node2->GetName();
	//	printf("            -- node2 name: %s\n",name2);
	if (strstr(name2,"VirtualDetector") != 0) {
	  node2->SetVisibility(0);
	  node2->SetVisDaughters(0);
	  node2->SetVisLeaves(0);
	}
	else if (strstr(name2,"MBSMother") != 0) {
	  node2->SetVisibility(0);
	  node2->SetVisDaughters(0);
	  node2->SetVisLeaves(0);
	}
	else if (strstr(name2,"VPSP") != 0) {
	  node2->SetVisibility(0);
	  node2->SetVisDaughters(0);
	  node2->SetVisLeaves(0);
	}
	else if (strstr(name2,"IFB") != 0) {
	  node2->SetVisibility(0);
	  node2->SetVisDaughters(0);
	  node2->SetVisLeaves(0);
	}
      }
    }
    else {
      node->SetVisibility(0);
      node->SetVisDaughters(0);
      node->SetVisLeaves(0);
    }
  }

  SetAbsorberColors();
  
  //  gm->GetVolume("InternalNeutronAbsorber3a")->SetLineColor(860); // default = 920

  //  gm->GetVolume("Foil")->SetLineColor(20);

  hall->Draw("ogl");
}

//-----------------------------------------------------------------------------
void DisplayMu2eGeometry::Help() {
  printf("dmg = new DisplayMu2eGeometry(\"/home/murat/figures/mu2e/gdml/mu2e_geometry_v6_1_4.gdml\")\n");
  printf("dmg->HideBuilding(1)\n");
  printf("dmg->gm->GetVolume(\"HallAir\")->Draw(\"ogl\")\n");
  printf("\n");
  printf("comment: TGeoManager::Import chokes on filenames like \"~/mu2e.gdml\"\n");
}
