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
#ifndef __STNTUPLE_DisplayMu2eGeometry__
#define __STNTUPLE_DisplayMu2eGeometry__

#include "TGeoVolume.h"
#include "TGeoManager.h"
#include "TString.h"

class DisplayMu2eGeometry {
public:
  TGeoManager* gm;

  TGeoNode*    fTop; 

  TGeoNode*    fDS1Vacuum; 
  TGeoNode*    fDS2Vacuum; 
  TGeoNode*    fDS3Vacuum; 

  TGeoNode*    fTS1Vacuum; 
  TGeoNode*    fTS2Vacuum; 
  TGeoNode*    fTS3Vacuum; 
  TGeoNode*    fTS4Vacuum; 
  TGeoNode*    fTS5Vacuum; 

  TGeoNode*    fSttMother; 
  TGeoNode*    fCalMother; 
  TGeoNode*    fTrkMother;
  TGeoNode*    fMbsMother;
  
  int          fTransp;
  int          fDefaultColor;

  DisplayMu2eGeometry(const char* Fn = "/home/murat/figures/mu2e/gdml/mu2e_geometry_v4_0_6.gdml", int OriginalColors = 0);
  ~DisplayMu2eGeometry();
  
  void SetRecursiveVisibility(TGeoVolume* Vol, int OnOff);
  void SetRecursiveVisibility(TGeoNode*   Vol, int OnOff);

  void SetDefaultColorTransp            ();
  void SetRecursiveColorTransp          (TGeoVolume* Vol , Int_t Color, Int_t Transp);
  void SetRecursiveColorTranspByName    (TGeoNode*   Vol , const char* Name   , Int_t Color, Int_t Transp);
  void SetRecursiveColorTranspByMaterial(TGeoNode*   Node, const char* MatName, Int_t Color, Int_t Transp);

  void SetRecursiveVisibilityColorTranspByNameAndMaterial(TGeoNode*   Top         ,
							  const char* Name        ,
							  const char* MatName     ,
							  int         Visibility  ,
							  int         Color       ,
							  int         Transparency);
    
  void SetRecursiveVisibilityByName    (TGeoNode* Node, const char* NamePattern, int OnOff);
  void SetRecursiveVisibilityByMaterial(TGeoNode* Node, const char* Material   , int OnOff);

  void HideBuilding(int OriginalColors);
  
				// Mu2e-specific - Node name starts with 'Pattern'
				// assume it is unique
  
  TGeoNode* FindNodeByName      (TGeoNode*   Top, const char* Name      );
  TGeoNode* FindNodeByVolumeName(TGeoNode*   Top, const char* VolumeName);
  TGeoNode* FindNodeByVolumeName(TGeoVolume* Top, const char* VolumeName);
    
  void DrawCRV();
  void DrawCalorimeter();
  void DrawCalorimeterDisk();
  void DrawExtShielding();
  void DrawDetectorSolenoid();
  void DrawDetectorSolenoidDev2();
  void DrawProductionTarget();
  void DrawStrawTracker();

  void SetAbsorberColors();
  void SetTrackerColors();
  void SetCalorimeterColors();
  
  static void Help();
};

#endif
