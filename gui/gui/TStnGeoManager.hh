#ifndef Stntuple_gui_TStnGeoManager_hh
#define Stntuple_gui_TStnGeoManager_hh

#include "TGeoVolume.h"
#include "TString.h"
#include "TGeoManager.h"

#include "Stntuple/gui/TEvdTracker.hh"

#include "Offline/DbTables/inc/TrkPanelMap.hh"
#include "Offline/TrackerConditions/inc/TrackerPanelMap.hh"

// namespace mu2e {
//   class TrkPanelMap;
//   class TrkPanelMap::Row;
//   class TrkPanelMapEntity;
//}

class  TEvdCrvSection;

class TStnGeoManager : public TNamed {
protected:

  static TStnGeoManager*  fgInstance;

public:
  enum { kNCrvSections = 22 };
  
  TGeoNode*       fTop; 

  TGeoNode*       fDs2Vacuum; 
  TGeoNode*       fDs3Vacuum; 
  TGeoNode*       fSttMother; 
  TGeoNode*       fCalMother; 
  TGeoNode*       fTrkMother;
  TGeoNode*       fMbsMother;

  TEvdCrvSection* fCrvSection[kNCrvSections];

  TGeoManager*    fGeoManager;

  stntuple::TEvdTracker*    fTracker;

  int             fTransp;

  TObjArray*      fListOfDetectors; 

  const mu2e::TrackerPanelMap* fTrackerPanelMap;
  
  float            fBField;          // by defautl , 1 T, but could be less.
                                     // Need to know to display straight cosmics

private:
  TStnGeoManager(const char* Name, const char* Fn, int OriginalColors);
  ~TStnGeoManager();
public:
  static TStnGeoManager* Instance(const char* Name = "StnGeoManager", const char* Fn = nullptr, int OriginalColors = 0);
//-----------------------------------------------------------------------------
// accessors
//-----------------------------------------------------------------------------
  float           BField() { return fBField; }
  TEvdCrvSection* CrvSection(int I) { return fCrvSection[I] ; }

					// here I explicitly assume that the name is known

  stntuple::TEvdTracker*    GetTracker() const {
    return (stntuple::TEvdTracker*) fListOfDetectors->FindObject("tracker");
  }
  
  TObject*        GetDetector(const char* Name) const { return fListOfDetectors->FindObject(Name); }
//-----------------------------------------------------------------------------
// setters
// from this point , TStnGeoManager owns detectors
//-----------------------------------------------------------------------------
  void AddDetector(TObject* Detector); //  { fListOfDetectors->Add(Detector); }

  void SetBField(float BField) { fBField = BField; }
  void SetRecursiveVisibility(TGeoVolume* Vol, int OnOff);
  void SetRecursiveVisibility(TGeoNode*   Vol, int OnOff);

  void SetDefaultColorTransp            (int Transp);
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

  void SetTrackerPanelMap(const mu2e::TrackerPanelMap* Map) { fTrackerPanelMap = Map; }

  void HideTsCoils (int OriginalColors);
  void HideDsCoils (int OriginalColors);
  void HideBuilding(int OriginalColors);
//-----------------------------------------------------------------------------
// Mu2e-specific - Node name starts with 'Pattern' assume it is unique
//-----------------------------------------------------------------------------
  TGeoNode* FindNodeByName      (TGeoNode*   Top, const char* Name      );
  TGeoNode* FindNodeByVolumeName(TGeoNode*   Top, const char* VolumeName);
  TGeoNode* FindNodeByVolumeName(TGeoVolume* Top, const char* VolumeName);
//-----------------------------------------------------------------------------
// other methods
//-----------------------------------------------------------------------------
  void DrawCRV();
  void DrawStrawTracker();
  void DrawCalorimeter();
  void DrawCalorimeterDisk();
  void DrawExtShielding();
  void DrawDetectorSolenoid();
  void DrawDetectorSolenoidDev2();
//-----------------------------------------------------------------------------
// geometry initialization
//-----------------------------------------------------------------------------
  int InitGeometry(); 
  int InitCalorimeterGeometry();
  int InitCrvGeometry();
  int InitTrackerGeometry();
//-----------------------------------------------------------------------------
// indexing
//-----------------------------------------------------------------------------
  const mu2e::TrkPanelMap::Row*  PanelMap(int PlaneID, int PanelID);
  
  ClassDef(TStnGeoManager,0)
};

#endif

