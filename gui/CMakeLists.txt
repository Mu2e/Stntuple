#
set (subpkg "gui"                    )
set (dict   "Stntuple_${subpkg}_dict")
set (lib    "Stntuple_${subpkg}"     )

include_directories(
  $ENV{ROOTSYS}/include/root
  ${CMAKE_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/${subpkg}
)

# link_directories($ENV{ROOTSYS}/lib)

set (include_files
  DisplayMu2eGeometry.hh TCalView.hh
  TCalVisNode.hh TCalVisNodeA.hh
  TCrvView.hh TCrvVisNode.hh
  TEvdCalSection.hh TEvdCluster.hh TEvdComboHit.hh TEvdCrvBar.hh TEvdCrvSection.hh
  TEvdCrystal.hh TEvdHelixSeed.hh TEvdHelixVisNode.hh TEvdMainFrame.hh
  TEvdPanel.hh TEvdPanelVisNode.hh TEvdPlane.hh
  TEvdSimParticle.hh TEvdStation.hh TEvdStrawHit.hh TEvdStraw.hh
  TEvdTracker.hh TEvdTimeCluster.hh TEvdTimeClusterVisNode.hh
  TEvdTrack.hh TEvdCosmicTrack.hh
  TEvdTrkStrawHit.hh TEventDisplayUtils.hh THeaderVisNode.hh TMcTruthVisNode.hh
  TMu2eBField.hh TEvdFrame.hh TStnGeoManager.hh TStnVisManager.hh
  TStnVisNode.hh
  TStnVisNodeA.hh
  TTrkVisNode.hh
  TEvdPanelVisNode.hh
)

ROOT_GENERATE_DICTIONARY(${dict} ${include_files} LINKDEF ${subpkg}/Linkdef.h)

set_property(SOURCE ${dict}.cxx PROPERTY generated)

set (source_files
  DisplayMu2eGeometry.cc
  TCalView.cc TCalVisNode.cc TCalVisNodeA.cc
  TCrvView.cc TCrvVisNode.cc
  TEvdCalSection.cc TEvdCluster.cc TEvdComboHit.cc TEvdCrvBar.cc TEvdCrvSection.cc
  TEvdCrystal.cc
  TEvdHelixSeed.cc TEvdHelixVisNode.cc
  TEvdMainFrame.cc
  TEvdPanel.cc TEvdPanelVisNode.cc 
  TEvdPlane.cc TEvdSimParticle.cc TEvdStation.cc TEvdStraw.cc TEvdStrawHit.cc
  TEvdTracker.cc
  TEvdTrkStrawHit.cc TEventDisplayUtils.cc THeaderVisNode.cc TMcTruthVisNode.cc
  TEvdTimeCluster.cc TEvdTimeClusterVisNode.cc
  TEvdTrack.cc TEvdCosmicTrack.cc
  TMu2eBField.cc TEvdFrame.cc TStnGeoManager.cc TStnVisManager.cc
  TStnVisNode.cc
  TStnVisNodeA.cc
  TTrkVisNode.cc
  TEvdPanelVisNode.cc
)

add_library(${lib} SHARED ${source_files} ${dict}.cxx)

target_link_libraries(${lib}
  Stntuple_print
  Stntuple_geom
  Stntuple_obj
  Stntuple_base
  
  Offline::RecoDataProducts
  Offline::MCDataProducts
  Offline::DataProducts
  Offline::CalorimeterGeom
  
  art::Framework_Core
  art::Framework_Services_Registry
  art::Framework_Principal
  art::Persistency_Common
  art::Persistency_Provenance
  art::Utilities

  ROOT::Physics ROOT::EG  ROOT::Eve ROOT::Graf ROOT::Graf3d ROOT::Gui ROOT::Gpad
  ROOT::Geom ROOT::Matrix ROOT::GenVector ROOT::Core
)

cet_add_compiler_flags(-g -O0)

install(TARGETS ${lib})
install_headers()
# install_source()
install(DIRECTORY ${subpkg} DESTINATION include/Stntuple)
