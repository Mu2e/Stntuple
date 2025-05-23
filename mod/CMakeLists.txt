#
set (subpkg "mod"                    )
set (dict   "Stntuple_${subpkg}_dict")
set (lib    "Stntuple_${subpkg}"     )

set (MuHitDisplayModule "Stntuple_MuHitDisplay_module" )
set (MuHitDisplayDict   "Stntuple_MuHitDisplay_module_dict" )

include_directories(
  $ENV{ROOTSYS}/include/root
  ${CMAKE_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/${subpkg}
)

# link_directories($ENV{ROOTSYS}/lib)

set (include_files
  # StntupleRecoAna_module.hh
  THistModule.hh TModule.hh
  TAnaRint.hh
)

ROOT_GENERATE_DICTIONARY(${MuHitDisplayDict} MuHitDisplay_module.hh
  LINKDEF ${subpkg}/MuHitDisplay_module_Linkdef.h
  OPTIONS -D__CLING__
)
set_property(SOURCE ${MuHitDisplayDict}.cxx PROPERTY generated)

ROOT_GENERATE_DICTIONARY(${dict} ${include_files}
  LINKDEF ${subpkg}/Linkdef.h
  OPTIONS -D__CLING__
)
set_property(SOURCE ${dict}.cxx PROPERTY generated)

cet_build_plugin(EventFilter art::module LIBRARIES REG
  Stntuple_mod
  Stntuple_gui
  Stntuple_print
  Stntuple_geom
  Stntuple_obj
  Stntuple_base
  art_root_io::TFileService_service
  artdaq-core-mu2e::Data
  ROOT::Hist
  ROOT::Tree
  ROOT::Core
  ROOT::RIO
  ROOT::Gui
  Offline::RecoDataProducts
)

cet_build_plugin(FillStntuple art::module LIBRARIES REG
  Stntuple_mod
  Stntuple_gui
  Stntuple_print
  Stntuple_geom
  Stntuple_obj
  Stntuple_base
  art_root_io::TFileService_service
  artdaq-core-mu2e::Data
  ROOT::Hist
  ROOT::Tree
  ROOT::Core
  ROOT::RIO
  ROOT::Gui
  Offline::RecoDataProducts
)

cet_build_plugin(InitStntuple art::module LIBRARIES REG
  Stntuple_mod
  Stntuple_gui
  Stntuple_print
  Stntuple_geom
  Stntuple_alg
  Stntuple_obj
  Stntuple_base
  art_root_io::TFileService_service
  artdaq-core-mu2e::Data
  ROOT::Hist
  ROOT::Tree
  ROOT::Core
  ROOT::RIO
  ROOT::Gui
  Offline::RecoDataProducts
)

add_library(${MuHitDisplayModule} SHARED MuHitDisplay_module.cc  ${MuHitDisplayDict}.cxx)
  
target_link_libraries(${MuHitDisplayModule} 
  PRIVATE
  Stntuple_mod
  Stntuple_gui
  Stntuple_print
  Stntuple_alg
  Stntuple_geom
  Stntuple_obj
  Stntuple_base
  Offline::RecoDataProducts
  art_root_io::TFileService_service
  artdaq-core-mu2e::Data
  ROOT::Hist
  ROOT::Tree
  ROOT::Core
  ROOT::RIO
  ROOT::Gui
)

cet_build_plugin(StntupleEventDump art::module LIBRARIES REG
  Stntuple_mod
  Stntuple_gui
  Stntuple_print
  Stntuple_geom
  Stntuple_obj
  Stntuple_base
  art_root_io::TFileService_service
  artdaq-core-mu2e::Data
  ROOT::Hist
  ROOT::Tree
  ROOT::Core
  ROOT::RIO
  ROOT::Gui
  Offline::RecoDataProducts
)

cet_build_plugin(StntupleMaker art::module LIBRARIES REG
  Stntuple_mod
  Stntuple_gui
  Stntuple_print
  Stntuple_geom
  Stntuple_obj
  Stntuple_base
  art_root_io::TFileService_service
  artdaq-core-mu2e::Data
  ROOT::Hist
  ROOT::Tree
  ROOT::Core
  ROOT::RIO
  ROOT::Gui
  Offline::RecoDataProducts
)

cet_build_plugin(StntupleRecoAna art::module LIBRARIES REG
  Stntuple_mod
  Stntuple_gui
  Stntuple_print
  Stntuple_geom
  Stntuple_obj
  Stntuple_base
  art_root_io::TFileService_service
  artdaq-core-mu2e::Data
  ROOT::Hist
  ROOT::Tree
  ROOT::Core
  ROOT::RIO
  ROOT::Gui
  Offline::RecoDataProducts
)

set (source_files
  autogenerated.cc
  addStntupleDataBlocks.cc 
  InitCalDataBlock.cc InitClusterBlock.cc InitCrvClusterBlock.cc InitCrvPulseBlock.cc
  InitGenpBlock.cc InitHeaderBlock.cc InitHelixBlock.cc InitPidBlock.cc InitSimpBlock.cc
  InitStepPointMCBlock.cc
  InitStrawHitBlock.cc InitTimeClusterBlock.cc
  InitTrackBlock.cc InitTrackBlock_KK.cc InitTrackSeedBlock.cc InitTrackStrawHitBlock.cc
  InitTriggerBlock.cc InitVDetDataBlock.cc
  StntupleModule.cc
  StntupleUtilities.cc TAnaRint.cc THistModule.cc TModule.cc
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

  ROOT::Physics ROOT::GenVector ROOT::Eve ROOT::Graf ROOT::Graf3d ROOT::Gui ROOT::Gpad
  ROOT::Geom ROOT::Matrix ROOT::Rint ROOT::Core
)

cet_add_compiler_flags(-g -O0)

install(TARGETS ${lib} ${MuHitDisplayModule})

# install_headers()
# install_source()

message( " --Stntuple/${subpkg}: CMAKE_INSTALL_PREFIX: " ${CMAKE_INSTALL_PREFIX})

install(DIRECTORY ${subpkg} DESTINATION include/Stntuple)
