#
set (subpkg "geom"                   )
set (dict   "Stntuple_${subpkg}_dict")
set (lib    "Stntuple_${subpkg}"     )

include_directories(
  $ENV{ROOTSYS}/include/root
  ${CMAKE_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/${subpkg}
)

# link_directories($ENV{ROOTSYS}/lib)

set (include_files
  stntuple_geom_tests.hh TCrvNumerology.hh
  TDiskCalorimeter.hh TDiskIndex.hh TDisk.hh TStnCrystal.hh
  TDiskCrystalMap.hh
  THexCrystalMap.hh TSqrCrystalMap.hh TSsqCrystalMap.hh 
  TStnTracker.hh
  TCrvChannelMap.hh
  TrkPanelMap.hh
  TCaloChannelMap.hh
)

ROOT_GENERATE_DICTIONARY(${dict} ${include_files} LINKDEF ${subpkg}/Linkdef.h)

set_property(SOURCE ${dict}.cxx PROPERTY generated)

set (source_files
  StnInitVirtualDetectors.cc stntuple_geom_tests.cc TCrvNumerology.cc
  TDiskCalorimeter.cc TDisk.cc TDiskCrystalMap.cc test_TDisk.cc
  THexCrystalMap.cc TSqrCrystalMap.cc TSsqCrystalMap.cc
  TStnCrystal.cc
  TStnTracker.cc
  TCrvChannelMap.cc
  TrkPanelMap.cc
  TCaloChannelMap.cc
)

add_library(${lib} SHARED ${source_files} ${dict}.cxx)

target_link_libraries(${lib}
  Stntuple_base
  fhiclcpp  cetlib
  ROOT::Physics ROOT::Matrix ROOT::Tree ROOT::Geom ROOT::Core)

cet_add_compiler_flags(-g -O0)

install(TARGETS ${lib})
# install_headers()
# install_source()
install(DIRECTORY ${subpkg} DESTINATION include/Stntuple)
