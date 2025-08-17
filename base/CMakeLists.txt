#
set (subpkg "base"                   )
set (dict   "Stntuple_${subpkg}_dict")
set (lib    "Stntuple_${subpkg}"     )

include_directories(
  $ENV{ROOTSYS}/include/root
  ${CMAKE_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/base
)

# link_directories($ENV{ROOTSYS}/lib)

ROOT_GENERATE_DICTIONARY(${dict}
  TBeamConditions.hh TBitset.hh TCalibData.hh TCalibManager.hh TCalibRunLists.hh
  TCalibRunRange.hh TCalibTable.hh TCdf2Filesets.hh TCdf2Files.hh TCdf2RunSectionRanges.hh
  TCdf2RunSections.hh THexIndex.hh TMatrix33.hh TMatrix55.hh TNamedHandle.hh TObjHandle.hh
  TOffscl.hh TSetRunMaps.hh TStnArrayF.hh TStnArrayI.hh TStnDataset.hh TStnFileset.hh TStnHexagon.hh
  TStnRunRecord.hh TStnShape.hh TStnSquare.hh TStnUtils.hh TStnView.hh TTriggerConditions.hh
  TUsedSet.hh TVisManager.hh TVisNode.hh TVisNodeA.hh
  LINKDEF base/Linkdef.h)

set_property(SOURCE ${dict}.cxx PROPERTY generated)

add_library(${lib} SHARED
  TBeamConditions.cc  TCalibRunLists.cc  TCdf2Filesets.cc
  TMatrix33.cc     TOffscl.cc      TStnDataset.cc    TStnShape.cc
  TTriggerConditions.cc TBitset.cc TCalibRunRange.cc TCdf2RunSectionRanges.cc
  TMatrix55.cc     TSetRunMaps.cc  TStnFileset.cc    TStnSquare.cc  TUsedSet.cc
  TCalibData.cc TCalibTable.cc     TCdf2RunSections.cc TNamedHandle.cc
  TStnArrayF.cc   TStnHexagon.cc    TStnUtils.cc   TVisManager.cc TCalibManager.cc
  TCdf2Files.cc THexIndex.cc  TObjHandle.cc TStnArrayI.cc   TStnRunRecord.cc
  TStnView.cc    TVisNode.cc TVisNodeA.cc
  ${dict}.cxx
)

target_link_libraries(${lib} ROOT::Physics ROOT::Core ROOT::Matrix ROOT::Tree ROOT::Graf ROOT::Gpad  ROOT::Geom)

cet_add_compiler_flags(-g -O0)

install(TARGETS Stntuple_base)
install_headers()
# install_source()
install(DIRECTORY ${subpkg} DESTINATION include/Stntuple)
