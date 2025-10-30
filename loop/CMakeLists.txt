#
set (subpkg "loop")
set (dict   "Stntuple_${subpkg}_dict")
set (lib    "Stntuple_${subpkg}")

include_directories(
  $ENV{ROOTSYS}/include/root
  ${CMAKE_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/${subpkg}
)

# link_directories($ENV{ROOTSYS}/lib)

set (include_files
  TDFCModule.hh TDuplicateFilterModule.hh THttpCatalogServer.hh
  TSamCatalogServer.hh TStnAna.hh TStnCatalog.hh TStnCatalogServer.hh
  TStnInputModule.hh TStnModule.hh TStnOutputModule.hh
  TStnRun2InputModule.hh TTxtCatalogServer.hh
)
if ($ENV{STNTUPLE_MC_GEN})
  append (include_files TStnGeneratorModule.hh)
endif()

ROOT_GENERATE_DICTIONARY(${dict} ${include_files} LINKDEF ${subpkg}/Linkdef.h)

set_property(SOURCE ${dict}.cxx PROPERTY generated)

set (source_files
  StnGenInitHeaderBlock.cc TDFCModule.cc TDuplicateFilterModule.cc
  THttpCatalogServer.cc TSamCatalogServer.cc TStnAna.cc TStnCatalog.cc TStnCatalogServer.cc
  TStnInputModule.cc TStnModule.cc TStnOutputModule.cc TStnRun2InputModule.cc
  TTxtCatalogServer.cc
)
if ($ENV{STNTUPLE_MC_GEN})
  append (source_files TStnGeneratorModule.cc StnGenInitGenpBlock.cc )
endif()

add_library(${lib} SHARED ${source_files}     ${dict}.cxx)
add_library(${CMAKE_PROJECT_NAME}::${subpkg} ALIAS ${lib})

message(" -- P.Murat CMAKE_PROJECT_NAME::subpkg ${CMAKE_PROJECT_NAME}::${subpkg}")

target_link_libraries(${lib}
  Stntuple_obj
  Stntuple_base
  ROOT::Physics ROOT::GenVector ROOT::Graf ROOT::Graf3d ROOT::Geom ROOT::Matrix ROOT::Core)

cet_add_compiler_flags(-g -O0)

install(TARGETS ${lib})

# install_headers()
# install_source()

install(DIRECTORY ${subpkg} DESTINATION include/Stntuple)
