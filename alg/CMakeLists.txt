#
SET (CMAKE_Fortran_COMPILER  gfortran)
ENABLE_LANGUAGE(Fortran)

set (subpkg "alg")
set (dict "Stntuple_${subpkg}_dict")
set (lib  "Stntuple_${subpkg}")

include_directories(
  $ENV{ROOTSYS}/include/root
  ${CMAKE_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/${subpkg}
)

# link_directories($ENV{ROOTSYS}/lib)

ROOT_GENERATE_DICTIONARY(${dict}
  smooth.hh TEmuLogLH.hh TStnTrackID.hh TStntuple.hh
  LINKDEF ${subpkg}/Linkdef.h)

set_property(SOURCE ${dict}.cxx PROPERTY generated)

add_library(${lib} SHARED
  InitPbarCommon.cc RMCSpectra.cc smooth.cc TEmuLogLH.cc TStnTrackID.cc TStntuple.cc TStntuple_Pbar_Striganov.cc
  pbar_common_address.f pbar.f 
  ${dict}.cxx
)

target_link_libraries(${lib} ROOT::Core ROOT::Matrix ROOT::Tree ROOT::Graf ROOT::Gpad)

cet_add_compiler_flags(-g -O0)

install(TARGETS ${lib})
# install_headers()
# install_source()
install(DIRECTORY ${subpkg} DESTINATION include/Stntuple)
