#
set (subpkg "stat")
set (dict   "Stntuple_${subpkg}_dict")
set (lib    "Stntuple_${subpkg}")

include_directories(
  $ENV{ROOTSYS}/include/root
  ${CMAKE_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/${subpkg}
)

# link_directories($ENV{ROOTSYS}/lib)

set (include_files
  crow_gardner.hh FCCalculator.hh model_t.hh Mu2e_model.hh
  parameter_t.hh pfixed_t.hh pflat_t.hh pgaus_t.hh plogn_t.hh ppoi_t.hh
  TBelt.hh TBeltLH.hh TFeldmanCousins.hh TKinLH.hh upper_limit.hh
)

ROOT_GENERATE_DICTIONARY(${dict} ${include_files} LINKDEF ${subpkg}/Linkdef.h)

set_property(SOURCE ${dict}.cxx PROPERTY generated)

set (source_files
  channel_t.cc crow_gardner.cc FCCalculator.cc Functions.cc model_t.cc Mu2e_model.cc
  parameter_t.cc pfixed_t.cc pflat_t.cc pgaus_t.cc plogn_t.cc Poisson_t.cc ppoi_t.cc
  TBelt.cc TBeltLH.cc TFeldmanCousins.cc TKinLH.cc upper_limit.cc var_t.cc
)

add_library(${lib} SHARED ${source_files} ${dict}.cxx)

target_link_libraries(${lib}
  Stntuple_val
  Stntuple_obj
  Stntuple_base
  ROOT::MathMore ROOT::GenVector ROOT::Graf ROOT::Graf3d ROOT::Geom ROOT::Matrix ROOT::Core)

cet_add_compiler_flags(-g -O0)

install(TARGETS ${lib})
install_headers()
# install_source()
