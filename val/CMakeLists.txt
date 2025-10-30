#
set (subpkg "val"                    )
set (dict   "Stntuple_${subpkg}_dict")
set (lib    "Stntuple_${subpkg}"     )

include_directories(
  $ENV{ROOTSYS}/include/root
  ${CMAKE_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/${subpkg}
)

# link_directories($ENV{ROOTSYS}/lib)

set (include_files
  hist_data.hh hist_file.hh stn_book.hh stn_catalog.hh stn_dataset.hh
  TBadFolder.hh TGoodFolder.hh THistComp.hh
  stntuple_val_functions.hh
)

ROOT_GENERATE_DICTIONARY(${dict} ${include_files} LINKDEF ${subpkg}/Linkdef.h)

set_property(SOURCE ${dict}.cxx PROPERTY generated)

set (source_files
  compare_stn_hist.cc create_eff_hist.cc create_ntuple.cc folder_to_dir.cc get_hist.cc
  hist_data.cc new_slide.cc overlay_2h.cc plot_data.cc plot_hist.cc stn_book.cc
  stn_catalog.cc stn_dataset.cc TBadFolder.cc TGoodFolder.cc THistComp.cc
)

add_library(${lib} SHARED ${source_files} ${dict}.cxx)

target_link_libraries(${lib}
  Stntuple_obj
  Stntuple_base
  ROOT::GenVector ROOT::Graf ROOT::Graf3d ROOT::Geom ROOT::Matrix ROOT::Core)

cet_add_compiler_flags(-g -O0)

install(TARGETS ${lib})
# install_headers()
# install_source()
install(DIRECTORY ${subpkg} DESTINATION include/Stntuple)
