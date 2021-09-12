# If dune-uggrid was found this module adds the dgf magic to config.h
# and makes add_dune_ug_flags available.
#
# .. cmake_function:: add_dune_ug_flags
#
#    .. cmake_param:: targets
#       :single:
#       :positional:
#       :required:
#
#       The targets to add the UG flags to.
#

# Add dgf magic to config.h and register flags
if(dune-uggrid_FOUND)
  dune_register_package_flags(COMPILE_DEFINITIONS "ENABLE_UG=1"
                              LIBRARIES "duneuggrid")

  dune_define_gridtype(GRID_CONFIG_H_BOTTOM GRIDTYPE UGGRID
    ASSERTION "GRIDDIM == WORLDDIM"
    DUNETYPE "Dune::UGGrid< dimgrid >"
    HEADERS "dune/grid/uggrid.hh" "dune/grid/io/file/dgfparser/dgfug.hh")
endif()

# Add flags to targets
function(add_dune_ug_flags _target)
  if(dune-uggrid_FOUND)
    target_link_libraries(${_target} PUBLIC duneuggrid)
  endif()
endfunction(add_dune_ug_flags)
