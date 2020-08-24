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
#    .. cmake_param:: SOURCE_ONLY
#       :option:
#
#       TODO doc me
#       old doc: indicates that the targets are source files.
#
#    .. cmake_param:: OBJECT
#       :option:
#
#       TODO doc me
#       old doc: indicates that the targets are object libraries.
#

# Add dgf magic to config.h and register flags
if(dune-uggrid_FOUND)
  dune_define_gridtype(GRID_CONFIG_H_BOTTOM GRIDTYPE UGGRID
    ASSERTION "GRIDDIM == WORLDDIM"
    DUNETYPE "Dune::UGGrid< dimgrid >"
    HEADERS "dune/grid/uggrid.hh" "dune/grid/io/file/dgfparser/dgfug.hh")
  set(HAVE_UG TRUE)
endif()

# Add flags to targets
function(add_dune_ug_flags)
  if(dune-uggrid_FOUND)
    foreach(_target ${ARGV})
      target_link_libraries(${_target} PUBLIC Dune::dune-uggrid ugL)
      # link dimension dependent targets
      foreach (dim ${UG_ENABLED_DIMENSIONS})
        target_link_libraries(${_target} PUBLIC ugS${dim})
      endforeach (dim)
      target_compile_definitions(${_target} PUBLIC "ENABLE_UG=1;${UG_DEFINITIONS}")
    endforeach()
  endif()
endfunction(add_dune_ug_flags)
