# Module providing convenience methods for compile binaries with psurface support.
#
# .. cmake_function:: add_dune_psurface_flags
#
#    .. cmake_param:: targets
#       :single:
#       :required:
#       :positional:
#
#       the targets to add the Psurface flags to.
#

# set HAVE_PSURFACE for config.h
set(HAVE_PSURFACE ${Psurface_FOUND})

# register all psurface related flags
if(Psurface_FOUND)
  dune_register_package_flags(LIBRARIES "Psurface::Psurface")

  dune_generate_pkg_config("psurface"
    NAME "Psurface"
    DESCRIPTION "Bijections between triangulated surfaces"
    URL "https://github.com/psurface/psurface"
    CFLAGS "-I${PSURFACE_INCLUDE_DIR}"
    LIBS "${PSURFACE_LIBRARY}")
  dune_add_pkg_config_requirement("psurface")
  dune_add_pkg_config_flags("-DHAVE_PSURFACE")
endif()

function(add_dune_psurface_flags _targets)
  if(Psurface_FOUND)
    foreach(_target ${_targets})
      target_link_libraries(${_target} PUBLIC Psurface::Psurface)
    endforeach(_target ${_targets})
  endif(Psurface_FOUND)
endfunction(add_dune_psurface_flags)
