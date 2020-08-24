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

include_guard(GLOBAL)

# set HAVE_PSURFACE for config.h
set(HAVE_PSURFACE ${Psurface_FOUND})

function(add_dune_psurface_flags _targets)
  if(Psurface_FOUND)
    foreach(_target ${_targets})
      target_link_libraries(${_target} PUBLIC ${PSURFACE_DUNE_LIBRARIES})
      target_compile_options(${_target} PUBLIC ${PSURFACE_DUNE_COMPILE_FLAGS})
    endforeach(_target ${_targets})
  endif(Psurface_FOUND)
endfunction(add_dune_psurface_flags)
