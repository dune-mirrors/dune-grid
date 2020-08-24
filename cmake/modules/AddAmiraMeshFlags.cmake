# Module providing convenience methods for compile binaries with AmiraMesh support.
#
# .. cmake_function:: add_dune_amiramesh_flags
#
#    .. cmake_param:: targets
#       :single:
#       :positional:
#       :required:
#
#       The targets to add the amiramesh flags to.
#

include_guard(GLOBAL)

# set HAVE_AMIRAMESH for config.h
set(HAVE_AMIRAMESH ${AmiraMesh_FOUND})

function(add_dune_amiramesh_flags _targets)
  if(AmiraMesh_FOUND)
    foreach(_target ${_targets})
      target_link_libraries(${_target} PUBLIC ${AMIRAMESH_LIBRARIES})
      target_compile_options(${_target} PUBLIC ${AMIRAMESH_COMPILE_FLAGS})
    endforeach(_target ${_targets})
  endif(AmiraMesh_FOUND)
endfunction(add_dune_amiramesh_flags)
