# OBJECT is needed since CMake 3.0.0
# but it conflicts with shared libs
set(_OBJECT_FLAG "OBJECT")
if(BUILD_SHARED_LIBS)
  set(_OBJECT_FLAG "")
endif()

if(ALBERTA_FOUND)
foreach(_dim ${ALBERTA_WORLD_DIMS})
  dune_add_library(dunealbertagrid_${_dim}d APPEND
    _DUNE_TARGET_OBJECTS:albertagrid_${_dim}d_
    ADD_LIBS ${DUNE_LIBS})
  add_dune_alberta_flags(dunealbertagrid_${_dim}d ${_OBJECT_FLAG} NO_LINK_DUNEALBERTAGRID GRIDDIM ${_dim})
  list(APPEND DUNE_ALBERTA_LIBS dunealbertagrid_${_dim}d)
endforeach(_dim "${ALBERTA_DIMS}")
endif(ALBERTA_FOUND)

install(FILES dunegridam2cmake.lib
  DESTINATION ${CMAKE_INSTALL_BINDIR}/../lib)
