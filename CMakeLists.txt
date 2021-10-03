cmake_minimum_required(VERSION 3.13)
project("dune-grid" C CXX)

if(NOT (dune-common_DIR OR dune-common_ROOT OR
        "${CMAKE_PREFIX_PATH}" MATCHES ".*dune-common.*"))
    string(REPLACE  ${PROJECT_NAME} dune-common dune-common_DIR
      ${PROJECT_BINARY_DIR})
endif()

#find dune-common and set the module path
find_package(dune-common REQUIRED)

# start a dune project with information from dune.module
dune_project(dune-grid)

if(Alberta_FOUND)
  foreach(_dim ${ALBERTA_WORLD_DIMS})
    dune_add_library(dunealbertagrid${_dim}d EXPORT_NAME albertagrid${_dim}d)
    target_link_libraries(dunealbertagrid${_dim}d PUBLIC Dune::grid Alberta::AlbertaGrid_${_dim}d)
    target_compile_definitions(dunealbertagrid${_dim}d PRIVATE ALBERTA_DIM=${_dim} ENABLE_ALBERTA=1)
  endforeach(_dim)
endif()

add_subdirectory("dune")
add_subdirectory("doc")
add_subdirectory("lib")
add_subdirectory("src")
add_subdirectory("cmake/modules")

# if Python bindings are enabled, include necessary sub directories.
if( DUNE_ENABLE_PYTHONBINDINGS )
  add_subdirectory("python")
  dune_python_install_package(PATH python)
endif()

# Use a  custom section to export path to the grid examples.
set(DUNE_CUSTOM_PKG_CONFIG_SECTION "#Export the directory with the grid example for downstream modules
if(${ProjectName}_INSTALLED)
  set(DUNE_GRID_EXAMPLE_GRIDS_PATH \"\${PACKAGE_PREFIX_DIR}/${CMAKE_INSTALL_DOCDIR}/grids/\")
else()
  set(DUNE_GRID_EXAMPLE_GRIDS_PATH \"\${PACKAGE_PREFIX_DIR}/doc/grids/\")
endif()
")
# finalize the dune project, e.g. generating config.h etc.
finalize_dune_project(GENERATE_CONFIG_H_CMAKE)
