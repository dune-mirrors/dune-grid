# SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
# SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

add_subdirectory(tutorial)
add_python_targets(grid
  __init__
  __main__
  _grids
  grid_generator
  core
  map
  datahandle
  mesh2dgf
)
dune_add_pybind11_module(NAME _grid)
set_property(TARGET _grid PROPERTY LINK_LIBRARIES Dune::Grid APPEND)
# pass the path to the example grids to the Python module as a compile definition
if(DUNE_GRID_EXAMPLE_GRIDS_PATH)
  target_compile_definitions(_grid PRIVATE DUNE_GRID_EXAMPLE_GRIDS_PATH=${DUNE_GRID_EXAMPLE_GRIDS_PATH})
else()
  target_compile_definitions(_grid PRIVATE DUNE_GRID_EXAMPLE_GRIDS_PATH=\"${PROJECT_SOURCE_DIR}/doc/grids/\")
endif()
if(SKBUILD)
  install(TARGETS _grid LIBRARY DESTINATION python/dune/grid)
endif()
