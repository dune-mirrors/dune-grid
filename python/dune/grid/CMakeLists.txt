add_python_targets(grid
  __init__
  __main__
  _grids
  grid_generator
  core
  map
  datahandle
)
dune_add_pybind11_module(NAME _grid)
target_link_libraries(_grid PUBLIC Dune::grid)
install(TARGETS _grid LIBRARY DESTINATION python/dune/grid)
