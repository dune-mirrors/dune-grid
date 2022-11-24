# SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
# SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

add_subdirectory(gridinfo-gmsh)

add_executable(periodic-grid periodic-grid.cc)
target_link_libraries(periodic-grid PRIVATE dunegrid)
add_dune_alberta_flags(periodic-grid WORLDDIM 2)