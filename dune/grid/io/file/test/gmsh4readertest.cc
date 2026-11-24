// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later

// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include <cmath>
#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/test/testsuite.hh>

#include <dune/grid/io/file/gmshreader.hh>
#include <dune/grid/io/file/gmsh/gmsh4reader.hh>

#include <dune/grid/test/gridcheck.hh>

#include <dune/grid/uggrid.hh>

using namespace Dune;

// Compare the number of entities of all codimensions and the sum of all vertex coordinates
template <class GridA, class GridB>
bool equalGrids (const GridA& gridA, const GridB& gridB)
{
  auto gvA = gridA.leafGridView();
  auto gvB = gridB.leafGridView();
  for (int codim = 0; codim <= GridA::dimension; ++codim)
    if (gvA.size(codim) != gvB.size(codim))
      return false;

  auto coordinateSum = [](const auto& gridView) {
    double sum = 0.0;
    for (const auto& v : vertices(gridView))
      for (const auto& x_i : v.geometry().corner(0))
        sum += x_i;
    return sum;
  };
  using std::abs;
  return abs(coordinateSum(gvA) - coordinateSum(gvB)) < 1.e-10;
}

// Read the ascii grid `name`-ascii.msh and compare it with the grids read from the
// binary files `name``suffix` for all given suffixes.
template <int dim>
void test (TestSuite& testSuite, const std::string& path, const std::string& name,
           const std::vector<std::string>& binarySuffixes)
{
  using Grid = UGGrid<dim>;
  using Reader = Impl::Gmsh::Gmsh4Reader<Grid>;
  using Reader32 = Impl::Gmsh::Gmsh4Reader<Grid, Impl::Gmsh::ContinuousGridCreator<Grid>, std::uint32_t>;

  std::cout << "Reading " << name << "-ascii.msh" << std::endl;
  auto ascii = Reader::createGridFromFile(path + name + "-ascii.msh");
  gridcheck(*ascii);

  for (const auto& suffix : binarySuffixes)
  {
    std::cout << "Reading " << name << suffix << std::endl;
    auto binary = Reader::createGridFromFile(path + name + suffix);
    gridcheck(*binary);
    testSuite.check(equalGrids(*ascii, *binary), name + suffix)
      << "Grid read from binary file differs from grid read from ascii file";

    // read with a size_type that differs from the data-size of most files
    auto binary32 = Reader32::createGridFromFile(path + name + suffix);
    testSuite.check(equalGrids(*ascii, *binary32), name + suffix + " [uint32_t]")
      << "Grid read with size_type uint32_t differs from grid read from ascii file";
  }
}

int main(int argc, char** argv)
{
  MPIHelper::instance(argc, argv);
  TestSuite testSuite;

  const std::string path = DUNE_GRID_EXAMPLE_GRIDS_PATH "gmsh/";

  // The *-datasize4* files are copies of the *-binary.msh files with all size_t
  // values stored with 4 bytes, as written by Gmsh on a 32-bit machine. The
  // *-bigendian file is additionally stored in big-endian byte order.
  test<2>(testSuite, path, "hybrid-testgrid-2d-v4",
          {"-binary.msh", "-binary-datasize4.msh"});
  test<3>(testSuite, path, "hybrid-testgrid-3d-v4",
          {"-binary.msh", "-binary-datasize4-bigendian.msh"});

  { // Read a grid with a $PhysicalNames section through the GmshReader interface
    std::cout << "Reading unitsquare-v4-physicalnames.msh" << std::endl;
    auto grid = GmshReader<UGGrid<2> >::read(path + "unitsquare-v4-physicalnames.msh");
    gridcheck(*grid);
    testSuite.check(grid->leafGridView().size(0) > 0, "unitsquare-v4-physicalnames.msh")
      << "Grid with physical names is empty";
  }

  return testSuite.exit();
}
