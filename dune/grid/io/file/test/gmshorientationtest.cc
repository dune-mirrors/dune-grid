// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later

// Test that all elements of grids read by the GmshReader have a positive
// orientation, i.e., a positive determinant of the Jacobian of the element
// geometry, for all element types in the hybrid test grids (triangles,
// quadrilaterals, tetrahedra, hexahedra, prisms, pyramids). The grids are read
// from files in Gmsh format version 2 and version 4 (ascii and binary).

#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include <cmath>
#include <iostream>
#include <map>
#include <string>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/test/testsuite.hh>
#include <dune/geometry/referenceelements.hh>

#include <dune/grid/io/file/gmshreader.hh>
#include <dune/grid/uggrid.hh>

using namespace Dune;

// Check the orientation of all elements and return the volume of the grid
template <int dim>
double checkOrientation (TestSuite& testSuite, const std::string& path, const std::string& filename)
{
  auto grid = GmshReader<UGGrid<dim> >::read(path + filename, /*verbose*/ false);

  std::map<GeometryType, int> numPositive, numNegative;
  double volume = 0.0;
  for (const auto& element : elements(grid->leafGridView())) {
    const auto geometry = element.geometry();
    const auto center = referenceElement(geometry).position(0, 0);
    if (geometry.jacobianTransposed(center).determinant() > 0.0)
      ++numPositive[element.type()];
    else
      ++numNegative[element.type()];
    volume += geometry.volume();
  }

  for (const auto& [type, n] : numPositive)
    std::cout << filename << ": " << n << " " << type << " with positive orientation" << std::endl;
  for (const auto& [type, n] : numNegative)
    testSuite.check(false, filename) << n << " " << type << " with negative orientation";

  return volume;
}

template <int dim>
void test (TestSuite& testSuite, const std::string& path, const std::string& name)
{
  const double volume = checkOrientation<dim>(testSuite, path, name + ".msh");
  for (std::string suffix : {"-v4-ascii.msh", "-v4-binary.msh"}) {
    const double volumeV4 = checkOrientation<dim>(testSuite, path, name + suffix);
    testSuite.check(std::abs(volume - volumeV4) < 1.e-10, name + suffix)
      << "volume " << volumeV4 << " differs from volume " << volume << " of the version 2 grid";
  }
}

int main (int argc, char** argv)
{
  MPIHelper::instance(argc, argv);
  TestSuite testSuite;

  const std::string path = DUNE_GRID_EXAMPLE_GRIDS_PATH "gmsh/";

  test<2>(testSuite, path, "hybrid-testgrid-2d");
  test<3>(testSuite, path, "hybrid-testgrid-3d");

  return testSuite.exit();
}
