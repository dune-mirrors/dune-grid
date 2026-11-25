// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

/** \file
 * \brief Test the element and boundary data provided by the GmshReader
 *
 * The documentation of `GmshReader::elementData()` and `GmshReader::boundaryData()`
 * states that the vectors map the element index and the boundary segment index to the
 * Gmsh physical entity. Actually, the vectors are indexed by the insertion index of the
 * grid factory. See https://gitlab.dune-project.org/core/dune-grid/-/issues/197
 *
 * The meshes used here:
 * - `partialboundary.msh`: The unit square split into the physical surfaces 1 (x < 0.5)
 *   and 2 (x > 0.5). Only 3 of the 8 boundary edges are in the file: the left side with
 *   physical tag 10 and the right half of the bottom side with physical tag 20. The
 *   element order in the file does not follow the geometry.
 * - `partialboundary-v4.msh`: The same mesh in the Gmsh 4.1 format.
 * - `interiorinterface.msh`: The unit square split along the diagonal into the physical
 *   surfaces 1 (above) and 2 (below). The diagonal is a physical curve with tag 30, i.e.,
 *   an interior interface, and the left side has the physical tag 10.
 */

#include <config.h>

#include <cmath>
#include <exception>
#include <iostream>
#include <string>
#include <vector>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/test/testsuite.hh>

#include <dune/grid/common/gridfactory.hh>
#include <dune/grid/common/mcmgmapper.hh>
#include <dune/grid/io/file/gmshreader.hh>

#if GMSH_UGGRID
#include <dune/grid/uggrid.hh>
#endif

#if GMSH_ALBERTAGRID
#include <dune/grid/albertagrid.hh>
#include <dune/grid/albertagrid/gridfactory.hh>
#endif

using namespace Dune;

const std::string path = std::string(DUNE_GRID_EXAMPLE_GRIDS_PATH) + "gmsh/";

constexpr auto readerOptions = Gmsh::ReaderOptions::insertBoundarySegments
                             | Gmsh::ReaderOptions::readElementData
                             | Gmsh::ReaderOptions::readBoundaryData;

// The physical tag of a boundary intersection of the mesh `partialboundary.msh`,
// or 0 if the boundary segment is not in the file
template <class Intersection>
int expectedBoundaryTag (const Intersection& intersection)
{
  auto center = intersection.geometry().center();
  if (std::abs(center[0]) < 1e-8)
    return 10;
  if (std::abs(center[1]) < 1e-8 && center[0] > 0.5)
    return 20;
  return 0;
}


// The status quo: the data is indexed by the insertion index of the factory.
// This works if the factory is available, but only for the inserted boundary segments.
template <class Grid>
TestSuite testInsertionIndexing (const std::string& fileName)
{
  TestSuite t("insertion indexing (" + fileName + ")");
  try {
    GridFactory<Grid> factory;
    GmshReader<Grid> reader(path + fileName, factory, readerOptions);
    auto grid = factory.createGrid();
    auto gridView = grid->leafGridView();

    t.check(reader.hasElementData(), "reader has element data");
    t.check(reader.hasBoundaryData(), "reader has boundary data");
    if (!reader.hasElementData() || !reader.hasBoundaryData())
      return t;
    const auto& elementData = reader.elementData();
    const auto& boundaryData = reader.boundaryData();

    std::size_t numInserted = 0;
    for (const auto& element : elements(gridView)) {
      int expected = element.geometry().center()[0] < 0.5 ? 1 : 2;
      auto index = factory.insertionIndex(element);
      t.check(index < elementData.size(), "element insertion index in range");
      if (index < elementData.size())
        t.check(elementData[index] == expected, "element data by insertion index")
          << "got " << elementData[index] << ", expected " << expected;

      for (const auto& intersection : intersections(gridView, element)) {
        if (!intersection.boundary() || !factory.wasInserted(intersection))
          continue;
        ++numInserted;
        auto bIndex = factory.insertionIndex(intersection);
        t.check(bIndex < boundaryData.size(), "boundary insertion index in range");
        if (bIndex < boundaryData.size())
          t.check(boundaryData[bIndex] == expectedBoundaryTag(intersection), "boundary data by insertion index")
            << "got " << boundaryData[bIndex] << ", expected " << expectedBoundaryTag(intersection);
      }
    }
    t.check(numInserted == 3, "number of inserted boundary segments")
      << numInserted << " inserted, expected 3";
  }
  catch (const std::exception& e) {
    t.check(false, "no exception") << e.what();
  }
  return t;
}


// The documented behavior: the data is indexed by the element index and by the
// boundary segment index of the grid.
template <class Grid>
TestSuite testDocumentedIndexing (const std::string& fileName)
{
  TestSuite t("documented indexing (" + fileName + ")");
  try {
    GmshReader<Grid> reader(path + fileName, readerOptions);
    auto grid = reader.createGrid();
    auto gridView = grid->leafGridView();
    MultipleCodimMultipleGeomTypeMapper mapper(gridView, mcmgElementLayout());

    const auto& elementData = reader.elementData();
    const auto& boundaryData = reader.boundaryData();

    t.check(elementData.size() == mapper.size(), "size of element data")
      << elementData.size() << " entries for " << mapper.size() << " elements";
    t.check(boundaryData.size() == grid->numBoundarySegments(), "size of boundary data")
      << boundaryData.size() << " entries for " << grid->numBoundarySegments() << " boundary segments";

    for (const auto& element : elements(gridView)) {
      int expected = element.geometry().center()[0] < 0.5 ? 1 : 2;
      auto index = mapper.index(element);
      t.check(index < elementData.size(), "element index in range");
      if (index < elementData.size())
        t.check(elementData[index] == expected, "element data by element index")
          << "got " << elementData[index] << ", expected " << expected;

      for (const auto& intersection : intersections(gridView, element)) {
        if (!intersection.boundary())
          continue;
        auto bIndex = intersection.boundarySegmentIndex();
        int expectedTag = expectedBoundaryTag(intersection);
        t.check(bIndex < boundaryData.size(), "boundary segment index in range")
          << "segment " << bIndex << " (expected tag " << expectedTag << ") is out of range";
        if (bIndex < boundaryData.size())
        {
          // Which value a segment without data should have is not specified yet,
          // but it must not be the tag of another segment
          if (expectedTag != 0)
            t.check(boundaryData[bIndex] == expectedTag, "boundary data by boundary segment index")
              << "segment " << bIndex << ": got " << boundaryData[bIndex] << ", expected " << expectedTag;
          else
            t.check(boundaryData[bIndex] != 10 && boundaryData[bIndex] != 20, "no data for segments not in the file")
              << "segment " << bIndex << " is not in the file, but got the tag " << boundaryData[bIndex];
        }
      }
    }
  }
  catch (const std::exception& e) {
    t.check(false, "no exception") << e.what();
  }
  return t;
}


// Codim-1 elements in the interior of the domain, e.g., on the interface between two
// physical surfaces, are not boundary segments and must not be inserted as such.
template <class Grid>
TestSuite testInteriorInterface ()
{
  TestSuite t("interior interface (interiorinterface.msh)");
  try {
    GridFactory<Grid> factory;
    GmshReader<Grid> reader(path + "interiorinterface.msh", factory, readerOptions);
    auto grid = factory.createGrid();
    auto gridView = grid->leafGridView();

    std::size_t numInserted = 0;
    for (const auto& element : elements(gridView))
      for (const auto& intersection : intersections(gridView, element))
        if (intersection.boundary() && factory.wasInserted(intersection))
          ++numInserted;

    t.check(numInserted == 1, "number of inserted boundary segments")
      << numInserted << " inserted, expected 1 (the left side)";
    t.check(reader.boundaryData().size() == numInserted, "every boundary data entry belongs to a boundary segment")
      << reader.boundaryData().size() << " entries for " << numInserted << " inserted boundary segments";
  }
  catch (const std::exception& e) {
    t.check(false, "no exception") << e.what();
  }
  return t;
}


int main (int argc, char** argv)
{
  MPIHelper::instance(argc, argv);
  TestSuite t;

#if GMSH_UGGRID
  using Grid = UGGrid<2>;
#elif GMSH_ALBERTAGRID
  using Grid = AlbertaGrid<2,2>;
#endif

  t.subTest(testInsertionIndexing<Grid>("partialboundary.msh"));
  t.subTest(testDocumentedIndexing<Grid>("partialboundary.msh"));
  t.subTest(testInteriorInterface<Grid>());

  // Gmsh version 4 files: element and boundary data are not read at all
  t.subTest(testInsertionIndexing<Grid>("partialboundary-v4.msh"));

  return t.exit();
}
