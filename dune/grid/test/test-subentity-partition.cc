// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

#include <dune-grid-config.hh>

#include <array>
#include <bitset>
#include <iostream>
#include <utility>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/test/testsuite.hh>
#include <dune/grid/common/entity.hh>
#include <dune/grid/geometrygrid.hh>
#include <dune/grid/identitygrid.hh>
#include <dune/grid/onedgrid.hh>
#include <dune/grid/yaspgrid.hh>
#if HAVE_DUNE_UGGRID
#include <dune/grid/uggrid.hh>
#include <dune/grid/utility/structuredgridfactory.hh>
#endif

// Exercise the actual default implementation even when the grid specializes it.
template<int codim, int dim, class Grid>
class DefaultPartitionEntity
  : public Dune::EntityDefaultImplementation<codim, dim, Grid, DefaultPartitionEntity>
{
  using Base = Dune::EntityDefaultImplementation<codim, dim, Grid, DefaultPartitionEntity>;
  using Entity = typename Grid::template Codim<0>::Entity;
  Entity entity_;

public:
  explicit DefaultPartitionEntity(const Entity& entity) : entity_(entity) {}

  using Base::partitionType;
  Dune::PartitionType partitionType() const { return entity_.partitionType(); }

  template<int cd>
  auto subEntity(int i) const { return entity_.template subEntity<cd>(i); }
};

template<class GridView>
void checkView(Dune::TestSuite& test, const GridView& view)
{
  using Grid = typename GridView::Grid;
  constexpr int dim = GridView::dimension;
  for (const auto& element : elements(view, Dune::Partitions::all))
  {
    const DefaultPartitionEntity<0, dim, Grid> fallback(element);
    test.check(element.partitionType(0, 0) == element.partitionType(), "element partition");
    test.check(fallback.partitionType() == element.partitionType(), "zero-argument overload");
    Dune::Hybrid::forEach(std::make_index_sequence<dim+1>{}, [&](auto cd) {
      for (unsigned int i = 0; i < element.subEntities(cd); ++i)
      {
        const auto expected = element.template subEntity<cd>(i).partitionType();
        test.check(element.partitionType(i, cd) == expected, "subentity partition")
          << "level=" << element.level() << " codim=" << cd << " i=" << i;
        test.check(fallback.partitionType(i, cd) == expected, "default subentity partition")
          << "level=" << element.level() << " codim=" << cd << " i=" << i;
      }
    });
    bool threw = false;
    try { fallback.partitionType(0, dim+1); }
    catch (const Dune::Exception&) { threw = true; }
    test.check(threw, "default rejects invalid codimension");
  }
}

template<class Grid>
void checkGrid(Dune::TestSuite& test, const Grid& grid)
{
  checkView(test, grid.leafGridView());
  for (int level = 0; level <= grid.maxLevel(); ++level)
    checkView(test, grid.levelGridView(level));
}

template<int dim>
struct IdentityCoordinates
  : Dune::AnalyticalCoordFunction<double, dim, dim, IdentityCoordinates<dim>>
{
  void evaluate(const Dune::FieldVector<double, dim>& x,
                Dune::FieldVector<double, dim>& y) const { y = x; }
};

template<int dim>
void checkYasp(Dune::TestSuite& test)
{
  using Grid = Dune::YaspGrid<dim>;
  std::array<int, dim> cells;
  cells.fill(8);
  for (bool periodic : {false, true})
    for (int overlap : {0, 1})
    {
      std::bitset<dim> periodic_directions;
      if (periodic)
        periodic_directions.set();
      Grid grid(Dune::FieldVector<double, dim>(1), cells, periodic_directions, overlap);
      grid.globalRefine(1);
      checkGrid(test, grid);
      Dune::IdentityGrid<Grid> identity(grid);
      checkGrid(test, identity);
      IdentityCoordinates<dim> coordinates;
      Dune::GeometryGrid<Grid, IdentityCoordinates<dim>> geometry(grid, coordinates);
      checkGrid(test, geometry);
    }
}

#if HAVE_DUNE_UGGRID
template<int dim>
void checkUG(Dune::TestSuite& test)
{
  using Grid = Dune::UGGrid<dim>;
  std::array<unsigned int, dim> cells;
  cells.fill(4);
  for (bool simplex : {false, true})
  {
    auto grid = simplex
      ? Dune::StructuredGridFactory<Grid>::createSimplexGrid(
          Dune::FieldVector<double, dim>(0), Dune::FieldVector<double, dim>(1), cells)
      : Dune::StructuredGridFactory<Grid>::createCubeGrid(
          Dune::FieldVector<double, dim>(0), Dune::FieldVector<double, dim>(1), cells);
    grid->loadBalance();
    grid->globalRefine(1);
    checkGrid(test, *grid);
  }
}
#endif

int main(int argc, char** argv)
{
  const auto& mpi = Dune::MPIHelper::instance(argc, argv);
  Dune::TestSuite test;
  try
  {
    checkYasp<1>(test);
    checkYasp<2>(test);
    checkYasp<3>(test);
    Dune::OneDGrid line(8, 0., 1.);
    line.globalRefine(1);
    checkGrid(test, line);
#if HAVE_DUNE_UGGRID
    checkUG<2>(test);
    checkUG<3>(test);
#endif
  }
  catch (const Dune::Exception& exception)
  {
    std::cerr << "Rank " << mpi.rank() << ": " << exception << '\n';
    return 1;
  }
  return test.exit();
}
