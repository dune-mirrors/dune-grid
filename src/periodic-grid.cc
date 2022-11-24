#include <config.h>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/grid/albertagrid.hh>
#include <dune/grid/albertagrid/albertareader.hh>
#include <dune/grid/io/file/gmshreader.hh>

int main(int argc, char** argv)
{
  using namespace Dune;
  MPIHelper::instance(argc, argv);

  using Grid = AlbertaGrid<2,2>;
  using Factory = GridFactory<Grid>;

  std::string filename = "periodic";
  if (argc > 1)
    filename = std::string(argv[1]);

  {
    // use AlbertaReader
    Factory factory;
    factory.insertFaceTransformation({{1.0,0.0},{0.0,1.0}}, {1.0,0.0});
    AlbertaReader<Grid> reader;
    reader.readGrid(filename + ".2d", factory);
    auto gridPtr = factory.createGrid();
  }

  {
    // use GmshReader
    Factory factory;
    factory.insertFaceTransformation({{1.0,0.0},{0.0,1.0}}, {1.0,0.0});
    GmshReader<Grid>::read(factory, filename + ".msh");
    auto gridPtr = factory.createGrid();
  }

}