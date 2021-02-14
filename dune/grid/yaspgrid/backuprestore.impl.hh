// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_YASPGRID_BACKUPRESTORE_IMPL_HH
#define DUNE_GRID_YASPGRID_BACKUPRESTORE_IMPL_HH

namespace Dune {

template<int dim, class Coordinates>
void BackupRestoreFacility<Dune::YaspGrid<dim, Coordinates> >::backup ( const Grid &grid, const std::string &filename )
{
  if (grid.comm().rank() == 0)
  {
    std::ofstream file(filename);
    if( file )
    {
      backup(grid,file);
      file.close();
    }
    else
      std::cerr << "ERROR: BackupRestoreFacility::backup: couldn't open file `" << filename << "'" << std::endl;
  }
}


template<int dim, class Coordinates>
void BackupRestoreFacility<Dune::YaspGrid<dim, Coordinates> >::backup ( const Grid &grid, std::ostream &stream )
{
  stream << "YaspGrid BackupRestore Format Version: " << YASPGRID_BACKUPRESTORE_FORMAT_VERSION << std::endl;
  stream << "Torus structure: ";
  for (int i=0; i<dim; i++)
    stream << grid.torus().dims(i) << " ";
  stream << std::endl << "Refinement level: " << grid.maxLevel() << std::endl;
  stream << "Periodicity: ";
  for (int i=0; i<dim; i++)
    stream << (grid.isPeriodic(i) ? "1 " : "0 ");
  stream << std::endl << "Overlap: " << grid.overlapSize(0,0) << std::endl;
  stream << "KeepPhysicalOverlap: ";
  for (typename Grid::YGridLevelIterator i=++grid.begin(); i != grid.end(); ++i)
    stream << (i->keepOverlap ? "1" : "0") << " ";
  stream << std::endl;
  stream << "Coarse Size: ";
  for (int i=0; i<dim; i++)
    stream << grid.levelSize(0,i) << " ";
  stream << std::endl;
  stream << "Meshsize: " ;
  for (int i=0; i<dim; i++)
    stream << grid.begin()->coords.meshsize(i,0) << " ";
  stream << std::endl;
  MaybeHaveOrigin<Coordinates>::writeOrigin(stream, grid.begin()->coords);
}


template<int dim, class Coordinates>
Dune::YaspGrid<dim, Coordinates>*
BackupRestoreFacility<Dune::YaspGrid<dim, Coordinates> >::restore (const std::string &filename, Comm comm)
{
  std::ifstream file(filename);
  if( file )
    return restore(file,comm);
  else
  {
    std::cerr << "ERROR: BackupRestoreFacility::restore: couldn't open file `" << filename << "'" << std::endl;
    return 0;
  }
}


template<int dim, class Coordinates>
Dune::YaspGrid<dim, Coordinates>*
BackupRestoreFacility<Dune::YaspGrid<dim, Coordinates> >::restore (std::istream &stream, Comm comm)
{
  std::string input;

  int version;
  stream >> input >> input >> input >> input;
  stream >> version;
  if (version != YASPGRID_BACKUPRESTORE_FORMAT_VERSION)
    DUNE_THROW(Dune::Exception, "Your YaspGrid backup file is written in an outdated format!");

  std::array<int,dim> torus_dims;
  stream >> input >> input;
  for (int i=0; i<dim; i++)
    stream >> torus_dims[i];

  int refinement;
  stream >> input >> input;
  stream >> refinement;

  std::bitset<dim> periodic;
  bool b;
  stream >> input;
  for (int i=0; i<dim; i++)
  {
    stream >> b;
    periodic[i] = b;
  }

  int overlap;
  stream >> input;
  stream >> overlap;

  std::vector<bool> physicalOverlapSize;
  physicalOverlapSize.resize(refinement);
  stream >> input;
  for (int i=0; i<refinement; ++i)
  {
    stream >> b;
    physicalOverlapSize[i] = b;
  }

  std::array<int,dim> coarseSize;
  stream >> input >> input;
  for (int i=0; i<dim; i++)
    stream >> coarseSize[i];

  Dune::FieldVector<ctype,dim> h;
  stream >>  input;
  for (int i=0; i<dim; i++)
    stream >> h[i];

  Dune::FieldVector<ctype,dim> origin;
  MaybeHaveOrigin<Coordinates>::readOrigin(stream, origin);

  // the constructor takes the upper right corner...
  Dune::FieldVector<ctype,dim> length(h);
  for (int i=0; i<dim; i++)
    length[i] *= coarseSize[i];

  YaspFixedSizePartitioner<dim> lb(torus_dims);

  Grid* grid = MaybeHaveOrigin<Coordinates>::createGrid(origin, length, coarseSize, periodic, overlap, comm, &lb);

  for (int i=0; i<refinement; ++i)
  {
    grid->refineOptions(physicalOverlapSize[i]);
    grid->globalRefine(1);
  }

  return grid;
}


template<int dim, class ctype>
void BackupRestoreFacility<YaspGrid<dim,TensorProductCoordinates<ctype,dim> > >::backup ( const Grid &grid, const std::string &filename )
{
  std::ostringstream filename_str;
  filename_str << filename << grid.comm().rank();
  std::ofstream file( filename_str.str() );
  if( file )
  {
    backup(grid,file);
    file.close();
  }
  else
    std::cerr << "ERROR: BackupRestoreFacility::backup: couldn't open file `" << filename_str.str() << "'" << std::endl;
}


template<int dim, class ctype>
void BackupRestoreFacility<YaspGrid<dim,TensorProductCoordinates<ctype,dim> > >::backup ( const Grid &grid, std::ostream &stream )
{
  stream << "YaspGrid BackupRestore Format Version: " << YASPGRID_BACKUPRESTORE_FORMAT_VERSION << std::endl;
  stream << "Torus structure: ";
  for (int i=0; i<dim; i++)
    stream << grid.torus().dims(i) << " ";
  stream << std::endl << "Refinement level: " << grid.maxLevel() << std::endl;
  stream << "Periodicity: ";
  for (int i=0; i<dim; i++)
    stream << (grid.isPeriodic(i) ? "1 " : "0 ");
  stream << std::endl << "Overlap: " << grid.overlapSize(0,0) << std::endl;
  stream << "KeepPhysicalOverlap: ";
  for (typename Grid::YGridLevelIterator i=++grid.begin(); i != grid.end(); ++i)
    stream << (i->keepOverlap ? "1" : "0") << " ";
  stream << std::endl;
  stream << "Coarse Size: ";
  for (int i=0; i<dim; i++)
    stream << grid.levelSize(0,i) << " ";
  stream << std::endl;

  grid.begin()->coords.print(stream);
}


template<int dim, class ctype>
YaspGrid<dim,TensorProductCoordinates<ctype,dim> >*
BackupRestoreFacility<YaspGrid<dim,TensorProductCoordinates<ctype,dim> > >::restore (const std::string &filename, Comm comm)
{
  std::ostringstream filename_str;
  filename_str << filename;
  filename_str << comm.rank();
  std::ifstream file(filename_str.str());
  if( file )
    return restore(file, comm);
  else
  {
    std::cerr << "ERROR: BackupRestoreFacility::restore: couldn't open file `" << filename_str.str() << "'" << std::endl;
    return 0;
  }
}


template<int dim, class ctype>
YaspGrid<dim,TensorProductCoordinates<ctype,dim> >*
BackupRestoreFacility<YaspGrid<dim,TensorProductCoordinates<ctype,dim> > >::restore (std::istream &stream, Comm comm)
{
  std::string input;

  int version;
  stream >> input >> input >> input >> input;
  stream >> version;
  if (version != YASPGRID_BACKUPRESTORE_FORMAT_VERSION)
    DUNE_THROW(Dune::Exception, "Your YaspGrid backup file is written in an outdated format!");

  std::array<int,dim> torus_dims;
  stream >> input >> input;
  for (int i=0; i<dim; i++)
    stream >> torus_dims[i];

  int refinement;
  stream >> input >> input;
  stream >> refinement;

  std::bitset<dim> periodic;
  bool b;
  stream >> input;
  for (int i=0; i<dim; i++)
  {
    stream >> b;
    periodic[i] = b;
  }

  int overlap;
  stream >> input;
  stream >> overlap;

  std::vector<bool> physicalOverlapSize;
  physicalOverlapSize.resize(refinement);
  stream >> input;
  for (int i=0; i<refinement; ++i)
  {
    stream >> b;
    physicalOverlapSize[i] = b;
  }


  std::array<int,dim> coarseSize;
  stream >> input >> input;
  for (int i=0; i<dim; i++)
    stream >> coarseSize[i];

  std::array<std::vector<ctype>,dim> coords;
  stream >> input >> input >> input >> input;
  for (int d=0; d<dim; d++)
  {
    stream >> input >> input;
    int size;
    stream >> size;
    stream >> input;
    ctype tmp;
    coords[d].resize(size);
    for (int i=0; i<size; i++)
    {
      stream >> tmp;
      coords[d][i] = tmp;
    }
  }

  YaspFixedSizePartitioner<dim> lb(torus_dims);
  Grid* grid = new Grid(coords, periodic, overlap, comm, coarseSize, &lb);

  for (int i=0; i<refinement; ++i)
  {
    grid->refineOptions(physicalOverlapSize[i]);
    grid->globalRefine(1);
  }

  return grid;
}

} // namespace Dune

#endif // #ifndef DUNE_GRID_YASPGRID_BACKUPRESTORE_IMPL_HH
