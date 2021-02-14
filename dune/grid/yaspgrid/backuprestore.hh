// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_YASPGRID_BACKUPRESTORE_HH
#define DUNE_GRID_YASPGRID_BACKUPRESTORE_HH

//- system headers
#include <fstream>

//- Dune headers
#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>
#include <dune/grid/common/backuprestore.hh>
#include <dune/grid/yaspgrid.hh>

// bump this version number up if you introduce any changes
// to the outout format of the YaspGrid BackupRestoreFacility.
#define YASPGRID_BACKUPRESTORE_FORMAT_VERSION 2

namespace Dune
{

  template<class Coordinates>
  struct MaybeHaveOrigin
  {

    template<class S>
    static void writeOrigin(S& s, const Coordinates& coord)
    {}

    template<class S>
    static void readOrigin(S& s, Dune::FieldVector<typename Coordinates::ctype,Coordinates::dimension>& coord)
    {}

    template<typename... A>
    static typename Dune::YaspGrid<Coordinates::dimension,Coordinates>* createGrid(
      const Dune::FieldVector<typename Coordinates::ctype,Coordinates::dimension>& lowerleft, A... args)
    {
      return new Dune::YaspGrid<Coordinates::dimension,Coordinates>(args...);
    }
  };

  template<class ctype, int dim>
  struct MaybeHaveOrigin<Dune::EquidistantOffsetCoordinates<ctype, dim> >
  {
    typedef typename Dune::EquidistantOffsetCoordinates<ctype, dim> Coordinates;

    template<class S>
    static void writeOrigin(S& s, const Coordinates& coord)
    {
      s << "Origin: ";
      for (int i=0; i<dim; i++)
        s << coord.origin(i) << " ";
      s << std::endl;
    }

    template<class S>
    static void readOrigin(S& s, Dune::FieldVector<ctype, dim>& coord)
    {
      std::string token;
      s >> token;
      for (int i=0; i<dim; i++)
        s >> coord[i];
    }

    template<typename... A>
    static typename Dune::YaspGrid<Coordinates::dimension,Coordinates>* createGrid(
      const Dune::FieldVector<typename Coordinates::ctype,Coordinates::dimension>& lowerleft,
      const Dune::FieldVector<typename Coordinates::ctype,Coordinates::dimension>& extension, A... args)
    {
      Dune::FieldVector<typename Coordinates::ctype,Coordinates::dimension> upperright(lowerleft);
      upperright += extension;
      return new Dune::YaspGrid<Coordinates::dimension,Coordinates>(lowerleft, upperright, args...);
    }
  };

  /** \copydoc Dune::BackupRestoreFacility */
  template<int dim, class Coordinates>
  struct BackupRestoreFacility<Dune::YaspGrid<dim, Coordinates> >
  {
    // type of grid
    typedef typename Dune::YaspGrid<dim, Coordinates> Grid;
    typedef typename Grid::ctype ctype;
    typedef typename Grid::Traits::CollectiveCommunication Comm;

    /** \copydoc Dune::BackupRestoreFacility::backup(grid,filename)  */
    static void backup ( const Grid &grid, const std::string &filename );

    /** \copydoc Dune::BackupRestoreFacility::backup(grid,stream)  */
    static void backup ( const Grid &grid, std::ostream &stream );

    /** \copydoc Dune::BackupRestoreFacility::restore(filename) */
    static Grid *restore (const std::string &filename, Comm comm = Comm());

    /** \copydoc Dune::BackupRestoreFacility::restore(stream) */
    static Grid *restore (std::istream &stream, Comm comm = Comm());
  };

  extern template struct BackupRestoreFacility<YaspGrid<2,EquidistantCoordinates<double,2> > >;
  extern template struct BackupRestoreFacility<YaspGrid<3,EquidistantCoordinates<double,3> > >;
  extern template struct BackupRestoreFacility<YaspGrid<2,EquidistantOffsetCoordinates<double,2> > >;
  extern template struct BackupRestoreFacility<YaspGrid<3,EquidistantOffsetCoordinates<double,3> > >;


  /** \copydoc Dune::BackupRestoreFacility */
  template<int dim, class ctype>
  struct BackupRestoreFacility<YaspGrid<dim,TensorProductCoordinates<ctype,dim> > >
  {
    // type of grid
    typedef YaspGrid<dim,TensorProductCoordinates<ctype,dim> > Grid;
    typedef typename Grid::Traits::CollectiveCommunication Comm;

    /** \copydoc Dune::BackupRestoreFacility::backup(grid,filename)  */
    static void backup ( const Grid &grid, const std::string &filename );

    /** \copydoc Dune::BackupRestoreFacility::backup(grid,stream)  */
    static void backup ( const Grid &grid, std::ostream &stream );

    /** \copydoc Dune::BackupRestoreFacility::restore(filename) */
    static Grid *restore (const std::string &filename, Comm comm = Comm());

    /** \copydoc Dune::BackupRestoreFacility::restore(stream) */
    static Grid *restore (std::istream &stream, Comm comm = Comm());
  };

  extern template struct BackupRestoreFacility<YaspGrid<2,TensorProductCoordinates<double,2> > >;
  extern template struct BackupRestoreFacility<YaspGrid<3,TensorProductCoordinates<double,3> > >;

} // namespace Dune

#include <dune/grid/yaspgrid/backuprestore.impl.hh>

#endif // #ifndef DUNE_GRID_YASPGRID_BACKUPRESTORE_HH
