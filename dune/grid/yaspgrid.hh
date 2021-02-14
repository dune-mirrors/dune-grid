// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_YASPGRID_HH
#define DUNE_GRID_YASPGRID_HH

#include <iostream>
#include <vector>
#include <algorithm>
#include <stack>
#include <type_traits>

// either include stdint.h or provide fallback for uint8_t
#if HAVE_STDINT_H
#include <stdint.h>
#else
typedef unsigned char uint8_t;
#endif

#include <dune/grid/common/backuprestore.hh>
#include <dune/grid/common/grid.hh>     // the grid base classes
#include <dune/grid/common/capabilities.hh> // the capabilities
#include <dune/common/hybridutilities.hh>
#include <dune/common/power.hh>
#include <dune/common/bigunsignedint.hh>
#include <dune/common/typetraits.hh>
#include <dune/common/reservedvector.hh>
#include <dune/common/parallel/communication.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/geometry/axisalignedcubegeometry.hh>
#include <dune/geometry/type.hh>
#include <dune/grid/common/indexidset.hh>
#include <dune/grid/common/datahandleif.hh>


#if HAVE_MPI
#include <dune/common/parallel/mpicommunication.hh>
#endif

/*! \file yaspgrid.hh
 * YaspGrid stands for yet another structured parallel grid.
 * It will implement the dune grid interface for structured grids
 * with arbitrary overlap, parallel features with two overlap
 * models, periodic boundaries and a fast implementation allowing on-the-fly computations.
 */

namespace Dune {

  /* some sizes for building global ids
   */
  const int yaspgrid_dim_bits = 24; // bits for encoding each dimension
  const int yaspgrid_level_bits = 5; // bits for encoding level number


  //************************************************************************
  // forward declaration of templates

  template<int dim, class Coordinates>                             class YaspGrid;
  template<int mydim, int cdim, class GridImp>  class YaspGeometry;
  template<int codim, int dim, class GridImp>   class YaspEntity;
  template<int codim, class GridImp>            class YaspEntitySeed;
  template<int codim, PartitionIteratorType pitype, class GridImp> class YaspLevelIterator;
  template<class GridImp>            class YaspIntersectionIterator;
  template<class GridImp>            class YaspIntersection;
  template<class GridImp>            class YaspHierarchicIterator;
  template<class GridImp, bool isLeafIndexSet>                     class YaspIndexSet;
  template<class GridImp>            class YaspGlobalIdSet;
  template<class GridImp>            class YaspPersistentContainerIndex;

} // namespace Dune

#include <dune/grid/yaspgrid/coordinates.hh>
#include <dune/grid/yaspgrid/torus.hh>
#include <dune/grid/yaspgrid/ygrid.hh>
#include <dune/grid/yaspgrid/yaspgridgeometry.hh>
#include <dune/grid/yaspgrid/yaspgridentity.hh>
#include <dune/grid/yaspgrid/yaspgridintersection.hh>
#include <dune/grid/yaspgrid/yaspgridintersectioniterator.hh>
#include <dune/grid/yaspgrid/yaspgridhierarchiciterator.hh>
#include <dune/grid/yaspgrid/yaspgridentityseed.hh>
#include <dune/grid/yaspgrid/yaspgridleveliterator.hh>
#include <dune/grid/yaspgrid/yaspgridindexsets.hh>
#include <dune/grid/yaspgrid/yaspgrididset.hh>
#include <dune/grid/yaspgrid/yaspgridpersistentcontainer.hh>

namespace Dune {

  template<int dim, class Coordinates>
  struct YaspGridFamily
  {
#if HAVE_MPI
    typedef CollectiveCommunication<MPI_Comm> CCType;
#else
    typedef CollectiveCommunication<No_Comm> CCType;
#endif

    typedef GridTraits<dim,                                     // dimension of the grid
        dim,                                                    // dimension of the world space
        Dune::YaspGrid<dim, Coordinates>,
        YaspGeometry,YaspEntity,
        YaspLevelIterator,                                      // type used for the level iterator
        YaspIntersection,              // leaf  intersection
        YaspIntersection,              // level intersection
        YaspIntersectionIterator,              // leaf  intersection iter
        YaspIntersectionIterator,              // level intersection iter
        YaspHierarchicIterator,
        YaspLevelIterator,                                      // type used for the leaf(!) iterator
        YaspIndexSet< const YaspGrid< dim, Coordinates >, false >,                  // level index set
        YaspIndexSet< const YaspGrid< dim, Coordinates >, true >,                  // leaf index set
        YaspGlobalIdSet<const YaspGrid<dim, Coordinates> >,
        bigunsignedint<dim*yaspgrid_dim_bits+yaspgrid_level_bits+dim>,
        YaspGlobalIdSet<const YaspGrid<dim, Coordinates> >,
        bigunsignedint<dim*yaspgrid_dim_bits+yaspgrid_level_bits+dim>,
        CCType,
        DefaultLevelGridViewTraits, DefaultLeafGridViewTraits,
        YaspEntitySeed>
    Traits;
  };

#ifndef DOXYGEN
  template<int dim, int codim>
  struct YaspCommunicateMeta {
    template<class G, class DataHandle>
    static void comm (const G& g, DataHandle& data, InterfaceType iftype, CommunicationDirection dir, int level)
    {
      if (data.contains(dim,codim))
      {
        g.template communicateCodim<DataHandle,codim>(data,iftype,dir,level);
      }
      YaspCommunicateMeta<dim,codim-1>::comm(g,data,iftype,dir,level);
    }
  };

  template<int dim>
  struct YaspCommunicateMeta<dim,0> {
    template<class G, class DataHandle>
    static void comm (const G& g, DataHandle& data, InterfaceType iftype, CommunicationDirection dir, int level)
    {
      if (data.contains(dim,0))
        g.template communicateCodim<DataHandle,0>(data,iftype,dir,level);
    }
  };
#endif

  //************************************************************************
  /*!
   * \brief [<em> provides \ref Dune::Grid </em>]
   * \brief Provides a distributed structured cube mesh.
   * \ingroup GridImplementations
   *
   * YaspGrid stands for yet another structured parallel grid.
   * It implements the dune grid interface for structured grids
   * with arbitrary overlap (including zero),
   * periodic boundaries, and a fast implementation allowing on-the-fly computations.
   *
   * YaspGrid supports three coordinate modes: \ref EquidistantCoordinates,
   * \ref EquidistantOffsetCoordinates, and \ref Dune::TensorProductCoordinates.
   *
   * \tparam dim The dimension of the grid and its surrounding world
   * \tparam Coordinates The coordinate mode of the grid.
   */
  template<int dim, class Coordinates = EquidistantCoordinates<double, dim> >
  class YaspGrid
    : public GridDefaultImplementation<dim,dim,typename Coordinates::ctype,YaspGridFamily<dim, Coordinates> >
  {

    template<int, PartitionIteratorType, typename>
    friend class YaspLevelIterator;

    template<typename>
    friend class YaspHierarchicIterator;

  protected:

  public:
    //! Type used for coordinates
    typedef typename Coordinates::ctype ctype;
#if HAVE_MPI
    typedef CollectiveCommunication<MPI_Comm> CollectiveCommunicationType;
#else
    typedef CollectiveCommunication<No_Comm> CollectiveCommunicationType;
#endif

#ifndef DOXYGEN
    typedef typename Dune::YGrid<Coordinates> YGrid;
    typedef typename Dune::YGridList<Coordinates>::Intersection Intersection;

    /** \brief A single grid level within a YaspGrid
     */
    struct YGridLevel {

      /** \brief Level number of this level grid */
      int level() const
      {
        return level_;
      }

      Coordinates coords;

      std::array<YGrid, dim+1> overlapfront;
      std::array<YGridComponent<Coordinates>, StaticPower<2,dim>::power> overlapfront_data;
      std::array<YGrid, dim+1> overlap;
      std::array<YGridComponent<Coordinates>, StaticPower<2,dim>::power> overlap_data;
      std::array<YGrid, dim+1> interiorborder;
      std::array<YGridComponent<Coordinates>, StaticPower<2,dim>::power> interiorborder_data;
      std::array<YGrid, dim+1> interior;
      std::array<YGridComponent<Coordinates>, StaticPower<2,dim>::power> interior_data;

      std::array<YGridList<Coordinates>,dim+1> send_overlapfront_overlapfront;
      std::array<std::deque<Intersection>, StaticPower<2,dim>::power>  send_overlapfront_overlapfront_data;
      std::array<YGridList<Coordinates>,dim+1> recv_overlapfront_overlapfront;
      std::array<std::deque<Intersection>, StaticPower<2,dim>::power>  recv_overlapfront_overlapfront_data;

      std::array<YGridList<Coordinates>,dim+1> send_overlap_overlapfront;
      std::array<std::deque<Intersection>, StaticPower<2,dim>::power>  send_overlap_overlapfront_data;
      std::array<YGridList<Coordinates>,dim+1> recv_overlapfront_overlap;
      std::array<std::deque<Intersection>, StaticPower<2,dim>::power>  recv_overlapfront_overlap_data;

      std::array<YGridList<Coordinates>,dim+1> send_interiorborder_interiorborder;
      std::array<std::deque<Intersection>, StaticPower<2,dim>::power>  send_interiorborder_interiorborder_data;
      std::array<YGridList<Coordinates>,dim+1> recv_interiorborder_interiorborder;
      std::array<std::deque<Intersection>, StaticPower<2,dim>::power>  recv_interiorborder_interiorborder_data;

      std::array<YGridList<Coordinates>,dim+1> send_interiorborder_overlapfront;
      std::array<std::deque<Intersection>, StaticPower<2,dim>::power>  send_interiorborder_overlapfront_data;
      std::array<YGridList<Coordinates>,dim+1> recv_overlapfront_interiorborder;
      std::array<std::deque<Intersection>, StaticPower<2,dim>::power>  recv_overlapfront_interiorborder_data;

      // general
      YaspGrid<dim,Coordinates>* mg;  // each grid level knows its multigrid
      int overlapSize;           // in mesh cells on this level
      bool keepOverlap;

      /** \brief The level number within the YaspGrid level hierarchy */
      int level_;
    };

    //! define types used for arguments
    typedef std::array<int, dim> iTupel;
    typedef FieldVector<ctype, dim> fTupel;

    // communication tag used by multigrid
    enum { tag = 17 };
#endif

    //! return reference to torus
    const Torus<CollectiveCommunicationType, dim>& torus () const
    {
      return _torus;
    }

    //! return number of cells on finest level in given direction on all processors
    int globalSize(int i) const
    {
      return levelSize(maxLevel(),i);
    }

    //! return number of cells on finest level on all processors
    iTupel globalSize() const
    {
      return levelSize(maxLevel());
    }

    //! return size of the grid (in cells) on level l in direction i
    int levelSize(int l, int i) const
    {
      return _coarseSize[i] * (1 << l);
    }

    //! return size vector of the grid (in cells) on level l
    iTupel levelSize(int l) const
    {
      iTupel s;
      for (int i=0; i<dim; ++i)
        s[i] = levelSize(l,i);
      return s;
    }

    //! return whether the grid is periodic in direction i
    bool isPeriodic(int i) const
    {
      return _periodic[i];
    }

    bool getRefineOption() const
    {
      return keep_ovlp;
    }

    //! Iterator over the grid levels
    typedef typename ReservedVector<YGridLevel,32>::const_iterator YGridLevelIterator;

    //! return iterator pointing to coarsest level
    YGridLevelIterator begin () const
    {
      return YGridLevelIterator(_levels,0);
    }

    //! return iterator pointing to given level
    YGridLevelIterator begin (int i) const
    {
      if (i<0 || i>maxLevel())
        DUNE_THROW(GridError, "level not existing");
      return YGridLevelIterator(_levels,i);
    }

    //! return iterator pointing to one past the finest level
    YGridLevelIterator end () const
    {
      return YGridLevelIterator(_levels,maxLevel()+1);
    }

    // static method to create the default load balance strategy
    static const YLoadBalanceDefault<dim>* defaultLoadbalancer()
    {
      static YLoadBalanceDefault<dim> lb;
      return & lb;
    }

  protected:
    /** \brief Make a new YGridLevel structure
     *
     * \param coords      the coordinate container
     * \param periodic    indicate periodicity for each direction
     * \param o_interior  origin of interior (non-overlapping) cell decomposition
     * \param overlap     to be used on this grid level
     */
    void makelevel (const Coordinates& coords, std::bitset<dim> periodic, iTupel o_interior, int overlap);

#ifndef DOXYGEN
    /** \brief special data structure to communicate ygrids
     * Historically, this was needed because Ygrids had virtual functions and
     * a communicated virtual function table pointer introduced a bug. After the
     * change to tensorproductgrid, the dynamic polymorphism was removed, still this
     * is kept because it allows to communicate ygrids, that only have index, but no
     * coordinate information. This is sufficient, because all communicated YGrids are
     * intersected with a local grid, which has coordinate information.
     */
    struct mpifriendly_ygrid {
      mpifriendly_ygrid ()
      {
        std::fill(origin.begin(), origin.end(), 0);
        std::fill(size.begin(), size.end(), 0);
      }
      mpifriendly_ygrid (const YGridComponent<Coordinates>& grid)
        : origin(grid.origin()), size(grid.size())
      {}
      iTupel origin;
      iTupel size;
    };
#endif

    /** \brief Construct list of intersections with neighboring processors
     *
     * \param recvgrid the grid stored in this processor
     * \param sendgrid the subgrid to be sent to neighboring processors
     * \param sendlist the deque to fill with send intersections
     * \param recvlist the deque to fill with recv intersections
     * \returns two lists: Intersections to be sent and Intersections to be received
     */
    void intersections(const YGridComponent<Coordinates>& sendgrid, const YGridComponent<Coordinates>& recvgrid,
                        std::deque<Intersection>& sendlist, std::deque<Intersection>& recvlist);

  protected:

    typedef const YaspGrid<dim,Coordinates> GridImp;

    void init();

    void boundarysegmentssize();

  public:

    // define the persistent index type
    typedef bigunsignedint<dim*yaspgrid_dim_bits+yaspgrid_level_bits+dim> PersistentIndexType;

    //! the GridFamily of this grid
    typedef YaspGridFamily<dim, Coordinates> GridFamily;
    // the Traits
    typedef typename YaspGridFamily<dim, Coordinates>::Traits Traits;

    // need for friend declarations in entity
    typedef YaspIndexSet<YaspGrid<dim, Coordinates>, false > LevelIndexSetType;
    typedef YaspIndexSet<YaspGrid<dim, Coordinates>, true > LeafIndexSetType;
    typedef YaspGlobalIdSet<YaspGrid<dim, Coordinates> > GlobalIdSetType;

    /** Standard constructor for a YaspGrid with a given Coordinates object
     *  @param coordinates Object that stores or computes the vertex coordinates
     *  @param periodic tells if direction is periodic or not
     *  @param overlap size of overlap on coarsest grid (same in all directions)
     *  @param comm the collective communication object for this grid. An MPI communicator can be given here.
     *  @param lb pointer to an overloaded YLoadBalance instance
     */
    YaspGrid (const Coordinates& coordinates,
              std::bitset<dim> periodic = std::bitset<dim>(0ULL),
              int overlap = 1,
              CollectiveCommunicationType comm = CollectiveCommunicationType(),
              const YLoadBalance<dim>* lb = defaultLoadbalancer());

    /** Standard constructor for an equidistant YaspGrid
     *  @param L extension of the domain
     *  @param s number of cells on coarse mesh in each direction
     *  @param periodic tells if direction is periodic or not
     *  @param overlap size of overlap on coarsest grid (same in all directions)
     *  @param comm the collective communication object for this grid. An MPI communicator can be given here.
     *  @param lb pointer to an overloaded YLoadBalance instance
     */
    YaspGrid (Dune::FieldVector<ctype, dim> L,
              std::array<int, dim> s,
              std::bitset<dim> periodic = std::bitset<dim>(0ULL),
              int overlap = 1,
              CollectiveCommunicationType comm = CollectiveCommunicationType(),
              const YLoadBalance<dim>* lb = defaultLoadbalancer());

    /** Constructor for an equidistant YaspGrid with non-trivial origin
     *  @param lowerleft Lower left corner of the domain
     *  @param upperright Upper right corner of the domain
     *  @param s number of cells on coarse mesh in each direction
     *  @param periodic tells if direction is periodic or not
     *  @param overlap size of overlap on coarsest grid (same in all directions)
     *  @param comm the collective communication object for this grid. An MPI communicator can be given here.
     *  @param lb pointer to an overloaded YLoadBalance instance
     */
    YaspGrid (Dune::FieldVector<ctype, dim> lowerleft,
              Dune::FieldVector<ctype, dim> upperright,
              std::array<int, dim> s,
              std::bitset<dim> periodic = std::bitset<dim>(0ULL),
              int overlap = 1,
              CollectiveCommunicationType comm = CollectiveCommunicationType(),
              const YLoadBalance<dim>* lb = defaultLoadbalancer());

    /** @brief Standard constructor for a tensorproduct YaspGrid
     *  @param coords coordinate vectors to be used for coarse grid
     *  @param periodic tells if direction is periodic or not
     *  @param overlap size of overlap on coarsest grid (same in all directions)
     *  @param comm the collective communication object for this grid. An MPI communicator can be given here.
     *  @param lb pointer to an overloaded YLoadBalance instance
     */
    YaspGrid (std::array<std::vector<ctype>, dim> coords,
              std::bitset<dim> periodic = std::bitset<dim>(0ULL),
              int overlap = 1,
              CollectiveCommunicationType comm = CollectiveCommunicationType(),
              const YLoadBalance<dim>* lb = defaultLoadbalancer());

  private:

    /** @brief Constructor for a tensorproduct YaspGrid with only coordinate
     *         information on this processor
     *  @param comm MPI communicator where this mesh is distributed to
     *  @param coords coordinate vectors to be used for coarse grid
     *  @param periodic tells if direction is periodic or not
     *  @param overlap size of overlap on coarsest grid (same in all directions)
     *  @param coarseSize the coarse size of the global grid
     *  @param lb pointer to an overloaded YLoadBalance instance
     *
     *  @warning The construction of overlapping coordinate ranges is
     *           an error-prone procedure. For this reason, it is kept private.
     *           You can safely use it through BackupRestoreFacility. All other
     *           use is not supported for the moment.
     */
    YaspGrid (std::array<std::vector<ctype>, dim> coords,
              std::bitset<dim> periodic,
              int overlap,
              CollectiveCommunicationType comm,
              std::array<int,dim> coarseSize,
              const YLoadBalance<dim>* lb = defaultLoadbalancer());

    // the backup restore facility needs to be able to use above constructor
    friend struct BackupRestoreFacility<YaspGrid<dim,Coordinates> >;

    // do not copy this class
    YaspGrid(const YaspGrid&);

  public:

    /*! Return maximum level defined in this grid. Levels are numbered
          0 ... maxlevel with 0 the coarsest level.
     */
    int maxLevel() const
    {
      return _levels.size()-1;
    }

    //! refine the grid refCount times.
    void globalRefine (int refCount);

    /**
       \brief set options for refinement
       @param keepPhysicalOverlap [true] keep the physical size of the overlap, [false] keep the number of cells in the overlap.  Default is [true].
     */
    void refineOptions (bool keepPhysicalOverlap)
    {
      keep_ovlp = keepPhysicalOverlap;
    }

    /** \brief Marks an entity to be refined/coarsened in a subsequent adapt.

       \param[in] refCount Number of subdivisions that should be applied. Negative value means coarsening.
       \param[in] e        Entity to Entity that should be refined

       \return true if Entity was marked, false otherwise.

       \note
          -  On yaspgrid marking one element will mark all other elements of the level as well
          -  If refCount is lower than refCount of a previous mark-call, nothing is changed
     */
    bool mark( int refCount, const typename Traits::template Codim<0>::Entity & e );

    /** \brief returns adaptation mark for given entity

       \param[in] e   Entity for which adaptation mark should be determined

       \return int adaptation mark, here the default value 0 is returned
     */
    int getMark ( const typename Traits::template Codim<0>::Entity &e ) const
    {
      return ( e.level() == maxLevel() ) ? adaptRefCount : 0;
    }

    //! map adapt to global refine
    bool adapt ();

    //! returns true, if the grid will be coarsened
    bool preAdapt ();

    //! clean up some markers
    void postAdapt();

    //! one past the end on this level
    template<int cd, PartitionIteratorType pitype>
    typename Traits::template Codim<cd>::template Partition<pitype>::LevelIterator lbegin (int level) const
    {
      return levelbegin<cd,pitype>(level);
    }

    //! Iterator to one past the last entity of given codim on level for partition type
    template<int cd, PartitionIteratorType pitype>
    typename Traits::template Codim<cd>::template Partition<pitype>::LevelIterator lend (int level) const
    {
      return levelend<cd,pitype>(level);
    }

    //! version without second template parameter for convenience
    template<int cd>
    typename Traits::template Codim<cd>::template Partition<All_Partition>::LevelIterator lbegin (int level) const
    {
      return levelbegin<cd,All_Partition>(level);
    }

    //! version without second template parameter for convenience
    template<int cd>
    typename Traits::template Codim<cd>::template Partition<All_Partition>::LevelIterator lend (int level) const
    {
      return levelend<cd,All_Partition>(level);
    }

    //! return LeafIterator which points to the first entity in maxLevel
    template<int cd, PartitionIteratorType pitype>
    typename Traits::template Codim<cd>::template Partition<pitype>::LeafIterator leafbegin () const
    {
      return levelbegin<cd,pitype>(maxLevel());
    }

    //! return LeafIterator which points behind the last entity in maxLevel
    template<int cd, PartitionIteratorType pitype>
    typename Traits::template Codim<cd>::template Partition<pitype>::LeafIterator leafend () const
    {
      return levelend<cd,pitype>(maxLevel());
    }

    //! return LeafIterator which points to the first entity in maxLevel
    template<int cd>
    typename Traits::template Codim<cd>::template Partition<All_Partition>::LeafIterator leafbegin () const
    {
      return levelbegin<cd,All_Partition>(maxLevel());
    }

    //! return LeafIterator which points behind the last entity in maxLevel
    template<int cd>
    typename Traits::template Codim<cd>::template Partition<All_Partition>::LeafIterator leafend () const
    {
      return levelend<cd,All_Partition>(maxLevel());
    }

    // \brief obtain Entity from EntitySeed. */
    template <typename Seed>
    typename Traits::template Codim<Seed::codimension>::Entity
    entity(const Seed& seed) const;

    //! return size (= distance in graph) of overlap region
    int overlapSize (int level, int codim) const
    {
      YGridLevelIterator g = begin(level);
      return g->overlapSize;
    }

    //! return size (= distance in graph) of overlap region
    int overlapSize (int codim) const
    {
      YGridLevelIterator g = begin(maxLevel());
      return g->overlapSize;
    }

    //! return size (= distance in graph) of ghost region
    int ghostSize (int level, int codim) const
    {
      return 0;
    }

    //! return size (= distance in graph) of ghost region
    int ghostSize (int codim) const
    {
      return 0;
    }

    //! number of entities per level and codim in this process
    int size (int level, int codim) const;

    //! number of leaf entities per codim in this process
    int size (int codim) const
    {
      return size(maxLevel(),codim);
    }

    //! number of entities per level and geometry type in this process
    int size (int level, GeometryType type) const
    {
      return (type.isCube()) ? size(level,dim-type.dim()) : 0;
    }

    //! number of leaf entities per geometry type in this process
    int size (GeometryType type) const
    {
      return size(maxLevel(),type);
    }

    //! \brief returns the number of boundary segments within the macro grid
    size_t numBoundarySegments () const
    {
      return nBSegments;
    }

    //! \brief returns the size of the physical domain
    const Dune::FieldVector<ctype, dim>& domainSize () const {
      return _L;
    }

    /*! The new communication interface

       communicate objects for all codims on a given level
     */
    template<class DataHandleImp, class DataType>
    void communicate (CommDataHandleIF<DataHandleImp,DataType> & data, InterfaceType iftype, CommunicationDirection dir, int level) const
    {
      YaspCommunicateMeta<dim,dim>::comm(*this,data,iftype,dir,level);
    }

    /*! The new communication interface

       communicate objects for all codims on the leaf grid
     */
    template<class DataHandleImp, class DataType>
    void communicate (CommDataHandleIF<DataHandleImp,DataType> & data, InterfaceType iftype, CommunicationDirection dir) const
    {
      YaspCommunicateMeta<dim,dim>::comm(*this,data,iftype,dir,this->maxLevel());
    }

    /*! The new communication interface

       communicate objects for one codim
     */
    template<class DataHandle, int codim>
    void communicateCodim (DataHandle& data, InterfaceType iftype, CommunicationDirection dir, int level) const;

    // The new index sets from DDM 11.07.2005
    const typename Traits::GlobalIdSet& globalIdSet() const
    {
      return theglobalidset;
    }

    const typename Traits::LocalIdSet& localIdSet() const
    {
      return theglobalidset;
    }

    const typename Traits::LevelIndexSet& levelIndexSet(int level) const
    {
      if (level<0 || level>maxLevel()) DUNE_THROW(RangeError, "level out of range");
      return *(indexsets[level]);
    }

    const typename Traits::LeafIndexSet& leafIndexSet() const
    {
      return leafIndexSet_;
    }

    /*! @brief return a collective communication object
     */
    const CollectiveCommunicationType& comm () const
    {
      return ccobj;
    }

  private:

    // number of boundary segments of the level 0 grid
    int nBSegments;

    // Index classes need access to the real entity
    friend class Dune::YaspIndexSet<const Dune::YaspGrid<dim, Coordinates>, true >;
    friend class Dune::YaspIndexSet<const Dune::YaspGrid<dim, Coordinates>, false >;
    friend class Dune::YaspGlobalIdSet<const Dune::YaspGrid<dim, Coordinates> >;
    friend class Dune::YaspPersistentContainerIndex<const Dune::YaspGrid<dim, Coordinates> >;

    friend class Dune::YaspIntersectionIterator<const Dune::YaspGrid<dim, Coordinates> >;
    friend class Dune::YaspIntersection<const Dune::YaspGrid<dim, Coordinates> >;
    friend class Dune::YaspEntity<0, dim, const Dune::YaspGrid<dim, Coordinates> >;

    template<int codim_, int dim_, class GridImp_, template<int,int,class> class EntityImp_>
    friend class Entity;

    template<class DT>
    class MessageBuffer {
    public:
      // Constructor
      MessageBuffer (DT *p)
      {
        a=p;
        i=0;
        j=0;
      }

      // write data to message buffer, acts like a stream !
      template<class Y>
      void write (const Y& data)
      {
        static_assert(( std::is_same<DT,Y>::value ), "DataType mismatch");
        a[i++] = data;
      }

      // read data from message buffer, acts like a stream !
      template<class Y>
      void read (Y& data) const
      {
        static_assert(( std::is_same<DT,Y>::value ), "DataType mismatch");
        data = a[j++];
      }

    private:
      DT *a;
      int i;
      mutable int j;
    };

    //! one past the end on this level
    template<int cd, PartitionIteratorType pitype>
    YaspLevelIterator<cd,pitype,GridImp> levelbegin (int level) const;

    //! Iterator to one past the last entity of given codim on level for partition type
    template<int cd, PartitionIteratorType pitype>
    YaspLevelIterator<cd,pitype,GridImp> levelend (int level) const;

    CollectiveCommunicationType ccobj;

    Torus<CollectiveCommunicationType,dim> _torus;

    std::vector< std::shared_ptr< YaspIndexSet<const YaspGrid<dim,Coordinates>, false > > > indexsets;
    YaspIndexSet<const YaspGrid<dim,Coordinates>, true> leafIndexSet_;
    YaspGlobalIdSet<const YaspGrid<dim,Coordinates> > theglobalidset;

    Dune::FieldVector<ctype, dim> _L;
    iTupel _s;
    std::bitset<dim> _periodic;
    iTupel _coarseSize;
    ReservedVector<YGridLevel,32> _levels;
    int _overlap;
    bool keep_ovlp;
    int adaptRefCount;
    bool adaptActive;
  };

  //! Explicit template instantiation

  extern template class YaspGrid<2, EquidistantCoordinates<double, 2>>;
  extern template class YaspGrid<3, EquidistantCoordinates<double, 3>>;
  extern template class YaspGrid<2, EquidistantOffsetCoordinates<double, 2>>;
  extern template class YaspGrid<3, EquidistantOffsetCoordinates<double, 3>>;
  extern template class YaspGrid<2, TensorProductCoordinates<double, 2>>;
  extern template class YaspGrid<3, TensorProductCoordinates<double, 3>>;


  //! Output operator for multigrids

  template <int d, class CC>
  std::ostream& operator<< (std::ostream& s, const YaspGrid<d,CC>& grid);

  namespace Capabilities
  {

    /** \struct hasEntity
       \ingroup YaspGrid
     */

    /** \struct hasBackupRestoreFacilities
       \ingroup YaspGrid
     */
    template<int dim, class Coordinates>
    struct hasBackupRestoreFacilities< YaspGrid<dim, Coordinates> >
    {
      static const bool v = true;
    };

    /** \brief YaspGrid has only one geometry type for codim 0 entities
       \ingroup YaspGrid
     */
    template<int dim, class Coordinates>
    struct hasSingleGeometryType< YaspGrid<dim, Coordinates> >
    {
      static const bool v = true;
      static const unsigned int topologyId = GeometryTypes::cube(dim).id();
    };

    /** \brief YaspGrid is a Cartesian grid
        \ingroup YaspGrid
     */
    template<int dim, class Coordinates>
    struct isCartesian< YaspGrid<dim, Coordinates> >
    {
      static const bool v = true;
    };

    /** \brief YaspGrid has entities for all codimensions
       \ingroup YaspGrid
     */
    template<int dim, class Coordinates, int codim>
    struct hasEntity< YaspGrid<dim, Coordinates>, codim>
    {
      static const bool v = true;
    };

    /**
     * \brief YaspGrid can iterate over all codimensions
     * \ingroup YaspGrid
     **/
    template<int dim, class Coordinates, int codim>
    struct hasEntityIterator<YaspGrid<dim, Coordinates>, codim>
    {
      static const bool v = true;
    };

    /** \brief YaspGrid can communicate on all codimensions
     *  \ingroup YaspGrid
     */
    template<int dim, int codim, class Coordinates>
    struct canCommunicate< YaspGrid< dim, Coordinates>, codim >
    {
      static const bool v = true;
    };

    /** \brief YaspGrid is levelwise conforming
       \ingroup YaspGrid
     */
    template<int dim, class Coordinates>
    struct isLevelwiseConforming< YaspGrid<dim, Coordinates> >
    {
      static const bool v = true;
    };

    /** \brief YaspGrid is leafwise conforming
       \ingroup YaspGrid
     */
    template<int dim, class Coordinates>
    struct isLeafwiseConforming< YaspGrid<dim, Coordinates> >
    {
      static const bool v = true;
    };

  }

} // end namespace

// Include the implementation of the template class
#include <dune/grid/yaspgrid/yaspgrid.impl.hh>
// Include the specialization of the StructuredGridFactory class for YaspGrid
#include <dune/grid/yaspgrid/structuredyaspgridfactory.hh>
// Include the specialization of the BackupRestoreFacility class for YaspGrid
#include <dune/grid/yaspgrid/backuprestore.hh>

#endif
