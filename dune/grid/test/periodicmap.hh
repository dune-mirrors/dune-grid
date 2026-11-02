// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_TEST_PERIODICMAP_HH
#define DUNE_GRID_TEST_PERIODICMAP_HH

#include <algorithm>
#include <array>
#include <bitset>
#include <cmath>
#include <cstddef>
#include <functional>
#include <map>
#include <utility>
#include <vector>

#include <dune/common/float_cmp.hh>
#include <dune/common/ftraits.hh>
#include <dune/common/hybridutilities.hh>
#include <dune/common/reservedvector.hh>

#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/type.hh>
#include <dune/geometry/typeindex.hh>

#include <dune/grid/common/capabilities.hh>
#include <dune/grid/common/gridenums.hh>
#include <dune/grid/common/rangegenerators.hh>
#include <dune/grid/test/coordinatetolerance.hh>

/** \file
 *  \brief Description of periodic grids for the grid checks
 *
 *  Some grids store several copies of the same entity, e.g. a periodic
 *  YaspGrid holds copies of the entities at the opposite side of the domain
 *  in its overlap. This also happens on a single process: the "neighbor"
 *  across the periodic boundary is the process itself. Such copies
 *  - have different geometries (they are translated by a period),
 *  - have different indices (each copy is a separate entity of the grid view),
 *  - but have the same (local and global) id, as they represent the same entity.
 *
 *  Many grid checks compare ids with entity positions and therefore need to
 *  know which positions are identified. The test describes this by a
 *  `PeriodicMap`: a function that maps each global coordinate to the unique
 *  representative of its point in the periodic domain. The requirements on
 *  such a map are listed in the documentation of `PeriodicMap`. The checks
 *  derive everything else from it: `PeriodicIdentification` groups the
 *  entities of one codimension of a grid view into classes of copies of the
 *  same entity.
 *
 *  Two standard implementations are provided:
 *  - `noPeriodicMap()`: the identity, i.e. the description of a non-periodic
 *    grid. It is the default argument of all checks, which then behave as
 *    without periodic information.
 *  - `periodicBoxMap()`: an axis-aligned box that is periodic in selected
 *    directions, e.g. the domain of a periodic YaspGrid.
 *
 *  Usage, e.g. for a YaspGrid on [lower,upper] that is periodic in direction 0:
 *  \code
 *  std::bitset<dim> periodic;
 *  periodic[0] = true;
 *  const Dune::GridCheck::PeriodicMap<Grid> periodicMap
 *    = Dune::GridCheck::periodicBoxMap(lower, upper, periodic);
 *  gridcheck(grid, periodicMap);
 *  checkIterators(grid.leafGridView(), periodicMap);
 *  checkPartitionType(grid.leafGridView(), periodicMap);
 *  \endcode
 *
 *  The map is purely geometric and the same on all processes, but what a
 *  process can conclude from it differs:
 *  - Copies stored on the same process ("local coincidences") occur on a
 *    single process, or in parallel if a process is its own neighbor across a
 *    periodic boundary. They are visible to the local checks (ids, iterators,
 *    partition types).
 *  - In parallel, the original of a copy in the overlap usually is stored on
 *    another process. The map moves the copy, but there is no local
 *    coincidence. Such identifications can only be verified by communication
 *    (see check_communication_correctness).
 *  Hence, the local checks use the classes of `PeriodicIdentification` and do
 *  not conclude anything from an entity being moved by the map alone.
 */

namespace Dune {
  namespace GridCheck {

    /** \brief Map a global coordinate to its unique representative in the periodic unit cell
     *
     *  A periodic map describes which global coordinates denote the same point
     *  of a periodic domain: two coordinates x and y denote the same point if
     *  and only if map(x) and map(y) coincide. The grid checks rely on the
     *  following requirements:
     *  - Unique representatives: the map is constant on the coordinates of one
     *    point and maps different points to different representatives.
     *    Coincidence is understood up to the tolerance of coordinateTolerance().
     *  - Idempotence: map(map(x)) == map(x), i.e. a representative is mapped to
     *    itself.
     *  - Defined for every coordinate the grid stores, in particular for the
     *    overlap outside of the periodic unit cell.
     *  - Pointwise and purely geometric: the map is applied to the centers and
     *    to each corner of the entities separately (see EntityPosition). Hence,
     *    the copies of an entity must be mapped point by point onto the same
     *    points. The map must not depend on entities, indices, or the process;
     *    it is the same function on all processes.
     *  - Robust to round-off: coordinates of the same point that differ by
     *    round-off must be mapped to representatives that coincide up to the
     *    tolerance. In particular, points on the "seam" of the unit cell must
     *    be mapped consistently, e.g. a point on the upper boundary onto the
     *    lower boundary (periodicBoxMap() snaps such points).
     *
     *  The map does not need to be a translation or continuous; any
     *  identification of points that satisfies these requirements can be
     *  described. The identity describes a non-periodic grid.
     *
     *  Use in the grid checks:
     *  - checkIdSet() (via gridcheck()): entities with the same id have the same
     *    mapped center; periodic copies on one level have the same id.
     *  - checkIterators() and checkIndexSet() (via gridcheck()): each id is
     *    found once for each periodic copy stored on this process.
     *  - checkPartitionType(): among the copies on one process, at most one is
     *    an interior entity.
     *  - checkCommunication(), check_communication_correctness(): a copy
     *    receives the data of its original; communicating entities have the
     *    same global id.
     *  The classes of copies are built by PeriodicIdentification.
     *
     *  \tparam GridOrGridView a grid or a grid view; only the global coordinate type is used
     *
     *  \see noPeriodicMap(), periodicBoxMap(), PeriodicIdentification
     */
    template<class GridOrGridView>
    using PeriodicMap = std::function<
      typename GridOrGridView::template Codim<0>::Geometry::GlobalCoordinate(
        typename GridOrGridView::template Codim<0>::Geometry::GlobalCoordinate)>;

    /** \brief The identity map, i.e. the description of a non-periodic grid
     *
     *  This is the default argument of all grid checks that take a periodic
     *  map; with it, the checks behave as without periodic information.
     *
     *  \see PeriodicMap
     */
    template<class GridOrGridView>
    PeriodicMap<GridOrGridView> noPeriodicMap ()
    {
      return [](auto x) { return x; };
    }

    /** \brief Container of points with associated values, searched up to an absolute tolerance
     *
     *  Periodic maps subtract periods, so representatives of the same point
     *  differ by round-off. Hence, points must be compared up to an absolute
     *  tolerance, see coordinateTolerance(). Approximate equality is not
     *  transitive, and no ordering derived from it (e.g. a lexicographic order
     *  that treats approximately equal components as equal) is a strict weak
     *  ordering, not even for well-separated points. Therefore, it cannot be
     *  used as comparator of std::map or std::sort.
     *
     *  Instead, the points are sorted into buckets: a point x belongs to the
     *  cell q of a regular grid of width w = 2*tolerance, with
     *  q_i = floor(x_i / w). Cells are integer tuples, which have an exact
     *  ordering. A stored point that equals a query point up to the tolerance
     *  (in each component) lies in the cell of the query point or in a
     *  neighboring cell. A query therefore searches the 3^n cells around the
     *  cell of the query point, n the number of components, and compares the
     *  candidates with coordinatesEqual() and an additional predicate only.
     *  The cost of a query is O(3^n times the number of points per cell).
     *
     *  \tparam Coordinate  coordinate type (a FieldVector)
     *  \tparam Value       value stored with each point
     */
    template<class Coordinate, class Value>
    class ApproximatePointMap
    {
      static constexpr int n = Coordinate::dimension;
      using Cell = std::array<long long, n>;

    public:
      using Field = typename FieldTraits<Coordinate>::field_type;

      /** \brief Construct an empty container
       *
       *  \param tolerance  absolute tolerance, see coordinateTolerance()
       */
      explicit ApproximatePointMap (Field tolerance)
        : tolerance_(tolerance), width_(2*tolerance)
      {
        // any positive width finds all points within the tolerance; a zero
        // tolerance means an exact comparison
        if (!(width_ > 0))
          width_ = 1;
      }

      /** \brief Find a stored point that equals x up to the tolerance and whose value satisfies pred
       *
       *  \param x     the query point
       *  \param pred  predicate on the stored values, e.g. a further comparison
       *  \returns a pointer to the value of such a point, or nullptr if there is
       *           none; it is invalidated by the next call of insert()
       */
      template<class Predicate>
      Value* find (const Coordinate& x, Predicate&& pred)
      {
        const Cell cell = cellOf(x);
        int numNeighbors = 1;
        for (int i = 0; i < n; ++i)
          numNeighbors *= 3;
        for (int k = 0; k < numNeighbors; ++k)
        {
          // the k-th neighbor cell, offsets -1, 0, +1 encoded in the digits of k in base 3
          Cell neighbor = cell;
          for (int i = 0, m = k; i < n; ++i, m /= 3)
            neighbor[i] += m % 3 - 1;
          const auto bucket = buckets_.find(neighbor);
          if (bucket == buckets_.end())
            continue;
          for (auto& [point, value] : bucket->second)
            if (coordinatesEqual(point, x, tolerance_) && pred(std::as_const(value)))
              return &value;
        }
        return nullptr;
      }

      /** \brief Insert the point x with the given value (equal points are not merged)
       *
       *  \param x      the point
       *  \param value  the value stored with x
       */
      void insert (const Coordinate& x, Value value)
      {
        buckets_[cellOf(x)].emplace_back(x, std::move(value));
      }

    private:
      Cell cellOf (const Coordinate& x) const
      {
        using std::floor;
        Cell cell;
        for (int i = 0; i < n; ++i)
          cell[i] = static_cast<long long>(floor(x[i] / width_));
        return cell;
      }

      Field tolerance_;
      Field width_;
      std::map<Cell, std::vector<std::pair<Coordinate, Value>>> buckets_;
    };

    /** \brief Whether the map moves the point x, i.e. whether x is a periodic copy
     *
     *  \param periodicMap  the periodic map, see PeriodicMap
     *  \param x            the point
     *  \param tolerance    absolute tolerance, see coordinateTolerance()
     */
    template<class Map, class Coordinate>
    bool isPeriodicCopy (const Map& periodicMap, const Coordinate& x,
                         typename FloatCmp::EpsilonType<Coordinate>::Type tolerance)
    {
      return !coordinatesEqual(periodicMap(x), x, tolerance);
    }

    /** \brief Periodic map for an axis-aligned box [lower,upper)
     *
     *  In each direction with `periodic[i] == true` the coordinate is mapped
     *  into the half-open interval [lower[i],upper[i]) by subtracting a
     *  multiple of the length upper[i]-lower[i]; the other directions are not
     *  changed. Points on the upper boundary, up to a tolerance of
     *  relativeCoordinateTolerance() (in the field type of Coordinate) times
     *  the length of the box, are mapped to the lower boundary, so that both
     *  sides of the seam have the same representative (see also
     *  coordinateTolerance()).
     *
     *  The returned std::function converts to PeriodicMap<Grid> if Coordinate
     *  is the global coordinate type of the grid. For an example see
     *  yaspPeriodicMap() in dune/grid/test/yasp/test-yaspgrid.hh.
     *
     *  \param lower     lower left corner of the box
     *  \param upper     upper right corner of the box
     *  \param periodic  the directions in which the box is periodic
     *
     *  \see PeriodicMap
     */
    template<class Coordinate, std::size_t dim>
    std::function<Coordinate(Coordinate)> periodicBoxMap (const Coordinate& lower, const Coordinate& upper,
                                                          std::bitset<dim> periodic)
    {
      using Field = typename FieldTraits<Coordinate>::field_type;
      const Field relTol = relativeCoordinateTolerance<Field>();
      return [=](Coordinate x) {
        using std::floor;
        for (std::size_t i = 0; i < dim; ++i)
        {
          if (!periodic[i])
            continue;
          const auto length = upper[i] - lower[i];
          auto t = (x[i] - lower[i]) / length;
          t -= floor(t);
          x[i] = lower[i] + t*length;
          // snap round-off at the upper boundary to the lower boundary
          if (coordinatesEqual(x[i], upper[i], relTol*length))
            x[i] = lower[i];
        }
        return x;
      };
    }

    /** \brief Position of an entity up to the periodic map: geometry type, mapped center, and mapped corners
     *
     *  Periodic copies of an entity have the same position. Neither the center
     *  nor the corners alone identify an entity:
     *  - different entities can share their center, e.g. a coarse face and the
     *    middle child of the refined face of its neighbor in a nonconforming
     *    leaf view,
     *  - different entities can share their set of mapped corners, e.g. the two
     *    cells of a grid with two cells in a periodic direction.
     *
     *  The corners are compared as sets (see equals()), so that the position
     *  does not depend on the numbering of the corners of the element the
     *  entity is taken from.
     *
     *  \tparam Coordinate  global coordinate type
     *  \tparam dim         dimension of the grid (bounds the number of corners by 2^dim)
     */
    template<class Coordinate, int dim>
    class EntityPosition
    {
    public:
      using Corners = ReservedVector<Coordinate, (1 << dim)>;
      using Field = typename FieldTraits<Coordinate>::field_type;

      /** \brief construct from the mapped center and the mapped corners
       *
       *  \param type     global geometry type index (GlobalGeometryTypeIndex)
       *  \param center   the mapped center
       *  \param corners  the mapped corners, in any order
       */
      EntityPosition (std::size_t type, const Coordinate& center, const Corners& corners)
        : type_(type), center_(center), corners_(corners)
      {}

      std::size_t type () const { return type_; }
      const Coordinate& center () const { return center_; }
      //! the mapped corners, in the order given to the constructor
      const Corners& corners () const { return corners_; }

      /** \brief Whether two positions are equal up to the tolerance
       *
       *  Positions are equal if they have the same type, their centers are
       *  equal, and their sets of corners are equal, i.e. they have the same
       *  number of corners and each corner of one position equals some corner
       *  of the other one (checked in both directions). Coordinates are
       *  compared with coordinatesEqual().
       *
       *  \param other      the position to compare with
       *  \param tolerance  absolute tolerance, see coordinateTolerance()
       */
      bool equals (const EntityPosition& other, Field tolerance) const
      {
        if (type_ != other.type_ || corners_.size() != other.corners_.size()
            || !coordinatesEqual(center_, other.center_, tolerance))
          return false;
        return containsCorners(other, tolerance) && other.containsCorners(*this, tolerance);
      }

    private:
      // whether each corner of other equals some corner of this position
      bool containsCorners (const EntityPosition& other, Field tolerance) const
      {
        return std::all_of(other.corners_.begin(), other.corners_.end(), [&](const Coordinate& x) {
          return std::any_of(corners_.begin(), corners_.end(), [&](const Coordinate& y) {
            return coordinatesEqual(x, y, tolerance);
          });
        });
      }

      std::size_t type_;
      Coordinate center_;
      Corners corners_;
    };

    /** \brief Group the entities of one codimension of a grid view into classes of local periodic copies
     *
     *  Entities are identified by their geometry type and their index in the
     *  index set of the grid view (this works for all grid views, a mapper
     *  would require a copy-assignable grid view). Two entities belong to the
     *  same class if they have the same EntityPosition, i.e., the same geometry
     *  type, and their centers and their corners are mapped to the same points
     *  by the periodic map. For the identity map, each class consists of a
     *  single entity unless the grid view contains two entities with the same
     *  geometry.
     *
     *  Only entities stored on this process are considered. A class with more
     *  than one member is a local coincidence. In parallel, the original of a
     *  copy is usually stored on another process, then the class of the copy
     *  has a single member.
     *
     *  Ghost entities are not identified: each ghost entity is a class of its
     *  own. Ghosts are virtual copies of entities of other processes with less
     *  consistency requirements; e.g., ALUGrid may store the same remote
     *  element (and its faces) several times as ghost, with different indices
     *  and ids. The grid interface does not forbid this.
     */
    template<class GridView>
    class PeriodicIdentification
    {
    public:
      using Index = typename GridView::IndexSet::IndexType;
      //! the global geometry type index and the index of an entity
      //! (not the GeometryType itself: GeometryTypes that compare equal can still
      //! differ in their topology id, and std::pair compares them via the id)
      using Key = std::pair<std::size_t, Index>;
      using Coordinate = typename GridView::template Codim<0>::Geometry::GlobalCoordinate;

      /** \brief Build the classes of the entities of one codimension
       *
       *  \param gridView     the grid view; must outlive this object (its index set is stored)
       *  \param codim        the codimension of the entities to group
       *  \param periodicMap  the periodic map, see PeriodicMap
       */
      PeriodicIdentification (const GridView& gridView, int codim, const PeriodicMap<GridView>& periodicMap)
        : indexSet_(gridView.indexSet())
      {
        constexpr int dim = GridView::dimension;
        const auto tolerance = coordinateTolerance(gridView);

        // group the entities by their position up to the periodic map,
        // searching the groups by the mapped center
        using Position = EntityPosition<Coordinate, dim>;
        std::vector<std::vector<Key>> groups;
        std::vector<Position> groupPositions;
        ApproximatePointMap<Coordinate, std::size_t> groupIds(tolerance);
        std::map<Key, bool> seen;

        for (const auto& element : elements(gridView))
        {
          const auto& geometry = element.geometry();
          const auto refElement = referenceElement<typename GridView::ctype, dim>(element.type());
          for (int i = 0; i < refElement.size(codim); ++i)
          {
            const Key key = subKey(element, i, codim);
            if (!seen.emplace(key, true).second)
              continue;

            // ghost entities are not identified, see the class documentation
            if (subPartitionType(element, i, codim) == GhostEntity)
            {
              addClass({ key });
              continue;
            }

            typename Position::Corners corners;
            for (int j = 0; j < refElement.size(i, codim, dim); ++j)
              corners.push_back(periodicMap(geometry.corner(refElement.subEntity(i, codim, j, dim))));
            const Position position(key.first, periodicMap(geometry.global(refElement.position(i, codim))),
                                    corners);

            const std::size_t* groupId = groupIds.find(position.center(), [&](std::size_t g) {
              return groupPositions[g].equals(position, tolerance);
            });
            if (groupId)
              groups[*groupId].push_back(key);
            else
            {
              groupIds.insert(position.center(), groups.size());
              groupPositions.push_back(position);
              groups.push_back({ key });
            }
          }
        }

        for (const auto& group : groups)
          addClass(group);
      }

      /** \brief the key (geometry type and index) of an entity
       *
       *  \param entity  an entity of the grid view
       */
      template<class Entity>
      Key key (const Entity& entity) const
      {
        return Key(GlobalGeometryTypeIndex::index(entity.type()), indexSet_.index(entity));
      }

      /** \brief the key (geometry type and index) of a subentity of an element
       *
       *  \param element  an element of the grid view
       *  \param i        number of the subentity within the element
       *  \param codim    codimension of the subentity
       */
      template<class Element>
      Key subKey (const Element& element, int i, int codim) const
      {
        const auto type = referenceElement<typename GridView::ctype, GridView::dimension>(element.type()).type(i, codim);
        return Key(GlobalGeometryTypeIndex::index(type), indexSet_.subIndex(element, i, codim));
      }

      //! the number of classes
      std::size_t size () const { return members_.size(); }

      /** \brief the class of the entity with the given key
       *
       *  \param key  the key of an entity of the grid view and codimension, see key() and subKey()
       */
      std::size_t classOf (const Key& key) const { return classOf_.at(key); }

      /** \brief the keys of all entities in a class
       *
       *  \param classId  the class, in [0, size())
       */
      const std::vector<Key>& members (std::size_t classId) const { return members_[classId]; }

      //! whether some class has more than one member, i.e. this process stores several copies of an entity
      bool hasCoincidences () const { return hasCoincidences_; }

    private:
      // the partition type of a subentity of an element, for a codim given at run time
      template<class Element>
      static PartitionType subPartitionType (const Element& element, int i, int codim)
      {
        PartitionType partitionType = InteriorEntity;
        Hybrid::switchCases(std::make_index_sequence<GridView::dimension+1>{}, codim, [&](auto cd) {
          if constexpr (Capabilities::hasEntity<typename GridView::Grid, cd>::v)
            partitionType = element.template subEntity<cd>(i).partitionType();
        });
        return partitionType;
      }

      void addClass (const std::vector<Key>& keys)
      {
        for (const auto& key : keys)
          classOf_[key] = members_.size();
        members_.push_back(keys);
        hasCoincidences_ = hasCoincidences_ || (keys.size() > 1);
      }

      const typename GridView::IndexSet& indexSet_;
      std::map<Key, std::size_t> classOf_;
      std::vector<std::vector<Key>> members_;
      bool hasCoincidences_ = false;
    };

  } // end namespace GridCheck
} // end namespace Dune

#endif // DUNE_GRID_TEST_PERIODICMAP_HH
