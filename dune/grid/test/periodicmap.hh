// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_TEST_PERIODICMAP_HH
#define DUNE_GRID_TEST_PERIODICMAP_HH

#include <bitset>
#include <cmath>
#include <cstddef>
#include <functional>
#include <map>
#include <utility>
#include <vector>

#include <dune/common/float_cmp.hh>

#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/type.hh>
#include <dune/geometry/typeindex.hh>

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

    /** \brief Lexicographic order on coordinates that treats approximately equal components as equal
     *
     *  Periodic maps subtract periods, so representatives of the same point
     *  differ by round-off. Components are compared with coordinatesEqual()
     *  and an absolute tolerance, see coordinateTolerance(). This is a strict
     *  weak ordering for point sets in which distinct points differ by much
     *  more than the tolerance in each differing component, which is the case
     *  for the entity centers of a mesh.
     */
    struct ApproxLess
    {
      double tolerance;

      template<class Coordinate>
      bool operator() (const Coordinate& a, const Coordinate& b) const
      {
        for (std::size_t i = 0; i < a.size(); ++i)
        {
          if (coordinatesEqual(a[i], b[i], tolerance))
            continue;
          return a[i] < b[i];
        }
        return false;
      }
    };

    /** \brief Whether the map moves the point x, i.e. whether x is a periodic copy
     *
     *  \param tolerance  absolute tolerance, see coordinateTolerance()
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
     *  changed. Points on the upper boundary, up to round-off of 1e-8 relative
     *  to the length of the box, are mapped to the lower boundary, so that both
     *  sides of the seam have the same representative.
     *
     *  The returned std::function converts to PeriodicMap<Grid> if Coordinate
     *  is the global coordinate type of the grid. For an example see
     *  yaspPeriodicMap() in dune/grid/test/yasp/test-yaspgrid.hh.
     *
     *  \see PeriodicMap
     */
    template<class Coordinate, std::size_t dim>
    std::function<Coordinate(Coordinate)> periodicBoxMap (const Coordinate& lower, const Coordinate& upper,
                                                          std::bitset<dim> periodic)
    {
      return [=](Coordinate x) {
        for (std::size_t i = 0; i < dim; ++i)
        {
          if (!periodic[i])
            continue;
          const auto length = upper[i] - lower[i];
          auto t = (x[i] - lower[i]) / length;
          t -= std::floor(t);
          // snap round-off at the upper boundary to the lower boundary
          if (1 - t <= 1e-8)
            t = 0;
          x[i] = lower[i] + t*length;
        }
        return x;
      };
    }

    /** \brief Group the entities of one codimension of a grid view into classes of local periodic copies
     *
     *  Entities are identified by their geometry type and their index in the
     *  index set of the grid view (this works for all grid views, a mapper
     *  would require a copy-assignable grid view). Two entities belong to the
     *  same class if they have the same geometry type, their centers are
     *  mapped to the same point by the periodic map, and at least one of the
     *  entities of this point is moved by the map, i.e. is a periodic copy.
     *  Entities that just share their center without periodicity (e.g. UGGrid
     *  stores different entities with the same center in its leaf view) are
     *  not identified. Hence, for the identity map each class consists of a
     *  single entity.
     *
     *  Only entities stored on this process are considered. A class with more
     *  than one member is a local coincidence. In parallel, the original of a
     *  copy is usually stored on another process, then the class of the copy
     *  has a single member.
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

      PeriodicIdentification (const GridView& gridView, int codim, const PeriodicMap<GridView>& periodicMap)
        : indexSet_(gridView.indexSet())
      {
        constexpr int dim = GridView::dimension;
        const auto tolerance = coordinateTolerance(gridView);

        // group the entities by geometry type and mapped center, and remember
        // whether a group contains a periodic copy
        struct Group { std::vector<Key> keys; bool hasCopy = false; };
        std::vector<Group> groups;
        std::map<std::size_t, std::map<Coordinate, std::size_t, ApproxLess>> groupIds;
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

            const Coordinate center = geometry.global(refElement.position(i, codim));
            const Coordinate mapped = periodicMap(center);
            auto& groupIdsOfType = groupIds.try_emplace(key.first, ApproxLess{ tolerance }).first->second;
            auto it = groupIdsOfType.find(mapped);
            if (it == groupIdsOfType.end())
            {
              it = groupIdsOfType.emplace(mapped, groups.size()).first;
              groups.emplace_back();
            }
            groups[it->second].keys.push_back(key);
            groups[it->second].hasCopy = groups[it->second].hasCopy || !coordinatesEqual(mapped, center, tolerance);
          }
        }

        // a group with a periodic copy is one class, otherwise each entity is its own class
        for (const auto& group : groups)
        {
          if (group.hasCopy)
            addClass(group.keys);
          else
            for (const auto& key : group.keys)
              addClass({ key });
        }
      }

      //! the key (geometry type and index) of an entity
      template<class Entity>
      Key key (const Entity& entity) const
      {
        return Key(GlobalGeometryTypeIndex::index(entity.type()), indexSet_.index(entity));
      }

      //! the key (geometry type and index) of a subentity of an element
      template<class Element>
      Key subKey (const Element& element, int i, int codim) const
      {
        const auto type = referenceElement<typename GridView::ctype, GridView::dimension>(element.type()).type(i, codim);
        return Key(GlobalGeometryTypeIndex::index(type), indexSet_.subIndex(element, i, codim));
      }

      //! the number of classes
      std::size_t size () const { return members_.size(); }

      //! the class of the entity with the given key
      std::size_t classOf (const Key& key) const { return classOf_.at(key); }

      //! the keys of all entities in a class
      const std::vector<Key>& members (std::size_t classId) const { return members_[classId]; }

      //! whether some class has more than one member, i.e. this process stores several copies of an entity
      bool hasCoincidences () const { return hasCoincidences_; }

    private:
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
