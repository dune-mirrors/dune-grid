// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_CONCEPTS_GRIDVIEW_HH
#define DUNE_GRID_CONCEPTS_GRIDVIEW_HH

#include <concepts>

#include <dune/grid/common/capabilities.hh>
#include <dune/grid/common/gridenums.hh>
#include <dune/grid/concepts/entityiterator.hh>
#include <dune/grid/concepts/geometry.hh>
#include <dune/grid/concepts/indexidset.hh>
#include <dune/grid/concepts/intersection.hh>
#include <dune/grid/concepts/intersectioniterator.hh>
#include <dune/grid/concepts/archetypes/datahandle.hh>

#include <dune/common/indices.hh>

namespace Dune::Concept {
namespace Impl {

  template<class ES, int codim, Dune::PartitionIteratorType partition,
    class Traits = typename ES::template Codim<codim>::template Partition<partition>>
  concept EntitySetPartition =
    EntityIterator<typename Traits::Iterator> &&
  requires(const ES es)
  {
    { es.template begin<codim,partition>() } -> std::convertible_to<typename Traits::Iterator>;
    { es.template end<codim,partition>()   } -> std::convertible_to<typename Traits::Iterator>;
  };

  template<class ES, int codim,
    class Traits = typename ES::template Codim<codim>>
  concept EntitySetCodim =
    Geometry<typename Traits::Geometry> &&
    Geometry<typename Traits::LocalGeometry> &&
    EntityIterator<typename Traits::Iterator> &&
  requires(const ES es)
  {
    { es.template begin<codim>() } -> std::convertible_to<typename Traits::Iterator>;
    { es.template end<codim>()   } -> std::convertible_to<typename Traits::Iterator>;

    requires (codim != 0) || requires(const typename Traits::Entity& entity)
    {
      { es.ibegin(entity) } -> std::convertible_to<typename ES::IntersectionIterator>;
      { es.iend(entity)   } -> std::convertible_to<typename ES::IntersectionIterator>;
    };
  };

  template<class ES, int codim>
  concept EntitySetAllPartitions = EntitySetCodim<ES,codim> &&
    EntitySetPartition<ES,codim,Dune::PartitionIteratorType::InteriorBorder_Partition> &&
    EntitySetPartition<ES,codim,Dune::PartitionIteratorType::Overlap_Partition> &&
    EntitySetPartition<ES,codim,Dune::PartitionIteratorType::OverlapFront_Partition> &&
    EntitySetPartition<ES,codim,Dune::PartitionIteratorType::All_Partition> &&
    EntitySetPartition<ES,codim,Dune::PartitionIteratorType::Ghost_Partition>;

  template<class GV, class Grid, int codim>
    requires Dune::Capabilities::hasEntityIterator<Grid,codim>::v
  void requireGridViewCodim()
    requires EntitySetAllPartitions<GV,codim> {}

  template<class GV, class Grid, int codim>
    requires (not Dune::Capabilities::hasEntityIterator<Grid,codim>::v)
  void requireGridViewCodim() {}

  template <class GV, int dim>
  concept GridViewAllCodims = requires(std::make_integer_sequence<int,dim+1> dims)
  {
    []<int... d>(std::integer_sequence<int,d...>) requires
      requires { (requireGridViewCodim<GV,typename GV::Grid,(dim-d)>(),...); } {} (dims);
  };

} // end namespace Impl

/**
 * @brief Model of a grid view
 * @ingroup GridConcepts
 * @details Dune::GridView is a template for this model
 */
template<class GV>
concept GridView = std::copyable<GV> &&
  IndexSet<typename GV::IndexSet> &&
  Intersection<typename GV::Intersection> &&
  IntersectionIterator<typename GV::IntersectionIterator> &&
requires(const GV gv, int codim, Dune::GeometryType type)
{
  typename GV::Traits;
  typename GV::ctype;
  { GV::conforming        } -> std::convertible_to<bool>;
  { GV::dimension         } -> std::convertible_to<int>;
  { GV::dimensionworld    } -> std::convertible_to<int>;
  { gv.grid()             } -> std::convertible_to<const typename GV::Grid&>;
  { gv.indexSet()         } -> std::convertible_to<const typename GV::IndexSet&>;
  { gv.size(codim)        } -> std::convertible_to<int>;
  { gv.size(type)         } -> std::convertible_to<int>;
  { gv.comm()             } -> std::convertible_to<typename GV::CollectiveCommunication>;
  { gv.overlapSize(codim) } -> std::convertible_to<int>;
  { gv.ghostSize(codim)   } -> std::convertible_to<int>;

  requires requires(Archetypes::CommDataHandle<std::byte>& handle,
                    InterfaceType iface, CommunicationDirection dir)
  {
    gv.communicate(handle, iface, dir);
  };
} && Impl::GridViewAllCodims<GV, GV::dimension>;

}  // end namespace Dune::Concept

#endif // DUNE_GRID_CONCEPTS_GRIDVIEW_HH
