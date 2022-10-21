#ifndef DUNE_GRID_CONCEPTS_GRIDVIEW_HH
#define DUNE_GRID_CONCEPTS_GRIDVIEW_HH

#include <dune/grid/concepts/entityset.hh>

#include <dune/grid/common/capabilities.hh>

#include <dune/common/indices.hh>

namespace Dune::Concept {

namespace Impl {

  template<class ES, int codim,  Dune::PartitionIteratorType partition>
  concept EntityPartitionSpan = requires(ES es) {
    requires EntityIterator<typename ES::template Codim<codim>::template Partition<partition>::Iterator>;
    { es.template begin<codim,partition>()  } -> std::convertible_to< typename ES::template Codim<codim>::template Partition<partition>::Iterator>;
    { es.template end<codim,partition>()    } -> std::convertible_to< typename ES::template Codim<codim>::template Partition<partition>::Iterator>;
  };


  template<class ES, int codim>
  concept EntitySetAllPartitions = requires(ES es) {
    requires EntitySet<ES,codim>;
    requires EntityPartitionSpan<ES,codim,Dune::PartitionIteratorType::InteriorBorder_Partition>;
    requires EntityPartitionSpan<ES,codim,Dune::PartitionIteratorType::Overlap_Partition>;
    requires EntityPartitionSpan<ES,codim,Dune::PartitionIteratorType::OverlapFront_Partition>;
    requires EntityPartitionSpan<ES,codim,Dune::PartitionIteratorType::All_Partition>;
    requires EntityPartitionSpan<ES,codim,Dune::PartitionIteratorType::Ghost_Partition>;
  };

  template<class ES, int codim>
  concept EntitySetAllPartitionsIfSupported
    = ((not Dune::Capabilities::hasEntityIterator<ES,codim>::value) || EntitySetAllPartitions<ES,codim>);

  template<class ES, int... codim>
    requires (EntitySetAllPartitionsIfSupported<ES,(codim+1)> &&...)
  void requireEntitySetAllPartitionsIfSupported(std::integer_sequence<int,codim...>) {};

}

/**
 * @brief Model of a grid view
 * @ingroup GridConcepts
 * @details Dune::GridView is a template for this model
 */
template<class GV>
concept GridView = requires {
  requireEntitySetAllPartitionsIfSupported<GV>(std::make_integer_sequence<int,GV::dimension>{});
};

}  // end namespace Dune::Concept

#endif // DUNE_GRID_CONCEPTS_GRIDVIEW_HH
