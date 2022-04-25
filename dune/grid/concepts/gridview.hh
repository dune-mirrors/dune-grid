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
  requires (not Dune::Capabilities::hasEntityIterator<ES,codim>::value) || EntitySetAllPartitions<ES,codim>
  void requireEntitySetAllPartitionsIfSupported();

}

/**
 * @brief Model of a grid view
 * @ingroup GridConcepts
 * @details Dune::GridView is a template for this model
 */
template<class GV>
concept GridView = requires(GV gv)
{
  requires Impl::EntitySetAllPartitions<GV,0>; // minimal requirement
  Dune::unpackIntegerSequence([](auto... i){
    (Impl::requireEntitySetAllPartitionsIfSupported<GV,i+1>(), ...);
  }, std::make_index_sequence<GV::dimension>{});
};

}  // end namespace Dune::Concept

#endif // DUNE_GRID_CONCEPTS_GRIDVIEW_HH
