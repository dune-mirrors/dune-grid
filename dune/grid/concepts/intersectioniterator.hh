#ifndef DUNE_GRID_CONCEPTS_INTERSECTION_ITERATOR_HH
#define DUNE_GRID_CONCEPTS_INTERSECTION_ITERATOR_HH

#include <dune/grid/concepts/intersection.hh>

#include <dune/common/concept.hh>

namespace Dune {
  namespace Concept
  {

#if DUNE_HAVE_CXX_CONCEPTS
namespace Concept {
    template<class I>
    concept IntersectionIterator = requires(I i)
    {
      requires Intersection<typename I::Intersection>;
      i++; // FIXME set type requirement
      ++i; // FIXME set type requirement
      *i; // FIXME set type requirement
      i.operator ->(); // FIXME set type requirement
      { i==i } -> Std::convertible_to<bool>;
      { i!=i } -> Std::convertible_to<bool>;
      requires Std::default_initializable<I>;
    };
}
#endif

    struct IntersectionIterator
    {
      template<class I>
      auto require(I&& i) -> decltype(
        requireConcept<Dune::Concept::Intersection, typename I::Intersection>(),
        i++, // FIXME set type requirement
        ++i, // FIXME set type requirement
        *i, // FIXME set type requirement
        i.operator ->(), // FIXME set type requirement
        requireConvertible<bool>(i==i),
        requireConvertible<bool>(i!=i),
        I{}, // default constructible
        I{i} // copy constructible
      );
    };

  }

  template <class I>
  constexpr void expectIntersectionIterator()
  {
#if DUNE_HAVE_CXX_CONCEPTS
    static_assert(Dune::Concept::Concept::IntersectionIterator<I>);
#else
    static_assert(models<Concept::IntersectionIterator, I>());
#endif
  }

}  // end namespace Dune

#endif