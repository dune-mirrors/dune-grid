// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_VIRTUALIZEDGRID_INDEXSETS_DEFINITION_HH
#define DUNE_VIRTUALIZEDGRID_INDEXSETS_DEFINITION_HH

/** \file
 * \brief The index and id sets for the VirtualizedGrid class
 */

#include <cstddef>
#include <utility>
#include <vector>

#include <dune/geometry/dimension.hh>
#include <dune/geometry/type.hh>
#include <dune/grid/virtualizedgrid/entity.hh>

namespace Dune {
namespace Impl {

template<class GridImp>
struct VirtualizedGridLevelIndexSetDefinition
{
  using Types = std::vector<GeometryType>;

  enum { dim = GridImp::dimension };

  template<int codim>
  struct InterfaceCodim
  {
    using Entity = typename GridImp::Traits::template Codim<codim>::Entity;

    virtual ~InterfaceCodim () = default;
    virtual int index (Codim<codim>, const Entity& e) const = 0;
    virtual int subIndex (Codim<codim>, const Entity& e, int i, int cd) const = 0;
    virtual bool contains (Codim<codim>, const Entity& e) const = 0;
  };

  template<int... codims>
  struct InterfaceImpl
      : virtual InterfaceCodim<codims>...
  {
    virtual ~InterfaceImpl () = default;
    virtual std::size_t size (int codim) const = 0;
    virtual std::size_t size (GeometryType type) const = 0;
    virtual Types types (int codim) const = 0;

    using InterfaceCodim<codims>::index...;
    using InterfaceCodim<codims>::subIndex...;
    using InterfaceCodim<codims>::contains...;
  };

  template<class Seq> struct Interface_t;
  template<int... codims>
  struct Interface_t<std::integer_sequence<int,codims...>> {
    using type = InterfaceImpl<codims...>;
  };
  using Interface = typename Interface_t<std::make_integer_sequence<int,dim+1>>::type;


  template<class Derived, class Impl, int codim>
  struct ImplementationCodim
    : virtual InterfaceCodim<codim>
  {
    using Entity = typename GridImp::Traits::template Codim<codim>::Entity;
    using EntityImpl = typename Impl::Wrapped::template Codim<codim>::Entity;
    using EntityModel = typename VirtualizedGridEntity<codim,dim,GridImp>::template Implementation<const EntityImpl>;
    // using EntityModel = typename VirtualizedGridEntity<codim,dim,GridImp>::template WrapperImplementation<EntityImpl>;

    int index (Codim<codim>, const Entity& e) const final {
      return derived().get().index(upcast<EntityModel>(e));
    }
    int subIndex (Codim<codim>, const Entity& e, int i, int cd) const final {
      return derived().get().template subIndex<codim>(upcast<EntityModel>(e), i, cd);
    }
    bool contains (Codim<codim>, const Entity& e) const final {
      return derived().get().contains(upcast<EntityModel>(e));
    }

  private:
    const Derived& derived () const { return static_cast<const Derived&>(*this); }
  };

  template<class Impl, int... codims>
  struct ModelImpl
    : public Impl
    , public ImplementationCodim<ModelImpl<Impl,codims...>, Impl, codims>...
  {
    using Impl::Impl;

    std::size_t size (int codim) const final { return this->get().size(codim); }
    std::size_t size (GeometryType type) const final { return this->get().size(type); }
    Types types (int codim) const final { return this->get().types(codim); }
  };

  template<class I, class Seq> struct Model_t;
  template<class I, int... codims>
  struct Model_t<I,std::integer_sequence<int,codims...>> {
    using type = ModelImpl<I,codims...>;
  };
  template<class I>
  using Model = typename Model_t<I,std::make_integer_sequence<int,dim+1>>::type;

  using Base = Dune::TypeErasure::TypeErasureBase<Interface,Model>;
};

} // end namespace Impl
} // end namespace Dune


#endif // DUNE_VIRTUALIZEDGRID_INDEXSETS_DEFINITION_HH
