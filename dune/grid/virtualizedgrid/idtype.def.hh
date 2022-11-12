// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_VIRTUALIZED_GRID_IDTYPE_DEFINITION_HH
#define DUNE_VIRTUALIZED_GRID_IDTYPE_DEFINITION_HH

#include <cstddef>
#include <sstream>
#include <string>

#include <dune/common/visibility.hh>
#include <dune/grid/virtualizedgrid/typeerasure/typeerasure.hh>

namespace Dune {
namespace Impl {

struct VirtualizedGridIdTypeDefinition
{
  struct Interface
  {
    virtual ~Interface () = default;
    virtual bool operator== (const Interface& other) const = 0;
    virtual bool operator!= (const Interface& other) const = 0;
    virtual bool operator< (const Interface& other) const = 0;
    virtual bool operator<= (const Interface& other) const = 0;
    virtual std::size_t hash () const = 0;
    virtual std::string str() const = 0;
  };

  template<class Impl>
  struct Model
    : public Impl
  {
    using Impl::Impl;

    bool operator== (const Interface& other) const final
    {
      return this->get() == static_cast<const Model<Impl>&>(other).get();
    }

    bool operator!= (const Interface& other) const final
    {
      return this->get() != static_cast<const Model<Impl>&>(other).get();
    }

    bool operator< (const Interface& other) const final
    {
      return this->get() < static_cast<const Model<Impl>&>(other).get();
    }

    bool operator<= (const Interface& other) const final
    {
      return this->get() <= static_cast<const Model<Impl>&>(other).get();
    }

    std::string str () const final
    {
      std::stringstream ss;
      ss << this->get() << std::endl;
      return ss.str();
    }

    std::size_t hash () const final
    {
      return std::hash<typename Impl::Wrapped>()(this->get());
    }
  };

  using Base = Dune::TypeErasure::TypeErasureBase<Interface,Model>;
};

} // end namespace Impl
} // end namespace Dune


#endif  // DUNE_VIRTUALIZED_GRID_IDTYPE_DEFINITION_HH
