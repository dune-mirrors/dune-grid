// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_VIRTUALIZED_GRID_IDTYPE_HH
#define DUNE_VIRTUALIZED_GRID_IDTYPE_HH

/**
 * \file
 * \brief The VirtualizedGridIdType class
 */

#include <dune/grid/virtualizedgrid/idtype.def.hh>

namespace Dune {


  /**
   * \brief The IdType class provides a virtualized id type.
   * \ingroup VirtualizedGrid
   *
   */
  class VirtualizedGridIdType
      : public Impl::VirtualizedGridIdTypeDefinition::Base
  {
    using TypeErasureBase = Impl::VirtualizedGridIdTypeDefinition::Base;

  public:

    VirtualizedGridIdType() = default;

    template<class Impl, disableCopyMove<VirtualizedGridIdType,Impl> = 0>
    VirtualizedGridIdType (Impl&& impl)
      : TypeErasureBase(std::forward<Impl>(impl))
    {}

    bool operator== (const VirtualizedGridIdType& other) const
    {
      return this->asInterface().operator==(other.asInterface());
    }

    bool operator!= (const VirtualizedGridIdType& other) const
    {
      return this->asInterface().operator!=(other.asInterface());
    }

    bool operator< (const VirtualizedGridIdType& other) const
    {
      return this->asInterface().operator<(other.asInterface());
    }

    bool operator<= (const VirtualizedGridIdType& other) const
    {
      return this->asInterface().operator<=(other.asInterface());
    }

    std::string str () const
    {
      return this->asInterface().str();
    }

    std::size_t hash () const
    {
      return this->asInterface().hash();
    }
  };

  inline std::ostream &operator<< ( std::ostream &out, const VirtualizedGridIdType &idtype )
  {
    return out << idtype.str();
  }

} // namespace Dune


namespace std
{
  template <> struct hash<Dune::VirtualizedGridIdType>
  {
    size_t operator()(const Dune::VirtualizedGridIdType& x) const
    {
      return x.hash();
    }
  };
}

#endif  // #define DUNE_VIRTUALIZED_GRID_IDTYPE_HH
