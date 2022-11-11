// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_VIRTUALIZED_GRID_ENTITY_SEED_HH
#define DUNE_VIRTUALIZED_GRID_ENTITY_SEED_HH

/**
 * \file
 * \brief The VirtualizedGridEntitySeed class
 */


namespace Dune {


  /**
   * \brief The EntitySeed class provides the minimal information needed to restore an Entity using the grid.
   * \ingroup VirtualizedGrid
   *
   */
  template<int codim, class GridImp>
  class VirtualizedGridEntitySeed
  {
  protected:

    // Entity type of the grid
    typedef typename GridImp::Traits::template Codim<codim>::Entity Entity;

  public:
    // VIRTUALIZATION BEGIN
    struct Interface
    {
      virtual ~Interface () = default;
      virtual Interface *clone () const = 0;
      virtual bool isValid () const = 0;
    };

    template< class I >
    struct DUNE_PRIVATE Implementation final
      : public Interface
    {
      Implementation ( I&& i ) : impl_( std::forward<I>(i) ) {}
      Implementation *clone() const override { return new Implementation( *this ); }
      bool isValid() const override { return impl().isValid(); }

      const auto &impl () const { return impl_; }
      auto &impl () { return impl_; }

    private:
      I impl_;
    };
    // VIRTUALIZATION END

    enum {codimension = codim};

    /**
     * \brief Construct an empty (i.e. isValid() == false) seed.
     */
    VirtualizedGridEntitySeed()
    {}

    /**
     * \brief Create EntitySeed from implementation entity
     */
    template< class ImplEntitySeed >
    VirtualizedGridEntitySeed(ImplEntitySeed&& implEntitySeed)
    : impl_( new Implementation<ImplEntitySeed>( std::forward<ImplEntitySeed>(implEntitySeed) ) )
    {}

    VirtualizedGridEntitySeed(const VirtualizedGridEntitySeed& other)
    : impl_( other.impl_ ? other.impl_->clone() : nullptr )
    {}

    VirtualizedGridEntitySeed ( VirtualizedGridEntitySeed && ) = default;

    VirtualizedGridEntitySeed& operator=(const VirtualizedGridEntitySeed& other)
    {
      impl_.reset( other.impl_ ? other.impl_->clone() : nullptr );
      return *this;
    }

    VirtualizedGridEntitySeed& operator=( VirtualizedGridEntitySeed&& ) = default;

    /**
     * \brief Check whether it is safe to create an Entity from this Seed
     */
    bool isValid() const
    {
      return impl_ && impl_->isValid();
    }

    std::unique_ptr<Interface> impl_;
  };

} // namespace Dune

#endif  // #define DUNE_VIRTUALIZED_GRID_ENTITY_SEED_HH