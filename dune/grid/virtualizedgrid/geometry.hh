// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_VIRTUALIZEDGRIDGEOMETRY_HH
#define DUNE_VIRTUALIZEDGRIDGEOMETRY_HH

/** \file
 * \brief The VirtualizedGridGeometry class and its specializations
 */

#include <dune/common/fmatrix.hh>
#include <dune/common/typetraits.hh>
#include <dune/grid/common/geometry.hh>
#include <dune/grid/virtualizedgrid/geometry.def.hh>

namespace Dune {

  template<int mydim, int coorddim, class GridImp>
  class VirtualizedGridGeometry
    : public GeometryDefaultImplementation <mydim, coorddim, GridImp, VirtualizedGridGeometry>
    , public Impl::VirtualizedGridGeometryDefinition<mydim,coorddim,GridImp>::Base
  {
    using TypeErasureBase = typename Impl::VirtualizedGridGeometryDefinition<mydim,coorddim,GridImp>::Base;

  public:
    typedef typename GridImp::ctype ctype;
    typedef ctype Volume;
    typedef FieldMatrix< ctype, mydim, coorddim > JacobianTransposed;
    typedef FieldMatrix< ctype, coorddim, mydim > JacobianInverseTransposed;
    typedef FieldVector< ctype, coorddim > GlobalCoordinate;

    struct Cache
    {
      template<class Geometry>
      void update (const Geometry& geometry)
      {
        type_ = geometry.type();
        affine_ = geometry.affine();
        corners_ = geometry.corners();

        corner_.resize(corners_);
        for (int i = 0; i < corners_; ++i)
          corner_[i] = geometry.corner(i);

        volume_ = geometry.volume();
        center_ = geometry.center();
      }

      //! return the element type identifier
      GeometryType type () const {
        return type_;
      }

      //! return wether we have an affine mapping
      bool affine () const {
        return affine_;
      }

      //! return the number of corners of this element. Corners are numbered 0...n-1
      int corners () const {
        return corners_;
      }

      //! access to coordinates of corners. Index is the number of the corner
      const FieldVector<ctype, coorddim>& corner (int i) const {
        return corner_[i];
      }

      //! return volume of the geometry
      Volume volume () const {
        return volume_;
      }

      //! return center of the geometry
      GlobalCoordinate center () const {
        return center_;
      }

    private:
      GeometryType type_{};
      bool affine_{};
      int corners_{};
      std::vector<FieldVector<ctype, coorddim>> corner_{};
      Volume volume_{};
      GlobalCoordinate center_{};
    };

  public:

    /** constructor from host geometry
     */
    template<class Impl>
    VirtualizedGridGeometry(Impl&& impl)
      : TypeErasureBase(std::forward<Impl>(impl))
    {
      update(impl);
    }

    /** \brief Return the element type identifier
     */
    GeometryType type () const {
      return cache().type();
    }

    // return wether we have an affine mapping
    bool affine () const {
      return cache().affine();
    }

    //! return the number of corners of this element. Corners are numbered 0...n-1
    int corners () const {
      return cache().corners();
    }

#ifndef DUNE_VIRTUALIZEDGRID_NO_CACHE
    //! return volume of the geometry
    Volume volume () const {
      return cache().volume();
    }

    //! return center of the geometry
    GlobalCoordinate center () const {
      return cache().center();
    }
#endif

    //! access to coordinates of corners. Index is the number of the corner
    FieldVector<ctype, coorddim> corner (int i) const {
      return cache().corner(i);
    }

    /** \brief Maps a local coordinate within reference element to
     * global coordinate in element  */
    FieldVector<ctype, coorddim> global (const FieldVector<ctype, mydim>& local) const {
      return this->asInterface().global(local);
    }

    /** \brief Return the transposed of the Jacobian
     */
    JacobianTransposed jacobianTransposed (const FieldVector<ctype, mydim>& local) const {
      return this->asInterface().jacobianTransposed(local);
    }

    /** \brief Maps a global coordinate within the element to a
     * local coordinate in its reference element */
    FieldVector<ctype, mydim> local (const FieldVector<ctype, coorddim>& global) const {
      return this->asInterface().local(global);
    }

    ctype integrationElement (const FieldVector<ctype, mydim>& local) const {
      return this->asInterface().integrationElement(local);
    }

    //! The Jacobian matrix of the mapping from the reference element to this element
    JacobianInverseTransposed jacobianInverseTransposed (const FieldVector<ctype, mydim>& local) const {
      return this->asInterface().jacobianInverseTransposed(local);
    }

    const auto& cache() const {
#ifndef DUNE_VIRTUALIZEDGRID_NO_CACHE
      return cache_;
#else
      return this->asInterface();
#endif
    }

    template<class Grid>
    void update (const Grid& grid) {
#ifndef DUNE_VIRTUALIZEDGRID_NO_CACHE
      cache_.update(grid);
#endif
    }

  private:
    Cache cache_;
  };

}  // namespace Dune

#endif
