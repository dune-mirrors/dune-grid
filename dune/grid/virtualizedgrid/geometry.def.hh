// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_VIRTUALIZEDGRID_GEOMETRY_DEFINITION_HH
#define DUNE_GRID_VIRTUALIZEDGRID_GEOMETRY_DEFINITION_HH

/** \file
 * \brief The VirtualizedGridGeometry class and its specializations
 */

#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>
#include <dune/geometry/type.hh>
#include <dune/grid/virtualizedgrid/typeerasure/typeerasure.hh>

namespace Dune {
namespace Impl {

template<int mydim, int coorddim, class GridImp>
struct VirtualizedGridGeometryDefinition
{
  using ctype = typename GridImp::ctype;
  using Volume = ctype;
  using JacobianTransposed = FieldMatrix< ctype, mydim, coorddim >;
  using JacobianInverseTransposed = FieldMatrix< ctype, coorddim, mydim >;
  using GlobalCoordinate = FieldVector< ctype, coorddim >;

  struct Interface
  {
    virtual ~Interface () = default;
    virtual GeometryType type () const = 0;
    virtual bool affine() const = 0;
    virtual int corners () const = 0;
    virtual FieldVector<ctype, coorddim> corner (int i) const = 0;
    virtual FieldVector<ctype, coorddim> global (const FieldVector<ctype, mydim>& local) const = 0;
    virtual JacobianTransposed jacobianTransposed (const FieldVector<ctype, mydim>& local) const = 0;
    virtual FieldVector<ctype, mydim> local (const FieldVector<ctype, coorddim>& global) const = 0;
    virtual ctype integrationElement (const FieldVector<ctype, mydim>& local) const = 0;
    virtual JacobianInverseTransposed jacobianInverseTransposed (const FieldVector<ctype, mydim>& local) const = 0;
  };

  template<class Impl>
  struct Model
    : public Impl
  {
    using Impl::Impl;

    GeometryType type () const final { return this->get().type(); }
    bool affine() const final { return this->get().affine(); }
    int corners () const final { return this->get().corners(); }
    FieldVector<ctype, coorddim> corner (int i) const final { return this->get().corner(i); }
    FieldVector<ctype, coorddim> global (const FieldVector<ctype, mydim>& local) const final { return this->get().global(local); }
    JacobianTransposed jacobianTransposed (const FieldVector<ctype, mydim>& local) const final { return this->get().jacobianTransposed(local); }
    FieldVector<ctype, mydim> local (const FieldVector<ctype, coorddim>& global) const final { return this->get().local(global); }
    ctype integrationElement (const FieldVector<ctype, mydim>& local) const final { return this->get().integrationElement(local); }
    JacobianInverseTransposed jacobianInverseTransposed (const FieldVector<ctype, mydim>& local) const final { return this->get().jacobianInverseTransposed(local); }
  };

  using Base = Dune::TypeErasure::TypeErasureBase<Interface,Model>;
};

} // end namespace Impl
} // end namespace Dune

#endif // DUNE_GRID_VIRTUALIZEDGRID_GEOMETRY_DEFINITION_HH
