// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_CONCEPTS_ARCHETYPES_GEOMETRY_HH
#define DUNE_GRID_CONCEPTS_ARCHETYPES_GEOMETRY_HH

#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>
#include <dune/geometry/type.hh>

namespace Dune::Concept::Archetypes {

struct ReferenceElement {};

template <int mydim, int cdim = mydim>
struct Geometry
{
  static constexpr int mydimension = mydim;
  static constexpr int coorddimension = cdim;

  using ctype = double;
  using LocalCoordinate = Dune::FieldVector<ctype, mydim>;
  using GlobalCoordinate = Dune::FieldVector<ctype, cdim>;
  using Volume = ctype;
  using Jacobian = Dune::FieldMatrix<ctype, cdim, mydim>;
  using JacobianTransposed = Dune::FieldMatrix<ctype, mydim, cdim>;
  using JacobianInverse = Dune::FieldMatrix<ctype, mydim, cdim>;
  using JacobianInverseTransposed = Dune::FieldMatrix<ctype, cdim, mydim>;

  Dune::GeometryType type () const { return Dune::GeometryTypes::none(mydim); }
  bool affine () const { return false; }
  int corners() const { return 0; }

  GlobalCoordinate corner(int i) const { return {}; }
  GlobalCoordinate global(const LocalCoordinate&) const { return {}; }
  LocalCoordinate local(const GlobalCoordinate&) const { return {}; }
  GlobalCoordinate center() const { return {}; }

  Volume integrationElement(const LocalCoordinate&) const { return 0; }
  Volume volume() const { return 0; }

  Jacobian jacobian(const LocalCoordinate&) const { return {}; }
  JacobianTransposed jacobianTransposed(const LocalCoordinate&) const { return {}; }
  JacobianInverse jacobianInverse(const LocalCoordinate&) const { return {}; }
  JacobianInverseTransposed jacobianInverseTransposed(const LocalCoordinate&) const { return {}; }

  friend Archetypes::ReferenceElement referenceElement(const Geometry&) { return Archetypes::ReferenceElement{}; }
};

} // end namespace Dune::Concept::Archetypes


#endif // DUNE_GRID_CONCEPTS_ARCHETYPES_GEOMETRY_HH
