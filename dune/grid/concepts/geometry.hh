// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_CONCEPTS_GEOMETRY_HH
#define DUNE_GRID_CONCEPTS_GEOMETRY_HH

#include <concepts>

#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>
#include <dune/geometry/type.hh>
#include <dune/grid/concepts/archetypes/geometry.hh>

namespace Dune::Concept {

template<class R>
concept ReferenceElement = true;

static_assert(ReferenceElement< Archetypes::ReferenceElement >);


/**
 * @brief Model of a geometry object
 * @ingroup GridConcepts
 * @details Dune::Geometry is a template for this model
 */
template<class G>
concept Geometry = requires(const G g, typename G::GlobalCoordinate global, typename G::LocalCoordinate local)
{
  typename G::ctype;
  { G::mydimension                     } -> std::convertible_to<int>;
  { G::coorddimension                  } -> std::convertible_to<int>;
  { g.type()                           } -> std::same_as<Dune::GeometryType>;
  { g.affine()                         } -> std::convertible_to<bool>;
  { g.corners()                        } -> std::convertible_to<int>;
  { g.corner(/*i*/ int{})              } -> std::same_as<typename G::GlobalCoordinate>;
  { g.global(local)                    } -> std::same_as<typename G::GlobalCoordinate>;
  { g.local(global)                    } -> std::same_as<typename G::LocalCoordinate>;
  { g.integrationElement(local)        } -> std::same_as<typename G::Volume>;
  { g.volume()                         } -> std::same_as<typename G::Volume>;
  { g.center()                         } -> std::same_as<typename G::GlobalCoordinate>;
  { g.jacobian(local)                  } -> std::same_as<typename G::Jacobian>;
  { g.jacobianInverse(local)           } -> std::same_as<typename G::JacobianInverse>;
  { g.jacobianTransposed(local)        } -> std::same_as<typename G::JacobianTransposed>;
  { g.jacobianInverseTransposed(local) } -> std::same_as<typename G::JacobianInverseTransposed>;
  { referenceElement(g)                } -> ReferenceElement;
};

static_assert(Geometry< Archetypes::Geometry<2,3> >);

} // end namespace Dune::Concept

#endif // DUNE_GRID_CONCEPTS_GEOMETRY_HH
