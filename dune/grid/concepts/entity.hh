// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_CONCEPT_ENTITY_HH
#define DUNE_GRID_CONCEPT_ENTITY_HH

#include <concepts>
#include <utility>

#include <dune/geometry/type.hh>
#include <dune/grid/common/gridenums.hh>
#include <dune/grid/concepts/entity.hh>
#include <dune/grid/concepts/geometry.hh>
#include <dune/grid/concepts/archetypes/entity.hh>

namespace Dune::Concept {

/**
 * @brief Model of an entity seed
 * @ingroup GridConcepts
 * @details Dune::EntitySeed is a template for this model
 */
template<class S>
concept EntitySeed = std::semiregular<S> && requires(const S seed)
{
  { S::codimension  } -> std::convertible_to<int>;
  { seed.isValid()  } -> std::convertible_to<bool>;
};

static_assert(EntitySeed< Archetypes::EntitySeed<0> >);


/**
 * @brief Model of a grid entity for any codimension
 * @ingroup GridConcepts
 * @details Dune::Entity is a template for this model.
 */
template<class E>
concept EntityGeneral = std::regular<E> &&
  Geometry<typename E::Geometry> &&
  EntitySeed<typename E::EntitySeed> &&
requires(const E e, unsigned int codim)
{
  requires E::mydimension == (E::dimension - E::codimension);
  { e.level()            } -> std::convertible_to<int>;
  { e.partitionType()    } -> std::convertible_to<Dune::PartitionType>;
  { e.geometry()         } -> std::convertible_to<typename E::Geometry>;
  { e.type()             } -> std::convertible_to<Dune::GeometryType>;
  { e.subEntities(codim) } -> std::convertible_to<unsigned int>;
  { e.seed()             } -> std::convertible_to<typename E::EntitySeed>;
};

static_assert(EntityGeneral< Archetypes::Entity<2,0> >);


namespace Impl {

  template<class E, int codim>
  concept EntityCodimExtended = requires(const E e, int subEntity)
  {
    { e.template subEntity<codim>(subEntity) } -> EntityGeneral;
  };

  template<class E, int first, int last>
  concept EntityAllCodimsExtended = requires(std::make_integer_sequence<int,last-first> codims)
  {
    []<int... c>(std::integer_sequence<int,c...>)
        requires (EntityCodimExtended<E,(first+c)> &&...) {} (codims);
  };

} // end namespace Impl

/**
 * @brief Model of a grid entity with extended requirements for codimension 0
 * @ingroup GridConcepts
 * @details Dune::Entity of codimension 0 is a template for this model.
 */
template<class E>
concept EntityExtended = EntityGeneral<E> &&
  Geometry<typename E::LocalGeometry> &&
requires(const E e, int maxLevel)
{
  requires (E::codimension == 0);
  { e.father()                   } -> std::convertible_to<E>;
  { e.hasFather()                } -> std::convertible_to<bool>;
  { e.isLeaf()                   } -> std::convertible_to<bool>;
  { e.isRegular()                } -> std::convertible_to<bool>;
  { e.geometryInFather()         } -> std::convertible_to<typename E::LocalGeometry>;
  { e.hbegin(maxLevel)           } -> std::convertible_to<typename E::HierarchicIterator>;
  { e.hend(maxLevel)             } -> std::convertible_to<typename E::HierarchicIterator>;
  { e.isNew()                    } -> std::convertible_to<bool>;
  { e.mightVanish()              } -> std::convertible_to<bool>;
  { e.hasBoundaryIntersections() } -> std::convertible_to<bool>;

  requires std::same_as<E, typename E::template Codim<0>::Entity>;
} &&
Impl::EntityCodimExtended<E,0> &&
Impl::EntityAllCodimsExtended<E,1,E::dimension+1>;

/**
 * @brief Model of a grid entity
 * @ingroup GridConcepts
 * @details Codimension 0 entities are required to have an extended interface.
 * Dune::Entity is a template for this model.
 */
template<class E>
concept Entity = EntityGeneral<E> && ((E::codimension != 0) || EntityExtended<E>);

} // end namespace Dune::Concept


#endif // DUNE_GRID_CONCEPT_ENTITY_HH
