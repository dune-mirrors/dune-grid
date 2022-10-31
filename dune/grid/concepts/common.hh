// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_CONCEPTS_COMMON_HH
#define DUNE_GRID_CONCEPTS_COMMON_HH

#include <concepts>
#include <cstddef>

namespace Dune::Concept {

template<class B>
concept BooleanTestable = std::convertible_to<B,bool> &&
requires(B b)
{
  { !b } -> std::convertible_to<bool>;
};

template<class I>
concept Integer = std::integral<I> &&
  !std::same_as<I,bool> &&
  !std::same_as<I,char> &&
  !std::same_as<I,unsigned char> &&
  !std::same_as<I,char8_t> &&
  !std::same_as<I,char16_t> &&
  !std::same_as<I,char32_t> &&
  !std::same_as<I,wchar_t>;

} // end namespace Dune::Concept


#endif // DUNE_GRID_CONCEPTS_COMMON_HH
