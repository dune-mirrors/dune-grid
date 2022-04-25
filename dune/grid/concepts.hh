#ifndef DUNE_GRID_CONCEPTS_HH
#define DUNE_GRID_CONCEPTS_HH

/**
 * This file contains a convenience definition that checks if concepts are available
 * If DUNE_GRID_HAVE_CONCEPTS is true, the dune-grid concepts are available and
 * have been included.
 */

#if DUNE_HAVE_CXX_UNEVALUATED_CONTEXT_LAMBDA && __cpp_concepts > 201907L && __cpp_lib_concepts > 202002L

// Include all concept headers
#include <dune/grid/concepts/entity.hh>
#include <dune/grid/concepts/entityiterator.hh>
#include <dune/grid/concepts/intersection.hh>
#include <dune/grid/concepts/intersectioniterator.hh>
#include <dune/grid/concepts/geometry.hh>
#include <dune/grid/concepts/indexset.hh>
#include <dune/grid/concepts/gridview.hh>
#include <dune/grid/concepts/grid.hh>

//! Grid concepts are availalbe
#define DUNE_GRID_HAVE_CONCEPTS 1

#else

//! Grid concepts are not availalbe
#define DUNE_GRID_HAVE_CONCEPTS 0

#endif

#endif // DUNE_GRID_CONCEPTS_HH
