// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#if HAVE_CONFIG_H
#include <config.h>
#endif

#include <dune/grid/yaspgrid.hh>

namespace Dune {

template class YaspGrid<2, EquidistantCoordinates<double, 2>>;
template class YaspGrid<3, EquidistantCoordinates<double, 3>>;
template class YaspGrid<2, EquidistantOffsetCoordinates<double, 2>>;
template class YaspGrid<3, EquidistantOffsetCoordinates<double, 3>>;
template class YaspGrid<2, TensorProductCoordinates<double, 2>>;
template class YaspGrid<3, TensorProductCoordinates<double, 3>>;

} // end namespace Dune
