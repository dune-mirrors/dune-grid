set(HEADERS nulliteratorfactory.hh
  onedgridentity.hh
  onedgridentityseed.hh
  onedgridfactory.hh
  onedgridhieriterator.hh
  onedgridindexsets.hh
  onedgridleafiterator.hh
  onedgridleveliterator.hh
  onedgridlist.hh
  onedgridintersections.hh
  onedgridintersectioniterators.hh
  onedgridviews.hh)

install(FILES ${HEADERS} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/grid/onedgrid/)

# Some of the OnedGrid headers do not pass make headercheck,
# and they don't need to necessarily. Remove this line when headers are fixed
exclude_from_headercheck( onedgridleafiterator.hh onedgridleveliterator.hh onedgridhieriterator.hh )

dune_add_library(onedgrid OBJECT
  SOURCES
  onedgrid.cc
  onedgridfactory.cc
  nulliteratorfactory.cc)
