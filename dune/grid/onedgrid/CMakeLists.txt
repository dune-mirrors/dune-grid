set(HEADERS nulliteratorfactory.hh
  onedgridentity.hh
  onedgridentitypointer.hh
  onedgridentityseed.hh
  onedgridfactory.hh
  onedgridgeometry.hh
  onedgridhieriterator.hh
  onedgridindexsets.hh
  onedgridleafiterator.hh
  onedgridleveliterator.hh
  onedgridlist.hh
  onedgridintersections.hh
  onedgridintersectioniterators.hh)

install(FILES ${HEADERS} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/grid/onedgrid/)

dune_library_add_sources(
  dunegrid
  SOURCES
  onedgrid.cc
  onedgridfactory.cc
  nulliteratorfactory.cc
  )
