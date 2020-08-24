include(GNUInstallDirs)

install(FILES
  AddAlbertaFlags.cmake
  AddAmiraMeshFlags.cmake
  AddPsurfaceFlags.cmake
  AddUGFlags.cmake
  DuneGridMacros.cmake
  FindAlberta.cmake
  FindAmiraMesh.cmake
  FindPsurface.cmake
  GridType.cmake
  DESTINATION ${DUNE_INSTALL_MODULEDIR})
