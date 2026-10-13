.. toctree::
   :maxdepth: 2
   :caption: Contents:

   realgrid
   entity_0
   indexset

dune-grid Python API documentation
==================================

**dune-grid** is a module in the Dune ecosystem of finite element software.
dune-grid provides an interface to the grids used for the finite element
method, and several implementations of this interface.

.. note::

   Some parts of the Dune Python API are compiled Just-In-Time, i.e., only
   when you need them. When such parts appear in this documentation,
   their names may contain funky hash strings such as
   "dune.generated.hierarchicalgrid_966e2a5c8356c5b278ccd3acad180f0a".
   Do not let that disturb you.

.. autofunction:: dune.grid.yaspGrid

.. autofunction:: dune.grid.GridFunction
