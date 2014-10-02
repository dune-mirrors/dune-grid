// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_YASPGRIDLEVELITERATOR_HH
#define DUNE_GRID_YASPGRIDLEVELITERATOR_HH

/** \file
 * \brief The YaspLevelIterator class
 */

#include <cstddef>
#include <iterator>

namespace Dune {


  /** \brief Iterates over entities of one grid level
   */
  template<int codim, PartitionIteratorType pitype, class GridImp>
  class YaspLevelIterator :
    public YaspEntityPointer<codim,GridImp>
  {
    //! know your own dimension
    enum { dim=GridImp::dimension };
    //! know your own dimension of world
    enum { dimworld=GridImp::dimensionworld };
    typedef typename GridImp::ctype ctype;
  public:
    typedef typename GridImp::template Codim<codim>::Entity Entity;
    typedef typename GridImp::YGridLevelIterator YGLI;
    typedef typename GridImp::YGrid::Iterator I;
    typedef std::random_access_iterator_tag iterator_category;

    //! default constructor
    YaspLevelIterator ()
    {}

    //! constructor
    YaspLevelIterator (const GridImp * yg, const YGLI & g, const I& it) :
      YaspEntityPointer<codim,GridImp>(yg,g,it) {}

    //! copy constructor
    YaspLevelIterator (const YaspLevelIterator& i) :
      YaspEntityPointer<codim,GridImp>(i) {}

    //! increment
    void increment()
    {
      ++(this->_it);
    }

    //! decrement
    void decrement()
    {
      --(this->_it);
    }

    //! advance iterator by given amount, forward or backward
    void advance(std::ptrdiff_t distance)
    {
      this->_it += distance;
    }

    //! compute distance from this iterator to another
    int distanceTo(const YaspLevelIterator &other) const
    {
      return other._it - this->_it;
    }
  };

}

#endif   // DUNE_GRID_YASPGRIDLEVELITERATOR_HH
