// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_ENTITYITERATOR_HH
#define DUNE_GRID_ENTITYITERATOR_HH

#include <cstddef>
#include <iterator>

#include <dune/common/dotproduct.hh> // for AlwaysVoid

#include <dune/grid/common/entitypointer.hh>

namespace Dune
{

  //! Determine iterator category for an entity interator implementation
  /**
   * This handles the default case when the iterator implementation does not
   * specify a category; we just assume category std::forward_iterator_tag.
   *
   * \note IteratorImp must declare the iterator category as member type
   *       IteratorImp::iterator_category (or not specify a category, in which
   *       case the implementation is axpected to implement a forward
   *       iterator).  In particular, since IteratorImp is not in itself an
   *       iterator, this class does not use std::iterator_traits.
   */
  template< class IteratorImp, class = void >
  struct EntityIteratorCategory
  {
    //! The iterator category of IteratorImp
    typedef std::forward_iterator_tag type;
  };

  //! Determine iterator tag for an entity interator implementation
  /**
   * This handles the case when the iterator implementation does specify a
   * category; we just export that category.
   */
  template< class IteratorImp >
  struct EntityIteratorCategory<
    IteratorImp,
    typename AlwaysVoid<typename IteratorImp::iterator_category>::type>
  {
    //! The iterator category of IteratorImp
    typedef typename IteratorImp::iterator_category type;
  };

  //! Base class for the EntityIterators facade class
  /**
   * Each specialization of this template implements its own level of iterator
   * functionality.  EntityIterator derives from EntityIteratorBase<...,
   * EntityIteratorCategory<IteratorImp>::type>.  The specialization for
   * random access iterators derives from the specialization for forward
   * iterators and does not implement forward iterator functionality itself.
   *
   * This general template implements the methods to make EntityIterator a
   * forward iterator.  The implementation IteratorImp must
   * - be default constructible,
   * - be copy constructible,
   * - be copy assignable,
   * - provide the method increment() to advance the iterator, and
   * - must provide the methods equals() and dereference(), which are
   *   already required by the entity pointer.
   * .
   *
   * \tparam codim       Codimention to iterate over.
   * \tparam Grid        Type of the grid.
   * \tparam IteratorImp Class implementing the iterator (see above for
   *                     requirements).
   * \tparam MostDerived Type of the final EntityIterator.  This is used to
   *                     return the correct type from iterator interface
   *                     methods.
   * \tparam category    Iterator category.  This is used to select the
   *                     correct specialization.
   */
  template< int codim, class Grid, class IteratorImp, class MostDerived,
            class category >
  class EntityIteratorBase
    : public EntityPointer< Grid, IteratorImp >
  {
    static_assert(std::is_base_of<std::forward_iterator_tag, category>::value,
                  "The iterator_category member of IteratorImp must be at "
                  "least std::forward_iterator_tag, or it must not exist, in "
                  "which case std::forward_iterator_tag is assumed.");

    typedef EntityPointer< Grid, IteratorImp > Base;

  protected:
    using Base::realIterator;

  public:
    //! export type of entity
    typedef typename Grid::template Codim< codim >::Entity Entity;

    //! export type of difference (for std::iterator_traits)
    typedef std::ptrdiff_t difference_type;
    //! export type of values (for std::iterator_traits)
    typedef const Entity value_type;
    //! export type of pointers (for std::iterator_traits)
    typedef value_type *pointer;
    //! export type of references (for std::iterator_traits)
    typedef value_type &reference;
    //! export iterator category (for std::iterator_traits)
    /**
     * This base class implements the forward iterator interface, so we
     * unconditionally set this to std::forward_iterator_tag.  If a derived
     * class implements additional functionality, it should adjust this member
     * type accordingly.
     */
    typedef std::forward_iterator_tag iterator_category;

    /** \brief prefix increment operator */
    MostDerived &operator++ ()
    {
      realIterator.increment();
      return static_cast<MostDerived&>(*this);
    }

    /** \brief postfix increment operator */
    MostDerived operator++ (int)
    {
      MostDerived tmp(static_cast<MostDerived&>(*this));
      realIterator.increment();
      return tmp;
    }

    /** \brief default construct (undefined) iterator */
    EntityIteratorBase ( )
    {}

    /** \name Implementor's interface
     *  \{
     */

    /** \brief copy constructor from implementaton */
    EntityIteratorBase ( const IteratorImp &imp )
      : Base( imp )
    {}

    /** \} */
  };

  //! Base class for the EntityIterators facade class
  /**
   * This specialization implements the methods to make EntityIterator a
   * random access iterator.  The implementation IteratorImp must
   * - fulfill all the requirements for EntityIteratorBase<codim, Grid,
   *   IteratorImp, MostDerived, std::forward_iterator_tag>,
   * - provide the method advance(n) to advance the iterator by an arbitrary
   *   (positive or negative) amount, and
   * - provide the method distanceTo(otherIT) to compute the distance
   *   between two iterators.
   * .
   *
   * \tparam codim       Codimention to iterate over.
   * \tparam Grid        Type of the grid.
   * \tparam IteratorImp Class implementing the iterator (see above for
   *                     requirements).
   * \tparam MostDerived Type of the final EntityIterator.  This is used to
   *                     return the correct type from iterator interface
   *                     methods.
   *
   * This class provides the forward iterator functionality by deriving from
   * EntityIteratorBase<codim, Grid, IteratorImp, MostDerived,
   * std::forward_iterator_tag>, only the more advanced functionality is
   * implemented here.
   *
   * \note Subscription (operator[]()) is not implemented.  A reasonable
   *       approximation needs Entities as first class objects.
   */
  template< int codim, class Grid, class IteratorImp, class MostDerived>
  class EntityIteratorBase< codim, Grid, IteratorImp, MostDerived,
                            std::random_access_iterator_tag >
    : public EntityIteratorBase< codim, Grid, IteratorImp, MostDerived,
                                 std::forward_iterator_tag >
  {
    typedef EntityIteratorBase< codim, Grid, IteratorImp, MostDerived,
                                std::forward_iterator_tag > Base;

  protected:
    using Base::realIterator;

  public:

    //! export type of difference (for std::iterator_traits)
    typedef typename Base::difference_type difference_type;
    //! export iterator category (for std::iterator_traits)
    typedef std::random_access_iterator_tag iterator_category;

    /** \brief prefix decrement operator */
    MostDerived &operator-- ()
    {
      realIterator.advance(-1);
      return static_cast<MostDerived&>(*this);
    }

    /** \brief operator+= */
    MostDerived &operator+= (difference_type n)
    {
      realIterator.advance(n);
      return static_cast<MostDerived&>(*this);
    }

    /** \brief operator-= */
    MostDerived &operator-= (difference_type n)
    {
      realIterator.advance(-n);
      return static_cast<MostDerived&>(*this);
    }

    /** \brief operator+ */
    MostDerived operator+ (difference_type n) const
    {
      MostDerived tmp(static_cast<const MostDerived&>(*this));
      tmp += n;
      return tmp;
    }

    /** \brief operator- */
    MostDerived operator- (difference_type n) const
    {
      MostDerived tmp(static_cast<const MostDerived&>(*this));
      tmp -= n;
      return tmp;
    }

    /** \brief operator- */
    difference_type operator- (const EntityIteratorBase &other) const
    {
      return realIterator.distanceTo(other.realIterator);
    }

    // Can't really implement this: would need Entities as first class objects
    // to do a reasonable approximation, or a proxy that implements all entity
    // methods and converts to Entity&.
    //
    // /** \brief operator[] */
    // const Entity operator[] (difference_type n) const
    // {
    //   return *(*this+n);
    // }

    /** \brief operator< */
    bool operator< (const EntityIteratorBase &other) const
    {
      return realIterator.distanceTo(other.realIterator) < 0;
    }

    /** \brief operator> */
    bool operator> (const EntityIteratorBase &other) const
    {
      return realIterator.distanceTo(other.realIterator) > 0;
    }

    /** \brief operator< */
    bool operator<= (const EntityIteratorBase &other) const
    {
      return realIterator.distanceTo(other.realIterator) <= 0;
    }

    /** \brief operator> */
    bool operator>= (const EntityIteratorBase &other) const
    {
      return realIterator.distanceTo(other.realIterator) >= 0;
    }

    /** \brief default construct (undefined) iterator */
    EntityIteratorBase ( )
    {}

    /** \name Implementor's interface
     *  \{
     */

    /** \brief copy constructor from implementaton */
    EntityIteratorBase ( const IteratorImp &imp )
      : Base( imp )
    {}

    /** \} */
  };

  /** \class EntityIterator
   *  \brief interface class for an iterator over grid entities
   *  \ingroup GIEntityPointer
   *
   *  An entity iterator is an iterator over a subset of entities within a
   *  hierarchical grid. It is an extension of the Dune::EntityPointer
   *  interface.
   *
   *  Examples of entity iterators are:
   *  - iterators over the leaf level (LeafGridView::Iterator)
   *  - iterators over a grid level (LevelGridView::Iterator)
   *  - iterators over the children of an entity (Grid::HierarchicIterator)
   *  .
   *
   *  See also the documentation of Dune::EntityPointer.
   *
   *  \tparam  codim        codimension of entities this iterator walks over
   *  \tparam  Grid         type of the grid implementation
   *  \tparam  IteratorImp  type of the iterator implementation
   *
   *  \note For grid implementors: your iterator implementation must provide
   *        at least enough functionality to implement a forward iterator (see
   *        EntityIteratorBase<codim, Grid, IteratorImp, MostDerived,
   *        std::forward_iterator_tag> for the details).  If it provides
   *        enough functionality to implement a random access iterator (see
   *        EntityIteratorBase<codim, Grid, IteratorImp, MostDerived,
   *        std::random_access_iterator_tag>), you should give it a member
   *        \code typedef std::random_access_iterator_tag
   *        iterator_category;\endcode.  Setting the iterator category to any
   *        iterator tag other than std::random_access_iterator_tag should not
   *        hurt, as long as it is at least std::forward_iterator_tag.
   */
  template< int codim, class Grid, class IteratorImp>
  class EntityIterator
    : public EntityIteratorBase< codim, Grid, IteratorImp,
                                 EntityIterator<codim, Grid, IteratorImp>,
                                 typename EntityIteratorCategory<IteratorImp>::type>
  {
    typedef EntityIteratorBase< codim, Grid, IteratorImp, EntityIterator,
                                typename EntityIteratorCategory<IteratorImp>::type> Base;

  public:

    /** \brief default construct (undefined) iterator */
    EntityIterator ()
    {}

    /** \name Implementor's interface
     *  \{
     */

    /** \brief copy constructor from implementaton */
    EntityIterator ( const IteratorImp &imp )
      : Base( imp )
    {}

    /** \} */
  };

} // namespace Dune

#endif // #ifndef DUNE_GRID_ENTITYITERATOR_HH
