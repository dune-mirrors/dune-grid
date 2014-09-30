// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_GRID_TEST_CHECKRANDOMACCESSITERATORS_HH
#define DUNE_GRID_TEST_CHECKRANDOMACCESSITERATORS_HH

#include <iostream>
#include <iterator>
#include <ostream>

#include <dune/common/forloop.hh>
#include <dune/common/typetraits.hh>
#include <dune/common/test/iteratortest.hh>

#include <dune/grid/common/capabilities.hh>

// This is a copy of
// dune/common/test/iteratortest.hh:testRandomAccessIterator(), with the test
// for subscription turned into a test for *(iterator+value)
template<class Iter, class Opt>
int testRandomAccessIteratorNoSubscription(Iter begin, Iter end, Opt opt){
  int ret=testBidirectionalIterator(begin, end, opt);

  typename Iter::difference_type size = end-begin;

  srand(300);

  int no= (size>10) ? 10 : size;

  for(int i=0; i < no; i++)
  {
    int index = static_cast<int>(size*(rand()/(RAND_MAX+1.0)));
    opt(*(begin+index));
  }

  // Test the less than operator
  if(begin != end &&!( begin<end))
  {
    std::cerr<<"! (begin()<end())"<<std::endl;
    ret++;
  }

  if(begin != end) {
    if(begin-end >= 0) {
      std::cerr<<"begin!=end, but begin-end >= 0!"<<std::endl;
      ret++;
    }
    if(end-begin <= 0) {
      std::cerr<<"begin!=end, but end-begin <= 0!"<<std::endl;
      ret++;
    }
  }

  for(int i=0; i < no; i++)
  {
    int index = static_cast<int>(size*(rand()/(RAND_MAX+1.0)));
    Iter rand(begin), test(begin), res;
    rand+=index;

    if((res=begin+index) != rand)
    {
      std::cerr << " i+n should have the result i+=n, where i is the "
                <<"iterator and n is the difference type!" <<std::endl;
      ret++;
    }
    for(int i=0; i< index; i++) ++test;

    if(test != rand)
    {
      std::cerr << "i+=n should have the same result as applying the"
                << "increment ooperator n times!"<< std::cerr;
      ret++;
    }

    rand=end, test=end;
    rand-=index;


    if((end-index) != rand)
    {
      std::cerr << " i-n should have the result i-=n, where i is the "
                <<"iterator and n is the difference type!" <<std::endl;
      ret++;
    }
    for(int i=0; i< index; i++) --test;

    if(test != rand)
    {
      std::cerr << "i+=n should have the same result as applying the"
                << "increment ooperator n times!"<< std::cerr;
      ret++;
    }
  }

  for(int i=0; i < no; i++)
  {
    Iter iter1 = begin+static_cast<int>(size*(rand()/(RAND_MAX+1.0)));
    Iter iter2 = begin+static_cast<int>(size*(rand()/(RAND_MAX+1.0)));
    typename Iter::difference_type diff = iter2 -iter1;
    if((iter1+diff)!=iter2) {
      std::cerr<< "i+(j-i) = j should hold, where i,j are iterators!"<<std::endl;
      ret++;
    }
  }

  return ret;
}

struct NoOp {
  template<class T>
  void operator()(const T &)
  {}
};

// CheckCodimIterators
// -------------------

template< class GridView, int codim, class = void >
struct RACheckCodimIterators
{
  static void apply ( const GridView &gv, bool &haveRandomAccess,
                      bool &result )
  { }
};

template< class GridView, int codim >
struct RACheckCodimIterators<
  GridView, codim,
  typename Dune::enable_if<
    Dune::Capabilities::hasEntity<
      typename GridView::Grid, codim
      >::v &&
    Dune::IsBaseOf<
      std::random_access_iterator_tag,
      typename std::iterator_traits<
        typename GridView::template Codim<codim>::Iterator
        >::iterator_category
      >::value
    >::type>
{
  static void apply ( const GridView &gv, bool &haveRandomAccess,
                      bool &result )
  {
    haveRandomAccess = true;

    std::cout << "Checking random-access iterators for codim " << codim
              << std::endl;

    typedef typename GridView::template Codim<codim>::Iterator Iterator;
    typedef typename std::iterator_traits<Iterator>::difference_type
      Difference;

    Difference size = gv.size(codim);
    if(gv.template end<codim>() - gv.template begin<codim>() != size)
    {
      std::cerr << "Error: Codim " << codim << ": end() - begin() != size()."
                << std::endl;
      result = false;
    };

    if(gv.template begin<codim>() + size != gv.template end<codim>())
    {
      std::cerr << "Error: Codim " << codim << ": begin() + size() != end()."
                << std::endl;
      result = false;
    };

    if(gv.template end<codim>() - size != gv.template begin<codim>())
    {
      std::cerr << "Error: Codim " << codim << ": end() - size() != begin()."
                << std::endl;
      result = false;
    };

    {
      Iterator it = gv.template begin<codim>();
      it += size;
      if(it != gv.template end<codim>())
      {
        std::cerr << "Error: Codim " << codim << ": it = begin(); it += "
                  << "size() != end()." << std::endl;
        result = false;
      };
    }

    {
      Iterator it = gv.template end<codim>();
      it -= size;
      if(it != gv.template begin<codim>())
      {
        std::cerr << "Error: Codim " << codim << ": it = end(); it -= size() "
                  << "!= begin()." << std::endl;
        result = false;
      };
    }

    if(testRandomAccessIteratorNoSubscription(gv.template begin<codim>(),
                                              gv.template end<codim>(),
                                              NoOp()))
      result = false;
  }
};


// CheckIterators
// --------------

template< class GridView >
class RACheckIterators
{
  template< int codim >
  struct CheckCodim;

public:
  static void apply ( const GridView &gv, bool &haveRandomAccess,
                      bool &result )
  {
    std::cout << "Checking iterators for higher codimension..." << std::endl;
    Dune::ForLoop< CheckCodim, 0, GridView::dimension >::
      apply( gv, haveRandomAccess, result );
  }
};


// CheckIterators::CheckCodim
// --------------------------

template< class GridView >
template< int codim >
struct RACheckIterators< GridView >::CheckCodim
{
  static void apply ( const GridView &gv, bool &haveRandomAccess,
                      bool &result )
  {
    RACheckCodimIterators< GridView, codim >::
      apply( gv, haveRandomAccess, result );
  }
};


template< class GridView >
void checkRandomAccessIterators ( const GridView &gv )
{
  bool haveRandomAccess = false;
  bool result = true;
  RACheckIterators< GridView >::apply( gv, haveRandomAccess, result );

  if(!haveRandomAccess)
  {
    std::cerr << "Error: None of the iterators for any codim declare "
              << "themselves as random access" << std::endl;
    result = false;
  }

  if(!result)
    DUNE_THROW(Dune::Exception, "checkRandomAccessIterators() failed.");
}

#endif // DUNE_GRID_TEST_CHECKRANDOMACCESSITERATORS_HH
