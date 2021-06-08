// Copyright (C) 2021 Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the Eclipse Public License.

#ifndef MC__SPOLYMON_H
#define MC__SPOLYMON_H

#include <iostream>
#include <iomanip>
//#include <typeinfo>
#include <sstream>
#include <string>
//#include <vector>
//#include <set>
#include <map>
//#include <stdarg.h>
#include <cassert>
//#include <climits>
//#include <limits>
//#include <stdlib.h>
//#include <complex>
//#include <numeric>

//#include "mcfunc.hpp"
//#include "mclapack.hpp"
//#include "mcop.hpp"

namespace mc
{

//! @brief C++ class for sparse monomial storage and manipulation
////////////////////////////////////////////////////////////////////////
//! mc::SPolyMon is a C++ class for sparse monomial storage and
//! manipulation.
////////////////////////////////////////////////////////////////////////
struct SPolyMon
////////////////////////////////////////////////////////////////////////
{
  //! @brief Monomial total order
  unsigned tord;

  //! @brief Monomial variables and partial orders 
  std::map< unsigned, unsigned > expr;

  //! @brief Constructor of constant monomial
  SPolyMon
    ()
    : tord( 0 )
    {}

  //! @brief Constructor of monomial for variable with index <a>i</a>
  SPolyMon
    ( unsigned const i )
    : tord( 1 )
    { expr.insert( std::make_pair( i, 1 ) ); }

  //! @brief Copy constructor of monomial
  SPolyMon
    ( unsigned const tord_, std::map< unsigned, unsigned > const& expr_ )
    : tord( tord_ ), expr( expr_ )
    {}
    
  //! @brief Append monomial to output stream <a>out</a>
  std::string display
    ( int const& basis )
    const;

  //! @brief Greatest common exponent of terms in monomial
  unsigned int gcexp
    ()
    const;

  //! @brief Least exponent among all terms in monomial
  unsigned int lexp
    ()
    const;

  //! @brief Greatest exponent among all terms in monomial
  unsigned int gexp
    ()
    const;

  //! @brief Test for intersection
  bool inter
    ( SPolyMon const& Mon )
    const;

  //! @brief Test for proper subset
  bool subset
    ( SPolyMon const& Mon )
    const;
    
  //! @brief Test for subset
  bool subseteq
    ( SPolyMon const& Mon )
    const;
    
  //! @brief Overloaded operator '+=' for monomial
  SPolyMon& operator+=
    ( SPolyMon const& mon );

  //! @brief Overloaded operator '-=' for monomial
  SPolyMon& operator-=
    ( SPolyMon const& mon );

  //! @brief Overloaded operator '*=' for monomial
  SPolyMon& operator*=
    ( unsigned const& factor );

  //! @brief Overloaded operator '/=' for monomial
  SPolyMon& operator/=
    ( unsigned const& factor );

  //! @brief Overloaded operator '[]' for extracting from monomial
  SPolyMon operator[]
    ( unsigned const& ivar )
    const;

  //! @brief Exceptions of mc::SPolyMon
  class Exceptions
  {
   public:
    //! @brief Enumeration type for mc::Interval exceptions
    enum TYPE{
      SUB=1,    //!< Subtraction of a monomial that is not a proper subset
      DIV,	    //!< Division by a factor greater than the greatest common exponent
    };
    //! @brief Constructor for error flag <a>ierr</a>
    Exceptions( TYPE ierr ) : _ierr( ierr ){}
    //! @brief Return error flag
    int ierr(){ return _ierr; }
    //! @brief Return error description
    std::string what(){
      switch( _ierr ){
      case SUB:
        return "mc::SPolyMon\t Subtraction of a monomial that is not a proper subset";
      case DIV:
        return "mc::SPolyMon\t Division by a factor greater than the greatest common exponent";
      }
      return "mc::SPolyMon\t Undocumented error";
    }
   private:
    TYPE _ierr;
  };
};

//inline std::ostream&
//operator<<
//( std::ostream& out, SPolyMon const& mon )
//{
//  // Sparse multivariate polynomial
//  if( mon.expr.empty() )  out << "1";
//  for( auto&& ie=mon.expr.begin(); ie!=mon.expr.end(); ++ie ){
//    if( ie != mon.expr.begin() ) out << "·";
//    out << "T" << ie->second << "[" << ie->first << "]";
//  }
//  return out;
//}

inline std::string
SPolyMon::display
( int const& basis )
const
{
  std::ostringstream out;
  // Sparse monomial
  if( expr.empty() )  out << "1";
  for( auto ie=expr.begin(); ie!=expr.end(); ++ie ){
    if( ie != expr.begin() ) out << "·";
    switch( basis ){
     case 0:
      out << "[" << ie->first << "]";
      if( ie->second > 1 ) out << "^" << ie->second;
      break;
     default:
      out << "T" << ie->second << "[" << ie->first << "]";
      break;
    }
  }
  return out.str();
}

inline SPolyMon
SPolyMon::operator[]
( unsigned const& ivar )
const
{
  auto itvar = expr.find( ivar );
  if( itvar == expr.end() ) return SPolyMon();
  return SPolyMon( itvar->second, {*itvar} );
}

inline unsigned int
SPolyMon::gcexp
()
const
{
  if( expr.size() <= 1 ) return tord;
  auto&& it = expr.cbegin();
  unsigned gce = it->second;
  for( ++it; gce>1 && it!=expr.cend(); ++it )
    gce = std::gcd( it->second, gce );
  return gce;
}

inline unsigned int
SPolyMon::gexp
()
const
{
  if( expr.size() <= 1 ) return tord;
  auto&& it = expr.cbegin();
  unsigned ge = it->second;
  for( ++it; it!=expr.cend(); ++it )
    if( ge < it->second ) ge = it->second;
  return ge;
}

inline unsigned int
SPolyMon::lexp
()
const
{
  if( expr.size() <= 1 ) return tord;
  auto&& it = expr.cbegin();
  unsigned le = it->second;
  for( ++it; it!=expr.cend(); ++it )
    if( le > it->second ) le = it->second;
  return le;
}

inline bool
SPolyMon::inter
( SPolyMon const& Mon )
const
{
  for( auto&& it1=expr.begin(); it1!=expr.end(); ++it1 )
    // Monomials share a common variable (regardless of exponent)
    if( Mon.expr.find( it1->first ) != Mon.expr.end() ) return true;
  return false;
}

inline bool
SPolyMon::subset
( SPolyMon const& Mon )
const
{
  // Total order must be larger for Mon than *this
  if( tord >= Mon.tord ) return false;
  for( auto&& it1=expr.begin(); it1!=expr.end(); ++it1 ){
    auto&& it2 = Mon.expr.find( it1->first );
    // Mon must comprise (at least) the same variables as *this, and the order
    // of each variable in Mon must have a larger order than in *this
    if( it2 == Mon.expr.end() || it1->second > it2->second ) return false;
  }
  return true;
}

inline bool
SPolyMon::subseteq
( SPolyMon const& Mon )
const
{
  // Total order must be larger or equal for Mon than *this
  if( tord > Mon.tord ) return false;
  for( auto&& it1=expr.begin(); it1!=expr.end(); ++it1 ){
    auto&& it2 = Mon.expr.find( it1->first );
    // Mon must comprise (at least) the same variables as *this, and the order
    // of each variable in Mon must have a larger order than in *this
    if( it2 == Mon.expr.end() || it1->second > it2->second ) return false;
  }
  return true;
}

inline SPolyMon&
SPolyMon::operator*=
( unsigned const& factor )
{
  // Return if factor is unity
  if( factor == 1 )
    return *this;

  // divide monomial partial orders by factor
  tord = 0;
  for( auto& varpow : expr ){
    varpow.second *= factor;
    tord += varpow.second;
  }
  return *this;
}

inline SPolyMon
operator*
( SPolyMon const& Mon, unsigned const& factor )
{
  SPolyMon Mon2( Mon );
  return( Mon2 *= factor );
}

inline SPolyMon&
SPolyMon::operator/=
( unsigned const& factor )
{
  // Return if factor is unity
  if( factor == 1 )
    return *this;

  // factor may not be greater than the least exponent
  if( lexp() < factor )
    throw typename SPolyMon::Exceptions( SPolyMon::Exceptions::DIV );

  // divide monomial partial orders by factor
  tord = 0;
  for( auto& varpow : expr ){
    varpow.second /= factor;
    tord += varpow.second;
  }
  return *this;
}

inline SPolyMon
operator/
( SPolyMon const& Mon, unsigned const& factor )
{
  SPolyMon Mon2( Mon );
  return( Mon2 /= factor );
}

inline SPolyMon&
SPolyMon::operator+=
( SPolyMon const& mon )
{
  for( auto const& varpow : mon.expr ){
    auto&& it = expr.insert( varpow );
    // If element from mon was not inserted, increment existing variable order
    if( !it.second ) it.first->second += varpow.second;
  }
  tord += mon.tord;
  return *this;
}

inline SPolyMon
operator+
( SPolyMon const& Mon1, SPolyMon const& Mon2 )
{
  // Create a copy of Mon1 and add terms with Mon2 
  SPolyMon Mon3( Mon1 );
  for( auto&& it2=Mon2.expr.begin(); it2!=Mon2.expr.end(); ++it2 ){
    auto&& it3 = Mon3.expr.insert( *it2 );
    // If element from Mon2 was not inserted, increment existing variable order
    if( !it3.second ) it3.first->second += it2->second;
  }
  Mon3.tord += Mon2.tord;
  return Mon3;
}

inline SPolyMon&
SPolyMon::operator-=
( SPolyMon const& mon )
{
  // mon must be a proper subset
  if( !mon.subseteq( *this ) )
    throw Exceptions( Exceptions::SUB );

  // Return constant monomial if mon is identical
  if( mon.tord == tord ){
    expr.clear();
    tord = 0;
    return *this;
  }

  // Cancel the common terms with mon
  for( auto& varpow : mon.expr ){
    auto&& it = expr.find( varpow.first );
    assert( it != mon.expr.end() );
    if( it->second == varpow.second )
      expr.erase( it );
    else
      it->second -= varpow.second;
  }
  tord -= mon.tord;
  return *this;
}

inline SPolyMon
operator-
( SPolyMon const& Mon1, SPolyMon const& Mon2 )
{
  // Mon2 must be a proper subset of Mon1 
  if( !Mon2.subseteq( Mon1 ) )
    throw typename SPolyMon::Exceptions( SPolyMon::Exceptions::SUB );

  // Return constant monomial if Mon1 and Mon2 are identical 
  if( Mon1.tord == Mon2.tord )
    return SPolyMon();

  // Create a copy of Mon1 and cancel the common terms with Mon2 
  SPolyMon Mon3( Mon1 );
  for( auto&& it2=Mon2.expr.begin(); it2!=Mon2.expr.end(); ++it2 ){
    auto&& it3 = Mon3.expr.find( it2->first );
    assert( it3 != Mon3.expr.end() );
    if( it3->second == it2->second )
      Mon3.expr.erase( it3 );
    else
      it3->second -= it2->second;
  }
  Mon3.tord -= Mon2.tord;
  return Mon3;
}

//! @brief C++ structure for ordering of monomials in graded lexicographic order (grlex)
struct lt_SPolyMon
{
  bool operator()
    ( SPolyMon const& Mon1, SPolyMon const& Mon2 ) const
    {
      // Order monomials based on their total order first
      if( Mon1.tord < Mon2.tord ) return true;
      if( Mon1.tord > Mon2.tord ) return false;
      // Account for the case of an empty list
      if( Mon2.expr.empty() ) return false;
      if( Mon1.expr.empty() ) return true;
      // Order in graded lexicographic order next
      for( auto&& it1=Mon1.expr.begin(),&& it2=Mon2.expr.begin(); it1!=Mon1.expr.end(); ++it1, ++it2 ){
        if( it1->first < it2->first ) return true;
        if( it2->first < it1->first ) return false;
        if( it1->second > it2->second ) return true;
        if( it1->second < it2->second ) return false;
      }
      return false;
    }
};

//! @brief C++ structure for ordering of monomials in graded lexicographic order (grlex)
struct lt_pSPolyMon
{
  bool operator()
    ( SPolyMon const* Mon1, SPolyMon const* Mon2 ) const
    {
      assert( Mon1 && Mon2 );
      return lt_SPolyMon()( *Mon1, *Mon2 );
    }
};

} // namespace mc

#endif

