// Copyright (C) 2022 Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the Eclipse Public License.

#ifndef MC__SMON_H
#define MC__SMON_H

#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <map>
#include <cassert>
#include <functional>
#include <numeric>

namespace mc
{

//! @brief C++ class for sparse monomial storage and manipulation
////////////////////////////////////////////////////////////////////////
//! mc::SMon is a C++ class for sparse monomial storage and
//! manipulation.
////////////////////////////////////////////////////////////////////////
template <typename KEY=unsigned, typename COMP=std::less<unsigned>>
struct SMon
////////////////////////////////////////////////////////////////////////
{
  //! @brief Monomial total order
  unsigned tord;

  //! @brief Monomial variables and partial orders 
  std::map<KEY,unsigned,COMP> expr;

  //! @brief Constructor of constant monomial
  SMon
    ()
    : tord( 0 )
    {}

  //! @brief Constructor of monomial for variable with index <a>i</a>
  SMon
    ( KEY const& var, const unsigned ord=1 )
    : tord( ord )
    { if( tord ) expr.insert( std::make_pair( var, ord ) ); }

  //! @brief Copy constructor of monomial
  SMon
    ( unsigned const tord_, std::map<KEY,unsigned,COMP> const& expr_ )
    : tord( tord_ ), expr( expr_ )
    {}

  //! @brief Copy constructor of monomial with conversion
  template <typename OTHERKEY, typename OTHERCOMP>
  SMon
    ( unsigned const tord_, std::map<OTHERKEY,unsigned,OTHERCOMP> const& expr_,
      std::map<OTHERKEY,KEY,OTHERCOMP>& match )
    : tord( tord_ )
    { for( auto const& [key,ord] : expr_ ) expr[match[key]] = ord; }

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
    ( SMon<KEY,COMP> const& Mon )
    const;

  //! @brief Test for proper subset
  bool subset
    ( SMon<KEY,COMP> const& Mon )
    const;
    
  //! @brief Test for subset
  bool subseteq
    ( SMon<KEY,COMP> const& Mon )
    const;
    
  //! @brief Overloaded operator '+=' for monomial
  SMon<KEY,COMP>& operator+=
    ( SMon<KEY,COMP> const& mon );

  //! @brief Overloaded operator '-=' for monomial
  SMon<KEY,COMP>& operator-=
    ( SMon<KEY,COMP> const& mon );

  //! @brief Overloaded operator '*=' for monomial
  SMon<KEY,COMP>& operator*=
    ( unsigned const& factor );

  //! @brief Overloaded operator '/=' for monomial
  SMon<KEY,COMP>& operator/=
    ( unsigned const& factor );

  //! @brief Overloaded operator '[]' for extracting from monomial
  SMon<KEY,COMP> operator[]
    ( KEY const& var )
    const;

  //! @brief Exceptions of mc::SMon
  class Exceptions
  {
   public:
    //! @brief Enumeration type for mc::Interval exceptions
    enum TYPE{
      SUB=1,   //!< Subtraction of a monomial that is not a proper subset
      DIV,	//!< Division by a factor greater than the greatest common exponent
      CONV	//!< Conversion with different variable indexing failed
    };
    //! @brief Constructor for error flag <a>ierr</a>
    Exceptions( TYPE ierr ) : _ierr( ierr ){}
    //! @brief Return error flag
    int ierr(){ return _ierr; }
    //! @brief Return error description
    std::string what(){
      switch( _ierr ){
      case SUB:
        return "mc::SMon\t Subtraction of a monomial that is not a proper subset";
      case DIV:
        return "mc::SMon\t Division by a factor greater than the greatest common exponent";
      case CONV:
        return "mc::SMon\t Conversion with different variable indexing failed";
      }
      return "mc::SMon\t Undocumented error";
    }
   private:
    TYPE _ierr;
  };

  //! @brief Overloads turning either a reference or a pointer into a pointer
  template<typename T>
  static T* ptr
   ( T& obj )
   { return &obj; }
  template<typename T>
  static T* ptr
   ( T* obj )
   { return obj; }
};

template <typename KEY, typename COMP>
inline std::string
SMon<KEY,COMP>::display
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
      out << "[" << *ptr(ie->first) << "]";
      if( ie->second > 1 ) out << "^" << ie->second;
      break;
     default:
      out << "T" << ie->second << "[" << *ptr(ie->first) << "]";
      break;
    }
  }
  return out.str();
}

template <typename KEY, typename COMP>
inline SMon<KEY,COMP>
SMon<KEY,COMP>::operator[]
( KEY const& ivar )
const
{
  auto itvar = expr.find( ivar );
  if( itvar == expr.end() ) return SMon<KEY,COMP>();
  return SMon<KEY,COMP>( itvar->second, {*itvar} );
}

template <typename KEY, typename COMP>
inline unsigned int
SMon<KEY,COMP>::gcexp
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

template <typename KEY, typename COMP>
inline unsigned int
SMon<KEY,COMP>::gexp
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

template <typename KEY, typename COMP>
inline unsigned int
SMon<KEY,COMP>::lexp
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

template <typename KEY, typename COMP>
inline bool
SMon<KEY,COMP>::inter
( SMon<KEY,COMP> const& Mon )
const
{
  for( auto&& it1=expr.begin(); it1!=expr.end(); ++it1 )
    // Monomials share a common variable (regardless of exponent)
    if( Mon.expr.find( it1->first ) != Mon.expr.end() ) return true;
  return false;
}

template <typename KEY, typename COMP>
inline bool
SMon<KEY,COMP>::subset
( SMon<KEY,COMP> const& Mon )
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

template <typename KEY, typename COMP>
inline bool
SMon<KEY,COMP>::subseteq
( SMon<KEY,COMP> const& Mon )
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

template <typename KEY, typename COMP>
inline SMon<KEY,COMP>&
SMon<KEY,COMP>::operator*=
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

template <typename KEY, typename COMP>
inline SMon<KEY,COMP>
operator*
( SMon<KEY,COMP> const& Mon, unsigned const& factor )
{
  SMon<KEY,COMP> Mon2( Mon );
  return( Mon2 *= factor );
}

template <typename KEY, typename COMP>
inline SMon<KEY,COMP>&
SMon<KEY,COMP>::operator/=
( unsigned const& factor )
{
  // Return if factor is unity
  if( factor == 1 )
    return *this;

  // factor may not be greater than the least exponent
  if( lexp() < factor )
    throw typename SMon<KEY,COMP>::Exceptions( SMon<KEY,COMP>::Exceptions::DIV );

  // divide monomial partial orders by factor
  tord = 0;
  for( auto& varpow : expr ){
    varpow.second /= factor;
    tord += varpow.second;
  }
  return *this;
}

template <typename KEY, typename COMP>
inline SMon<KEY,COMP>
operator/
( SMon<KEY,COMP> const& Mon, unsigned const& factor )
{
  SMon<KEY,COMP> Mon2( Mon );
  return( Mon2 /= factor );
}

template <typename KEY, typename COMP>
inline SMon<KEY,COMP>&
SMon<KEY,COMP>::operator+=
( SMon<KEY,COMP> const& mon )
{
  for( auto const& varpow : mon.expr ){
    auto&& it = expr.insert( varpow );
    // If element from mon was not inserted, increment existing variable order
    if( !it.second ) it.first->second += varpow.second;
  }
  tord += mon.tord;
  return *this;
}

template <typename KEY, typename COMP>
inline SMon<KEY,COMP>
operator+
( SMon<KEY,COMP> const& Mon1, SMon<KEY,COMP> const& Mon2 )
{
  // Create a copy of Mon1 and add terms with Mon2 
  SMon<KEY,COMP> Mon3( Mon1 );
  return( Mon3 += Mon2 );
//  for( auto&& it2=Mon2.expr.begin(); it2!=Mon2.expr.end(); ++it2 ){
//    auto&& it3 = Mon3.expr.insert( *it2 );
//    // If element from Mon2 was not inserted, increment existing variable order
//    if( !it3.second ) it3.first->second += it2->second;
//  }
//  Mon3.tord += Mon2.tord;
//  return Mon3;
}

template <typename KEY, typename COMP>
inline SMon<KEY,COMP>&
SMon<KEY,COMP>::operator-=
( SMon<KEY,COMP> const& mon )
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

template <typename KEY, typename COMP>
inline SMon<KEY,COMP>
operator-
( SMon<KEY,COMP> const& Mon1, SMon<KEY,COMP> const& Mon2 )
{
  // Create a copy of Mon1 and add terms with Mon2 
  SMon<KEY,COMP> Mon3( Mon1 );
  return( Mon3 -= Mon2 );
//  // Mon2 must be a proper subset of Mon1 
//  if( !Mon2.subseteq( Mon1 ) )
//    throw typename SMon<KEY,COMP>::Exceptions( SMon<KEY,COMP>::Exceptions::SUB );

//  // Return constant monomial if Mon1 and Mon2 are identical 
//  if( Mon1.tord == Mon2.tord )
//    return SMon();

//  // Create a copy of Mon1 and cancel the common terms with Mon2 
//  SMon Mon3( Mon1 );
//  for( auto&& it2=Mon2.expr.begin(); it2!=Mon2.expr.end(); ++it2 ){
//    auto&& it3 = Mon3.expr.find( it2->first );
//    assert( it3 != Mon3.expr.end() );
//    if( it3->second == it2->second )
//      Mon3.expr.erase( it3 );
//    else
//      it3->second -= it2->second;
//  }
//  Mon3.tord -= Mon2.tord;
//  return Mon3;
}

//! @brief C++ structure for ordering of monomials in graded lexicographic order (grlex)
template <typename COMP=std::less<unsigned>>
struct lt_SMon
{
  template <typename KEY>
  bool operator()
    ( SMon<KEY,COMP> const& Mon1, SMon<KEY,COMP> const& Mon2 )
    const
    {
      // Order monomials based on their total order first
      if( Mon1.tord < Mon2.tord ) return true;
      if( Mon1.tord > Mon2.tord ) return false;
      // Account for the case of an empty list
      if( Mon2.expr.empty() ) return false;
      if( Mon1.expr.empty() ) return true;
      // Order in graded lexicographic order next
      for( auto&& it1=Mon1.expr.begin(),&& it2=Mon2.expr.begin(); it1!=Mon1.expr.end(); ++it1, ++it2 ){
        if( COMP()( it1->first, it2->first ) ) return true;
        if( COMP()( it2->first, it1->first ) ) return false;
        if( it1->second > it2->second ) return true;
        if( it1->second < it2->second ) return false;
      }
      return false;
    }
};

//! @brief C++ structure for ordering of monomials in graded lexicographic order (grlex)
template <typename COMP=std::less<unsigned>>
struct lt_pSMon
{
  template <typename KEY>
  bool operator()
    ( SMon<KEY,COMP> const* Mon1, SMon<KEY,COMP> const* Mon2 )
    const
    {
      assert( Mon1 && Mon2 );
      return lt_SMon<COMP>()( *Mon1, *Mon2 );
    }
};

} // namespace mc

#endif

