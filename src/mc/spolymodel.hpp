// Copyright (C) 2009-2016 Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the Eclipse Public License.

#ifndef MC__SPOLYMODEL_H
#define MC__SPOLYMODEL_H

#include <iostream>
#include <iomanip>
#include <typeinfo>
#include <sstream>
#include <string>
#include <vector>
#include <set>
#include <map>
#include <stdarg.h>
#include <cassert>
#include <climits>
#include <limits>
#include <stdlib.h>
#include <complex>
#include <numeric>

#include "mcfunc.hpp"
#include "mclapack.hpp"
#include "mcop.hpp"

#undef  MC__SPOLYMODEL_DEBUG
#define MC__SPOLYMODEL_DEBUG_POLYBOUND
#define MC__SPOLYMODEL_CHECK

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

//! @brief C++ base class for the computation of sparse polynomial models for factorable functions: Variable
////////////////////////////////////////////////////////////////////////
//! mc::SPolyVar is a C++ base class for definition of polynomial model
//! variables in sparse format.
////////////////////////////////////////////////////////////////////////
template <typename T>
class SPolyVar
////////////////////////////////////////////////////////////////////////
{
public:

  // Monomial representation: <total order, <variable index, order>>
  typedef std::map< SPolyMon, double, lt_SPolyMon > t_coefmon;
  typedef std::set< unsigned > t_var;
  
protected:

  //!brief Unit ball in T arithmetic
  static T _TOne;

  //! @brief Set of participating variables in polynomial expression
  t_var _ndxvar;

  //! @brief Array of size <tt>_nmon</tt> with monomial coefficients of variable
  t_coefmon _coefmon;

  //! @brief Remainder bound of variable
  T _bndrem;

  //! @brief Pointer to variable bound in underlying T arithmetic (possibly NULL if not computed)
  mutable T* _bndT;

  //! @brief Pointer to polynomial bound of variable (possibly NULL if not available)
  mutable T* _bndpol;

  //! @brief Initialize private/protected members of model variable
  void _init
    ();

  //! @brief Reinitialize private/protected members of model variable
  void _reinit
    ();

  //! @brief Clean up private/protected members of variable
  void _cleanup
    ();

  //! @brief Set variable equal to <a>var</a>
  SPolyVar<T>& _set
    ( const SPolyVar<T>&var );

  //! @brief Center remainder error term <tt>_bndrem</tt>
  virtual void _center();

  //! @brief Set variable bound in unerlying T arithmetic
  virtual void _set_bndT
    ( const T&bndT );

  //! @brief Set variable bound in unerlying T arithmetic
  virtual void _set_bndT
    ( const T*bndT );

  //! @brief Unset variable bound in underlying T arithmetic
  virtual void _unset_bndT
    ();

  //! @brief Set polynomial bound in variable as <tt>bndpol</tt>
  virtual void _set_bndpol
    ( const T&bndpol );

  //! @brief Set polynomial bound in variable as <tt>bndpol</tt>
  virtual void _set_bndpol
    ( const T*bndpol );

  //! @brief Unset polynomial bound in variable
  virtual void _unset_bndpol
    ();

  //! @brief Polynomial range bounder using specified bounder <a>type</a> (pure virtual)
  virtual T _polybound
    ( const int type )
    const
    = 0;
 
  //! @brief Polynomial range bounder using default bounder (pure virtual)
  virtual T _polybound
    ()
    const
    = 0;

public:
  /** @addtogroup POLYMOD Polynomial Model Arithmetic for Factorable Functions
   *  @{
   */
  //! @brief Constructor of variable linked to polynomial model environment <a>env</a>
  SPolyVar
    ()
    { _init(); }

  //! @brief Copy constructor of variable
  SPolyVar
    ( SPolyVar<T> const& var )
    { _init(); _set( var ); }

  //! @brief Destructor of variable
  virtual ~SPolyVar()
    { delete _bndpol; delete _bndT; }

  //! @brief Set multivariate polynomial coefficients in variable as <tt>coefmon</tt>
  virtual SPolyVar<T>& set
    ( t_coefmon& coefmon )
    { _coefmon = coefmon; _unset_bndT(); _unset_bndpol();
      return *this; } // this is assuming the same order and number of variables

  //! @brief Set remainder term in variable as <tt>bndrem</tt>
  virtual SPolyVar<T>& set
    ( T const& bndrem )
    { _bndrem = bndrem; return *this; }

  //! @brief Set multivariate polynomial coefficients and remainder term equal to those in variable <tt>var</tt>, possibly defined in another polynomial model environment with fewer variables or with a different expansion order. Coefficients involving other variables or higher order are initialized to 0 if <tt>reset=true</tt> (default), otherwise they are left unmodified. Higher-order terms in TV are bounded and added to the remainder bound.
  virtual SPolyVar<T>& set
    ( SPolyVar<T> const& var, bool const reset=true );
/*
  //! @brief Copy multivariate polynomial coefficients from current variable into variable <tt>var</tt>, possibly defined in another polynomial model environment with less variables or with a lower expansion order. Copied coefficients are reset to 0 in current Taylor variable if <tt>reset=true</tt>, otherwise they are left unmodified (default).
  virtual SPolyVar<T>& get
    ( SPolyVar<T>&var, const bool reset=false );
*/
  //! @brief Compute bound on variable using bounder <a>type</a>
  virtual T bound
    ( int const type ) const
    { if( !_bndT ) return _polybound(type) + _bndrem;
      else{ T bndT; return Op<T>::inter( bndT, _polybound(type) + _bndrem, *_bndT )?
                           bndT: _polybound(type) + _bndrem; } }

  //! @brief Retreive bound on variable using default bounder
  T bound() const
    { if( !_bndT ) return bndpol() + _bndrem;
      else{ T bndT; return Op<T>::inter( bndT, bndpol() + _bndrem, *_bndT )?
                           bndT: bndpol() + _bndrem; } }

  //! @brief Retreive bound on multivariate polynomial using default bounder
  T bndpol() const
    { if( !_bndpol ) _bndpol = new T( _polybound() );
      return *_bndpol; }

  //! @brief Retreive bound on all terms with (total) order <tt>minord</tt> or higher in polynomial part
  virtual T bndord
    ( const unsigned minord )
    const
    = 0;

  //! @brief Order of polynomial model
  virtual unsigned maxord
    ()
    const
    = 0;

  //! @brief Max order of participating variables in polynomial variable
  unsigned nord
    ()
    const
    { //std::cout << "NORD = " << _coefmon.crbegin()->first.tord << std::endl; 
      return _coefmon.crbegin()->first.tord; }

  //! @brief Number of participating variables in polynomial variable
  unsigned nvar
    ()
    const
    { return _ndxvar.size(); }

  //! @brief Total number of monomial terms in polynomial variable
  unsigned nmon
    ()
    const
    { return _coefmon.size(); }

  //! @brief Return remainder term of variable
  T remainder
    ()
    const
    { return _bndrem; }

  //! @brief Shortcut to mc::SPolyVar::bound
  virtual T B
    ( const int type )
    const
    { return bound( type ); }

  //! @brief Shortcut to mc::SPolyVar::bound
  virtual T B
    ()
    const
    { return bound(); }

  //! @brief Shortcut to mc::SPolyVar::remainder
  T R
    ()
    const
    { return remainder(); }

  //! @brief Get const map of monomial coefficients
  t_coefmon const& coefmon
    ()
    const
    { return _coefmon; }

  //! @brief Get map of monomial coefficients
  t_coefmon& coefmon
    ()
    { return _coefmon; }

  //! @brief Get const set of participating variables
  t_var const& ndxvar
    ()
    const
    { return _ndxvar; }

  //! @brief Get set of participating variables
  t_var& ndxvar
    ()
    { return _ndxvar; }

  //! @brief Overloaded operator '=' for polynomial model variables
  virtual SPolyVar<T>& operator=
    ( SPolyVar<T> const& var )
    { _set( var ); return *this; }

  std::string display
    ( typename SPolyVar::t_coefmon coefmon, int const& basis )
    const;
  /** @} */
};

///////////////////////////////// SPolyVar /////////////////////////////////////

template <typename T>
inline
T SPolyVar<T>::_TOne
  = 2.*Op<T>::zeroone()-1.;

template <typename T>
inline
void
SPolyVar<T>::_init
()
{
  _bndpol = 0;
  _bndT = 0;
  _bndrem = 0.;
  return;
}

template <typename T> inline void
SPolyVar<T>::_cleanup
()
{
  delete _bndpol; _bndpol = 0;
  delete _bndT;   _bndT = 0;
  _ndxvar.clear();
  _coefmon.clear();
}

template <typename T> inline void
SPolyVar<T>::_reinit
()
{
  _cleanup();
  _init();
}

template <typename T> inline SPolyVar<T>&
SPolyVar<T>::_set
( const SPolyVar<T>&var )
{
  // Same SPolyVar?
  if( this == &var ) return *this;

  // Set coefficients and remainder
  _ndxvar  = var._ndxvar;
  _coefmon = var._coefmon;
  _bndrem  = var._bndrem;

  // Set polynomial bound
  _set_bndpol( var._bndpol );

  // Set underlying variable bound
  _set_bndT( var._bndT );

  return *this;
}

template <typename T> inline SPolyVar<T>&
SPolyVar<T>::set
( const SPolyVar<T>& var, const bool reset )
{
  if( reset ){
    _ndxvar.clear();
    _coefmon.clear();
  }
  auto&& itord = var._coefmon.upper_bound( SPolyMon( maxord()+1, std::map<unsigned,unsigned>() ) );
  for( auto&& it=var._coefmon.begin(); it!=itord; ++it ){
    for( auto&& [var,order] : it->first.expr ) _ndxvar.insert( var );
    _coefmon.insert( *it );
  }
  _bndrem  = var._bndrem + var.bndord(maxord()+1);  // var may have a higher order than *this
  _unset_bndT();
  _unset_bndpol();

  return *this;
}
/*
template <typename T> inline SPolyVar<T>&
SPolyVar<T>::get
( SPolyVar<T>& var, const bool reset )
{

  if( !_PM ){
    var._ndxmon.clear();
    if( var._PM && var._PM->sparse() ) var._ndxmon.insert(0);
    var._coefmon[0] = _coefmon[0];
    // Reset monomial coefficients to 0
    if( reset ) _coefmon[0] = 0.;
    return *this;
  }
  if( !var._PM || nvar() < var.nvar() || nord() < var.nord() ) // looks fishy...
    return *this;

  // Copy monomial coefficients from *this into var
  unsigned*iexp = new unsigned[nvar()];
  // Case: sparse representation
  for( auto it=_ndxmon.begin(); it!=_ndxmon.end() && *it<_posord(nord()+1); ++it ){
    for( unsigned ivar=0; ivar<nvar(); ivar++ )
      iexp[ivar] = ( ivar<var.nvar()? var._expmon(*it)[ivar]: 0 );
    //unsigned imon = _loc_expmon(iexp);
    auto pmon = _ndxmon.insert( _loc_expmon(iexp) );
    var._coefmon[*it] = _coefmon[*(pmon.first)];
    // Reset monomial coefficients to 0
    if( reset ) _coefmon[*(pmon.first)] = 0.;
  }
  // Case: dense representation
  for( unsigned jmon=0; _ndxmon.empty() && jmon<var.nmon(); jmon++ ){
    for( unsigned ivar=0; ivar<nvar(); ivar++ )
      iexp[ivar] = ( ivar<var.nvar()? var._expmon(jmon)[ivar]: 0 );
    var._coefmon[jmon] = _coefmon[_loc_expmon(iexp)];
    // Reset monomial coefficients to 0
    if( reset ) _coefmon[_loc_expmon(iexp)] = 0.;
  }
  delete[] iexp;
  *var._bndrem = 0.;
  var._unset_bndT();
  var._unset_bndpol();
  var._bndord_uptd = false;

  if( reset ){
    _unset_bndT();
    _unset_bndpol();
    _bndord_uptd = false;
  }

  return *this;
}
*/
template <typename T> inline void
SPolyVar<T>::_unset_bndpol
()
{
  delete _bndpol;
  _bndpol = 0;
}

template <typename T> inline void
SPolyVar<T>::_set_bndpol
( const T*bndpol )
{
  if( !bndpol ){
    if( _bndpol ) delete _bndpol;
    _bndpol = 0;
  }
  else if( !_bndpol )
    _bndpol = new T( *bndpol );
  else
    *_bndpol = *bndpol;
}

template <typename T> inline void
SPolyVar<T>::_set_bndpol
( const T&bndpol )
{
  if( !_bndpol )
    _bndpol = new T( bndpol );
  else
    *_bndpol = bndpol;
}

template <typename T> inline void
SPolyVar<T>::_unset_bndT
()
{
  delete _bndT;
  _bndT = 0;
}

template <typename T> inline void
SPolyVar<T>::_set_bndT
( const T*bndT )
{
  if( !bndT ){
    if( _bndT ) delete _bndT;
    _bndT = 0;
  }
  else if( !_bndT )
    _bndT = new T( *bndT );
  else
    *_bndT = *bndT;
}

template <typename T> inline void
SPolyVar<T>::_set_bndT
( const T&bndT )
{
  if( !_bndT )
    _bndT = new T( bndT );
  else
    *_bndT = bndT;
}

template <typename T> inline void
SPolyVar<T>::_center()
{
  const double remmid = Op<T>::mid(_bndrem);
  if( remmid == 0. ) return;
  if( _coefmon.empty() || _coefmon.begin()->first.tord ) 
    _coefmon.insert( std::make_pair( SPolyMon(), remmid ) );
  else
    _coefmon.begin()->second += remmid;
  _bndrem -= remmid;
  if( _bndpol ) *_bndpol += remmid;
}

template <typename T> inline std::string
SPolyVar<T>::display
( typename SPolyVar::t_coefmon coefmon, int const& basis )
const
{
  unsigned IDISP = 5;
  std::ostringstream out;
  out << std::endl << std::scientific << std::setprecision(IDISP)
      << std::right;

  // Sparse multivariate polynomial
  for( auto const& [mon,coef] : coefmon )
    out << std::right << std::setw(IDISP+7) << coef << "  " << std::setw(2) << mon.tord << "  "// << "   "
        << mon.display( basis ) << std::endl;

  return out.str();
}

} // namespace mc

#endif

