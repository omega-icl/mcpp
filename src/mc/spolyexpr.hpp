// Copyright (C) 2017 Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the Eclipse Public License.

#ifndef MC__SPOLYEXPR_H
#define MC__SPOLYEXPR_H

#include <set>
#include <map>
#include <numeric>

#include "ffunc.hpp"

#undef  MC__SPOLYEXPR_DEBUG_SPROD

namespace mc
{
class SPolyExpr;

//! @brief C++ structure for defining monomials in multivariate polynomials
struct FFMon
{
  //! @brief Monomial total order
  unsigned tord;

  //! @brief Monomial variables and partial orders 
  std::map< FFVar const*, unsigned > expr;

  //! @brief Constructor of constant monomial
  FFMon
    ()
    : tord( 0 )
    {}

  //! @brief Constructor of monomial for variable <a>x</a>
  FFMon
    ( FFVar const& x )
    : tord( 1 )
    { expr.insert( std::make_pair( &x, 1 ) ); }

  //! @brief Copy constructor of monomial
  FFMon
    ( unsigned const tord_, std::map< FFVar const*, unsigned > const& expr_ )
    : tord( tord_ ), expr( expr_ )
    {}

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

  //! @brief Overloaded operator '+=' for monomial
  FFMon& operator+=
    ( FFMon const& mon );

  //! @brief Overloaded operator '-=' for monomial
  FFMon& operator-=
    ( FFMon const& mon );

  //! @brief Overloaded operator '/=' for monomial
  FFMon& operator/=
    ( unsigned const& factor );
    
  //! @brief Exceptions of mc::FFMon
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
        return "mc::FFMon\t Subtraction of a monomial that is not a proper subset";
      case DIV:
        return "mc::FFMon\t Division by a factor greater than the greatest common exponent";
      }
      return "mc::FFMon\t Undocumented error";
    }
   private:
    TYPE _ierr;
  };
};

inline unsigned int
FFMon::gcexp
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
FFMon::gexp
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
FFMon::lexp
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
operator<
( FFMon const& Mon1, FFMon const& Mon2 )
{
  // Total order must be larger for Mon2
  if( Mon1.tord >= Mon2.tord ) return false;
  for( auto it1=Mon1.expr.begin(); it1!=Mon1.expr.end(); ++it1 ){
    auto it2 = Mon2.expr.find( it1->first );
    // Mon2 must comprise (at least) the same variables as Mon1, and the order
    // of each variable in Mon2 must have a larger order than in Mon1
    if( it2 == Mon2.expr.end() || it1->second > it2->second ) return false;
  }
  return true;
}

inline bool
operator<=
( FFMon const& Mon1, FFMon const& Mon2 )
{
  // Total order must be larger or equal for Mon2
  if( Mon1.tord > Mon2.tord ) return false;
  for( auto it1=Mon1.expr.begin(); it1!=Mon1.expr.end(); ++it1 ){
    auto it2 = Mon2.expr.find( it1->first );
    // Mon2 must comprise (at least) the same variables as Mon1, and the order
    // of each variable in Mon2 must have a larger order than in Mon1
    if( it2 == Mon2.expr.end() || it1->second > it2->second ) return false;
  }
  return true;
}

inline FFMon&
FFMon::operator/=
( unsigned const& factor )
{
  // Return if factor is unity
  if( factor == 1 )
    return *this;

  // factor may not be greater than the least exponent
  if( lexp() < factor )
    throw typename FFMon::Exceptions( FFMon::Exceptions::DIV );

  // divide monomial partial orders by factor
  tord = 0;
  for( auto&& varpow : expr ){
    varpow.second /= factor;
    tord += varpow.second;
  }
  return *this;
}

inline FFMon
operator/
( FFMon const& Mon, unsigned const& factor )
{
  FFMon Mon2( Mon );
  return( Mon2 /= factor );
}

inline FFMon&
FFMon::operator+=
( FFMon const& mon )
{
  for( auto&& varpow : mon.expr ){
    auto it = expr.insert( varpow );
    // If element from mon was not inserted, increment existing variable order
    if( !it.second ) it.first->second += varpow.second;
  }
  tord += mon.tord;
  return *this;
}

inline FFMon
operator+
( FFMon const& Mon1, FFMon const& Mon2 )
{
  // Create a copy of Mon1 and add terms with Mon2 
  FFMon Mon3( Mon1 );
  for( auto it2=Mon2.expr.begin(); it2!=Mon2.expr.end(); ++it2 ){
    auto it3 = Mon3.expr.insert( *it2 );
    // If element from Mon2 was not inserted, increment existing variable order
    if( !it3.second ) it3.first->second += it2->second;
  }
  Mon3.tord += Mon2.tord;
  return Mon3;
}

inline FFMon&
FFMon::operator-=
( FFMon const& mon )
{
  // mon must be a proper subset
  if( !( mon <= *this ) )
    throw Exceptions( Exceptions::SUB );

  // Return constant monomial if mon is identical
  if( mon.tord == tord ){
    expr.clear();
    tord = 0;
    return *this;
  }

  // Cancel the common terms with mon
  for( auto&& varpow : mon.expr ){
    auto it = expr.find( varpow.first );
    assert( it != mon.expr.end() );
    if( it->second == varpow.second )
      expr.erase( it );
    else
      it->second -= varpow.second;
  }
  tord -= mon.tord;
  return *this;
}

inline FFMon
operator-
( FFMon const& Mon1, FFMon const& Mon2 )
{
  // Mon2 must be a proper subset of Mon1 
  if( !( Mon2 <= Mon1 ) )
    throw typename FFMon::Exceptions( FFMon::Exceptions::SUB );

  // Return constant monomial if Mon1 and Mon2 are identical 
  if( Mon1.tord == Mon2.tord )
    return FFMon();

  // Create a copy of Mon1 and cancel the common terms with Mon2 
  FFMon Mon3( Mon1 );
  for( auto it2=Mon2.expr.begin(); it2!=Mon2.expr.end(); ++it2 ){
    auto it3 = Mon3.expr.find( it2->first );
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
struct lt_FFMon
{
  bool operator()
    ( FFMon const& Mon1, FFMon const& Mon2 ) const
    {
      // Order monomials based on their total order first
      if( Mon1.tord < Mon2.tord ) return true;
      if( Mon1.tord > Mon2.tord ) return false;
      // Order in graded lexicographic order next
      for( auto it1=Mon1.expr.begin(), it2=Mon2.expr.begin(); it1!=Mon1.expr.end(); ++it1, ++it2 ){
        if( lt_FFVar()( it1->first, it2->first ) ) return true;
        if( lt_FFVar()( it2->first, it1->first ) ) return false;
        if( it1->second > it2->second ) return true;
        if( it1->second < it2->second ) return false;
      }
      return false;
    }
};

//! @brief C++ structure for ordering of monomials in graded lexicographic order (grlex)
struct lt_pFFMon
{
  bool operator()
    ( FFMon const* Mon1, FFMon const* Mon2 ) const
    {
      assert( Mon1 && Mon2 );
      return lt_FFMon()( *Mon1, *Mon2 );
    }
};

//! @brief C++ class for sparse polynomial representation and arithmetic
////////////////////////////////////////////////////////////////////////
//! mc::SPolyExpr is a C++ class for sparse polynomial representation
//! and arithmetic
////////////////////////////////////////////////////////////////////////
class SPolyExpr
////////////////////////////////////////////////////////////////////////
{
public:

  // Monomial representation: FFMon := <total order, <DAG variable ptr, order>>
  typedef std::map< FFMon, double, lt_FFMon > t_ffpoly;
  typedef std::set< FFVar const*, lt_FFVar > t_ffvar;

  // Friends for arithmetic operations
  friend SPolyExpr sqr ( SPolyExpr const& );
  friend SPolyExpr operator- ( SPolyExpr const& );
  friend std::ostream& operator<< ( std::ostream &, SPolyExpr const& );

protected:

  //! @brief Map of monomial terms in polynomial expression
  t_ffpoly _mapmon;

  //! @brief Set of participating variables in polynomial expression
  t_ffvar _setvar;

  //! @brief Initialize private/protected members of model variable
  void _init
    ();

  //! @brief Reinitialize private/protected members of model variable
  void _reinit
    ();

  //! @brief Clean up private/protected members of variable
  void _cleanup
    ();

  //! @brief Set polynomial expression equal to <a>spoly</a>
  SPolyExpr& _set
    ( SPolyExpr const& spoly );

  //! @brief Set polynomial expression equal to constant <a>d</a>
  SPolyExpr& _set
    ( double const d );

  //! @brief Set polynomial expression equal to variable <a>x</a>
  SPolyExpr& _set
    ( FFVar const& x );

  //! @brief Set polynomial expression equal to monomial <a>mon</a>
  SPolyExpr& _set
    ( std::pair< FFMon, double > const& mon );

  //! @brief Clean sparse polynomial by removing zero entries
  void _clean
    ( t_ffpoly & spoly )
    const;

  //! @brief Build univariate sparse polynomial (with sparse polynomial coefficients)
  void _svec1D
    ( typename t_ffvar::const_iterator itvar, std::pair<FFMon,double> const& mon,
      std::map<unsigned,t_ffpoly> & mapspoly )
    const;

  //! @brief Recursive product of univariate sparse polynomials
  void _sprod1D
    ( const std::map<unsigned,t_ffpoly>& sp1map,
      const std::map<unsigned,t_ffpoly>& sp2map,
      t_ffpoly& mapmon, typename t_ffvar::const_iterator itvar )
    const;

  //! @brief Scaling of univariate sparse polynomial
  void _sscal1D
    ( const t_ffpoly& spoly, const double& dscal, t_ffpoly& spscal )
    const;

  //! @brief Lifting of sparse polynomial coefficient maps
  void _slift1D
    ( const t_ffpoly& spoly, const double& dscal, t_ffpoly& splift )
    const;

  //! @brief Lifting of sparse polynomial coefficient maps
  void _slift1D
    ( const t_ffpoly& spoly, const double& dscal, t_ffpoly& splift,
      typename t_ffvar::const_iterator itvar, const unsigned iord )
    const;

  //! @brief Display of recursive univariate sparse polynomial
  void _sdisp1D
    ( const t_ffpoly&spoly, const std::string&name, std::ostream&osos=std::cout )
    const;

  //! @brief Display of recursive univariate sparse polynomial
  void _sdisp1D
    ( const std::map<unsigned,t_ffpoly>&mapspoly, typename t_ffvar::const_iterator itvar, 
      const std::string&name, std::ostream&osos=std::cout )
    const;

public:

  //! @brief Default constructor of sparse polynomial expression
  SPolyExpr
    ()
    { _init(); }

  //! @brief Constructor of sparse polynomial expression as constant
  SPolyExpr
    ( double const d )
    { _init(); _set( d ); }

  //! @brief Constructor of sparse polynomial expression as variable
  SPolyExpr
    ( FFVar const& x )
    { _init(); _set( x ); }

  //! @brief Constructor of sparse polynomial expression as monomial
  SPolyExpr
    ( std::pair< FFMon, double > const& mon )
    { _init(); _set( mon ); }

  //! @brief Copy constructor of sparse polynomial expression
  SPolyExpr
    ( SPolyExpr const& spoly )
    { _init(); _set( spoly ); }

  //! @brief Destructor of sparse polynomial expression
  virtual ~SPolyExpr()
    {}

  //! @brief Insert sparse polynomial into DAG
  FFVar insert
    ( FFGraph*dag )
    const;

  //! @brief Extract univariate polynpomials from a multivariaate polynomial
  std::map< FFVar const*, SPolyExpr, lt_FFVar > extract_univariate
    ( unsigned const ordmin=3 )
    const;

  //! @brief Total number of monomial terms in polynomial variable
  unsigned nmon
    ()
    const
    { return _mapmon.size(); };

  //! @brief Get const map of monomial coefficients
  const t_ffpoly& mapmon
    ()
    const
    { return _mapmon; }

  //! @brief Get map of monomial coefficients
  t_ffpoly& mapmon
    ()
    { return _mapmon; }

  //! @brief Overloaded operator '=' for sparse polynomial
  SPolyExpr& operator=
    ( SPolyExpr const& spoly )
    { _set( spoly ); return *this; }

  //! @brief Overloaded operator '=' for constant 
  SPolyExpr& operator=
    ( double const d )
    { _set( d ); return *this; }

  //! @brief Overloaded operator '=' for variable
  SPolyExpr& operator=
    ( FFVar const& x )
    { _set( x ); return *this; }

  //! @brief Overloaded operator '=' for monomial
  SPolyExpr& operator=
    ( std::pair< FFMon, double > const& mon )
    { _set( mon ); return *this; }

  //! @brief Overloaded operator '+=' for sparse polynomial
  SPolyExpr& operator+=
    ( const SPolyExpr& spoly );

  //! @brief Overloaded operator '+=' for monomial
  SPolyExpr& operator+=
    ( std::pair< FFMon, double > const& mon );

  //! @brief Overloaded operator '-=' for sparse polynomial
  SPolyExpr& operator-=
    ( const SPolyExpr& spoly );

  //! @brief Overloaded operator '*=' for sparse polynomial
  SPolyExpr& operator*=
    ( const SPolyExpr& spoly );

  //! @brief Overloaded operator '*=' for real scalar
  SPolyExpr& operator*=
    ( const double d );

  //! @brief Overloaded operator '/=' for sparse polynomial
  SPolyExpr& operator/=
    ( const SPolyExpr& spoly );

  //! @brief Overloaded operator '/=' for real denominator
  SPolyExpr& operator/=
    ( const double d );

  //! @brief Exceptions of mc::SPolyExpr
  class Exceptions
  {
   public:
    //! @brief Enumeration type for SPolyExpr exception handling
    enum TYPE{
      DIVZERO = 1,    //!< Scalar division by zero
      DIVPOLY = 2,    //!< Division between two polynomials
      INTERNAL = -33  //!< Internal error
    };
    //! @brief Constructor for error <a>ierr</a>
    Exceptions( TYPE ierr ) : _ierr( ierr ){}
    //! @brief Error flag
    int ierr(){ return _ierr; }
    //! @brief Error description
    std::string what(){
      switch( _ierr ){
      case DIVZERO:
        return "mc::SPolyExpr\t Scalar division by zero";
      case DIVPOLY:
        return "mc::SPolyExpr\t Division between two polynomials";
      case INTERNAL:
      default:
        return "mc::SPolyExpr\t Internal error";
      }
    }
   private:
    TYPE _ierr;
  };

  //! @brief Options of mc::SPolyExpr
  static struct Options
  {
    //! @brief Constructor
    Options():
      BASIS(MONOM), REMZERO(true), DISPLEN(5)
      {}
    //! @brief Assignment of mc::SPolyExpr::Options
    Options& operator=
      ( Options& opt ){
        BASIS   = opt.BASIS;
        REMZERO = opt.REMZERO;
        DISPLEN = opt.DISPLEN;
        return *this;
      }
    //! @brief Available basis representations
    enum BASIS_TYPE{
      MONOM=0,	//!< Monomial basis
      CHEB	    //!< Chebyshev basis
    };
    //! @brief Basis representation of sparse polynomial
    BASIS_TYPE BASIS;
    //! @brief Whether to remove zeros entries from sparse polynomials
    bool REMZERO;
    //! @brief Number of digits in output stream for sparse polynomial coefficients
    unsigned DISPLEN;
  } options;
};

inline SPolyExpr::Options SPolyExpr::options;

//////////////////////////////// SPolyExpr ////////////////////////////////////

inline void
SPolyExpr::_init
()
{
//  _dag = 0;
  return;
}

inline void
SPolyExpr::_cleanup
()
{
  _mapmon.clear();
  _setvar.clear();
}

inline void
SPolyExpr::_reinit
()
{
  _cleanup();
  _init();
}

inline SPolyExpr&
SPolyExpr::_set
( double const d )
{
  _reinit();
  if( !isequal( d, 0. ) )
    _mapmon.insert( std::make_pair( FFMon(), d ) );
  return *this;
}

inline SPolyExpr&
SPolyExpr::_set
( FFVar const& x )
{
  if( x.cst() ) return _set( x.num().val() );
  _reinit();
  _mapmon.insert( std::make_pair( FFMon( x ), 1. ) );
  _setvar.insert( &x );
  return *this;
}

inline SPolyExpr&
SPolyExpr::_set
( std::pair< FFMon, double > const& mon )
{
  _cleanup();
  _mapmon.insert( mon );
  for( auto && var : mon.first.expr )
    _setvar.insert( var.first );
  return *this;
}

inline SPolyExpr&
SPolyExpr::_set
( const SPolyExpr&spoly )
{
  if( this == &spoly ) return *this;
  _mapmon = spoly._mapmon;
  _setvar = spoly._setvar;
  return *this;
}

inline void
SPolyExpr::_clean
( t_ffpoly& spoly )
const
{
  if( !options.REMZERO ) return;
  for( auto itmon = spoly.begin(); itmon!=spoly.end(); )
    if( isequal( itmon->second, 0. ) )
      itmon = spoly.erase( itmon );
    else
      ++itmon;
  return;
}

inline std::map< FFVar const*, SPolyExpr, lt_FFVar >
SPolyExpr::extract_univariate
( unsigned const ordmin )
const
{
  std::map< FFVar const*, SPolyExpr, lt_FFVar > dec;

  // Iterate through list of monomials by decreasing order
  for( auto it = _mapmon.crbegin(); it != _mapmon.crend(); ++it ){
    FFMon const& mon = it->first;
    FFVar const* pvar = 0;

    // Isolate monomials that depend on a single variable
    if( mon.expr.size() == 1 ){
      pvar = mon.expr.cbegin()->first;
      auto jt = dec.find( pvar );
      // Add univariate term if a higher-order univariate term is present already
      if( jt != dec.end() ){
        jt->second += *it;
        continue;
      }
      // Create new univariate expression if order is sufficiently large
      else if( mon.tord >= ordmin ){
        dec.insert( std::make_pair( pvar, *it ) ); // NEED CONSTRUCTOR FOR MONOMIALS!
        continue;
      }
    }

    // Otherwise, retain term in the multivariate polynomial
    auto pt = dec.insert( std::make_pair( pvar, *it ) ); // NEED CONSTRUCTOR FOR MONOMIALS!
    if( !pt.second )
      pt.first->second += *it;
  }
  
  return dec;
}

inline FFVar
SPolyExpr::insert
( FFGraph*dag )
const
{
  FFVar var = 0.;
  for( auto it=_mapmon.begin(); it!=_mapmon.end(); ++it ){
    FFVar prodmon = it->second;
    for( auto ie=it->first.expr.begin(); ie!=it->first.expr.end(); ++ie ){
      auto itVar = dag->Vars().find( const_cast<FFVar*>(ie->first) );
      if( itVar == dag->Vars().end() )
        throw typename SPolyExpr::Exceptions( SPolyExpr::Exceptions::INTERNAL );
      switch( options.BASIS ){
       case Options::MONOM:
        prodmon *= pow( **itVar, (int)ie->second );
        break;
       case Options::CHEB:
        prodmon *= cheb( **itVar, ie->second );
        break;
      }
    }
    if( it->first.expr.empty() ) var = prodmon;
    else                         var+= prodmon;
  }
  // ADD OPTION TO RETURN A PRODMON OR MONOM TERM?
  return var;
}

inline std::ostream&
operator<<
( std::ostream& out, FFMon const& mon )
{
  // Sparse multivariate polynomial
  if( mon.expr.empty() )  out << "1";
  for( auto ie=mon.expr.begin(); ie!=mon.expr.end(); ++ie ){
    if( ie != mon.expr.begin() ) out << "·";
    switch( SPolyExpr::options.BASIS ){
     case SPolyExpr::Options::MONOM:
      out << *ie->first;
      if( ie->second > 1 ) out << "^" << ie->second;
      break;
     case SPolyExpr::Options::CHEB:
      out << "T" << ie->second << "[" << *ie->first << "]";
      break;
    }
  }
  return out;
}

inline std::ostream&
operator<<
( std::ostream& out, SPolyExpr const& spoly )
{
  const unsigned DISPLEN = SPolyExpr::options.DISPLEN;
  out << std::endl << std::scientific << std::setprecision(DISPLEN)
      << std::right;

  // Sparse multivariate polynomial
  for( auto it=spoly._mapmon.begin(); it!=spoly._mapmon.end(); ++it )
    out << "  " << std::right << std::setw(DISPLEN+7) << it->second << "   " << it->first
        << std::endl;

  return out;
}

inline SPolyExpr
operator+
( const SPolyExpr&spoly )
{
  return spoly;
}

inline SPolyExpr&
SPolyExpr::operator+=
( const SPolyExpr&spoly )
{
  for( auto it=spoly._mapmon.begin(); it!=spoly._mapmon.end(); ++it ){
    // No warm-start for insert unfortunately...
    auto pt = _mapmon.insert( *it );
    if( !pt.second ) pt.first->second += it->second;
  }

  _clean( _mapmon );
  _setvar.insert( spoly._setvar.begin(), spoly._setvar.end() );

  return *this;
}

inline SPolyExpr&
SPolyExpr::operator+=
( std::pair< FFMon, double > const& mon )
{
  // No warm-start for insert unfortunately...
  auto pt = _mapmon.insert( mon );
  if( !pt.second ) pt.first->second += mon.second;
  
  _clean( _mapmon );
  for( auto && var : mon.first.expr )
    _setvar.insert( var.first );

  return *this;
}

inline SPolyExpr
operator+
( const SPolyExpr&spoly1, const SPolyExpr&spoly2 )
{
  if( spoly1.nmon() >= spoly2.nmon() ){
    SPolyExpr spoly3( spoly1 );
    spoly3 += spoly2;
    return spoly3;
  }
  SPolyExpr spoly3( spoly2 );
  spoly3 += spoly1;
  return spoly3;
}

inline SPolyExpr
operator-
( const SPolyExpr&spoly )
{
  SPolyExpr spoly2( spoly );
  for( auto it=spoly2._mapmon.begin(); it!=spoly2._mapmon.end(); ++it )
    it->second = -it->second;
  spoly2._setvar.insert( spoly._setvar.begin(), spoly._setvar.end() );
  return spoly2;
}

inline SPolyExpr&
SPolyExpr::operator-=
( const SPolyExpr&spoly )
{
  for( auto it=spoly._mapmon.begin(); it!=spoly._mapmon.end(); ++it ){
    // No warm-start for insert unfortunately...
    auto pt = _mapmon.insert( *it );
    if( pt.second ) pt.first->second = -it->second;
    else            pt.first->second -= it->second;
  }
  _clean( _mapmon );
  _setvar.insert( spoly._setvar.begin(), spoly._setvar.end() );
  return *this;
}

inline SPolyExpr
operator-
( const SPolyExpr&spoly1, const SPolyExpr&spoly2 )
{
  if( spoly1.nmon() >= spoly2.nmon() ){
    SPolyExpr spoly3( spoly1 );
    spoly3 -= spoly2;
    return spoly3;
  }
  SPolyExpr spoly3( -spoly2 );
  spoly3 += spoly1;
  return spoly3;
}

inline SPolyExpr&
SPolyExpr::operator/=
( const double d )
{
  if( d == 0. ) throw typename SPolyExpr::Exceptions( SPolyExpr::Exceptions::DIVZERO );
  if( d == 1. ) return *this;
  for( auto&& mon : _mapmon )
    mon.second /= d;
  return *this;
}

inline SPolyExpr&
SPolyExpr::operator/=
( const SPolyExpr&spoly )
{
  throw typename SPolyExpr::Exceptions( SPolyExpr::Exceptions::DIVPOLY );
}

inline SPolyExpr&
SPolyExpr::operator*=
( const double d )
{
  if( d == 0. ){ *this = 0; return *this; };
  if( d == 1. ) return *this;
  for( auto&& mon : _mapmon )
    mon.second /= d;
  return *this;
}

inline SPolyExpr&
SPolyExpr::operator*=
( const SPolyExpr&spoly )
{
  // Consolidate set of participating variables in product term
  _setvar.insert( spoly._setvar.begin(), spoly._setvar.end() );

  // Construct vectors of coefficients for first partipating variable
  auto itvar = _setvar.begin();
  std::map<unsigned,t_ffpoly> sp1map, sp2map;
  for( auto itmon=_mapmon.begin(); itmon!=_mapmon.end(); ++itmon )
    _svec1D( itvar, *itmon, sp1map );
#ifdef MC__SPOLYEXPR_DEBUG_SPROD
  _sdisp1D( sp1map, itvar, "Poly #1: " );
#endif
  for( auto itmon=spoly._mapmon.begin(); itmon!=spoly._mapmon.end(); ++itmon )
    _svec1D( itvar, *itmon, sp2map );
#ifdef MC__SPOLYEXPR_DEBUG_SPROD
  _sdisp1D( sp2map, itvar, "Poly #2: " );
#endif

  // Call recursive product of univariate Chebyshev polynomials
  _sprod1D( sp1map, sp2map, _mapmon, itvar );
  return *this;
}

inline void
SPolyExpr::_sprod1D
( const std::map<unsigned,t_ffpoly>& sp1map,
  const std::map<unsigned,t_ffpoly>& sp2map,
  t_ffpoly& mapmon, typename t_ffvar::const_iterator itvar )
const
{
  // construct product matrix of polynomial coefficients
  auto itvarnext = itvar; if( !_setvar.empty() ) ++itvarnext;
  std::map<std::pair<unsigned,unsigned>,t_ffpoly> sp12map;
  auto it1map=sp1map.begin();
  for( unsigned imon1=0; it1map!=sp1map.end(); ++it1map, imon1++ ){
    // empty monomial in sp1
    if( it1map->second.empty() ) continue; 
    // constant monomial in sp1
    if( it1map->second.size() == 1 && !it1map->second.begin()->first.tord ){
      auto it2map=sp2map.begin();
      for( unsigned imon2=0;  it2map!=sp2map.end(); ++it2map, imon2++ ){
        auto pmon = sp12map.insert( std::make_pair( std::make_pair(it1map->first,it2map->first), t_ffpoly() ) );
        assert( pmon.second ); // map is initially empty
        _sscal1D( it2map->second, it1map->second.begin()->second, pmon.first->second );
      }
      continue;
    }
    // general monomial in sp1
    auto it2map=sp2map.begin();
    for( unsigned imon2=0;  it2map!=sp2map.end(); ++it2map, imon2++ ){
      // empty monomial in sp2
      if( it2map->second.empty() ) continue; // no term
      // constant monomial in sp2
      if( it2map->second.size() == 1 && !it2map->second.begin()->first.tord ){
        auto pmon = sp12map.insert( std::make_pair( std::make_pair(it1map->first,it2map->first), t_ffpoly() ) );
        assert( pmon.second ); // map is initially empty
        _sscal1D( it1map->second, it2map->second.begin()->second, pmon.first->second );
        continue;
      }
#ifdef MC__SPOLYEXPR_DEBUG_SPROD
      std::cout << "Term (" << imon1 << "," << imon2 << "):\n";
#endif
      std::map<unsigned,t_ffpoly> sp11map, sp22map;
      for( auto itmon=it1map->second.begin(); itmon!=it1map->second.end(); ++itmon )
        _svec1D( itvarnext, *itmon, sp11map );
#ifdef MC__SPOLYEXPR_DEBUG_SPROD
      _sdisp1D( sp11map, itvarnext, "Poly #1: " );
#endif
      for( auto itmon=it2map->second.begin(); itmon!=it2map->second.end(); ++itmon )
        _svec1D( itvarnext, *itmon, sp22map );
#ifdef MC__SPOLYEXPR_DEBUG_SPROD
      _sdisp1D( sp22map, itvarnext, "Poly #2: " );
#endif
      auto pmon = sp12map.insert( std::make_pair( std::make_pair(it1map->first,it2map->first), t_ffpoly() ) );
      assert( pmon.second ); // map is initially empty
      _sprod1D( sp11map, sp22map, pmon.first->second, itvarnext );
    }
  }

  // construct 1D product result and augment remainder as appropriate
  mapmon.clear();
  switch( options.BASIS ){
   case Options::MONOM: default:
    for( auto it12 = sp12map.begin(); it12!=sp12map.end(); ++it12 ){
      // Product involving two constant terms
      if( !it12->first.first && !it12->first.second )
        _slift1D( it12->second, 1., mapmon );
      // Product between non-constant monomial basis functions
      else
        _slift1D( it12->second, 1., mapmon, itvar, it12->first.first + it12->first.second );
    }
    break;
   case Options::CHEB:
    for( auto it12 = sp12map.begin(); it12!=sp12map.end(); ++it12 ){
      // Product involving two constant terms
      if( !it12->first.first && !it12->first.second )
        _slift1D( it12->second, 1., mapmon );
      // Product involving exactly one constant term
      else if( !it12->first.first || !it12->first.second )
        _slift1D( it12->second, 1., mapmon, itvar, it12->first.first + it12->first.second );
      // Product between non-constant Chebyshev basis functions
      else{
        _slift1D( it12->second, .5, mapmon, itvar, it12->first.first + it12->first.second );
        if( it12->first.first == it12->first.second )
          _slift1D( it12->second, .5, mapmon );
        else
          _slift1D( it12->second, .5, mapmon, itvar, it12->first.first > it12->first.second?
                    it12->first.first - it12->first.second: it12->first.second - it12->first.first );
      }
    }
  }

  // clean sparse polynomial
  _clean( mapmon );
#ifdef MC__SPOLYEXPR_DEBUG_SPROD
  _sdisp1D( mapmon, "Prod: " );
#endif
}

inline void
SPolyExpr::_sscal1D
( const t_ffpoly& spoly, const double& dscal, t_ffpoly& spscal )
const
{
  if( isequal(dscal,0.) ) return;
  spscal = spoly;
  if( isequal(dscal,1.) ) return;
  for( auto itmon=spscal.begin(); itmon!=spscal.end(); ++itmon )
    itmon->second *= dscal;
  return;
}

inline void
SPolyExpr::_slift1D
( const t_ffpoly& spoly, const double& dscal, t_ffpoly& splift )
const
{
  if( isequal(dscal,0.) ) return;
  for( auto itmon=spoly.begin(); itmon!=spoly.end(); ++itmon ){
    auto pmon = splift.insert( *itmon );
    if( pmon.second ){ if( !isequal(dscal,1.) ) pmon.first->second *= dscal; }
    else{ pmon.first->second += isequal(dscal,1.)? itmon->second: itmon->second * dscal; }
  }
}

inline void
SPolyExpr::_slift1D
( const t_ffpoly& spoly, const double& dscal, t_ffpoly& splift,
 typename t_ffvar::const_iterator itvar, const unsigned iord )
const
{
  for( auto itmon=spoly.begin(); itmon!=spoly.end(); ++itmon ){
    FFMon mon = itmon->first;
    mon.tord += iord;
    mon.expr.insert( mon.expr.begin(), std::make_pair( *itvar, iord ) );
    auto pmon = splift.insert( std::make_pair( mon, isequal(dscal,1.)? itmon->second: itmon->second*dscal ) );
    if( !pmon.second ) pmon.first->second += isequal(dscal,1.)? itmon->second: itmon->second*dscal;
  }
}

inline void
SPolyExpr::_svec1D
( typename t_ffvar::const_iterator itvar, const std::pair<FFMon,double>&mon,
  std::map<unsigned,t_ffpoly>&mapspoly )
const
{
  auto ie = mon.first.expr.begin();
  // no dependence on variable #itvar 
  if( !mon.first.tord || (ie->first != *itvar && *ie->first != **itvar) )
    mapspoly[ 0 ].insert( mon );
  // dependence on variable #itvar of order ie
  else
    //mapspoly[ ie->second ].insert( std::make_pair( std::make_pair( mon.first.tord-ie->second,
    //  std::map<const FFVar*,unsigned>( ++mon.first.expr.begin(),mon.first.expr.end() ) ),
    //  mon.second ) );
    mapspoly[ ie->second ].insert( std::make_pair( FFMon( mon.first.tord - ie->second,
      std::map<const FFVar*,unsigned>( ++mon.first.expr.begin(),mon.first.expr.end() ) ),
      mon.second ) );
}

inline void
SPolyExpr::_sdisp1D
( const t_ffpoly&spoly, const std::string&name, std::ostream&os )
const
{
  os << name;
  for( auto itmon=spoly.begin(); itmon!=spoly.end(); ++itmon ){
    if( itmon != spoly.begin() ) os << " + ";
    os << itmon->second;
    for( auto ie=itmon->first.expr.begin(); ie!=itmon->first.expr.end(); ++ie )
      os << "·T" << ie->second << "[" << *ie->first << "]";
  }
  os << std::endl;
}

inline void
SPolyExpr::_sdisp1D
( const std::map<unsigned,t_ffpoly>&mapspoly, typename t_ffvar::const_iterator itvar, 
  const std::string&name, std::ostream&os )
const
{
  os << name;
  for( auto itpoly=mapspoly.begin(); itpoly!=mapspoly.end(); ++itpoly ){
    if( itpoly != mapspoly.begin() ) os << " + ";
    os << "T" << itpoly->first << "[" << **itvar << "] · { ";
    for( auto itmon=itpoly->second.begin(); itmon!=itpoly->second.end(); ++itmon ){
      if( itmon != itpoly->second.begin() ) os << " + ";
      os << itmon->second;
      for( auto ie=itmon->first.expr.begin(); ie!=itmon->first.expr.end(); ++ie )
        os << "·T" << ie->second << "[" << *ie->first << "]";
    }
    os << " }";
  }
  os << std::endl;
}

inline SPolyExpr
operator*
( const SPolyExpr&spoly1, const SPolyExpr&spoly2 )
{
  if( spoly1.nmon() >= spoly2.nmon() ){
    SPolyExpr spoly3( spoly1 );
    spoly3 *= spoly2;
    return spoly3;
  }
  SPolyExpr spoly3( spoly2 );
  spoly3 *= spoly1;
  return spoly3;
}

inline SPolyExpr
sqr
( const SPolyExpr&spoly )
{
  // Construct vectors of coefficients for first partipating variable
  auto itvar = spoly._setvar.begin();
  std::map<unsigned,SPolyExpr::t_ffpoly> spmap;
  for( auto itmon=spoly._mapmon.begin(); itmon!=spoly._mapmon.end(); ++itmon )
    spoly._svec1D( itvar, *itmon, spmap );
#ifdef MC__SPOLYEXPR_DEBUG_SPROD
  spoly._sdisp1D( spmap, itvar, "Poly #1: " );
#endif

  // Call recursive product of univariate Chebyshev polynomials
  SPolyExpr spsqr;
  spoly._sprod1D( spmap, spmap, spsqr._mapmon, itvar );
  spsqr._setvar = spoly._setvar;
  return spsqr;
}

inline SPolyExpr
pow
( const SPolyExpr&spoly, const unsigned n )
{
  switch( n ){
   case 0:  return 1.;
   case 1:  return spoly;
   case 2:  return sqr( spoly );
   default: return n%2 ? sqr( pow( spoly, n/2 ) ) * spoly : sqr( pow( spoly, n/2 ) );
  }
}

inline SPolyExpr
cheb
( const SPolyExpr&spoly, const unsigned n )
{
  switch( n ){
   case 0:  return 1.;
   case 1:  return spoly;
   case 2:  return 2 * sqr( spoly ) - 1;
   default: return 2 * spoly * cheb( spoly, n-1 ) - cheb( spoly, n-2 );
  }
}

inline SPolyExpr
prod
( const unsigned int npoly, const SPolyExpr*ppoly )
{
  switch( npoly ){
   case 0:  return 1.;
   case 1:  return ppoly[0];
   default: return ppoly[0] * prod( npoly-1, ppoly+1 );
  }
}

inline SPolyExpr
monom
( const unsigned int npoly, const SPolyExpr*ppoly, const unsigned*k, const bool chebbasis=false )
{
  switch( npoly ){
   case 0:  return 1.;
   case 1:  return chebbasis? cheb( ppoly[0], k[0] ): pow( ppoly[0], k[0] );
   default: return ( chebbasis? cheb( ppoly[0], k[0] ): pow( ppoly[0], k[0] ) ) * monom( npoly-1, ppoly+1, k+1 );
  }
}

inline SPolyExpr
operator/
( const SPolyExpr&spoly1, const SPolyExpr&spoly2 )
{
  SPolyExpr spoly3( spoly1 );
  return spoly3 /= spoly2;
}

} // namespace mc

//#include "mcop.hpp"

namespace mc
{

//! @brief Specialization of the structure mc::Op to allow usage of the type mc::SPolyExpr for DAG evaluation or as a template parameter in other MC++ classes
template <> struct Op<mc::SPolyExpr>
{
  typedef mc::SPolyExpr T;
  static T point( const double c ) { return T(c); }
  static T zeroone() { throw std::runtime_error("operation not permitted"); }
  static void I(T& x, const T&y) { x = y; }
  static double l(const T& x) { throw std::runtime_error("operation not permitted"); }
  static double u(const T& x) { throw std::runtime_error("operation not permitted"); }
  static double abs (const T& x) { throw std::runtime_error("operation not permitted");  }
  static double mid (const T& x) { throw std::runtime_error("operation not permitted");  }
  static double diam(const T& x) { throw std::runtime_error("operation not permitted"); }
  static T inv (const T& x) { throw std::runtime_error("operation not permitted");  }
  static T sqr (const T& x) { return mc::sqr(x);  }
  static T sqrt(const T& x) { throw std::runtime_error("operation not permitted"); }
  static T exp (const T& x) { throw std::runtime_error("operation not permitted");  }
  static T log (const T& x) { throw std::runtime_error("operation not permitted");  }
  static T xlog(const T& x) { throw std::runtime_error("operation not permitted"); }
  static T lmtd(const T& x, const T& y) { throw std::runtime_error("operation not permitted"); }
  static T rlmtd(const T& x, const T& y) { throw std::runtime_error("operation not permitted"); }
  static T fabs(const T& x) { throw std::runtime_error("operation not permitted"); }
  static T sin (const T& x) { throw std::runtime_error("operation not permitted");  }
  static T cos (const T& x) { throw std::runtime_error("operation not permitted");  }
  static T tan (const T& x) { throw std::runtime_error("operation not permitted");  }
  static T asin(const T& x) { throw std::runtime_error("operation not permitted"); }
  static T acos(const T& x) { throw std::runtime_error("operation not permitted"); }
  static T atan(const T& x) { throw std::runtime_error("operation not permitted"); }
  static T sinh(const T& x) { throw std::runtime_error("operation not permitted"); }
  static T cosh(const T& x) { throw std::runtime_error("operation not permitted"); }
  static T tanh(const T& x) { throw std::runtime_error("operation not permitted"); }
  static T erf (const T& x) { throw std::runtime_error("operation not permitted");  }
  static T erfc(const T& x) { throw std::runtime_error("operation not permitted"); }
  static T fstep(const T& x) { throw std::runtime_error("operation not permitted"); }
  static T bstep(const T& x) { throw std::runtime_error("operation not permitted"); }
  static T min (const T& x, const T& y) { throw std::runtime_error("operation not permitted");  }
  static T max (const T& x, const T& y) { throw std::runtime_error("operation not permitted");  }
  static T arh (const T& x, const double k) { throw std::runtime_error("operation not permitted"); }
  template <typename EXP> static T pow(const T& x, const EXP& y) { return mc::pow(x,y); }
  static T cheb (const T& x, const unsigned n) { return mc::cheb(x,n); }
  static T prod (const unsigned int n, const T* x) { return mc::prod(n,x); }
  static T monom(const unsigned int n, const T* x, const unsigned* k, const bool cheb=false) { return mc::monom(n,x,k,cheb); }
  static T hull(const T& x, const T& y) { throw std::runtime_error("operation not permitted"); }
  static bool inter(T& xIy, const T& x, const T& y) { throw std::runtime_error("operation not permitted"); }
  static bool eq(const T& x, const T& y) { throw std::runtime_error("operation not permitted"); }
  static bool ne(const T& x, const T& y) { throw std::runtime_error("operation not permitted"); }
  static bool lt(const T& x, const T& y) { throw std::runtime_error("operation not permitted");  }
  static bool le(const T& x, const T& y) { throw std::runtime_error("operation not permitted"); }
  static bool gt(const T& x, const T& y) { throw std::runtime_error("operation not permitted");  }
  static bool ge(const T& x, const T& y) { throw std::runtime_error("operation not permitted"); }
};

} // namespace mc

#endif

