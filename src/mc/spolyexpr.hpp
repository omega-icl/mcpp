// Copyright (C) 2017 Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the Eclipse Public License.

#ifndef MC__SPOLYEXPR_H
#define MC__SPOLYEXPR_H

#include <set>
#include <map>

#include "ffunc.hpp"

#undef  MC__SPOLYEXPR_DEBUG_SPROD

namespace mc
{
//! @brief C++ structure for ordering of monomial in graded lexicographic order (grlex)
struct lt_grlex
{
  typedef std::pair< unsigned, std::map< const FFVar*, unsigned > > t_ffmon;

  bool operator()
    ( const t_ffmon&Mon1, const t_ffmon&Mon2 ) const
    {
      // Order monomials based on their total order first
      if( Mon1.first < Mon2.first ) return true;
      if( Mon1.first > Mon2.first ) return false;
      // Account for the case of an empty list
      if( Mon1.second.empty() ) return true;
      if( Mon2.second.empty() ) return false;
      // Order in graded lexicographic order next
      for( auto it1=Mon1.second.begin(), it2=Mon2.second.begin(); it1!=Mon1.second.end(); ++it1, ++it2 ){
        if( lt_FFVar()( it1->first, it2->first ) ) return true;
        if( lt_FFVar()( it2->first, it1->first ) ) return false;
        if( it1->second > it2->second ) return true;
        if( it1->second < it2->second ) return false;
      }
      return false;
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

  // Monomial representation: <total order, <DAG variable ptr, order>>
  typedef std::pair< unsigned, std::map< const FFVar*, unsigned > > t_ffmon;
  typedef std::map< t_ffmon, double > t_ffpoly;
  typedef std::set< const FFVar* > t_ffvar;

  // Friends for arithmetic operations
  friend SPolyExpr sqr ( const SPolyExpr& );
  friend SPolyExpr operator- ( const SPolyExpr& );
  friend std::ostream& operator<< ( std::ostream&, const SPolyExpr& );

protected:
//  //! @brief Pointer to underlying factorable function DAG - _dag := NULL for variable identifier NOREF
//  mutable FFGraph *_dag;

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
    ( const SPolyExpr& spoly );

  //! @brief Set polynomial expression equal to constant <a>d</a>
  SPolyExpr& _set
    ( const double d );

  //! @brief Set polynomial expression equal to variable <a>x</a>
  SPolyExpr& _set
    ( const FFVar& x );

  //! @brief Clean sparse polynomial by removing zero entries
  void _clean
    ( t_ffpoly& spoly )
    const;

  //! @brief Build univariate sparse polynomial (with sparse polynomial coefficients)
  void _svec1D
    ( typename t_ffvar::const_iterator itvar, const std::pair<t_ffmon,double>&mon,
      std::map<unsigned,t_ffpoly>&mapspoly )
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
  /** @ingroup FFunc
   *  @{
   */
  //! @brief Default constructor of sparse polynomial expression
  SPolyExpr
    ()
    { _init(); }

  //! @brief Constructor of sparse polynomial expression as constant
  SPolyExpr
    ( const double d )
    { _init(); _set( d ); }

  //! @brief Constructor of sparse polynomial expression as variable
  SPolyExpr
    ( const FFVar& x )
    { _init(); _set( x ); }

  //! @brief Copy constructor of sparse polynomial expression
  SPolyExpr
    ( const SPolyExpr& spoly )
    { _init(); _set( spoly ); }

  //! @brief Destructor of sparse polynomial expression
  virtual ~SPolyExpr()
    {}

//  //! @brief Set multivariate polynomial coefficients in variable as <tt>coefmon</tt>
//  SPolyExpr& set
//    ( t_ffpoly& mapmon )
//    { _mapmon = mapmon;
//      return *this; } // this is assuming the same order and number of variables

  //! @brief Insert sparse polynomial into DAG
  FFVar insert
    ( FFGraph*dag );

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
    ( const SPolyExpr& spoly )
    { _set( spoly ); return *this; }

  //! @brief Overloaded operator '=' for constant 
  SPolyExpr& operator=
    ( const double d )
    { _set( d ); return *this; }

  //! @brief Overloaded operator '=' for variable
  SPolyExpr& operator=
    ( const FFVar& x )
    { _set( x ); return *this; }

  //! @brief Overloaded operator '+=' for sparse polynomial
  SPolyExpr& operator+=
    ( const SPolyExpr& spoly );

  //! @brief Overloaded operator '-=' for sparse polynomial
  SPolyExpr& operator-=
    ( const SPolyExpr& spoly );

  //! @brief Overloaded operator '*=' for sparse polynomial
  SPolyExpr& operator*=
    ( const SPolyExpr& spoly );

  //! @brief Overloaded operator '*=' for real scalar
  SPolyExpr& operator*=
    ( const double d );

  //! @brief Overloaded operator '/=' for real denominator
  SPolyExpr& operator/=
    ( const double d );

  //! @brief Exceptions of mc::SPolyExpr
  class Exceptions
  {
   public:
    //! @brief Enumeration type for SPolyExpr exception handling
    enum TYPE{
//      DAG=0,          //!< Operation between sparse polynomials linked to different DAGs
      DIVZERO = 1,    //!< Scalar division by zero
      INTERNAL = -33  //!< Internal error
    };
    //! @brief Constructor for error <a>ierr</a>
    Exceptions( TYPE ierr ) : _ierr( ierr ){}
    //! @brief Error flag
    int ierr(){ return _ierr; }
    //! @brief Error description
    std::string what(){
      switch( _ierr ){
//      case DAG:
//        return "mc::SPolyExpr\t Operation between sparse polynomials linked to different DAGs is not allowed";
      case DIVZERO:
        return "mc::SPolyExpr\t Scalar division by zero";
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
  /** @} */
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
( const double d )
{
  _reinit();
  if( !isequal( d, 0. ) )
    _mapmon.insert( std::make_pair( std::make_pair(0,std::map<const FFVar*,unsigned>()), d ) );
  return *this;
}

inline SPolyExpr&
SPolyExpr::_set
( const FFVar&x )
{
  if( x.cst() ) return _set( x.num().val() );
  _reinit();
//  _dag = x.dag();
  std::map<const FFVar*,unsigned> mon; mon.insert( std::make_pair(&x,1) );
  _mapmon.insert( std::make_pair( std::make_pair(1,mon), 1. ) );
  _setvar.insert( &x );
  return *this;
}

inline SPolyExpr&
SPolyExpr::_set
( const SPolyExpr&spoly )
{
  if( this == &spoly ) return *this;
//  _dag = spoly._dag;
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

inline FFVar
SPolyExpr::insert
( FFGraph*dag )
{
  FFVar var = 0.;
  for( auto it=_mapmon.begin(); it!=_mapmon.end(); ++it ){
    FFVar prodmon = it->second;
    for( auto ie=it->first.second.begin(); ie!=it->first.second.end(); ++ie ){
      auto itVar = dag->Vars().find( const_cast<FFVar*>(ie->first) );
      if( itVar == dag->Vars().end() )
        throw typename SPolyExpr::Exceptions( SPolyExpr::Exceptions::INTERNAL );
      switch( SPolyExpr::options.BASIS ){
       case SPolyExpr::Options::MONOM:
        prodmon *= pow( **itVar, (int)ie->second );
        break;
       case SPolyExpr::Options::CHEB:
        prodmon *= cheb( **itVar, ie->second );
        break;
      }
    }
    if( it->first.second.empty() ) var = prodmon;
    else                           var+= prodmon;
  }
  // ADD OPTION TO RETURN A PRODMON OR MONOM TERM?
  return var;
}

inline std::ostream&
operator<<
( std::ostream&out, const SPolyExpr&spoly )
{
  const unsigned DISPLEN = SPolyExpr::options.DISPLEN;
  out << std::endl << std::scientific << std::setprecision(DISPLEN)
      << std::right;

  // Sparse multivariate polynomial
  for( auto it=spoly._mapmon.begin(); it!=spoly._mapmon.end(); ++it ){
    out << std::right << std::setw(DISPLEN+7) << it->second << "   ";
    for( auto ie=it->first.second.begin(); ie!=it->first.second.end(); ++ie ){
      if( ie!=it->first.second.begin() ) out << "·";
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
    out << std::endl;
  }

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
//  if( _dag && spoly._dag && _dag != spoly._dag )
//    throw typename SPolyExpr::Exceptions( SPolyExpr::Exceptions::DAG );

  for( auto it=spoly._mapmon.begin(); it!=spoly._mapmon.end(); ++it ){
    // No warm-start for insert unfortunately...
    auto pt = _mapmon.insert( *it );
    if( !pt.second ) pt.first->second += it->second;
  }
  _clean( _mapmon );
  _setvar.insert( spoly._setvar.begin(), spoly._setvar.end() );
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
//  if( _dag && spoly._dag && _dag != spoly._dag )
//    throw typename SPolyExpr::Exceptions( SPolyExpr::Exceptions::DAG );

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
//  if( _dag && spoly._dag && _dag != spoly._dag )
//    throw typename SPolyExpr::Exceptions( SPolyExpr::Exceptions::DAG );

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
    if( it1map->second.size() == 1 && !it1map->second.begin()->first.first ){
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
      if( it2map->second.size() == 1 && !it2map->second.begin()->first.first ){
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
    t_ffmon ffmon = itmon->first;
    ffmon.first += iord;
    ffmon.second.insert( ffmon.second.begin(), std::make_pair( *itvar, iord ) );
    auto pmon = splift.insert( std::make_pair( ffmon, isequal(dscal,1.)? itmon->second: itmon->second*dscal ) );
    if( !pmon.second ) pmon.first->second += isequal(dscal,1.)? itmon->second: itmon->second*dscal;
  }
}

inline void
SPolyExpr::_svec1D
( typename t_ffvar::const_iterator itvar, const std::pair<t_ffmon,double>&mon,
  std::map<unsigned,t_ffpoly>&mapspoly )
const
{
  auto ie = mon.first.second.begin();
  // no dependence on variable #itvar 
  if( !mon.first.first || *ie->first != **itvar )
    mapspoly[ 0 ].insert( mon );
  // dependence on variable #itvar of order ie
  else
    mapspoly[ ie->second ].insert( std::make_pair( std::make_pair( mon.first.first-ie->second,
      std::map<const FFVar*,unsigned>( ++mon.first.second.begin(),mon.first.second.end() ) ),
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
    for( auto ie=itmon->first.second.begin(); ie!=itmon->first.second.end(); ++ie )
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
      for( auto ie=itmon->first.second.begin(); ie!=itmon->first.second.end(); ++ie )
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

} // namespace mc

#endif
