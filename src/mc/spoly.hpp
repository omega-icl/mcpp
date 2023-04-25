// Copyright (C) 2022 Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the Eclipse Public License.

#ifndef MC__SPOLY_H
#define MC__SPOLY_H

#include <vector>
#include <set>
#include "smon.hpp"
#include "mcfunc.hpp"

#undef  MC__SPOLY_DEBUG_SPROD

namespace mc
{

//! @brief C++ class for sparse polynomial representation and arithmetic
////////////////////////////////////////////////////////////////////////
//! mc::SPoly is a C++ class for sparse polynomial representation
//! and arithmetic
////////////////////////////////////////////////////////////////////////
template <typename KEY=unsigned, typename COMP=std::less<unsigned>>
class SPoly
////////////////////////////////////////////////////////////////////////
{
public:

  // Monomial representation: SMon := <total order, <variable, order>>
  typedef SMon< KEY, COMP > t_mon;
  typedef lt_SMon< COMP > lt_mon;
  typedef std::map< t_mon, double, lt_mon > t_poly;
  typedef std::set< KEY, COMP > t_var;

  // Friends for arithmetic operations
  template <typename K, typename C> friend SPoly<K,C> sqr ( SPoly<K,C> const& );
  template <typename K, typename C> friend SPoly<K,C> operator- ( SPoly<K,C> const& );
  template <typename K, typename C> friend std::ostream& operator<< ( std::ostream&, SPoly<K,C> const& );

protected:

  //! @brief Map of monomial terms in polynomial expression
  t_poly _mapmon;

  //! @brief Set of participating variables in polynomial expression
  t_var _setvar;

  //! @brief Initialize private/protected members of model variable
  void _init
    ();

  //! @brief Reinitialize private/protected members of model variable
  void _reinit
    ();

  //! @brief Clean up private/protected members of variable
  void _cleanup
    ();

  //! @brief Set polynomial expression equal to constant <a>d</a>
  SPoly<KEY,COMP>& _set
    ( double const d );

  //! @brief Set polynomial expression equal to variable <a>x</a>
  SPoly<KEY,COMP>& _set
    ( KEY const& x );

  //! @brief Set polynomial expression equal to monomial <a>mon</a>
  SPoly<KEY,COMP>& _set
    ( std::pair<t_mon,double> const& mon );

  //! @brief Set polynomial expression equal to monomial map <a>mapmon</a>
  SPoly<KEY,COMP>& _set
    ( t_poly const& mapmon );

  //! @brief Set polynomial expression equal to <a>spoly</a>
  SPoly<KEY,COMP>& _set
    ( SPoly<KEY,COMP> const& spoly );
    
  //! @brief Clean sparse polynomial by removing zero entries
  void _clean
    ( t_poly & spoly )
    const;

  //! @brief Recursive product of univariate sparse polynomials
  void _sprod1D
    ( std::map<unsigned,t_poly> const& sp1map,
      std::map<unsigned,t_poly> const& sp2map,
      t_poly& mapmon, t_var const& setvar,
      typename t_var::const_iterator itvar )
    const;

  //! @brief Display of recursive univariate sparse polynomial
  void _sdisp1D
    ( t_poly const& spoly, std::string const& name, std::ostream& os=std::cout )
    const;

  //! @brief Display of recursive univariate sparse polynomial
  void _sdisp1D
    ( std::map<unsigned,t_poly> const& mapspoly,
      typename t_var::const_iterator itvar, 
      std::string const& name, std::ostream& os=std::cout )
    const;

  //! @brief Build univariate sparse polynomial (with sparse polynomial coefficients)
  static void _svec1D
    ( typename t_var::const_iterator itvar,
      std::pair<t_mon,double> const& mon,
      std::map<unsigned,t_poly>& mapspoly );

  //! @brief Scaling of univariate sparse polynomial
  static void _sscal1D
    ( t_poly const& spoly, double const& dscal, t_poly& spscal );

  //! @brief Lifting of sparse polynomial coefficient maps
  static void _slift1D
    ( t_poly const& spoly, double const& dscal, t_poly& splift );

  //! @brief Lifting of sparse polynomial coefficient maps
  static void _slift1D
    ( t_poly const& spoly, double const& dscal, t_poly& splift,
      typename t_var::const_iterator itvar, unsigned const iord );
      
  //! @brief Convert univariate polynomial from Chebyshev to monomial basis
  void _convert
    ( t_poly& mapmon, t_var const& setvar, int const BASIS )
    const;

  //! @brief Build univariate sparse polynomial (with sparse polynomial coefficients)
  static void _svec1Dfull
    ( typename t_var::const_iterator itvar,
      std::pair<t_mon,double> const& mon,
      std::map<unsigned,t_poly>& mapspoly );

  //! @brief Convert univariate polynomial from monomial to Chebyshev basis
  static void _pow2cheb
    ( std::map<unsigned,t_poly>& mapspoly );
    
  //! @brief Convert univariate polynomial from Chebyshev to monomial basis
  static void _cheb2pow
    ( std::map<unsigned,t_poly>& mapspoly );

public:

  //! @brief Default constructor of sparse polynomial expression
  SPoly
    ()
    { _init(); }

  //! @brief Constructor of sparse polynomial expression as constant
  SPoly
    ( double const d )
    { _init(); _set( d ); }

  //! @brief Constructor of sparse polynomial expression as variable
  //SPoly
  //  ( KEY const& x )
  //  { _init(); _set( x ); }

  //! @brief Constructor of sparse polynomial expression as monomial
  SPoly
    ( std::pair< t_mon, double > const& mon )
    { _init(); _set( mon ); }

  //! @brief Constructor of sparse polynomial expression as monomial map
  SPoly
    ( t_poly const& mapmon )
    { _init(); _set( mapmon ); }

  //! @brief Copy constructor of sparse polynomial expression
  SPoly
    ( SPoly<KEY,COMP> const& spoly )
    { _init(); _set( spoly ); }

  //! @brief Destructor of sparse polynomial expression
  virtual ~SPoly()
    {}

  //! @brief Extract univariate polynomials from a multivariate polynomial
  //std::map< KEY, SPoly<KEY,COMP>, COMP > extract_univariate
  //  ( unsigned const ordmin=3 )
  //  const;

  //! @brief Display sparse coefficient map
  std::string display
    ( t_poly spoly, int const& BASIS=options.BASIS, int const& DISPLEN=options.DISPLEN )
    const;

  //! @brief Convert coefficient map to desired BASIS
  SPoly<KEY,COMP>& convert
    ( int const BASIS );

  //! @brief Maximal degree of monomial terms in polynomial variable
  unsigned maxord
    ()
    const
    { return _mapmon.empty()? 0: _mapmon.crbegin()->first.tord; }

  //! @brief Minimal degree of monomial terms in polynomial variable
  unsigned minord
    ()
    const
    { return _mapmon.empty()? 0: _mapmon.cbegin()->first.tord; }

  //! @brief Minimal degree of variable monomial terms in polynomial variable
  unsigned minord
    ( KEY const& x )
    const
    { if( _mapmon.empty() || !minord() ) return 0;
      unsigned ord = maxord();
      for( auto const& [mon,coef] : mapmon ){
        unsigned const exp = mon.exp(x);
        if( exp < ord ) ord = exp;
      } }

  //! @brief Total number of monomial terms in polynomial variable
  unsigned nmon
    ()
    const
    { return _mapmon.size(); };

  //! @brief Get coefficient of monomial mon
  double coefmon
    ( t_mon const& mon )
    const
    { auto it = _mapmon.find( mon ); 
      return( it != _mapmon.end()? it->second: 0e0 ); }

  //! @brief Get const map of monomial coefficients
  t_poly const& mapmon
    ()
    const
    { return _mapmon; }

  //! @brief Get map of monomial coefficients
  t_poly& mapmon
    ()
    { return _mapmon; }

  //! @brief Get const map of monomial coefficients
  t_var const& setvar
    ()
    const
    { return _setvar; }

  //! @brief Get map of monomial coefficients
  t_var& setvar
    ()
    { return _setvar; }

  //! @brief Set variable key
  SPoly<KEY,COMP>& var
    ( KEY const& x )
    { _set( x ); return *this; }

  //! @brief Overloaded operator '=' for sparse polynomial
  SPoly<KEY,COMP>& operator=
    ( SPoly<KEY,COMP> const& spoly )
    { _set( spoly ); return *this; }

  //! @brief Overloaded operator '=' for constant 
  SPoly<KEY,COMP>& operator=
    ( double const& d )
    { _set( d ); return *this; }

  //! @brief Overloaded operator '=' for monomial
  SPoly<KEY,COMP>& operator=
    ( std::pair< t_mon, double > const& mon )
    { _set( mon ); return *this; }

  //! @brief Overloaded operator '=' for monomial map
  SPoly<KEY,COMP>& operator=
    ( t_poly const& mapmon )
    { _set( mapmon ); return *this; }

  //! @brief Overloaded operator '+=' for constant
  SPoly<KEY,COMP>& operator+=
    ( double const& d );

  //! @brief Overloaded operator '+=' for sparse polynomial
  SPoly<KEY,COMP>& operator+=
    ( SPoly<KEY,COMP> const& spoly );

  //! @brief Overloaded operator '+=' for monomial
  SPoly<KEY,COMP>& operator+=
    ( std::pair< t_mon, double > const& mon );

  //! @brief Overloaded operator '-=' for constant
  SPoly<KEY,COMP>& operator-=
    ( double const& d );

  //! @brief Overloaded operator '-=' for sparse polynomial
  SPoly<KEY,COMP>& operator-=
    ( SPoly<KEY,COMP> const& spoly );

  //! @brief Overloaded operator '*=' for sparse polynomial
  SPoly<KEY,COMP>& operator*=
    ( SPoly<KEY,COMP> const& spoly );

  //! @brief Overloaded operator '*=' for constant
  SPoly<KEY,COMP>& operator*=
    ( double const& d );

  //! @brief Overloaded operator '/=' for sparse polynomial
  SPoly<KEY,COMP>& operator/=
    ( SPoly<KEY,COMP> const& spoly );

  //! @brief Overloaded operator '/=' for constant
  SPoly<KEY,COMP>& operator/=
    ( double const& d );

  //! @brief Exceptions of mc::SPoly
  class Exceptions
  {
   public:
    //! @brief Enumeration type for SPoly exception handling
    enum TYPE{
      DIVZERO = 1,    //!< Scalar division by zero
      DIVPOLY = 2,    //!< Division between two polynomials
      //WRONGBASIS = 3, //!< Operation between two polynomials in different bases
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
        return "mc::SPoly\t Scalar division by zero";
      case DIVPOLY:
        return "mc::SPoly\t Division between two polynomials";
      //case WRONGBASIS:
      //  return "mc::SPoly\t Operation between two polynomials in different bases";
      case INTERNAL:
      default:
        return "mc::SPoly\t Internal error";
      }
    }
   private:
    TYPE _ierr;
  };

  //! @brief Options of mc::SPoly
  static struct Options
  {
    //! @brief Constructor
    Options():
      BASIS(MONOM), REMZERO(true), DISPLEN(5)
      {}
    //! @brief Assignment of mc::SPoly::Options
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
    int BASIS;
    //! @brief Whether to remove zeros entries from sparse polynomials
    bool REMZERO;
    //! @brief Number of digits in output stream for sparse polynomial coefficients
    unsigned DISPLEN;
  } options;
};

template <typename KEY, typename COMP>
inline typename SPoly<KEY,COMP>::Options SPoly<KEY,COMP>::options;

//////////////////////////////// SPoly ////////////////////////////////////

template <typename KEY, typename COMP>
inline void
SPoly<KEY,COMP>::_init
()
{
  return;
}

template <typename KEY, typename COMP>
inline void
SPoly<KEY,COMP>::_cleanup
()
{
  _mapmon.clear();
  _setvar.clear();
}

template <typename KEY, typename COMP>
inline void
SPoly<KEY,COMP>::_reinit
()
{
  _cleanup();
  _init();
}

template <typename KEY, typename COMP>
inline SPoly<KEY,COMP>&
SPoly<KEY,COMP>::_set
( double const d )
{
  _reinit();
  if( !isequal( d, 0. ) )
    _mapmon.insert( std::make_pair( t_mon(), d ) );
  return *this;
}

template <typename KEY, typename COMP>
inline SPoly<KEY,COMP>&
SPoly<KEY,COMP>::_set
( KEY const& x )
{
  _reinit();
  _mapmon.insert( std::make_pair( t_mon( x ), 1. ) );
  _setvar.insert( x );
  return *this;
}

template <typename KEY, typename COMP>
inline SPoly<KEY,COMP>&
SPoly<KEY,COMP>::_set
( std::pair<t_mon,double> const& mon )
{
  _cleanup();
  _mapmon.insert( mon );
  for( auto && var : mon.first.expr )
    _setvar.insert( var.first );
  return *this;
}

template <typename KEY, typename COMP>
inline SPoly<KEY,COMP>&
SPoly<KEY,COMP>::_set
( t_poly const& mapmon )
{
  _cleanup();
  _mapmon = mapmon;
  for( auto const& [mon,coef] : mapmon )
    for( auto && var : mon.expr )
      _setvar.insert( var.first );
  return *this;
}

template <typename KEY, typename COMP>
inline SPoly<KEY,COMP>&
SPoly<KEY,COMP>::_set
( SPoly<KEY,COMP> const& spoly )
{
  if( this == &spoly ) return *this;
  _mapmon = spoly._mapmon;
  _setvar = spoly._setvar;
  return *this;
}

template <typename KEY, typename COMP>
inline void
SPoly<KEY,COMP>::_clean
( t_poly& spoly )
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

template <typename KEY, typename COMP>
inline void
SPoly<KEY,COMP>::_pow2cheb
( std::map<unsigned,t_poly>& spmap )
{
  // Converts COEFF in monomial basis to the Chebyshev basis.
  // Based on SCONMC function (http://www.netlib.org/math/MATH77/dconmc.f, http://www.netlib.org/math/docpdf/ch11-03.pdf)
  int const N = spmap.crbegin()->first;
  if( N <= 1 ) return;
  // TP = .5D0**(N-1)
  // COEFF(N) = TP * COEFF(N)
  // COEFF(N-1) = TP * COEFF(N-1)
  double TP = std::pow(0.5,N-1);
  for( auto& [mon,coef] : spmap[N] )   coef *= TP;
  for( auto& [mon,coef] : spmap[N-1] ) coef *= TP;
  //    do 20 J = N-2, 0, -1
  //      TP = 2.D0 * TP
  //      COEFF(J) = TP * COEFF(J)
  //      COEFF(J+1) = 2.D0 * COEFF(J+1)
  for( int J=N-2; J>=0; --J ){
    TP *= 2e0;
    for( auto& [mon,coef] : spmap[J] )   coef *= TP;
    for( auto& [mon,coef] : spmap[J+1] ) coef *= 2e0;
    //    do 10 I = J, N-2
    //      COEFF(I) = COEFF(I) + COEFF(I+2)
    // 10    continue
    for( int I=J; I<=N-2; ++I ){
      for( auto& [mon,coef] : spmap[I+2] ){
        auto [itmon,ins] = spmap[I].insert( std::make_pair( mon, coef ) );
        if( !ins ) itmon->second += coef;
        if( isequal( itmon->second, 0. ) ) spmap[I].erase( itmon );
      }
    }
  // 20 continue
  }
  //    return
  //    end  
}

template <typename KEY, typename COMP>
inline void
SPoly<KEY,COMP>::_cheb2pow
( std::map<unsigned,t_poly>& spmap )
{
  // Converts COEFF in Chebyshev basis to the monomial basis.
  // Based on SCONCM function (http://www.netlib.org/math/MATH77/dconcm.f, http://www.netlib.org/math/docpdf/ch11-03.pdf)
  int const N = spmap.crbegin()->first;
  if( N <= 1 ) return;
  //    TP = 1.D0
  double TP = 1e0;
  //    do 20 J = 0, N-2
  for( int J=0; J<=N-2; ++J ){
  //       do 10 I = N-2, J, -1
  //          COEFF(I) = COEFF(I) - COEFF(I+2)
  // 10    continue
    for( int I=N-2; I>=J; --I ){
      for( auto& [mon,coef] : spmap[I+2] ){
        auto [itmon,ins] = spmap[I].insert( std::make_pair( mon, -coef ) );
        if( !ins ) itmon->second -= coef;
        if( isequal( itmon->second, 0. ) ) spmap[I].erase( itmon );
      }
    }
  //       COEFF(J+1) = .5D0 * COEFF(J+1)
  //       COEFF(J) = TP * COEFF(J)
  //       TP = 2.D0 * TP
  // 20 continue
    for( auto& [mon,coef] : spmap[J+1] ) coef /= 2e0;
    for( auto& [mon,coef] : spmap[J] )   coef *= TP;
    TP *= 2e0;
  }
  //    COEFF(N) = TP * COEFF(N)
  //    COEFF(N-1) = TP * COEFF(N-1)
  for( auto& [mon,coef] : spmap[N] )   coef *= TP;
  for( auto& [mon,coef] : spmap[N-1] ) coef *= TP;
  //    return
  //    end  
}

template <typename KEY, typename COMP>
inline void
SPoly<KEY,COMP>::_svec1Dfull
( typename t_var::const_iterator itvar,
  std::pair<t_mon,double> const& coefmon,
  std::map<unsigned,t_poly>& mapspoly )
{
  auto& [mon,coef] = coefmon;
  auto ie = mon.expr.find( *itvar );
  if( ie == mon.expr.end() ) // no dependence on variable *itvar 
    mapspoly[ 0 ].insert( coefmon );
  else{ // dependence on variable *itvar of order iord
    auto& [ivar,iord] = *ie;
    t_mon monmod( mon.tord - iord, mon.expr );
    monmod.expr.erase( ivar ); // remove *itvar entry
    mapspoly[ iord ].insert( std::make_pair( monmod, coef ) );
  }
}

template <typename KEY, typename COMP>
inline SPoly<KEY,COMP>&
SPoly<KEY,COMP>::convert
( int const BASIS )
{
  if( BASIS != options.BASIS) _convert( _mapmon, _setvar, BASIS );
  return *this;
}

template <typename KEY, typename COMP>
inline void
SPoly<KEY,COMP>::_convert
( t_poly& mapmon, t_var const& setvar, int const BASIS )
const
{
  // Constant or linear polynomial - nothing to do
  if( mapmon.empty() || setvar.empty() || mapmon.crbegin()->first.tord <= 1 ) return;

  // Apply conversion variable by variable
  for( auto itvar=setvar.cbegin(); itvar!=setvar.cend(); ++itvar ){
    // Isolate factors of variable *itvar and convert to desired BASIS
    std::map<unsigned,t_poly> spmap;
    for( auto const& mon : mapmon ) _svec1Dfull( itvar, mon, spmap );
#ifdef MC__SPOLY_DEBUG_CONVERT
    _sdisp1D( spmap, itvar, "before:" );
#endif
    switch( BASIS ){
      case Options::CHEB:  _pow2cheb( spmap ); break;
      case Options::MONOM: _cheb2pow( spmap ); break;
    }
#ifdef MC__SPOLY_DEBUG_CONVERT
    _sdisp1D( spmap, itvar, "after:" );
#endif
    // Merge back into multivariate polynomial
    mapmon.clear();
    for( auto&& [iord,spoly] : spmap ){
      if( !iord ){
        mapmon = spmap[0];
        continue;
      }
      for( auto&& [mon,coef] : spoly ){
        t_mon monmod( mon.tord + iord, mon.expr );
        monmod.expr[ *itvar ] = iord;
        mapmon[ monmod ] = coef;
      }
    }
#ifdef MC__SPOLY_DEBUG_CONVERT
    std::cout << display( mapmon, options.BASIS, options.DISPLEN );
#endif
    //_clean( mapmon );
  }
}

//template <typename KEY, typename COMP>
//inline std::map< KEY, SPoly<KEY,COMP>, COMP >
//SPoly::extract_univariate
//( unsigned const ordmin )
//const
//{
//  std::map< KEY, SPoly<KEY,COMP>, COMP > dec;

//  // Iterate through list of monomials by decreasing order
//  for( auto it = _mapmon.crbegin(); it != _mapmon.crend(); ++it ){
//    auto const& mon = it->first;
//    KEY pvar = 0;

//    // Isolate monomials that depend on a single variable
//    if( mon.expr.size() == 1 ){
//      pvar = mon.expr.cbegin()->first;
//      auto jt = dec.find( pvar );
//      // Add univariate term if a higher-order univariate term is present already
//      if( jt != dec.end() ){
//        jt->second += *it;
//        continue;
//      }
//      // Create new univariate expression if order is sufficiently large
//      else if( mon.tord >= ordmin ){
//        dec.insert( std::make_pair( pvar, *it ) ); // NEED CONSTRUCTOR FOR MONOMIALS!
//        continue;
//      }
//    }

//    // Otherwise, retain term in the multivariate polynomial
//    auto pt = dec.insert( std::make_pair( pvar, *it ) ); // NEED CONSTRUCTOR FOR MONOMIALS!
//    if( !pt.second )
//      pt.first->second += *it;
//  }
//  
//  return dec;
//}

//inline FFVar
//SPoly::insert
//( FFGraph*dag )
//const
//{
//  FFVar var = 0.;
//  for( auto it=_mapmon.begin(); it!=_mapmon.end(); ++it ){
//    FFVar prodmon = it->second;
//    for( auto ie=it->first.expr.begin(); ie!=it->first.expr.end(); ++ie ){
//      auto itVar = dag->Vars().find( const_cast<FFVar*>(ie->first) );
//      if( itVar == dag->Vars().end() )
//        throw typename SPoly::Exceptions( SPoly::Exceptions::INTERNAL );
//      switch( options.BASIS ){
//       case Options::MONOM:
//        prodmon *= pow( **itVar, (int)ie->second );
//        break;
//       case Options::CHEB:
//        prodmon *= cheb( **itVar, ie->second );
//        break;
//      }
//    }
//    if( it->first.expr.empty() ) var = prodmon;
//    else                         var+= prodmon;
//  }
//  // ADD OPTION TO RETURN A PRODMON OR MONOM TERM?
//  return var;
//}

template <typename KEY, typename COMP>
inline
std::string
SPoly<KEY,COMP>::display
( t_poly mapmon, int const& BASIS, int const& DISPLEN )
const
{
  std::ostringstream out;
  out << std::endl << std::scientific << std::setprecision(DISPLEN);
  for( auto const& [mon,coef] : mapmon )
    out << std::right << std::setw(DISPLEN+7) << coef << "  " << std::setw(2) << mon.tord << "  "
        << mon.display( BASIS ) << std::endl;
  return out.str();
}

template <typename KEY, typename COMP>
inline std::ostream&
operator<<
( std::ostream& out, SPoly<KEY,COMP> const& spoly )
{
  return out << spoly.display( spoly._mapmon, spoly.options.BASIS, spoly.options.DISPLEN );
}

template <typename KEY, typename COMP>
inline SPoly<KEY,COMP>
operator+
( SPoly<KEY,COMP> const& spoly )
{
  return spoly;
}

template <typename KEY, typename COMP>
inline SPoly<KEY,COMP>&
SPoly<KEY,COMP>::operator+=
( SPoly<KEY,COMP> const& spoly )
{
  _setvar.insert( spoly._setvar.begin(), spoly._setvar.end() );
  for( auto&& [mon,coef] : spoly._mapmon ){
    // No warm-start with map::insert unfortunately...
    auto [itmon,ins] = _mapmon.insert( std::make_pair( mon, coef ) );
    if( !ins ) itmon->second += coef;
  }
  _clean( _mapmon );
  return *this;
}

template <typename KEY, typename COMP>
inline SPoly<KEY,COMP>&
SPoly<KEY,COMP>::operator+=
( double const& d )
{
  if( isequal( d, 0. ) ) return *this;
  auto [itmon,ins] = _mapmon.insert( std::make_pair( t_mon(), d ) );
  if( !ins ) itmon->second += d;
  _clean( _mapmon );
  return *this;
}

template <typename KEY, typename COMP>
inline SPoly<KEY,COMP>
operator+
( SPoly<KEY,COMP> const& spoly1, double const d )
{
  SPoly<KEY,COMP> spoly2( spoly1 );
  spoly2 += d;
  return spoly2;
}

template <typename KEY, typename COMP>
inline SPoly<KEY,COMP>
operator+
( double const d, SPoly<KEY,COMP> const& spoly1 )
{
  SPoly<KEY,COMP> spoly2( spoly1 );
  spoly2 += d;
  return spoly2;
}

template <typename KEY, typename COMP>
inline SPoly<KEY,COMP>&
SPoly<KEY,COMP>::operator+=
( std::pair<t_mon,double> const& mon )
{
  // No warm-start for insert unfortunately...
  for( auto && [ivar,iord] : mon.first.expr ) _setvar.insert( ivar );
  auto [itmon,ins] = _mapmon.insert( mon );
  if( !ins ) itmon->second += mon.second;
  _clean( _mapmon );
  return *this;
}

template <typename KEY, typename COMP>
inline SPoly<KEY,COMP>
operator+
( SPoly<KEY,COMP> const& spoly1, SPoly<KEY,COMP> const& spoly2 )
{
  if( spoly1.nmon() >= spoly2.nmon() ){
    SPoly<KEY,COMP> spoly3( spoly1 );
    spoly3 += spoly2;
    return spoly3;
  }
  SPoly<KEY,COMP> spoly3( spoly2 );
  spoly3 += spoly1;
  return spoly3;
}

template <typename KEY, typename COMP>
inline SPoly<KEY,COMP>
operator-
( SPoly<KEY,COMP> const& spoly )
{
  SPoly<KEY,COMP> spoly2( spoly );
  for( auto&& [mon,coef] : spoly2._mapmon ) coef *= -1;
  return spoly2;
}

template <typename KEY, typename COMP>
inline SPoly<KEY,COMP>&
SPoly<KEY,COMP>::operator-=
( double const& d )
{
  if( isequal( d, 0. ) ) return *this;
  auto [itmon,ins] = _mapmon.insert( std::make_pair( t_mon(), -d ) );
  if( !ins ) itmon->second -= d;
  _clean( _mapmon );
  return *this;
}

template <typename KEY, typename COMP>
inline SPoly<KEY,COMP>
operator-
( SPoly<KEY,COMP> const& spoly1, double const d )
{
  SPoly<KEY,COMP> spoly2( spoly1 );
  spoly2 -= d;
  return spoly2;
}

template <typename KEY, typename COMP>
inline SPoly<KEY,COMP>
operator-
( double const d, SPoly<KEY,COMP> const& spoly1 )
{
  SPoly<KEY,COMP> spoly2( -spoly1 );
  spoly2 += d;
  return spoly2;
}

template <typename KEY, typename COMP>
inline SPoly<KEY,COMP>&
SPoly<KEY,COMP>::operator-=
( SPoly<KEY,COMP> const& spoly )
{
  _setvar.insert( spoly._setvar.begin(), spoly._setvar.end() );
  for( auto&& [mon,coef] : spoly._mapmon ){
    // No warm-start with map::insert unfortunately...
    auto [itmon,ins] = _mapmon.insert( std::make_pair( mon, -coef ) );
    if( !ins ) itmon->second -= coef;
  }
  _clean( _mapmon );
  return *this;
}

template <typename KEY, typename COMP>
inline SPoly<KEY,COMP>
operator-
( SPoly<KEY,COMP> const& spoly1, SPoly<KEY,COMP> const& spoly2 )
{
  if( spoly1.nmon() >= spoly2.nmon() ){
    SPoly<KEY,COMP> spoly3( spoly1 );
    spoly3 -= spoly2;
    return spoly3;
  }
  SPoly<KEY,COMP> spoly3( -spoly2 );
  spoly3 += spoly1;
  return spoly3;
}

template <typename KEY, typename COMP>
inline SPoly<KEY,COMP>&
SPoly<KEY,COMP>::operator/=
( double const& d )
{
  if( isequal( d, 0. ) ) throw typename SPoly<KEY,COMP>::Exceptions( SPoly<KEY,COMP>::Exceptions::DIVZERO );
  if( isequal( d, 1. ) ) return *this;
  for( auto&& [mon,coef] : _mapmon ) coef /= d;
  return *this;
}

template <typename KEY, typename COMP>
inline SPoly<KEY,COMP>&
SPoly<KEY,COMP>::operator/=
( SPoly<KEY,COMP> const& spoly )
{
  if( !spoly._mapmon.size() ) throw typename SPoly<KEY,COMP>::Exceptions( SPoly<KEY,COMP>::Exceptions::DIVZERO );
  if( spoly._mapmon.size() == 1 && !spoly._mapmon.cbegin()->first.tord ) return operator/=( spoly._mapmon.cbegin()->second );
  throw typename SPoly<KEY,COMP>::Exceptions( SPoly<KEY,COMP>::Exceptions::DIVPOLY );
}

template <typename KEY, typename COMP>
inline SPoly<KEY,COMP>
operator/
( double const& d, SPoly<KEY,COMP> const& spoly1 )
{
  SPoly<KEY,COMP> spoly2( d );
  return spoly2 /= spoly1;
}

template <typename KEY, typename COMP>
inline SPoly<KEY,COMP>
operator/
( SPoly<KEY,COMP> const& spoly1, SPoly<KEY,COMP> const& spoly2 )
{
  SPoly<KEY,COMP> spoly3( spoly1 );
  return spoly3 /= spoly2;
}

template <typename KEY, typename COMP>
inline SPoly<KEY,COMP>&
SPoly<KEY,COMP>::operator*=
( double const& d )
{
  if( isequal( d, 0. ) ){ *this = 0; return *this; };
  if( isequal( d, 1. ) ) return *this;
  for( auto&& [mon,coef] : _mapmon ) coef *= d;
  return *this;
}

template <typename KEY, typename COMP>
inline SPoly<KEY,COMP>
operator*
( SPoly<KEY,COMP> const& spoly1, double const d )
{
  SPoly<KEY,COMP> spoly2( spoly1 );
  spoly2 *= d;
  return spoly2;
}

template <typename KEY, typename COMP>
inline SPoly<KEY,COMP>
operator*
( double const d, SPoly<KEY,COMP> const& spoly1 )
{
  SPoly<KEY,COMP> spoly2( spoly1 );
  spoly2 *= d;
  return spoly2;
}

template <typename KEY, typename COMP>
inline SPoly<KEY,COMP>&
SPoly<KEY,COMP>::operator*=
( SPoly<KEY,COMP> const& spoly )
{
  // Construct vectors of coefficients for first participating variable
  _setvar.insert( spoly._setvar.begin(), spoly._setvar.end() );
  auto itvar = _setvar.begin();
  std::map<unsigned,t_poly> sp1map, sp2map;
  for( auto&& mon : _mapmon ) _svec1D( itvar, mon, sp1map );
#ifdef MC__SPOLY_DEBUG_SPROD
  _sdisp1D( sp1map, itvar, "Var1: " );
#endif
  for( auto&& mon : spoly._mapmon ) _svec1D( itvar, mon, sp2map );
#ifdef MC__SPOLY_DEBUG_SPROD
  _sdisp1D( sp2map, itvar, "Var2: " );
#endif

  // Call recursive product of univariate Chebyshev polynomials
  _sprod1D( sp1map, sp2map, _mapmon, _setvar, itvar );
  return *this;
}

template <typename KEY, typename COMP>
inline SPoly<KEY,COMP>
operator*
( SPoly<KEY,COMP> const& spoly1, SPoly<KEY,COMP> const& spoly2 )
{
  if( spoly1.nmon() >= spoly2.nmon() ){
    SPoly<KEY,COMP> spoly3( spoly1 );
    spoly3 *= spoly2;
    return spoly3;
  }
  SPoly<KEY,COMP> spoly3( spoly2 );
  spoly3 *= spoly1;
  return spoly3;
}

template <typename KEY, typename COMP>
inline SPoly<KEY,COMP>
sqr
( SPoly<KEY,COMP> const& spoly )
{
 // Construct coefficient maps for first participating variable
  auto itvar = spoly._setvar.begin();
  std::map<unsigned,typename SPoly<KEY,COMP>::t_poly> spmap;
  for( auto&& mon : spoly._mapmon ) spoly._svec1D( itvar, mon, spmap );
#ifdef MC__SPOLY_DEBUG_SQR
  spoly._sdisp1D( spmap, itvar, "Var: " );
#endif

  // Call recursive product of univariate Chebyshev polynomials
  SPoly<KEY,COMP> spsqr;
  spsqr._setvar = spoly._setvar;
  spoly._sprod1D( spmap, spmap, spsqr._mapmon, spsqr._setvar, itvar );
  return spsqr;
}

template <typename KEY, typename COMP>
inline void
SPoly<KEY,COMP>::_sprod1D
( std::map<unsigned,t_poly> const& sp1map,
  std::map<unsigned,t_poly> const& sp2map,
  t_poly& mapmon, t_var const& setvar,
  typename t_var::const_iterator itvar )
const
{
  // construct product matrix of polynomial coefficients
  auto itvarnext = itvar;
  if( !setvar.empty() ) std::advance( itvarnext, 1 );
  std::map<std::pair<unsigned,unsigned>,t_poly> sp12map;
  for( auto& [ndx1,poly1] : sp1map ){
    // empty monomial in sp1
    if( poly1.empty() ) continue; 
    // constant monomial in sp1
    if( poly1.size() == 1 && !poly1.begin()->first.tord ){
      for( auto& [ndx2,poly2] : sp2map ){
        auto [it12,ins] = sp12map.insert( std::make_pair( std::make_pair(ndx1,ndx2), t_poly() ) );
        assert( ins ); // map is initially empty
        _sscal1D( poly2, poly1.begin()->second, it12->second );
      }
      continue;
    }
    // general monomial in sp1
    for( auto& [ndx2,poly2] : sp2map ){
      // empty monomial in sp2
      if( poly2.empty() ) continue; // no term
      // constant monomial in sp2
      if( poly2.size() == 1 && !poly2.begin()->first.tord ){
        auto [it12,ins] = sp12map.insert( std::make_pair( std::make_pair(ndx1,ndx2), t_poly() ) );
        assert( ins ); // map is initially empty
        _sscal1D( poly1, poly2.begin()->second, it12->second );
        continue;
      }
#ifdef MC__SPOLY_DEBUG_SPROD
      std::cout << "Term (" << ndx1 << "," << ndx2 << "):\n";
#endif
      std::map<unsigned,t_poly> sp11map, sp22map;
      for( auto& mon1 : poly1 )
        _svec1D( itvarnext, mon1, sp11map );
#ifdef MC__SPOLYEXPR_DEBUG_SPROD
      _sdisp1D( sp11map, itvarnext, "Poly #1: " );
#endif
      for( auto& mon2 : poly2 )
        _svec1D( itvarnext, mon2, sp22map );
#ifdef MC__SPOLYEXPR_DEBUG_SPROD
      _sdisp1D( sp22map, itvarnext, "Poly #2: " );
#endif
      auto [it12,ins] = sp12map.insert( std::make_pair( std::make_pair(ndx1,ndx2), t_poly() ) );
      assert( ins ); // map is initially empty
      _sprod1D( sp11map, sp22map, it12->second, setvar, itvarnext );
    }
  }

  // construct 1D product result and augment remainder as appropriate
  mapmon.clear();
  for( auto const& [ndx12,poly12] : sp12map ){
    auto const& [ndx1,ndx2] = ndx12;
    // Product involving two constant terms
    if( !ndx1 && !ndx2 )
      _slift1D( poly12, 1., mapmon );
    // Product involving exactly one constant term
    else if( !ndx1 || !ndx2 )
      _slift1D( poly12, 1., mapmon, itvar, ndx1+ndx2 );
    // Product between non-constant terms
    else{
      switch( options.BASIS ){
        // Chebyshev basis functions
        case Options::CHEB:
          _slift1D( poly12, .5, mapmon, itvar, ndx1+ndx2 );
          if( ndx1 == ndx2 )
            _slift1D( poly12, .5, mapmon );
          else if( ndx1 > ndx2 )
            _slift1D( poly12, .5, mapmon, itvar, ndx1-ndx2 );
          else
            _slift1D( poly12, .5, mapmon, itvar, ndx2-ndx1 );
          break;
        // Power basis functions
        case Options::MONOM:
          _slift1D( poly12, 1., mapmon, itvar, ndx1+ndx2 );
          break;
      }
    }
  }

  // clean sparse polynomial
  _clean( mapmon );
#ifdef MC__SPOLY_DEBUG_SPROD
  _sdisp1D( mapmon, "Prod: " );
#endif
}

template <typename KEY, typename COMP>
inline void
SPoly<KEY,COMP>::_sscal1D
( t_poly const& spoly, double const& dscal, t_poly& spscal )

{
  if( isequal(dscal,0.) ) return;
  spscal = spoly;
  if( isequal(dscal,1.) ) return;
  for( auto& [mon,coef] : spscal )
    coef *= dscal;
  return;
}

template <typename KEY, typename COMP>
inline void
SPoly<KEY,COMP>::_slift1D
( t_poly const& spoly, double const& dscal, t_poly& splift )
{
  if( isequal(dscal,0.) ) return;
  for( auto const& [mon,coef] : spoly ){
    auto [itmon,ins] = splift.insert( std::make_pair(mon,coef*dscal) );
    if( !ins ) itmon->second += coef*dscal;
  }
}

template <typename KEY, typename COMP>
inline void
SPoly<KEY,COMP>::_slift1D
( t_poly const& spoly, double const& dscal, t_poly& splift,
 typename t_var::const_iterator itvar, unsigned const ndxord )
{
  for( auto const& [mon,coef] : spoly ){
    auto monlift = mon; // local copy for modification
    monlift.tord += ndxord;
    assert( monlift.expr.insert( std::make_pair( *itvar, ndxord ) ).second );
    auto [itmon,ins] = splift.insert( std::make_pair( monlift, coef*dscal ) );
    if( !ins ) itmon->second += coef*dscal;
  }
}

template <typename KEY, typename COMP>
inline void
SPoly<KEY,COMP>::_svec1D
( typename t_var::const_iterator itvar,
  std::pair<t_mon,double> const& mon,
  std::map<unsigned,t_poly>& mapspoly )
{
  auto const& [ivar,iord] = *mon.first.expr.begin();
  if( !mon.first.tord || ivar != *itvar ) // no dependence on variable *itvar 
    mapspoly[ 0 ].insert( mon );
  else // dependence on variable *itvar of order iord
    mapspoly[ iord ].insert( std::make_pair( t_mon( mon.first.tord - iord,
      std::map<KEY,unsigned,COMP>( ++mon.first.expr.begin(), mon.first.expr.end() ) ),
      mon.second ) );
}

template <typename KEY, typename COMP>
inline void
SPoly<KEY,COMP>::_sdisp1D
( t_poly const& spoly, std::string const& name, std::ostream& os )
const
{
  os << name;
  for( auto itmon=spoly.begin(); itmon!=spoly.end(); ++itmon ){
    if( itmon != spoly.begin() ) os << " + ";
    os << itmon->second;
    for( auto const& [ivar,iord] : itmon->first.expr )
      os << "·" << t_mon( ivar, iord ).display( options.BASIS );
  }
}

template <typename KEY, typename COMP>
inline void
SPoly<KEY,COMP>::_sdisp1D
( std::map<unsigned,t_poly> const& mapspoly, typename t_var::const_iterator itvar, 
  std::string const& name, std::ostream& os )
const
{
  os << name;
  bool first = true;
  for( auto const& [iord,spoly] : mapspoly ){
    if( !first ) os << " + " << t_mon( *itvar, iord ).display( options.BASIS ) << " ·";
    os << " { ";
    _sdisp1D( spoly, "", os );
    os << " }";
    first = false;
  }
  os << std::endl;
}

template <typename KEY, typename COMP>
inline SPoly<KEY,COMP>
pow
( SPoly<KEY,COMP> const& spoly, unsigned const n )
{
  switch( n ){
   case 0:  return 1.;
   case 1:  return spoly;
   case 2:  return sqr( spoly );
   default: return n%2 ? sqr( pow( spoly, n/2 ) ) * spoly : sqr( pow( spoly, n/2 ) );
  }
}

template <typename KEY, typename COMP>
inline SPoly<KEY,COMP>
cheb
( SPoly<KEY,COMP> const& spoly, unsigned const n )
{
  switch( n ){
   case 0:  return 1.;
   case 1:  return spoly;
   case 2:  return 2 * sqr( spoly ) - 1;
   default: return 2 * spoly * cheb( spoly, n-1 ) - cheb( spoly, n-2 );
  }
}

template <typename KEY, typename COMP>
inline SPoly<KEY,COMP>
prod
( unsigned int const npoly, SPoly<KEY,COMP> const* ppoly )
{
  switch( npoly ){
   case 0:  return 1.;
   case 1:  return ppoly[0];
   default: return ppoly[0] * prod( npoly-1, ppoly+1 );
  }
}

template <typename KEY, typename COMP>
inline SPoly<KEY,COMP>
monom
( unsigned int const npoly, SPoly<KEY,COMP> const* ppoly, unsigned const* k,
  bool const chebbasis=false )
{
  switch( npoly ){
   case 0:  return 1.;
   case 1:  return chebbasis? cheb( ppoly[0], k[0] ): pow( ppoly[0], k[0] );
   default: return ( chebbasis? cheb( ppoly[0], k[0] ): pow( ppoly[0], k[0] ) ) * monom( npoly-1, ppoly+1, k+1 );
  }
}

} // namespace mc

#include "mcop.hpp"

namespace mc
{

//! @brief Specialization of the structure mc::Op to allow usage of the type mc::SPoly for DAG evaluation or as a template parameter in other MC++ classes
template <typename KEY, typename COMP> struct Op< mc::SPoly<KEY,COMP> >
{
  typedef mc::SPoly<KEY,COMP> T;
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

