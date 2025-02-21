// Copyright (C) Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the Eclipse Public License.

/*!
\page page_FFEXPR String Expression of Factorable Functions
\author Benoit C. Chachuat
\version 2.0
\date 2024
\bug No known bugs.

mc::FFExpr is a C++ class that constructs strings representing mathematical expressions in factorable functions. It relies on the operator overloading and function overloading mechanisms of C++. The overloaded operators are: `+', `-', `*', and `/'; the overloaded functions include: `exp', `log', `sqr', `pow', `cheb', `sqrt', `fabs', `xlog', `min', `max', `cos', `sin', `tan', `acos', `asin', `atan', `cosh', `sinh', `tanh'.

This class may be used in combination with mc::FFGraph::eval for automatic string representation of subgraphs in DAGs.

\section sec_FFExprOpt Available Options for String Expression

The class mc::FFExpr has a static public member called mc::FFExpr::options that can be used to set/modify the options. For instance:

\code
      FFExpr::options.LANG = FFExpr::Options::GAMS;
\endcode

The available options are the following:

<TABLE border="1">
<CAPTION><EM>Options in mc::FFExpr::Options: name, type and description</EM></CAPTION>
     <TR><TH><b>Name</b>  <TD><b>Type</b><TD><b>Default</b>
         <TD><b>Description</b>
     <TR><TH><tt>LANG</tt> <TD><tt>FFExpr::Options::TYPE</tt> <TD> FFExpr::Options::DAG <TD>Defines the language for the string expressions
</TABLE>


\section sec_FFExprErr Errors Encountered during String Expression

Errors are managed based on the exception handling mechanism of the C++ language. Each time an error is encountered, a class object of type FFExpr::Exceptions is thrown, which contains the type of error. It is the user's responsibility to test whether an exception was thrown during a calculation, and then make the appropriate changes. Should an exception be thrown and not caught by the calling program, the execution will stop.

Possible errors encountered in determining the structure of a factorable function are:

<TABLE border="1">
<CAPTION><EM>Errors during invertible structure detection</EM></CAPTION>
     <TR><TH><b>Number</b> <TD><b>Description</b>
     <TR><TH><tt>-1</tt> <TD>Internal error
</TABLE>
*/

#ifndef MC__FFEXPR_HPP
#define MC__FFEXPR_HPP

#include <iostream>
#include <sstream>
#include <map>

#include "ffunc.hpp"

namespace mc
{

//! @brief C++ class for string expression of factorable functions
////////////////////////////////////////////////////////////////////////
//! mc::FFExpr is a C++ class for string expression of factorable
//! functions
////////////////////////////////////////////////////////////////////////
class FFExpr
////////////////////////////////////////////////////////////////////////
{
private:

  //! @brief Parent precedence. 0: VAR/CST/UNIV, 1:MULT/DIV, 2:ADD/SUB
  unsigned _prec;

  //! @brief String expression
  mutable std::ostringstream _ostr;

  //! @brief Real values
  static std::string _d2s
    ( double const& c )
    { std::ostringstream ostr;
      ostr << std::setprecision(options.DISPLEN);
      ostr << c;
      return ostr.str(); }

public:

  /** @defgroup FFExpr String Expression of Factorable Functions
   *  @{
   */
  //! @brief Exceptions of mc::FFExpr
  class Exceptions
  {
  public:
    //! @brief Enumeration type for FFExpr exception handling
    enum TYPE{
      UNDEF=-33, //!< Internal error
    };
    //! @brief Constructor for error <a>ierr</a>
    Exceptions( TYPE ierr ) : _ierr( ierr ){}

    //! @brief Inline function returning the error flag
    int ierr(){ return _ierr; }
  private:
    TYPE _ierr;
    //! @brief Error description
    std::string what(){
      switch( _ierr ){
      case UNDEF:
        return "mc::FFExpr\t Unavailable feature";
      default:
        return "mc::FFExpr\t Undocumented error";
      }
    }
  };
  
  //! @brief Options of FFExpr
  static struct Options
  {
    //! @brief String expression language type
    enum TYPE{
      DAG=0,  //!< DAG expressions
      GAMS    //!< GAMS language
    };
    //! @brief Constructor
    Options():
      LANG( DAG ), DISPLEN(14)
      {}
    //! @brief Set of allowed invertible operations
    TYPE LANG;
    //! @brief Number of digits in output stream for sparse polynomial coefficients
    unsigned DISPLEN;
  } options;

  //! @brief Default constructor (needed to declare arrays of FFExpr class)
  FFExpr
    ( double const& c=0. )
    : _prec( 0 )
    { if( c > 0 )
        _ostr << _d2s(c);
      else if( c < 0 )
        _ostr << "(" << _d2s(c) << ")"; }

  //! @brief Constructor from FFVar
  FFExpr
    ( FFVar const& X )
    : _prec( 0 )
    { switch( options.LANG ){
       case Options::DAG:
       case Options::GAMS:
         _ostr << X.name();
	 break;
       default:
         throw Exceptions( Exceptions::UNDEF );
      }
    }

  //! @brief Constructor from string
  FFExpr
    ( std::string const& X )
    : _prec( 0 )
    { switch( options.LANG ){
       case Options::DAG:
       case Options::GAMS:
         _ostr << X;
	 break;
       default:
         throw Exceptions( Exceptions::UNDEF );
      }
    }

  //! @brief Copy constructor
  FFExpr
    ( FFExpr const& E )
    : _prec( E._prec )
    { _ostr << E._ostr.str(); }
    
  //! @brief Destructor
  ~FFExpr()
    {}

  //! @brief Set variable to DAG FFVar
  FFExpr& set
    ( FFVar const& X )
    { _prec = 0;
      std::ostringstream otmp;
      _ostr.swap(otmp);
      switch( options.LANG ){
       case Options::DAG:
       case Options::GAMS:
         //std::cout << "FFExpr::set: " << X.name() << std::endl;
         _ostr << X.name();
	 break;
       default:
         throw Exceptions( Exceptions::UNDEF );
      }
      return *this; }

  //! @brief Set variable to constant
  FFExpr& set
    ( double const& c )
    { _prec = 0;
      std::ostringstream otmp;
      _ostr.swap(otmp);
      if( c >= 0 )
        _ostr << _d2s(c);
      else if( c < 0 )
        _ostr << "(" << _d2s(c) << ")";
      return *this; }

  //! @brief Retrieve string expression
  std::ostringstream const& ostr
    ()
    const
    { return _ostr; }

  //! @brief Retrieve or set string expression
  std::ostringstream& ostr
    ()
    { return _ostr; }

  //! @brief Compose string expression
  static FFExpr compose
    ( std::string const& UNIV, FFExpr const& E );
  static FFExpr compose
    ( std::string const& UNIV, FFExpr const& E1, FFExpr const& E2 );
  static FFExpr compose
    ( std::string const& UNIV, unsigned const n, FFExpr const* E );
  static FFExpr compose
    ( std::string const& UNIV, FFExpr const& E, int const n );
  static FFExpr compose
    ( std::string const& UNIV, FFExpr const& E, double const& d );

  //! @brief Retrieve string expression for DAG dependent
  static FFExpr dep
    ( FFVar const& dep );
  //! @brief Retrieve string expression for subgraph
  static std::vector<FFExpr> subgraph
    ( FFGraph* dag, FFSubgraph& sg );
  /** @} */

  // other operator overloadings (inlined)
  FFExpr& operator=
    ( double const& c )
    { return set( c ); }
      
  FFExpr& operator=
    ( FFVar const& X )
    { return set( X ); }
      
  FFExpr& operator=
    ( FFExpr const& E )
    { if( this != &E ){
        _prec = E._prec;
        std::ostringstream otmp;
        _ostr.swap(otmp);
        _ostr << E._ostr.str();
      }
      return *this; }

  FFExpr& operator=
    ( FFExpr&& E )
    { if( this != &E ){
        _prec = E._prec;
        _ostr.swap( E._ostr );
      }
      return *this; }

  FFExpr& operator+=
    ( double const& c )
    { if( c == 0. )        return *this;
      if( !_ostr.tellp() ) return( *this = c );
      _prec = 2;
      if( c > 0 )
        _ostr << " + " << _d2s(c);
      else
        _ostr << " - " << _d2s(-c);
      return *this; }

  FFExpr& operator+=
    ( FFExpr const& E )
    { if( this == &E )       return operator*=( 2 );
      if( !E._ostr.tellp() ) return *this;
      if( !_ostr.tellp() )   return( *this = E );
      _prec = 2;
      _ostr << " + " << E._ostr.str();
      return *this; }

  FFExpr& operator-=
    ( double const& c )
    { if( c == 0. )        return *this;
      if( !_ostr.tellp() ) return( *this = -c );
      _prec = 2;
      if( c > 0 )
        _ostr << " - " << _d2s(c);
      else
        _ostr << " + " << _d2s(-c);
      return *this; }

  FFExpr& operator-=
    ( FFExpr const& E )
    { if( this == &E )       return( *this = 0 );
      if( !E._ostr.tellp() ) return *this;
      _prec = _ostr.tellp()? 2: 0;
      if( !_prec ) _ostr << "(";
      _ostr << " - ";
      if( E._prec > 1 ) _ostr << "( ";
      _ostr << E._ostr.str();
      if( E._prec > 1 ) _ostr << " )";
      if( !_prec ) _ostr << " )";
      return *this; }

  FFExpr& operator*=
    ( double const& c )
    { if( c == 0. || !_ostr.tellp() ) return( *this = 0 );
      if( c == 1. )                   return *this;
      if( _prec > 0 ){
        std::ostringstream otmp;
	otmp << "( " << _ostr.str() << " )";
        _ostr.swap(otmp);
      }
      if( c > 0 )
        _ostr << " * " << _d2s(c);
      else
        _ostr << " * (" << _d2s(c) << ")";
      _prec = 1;
      return *this; }

  FFExpr& operator*=
    ( FFExpr const& E )
    { if( !E._ostr.tellp() || !_ostr.tellp() ) return( *this = 0 );
      if( _prec > 0 ){
        std::ostringstream otmp;
	otmp << "( " << _ostr.str() << " )";
        _ostr.swap(otmp);
      }
      _ostr << " * ";
      if( E._prec > 1 ) _ostr << "( ";
      _ostr << E._ostr.str();
      if( E._prec > 1 ) _ostr << " )";
      _prec = 1;
      return *this; }

  FFExpr& operator/=
    ( double const& c )
    { if( !_ostr.tellp() ) return( *this = 0 );
      if( c == 1. )        return *this;
      if( _prec > 0 ){
        std::ostringstream otmp;
	otmp << "( " << _ostr.str() << " )";
        _ostr.swap(otmp);
      }
      if( c >= 0 )
        _ostr << " / " << _d2s(c);
      else
        _ostr << " / (" << _d2s(c) << ")";
      _prec = 1;
      return *this; }

  FFExpr& operator/=
    ( FFExpr const& E )
    { if( this == &E ) return( *this = 1 );
      if( !E._ostr.tellp() ) return operator/=( 0 );
      if( !_ostr.tellp() )   return( *this = 0 );
      if( _prec > 0 ){
        std::ostringstream otmp;
	otmp << "( " << _ostr.str() << " )";
        _ostr.swap(otmp);
      }
      _ostr << " / ";
      if( E._prec > 0 ) _ostr << "( ";
      _ostr << E._ostr.str();
      if( E._prec > 0 ) _ostr << " )";
      _prec = 1;
      return *this; }
};

////////////////////////////////////////////////////////////////////////

inline FFExpr::Options FFExpr::options;

inline
std::ostream&
operator<<
( std::ostream& out, FFExpr const& E )
{
  return out << E.ostr().str();
}

inline
FFExpr
operator+
( FFExpr const& E )
{
  return E;
}

inline FFExpr
operator-
( FFExpr const& E )
{
  FFExpr _E;
  return _E -= E;
}

inline FFExpr
operator+
( double const& c, FFExpr const& E )
{
  FFExpr _E( E );
  return _E += c;
}

inline FFExpr
operator+
( FFExpr const& E, double const& c )
{
  FFExpr _E( E );
  return _E += c;
}

inline FFExpr
operator+
( FFExpr const& E1, FFExpr const& E2 )
{
  FFExpr _E( E1 );
  return _E += E2;
}

inline FFExpr
sum
( unsigned const n, FFExpr const* E )
{
  switch( n ){
   case 0:  return 0.;
   case 1:  return E[0];
   case 2:  return E[0] + E[1];
   default: return E[0] + sum( n-1, E+1 );
  }
}

inline FFExpr
operator-
( double const& c, FFExpr const& E )
{
  FFExpr _E( c );
  return _E -= E;
}

inline FFExpr
operator-
( FFExpr const& E, double const& c )
{
  FFExpr _E( E );
  return _E -= c;
}

inline FFExpr
operator-
( FFExpr const& E1, FFExpr const& E2 )
{
  FFExpr _E( E1 );
  return _E -= E2;
}

inline FFExpr
operator*
( double const& c, FFExpr const& E )
{
  FFExpr _E( E );
  return _E *= c;
}

inline FFExpr
operator*
( FFExpr const& E, double const& c )
{
  FFExpr _E( E );
  return _E *= c;
}

inline FFExpr
sqr
( FFExpr const& E )
{
  switch( FFExpr::options.LANG ){
   case FFExpr::Options::DAG:
    return FFExpr::compose( FFOp(FFOp::SQR).name(), E );
   case FFExpr::Options::GAMS:
    return FFExpr::compose( "POWER", E, 2 );
   default:
    throw typename FFExpr::Exceptions( FFExpr::Exceptions::UNDEF );
  }
}

inline FFExpr
operator*
( FFExpr const& E1, FFExpr const& E2 )
{
  if( &E1 == &E2 )
    return sqr( E1 );

  FFExpr _E( E1 );
  return _E *= E2;
}

inline
std::vector<FFExpr>
FFExpr::subgraph
( FFGraph* dag, FFSubgraph& sg )
{
  std::vector<FFVar> Var;
  Var.reserve( sg.v_indep.size() );
  std::vector<FFExpr> EVar;
  EVar.reserve( sg.v_indep.size() );
  for( auto const& v : sg.v_indep ){
    Var.push_back( *v );
    EVar.push_back( FFExpr( *v ) );
  }

  std::vector<FFVar> Dep;
  Dep.reserve( sg.v_dep.size() );
  std::vector<FFExpr> EDep( sg.v_dep.size() );
  for( auto const& v : sg.v_dep ){
    Dep.push_back( *v );
  }

  dag->eval( sg, Dep.size(), Dep.data(), EDep.data(), Var.size(), Var.data(), EVar.data() );
  return EDep;
}

inline
FFExpr
FFExpr::dep
( FFVar const& vdep )
{
  if( !vdep.dag() ) return FFExpr();
  FFSubgraph&& sg = vdep.dag()->subgraph( 1, &vdep );
  std::vector<FFVar> vvar;
  vvar.reserve( sg.v_indep.size() );
  std::vector<FFExpr> evar;
  evar.reserve( sg.v_indep.size() );
  for( auto const& v : sg.v_indep ){
    vvar.push_back( *v );
    evar.push_back( FFExpr( *v ) );
  }

  FFExpr edep;
  dynamic_cast<FFGraph*>(vdep.dag())->eval( sg, 1, &vdep, &edep, vvar.size(), vvar.data(), evar.data() );
  return edep;
}

inline
FFExpr
FFExpr::compose
( std::string const& UNIV, FFExpr const& E )
{
  assert( E._ostr.tellp() );
  FFExpr _E; // sets _E._prec = 0 by default
  _E.ostr() << UNIV << "( " << E.ostr().str() << " )";
  return _E;
}

inline
FFExpr
FFExpr::compose
( std::string const& UNIV, FFExpr const& E1, FFExpr const& E2 )
{
  assert( E1._ostr.tellp() );
  assert( E2._ostr.tellp() );
  FFExpr _E; // sets _E._prec = 0 by default
  _E.ostr() << UNIV << "( " << E1.ostr().str() << ", " << E2.ostr().str() << " )";
  return _E;
}

inline
FFExpr
FFExpr::compose
( std::string const& UNIV, unsigned const n, FFExpr const* E )
{
  FFExpr _E; // sets _E._prec = 0 by default
  _E.ostr() << UNIV << "( ";
  for( unsigned i=0; i<n; ++i ){
    assert( E[i]._ostr.tellp() );
    if( i ) _E.ostr() << ", ";
    _E.ostr() << E[i].ostr().str();
  }
  _E.ostr() << " )";
  return _E;
}

inline
FFExpr
FFExpr::compose
( std::string const& UNIV, FFExpr const& E, int const n )
{
  assert( E._ostr.tellp() );
  FFExpr _E; // sets _E._prec = 0 by default
  _E.ostr() << UNIV << "( " << E.ostr().str() << ", " << n << " )";
  return _E;
}

inline
FFExpr
FFExpr::compose
( std::string const& UNIV, FFExpr const& E, double const& d )
{
  assert( E._ostr.tellp() );
  FFExpr _E; // sets _E._prec = 0 by default
  _E.ostr() << UNIV << "( " << E.ostr().str() << ", " << _d2s(d) << " )";
  return _E;
}

inline FFExpr
operator/
( double const& c, FFExpr const& E )
{
  FFExpr _E( c );
  return _E /= E;
}

inline FFExpr
inv
( FFExpr const& E )
{
  switch( FFExpr::options.LANG ){
   case FFExpr::Options::DAG:
    return FFExpr::compose( FFOp(FFOp::INV).name(), E );
   case FFExpr::Options::GAMS:
    return( 1. / E );
   default:
    throw typename FFExpr::Exceptions( FFExpr::Exceptions::UNDEF );
  }
}

inline FFExpr
prod
( unsigned const n, FFExpr const* E )
{
  switch( FFExpr::options.LANG ){
   case FFExpr::Options::DAG:
    return FFExpr::compose( FFOp(FFOp::PROD).name(), n, E );
   case FFExpr::Options::GAMS:
    switch( n ){
     case 0:  return 0.;
     case 1:  return E[0];
     case 2:  return E[0] * E[1];
     default: return E[0] + prod( n-1, E+1 );
    }
   default:
    throw typename FFExpr::Exceptions( FFExpr::Exceptions::UNDEF );
  }
}

inline FFExpr
pow
( FFExpr const& E, const int n )
{
  switch( FFExpr::options.LANG ){
   case FFExpr::Options::DAG:
    return FFExpr::compose( FFOp(FFOp::IPOW).name(), E, n );
   case FFExpr::Options::GAMS:
    if( n == 0  ) return FFExpr( 1 );
    if( n == 1  ) return E;
    if( n == 2  ) return sqr( E );
    if( n == -1 ) return inv( E );
    if( n <  -1 ) return inv( pow( E, -n ) );
    return FFExpr::compose( "POWER", E, n );
   default:
    throw typename FFExpr::Exceptions( FFExpr::Exceptions::UNDEF );
  }
}

inline FFExpr
cheb
( FFExpr const& E, const unsigned n )
{
  switch( FFExpr::options.LANG ){
   case FFExpr::Options::DAG:
    return FFExpr::compose( FFOp(FFOp::CHEB).name(), E, (int)n );
   case FFExpr::Options::GAMS:
    if( n == 0  ) return FFExpr( 1 );
    if( n == 1  ) return E;
    if( n == 2  ) return 2*sqr( E )-1;
    if( n == 3  ) return (4*sqr( E )-3.)*E;
    return 2.*E*cheb(E,n-1)-cheb(E,n-2);
   default:
    throw typename FFExpr::Exceptions( FFExpr::Exceptions::UNDEF );
  }
}

inline FFExpr
monom
( unsigned const n, FFExpr const* E, unsigned const* k )
{
  switch( n ){
   case 0:  return 0.;
   case 1:  return pow( E[0], (int)k[0] );
   default: return pow( E[0], (int)k[0] ) * monom( n-1, E+1, k+1 );
  }
}

inline FFExpr
operator/
( FFExpr const& E, double const& c )
{
  FFExpr _E( E );
  return _E /= c;
}

inline FFExpr
operator/
( FFExpr const& E1, FFExpr const& E2 )
{
  if( &E1 == &E2 )
    return( 1 );

  FFExpr _E( E1 );
  return _E /= E2;
}

inline FFExpr
sqrt
( FFExpr const& E )
{
  switch( FFExpr::options.LANG ){
   case FFExpr::Options::DAG:
    return FFExpr::compose( FFOp(FFOp::SQRT).name(), E );
   case FFExpr::Options::GAMS:
    return FFExpr::compose( "SQRT", E );
   default:
    throw typename FFExpr::Exceptions( FFExpr::Exceptions::UNDEF );
  }
}

inline FFExpr
exp
( FFExpr const& E )
{
  switch( FFExpr::options.LANG ){
   case FFExpr::Options::DAG:
    return FFExpr::compose( FFOp(FFOp::EXP).name(), E );
   case FFExpr::Options::GAMS:
    return FFExpr::compose( "EXP", E );
   default:
    throw typename FFExpr::Exceptions( FFExpr::Exceptions::UNDEF );
  }
}

inline FFExpr
log
( FFExpr const& E )
{
  switch( FFExpr::options.LANG ){
   case FFExpr::Options::DAG:
    return FFExpr::compose( FFOp(FFOp::LOG).name(), E );
   case FFExpr::Options::GAMS:
    return FFExpr::compose( "LOG", E );
   default:
    throw typename FFExpr::Exceptions( FFExpr::Exceptions::UNDEF );
  }
}

inline FFExpr
xlog
( FFExpr const& E )
{
  switch( FFExpr::options.LANG ){
   case FFExpr::Options::DAG:
    return FFExpr::compose( FFOp(FFOp::XLOG).name(), E );
   case FFExpr::Options::GAMS:
    return( - FFExpr::compose( "ENTROPY", E ) );
   default:
    throw typename FFExpr::Exceptions( FFExpr::Exceptions::UNDEF );
  }
}

inline FFExpr
lmtd
( FFExpr const& E1, FFExpr const& E2 )
{
  if( &E1 == &E2 )
    return( E1 );
  return( ( log( E1 ) - log( E2 ) ) / ( E1 - E2 ) );
}

inline FFExpr
rlmtd
( FFExpr const& E1, FFExpr const& E2 )
{
  if( &E1 == &E2 )
    return( 1 / E1 );
  return( ( E1 - E2 ) / ( log( E1 ) - log( E2 ) ) );
}

inline FFExpr
erf
( FFExpr const& E )
{
  switch( FFExpr::options.LANG ){
   case FFExpr::Options::DAG:
    return FFExpr::compose( FFOp(FFOp::ERF).name(), E );
   case FFExpr::Options::GAMS:
    return FFExpr::compose( "ERRORF", E );
   default:
    throw typename FFExpr::Exceptions( FFExpr::Exceptions::UNDEF );
  }
}

inline FFExpr
erfc
( FFExpr const& E )
{
  return( 1 - erf( E ) );
}

inline FFExpr
fabs
( FFExpr const& E )
{
  switch( FFExpr::options.LANG ){
   case FFExpr::Options::DAG:
    return FFExpr::compose( FFOp(FFOp::FABS).name(), E );
   case FFExpr::Options::GAMS:
    return FFExpr::compose( "ABS", E );
   default:
    throw typename FFExpr::Exceptions( FFExpr::Exceptions::UNDEF );
  }
}

inline FFExpr
pow
( FFExpr const& E, double const& r )
{
  switch( FFExpr::options.LANG ){
   case FFExpr::Options::DAG:
    return FFExpr::compose( FFOp(FFOp::DPOW).name(), E, r );
   case FFExpr::Options::GAMS:
    if( r == 0.  ) return FFExpr( 1 );
    if( r == 1.  ) return E;
    if( r == 2.  ) return sqr( E );
    if( r == -1. ) return inv( E );
    return FFExpr::compose( "RPOWER", E, r );
   default:
    throw typename FFExpr::Exceptions( FFExpr::Exceptions::UNDEF );
  }
}

inline FFExpr
pow
( FFExpr const& E1, FFExpr const& E2 )
{
  switch( FFExpr::options.LANG ){
   case FFExpr::Options::DAG:
    return FFExpr::compose( FFOp(FFOp::DPOW).name(), E1, E2 );
   case FFExpr::Options::GAMS:
    return( exp( E2 * log( E1 ) ) );
   default:
    throw typename FFExpr::Exceptions( FFExpr::Exceptions::UNDEF );
  }
}

inline FFExpr
min
( FFExpr const& E1, FFExpr const& E2 )
{
  switch( FFExpr::options.LANG ){
   case FFExpr::Options::DAG:
    return FFExpr::compose( FFOp(FFOp::MINF).name(), E1, E2 );
   case FFExpr::Options::GAMS:
    return FFExpr::compose( "MIN", E1, E2 );
   default:
    throw typename FFExpr::Exceptions( FFExpr::Exceptions::UNDEF );
  }
}

inline FFExpr
max
( FFExpr const& E1, FFExpr const& E2 )
{
  switch( FFExpr::options.LANG ){
   case FFExpr::Options::DAG:
    return FFExpr::compose( FFOp(FFOp::MAXF).name(), E1, E2 );
   case FFExpr::Options::GAMS:
    return FFExpr::compose( "MAX", E1, E2 );
   default:
    throw typename FFExpr::Exceptions( FFExpr::Exceptions::UNDEF );
  }
}

inline FFExpr
min
( const unsigned n, FFExpr const* E )
{
  switch( FFExpr::options.LANG ){
   case FFExpr::Options::DAG:
    return FFExpr::compose( FFOp(FFOp::MINF).name(), n, E );
   case FFExpr::Options::GAMS:
    return FFExpr::compose( "MIN", n, E );
   default:
    throw typename FFExpr::Exceptions( FFExpr::Exceptions::UNDEF );
  }
}

inline FFExpr
max
( const unsigned n, FFExpr const* E )
{
  switch( FFExpr::options.LANG ){
   case FFExpr::Options::DAG:
    return FFExpr::compose( FFOp(FFOp::MAXF).name(), n, E );
   case FFExpr::Options::GAMS:
    return FFExpr::compose( "MAX", n, E );
   default:
    throw typename FFExpr::Exceptions( FFExpr::Exceptions::UNDEF );
  }
}

inline FFExpr
cos
( FFExpr const& E )
{
  switch( FFExpr::options.LANG ){
   case FFExpr::Options::DAG:
    return FFExpr::compose( FFOp(FFOp::COS).name(), E );
   case FFExpr::Options::GAMS:
    return FFExpr::compose( "COS", E );
   default:
    throw typename FFExpr::Exceptions( FFExpr::Exceptions::UNDEF );
  }
}

inline FFExpr
sin
( FFExpr const& E )
{
  switch( FFExpr::options.LANG ){
   case FFExpr::Options::DAG:
    return FFExpr::compose( FFOp(FFOp::SIN).name(), E );
   case FFExpr::Options::GAMS:
    return FFExpr::compose( "SIN", E );
   default:
    throw typename FFExpr::Exceptions( FFExpr::Exceptions::UNDEF );
  }
}

inline FFExpr
tan
( FFExpr const& E )
{
  switch( FFExpr::options.LANG ){
   case FFExpr::Options::DAG:
    return FFExpr::compose( FFOp(FFOp::TAN).name(), E );
   case FFExpr::Options::GAMS:
    return FFExpr::compose( "TAN", E );
   default:
    throw typename FFExpr::Exceptions( FFExpr::Exceptions::UNDEF );
  }
}

inline FFExpr
acos
( FFExpr const& E )
{
  switch( FFExpr::options.LANG ){
   case FFExpr::Options::DAG:
    return FFExpr::compose( FFOp(FFOp::ACOS).name(), E );
   case FFExpr::Options::GAMS:
    return FFExpr::compose( "ARCCOS", E );
   default:
    throw typename FFExpr::Exceptions( FFExpr::Exceptions::UNDEF );
  }
}

inline FFExpr
asin
( FFExpr const& E )
{
  switch( FFExpr::options.LANG ){
   case FFExpr::Options::DAG:
    return FFExpr::compose( FFOp(FFOp::ASIN).name(), E );
   case FFExpr::Options::GAMS:
    return FFExpr::compose( "ARCSIN", E );
   default:
    throw typename FFExpr::Exceptions( FFExpr::Exceptions::UNDEF );
  }
}

inline FFExpr
atan
( FFExpr const& E )
{
  switch( FFExpr::options.LANG ){
   case FFExpr::Options::DAG:
    return FFExpr::compose( FFOp(FFOp::ATAN).name(), E );
   case FFExpr::Options::GAMS:
    return FFExpr::compose( "ARCTAN", E );
   default:
    throw typename FFExpr::Exceptions( FFExpr::Exceptions::UNDEF );
  }
}

inline FFExpr
cosh
( FFExpr const& E )
{
  switch( FFExpr::options.LANG ){
   case FFExpr::Options::DAG:
    return FFExpr::compose( FFOp(FFOp::COSH).name(), E );
   case FFExpr::Options::GAMS:
    return FFExpr::compose( "COSH", E );
   default:
    throw typename FFExpr::Exceptions( FFExpr::Exceptions::UNDEF );
  }
}

inline FFExpr
sinh
( FFExpr const& E )
{
  switch( FFExpr::options.LANG ){
   case FFExpr::Options::DAG:
    return FFExpr::compose( FFOp(FFOp::SINH).name(), E );
   case FFExpr::Options::GAMS:
    return FFExpr::compose( "SINH", E );
   default:
    throw typename FFExpr::Exceptions( FFExpr::Exceptions::UNDEF );
  }
}

inline FFExpr
tanh
( FFExpr const& E )
{
  switch( FFExpr::options.LANG ){
   case FFExpr::Options::DAG:
    return FFExpr::compose( FFOp(FFOp::TANH).name(), E );
   case FFExpr::Options::GAMS:
    return FFExpr::compose( "TANH", E );
   default:
     throw typename FFExpr::Exceptions( FFExpr::Exceptions::UNDEF );
  }
}

inline FFExpr
fstep
( FFExpr const& E )
{
  switch( FFExpr::options.LANG ){
   case FFExpr::Options::DAG:
    return FFExpr::compose( FFOp(FFOp::FSTEP).name(), E );
   case FFExpr::Options::GAMS:
   default:
     throw typename FFExpr::Exceptions( FFExpr::Exceptions::UNDEF );
  }
}

inline FFExpr
bstep
( FFExpr const& E )
{
  return( fstep( -E ) );
}

} // namespace mc

namespace mc
{

//! @brief C++ structure for specialization of the mc::Op templated structure for use of mc::FFEpr in DAG evaluation or as template parameter in other MC++ types
template <> struct Op< mc::FFExpr >
{
  typedef mc::FFExpr FFE;
  static FFE point( const double c ) { return FFE(c); }
  static FFE zeroone() { throw typename FFExpr::Exceptions( FFExpr::Exceptions::UNDEF ); }
  static void I(FFE& x, const FFE&y) { x = y; }
  static double l(const FFE& x) { throw typename FFExpr::Exceptions( FFExpr::Exceptions::UNDEF ); }
  static double u(const FFE& x) { throw typename FFExpr::Exceptions( FFExpr::Exceptions::UNDEF ); }
  static double abs (const FFE& x) { throw typename FFExpr::Exceptions( FFExpr::Exceptions::UNDEF ); }
  static double mid (const FFE& x) { throw typename FFExpr::Exceptions( FFExpr::Exceptions::UNDEF ); }
  static double diam(const FFE& x) { throw typename FFExpr::Exceptions( FFExpr::Exceptions::UNDEF ); }
  static FFE inv (const FFE& x) { return mc::inv(x);  }
  static FFE sqr (const FFE& x) { return mc::sqr(x);  }
  static FFE sqrt(const FFE& x) { return mc::sqrt(x); }
  static FFE exp (const FFE& x) { return mc::exp(x);  }
  static FFE log (const FFE& x) { return mc::log(x);  }
  static FFE xlog(const FFE& x) { return mc::xlog(x); }
  static FFE lmtd(const FFE& x, const FFE& y) { return mc::lmtd(x,y); }
  static FFE rlmtd(const FFE& x, const FFE& y) { return mc::rlmtd(x,y); }
  static FFE fabs(const FFE& x) { return mc::fabs(x); }
  static FFE sin (const FFE& x) { return mc::sin(x);  }
  static FFE cos (const FFE& x) { return mc::cos(x);  }
  static FFE tan (const FFE& x) { return mc::tan(x);  }
  static FFE asin (const FFE& x) { return mc::asin(x);  }
  static FFE acos (const FFE& x) { return mc::acos(x);  }
  static FFE atan (const FFE& x) { return mc::atan(x);  }
  static FFE sinh (const FFE& x) { return mc::sinh(x);  }
  static FFE cosh (const FFE& x) { return mc::cosh(x);  }
  static FFE tanh (const FFE& x) { return mc::tanh(x);  }
  static FFE erf (const FFE& x) { return mc::erf(x);  }
  static FFE erfc (const FFE& x) { return mc::erfc(x);  }
  static FFE fstep(const FFE& x) { return mc::fstep(x); }
  static FFE bstep(const FFE& x) { return mc::bstep(x); }
  static FFE hull(const FFE& x, const FFE& y) { throw typename FFExpr::Exceptions( FFExpr::Exceptions::UNDEF ); }
  static FFE min (const FFE& x, const FFE& y) { return mc::min(x,y); }
  static FFE max (const FFE& x, const FFE& y) { return mc::max(x,y); }
  static FFE arh (const FFE& x, const double k) { throw typename FFExpr::Exceptions( FFExpr::Exceptions::UNDEF ); }
  template <typename X, typename Y> static FFE pow(const X& x, const Y& y) { return mc::pow(x,y); }
  static FFE cheb(const FFE& x, const unsigned n) { return mc::cheb(x,n); }
  static FFE prod (const unsigned n, const FFE* x) { return mc::prod(n,x); }
  static FFE monom (const unsigned n, const FFE* x, const unsigned* k) { return mc::monom(n,x,k); }
  static bool inter(FFE& xIy, const FFE& x, const FFE& y) { throw typename FFExpr::Exceptions( FFExpr::Exceptions::UNDEF ); }
  static bool eq(const FFE& x, const FFE& y) { throw typename FFExpr::Exceptions( FFExpr::Exceptions::UNDEF ); }
  static bool ne(const FFE& x, const FFE& y) { throw typename FFExpr::Exceptions( FFExpr::Exceptions::UNDEF ); }
  static bool lt(const FFE& x, const FFE& y) { throw typename FFExpr::Exceptions( FFExpr::Exceptions::UNDEF ); }
  static bool le(const FFE& x, const FFE& y) { throw typename FFExpr::Exceptions( FFExpr::Exceptions::UNDEF ); }
  static bool gt(const FFE& x, const FFE& y) { throw typename FFExpr::Exceptions( FFExpr::Exceptions::UNDEF ); }
  static bool ge(const FFE& x, const FFE& y) { throw typename FFExpr::Exceptions( FFExpr::Exceptions::UNDEF ); }
};

} // namespace mc

#endif
