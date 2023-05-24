// Copyright (C) 2022 Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the Eclipse Public License.

/*!
\page page_FFINV Invertible Structure Detection in Factorable Functions
\author Benoit C. Chachuat
\version 2.0
\date 2022
\bug No known bugs.

mc::FFInv is a C++ class that determines the structure of mathematical expressions with a view to inverting them, that is, express one of the participating variables as a function of the rest. It relies on the operator overloading and function overloading mechanisms of C++. The overloaded operators are: `+', `-', `*', and `/'; the overloaded functions include: `exp', `log', `sqr', `pow', `cheb', `sqrt', `fabs', `xlog', `min', `max', `cos', `sin', `tan', `acos', `asin', `atan', `cosh', `sinh', `tanh'.

This class may be used in combination with mc::FFDep for detecting the type of dependencies in factorable functions.


\section sec_FFInvEval How Do I Determine the Invertibility of a Factorable Function?

Suppose you are given 4 variables \f$x_0,\ldots,x_3\f$ and want to analyze the structure of the vector function:
\f{eqnarray*}
  {\bf f}({\bf x}) = \left(\begin{array}{c} f_0({\bf x})\\ f_1({\bf x})\end{array}\right) = \left(\begin{array}{c} \displaystyle x_2 x_3+\frac{x_0}{x_2}\\x_0^2[\exp(x_2)+x_3]+x_1 \end{array}\right)
\f}

First, define the variables \f$x_0\ldots,x_3\f$ as

\code
      const int NX = 4;
      mc::FFInv X[NX];
      for( int i=0; i<NX; i++ ) X[i].indep(i);
\endcode

The first line declares <tt>X</tt> as an array of mc::FFInv class objects, then the second line defines X[0] ... X[NX-1] as independent variables with indices 0 ... NX-1, respectively.

Once the independent variables \f${\bf x}\f$ have been defined, the invertible structure of \f${\bf f}({\bf x})\f$ can be detected using the mc::FFInv arithmetic as follows:

\code
      const int NF = 2;
      mc::FFInv F[NF] = { X[2] * X[3] + X[0] / X[2],
                          sqr( X[0] ) * ( exp( X[2] ) + X[3] ) + X[1] };
\endcode

Retrieve information about the invertible structure of \f$f_1\f$ and \f$f_2\f$ using

\code
      auto&& F0_inv = F[0].inv();
      auto&& F1_inv = F[1].inv();
\endcode

or display the invertible structure as

\code
      std::cout << "Invertible structure of F[0]: " << F[0] << std::endl;
      std::cout << "Invertible structure of F[1]: " << F[1] << std::endl;
\endcode

The corresponding output is

\verbatim
      Invertible structure of F[0]: { 0S 2U 3S }
      Invertible structure of F[1]: { 0U 1L 2N 3S }
\endverbatim

These results indicate that X[0], X[2] and X[3] participate in F[0], but not X[1]: X[0] and X[3] are separably linear (S), but the invertibility of X[2] is undetermined (U) due to multiple occurences.

Likewise, all four variables X[0], X[1], X[2] and X[3] participate in F[1]: the invertibility of X[0] is undetermined &ndash; since IPOW is not allowed as invertible in the option set FFInv::Options::INVOP by default); X[1] participates linearly and is thus invertible; X[2] is nonlinearly invertible (N) &ndash; insofar as EXP is indeed allowed as invertible in the option set FFInv::Options::INVOP; and X[3] is separably linear (S).


\section sec_FFInvOpt Available Options for Invertible Structure Detection

The class mc::FFInv has a static public member called mc::FFInv::options that can be used to set/modify the options. For instance:

\code
      FFInv::options.INVOPT = { FFInv::Options::INV, FFInv::Options::EXP };
\endcode

The available options are the following:

<TABLE border="1">
<CAPTION><EM>Options in mc::FFInv::Options: name, type and description</EM></CAPTION>
     <TR><TH><b>Name</b>  <TD><b>Type</b><TD><b>Default</b>
         <TD><b>Description</b>
     <TR><TH><tt>INVOPT</tt> <TD><tt>std::set<NLINV></tt> <TD> { FFInv::Options::INV, FFInv::Options::SQRT, FFInv::Options::EXP, FFInv::Options::LOG, FFInv::Options::RPOW }
         <TD>Defines the set of nonlinear invertible operations - FFInv::Options::IPOW not included by default.
</TABLE>


\section sec_FFInvErr Errors Encountered during Invertible Structure Detection

Errors are managed based on the exception handling mechanism of the C++ language. Each time an error is encountered, a class object of type FFInv::Exceptions is thrown, which contains the type of error. It is the user's responsibility to test whether an exception was thrown during a calculation, and then make the appropriate changes. Should an exception be thrown and not caught by the calling program, the execution will stop.

Possible errors encountered in determining the structure of a factorable function are:

<TABLE border="1">
<CAPTION><EM>Errors during invertible structure detection</EM></CAPTION>
     <TR><TH><b>Number</b> <TD><b>Description</b>
     <TR><TH><tt>-1</tt> <TD>Internal error
     <TR><TH><tt>-33</tt> <TD>Call to unavailable feature
</TABLE>
*/

#ifndef MC__FFINV_HPP
#define MC__FFINV_HPP

#include <iostream>
#include <map>

namespace mc
{

//! @brief C++ class for invertible structure detection in factorable functions
////////////////////////////////////////////////////////////////////////
//! mc::FFInv is a C++ class for invertible structure detection in
//! factorable functions
////////////////////////////////////////////////////////////////////////
class FFInv
////////////////////////////////////////////////////////////////////////
{
  // friends of class FFInv for operator and function overloading
  friend FFInv operator+  ( const FFInv& );
  friend FFInv operator+  ( const FFInv&, const FFInv& );
  friend FFInv operator+  ( const double&, const FFInv& );
  friend FFInv operator+  ( const FFInv&, double const& );
  friend FFInv operator-  ( const FFInv& );
  friend FFInv operator-  ( const FFInv&, const FFInv& );
  friend FFInv operator-  ( const double&, const FFInv& );
  friend FFInv operator-  ( const FFInv&, double const& );
  friend FFInv operator*  ( const FFInv&, const FFInv& );
  friend FFInv operator*  ( const FFInv&, double const& );
  friend FFInv operator*  ( const double&, const FFInv& );
  friend FFInv operator/  ( const FFInv&, const FFInv& );
  friend FFInv operator/  ( const FFInv&, double const& );
  friend FFInv operator/  ( const double&, const FFInv& );
  friend std::ostream& operator<< ( std::ostream&, const FFInv& );
  friend bool operator==  ( const FFInv&, const FFInv& );
  friend bool operator!=  ( const FFInv&, const FFInv& );
  friend bool operator<=  ( const FFInv&, const FFInv& );
  friend bool operator>=  ( const FFInv&, const FFInv& );
  friend bool operator<   ( const FFInv&, const FFInv& );
  friend bool operator>   ( const FFInv&, const FFInv& );
  friend FFInv inv   ( const FFInv& );
  friend FFInv sqr   ( const FFInv& );
  friend FFInv exp   ( const FFInv& );
  friend FFInv log   ( const FFInv& );
  friend FFInv xlog  ( const FFInv& );
  friend FFInv lmtd  ( const FFInv&, const FFInv& );
  friend FFInv rlmtd ( const FFInv&, const FFInv& );
  friend FFInv cos   ( const FFInv& );
  friend FFInv sin   ( const FFInv& );
  friend FFInv tan   ( const FFInv& );
  friend FFInv acos  ( const FFInv& );
  friend FFInv asin  ( const FFInv& );
  friend FFInv atan  ( const FFInv& );
  friend FFInv cosh  ( const FFInv& );
  friend FFInv sinh  ( const FFInv& );
  friend FFInv tanh  ( const FFInv& );
  friend FFInv fabs  ( const FFInv& );
  friend FFInv sqrt  ( const FFInv& );
  friend FFInv erf   ( const FFInv& );
  friend FFInv erfc  ( const FFInv& );
  friend FFInv fstep ( const FFInv& );
  friend FFInv bstep ( const FFInv& );
  friend FFInv cheb  ( const FFInv&, const unsigned );
  friend FFInv pow   ( const FFInv&, const int );
  friend FFInv pow   ( const FFInv&, double const& );
  friend FFInv pow   ( const FFInv&, const FFInv& );
  friend FFInv min   ( const FFInv&, const FFInv& );
  friend FFInv max   ( const FFInv&, const FFInv& );
  friend FFInv inter ( const FFInv&, const FFInv& );
  friend FFInv min   ( const unsigned int, const FFInv* );
  friend FFInv max   ( const unsigned int, const FFInv* );
  friend FFInv sum   ( const unsigned int, const FFInv* );
  friend FFInv prod  ( const unsigned int, const FFInv* );
  friend FFInv monom ( const unsigned int, const FFInv*, const unsigned* );

public:

  /** @defgroup FFInv Invertible Structure Detection for Factorable Functions
   *  @{
   */
  //! @brief Exceptions of mc::FFInv
  class Exceptions
  {
  public:
    //! @brief Enumeration type for FFInv exception handling
    enum TYPE{
      INTERN=-1, //!< Internal error
      UNDEF=-33	  //!< Error due to calling an unavailable feature
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
        return "mc::FFInv\t Unavailable feature";
      case INTERN:
        return "mc::FFInv\t Internal error";
      default:
        return "mc::FFInv\t Undocumented error";
      }
    }
  };
  
  //! @brief Options of FFInv
  static struct Options
  {
    //! @brief Invertible nonlinear operations
    enum NLINV{
      INV=0, //!< Inverse
      SQRT,  //!< Square-root
      EXP,   //!< Exponential
      LOG,   //!< Logarithm
      IPOW,  //!< Integer power
      RPOW   //!< Real power
    };
    //! @brief Constructor
    Options():
      INVOP( {INV,SQRT,EXP,LOG,RPOW} )
      {}
    //! @brief Set of allowed invertible operations
    std::set<NLINV> INVOP;
  } options;

  //! @brief Invertible type
  enum TYPE{
    L=0, //!< Linear
    S,   //!< Separably linear
    N,   //!< Separably nonlinear
    U    //!< Undetermined
  };

  //! @brief Typedef for invertibility map
  typedef std::map<int,TYPE> t_FFINV;

  //! @brief Default constructor (needed to declare arrays of FFInv class)
  FFInv
    ( double const& c=0. )
    {}
  //! @brief Copy constructor
  FFInv
    ( FFInv const& S ):
    _inv( S._inv )
    {}
  //! @brief Destructor
  ~FFInv()
    {}

  //! @brief Set independent variable with index <a>ind</a>
  FFInv& indep
    ( int const ind )
    {
      _inv.clear();
      _inv.insert( std::make_pair( ind, TYPE::L ) );
      return *this;
    }

  //! @brief Determine if current expression is invertible in the variable of index <a>ind</a>
  std::pair<bool,TYPE> inv
    ( int const ind )
    {
      auto it = _inv.find( ind );
      return( it == _inv.end()? std::make_pair( false, TYPE::U )
                              : std::make_pair( true,  it->second ) );
    }

  //! @brief Combine with the invertibility sets of another variable
  FFInv& combine
    ( FFInv const& S, TYPE const& invmin );
  //! @brief Combine the invertibility sets of two variables
  static FFInv combine
    ( FFInv const& S1, FFInv const& S2, TYPE const& invmin );
  //! @brief Update invertibility type of variables
  FFInv& update
    ( TYPE const& invmin );
  //! @brief Copy variables and update invertibility type
  static FFInv copy
    ( FFInv const& S, TYPE const& invmin );

  //! @brief Return invertibility map
  const t_FFINV& inv
    ()
    const
    { return _inv; }
  t_FFINV& inv
    ()
    { return _inv; }
  /** @} */

private:
  //! @brief Invertibility map
  t_FFINV _inv;

  //! @brief Check invertibility of operation op
  static TYPE _nlop
    ( Options::NLINV const& op )
    { return( options.INVOP.find( op ) != options.INVOP.end() ? TYPE::N : TYPE::U ); }

public:
  // other operator overloadings (inlined)
  FFInv& operator=
    ( double const& c )
    { _inv.clear();
      return *this; }
  FFInv& operator=
    ( FFInv const& S )
    { if( this != &S )
        _inv = S._inv;
      return *this; }
  FFInv& operator+=
    ( double const& c )
    { return *this; }
  FFInv& operator+=
    ( FFInv const& S )
    { return combine( S, TYPE::L ); }
  FFInv& operator-=
    ( double const& c )
    { return *this; }
  FFInv& operator-=
    ( FFInv const& S )
    { return combine( S, TYPE::L ); }
  FFInv& operator*=
    ( double const& c )
    { return *this; }
  FFInv& operator*=
    ( FFInv const& S );
  FFInv& operator/=
    ( double const& c )
    { return *this; }
  FFInv& operator/=
    ( FFInv const& S );
};

////////////////////////////////////////////////////////////////////////

inline FFInv::Options FFInv::options;

inline FFInv&
FFInv::update
( TYPE const& invmin )
{
  for( auto& [var,inv] : _inv )
    if( inv < invmin ) inv = invmin;
  return *this;
}

inline FFInv
FFInv::copy
( FFInv const& S, TYPE const& invmin )
{
  FFInv S2( S );
  return( invmin > TYPE::L? S2.update( invmin ): S2 );
}

inline FFInv&
FFInv::combine
( FFInv const& S, TYPE const& invmin )
{
  for( auto const& [varS,invS] : S._inv ){
    auto ins = _inv.insert( std::make_pair(varS,invS) );
    if( ins.second ) continue;
    auto& inv = ins.first->second;
    if( inv <= TYPE::L && invS <= TYPE::L )
      inv = TYPE::L;
    else if( inv <= TYPE::S && invS <= TYPE::S )
      inv = TYPE::S;
    else
      inv = TYPE::U;
  }

  return( invmin > TYPE::L? update( invmin ): *this );
}

inline FFInv
FFInv::combine
( FFInv const& S1, FFInv const& S2, TYPE const& invmin )
{
  FFInv S3( S1 );
  return S3.combine( S2, invmin );
}

inline std::ostream&
operator<<
( std::ostream& out, FFInv const& S )
{
  out << "{ ";
  auto iti = S._inv.begin();
  for( ; iti != S._inv.end(); ++iti ){
    out << iti->first;
    switch( iti->second ){
     case FFInv::TYPE::L: out << "L"; break;
     case FFInv::TYPE::S: out << "S"; break;
     case FFInv::TYPE::N: out << "N"; break;
     case FFInv::TYPE::U: out << "U"; break;
    }
    out << " ";
  }
  out << "}";
  return out;
}

inline FFInv
operator+
( FFInv const& S )
{
  return S;
}

inline FFInv
operator-
( FFInv const& S )
{
  return S;
}

inline FFInv
operator+
( double const& c, FFInv const& S )
{
  return S;
}

inline FFInv
operator+
( FFInv const& S, double const& c )
{
  return S;
}

inline FFInv
operator+
( FFInv const& S1, FFInv const& S2 )
{
  if( S1._inv.empty() ) return S2;
  if( S2._inv.empty() ) return S1;
  return FFInv::combine( S1, S2, FFInv::TYPE::L );
}

inline FFInv
sum
( const unsigned int n, const FFInv*S )
{
  switch( n ){
   case 0:  return 0.;
   case 1:  return S[0];
   case 2:  return S[0] + S[1];
   default: return S[0] + sum( n-1, S+1 );
  }
}

inline FFInv
operator-
( double const& c, FFInv const& S )
{
  return S;
}

inline FFInv
operator-
( FFInv const& S, double const& c )
{
  return S;
}

inline FFInv
operator-
( FFInv const& S1, FFInv const& S2 )
{
  if( S1._inv.empty() ) return S2;
  if( S2._inv.empty() ) return S1;
  return FFInv::combine( S1, S2, FFInv::TYPE::L );
}

inline FFInv
operator*
( double const& c, FFInv const& S )
{
  return S;
}

inline FFInv
operator*
( FFInv const& S, double const& c )
{
  return S;
}

inline FFInv&
FFInv::operator*=
( FFInv const& S )
{
  // Don't alter *this if S has no dependencies
  if( S._inv.empty() )
    return *this;

  // Copy S if *this has no dependencies
  if( _inv.empty() )
    return combine( S, TYPE::L );

  // Set invertibility rules in product term
  for( auto const& [varS,invS] : S._inv ){
    auto ins = _inv.insert( std::make_pair(varS,invS) );
    if( ins.second ) continue;
    ins.first->second = TYPE::U;
  }
  for( auto& [var,inv] : _inv )
    if( inv < TYPE::S ) inv = TYPE::S;

  return *this;
}

inline FFInv
operator*
( FFInv const& S1, FFInv const& S2 )
{
  if( &S1 == &S2 )
    return sqr( S1 );

  FFInv S3( S1 );
  return S3 *= S2 ;
}

inline FFInv
sqr
( FFInv const& S )
{
  return FFInv::copy( S, FFInv::_nlop( FFInv::Options::IPOW ) );
//  
//  // Don't alter *this if S has no dependencies
//  if( S._inv.empty() )
//    return S;
//  FFInv S2( S );

//  // Set invertibility rules in squared term
//  for( auto& [varS2,invS2] : S2._inv )
//    invS2 = FFInv::_nlop( FFInv::Options::IPOW );
//  return( invmin > TYPE::L? S2.update( invmin ): S2 );
//  
//  return S2;
}

inline FFInv
prod
( const unsigned int n, const FFInv*S )
{
  switch( n ){
   case 0:  return 0.;
   case 1:  return S[0];
   case 2:  return S[0] * S[1];
   default: return S[0] * prod( n-1, S+1 );
  }
}

inline FFInv
monom
( const unsigned int n, const FFInv*S, const unsigned*k )
{
  switch( n ){
   case 0:  return 0.;
   case 1:  return pow( S[0], (int)k[0] );
   default: return pow( S[0], (int)k[0] ) * monom( n-1, S+1, k+1 );
  }
}

inline FFInv
inv
( FFInv const& S )
{
  return FFInv::copy( S, FFInv::_nlop( FFInv::Options::INV ) );
}


inline FFInv
operator/
( FFInv const& S, double const& c )
{
  return S;
}

inline FFInv
operator/
( double const& c, FFInv const& S )
{
  return mc::inv( S );
}

inline FFInv&
FFInv::operator/=
( FFInv const& S )
{
  if( S._inv.empty() ) return *this;
  return operator*=( mc::inv(S) );
}

inline FFInv
operator/
( FFInv const& S1, FFInv const& S2 )
{
  FFInv S3( S1 );
  return S3 /= S2 ;
}

inline FFInv
sqrt
( FFInv const& S )
{
  return FFInv::copy( S, FFInv::_nlop( FFInv::Options::SQRT ) );
}

inline FFInv
exp
( FFInv const& S )
{
  return FFInv::copy( S, FFInv::_nlop( FFInv::Options::EXP ) );
}

inline FFInv
arh
( FFInv const& S, double const& a )
{
  return FFInv::copy( S, FFInv::TYPE::U );
}

inline FFInv
log
( FFInv const& S )
{
  return FFInv::copy( S, FFInv::_nlop( FFInv::Options::LOG ) );
}

inline FFInv
xlog
( FFInv const& S )
{
  return FFInv::copy( S, FFInv::TYPE::U );
}

inline FFInv
lmtd
( FFInv const& S1, FFInv const& S2 )
{
  return FFInv::combine( S1, S2, FFInv::TYPE::U );
}

inline FFInv
rlmtd
( FFInv const& S1, FFInv const& S2 )
{
  return FFInv::combine( S1, S2, FFInv::TYPE::U );
}

inline FFInv
erf
( FFInv const& S )
{
  return FFInv::copy( S, FFInv::TYPE::U );
}

inline FFInv
erfc
( FFInv const& S )
{
  return FFInv::copy( S, FFInv::TYPE::U );
}

inline FFInv
fstep
( FFInv const& S )
{
  return FFInv::copy( S, FFInv::TYPE::U );
}

inline FFInv
bstep
( FFInv const& S )
{
  return FFInv::copy( S, FFInv::TYPE::U );
}

inline FFInv
fabs
( FFInv const& S )
{
  return FFInv::copy( S, FFInv::TYPE::U );
}

inline FFInv
cheb
( FFInv const& S, const unsigned n )
{
  if( n == 0 ) return FFInv();
  if( n == 1 ) return S;
  if( n == 2 ) return sqr(S);
  return FFInv::copy( S, FFInv::TYPE::U );
}

inline FFInv
pow
( FFInv const& S, const int n )
{
  if( n == 0  ) return FFInv();
  if( n == 1  ) return S;
  if( n == 2  ) return sqr(S);
  if( n == -1 ) return inv(S);
  if( n <  -1 ) return inv( pow(S,-n) );
  return FFInv::copy( S, FFInv::_nlop( FFInv::Options::IPOW ) );
}

inline FFInv
pow
( FFInv const& S, double const& a )
{
  if( a == 0.  ) return FFInv();
  if( a == 1.  ) return S;
  if( a == 2.  ) return sqr(S);
  if( a == -1. ) return inv(S);
  if( a == std::rint(a) && a > 2.  ) return FFInv::copy( S, FFInv::TYPE::U );
  if( a == std::rint(a) && a < -1. ) return FFInv::copy( S, FFInv::TYPE::U );
  return FFInv::copy( S, FFInv::_nlop( FFInv::Options::RPOW ) );
}

inline FFInv
pow
( FFInv const& S1, FFInv const& S2 )
{
  return FFInv::combine( S1, S2, FFInv::TYPE::U );
}

inline FFInv
min
( FFInv const& S1, FFInv const& S2 )
{
  return FFInv::combine( S1, S2, FFInv::TYPE::U );
}

inline FFInv
max
( FFInv const& S1, FFInv const& S2 )
{
  return FFInv::combine( S1, S2, FFInv::TYPE::U );
}

inline FFInv
min
( const unsigned int n, const FFInv*S )
{
  switch( n ){
   case 0:  return FFInv();
   case 1:  return S[0];
   case 2:  return min( S[0], S[1] );
   default: return min( S[0], min( n-1, S+1 ) );
  }
}

inline FFInv
max
( const unsigned int n, const FFInv*S )
{
  switch( n ){
   case 0:  return FFInv();
   case 1:  return S[0];
   case 2:  return max( S[0], S[1] );
   default: return max( S[0], max( n-1, S+1 ) );
  }
}

inline FFInv
cos
( FFInv const& S )
{
  return FFInv::copy( S, FFInv::TYPE::U );
}

inline FFInv
sin
( FFInv const& S )
{
  return FFInv::copy( S, FFInv::TYPE::U );
}

inline FFInv
tan
( FFInv const& S )
{
  return FFInv::copy( S, FFInv::TYPE::U );
}

inline FFInv
acos
( FFInv const& S )
{
  return FFInv::copy( S, FFInv::TYPE::U );
}

inline FFInv
asin
( FFInv const& S )
{
  return FFInv::copy( S, FFInv::TYPE::U );
}

inline FFInv
atan
( FFInv const& S )
{
  return FFInv::copy( S, FFInv::TYPE::U );
}

inline FFInv
cosh
( FFInv const& S )
{
  return FFInv::copy( S, FFInv::TYPE::U );
}

inline FFInv
sinh
( FFInv const& S )
{
  return FFInv::copy( S, FFInv::TYPE::U );
}

inline FFInv
tanh
( FFInv const& S )
{
  return FFInv::copy( S, FFInv::TYPE::U );
}

inline FFInv
inter
( FFInv const& S1, FFInv const& S2 )
{
  return FFInv::combine( S1, S2, FFInv::TYPE::U );
}

inline bool
operator==
( FFInv const& S1, FFInv const& S2 )
{
  return( S1._inv == S2._inv );
}

inline bool
operator!=
( FFInv const& S1, FFInv const& S2 )
{
  return( S1._inv != S2._inv );
}

inline bool
operator<=
( FFInv const& S1, FFInv const& S2 )
{
  return( S1._inv <= S2._inv );
}

inline bool
operator>=
( FFInv const& S1, FFInv const& S2 )
{
  return( S1._inv >= S2._inv );
}

inline bool
operator<
( FFInv const& S1, FFInv const& S2 )
{
  return( S1._inv < S2._inv );
}

inline bool
operator>
( FFInv const& S1, FFInv const& S2 )
{
  return( S1._inv > S2._inv );
}

} // namespace mc

namespace mc
{

//! @brief Specialization of the structure mc::Op to allow usage of the type mc::Interval for DAG evaluation or as a template parameter in other MC++ classes
template <> struct Op< mc::FFInv >
{
  typedef mc::FFInv FV;
  static FV point( const double c ) { return FV(c); }
  static FV zeroone() { throw typename FFInv::Exceptions( FFInv::Exceptions::UNDEF ); }
  static void I(FV& x, const FV&y) { x = y; }
  static double l(const FV& x) { throw typename FFInv::Exceptions( FFInv::Exceptions::UNDEF ); }
  static double u(const FV& x) { throw typename FFInv::Exceptions( FFInv::Exceptions::UNDEF ); }
  static double abs (const FV& x) { throw typename FFInv::Exceptions( FFInv::Exceptions::UNDEF );  }
  static double mid (const FV& x) { throw typename FFInv::Exceptions( FFInv::Exceptions::UNDEF );  }
  static double diam(const FV& x) { throw typename FFInv::Exceptions( FFInv::Exceptions::UNDEF ); }
  static FV inv (const FV& x) { return mc::inv(x);  }
  static FV sqr (const FV& x) { return mc::sqr(x);  }
  static FV sqrt(const FV& x) { return mc::sqrt(x); }
  static FV exp (const FV& x) { return mc::exp(x);  }
  static FV log (const FV& x) { return mc::log(x);  }
  static FV xlog(const FV& x) { return mc::xlog(x); }
  static FV lmtd(const FV& x, const FV& y) { return mc::lmtd(x,y); }
  static FV rlmtd(const FV& x, const FV& y) { return mc::rlmtd(x,y); }
  static FV fabs(const FV& x) { return mc::fabs(x); }
  static FV sin (const FV& x) { return mc::sin(x);  }
  static FV cos (const FV& x) { return mc::cos(x);  }
  static FV tan (const FV& x) { return mc::tan(x);  }
  static FV asin(const FV& x) { return mc::asin(x); }
  static FV acos(const FV& x) { return mc::acos(x); }
  static FV atan(const FV& x) { return mc::atan(x); }
  static FV sinh(const FV& x) { return mc::sinh(x); }
  static FV cosh(const FV& x) { return mc::cosh(x); }
  static FV tanh(const FV& x) { return mc::tanh(x); }
  static FV erf (const FV& x) { return mc::erf(x);  }
  static FV erfc(const FV& x) { return mc::erfc(x); }
  static FV fstep(const FV& x) { return mc::fstep(x); }
  static FV bstep(const FV& x) { return mc::bstep(x); }
  static FV hull(const FV& x, const FV& y) { throw typename FFInv::Exceptions( FFInv::Exceptions::UNDEF ); }
  static FV min (const FV& x, const FV& y) { return mc::min(x,y);  }
  static FV max (const FV& x, const FV& y) { return mc::max(x,y);  }
  static FV arh (const FV& x, const double k) { return mc::exp(-k/x); }
  template <typename X, typename Y> static FV pow(const X& x, const Y& y) { return mc::pow(x,y); }
  static FV cheb(const FV& x, const unsigned n) { return mc::cheb(x,n); }
  static FV prod(const unsigned int n, const FV* x) { return mc::prod(n,x); }
  static FV monom(const unsigned int n, const FV* x, const unsigned* k) { return mc::monom(n,x,k); }
  static bool inter(FV& xIy, const FV& x, const FV& y) { xIy = mc::inter(x,y); return true; }
  static bool eq(const FV& x, const FV& y) { throw typename FFInv::Exceptions( FFInv::Exceptions::UNDEF ); }
  static bool ne(const FV& x, const FV& y) { throw typename FFInv::Exceptions( FFInv::Exceptions::UNDEF ); }
  static bool lt(const FV& x, const FV& y) { throw typename FFInv::Exceptions( FFInv::Exceptions::UNDEF ); }
  static bool le(const FV& x, const FV& y) { throw typename FFInv::Exceptions( FFInv::Exceptions::UNDEF ); }
  static bool gt(const FV& x, const FV& y) { throw typename FFInv::Exceptions( FFInv::Exceptions::UNDEF ); }
  static bool ge(const FV& x, const FV& y) { throw typename FFInv::Exceptions( FFInv::Exceptions::UNDEF ); }
};

} // namespace mc

#endif
