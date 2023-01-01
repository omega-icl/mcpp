// Copyright (C) 2017 Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the Eclipse Public License.

/*!
\page page_FFDEP Dependence Structure Detection for Factorable Functions
\author Benoit C. Chachuat
\version 1.0
\date 2017
\bug No known bugs.

mc::FFDep is a C++ class that determines the structure of mathematical expressions, namely their sparsity pattern and linearity, for a given set of participating variables. It relies on the operator overloading and function overloading mechanisms of C++. The overloaded operators are: `+', `-', `*', and `/'; the overloaded functions include: `exp', `log', `sqr', `pow', `cheb', `sqrt', `fabs', `xlog', `min', `max', `cos', `sin', `tan', `acos', `asin', `atan', `cosh', `sinh', `tanh'.


\section sec_FFDepEval How Do I Determine the Structure of a Factorable Function?

Suppose you are given 4 variables \f$x_0,\ldots,x_3\f$ and want to determine the sparsity pattern and linearity of the vector following function
\f{eqnarray*}
  {\bf f}({\bf x}) = \left(\begin{array}{c} f_0({\bf x})\\ f_1({\bf x})\end{array}\right) = \left(\begin{array}{c} \displaystyle x_2 x_3+\frac{x_0}{x_2}\\x_0^2[\exp(x_2)+x_3]+x_1 \end{array}\right)
\f}

First, define the variables \f$x_0\ldots,x_3\f$ as

\code
      const int NX = 4;
      mc::FFDep X[NX];
      for( int i=0; i<NX; i++ ) X[i].indep(i);
\endcode

Essentially, the first line means that <tt>X</tt> is an array of mc::FFDep class objects, and the second line defines X[0] ... X[NX-1] as independent variables with indices 0 ... NX-1, respectively.

Once the independent variables \f${\bf x}\f$ have been defined, determine the structure of \f${\bf f}({\bf x})\f$ simply as

\code
      const int NF = 2;
      mc::FFDep F[NF] = { X[2] * X[3] + X[0] / X[2],
                          sqr( X[0] ) * ( exp( X[2] ) + X[3] ) + X[1] };
\endcode

Retrieve the structure - both the sparsity pattern and the dependence type - of \f$f_1\f$ and \f$f_2\f$ as

\code
      auto&& F0_dep = F[0].dep();
      auto&& F1_dep = F[1].dep();
\endcode

You can also display the structure as

\code
      std::cout << "Variable dependence of F[0]: " << F[0] << std::endl;
      std::cout << "Variable dependence of F[1]: " << F[1] << std::endl;
\endcode

The corresponding output is

\verbatim
      Variable dependence of F[0]: { 0R 2R 3Q }
      Variable dependence of F[1]: { 0N 1L 2N 3P }
\endverbatim

which indicates that X[0], X[2] and X[3] participate in F[0], but not X[1], and that X[0] and X[2] participate in rational 'R' terms, whereas X[3] participates in quadratic 'Q' terms. Likewise, all four variables X[0], X[1], X[2] and X[3] participate in F[1], with X[0] and X[2] in nonlinea 'N' terms, X[3] in polynomial 'P' terms. and X[1] in linear terms.


\section sec_FFDepErr Errors Encountered in Determining the Structure of a Factorable Function?

Errors are managed based on the exception handling mechanism of the C++ language. Each time an error is encountered, a class object of type FFDep::Exceptions is thrown, which contains the type of error. It is the user's responsibility to test whether an exception was thrown during a calculation, and then make the appropriate changes. Should an exception be thrown and not caught by the calling program, the execution will stop.

Possible errors encountered in determining the structure of a factorable function are:

<TABLE border="1">
<CAPTION><EM>Errors during Structure Determination</EM></CAPTION>
     <TR><TH><b>Number</b> <TD><b>Description</b>
     <TR><TH><tt>-1</tt> <TD>Internal error
     <TR><TH><tt>-33</tt> <TD>Call to unavailable feature
</TABLE>
*/

#ifndef MC__FFDEP_HPP
#define MC__FFDEP_HPP

#include <iostream>
#include <map>

namespace mc
{

//! @brief C++ class for evaluation of the sparsity pattern of a factorable function
////////////////////////////////////////////////////////////////////////
//! mc::FFDep is a C++ class for evaluating the sparsity pattern of a
//! factorable function
////////////////////////////////////////////////////////////////////////
class FFDep
////////////////////////////////////////////////////////////////////////
{
  // friends of class FFDep for operator and function overloading
  friend FFDep operator+  ( const FFDep& );
  friend FFDep operator+  ( const FFDep&, const FFDep& );
  friend FFDep operator+  ( const double, const FFDep& );
  friend FFDep operator+  ( const FFDep&, const double );
  friend FFDep operator-  ( const FFDep& );
  friend FFDep operator-  ( const FFDep&, const FFDep& );
  friend FFDep operator-  ( const double, const FFDep& );
  friend FFDep operator-  ( const FFDep&, const double );
  friend FFDep operator*  ( const FFDep&, const FFDep& );
  friend FFDep operator*  ( const FFDep&, const double );
  friend FFDep operator*  ( const double, const FFDep& );
  friend FFDep operator/  ( const FFDep&, const FFDep& );
  friend FFDep operator/  ( const FFDep&, const double );
  friend FFDep operator/  ( const double, const FFDep& );
  friend std::ostream& operator<< ( std::ostream&, const FFDep& );
  friend bool operator==  ( const FFDep&, const FFDep& );
  friend bool operator!=  ( const FFDep&, const FFDep& );
  friend bool operator<=  ( const FFDep&, const FFDep& );
  friend bool operator>=  ( const FFDep&, const FFDep& );
  friend bool operator<   ( const FFDep&, const FFDep& );
  friend bool operator>   ( const FFDep&, const FFDep& );
  friend FFDep inv   ( const FFDep& );
  friend FFDep sqr   ( const FFDep& );
  friend FFDep exp   ( const FFDep& );
  friend FFDep log   ( const FFDep& );
  friend FFDep xlog  ( const FFDep& );
  friend FFDep lmtd  ( const FFDep&, const FFDep& );
  friend FFDep rlmtd ( const FFDep&, const FFDep& );
  friend FFDep cos   ( const FFDep& );
  friend FFDep sin   ( const FFDep& );
  friend FFDep tan   ( const FFDep& );
  friend FFDep acos  ( const FFDep& );
  friend FFDep asin  ( const FFDep& );
  friend FFDep atan  ( const FFDep& );
  friend FFDep cosh  ( const FFDep& );
  friend FFDep sinh  ( const FFDep& );
  friend FFDep tanh  ( const FFDep& );
  friend FFDep fabs  ( const FFDep& );
  friend FFDep sqrt  ( const FFDep& );
  friend FFDep erf   ( const FFDep& );
  friend FFDep erfc  ( const FFDep& );
  friend FFDep fstep ( const FFDep& );
  friend FFDep bstep ( const FFDep& );
  friend FFDep cheb  ( const FFDep&, const unsigned );
  friend FFDep pow   ( const FFDep&, const int );
  friend FFDep pow   ( const FFDep&, const double );
  friend FFDep pow   ( const FFDep&, const FFDep& );
  friend FFDep min   ( const FFDep&, const FFDep& );
  friend FFDep max   ( const FFDep&, const FFDep& );
  friend FFDep inter ( const FFDep&, const FFDep& );
  friend FFDep min   ( const unsigned int, const FFDep* );
  friend FFDep max   ( const unsigned int, const FFDep* );
  friend FFDep sum   ( const unsigned int, const FFDep* );
  friend FFDep prod  ( const unsigned int, const FFDep* );
  friend FFDep monom ( const unsigned int, const FFDep*, const unsigned* );

public:

  /** @defgroup FFDep Dependence Structure Detection for Factorable Functions
   *  @{
   */
  //! @brief Exceptions of mc::FFDep
  class Exceptions
  {
  public:
    //! @brief Enumeration type for FFDep exception handling
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
        return "mc::FFDep\t Unavailable feature";
      case INTERN:
        return "mc::FFDep\t Internal error";
      default:
        return "mc::FFDep\t Undocumented error";
      }
    }
  };

  //! @brief Dependence type
  enum TYPE{
    L=0, //!< Linear
    Q,   //!< Quadratic
    P,   //!< Polynomial
    R,   //!< Rational
    N    //!< General nonlinear
  };

  //! @brief Typedef for dependency map
  typedef std::map<int,TYPE> t_FFDep;

  //! @brief Default constructor (needed to declare arrays of FFDep class)
  FFDep
    ( const double c=0. )
    {}
  //! @brief Copy constructor
  FFDep
    ( const FFDep&S ):
    _dep(S._dep)
    {}
  //! @brief Destructor
  ~FFDep()
    {}

  //! @brief Set independent variable with index <a>ind</a>
  FFDep& indep
    ( const int ind )
    { _dep.clear();
      _dep.insert( std::make_pair(ind,TYPE::L) );
      return *this; }

  //! @brief Determine if current object is dependent on the variable of index <a>ind</a>
  std::pair<bool,TYPE> dep
    ( const int ind )
    { auto it = _dep.find(ind);
      return( it==_dep.end()? std::make_pair(false,TYPE::L):
                              std::make_pair(true,it->second) ); }

  //! @brief Combines with the dependency sets of another variable
  FFDep& combine
    ( const FFDep&S, const TYPE&dep );
  //! @brief Combines the dependency sets of two variables
  static FFDep combine
    ( const FFDep&S1, const FFDep&S2, const TYPE&dep );
  //! @brief Update type of dependent variables
  FFDep& update
    ( const TYPE&dep );
  //! @brief Copy dependent variables and update type
  static FFDep copy
    ( const FFDep&S, const TYPE&dep );
  //! @brief Return worst dependency
  TYPE worst() const;

  //! @brief Return dependency map
  const t_FFDep& dep() const
    { return _dep; }
  t_FFDep& dep()
    { return _dep; }
  /** @} */
  
private:
  //! @brief Dependency set
  t_FFDep _dep;

public:
  // other operator overloadings (inlined)
  FFDep& operator=
    ( const double c )
    { _dep.clear();
      return *this; }
  FFDep& operator=
    ( const FFDep&S )
    { if( this != &S ) _dep = S._dep;
      return *this; }
  FFDep& operator+=
    ( const double c )
    { return *this; }
  FFDep& operator+=
    ( const FFDep&S )
    { return combine( S, TYPE::L ); }
  FFDep& operator-=
    ( const double c )
    { return *this; }
  FFDep& operator-=
    ( const FFDep&S )
    { return combine( S, TYPE::L ); }
  FFDep& operator*=
    ( const double c )
    { return *this; }
  FFDep& operator*=
    ( const FFDep&S );
  FFDep& operator/=
    ( const double c )
    { return *this; }
  FFDep& operator/=
    ( const FFDep&S );
};

////////////////////////////////////////////////////////////////////////

inline FFDep::TYPE
FFDep::worst
() const
{
  auto it = _dep.begin();
  TYPE depw = TYPE::L;
  for( ; it != _dep.end(); ++it )
    if( it->second > depw ) depw = it->second;
  return depw;
}

inline FFDep&
FFDep::update
( const TYPE&dep )
{
  auto it = _dep.begin();
  for( ; it != _dep.end(); ++it )
    if( it->second < dep ) it->second = dep;
  return *this;
}

inline FFDep
FFDep::copy
( const FFDep&S, const TYPE&dep )
{
  FFDep S2( S );
  return S2.update( dep ); 
}

inline FFDep&
FFDep::combine
( const FFDep&S, const TYPE&dep )
{
  auto cit = S._dep.begin();
  for( ; cit != S._dep.end(); ++cit ){
    auto ins = _dep.insert( *cit );
    if( !ins.second && ins.first->second < cit->second ) ins.first->second = cit->second;
  }
  return( dep? update( dep ): *this );
}

inline FFDep
FFDep::combine
( const FFDep&S1, const FFDep&S2, const TYPE&dep )
{
  FFDep S3( S1 );
  return S3.combine( S2, dep );
}

inline std::ostream&
operator<<
( std::ostream&out, const FFDep&S)
{
  out << "{ ";
  auto iS = S._dep.begin();
  for( ; iS != S._dep.end(); ++iS ){
    out << iS->first;
    switch( iS->second ){
     case FFDep::TYPE::L: out << "L "; break;
     case FFDep::TYPE::Q: out << "Q "; break;
     case FFDep::TYPE::P: out << "P "; break;
     case FFDep::TYPE::R: out << "R "; break;
     case FFDep::TYPE::N: out << "N "; break;
    }
  }
  out << "}";
  return out;
}

inline FFDep
operator+
( const FFDep&S )
{
  return S;
}

inline FFDep
operator-
( const FFDep&S )
{
  return S;
}

inline FFDep
operator+
( const double c, const FFDep&S )
{
  return S;
}

inline FFDep
operator+
( const FFDep&S, const double c )
{
  return S;
}

inline FFDep
operator+
( const FFDep&S1, const FFDep&S2 )
{ 
  if( S1._dep.empty() ) return S2;
  if( S2._dep.empty() ) return S1;
  return FFDep::combine( S1, S2, FFDep::TYPE::L );
}

inline FFDep
sum
( const unsigned int n, const FFDep*S )
{
  switch( n ){
   case 0:  return 0.;
   case 1:  return S[0];
   case 2:  return S[0] + S[1];
   default: return S[0] + sum( n-1, S+1 );
  }
}

inline FFDep
operator-
( const double c, const FFDep&S )
{
  return S;
}

inline FFDep
operator-
( const FFDep&S, const double c )
{
  return S;
}

inline FFDep
operator-
( const FFDep&S1, const FFDep&S2 )
{
  if( S1._dep.empty() ) return S2;
  if( S2._dep.empty() ) return S1;
  return FFDep::combine( S1, S2, FFDep::TYPE::L );
}

inline FFDep
operator*
( const double c, const FFDep&S )
{
  return S;
}

inline FFDep
operator*
( const FFDep&S, const double c )
{
  return S;
}

inline FFDep&
FFDep::operator*=
( const FFDep&S )
{
  if( S._dep.empty() )
    return *this;

  if( _dep.empty() )
    return combine( S, TYPE::L );

  TYPE w = worst(), depw = S.worst(), wmax = ( w>depw? w: depw );

  if( wmax == TYPE::L )
    return combine( S, TYPE::Q );

  if( w == TYPE::L ){
    auto it = _dep.begin();
    TYPE wthres = ( depw>TYPE::P? depw: TYPE::P );
    for( ; it != _dep.end(); ++it ){
      //if( it->second < TYPE::P ) it->second = TYPE::P;
      if( it->second < wthres ) it->second = wthres;
    }
    auto cit = S._dep.cbegin();
    for( ; cit != S._dep.cend(); ++cit ){
      auto ins = _dep.insert( *cit );
      switch( ins.first->second ){
       case TYPE::L: ins.first->second = TYPE::Q; break;
       case TYPE::Q: ins.first->second = TYPE::P; break;
       case TYPE::P:
       case TYPE::R:
       case TYPE::N: break;
      }
    }
  }
  
  else if( depw == TYPE::L ){
    auto it = _dep.begin();
    for( ; it != _dep.end(); ++it ){
      switch( it->second ){
       case TYPE::L: it->second = TYPE::Q; break;
       case TYPE::Q: it->second = TYPE::P; break;
       case TYPE::P:
       case TYPE::R:
       case TYPE::N: break;
      }
    }
    auto cit = S._dep.cbegin();
    TYPE wthres = ( w>TYPE::P? w: TYPE::P );
    for( ; cit != S._dep.cend(); ++cit ){
      auto ins = _dep.insert( *cit );
      //if( ins.first->second < TYPE::P ) ins.first->second = TYPE::P;
      if( ins.first->second < wthres ) ins.first->second = wthres;
    }
    return *this;
  }

  else{
    auto it = _dep.begin();
    TYPE wthres = ( depw>TYPE::P? depw: TYPE::P );
    for( ; it != _dep.end(); ++it ){
      //if( it->second < TYPE::P ) it->second = TYPE::P;
      if( it->second < wthres ) it->second = wthres;
    }
    auto cit = S._dep.cbegin();
    wthres = ( w>TYPE::P? w: TYPE::P );
    for( ; cit != S._dep.cend(); ++cit ){
      auto ins = _dep.insert( *cit );
      //if( ins.first->second < TYPE::P ) ins.first->second = TYPE::P;
      if( ins.first->second < wthres ) ins.first->second = wthres;
    }
    //return combine( S, wmax>TYPE::P? wmax: TYPE::P ); 
    //return combine( S, TYPE::P ); 
  }

  return *this;
}

inline FFDep
operator*
( const FFDep&S1, const FFDep&S2 )
{
  FFDep S3( S1 );
  return S3 *= S2 ;
}

inline FFDep
sqr
( const FFDep&S )
{
  FFDep S2( S );
  return S2 *= S ;
}

inline FFDep
prod
( const unsigned int n, const FFDep*S )
{
  switch( n ){
   case 0:  return 0.;
   case 1:  return S[0];
   case 2:  return S[0] * S[1];
   default: return S[0] * prod( n-1, S+1 );
  }
}

inline FFDep
monom
( const unsigned int n, const FFDep*S, const unsigned*k )
{
  switch( n ){
   case 0:  return 0.;
   case 1:  return pow( S[0], (int)k[0] );
   default: return pow( S[0], (int)k[0] ) * monom( n-1, S+1, k+1 );
  }
}

inline FFDep
operator/
( const FFDep&S, const double c )
{
  return S;
}

inline FFDep
operator/
( const double c, const FFDep&S )
{
  return inv( S );
}

inline FFDep&
FFDep::operator/=
( const FFDep&S )
{
  if( S._dep.empty() ) return *this;
  return operator*=( inv(S) );
}

inline FFDep
operator/
( const FFDep&S1, const FFDep&S2 )
{
  FFDep S3( S1 );
  return S3 /= S2 ;
}

inline FFDep
inv
( const FFDep &S )
{
  return FFDep::copy( S, FFDep::TYPE::R );
}

inline FFDep
sqrt
( const FFDep&S )
{
  return FFDep::copy( S, FFDep::TYPE::N );
}

inline FFDep
exp
( const FFDep &S )
{
  return FFDep::copy( S, FFDep::TYPE::N );
}

inline FFDep
arh
( const FFDep &S, const double a )
{
  return FFDep::copy( S, FFDep::TYPE::N );
}

inline FFDep
log
( const FFDep &S )
{
  return FFDep::copy( S, FFDep::TYPE::N );
}

inline FFDep
xlog
( const FFDep&S )
{
  return FFDep::copy( S, FFDep::TYPE::N );
}

inline FFDep
lmtd
( const FFDep&S1, const FFDep&S2 )
{
  return FFDep::combine( S1, S2, FFDep::TYPE::N );
}

inline FFDep
rlmtd
( const FFDep&S1, const FFDep&S2 )
{
  return FFDep::combine( S1, S2, FFDep::TYPE::N );
}

inline FFDep
erf
( const FFDep &S )
{
  return FFDep::copy( S, FFDep::TYPE::N );
}

inline FFDep
erfc
( const FFDep &S )
{
  return FFDep::copy( S, FFDep::TYPE::N );
}

inline FFDep
fstep
( const FFDep &S )
{
  return FFDep::copy( S, FFDep::TYPE::N );
}

inline FFDep
bstep
( const FFDep &S )
{
  return FFDep::copy( S, FFDep::TYPE::N );
}

inline FFDep
fabs
( const FFDep&S )
{
  return FFDep::copy( S, FFDep::TYPE::N );
}

inline FFDep
cheb
( const FFDep&S, const unsigned n )
{
  if( n == 0 ) return FFDep();
  if( n == 1 ) return S;
  if( n == 2 ) return sqr(S);
  return FFDep::copy( S, FFDep::TYPE::P );
}

inline FFDep
pow
( const FFDep&S, const int n )
{
  if( n == 0 ) return FFDep();
  if( n == 1 ) return S;
  if( n == 2 ) return sqr(S);
  if( n < 0  ) return FFDep::copy( S, FFDep::TYPE::R );
  return FFDep::copy( S, FFDep::TYPE::P );
}

inline FFDep
pow
( const FFDep&S, const double a )
{
  if( a == 0. ) return FFDep();
  if( a == 1. ) return S;
  if( a == 2. ) return sqr(S);
  if( a > 1. && a == std::rint(a) ) return FFDep::copy( S, FFDep::TYPE::P );
  if( a < 0. && a == std::rint(a) ) return FFDep::copy( S, FFDep::TYPE::R );
  return FFDep::copy( S, FFDep::TYPE::N );
}

inline FFDep
pow
( const FFDep&S1, const FFDep&S2 )
{
  return FFDep::combine( S1, S2, FFDep::TYPE::N );
}

inline FFDep
min
( const FFDep&S1, const FFDep&S2 )
{
  return FFDep::combine( S1, S2, FFDep::TYPE::N );
}

inline FFDep
max
( const FFDep&S1, const FFDep&S2 )
{
  return FFDep::combine( S1, S2, FFDep::TYPE::N );
}

inline FFDep
min
( const unsigned int n, const FFDep*S )
{
  switch( n ){
   case 0:  return FFDep();
   case 1:  return S[0];
   case 2:  return min( S[0], S[1] );
   default: return min( S[0], min( n-1, S+1 ) );
  }
}

inline FFDep
max
( const unsigned int n, const FFDep*S )
{
  switch( n ){
   case 0:  return FFDep();
   case 1:  return S[0];
   case 2:  return max( S[0], S[1] );
   default: return max( S[0], max( n-1, S+1 ) );
  }
}

inline FFDep
cos
( const FFDep&S )
{
  return FFDep::copy( S, FFDep::TYPE::N );
}

inline FFDep
sin
( const FFDep &S )
{
  return FFDep::copy( S, FFDep::TYPE::N );
}

inline FFDep
tan
( const FFDep&S )
{
  return FFDep::copy( S, FFDep::TYPE::N );
}

inline FFDep
acos
( const FFDep &S )
{
  return FFDep::copy( S, FFDep::TYPE::N );
}

inline FFDep
asin
( const FFDep &S )
{
  return FFDep::copy( S, FFDep::TYPE::N );
}

inline FFDep
atan
( const FFDep &S )
{
  return FFDep::copy( S, FFDep::TYPE::N );
}

inline FFDep
cosh
( const FFDep &S )
{
  return FFDep::copy( S, FFDep::TYPE::N );
}

inline FFDep
sinh
( const FFDep &S )
{
  return FFDep::copy( S, FFDep::TYPE::N );
}

inline FFDep
tanh
( const FFDep &S )
{
  return FFDep::copy( S, FFDep::TYPE::N );
}

inline FFDep
inter
( const FFDep&S1, const FFDep&S2 )
{
  return FFDep::combine( S1, S2, FFDep::TYPE::N );
}

inline bool
operator==
( const FFDep&S1, const FFDep&S2 )
{
  return( S1._dep == S2._dep );
}

inline bool
operator!=
( const FFDep&S1, const FFDep&S2 )
{
  return( S1._dep != S2._dep );
}

inline bool
operator<=
( const FFDep&S1, const FFDep&S2 )
{
  return( S1._dep <= S2._dep );
}

inline bool
operator>=
( const FFDep&S1, const FFDep&S2 )
{
  return( S1._dep >= S2._dep );
}

inline bool
operator<
( const FFDep&S1, const FFDep&S2 )
{
  return( S1._dep < S2._dep );
}

inline bool
operator>
( const FFDep&S1, const FFDep&S2 )
{
  return( S1._dep > S2._dep );
}

} // namespace mc

namespace mc
{

//! @brief Specialization of the structure mc::Op to allow usage of the type mc::Interval for DAG evaluation or as a template parameter in other MC++ classes
template <> struct Op< mc::FFDep >
{
  typedef mc::FFDep FV;
  static FV point( const double c ) { return FV(c); }
  static FV zeroone() { throw typename FFDep::Exceptions( FFDep::Exceptions::UNDEF ); }
  static void I(FV& x, const FV&y) { x = y; }
  static double l(const FV& x) { throw typename FFDep::Exceptions( FFDep::Exceptions::UNDEF ); }
  static double u(const FV& x) { throw typename FFDep::Exceptions( FFDep::Exceptions::UNDEF ); }
  static double abs (const FV& x) { throw typename FFDep::Exceptions( FFDep::Exceptions::UNDEF );  }
  static double mid (const FV& x) { throw typename FFDep::Exceptions( FFDep::Exceptions::UNDEF );  }
  static double diam(const FV& x) { throw typename FFDep::Exceptions( FFDep::Exceptions::UNDEF ); }
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
  static FV hull(const FV& x, const FV& y) { throw typename FFDep::Exceptions( FFDep::Exceptions::UNDEF ); }
  static FV min (const FV& x, const FV& y) { return mc::min(x,y);  }
  static FV max (const FV& x, const FV& y) { return mc::max(x,y);  }
  static FV arh (const FV& x, const double k) { return mc::exp(-k/x); }
  template <typename X, typename Y> static FV pow(const X& x, const Y& y) { return mc::pow(x,y); }
  static FV cheb(const FV& x, const unsigned n) { return mc::cheb(x,n); }
  static FV prod(const unsigned int n, const FV* x) { return mc::prod(n,x); }
  static FV monom(const unsigned int n, const FV* x, const unsigned* k) { return mc::monom(n,x,k); }
  static bool inter(FV& xIy, const FV& x, const FV& y) { xIy = mc::inter(x,y); return true; }
  static bool eq(const FV& x, const FV& y) { throw typename FFDep::Exceptions( FFDep::Exceptions::UNDEF ); }
  static bool ne(const FV& x, const FV& y) { throw typename FFDep::Exceptions( FFDep::Exceptions::UNDEF ); }
  static bool lt(const FV& x, const FV& y) { throw typename FFDep::Exceptions( FFDep::Exceptions::UNDEF ); }
  static bool le(const FV& x, const FV& y) { throw typename FFDep::Exceptions( FFDep::Exceptions::UNDEF ); }
  static bool gt(const FV& x, const FV& y) { throw typename FFDep::Exceptions( FFDep::Exceptions::UNDEF ); }
  static bool ge(const FV& x, const FV& y) { throw typename FFDep::Exceptions( FFDep::Exceptions::UNDEF ); }
};

} // namespace mc

#endif
