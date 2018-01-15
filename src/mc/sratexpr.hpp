// Copyright (C) 2017 Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the Eclipse Public License.

#ifndef MC__SRATEXPR_H
#define MC__SRATEXPR_H

#include "spolyexpr.hpp"

#undef  MC__SPOLYMODEL_DEBUG
#define MC__SPOLYMODEL_DEBUG_POLYBOUND
#define MC__SPOLYMODEL_CHECK
#undef  MC__SPOLYEXPR_DEBUG_SPROD

namespace mc
{
//! @brief C++ class for reformulation of a factorable function in sparse rational form
////////////////////////////////////////////////////////////////////////
//! mc::SRatEnv is a C++ class for reformulation of a factorable
//! function in sparse rational form
////////////////////////////////////////////////////////////////////////
class SRatEnv
////////////////////////////////////////////////////////////////////////
{
};

//! @brief C++ class for sparse rational function representation and arithmetic
////////////////////////////////////////////////////////////////////////
//! mc::SRatExpr is a C++ class for sparse rational function
//! representation and arithmetic
////////////////////////////////////////////////////////////////////////
class SRatExpr
////////////////////////////////////////////////////////////////////////
{
private:
  //! @brief Pointer to sparse rational function environment
  SRatEnv *_env;

protected:
  // Denominator sparse polynomial
  SPolyExpr _denom;

  // numerator sparse polynomial
  SPolyExpr _numer;

  //! @brief Set rational expression equal to <a>srat</a>
  SRatExpr& _set
    ( const SRatExpr& srat );

  //! @brief Set rational expression equal to constant <a>d</a>
  SRatExpr& _set
    ( const double d );

  //! @brief Set rational expression equal to variable <a>x</a>
  SRatExpr& _set
    ( const FFVar& x );

public:
  /** @ingroup FFunc
   *  @{
   */
  //! @brief Default constructor of sparse rational expression
  SRatExpr
    ()
    {}

  //! @brief Constructor of sparse rational expression as constant
  SRatExpr
    ( const double d )
    { _set( d ); }

  //! @brief Constructor of sparse rational expression as variable
  SRatExpr
    ( const FFVar& x)
    { _set( x ); }

  //! @brief Constructor of sparse rational expression as variable
  SRatExpr
    ( const SPolyExpr& n, const SPolyExpr& d )
    { _numer = n; _denom = d; }

  //! @brief Copy constructor of sparse rational expression
  SRatExpr
    ( const SRatExpr& srat )
    { _set( srat ); }

  //! @brief Destructor of sparse rational expression
  virtual ~SRatExpr()
    {}

  //! @brief Overloaded operator '=' for sparse rational expression
  SRatExpr& operator=
    ( const SRatExpr& srat )
    { _set( srat ); return *this; }

  //! @brief Overloaded operator '=' for constant 
  SRatExpr& operator=
    ( const double d )
    { _set( d ); return *this; }

  //! @brief Overloaded operator '=' for variable
  SRatExpr& operator=
    ( const FFVar& x )
    { _set( x ); return *this; }

  //! @brief Overloaded operator '+=' for sparse rational function
  SRatExpr& operator+=
    ( const SRatExpr& srat );

  //! @brief Overloaded operator '-=' for sparse rational function
  SRatExpr& operator-=
    ( const SRatExpr& srat );

  //! @brief Overloaded operator '*=' for sparse rational function
  SRatExpr& operator*=
    ( const SRatExpr& srat );

  //! @brief Overloaded operator '/=' for sparse rational function
  SRatExpr& operator/=
    ( const SRatExpr& srat );

  // Denominator sparse polynomial
  const SPolyExpr& denom
    ()
    const
    { return _denom; };

  // numerator sparse polynomial
  const SPolyExpr& numer
    ()
    const
    { return _numer; };

  //! @brief Exceptions of mc::SRatExpr
  class Exceptions
  {
   public:
    //! @brief Enumeration type for SRatExpr exception handling
    enum TYPE{
      ENV=0,          //!< Operation between sparse rational expressions linked to different environments
      INTERNAL = -33  //!< Internal error
    };
    //! @brief Constructor for error <a>ierr</a>
    Exceptions( TYPE ierr ) : _ierr( ierr ){}
    //! @brief Error flag
    int ierr(){ return _ierr; }
    //! @brief Error description
    std::string what(){
      switch( _ierr ){
      case ENV:
        return "mc::SRatExpr\t Operation between sparse rational expressions linked to different environments is not allowed";
      case INTERNAL:
      default:
        return "mc::SRatExpr\t Internal error";
      }
    }
   private:
    TYPE _ierr;
  };

//  //! @brief Options of mc::SPolyExpr
//  static struct Options
//  {
//    //! @brief Constructor
//    Options():
//      BASIS(MONOM), REMZERO(true), DISPLEN(5)
//      {}
//    //! @brief Assignment of mc::SCModel::Options
//    Options& operator=
//      ( Options& opt ){
//        BASIS   = opt.BASIS;
//        DISPLEN = opt.DISPLEN;
//        return *this;
//      }
//    //! @brief Available basis representations
//    enum BASIS_TYPE{
//      MONOM=0,	//!< Monomial basis
//      CHEB	    //!< Chebyshev basis
//    };
//    //! @brief Basis representation of sparse polynomial
//    BASIS_TYPE BASIS;
//    //! @brief Whether to remove zeros entries from sparse polynomials
//    bool REMZERO;
//    //! @brief Number of digits in output stream for sparse polynomial coefficients
//    unsigned DISPLEN;
//  } options;
  /** @} */
};

//SPolyExpr::Options SPolyExpr::options;

///////////////////////////////// SRatExpr /////////////////////////////////////

inline SRatExpr&
SRatExpr::_set
( const double d )
{
  _numer = d;
  _denom = 1;
  return *this;
}

inline SRatExpr&
SRatExpr::_set
( const FFVar&x )
{
  _numer = x;
  _denom = 1;
  return *this;
}

inline SRatExpr&
SRatExpr::_set
( const SRatExpr&srat )
{
  if( this == &srat ) return *this;
  _numer = srat._numer;
  _denom = srat._denom;
  return *this;
}

inline std::ostream&
operator<<
( std::ostream&out, const SRatExpr&srat )
{
  out << "NUMERATOR:"   << srat.numer() << std::endl;
  out << "DENOMINATOR:" << srat.denom();
  return out;
}

inline SRatExpr
operator+
( const SRatExpr&srat )
{
  return srat;
}

inline SRatExpr&
SRatExpr::operator+=
( const SRatExpr&srat )
{
  _numer *= srat._denom;
  _numer += srat._numer * _denom;
  _denom *= srat._denom;
  return *this;
}

inline SRatExpr
operator+
( const SRatExpr&srat1, const SRatExpr&srat2 )
{
  SRatExpr srat3( srat1 );
  srat3 += srat2;
  return srat3;
}

inline SRatExpr
operator-
( const SRatExpr&srat )
{
  return SRatExpr( -srat.numer(), srat.denom() );
}

inline SRatExpr&
SRatExpr::operator-=
( const SRatExpr&srat )
{
  _numer *= srat._denom;
  _numer -= srat._numer * _denom;
  _denom *= srat._denom;
  return *this;
}

inline SRatExpr
operator-
( const SRatExpr&srat1, const SRatExpr&srat2 )
{
  SRatExpr srat3( srat1 );
  srat3 -= srat2;
  return srat3;
}

inline SRatExpr&
SRatExpr::operator*=
( const SRatExpr&srat )
{
  _numer *= srat._numer;
  _denom *= srat._denom;
  return *this;
}

inline SRatExpr
operator*
( const SRatExpr&srat1, const SRatExpr&srat2 )
{
  SRatExpr srat3( srat1 );
  srat3 *= srat2;
  return srat3;
}

inline SRatExpr
inv
( const SRatExpr&srat )
{
  return SRatExpr( srat.denom(), srat.numer() );
}

inline SRatExpr
sqr
( const SRatExpr&srat )
{
  return SRatExpr( sqr( srat.numer() ), sqr( srat.denom() ) );
}

inline SRatExpr
pow
( const SRatExpr&srat, const int n )
{
  if( n < 0 ) return pow( inv( srat ), -n );
  switch( n ){
   case 0:  return 1.;
   case 1:  return srat;
   case 2:  return sqr( srat );
   default: return SRatExpr( pow( srat.numer(), (unsigned)n ), pow( srat.denom(), (unsigned)n ) );
  }
}

inline SRatExpr&
SRatExpr::operator/=
( const SRatExpr&srat )
{
  _numer *= srat._denom;
  _denom *= srat._numer;
  return *this;
}

inline SRatExpr
operator/
( const SRatExpr&srat1, const SRatExpr&srat2 )
{
  SRatExpr srat3( srat1 );
  srat3 /= srat2;
  return srat3;
}

} // namespace mc

#endif

