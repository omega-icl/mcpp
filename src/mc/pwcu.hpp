// Copyright (C) Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the Eclipse Public License.

/*!
\page page_PWCU Piecewise constant univariate estimators on fixed partition
\author Yanlin Zha, Beno&icirc;t Chachuat

The class mc::PWCU provides an implementation of piecewise constant univariate estimators for use within a superposition relaxation (see \ref page_SUPREL). It support the definition of such estimators, their propagation through linear operations and compositions with convex/concave monotonic terms, and their bounding and evaluation.
*/

#ifndef MC__PWCU_HPP
#define MC__PWCU_HPP

#include <iostream>
#include <iomanip> 
#include <vector> 
#include <algorithm>

#include "mcop.hpp"
#include "mcfunc.hpp"

#define MC__PWCU_DEBUG
#define MC__PWCU_TRACE
#define MC__PWCU_CHECK
#undef  MC__PWCU_FULL_UPDATE

namespace mc
{
//! @brief C++ class for propagation of piecewise constant univariate estimators
////////////////////////////////////////////////////////////////////////
//! mc::PWCU is a C++ class for propagation of piecewise constant 
//! univariate estimators through factorable expressions.
////////////////////////////////////////////////////////////////////////
class PWCU
////////////////////////////////////////////////////////////////////////
{
 private:

  //! @brief initial abscissa
  double _xL;

  //! @brief final abscissa
  double _xU;

  //! @brief vector of y lower bounds
  std::vector<double> _yL;

  //! @brief vector of y upper bounds
  std::vector<double> _yU;

 public:

  //! @brief Options of mc::PWCU
  static struct Options
  {
    //! @brief Constructor
    Options():
      BKPTATOL( 1e2*DBL_EPSILON ),
      BKPTRTOL( 1e2*DBL_EPSILON ),
      DISPNUM( 5 )
      {}
    //! @brief Assignment
    Options& operator=
      ( Options const& opt )
      {
        BKPTATOL        = opt.BKPTATOL;
        BKPTRTOL        = opt.BKPTRTOL;
        DISPNUM         = opt.DISPNUM;
        return *this;
      }
    //! @brief Assignment
    void reset
      ()
      {
        BKPTATOL = 1e2*DBL_EPSILON;
        BKPTRTOL = 1e2*DBL_EPSILON;
        DISPNUM  = 5;
      }
    //! @brief Absolute tolerance in breakpoints - Default: 1e2*DBL_EPSILON
    double BKPTATOL;
    //! @brief Relative tolerance in breakpoints - Default: 1e2*DBL_EPSILON
    double BKPTRTOL;
    //! @brief Number of numerical digits displayed with << operator - Default: 5
    unsigned int DISPNUM;
  } options;

  //! @brief Exceptions of mc::PWCU
  class Exceptions
  {
   public:
    //! @brief Enumeration type for SElimVar exception handling
    enum TYPE{
      RANGE=0,        //!< Operation on variable with empty range
      SIZE,           //!< Inconsistent array or list sizes in estimator
      EXTRAPOL,       //!< Extrapolation outside of variable range
      DIV,	      //!< Division by zero
      INTERNAL=-1,    //!< Internal error
      UNDEF=-33       //!< Feature not yet implemented
    };
    //! @brief Constructor for error <a>ierr</a>
    Exceptions( TYPE ierr ) : _ierr( ierr ){}
    //! @brief Error flag
    int ierr(){ return _ierr; }
    //! @brief Error description
    std::string what(){
      switch( _ierr ){
      case RANGE:
        return "mc::PWCU\t Operation on variable with empty range";
      case SIZE:
        return "mc::PWCU\t Inconsistent array or list sizes in estimator";
      case EXTRAPOL:
        return "mc::PWCU\t Extrapolation outside of variable range";
      case DIV:
        return "mc::PWCU\t Division by zero";
      case UNDEF:
        return "mc::PWCU\t Feature not yet implemented";
      case INTERNAL:
      default:
        return "mc::PWCU\t Internal error";
      }
    }
   private:
    TYPE _ierr;
  };

  //! @brief Default constructor
  PWCU
    ()
    : _yL(), _yU()
    {}

  //! @brief Constructor of constant
  PWCU
    ( double const& xL, double const& xU, double const& y, size_t const& N )
    : _xL( xL ), _xU( xU ), _yL( N>1?N:1, y-DBL_EPSILON ), _yU( N>1?N:1, y+DBL_EPSILON )
    {}

  //! @brief Constructor of variable
  PWCU
    ( double xL, double xU, size_t const& N )
    : _xL( xL ), _xU( xU ), _yL( N>1? N: 1 ), _yU( N>1? N: 1 )
    {
      xU = (xU-xL)/_yL.size();
#ifdef MC__PWCU_CHECK
      if( xU <= 0. )
        throw Exceptions( Exceptions::RANGE );
#endif
      for( auto iyL=_yL.begin(), iyU=_yU.begin(); iyL!=_yL.end(); ++iyL, ++iyU ){
        *iyL = xL;
        xL += xU;
        *iyU = xL;
      }
    }

  //! @brief Copy constructor
  PWCU
    ( PWCU const& var )
    : _xL( var._xL ), _xU( var._xU ), _yL( var._yL ), _yU( var._yU )
    {}

  //! @brief Move constructor
  PWCU
    ( PWCU && var )
    : _xL( std::move(var._xL) ), _xU( std::move(var._xU) ),
      _yL( std::move(var._yL) ), _yU( std::move(var._yU) )
    {}

  //! @brief Set constant estimator
  PWCU& set
    ( double const& y )
    {
      _yL.assign( _yL.size(), y-DBL_EPSILON );
      _yU.assign( _yU.size(), y+DBL_EPSILON );
      return *this;
    }

  //! @brief Set constant estimator
  PWCU& set
    ( double const& xL, double const& xU, double const& y, size_t const& N )
    {
      _xL = xL;
      _xU = xU;
      _yL.assign( N>1?N:1, y-DBL_EPSILON );
      _yU.assign( N>1?N:1, y+DBL_EPSILON );

      return *this;
    }

  //! @brief Set variable estimator
  PWCU& set
    ( double xL, double xU, size_t const& N )
    {
      _xL = xL;
      _xU = xU;

      _yL.resize( N>1? N: 1 );
      _yU.resize( N>1? N: 1 );
      xU = (xU-xL)/_yL.size();
#ifdef MC__PWCU_CHECK
      if( xU <= 0. )
        throw Exceptions( Exceptions::RANGE );
#endif
      for( auto iyL=_yL.begin(), iyU=_yU.begin(); iyL!=_yL.end(); ++iyL, ++iyU ){
        *iyL = xL;
        xL += xU;
        *iyU = xL;
      }

      return *this;
    }

  //! @brief Evaluate lower bound at point
  double l
    ( double const& x )
    const
    {
#ifdef MC__PWCU_CHECK
      if( _yL.empty() )
        throw Exceptions( Exceptions::SIZE );
#endif
      double const dx = (_xU-_xL)/_yL.size();
      double xi = _xL;
      for( auto iyL=_yL.cbegin(); iyL!=_yL.cend(); ++iyL ){
        if( xi + dx < x ){
          xi += dx;
          continue;
        }
        return *iyL;
      }
      if( isequal( xi, x, options.BKPTATOL, options.BKPTRTOL ) )
        return *_yL.crbegin();
      throw Exceptions( Exceptions::EXTRAPOL );
    }

  //! @brief Evaluate upper bound at point
  double u
    ( double const& x )
    const
    {
#ifdef MC__PWCU_CHECK
      if( _yU.empty() )
        throw Exceptions( Exceptions::SIZE );
#endif
      double const dx = (_xU-_xL)/_yU.size();
      double xi = _xL;
      for( auto iyU=_yU.cbegin(); iyU!=_yU.cend(); ++iyU ){
        if( xi + dx < x ){
          xi += dx;
          continue;
        }
        return *iyU;
      }
      if( isequal( xi, x, options.BKPTATOL, options.BKPTRTOL ) )
        return *_yU.crbegin();
      throw Exceptions( Exceptions::EXTRAPOL );
    }

  //! @brief Evaluate lower range
  double l
    ()
    const
    {
#ifdef MC__PWCU_CHECK
      if( _yL.empty() )
        throw Exceptions( Exceptions::SIZE );
#endif
      return *std::min_element( _yL.cbegin(), _yL.cend() );
    }

  //! @brief Evaluate upper range
  double u
    ()
    const
    {
#ifdef MC__PWCU_CHECK
      if( _yU.empty() )
        throw Exceptions( Exceptions::SIZE );
#endif
      return *std::max_element( _yU.cbegin(), _yU.cend() );
    }

  //! @brief Evaluate range
  double w
    ()
    const
    {
      return u() - l();
    }

  //! @brief Display estimator
  std::ostream& display
    ( std::ostream& os=std::cout, size_t const dispnum=options.DISPNUM )
    const
    {
#ifdef MC__PWCU_CHECK
      if( _yL.empty() || _yL.size() != _yU.size() )
        throw Exceptions( Exceptions::SIZE );
#endif
      double const dx = (_xU-_xL)/_yL.size();
      double xi = _xL;
      os << "{" << std::scientific << std::setprecision(dispnum) << std::right
         << std::setw(dispnum+7) << xi;// << std::setw(dispnum+8) << yi;
      for( auto iyL=_yL.cbegin(), iyU=_yU.cbegin(); iyL!=_yL.cend(); ++iyL, ++iyU )
        os << " <" << std::setw(dispnum+7) << *iyL << " : " << std::setw(dispnum+7) << *iyU <<  "> " 
           << std::setw(dispnum+7) << (xi += dx);
      return os << " }";
    }

  //! @brief Retreive number of segments
  size_t size
    ()
    const
    {
#ifdef MC__PWLU_CHECK
      if( _yL.size() != _yU.size() )
        throw Exceptions( Exceptions::SIZE );
#endif
      return _yL.size();
    }

  //! @brief Retreive lower levels
  std::vector<double> const& yL
    ()
    const
    { return _yL; }

  //! @brief Retreive upper levels
  std::vector<double> const& yU
    ()
    const
    { return _yU; }

  //! @brief Retreive initial abscissa
  double xL
    ()
    const
    { return _xL; }

  //! @brief Retreive final abscissa
  double xU
    ()
    const
    { return _xU; }

  //! @brief Retreive/set lower levels
  std::vector<double>& yL
    ()
    { return _yL; }

  //! @brief Retreive/set upper levels
  std::vector<double>& yU
    ()
    { return _yU; }

  //! @brief Retreive/set initial abscissa
  double& xL
    ()
    { return _xL; }

  //! @brief Retreive/set final abscissa
  double& xU
    ()
    { return _xU; }

  PWCU& operator=
    ( PWCU const& var )
    {
      _xL = var._xL;
      _xU = var._xU;
      _yL = var._yL;
      _yU = var._yU;
      return *this;
    }

  PWCU& operator=
    ( PWCU && var )
    {
      _xL = std::move( var._xL );
      _xU = std::move( var._xU );
      _yL = std::move( var._yL );
      _yU = std::move( var._yU );
      return *this;
    }

  PWCU& operator+=
    ( double const& cst )
    {
      for( auto iyL=_yL.begin(), iyU=_yU.begin(); iyL!=_yL.end(); ++iyL, ++iyU ){
        *iyL += cst;
        *iyU += cst;
      }
      return *this;
    }

  PWCU& operator-=
    ( double const& cst )
    {
      for( auto iyL=_yL.begin(), iyU=_yU.begin(); iyL!=_yL.end(); ++iyL, ++iyU ){
        *iyL -= cst;
        *iyU -= cst;
      }
      return *this;
    }
 
  PWCU& operator*=
    ( double const& cst )
    {
      if( cst == 1. )
        return *this;
      for( auto iyL=_yL.begin(), iyU=_yU.begin(); iyL!=_yL.end(); ++iyL, ++iyU ){
        *iyL *= cst;
        *iyU *= cst;
      }
      if( cst < 0 )
        std::swap( _yL, _yU );
      return *this;
    }

  PWCU& operator/=
    ( double const& cst )
    {
      if( cst == 1. )
        return *this;
      for( auto iyL=_yL.begin(), iyU=_yU.begin(); iyL!=_yL.end(); ++iyL, ++iyU ){
        *iyL /= cst;
        *iyU /= cst;
      }
      if( cst < 0 )
        std::swap( _yL, _yU );
      return *this;
    }
           
  PWCU& operator+=
    ( PWCU const& var )
    {
#ifdef MC__PWCU_CHECK
      if( _yL.empty() || _yL.size() != var._yL.size() 
       || _yU.empty() || _yU.size() != var._yU.size() )
        throw Exceptions( Exceptions::SIZE );
#endif
      auto iyL1 = _yL.begin(),      iyU1 = _yU.begin();
      auto iyL2 = var._yL.cbegin(), iyU2 = var._yU.cbegin();
      for( ; iyL1 != _yL.end(); ++iyL1, ++iyU1, ++iyL2, ++iyU2 ){
        *iyL1 += *iyL2;
        *iyU1 += *iyU2;
      }
      return *this;
    }

  PWCU& operator-=
    ( PWCU const& var )
    {
#ifdef MC__PWCU_CHECK
      if( _yL.empty() || _yL.size() != var._yL.size() 
       || _yU.empty() || _yU.size() != var._yU.size() )
        throw Exceptions( Exceptions::SIZE );
#endif
      auto iyL1 = _yL.begin(),      iyU1 = _yU.begin();
      auto iyL2 = var._yL.cbegin(), iyU2 = var._yU.cbegin();
      for( ; iyL1 != _yL.end(); ++iyL1, ++iyU1, ++iyL2, ++iyU2 ){
        *iyL1 -= *iyU2;
        *iyU1 -= *iyL2;
      }
      return *this;
    }

  template <typename UNIV, typename DUNIV>
  PWCU& compose
    ( UNIV const& f, DUNIV const& df, bool const under, bool const cvx, bool const inc )
    {
      if( _yL.empty() || _yL.size() != _yU.size() )
        return *this;

      for( auto iyL=_yL.begin(), iyU=_yU.begin(); iyL!=_yL.end(); ++iyL, ++iyU ){
        *iyL = f( *iyL );
        *iyU = f( *iyU );
      }

      if( !inc )
        std::swap( _yL, _yU );

      return *this;
    }

  PWCU& min
    ( double const& c )
    {
      if( _yL.empty() || _yL.size() != _yU.size() )
        return *this;

      for( auto iyL=_yL.begin(), iyU=_yU.begin(); iyL!=_yL.end(); ++iyL, ++iyU ){
        if( *iyL > c ) *iyL = c - DBL_EPSILON; 
        if( *iyU > c ) *iyU = c; 
      }

      return *this;
    }

  PWCU& max
    ( double const& c )
    {
      if( _yL.empty() || _yL.size() != _yU.size() )
        return *this;

      for( auto iyL=_yL.begin(), iyU=_yU.begin(); iyL!=_yL.end(); ++iyL, ++iyU ){
        if( *iyL < c ) *iyL = c; 
        if( *iyU < c ) *iyU = c + DBL_EPSILON; 
      }

      return *this;
    }

  PWCU& reduce
    ( bool const under, size_t const nseg )
    {
      return *this;
    }

  PWCU& clean
    ( bool const under )
    {
      return *this;
    }
};

////////////////////////////////////////////////////////////////////////

inline PWCU::Options PWCU::options;

std::ostream& operator<<
( std::ostream& os, PWCU const& var )
{
  var.display( os );
  return os;
}

inline
PWCU operator+
( PWCU const& var )
{
  return var;
}

inline
PWCU operator+
( PWCU && var )
{
  PWCU res( std::move(var) );  
  return var;
}

inline
PWCU operator+
( PWCU const& var1, PWCU const& var2 )
{
  PWCU res( var1 );
  res += var2;
  return res;
}
/*
inline
PWCU operator+
( PWCU const& var1, PWCU && var2 )
{
  PWCU res( std::move(var2) );
  res += var1;
  return res;
}
*/
inline
PWCU operator+
( PWCU && var1, PWCU const& var2 )
{
  PWCU res( std::move(var1) );
  res += var2;
  return res;
}

inline
PWCU operator+
( PWCU const& var1, double const& cst2 )
{
  PWCU res( var1 );
  res += cst2;
  return res;
}

inline
PWCU operator+
( PWCU && var1, double const& cst2 )
{
  PWCU res( std::move(var1) );
  res += cst2;
  return res;
}

inline
PWCU operator+
( double const& cst1, PWCU const& var2 )
{
  PWCU res( var2 );
  res += cst1;
  return res;
}

inline
PWCU operator+
( double const& cst1, PWCU && var2 )
{
  PWCU res( std::move(var2) );
  res += cst1;
  return res;
}

inline
PWCU operator-
( PWCU const& var )
{
  PWCU res( var );
  return( res *= -1 );
}

inline
PWCU operator-
( PWCU && var )
{
  PWCU res( std::move(var) );
  return( res *= -1 );
}

inline
PWCU operator-
( PWCU const& var1, PWCU const& var2 )
{
  PWCU res( var1 );
  res -= var2;
  return res;
}
/*
inline
PWCU operator-
( PWCU const& var1, PWCU && var2 )
{
  PWCU res( std::move(-var2) );
  res += var1;  
  return res;
}
*/
inline
PWCU operator-
( PWCU && var1, PWCU const& var2 )
{
  PWCU res( std::move(var1) );
  res -= var2;
  return res;
}

inline
PWCU operator-
( PWCU const& var1, double const& cst2 )
{
  PWCU res( var1 );
  res += -cst2;
  return res;
}

inline
PWCU operator-
( PWCU && var1, double const& cst2 )
{
  PWCU res( std::move(var1) );
  res += -cst2;
  return res;
}

inline
PWCU operator-
( double const& cst1, PWCU const& var2 )
{
  PWCU res( -var2 );
  res += cst1;
  return res;
}

inline
PWCU operator-
( double const& cst1, PWCU && var2 )
{
  PWCU res( -std::move(var2) );
  res += cst1;
  return res;
}

inline
PWCU operator*
( PWCU const& var1, double const& cst2 )
{
  PWCU res( var1 );
  res *= cst2;
  return res;
}

inline
PWCU operator*
( PWCU && var1, double const& cst2 )
{
  PWCU res( std::move(var1) );
  res *= cst2;
  return res;
}

inline
PWCU operator*
( double const& cst1, PWCU const& var2 )
{
  PWCU res( var2 );
  res *= cst1;
  return res;
}

inline
PWCU operator*
( double const& cst1, PWCU && var2 )
{
  PWCU res( std::move(var2) );
  res *= cst1;
  return res;
}

inline
PWCU operator/
( PWCU const& var1, double const& cst2 )
{
  PWCU var3( var1 );
  var3 /= cst2;
  return var3;
}

inline
PWCU operator/
( PWCU && var1, double const& cst2 )
{
  PWCU var3( std::move(var1) );
  var3 /= cst2;
  return var3;
}

} // namespace mc

#endif
