// Copyright (C) Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the Eclipse Public License.

/*!
\page page_PWLU Piecewise linear univariate estimators on adaptive partition
\author Yanlin Zha, Beno&icirc;t Chachuat

The class mc::PWLU provides an implementation of piecewise linear univariate estimators for use within a superposition relaxation (see \ref page_SUPREL). It support the definition of such estimators, their propagation through linear operations and compositions with convex/concave monotonic terms, their simplification or reduction, and their bounding and evaluation.
*/

#ifndef MC__PWLU_HPP
#define MC__PWLU_HPP

#include <fstream>
#include <iostream>
#include <iomanip> 
#include <vector> 
#include <algorithm>
#include <cmath>
#include <limits>
#include <bitset>
#include <cassert>

#include "mcfunc.hpp"

//#define MC__PWLU_DEBUG
//#define MC__PWLU_TRACE
//#define MC__PWLU_CHECK
namespace mc
{
//! @brief C++ class for propagation of piecewise linear univariate estimators
////////////////////////////////////////////////////////////////////////
//! mc::PWLU is a C++ class for propagation of piecewise linear 
//! univariate estimators through factorable expressions.
////////////////////////////////////////////////////////////////////////
class PWLU
////////////////////////////////////////////////////////////////////////
{
private:

  //! @brief initial abscissa
  double _x0;

  //! @brief initial ordinate
  double _y0;

  //! @brief list of distances between breakpoints
  std::vector<double> _dx;
  
  //! @brief list of slopes between breakpoints
  std::vector<double> _dy;

public:

  //! @brief Options of mc::PWLU
  static struct Options
  {
    //! @brief Constructor
    Options()
      {
        reset();
      }
    //! @brief Assignment
    Options& operator=
      ( Options const& opt )
      {
        BKPTATOL        = opt.BKPTATOL;
        BKPTRTOL        = opt.BKPTRTOL;
        REDUCEMETH      = opt.REDUCEMETH;
        DISPNUM         = opt.DISPNUM;
        return *this;
      }
    //! @brief Reset to default options
    void reset
      ()
      {
        BKPTATOL   = 1e3*DBL_EPSILON;
        BKPTRTOL   = 1e3*DBL_EPSILON;
        REDUCEMETH = 0;
        DISPNUM    = 5;
      }
    //! @brief Absolute tolerance in breakpoints - Default: 1e3*DBL_EPSILON
    double BKPTATOL;
    //! @brief Relative tolerance in breakpoints - Default: 1e3*DBL_EPSILON
    double BKPTRTOL;
    //! @brief Method for breakpoint reduction - <=0: Standard library (default); >0: Tailored heap implementation over threshold
    int REDUCEMETH;
    //! @brief Number of numerical digits displayed with << operator - Default: 5
    unsigned int DISPNUM;
  } options;

  //! @brief Exceptions of mc::PWLU
  class Exceptions
  {
   public:
    //! @brief Enumeration type for PWLU exception handling
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
        return "mc::PWLU\t Operation on variable with empty range";
      case SIZE:
        return "mc::PWLU\t Inconsistent array or list sizes in estimator";
      case EXTRAPOL:
        return "mc::PWLU\t Extrapolation outside of variable range";
      case DIV:
        return "mc::PWLU\t Division by zero";
      case UNDEF:
        return "mc::PWLU\t Feature not yet implemented";
      case INTERNAL:
      default:
        return "mc::PWLU\t Internal error";
      }
    }
   private:
    TYPE _ierr;
  };

  //! @brief Default constructor
  PWLU
    ()
    : _x0(0.), _y0(0.), _dx(), _dy()
    {}

  //! @brief Constructor of constant
  PWLU
    ( double const& xL, double const& xU, double const& y )
    : _x0(xL), _y0(y), _dx({xU-xL}), _dy({0.})
    {
#ifdef MC__PWLU_CHECK
      if( xU <= xL )
        throw Exceptions( Exceptions::RANGE );
#endif
    }

  //! @brief Constructor of variable
  PWLU
    ( double const& xL, double const& xU, size_t const N=1 )
    : _x0(xL), _y0(xL), _dx(N>1?N:1, N>1?(xU-xL)/N:(xU-xL)), _dy( N>1?N:1, 1. )
    {
#ifdef MC__PWLU_CHECK
      if( xU < xL )
        throw Exceptions( Exceptions::RANGE );
#endif
    }

  //! @brief Constructor of variable
  PWLU
    ( double const& xL, std::vector<double> const& dx )
    : _x0(xL), _y0(xL), _dx(dx), _dy(dx.size(),1)
    {
#ifdef MC__PWLU_CHECK
      for( auto const& dxi : _dx )
        if( dxi <= 0. ) throw Exceptions( Exceptions::RANGE );
#endif
    }

  //! @brief Constructor of variable
  PWLU
    ( double const& xL, double const& xU, double const& yL, double const& yU )
    : _x0(xL), _y0(yL), _dx({xU-xL}), _dy({(yU-yL)/(xU-xL)})
    {
#ifdef MC__PWLU_CHECK
      if( xU < xL )
        throw Exceptions( Exceptions::RANGE );
#endif
    }

  //! @brief Constructor of variable
  PWLU
    ( double const& xL, double const& yL, std::vector<double> const& dx, std::vector<double> const& dy )
    : _x0(xL), _y0(yL), _dx(dx), _dy(dy)
    {
#ifdef MC__PWLU_CHECK
      if( _dx.size() != _dy.size() )
        throw Exceptions( Exceptions::SIZE );
      for( auto const& dxi : _dx )
        if( dxi <= 0. ) throw Exceptions( Exceptions::RANGE );
#endif
    }

  //! @brief Copy constructor
  PWLU
    ( PWLU const& var )
    : _x0(var._x0), _y0(var._y0), _dx(var._dx), _dy(var._dy)
    {}

  //! @brief Move constructor
  PWLU
    ( PWLU && var )
    : _x0(std::move(var._x0)), _y0(std::move(var._y0)), _dx(std::move(var._dx)), _dy(std::move(var._dy))
    {}

  //! @brief Set constant estimator
  PWLU& set
    ( double const& c )
    {
      _dx = { xU() - _x0 };
      _y0 = c;
      _dy = { 0. };
      return *this;
    }

  //! @brief Set constant estimator
  PWLU& set
    ( double const& xL, double const& xU, double const& y )
    {
#ifdef MC__PWLU_CHECK
      if( xU < xL )
        throw Exceptions( Exceptions::RANGE );
#endif
      _x0 = xL;
      _y0 = y;
      _dx = { xU-xL };
      _dy = { 0. };
      return *this;
    }

  //! @brief Set variable estimator
  PWLU& set
    ( double const& xL, double const& xU, size_t const N=1 )
    {
#ifdef MC__PWLU_CHECK
      if( xU < xL )
        throw Exceptions( Exceptions::RANGE );
#endif
      _x0 = xL;
      _y0 = xL;
      _dx = std::vector<double>( N>1?N:1, N>1?(xU-xL)/N:(xU-xL) );
      _dy = std::vector<double>( N>1?N:1, 1. );
      return *this;
    }

  //! @brief Set variable estimator
  PWLU& set
    ( double const& xL, std::vector<double> const& dx )
    {
#ifdef MC__PWLU_CHECK
      for( auto const& dxi : dx )
        if( dxi <= 0. ) throw Exceptions( Exceptions::RANGE );
#endif
      _x0 = xL;
      _y0 = xL;
      _dx = dx;
      _dy = std::vector<double>( _dx.size(), 1. );
      return *this;
    }

  //! @brief Set variable estimator
  PWLU& set
    ( double const& xL, double const& xU, double const& yL, double const& yU )
    {
#ifdef MC__PWLU_CHECK
      if( xU < xL )
        throw Exceptions( Exceptions::RANGE );
#endif
      _x0 = xL;
      _y0 = yL;
      _dx = { xU-xL };
      _dy = { (yU-yL)/(xU-xL) };
      return *this;
    }

  //! @brief Set variable estimator
  PWLU& set
    ( double const& xL, double const& yL, std::vector<double> const& dx, std::vector<double> const& dy )
    {
#ifdef MC__PWLU_CHECK
      if( dx.size() != dy.size() )
        throw Exceptions( Exceptions::SIZE );
      for( auto const& dxi : _dx )
        if( dxi <= 0. ) throw Exceptions( Exceptions::RANGE );
#endif
      _x0 = xL;
      _y0 = yL;
      _dx = dx;
      _dy = dy;
      return *this;
    }

  //! @brief Insert breakpoints
  PWLU& insert
    ( std::vector<double> const& x )
    {
      _dx.reserve( _dx.size() + x.size() );
      _dy.reserve( _dy.size() + x.size() );
      for( auto const& xi : x )
        insert( xi );
      return *this;
    }

  //! @brief Insert breakpoint
  PWLU& insert
    ( double const& x )
    {
#ifdef MC__PWLU_CHECK
      if( _dx.empty() )
        throw Exceptions( Exceptions::RANGE );
      if( _dx.size() != _dy.size() )
        throw Exceptions( Exceptions::SIZE );
#endif
      if( isequal( x, _x0, options.BKPTATOL, options.BKPTRTOL ) )
        return *this; // no insertion
      double xi = _x0;
      for( auto idx=_dx.begin(), idy=_dy.begin(); idx!=_dx.end(); ++idx, ++idy ){
        if( isequal( x, xi + *idx, options.BKPTATOL, options.BKPTRTOL ) )
          return *this; // no insertion
        else if( xi + *idx <= x ){
          xi += *idx;
          continue;
        }
        double dx = xi + *idx - x, dy = *idy;
        *idx -= dx;
        _dx.insert( ++idx, dx );
        _dy.insert( ++idy, dy );
        break;
      }
      return *this;
    }

  //! @brief Evaluate estimator at point
  double val
    ( double const& x )
    const
    {
#ifdef MC__PWLU_CHECK
      if( _dx.empty() )
        throw Exceptions( Exceptions::RANGE );
      if( _dx.size() != _dy.size() )
        throw Exceptions( Exceptions::SIZE );
#endif
      double xi = _x0, yi = _y0;
      for( auto idx=_dx.cbegin(), idy=_dy.cbegin(); idx!=_dx.cend(); ++idx, ++idy ){
        if( xi + *idx < x ){
          xi += *idx;
          yi += *idx * *idy;
          continue;
        }
        //std::cout << "x: " << x << "  xi: " << xi << "  yi: " << yi << "  y: " << yi + (x-xi) * *idy << std::endl;
        return yi + (x-xi) * *idy;
      }
      //std::cout << "x: " << x << "  xi: " << xi << "  yi: " << yi << std::endl;
      if( isequal( xi, x, options.BKPTATOL, options.BKPTRTOL ) )
        return yi;
      throw Exceptions( Exceptions::EXTRAPOL );
    }

  //! @brief Evaluate estimator upper bound at point
  double l
    ( double const& x )
    const
    {
      return val( x );
    }

  //! @brief Evaluate estimator upper bound at point
  double u
    ( double const& x )
    const
    {
      return val( x );
    }
    
  //! @brief Evaluate lower range
  double l
    ()
    const
    {
#ifdef MC__PWLU_CHECK
      if( _dx.empty() )
        throw Exceptions( Exceptions::RANGE );
      if( _dx.size() != _dy.size() )
        throw Exceptions( Exceptions::SIZE );
#endif
      double l = _y0, yi = _y0;
      for( auto idx=_dx.cbegin(), idy=_dy.cbegin(); idx!=_dx.cend(); ++idx, ++idy )
        if( (yi += *idx * *idy) < l ) l = yi;
      return l;
    }

  //! @brief Evaluate upper range
  double u
    ()
    const
    {
#ifdef MC__PWLU_CHECK
      if( _dx.empty() )
        throw Exceptions( Exceptions::RANGE );
      if( _dx.size() != _dy.size() )
        throw Exceptions( Exceptions::SIZE );
#endif
      double u = _y0, yi = _y0;
      for( auto idx=_dx.cbegin(), idy=_dy.cbegin(); idx!=_dx.cend(); ++idx, ++idy )
        if( (yi += *idx * *idy) > u ) u = yi;
      return u;
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
#ifdef MC__PWLU_CHECK
      if( _dx.size() != _dy.size() )
        throw Exceptions( Exceptions::SIZE );
#endif
      double xi = _x0, yi = _y0;
      os << "{" << std::scientific << std::setprecision(dispnum) << std::right
         << std::setw(dispnum+8) << xi << ":" << std::setw(dispnum+8) << yi;
      for( auto idx=_dx.cbegin(), idy=_dy.cbegin(); idx!=_dx.cend(); ++idx, ++idy )
        os << ", " << std::setw(dispnum+8) << (xi += *idx) << ":" << std::setw(dispnum+8) << (yi += *idx * *idy);
      return os << " }";
    }

  //! @brief Retreive number of segments
  size_t size
    ()
    const
    {
#ifdef MC__PWLU_CHECK
      if( _dx.size() != _dy.size() )
        throw Exceptions( Exceptions::SIZE );
#endif
      return _dx.size();
    }

  //! @brief Retreive distances between breakpoints
  std::vector<double> const& dx
    ()
    const
    { return _dx; }

  //! @brief Retreive slopes between breakpoints
  std::vector<double> const& dy
    ()
    const
    { return _dy; }

  //! @brief Retreive lower abscissa
  double xL
    ()
    const
    { return _x0; }

  //! @brief Retreive upper abscissa
  double xU
    ()
    const
    {
#ifdef MC__PWLU_CHECK
      if( _dx.empty() )
        throw Exceptions( Exceptions::RANGE );
      if( _dx.size() != _dy.size() )
        throw Exceptions( Exceptions::SIZE );
#endif
      double xi = _x0;
      for( auto idx=_dx.cbegin(); idx!=_dx.cend(); ++idx )
        xi += *idx;
      return xi;
    }

  //! @brief Retreive initial ordinate
  double yL
    ()
    const
    { return _y0; }

  //! @brief Retreive/set distances between breakpoints
  std::vector<double>& dx
    ()
    { return _dx; }

  //! @brief Retreive/set slopes between breakpoints
  std::vector<double>& dy
    ()
    { return _dy; }

  //! @brief Retreive/set lower abscissa
  double& xL
    ()
    { return _x0; }

  //! @brief Retreive/set initial ordinate
  double& yL
    ()
    { return _y0; }

  PWLU& operator=
    ( PWLU const& var )
    {
      _x0 = var._x0;
      _y0 = var._y0;
      _dx = var._dx;
      _dy = var._dy;
      return *this;
    }

  PWLU& operator=
    ( PWLU && var )
    {
      _x0 = std::move( var._x0 );
      _y0 = std::move( var._y0 );
      _dx = std::move( var._dx );
      _dy = std::move( var._dy );
      return *this;
    }

  PWLU& operator+=
    ( double const& cst )
    {
      _y0 += cst;
      return *this;
    }

  PWLU& operator-=
    ( double const& cst )
    {
      _y0 -= cst;
      return *this;
    }
 
  PWLU& operator*=
    ( double const& cst )
    {
      if( cst == 1. )
        return *this;
      _y0 *= cst;
      for( auto& dyi : _dy )
        dyi *= cst;
      return *this;
    }

  PWLU& operator/=
    ( double const& cst )
    {
      if( cst == 1. )
        return *this;
      _y0 /= cst;
      for( auto& dyi : _dy )
        dyi /= cst;
      return *this;
    }
           
  PWLU& operator+=
    ( PWLU const& var )
    {
#ifdef MC__PWLU_CHECK
      if( _dx.empty() || var._dx.empty() )
        throw Exceptions( Exceptions::RANGE );
      if( _dx.size() != _dy.size() || var._dx.size() != var._dy.size() )
        throw Exceptions( Exceptions::SIZE );
#endif

      if( !isequal( _x0, var._x0 ) )
        throw Exceptions( Exceptions::RANGE );

      _dx.reserve( _dx.size() + var._dx.size() );
      _dy.reserve( _dy.size() + var._dy.size() );

      _y0 += var._y0;
      double d2 = 0.;
      auto idx1 = _dx.begin(),      idy1 = _dy.begin();
      auto idx2 = var._dx.cbegin(), idy2 = var._dy.cbegin();
      for( ; idx1 != _dx.end() && idx2 != var._dx.cend(); ++idx1, ++idy1 ){
        //std::cout << std::scientific << std::setprecision(14) << "*idx1=" << *idx1 << "  *idx2-d2=" << *idx2-d2 << std::endl;

        if( isequal( *idx1, *idx2-d2, options.BKPTATOL, options.BKPTRTOL ) ){
          //std::cout << "case 1\n";
          *idy1 += *idy2; // update slope
          ++idx2;
          ++idy2;
          d2 = 0.;
        }

        else if( *idx1 > *idx2-d2 ){
          //std::cout << "case 2\n";
          auto idx11 = idx1, idy11 = idy1;
          idx11 = _dx.insert( ++idx11, *idx1 - *idx2 + d2 );
          idy11 = _dy.insert( ++idy11, *idy1 );
          *idx1 -= *idx11;
          *idy1 += *idy2; // update slope
          ++idx2;
          ++idy2;
          d2 = 0.;
        }

        else{
          //std::cout << "case 3:\n";
          *idy1 += *idy2; // update slope
          d2 += *idx1;
        }
      }

      if( idx1 != _dx.end() ){//&& isequal( *idx1, 0., options.BKPTATOL, options.BKPTRTOL ) ){
#ifdef MC__PWLU_CHECK
        std::cout << " xU = " << xU() << "  var.xU = " << var.xU() << std::endl;
#endif
        *idy1 += *(--idy2);
        ++idx1;
        ++idy1;
        ++idy2;
      }

      if( idx1 != _dx.end() ){//|| idx2 != var._dx.cend() ){
#ifdef MC__PWLU_CHECK
        if( idx1 != _dx.end() )
          std::cout << " dx = " << *idx1 << std::endl;
        if( idx2 != var._dx.cend() ){
          display( std::cout );
          var.display( std::cout );
          std::cout << " var.dx = " << *idx2 << std::endl;
        }
        std::cout << " xU = " << xU() << "  var.xU = " << var.xU() << std::endl;
#endif
        throw Exceptions( Exceptions::RANGE );
      }

      return *this;
    }

  PWLU& operator-=
    ( PWLU const& var )
    {
#ifdef MC__PWLU_CHECK
      if( _dx.empty() || var._dx.empty() )
        throw Exceptions( Exceptions::RANGE );
      if( _dx.size() != _dy.size() || var._dx.size() != var._dy.size() )
        throw Exceptions( Exceptions::SIZE );
#endif

      if( !isequal( _x0, var._x0 ) )
        throw Exceptions( Exceptions::RANGE );

      _dx.reserve( _dx.size() + var._dx.size() );
      _dy.reserve( _dy.size() + var._dy.size() );

      _y0 -= var._y0;
      double d2 = 0.;
      auto idx1 = _dx.begin(),      idy1 = _dy.begin();
      auto idx2 = var._dx.cbegin(), idy2 = var._dy.cbegin();
      for( ; idx1 != _dx.end() && idx2 != var._dx.cend(); ++idx1, ++idy1 ){

        if( isequal( *idx1, *idx2-d2, options.BKPTATOL, options.BKPTRTOL ) ){
          *idy1 -= *idy2; // update slope
          ++idx2;
          ++idy2;
          d2 = 0.;
        }

        else if( *idx1 > *idx2-d2 ){
          auto idx11 = idx1, idy11 = idy1;
          idx11 = _dx.insert( ++idx11, *idx1 - *idx2 + d2 );
          idy11 = _dy.insert( ++idy11, *idy1 );
          *idx1 -= *idx11;
          *idy1 -= *idy2; // update slope
          ++idx2;
          ++idy2;
          d2 = 0.;
        }

        else{
          *idy1 -= *idy2; // update slope
          d2 += *idx1;
        }
      }

      if( idx1 != _dx.end() ){//&& isequal( *idx1, 0., options.BKPTATOL, options.BKPTRTOL ) ){
#ifdef MC__PWLU_CHECK
        std::cout << " xU = " << xU() << "  var.xU = " << var.xU() << std::endl;
#endif
        *idy1 -= *(--idy2);
        ++idx1;
        ++idy1;
        ++idy2;
      }

      if( idx1 != _dx.end() ){//|| idx2 != var._dx.cend() ){
#ifdef MC__PWLU_CHECK
        if( idx1 != _dx.end() )
          std::cout << " dx = " << *idx1 << std::endl;
        if( idx2 != var._dx.cend() ){
          var.display( std::cout );
          std::cout << " var.dx = " << *idx2 << std::endl;
        }
        std::cout << " xU = " << xU() << "  var.xU = " << var.xU() << std::endl;
#endif
        throw Exceptions( Exceptions::RANGE );
      }
      
      return *this;
    }

  template <typename UNIV, typename DUNIV>
  PWLU& compose
    ( UNIV const& f, DUNIV const& df, bool const under, bool const cvx, bool const inc=true )
    {
      if( cvx && !under || !cvx && under ){
        double f1 = f( _y0 ), f2;
        double yi = _y0;
        _y0 = f1;
        //std::cout << yi << ": " << f1 << std::endl;
        for( auto idx=_dx.begin(), idy=_dy.begin(); idx!=_dx.end(); ++idx, ++idy ){
          yi += *idx * *idy;
          f2 = f( yi );
          *idy = ( f2 - f1 ) / *idx; 
          //std::cout << yi << ": " << f2 << "  " << *idy << std::endl;
          std::swap( f1, f2 );
        }
      }

      else{
        _dx.reserve( 2*_dx.size() );
        _dy.reserve( 2*_dy.size() );

        double f1  = f( _y0 ), f2;
        double df1 = df( _y0 ), df2;
        double yi = _y0;
        _y0 = f1;
        //std::cout << yi << ": " << f1 << std::endl;
        for( auto idx=_dx.begin(), idy=_dy.begin(); idx!=_dx.end(); ++idx, ++idy ){
          yi += *idx * *idy;
          f2  = f( yi );
          df2 = df( yi );

          //if( isequal( f1 + 0.5 * df1 * *idy * *idx, f2 - 0.5 * df2 * *idy * *idx, options.BKPTATOL, options.BKPTRTOL ) ){
          if( isequal( f1, f2 - df2 * *idy * *idx, options.BKPTATOL, options.BKPTRTOL )
           && isequal( f1 + df1 * *idy * *idx, f2, options.BKPTATOL, options.BKPTRTOL ) ){
          //if( isequal( df1 * *idy * *idx, df2 * *idy * *idx, options.BKPTATOL, options.BKPTRTOL ) ){
            *idy = ( f2 - f1 ) / *idx; // use secant approximation
          }
          //double denom = *idy * ( df2 - df1 ); 
          //if( isequal( denom, 0., options.BKPTATOL, options.BKPTRTOL ) ){
          //  *idy = ( f2 - f1 ) / *idx;
          //}
          else{
            double dx = ( f2 - f1 - *idx * *idy * df1 ) / ( *idy * ( df2 - df1 ) );
            if( dx < 0. || dx > *idx ){
#ifdef MC__PWLU_CHECK
              std::cout << " dx = " << dx << " not in [0," << *idx << "]" << std::endl;
#endif
              *idy = ( f2 - f1 ) / *idx; // use secant approximation
            }
            else{
              auto idx0 = idx, idy0 = idy;
              idx = _dx.insert( ++idx, dx );
              idy = _dy.insert( ++idy, *idy0 * df2 );
              *idx0 -= *idx;
              *idy0 *= df1;
            }
          }
          std::swap( f1, f2 );
          std::swap( df1, df2 );
        }
      }

      return *this;
    }

  PWLU& max
    ( double const& c )
    {
      if( _dx.empty() )
        return *this;

#ifdef MC__PWLU_CHECK
      if( _dx.size() != _dy.size() )
        throw Exceptions( Exceptions::SIZE );
#endif
      _dx.reserve( 2*_dx.size() );
      _dy.reserve( 2*_dy.size() );

      double y1 = _y0, y2;
      if( _y0 < c ) _y0 = c;
      for( auto idx=_dx.begin(), idy=_dy.begin(); idx!=_dx.end(); ++idx, ++idy ){
        y2 = y1 + *idx * *idy;
        if( ( y1 < c || isequal( y1, c, options.BKPTATOL, options.BKPTRTOL ) )
         && ( y2 < c || isequal( y2, c, options.BKPTATOL, options.BKPTRTOL ) ) ){
          *idy = 0.;
        }
        else if( y1 < c ){
          // y2 - (*idy)*dx2 = c -> dx2 = (y2-c)/(*idy)
          auto idx0 = idx, idy0 = idy;
          idx = _dx.insert( ++idx, ( y2 - c ) / *idy0 );
          idy = _dy.insert( ++idy, *idy0 );
          *idx0 -= *idx;
          *idy0 = 0.;
        }
        else if( y2 < c ){
          // y2 - (*idy)*dx2 = c -> dx2 = (y2-c)/(*idy)
          auto idx0 = idx, idy0 = idy;
          idx = _dx.insert( ++idx, ( y2 - c ) / *idy0 );
          idy = _dy.insert( ++idy, 0. );
          *idx0 -= *idx;
        }
        std::swap( y1, y2 );
      }
      return *this;
    }

  PWLU& min
    ( double const& c )
    {
      if( _dx.empty() )
        return *this;

#ifdef MC__PWLU_CHECK
      if( _dx.size() != _dy.size() )
        throw Exceptions( Exceptions::SIZE );
#endif
      _dx.reserve( 2*_dx.size() );
      _dy.reserve( 2*_dy.size() );

      double y1 = _y0, y2;
      if( _y0 > c ) _y0 = c;
      for( auto idx=_dx.begin(), idy=_dy.begin(); idx!=_dx.end(); ++idx, ++idy ){
        y2 = y1 + *idx * *idy;
        if( ( y1 > c || isequal( y1, c, options.BKPTATOL, options.BKPTRTOL ) )
         && ( y2 > c || isequal( y2, c, options.BKPTATOL, options.BKPTRTOL ) ) ){
          *idy = 0.;
        }
        else if( y1 > c ){
          // y2 - (*idy)*dx2 = c -> dx2 = (y2-c)/(*idy)
          auto idx0 = idx, idy0 = idy;
          idx = _dx.insert( ++idx, ( y2 - c ) / *idy0 );
          idy = _dy.insert( ++idy, *idy0 );
          *idx0 -= *idx;
          *idy0 = 0.;
        }
        else if( y2 > c ){
          // y2 - (*idy)*dx2 = c -> dx2 = (y2-c)/(*idy)
          auto idx0 = idx, idy0 = idy;
          idx = _dx.insert( ++idx, ( y2 - c ) / *idy0 );
          idy = _dy.insert( ++idy, 0. );
          *idx0 -= *idx;
        }
        std::swap( y1, y2 );
      }
      return *this;
    }
    
  PWLU& clean
    ( bool const under )
    {
      // Initial stage
      auto idx=_dx.begin(), idy=_dy.begin();
      for( ; ; ){
        if( isequal( *idx, 0., options.BKPTATOL, options.BKPTRTOL ) ){
          double dx = *idx, dy = *idy;
          idx = _dx.erase( idx );
          idy = _dy.erase( idy );
          if( ( under && dy > *idy ) || ( !under && dy < *idy ) )
            _average( *idx, *idy, dx, dy);
          else{
            *idx += dx;
            _y0  += dx * ( dy - *idy );
          }
          continue;
        }
        break;
      }
        
      // Sunsequent stages
      for( ; ; ){
        if( isequal( *idx, 0., options.BKPTATOL, options.BKPTRTOL ) ){
          double dx = *idx, dy = *idy;
          idx = _dx.erase( idx );
          idy = _dy.erase( idy );
          // Final stage
          if( idy == _dy.end() ){
            --idx; 
            --idy;
            if( ( under && dy < *idy ) || ( !under && dy > *idy ) )
              _average( *idx, *idy, dx, dy);
            else
              *idx += dx;
            break;
          }
          // Intermediate stage
          if( ( under && dy > *idy ) || ( !under && dy < *idy ) )
            _average( *idx, *idy, dx, dy);
          else{
            --idx;
            --idy;
            if( ( under && *idy < dy ) || ( !under && *idy > dy ) ){
              *idx += dx;
              double oy = dx *( dy - *idy );
              ++idx;
              ++idy;
              *(idy) += oy / *idx;
            }
            else{
              _average( *idx, *idy, dx, dy);
              ++idx;
              ++idy;
            }
          }
          continue;
        }
        else{
          ++idx;
          ++idy;
          if( idy == _dy.end() ) break;
        }
      }
      
      return *this;
    }

  PWLU& merge
    ( bool const under )
    {
      auto idx = _dx.begin(), idy = _dy.begin();
      for( ; ; ){
        auto idx2 = idx, idy2 = idy;
        ++idx2;
        ++idy2;
        if( idy2 == _dy.end() ) break;
        if( isequal( *idy, *idy2, options.BKPTATOL, options.BKPTRTOL ) ){
          *idx += *idx2;
          if( ( under && *idy > *idy2 ) || ( !under && *idy < *idy2 ) )
            *idy += *idx2 / *idx * ( *idy2 - *idy );
          _dx.erase( idx2 );
          _dy.erase( idy2 );
        }
        else{
          ++idx;
          ++idy;
          if( idy == _dy.end() ) break;
        }
      }
      return *this;
    }

  PWLU& reduce
    ( bool const under, size_t const nseg ){
      clean( under );
      if( !nseg || _dx.size() <= nseg || _dx.size() < 3 )
        return *this;

      // use standard library methods std::erase with std::min_element
      if( options.REDUCEMETH <= 0 || _dx.size() <= nseg + options.REDUCEMETH )
        return _reduce_std( under, nseg );

      // use tailored heap-based implementation
      return _reduce_heap( under, nseg );
    }

private:

  static void _average
    ( double& dx, double& s, double const& dx0, double const& s0 )
    {
      s  *= dx;
      s  += (s0 * dx0);
      dx += dx0;     
      s /= dx;
    }    

  static void _loss
    ( bool const isUnder, double const& slpDiff1, double const& slpDiff2,
      double const& incA0, double const& slpA0, 
      double const& incA1, double const& slpA1,
      double const& incA2, double const& slpA2,
      std::tuple<double,int,double,double>& loss );

  static double _area
    ( double const& inc1, double const& inc2, double const& slpdiff)
    {
      return inc1*inc2*slpdiff;
    }

  static std::pair<double,double> _intersect
    ( double const& x1, double const& y1, double const& x2, double const& y2,
      double const& x3, double const& y3, double const& x4, double const& y4 )
    {
      double const d = (y1-y2)*(x4-x3)+(y4-y3)*(x2-x1);
      // Avoid division by zero
      double x = !isequal( d, 0. )? ((y1*x2-y2*x1)*(x4-x3)+(y4*x3-y3*x4)*(x2-x1))/d: 0.5*(x2+x3);
      if( x < x2 || x > x3 ) x = 0.5*(x2+x3);
      return std::make_pair( x, (y1*(x2-x)+y2*(x-x1))/(x2-x1) );
    }

  static void _underestimate_in_window3
    ( const std::vector<double> & x, const std::vector<double> & y, const unsigned int pos, 
      const double p1AreaLoss, const double p2AreaLoss, std::vector<double> & data, const bool isUnder );

  static bool _cmp_indexed_areaLoss
    ( const std::pair<double, int>& a, const std::pair<double, int>& b );

  static std::pair<double, int> _heap_pop
    ( std::vector<std::pair<double, int>> & heap, std::vector<int> & hind );

  static void _modify_heap
    ( std::vector<std::pair<double, int>> & heap, const int pos, std::vector<int> & hind );

  static void _heap_bubble_up
    ( std::vector<std::pair<double, int>> & heap, int pos, std::vector<int> & hind );

  static void _heap_bubble_down
    ( std::vector<std::pair<double, int>> & heap, int pos, std::vector<int> & hind );

  PWLU& _reduce_heap
    ( bool const under, size_t const nseg );

  PWLU& _reduce_std
    ( bool const under, size_t const nseg );

};

////////////////////////////////////////////////////////////////////////

inline PWLU::Options PWLU::options;

inline
PWLU& PWLU::_reduce_std
( bool const under, size_t const nseg )
{
#ifdef MC__PWLU_CHECK
  if( _dx.size() != _dy.size() )
    throw Exceptions( Exceptions::SIZE );
#endif

  // Remove any degenerate segments, reporting on adjacent segments as necessary
  clean( under );
  if( !nseg || _dx.size() <= nseg || _dx.size() < 3 )
    return *this;  

#ifdef MC__PWLU_DEBUG
  std::cout << std::scientific << std::setprecision(16);
  std::cout << "INITIAL VARIABLE AFTER CLEAN:" << std::endl;
  std::cout << "under:" << " " << under << std::endl;
  std::cout << "x0:" << " " << _x0;
  std::cout << std::endl << "dx:";
  for( auto dxi : _dx ) std::cout << " " << dxi;
  std::cout << std::endl << "y0:" << " " << _y0;
  std::cout << std::endl << "dy:";
  for( auto dyi : _dy ) std::cout << " " << dyi;
  std::cout << std::endl;
#endif

  // Compute minimal loss for each 3-segment combination
  static thread_local std::vector<std::tuple<double,int,double,double>> vloss;  
  vloss.resize( _dx.size()-2 );  

  static thread_local double lcv[3];
  double *lcv0 = lcv, *lcv2 = lcv0+1, *lcv3 = lcv2+1;

  auto idx = _dx.begin(), idy = _dy.begin();
  double *incA0 = &(*idx), *incA1 = &(*(++idx)), *incA2 = &(*(++idx));
  double *slpA0 = &(*idy), *slpA1 = &(*(++idy)), *slpA2 = &(*(++idy));
  lcv[1]  = under? (*slpA0)-(*slpA1):(*slpA1)-(*slpA0); 
  lcv[2]  = under? (*slpA1)-(*slpA2):(*slpA2)-(*slpA1);

  for( auto ploss = vloss.data(); ; ploss++ ){
    _loss( under, *lcv2, *lcv3, *incA0, *slpA0, *incA1, *slpA1, *incA2, *slpA2, *ploss );
#ifdef MC__PWLU_DEBUG
    std::cout << "loss: " << std::get<0>( *ploss ) << std::endl;
#endif
    
    if( ++idx == _dx.end() ) break;
    lcv0 = lcv2, lcv2 = lcv3, lcv3 = lcv0;
    incA0 = incA1; incA1 = incA2; incA2 = &(*idx);
    slpA0 = slpA1; slpA1 = slpA2; slpA2 = &(*(++idy));   
    *lcv3 = under? (*slpA1)-(*slpA2):(*slpA2)-(*slpA1); 
  }

  // Main reduction loop
  for( ; ; ){

#ifdef MC__PWLU_DEBUG
    std::cout << "dx:";
    for( auto dxi : _dx ) std::cout << " " << dxi;
    std::cout << std::endl << "dy:";
    for( auto dyi : _dy ) std::cout << " " << dyi;
    std::cout << std::endl;

    std::cout << "Candidate 3-segments:\n";
    for( auto itloss=vloss.cbegin(); itloss!=vloss.cend(); ++itloss )
      std::cout << std::get<0>(*itloss) << " "
                << std::get<1>(*itloss) << " "
                << std::get<2>(*itloss) << " "
                << std::get<3>(*itloss) << " " << std::endl;
#endif

    // Determine minimal element
    auto ilossopt = std::min_element(
      vloss.cbegin(), vloss.cend(),
      [=]( std::tuple<double,int,double,double> const& a, std::tuple<double,int,double,double> const& b )
      { return std::get<0>(a) < std::get<0>(b); }
    );

    size_t pos = std::distance( vloss.cbegin(), ilossopt );
    auto idxopt = _dx.begin() + pos, idyopt = _dy.begin() + pos;

#ifdef MC__PWLU_DEBUG
    std::cout << "Minimum-segments at" << pos << " :\n";
    std::cout << std::get<0>(*ilossopt) << " "
              << std::get<1>(*ilossopt) << " "
              << std::get<2>(*ilossopt) << " "
              << std::get<3>(*ilossopt) << " " << std::endl;
#endif

    // Remove minimal element
    switch (std::get<1>(*ilossopt)) {
        case 0:       // remove two vtx then add one vtx == remove mid seg
            *idxopt = std::get<2>(*ilossopt);
            idxopt = _dx.erase( ++idxopt );
            idyopt = _dy.erase( ++idyopt );
            *idxopt = std::get<3>(*ilossopt);   
            break;
        case 2:        // remove one vtx == remove mid seg and update one seg
            std::advance( idxopt, 1 );
            std::advance( idyopt, 1 );
        case 1:       // remove one vtx == remove mid seg and update one seg
            idxopt = _dx.erase( idxopt );
            idyopt = _dy.erase( idyopt );
            *idxopt = std::get<3>(*ilossopt);
            *idyopt = std::get<2>(*ilossopt);  
            break;
        default:
            std::cout << "error in reduce" << std::endl;
            break;
    }

    if( _dx.size() <= nseg || _dx.size() < 3 ){
#ifdef MC__PWLU_DEBUG
      std::cout << "dx:";
      for( auto dxi : _dx ) std::cout << " " << dxi;
      std::cout << std::endl << "dy:";
      for( auto dyi : _dy ) std::cout << " " << dyi;
      std::cout << std::endl;
#endif
      return *this;
    }

    // Update minimal loss
    vloss.erase( ilossopt );
    auto itloss = vloss.begin();

    idx  = _dx.begin(),     idy  = _dy.begin();
    size_t p = 0;
    if( pos > 2 ){
      p = pos - 2;
      idx += p;    
      idy += p;
      itloss += p;
    }

    incA0 = &(*idx); incA1 = &(*(++idx)); incA2 = &(*(++idx));
    slpA0 = &(*idy); slpA1 = &(*(++idy)); slpA2 = &(*(++idy));

    lcv[1]  = under? (*slpA0)-(*slpA1):(*slpA1)-(*slpA0); 
    lcv[2]  = under? (*slpA1)-(*slpA2):(*slpA2)-(*slpA1);
    _loss( under, lcv[1], lcv[2], *incA0, *slpA0, *incA1, *slpA1, *incA2, *slpA2, *itloss );
#ifdef MC__PWLU_DEBUG
    std::cout << std::get<0>( *itloss ) << std::endl;
#endif
    if( ++idx==_dx.end() || ++p>pos+2 ) continue;

    incA0 = &(*idx);
    slpA0 = &(*(++idy));
    lcv[0] = under? (*slpA2)-(*slpA0):(*slpA0)-(*slpA2); 
    _loss( under, lcv[2], lcv[0], *incA1, *slpA1, *incA2, *slpA2, *incA0, *slpA0,*(++itloss) );    
#ifdef MC__PWLU_DEBUG
    std::cout << std::get<0>( *itloss ) << std::endl;
#endif
    if( ++idx==_dx.end() || ++p>pos+2 ) continue;

    incA1 = &(*idx);
    slpA1 = &(*(++idy));
    lcv[1] = under? (*slpA0)-(*slpA1):(*slpA1)-(*slpA0); 
    _loss( under, lcv[0], lcv[1],  *incA2, *slpA2, *incA0, *slpA0, *incA1, *slpA1, *(++itloss) );      
#ifdef MC__PWLU_DEBUG
    std::cout << std::get<0>( *itloss ) << std::endl;
#endif
    if( ++idx==_dx.end() || ++p>pos+2 ) continue;

    incA2 = &(*idx);
    slpA2 = &(*(++idy));
    lcv[2] = under? (*slpA1)-(*slpA2):(*slpA2)-(*slpA1); 
    _loss( under, lcv[1], lcv[2], *incA0, *slpA0, *incA1, *slpA1, *incA2, *slpA2, *(++itloss) );      
#ifdef MC__PWLU_DEBUG
    std::cout << std::get<0>( *itloss ) << std::endl;
#endif
    if( ++idx==_dx.end() || ++p>pos+2 ) continue;

    incA0 = &(*idx);
    slpA0 = &(*(++idy));
    lcv[0] = under? (*slpA2)-(*slpA0):(*slpA0)-(*slpA2); 
    _loss( under, lcv[2], lcv[0], *incA1, *slpA1, *incA2, *slpA2, *incA0, *slpA0, *(++itloss) );      
#ifdef MC__PWLU_DEBUG
    std::cout << std::get<0>( *itloss ) << std::endl;
#endif
    if( ++idx==_dx.end() || ++p>pos+2 ) continue;
  }

  return *this;
}

inline
void PWLU::_loss
( bool const isUnder, double const& slpDiff1, double const& slpDiff2,
  double const& incA0, double const& slpA0, 
  double const& incA1, double const& slpA1,
  double const& incA2, double const& slpA2,
  std::tuple<double,int,double,double>& loss )
{
  if( slpDiff1 >= 0. && slpDiff2 >= 0. ){ // ccv ccv 
    //std::cerr << "case 1" << std::endl;
    double const p1AreaLoss = slpDiff1 * incA0 * incA1;
    double const p2AreaLoss = slpDiff2 * incA1 * incA2;

    if( p1AreaLoss < p2AreaLoss ){
      double const deltaXL = incA1 + incA0;
      double const ratioL = incA0/deltaXL;
      double const slpXL = slpA0*ratioL + slpA1*(1 - ratioL);
      loss = std::make_tuple( p1AreaLoss, 1, slpXL, deltaXL);
    }

    else{
      double const deltaXR = incA1 + incA2;
      double const ratioL = incA1/deltaXR;
      double const slpXR = slpA1*ratioL + slpA2*(1 - ratioL);            
      loss = std::make_tuple( p2AreaLoss, 2, slpXR, deltaXR);
    }
  }        

  else if( slpDiff2 >= 0. ){ // cvx ccv
    //std::cerr << "case 2" << std::endl;
    double const deltaXR = incA1 + incA2;
    double const ratioL = incA1/deltaXR;
    double const ratioR = 1 - ratioL;    
    double const slpXR = slpA1*ratioL + slpA2*ratioR;
    double slpDiff3 = isUnder?slpA0 - slpXR:slpXR - slpA0;        
    
    if( slpDiff3 > 1e-14 ){
      //std::cerr << "case 2.1" << std::endl;
      double const _tmpL = ((slpA1-slpA0)*incA1)/(slpA0-slpA2);
      double const _tmpR = incA2 - _tmpL;
      if( _tmpR > options.BKPTATOL ){
        loss = std::make_tuple( slpDiff2*incA1*_tmpL, 0, incA0 + incA1 + _tmpL, _tmpR );
        return;
      }
    }
    loss = std::make_tuple( slpDiff2*incA1*incA2, 2, slpXR, deltaXR);
  }

  else if( slpDiff1 >= 0. ){ // ccv cvx
    //std::cerr << "case 3" << std::endl;    
    double const deltaXL = incA1 + incA0;
    double const ratioL = incA0/deltaXL;
    double const ratioR = 1 - ratioL;
    double const slpXL = slpA0*ratioL + slpA1*ratioR;
    double slpDiff3 = isUnder?slpXL - slpA2:slpA2 - slpXL;  
    
    if( slpDiff3 > 1e-14 ){
      //std::cerr << "case 3.1" << std::endl;
      double const _tmpR = ((slpA1-slpA2)*incA1)/(slpA2-slpA0);
      double const _tmpL = incA0 - _tmpR;     
      if( _tmpL > options.BKPTATOL ){
        loss = std::make_tuple( (slpA0-slpA1)*incA1*_tmpR, 0, _tmpL, _tmpR + incA1 + incA2 );
        return;
      }                
    
    }
    loss = std::make_tuple( slpDiff1*incA0*incA1, 1, slpXL, deltaXL );    
  }

  else{
    //std::cerr << "case 4" << std::endl;
    double slpDiff3 = isUnder?slpA0-slpA2:slpA2-slpA0;  
    double const ratioR = slpDiff1/slpDiff3;
    double const ratioL = slpDiff2/slpDiff3;
    double const deltaXR = ratioR*incA1;
    double const deltaXL = ratioL*incA1;  
    loss = std::make_tuple( -slpDiff3*deltaXL*deltaXR, 0, deltaXL + incA0, deltaXR + incA2 );
  }

  return;
}

//! static member function 
//! @brief Linear interpolation for aligning a univariate linear estimator with a univariate piecewise linear estimator w.r.t. the same variable
//! @param[in] pos is
//! @param[in] p1AreaLoss is 
//! @param[in] p2AreaLoss is t
//! @param[out] data 
inline
void PWLU::_underestimate_in_window3
( const std::vector<double>& incA, const std::vector<double>& slpA, const unsigned int pos, 
  const double slpDiff1, const double slpDiff2, std::vector<double>& data, const bool isUnder )
{

  double & deltaXL  = data[0];
  double & deltaXR  = data[1];
  double & slpXL = data[2];     
  double & slpXR = data[3];     
  double & aLoss = data[4]; 

  // T slpDiff1 = isUnder?(slpA[0]-slpA[1]):(slpA[1]-slpA[0]); 
  // T slpDiff2 = isUnder?(slpA[1]-slpA[2]):(slpA[2]-slpA[1]);
  unsigned int pos1(pos+1),pos2(pos+2);
  if (slpDiff1 >= 0. && slpDiff2 >= 0.){ // ccv ccv 
    //std::cerr << "case 1" << std::endl;

    double const p1AreaLoss = slpDiff1*incA[pos]*incA[pos1];
    double const p2AreaLoss = slpDiff2*incA[pos1]*incA[pos2];
    if (p1AreaLoss < p2AreaLoss){
        deltaXL = incA[pos1] + incA[pos];
        deltaXR = incA[pos2];
        double const ratioL = incA[pos]/deltaXL;
        double const ratioR = 1 - ratioL;
        slpXL = slpA[pos]*ratioL + slpA[pos1]*ratioR;
        slpXR = slpA[pos2];            
        aLoss = p1AreaLoss;
    }
    else{
        deltaXL = incA[pos];
        deltaXR = incA[pos1] + incA[pos2];
        double const ratioL = incA[pos1]/deltaXR;
        double const ratioR = 1 - ratioL;
        slpXL = slpA[pos];
        slpXR = slpA[pos1]*ratioL + slpA[pos2]*ratioR;            
        aLoss = p2AreaLoss;
    }
  }        
  else if (slpDiff2 >= 0.){ // cvx ccv
    //std::cerr << "case 2" << std::endl;

    deltaXL = incA[pos];
    deltaXR = incA[pos1] + incA[pos2];
    double const ratioL = incA[pos1]/deltaXR;
    double const ratioR = 1 - ratioL;    
    slpXL = slpA[pos];
    slpXR = slpA[pos1]*ratioL + slpA[pos2]*ratioR;
    aLoss = slpDiff2*incA[pos1]*incA[pos2];
    double slpDiff3 = isUnder?slpA[pos] - slpXR:slpXR - slpA[pos];        

    if (slpDiff3 > 1e-14){
      //std::cerr << "case 2.1" << std::endl;
      double const _tmpL = ((slpA[pos1]-slpA[pos])*incA[pos1])/(slpA[pos]-slpA[pos2]);
      double const _tmpR = incA[pos2] - _tmpL;
      if( _tmpR > options.BKPTATOL ){
        deltaXR = _tmpR; 
        slpXR = slpA[pos2];
        aLoss = slpDiff2*incA[pos1]*_tmpL; 
        deltaXL = incA[pos] + incA[pos1] + _tmpL;
      }
    }
  }
  else if (slpDiff1 >= 0.){ // ccv cvx
    //std::cerr << "case 3" << std::endl;

    deltaXL = incA[pos1] + incA[pos];
    deltaXR = incA[pos2];
    double const ratioL = incA[pos]/deltaXL;
    double const ratioR = 1 - ratioL;
    slpXL = slpA[pos]*ratioL + slpA[pos1]*ratioR;
    slpXR = slpA[pos2];            
    aLoss = slpDiff1*incA[pos]*incA[pos1];
    double slpDiff3 = isUnder?slpXL - slpA[pos2]:slpA[pos2] - slpXL;  

    if (slpDiff3 > 1e-14){
      //std::cerr << "case 3.1" << std::endl;
      double const _tmpR = ((slpA[pos1]-slpA[pos2])*incA[pos1])/(slpA[pos2]-slpA[pos]);
      // if(_tmpR < 0){
      //   std::cout << isUnder << std::endl;
      //   std::cout << _tmpR << std::endl;
      //   std::cout << slpA[pos] << ", " << slpA[pos1] << ", " << slpA[pos2] << std::endl;
      //   std::cout << incA[pos] << ", " << incA[pos1] << ", " << incA[pos2] << std::endl;
      //   assert(_tmpR >=0.);
      // }
      double const _tmpL = incA[pos] - _tmpR;     
      if( _tmpL > options.BKPTATOL ){
        deltaXL = _tmpL;
        slpXL = slpA[pos];
        aLoss = (slpA[pos]-slpA[pos1])*incA[pos1]*_tmpR; 
        deltaXR = _tmpR + incA[pos1] + incA[pos2];
      }                

    }
  }
  else{
    //std::cerr << "case 4" << std::endl;
    double slpDiff3 = isUnder?slpA[pos]-slpA[pos2]:slpA[pos2]-slpA[pos];  
    double const ratioR = slpDiff1/slpDiff3;
    double const ratioL = slpDiff2/slpDiff3;
    
    slpXL = slpA[pos];
    slpXR = slpA[pos2];      
    deltaXR = ratioR*incA[pos1];
    deltaXL = ratioL*incA[pos1];  
    aLoss = -slpDiff3*deltaXL*deltaXR; 
    deltaXR = deltaXR + incA[pos2];
    deltaXL = deltaXL + incA[pos];
  }
#ifdef MC__CPWLUSG_DEBUG    
    if(deltaXR < 1e-15){
      std::cerr << deltaXR << std::endl;
      assert(deltaXR >= 1e-15);
    }
    if(deltaXL < 1e-15){
      std::cerr << deltaXL << std::endl;
      assert(deltaXL >= 1e-15);
    }    

    assert(std::fabs(slpXL) <= 50);
    assert(std::fabs(slpXR) <= 50);  
#endif 
}

//! @brief Linear interpolation for aligning a univariate linear estimator with a univariate piecewise linear estimator w.r.t. the same variable
//! @param[in] nbpsMax is
//! @param[in] isUnder is 
//! @param[out] first 
//! @param[out] second
inline
PWLU& PWLU::_reduce_heap
( bool const isUnder, size_t const nseg )
{
  // std::cout <<  "        in reduce" << std::endl;
  std::vector<double> & incA = _dx;
  std::vector<double> & slpA = _dy;
  const unsigned int n = _dx.size() + 1;
  const unsigned int k = nseg;

  // std::cout <<  "          max pieces = " << k << std::endl;
  // std::cout <<  "          current number of vertices = " << n << std::endl;

  // double const xRightEndpoint = first.back();
  // first[n-1] = 0.;
  // for (unsigned int i = 0; i < n-1; i++)
  //   first[n-1] += first[i];
  // first[n-1] = xRightEndpoint - first[n-1];  

#ifdef MC__CPWLUSG_DEBUG    
          for (unsigned int jj = 1; jj < first.size(); jj++)
            if(first[jj] < 1e-15){std::cout << jj << " in " << first.size() << first[jj] << std::endl; assert(first[jj] >= 1e-15);};
#endif

  double slpDiff1;// = isUnder?(slpA[0]-slpA[1]):(slpA[1]-slpA[0]); 
  double slpDiff2;// = isUnder?(slpA[1]-slpA[2]):(slpA[2]-slpA[1]);

  std::vector<double> data(5);
  std::vector<double> deltaXL(n-3);
  std::vector<double> deltaXR(n-3);
  std::vector<double> slpXL(n-3);
  std::vector<double> slpXR(n-3);  
  std::vector<std::pair<double, int>> areaLoss(n-3);

  std::vector<double> yval(n);
  yval[0] = _y0;//slpA[0];

  double yMinOrMax = _y0;//slpA[0];

  slpDiff1 = isUnder?(slpA[0]-slpA[1]):(slpA[1]-slpA[0]);//isUnder?(slpA[1]-slpA[2]):(slpA[2]-slpA[1]); 
  for(unsigned int j = 0; j < n - 3; j++){       
    slpDiff2 = isUnder?(slpA[j+1]-slpA[j+2]):(slpA[j+2]-slpA[j+1]);
    _underestimate_in_window3(_dx,_dy,j,slpDiff1,slpDiff2,data,isUnder);
    deltaXL[j] = data[0];
    deltaXR[j] = data[1];
    slpXL[j] = data[2];
    slpXR[j] = data[3];

#ifdef MC__CPWLUSG_DEBUG            
    assert(std::fabs(slpXL[j-1]) <= 50);
    assert(std::fabs(slpXR[j-1]) <= 50);   
#endif       

    areaLoss[j] = std::make_pair(data[4],j);// AreaLoss(aLoss,j)
    //std::cout << "j: " << j << " area: " << data[2] << std::endl;
    yval[j+1] = yval[j] + slpA[j]*incA[j];
    // std::cout << "j: " << j << " yval[j-1] " << yval[j-1] 
    //           << " yval[j] " << yval[j]  
    //           << " slpA[j] " << slpA[j]
    //           << " incA[j] " << incA[j]
    //           << std::endl;
    yMinOrMax = isUnder?std::min(yMinOrMax,yval[j+1]):std::max(yMinOrMax,yval[j+1]);
    slpDiff1 = slpDiff2;
  }


  yval[n-2] = slpA[n-3]*incA[n-3]+yval[n-3];  
  yval[n-1] = slpA[n-2]*incA[n-2]+yval[n-2];

  // std::cout << "          input vertices" << std::endl;
  // for (unsigned int jj = 0; jj < yval.size(); jj++)
  //   std::cout <<"              " << jj << " : " << yval[jj] << "  - " << incA[jj] << std::endl; 
  // for (unsigned int jj = 0; jj < incA.size(); jj++)
  //   assert(incA[jj] > 0);
  
#ifdef MC__CPWLUSG_DEBUG    
          for (unsigned int jj = 1; jj < yval.size(); jj++)
            if(std::fabs(yval[jj]) > 1000){std::cout << jj << " in " << yval.size() << " "<< yval[jj] << std::endl; assert(std::fabs(yval[jj]) < 1000);};
#endif  

  //double lbnd;
  //double ubnd;

  if(isUnder){
    yMinOrMax = std::min(yMinOrMax,yval[n-2]);      
    yMinOrMax = std::min(yMinOrMax,yval[n-1]);    
    // if(!_lbnd.second) {_lbnd.second = true;_lbnd.first = yMinOrMax;};
    // _ubnd.second = false;
    //lbnd = yMinOrMax;
  }
  else{
    yMinOrMax = std::max(yMinOrMax,yval[n-2]);      
    yMinOrMax = std::max(yMinOrMax,yval[n-1]);
    // if(!_ubnd.second) {_ubnd.second = true;_ubnd.first = yMinOrMax;};
    // _lbnd.second = false;
    //ubnd = yMinOrMax;
  }
  
  // std::cout << "yMinOrMax " << yMinOrMax << std::endl;

  std::make_heap(areaLoss.begin(), areaLoss.end(), _cmp_indexed_areaLoss);

  std::vector<int> hind(n-3);       // mapping each window index to the position in the heap areaLoss
  for (unsigned int i = 0; i < n-3; i++) {
      hind[areaLoss[i].second] = (int) i;
  }        
    

  std::vector<int> nxtInd(n); // = [i for i in range(1,n+1)]
  std::vector<int> prvInd(n); // = [i for i in range(-1,n-1,1)]
  for (int i = 0; i < ((int) n); i++) {
      prvInd[i] = i - 1;  
      nxtInd[i] = i + 1;        
  }                 

  auto & aMin = areaLoss[0]; // bind the reference to the top of the heap, heappop(areaLoss,hind) in python 
    // std::cout << "aMin " << aMin.first << " at " << aMin.second << std::endl;
    // std::cout << "areaLoss len" << areaLoss.size() << std::endl;
    // std::cout << "aMin1 " << areaLoss[1].first << " at " << areaLoss[1].second << std::endl;
  //unsigned int ctr = areaLoss.size();

  // if (isUnder){
  //   while (yval[aMin.second] + slpXL[aMin.second]*deltaXL[aMin.second] + 1e-14 < yMinOrMax){
  //     areaLoss[0].first = DBL_MAX;
  //     _modify_heap(areaLoss,0,hind);
  //     ctr --;
  //     if(ctr == 0 ){
  //       // _lbnd.first = yval[aMin.second] + slpXL[aMin.second]*deltaXL[aMin.second]; 
  //       // yMinOrMax = _lbnd.first; 
  //       yMinOrMax = yval[aMin.second] + slpXL[aMin.second]*deltaXL[aMin.second]; 
  //       break;
  //     }      
  //     // std::cout << "aMin " << aMin.first << " at " << aMin.second << std::endl;
  //     // if(areaLoss.size() <= 2) {
  //     //   std::cout << " here we only have 1-2 choices" << std::endl;
  //     //   throw typename PWLU::Exceptions( PWLU::Exceptions::UNDEF );
  //     // }        
  //     // aMin = areaLoss[0]; // as the reference binds to that, we do not need assignment
  //   }
  // }
  // else{
  //   while (yval[aMin.second] + slpXL[aMin.second]*deltaXL[aMin.second] + 1e-14 > yMinOrMax){
  //     areaLoss[0].first = DBL_MAX;
  //     _modify_heap(areaLoss,0,hind);
  //     ctr --;
  //     if(ctr == 0 ){
  //       // _ubnd.first = yval[aMin.second] + slpXL[aMin.second]*deltaXL[aMin.second]; 
  //       // yMinOrMax = _ubnd.first;
  //       yMinOrMax = yval[aMin.second] + slpXL[aMin.second]*deltaXL[aMin.second]; 
  //       break;
  //     };
  //     // aMin = areaLoss[0]; // as the reference binds to that, we do not need assignment
  //   }       
  // }    

  //std::cout <<  "            aMin.second " << aMin.second << " with " << std::scientific << std::setprecision(14) << aMin.first << std::endl;    
  // updating the indexes so that we can skip the second element in the window, and 
  //   update the third one with the computed one
  // std::cout << "update the third one with the computed one" << std::endl;
  int pos = aMin.second;
  int id_1 = nxtInd[aMin.second];
  int id   = nxtInd[id_1];
  prvInd[id]   = aMin.second;
  nxtInd[aMin.second] = id;
 
  // std::cout << "        pos = "<< pos << std::endl;

  incA[pos] = deltaXL[aMin.second];
  slpA[pos] = slpXL[aMin.second];
  incA[id]  = deltaXR[aMin.second];
  slpA[id]  = slpXR[aMin.second];
  yval[id] = yval[pos] + slpXL[aMin.second]*deltaXL[aMin.second];

#ifdef MC__CPWLUSG_DEBUG    
    assert(incA[pos] >= 1e-15);
    assert(incA[id] >= 1e-15);    
    assert(std::fabs(slpA[pos]) <= 50);
    assert(std::fabs(slpA[id]) <= 50);  
#endif

  //T accAreaLoss = aMin.first;
  // std::cout << "n-k-1 " << n-k-1 << std::endl;    
  for(unsigned int i = 1; i < n-k-1; i ++ ){
    // std::cout << "i =  " << i << std::endl; 
    _heap_pop(areaLoss,hind);
    const int idx1 = pos;
    //const int idx2 = id; 
    const int idx3 = nxtInd[id]; 

    if (pos > 0){
      const int idx0 = prvInd[idx1];

      const std::vector<double> incIn{incA[idx0],incA[pos],incA[id]}; 
      const std::vector<double> slpIn{slpA[idx0],slpA[pos],slpA[id]};

      slpDiff1 = isUnder?(slpIn[0]-slpIn[1]):(slpIn[1]-slpIn[0]);       
      slpDiff2 = isUnder?(slpIn[1]-slpIn[2]):(slpIn[2]-slpIn[1]); 
      _underestimate_in_window3(incIn,slpIn,0,slpDiff1,slpDiff2,data,isUnder);
      areaLoss[hind[idx0]].first = data[4]; // aLoss;
      deltaXL[idx0] = data[0];
      deltaXR[idx0] = data[1];
      slpXL[idx0] = data[2];
      slpXR[idx0] = data[3];    
      _modify_heap(areaLoss,hind[idx0],hind);
      // std::cout << "data[0] " << data[0] << std::endl;
      // std::cout << "data[1] " << data[1] << std::endl;     
      // std::cout << "data[2] " << data[2] << std::endl;            
      if (idx0 > 0){

        const int idx_1 = prvInd[idx0];
        const std::vector<double> incIn{incA[idx_1],incA[idx0],incA[pos]}; 
        const std::vector<double> slpIn{slpA[idx_1],slpA[idx0],slpA[pos]};                

        slpDiff2 = slpDiff1;
        slpDiff1 = isUnder?(slpIn[0]-slpIn[1]):(slpIn[1]-slpIn[0]);
        _underestimate_in_window3(incIn,slpIn,0,slpDiff1,slpDiff2,data,isUnder);          
        areaLoss[hind[idx_1]].first = data[4];
        deltaXL[idx_1] = data[0];
        deltaXR[idx_1] = data[1];
        slpXL[idx_1] = data[2];
        slpXR[idx_1] = data[3];  
        _modify_heap(areaLoss,hind[idx_1],hind);
        // std::cout << "data[0] " << data[0] << std::endl;
        // std::cout << "data[1] " << data[1] << std::endl;     
        // std::cout << "data[2] " << data[2] << std::endl;                
      }
    } 
    if (idx3 <= ((int )n - 2)){  
    //  the index of the right endpoint is n-1

      const int idx4 = nxtInd[idx3];
      const std::vector<double> incIn{incA[pos],incA[id],incA[idx3]}; 
      const std::vector<double> slpIn{slpA[pos],slpA[id],slpA[idx3]};     
      slpDiff1 = isUnder?(slpIn[0]-slpIn[1]):(slpIn[1]-slpIn[0]);
      slpDiff2 = isUnder?(slpIn[1]-slpIn[2]):(slpIn[2]-slpIn[1]);
      _underestimate_in_window3(incIn,slpIn,0,slpDiff1,slpDiff2,data,isUnder);                          
      areaLoss[hind[id_1]].first = data[4];
      areaLoss[hind[id_1]].second = pos;
      deltaXL[pos] = data[0];
      deltaXR[pos] = data[1];
      slpXL[pos] = data[2];
      slpXR[pos] = data[3];  
      //hind[idx1] = hind[id_1]
      _modify_heap(areaLoss,hind[id_1],hind);

      // std::cout << "data[0] " << data[0] << std::endl;
      // std::cout << "data[1] " << data[1] << std::endl;     
      // std::cout << "data[2] " << data[2] << std::endl;    
      if (idx4 <= ((int )n - 2)){


        const std::vector<double> incIn{incA[id],incA[idx3],incA[idx4]}; 
        const std::vector<double> slpIn{slpA[id],slpA[idx3],slpA[idx4]};                     
        slpDiff1 = slpDiff2;
        slpDiff2 = isUnder?(slpIn[1]-slpIn[2]):(slpIn[2]-slpIn[1]);
        _underestimate_in_window3(incIn,slpIn,0,slpDiff1,slpDiff2,data,isUnder);    

        areaLoss[hind[id]].first = data[4]; // aLoss
        deltaXL[id] = data[0];
        deltaXR[id] = data[1];
        slpXL[id] = data[2];
        slpXR[id] = data[3]; 
        _modify_heap(areaLoss,hind[id],hind);
        // std::cout << "data[0] " << data[0] << std::endl;
        // std::cout << "data[1] " << data[1] << std::endl;     
        // std::cout << "data[2] " << data[2] << std::endl;                      
      }
    } 

    // std::cout << "aMin " << areaLoss[0].first << " at " << areaLoss[0].second << std::endl;
    // std::cout << "aMin1 " << areaLoss[1].first << " at " << areaLoss[1].second << std::endl;


    // // if(aMin.first != areaLoss[0].first) std::cout << "error" << std::endl; // test to make sure the reference binds to the heap top
    // // aMin = areaLoss[0]; // as the reference binds to that, we do not need assignment, heappop(areaLoss,hind) in python
    // unsigned int ctr = areaLoss.size();
    // if (isUnder){
    //   while (yval[aMin.second] + slpXL[aMin.second]*deltaXL[aMin.second] + 1e-14 < yMinOrMax){
    //     areaLoss[0].first = DBL_MAX;
    //     _modify_heap(areaLoss,0,hind);
    //     ctr --;
    //     if(ctr == 0 ){
    //       // _lbnd.first = yval[aMin.second] + slpXL[aMin.second]*deltaXL[aMin.second]; 
    //       // yMinOrMax = _lbnd.first; 
    //       yMinOrMax = yval[aMin.second] + slpXL[aMin.second]*deltaXL[aMin.second]; 
    //       break;
    //     }      
    //     // std::cout << "aMin " << aMin.first << " at " << aMin.second << std::endl;
    //     // if(areaLoss.size() <= 2) {
    //     //   std::cout << " here we only have 1-2 choices" << std::endl;
    //     //   throw typename PWLU::Exceptions( PWLU::Exceptions::UNDEF );
    //     //}        
    //     //aMin = areaLoss[0]; // as the reference binds to that, we do not need assignment
    //   }
    // }
    // else{
    //   while (yval[aMin.second] + slpXL[aMin.second]*deltaXL[aMin.second] + 1e-14 > yMinOrMax){
    //     areaLoss[0].first = DBL_MAX;
    //     _modify_heap(areaLoss,0,hind);
    //     ctr --;
    //     if(ctr == 0 ){
    //       // _ubnd.first = yval[aMin.second] + slpXL[aMin.second]*deltaXL[aMin.second]; 
    //       // yMinOrMax = _ubnd.first;
    //       yMinOrMax = yval[aMin.second] + slpXL[aMin.second]*deltaXL[aMin.second];           
    //       break;
    //     };
    //     //aMin = areaLoss[0]; // as the reference binds to that, we do not need assignment
    //   }       
    // }    
    // //std::cout <<  "              aMin.second " << aMin.second << " with " << std::scientific << std::setprecision(14) << aMin.first << std::endl;

    pos = aMin.second;
    id_1 = nxtInd[aMin.second];
    id   = nxtInd[id_1];
    prvInd[id]   = aMin.second;
    nxtInd[aMin.second] = id;

    incA[pos] = deltaXL[aMin.second];
    slpA[pos] = slpXL[aMin.second];
    incA[id]  = deltaXR[aMin.second];
    slpA[id]  = slpXR[aMin.second];
    yval[id] = yval[pos] + slpXL[aMin.second]*deltaXL[aMin.second];

#ifdef MC__CPWLUSG_DEBUG    
    assert(incA[pos+1] >= 1e-15);
    assert(incA[id+1] >= 1e-15);    

    assert(std::fabs(slpA[pos+1]) <= 50);
    assert(std::fabs(slpA[id+1]) <= 50);  
#endif    
    //accAreaLoss = accAreaLoss + aMin.first;
  }



  // std::vector<double> incC(k); 
  // std::vector<double> slpC(k);  

  unsigned int i = 0;
  for(unsigned int j = 0; j < k; j++){
    incA[j] = incA[i];
    slpA[j] = slpA[i];
#ifdef MC__CPWLUSG_DEBUG        
    assert(i+1 < first.size()); 
#endif     
    // if(j==k-1){
    //   std::cout << i+1 << " , " << first.size() << std::endl;
    // }
    // std::cout << "    i = " << i << std::endl;
    i = nxtInd[i];
  }
    
  // std::cout << " , " << incA.size() << std::endl;  
  // std::cout << " , " << _dy.size() << std::endl;
  _dx.resize(k);
  _dy.resize(k);
  // _dx.swap(incC);
  // _dy.swap(slpC);


  // std::cout << "          output vertices" << std::endl;
  // std::cout <<"              0" << " : " << yval[0] << std::endl; 
  // for (unsigned int jj = 0; jj < _dx.size(); jj++){
  //   yval[jj+1] = yval[jj] + _dx[jj]*_dy[jj];
  //   std::cout <<"              " << jj+1 << " : " << yval[jj+1] << "  - " << _dx[jj] << std::endl; 
  // }

  // for (unsigned int jj = 0; jj < _dx.size(); jj++)
  //   assert(_dx[jj] > 0);

#ifdef MC__CPWLUSG_DEBUG    
          for (unsigned int jj = 1; jj < first.size(); jj++)
            if(first[jj] < 1e-15){std::cout  << jj << " in " << first.size() << "\n"   << first[jj] << std::endl; assert(first[jj] >= 1e-15);};
#endif

  // first.back() = xRightEndpoint;

#ifdef MC__CPWLUSG_DEBUG              
  for (unsigned int jj = 1; jj < second.size(); jj++)
    if(std::fabs(second[jj]) > 50) {
      std::cout << jj << " in " << second.size() << std::endl;
      std::cout << second[jj] << std::endl;
      assert(std::fabs(second[jj]) <= 50);
    }  
#endif
  //std::cout <<  "leave segment reduction" << std::endl;

  return *this;
}  

inline
bool PWLU::_cmp_indexed_areaLoss
( const std::pair<double, int>& a, const std::pair<double, int>& b )
{
  return a.first > b.first;
}

inline
std::pair<double, int> PWLU::_heap_pop
( std::vector<std::pair<double, int>>& heap, std::vector<int>& hind )
{
  // Pop the smallest item off the heap, maintaining the heap invariant."""
  auto lastelt =  heap[heap.size()-1];
  heap.pop_back();    // raises appropriate IndexError if heap is empty
  if (!heap.empty()){
    auto returnitem = heap[0];
    hind[heap[0].second] = -1;
    heap[0] = lastelt;
    hind[lastelt.second] = 0;
    _heap_bubble_down(heap, 0, hind);
    return returnitem;
  }
  hind[lastelt.second] = -1;
  return lastelt;
}

inline 
void PWLU::_modify_heap
( std::vector<std::pair<double, int>>& heap, const int pos, std::vector<int>& hind )
{
  if (pos == 0){
    _heap_bubble_down(heap, pos, hind);
    //print("pos",pos,"heap[pos].al",heap[pos].al)
  }
  else{
    unsigned int parentpos = (pos - 1) >> 1;
    // print("pos",pos,"parentpos",parentpos)
    // print("heap[pos]",heap[pos].al,"heap[parentpos]",heap[parentpos].al)
    if (heap[pos].first < heap[parentpos].first){
      std::swap(heap[pos],heap[parentpos]);
      // print("heap[pos]",heap[pos],"heap[parentpos]",heap[parentpos])
      hind[heap[pos].second] = pos;
      hind[heap[parentpos].second] = parentpos;
      _heap_bubble_up(heap, parentpos, hind);
    }
    else
      _heap_bubble_down(heap, pos, hind);
  }

}

inline 
void PWLU::_heap_bubble_up
(std::vector<std::pair<double, int>> & heap, int pos, std::vector<int> & hind)
{
    auto newitem = heap[pos];
    // Follow the path to the root, moving parents down until finding a place
    // newitem fits.
    while (pos > 0){
        unsigned int parentpos = (pos - 1) >> 1;
        auto parent = heap[parentpos];
        if (newitem.first < parent.first){
            heap[pos] = parent;
            hind[parent.second] = pos;
            pos = parentpos;
            continue;
        }
        break;
    }
    heap[pos] = newitem;
    hind[newitem.second] = pos;
}

inline
void PWLU::_heap_bubble_down
( std::vector<std::pair<double, int>>& heap, int pos, std::vector<int>& hind )
{
    // print("bubble_down")
    int endpos = heap.size();
    //int startpos = pos;
    auto newitem = heap[pos];
    // print("endpos ", endpos,"startpos ",startpos,"newitem ",newitem)
    //  Bubble up the smaller child until hitting a leaf.
    int childpos = 2*pos + 1;    // leftmost child position
    while (childpos < endpos){
        // Set childpos to index of smaller child.
        int rightpos = childpos + 1;
        // print("left: heap[childpos] ", heap[childpos], "heap[rightpos]", heap[rightpos])
        if (rightpos < endpos && (heap[childpos].first >= heap[rightpos].first)){
            childpos = rightpos;
        // Move the smaller child up.
        }
        if (newitem.first > heap[childpos].first){
            heap[pos] = heap[childpos];
            hind[heap[pos].second] = pos;
            pos = childpos;
            childpos = 2*pos + 1;
        }
        else break;
    }
    // The leaf at pos is empty now.  Put newitem there, and bubble it up
    // to its final resting place (by sifting its parents down).
    heap[pos] = newitem;
    hind[newitem.second] = pos;
}

inline
std::ostream& operator<<
( std::ostream& os, PWLU const& var )
{
  var.display( os );
  return os;
}

inline
PWLU operator+
( PWLU const& var )
{
  return var;
}

inline
PWLU operator+
( PWLU && var )
{
  PWLU res( std::move(var) );  
  return var;
}

inline
PWLU operator+
( PWLU const& var1, PWLU const& var2 )
{
  if( var1.dx().size() < var2.dx().size() ){
    PWLU res( var2 );
    res += var1;
    return res;
  }

  PWLU res( var1 );
  res += var2;
  return res;
}
/*
inline
PWLU operator+
( PWLU const& var1, PWLU && var2 )
{
  PWLU res( std::move(var2) );
  res += var1;
  return res;
}
*/
inline
PWLU operator+
( PWLU && var1, PWLU const& var2 )
{
  PWLU res( std::move(var1) );
  res += var2;
  return res;
}

inline
PWLU operator+
( PWLU const& var1, double const& cst2 )
{
  PWLU res( var1 );
  res += cst2;
  return res;
}

inline
PWLU operator+
( PWLU && var1, double const& cst2 )
{
  PWLU res( std::move(var1) );
  res += cst2;
  return res;
}

inline
PWLU operator+
( double const& cst1, PWLU const& var2 )
{
  PWLU res( var2 );
  res += cst1;
  return res;
}

inline
PWLU operator+
( double const& cst1, PWLU && var2 )
{
  PWLU res( std::move(var2) );
  res += cst1;
  return res;
}

inline
PWLU operator-
( PWLU const& var )
{
  PWLU res( var );
  return( res *= -1 );
}

inline
PWLU operator-
( PWLU && var )
{
  PWLU res( std::move(var) );
  return( res *= -1 );
}

inline
PWLU operator-
( PWLU const& var1, PWLU const& var2 )
{
  PWLU res( var1 );
  res -= var2;
  return res;
}
/*
inline
PWLU operator-
( PWLU const& var1, PWLU && var2 )
{
  PWLU res( std::move(-var2) );
  res += var1;  
  return res;
}
*/
inline
PWLU operator-
( PWLU && var1, PWLU const& var2 )
{
  PWLU res( std::move(var1) );
  res -= var2;
  return res;
}

inline
PWLU operator-
( PWLU const& var1, double const& cst2 )
{
  PWLU res( var1 );
  res += -cst2;
  return res;
}

inline
PWLU operator-
( PWLU && var1, double const& cst2 )
{
  PWLU res( std::move(var1) );
  res += -cst2;
  return res;
}

inline
PWLU operator-
( double const& cst1, PWLU const& var2 )
{
  PWLU res( -var2 );
  res += cst1;
  return res;
}

inline
PWLU operator-
( double const& cst1, PWLU && var2 )
{
  PWLU res( -std::move(var2) );
  res += cst1;
  return res;
}

inline
PWLU operator*
( PWLU const& var1, double const& cst2 )
{
  PWLU res( var1 );
  res *= cst2;
  return res;
}

inline
PWLU operator*
( PWLU && var1, double const& cst2 )
{
  PWLU res( std::move(var1) );
  res *= cst2;
  return res;
}

inline
PWLU operator*
( double const& cst1, PWLU const& var2 )
{
  PWLU res( var2 );
  res *= cst1;
  return res;
}

inline
PWLU operator*
( double const& cst1, PWLU && var2 )
{
  PWLU res( std::move(var2) );
  res *= cst1;
  return res;
}

inline
PWLU operator/
( PWLU const& var1, double const& cst2 )
{
  PWLU var3( var1 );
  var3 /= cst2;
  return var3;
}

inline
PWLU operator/
( PWLU && var1, double const& cst2 )
{
  PWLU var3( std::move(var1) );
  var3 /= cst2;
  return var3;
}

} // namespace mc

#endif
