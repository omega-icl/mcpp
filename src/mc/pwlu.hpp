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

#include <iostream>
#include <iomanip> 
#include <vector> 
#include <algorithm>

#include "mcfunc.hpp"

//#define MC__PWLU_DEBUG
//#define MC__PWLU_TRACE
//#define MC__PWLU_CHECK
#undef  MC__PWLU_FULL_UPDATE

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
    Options():
      BKPTATOL( 1e2*DBL_EPSILON ),
      BKPTRTOL( 1e2*DBL_EPSILON ),
      LOSSRTOL( 1e-2 ),
      DISPNUM( 5 )
      {}
    //! @brief Assignment
    Options& operator=
      ( Options const& opt )
      {
        BKPTATOL        = opt.BKPTATOL;
        BKPTRTOL        = opt.BKPTRTOL;
        LOSSRTOL        = opt.LOSSRTOL;
        DISPNUM         = opt.DISPNUM;
        return *this;
      }
    //! @brief Assignment
    void reset
      ()
      {
        BKPTATOL = 1e2*DBL_EPSILON;
        BKPTRTOL = 1e2*DBL_EPSILON;
        LOSSRTOL = 1e-2;
        DISPNUM  = 5;
      }
    //! @brief Absolute tolerance in breakpoints - Default: 1e2*DBL_EPSILON
    double BKPTATOL;
    //! @brief Relative tolerance in breakpoints - Default: 1e2*DBL_EPSILON
    double BKPTRTOL;
    //! @brief Relative tolerance in reduction loss - Default: 1e-2
    double LOSSRTOL;
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
    ( bool const under, size_t const nseg )
    {
      //if( !nseg )
      //  return *this;
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
      static thread_local std::vector<std::tuple<double,int,int,double,double>> vloss;
      vloss.resize( _dx.size()-2 );

      static thread_local double x[5], y[5], lcv[3];
      double *x0 = x, *x1 = x0+1, *x2 = x1+1, *x3 = x2+1, *x4 = x3+1;
      double *y0 = y, *y1 = y0+1, *y2 = y1+1, *y3 = y2+1, *y4 = y3+1;
      double *lcv0 = lcv, *lcv2 = lcv0+1, *lcv3 = lcv2+1;
      
      auto idx = _dx.begin(), idy = _dy.begin();
      x[1] = _x0,             y[1] = _y0;
      x[2] = x[1] + *idx,     y[2] = y[1] + *idy * *idx,         lcv[1]  = *idy;
      x[3] = x[2] + *(++idx), y[3] = y[2] + *(++idy) * *idx,     lcv[1] -= *idy; lcv[2]  = *idy; 
      x[4] = x[3] + *(++idx), y[4] = y[3] + *(++idy) * *idx;                     lcv[2] -= *idy;

      double totalloss( 0. );
      for( auto ploss = vloss.data(); ; ploss++ ){
        _loss( under, *lcv2<0, *lcv3<0, *x1, *y1, *x2, *y2, *x3, *y3, *x4, *y4, *ploss );
#ifdef MC__PWLU_DEBUG
        std::cout << "loss: " << std::get<0>( *ploss ) << std::endl;
#endif
        totalloss += std::get<0>( *ploss );
        
        if( ++idx == _dx.end() ) break;
        lcv0 = lcv2, lcv2 = lcv3, lcv3 = lcv0;
        x0 = x1; x1 = x2; x2 = x3; x3 = x4; x4 = x0; *x4 = *x3 + *idx;            *lcv3  = *idy;
        y0 = y1; y1 = y2; y2 = y3; y3 = y4; y4 = y0; *y4 = *y3 + *(++idy) * *idx; *lcv3 -= *idy;
      }
      totalloss *= options.LOSSRTOL;
      
      // Main reduction loop
      for( ; ; ){

#ifdef MC__PWLU_DEBUG
        std::cout << "dx:";
        for( auto dxi : _dx ) std::cout << " " << dxi;
        std::cout << std::endl << "dy:";
        for( auto dyi : _dy ) std::cout << " " << dyi;
        std::cout << std::endl;

        std::cout << "Candidate 3-segments:\n";
        //for( auto itloss=loss.cbegin(); itloss!=loss.cend(); ++itloss )
        for( auto itloss=vloss.cbegin(); itloss!=vloss.cend(); ++itloss )
          std::cout << std::get<0>(**itloss) << " "
                    << std::get<1>(**itloss) << " "
                    << std::get<2>(**itloss) << " "
                    << std::get<3>(**itloss) << " "
                    << std::get<4>(**itloss) << " " << std::endl;
#endif

        // Determine minimal element
        auto ilossopt = std::min_element(
          vloss.cbegin(), vloss.cend(),
          [=]( std::tuple<double,int,int,double,double> const& a, std::tuple<double,int,int,double,double> const& b )
          { // Below threshold eliminate breakpoints that create smaller gaps
            auto const& fmax = [=]( const double& x1, const double& x2 ){ return x1>x2? x1: x2; };
            if( std::get<0>(a) < totalloss && std::get<0>(b) < totalloss )
              return fmax( std::get<3>(a),std::get<4>(a) ) < fmax( std::get<3>(b),std::get<4>(b) );
            // Above threshold eliminate breakpoints according to area loss
            return std::get<0>(a) < std::get<0>(b); }
        );

        auto idxopt = _dx.begin(), idyopt = _dy.begin();
        size_t pos = std::distance( vloss.cbegin(), ilossopt );
        if( pos ){
          std::advance( idxopt, pos );
          std::advance( idyopt, pos );
        }

        // Remove minimal element
        if( !std::get<2>(*ilossopt) ){
          std::advance( idxopt, std::get<1>(*ilossopt)-1 );
          std::advance( idyopt, std::get<1>(*ilossopt)-1 );
          double dx = *idxopt, dy = *idyopt;
          idxopt = _dx.erase( idxopt );
          idyopt = _dy.erase( idyopt );
          _average( *idxopt, *idyopt, dx, dy );
        }
        else{
          *idxopt = std::get<3>(*ilossopt);
          idxopt = _dx.erase( ++idxopt );
          idyopt = _dy.erase( ++idyopt );
          *idxopt = std::get<4>(*ilossopt);        
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

#ifdef MC__PWLU_FULL_UPDATE
        x0 = x, x1 = x0+1, x2 = x1+1, x3 = x2+1, x4 = x3+1;
        y0 = y, y1 = y0+1, y2 = y1+1, y3 = y2+1, y4 = y3+1;
        lcv0 = lcv, lcv2 = lcv0+1, lcv3 = lcv2+1;
      
        x[1] = _x0,             y[1] = _y0;
        x[2] = x[1] + *idx,     y[2] = y[1] + *idy * *idx,         lcv[1]  = *idy;
        x[3] = x[2] + *(++idx), y[3] = y[2] + *(++idy) * *idx,     lcv[1] -= *idy; lcv[2]  = *idy; 
        x[4] = x[3] + *(++idx), y[4] = y[3] + *(++idy) * *idx;                     lcv[2] -= *idy;

        for( ; itloss!=vloss.end(); ++itloss ){
          _loss( under, *lcv2<0, *lcv3<0, *x1, *y1, *x2, *y2, *x3, *y3, *x4, *y4, *itloss );
#ifdef MC__PWLU_DEBUG
          std::cout << std::get<0>( *itloss ) << std::endl;
#endif
          if( ++idx == _dx.end() ) break;
          lcv0 = lcv2, lcv2 = lcv3, lcv3 = lcv0;
          x0 = x1; x1 = x2; x2 = x3; x3 = x4; x4 = x0; *x4 = *x3 + *idx;            *lcv3  = *idy;
          y0 = y1; y1 = y2; y2 = y3; y3 = y4; y4 = y0; *y4 = *y3 + *(++idy) * *idx; *lcv3 -= *idy;
        }

#else
        size_t p = 0;
        x[1] = _x0,             y[1] = _y0;
        for( ; p+2<pos; ++p, ++idx, ++idy, ++itloss )
          x[1] += *idx, y[1] += *idy * *idx;
        x[2] = x[1] + *idx,     y[2] = y[1] + *idy * *idx,         lcv[1]  = *idy;
        x[3] = x[2] + *(++idx), y[3] = y[2] + *(++idy) * *idx,     lcv[1] -= *idy; lcv[2]  = *idy; 
        x[4] = x[3] + *(++idx), y[4] = y[3] + *(++idy) * *idx;                     lcv[2] -= *idy;
        _loss( under, lcv[1]<0, lcv[2]<0, x[1], y[1], x[2], y[2], x[3], y[3], x[4], y[4], *itloss );      
#ifdef MC__PWLU_DEBUG
        std::cout << std::get<0>( *itloss ) << std::endl;
#endif
        if( ++idx==_dx.end() || ++p>pos+2 ) continue;

        x[0] = x[4] + *idx;            lcv[0]  = *idy;
        y[0] = y[4] + *(++idy) * *idx; lcv[0] -= *idy;
        _loss( under, lcv[2]<0, lcv[0]<0, x[2], y[2], x[3], y[3], x[4], y[4], x[0], y[0], *(++itloss) );
#ifdef MC__PWLU_DEBUG
        std::cout << std::get<0>( *itloss ) << std::endl;
#endif
        if( ++idx==_dx.end() || ++p>pos+2 ) continue;

        x[1] = x[0] + *idx;            lcv[1]  = *idy;
        y[1] = y[0] + *(++idy) * *idx; lcv[1] -= *idy;
        _loss( under, lcv[0]<0, lcv[1]<0, x[3], y[3], x[4], y[4], x[0], y[0], x[1], y[1], *(++itloss) );
#ifdef MC__PWLU_DEBUG
        std::cout << std::get<0>( *itloss ) << std::endl;
#endif
        if( ++idx==_dx.end() || ++p>pos+2 ) continue;

        x[2] = x[1] + *idx;            lcv[2]  = *idy;
        y[2] = y[1] + *(++idy) * *idx; lcv[2] -= *idy;
        _loss( under, lcv[1]<0, lcv[2]<0, x[4], y[4], x[0], y[0], x[1], y[1], x[2], y[2], *(++itloss) );
#ifdef MC__PWLU_DEBUG
        std::cout << std::get<0>( *itloss ) << std::endl;
#endif
        if( ++idx==_dx.end() || ++p>pos+2 ) continue;

        x[3] = x[2] + *idx;            lcv[0]  = *idy;
        y[3] = y[2] + *(++idy) * *idx; lcv[0] -= *idy;
        _loss( under, lcv[2]<0, lcv[0]<0, x[0], y[0], x[1], y[1], x[2], y[2], x[3], y[3], *(++itloss) );
#ifdef MC__PWLU_DEBUG
        std::cout << std::get<0>( *itloss ) << std::endl;
#endif
        if( ++idx==_dx.end() || ++p>pos+2 ) continue;
#endif
      }

      return *this;
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
    ( bool const under, bool const lcv2, bool const lcv3,
      double const& x1, double const& y1, double const& x2, double const& y2,
      double const& x3, double const& y3, double const& x4, double const& y4,
      std::tuple<double,int,int,double,double>& loss )
    {
      if( ( under && lcv2 && lcv3 ) || ( !under && !lcv2 && !lcv3 ) ){
#ifdef MC__PWLU_DEBUG
        std::cout << "remove (" << x2 << "," << y2 << ")" << std::endl;
        std::cout << "remove (" << x3 << "," << y3 << ")" << std::endl;
#endif
        // cross point of (x1,y1)-(x2,y2) and (x3,y3)-(x4,y4) lines
        auto const& [x5,y5] = _intersect( x1,y1, x2,y2, x3,y3, x4,y4 );
#ifdef MC__PWLU_DEBUG
        std::cout << "insert (" << x5 << "," << y5 << ")" << std::endl;
#endif
        loss = std::make_tuple( _area( x2,y2, x3,y3, x5,y5 ), 1, 2, x5-x1, x4-x5 );
        return;
      }

      else if( under && lcv2 && !lcv3 ){
        if( (y2-y1)*(x4-x1) <= (y4-y1)*(x2-x1) + options.BKPTATOL ){
#ifdef MC__PWLU_DEBUG
          std::cout << "remove (" << x3 << "," << y3 << ")" << std::endl;
#endif
          // keep (x2,y2)
          loss = std::make_tuple( _area( x2,y2, x3,y3, x4,y4 ), 2, 0, x2-x1, x4-x2 );
        }
        else{
#ifdef MC__PWLU_DEBUG
          std::cout << "remove (" << x2 << "," << y2 << ")" << std::endl;
          std::cout << "remove (" << x3 << "," << y3 << ")" << std::endl;
#endif
          // cross point of (x1,y1)-(x2,y2) and (x3,y3)-(x4,y4) lines
          auto const& [x5,y5] = _intersect( x1,y1, x2,y2, x3,y3, x4,y4 );
#ifdef MC__PWLU_DEBUG
          std::cout << "insert (" << x5 << "," << y5 << ")" << std::endl;
#endif
          loss = std::make_tuple( _area( x2,y2, x3,y3, x5,y5 ), 1, 2, x5-x1, x4-x5 );
        }
        return;
      }

      else if( !under && !lcv2 && lcv3 ){
        if( (y2-y1)*(x4-x1) + options.BKPTATOL >= (y4-y1)*(x2-x1) ){
#ifdef MC__PWLU_DEBUG
          std::cout << "remove (" << x3 << "," << y3 << ")" << std::endl;
#endif
          // keep (x2,y2)
          loss = std::make_tuple( _area( x2,y2, x3,y3, x4,y4 ), 2, 0, x2-x1, x4-x2 );
        }
        else{
#ifdef MC__PWLU_DEBUG
          std::cout << "remove (" << x2 << "," << y2 << ")" << std::endl;
          std::cout << "remove (" << x3 << "," << y3 << ")" << std::endl;
#endif
          // cross point of (x1,y1)-(x2,y2) and (x3,y3)-(x4,y4) lines
          auto const& [x5,y5] = _intersect( x1,y1, x2,y2, x3,y3, x4,y4 );
#ifdef MC__PWLU_DEBUG
          std::cout << "insert (" << x5 << "," << y5 << ")" << std::endl;
#endif
          loss = std::make_tuple( _area( x2,y2, x3,y3, x5,y5 ), 1, 2, x5-x1, x4-x5 );
        }
        return;
      }

      else if( under && !lcv2 && lcv3 ){
        if( (y4-y3)*(x4-x1) + options.BKPTATOL >= (y4-y1)*(x4-x3) ){
#ifdef MC__PWLU_DEBUG
          std::cout << "remove (" << x2 << "," << y2 << ")" << std::endl;
#endif
          // keep (x3,y3)
          loss = std::make_tuple( _area( x1,y1, x2,y2, x3,y3 ), 1, 0, x3-x1, x4-x3 );
        }
        else{
#ifdef MC__PWLU_DEBUG
          std::cout << "remove (" << x2 << "," << y2 << ")" << std::endl;
          std::cout << "remove (" << x3 << "," << y3 << ")" << std::endl;
#endif
          // cross point of (x1,y1)-(x2,y2) and (x3,y3)-(x4,y4) lines
          auto const& [x5,y5] = _intersect( x1,y1, x2,y2, x3,y3, x4,y4 );
#ifdef MC__PWLU_DEBUG
          std::cout << "insert (" << x5 << "," << y5 << ")" << std::endl;
#endif
          loss = std::make_tuple( _area( x2,y2, x3,y3, x5,y5 ), 1, 2, x5-x1, x4-x5 );
        }
        return;
      }

      else if( !under && lcv2 && !lcv3 ){
        if( (y4-y3)*(x4-x1) <= (y4-y1)*(x4-x3) + options.BKPTATOL ){
#ifdef MC__PWLU_DEBUG
          std::cout << "remove (" << x2 << "," << y2 << ")" << std::endl;
#endif
          // keep (x3,y3)
          loss = std::make_tuple( _area( x1,y1, x2,y2, x3,y3 ), 1, 0, x3-x1, x4-x3 );
        }
        else{
#ifdef MC__PWLU_DEBUG
          std::cout << "remove (" << x2 << "," << y2 << ")" << std::endl;
          std::cout << "remove (" << x3 << "," << y3 << ")" << std::endl;
#endif
          // cross point of (x1,y1)-(x2,y2) and (x3,y3)-(x4,y4) lines
          auto const& [x5,y5] = _intersect( x1,y1, x2,y2, x3,y3, x4,y4 );
#ifdef MC__PWLU_DEBUG
          std::cout << "insert (" << x5 << "," << y5 << ")" << std::endl;
#endif
          loss = std::make_tuple( _area( x2,y2, x3,y3, x5,y5 ), 1, 2, x5-x1, x4-x5 );
        }
        return;
      }

      // keep either (x2,y2) or (x3,y3)
      double const& A1 = _area( x1,y1, x2,y2, x3,y3 );
      double const& A2 = _area( x2,y2, x3,y3, x4,y4 );
      if( A1 < A2 ){
#ifdef MC__PWLU_DEBUG
        std::cout << "remove (" << x2 << "," << y2 << ")" << std::endl;
#endif
        loss = std::make_tuple( A1, 1, 0, x3-x1, x4-x3 );
      }
      else{
#ifdef MC__PWLU_DEBUG
        std::cout << "remove (" << x3 << "," << y3 << ")" << std::endl;
#endif
        loss = std::make_tuple( A2, 2, 0, x2-x1, x4-x2 );
      }
      return;
    }

  static double _area
    ( double const& x1, double const& y1, double const& x2, double const& y2,
      double const& x3, double const& y3 )
    {
      return std::fabs((x1-x3)*(y2-y1)-(x1-x2)*(y3-y1));
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
};

////////////////////////////////////////////////////////////////////////

inline PWLU::Options PWLU::options;

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
