// Copyright (C) Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the Eclipse Public License.

/*!
\page page_PWCU Piecewise constant and linear univariate estimators on fixed
partition \author Yanlin Zha, Beno&icirc;t Chachuat

The class mc::PWCU provides an implementation of piecewise constant univariate
estimators for use within a superposition relaxation (see \ref page_SUPREL). It
support the definition of such estimators, their propagation through linear
operations and compositions with convex/concave monotonic terms, and their
bounding and evaluation. It also (optionally) propagates slopes to enable
piecewise linear estimators that are quadratically convergent.
*/

#ifndef MC__PWCU_HPP
#define MC__PWCU_HPP

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <tuple>
#include <vector>

#include "mcfunc.hpp"

#undef MC__PWCU_DEBUG
#undef MC__PWCU_TRACE
#undef MC__PWCU_CHECK

namespace mc
{
//! @brief C++ class for propagation of piecewise constant univariate
//! estimators, supplemented with slopes
////////////////////////////////////////////////////////////////////////
//! mc::PWCU is a C++ class for propagation of piecewise constant
//! univariate estimators, supplemented with slopes, through factorable
//! expressions.
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

  //! @brief vector of y reference points
  std::vector<double> _y0;

  //! @brief vector of lower slopes at start
  std::vector<double> _sLL;

  //! @brief vector of lower slopes at end
  std::vector<double> _sLU;

  //! @brief vector of upper slopes at start
  std::vector<double> _sUL;

  //! @brief vector of upper slopes at end
  std::vector<double> _sUU;

  double
  _tighten(bool const under, size_t const i, size_t const N) const
  {
#ifdef MC__PWCU_CHECK
    if (i >= N || !options.SLOPEUSE) throw Exceptions(Exceptions::INTERNAL);
#endif

    double const dx = (_xU - _xL) / N;
    if (under)
    {
      if (_y0[i] <= _y0[i + 1])
      {
        if (_sLL[i] >= 0)
          return _y0[i];
        else if (isequal(_sLL[i], _sLU[i], options.BKPTATOL, options.BKPTRTOL))
          return std::max(_y0[i], _yL[i]);
        double const yint = _y0[i] + _sLL[i] / (_sLL[i] - _sLU[i]) *
                                         (_y0[i + 1] - _y0[i] - _sLU[i] * dx);
#ifdef MC__PWCU_CHECK
        if (yint >
            _y0[i] + options.BKPTATOL + std::fabs(_y0[i]) * options.BKPTRTOL)
        {
          std::cout << "yint = " << yint << ">? _y0[" << i << "] = " << _y0[i]
                    << std::endl;
          throw Exceptions(Exceptions::INTERNAL);
        }
#endif
        return std::max(yint, _yL[i]);  // intersect with constant bound
      }
      /*else if( _y0[i] > _y0[i+1] ){*/
      if (_sLU[i] <= 0)
        return _y0[i + 1];
      else if (isequal(_sLL[i], _sLU[i], options.BKPTATOL, options.BKPTRTOL))
        return std::max(_y0[i + 1], _yL[i]);
      double const yint = _y0[i] + _sLL[i] / (_sLL[i] - _sLU[i]) *
                                       (_y0[i + 1] - _y0[i] - _sLU[i] * dx);
#ifdef MC__PWCU_CHECK
      if (yint > _y0[i + 1] + options.BKPTATOL +
                     std::fabs(_y0[i + 1]) * options.BKPTRTOL)
      {
        std::cout << "yint = " << yint << ">? _y0[" << i + 1
                  << "] = " << _y0[i + 1] << std::endl;
        throw Exceptions(Exceptions::INTERNAL);
      }
#endif
      return std::max(yint, _yL[i]);  // intersect with constant bound
      /*}*/
    }

    /*else if( !under ){*/
    if (_y0[i] <= _y0[i + 1])
    {
      if (_sUU[i] >= 0)
        return _y0[i + 1];
      else if (isequal(_sUL[i], _sUU[i], options.BKPTATOL, options.BKPTRTOL))
        return std::min(_y0[i + 1], _yU[i]);
      double const yint = _y0[i] + _sUL[i] / (_sUL[i] - _sUU[i]) *
                                       (_y0[i + 1] - _y0[i] - _sUU[i] * dx);
#ifdef MC__PWCU_CHECK
      if (yint < _y0[i + 1] - options.BKPTATOL -
                     std::fabs(_y0[i + 1]) * options.BKPTRTOL)
      {
        std::cout << "yint = " << yint << "<? _y0[" << i + 1
                  << "] = " << _y0[i + 1] << std::endl;
        throw Exceptions(Exceptions::INTERNAL);
      }
#endif
      return std::min(yint, _yU[i]);  // intersect with constant bound
    }
    /*else if( _y0[i] > _y0[i+1] ){*/
    if (_sUL[i] <= 0)
      return _y0[i];
    else if (isequal(_sUL[i], _sUU[i], options.BKPTATOL, options.BKPTRTOL))
      return std::min(_y0[i], _yU[i]);
    double const yint = _y0[i] + _sUL[i] / (_sUL[i] - _sUU[i]) *
                                     (_y0[i + 1] - _y0[i] - _sUU[i] * dx);
#ifdef MC__PWCU_CHECK
    if (yint < _y0[i] - options.BKPTATOL - std::fabs(_y0[i]) * options.BKPTRTOL)
    {
      std::cout << "yint = " << yint << "<? _y0[" << i << "] = " << _y0[i]
                << std::endl;
      throw Exceptions(Exceptions::INTERNAL);
    }
#endif
    return std::min(yint, _yU[i]);  // intersect with constant bound
                                    /*}*/
    /*}*/
  }

  // cvx=0: concave  cvx=1: convex
  template <typename UNIV, typename DUNIV>
  PWCU&
  _compose(UNIV const& f, DUNIV const& df, int const cvx, bool const inc,
           double const& xmid = 0.0)
  {
    if (_yL.empty() || _yL.size() != _yU.size()) return *this;

    size_t const N = _yL.size();

    // Bounding valid regardless of convexity
    if (options.SLOPEUSE)
    {
#ifdef MC__PWCU_CHECK
      if (_y0.size() != N + 1 || _sLL.size() != N || _sLU.size() != N ||
          _sUL.size() != N || _sUU.size() != N)
        throw Exceptions(Exceptions::SLOPE);
#endif

      for (unsigned i = 0; i < N; ++i)
      {
        double sL = df(_yL[i]);
        double sU = df(_yU[i]);
        if (cvx > 1 && xmid > _yL[i] && xmid < _yU[i])
        {
          // double const sM = df( xmid );
          // std::cout << "{ " << _yL[i] << " " << xmid << " " << _yU[i] << " }
          // >> "; std::cout << "{ " << sL << " " << sM << " " << sU << " } >> [
          // ";
          std::tie(sL, sU) = extreme(sL, sU, df(xmid));
          // std::cout << sL << " " << sU << " ]\n";
        }

        double const sLL_sL = _sLL[i] * sL, sLL_sU = _sLL[i] * sU,
                     sUL_sL = _sUL[i] * sL, sUL_sU = _sUL[i] * sU;
        _sLL[i] = std::min(std::min(sLL_sL, sLL_sU), std::min(sUL_sL, sUL_sU));
        _sUL[i] = std::max(std::max(sLL_sL, sLL_sU), std::max(sUL_sL, sUL_sU));

        double const sLU_sL = _sLU[i] * sL, sLU_sU = _sLU[i] * sU,
                     sUU_sL = _sUU[i] * sL, sUU_sU = _sUU[i] * sU;
        _sLU[i] = std::max(std::max(sLU_sL, sLU_sU), std::max(sUU_sL, sUU_sU));
        _sUU[i] = std::min(std::min(sLU_sL, sLU_sU), std::min(sUU_sL, sUU_sU));

        _y0[i] = f(_y0[i]);
        _yL[i] = f(_yL[i]);
        _yU[i] = f(_yU[i]);
      }

      // Update endpoint
      _y0[N] = f(_y0[N]);
    }

    else
    {
      for (unsigned i = 0; i < N; ++i)
      {
        _yL[i] = f(_yL[i]);
        _yU[i] = f(_yU[i]);
      }
    }

    if (!inc)
    {
      std::swap(_yL, _yU);
    }

    // Backpropagate slopes to constant bounds
    if (options.SLOPEUSE > 1)
    {
      for (unsigned i = 0; i < N; ++i)
      {
        _yL[i] = _tighten(1, i, N);
        _yU[i] = _tighten(0, i, N);
      }
    }

    return *this;
  }
  /*
    // cvx=0: concave  cvx=1: convex
    template <typename UNIV, typename DUNIV>
    PWCU& _compose
      ( UNIV const& f, DUNIV const& df, int const cvx, bool const inc )
      {
        if( _yL.empty() || _yL.size() != _yU.size() )
          return *this;

        size_t const N = _yL.size();

        // Bounding valid regardless of convexity
        if( options.SLOPEUSE ){
  #ifdef MC__PWCU_CHECK
          if( _y0.size() != N+1 || _sLL.size() != N || _sLU.size() != N ||
  _sUL.size() != N || _sUU.size() != N ) throw Exceptions( Exceptions::SLOPE );
  #endif
          //double sL = df( _y0[0] ), sU, sM;

          for( unsigned i=0; i<N; ++i ){
            double const sL = df( _yL[i] );
            double const sU = df( _yU[i] );
            _yL[i] = f( _yL[i] );
            _yU[i] = f( _yU[i] );

            double const sLL_sL = _sLL[i]*sL, sLL_sU = _sLL[i]*sU, sUL_sL =
  _sUL[i]*sL, sUL_sU = _sUL[i]*sU; _sLL[i] = std::min( std::min( sLL_sL, sLL_sU
  ), std::min( sUL_sL, sUL_sU ) ); _sUL[i] = std::max( std::max( sLL_sL, sLL_sU
  ), std::max( sUL_sL, sUL_sU ) );

            double const sLU_sL = _sLU[i]*sL, sLU_sU = _sLU[i]*sU, sUU_sL =
  _sUU[i]*sL, sUU_sU = _sUU[i]*sU; _sLU[i] = std::max( std::max( sLU_sL, sLU_sU
  ), std::max( sUU_sL, sUU_sU ) ); _sUU[i] = std::min( std::min( sLU_sL, sLU_sU
  ), std::min( sUU_sL, sUU_sU ) );

            //double const dy0 = _y0[i+1] - _y0[i];
            //sM = ( std::fabs( dy0 ) > options.BKPTATOL?
            //       ( f( _y0[i+1] ) - f( _y0[i] ) ) / dy0:
            //       df( 0.5 * ( _y0[i+1] + _y0[i] ) )      );
            //sU = df( _y0[i+1] );

            //double const sLL_sL = _sLL[i]*sL, sLL_sM = _sLL[i]*sM, sUL_sL =
  _sUL[i]*sL, sUL_sM = _sUL[i]*sM;
            //_sLL[i] = std::min( std::min( sLL_sL, sLL_sM ), std::min( sUL_sL,
  sUL_sM ) );
            //_sUL[i] = std::max( std::max( sLL_sL, sLL_sM ), std::max( sUL_sL,
  sUL_sM ) );

            //double const sLU_sU = _sLU[i]*sU, sLU_sM = _sLU[i]*sM, sUU_sU =
  _sUU[i]*sU, sUU_sM = _sUU[i]*sM;
            //_sLU[i] = std::max( std::max( sLU_sU, sLU_sM ), std::max( sUU_sU,
  sUU_sM ) );
            //_sUU[i] = std::min( std::min( sLU_sU, sLU_sM ), std::min( sUU_sU,
  sUU_sM ) );

            _y0[i] = f( _y0[i] );
            //std::swap( sL, sU );
          }

          // Update endpoint
          _y0[N] = f( _y0[N] );
        }

        else{
          for( unsigned i=0; i<N; ++i ){
            _yL[i] = f( _yL[i] );
            _yU[i] = f( _yU[i] );
          }
        }

        if( !inc ){
          std::swap( _yL,  _yU );
          //std::swap( _sLL, _sUL );
          //std::swap( _sLU, _sUU );
        }

        // Backpropagate slopes to constant bounds
        if( options.SLOPEUSE > 1 ){
          for( unsigned i=0; i<N; ++i ){
            _yL[i] = _tighten( 1, i, N );
            _yU[i] = _tighten( 0, i, N );
          }
        }

        return *this;
      }

    // cvx=0: concave  cvx=1: convex
    template <typename UNIV, typename DUNIV>
    PWCU& _compose
      ( UNIV const& f, DUNIV const& df, int const cvx, bool const inc )
      {
        if( _yL.empty() || _yL.size() != _yU.size() )
          return *this;

        size_t const N = _yL.size();

        // Bounding valid regardless of convexity
        if( options.SLOPEUSE ){
  #ifdef MC__PWCU_CHECK
          if( _y0.size() != N+1 || _sLL.size() != N || _sLU.size() != N ||
  _sUL.size() != N || _sUU.size() != N ) throw Exceptions( Exceptions::SLOPE );
  #endif

          for( unsigned i=0; i<N; ++i ){
            _yL[i] = f( _yL[i] );
            _yU[i] = f( _yU[i] );
            // Convex increasing / concave decreasing univariate
            if( (cvx && inc) || (!cvx && !inc) ){
              _sLL[i] *= df( _y0[i] );
              _sLU[i] *= df( _y0[i+1] );
              double const dy0 = _y0[i+1] - _y0[i];
              double const sU = ( std::fabs( dy0 ) > options.BKPTATOL?
                                  ( f( _y0[i+1] ) - f( _y0[i] ) ) / dy0:
                                  df( 0.5 * ( _y0[i+1] + _y0[i] ) ) );
              _sUL[i] *= sU;
              _sUU[i] *= sU;
            }
            // Concave increasing / convex decreasing univariate
            else{
              _sUL[i] *= df( _y0[i] );
              _sUU[i] *= df( _y0[i+1] );
              double const dy0 = _y0[i+1] - _y0[i];
              double const sL = ( std::fabs( dy0 ) > options.BKPTATOL?
                                  ( f( _y0[i+1] ) - f( _y0[i] ) ) / dy0:
                                  df( 0.5 * ( _y0[i+1] + _y0[i] ) ) );
              _sLL[i] *= sL;
              _sLU[i] *= sL;
            }
            _y0[i] = f( _y0[i] );
          }
          // Update endpoint
          _y0[N] = f( _y0[N] );
        }

        else{
          for( unsigned i=0; i<N; ++i ){
            _yL[i] = f( _yL[i] );
            _yU[i] = f( _yU[i] );
          }
        }

        if( !inc ){
          std::swap( _yL,  _yU );
          std::swap( _sLL, _sUL );
          std::swap( _sLU, _sUU );
        }

        // Backpropagate slopes to constant bounds
        if( options.SLOPEUSE > 1 ){
          for( unsigned i=0; i<N; ++i ){
            _yL[i] = _tighten( 1, i, N );
            _yU[i] = _tighten( 0, i, N );
          }
        }

        return *this;
      }
  */
  PWCU&
  _neg()
  {
    if (_yL.empty() || _yL.size() != _yU.size()) return *this;

    size_t const N = _yL.size();

    // Bounding valid regardless of convexity
    if (options.SLOPEUSE)
    {
#ifdef MC__PWCU_CHECK
      if (_y0.size() != N + 1 || _sLL.size() != N || _sLU.size() != N ||
          _sUL.size() != N || _sUU.size() != N)
        throw Exceptions(Exceptions::SLOPE);
#endif
      for (unsigned i = 0; i < N; ++i)
      {
        _yL[i]  = -_yL[i];
        _yU[i]  = -_yU[i];
        _y0[i]  = -_y0[i];
        _sLL[i] = -_sLL[i];
        _sLU[i] = -_sLU[i];
        _sUL[i] = -_sUL[i];
        _sUU[i] = -_sUU[i];
      }
      _y0[N] = -_y0[N];

      std::swap(_yL, _yU);
      std::swap(_sLL, _sUL);
      std::swap(_sLU, _sUU);
    }

    else
    {
      for (unsigned i = 0; i < N; ++i)
      {
        _yL[i] = -_yL[i];
        _yU[i] = -_yU[i];
      }

      std::swap(_yL, _yU);
    }

    return *this;
  }

 public:
  //! @brief Options of mc::PWCU
  static struct Options
  {
    //! @brief Constructor
    Options() { reset(); }
    //! @brief Assignment
    Options&
    operator=(Options const& opt)
    {
      BKPTATOL = opt.BKPTATOL;
      BKPTRTOL = opt.BKPTRTOL;
      DISPNUM  = opt.DISPNUM;
      SLOPEUSE = opt.SLOPEUSE;
      return *this;
    }
    //! @brief Assignment
    void
    reset()
    {
      BKPTATOL = 1e2 * DBL_EPSILON;
      BKPTRTOL = 1e2 * DBL_EPSILON;
      DISPNUM  = 5;
      SLOPEUSE = 2;
    }
    //! @brief Absolute tolerance in breakpoints - Default: 1e2*DBL_EPSILON
    double BKPTATOL;
    //! @brief Relative tolerance in breakpoints - Default: 1e2*DBL_EPSILON
    double BKPTRTOL;
    //! @brief Number of numerical digits displayed with << operator - Default:
    //! 5
    unsigned DISPNUM;
    //! @brief Enable propagation of slopes - Default: 2
    unsigned SLOPEUSE;
  } options;

  //! @brief Exceptions of mc::PWCU
  class Exceptions
  {
   public:
    //! @brief Enumeration type for SElimVar exception handling
    enum TYPE
    {
      RANGE = 0,      //!< Operation on variable with empty range
      SIZE,           //!< Inconsistent vector size in estimator
      SLOPE,          //!< Inconsistent slope vector in estimator
      EXTRAPOL,       //!< Extrapolation outside of variable range
      DIV,            //!< Division by zero
      INTERNAL = -1,  //!< Internal error
      UNDEF    = -33  //!< Feature not yet implemented
    };
    //! @brief Constructor for error <a>ierr</a>
    Exceptions(TYPE ierr) : _ierr(ierr) {}
    //! @brief Error flag
    int
    ierr()
    {
      return _ierr;
    }
    //! @brief Error description
    std::string
    what()
    {
      switch (_ierr)
      {
        case RANGE:
          return "mc::PWCU\t Operation on variable with empty range";
        case SIZE:
          return "mc::PWCU\t Inconsistent vector size in estimator";
        case SLOPE:
          return "mc::PWCU\t Inconsistent slope vector in estimator";
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
  PWCU() : _yL(), _yU(), _y0(), _sLL(), _sLU(), _sUL(), _sUU() {}

  //! @brief Constructor of constant
  PWCU(double const& xL, double const& xU, double const& y, size_t const& N)
      : _xL(xL),
        _xU(xU),
        _yL(N > 1 ? N : 1, y - DBL_EPSILON),
        _yU(N > 1 ? N : 1, y + DBL_EPSILON)
  {
    if (options.SLOPEUSE)
    {
      size_t const n = N > 1 ? N : 1;
      _y0.assign(n + 1, y);
      _sLL.assign(n, 0.0);
      _sLU.assign(n, 0.0);
      _sUL.assign(n, 0.0);
      _sUU.assign(n, 0.0);
    }
  }

  //! @brief Constructor of variable
  PWCU(double xL, double xU, size_t const& N)
      : _xL(xL), _xU(xU), _yL(N > 1 ? N : 1), _yU(N > 1 ? N : 1)
  {
    size_t const n = N > 1 ? N : 1;
    xU             = (xU - xL) / n;
#ifdef MC__PWCU_CHECK
    if (xU <= 0.) throw Exceptions(Exceptions::RANGE);
#endif

    if (options.SLOPEUSE)
    {
      _y0.resize(n + 1);
      _sLL.assign(n, 1.0);
      _sLU.assign(n, 1.0);
      _sUL.assign(n, 1.0);
      _sUU.assign(n, 1.0);
    }

    for (size_t i = 0; i < n; ++i)
    {
      if (options.SLOPEUSE) _y0[i] = xL;
      _yL[i] = xL;
      xL += xU;
      _yU[i] = xL;
    }
    if (options.SLOPEUSE) _y0[n] = xL;
  }

  //! @brief Copy constructor
  PWCU(PWCU const& var)
      : _xL(var._xL),
        _xU(var._xU),
        _yL(var._yL),
        _yU(var._yU),
        _y0(var._y0),
        _sLL(var._sLL),
        _sLU(var._sLU),
        _sUL(var._sUL),
        _sUU(var._sUU)
  {
  }

  //! @brief Move constructor
  PWCU(PWCU&& var)
      : _xL(std::move(var._xL)),
        _xU(std::move(var._xU)),
        _yL(std::move(var._yL)),
        _yU(std::move(var._yU)),
        _y0(std::move(var._y0)),
        _sLL(std::move(var._sLL)),
        _sLU(std::move(var._sLU)),
        _sUL(std::move(var._sUL)),
        _sUU(std::move(var._sUU))
  {
  }

  //! @brief Set constant estimator
  PWCU&
  set(double const& y)
  {
    _yL.assign(_yL.size(), y - DBL_EPSILON);
    _yU.assign(_yU.size(), y + DBL_EPSILON);

    if (options.SLOPEUSE)
    {
      _y0.assign(_y0.size(), y);
      _sLL.assign(_yL.size(), 0.0);
      _sLU.assign(_yL.size(), 0.0);
      _sUL.assign(_yU.size(), 0.0);
      _sUU.assign(_yU.size(), 0.0);
    }
    else
    {
      _y0.clear();
      _sLL.clear();
      _sLU.clear();
      _sUL.clear();
      _sUU.clear();
    }

    return *this;
  }

  //! @brief Set constant estimator
  PWCU&
  set(double const& xL, double const& xU, double const& y, size_t const& N)
  {
    _xL = xL;
    _xU = xU;

    size_t const n = N > 1 ? N : 1;
    _yL.assign(n, y - DBL_EPSILON);
    _yU.assign(n, y + DBL_EPSILON);

    if (options.SLOPEUSE)
    {
      _y0.assign(n + 1, y);
      _sLL.assign(n, 0.0);
      _sLU.assign(n, 0.0);
      _sUL.assign(n, 0.0);
      _sUU.assign(n, 0.0);
    }
    else
    {
      _y0.clear();
      _sLL.clear();
      _sLU.clear();
      _sUL.clear();
      _sUU.clear();
    }

    return *this;
  }

  //! @brief Set variable estimator
  PWCU&
  set(double xL, double xU, size_t const N)
  {
    _xL = xL;
    _xU = xU;

    size_t const n = N > 1 ? N : 1;
    _yL.resize(n);
    _yU.resize(n);

    if (options.SLOPEUSE)
    {
      _y0.resize(n + 1);
      _sLL.assign(n, 1.0);
      _sLU.assign(n, 1.0);
      _sUL.assign(n, 1.0);
      _sUU.assign(n, 1.0);
    }
    else
    {
      _y0.clear();
      _sLL.clear();
      _sLU.clear();
      _sUL.clear();
      _sUU.clear();
    }

    xU = (xU - xL) / n;
#ifdef MC__PWCU_CHECK
    if (xU <= 0.) throw Exceptions(Exceptions::RANGE);
#endif
    for (size_t i = 0; i < n; ++i)
    {
      if (options.SLOPEUSE) _y0[i] = xL;
      _yL[i] = xL;
      xL += xU;
      _yU[i] = xL;
    }
    if (options.SLOPEUSE) _y0[n] = xL;

    return *this;
  }

  //! @brief Evaluate lower bound at point
  double
  l(double const& x, unsigned const opt = 1) const
  {
#ifdef MC__PWCU_CHECK
    if (_yL.empty()) throw Exceptions(Exceptions::SIZE);
#endif
    size_t const N  = _yL.size();
    double const dx = (_xU - _xL) / N;
    double xi       = _xL;

    if (opt && options.SLOPEUSE)
    {
#ifdef MC__PWCU_CHECK
      if (_sLL.empty() || _sLU.empty()) throw Exceptions(Exceptions::SLOPE);
#endif
      for (unsigned i = 0; i < N; ++i)
      {
        if (xi + dx < x)
        {
          xi += dx;
          continue;
        }
        return std::max(_yL[i],
                        std::max(_y0[i] + _sLL[i] * (x - xi),
                                 _y0[i + 1] + _sLU[i] * (x - (xi + dx))));
      }
      if (isequal(xi, x, options.BKPTATOL, options.BKPTRTOL)) return _y0[N];
    }

    else
    {
      for (unsigned i = 0; i < N; ++i)
      {
        if (xi + dx < x)
        {
          xi += dx;
          continue;
        }
        return _yL[i];
      }
      if (isequal(xi, x, options.BKPTATOL, options.BKPTRTOL)) return _yL.back();
    }

    throw Exceptions(Exceptions::EXTRAPOL);
  }

  //! @brief Evaluate upper bound at point
  double
  u(double const& x, unsigned const opt = 1) const
  {
#ifdef MC__PWCU_CHECK
    if (_yU.empty()) throw Exceptions(Exceptions::SIZE);
#endif
    size_t const N  = _yL.size();
    double const dx = (_xU - _xL) / N;
    double xi       = _xL;

    if (opt && options.SLOPEUSE)
    {
#ifdef MC__PWCU_CHECK
      if (_sUL.empty() || _sUU.empty()) throw Exceptions(Exceptions::SLOPE);
#endif
      for (unsigned i = 0; i < N; ++i)
      {
        if (xi + dx < x)
        {
          xi += dx;
          continue;
        }
        return std::min(_yU[i],
                        std::min(_y0[i] + _sUL[i] * (x - xi),
                                 _y0[i + 1] + _sUU[i] * (x - (xi + dx))));
      }
      if (isequal(xi, x, options.BKPTATOL, options.BKPTRTOL)) return _y0[N];
    }

    else
    {
      for (unsigned i = 0; i < N; ++i)
      {
        if (xi + dx < x)
        {
          xi += dx;
          continue;
        }
        return _yU[i];
      }
      if (isequal(xi, x, options.BKPTATOL, options.BKPTRTOL)) return _yU.back();
    }

    throw Exceptions(Exceptions::EXTRAPOL);
  }

  //! @brief Evaluate lower range
  double
  l() const
  {
#ifdef MC__PWCU_CHECK
    if (_yL.empty()) throw Exceptions(Exceptions::SIZE);
#endif
    return *std::min_element(_yL.cbegin(), _yL.cend());
  }

  //! @brief Evaluate upper range
  double
  u() const
  {
#ifdef MC__PWCU_CHECK
    if (_yU.empty()) throw Exceptions(Exceptions::SIZE);
#endif
    return *std::max_element(_yU.cbegin(), _yU.cend());
  }

  //! @brief Evaluate range
  double
  w() const
  {
    return std::max(u() - l(), 0.);
  }

  //! @brief Display estimator
  std::ostream&
  display(std::ostream& os     = std::cout,
          size_t const dispnum = options.DISPNUM) const
  {
#ifdef MC__PWCU_CHECK
    if (_yL.empty() || _yL.size() != _yU.size())
      throw Exceptions(Exceptions::SIZE);
#endif
    size_t const N  = _yL.size();
    double const dx = (_xU - _xL) / N;
    double xi       = _xL;

    if (options.SLOPEUSE)
    {
#ifdef MC__PWCU_CHECK
      if (_y0.size() != N + 1 || _sLL.size() != N || _sLU.size() != N ||
          _sUL.size() != N || _sUU.size() != N)
        throw Exceptions(Exceptions::SLOPE);
#endif
      os << "{" << std::scientific << std::setprecision(dispnum) << std::right
         << std::setw(dispnum + 8) << _xL << ":" << std::setw(dispnum + 8)
         << _y0[0];
      for (unsigned i = 1; i <= N; ++i)
        os << ", " << std::setw(dispnum + 8) << (xi += dx) << ": "
           << std::setw(dispnum + 8) << _y0[i];
      return os  //<< ", " << std::setw(dispnum+8) << _xU << ":" <<
                 // std::setw(dispnum+8) << _y0[N]
             << " }";
    }

    else
    {
      os << "{" << std::scientific << std::setprecision(dispnum) << std::right
         << std::setw(dispnum + 7) << xi;
      for (unsigned i = 0; i < N; ++i)
        os << " <" << std::setw(dispnum + 7) << _yL[i] << " : "
           << std::setw(dispnum + 7) << _yU[i] << "> " << std::setw(dispnum + 7)
           << (xi += dx);
    }

    return os << " }";
  }

  //! @brief Retreive number of segments
  size_t
  size() const
  {
#ifdef MC__PWCU_CHECK
    if (_yL.size() != _yU.size()) throw Exceptions(Exceptions::SIZE);
#endif
    return _yL.size();
  }

  //! @brief Retreive lower levels
  std::vector<double> const&
  yL() const
  {
    return _yL;
  }

  //! @brief Retreive upper levels
  std::vector<double> const&
  yU() const
  {
    return _yU;
  }

  //! @brief Retreive knot values
  std::vector<double> const&
  y0() const
  {
    return _y0;
  }

  //! @brief Retreive lower slopes at start
  std::vector<double> const&
  sLL() const
  {
    return _sLL;
  }

  //! @brief Retreive lower slopes at end
  std::vector<double> const&
  sLU() const
  {
    return _sLU;
  }

  //! @brief Retreive upper slopes at start
  std::vector<double> const&
  sUL() const
  {
    return _sUL;
  }

  //! @brief Retreive upper slopes at end
  std::vector<double> const&
  sUU() const
  {
    return _sUU;
  }

  //! @brief Retreive initial abscissa
  double
  xL() const
  {
    return _xL;
  }

  //! @brief Retreive final abscissa
  double
  xU() const
  {
    return _xU;
  }

  //! @brief Retreive/set lower levels
  std::vector<double>&
  yL()
  {
    return _yL;
  }

  //! @brief Retreive/set upper levels
  std::vector<double>&
  yU()
  {
    return _yU;
  }

  //! @brief Retreive/set knot values
  std::vector<double>&
  y0()
  {
    return _y0;
  }

  //! @brief Retreive/set lower slopes at start
  std::vector<double>&
  sLL()
  {
    return _sLL;
  }

  //! @brief Retreive/set lower slopes at end
  std::vector<double>&
  sLU()
  {
    return _sLU;
  }

  //! @brief Retreive/set upper slopes at start
  std::vector<double>&
  sUL()
  {
    return _sUL;
  }

  //! @brief Retreive/set upper slopes at end
  std::vector<double>&
  sUU()
  {
    return _sUU;
  }

  //! @brief Retreive/set initial abscissa
  double&
  xL()
  {
    return _xL;
  }

  //! @brief Retreive/set final abscissa
  double&
  xU()
  {
    return _xU;
  }

  PWCU&
  operator=(PWCU const& var)
  {
    _xL = var._xL;
    _xU = var._xU;
    _yL = var._yL;
    _yU = var._yU;

    if (options.SLOPEUSE)
    {
      _y0  = var._y0;
      _sLL = var._sLL;
      _sLU = var._sLU;
      _sUL = var._sUL;
      _sUU = var._sUU;
    }

    return *this;
  }

  PWCU&
  operator=(PWCU&& var)
  {
    _xL = std::move(var._xL);
    _xU = std::move(var._xU);
    _yL = std::move(var._yL);
    _yU = std::move(var._yU);

    if (options.SLOPEUSE)
    {
      _y0  = std::move(var._y0);
      _sLL = std::move(var._sLL);
      _sLU = std::move(var._sLU);
      _sUL = std::move(var._sUL);
      _sUU = std::move(var._sUU);
    }

    return *this;
  }

  PWCU&
  operator+=(double const& cst)
  {
    if (cst == 0.) return *this;

    size_t const N = _yL.size();

    if (options.SLOPEUSE)
    {
#ifdef MC__PWCU_CHECK
      if (_y0.size() != N + 1 || _sLL.size() != N || _sLU.size() != N ||
          _sUL.size() != N || _sUU.size() != N)
        throw Exceptions(Exceptions::SLOPE);
#endif
      for (size_t i = 0; i < N; ++i)
      {
        _yL[i] += cst;
        _yU[i] += cst;
        _y0[i] += cst;
      }
      _y0[N] += cst;
    }

    else
    {
      for (size_t i = 0; i < N; ++i)
      {
        _yL[i] += cst;
        _yU[i] += cst;
      }
    }

    return *this;
  }

  PWCU&
  operator-=(double const& cst)
  {
    if (cst == 0.) return *this;

    size_t const N = _yL.size();

    if (options.SLOPEUSE)
    {
#ifdef MC__PWCU_CHECK
      if (_y0.size() != N + 1 || _sLL.size() != N || _sLU.size() != N ||
          _sUL.size() != N || _sUU.size() != N)
        throw Exceptions(Exceptions::SLOPE);
#endif
      for (size_t i = 0; i < N; ++i)
      {
        _yL[i] -= cst;
        _yU[i] -= cst;
        _y0[i] -= cst;
      }
      _y0[N] -= cst;
    }

    else
    {
      for (size_t i = 0; i < N; ++i)
      {
        _yL[i] -= cst;
        _yU[i] -= cst;
      }
    }

    return *this;
  }

  PWCU&
  operator*=(double const& cst)
  {
    if (cst == 1.) return *this;

    size_t const N = _yL.size();

    if (options.SLOPEUSE)
    {
#ifdef MC__PWCU_CHECK
      if (_y0.size() != N + 1 || _sLL.size() != N || _sLU.size() != N ||
          _sUL.size() != N || _sUU.size() != N)
        throw Exceptions(Exceptions::SLOPE);
#endif
      for (size_t i = 0; i < N; ++i)
      {
        _yL[i] *= cst;
        _yU[i] *= cst;
        _y0[i] *= cst;
        _sLL[i] *= cst;
        _sLU[i] *= cst;
        _sUL[i] *= cst;
        _sUU[i] *= cst;
      }
      _y0[N] *= cst;
      if (cst < 0.)
      {
        std::swap(_yL, _yU);
        std::swap(_sLL, _sUL);
        std::swap(_sLU, _sUU);
      }
    }

    else
    {
      for (size_t i = 0; i < N; ++i)
      {
        _yL[i] *= cst;
        _yU[i] *= cst;
      }
      if (cst < 0.) std::swap(_yL, _yU);
    }

    return *this;
  }

  PWCU&
  operator/=(double const& cst)
  {
    if (cst == 1.)
      return *this;
    else if (cst == 0.)
      throw Exceptions(Exceptions::DIV);

    size_t const N = _yL.size();

    if (options.SLOPEUSE)
    {
#ifdef MC__PWCU_CHECK
      if (_y0.size() != N + 1 || _sLL.size() != N || _sLU.size() != N ||
          _sUL.size() != N || _sUU.size() != N)
        throw Exceptions(Exceptions::SLOPE);
#endif
      for (size_t i = 0; i < N; ++i)
      {
        _yL[i] /= cst;
        _yU[i] /= cst;
        _y0[i] /= cst;
        _sLL[i] /= cst;
        _sLU[i] /= cst;
        _sUL[i] /= cst;
        _sUU[i] /= cst;
      }
      _y0[N] /= cst;
      if (cst < 0.)
      {
        std::swap(_yL, _yU);
        std::swap(_sLL, _sUL);
        std::swap(_sLU, _sUU);
      }
    }

    else
    {
      for (size_t i = 0; i < N; ++i)
      {
        _yL[i] /= cst;
        _yU[i] /= cst;
      }
      if (cst < 0.) std::swap(_yL, _yU);
    }

    return *this;
  }

  PWCU&
  operator+=(PWCU const& var)
  {
#ifdef MC__PWCU_CHECK
    if (_yL.empty() || _yL.size() != var._yL.size() || _yU.empty() ||
        _yU.size() != var._yU.size())
      throw Exceptions(Exceptions::SIZE);
#endif

    size_t const N = _yL.size();

    if (options.SLOPEUSE)
    {
#ifdef MC__PWCU_CHECK
      if (_y0.size() != N + 1 || _sLL.size() != N || _sLU.size() != N ||
          _sUL.size() != N || _sUU.size() != N)
        throw Exceptions(Exceptions::SLOPE);
#endif
      for (size_t i = 0; i < N; ++i)
      {
        _yL[i] += var._yL[i];
        _yU[i] += var._yU[i];
        _y0[i] += var._y0[i];
        _sLL[i] += var._sLL[i];
        _sLU[i] += var._sLU[i];
        _sUL[i] += var._sUL[i];
        _sUU[i] += var._sUU[i];
      }
      _y0[N] += var._y0[N];
    }

    else
    {
      for (size_t i = 0; i < N; ++i)
      {
        _yL[i] += var._yL[i];
        _yU[i] += var._yU[i];
      }
    }

    // Backpropagate slopes to constant bounds
    if (options.SLOPEUSE > 1)
    {
      for (size_t i = 0; i < N; ++i)
      {
        _yL[i] = _tighten(1, i, N);
        _yU[i] = _tighten(0, i, N);
      }
    }

    return *this;
  }

  PWCU&
  operator-=(PWCU const& var)
  {
#ifdef MC__PWCU_CHECK
    if (_yL.empty() || _yL.size() != var._yL.size() || _yU.empty() ||
        _yU.size() != var._yU.size())
      throw Exceptions(Exceptions::SIZE);
#endif

    size_t const N = _yL.size();

    if (options.SLOPEUSE)
    {
#ifdef MC__PWCU_CHECK
      if (_y0.size() != N + 1 || _sLL.size() != N || _sLU.size() != N ||
          _sUL.size() != N || _sUU.size() != N)
        throw Exceptions(Exceptions::SLOPE);
#endif
      for (size_t i = 0; i < N; ++i)
      {
        _yL[i] -= var._yU[i];
        _yU[i] -= var._yL[i];
        _y0[i] -= var._y0[i];
        _sLL[i] -= var._sUL[i];
        _sLU[i] -= var._sUU[i];
        _sUL[i] -= var._sLL[i];
        _sUU[i] -= var._sLU[i];
      }
      _y0[N] -= var._y0[N];
    }

    else
    {
      for (size_t i = 0; i < N; ++i)
      {
        _yL[i] -= var._yU[i];
        _yU[i] -= var._yL[i];
      }
    }

    // Backpropagate slopes to constant bounds
    if (options.SLOPEUSE > 1)
    {
      for (size_t i = 0; i < N; ++i)
      {
        _yL[i] = _tighten(1, i, N);
        _yU[i] = _tighten(0, i, N);
      }
    }

    return *this;
  }

  PWCU&
  neg()
  {
    return _neg();
  }

  // cvx=0: concave  cvx=1: convex  cvx=2: concavoconvex  cvx=3: convexoconcave
  // cvx=4: concave nonotonic  cvx=5: convex non-monotonic
  template <typename UNIV, typename DUNIV>
  PWCU&
  compose(UNIV const& f, DUNIV const& df, bool const under, int const cvx,
          bool const inc = 0, double const& xmid = 0.0)
  {
    if (_yL.empty() || _yL.size() != _yU.size()) return *this;

    if (cvx < 0 || cvx > 5) throw Exceptions(Exceptions::INTERNAL);

    return this->_compose(f, df, cvx, inc, xmid);
    /*
          // univariate is either convex or concave
          if( cvx == 0 || cvx == 1 )
            return this->_compose( f, df, cvx, inc );

          // univariate is either concavoconvex or convexoconcave with
       inflection at xmid else if( cvx == 2 || cvx == 3 ){ double const fmid =
       f( xmid ), dfmid = df( xmid ); auto const& fcv  = [&]( const double& x )
                                  { double const z = f(x), t =
       fmid+dfmid*(x-xmid); return z>t?z:t; }; auto const& dfcv = [&]( const
       double& x ) { double const z = f(x), t = fmid+dfmid*(x-xmid); return
       z>t?df(x):dfmid; }; auto const& fcc  = [&]( const double& x ) { double
       const z = f(x)-fmid, t = dfmid*(x-xmid); return z<t?z-t:0; }; auto const&
       dfcc = [&]( const double& x ) { double const z = f(x)-fmid, t =
       dfmid*(x-xmid); return z<t?df(x)-dfmid:0; }; PWCU copy( *this ); if( cvx
       == 2 ) return this->_compose( fcv, dfcv, 1, inc ) += copy._compose( fcc,
       dfcc, 0, 1 ); else return this->_compose( fcv, dfcv, 1, inc ) +=
       copy._compose( fcc, dfcc, 0, 0 );
          }

          // univariate is either convex or concave non-monotonic with optimum
       at xmid else if( cvx == 4 || cvx == 5 ){ double const fmid = f( xmid );
            auto const& fr  = [&]( const double& x )
                                 { return x>xmid?f(x):fmid; };
            auto const& dfr = [&]( const double& x )
                                 { return x>xmid?df(x):0; };
            auto const& fl  = [&]( const double& x )
                                 { return x<xmid?f(x)-fmid:0; };
            auto const& dfl = [&]( const double& x )
                                 { return x<xmid?df(x):0; };

            PWCU copy( *this );
            if( cvx == 4 )
              return this->_compose( fr, dfr, 0, 0 ) += copy._compose( fl, dfl,
       0, 1 ); else return this->_compose( fr, dfr, 1, 1 ) += copy._compose( fl,
       dfl, 1, 0 );
          }

          throw Exceptions( Exceptions::INTERNAL );
    */
  }

  PWCU&
  min(double const& c)
  {
    if (_yL.empty() || _yL.size() != _yU.size()) return *this;

    auto const& fmin  = [&](const double& x) { return x < c ? x : c; };
    auto const& dfmin = [&](const double& x) { return x < c ? 1 : 0; };

    return this->_compose(fmin, dfmin, 0, 1);
  }

  PWCU&
  max(double const& c)
  {
    if (_yL.empty() || _yL.size() != _yU.size()) return *this;

    auto const& fmax  = [&](const double& x) { return x > c ? x : c; };
    auto const& dfmax = [&](const double& x) { return x > c ? 1 : 0; };

    return this->_compose(fmax, dfmax, 1, 1);
  }

  PWCU&
  reduce(bool const under, size_t const nseg)
  {
    return *this;
  }

  PWCU&
  clean(bool const under)
  {
    return *this;
  }
};

////////////////////////////////////////////////////////////////////////

inline PWCU::Options PWCU::options;

inline std::ostream&
operator<<(std::ostream& os, PWCU const& var)
{
  var.display(os);
  return os;
}

inline PWCU
operator+(PWCU const& var)
{
  return var;
}

inline PWCU
operator+(PWCU&& var)
{
  PWCU res(std::move(var));
  return var;
}

inline PWCU
operator+(PWCU const& var1, PWCU const& var2)
{
  PWCU res(var1);
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
inline PWCU
operator+(PWCU&& var1, PWCU const& var2)
{
  PWCU res(std::move(var1));
  res += var2;
  return res;
}

inline PWCU
operator+(PWCU const& var1, double const& cst2)
{
  PWCU res(var1);
  res += cst2;
  return res;
}

inline PWCU
operator+(PWCU&& var1, double const& cst2)
{
  PWCU res(std::move(var1));
  res += cst2;
  return res;
}

inline PWCU
operator+(double const& cst1, PWCU const& var2)
{
  PWCU res(var2);
  res += cst1;
  return res;
}

inline PWCU
operator+(double const& cst1, PWCU&& var2)
{
  PWCU res(std::move(var2));
  res += cst1;
  return res;
}

inline PWCU
operator-(PWCU const& var)
{
  PWCU res(var);
  return (res *= -1);
}

inline PWCU
operator-(PWCU&& var)
{
  PWCU res(std::move(var));
  return (res *= -1);
}

inline PWCU
operator-(PWCU const& var1, PWCU const& var2)
{
  PWCU res(var1);
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
inline PWCU
operator-(PWCU&& var1, PWCU const& var2)
{
  PWCU res(std::move(var1));
  res -= var2;
  return res;
}

inline PWCU
operator-(PWCU const& var1, double const& cst2)
{
  PWCU res(var1);
  res += -cst2;
  return res;
}

inline PWCU
operator-(PWCU&& var1, double const& cst2)
{
  PWCU res(std::move(var1));
  res += -cst2;
  return res;
}

inline PWCU
operator-(double const& cst1, PWCU const& var2)
{
  PWCU res(-var2);
  res += cst1;
  return res;
}

inline PWCU
operator-(double const& cst1, PWCU&& var2)
{
  PWCU res(-std::move(var2));
  res += cst1;
  return res;
}

inline PWCU
operator*(PWCU const& var1, double const& cst2)
{
  PWCU res(var1);
  res *= cst2;
  return res;
}

inline PWCU
operator*(PWCU&& var1, double const& cst2)
{
  PWCU res(std::move(var1));
  res *= cst2;
  return res;
}

inline PWCU
operator*(double const& cst1, PWCU const& var2)
{
  PWCU res(var2);
  res *= cst1;
  return res;
}

inline PWCU
operator*(double const& cst1, PWCU&& var2)
{
  PWCU res(std::move(var2));
  res *= cst1;
  return res;
}

inline PWCU
operator/(PWCU const& var1, double const& cst2)
{
  PWCU var3(var1);
  var3 /= cst2;
  return var3;
}

inline PWCU
operator/(PWCU&& var1, double const& cst2)
{
  PWCU var3(std::move(var1));
  var3 /= cst2;
  return var3;
}

}  // namespace mc

#endif
