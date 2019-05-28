// Copyright (C) 2009-2019 Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the Eclipse Public License.

#ifndef MC__MCBOOST_HPP
#define MC__MCBOOST_HPP

#include <boost/numeric/interval.hpp>
#include "mcfunc.hpp"
#include "fadbad.h"

//namespace mc
//{
//typedef boost::numeric::interval_lib::save_state<boost::numeric::interval_lib::rounded_transc_opp<double>> T_boost_round;
//typedef boost::numeric::interval_lib::checking_strict<double> T_boost_check;
//typedef boost::numeric::interval<double,boost::numeric::interval_lib::policies<T_boost_round,T_boost_check>> T_policy;
//}

namespace fadbad
{
//! @brief Specialization of the structure fadbad::Op for use of the type boost::numeric::interval_lib::interval<double> of the <A href="https://www.boost.org/doc/libs/1_68_0/libs/numeric/interval/doc/interval.htm">Boost Interval Arithmetic Library</A> as a template parameter of the classes fadbad::F, fadbad::B and fadbad::T of FADBAD++
template <> struct Op< boost::numeric::interval < double,boost::numeric::interval_lib::policies<
                       boost::numeric::interval_lib::save_state<boost::numeric::interval_lib::rounded_transc_opp<double>>,
                       boost::numeric::interval_lib::checking_strict<double>
                     > > >
{
  typedef double Base;
  typedef boost::numeric::interval_lib::save_state<boost::numeric::interval_lib::rounded_transc_opp<double>> R;
  typedef boost::numeric::interval_lib::checking_strict<double> P;
  typedef boost::numeric::interval<double,boost::numeric::interval_lib::policies<R,P>> T;
  //typedef boost::numeric::interval<double,T_boost_policy> T;
  static Base myInteger( const int i ) { return Base(i); }
  static Base myZero() { return myInteger(0); }
  static Base myOne() { return myInteger(1);}
  static Base myTwo() { return myInteger(2); }
  static double myPI() { return mc::PI; }
  static T myPos( const T& x ) { return  x; }
  static T myNeg( const T& x ) { return -x; }
  template <typename U> static T& myCadd( T& x, const U& y ) { return x+=y; }
  template <typename U> static T& myCsub( T& x, const U& y ) { return x-=y; }
  template <typename U> static T& myCmul( T& x, const U& y ) { return x*=y; }
  template <typename U> static T& myCdiv( T& x, const U& y ) { return x/=y; }
  static T myInv( const T& x ) { return T(1.)/x; }
  static T mySqr( const T& x ) { return boost::numeric::square(x); }
  static T myPow( const T& x, const int n ) { return boost::numeric::pow(x,n); }
  template <typename X, typename Y> static T pow(const X& x, const Y& y) {
    return boost::numeric::exp( y * boost::numeric::log(x) );
  }
  static T mySqrt( const T& x ) { return boost::numeric::sqrt(x); }
  static T myLog( const T& x ) { return boost::numeric::log(x); }
  static T myExp( const T& x ) { return boost::numeric::exp(x); }
  static T mySin( const T& x ) { return boost::numeric::sin( x ); }
  static T myCos( const T& x ) { return boost::numeric::cos( x ); }
  static T myTan( const T& x ) { return boost::numeric::tan( x ); }
  static T myAsin( const T& x ) { return boost::numeric::asin( x ); }
  static T myAcos( const T& x ) { return boost::numeric::acos( x ); }
  static T myAtan( const T& x ) { return boost::numeric::atan( x ); }
  static T mySinh( const T& x ) { return boost::numeric::sinh( x ); }
  static T myCosh( const T& x ) { return boost::numeric::cosh( x ); }
  static T myTanh( const T& x ) { return boost::numeric::tanh( x ); }
  static bool myEq(const T& x, const T& y) { return boost::numeric::equal( x, y ); }
  static bool myNe(const T& x, const T& y) { return !boost::numeric::equal( x, y ); }
  static bool myLt(const T& x, const T& y) { return boost::numeric::proper_subset( x, y );  }
  static bool myLe(const T& x, const T& y) { return boost::numeric::subset( x, y ); }
  static bool myGt(const T& x, const T& y) { return boost::numeric::proper_subset( y, x );  }
  static bool myGe(const T& x, const T& y) { return boost::numeric::subset( y, x ); }
};

} // end namespace fadbad

#include "mcop.hpp"

namespace mc
{
//! @brief Specialization of the structure mc::Op for use of the type boost:numeric:interval_lib::interval<double> of the <A href="https://www.boost.org/doc/libs/1_68_0/libs/numeric/interval/doc/interval.htm">Boost Interval Arithmetic Library</A> as a template parameter in other MC++ types
//template <> struct Op< boost::numeric::interval<double> >
//{
//  typedef boost::numeric::interval<double> T;
template <> struct Op< boost::numeric::interval < double,boost::numeric::interval_lib::policies<
                       boost::numeric::interval_lib::save_state<boost::numeric::interval_lib::rounded_transc_opp<double>>,
                       boost::numeric::interval_lib::checking_strict<double>
                     > > >
{
  typedef boost::numeric::interval_lib::save_state<boost::numeric::interval_lib::rounded_transc_opp<double>> R;
  typedef boost::numeric::interval_lib::checking_strict<double> P;
  typedef boost::numeric::interval<double,boost::numeric::interval_lib::policies<R,P>> T;
  static T point( const double c ) { return T( c ); }
  static T zeroone() { return T( 0., 1. ); }
  static void I(T& x, const T& y) { x = y; }
  static double l(const T& x) { return boost::numeric::lower(x); }
  static double u(const T& x) { return boost::numeric::upper(x); }
  static double abs (const T& x) { return boost::numeric::norm(x);  }
  static double mid (const T& x) { return boost::numeric::median(x);  }
  static double diam(const T& x) { return boost::numeric::width(x); }
  static T inv (const T& x) { return boost::numeric::interval_lib::multiplicative_inverse(x); }
  static T sqr (const T& x) { return boost::numeric::square(x);  }
  static T sqrt(const T& x) { return boost::numeric::sqrt(x); }
  static T exp (const T& x) { return boost::numeric::exp(x);  }
  static T log (const T& x) { return boost::numeric::log(x);  }
  static T xlog0(const T& x) { return boost::numeric::log(x)*x; }
  static T xlog (const T& x) {
    T zmin = boost::numeric::exp( T(-1.) );
    T mono = boost::numeric::hull( xlog0( boost::numeric::lower(x) ), xlog0( boost::numeric::upper(x) ) );
    return overlap( x, zmin )? boost::numeric::hull( mono, xlog0(zmin) ): mono; 
  }
  static T lmtd0(const T& x,const T& y){ return x == y? x: (x-y)/(boost::numeric::log(x)-boost::numeric::log(y)); }
  static T lmtd (const T& x,const T& y){
    T vmin = lmtd0( boost::numeric::lower(x), boost::numeric::lower(y) );
    T vmax = lmtd0( boost::numeric::upper(x), boost::numeric::upper(y) );
    return boost::numeric::hull( vmin, vmax );
  }
  static T rlmtd0(const T& x,const T& y){ return x == y? inv(x): (boost::numeric::log(x)-boost::numeric::log(y))/(x-y); }
  static T rlmtd (const T& x,const T& y){
    T vmin = rlmtd0( boost::numeric::upper(x), boost::numeric::upper(y) );
    T vmax = rlmtd0( boost::numeric::lower(x), boost::numeric::lower(y) );
    return boost::numeric::hull( vmin, vmax );
  }
  static T fabs(const T& x) { return boost::numeric::abs(x); }
  static T sin (const T& x) { return boost::numeric::sin(x);  }
  static T cos (const T& x) { return boost::numeric::cos(x);  }
  static T tan (const T& x) { return boost::numeric::tan(x);  }
  static T asin(const T& x) { return boost::numeric::asin(x); }
  static T acos(const T& x) { return boost::numeric::acos(x); }
  static T atan(const T& x) { return boost::numeric::atan(x); }
  static T sinh(const T& x) { return boost::numeric::sinh(x); }
  static T cosh(const T& x) { return boost::numeric::cosh(x); }
  static T tanh(const T& x) { return boost::numeric::tanh(x); }
  static T erf (const T& x) { 
    return T( std::erf( boost::numeric::lower(x) ), std::erf( boost::numeric::upper(x) ) );
  }
  static T erfc(const T& x) {
    return T( std::erfc( boost::numeric::upper(x) ), std::erfc( boost::numeric::lower(x) ) );
  }
  static T fstep(const T& x) {
    if( boost::numeric::upper(x) < 0 ) return T(0); 
    if( boost::numeric::lower(x) > 0 ) return T(1); 
    return T(0,1);
  }
  static T bstep(const T& x) {
    if( boost::numeric::lower(x) > 0 ) return T(0); 
    if( boost::numeric::upper(x) < 0 ) return T(1); 
    return T(0,1);
  }
  static T hull(const T& x, const T& y) { return boost::numeric::hull( x, y ); }
  static T min (const T& x, const T& y) { return boost::numeric::min( x, y ); }
  static T max (const T& x, const T& y) { return boost::numeric::max( x, y ); }
  static T arh (const T& x, const double k) { return boost::numeric::exp( -x / k ); }
  static T cheb (const T& x, const unsigned n) { return T(-1.,1.); }
  static T pow(const T& x, const int n) { return boost::numeric::pow( x, n ); }
  template <typename X, typename Y> static T pow(const X& x, const Y& y) {
    return boost::numeric::exp( y * boost::numeric::log(x) );
  }
  static T prod (const unsigned int n, const T* x) { return n? x[0] * prod(n-1, x+1): 1.; }
  static T monom (const unsigned int n, const T* x, const unsigned* k) {
    return n? boost::numeric::pow(x[0], k[0]) * monom(n-1, x+1, k+1): 1.;
  }
  static bool inter(T& xIy, const T& x, const T& y) {
    xIy = boost::numeric::intersect(x,y);
    return !boost::numeric::empty( xIy );
  }  
  static bool eq(const T& x, const T& y) { return boost::numeric::equal( x, y ); }
  static bool ne(const T& x, const T& y) { return !boost::numeric::equal( x, y ); }
  static bool lt(const T& x, const T& y) { return boost::numeric::proper_subset( x, y );  }
  static bool le(const T& x, const T& y) { return boost::numeric::subset( x, y ); }
  static bool gt(const T& x, const T& y) { return boost::numeric::proper_subset( y, x );  }
  static bool ge(const T& x, const T& y) { return boost::numeric::subset( y, x ); }
};

} // namespace mc

#include <iostream>
#include <iomanip>

namespace boost{
namespace numeric{

template <typename POLICY> 
inline std::ostream&
operator<<
( std::ostream&os, const boost::numeric::interval<double,POLICY>&I)
{
  const int DISPLAY_DIGITS = 16;
  os << std::right << std::scientific << std::setprecision(DISPLAY_DIGITS);
  os << "[ "  << std::setw(DISPLAY_DIGITS+7) << I.lower()
     << " : " << std::setw(DISPLAY_DIGITS+7) << I.upper() << " ]";
  return os;
}

} // namespace numeric
} // namespace boost
#endif
