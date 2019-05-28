// Copyright (C) 2009-2019 Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the Eclipse Public License.

#ifndef MC__MCPROFIL_HPP
#define MC__MCPROFIL_HPP

#include "mcfunc.hpp"
#include <Interval.h>
#include <Functions.h>
#include <Constants.h>

#include "fadbad.h"

namespace fadbad
{

//! @brief Specialization of the structure fadbad::Op for use of the type ::INTERVAL of <A href="http://www.ti3.tu-harburg.de/Software/PROFILEnglisch.html">PROFIL</A> as a template parameter of the classes fadbad::F, fadbad::B and fadbad::T of FADBAD++
template <> struct Op< ::INTERVAL >
{
  typedef double Base;
  typedef ::INTERVAL T;
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
  static T mySqr( const T& x ) { return ::Sqr(x); }
  static T myPow( const T& x, const int n ) {
    return( (n>=3&&n%2)? T(::Power(Inf(x),n),::Power(Sup(x),n)): ::Power(x,n) );
  }
  template <typename X, typename Y> static T myPow( const X& x, const Y& y ) { return ::Power(x,y); }
  //static T myCheb( const T& x, const unsigned n ) { return T(-1.,1.); }
  static T mySqrt( const T& x ) { return ::Sqrt( x ); }
  static T myLog( const T& x ) { return ::Log( x ); }
  static T myExp( const T& x ) { return ::Exp( x ); }
  static T mySin( const T& x ) { return ::Sin( x ); }
  static T myCos( const T& x ) { return ::Cos( x ); }
  static T myTan( const T& x ) { return ::Tan( x ); }
  static T myAsin( const T& x ) { return ::ArcSin( x ); }
  static T myAcos( const T& x ) { return ::ArcCos( x ); }
  static T myAtan( const T& x ) { return ::ArcTan( x ); }
  static T mySinh( const T& x ) { return ::Sinh( x ); }
  static T myCosh( const T& x ) { return ::Cosh( x ); }
  static T myTanh( const T& x ) { return ::Tanh( x ); }
  static bool myEq( const T& x, const T& y ) { return x == y; }
  static bool myNe( const T& x, const T& y ) { return x != y; }
  static bool myLt( const T& x, const T& y ) { return x <  y; }
  static bool myLe( const T& x, const T& y ) { return x <= y; }
  static bool myGt( const T& x, const T& y ) { return y <  x; }
  static bool myGe( const T& x, const T& y ) { return y <= x; }
};

} // end namespace fadbad

#include "mcop.hpp"

namespace mc
{

//! @brief Specialization of the structure mc::Op for use of the type ::INTERVAL of <A href="http://www.ti3.tu-harburg.de/Software/PROFILEnglisch.html">PROFIL</A> as a template parameter in the other MC++ types
template <> struct Op< ::INTERVAL >
{
  typedef ::INTERVAL T;
  static T point( const double c ) { return T(c); }
  static T zeroone() { return T(0.,1.); }
  static void I(T& x, const T& y) { x = y; }
  static double l(const T& x) { return ::Inf(x); }
  static double u(const T& x) { return ::Sup(x); }
  static double abs (const T& x) { return ::Abs(x);  }
  static double mid (const T& x) { return ::Mid(x);  }
  static double diam(const T& x) { return ::Diam(x); }
  static T inv (const T& x) { return T(1.)/x;  }
  static T sqr (const T& x) { return ::Sqr(x);  }
  static T sqrt(const T& x) {
    //if( ::Inf(x) < 0. ) throw std::runtime_error("negative square root in PROFIL");
    return ::Sqrt(x);
  }
  static T exp (const T& x) { return ::Exp(x);  }
  static T log (const T& x) { return ::Log(x);  }
  static T xlog0(const T& x) { return ::Log(x)*x; }
  static T xlog (const T& x) {
    T zmin = ::Exp(-1);
    T mono = ::Hull( xlog0( ::Inf(x) ), xlog0( ::Sup(x) ) );
    return ::Inf(x)>=::Sup(zmin)||::Sup(x)<=::Inf(zmin)? ::Hull( mono, xlog0(zmin) ): mono; 
  }
  static T lmtd0(const T& x,const T& y){ return x==y? x: (x-y)/(::Log(x)-::Log(y)); }
  static T lmtd (const T& x,const T& y){
    T vmin = lmtd0( ::Inf(x), ::Inf(y) );
    T vmax = lmtd0( ::Sup(x), ::Sup(y) );
    return ::Hull( vmin, vmax );
  }
  static T rlmtd0(const T& x,const T& y){ return x==y? inv(x): (::Log(x)-::Log(y))/(x-y); }
  static T rlmtd (const T& x,const T& y){
    T vmin = rlmtd0( ::Sup(x), ::Sup(y) );
    T vmax = rlmtd0( ::Inf(x), ::Inf(y) );
    return ::Hull( vmin, vmax );
  }
  static T fabs(const T& x) { return T(::Pred(mc::mid(::Inf(x),::Sup(x),0.)),::Succ(::Abs(x))); }
  static T sin (const T& x) { return ::Sin(x);  }
  static T cos (const T& x) { return ::Cos(x);  }
  static T tan (const T& x) { return ::Tan(x);  }
  static T asin(const T& x) { return ::ArcSin(x); }
  static T acos(const T& x) { return ::ArcCos(x); }
  static T atan(const T& x) { return ::ArcTan(x); }
  static T sinh(const T& x) { return ::Sinh(x); }
  static T cosh(const T& x) { return ::Cosh(x); }
  static T tanh(const T& x) { return ::Tanh(x); }
  static T erf (const T& x) {
    return T( ::Pred(std::erf( ::Inf(x) )), ::Succ(std::erf( ::Sup(x) )) );
  }
  static T erfc(const T& x) {
    return T( ::Pred(std::erfc( ::Sup(x) )), ::Succ(std::erfc( ::Inf(x) )) );
  }
  static T fstep(const T& x) {
    if( ::Sup(x) < 0 ) return T(0); 
    if( ::Inf(x) > 0 ) return T(1); 
    return T(0,1);
  }
  static T bstep(const T& x) {
    if( ::Inf(x) > 0 ) return T(0); 
    if( ::Sup(x) < 0 ) return T(1); 
    return T(0,1);
  }
  static T hull(const T& x, const T& y) { return ::Hull(x,y); }
  static T min (const T& x, const T& y) {
    return T( ::Pred(std::min(::Inf(x),::Inf(y))), ::Succ(std::min(::Sup(x),::Sup(y))) );
  }
  static T max (const T& x, const T& y) {
    return T( ::Pred(std::max(::Inf(x),::Inf(y))), ::Succ(std::max(::Sup(x),::Sup(y))) );
  }
  static T arh (const T& x, const double k) { return ::Exp( -x / k ); }
  static T pow(const T& x, const int n) { 
    return (n>=3&&n%2)? T(::Power(Inf(x),n),::Power(Sup(x),n)): ::Power(x,n);
  }
  template <typename X, typename Y> static T pow(const X& x, const Y& y) { return ::Power(x,y); }
  static T cheb (const T& x, const unsigned n) { return T(-1.,1.); }
  static T prod (const unsigned int n, const T* x) {
    return n? x[0] * prod(n-1, x+1): 1.;
  }
  static T monom (const unsigned int n, const T* x, const unsigned* k) {
    return n? ::Power(x[0], k[0]) * monom(n-1, x+1, k+1): 1.;
  }
  static bool inter(T& xIy, const T& x, const T& y) { return ::Intersection(xIy,x,y); }
  static bool eq(const T& x, const T& y) { return x==y; }
  static bool ne(const T& x, const T& y) { return x!=y; }
  static bool lt(const T& x, const T& y) { return x<y;  }
  static bool le(const T& x, const T& y) { return x<=y; }
  static bool gt(const T& x, const T& y) { return y<x;  }
  static bool ge(const T& x, const T& y) { return y<=x; }
};

} // namespace mc

#endif
