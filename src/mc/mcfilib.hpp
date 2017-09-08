// Copyright (C) 2009-2017 Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the Eclipse Public License.

#ifndef MC__MCFILIB_HPP
#define MC__MCFILIB_HPP

#include "mcop.hpp"
#include "mcfunc.hpp"
#include "interval/interval.hpp"

namespace mc
{

//! @brief Specialization of the structure mc::Op for use of the type filib::interval<double> of <A href="http://www.math.uni-wuppertal.de/~xsc/software/filib.html">FILIB++</A> as a template parameter in other MC++ types
template <> struct Op< filib::interval<double> >
{
  typedef filib::interval<double> T;
  static T point( const double c ) { return T(c); }
  static T zeroone() { return T(0.,1.); }
  static void I(T& x, const T& y) { x = y; }
  static double l(const T& x) { return filib::inf(x); }
  static double u(const T& x) { return filib::sup(x); }
  static double abs (const T& x) { return filib::mag(x);  }
  static double mid (const T& x) { return filib::mid(x);  }
  static double diam(const T& x) { return filib::diam(x); }
  static T inv (const T& x) { return T(1.)/x;  }
  static T sqr (const T& x) { return filib::sqr(x);  }
  static T sqrt(const T& x) { return filib::sqrt(x); }
  static T exp (const T& x) { return filib::exp(x);  }
  static T log (const T& x) { return filib::log(x);  }
  static T xlog(const T& x) { return filib::log(x)*x; }
  // !!THE RESULT IS NOT VERIFIED!!
  static T lmtd(const T& x,const T& y) { return T(mc::lmtd(filib::inf(x),filib::inf(y)),mc::lmtd(filib::sup(x),filib::sup(y)) ) ; }
  static T rlmtd(const T& x,const T& y) { return T(mc::rlmtd(filib::sup(x),filib::sup(y)),mc::rlmtd(filib::inf(x),filib::inf(y)) ) ; }
  static T fabs(const T& x) { return filib::abs(x); }
  static T sin (const T& x) { return filib::sin(x);  }
  static T cos (const T& x) { return filib::cos(x);  }
  static T tan (const T& x) { return filib::tan(x);  }
  static T asin(const T& x) { return filib::asin(x); }
  static T acos(const T& x) { return filib::acos(x); }
  static T atan(const T& x) { return filib::atan(x); }
  static T sinh(const T& x) { return filib::sinh(x); }
  static T cosh(const T& x) { return filib::cosh(x); }
  static T tanh(const T& x) { return filib::tanh(x); }
  static T erf (const T& x) { throw std::runtime_error("operation not permitted"); }
  static T erfc(const T& x) { throw std::runtime_error("operation not permitted"); }
  static T fstep(const T& x) { throw std::runtime_error("operation not permitted"); }
  static T bstep(const T& x) { throw std::runtime_error("operation not permitted"); }
  static T hull(const T& x, const T& y) { return x.hull(y); }
  static T min (const T& x, const T& y) { return x.imin(y); }
  static T max (const T& x, const T& y) { return x.imax(y); }
  static T arh (const T& x, const double k) { return filib::exp(-x/k); }
  static T cheb (const T& x, const unsigned n) { return T(-1.,1.); }
  template <typename X> static T pow(const X& x, const int n) { return filib::power(x,n); }
  template <typename X, typename Y> static T pow(const X& x, const Y& y) { return filib::pow(x,y); }
  static T monom (const unsigned int n, const T* x, const unsigned* k) { return n? filib::power(x[0], k[0]) * monom(n-1, x+1, k+1): 1.; }
  static bool inter(T& xIy, const T& x, const T& y) { xIy = x.intersect(y); return true; }  
  static bool eq(const T& x, const T& y) { return x.seq(y); }
  static bool ne(const T& x, const T& y) { return x.sne(y); }
  static bool lt(const T& x, const T& y) { return x.slt(y); }
  static bool le(const T& x, const T& y) { return x.sle(y); }
  static bool gt(const T& x, const T& y) { return x.sgt(y); }
  static bool ge(const T& x, const T& y) { return x.sge(y); }
};

//! @brief Specialization of the structure mc::Op for use of the type filib::interval<double,filib::native_switched,filib::i_mode_extended> of <A href="http://www.math.uni-wuppertal.de/~xsc/software/filib.html">FILIB++</A> as a template parameter in other MC++ types
template <> struct Op< filib::interval<double,filib::native_switched,filib::i_mode_extended> >
{
  typedef filib::interval<double,filib::native_switched,filib::i_mode_extended> T;
  static T point( const double c ) { return T(c); }
  static T zeroone() { return T(0.,1.); }
  static void I(T& x, const T& y) { x = y; }
  static double l(const T& x) { return filib::inf(x); }
  static double u(const T& x) { return filib::sup(x); }
  static double abs (const T& x) { return filib::mag(x);  }
  static double mid (const T& x) { return filib::mid(x);  }
  static double diam(const T& x) { return filib::diam(x); }
  static T inv (const T& x) { return T(1.)/x;  }
  static T sqr (const T& x) { return filib::sqr(x);  }
  static T sqrt(const T& x) { return filib::sqrt(x); }
  static T exp (const T& x) { return filib::exp(x);  }
  static T log (const T& x) { return filib::log(x);  }
  static T xlog(const T& x) { return filib::log(x)*x; }
  // !!THE RESULT IS NOT VERIFIED!!
  static T lmtd(const T& x,const T& y) { return T(mc::lmtd(filib::inf(x),filib::inf(y)),mc::lmtd(filib::sup(x),filib::sup(y)) ) ; }
  static T rlmtd(const T& x,const T& y) { return T(mc::rlmtd(filib::sup(x),filib::sup(y)),mc::rlmtd(filib::inf(x),filib::inf(y)) ) ; }
  static T fabs(const T& x) { return filib::abs(x); }
  static T sin (const T& x) { return filib::sin(x);  }
  static T cos (const T& x) { return filib::cos(x);  }
  static T tan (const T& x) { return filib::tan(x);  }
  static T asin(const T& x) { return filib::asin(x); }
  static T acos(const T& x) { return filib::acos(x); }
  static T atan(const T& x) { return filib::atan(x); }
  static T sinh(const T& x) { return filib::sinh(x); }
  static T cosh(const T& x) { return filib::cosh(x); }
  static T tanh(const T& x) { return filib::tanh(x); }
  static T hull(const T& x, const T& y) { return x.hull(y); }
  static T min (const T& x, const T& y) { return x.imin(y); }
  static T max (const T& x, const T& y) { return x.imax(y); }
  static T arh (const T& x, const double k) { return filib::exp(-x/k); }
  static T cheb (const T& x, const unsigned n) { return T(-1.,1.); }
  template <typename X> static T pow(const X& x, const int n) { return filib::power(x,n); }
  template <typename X, typename Y> static T pow(const X& x, const Y& y) { return filib::pow(x,y); }
  static T monom (const unsigned int n, const T* x, const unsigned* k) { return n? filib::power(x[0], k[0]) * monom(n-1, x+1, k+1): 1.; }
  static bool inter(T& xIy, const T& x, const T& y) { xIy = x.intersect(y); return !xIy.isEmpty(); }  
  static bool eq(const T& x, const T& y) { return x.seq(y); }
  static bool ne(const T& x, const T& y) { return x.sne(y); }
  static bool lt(const T& x, const T& y) { return x.slt(y); }
  static bool le(const T& x, const T& y) { return x.sle(y); }
  static bool gt(const T& x, const T& y) { return x.sgt(y); }
  static bool ge(const T& x, const T& y) { return x.sge(y); }
};

} // namespace mc

#endif
