// Copyright (C) 2009-2013 Benoit Chachuat, Imperial College London.
// All Rights Reserved.

#ifndef MC__MCBADIFF_HPP
#define MC__MCBADIFF_HPP

#include "fadbad.h"
#include "badiff.h"
#include "mcop.hpp"

namespace mc
{

//! @brief C++ structure for specialization of the mc::Op templated structure to allow usage of the FADBAD type fadbad::B inside other MC++ type, e.g. mc::McCormick
template <> template<typename U> struct Op< fadbad::B<U> >
{
  typedef fadbad::B<U> TU;
  static TU point( const double c ) { throw std::runtime_error("mc::Op<fadbad::B<U>>::point -- operation not permitted"); }
  static TU zeroone() { throw std::runtime_error("mc::Op<fadbad::B<U>>::zeroone -- operation not permitted"); }
  static void I(TU& x, const TU&y) { x = y; }
  static double l(const TU& x) { throw std::runtime_error("mc::Op<fadbad::B<U>>::l -- operation not permitted"); }
  static double u(const TU& x) { throw std::runtime_error("mc::Op<fadbad::B<U>>::u -- operation not permitted"); }
  static double abs (const TU& x) { throw std::runtime_error("mc::Op<fadbad::B<U>>::abs -- operation not permitted"); }
  static double mid (const TU& x) { throw std::runtime_error("mc::Op<fadbad::B<U>>::mid -- operation not permitted"); }
  static double diam(const TU& x) { throw std::runtime_error("mc::Op<fadbad::B<U>>::diam -- operation not permitted"); }
  static TU inv (const TU& x) { return 1./x;  }
  static TU sqr (const TU& x) { return fadbad::sqr(x);  }
  static TU sqrt(const TU& x) { return fadbad::sqrt(x); }
  static TU log (const TU& x) { return fadbad::log(x);  }
  static TU xlog(const TU& x) { return x*fadbad::log(x); }
  static TU fabs(const TU& x) { throw std::runtime_error("mc::Op<fadbad::B<U>>::fabs -- operation not permitted"); }
  static TU exp (const TU& x) { return fadbad::exp(x);  }
  static TU sin (const TU& x) { return fadbad::sin(x);  }
  static TU cos (const TU& x) { return fadbad::cos(x);  }
  static TU tan (const TU& x) { return fadbad::tan(x);  }
  static TU asin(const TU& x) { return fadbad::asin(x); }
  static TU acos(const TU& x) { return fadbad::acos(x); }
  static TU atan(const TU& x) { return fadbad::atan(x); }
  static TU erf (const TU& x) { throw std::runtime_error("mc::Op<fadbad::B<U>>::erf -- operation not permitted"); }
  static TU erfc(const TU& x) { throw std::runtime_error("mc::Op<fadbad::B<U>>::erfc -- operation not permitted"); }
  static TU fstep(const TU& x) { throw std::runtime_error("mc::Op<fadbad::B<U>>::fstep -- operation not permitted"); }
  static TU bstep(const TU& x) { throw std::runtime_error("mc::Op<fadbad::B<U>>::bstep -- operation not permitted"); }
  static TU hull(const TU& x, const TU& y) { throw std::runtime_error("mc::Op<fadbad::B<U>>::hull -- operation not permitted"); }
  static TU min (const TU& x, const TU& y) { throw std::runtime_error("mc::Op<fadbad::B<U>>::min -- operation not permitted"); }
  static TU max (const TU& x, const TU& y) { throw std::runtime_error("mc::Op<fadbad::B<U>>::max -- operation not permitted"); }
  static TU arh (const TU& x, const double k) { return fadbad::exp(-k/x); }
  static TU cheb(const TU& x, const unsigned n) { throw std::runtime_error("mc::Op<fadbad::B<U>>::max -- operation not permitted"); }
  template <typename X, typename Y> static TU pow(const X& x, const Y& y) { return fadbad::pow(x,y); }
  static TU monomial (const unsigned int n, const U* x, const int* k) { throw std::runtime_error("mc::Op<fadbad::B<U>>::monomial -- operation not permitted"); }
  static bool inter(TU& xIy, const TU& x, const TU& y) { throw std::runtime_error("mc::Op<fadbad::B<U>>::inter -- operation not permitted"); }
  static bool eq(const TU& x, const TU& y) { throw std::runtime_error("mc::Op<fadbad::B<U>>::eq -- operation not permitted"); }
  static bool ne(const TU& x, const TU& y) { throw std::runtime_error("mc::Op<fadbad::B<U>>::ne -- operation not permitted"); }
  static bool lt(const TU& x, const TU& y) { throw std::runtime_error("mc::Op<fadbad::B<U>>::lt -- operation not permitted"); }
  static bool le(const TU& x, const TU& y) { throw std::runtime_error("mc::Op<fadbad::B<U>>::le -- operation not permitted"); }
  static bool gt(const TU& x, const TU& y) { throw std::runtime_error("mc::Op<fadbad::B<U>>::gt -- operation not permitted"); }
  static bool ge(const TU& x, const TU& y) { throw std::runtime_error("mc::Op<fadbad::B<U>>::ge -- operation not permitted"); }
};

} // namespace mc

#endif
