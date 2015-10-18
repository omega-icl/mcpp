// Copyright (C) 2009-2013 Benoit Chachuat, Imperial College London.
// All Rights Reserved.

#ifndef MC__MCFADIFF_HPP
#define MC__MCFADIFF_HPP

#include "fadbad.h"
#include "fadiff.h"

namespace fadbad
{
template <typename T, unsigned int N>
INLINE2 FTypeName<T,N> pow2(const FTypeName<T,N>& a, const int b)
{
	FTypeName<T,N> c(Op<T>::myPow(a.val(),b));
	if (!a.depend()) return c;
	T tmp(b*Op<T>::myPow(a.val(),b-1));
	c.setDepend(a);
	for(unsigned int i=0;i<N;++i) c[i]=tmp*a[i];
	return c;
}
template <typename T >
INLINE2 FTypeName<T,0> pow2(const FTypeName<T,0>& a, const int b)
{
	FTypeName<T,0> c(Op<T>::myPow(a.val(),b));
	if (!a.depend()) return c;
	T tmp(b*Op<T>::myPow(a.val(),b-1));
	c.setDepend(a);
	for(unsigned int i=0;i<c.size();++i) c[i]=tmp*a[i];
	return c;
}
template <typename T, unsigned int N>
INLINE2 FTypeName<T,N> cheb(const FTypeName<T,N>& a, const unsigned b)
{
	FTypeName<T,N> c(Op<T>::myCheb(a.val(),b));
	if (!a.depend()) return c;
        // dTn/dx = n*Un-1(x)
        // Un-1(x) = 2*(T1(x)+T3(x)+...+Tn-1(x)) if n even
        //           2*(T0(x)+T2(x)+...+Tn-1(x))-1 if n odd
	T tmp(0.);
        if( b%2 ){ // odd case
          for( unsigned k=0; k<b; k+=2 ) tmp += Op<T>::myCheb(a.val(),k);
          tmp *= 2.; tmp -= 1.;
        }
        else{ // even case
          for( unsigned k=1; k<b; k+=2 ) tmp += Op<T>::myCheb(a.val(),k);
          tmp *= 2.;
        }
	c.setDepend(a);
	for(unsigned int i=0;i<N;++i) c[i]=tmp*a[i];
	return c;
}
template <typename T>
INLINE2 FTypeName<T,0> cheb(const FTypeName<T,0>& a, const unsigned b)
{
	FTypeName<T,0> c(Op<T>::myCheb(a.val(),b));
	if (!a.depend()) return c;
        // dTn/dx = n*Un-1(x)
        // Un-1(x) = 2*(T1(x)+T3(x)+...+Tn-1(x)) if n even
        //           2*(T0(x)+T2(x)+...+Tn-1(x))-1 if n odd
	T tmp(0.);
        if( b%2 ){ // odd case
          for( unsigned k=0; k<b; k+=2 ) tmp += Op<T>::myCheb(a.val(),k);
          tmp *= 2.; tmp -= 1.;
        }
        else{ // even case
          for( unsigned k=1; k<b; k+=2 ) tmp += Op<T>::myCheb(a.val(),k);
          tmp *= 2.;
        }
	c.setDepend(a);
	for(unsigned int i=0;i<c.size();++i) c[i]=tmp*a[i];
	return c;
}
} // end namespace fadbad

#include "mcop.hpp"

namespace mc
{

//! @brief C++ structure for specialization of the mc::Op templated structure to allow usage of the FADBAD type fadbad::F inside other MC++ type, e.g. mc::McCormick
template <> template<typename U> struct Op< fadbad::F<U> >
{
  typedef fadbad::F<U> TU;
  static TU point( const double c ) { throw std::runtime_error("mc::Op<fadbad::F<U>>::point -- operation not permitted"); }
  static TU zeroone() { throw std::runtime_error("mc::Op<fadbad::F<U>>::zeroone -- operation not permitted"); }
  static void I(TU& x, const TU&y) { x = y; }
  static double l(const TU& x) { throw std::runtime_error("mc::Op<fadbad::F<U>>::l -- operation not permitted"); }
  static double u(const TU& x) { throw std::runtime_error("mc::Op<fadbad::F<U>>::u -- operation not permitted"); }
  static double abs (const TU& x) { throw std::runtime_error("mc::Op<fadbad::F<U>>::abs -- operation not permitted"); }
  static double mid (const TU& x) { throw std::runtime_error("mc::Op<fadbad::F<U>>::mid -- operation not permitted"); }
  static double diam(const TU& x) { throw std::runtime_error("mc::Op<fadbad::F<U>>::diam -- operation not permitted"); }
  static TU inv (const TU& x) { return 1./x;  }
  static TU sqr (const TU& x) { return fadbad::sqr(x);  }
  static TU sqrt(const TU& x) { return fadbad::sqrt(x); }
  static TU log (const TU& x) { return fadbad::log(x);  }
  static TU xlog(const TU& x) { return x*fadbad::log(x); }
  static TU fabs(const TU& x) { throw std::runtime_error("mc::Op<fadbad::F<U>>::fabs -- operation not permitted"); }
  static TU exp (const TU& x) { return fadbad::exp(x);  }
  static TU sin (const TU& x) { return fadbad::sin(x);  }
  static TU cos (const TU& x) { return fadbad::cos(x);  }
  static TU tan (const TU& x) { return fadbad::tan(x);  }
  static TU asin(const TU& x) { return fadbad::asin(x); }
  static TU acos(const TU& x) { return fadbad::acos(x); }
  static TU atan(const TU& x) { return fadbad::atan(x); }
  static TU erf (const TU& x) { throw std::runtime_error("mc::Op<fadbad::F<U>>::erf -- operation not permitted"); }
  static TU erfc(const TU& x) { throw std::runtime_error("mc::Op<fadbad::F<U>>::erfc -- operation not permitted"); }
  static TU fstep(const TU& x) { throw std::runtime_error("mc::Op<fadbad::F<U>>::fstep -- operation not permitted"); }
  static TU bstep(const TU& x) { throw std::runtime_error("mc::Op<fadbad::F<U>>::bstep -- operation not permitted"); }
  static TU hull(const TU& x, const TU& y) { throw std::runtime_error("mc::Op<fadbad::F<U>>::hull -- operation not permitted"); }
  static TU min (const TU& x, const TU& y) { throw std::runtime_error("mc::Op<fadbad::F<U>>::min -- operation not permitted"); }
  static TU max (const TU& x, const TU& y) { throw std::runtime_error("mc::Op<fadbad::F<U>>::max -- operation not permitted"); }
  static TU arh (const TU& x, const double k) { return fadbad::exp(-k/x); }
  static TU cheb(const TU& x, const unsigned n) { return fadbad::cheb(x,n); }
  template <typename X, typename Y> static TU pow(const X& x, const Y& y) { return fadbad::pow(x,y); }
  static TU monomial (const unsigned int n, const U* x, const int* k) { throw std::runtime_error("mc::Op<T<U>>::monomial -- operation not permitted"); }
  static bool inter(TU& xIy, const TU& x, const TU& y) { throw std::runtime_error("mc::Op<fadbad::F<U>>::inter -- operation not permitted"); }
  static bool eq(const TU& x, const TU& y) { throw std::runtime_error("mc::Op<fadbad::F<U>>::eq -- operation not permitted"); }
  static bool ne(const TU& x, const TU& y) { throw std::runtime_error("mc::Op<fadbad::F<U>>::ne -- operation not permitted"); }
  static bool lt(const TU& x, const TU& y) { throw std::runtime_error("mc::Op<fadbad::F<U>>::lt -- operation not permitted"); }
  static bool le(const TU& x, const TU& y) { throw std::runtime_error("mc::Op<fadbad::F<U>>::le -- operation not permitted"); }
  static bool gt(const TU& x, const TU& y) { throw std::runtime_error("mc::Op<fadbad::F<U>>::gt -- operation not permitted"); }
  static bool ge(const TU& x, const TU& y) { throw std::runtime_error("mc::Op<fadbad::F<U>>::ge -- operation not permitted"); }
};

} // namespace mc

#endif
