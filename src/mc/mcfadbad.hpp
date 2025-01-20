// Copyright (C) Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the Eclipse Public License.

#ifndef MC__MCFADBAD_HPP
#define MC__MCFADBAD_HPP

#include "mcop.hpp"
#include "fadiff.h"

namespace fadbad
{

template <typename T, unsigned int N>
INLINE2 FTypeName<T,N> pow2(const FTypeName<T,N>& a, const int b)
{
	FTypeName<T,N> c(Op<T>::myPow(a.val(),b));
	if (!a.depend()) return c;
	T tmp(Op<T>::myPow(a.val(),b-1)*Op<T>::myInteger(b));
	c.setDepend(a);
	for(unsigned int i=0;i<N;++i) c[i]=tmp*a[i];
	return c;
}

template <typename T >
INLINE2 FTypeName<T,0> pow2(const FTypeName<T,0>& a, const int b)
{
	FTypeName<T,0> c(Op<T>::myPow(a.val(),b));
	if (!a.depend()) return c;
	T tmp(Op<T>::myPow(a.val(),b-1)*Op<T>::myInteger(b));
	c.setDepend(a);
	for(unsigned int i=0;i<c.size();++i) c[i]=tmp*a[i];
	return c;
}
/*
template <typename T>
INLINE2 T cheb(const T& a, const unsigned b)
{
  switch( b ){
    case 0: return Op<T>::myOne();
    case 1: return a;
    default: return Op<T>::myTwo() * a * cheb(a,b-1) - cheb(a,b-2); }
}
*/
template <typename T, unsigned int N>
INLINE2 FTypeName<T,N> cheb(const FTypeName<T,N>& a, const unsigned b)
{
	FTypeName<T,N> c(mc::Op<T>::cheb(a.val(),b));
	if (!a.depend()) return c;

	T tmp( b%2? 0.5: 0. );
        for( int j=b-1; j>0; j-=2 )
          tmp += mc::Op<T>::cheb(a.val(),j);
        tmp *= 2*(double)b;

	c.setDepend(a);
	for(unsigned int i=0;i<N;++i) c[i]=tmp*a[i];
	return c;
}

template <typename T>
INLINE2 FTypeName<T,0> cheb(const FTypeName<T,0>& a, const unsigned b)
{
	FTypeName<T,0> c(mc::Op<T>::cheb(a.val(),b));
	if (!a.depend()) return c;

	T tmp( b%2? 0.5: 0. );
        for( int j=b-1; j>0; j-=2 )
          tmp += mc::Op<T>::cheb(a.val(),j);
        tmp *= 2*(double)b;

	c.setDepend(a);
	for(unsigned int i=0;i<c.size();++i) c[i]=tmp*a[i];
	return c;
}
/*
//@AVT.SVT: 01.06.2017
template <typename T, unsigned int N>
INLINE2 FTypeName<T,N> xlog (const FTypeName<T,N>& a)
{
    FTypeName<T,N> c(a*log(a));
    return c;
}
template <typename T>
INLINE2 FTypeName<T,0> xlog (const FTypeName<T,0>& a)
{
	FTypeName<T,0> c(a*log(a));
	return c;
}
*/
//@AVT.SVT: 07.06.2017
template <typename T, unsigned int N>
INLINE2 FTypeName<T,N> lmtd (const FTypeName<T,N>& a, const FTypeName<T,N>& b)
{   
    if(Op<T>::myEq(a.val(), b.val())){
      FTypeName<T,N> c(a.val());
      c.setDepend(a,b);
      for(unsigned int i=0;i<N;++i) c[i]=0.5*a[i]+0.5*b[i];
      return c;
    }
    FTypeName<T,N> c((a-b)/(mc::Op<T>::log(a)-mc::Op<T>::log(b)));
    return c;
}
template <typename T>
INLINE2 FTypeName<T,0> lmtd (const FTypeName<T,0>& a, const FTypeName<T,0>& b)
{
    if(Op<T>::myEq(a.val(), b.val())){
      FTypeName<T,0> c(a.val());
      c.setDepend(a,b);
      for(unsigned int i=0;i<c.size();++i) c[i]=0.5*a[i]+0.5*b[i];
      return c;
    }
    FTypeName<T,0> c((a-b)/(mc::Op<T>::log(a)-mc::Op<T>::log(b)));
    return c;
}
//@AVT.SVT: 08.06.2017
template <typename T, unsigned int N>
INLINE2 FTypeName<T,N> rlmtd (const FTypeName<T,N>& a, const FTypeName<T,N>& b)
{   
    if(Op<T>::myEq(a.val(), b.val())){
      FTypeName<T,N> c(1./a.val());
      if (!a.depend()) return c;
      c.setDepend(a,b);
      for(unsigned int i=0;i<N;++i) c[i]=-a[i]/(2.*mc::Op<T>::sqr(a.val()))-b[i]/(2.*mc::Op<T>::sqr(b.val()));
      return c;
    }
    FTypeName<T,N> c((mc::Op<T>::log(a)-mc::Op<T>::log(b))/(a-b));
    return c;
}
template <typename T>
INLINE2 FTypeName<T,0> rlmtd (const FTypeName<T,0>& a, const FTypeName<T,0>& b)
{
    if(Op<T>::myEq(a.val(), b.val())){
      FTypeName<T,0> c(1./a.val());
      if (!a.depend()) return c;
      c.setDepend(a,b);
      for(unsigned int i=0;i<c.size();++i) c[i]=-a[i]/(2.*mc::Op<T>::sqr(a.val()))-b[i]/(2.*mc::Op<T>::sqr(b.val()));
      return c;
    }
    FTypeName<T,0> c((mc::Op<T>::log(a)-mc::Op<T>::log(b))/(a-b));
    return c;
}
//@ICL: 04.01.2024
template <typename T, unsigned int N>
INLINE2 FTypeName<T,N> fabs(const FTypeName<T,N>& a)
{
	FTypeName<T,N> c(mc::Op<T>::fabs(a.val()));
	if (!a.depend()) return c;
	T tmp(mc::Op<T>::fstep(a.val())*Op<T>::myInteger(2)-Op<T>::myInteger(1));
	c.setDepend(a);
	for(unsigned int i=0;i<N;++i) c[i]=tmp*a[i];
	return c;
}

template <typename T >
INLINE2 FTypeName<T,0> fabs(const FTypeName<T,0>& a)
{
	FTypeName<T,0> c(mc::Op<T>::fabs(a.val()));
	if (!a.depend()) return c;
	T tmp(mc::Op<T>::fstep(a.val())*Op<T>::myInteger(2)-Op<T>::myInteger(1));
	c.setDepend(a);
	for(unsigned int i=0;i<c.size();++i) c[i]=tmp*a[i];
	return c;
}
template <typename T, unsigned int N>
INLINE2 FTypeName<T,N> min (const FTypeName<T,N>& a, const FTypeName<T,N>& b)
{
    //std::cout << "fadbad::min\n";
    FTypeName<T,N> c(0.5*(a+b-mc::Op<T>::fabs(a-b)));
    return c;
}
template <typename T>
INLINE2 FTypeName<T,0> min (const FTypeName<T,0>& a, const FTypeName<T,0>& b)
{
    //std::cout << "fadbad::min\n";
    FTypeName<T,0> c(0.5*(a+b-mc::Op<T>::fabs(a-b)));
    return c;
}
template <typename T, unsigned int N>
INLINE2 FTypeName<T,N> max (const FTypeName<T,N>& a, const FTypeName<T,N>& b)
{
    //std::cout << "fadbad::max\n";
    FTypeName<T,N> c(0.5*(a+b+mc::Op<T>::fabs(a-b)));
    return c;
}
template <typename T>
INLINE2 FTypeName<T,0> max (const FTypeName<T,0>& a, const FTypeName<T,0>& b)
{
    //std::cout << "fadbad::max\n";
    FTypeName<T,0> c(0.5*(a+b+mc::Op<T>::fabs(a-b)));
    return c;
}

} // end namespace fadbad

#include "tadiff.h"

namespace fadbad
{

template <typename U, int N>
struct TTypeNamePOW3 : public UnTTypeNameHV<U,N>
{
	int m_b;
	TTypeNamePOW3(const U& val, TTypeNameHV<U,N>* pOp, const int b):UnTTypeNameHV<U,N>(val,pOp),m_b(b){}
	TTypeNamePOW3(TTypeNameHV<U,N>* pOp, const int b):UnTTypeNameHV<U,N>(pOp),m_b(b){}
	unsigned int eval(const unsigned int k)
	{
		unsigned int l=this->opEval(k);
                if( m_b==0 ){
 			if (0==this->length()) { this->val(0)=Op<U>::myOne(); this->length()=1; }
			for(unsigned int i=this->length();i<l;++i) { this->val(i)=Op<U>::myZero(); }
                }
                else if( m_b==1 ){
			for(unsigned int i=this->length();i<l;++i) { this->val(i)=this->opVal(i); }
		}
                else if( m_b==2 ){
			if (0==this->length()) { this->val(0)=Op<U>::mySqr(this->opVal(0)); this->length()=1; }
			for(unsigned int i=this->length();i<l;++i)
			{
				this->val(i)=Op<U>::myZero();
				unsigned int m=(i+1)/2;
				for(unsigned int j=0;j<m;++j) Op<U>::myCadd(this->val(i), this->opVal(i-j)*this->opVal(j));
				Op<U>::myCmul(this->val(i), Op<U>::myTwo());
				if (0==i%2) Op<U>::myCadd(this->val(i), Op<U>::mySqr(this->opVal(m)));
			}
		}
		else if( m_b==3 ){
			if (0==this->length()) { this->val(0)=Op<U>::myPow(this->opVal(0),m_b); this->length()=1; }
			if (1<l && 1==this->length() ) { this->val(1)=Op<U>::myPow(this->opVal(0),m_b-1)
				*this->opVal(1)*Op<U>::myInteger(m_b); this->length()=2; }
			if (2<l && 2==this->length() ) { this->val(2)=Op<U>::myPow(this->opVal(0),m_b-2)
				*( this->opVal(0)*this->opVal(2) + Op<U>::myInteger(m_b-1)*Op<U>::mySqr(this->opVal(1)) )
				*Op<U>::myInteger(m_b); this->length()=3; }
			for(unsigned int i=this->length();i<l;++i)
                        {
                                this->val(i)=Op<U>::myZero();
				unsigned int m=(i+1)/2;
				for(unsigned int j=0;j<m;++j) Op<U>::myCadd(this->val(i), this->opVal(i-j)*this->opVal(j));
				Op<U>::myCmul(this->val(i), Op<U>::myTwo());
				if (0==i%2) Op<U>::myCadd(this->val(i), Op<U>::mySqr(this->opVal(m)));                        
                        }
			for(unsigned int i=l-1; i>=this->length();--i)
			{
				Op<U>::myCmul(this->val(i), this->opVal(0));
				for(unsigned int j=1;j<=i;++j) Op<U>::myCadd(this->val(i), this->val(i-j)*this->opVal(j));
			}
		}
                else{                       
			if (0==this->length()) { this->val(0)=Op<U>::myPow(this->opVal(0),m_b); this->length()=1; }
			for(unsigned int i=this->length();i<l;++i)
			{
				this->val(i)=Op<U>::myZero();
				for(unsigned int j=0;j<i;++j)
					Op<U>::myCadd(this->val(i), ( m_b - (m_b+Op<U>::myOne()) * Op<U>::myInteger(j) /
						Op<U>::myInteger(i) )*this->opVal(i-j)*this->val(j));
			}
                }
		return this->length()=l;
	}
private:
	void operator=(const TTypeNamePOW3<U,N>&){} // not allowed
};

template <typename U, int N>
TTypeName<U,N> pow(const TTypeName<U,N>& val, const int b)
{
	TTypeNameHV<U,N>* pHV=val.length()>0 ?
		new TTypeNamePOW3<U,N>(Op<U>::myPow(val.val(),b), val.getTTypeNameHV(), b):
		new TTypeNamePOW3<U,N>(val.getTTypeNameHV(), b);
	return TTypeName<U,N>(pHV);
}

} // end namespace fadbad

#include "badiff.h"

namespace mc
{

//! @brief C++ structure for specialization of the mc::Op templated structure for use of the FADBAD type fadbad::F inside other MC++ types
template< typename U > struct Op< fadbad::F<U> >
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
  static TU exp (const TU& x) { return fadbad::exp(x);  }
  static TU log (const TU& x) { return fadbad::log(x);  }
  static TU xlog(const TU& x) { return x*fadbad::log(x); }
  static TU lmtd(const TU& x, const TU& y) { return fadbad::lmtd(x,y); }
  static TU rlmtd(const TU& x, const TU& y) { return fadbad::rlmtd(x,y); }
  static TU fabs(const TU& x) { return fadbad::fabs(x); }
  static TU sin (const TU& x) { return fadbad::sin(x);  }
  static TU cos (const TU& x) { return fadbad::cos(x);  }
  static TU tan (const TU& x) { return fadbad::tan(x);  }
  static TU asin(const TU& x) { return fadbad::asin(x); }
  static TU acos(const TU& x) { return fadbad::acos(x); }
  static TU atan(const TU& x) { return fadbad::atan(x); }
  static TU sinh(const TU& x) { return fadbad::sinh(x); }
  static TU cosh(const TU& x) { return fadbad::cosh(x); }
  static TU tanh(const TU& x) { return fadbad::tanh(x); }
  static TU erf (const TU& x) { throw std::runtime_error("mc::Op<fadbad::F<U>>::erf -- operation not permitted"); }
  static TU erfc(const TU& x) { throw std::runtime_error("mc::Op<fadbad::F<U>>::erfc -- operation not permitted"); }
  static TU fstep(const TU& x) { throw std::runtime_error("mc::Op<fadbad::F<U>>::fstep -- operation not permitted"); }
  static TU bstep(const TU& x) { throw std::runtime_error("mc::Op<fadbad::F<U>>::bstep -- operation not permitted"); }
  static TU hull(const TU& x, const TU& y) { throw std::runtime_error("mc::Op<fadbad::F<U>>::hull -- operation not permitted"); }
  static TU min (const TU& x, const TU& y) { return 0.5*(x+y-fadbad::fabs(x-y)); }
  static TU max (const TU& x, const TU& y) { return 0.5*(x+y+fadbad::fabs(x-y)); }
  static TU arh (const TU& x, const double k) { return fadbad::exp(-k/x); }
  template <typename X, typename Y> static TU pow(const X& x, const Y& y) { return fadbad::pow(x,y); }
  static TU cheb(const TU& x, const unsigned n) { return fadbad::cheb(x,n); }
  static TU prod (const unsigned n, const TU* x) { switch( n ){ case 0: return 1.; case 1: return x[0]; default: return x[0]*prod(n-1,x+1); } }
  static TU monom (const unsigned n, const TU* x, const unsigned* k) { switch( n ){ case 0: return 1.; case 1: return pow(x[0],(int)k[0]); default: return pow(x[0],(int)k[0])*monom(n-1,x+1,k+1); } }
  static bool inter(TU& xIy, const TU& x, const TU& y) { xIy = x; return true; }
  static bool eq(const TU& x, const TU& y) { throw std::runtime_error("mc::Op<fadbad::F<U>>::eq -- operation not permitted"); }
  static bool ne(const TU& x, const TU& y) { throw std::runtime_error("mc::Op<fadbad::F<U>>::ne -- operation not permitted"); }
  static bool lt(const TU& x, const TU& y) { throw std::runtime_error("mc::Op<fadbad::F<U>>::lt -- operation not permitted"); }
  static bool le(const TU& x, const TU& y) { throw std::runtime_error("mc::Op<fadbad::F<U>>::le -- operation not permitted"); }
  static bool gt(const TU& x, const TU& y) { throw std::runtime_error("mc::Op<fadbad::F<U>>::gt -- operation not permitted"); }
  static bool ge(const TU& x, const TU& y) { throw std::runtime_error("mc::Op<fadbad::F<U>>::ge -- operation not permitted"); }
};

//! @brief C++ structure for specialization of the mc::Op templated structure for use of the FADBAD type fadbad::B inside other MC++ types
template< typename U > struct Op< fadbad::B<U> >
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
  static TU exp (const TU& x) { return fadbad::exp(x);  }
  static TU log (const TU& x) { return fadbad::log(x);  }
  static TU xlog(const TU& x) { return x*fadbad::log(x); }
  static TU lmtd(const TU& x, const TU& y) { throw std::runtime_error("mc::Op<fadbad::B<U>>::lmtd -- operation not permitted"); }
  static TU rlmtd(const TU& x, const TU& y) { throw std::runtime_error("mc::Op<fadbad::B<U>>::rlmtd -- operation not permitted"); }
  static TU fabs(const TU& x) { throw std::runtime_error("mc::Op<fadbad::B<U>>::fabs -- operation not permitted"); }
  static TU sin (const TU& x) { return fadbad::sin(x);  }
  static TU cos (const TU& x) { return fadbad::cos(x);  }
  static TU tan (const TU& x) { return fadbad::tan(x);  }
  static TU asin(const TU& x) { return fadbad::asin(x); }
  static TU acos(const TU& x) { return fadbad::acos(x); }
  static TU atan(const TU& x) { return fadbad::atan(x); }
  static TU sinh(const TU& x) { return fadbad::sinh(x); }
  static TU cosh(const TU& x) { return fadbad::cosh(x); }
  static TU tanh(const TU& x) { return fadbad::tanh(x); }
  static TU erf (const TU& x) { throw std::runtime_error("mc::Op<fadbad::B<U>>::erf -- operation not permitted"); }
  static TU erfc(const TU& x) { throw std::runtime_error("mc::Op<fadbad::B<U>>::erfc -- operation not permitted"); }
  static TU fstep(const TU& x) { throw std::runtime_error("mc::Op<fadbad::B<U>>::fstep -- operation not permitted"); }
  static TU bstep(const TU& x) { throw std::runtime_error("mc::Op<fadbad::B<U>>::bstep -- operation not permitted"); }
  static TU hull(const TU& x, const TU& y) { throw std::runtime_error("mc::Op<fadbad::B<U>>::hull -- operation not permitted"); }
  static TU min (const TU& x, const TU& y) { throw std::runtime_error("mc::Op<fadbad::B<U>>::min -- operation not permitted"); }
  static TU max (const TU& x, const TU& y) { throw std::runtime_error("mc::Op<fadbad::B<U>>::max -- operation not permitted"); }
  static TU arh (const TU& x, const double k) { return fadbad::exp(-k/x); }
  template <typename X, typename Y> static TU pow(const X& x, const Y& y) { return fadbad::pow(x,y); }
  static TU cheb(const TU& x, const unsigned n) { throw std::runtime_error("mc::Op<fadbad::B<U>>::cheb -- operation not permitted"); }
  static TU prod (const unsigned n, const TU* x) { switch( n ){ case 0: return 1.; case 1: return x[0]; default: return x[0]*prod(n-1,x+1); } }
  static TU monom (const unsigned n, const TU* x, const unsigned* k) { switch( n ){ case 0: return 1.; case 1: return pow(x[0],(int)k[0]); default: return pow(x[0],(int)k[0])*monom(n-1,x+1,k+1); } }
  static bool inter(TU& xIy, const TU& x, const TU& y) { xIy = x; return true; }
  static bool eq(const TU& x, const TU& y) { throw std::runtime_error("mc::Op<fadbad::B<U>>::eq -- operation not permitted"); }
  static bool ne(const TU& x, const TU& y) { throw std::runtime_error("mc::Op<fadbad::B<U>>::ne -- operation not permitted"); }
  static bool lt(const TU& x, const TU& y) { throw std::runtime_error("mc::Op<fadbad::B<U>>::lt -- operation not permitted"); }
  static bool le(const TU& x, const TU& y) { throw std::runtime_error("mc::Op<fadbad::B<U>>::le -- operation not permitted"); }
  static bool gt(const TU& x, const TU& y) { throw std::runtime_error("mc::Op<fadbad::B<U>>::gt -- operation not permitted"); }
  static bool ge(const TU& x, const TU& y) { throw std::runtime_error("mc::Op<fadbad::B<U>>::ge -- operation not permitted"); }
};

//! @brief C++ structure for specialization of the mc::Op templated structure for use of the FADBAD type fadbad::T inside other MC++ types
template< typename U > struct Op< fadbad::T<U> >
{
  typedef fadbad::T<U> TU;
  static TU point( const double c ) { throw std::runtime_error("mc::Op<fadbad::T<U>>::point -- operation not permitted"); }
  static TU zeroone() { throw std::runtime_error("mc::Op<fadbad::T<U>>::zeroone -- operation not permitted"); }
  static void I(TU& x, const TU&y) { x = y; }
  static double l(const TU& x) { throw std::runtime_error("mc::Op<fadbad::T<U>>::l -- operation not permitted"); }
  static double u(const TU& x) { throw std::runtime_error("mc::Op<fadbad::T<U>>::u -- operation not permitted"); }
  static double abs (const TU& x) { throw std::runtime_error("mc::Op<fadbad::T<U>>::abs -- operation not permitted"); }
  static double mid (const TU& x) { throw std::runtime_error("mc::Op<fadbad::T<U>>::mid -- operation not permitted"); }
  static double diam(const TU& x) { throw std::runtime_error("mc::Op<fadbad::T<U>>::diam -- operation not permitted"); }
  static TU inv (const TU& x) { return 1./x;  }
  static TU sqr (const TU& x) { return fadbad::sqr(x);  }
  static TU sqrt(const TU& x) { return fadbad::sqrt(x); }
  static TU exp (const TU& x) { return fadbad::exp(x);  }
  static TU log (const TU& x) { return fadbad::log(x);  }
  static TU xlog(const TU& x) { return x*fadbad::log(x); }
  static TU lmtd(const TU& x, const TU& y) { throw std::runtime_error("mc::Op<fadbad::T<U>>::lmtd -- operation not permitted"); }  
  static TU rlmtd(const TU& x, const TU& y) { throw std::runtime_error("mc::Op<fadbad::T<U>>::rlmtd -- operation not permitted"); }
  static TU fabs(const TU& x) { throw std::runtime_error("mc::Op<fadbad::T<U>>::fabs -- operation not permitted"); }
  static TU sin (const TU& x) { return fadbad::sin(x);  }
  static TU cos (const TU& x) { return fadbad::cos(x);  }
  static TU tan (const TU& x) { return fadbad::tan(x);  }
  static TU asin(const TU& x) { return fadbad::asin(x); }
  static TU acos(const TU& x) { return fadbad::acos(x); }
  static TU atan(const TU& x) { return fadbad::atan(x); }
  static TU sinh(const TU& x) { return fadbad::sinh(x); }
  static TU cosh(const TU& x) { return fadbad::cosh(x); }
  static TU tanh(const TU& x) { return fadbad::tanh(x); }
  static TU erf (const TU& x) { throw std::runtime_error("mc::Op<fadbad::T<U>>::erf -- operation not permitted"); }
  static TU erfc(const TU& x) { throw std::runtime_error("mc::Op<fadbad::T<U>>::erfc -- operation not permitted"); }
  static TU fstep(const TU& x) { throw std::runtime_error("mc::Op<fadbad::T<U>>::fstep -- operation not permitted"); }
  static TU bstep(const TU& x) { throw std::runtime_error("mc::Op<fadbad::T<U>>::bstep -- operation not permitted"); }
  static TU hull(const TU& x, const TU& y) { throw std::runtime_error("mc::Op<fadbad::T<U>>::hull -- operation not permitted"); }
  static TU min (const TU& x, const TU& y) { throw std::runtime_error("mc::Op<fadbad::T<U>>::min -- operation not permitted"); }
  static TU max (const TU& x, const TU& y) { throw std::runtime_error("mc::Op<fadbad::T<U>>::max -- operation not permitted"); }
  static TU arh (const TU& x, const double k) { return fadbad::exp(-k/x); }
  template <typename X, typename Y> static TU pow(const X& x, const Y& y) { return fadbad::pow(x,y); }
  static TU cheb(const TU& x, const unsigned n) { throw std::runtime_error("mc::Op<fadbad::T<U>>::cheb -- operation not permitted"); }
  static TU prod (const unsigned n, const TU* x) { switch( n ){ case 0: return 1.; case 1: return x[0]; default: return x[0]*prod(n-1,x+1); } }
  static TU monom (const unsigned n, const TU* x, const unsigned* k) { switch( n ){ case 0: return 1.; case 1: return pow(x[0],(int)k[0]); default: return pow(x[0],(int)k[0])*monom(n-1,x+1,k+1); } }
  static bool inter(TU& xIy, const TU& x, const TU& y) { xIy = x; return true; }
  static bool eq(const TU& x, const TU& y) { throw std::runtime_error("mc::Op<fadbad::T<U>>::eq -- operation not permitted"); }
  static bool ne(const TU& x, const TU& y) { throw std::runtime_error("mc::Op<fadbad::T<U>>::ne -- operation not permitted"); }
  static bool lt(const TU& x, const TU& y) { throw std::runtime_error("mc::Op<fadbad::T<U>>::lt -- operation not permitted"); }
  static bool le(const TU& x, const TU& y) { throw std::runtime_error("mc::Op<fadbad::T<U>>::le -- operation not permitted"); }
  static bool gt(const TU& x, const TU& y) { throw std::runtime_error("mc::Op<fadbad::T<U>>::gt -- operation not permitted"); }
  static bool ge(const TU& x, const TU& y) { throw std::runtime_error("mc::Op<fadbad::T<U>>::ge -- operation not permitted"); }
};

} // namespace mc

#endif
