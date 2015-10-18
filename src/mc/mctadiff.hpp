// Copyright (C) 2009-2013 Benoit Chachuat, Imperial College London.
// All Rights Reserved.

#ifndef MC__MCTADIFF_HPP
#define MC__MCTADIFF_HPP

#include "fadbad.h"
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

#include "mcop.hpp"

namespace mc
{

//! @brief C++ structure for specialization of the mc::Op templated structure to allow usage of the FADBAD type fadbad::T inside other MC++ type, e.g. mc::McCormick
template <> template<typename U> struct Op< fadbad::T<U> >
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
  static TU log (const TU& x) { return fadbad::log(x);  }
  static TU xlog(const TU& x) { return x*fadbad::log(x); }
  static TU fabs(const TU& x) { throw std::runtime_error("mc::Op<fadbad::T<U>>::fabs -- operation not permitted"); }
  static TU exp (const TU& x) { return fadbad::exp(x);  }
  static TU sin (const TU& x) { return fadbad::sin(x);  }
  static TU cos (const TU& x) { return fadbad::cos(x);  }
  static TU tan (const TU& x) { return fadbad::tan(x);  }
  static TU asin(const TU& x) { return fadbad::asin(x); }
  static TU acos(const TU& x) { return fadbad::acos(x); }
  static TU atan(const TU& x) { return fadbad::atan(x); }
  static TU erf (const TU& x) { throw std::runtime_error("mc::Op<fadbad::T<U>>::erf -- operation not permitted"); }
  static TU erfc(const TU& x) { throw std::runtime_error("mc::Op<fadbad::T<U>>::erfc -- operation not permitted"); }
  static TU fstep(const TU& x) { throw std::runtime_error("mc::Op<fadbad::T<U>>::fstep -- operation not permitted"); }
  static TU bstep(const TU& x) { throw std::runtime_error("mc::Op<fadbad::T<U>>::bstep -- operation not permitted"); }
  static TU hull(const TU& x, const TU& y) { throw std::runtime_error("mc::Op<fadbad::T<U>>::hull -- operation not permitted"); }
  static TU min (const TU& x, const TU& y) { throw std::runtime_error("mc::Op<fadbad::T<U>>::min -- operation not permitted"); }
  static TU max (const TU& x, const TU& y) { throw std::runtime_error("mc::Op<fadbad::T<U>>::max -- operation not permitted"); }
  static TU arh (const TU& x, const double k) { return fadbad::exp(-k/x); }
  static TU cheb(const TU& x, const unsigned n) { throw std::runtime_error("mc::Op<fadbad::B<U>>::max -- operation not permitted"); }
  template <typename X, typename Y> static TU pow(const X& x, const Y& y) { return fadbad::pow(x,y); }
  static TU monomial (const unsigned int n, const U* x, const int* k) { throw std::runtime_error("mc::Op<fadbad::T<U>>::monomial -- operation not permitted"); }
  static bool inter(TU& xIy, const TU& x, const TU& y) { throw std::runtime_error("mc::Op<fadbad::T<U>>::inter -- operation not permitted"); }
  static bool eq(const TU& x, const TU& y) { throw std::runtime_error("mc::Op<fadbad::T<U>>::eq -- operation not permitted"); }
  static bool ne(const TU& x, const TU& y) { throw std::runtime_error("mc::Op<fadbad::T<U>>::ne -- operation not permitted"); }
  static bool lt(const TU& x, const TU& y) { throw std::runtime_error("mc::Op<fadbad::T<U>>::lt -- operation not permitted"); }
  static bool le(const TU& x, const TU& y) { throw std::runtime_error("mc::Op<fadbad::T<U>>::le -- operation not permitted"); }
  static bool gt(const TU& x, const TU& y) { throw std::runtime_error("mc::Op<fadbad::T<U>>::gt -- operation not permitted"); }
  static bool ge(const TU& x, const TU& y) { throw std::runtime_error("mc::Op<fadbad::T<U>>::ge -- operation not permitted"); }
};

} // namespace mc

#endif
