// Copyright (C) Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the Eclipse Public License.

#ifndef MC__MCFADBAD_HPP
#define MC__MCFADBAD_HPP

#include "fadiff.h"
#include "mcop.hpp"

namespace fadbad
{

template <typename T, unsigned int N>
INLINE2 FTypeName<T, N>
pow2(const FTypeName<T, N>& a, const int b)
{
  FTypeName<T, N> c(Op<T>::myPow(a.val(), b));
  if (!a.depend()) return c;
  T tmp(Op<T>::myPow(a.val(), b - 1) * Op<T>::myInteger(b));
  c.setDepend(a);
  for (unsigned int i = 0; i < N; ++i) c[i] = tmp * a[i];
  return c;
}
template <typename T>
INLINE2 FTypeName<T, 0>
pow2(const FTypeName<T, 0>& a, const int b)
{
  FTypeName<T, 0> c(Op<T>::myPow(a.val(), b));
  if (!a.depend()) return c;
  T tmp(Op<T>::myPow(a.val(), b - 1) * Op<T>::myInteger(b));
  c.setDepend(a);
  for (unsigned int i = 0; i < c.size(); ++i) c[i] = tmp * a[i];
  return c;
}

template <typename T, unsigned int N>
INLINE2 FTypeName<T, N>
cheb(const FTypeName<T, N>& a, const unsigned b)
{
  FTypeName<T, N> c(mc::Op<T>::cheb(a.val(), b));
  if (!a.depend()) return c;

  T tmp(b % 2 ? 0.5 : 0.);
  for (int j = b - 1; j > 0; j -= 2) tmp += mc::Op<T>::cheb(a.val(), j);
  tmp *= 2 * (double)b;

  c.setDepend(a);
  for (unsigned int i = 0; i < N; ++i) c[i] = tmp * a[i];
  return c;
}
template <typename T>
INLINE2 FTypeName<T, 0>
cheb(const FTypeName<T, 0>& a, const unsigned b)
{
  FTypeName<T, 0> c(mc::Op<T>::cheb(a.val(), b));
  if (!a.depend()) return c;

  T tmp(b % 2 ? 0.5 : 0.);
  for (int j = b - 1; j > 0; j -= 2) tmp += mc::Op<T>::cheb(a.val(), j);
  tmp *= 2 * (double)b;

  c.setDepend(a);
  for (unsigned int i = 0; i < c.size(); ++i) c[i] = tmp * a[i];
  return c;
}

//@AVT.SVT: 01.06.2017
//@ICL: 26.04.2026
template <typename T, unsigned int N>
INLINE2 FTypeName<T, N>
xlog(const FTypeName<T, N>& a)
{
  // FTypeName<T,N> c(a*log(a));
  FTypeName<T, N> c(mc::Op<T>::xlog(a.val()));
  if (!a.depend()) return c;
  T tmp(Op<T>::myLog(a.val()) + Op<T>::myOne());
  c.setDepend(a);
  for (unsigned int i = 0; i < N; ++i) c[i] = tmp * a[i];
  return c;
}
template <typename T>
INLINE2 FTypeName<T, 0>
xlog(const FTypeName<T, 0>& a)
{
  // FTypeName<T,0> c(a*mc::Op<T>::log(a));
  FTypeName<T, 0> c(mc::Op<T>::xlog(a.val()));
  if (!a.depend()) return c;
  T tmp(Op<T>::myLog(a.val()) + Op<T>::myOne());
  c.setDepend(a);
  for (unsigned int i = 0; i < c.size(); ++i) c[i] = tmp * a[i];
  return c;
}

//@AVT.SVT: 07.06.2017
template <typename T, unsigned int N>
INLINE2 FTypeName<T, N>
lmtd(const FTypeName<T, N>& a, const FTypeName<T, N>& b)
{
  if (Op<T>::myEq(a.val(), b.val()))
  {
    FTypeName<T, N> c(a.val());
    c.setDepend(a, b);
    for (unsigned int i = 0; i < N; ++i) c[i] = 0.5 * a[i] + 0.5 * b[i];
    return c;
  }
  FTypeName<T, N> c((a - b) / (mc::Op<T>::log(a) - mc::Op<T>::log(b)));
  return c;
}
template <typename T>
INLINE2 FTypeName<T, 0>
lmtd(const FTypeName<T, 0>& a, const FTypeName<T, 0>& b)
{
  if (Op<T>::myEq(a.val(), b.val()))
  {
    FTypeName<T, 0> c(a.val());
    c.setDepend(a, b);
    for (unsigned int i = 0; i < c.size(); ++i) c[i] = 0.5 * a[i] + 0.5 * b[i];
    return c;
  }
  FTypeName<T, 0> c((a - b) / (mc::Op<T>::log(a) - mc::Op<T>::log(b)));
  return c;
}

//@AVT.SVT: 08.06.2017
template <typename T, unsigned int N>
INLINE2 FTypeName<T, N>
rlmtd(const FTypeName<T, N>& a, const FTypeName<T, N>& b)
{
  if (Op<T>::myEq(a.val(), b.val()))
  {
    FTypeName<T, N> c(1. / a.val());
    if (!a.depend()) return c;
    c.setDepend(a, b);
    for (unsigned int i = 0; i < N; ++i)
      c[i] = -a[i] / (2. * mc::Op<T>::sqr(a.val())) -
             b[i] / (2. * mc::Op<T>::sqr(b.val()));
    return c;
  }
  FTypeName<T, N> c((mc::Op<T>::log(a) - mc::Op<T>::log(b)) / (a - b));
  return c;
}
template <typename T>
INLINE2 FTypeName<T, 0>
rlmtd(const FTypeName<T, 0>& a, const FTypeName<T, 0>& b)
{
  if (Op<T>::myEq(a.val(), b.val()))
  {
    FTypeName<T, 0> c(1. / a.val());
    if (!a.depend()) return c;
    c.setDepend(a, b);
    for (unsigned int i = 0; i < c.size(); ++i)
      c[i] = -a[i] / (2. * mc::Op<T>::sqr(a.val())) -
             b[i] / (2. * mc::Op<T>::sqr(b.val()));
    return c;
  }
  FTypeName<T, 0> c((mc::Op<T>::log(a) - mc::Op<T>::log(b)) / (a - b));
  return c;
}

//@ICL: 04.01.2024
template <typename T, unsigned int N>
INLINE2 FTypeName<T, N>
fabs(const FTypeName<T, N>& a)
{
  FTypeName<T, N> c(mc::Op<T>::fabs(a.val()));
  if (!a.depend()) return c;
  T tmp(mc::Op<T>::fstep(a.val()) * Op<T>::myInteger(2) - Op<T>::myInteger(1));
  c.setDepend(a);
  for (unsigned int i = 0; i < N; ++i) c[i] = tmp * a[i];
  return c;
}
template <typename T>
INLINE2 FTypeName<T, 0>
fabs(const FTypeName<T, 0>& a)
{
  FTypeName<T, 0> c(mc::Op<T>::fabs(a.val()));
  if (!a.depend()) return c;
  T tmp(mc::Op<T>::fstep(a.val()) * Op<T>::myInteger(2) - Op<T>::myInteger(1));
  c.setDepend(a);
  for (unsigned int i = 0; i < c.size(); ++i) c[i] = tmp * a[i];
  return c;
}

//@ICL: 26.03.2026
template <typename T, typename U, unsigned int N>
INLINE2 FTypeName<T, N>
max(const FTypeName<T, N>& a, const U& b)
{
  FTypeName<T, N> c(mc::Op<T>::max(a.val(), b));
  if (!a.depend()) return c;
  T tmp(mc::Op<T>::fstep(a.val() - b));
  c.setDepend(a);
  for (unsigned int i = 0; i < N; ++i) c[i] = tmp * a[i];
  return c;
}
template <typename T, typename U>
INLINE2 FTypeName<T, 0>
max(const FTypeName<T, 0>& a, const U& b)
{
  FTypeName<T, 0> c(mc::Op<T>::max(a.val(), b));
  if (!a.depend()) return c;
  T tmp(mc::Op<T>::fstep(a.val() - b));
  c.setDepend(a);
  for (unsigned int i = 0; i < c.size(); ++i) c[i] = tmp * a[i];
  return c;
}

template <typename T, unsigned int N>
INLINE2 FTypeName<T, N>
max(const FTypeName<T, N>& a, const FTypeName<T, N>& b)
{
  FTypeName<T, N> c(mc::Op<T>::max(a.val(), b.val()));
  if (!a.depend() && !b.depend()) return c;
  c.setDepend(a, b);
  T tmp(mc::Op<T>::fstep(b.val() - a.val()));
  for (unsigned int i = 0; i < N; ++i)
    c[i] = (Op<T>::myInteger(1) - tmp) * a[i] + tmp * b[i];
  return c;
}
template <typename T>
INLINE2 FTypeName<T, 0>
max(const FTypeName<T, 0>& a, const FTypeName<T, 0>& b)
{
  FTypeName<T, 0> c(mc::Op<T>::max(a.val(), b.val()));
  if (!a.depend() && !b.depend()) return c;
  c.setDepend(a, b);
  T tmp(mc::Op<T>::fstep(b.val() - a.val()));
  for (unsigned int i = 0; i < c.size(); ++i)
    c[i] = (Op<T>::myInteger(1) - tmp) * a[i] + tmp * b[i];
  return c;
}

// template <typename T, unsigned int N>
// INLINE2 FTypeName<T,N> max(const FTypeName<T,N>& a, const FTypeName<T,N>& b)
//{
//         return 0.5*(a+b+mc::Op<FTypeName<T,N>>::fabs(a-b));
// }
// template <typename T>
// INLINE2 FTypeName<T,0> max(const FTypeName<T,0>& a, const FTypeName<T,0>& b)
//{
//         return 0.5*(a+b+mc::Op<FTypeName<T,0>>::fabs(a-b));
// }

template <typename T, typename U, unsigned int N>
INLINE2 FTypeName<T, N>
min(const FTypeName<T, N>& a, const U& b)
{
  FTypeName<T, N> c(mc::Op<T>::min(a.val(), b));
  if (!a.depend()) return c;
  T tmp(mc::Op<T>::fstep(b - a.val()));
  c.setDepend(a);
  for (unsigned int i = 0; i < N; ++i) c[i] = tmp * a[i];
  return c;
}
template <typename T, typename U>
INLINE2 FTypeName<T, 0>
min(const FTypeName<T, 0>& a, const U& b)
{
  FTypeName<T, 0> c(mc::Op<T>::min(a.val(), b));
  if (!a.depend()) return c;
  T tmp(mc::Op<T>::fstep(b - a.val()));
  c.setDepend(a);
  for (unsigned int i = 0; i < c.size(); ++i) c[i] = tmp * a[i];
  return c;
}

template <typename T, unsigned int N>
INLINE2 FTypeName<T, N>
min(const FTypeName<T, N>& a, const FTypeName<T, N>& b)
{
  FTypeName<T, N> c(mc::Op<T>::min(a.val(), b.val()));
  if (!a.depend() && !b.depend()) return c;
  c.setDepend(a, b);
  T tmp(mc::Op<T>::fstep(a.val() - b.val()));
  for (unsigned int i = 0; i < N; ++i)
    c[i] = (Op<T>::myInteger(1) - tmp) * a[i] + tmp * b[i];
  return c;
}
template <typename T>
INLINE2 FTypeName<T, 0>
min(const FTypeName<T, 0>& a, const FTypeName<T, 0>& b)
{
  FTypeName<T, 0> c(mc::Op<T>::min(a.val(), b.val()));
  if (!a.depend() && !b.depend()) return c;
  c.setDepend(a, b);
  T tmp(mc::Op<T>::fstep(a.val() - b.val()));
  for (unsigned int i = 0; i < c.size(); ++i)
    c[i] = (Op<T>::myInteger(1) - tmp) * a[i] + tmp * b[i];
  return c;
}

// template <typename T, unsigned int N>
// INLINE2 FTypeName<T,N> min (const FTypeName<T,N>& a, const FTypeName<T,N>& b)
//{
//         return 0.5*(a+b-mc::Op<FTypeName<T,N>>::fabs(a-b));
// }
// template <typename T>
// INLINE2 FTypeName<T,0> min (const FTypeName<T,0>& a, const FTypeName<T,0>& b)
//{
//         return 0.5*(a+b-mc::Op<FTypeName<T,0>>::fabs(a-b));
// }

//@ICL: 08.12.2025
template <typename T, unsigned int N>
INLINE2 FTypeName<T, N>
erf(const FTypeName<T, N>& a)
{
  FTypeName<T, N> c(mc::Op<T>::erf(a.val()));
  if (!a.depend()) return c;
  T tmp(2. / std::sqrt(mc::PI) * mc::Op<T>::exp(-mc::Op<T>::sqr(a.val())));
  c.setDepend(a);
  for (unsigned int i = 0; i < N; ++i) c[i] = tmp * a[i];
  return c;
}
template <typename T>
INLINE2 FTypeName<T, 0>
erf(const FTypeName<T, 0>& a)
{
  FTypeName<T, 0> c(mc::Op<T>::erf(a.val()));
  if (!a.depend()) return c;
  T tmp(2. / std::sqrt(mc::PI) * mc::Op<T>::exp(-mc::Op<T>::sqr(a.val())));
  c.setDepend(a);
  for (unsigned int i = 0; i < c.size(); ++i) c[i] = tmp * a[i];
  return c;
}

}  // end namespace fadbad

#include "badiff.h"

namespace fadbad
{
//@ICL: 26.04.2026
template <typename U>
struct BTypeNameABS : public UnBTypeNameHV<U>
{
  BTypeNameABS(const U& val, BTypeNameHV<U>* pOp) : UnBTypeNameHV<U>(val, pOp)
  {
  }
  virtual void
  propagate(typename Derivatives<U>::RecycleBin& bin)
  {
    U tmp(mc::Op<U>::fstep(this->op()->val()) * Op<U>::myInteger(2) -
          Op<U>::myInteger(1));
    this->op()->add(bin, tmp, this->m_derivatives);
  }

 private:
  void
  operator=(const BTypeNameABS<U>&)
  {
  }  // not allowed
};
template <typename U>
BTypeName<U>
fabs(const BTypeName<U>& x)
{
  return BTypeName<U>(static_cast<BTypeNameHV<U>*>(
      new BTypeNameABS<U>(mc::Op<U>::fabs(x.val()), x.getBTypeNameHV())));
}

template <typename U, typename V>
struct BTypeNameMAX2 : public UnBTypeNameHV<U>
{
  const V m_c;
  BTypeNameMAX2(const U& val, BTypeNameHV<U>* pOp, const V& c)
      : UnBTypeNameHV<U>(val, pOp), m_c(c)
  {
  }
  virtual void
  propagate(typename Derivatives<U>::RecycleBin& bin)
  {
    U tmp(mc::Op<U>::fstep(this->op()->val() - m_c));
    this->op()->add(bin, tmp, this->m_derivatives);
  }

 private:
  void
  operator=(const BTypeNameMAX2<U, V>&)
  {
  }  // not allowed
};
template <typename U, typename V>
BTypeName<U>
max(const BTypeName<U>& x, const V& c)
{
  return BTypeName<U>(static_cast<BTypeNameHV<U>*>(new BTypeNameMAX2<U, V>(
      mc::Op<U>::max(x.val(), c), x.getBTypeNameHV(), c)));
}
template <typename U, typename V>
BTypeName<U>
max(const V& c, const BTypeName<U>& x)
{
  return max(x, c);
}

// template <typename U>
// BTypeName<U> max(const BTypeName<U>& a, const BTypeName<U>& b)
//{
//         return 0.5*(mc::Op<BTypeName<U>>::fabs(a-b)+a+b);
// }
template <typename U>
struct BTypeNameMAX : public BinBTypeNameHV<U>
{
  BTypeNameMAX(const U& val, BTypeNameHV<U>* pOp1, BTypeNameHV<U>* pOp2)
      : BinBTypeNameHV<U>(val, pOp1, pOp2)
  {
  }
  virtual void
  propagate(typename Derivatives<U>::RecycleBin& bin)
  {
    U tmp(mc::Op<U>::fstep(this->op2()->val() - this->op1()->val()));
    this->op1()->add(bin, Op<U>::myInteger(1) - tmp, this->m_derivatives);
    this->op2()->add(bin, tmp, this->m_derivatives);
  }

 private:
  void
  operator=(const BTypeNameMAX<U>&)
  {
  }  // not allowed
};
template <typename U>
BTypeName<U>
max(const BTypeName<U>& x, const BTypeName<U>& y)
{
  return BTypeName<U>(static_cast<BTypeNameHV<U>*>(
      new BTypeNameMAX<U>(mc::Op<U>::max(x.val(), y.val()), x.getBTypeNameHV(),
                          y.getBTypeNameHV())));
}

template <typename U, typename V>
struct BTypeNameMIN2 : public UnBTypeNameHV<U>
{
  const V m_c;
  BTypeNameMIN2(const U& val, BTypeNameHV<U>* pOp, const V& c)
      : UnBTypeNameHV<U>(val, pOp), m_c(c)
  {
  }
  virtual void
  propagate(typename Derivatives<U>::RecycleBin& bin)
  {
    U tmp(mc::Op<U>::fstep(m_c - this->op()->val()));
    this->op()->add(bin, tmp, this->m_derivatives);
  }

 private:
  void
  operator=(const BTypeNameMIN2<U, V>&)
  {
  }  // not allowed
};
template <typename U, typename V>
BTypeName<U>
min(const BTypeName<U>& x, const V& c)
{
  return BTypeName<U>(static_cast<BTypeNameHV<U>*>(new BTypeNameMIN2<U, V>(
      mc::Op<U>::min(x.val(), c), x.getBTypeNameHV(), c)));
}
template <typename U, typename V>
BTypeName<U>
min(const V& c, const BTypeName<U>& x)
{
  return min(x, c);
}

// template <typename U>
// BTypeName<U> min(const BTypeName<U>& a, const BTypeName<U>& b)
//{
//         return 0.5*(a+b-mc::Op<BTypeName<U>>::fabs(a-b));
// }
template <typename U>
struct BTypeNameMIN : public BinBTypeNameHV<U>
{
  BTypeNameMIN(const U& val, BTypeNameHV<U>* pOp1, BTypeNameHV<U>* pOp2)
      : BinBTypeNameHV<U>(val, pOp1, pOp2)
  {
  }
  virtual void
  propagate(typename Derivatives<U>::RecycleBin& bin)
  {
    U tmp(mc::Op<U>::fstep(this->op1()->val() - this->op2()->val()));
    this->op1()->add(bin, Op<U>::myInteger(1) - tmp, this->m_derivatives);
    this->op2()->add(bin, tmp, this->m_derivatives);
  }

 private:
  void
  operator=(const BTypeNameMIN<U>&)
  {
  }  // not allowed
};
template <typename U>
BTypeName<U>
min(const BTypeName<U>& x, const BTypeName<U>& y)
{
  return BTypeName<U>(static_cast<BTypeNameHV<U>*>(
      new BTypeNameMIN<U>(mc::Op<U>::min(x.val(), y.val()), x.getBTypeNameHV(),
                          y.getBTypeNameHV())));
}

template <typename U>
struct BTypeNameXLOG : public UnBTypeNameHV<U>
{
  BTypeNameXLOG(const U& val, BTypeNameHV<U>* pOp) : UnBTypeNameHV<U>(val, pOp)
  {
  }
  virtual void
  propagate(typename Derivatives<U>::RecycleBin& bin)
  {
    U tmp(mc::Op<U>::log(this->op()->val()) + Op<U>::myInteger(1));
    this->op()->add(bin, tmp, this->m_derivatives);
  }

 private:
  void
  operator=(const BTypeNameXLOG<U>&)
  {
  }  // not allowed
};
template <typename U>
BTypeName<U>
xlog(const BTypeName<U>& x)
{
  return BTypeName<U>(static_cast<BTypeNameHV<U>*>(
      new BTypeNameXLOG<U>(mc::Op<U>::xlog(x.val()), x.getBTypeNameHV())));
}

template <typename U>
struct BTypeNameERF : public UnBTypeNameHV<U>
{
  BTypeNameERF(const U& val, BTypeNameHV<U>* pOp) : UnBTypeNameHV<U>(val, pOp)
  {
  }
  virtual void
  propagate(typename Derivatives<U>::RecycleBin& bin)
  {
    U tmp(2. / std::sqrt(mc::PI) *
          mc::Op<U>::exp(-mc::Op<U>::sqr(this->op()->val())));
    this->op()->add(bin, tmp, this->m_derivatives);
  }

 private:
  void
  operator=(const BTypeNameERF<U>&)
  {
  }  // not allowed
};
template <typename U>
BTypeName<U>
erf(const BTypeName<U>& x)
{
  return BTypeName<U>(static_cast<BTypeNameHV<U>*>(
      new BTypeNameERF<U>(mc::Op<U>::erf(x.val()), x.getBTypeNameHV())));
}

template <typename U>
struct BTypeNameCHEB : public UnBTypeNameHV<U>
{
  const unsigned _n;
  BTypeNameCHEB(const U& val, BTypeNameHV<U>* pOp, const unsigned n)
      : UnBTypeNameHV<U>(val, pOp), _n(n)
  {
  }
  virtual void
  propagate(typename Derivatives<U>::RecycleBin& bin)
  {
    U tmp(_n % 2 ? 0.5 : 0.);
    for (int j = _n - 1; j > 0; j -= 2)
      tmp += mc::Op<U>::cheb(this->op()->val(), j);
    tmp *= 2 * (double)_n;
    this->op()->add(bin, tmp, this->m_derivatives);
  }

 private:
  void
  operator=(const BTypeNameCHEB<U>&)
  {
  }  // not allowed
};
template <typename U>
BTypeName<U>
cheb(const BTypeName<U>& x, const unsigned n)
{
  return BTypeName<U>(static_cast<BTypeNameHV<U>*>(new BTypeNameCHEB<U>(
      mc::Op<U>::cheb(x.val(), n), x.getBTypeNameHV(), n)));
}

}  // namespace fadbad

#include "tadiff.h"

namespace fadbad
{

template <typename U, int N>
struct TTypeNamePOW3 : public UnTTypeNameHV<U, N>
{
  int m_b;
  TTypeNamePOW3(const U& val, TTypeNameHV<U, N>* pOp, const int b)
      : UnTTypeNameHV<U, N>(val, pOp), m_b(b)
  {
  }
  TTypeNamePOW3(TTypeNameHV<U, N>* pOp, const int b)
      : UnTTypeNameHV<U, N>(pOp), m_b(b)
  {
  }
  unsigned int
  eval(const unsigned int k)
  {
    unsigned int l = this->opEval(k);
    if (m_b == 0)
    {
      if (0 == this->length())
      {
        this->val(0)   = Op<U>::myOne();
        this->length() = 1;
      }
      for (unsigned int i = this->length(); i < l; ++i)
      {
        this->val(i) = Op<U>::myZero();
      }
    }
    else if (m_b == 1)
    {
      for (unsigned int i = this->length(); i < l; ++i)
      {
        this->val(i) = this->opVal(i);
      }
    }
    else if (m_b == 2)
    {
      if (0 == this->length())
      {
        this->val(0)   = Op<U>::mySqr(this->opVal(0));
        this->length() = 1;
      }
      for (unsigned int i = this->length(); i < l; ++i)
      {
        this->val(i)   = Op<U>::myZero();
        unsigned int m = (i + 1) / 2;
        for (unsigned int j = 0; j < m; ++j)
          Op<U>::myCadd(this->val(i), this->opVal(i - j) * this->opVal(j));
        Op<U>::myCmul(this->val(i), Op<U>::myTwo());
        if (0 == i % 2)
          Op<U>::myCadd(this->val(i), Op<U>::mySqr(this->opVal(m)));
      }
    }
    else if (m_b == 3)
    {
      if (0 == this->length())
      {
        this->val(0)   = Op<U>::myPow(this->opVal(0), m_b);
        this->length() = 1;
      }
      if (1 < l && 1 == this->length())
      {
        this->val(1) = Op<U>::myPow(this->opVal(0), m_b - 1) * this->opVal(1) *
                       Op<U>::myInteger(m_b);
        this->length() = 2;
      }
      if (2 < l && 2 == this->length())
      {
        this->val(2) =
            Op<U>::myPow(this->opVal(0), m_b - 2) *
            (this->opVal(0) * this->opVal(2) +
             Op<U>::myInteger(m_b - 1) * Op<U>::mySqr(this->opVal(1))) *
            Op<U>::myInteger(m_b);
        this->length() = 3;
      }
      for (unsigned int i = this->length(); i < l; ++i)
      {
        this->val(i)   = Op<U>::myZero();
        unsigned int m = (i + 1) / 2;
        for (unsigned int j = 0; j < m; ++j)
          Op<U>::myCadd(this->val(i), this->opVal(i - j) * this->opVal(j));
        Op<U>::myCmul(this->val(i), Op<U>::myTwo());
        if (0 == i % 2)
          Op<U>::myCadd(this->val(i), Op<U>::mySqr(this->opVal(m)));
      }
      for (unsigned int i = l - 1; i >= this->length(); --i)
      {
        Op<U>::myCmul(this->val(i), this->opVal(0));
        for (unsigned int j = 1; j <= i; ++j)
          Op<U>::myCadd(this->val(i), this->val(i - j) * this->opVal(j));
      }
    }
    else
    {
      if (0 == this->length())
      {
        this->val(0)   = Op<U>::myPow(this->opVal(0), m_b);
        this->length() = 1;
      }
      for (unsigned int i = this->length(); i < l; ++i)
      {
        this->val(i) = Op<U>::myZero();
        for (unsigned int j = 0; j < i; ++j)
          Op<U>::myCadd(this->val(i),
                        (m_b - (m_b + Op<U>::myOne()) * Op<U>::myInteger(j) /
                                   Op<U>::myInteger(i)) *
                            this->opVal(i - j) * this->val(j));
      }
    }
    return this->length() = l;
  }

 private:
  void
  operator=(const TTypeNamePOW3<U, N>&)
  {
  }  // not allowed
};

template <typename U, int N>
TTypeName<U, N>
pow(const TTypeName<U, N>& val, const int b)
{
  TTypeNameHV<U, N>* pHV =
      val.length() > 0 ? new TTypeNamePOW3<U, N>(Op<U>::myPow(val.val(), b),
                                                 val.getTTypeNameHV(), b)
                       : new TTypeNamePOW3<U, N>(val.getTTypeNameHV(), b);
  return TTypeName<U, N>(pHV);
}

}  // end namespace fadbad

namespace mc
{

//! @brief C++ structure for specialization of the mc::Op templated structure
//! for use of the FADBAD type fadbad::F inside other MC++ types
template <typename U>
struct Op<fadbad::F<U> >
{
  typedef fadbad::F<U> TU;
  static TU
  point(const double c)
  {
    throw std::runtime_error(
        "mc::Op<fadbad::F<U>>::point -- operation not permitted");
  }
  static TU
  zeroone()
  {
    throw std::runtime_error(
        "mc::Op<fadbad::F<U>>::zeroone -- operation not permitted");
  }
  static void
  I(TU& x, const TU& y)
  {
    x = y;
  }
  static double
  l(const TU& x)
  {
    throw std::runtime_error(
        "mc::Op<fadbad::F<U>>::l -- operation not permitted");
  }
  static double
  u(const TU& x)
  {
    throw std::runtime_error(
        "mc::Op<fadbad::F<U>>::u -- operation not permitted");
  }
  static double
  abs(const TU& x)
  {
    throw std::runtime_error(
        "mc::Op<fadbad::F<U>>::abs -- operation not permitted");
  }
  static double
  mid(const TU& x)
  {
    throw std::runtime_error(
        "mc::Op<fadbad::F<U>>::mid -- operation not permitted");
  }
  static double
  diam(const TU& x)
  {
    throw std::runtime_error(
        "mc::Op<fadbad::F<U>>::diam -- operation not permitted");
  }
  static TU
  inv(const TU& x)
  {
    return 1. / x;
  }
  static TU
  sqr(const TU& x)
  {
    return fadbad::sqr(x);
  }
  static TU
  sqrt(const TU& x)
  {
    return fadbad::sqrt(x);
  }
  static TU
  exp(const TU& x)
  {
    return fadbad::exp(x);
  }
  static TU
  log(const TU& x)
  {
    return fadbad::log(x);
  }
  static TU
  xlog(const TU& x)
  {
    return fadbad::xlog(x);
  }
  static TU
  lmtd(const TU& x, const TU& y)
  {
    return fadbad::lmtd(x, y);
  }
  static TU
  rlmtd(const TU& x, const TU& y)
  {
    return fadbad::rlmtd(x, y);
  }
  static TU
  fabs(const TU& x)
  {
    return fadbad::fabs(x);
  }
  static TU
  sin(const TU& x)
  {
    return fadbad::sin(x);
  }
  static TU
  cos(const TU& x)
  {
    return fadbad::cos(x);
  }
  static TU
  tan(const TU& x)
  {
    return fadbad::tan(x);
  }
  static TU
  asin(const TU& x)
  {
    return fadbad::asin(x);
  }
  static TU
  acos(const TU& x)
  {
    return fadbad::acos(x);
  }
  static TU
  atan(const TU& x)
  {
    return fadbad::atan(x);
  }
  static TU
  sinh(const TU& x)
  {
    return fadbad::sinh(x);
  }
  static TU
  cosh(const TU& x)
  {
    return fadbad::cosh(x);
  }
  static TU
  tanh(const TU& x)
  {
    return fadbad::tanh(x);
  }
  static TU
  erf(const TU& x)
  {
    return fadbad::erf(x);
  }
  static TU
  erfc(const TU& x)
  {
    return 1. - fadbad::erf(x);
  }
  static TU
  fstep(const TU& x)
  {
    throw std::runtime_error(
        "mc::Op<fadbad::F<U>>::fstep -- operation not permitted");
  }
  static TU
  bstep(const TU& x)
  {
    throw std::runtime_error(
        "mc::Op<fadbad::F<U>>::bstep -- operation not permitted");
  }
  static TU
  hull(const TU& x, const TU& y)
  {
    throw std::runtime_error(
        "mc::Op<fadbad::F<U>>::hull -- operation not permitted");
  }
  template <typename Y>
  static TU
  min(const TU& x, const Y& y)
  {
    return fadbad::min(x, y);
  }
  template <typename Y>
  static TU
  max(const TU& x, const Y& y)
  {
    return fadbad::max(x, y);
  }
  //  static TU min (const TU& x, const TU& y) { return
  //  0.5*(x+y-fadbad::fabs(x-y)); } static TU max (const TU& x, const TU& y) {
  //  return 0.5*(x+y+fadbad::fabs(x-y)); }
  static TU
  arh(const TU& x, const double k)
  {
    return fadbad::exp(-k / x);
  }
  template <typename X, typename Y>
  static TU
  pow(const X& x, const Y& y)
  {
    return fadbad::pow(x, y);
  }
  static TU
  cheb(const TU& x, const unsigned n)
  {
    return fadbad::cheb(x, n);
  }
  static TU
  prod(const unsigned n, const TU* x)
  {
    switch (n)
    {
      case 0:
        return 1.;
      case 1:
        return x[0];
      default:
        return x[0] * prod(n - 1, x + 1);
    }
  }
  static TU
  monom(const unsigned n, const TU* x, const unsigned* k)
  {
    switch (n)
    {
      case 0:
        return 1.;
      case 1:
        return pow(x[0], (int)k[0]);
      default:
        return pow(x[0], (int)k[0]) * monom(n - 1, x + 1, k + 1);
    }
  }
  static bool
  inter(TU& xIy, const TU& x, const TU& y)
  {
    xIy = x;
    return true;
  }
  static bool
  eq(const TU& x, const TU& y)
  {
    throw std::runtime_error(
        "mc::Op<fadbad::F<U>>::eq -- operation not permitted");
  }
  static bool
  ne(const TU& x, const TU& y)
  {
    throw std::runtime_error(
        "mc::Op<fadbad::F<U>>::ne -- operation not permitted");
  }
  static bool
  lt(const TU& x, const TU& y)
  {
    throw std::runtime_error(
        "mc::Op<fadbad::F<U>>::lt -- operation not permitted");
  }
  static bool
  le(const TU& x, const TU& y)
  {
    throw std::runtime_error(
        "mc::Op<fadbad::F<U>>::le -- operation not permitted");
  }
  static bool
  gt(const TU& x, const TU& y)
  {
    throw std::runtime_error(
        "mc::Op<fadbad::F<U>>::gt -- operation not permitted");
  }
  static bool
  ge(const TU& x, const TU& y)
  {
    throw std::runtime_error(
        "mc::Op<fadbad::F<U>>::ge -- operation not permitted");
  }
};

//! @brief C++ structure for specialization of the mc::Op templated structure
//! for use of the FADBAD type fadbad::B inside other MC++ types
template <typename U>
struct Op<fadbad::B<U> >
{
  typedef fadbad::B<U> TU;
  static TU
  point(const double c)
  {
    throw std::runtime_error(
        "mc::Op<fadbad::B<U>>::point -- operation not permitted");
  }
  static TU
  zeroone()
  {
    throw std::runtime_error(
        "mc::Op<fadbad::B<U>>::zeroone -- operation not permitted");
  }
  static void
  I(TU& x, const TU& y)
  {
    x = y;
  }
  static double
  l(const TU& x)
  {
    throw std::runtime_error(
        "mc::Op<fadbad::B<U>>::l -- operation not permitted");
  }
  static double
  u(const TU& x)
  {
    throw std::runtime_error(
        "mc::Op<fadbad::B<U>>::u -- operation not permitted");
  }
  static double
  abs(const TU& x)
  {
    throw std::runtime_error(
        "mc::Op<fadbad::B<U>>::abs -- operation not permitted");
  }
  static double
  mid(const TU& x)
  {
    throw std::runtime_error(
        "mc::Op<fadbad::B<U>>::mid -- operation not permitted");
  }
  static double
  diam(const TU& x)
  {
    throw std::runtime_error(
        "mc::Op<fadbad::B<U>>::diam -- operation not permitted");
  }
  static TU
  inv(const TU& x)
  {
    return 1. / x;
  }
  static TU
  sqr(const TU& x)
  {
    return fadbad::sqr(x);
  }
  static TU
  sqrt(const TU& x)
  {
    return fadbad::sqrt(x);
  }
  static TU
  exp(const TU& x)
  {
    return fadbad::exp(x);
  }
  static TU
  log(const TU& x)
  {
    return fadbad::log(x);
  }
  static TU
  xlog(const TU& x)
  {
    return fadbad::xlog(x);
  }
  static TU
  lmtd(const TU& x, const TU& y)
  {
    throw std::runtime_error(
        "mc::Op<fadbad::B<U>>::lmtd -- operation not permitted");
  }
  static TU
  rlmtd(const TU& x, const TU& y)
  {
    throw std::runtime_error(
        "mc::Op<fadbad::B<U>>::rlmtd -- operation not permitted");
  }
  static TU
  fabs(const TU& x)
  {
    return fadbad::fabs(x);
  }
  static TU
  sin(const TU& x)
  {
    return fadbad::sin(x);
  }
  static TU
  cos(const TU& x)
  {
    return fadbad::cos(x);
  }
  static TU
  tan(const TU& x)
  {
    return fadbad::tan(x);
  }
  static TU
  asin(const TU& x)
  {
    return fadbad::asin(x);
  }
  static TU
  acos(const TU& x)
  {
    return fadbad::acos(x);
  }
  static TU
  atan(const TU& x)
  {
    return fadbad::atan(x);
  }
  static TU
  sinh(const TU& x)
  {
    return fadbad::sinh(x);
  }
  static TU
  cosh(const TU& x)
  {
    return fadbad::cosh(x);
  }
  static TU
  tanh(const TU& x)
  {
    return fadbad::tanh(x);
  }
  static TU
  erf(const TU& x)
  {
    return fadbad::erf(x);
  }
  static TU
  erfc(const TU& x)
  {
    return 1. - fadbad::erf(x);
  }
  static TU
  fstep(const TU& x)
  {
    throw std::runtime_error(
        "mc::Op<fadbad::B<U>>::fstep -- operation not permitted");
  }
  static TU
  bstep(const TU& x)
  {
    throw std::runtime_error(
        "mc::Op<fadbad::B<U>>::bstep -- operation not permitted");
  }
  static TU
  hull(const TU& x, const TU& y)
  {
    throw std::runtime_error(
        "mc::Op<fadbad::B<U>>::hull -- operation not permitted");
  }
  template <typename Y>
  static TU
  min(const TU& x, const Y& y)
  {
    return fadbad::min(x, y);
  }
  template <typename Y>
  static TU
  max(const TU& x, const Y& y)
  {
    return fadbad::max(x, y);
  }
  static TU
  arh(const TU& x, const double k)
  {
    return fadbad::exp(-k / x);
  }
  template <typename X, typename Y>
  static TU
  pow(const X& x, const Y& y)
  {
    return fadbad::pow(x, y);
  }
  static TU
  cheb(const TU& x, const unsigned n)
  {
    return fadbad::cheb(x, n);
  }
  static TU
  prod(const unsigned n, const TU* x)
  {
    switch (n)
    {
      case 0:
        return 1.;
      case 1:
        return x[0];
      default:
        return x[0] * prod(n - 1, x + 1);
    }
  }
  static TU
  monom(const unsigned n, const TU* x, const unsigned* k)
  {
    switch (n)
    {
      case 0:
        return 1.;
      case 1:
        return pow(x[0], (int)k[0]);
      default:
        return pow(x[0], (int)k[0]) * monom(n - 1, x + 1, k + 1);
    }
  }
  static bool
  inter(TU& xIy, const TU& x, const TU& y)
  {
    xIy = x;
    return true;
  }
  static bool
  eq(const TU& x, const TU& y)
  {
    throw std::runtime_error(
        "mc::Op<fadbad::B<U>>::eq -- operation not permitted");
  }
  static bool
  ne(const TU& x, const TU& y)
  {
    throw std::runtime_error(
        "mc::Op<fadbad::B<U>>::ne -- operation not permitted");
  }
  static bool
  lt(const TU& x, const TU& y)
  {
    throw std::runtime_error(
        "mc::Op<fadbad::B<U>>::lt -- operation not permitted");
  }
  static bool
  le(const TU& x, const TU& y)
  {
    throw std::runtime_error(
        "mc::Op<fadbad::B<U>>::le -- operation not permitted");
  }
  static bool
  gt(const TU& x, const TU& y)
  {
    throw std::runtime_error(
        "mc::Op<fadbad::B<U>>::gt -- operation not permitted");
  }
  static bool
  ge(const TU& x, const TU& y)
  {
    throw std::runtime_error(
        "mc::Op<fadbad::B<U>>::ge -- operation not permitted");
  }
};

//! @brief C++ structure for specialization of the mc::Op templated structure
//! for use of the FADBAD type fadbad::T inside other MC++ types
template <typename U>
struct Op<fadbad::T<U> >
{
  typedef fadbad::T<U> TU;
  static TU
  point(const double c)
  {
    throw std::runtime_error(
        "mc::Op<fadbad::T<U>>::point -- operation not permitted");
  }
  static TU
  zeroone()
  {
    throw std::runtime_error(
        "mc::Op<fadbad::T<U>>::zeroone -- operation not permitted");
  }
  static void
  I(TU& x, const TU& y)
  {
    x = y;
  }
  static double
  l(const TU& x)
  {
    throw std::runtime_error(
        "mc::Op<fadbad::T<U>>::l -- operation not permitted");
  }
  static double
  u(const TU& x)
  {
    throw std::runtime_error(
        "mc::Op<fadbad::T<U>>::u -- operation not permitted");
  }
  static double
  abs(const TU& x)
  {
    throw std::runtime_error(
        "mc::Op<fadbad::T<U>>::abs -- operation not permitted");
  }
  static double
  mid(const TU& x)
  {
    throw std::runtime_error(
        "mc::Op<fadbad::T<U>>::mid -- operation not permitted");
  }
  static double
  diam(const TU& x)
  {
    throw std::runtime_error(
        "mc::Op<fadbad::T<U>>::diam -- operation not permitted");
  }
  static TU
  inv(const TU& x)
  {
    return 1. / x;
  }
  static TU
  sqr(const TU& x)
  {
    return fadbad::sqr(x);
  }
  static TU
  sqrt(const TU& x)
  {
    return fadbad::sqrt(x);
  }
  static TU
  exp(const TU& x)
  {
    return fadbad::exp(x);
  }
  static TU
  log(const TU& x)
  {
    return fadbad::log(x);
  }
  static TU
  xlog(const TU& x)
  {
    return x * fadbad::log(x);
  }
  static TU
  lmtd(const TU& x, const TU& y)
  {
    throw std::runtime_error(
        "mc::Op<fadbad::T<U>>::lmtd -- operation not permitted");
  }
  static TU
  rlmtd(const TU& x, const TU& y)
  {
    throw std::runtime_error(
        "mc::Op<fadbad::T<U>>::rlmtd -- operation not permitted");
  }
  static TU
  fabs(const TU& x)
  {
    throw std::runtime_error(
        "mc::Op<fadbad::T<U>>::fabs -- operation not permitted");
  }
  static TU
  sin(const TU& x)
  {
    return fadbad::sin(x);
  }
  static TU
  cos(const TU& x)
  {
    return fadbad::cos(x);
  }
  static TU
  tan(const TU& x)
  {
    return fadbad::tan(x);
  }
  static TU
  asin(const TU& x)
  {
    return fadbad::asin(x);
  }
  static TU
  acos(const TU& x)
  {
    return fadbad::acos(x);
  }
  static TU
  atan(const TU& x)
  {
    return fadbad::atan(x);
  }
  static TU
  sinh(const TU& x)
  {
    return fadbad::sinh(x);
  }
  static TU
  cosh(const TU& x)
  {
    return fadbad::cosh(x);
  }
  static TU
  tanh(const TU& x)
  {
    return fadbad::tanh(x);
  }
  static TU
  erf(const TU& x)
  {
    throw std::runtime_error(
        "mc::Op<fadbad::T<U>>::erf -- operation not permitted");
  }
  static TU
  erfc(const TU& x)
  {
    throw std::runtime_error(
        "mc::Op<fadbad::T<U>>::erfc -- operation not permitted");
  }
  static TU
  fstep(const TU& x)
  {
    throw std::runtime_error(
        "mc::Op<fadbad::T<U>>::fstep -- operation not permitted");
  }
  static TU
  bstep(const TU& x)
  {
    throw std::runtime_error(
        "mc::Op<fadbad::T<U>>::bstep -- operation not permitted");
  }
  static TU
  hull(const TU& x, const TU& y)
  {
    throw std::runtime_error(
        "mc::Op<fadbad::T<U>>::hull -- operation not permitted");
  }
  static TU
  min(const TU& x, const TU& y)
  {
    throw std::runtime_error(
        "mc::Op<fadbad::T<U>>::min -- operation not permitted");
  }
  static TU
  max(const TU& x, const TU& y)
  {
    throw std::runtime_error(
        "mc::Op<fadbad::T<U>>::max -- operation not permitted");
  }
  static TU
  arh(const TU& x, const double k)
  {
    return fadbad::exp(-k / x);
  }
  template <typename X, typename Y>
  static TU
  pow(const X& x, const Y& y)
  {
    return fadbad::pow(x, y);
  }
  static TU
  cheb(const TU& x, const unsigned n)
  {
    throw std::runtime_error(
        "mc::Op<fadbad::T<U>>::cheb -- operation not permitted");
  }
  static TU
  prod(const unsigned n, const TU* x)
  {
    switch (n)
    {
      case 0:
        return 1.;
      case 1:
        return x[0];
      default:
        return x[0] * prod(n - 1, x + 1);
    }
  }
  static TU
  monom(const unsigned n, const TU* x, const unsigned* k)
  {
    switch (n)
    {
      case 0:
        return 1.;
      case 1:
        return pow(x[0], (int)k[0]);
      default:
        return pow(x[0], (int)k[0]) * monom(n - 1, x + 1, k + 1);
    }
  }
  static bool
  inter(TU& xIy, const TU& x, const TU& y)
  {
    xIy = x;
    return true;
  }
  static bool
  eq(const TU& x, const TU& y)
  {
    throw std::runtime_error(
        "mc::Op<fadbad::T<U>>::eq -- operation not permitted");
  }
  static bool
  ne(const TU& x, const TU& y)
  {
    throw std::runtime_error(
        "mc::Op<fadbad::T<U>>::ne -- operation not permitted");
  }
  static bool
  lt(const TU& x, const TU& y)
  {
    throw std::runtime_error(
        "mc::Op<fadbad::T<U>>::lt -- operation not permitted");
  }
  static bool
  le(const TU& x, const TU& y)
  {
    throw std::runtime_error(
        "mc::Op<fadbad::T<U>>::le -- operation not permitted");
  }
  static bool
  gt(const TU& x, const TU& y)
  {
    throw std::runtime_error(
        "mc::Op<fadbad::T<U>>::gt -- operation not permitted");
  }
  static bool
  ge(const TU& x, const TU& y)
  {
    throw std::runtime_error(
        "mc::Op<fadbad::T<U>>::ge -- operation not permitted");
  }
};

}  // namespace mc

#endif
