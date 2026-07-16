// Copyright (C) Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// fdiff.hpp -- Independent (clean-room) forward-mode automatic differentiation
// types, a license-clean replacement for fadbad::F of FADBAD++ (fadiff.h). The
// derivative rules are the standard chain rule expressed entirely through the
// mc::Op<T> abstraction, so any evaluation type T that specialises mc::Op
// (double, mc::Interval, mc::McCormick, and the forward types themselves for
// higher-order derivatives) can be used.
//
//   * mc::F<T>   -- dense gradient stored as a contiguous vector of length N
//   * mc::SF<T>  -- sparse gradient stored as (seed index, derivative) pairs
//
// Both are a single implementation, Fimpl<T,MODE>, parameterised by the
// gradient representation (MODE 0 = dense std::vector, MODE 1 = sparse
// std::map), exposed through the alias templates below, with one mc::Op
// specialisation and one set of elementary functions serving both -- the dense
// and sparse variants share every derivative rule and differ only in the two
// gradient-combination primitives (_chain / _chain2) and the storage
// accessors, selected with `if constexpr`.
//
// The public member interface (val, x, size, dimension, depend, operator[],
// deriv, d, diff, setDepend, compound operators, begin/end/gradient) matches
// fadbad::F / the previous mc::F & mc::SF so existing call sites keep working.

#ifndef MC__FDIFF_HPP
#define MC__FDIFF_HPP

#include <vector>
#include <map>
#include <cstddef>
#include <cmath>
#include <stdexcept>
#include <algorithm>

#include "mcop.hpp"

namespace mc
{

// MODE 0 : dense gradient (std::vector<T>), empty == independent constant
// MODE 1 : sparse gradient (std::map<unsigned,T>) keyed by seed index
template <typename T, int MODE>
class Fimpl
{
public:
  typedef T UnderlyingType;
  typedef T value_type;
  typedef typename std::conditional<MODE==0, std::vector<T>,
                                    std::map<unsigned,T> >::type container;
  typedef typename container::iterator       iterator;
  typedef typename container::const_iterator const_iterator;

private:
  T         _v;        //!< function value
  container _g;        //!< gradient (dense: length N, empty==const; sparse: idx->deriv)
  unsigned  _n = 0;    //!< sparse seeding dimension (dense ignores; derives from _g)

  static const T& _zeroref() { static const T z = T(0.); return z; }

  //! @brief in-place gradient scaling g *= fac  (used by scalar *=)
  void _scale_inplace( const T& fac )
  { if constexpr( MODE==0 ){ for(auto& g:_g) g*=fac; } else { for(auto& kv:_g) kv.second*=fac; } }
  //! @brief in-place gradient division g /= d  (used by scalar /=)
  void _div_inplace( const T& d )
  { if constexpr( MODE==0 ){ for(auto& g:_g) g/=d; } else { for(auto& kv:_g) kv.second/=d; } }

  //! @brief (sparse only) g = fa*A + fb*B via a linear sorted merge with
  //! end-hinted inserts -- amortised O(1) per key, no per-key find().
  //! Only instantiated for MODE==1 (the dense _chain2 branch never calls it).
  static container _combine( const container& A, const T& fa,
                             const container& B, const T& fb )
  {
    container R;
    if constexpr( MODE==1 ){
      auto ia=A.begin(), ib=B.begin(); const auto ea=A.end(), eb=B.end();
      while( ia!=ea && ib!=eb ){
        if( ia->first < ib->first ){ R.emplace_hint(R.end(), ia->first, fa*ia->second); ++ia; }
        else if( ib->first < ia->first ){ R.emplace_hint(R.end(), ib->first, fb*ib->second); ++ib; }
        else{ R.emplace_hint(R.end(), ia->first, fa*ia->second + fb*ib->second); ++ia; ++ib; }
      }
      for( ; ia!=ea; ++ia ) R.emplace_hint(R.end(), ia->first, fa*ia->second);
      for( ; ib!=eb; ++ib ) R.emplace_hint(R.end(), ib->first, fb*ib->second);
    }
    return R;
  }

public:
  //-------------------------------------------------------- construction ----//
  Fimpl(): _v(), _g(), _n(0) {}
  template <typename V> Fimpl( const V& val ): _v(val), _g(), _n(0) {}
  Fimpl( const Fimpl& ) = default;
  Fimpl( Fimpl&& ) noexcept = default;
  template <typename V> Fimpl& operator=( const V& val ) { _v=val; _g.clear(); _n=0; return *this; }
  Fimpl& operator=( const Fimpl& ) = default;
  Fimpl& operator=( Fimpl&& ) noexcept = default;

  //------------------------------------------------------------ accessors ---//
  const T& val() const { return _v; }
  T&       x()         { return _v; }
  const T& x()   const { return _v; }

  unsigned size()      const { return static_cast<unsigned>(_g.size()); }        //!< dense: dim; sparse: nnz
  unsigned dimension() const { if constexpr(MODE==0) return static_cast<unsigned>(_g.size()); else return _n; }
  bool     depend()    const { return !_g.empty(); }

  const container& gradient() const { return _g; }
  iterator       begin()       { return _g.begin(); }
  iterator       end()         { return _g.end();   }
  const_iterator begin() const { return _g.begin(); }
  const_iterator end()   const { return _g.end();   }

  const T& deriv( const unsigned i ) const
  {
    if constexpr( MODE==0 ){ return i<_g.size()? _g[i] : _zeroref(); }
    else { auto it=_g.find(i); return it==_g.end()? _zeroref() : it->second; }
  }
  T& d( const unsigned i )
  {
    if constexpr( MODE==0 ){ if(i<_g.size()) return _g[i]; static thread_local T z; z=T(0.); return z; }
    else return _g[i];
  }
  const T& operator[]( const unsigned i ) const { if constexpr(MODE==0) return _g[i]; else return deriv(i); }
  T&       operator[]( const unsigned i )       { return _g[i]; }

  //-------------------------------------------------------------- seeding ---//
  T& diff( const unsigned idx, const unsigned N )
  {
    if constexpr( MODE==0 ){
      if( _g.empty() ) _g.assign( N, T(0.) );
      else{ if( _g.size()!=N ) throw std::runtime_error("mc::F::diff -- inconsistent gradient dimension");
            std::fill( _g.begin(), _g.end(), T(0.) ); }
      _g[idx]=T(1.); return _g[idx];
    }
    else{ _n=N; _g.clear(); auto res=_g.emplace(idx,T(1.)); return res.first->second; }
  }

  //-------------------------------------------------- dependency management -//
  void setDepend( const Fimpl& v )
  {
    if constexpr( MODE==0 ){
      if( _g.empty() ) _g.assign( v._g.size(), T(0.) );
      else if( _g.size()!=v._g.size() ) throw std::runtime_error("mc::F::setDepend -- inconsistent gradient dimension");
    }
    else if( !_n ) _n=v._n;
  }
  void setDepend( const Fimpl& v1, const Fimpl& v2 )
  {
    if constexpr( MODE==0 ){
      if( v1._g.size()!=v2._g.size() ) throw std::runtime_error("mc::F::setDepend -- inconsistent gradient dimension");
      if( _g.empty() ) _g.assign( v1._g.size(), T(0.) );
      else if( _g.size()!=v1._g.size() ) throw std::runtime_error("mc::F::setDepend -- inconsistent gradient dimension");
    }
    else if( !_n ) _n=std::max(v1._n,v2._n);
  }

  //------------------------------------------------------- chain-rule core --//
  //! @brief unary chain rule: result = ( fval, fac * grad )
  Fimpl _chain( const T& fval, const T& fac ) const
  {
    Fimpl c( fval );
    if constexpr( MODE==0 ){
      if( _g.empty() ) return c;
      c._g.reserve( _g.size() );
      for( auto const& gi:_g ) c._g.emplace_back( fac*gi );
    }
    else{
      c._n=_n;
      for( auto const& kv:_g ) c._g.emplace_hint( c._g.end(), kv.first, fac*kv.second );
    }
    return c;
  }
  //! @brief binary chain rule: result = ( fval, fa*a.grad + fb*b.grad ); a
  //! missing (empty) gradient is treated as identically zero.
  static Fimpl _chain2( const T& fval, const Fimpl& a, const T& fa,
                                       const Fimpl& b, const T& fb )
  {
    Fimpl c( fval );
    if constexpr( MODE==0 ){
      const std::size_t na=a._g.size(), nb=b._g.size();
      if( na && nb ){
        if( na!=nb ) throw std::runtime_error("mc::F::_chain2 -- inconsistent gradient dimension");
        c._g.reserve(na); for(std::size_t i=0;i<na;++i) c._g.emplace_back( fa*a._g[i] + fb*b._g[i] );
      }
      else if( na ){ c._g.reserve(na); for(auto const& gi:a._g) c._g.emplace_back( fa*gi ); }
      else if( nb ){ c._g.reserve(nb); for(auto const& gi:b._g) c._g.emplace_back( fb*gi ); }
    }
    else{ c._n=std::max(a._n,b._n); c._g=_combine(a._g,fa,b._g,fb); }
    return c;
  }
  //! @brief division a/b keeping the common-denominator form (a_i - (a/b) b_i)/b,
  //! which stays tighter than the split a_i/b - (a/b^2) b_i under the
  //! interval/McCormick dependency problem; used for both storage modes.
  static Fimpl _divide( const Fimpl& a, const Fimpl& b )
  {
    const T cv = a._v / b._v;
    if( !a.depend() && !b.depend() ) return Fimpl( cv );
    if( !b.depend() ){ Fimpl c(a); c /= b._v; return c; }         // a_i / b
    Fimpl c( cv );
    if constexpr( MODE==0 ){
      if( a.depend() ){
        const std::size_t na=a._g.size(), nb=b._g.size();
        if( na!=nb ) throw std::runtime_error("mc::F::_divide -- inconsistent gradient dimension");
        c._g.reserve(na); for(std::size_t i=0;i<na;++i) c._g.emplace_back( (a._g[i] - cv*b._g[i]) / b._v );
      }
      else{ c._g.reserve(b._g.size()); for(auto const& gi:b._g) c._g.emplace_back( -(cv*gi)/b._v ); }
    }
    else{
      c._n=std::max(a._n,b._n);
      if( a.depend() ){
        auto ia=a._g.begin(), ib=b._g.begin(); const auto ea=a._g.end(), eb=b._g.end();
        while( ia!=ea && ib!=eb ){
          if( ia->first<ib->first ){ c._g.emplace_hint(c._g.end(), ia->first, ia->second/b._v); ++ia; }
          else if( ib->first<ia->first ){ c._g.emplace_hint(c._g.end(), ib->first, -(cv*ib->second)/b._v); ++ib; }
          else{ c._g.emplace_hint(c._g.end(), ia->first, (ia->second - cv*ib->second)/b._v); ++ia; ++ib; }
        }
        for( ; ia!=ea; ++ia ) c._g.emplace_hint(c._g.end(), ia->first, ia->second/b._v);
        for( ; ib!=eb; ++ib ) c._g.emplace_hint(c._g.end(), ib->first, -(cv*ib->second)/b._v);
      }
      else for(auto const& kv:b._g) c._g.emplace_hint(c._g.end(), kv.first, -(cv*kv.second)/b._v);
    }
    return c;
  }

  //---------------------------------------------------- compound (Fimpl) ----//
  Fimpl& operator+=( const Fimpl& o )
  {
    _v += o._v;
    if( !o.depend() ) return *this;
    if constexpr( MODE==0 ){
      if( depend() ) for(std::size_t i=0;i<_g.size();++i) _g[i]+=o._g[i];
      else _g=o._g;
    }
    else{
      auto hint=_g.begin();
      for( auto const& kv:o._g ){ hint=_g.lower_bound(kv.first);
        if( hint!=_g.end() && hint->first==kv.first ) hint->second=hint->second+kv.second;
        else _g.emplace_hint(hint,kv.first,kv.second); }
      _n=std::max(_n,o._n);
    }
    return *this;
  }
  Fimpl& operator-=( const Fimpl& o )
  {
    _v -= o._v;
    if( !o.depend() ) return *this;
    if constexpr( MODE==0 ){
      if( depend() ) for(std::size_t i=0;i<_g.size();++i) _g[i]-=o._g[i];
      else{ _g.reserve(o._g.size()); for(auto const& g:o._g) _g.emplace_back(-g); }
    }
    else{
      auto hint=_g.begin();
      for( auto const& kv:o._g ){ hint=_g.lower_bound(kv.first);
        if( hint!=_g.end() && hint->first==kv.first ) hint->second=hint->second-kv.second;
        else _g.emplace_hint(hint,kv.first,-kv.second); }
      _n=std::max(_n,o._n);
    }
    return *this;
  }
  Fimpl& operator*=( const Fimpl& o ) { *this = (*this) * o; return *this; }
  Fimpl& operator/=( const Fimpl& o ) { *this = (*this) / o; return *this; }
  template <typename V> Fimpl& operator+=( const V& c ) { _v += c; return *this; }
  template <typename V> Fimpl& operator-=( const V& c ) { _v -= c; return *this; }
  template <typename V> Fimpl& operator*=( const V& c ) { _v *= c; _scale_inplace(T(c)); return *this; }
  template <typename V> Fimpl& operator/=( const V& c ) { _v /= c; _div_inplace(T(c)); return *this; }
};

// ============================= arithmetic operators =========================
template <typename T,int M> inline Fimpl<T,M> operator+( const Fimpl<T,M>& a, const Fimpl<T,M>& b )
{ Fimpl<T,M> c(a); c+=b; return c; }
template <typename T,int M> inline Fimpl<T,M> operator-( const Fimpl<T,M>& a, const Fimpl<T,M>& b )
{ Fimpl<T,M> c(a); c-=b; return c; }
template <typename T,int M> inline Fimpl<T,M> operator*( const Fimpl<T,M>& a, const Fimpl<T,M>& b )
{
  if( !a.depend() && !b.depend() ) return Fimpl<T,M>( a.val()*b.val() );
  if( !b.depend() ){ Fimpl<T,M> c(a); c*=b.val(); return c; }
  if( !a.depend() ){ Fimpl<T,M> c(b); c*=a.val(); return c; }
  return Fimpl<T,M>::_chain2( a.val()*b.val(), a, b.val(), b, a.val() );
}
template <typename T,int M> inline Fimpl<T,M> operator/( const Fimpl<T,M>& a, const Fimpl<T,M>& b )
{ return Fimpl<T,M>::_divide( a, b ); }

template <typename T,int M,typename U> inline Fimpl<T,M> operator+( const Fimpl<T,M>& a, const U& b ) { Fimpl<T,M> c(a); c.x()+=b; return c; }
template <typename T,int M,typename U> inline Fimpl<T,M> operator+( const U& a, const Fimpl<T,M>& b ) { Fimpl<T,M> c(b); c.x()+=a; return c; }
template <typename T,int M,typename U> inline Fimpl<T,M> operator-( const Fimpl<T,M>& a, const U& b ) { Fimpl<T,M> c(a); c.x()-=b; return c; }
template <typename T,int M,typename U> inline Fimpl<T,M> operator-( const U& a, const Fimpl<T,M>& b ) { return b._chain( a - b.val(), T(-1.) ); }
template <typename T,int M,typename U> inline Fimpl<T,M> operator*( const Fimpl<T,M>& a, const U& b ) { Fimpl<T,M> c(a); c*=b; return c; }
template <typename T,int M,typename U> inline Fimpl<T,M> operator*( const U& a, const Fimpl<T,M>& b ) { Fimpl<T,M> c(b); c*=a; return c; }
template <typename T,int M,typename U> inline Fimpl<T,M> operator/( const Fimpl<T,M>& a, const U& b ) { Fimpl<T,M> c(a); c/=b; return c; }
template <typename T,int M,typename U> inline Fimpl<T,M> operator/( const U& a, const Fimpl<T,M>& b )
{ T cv=a/b.val(); return b._chain( cv, -( cv/b.val() ) ); }

template <typename T,int M> inline Fimpl<T,M> operator+( const Fimpl<T,M>& a ) { return a; }
template <typename T,int M> inline Fimpl<T,M> operator-( const Fimpl<T,M>& a ) { return a._chain( -a.val(), T(-1.) ); }

// ============================= elementary functions =========================
// Each rule: value + derivative multiplier via mc::Op<T>, propagated with the
// unary (_chain) or binary (_chain2) chain rule -- identical for dense/sparse.
template <typename T,int M> inline Fimpl<T,M> inv ( const Fimpl<T,M>& a ) { T fv=mc::Op<T>::inv(a.val());  return a._chain( fv, -mc::Op<T>::sqr(fv) ); }
template <typename T,int M> inline Fimpl<T,M> sqr ( const Fimpl<T,M>& a ) { return a._chain( mc::Op<T>::sqr(a.val()), 2.*a.val() ); }
template <typename T,int M> inline Fimpl<T,M> sqrt( const Fimpl<T,M>& a ) { T fv=mc::Op<T>::sqrt(a.val()); return a._chain( fv, mc::Op<T>::inv(2.*fv) ); }
template <typename T,int M> inline Fimpl<T,M> exp ( const Fimpl<T,M>& a ) { T fv=mc::Op<T>::exp(a.val());  return a._chain( fv, fv ); }
template <typename T,int M> inline Fimpl<T,M> log ( const Fimpl<T,M>& a ) { return a._chain( mc::Op<T>::log(a.val()), mc::Op<T>::inv(a.val()) ); }
template <typename T,int M> inline Fimpl<T,M> xlog( const Fimpl<T,M>& a ) { return a._chain( mc::Op<T>::xlog(a.val()), mc::Op<T>::log(a.val())+1. ); }
template <typename T,int M> inline Fimpl<T,M> sin ( const Fimpl<T,M>& a ) { return a._chain( mc::Op<T>::sin(a.val()), mc::Op<T>::cos(a.val()) ); }
template <typename T,int M> inline Fimpl<T,M> cos ( const Fimpl<T,M>& a ) { return a._chain( mc::Op<T>::cos(a.val()), -mc::Op<T>::sin(a.val()) ); }
template <typename T,int M> inline Fimpl<T,M> tan ( const Fimpl<T,M>& a ) { T fv=mc::Op<T>::tan(a.val());  return a._chain( fv, 1.+mc::Op<T>::sqr(fv) ); }
template <typename T,int M> inline Fimpl<T,M> asin( const Fimpl<T,M>& a ) { return a._chain( mc::Op<T>::asin(a.val()),  mc::Op<T>::inv(mc::Op<T>::sqrt(1.-mc::Op<T>::sqr(a.val()))) ); }
template <typename T,int M> inline Fimpl<T,M> acos( const Fimpl<T,M>& a ) { return a._chain( mc::Op<T>::acos(a.val()), -mc::Op<T>::inv(mc::Op<T>::sqrt(1.-mc::Op<T>::sqr(a.val()))) ); }
template <typename T,int M> inline Fimpl<T,M> atan( const Fimpl<T,M>& a ) { return a._chain( mc::Op<T>::atan(a.val()),  mc::Op<T>::inv(1.+mc::Op<T>::sqr(a.val())) ); }
template <typename T,int M> inline Fimpl<T,M> sinh( const Fimpl<T,M>& a ) { return a._chain( mc::Op<T>::sinh(a.val()), mc::Op<T>::cosh(a.val()) ); }
template <typename T,int M> inline Fimpl<T,M> cosh( const Fimpl<T,M>& a ) { return a._chain( mc::Op<T>::cosh(a.val()), mc::Op<T>::sinh(a.val()) ); }
template <typename T,int M> inline Fimpl<T,M> tanh( const Fimpl<T,M>& a ) { T fv=mc::Op<T>::tanh(a.val()); return a._chain( fv, 1.-mc::Op<T>::sqr(fv) ); }
template <typename T,int M> inline Fimpl<T,M> fabs( const Fimpl<T,M>& a ) { return a._chain( mc::Op<T>::fabs(a.val()), 2.*mc::Op<T>::fstep(a.val())-1. ); }

template <typename T,int M> inline Fimpl<T,M> erf( const Fimpl<T,M>& a )
{ static const double c=2./std::sqrt(mc::PI); return a._chain( mc::Op<T>::erf(a.val()), c*mc::Op<T>::exp(-mc::Op<T>::sqr(a.val())) ); }
template <typename T,int M> inline Fimpl<T,M> erfc( const Fimpl<T,M>& a ) { return 1. - erf(a); }

template <typename T,int M> inline Fimpl<T,M> cheb( const Fimpl<T,M>& a, const unsigned n )
{
  T fv = mc::Op<T>::cheb( a.val(), n );
  // T_n'(x) = 2n * sum_{j=n-1,n-3,...} T_j(x)  (+ n T_0 for odd n), accumulated
  // in one O(n) Chebyshev recurrence rather than O(n^2) re-evaluation.
  const unsigned par = (n-1) & 1u;
  T fac( n & 1u ? 0.5 : 0. );
  if( n >= 2 ){
    const T& x = a.val();
    T Tkm2=T(1.), Tkm1=x;
    if( (1u & 1u)==par ) fac += Tkm1;
    for( unsigned k=2; k<n; ++k ){ T Tk=2.*x*Tkm1-Tkm2; if( (k&1u)==par ) fac+=Tk; Tkm2=Tkm1; Tkm1=Tk; }
  }
  fac *= 2.*(double)n;
  return a._chain( fv, fac );
}

template <typename T,int M,typename U> inline Fimpl<T,M> pow( const Fimpl<T,M>& a, const U& b )
{ T fv=mc::Op<T>::pow(a.val(),b); return a._chain( fv, b*mc::Op<T>::pow(a.val(),b-1.) ); }
template <typename T,int M> inline Fimpl<T,M> pow( const Fimpl<T,M>& a, const int b )
{ T fv=mc::Op<T>::pow(a.val(),b); return a._chain( fv, (double)b*mc::Op<T>::pow(a.val(),b-1) ); }
template <typename T,int M,typename U> inline Fimpl<T,M> pow( const U& a, const Fimpl<T,M>& b )
{ T fv=mc::Op<T>::pow(a,b.val()); return b._chain( fv, fv*std::log((double)a) ); }
template <typename T,int M> inline Fimpl<T,M> pow( const Fimpl<T,M>& a, const Fimpl<T,M>& b )
{
  if( !b.depend() ) return pow( a, b.val() );
  if( !a.depend() ) return pow( a.val(), b );
  T fv = mc::Op<T>::pow( a.val(), b.val() );
  return Fimpl<T,M>::_chain2( fv, a, b.val()*mc::Op<T>::pow(a.val(),b.val()-1.), b, fv*mc::Op<T>::log(a.val()) );
}

template <typename T,int M,typename U> inline Fimpl<T,M> max( const Fimpl<T,M>& a, const U& b ) { return a._chain( mc::Op<T>::max(a.val(),b), mc::Op<T>::fstep(a.val()-b) ); }
template <typename T,int M,typename U> inline Fimpl<T,M> max( const U& a, const Fimpl<T,M>& b ) { return max(b,a); }
template <typename T,int M> inline Fimpl<T,M> max( const Fimpl<T,M>& a, const Fimpl<T,M>& b )
{ T w=mc::Op<T>::fstep(b.val()-a.val()); return Fimpl<T,M>::_chain2( mc::Op<T>::max(a.val(),b.val()), a, 1.-w, b, w ); }
template <typename T,int M,typename U> inline Fimpl<T,M> min( const Fimpl<T,M>& a, const U& b ) { return a._chain( mc::Op<T>::min(a.val(),b), mc::Op<T>::fstep(b-a.val()) ); }
template <typename T,int M,typename U> inline Fimpl<T,M> min( const U& a, const Fimpl<T,M>& b ) { return min(b,a); }
template <typename T,int M> inline Fimpl<T,M> min( const Fimpl<T,M>& a, const Fimpl<T,M>& b )
{ T w=mc::Op<T>::fstep(a.val()-b.val()); return Fimpl<T,M>::_chain2( mc::Op<T>::min(a.val(),b.val()), a, 1.-w, b, w ); }

template <typename T,int M> inline Fimpl<T,M> lmtd( const Fimpl<T,M>& a, const Fimpl<T,M>& b )
{
  if( mc::Op<T>::eq( a.val(), b.val() ) ) return Fimpl<T,M>::_chain2( a.val(), a, T(0.5), b, T(0.5) );
  return ( a - b ) / ( log(a) - log(b) );
}
template <typename T,int M> inline Fimpl<T,M> rlmtd( const Fimpl<T,M>& a, const Fimpl<T,M>& b )
{
  if( mc::Op<T>::eq( a.val(), b.val() ) )
    return Fimpl<T,M>::_chain2( mc::Op<T>::inv(a.val()), a, -mc::Op<T>::inv(2.*mc::Op<T>::sqr(a.val())),
                                                         b, -mc::Op<T>::inv(2.*mc::Op<T>::sqr(b.val())) );
  return ( log(a) - log(b) ) / ( a - b );
}
template <typename T,int M> inline Fimpl<T,M> arh( const Fimpl<T,M>& a, const double k ) { return exp( -k/a ); }

// ============================== comparisons =================================
template <typename T,int M> inline bool operator==( const Fimpl<T,M>& a, const Fimpl<T,M>& b ) { return a.val()==b.val(); }
template <typename T,int M> inline bool operator!=( const Fimpl<T,M>& a, const Fimpl<T,M>& b ) { return a.val()!=b.val(); }
template <typename T,int M> inline bool operator< ( const Fimpl<T,M>& a, const Fimpl<T,M>& b ) { return a.val()< b.val(); }
template <typename T,int M> inline bool operator<=( const Fimpl<T,M>& a, const Fimpl<T,M>& b ) { return a.val()<=b.val(); }
template <typename T,int M> inline bool operator> ( const Fimpl<T,M>& a, const Fimpl<T,M>& b ) { return a.val()> b.val(); }
template <typename T,int M> inline bool operator>=( const Fimpl<T,M>& a, const Fimpl<T,M>& b ) { return a.val()>=b.val(); }
template <typename T,int M,typename U> inline bool operator==( const Fimpl<T,M>& a, const U& b ) { return a.val()==b; }
template <typename T,int M,typename U> inline bool operator==( const U& a, const Fimpl<T,M>& b ) { return a==b.val(); }
template <typename T,int M,typename U> inline bool operator!=( const Fimpl<T,M>& a, const U& b ) { return a.val()!=b; }
template <typename T,int M,typename U> inline bool operator!=( const U& a, const Fimpl<T,M>& b ) { return a!=b.val(); }
template <typename T,int M,typename U> inline bool operator< ( const Fimpl<T,M>& a, const U& b ) { return a.val()< b; }
template <typename T,int M,typename U> inline bool operator< ( const U& a, const Fimpl<T,M>& b ) { return a< b.val(); }
template <typename T,int M,typename U> inline bool operator<=( const Fimpl<T,M>& a, const U& b ) { return a.val()<=b; }
template <typename T,int M,typename U> inline bool operator<=( const U& a, const Fimpl<T,M>& b ) { return a<=b.val(); }
template <typename T,int M,typename U> inline bool operator> ( const Fimpl<T,M>& a, const U& b ) { return a.val()> b; }
template <typename T,int M,typename U> inline bool operator> ( const U& a, const Fimpl<T,M>& b ) { return a> b.val(); }
template <typename T,int M,typename U> inline bool operator>=( const Fimpl<T,M>& a, const U& b ) { return a.val()>=b; }
template <typename T,int M,typename U> inline bool operator>=( const U& a, const Fimpl<T,M>& b ) { return a>=b.val(); }

// =============================== public aliases =============================
template <typename T> using F  = mc::Fimpl<T,0>;   //!< dense forward-mode type
template <typename T> using SF = mc::Fimpl<T,1>;   //!< sparse forward-mode type

// =================== mc::Op specialization (covers F and SF) =================
template <typename U,int M> struct Op< mc::Fimpl<U,M> >
{
  typedef mc::Fimpl<U,M> TU;
  static TU point( const double c ) { throw std::runtime_error("mc::Op<mc::F<U>>::point -- operation not permitted"); }
  static TU zeroone() { throw std::runtime_error("mc::Op<mc::F<U>>::zeroone -- operation not permitted"); }
  static void I(TU& x, const TU& y) { x = y; }
  static double l(const TU& x) { throw std::runtime_error("mc::Op<mc::F<U>>::l -- operation not permitted"); }
  static double u(const TU& x) { throw std::runtime_error("mc::Op<mc::F<U>>::u -- operation not permitted"); }
  static double abs (const TU& x) { throw std::runtime_error("mc::Op<mc::F<U>>::abs -- operation not permitted"); }
  static double mid (const TU& x) { throw std::runtime_error("mc::Op<mc::F<U>>::mid -- operation not permitted"); }
  static double diam(const TU& x) { throw std::runtime_error("mc::Op<mc::F<U>>::diam -- operation not permitted"); }
  static TU inv (const TU& x) { return mc::inv(x);  }
  static TU sqr (const TU& x) { return mc::sqr(x);  }
  static TU sqrt(const TU& x) { return mc::sqrt(x); }
  static TU exp (const TU& x) { return mc::exp(x);  }
  static TU log (const TU& x) { return mc::log(x);  }
  static TU xlog(const TU& x) { return mc::xlog(x); }
  static TU lmtd(const TU& x, const TU& y) { return mc::lmtd(x,y); }
  static TU rlmtd(const TU& x, const TU& y) { return mc::rlmtd(x,y); }
  static TU fabs(const TU& x) { return mc::fabs(x); }
  static TU sin (const TU& x) { return mc::sin(x);  }
  static TU cos (const TU& x) { return mc::cos(x);  }
  static TU tan (const TU& x) { return mc::tan(x);  }
  static TU asin(const TU& x) { return mc::asin(x); }
  static TU acos(const TU& x) { return mc::acos(x); }
  static TU atan(const TU& x) { return mc::atan(x); }
  static TU sinh(const TU& x) { return mc::sinh(x); }
  static TU cosh(const TU& x) { return mc::cosh(x); }
  static TU tanh(const TU& x) { return mc::tanh(x); }
  static TU erf (const TU& x) { return mc::erf(x);  }
  static TU erfc(const TU& x) { return mc::erfc(x); }
  static TU fstep(const TU& x) { throw std::runtime_error("mc::Op<mc::F<U>>::fstep -- operation not permitted"); }
  static TU bstep(const TU& x) { throw std::runtime_error("mc::Op<mc::F<U>>::bstep -- operation not permitted"); }
  static TU hull(const TU& x, const TU& y) { throw std::runtime_error("mc::Op<mc::F<U>>::hull -- operation not permitted"); }
  template <typename Y> static TU min(const TU& x, const Y& y) { return mc::min(x,y); }
  template <typename Y> static TU max(const TU& x, const Y& y) { return mc::max(x,y); }
  static TU arh (const TU& x, const double k) { return mc::arh(x,k); }
  template <typename X, typename Y> static TU pow(const X& x, const Y& y) { return mc::pow(x,y); }
  static TU cheb(const TU& x, const unsigned n) { return mc::cheb(x,n); }
  static TU prod (const unsigned n, const TU* x) { switch(n){ case 0: return TU(1.); case 1: return x[0]; default: return x[0]*prod(n-1,x+1); } }
  static TU monom (const unsigned n, const TU* x, const unsigned* k) { switch(n){ case 0: return TU(1.); case 1: return pow(x[0],(int)k[0]); default: return pow(x[0],(int)k[0])*monom(n-1,x+1,k+1); } }
  static bool inter(TU& xIy, const TU& x, const TU& y) { xIy = x; return true; }
  static bool eq(const TU& x, const TU& y) { throw std::runtime_error("mc::Op<mc::F<U>>::eq -- operation not permitted"); }
  static bool ne(const TU& x, const TU& y) { throw std::runtime_error("mc::Op<mc::F<U>>::ne -- operation not permitted"); }
  static bool lt(const TU& x, const TU& y) { throw std::runtime_error("mc::Op<mc::F<U>>::lt -- operation not permitted"); }
  static bool le(const TU& x, const TU& y) { throw std::runtime_error("mc::Op<mc::F<U>>::le -- operation not permitted"); }
  static bool gt(const TU& x, const TU& y) { throw std::runtime_error("mc::Op<mc::F<U>>::gt -- operation not permitted"); }
  static bool ge(const TU& x, const TU& y) { throw std::runtime_error("mc::Op<mc::F<U>>::ge -- operation not permitted"); }
};

} // namespace mc

#endif // MC__FDIFF_HPP
