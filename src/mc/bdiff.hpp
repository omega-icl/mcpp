// bdiff.hpp -- reverse-mode (adjoint) automatic differentiation for MC++.
//
// Provides mc::B<T>  : dense reverse-mode type (adjoint per output seed stored
//                      in a std::vector), and
//          mc::SB<T> : sparse reverse-mode type (adjoint stored as
//                      (output-seed index -> derivative) pairs in a std::map).
//
// Both are templated on the evaluation type T and agree with mc::Op<T>, and
// expose the same public interface as fadbad::B (val/x, diff, d/deriv,
// compound and free operators, elementary functions).  The implementation is
// an independent, from-scratch reverse-mode design: each operation records a
// node in a reference-counted (std::shared_ptr) computation DAG together with
// the local partial derivatives of that operation; seeding an output with
// diff(idx,n) runs an explicit reverse-topological accumulation sweep so that
// d(i) on any input returns d(output_i)/d(input).  It does not reproduce the
// internal machinery of FADBAD++ -- only the public user interface.
//
// The local partial-derivative forms mirror those in fdiff.hpp (mc::F), which
// were chosen/verified to be favorable for interval and McCormick enclosure
// propagation (e.g. d(a/b)/db = -(a/b)/b rather than -a/b^2, inv' = -sqr).

#ifndef MC__BDIFF_HPP
#define MC__BDIFF_HPP

#include <memory>
#include <vector>
#include <map>
#include <stdexcept>
#include <cmath>
#include "mcop.hpp"

namespace mc
{

namespace bdiff_detail
{

// MODE = 0 : dense adjoint (std::vector<T>), indexed by output seed
// MODE = 1 : sparse adjoint (std::map<unsigned,T>), keyed by output seed
template <typename T, int MODE>
class Bimpl
{
  // dense adjoint keeps slot 0 inline (small-buffer): the common single-output
  // gradient then needs no per-node heap allocation at all.
  struct SBOAdj { T a0{}; std::vector<T> rest; };
  typedef typename std::conditional<MODE==0, SBOAdj,
                                    std::map<unsigned,T> >::type AdjT;

  struct Node
  {
    T                                                   _v;        // value
    std::vector< std::pair<std::shared_ptr<Node>, T> >  _edge;     // (operand, d(this)/d(operand))
    AdjT                                                _adj;      // accumulated adjoint(s)
    unsigned                                            _mark = 0; // DFS epoch stamp
    explicit Node( const T& v ) : _v(v) {}
  };

  std::shared_ptr<Node> _n;

  explicit Bimpl( const std::shared_ptr<Node>& nd ) : _n(nd) {}

  static unsigned& _epoch() { static thread_local unsigned e = 0; return e; }

  // ---- adjoint-storage helpers (dense SBO / sparse map via if constexpr) ----
  static T& _slot( AdjT& adj, unsigned idx )
  {
    if constexpr( MODE==0 ){ if( idx==0 ) return adj.a0;
      if( adj.rest.size()<idx ) adj.rest.resize( idx, T(0.) ); return adj.rest[idx-1]; }
    else return adj[idx];
  }
  static void _reset( AdjT& adj, unsigned idx )
  {
    if constexpr( MODE==0 ){ if( idx==0 ) adj.a0 = T(0.);
      else if( idx-1<adj.rest.size() ) adj.rest[idx-1] = T(0.); }
    else adj.erase(idx);
  }
  static void _accum( AdjT& adj, unsigned idx, const T& v )
  {
    if constexpr( MODE==0 ){ _slot(adj,idx) += v; }
    else { auto it=adj.find(idx); if( it==adj.end() ) adj.emplace(idx,v); else it->second += v; }
  }
  static T _read( const AdjT& adj, unsigned idx )
  {
    if constexpr( MODE==0 ){ if( idx==0 ) return adj.a0;
      return idx-1<adj.rest.size()? adj.rest[idx-1] : T(0.); }
    else { auto it=adj.find(idx); return it==adj.end()? T(0.) : it->second; }
  }

  // ---- reverse-topological accumulation sweep for output seed idx ----
  static void _backward( const std::shared_ptr<Node>& out, unsigned idx, unsigned /*n*/ )
  {
    // reusable scratch (single-threaded); cleared, not freed, between sweeps
    static thread_local std::vector<Node*> order;
    static thread_local std::vector< std::pair<Node*,std::size_t> > stack;
    order.clear(); stack.clear();
    const unsigned e = ++_epoch();
    // iterative post-order DFS (operands before self); epoch stamp = visited set
    stack.emplace_back( out.get(), 0 );
    out->_mark = e;
    while( !stack.empty() ){
      auto& top = stack.back();
      if( top.second < top.first->_edge.size() ){
        Node* ch = top.first->_edge[top.second].first.get();
        ++top.second;
        if( ch->_mark != e ){ ch->_mark = e; stack.emplace_back( ch, 0 ); }
      }
      else{ order.push_back( top.first ); stack.pop_back(); }
    }
    // clear slot idx on every reachable node, then seed the output with 1
    for( Node* p : order ) _reset( p->_adj, idx );
    _slot( out->_adj, idx ) = T(1.);
    // process in reverse post-order (self before operands => topological)
    for( auto it=order.rbegin(); it!=order.rend(); ++it ){
      Node* p = *it;
      if constexpr( MODE==0 ){
        const T a = _read( p->_adj, idx );
        for( auto& ed : p->_edge ) _accum( ed.first->_adj, idx, ed.second * a );
      }
      else{
        auto ai = p->_adj.find(idx);
        if( ai==p->_adj.end() ) continue;          // no adjoint reached this node
        const T a = ai->second;
        for( auto& ed : p->_edge ) _accum( ed.first->_adj, idx, ed.second * a );
      }
    }
  }

public:
  // ---------------- construction / assignment (handle semantics) ------------
  Bimpl() : _n( std::make_shared<Node>( T() ) ) {}
  template <typename V> Bimpl( const V& v ) : _n( std::make_shared<Node>( T(v) ) ) {}
  Bimpl( const Bimpl& ) = default;                 // shares the node
  Bimpl& operator=( const Bimpl& ) = default;      // rebinds to shared node
  template <typename V> Bimpl& operator=( const V& v )
  { _n = std::make_shared<Node>( T(v) ); return *this; }

  // ---------------- value access --------------------------------------------
  const T& val() const { return _n->_v; }
  T&       val()       { return _n->_v; }
  const T& x()   const { return _n->_v; }
  T&       x()         { return _n->_v; }

  // ---------------- derivative access / seeding -----------------------------
  //! @brief seed this node as output idx of n and back-propagate; returns seed slot
  T& diff( const unsigned idx, const unsigned n )
  { _backward( _n, idx, n ); return _slot( _n->_adj, idx ); }
  //! @brief d(output i)/d(this) after propagation (0 if this does not affect output i)
  T d    ( const unsigned i ) const { return _read( _n->_adj, i ); }
  T deriv( const unsigned i ) const { return _read( _n->_adj, i ); }

  // ---------------- node builders (used by the free functions) --------------
  static Bimpl _unary( const T& fv, const Bimpl& a, const T& fa )
  {
    auto nd = std::make_shared<Node>( fv );
    nd->_edge.reserve( 1 );
    nd->_edge.emplace_back( a._n, fa );
    return Bimpl( nd );
  }
  static Bimpl _binary( const T& fv, const Bimpl& a, const T& fa,
                                     const Bimpl& b, const T& fb )
  {
    auto nd = std::make_shared<Node>( fv );
    nd->_edge.reserve( 2 );
    nd->_edge.emplace_back( a._n, fa );
    nd->_edge.emplace_back( b._n, fb );
    return Bimpl( nd );
  }

  // ---------------- compound assignment -------------------------------------
  Bimpl& operator+=( const Bimpl& b ){ *this = *this + b; return *this; }
  Bimpl& operator-=( const Bimpl& b ){ *this = *this - b; return *this; }
  Bimpl& operator*=( const Bimpl& b ){ *this = *this * b; return *this; }
  Bimpl& operator/=( const Bimpl& b ){ *this = *this / b; return *this; }
  template <typename V> Bimpl& operator+=( const V& s ){ *this = *this + s; return *this; }
  template <typename V> Bimpl& operator-=( const V& s ){ *this = *this - s; return *this; }
  template <typename V> Bimpl& operator*=( const V& s ){ *this = *this * s; return *this; }
  template <typename V> Bimpl& operator/=( const V& s ){ *this = *this / s; return *this; }
};

// ============================ arithmetic operators ==========================
template <typename T, int M> inline Bimpl<T,M> operator+( const Bimpl<T,M>& a ) { return a; }
template <typename T, int M> inline Bimpl<T,M> operator-( const Bimpl<T,M>& a )
{ return Bimpl<T,M>::_unary( -a.val(), a, T(-1.) ); }

template <typename T, int M> inline Bimpl<T,M> operator+( const Bimpl<T,M>& a, const Bimpl<T,M>& b )
{ return Bimpl<T,M>::_binary( a.val()+b.val(), a, T(1.), b, T(1.) ); }
template <typename T, int M> inline Bimpl<T,M> operator-( const Bimpl<T,M>& a, const Bimpl<T,M>& b )
{ return Bimpl<T,M>::_binary( a.val()-b.val(), a, T(1.), b, T(-1.) ); }
template <typename T, int M> inline Bimpl<T,M> operator*( const Bimpl<T,M>& a, const Bimpl<T,M>& b )
{ return Bimpl<T,M>::_binary( a.val()*b.val(), a, b.val(), b, a.val() ); }
template <typename T, int M> inline Bimpl<T,M> operator/( const Bimpl<T,M>& a, const Bimpl<T,M>& b )
{
  const T fv = a.val()/b.val();
  // d/da = 1/b ; d/db = -(a/b)/b  (McCormick-favorable form, see fdiff.hpp)
  return Bimpl<T,M>::_binary( fv, a, mc::Op<T>::inv(b.val()), b, -( fv/b.val() ) );
}

template <typename T, int M, typename V> inline Bimpl<T,M> operator+( const Bimpl<T,M>& a, const V& s )
{ return Bimpl<T,M>::_unary( a.val()+s, a, T(1.) ); }
template <typename T, int M, typename V> inline Bimpl<T,M> operator+( const V& s, const Bimpl<T,M>& a )
{ return Bimpl<T,M>::_unary( s+a.val(), a, T(1.) ); }
template <typename T, int M, typename V> inline Bimpl<T,M> operator-( const Bimpl<T,M>& a, const V& s )
{ return Bimpl<T,M>::_unary( a.val()-s, a, T(1.) ); }
template <typename T, int M, typename V> inline Bimpl<T,M> operator-( const V& s, const Bimpl<T,M>& a )
{ return Bimpl<T,M>::_unary( s-a.val(), a, T(-1.) ); }
template <typename T, int M, typename V> inline Bimpl<T,M> operator*( const Bimpl<T,M>& a, const V& s )
{ return Bimpl<T,M>::_unary( a.val()*s, a, T(s) ); }
template <typename T, int M, typename V> inline Bimpl<T,M> operator*( const V& s, const Bimpl<T,M>& a )
{ return Bimpl<T,M>::_unary( s*a.val(), a, T(s) ); }
template <typename T, int M, typename V> inline Bimpl<T,M> operator/( const Bimpl<T,M>& a, const V& s )
{ return Bimpl<T,M>::_unary( a.val()/s, a, mc::Op<T>::inv( T(s) ) ); }
template <typename T, int M, typename V> inline Bimpl<T,M> operator/( const V& s, const Bimpl<T,M>& a )
{
  const T fv = s/a.val();
  return Bimpl<T,M>::_unary( fv, a, -( fv/a.val() ) );   // -(s/a)/a
}

// ============================ comparison operators ==========================
template <typename T, int M> inline bool operator==( const Bimpl<T,M>& a, const Bimpl<T,M>& b ){ return mc::Op<T>::eq(a.val(),b.val()); }
template <typename T, int M> inline bool operator!=( const Bimpl<T,M>& a, const Bimpl<T,M>& b ){ return mc::Op<T>::ne(a.val(),b.val()); }
template <typename T, int M> inline bool operator<=( const Bimpl<T,M>& a, const Bimpl<T,M>& b ){ return mc::Op<T>::le(a.val(),b.val()); }
template <typename T, int M> inline bool operator>=( const Bimpl<T,M>& a, const Bimpl<T,M>& b ){ return mc::Op<T>::ge(a.val(),b.val()); }
template <typename T, int M> inline bool operator<( const Bimpl<T,M>& a, const Bimpl<T,M>& b ){ return mc::Op<T>::lt(a.val(),b.val()); }
template <typename T, int M> inline bool operator>( const Bimpl<T,M>& a, const Bimpl<T,M>& b ){ return mc::Op<T>::gt(a.val(),b.val()); }

// ============================ elementary functions ==========================
template <typename T, int M> inline Bimpl<T,M> inv( const Bimpl<T,M>& a )
{ T fv=mc::Op<T>::inv(a.val()); return Bimpl<T,M>::_unary( fv, a, -mc::Op<T>::sqr(fv) ); }
template <typename T, int M> inline Bimpl<T,M> sqr( const Bimpl<T,M>& a )
{ return Bimpl<T,M>::_unary( mc::Op<T>::sqr(a.val()), a, 2.*a.val() ); }
template <typename T, int M> inline Bimpl<T,M> sqrt( const Bimpl<T,M>& a )
{ T fv=mc::Op<T>::sqrt(a.val()); return Bimpl<T,M>::_unary( fv, a, mc::Op<T>::inv(2.*fv) ); }
template <typename T, int M> inline Bimpl<T,M> exp( const Bimpl<T,M>& a )
{ T fv=mc::Op<T>::exp(a.val()); return Bimpl<T,M>::_unary( fv, a, fv ); }
template <typename T, int M> inline Bimpl<T,M> log( const Bimpl<T,M>& a )
{ return Bimpl<T,M>::_unary( mc::Op<T>::log(a.val()), a, mc::Op<T>::inv(a.val()) ); }
template <typename T, int M> inline Bimpl<T,M> xlog( const Bimpl<T,M>& a )
{ return Bimpl<T,M>::_unary( mc::Op<T>::xlog(a.val()), a, mc::Op<T>::log(a.val())+1. ); }
template <typename T, int M> inline Bimpl<T,M> sin( const Bimpl<T,M>& a )
{ return Bimpl<T,M>::_unary( mc::Op<T>::sin(a.val()), a, mc::Op<T>::cos(a.val()) ); }
template <typename T, int M> inline Bimpl<T,M> cos( const Bimpl<T,M>& a )
{ return Bimpl<T,M>::_unary( mc::Op<T>::cos(a.val()), a, -mc::Op<T>::sin(a.val()) ); }
template <typename T, int M> inline Bimpl<T,M> tan( const Bimpl<T,M>& a )
{ T fv=mc::Op<T>::tan(a.val()); return Bimpl<T,M>::_unary( fv, a, 1.+mc::Op<T>::sqr(fv) ); }
template <typename T, int M> inline Bimpl<T,M> asin( const Bimpl<T,M>& a )
{ return Bimpl<T,M>::_unary( mc::Op<T>::asin(a.val()), a, mc::Op<T>::inv(mc::Op<T>::sqrt(1.-mc::Op<T>::sqr(a.val()))) ); }
template <typename T, int M> inline Bimpl<T,M> acos( const Bimpl<T,M>& a )
{ return Bimpl<T,M>::_unary( mc::Op<T>::acos(a.val()), a, -mc::Op<T>::inv(mc::Op<T>::sqrt(1.-mc::Op<T>::sqr(a.val()))) ); }
template <typename T, int M> inline Bimpl<T,M> atan( const Bimpl<T,M>& a )
{ return Bimpl<T,M>::_unary( mc::Op<T>::atan(a.val()), a, mc::Op<T>::inv(1.+mc::Op<T>::sqr(a.val())) ); }
template <typename T, int M> inline Bimpl<T,M> sinh( const Bimpl<T,M>& a )
{ return Bimpl<T,M>::_unary( mc::Op<T>::sinh(a.val()), a, mc::Op<T>::cosh(a.val()) ); }
template <typename T, int M> inline Bimpl<T,M> cosh( const Bimpl<T,M>& a )
{ return Bimpl<T,M>::_unary( mc::Op<T>::cosh(a.val()), a, mc::Op<T>::sinh(a.val()) ); }
template <typename T, int M> inline Bimpl<T,M> tanh( const Bimpl<T,M>& a )
{ T fv=mc::Op<T>::tanh(a.val()); return Bimpl<T,M>::_unary( fv, a, 1.-mc::Op<T>::sqr(fv) ); }
template <typename T, int M> inline Bimpl<T,M> fabs( const Bimpl<T,M>& a )
{ return Bimpl<T,M>::_unary( mc::Op<T>::fabs(a.val()), a, 2.*mc::Op<T>::fstep(a.val())-1. ); }
template <typename T, int M> inline Bimpl<T,M> erf( const Bimpl<T,M>& a )
{ static const double c=2./std::sqrt(mc::PI);
  return Bimpl<T,M>::_unary( mc::Op<T>::erf(a.val()), a, c*mc::Op<T>::exp(-mc::Op<T>::sqr(a.val())) ); }
template <typename T, int M> inline Bimpl<T,M> erfc( const Bimpl<T,M>& a )
{ static const double c=2./std::sqrt(mc::PI);
  return Bimpl<T,M>::_unary( mc::Op<T>::erfc(a.val()), a, -c*mc::Op<T>::exp(-mc::Op<T>::sqr(a.val())) ); }

template <typename T, int M> inline Bimpl<T,M> cheb( const Bimpl<T,M>& a, const unsigned n )
{
  T fv = mc::Op<T>::cheb( a.val(), n );
  // T_n'(x) accumulated via a single O(n) Chebyshev recurrence (see fdiff.hpp)
  const unsigned par = (n-1) & 1u;
  T fac( n & 1u ? 0.5 : 0. );
  if( n >= 2 ){
    const T& x = a.val();
    T Tkm2 = T(1.), Tkm1 = x;
    if( (1u & 1u) == par ) fac += Tkm1;
    for( unsigned k=2; k<n; ++k ){
      T Tk = 2.*x*Tkm1 - Tkm2;
      if( (k & 1u) == par ) fac += Tk;
      Tkm2 = Tkm1; Tkm1 = Tk;
    }
  }
  fac *= 2.*(double)n;
  return Bimpl<T,M>::_unary( fv, a, fac );
}

// ---- powers ----
template <typename T, int M> inline Bimpl<T,M> pow( const Bimpl<T,M>& a, const int n )
{ return Bimpl<T,M>::_unary( mc::Op<T>::pow(a.val(),n), a, (double)n*mc::Op<T>::pow(a.val(),n-1) ); }
template <typename T, int M> inline Bimpl<T,M> pow( const Bimpl<T,M>& a, const double r )
{ return Bimpl<T,M>::_unary( mc::Op<T>::pow(a.val(),r), a, r*mc::Op<T>::pow(a.val(),r-1.) ); }
template <typename T, int M> inline Bimpl<T,M> pow( const Bimpl<T,M>& a, const Bimpl<T,M>& b )
{
  T fv = mc::Op<T>::pow( a.val(), b.val() );
  T fa = b.val() * mc::Op<T>::pow( a.val(), b.val()-1. );
  T fb = fv * mc::Op<T>::log( a.val() );
  return Bimpl<T,M>::_binary( fv, a, fa, b, fb );
}
template <typename T, int M, typename V> inline Bimpl<T,M> pow( const V& s, const Bimpl<T,M>& a )
{ T fv=mc::Op<T>::pow(T(s),a.val()); return Bimpl<T,M>::_unary( fv, a, fv*mc::Op<T>::log(T(s)) ); }

// ---- min / max ----
template <typename T, int M> inline Bimpl<T,M> min( const Bimpl<T,M>& a, const Bimpl<T,M>& b )
{ T w=mc::Op<T>::fstep(a.val()-b.val()); return Bimpl<T,M>::_binary( mc::Op<T>::min(a.val(),b.val()), a, 1.-w, b, w ); }
template <typename T, int M> inline Bimpl<T,M> max( const Bimpl<T,M>& a, const Bimpl<T,M>& b )
{ T w=mc::Op<T>::fstep(b.val()-a.val()); return Bimpl<T,M>::_binary( mc::Op<T>::max(a.val(),b.val()), a, 1.-w, b, w ); }
template <typename T, int M, typename V> inline Bimpl<T,M> min( const Bimpl<T,M>& a, const V& s )
{ return Bimpl<T,M>::_unary( mc::Op<T>::min(a.val(),T(s)), a, 1.-mc::Op<T>::fstep(a.val()-s) ); }
template <typename T, int M, typename V> inline Bimpl<T,M> min( const V& s, const Bimpl<T,M>& a )
{ return min( a, s ); }
template <typename T, int M, typename V> inline Bimpl<T,M> max( const Bimpl<T,M>& a, const V& s )
{ return Bimpl<T,M>::_unary( mc::Op<T>::max(a.val(),T(s)), a, mc::Op<T>::fstep(a.val()-s) ); }
template <typename T, int M, typename V> inline Bimpl<T,M> max( const V& s, const Bimpl<T,M>& a )
{ return max( a, s ); }

// ---- lmtd / rlmtd / arh (built from primitives where non-singular) ----
template <typename T, int M> inline Bimpl<T,M> lmtd( const Bimpl<T,M>& a, const Bimpl<T,M>& b )
{
  if( mc::Op<T>::eq( a.val(), b.val() ) )                       // removable singularity: value a, d = 1/2 each
    return Bimpl<T,M>::_binary( a.val(), a, T(0.5), b, T(0.5) );
  return ( a - b ) / ( log(a) - log(b) );                       // graph handles derivatives
}
template <typename T, int M> inline Bimpl<T,M> rlmtd( const Bimpl<T,M>& a, const Bimpl<T,M>& b )
{
  if( mc::Op<T>::eq( a.val(), b.val() ) )
    return Bimpl<T,M>::_binary( mc::Op<T>::inv(a.val()),
                                a, -mc::Op<T>::inv( 2.*mc::Op<T>::sqr(a.val()) ),
                                b, -mc::Op<T>::inv( 2.*mc::Op<T>::sqr(b.val()) ) );
  return ( log(a) - log(b) ) / ( a - b );
}
template <typename T, int M> inline Bimpl<T,M> arh( const Bimpl<T,M>& a, const double k )
{ return exp( ( -k )/a ); }

} // namespace bdiff_detail

// ============================ public type aliases ===========================
template <typename T> using B  = bdiff_detail::Bimpl<T,0>;   //!< dense reverse-mode type
template <typename T> using SB = bdiff_detail::Bimpl<T,1>;   //!< sparse reverse-mode type

// =================== mc::Op specialization (covers B and SB) =================
template <typename U, int M> struct Op< bdiff_detail::Bimpl<U,M> >
{
  typedef bdiff_detail::Bimpl<U,M> TU;
  static TU point( const double c ) { throw std::runtime_error("mc::Op<mc::B<U>>::point -- operation not permitted"); }
  static TU zeroone() { throw std::runtime_error("mc::Op<mc::B<U>>::zeroone -- operation not permitted"); }
  static void I(TU& x, const TU& y) { x = y; }
  static double l(const TU& x) { throw std::runtime_error("mc::Op<mc::B<U>>::l -- operation not permitted"); }
  static double u(const TU& x) { throw std::runtime_error("mc::Op<mc::B<U>>::u -- operation not permitted"); }
  static double abs (const TU& x) { throw std::runtime_error("mc::Op<mc::B<U>>::abs -- operation not permitted"); }
  static double mid (const TU& x) { throw std::runtime_error("mc::Op<mc::B<U>>::mid -- operation not permitted"); }
  static double diam(const TU& x) { throw std::runtime_error("mc::Op<mc::B<U>>::diam -- operation not permitted"); }
  static TU inv (const TU& x) { return bdiff_detail::inv(x);  }
  static TU sqr (const TU& x) { return bdiff_detail::sqr(x);  }
  static TU sqrt(const TU& x) { return bdiff_detail::sqrt(x); }
  static TU exp (const TU& x) { return bdiff_detail::exp(x);  }
  static TU log (const TU& x) { return bdiff_detail::log(x);  }
  static TU xlog(const TU& x) { return bdiff_detail::xlog(x); }
  static TU lmtd(const TU& x, const TU& y) { return bdiff_detail::lmtd(x,y); }
  static TU rlmtd(const TU& x, const TU& y) { return bdiff_detail::rlmtd(x,y); }
  static TU fabs(const TU& x) { return bdiff_detail::fabs(x); }
  static TU sin (const TU& x) { return bdiff_detail::sin(x);  }
  static TU cos (const TU& x) { return bdiff_detail::cos(x);  }
  static TU tan (const TU& x) { return bdiff_detail::tan(x);  }
  static TU asin(const TU& x) { return bdiff_detail::asin(x); }
  static TU acos(const TU& x) { return bdiff_detail::acos(x); }
  static TU atan(const TU& x) { return bdiff_detail::atan(x); }
  static TU sinh(const TU& x) { return bdiff_detail::sinh(x); }
  static TU cosh(const TU& x) { return bdiff_detail::cosh(x); }
  static TU tanh(const TU& x) { return bdiff_detail::tanh(x); }
  static TU erf (const TU& x) { return bdiff_detail::erf(x);  }
  static TU erfc(const TU& x) { return bdiff_detail::erfc(x); }
  static TU fstep(const TU& x) { throw std::runtime_error("mc::Op<mc::B<U>>::fstep -- operation not permitted"); }
  static TU bstep(const TU& x) { throw std::runtime_error("mc::Op<mc::B<U>>::bstep -- operation not permitted"); }
  static TU hull(const TU& x, const TU& y) { throw std::runtime_error("mc::Op<mc::B<U>>::hull -- operation not permitted"); }
  template <typename Y> static TU min(const TU& x, const Y& y) { return bdiff_detail::min(x,y); }
  template <typename Y> static TU max(const TU& x, const Y& y) { return bdiff_detail::max(x,y); }
  static TU arh (const TU& x, const double k) { return bdiff_detail::arh(x,k); }
  template <typename X, typename Y> static TU pow(const X& x, const Y& y) { return bdiff_detail::pow(x,y); }
  static TU cheb(const TU& x, const unsigned n) { return bdiff_detail::cheb(x,n); }
  static TU prod (const unsigned n, const TU* x) { switch( n ){ case 0: return TU(1.); case 1: return x[0]; default: return x[0]*prod(n-1,x+1); } }
  static TU monom (const unsigned n, const TU* x, const unsigned* k) { switch( n ){ case 0: return TU(1.); case 1: return pow(x[0],(int)k[0]); default: return pow(x[0],(int)k[0])*monom(n-1,x+1,k+1); } }
  static bool inter(TU& xIy, const TU& x, const TU& y) { xIy = x; return true; }
  static bool eq(const TU& x, const TU& y) { throw std::runtime_error("mc::Op<mc::B<U>>::eq -- operation not permitted"); }
  static bool ne(const TU& x, const TU& y) { throw std::runtime_error("mc::Op<mc::B<U>>::ne -- operation not permitted"); }
  static bool lt(const TU& x, const TU& y) { throw std::runtime_error("mc::Op<mc::B<U>>::lt -- operation not permitted"); }
  static bool le(const TU& x, const TU& y) { throw std::runtime_error("mc::Op<mc::B<U>>::le -- operation not permitted"); }
  static bool gt(const TU& x, const TU& y) { throw std::runtime_error("mc::Op<mc::B<U>>::gt -- operation not permitted"); }
  static bool ge(const TU& x, const TU& y) { throw std::runtime_error("mc::Op<mc::B<U>>::ge -- operation not permitted"); }
};

} // namespace mc

#endif // MC__BDIFF_HPP
