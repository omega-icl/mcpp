// tdiff.hpp -- Taylor-mode automatic differentiation for MC++.
//
// Provides mc::T<U> : forward Taylor-coefficient type templated on the
// evaluation type U, exposing the same public interface as fadbad::T
// (operator[] to read/write the k-th Taylor coefficient, eval(q) to compute
// coefficients through order q, length(), reset(), assignment from a value)
// and agreeing with mc::Op<U>.
//
// A T<U> carries a truncated Taylor series in a scalar parameter t,
//     x(t) = sum_k x[k] t^k,   x[k] = (1/k!) d^k x/dt^k (0),
// with the independent variable seeded as x[0]=value, x[1]=direction.  Each
// operation records a node in a reference-counted (std::shared_ptr) DAG;
// eval(q) lazily computes and *caches* coefficients through order q using the
// standard Taylor recurrences (Cauchy products for * and /, and ODE-based
// recurrences for the elementary functions).  Because coefficients are cached
// and only ever extended, an outer driver may feed higher-order coefficients
// back into variable nodes between eval() calls (as in FFGraph::TAD) without
// invalidating the coefficients already computed.
//
// This is an independent implementation of the public interface only; it does
// not reproduce the internal machinery of FADBAD++.  The set of supported
// functions matches mc::Op<fadbad::T<U>> (inv, sqr, sqrt, exp, log, xlog, sin,
// cos, tan, asin, acos, atan, sinh, cosh, tanh, arh, pow); fabs, erf, min/max,
// cheb, etc. are not defined for Taylor mode.

#ifndef MC__TDIFF_HPP
#define MC__TDIFF_HPP

#include <memory>
#include <vector>
#include <stdexcept>
#include "mcop.hpp"

namespace mc
{

namespace tdiff_detail
{

template <typename U>
class Timpl
{
  enum Kind { LEAF, NEG, ADD, SUB, MUL, DIV,
              SADD, SSUB, RSUB, SMUL, SDIV, RDIV,
              SQR, SQRT, EXP,
              LOG, ASIN, ACOS, ATAN,        // integrate:  z' = r * x',  r = _ops[1]
              SIN, COS, SINH, COSH, TAN, TANH };

  struct Node
  {
    Kind                                _k;
    std::vector< std::shared_ptr<Node> > _ops;
    double                              _s = 0.;   // scalar parameter (scalar ops)
    std::vector<U>                      _c;         // cached Taylor coefficients
    std::vector<U>                      _c2;        // companion series (sin/cos/... )
    std::vector<Node*>                  _topo;      // cached eval order (root nodes only)
    unsigned                            _mark = 0;  // DFS epoch stamp (topo build)
    explicit Node( Kind k ) : _k(k) {}
  };

  std::shared_ptr<Node> _n;
  explicit Timpl( const std::shared_ptr<Node>& nd ) : _n(nd) {}

  static unsigned& _epoch() { static thread_local unsigned e = 0; return e; }

  static std::shared_ptr<Node> _leaf( const U& v )
  { auto nd=std::make_shared<Node>(LEAF); nd->_c.assign(1,v); return nd; }

  // ---- compute coefficient m of an operation node (operands ready to m) ----
  static void _compute( Node* p, unsigned m )
  {
    const auto& A = p->_ops[0]->_c;
    switch( p->_k ){
      case NEG:  p->_c.push_back( -A[m] ); return;
      case ADD:  p->_c.push_back( A[m] + p->_ops[1]->_c[m] ); return;
      case SUB:  p->_c.push_back( A[m] - p->_ops[1]->_c[m] ); return;
      case MUL:{ const auto& B=p->_ops[1]->_c; U s=A[0]*B[m]; for(unsigned j=1;j<=m;++j) s=s+A[j]*B[m-j]; p->_c.push_back(s); return; }
      case DIV:{ const auto& B=p->_ops[1]->_c; U s=A[m]; for(unsigned j=0;j<m;++j) s=s-p->_c[j]*B[m-j]; p->_c.push_back( s/B[0] ); return; }
      case SADD: p->_c.push_back( m==0? A[0]+p->_s : A[m] ); return;
      case SSUB: p->_c.push_back( m==0? A[0]-p->_s : A[m] ); return;
      case RSUB: p->_c.push_back( m==0? p->_s-A[0] : -A[m] ); return;
      case SMUL: p->_c.push_back( p->_s * A[m] ); return;
      case SDIV: p->_c.push_back( A[m] / p->_s ); return;
      case RDIV:{ U s = (m==0? U(p->_s) : U(0.)); for(unsigned j=0;j<m;++j) s=s-p->_c[j]*A[m-j]; p->_c.push_back( s/A[0] ); return; }
      case SQR:{ U s=A[0]*A[m]; for(unsigned j=1;j<=m;++j) s=s+A[j]*A[m-j]; p->_c.push_back(s); return; }
      case SQRT:{
        if(m==0){ p->_c.push_back( mc::Op<U>::sqrt(A[0]) ); return; }
        U s=A[m]; for(unsigned j=1;j<m;++j) s=s-p->_c[j]*p->_c[m-j];
        p->_c.push_back( s/(2.*p->_c[0]) ); return; }
      case EXP:{
        if(m==0){ p->_c.push_back( mc::Op<U>::exp(A[0]) ); return; }
        U s=(double)1*A[1]*p->_c[m-1]; for(unsigned j=2;j<=m;++j) s=s+(double)j*A[j]*p->_c[m-j];
        p->_c.push_back( s/(double)m ); return; }
      case LOG:  case ASIN: case ACOS: case ATAN:{        // z0 = f(x0); z_m = (1/m) sum_{j=1}^m j x_j r_{m-j}
        if(m==0){
          U z0 = p->_k==LOG?  mc::Op<U>::log (A[0])
               : p->_k==ASIN? mc::Op<U>::asin(A[0])
               : p->_k==ACOS? mc::Op<U>::acos(A[0])
               :              mc::Op<U>::atan(A[0]);
          p->_c.push_back( z0 ); return; }
        const auto& R=p->_ops[1]->_c;
        U s=(double)1*A[1]*R[m-1]; for(unsigned j=2;j<=m;++j) s=s+(double)j*A[j]*R[m-j];
        p->_c.push_back( s/(double)m ); return; }
      case SIN: case COS:{     // primary + companion (cos/sin resp.) in _c2
        if(m==0){
          if(p->_k==SIN){ p->_c.push_back(mc::Op<U>::sin(A[0])); p->_c2.push_back(mc::Op<U>::cos(A[0])); }
          else          { p->_c.push_back(mc::Op<U>::cos(A[0])); p->_c2.push_back(mc::Op<U>::sin(A[0])); }
          return; }
        // sin_m = (1/m) sum j x_j cos_{m-j} ; cos_m = -(1/m) sum j x_j sin_{m-j}
        U ssin=(double)0, scos=(double)0;
        const std::vector<U>& S = (p->_k==SIN)? p->_c  : p->_c2;   // sin series
        const std::vector<U>& C = (p->_k==SIN)? p->_c2 : p->_c;    // cos series
        for(unsigned j=1;j<=m;++j){ ssin=ssin+(double)j*A[j]*C[m-j]; scos=scos+(double)j*A[j]*S[m-j]; }
        U sin_m =  ssin/(double)m, cos_m = -scos/(double)m;
        if(p->_k==SIN){ p->_c.push_back(sin_m); p->_c2.push_back(cos_m); }
        else          { p->_c.push_back(cos_m); p->_c2.push_back(sin_m); }
        return; }
      case SINH: case COSH:{
        if(m==0){
          if(p->_k==SINH){ p->_c.push_back(mc::Op<U>::sinh(A[0])); p->_c2.push_back(mc::Op<U>::cosh(A[0])); }
          else           { p->_c.push_back(mc::Op<U>::cosh(A[0])); p->_c2.push_back(mc::Op<U>::sinh(A[0])); }
          return; }
        U ssinh=(double)0, scosh=(double)0;
        const std::vector<U>& SH=(p->_k==SINH)? p->_c : p->_c2;
        const std::vector<U>& CH=(p->_k==SINH)? p->_c2: p->_c;
        for(unsigned j=1;j<=m;++j){ ssinh=ssinh+(double)j*A[j]*CH[m-j]; scosh=scosh+(double)j*A[j]*SH[m-j]; }
        U sinh_m=ssinh/(double)m, cosh_m=scosh/(double)m;
        if(p->_k==SINH){ p->_c.push_back(sinh_m); p->_c2.push_back(cosh_m); }
        else           { p->_c.push_back(cosh_m); p->_c2.push_back(sinh_m); }
        return; }
      case TAN: case TANH:{    // primary z ; companion w = 1 +/- z^2  in _c2
        if(m==0){
          U z0 = (p->_k==TAN)? mc::Op<U>::tan(A[0]) : mc::Op<U>::tanh(A[0]);
          p->_c.push_back( z0 );
          p->_c2.push_back( (p->_k==TAN)? 1.+z0*z0 : 1.-z0*z0 );
          return; }
        // z_m = (1/m) sum_{j=1}^m j x_j w_{m-j}
        U s=(double)1*A[1]*p->_c2[m-1]; for(unsigned j=2;j<=m;++j) s=s+(double)j*A[j]*p->_c2[m-j];
        U z_m = s/(double)m; p->_c.push_back( z_m );
        // w_m = +/-( 2 z_0 z_m + sum_{j=1}^{m-1} z_j z_{m-j} )   (w = 1 +/- z^2)
        U zz = 2.*p->_c[0]*z_m; for(unsigned j=1;j<m;++j) zz=zz+p->_c[j]*p->_c[m-j];
        p->_c2.push_back( (p->_k==TAN)? zz : -zz );
        return; }
      default: throw std::runtime_error("mc::T: unhandled operation");
    }
  }

  // build a post-order (sources-first) evaluation order for the DAG rooted at
  // root, cached in root->_topo so repeated eval() calls reuse it
  static void _build_topo( Node* root )
  {
    root->_topo.clear();
    static thread_local std::vector< std::pair<Node*,std::size_t> > stack;
    stack.clear();
    const unsigned e = ++_epoch();
    stack.emplace_back( root, 0 ); root->_mark = e;
    while( !stack.empty() ){
      auto& top = stack.back();
      if( top.second < top.first->_ops.size() ){
        Node* ch = top.first->_ops[top.second].get(); ++top.second;
        if( ch->_mark != e ){ ch->_mark = e; stack.emplace_back( ch, 0 ); }
      }
      else{ root->_topo.push_back( top.first ); stack.pop_back(); }
    }
  }

  // extend every node's cached coefficients through order k, in topological
  // order (operands before dependents), so each recurrence sees ready operands
  static void _grow( const std::shared_ptr<Node>& root, unsigned k )
  {
    if( root->_k==LEAF ){ if( root->_c.size()<=k ) root->_c.resize( k+1, U(0.) ); return; }
    if( root->_topo.empty() ) _build_topo( root.get() );
    for( Node* p : root->_topo ){
      if( p->_k==LEAF ){ if( p->_c.size()<=k ) p->_c.resize( k+1, U(0.) ); continue; }
      p->_c.reserve( k+1 );
      for( unsigned m=p->_c.size(); m<=k; ++m ) _compute( p, m );
    }
  }

  static Timpl _un ( Kind k, const Timpl& a )
  { auto nd=std::make_shared<Node>(k); nd->_ops.reserve(1); nd->_ops.push_back(a._n); return Timpl(nd); }
  static Timpl _un ( Kind k, const Timpl& a, double s )
  { auto nd=std::make_shared<Node>(k); nd->_ops.reserve(1); nd->_ops.push_back(a._n); nd->_s=s; return Timpl(nd); }
  static Timpl _bin( Kind k, const Timpl& a, const Timpl& b )
  { auto nd=std::make_shared<Node>(k); nd->_ops.reserve(2); nd->_ops.push_back(a._n); nd->_ops.push_back(b._n); return Timpl(nd); }
  static Timpl _integ( Kind k, const Timpl& a, const Timpl& r )  // LOG/ASIN/ACOS/ATAN
  { auto nd=std::make_shared<Node>(k); nd->_ops.reserve(2); nd->_ops.push_back(a._n); nd->_ops.push_back(r._n); return Timpl(nd); }

  template <typename V,int> friend struct OpHelper; // (unused hook)

public:
  // -------------------- construction / assignment --------------------------
  Timpl() : _n( std::make_shared<Node>(LEAF) ) {}                       // empty (length 0)
  Timpl( const U& v ) : _n( _leaf(v) ) {}
  template <typename V> Timpl( const V& v ) : _n( _leaf( U(v) ) ) {}
  Timpl( const Timpl& ) = default;                                     // shares node
  Timpl& operator=( const Timpl& ) = default;                          // rebinds
  Timpl& operator=( const U& v ) { _n = _leaf(v); return *this; }
  template <typename V> Timpl& operator=( const V& v ) { _n = _leaf(U(v)); return *this; }

  // -------------------- coefficient access / evaluation --------------------
  //! @brief read/write the k-th Taylor coefficient (resizes with zeros on a leaf/growing access)
  U&       operator[]( const unsigned k )       { if(_n->_c.size()<=k) _n->_c.resize(k+1,U(0.)); return _n->_c[k]; }
  const U& operator[]( const unsigned k ) const { return _n->_c[k]; }
  //! @brief compute Taylor coefficients through order q; returns q
  unsigned eval( const unsigned q ) { _grow(_n,q); return q; }
  //! @brief number of coefficients currently stored
  unsigned length() const { return (unsigned)_n->_c.size(); }
  //! @brief value (0th Taylor coefficient)
  const U& val() const { return _n->_c[0]; }
  U&       val()       { if(_n->_c.empty()) _n->_c.resize(1,U(0.)); return _n->_c[0]; }
  //! @brief drop cached coefficients of operation nodes so a re-eval recomputes them
  void reset()
  {
    std::vector<Node*> stk{ _n.get() };
    while(!stk.empty()){ Node* p=stk.back(); stk.pop_back();
      if(p->_k!=LEAF){ p->_c.clear(); p->_c2.clear(); }
      for(auto& o:p->_ops) stk.push_back(o.get());
    }
  }

  // -------------------- node builders exposed to free functions ------------
  static Timpl _unary ( Kind k, const Timpl& a )              { return _un(k,a); }
  static Timpl _sunary( Kind k, const Timpl& a, double s )    { return _un(k,a,s); }
  static Timpl _binary( Kind k, const Timpl& a, const Timpl& b){ return _bin(k,a,b); }
  static Timpl _integrate( Kind k, const Timpl& a, const Timpl& r ){ return _integ(k,a,r); }
  // Kind tags visible to the (friend) free functions:
  static constexpr Kind K_NEG=NEG,  K_ADD=ADD,  K_SUB=SUB,  K_MUL=MUL,  K_DIV=DIV,
    K_SADD=SADD, K_SSUB=SSUB, K_RSUB=RSUB, K_SMUL=SMUL, K_SDIV=SDIV, K_RDIV=RDIV,
    K_SQR=SQR, K_SQRT=SQRT, K_EXP=EXP, K_LOG=LOG, K_ASIN=ASIN, K_ACOS=ACOS, K_ATAN=ATAN,
    K_SIN=SIN, K_COS=COS, K_SINH=SINH, K_COSH=COSH, K_TAN=TAN, K_TANH=TANH;

  // compound assignment
  Timpl& operator+=( const Timpl& b ){ *this=*this+b; return *this; }
  Timpl& operator-=( const Timpl& b ){ *this=*this-b; return *this; }
  Timpl& operator*=( const Timpl& b ){ *this=*this*b; return *this; }
  Timpl& operator/=( const Timpl& b ){ *this=*this/b; return *this; }
  template <typename V> Timpl& operator+=( const V& s ){ *this=*this+s; return *this; }
  template <typename V> Timpl& operator-=( const V& s ){ *this=*this-s; return *this; }
  template <typename V> Timpl& operator*=( const V& s ){ *this=*this*s; return *this; }
  template <typename V> Timpl& operator/=( const V& s ){ *this=*this/s; return *this; }
};

// ============================ arithmetic operators ==========================
template <typename U> inline Timpl<U> operator+( const Timpl<U>& a ){ return a; }
template <typename U> inline Timpl<U> operator-( const Timpl<U>& a ){ return Timpl<U>::_unary(Timpl<U>::K_NEG,a); }
template <typename U> inline Timpl<U> operator+( const Timpl<U>& a, const Timpl<U>& b ){ return Timpl<U>::_binary(Timpl<U>::K_ADD,a,b); }
template <typename U> inline Timpl<U> operator-( const Timpl<U>& a, const Timpl<U>& b ){ return Timpl<U>::_binary(Timpl<U>::K_SUB,a,b); }
template <typename U> inline Timpl<U> operator*( const Timpl<U>& a, const Timpl<U>& b ){ return Timpl<U>::_binary(Timpl<U>::K_MUL,a,b); }
template <typename U> inline Timpl<U> operator/( const Timpl<U>& a, const Timpl<U>& b ){ return Timpl<U>::_binary(Timpl<U>::K_DIV,a,b); }

template <typename U,typename V> inline Timpl<U> operator+( const Timpl<U>& a, const V& s ){ return Timpl<U>::_sunary(Timpl<U>::K_SADD,a,(double)s); }
template <typename U,typename V> inline Timpl<U> operator+( const V& s, const Timpl<U>& a ){ return Timpl<U>::_sunary(Timpl<U>::K_SADD,a,(double)s); }
template <typename U,typename V> inline Timpl<U> operator-( const Timpl<U>& a, const V& s ){ return Timpl<U>::_sunary(Timpl<U>::K_SSUB,a,(double)s); }
template <typename U,typename V> inline Timpl<U> operator-( const V& s, const Timpl<U>& a ){ return Timpl<U>::_sunary(Timpl<U>::K_RSUB,a,(double)s); }
template <typename U,typename V> inline Timpl<U> operator*( const Timpl<U>& a, const V& s ){ return Timpl<U>::_sunary(Timpl<U>::K_SMUL,a,(double)s); }
template <typename U,typename V> inline Timpl<U> operator*( const V& s, const Timpl<U>& a ){ return Timpl<U>::_sunary(Timpl<U>::K_SMUL,a,(double)s); }
template <typename U,typename V> inline Timpl<U> operator/( const Timpl<U>& a, const V& s ){ return Timpl<U>::_sunary(Timpl<U>::K_SDIV,a,(double)s); }
template <typename U,typename V> inline Timpl<U> operator/( const V& s, const Timpl<U>& a ){ return Timpl<U>::_sunary(Timpl<U>::K_RDIV,a,(double)s); }

// ============================ elementary functions ==========================
template <typename U> inline Timpl<U> sqr ( const Timpl<U>& a ){ return Timpl<U>::_unary(Timpl<U>::K_SQR ,a); }
template <typename U> inline Timpl<U> sqrt( const Timpl<U>& a ){ return Timpl<U>::_unary(Timpl<U>::K_SQRT,a); }
template <typename U> inline Timpl<U> exp ( const Timpl<U>& a ){ return Timpl<U>::_unary(Timpl<U>::K_EXP ,a); }
template <typename U> inline Timpl<U> inv ( const Timpl<U>& a ){ return Timpl<U>::_sunary(Timpl<U>::K_RDIV,a,1.); }
template <typename U> inline Timpl<U> sin ( const Timpl<U>& a ){ return Timpl<U>::_unary(Timpl<U>::K_SIN ,a); }
template <typename U> inline Timpl<U> cos ( const Timpl<U>& a ){ return Timpl<U>::_unary(Timpl<U>::K_COS ,a); }
template <typename U> inline Timpl<U> sinh( const Timpl<U>& a ){ return Timpl<U>::_unary(Timpl<U>::K_SINH,a); }
template <typename U> inline Timpl<U> cosh( const Timpl<U>& a ){ return Timpl<U>::_unary(Timpl<U>::K_COSH,a); }
template <typename U> inline Timpl<U> tan ( const Timpl<U>& a ){ return Timpl<U>::_unary(Timpl<U>::K_TAN ,a); }
template <typename U> inline Timpl<U> tanh( const Timpl<U>& a ){ return Timpl<U>::_unary(Timpl<U>::K_TANH,a); }
// integrate-form: z' = r x'
template <typename U> inline Timpl<U> log ( const Timpl<U>& a ){ return Timpl<U>::_integrate(Timpl<U>::K_LOG , a, inv(a) ); }
template <typename U> inline Timpl<U> asin( const Timpl<U>& a ){ return Timpl<U>::_integrate(Timpl<U>::K_ASIN, a, inv(sqrt(1.-sqr(a))) ); }
template <typename U> inline Timpl<U> acos( const Timpl<U>& a ){ return Timpl<U>::_integrate(Timpl<U>::K_ACOS, a, -inv(sqrt(1.-sqr(a))) ); }
template <typename U> inline Timpl<U> atan( const Timpl<U>& a ){ return Timpl<U>::_integrate(Timpl<U>::K_ATAN, a, inv(1.+sqr(a)) ); }
// compositions
template <typename U> inline Timpl<U> xlog( const Timpl<U>& a ){ return a*log(a); }
template <typename U> inline Timpl<U> arh ( const Timpl<U>& a, const double k ){ return exp( (-k)*inv(a) ); }

// powers
template <typename U> inline Timpl<U> pow( const Timpl<U>& a, const int n )
{
  if( n==0 ) return Timpl<U>( U(1.) );
  if( n<0  ) return inv( pow(a,-n) );
  Timpl<U> r( U(1.) ), b(a); int e=n;
  while( e ){ if(e&1) r=r*b; e>>=1; if(e) b=b*b; }
  return r;
}
template <typename U> inline Timpl<U> pow( const Timpl<U>& a, const double r ){ return exp( r*log(a) ); }
template <typename U> inline Timpl<U> pow( const Timpl<U>& a, const Timpl<U>& b ){ return exp( b*log(a) ); }
template <typename U,typename V> inline Timpl<U> pow( const V& s, const Timpl<U>& a ){ return exp( a * Timpl<U>( mc::Op<U>::log(U(s)) ) ); }

} // namespace tdiff_detail

// ============================ public type alias =============================
template <typename U> using T = tdiff_detail::Timpl<U>;   //!< Taylor-mode type

// =================== mc::Op specialization for mc::T<U> =====================
template <typename U> struct Op< tdiff_detail::Timpl<U> >
{
  typedef tdiff_detail::Timpl<U> TU;
  static TU point( const double c ) { throw std::runtime_error("mc::Op<mc::T<U>>::point -- operation not permitted"); }
  static TU zeroone() { throw std::runtime_error("mc::Op<mc::T<U>>::zeroone -- operation not permitted"); }
  static void I(TU& x, const TU& y) { x = y; }
  static double l(const TU& x) { throw std::runtime_error("mc::Op<mc::T<U>>::l -- operation not permitted"); }
  static double u(const TU& x) { throw std::runtime_error("mc::Op<mc::T<U>>::u -- operation not permitted"); }
  static double abs (const TU& x) { throw std::runtime_error("mc::Op<mc::T<U>>::abs -- operation not permitted"); }
  static double mid (const TU& x) { throw std::runtime_error("mc::Op<mc::T<U>>::mid -- operation not permitted"); }
  static double diam(const TU& x) { throw std::runtime_error("mc::Op<mc::T<U>>::diam -- operation not permitted"); }
  static TU inv (const TU& x) { return tdiff_detail::inv(x);  }
  static TU sqr (const TU& x) { return tdiff_detail::sqr(x);  }
  static TU sqrt(const TU& x) { return tdiff_detail::sqrt(x); }
  static TU exp (const TU& x) { return tdiff_detail::exp(x);  }
  static TU log (const TU& x) { return tdiff_detail::log(x);  }
  static TU xlog(const TU& x) { return tdiff_detail::xlog(x); }
  static TU lmtd(const TU& x, const TU& y) { throw std::runtime_error("mc::Op<mc::T<U>>::lmtd -- operation not permitted"); }
  static TU rlmtd(const TU& x, const TU& y) { throw std::runtime_error("mc::Op<mc::T<U>>::rlmtd -- operation not permitted"); }
  static TU fabs(const TU& x) { throw std::runtime_error("mc::Op<mc::T<U>>::fabs -- operation not permitted"); }
  static TU sin (const TU& x) { return tdiff_detail::sin(x);  }
  static TU cos (const TU& x) { return tdiff_detail::cos(x);  }
  static TU tan (const TU& x) { return tdiff_detail::tan(x);  }
  static TU asin(const TU& x) { return tdiff_detail::asin(x); }
  static TU acos(const TU& x) { return tdiff_detail::acos(x); }
  static TU atan(const TU& x) { return tdiff_detail::atan(x); }
  static TU sinh(const TU& x) { return tdiff_detail::sinh(x); }
  static TU cosh(const TU& x) { return tdiff_detail::cosh(x); }
  static TU tanh(const TU& x) { return tdiff_detail::tanh(x); }
  static TU erf (const TU& x) { throw std::runtime_error("mc::Op<mc::T<U>>::erf -- operation not permitted"); }
  static TU erfc(const TU& x) { throw std::runtime_error("mc::Op<mc::T<U>>::erfc -- operation not permitted"); }
  static TU fstep(const TU& x) { throw std::runtime_error("mc::Op<mc::T<U>>::fstep -- operation not permitted"); }
  static TU bstep(const TU& x) { throw std::runtime_error("mc::Op<mc::T<U>>::bstep -- operation not permitted"); }
  static TU hull(const TU& x, const TU& y) { throw std::runtime_error("mc::Op<mc::T<U>>::hull -- operation not permitted"); }
  static TU min (const TU& x, const TU& y) { throw std::runtime_error("mc::Op<mc::T<U>>::min -- operation not permitted"); }
  static TU max (const TU& x, const TU& y) { throw std::runtime_error("mc::Op<mc::T<U>>::max -- operation not permitted"); }
  static TU arh (const TU& x, const double k) { return tdiff_detail::arh(x,k); }
  template <typename X, typename Y> static TU pow(const X& x, const Y& y) { return tdiff_detail::pow(x,y); }
  static TU cheb(const TU& x, const unsigned n) { throw std::runtime_error("mc::Op<mc::T<U>>::cheb -- operation not permitted"); }
  static TU prod (const unsigned n, const TU* x) { switch( n ){ case 0: return TU(U(1.)); case 1: return x[0]; default: return x[0]*prod(n-1,x+1); } }
  static TU monom (const unsigned n, const TU* x, const unsigned* k) { switch( n ){ case 0: return TU(U(1.)); case 1: return pow(x[0],(int)k[0]); default: return pow(x[0],(int)k[0])*monom(n-1,x+1,k+1); } }
  static bool inter(TU& xIy, const TU& x, const TU& y) { xIy = x; return true; }
  static bool eq(const TU& x, const TU& y) { throw std::runtime_error("mc::Op<mc::T<U>>::eq -- operation not permitted"); }
  static bool ne(const TU& x, const TU& y) { throw std::runtime_error("mc::Op<mc::T<U>>::ne -- operation not permitted"); }
  static bool lt(const TU& x, const TU& y) { throw std::runtime_error("mc::Op<mc::T<U>>::lt -- operation not permitted"); }
  static bool le(const TU& x, const TU& y) { throw std::runtime_error("mc::Op<mc::T<U>>::le -- operation not permitted"); }
  static bool gt(const TU& x, const TU& y) { throw std::runtime_error("mc::Op<mc::T<U>>::gt -- operation not permitted"); }
  static bool ge(const TU& x, const TU& y) { throw std::runtime_error("mc::Op<mc::T<U>>::ge -- operation not permitted"); }
};

} // namespace mc

#endif // MC__TDIFF_HPP
