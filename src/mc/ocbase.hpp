// Copyright (C) Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the Eclipse Public License.

#ifndef MC__OCBASE_HPP
#define MC__OCBASE_HPP

#include <type_traits>
#include <cassert>
#include <cmath>
#include <cfloat>
#include <limits>
#include <numeric>
#include <algorithm>
#include <optional>
#include <functional>
#include <utility>
#include <vector>
#include <set>
#include <map>
#include <tuple>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>

#include <armadillo>

#include "mcfunc.hpp" // for mc::PI
#include "ffexpr.hpp"
#include "ffdep.hpp"
#include "slift.hpp"
#include "smon.hpp"

namespace mc
{

//! @brief C++ class for computing quadrature nodes for various schemes and interpolating polynomials
////////////////////////////////////////////////////////////////////////
//! mc::BASE_OC is a C++ class for computing nodes for Gauss-Legendre,
//! Gauss-Legendre-Radau, Gauss-Legendre-Lobatto, and
//! Chebyshev-Gauss-Lobatto quadrature and evaluating Lagrange
//! interpolating polynomials between such nodes
////////////////////////////////////////////////////////////////////////
class BASE_OC
////////////////////////////////////////////////////////////////////////
{
public:

  //! @brief Default constructor
  BASE_OC
    ()
    {}

  //! @brief Destructor
  virtual ~BASE_OC
    ()
    {}

  //! @brief retrieve Gauss-Legendre quadrature nodes
  std::vector<double> const& lgnodes
    ()
    const
    { return _x; }

  //! @brief retrieve Gauss-Legendre quadrature nodes rescaled in [a,b]
  std::vector<double> lgnodes
    ( double const& a, double const& b )
    const
    {
      auto x = _x;
      double r = (b-a)/2., m=(a+b)/2.;
      for( auto& xi : x ) (xi *= r) += m;
      return x;
    }

  //! @brief retrieve Gauss-Legendre quadrature weights
  std::vector<double> const& lgweights
    ()
    const
    { return _w; }

  //! @brief generate <a>N</a> quadrature nodes in range [a,b] for: Legendre-Gauss (TYPE=0), Legendre-Gauss-Radau (TYPE=1 or -1), Legendre-Gauss-Lobatto (TYPE=2), or Chebyshev-Gauss-Lobatto (TYPE=3)
  bool set_lgnodes
    ( int const TYPE, size_t const N, double const& a=-1., double const& b=1.,
      double const& TOL=1E-7 );

  //! @brief set Lagrange evaluator for the current nodes by precomputing weights, so we can reuse them for many different y vectors to reduce the complexity
  bool set_lagrange
    ( double const& TOL=DBL_EPSILON );

  //! @brief Precompute and cache nodal differentiation matrices for Lagrange basis up to order MAXORD (enables fast evaluation when t matches a node)
  bool set_lagrange_ops
    ( int const MAXORD=2 );

  //! @brief Templated Lagrange evaluator of order ORD >= 0 at point t for given data yin 
  template <typename U>
  bool eval_lagrange
    ( U& yout, double const& t, int const ORD, std::vector<U> const& yin,
      double const& TOL=DBL_EPSILON, double const& RESCALE=1. )
    const;

  //! @brief Templated Lagrange evaluator of order ORD >= 0 at point t for strided/blocked data yin; STRIDE is the product of the dimensions before the active one; OUTER is the product of the dimensions after the active one; and the active dimension length is _x.size()
  template <typename U>
  bool eval_lagrange
    ( std::vector<U>& yout, double const& t, int const ORD, std::vector<U> const& yin,
      size_t const STRIDE, size_t const OUTER, double const& TOL=DBL_EPSILON,
      double const& RESCALE=1. )
    const;

  //! @brief Templated Lagrange evaluator of order ORD >= 0 at all nodes for given data yin
  template <typename U>
  bool eval_lagrange
    ( std::vector<U>& yout, int const ORD, std::vector<U> const& yin, double const& RESCALE=1. )
    const;

  //! @brief Templated Lagrange evaluator of order ORD >= 0 for strided/blocked data yin; STRIDE is the product of the dimensions before the active one; OUTER is the product of the dimensions after the active one; and the active dimension length is _x.size()
  template <typename U>
  bool eval_lagrange
    ( std::vector<U>& yout, int const ORD, std::vector<U> const& yin,
      size_t const STRIDE, size_t const OUTER, double const& RESCALE=1. )
    const;

  //! @brief Templated weight/gradient of Lagrange interpolant (or its ORD-th derivative) w.r.t. coefficients y:  p^(ORD)(t) = dot(w,y)
  template <typename U>
  bool w_eval_lagrange
    ( std::vector<U>& w, double const& t, int const ORD,
      double const& TOL=DBL_EPSILON )
    const;

  //! @brief Templated convenience: compute value and weight/gradient simultaneously
  template <typename U>
  bool eval_lagrange_grad
    ( U& val, std::vector<U>& w, double const& t, int const ORD,
      std::vector<U> const& y, double const& TOL=DBL_EPSILON )
    const;
    
  //! @brief Templated integral evaluator for the Lagrange interpolant over [a,b]
  template <typename U>
  bool quad_lagrange
    ( U& yout, std::vector<U> const& yin, double const& RESCALE=1. )
    const;

  //! @brief Templated integral evaluator for the Lagrange interpolant over [a,b]; STRIDE is the product of the dimensions before the active one; OUTER is the product of the dimensions after the active one; and the active dimension length is _x.size()
  template <typename U>
  bool quad_lagrange
    ( std::vector<U>& yout, std::vector<U> const& yin,
      size_t const STRIDE=1, size_t const OUTER=1, double const& RESCALE=1. )
    const;

  //! @brief Templated quadrature weight/gradient w.r.t. coefficients y:  ∫ p(t) dt = dot(wq,y)
  template <typename U>
  bool w_quad_lagrange
    ( std::vector<U>& wq )
    const;

  //! @brief Templated convenience: compute integral and quadrature weight/gradient simultaneously
  template <typename U>
  bool quad_lagrange_grad
    ( U& val, std::vector<U>& wq, std::vector<U> const& y )
    const;

protected:

  // storage for Legendre polynomial computation
  std::vector<double>             _P;

  // storage for quadrature nodes
  std::vector<double>             _x;

  // storage for quadrature weights
  std::vector<double>             _w;

  // storage for Lagrange weights: w_i = 1 / Π_{j≠i} (x_i - x_j)
  std::vector<double>             _Lw;

  // storage for Lagrange weights: c_i = Σ_{m≠i} 1/(x_i - x_m)
  std::vector<double>             _Lc;

  // Cached nodal differentiation matrices D^(m) (Armadillo), built on demand or via set_lagrange_ops
  mutable std::vector<arma::mat>  _D;
  mutable int                     _Dmax = -1;

  // storage for quadrature weights in Lagrange evaluation
  mutable std::vector<double>     _weval;

  // Ensure cached nodal differentiation matrices are available up to order MAXORD
  bool _ensure_lagrange_ops
    ( int const MAXORD )
    const;

  // generate <a>N</a> Legendre-Gauss collocation points between -1 and 1
  void _lgnodes
    ( size_t const N, double const& TOL );
  // generate <a>N</a> Legendre-Gauss-Radau collocation points between -1 and 1
  void _lgrnodes
    ( size_t const N, double const& TOL );
  // generate <a>N</a> Legendre-Gauss-Lobatto collocation points between -1 and 1
  void _lglnodes
    ( size_t const N, double const& TOL );

  // generate <a>N</a> Chebyshev-Gauss-Lobatto collocation points between -1 and 1
  void _cglnodes
    ( size_t const N );
    
  //! @brief Templated Lagrange evaluator at point t for strided/blocked data yin; STRIDE is the product of the dimensions before the active one; OUTER is the product of the dimensions after the active one; and the active dimension length is _x.size()
  template <typename U>
  bool _eval_lagrange
    ( U* yout, double const& t, U const* yin, size_t const STRIDE, size_t const OUTER,
      bool const sub, double const& TOL=DBL_EPSILON )
    const;

  //! @brief Like _eval_lagrange(yout,t,yin,...) but evaluates the ORD-th derivative
  //! of the Lagrange interpolant at point t, rescaled by RESCALE.
  //! Overwrites (sub=false) or subtracts from (sub=true) existing yout content.
  //! Use RESCALE = pow(2/w_elem, ORD) to convert from reference-interval derivative
  //! to physical-space derivative.
  template <typename U>
  bool _eval_lagrange
    ( U* yout, double const& t, int const ORD, U const* yin, size_t const STRIDE,
      size_t const OUTER, bool const sub, double const& TOL=DBL_EPSILON,
      double const& RESCALE=1. )
    const;
};

class OCBase;
class OCEnv;
class FFIntegral;
class FFPartial;

//! @brief Definition of domain and options for orthogonal collocation of distributed subexpressions
////////////////////////////////////////////////////////////////////////
//! mc::OCDom is a C++ class defining the domain and options for
//! orthogonal collocation of distributed subexpressions
////////////////////////////////////////////////////////////////////////
class OCDom
////////////////////////////////////////////////////////////////////////
: public virtual BASE_OC
{
  friend OCEnv;

public:

  //! @brief Exceptions
  class Exceptions
  {
  public:
    //! @brief Enumeration type for exception handling
    enum TYPE{
      BOUNDS = 1,	//!< Invalid domain bounds
      ELEMENTS,		//!< Invalid number of elements
      NODES,		//!< Invalid number of collocation points
      UNDEF = -33 	//!< Undefined
    };
    //! @brief Constructor for error <a>ierr</a>
    Exceptions( TYPE ierr=UNDEF ) : _ierr( ierr ){}
    //! @brief Inline function returning the error flag
    int ierr(){ return _ierr; }
    //! @brief Error description
    std::string what(){
      switch( _ierr ){
      case BOUNDS:
        return "OCDom::Exceptions  Invalid domain bounds";
      case ELEMENTS:
        return "OCDom::Exceptions  Invalid number of elements";
      case NODES:
        return "OCDom::Exceptions  Invalid number of collocation points";
      case UNDEF:
      default:
        return "OCDom::Exceptions  Undocumented exception";
      }
    }
  private:
    TYPE _ierr;
  };

  //! @brief Enumeration for equation domain: [a,b]:0, a:-1, b:-2, (a,b]=[a,b]-a:1, [a,b)=[a,b]-b:2, (a,b)=[a,b]-a-b:3
  enum DOM{
    ALL =  0, //!< [a,b]
    LB  = -1, //!< a
    UB  = -2  //!< b
  };

  //! @brief Enumeration for collocation type
  enum TYPE{
    LG  = 0, //!< Legendre-Gauss
    LGR = 1, //!< Legendre-Gauss-Radau, left endpoint included
    LGL = 2, //!< Legendre-Gauss-Lobatto, both endpoints included
    CGL = 3  //!< Chebyshev-Gauss-Lobatto, both endpoints included
  };

  //! @brief lower range of distributed domain
  double               lo_dom = 0.;

  //! @brief upper range of distributed domain
  double               up_dom = 0.;

  //! @brief type of collocation scheme
  TYPE                 type = LGR;

  //! @brief number of finite elements
  size_t               n_elem = 0;

  //! @brief representative finite-element width. For uniform grids this is the element width; for nonuniform grids this is the arithmetic mean and element-wise widths in elem_len are authoritative.
  double               w_elem = 0.;

  //! @brief number of collocation points on each element
  size_t               n_node = 0;

  //! @brief finite-element boundary coordinates, size n_elem+1
  std::vector<double>  elem_bnd;

  //! @brief finite-element lengths, size n_elem
  std::vector<double>  elem_len;

  //! @brief Default constructor
  OCDom
    ()
    {}

  //! @brief Uniform finite-element constructor
  OCDom
    ( double const& lo_dom, double const& up_dom,
      size_t n_elem, TYPE const& type, size_t n_node )
    : type   ( type   ),
      n_node ( n_node )
    {
      _set_uniform( lo_dom, up_dom, n_elem );
      if( !n_node ) throw Exceptions::NODES;
    }

  //! @brief Nonuniform finite-element constructor from element boundaries
  OCDom
    ( std::vector<double> const& elem_bnd, TYPE const& type, size_t n_node )
    : type   ( type   ),
      n_node ( n_node )
    {
      _set_boundaries( elem_bnd );
      if( !n_node ) throw Exceptions::NODES;
    }

  //! @brief Nonuniform finite-element constructor from lower bound and element lengths
  OCDom
    ( double const& lo_dom, std::vector<double> const& elem_len,
      TYPE const& type, size_t n_node )
    : type   ( type   ),
      n_node ( n_node )
    {
      _set_lengths( lo_dom, elem_len );
      if( !n_node ) throw Exceptions::NODES;
    }

  //! @brief Copy constructor
  OCDom
    ( OCDom const& other )
    : BASE_OC  ( other          ),
      lo_dom   ( other.lo_dom   ),
      up_dom   ( other.up_dom   ),
      type     ( other.type     ),
      n_elem   ( other.n_elem   ),
      w_elem   ( other.w_elem   ),
      n_node   ( other.n_node   ),
      elem_bnd ( other.elem_bnd ),
      elem_len ( other.elem_len )
    {
      _validate();
    }

  //! @brief Copy assignment
  OCDom& operator=
    ( OCDom const& other )
    {
      if( this == &other ) return *this;
      BASE_OC::operator=( other );
      lo_dom   = other.lo_dom;
      up_dom   = other.up_dom;
      type     = other.type;
      n_elem   = other.n_elem;
      w_elem   = other.w_elem;
      n_node   = other.n_node;
      elem_bnd = other.elem_bnd;
      elem_len = other.elem_len;
      _validate();
      return *this;
    }

  //! @brief True if all finite elements have the same length within tolerance
  bool uniform
    ()
    const
    {
      if( elem_len.empty() ) return true;
      double const w0 = elem_len.front();
      double const tol = 64. * DBL_EPSILON * std::max( 1., std::fabs( w0 ) );
      for( double const w : elem_len ) if( std::fabs( w - w0 ) > tol ) return false;
      return true;
    }

  //! @brief Lower bound of finite element <a>iel</a>
  double elem_lo
    ( size_t const iel )
    const
    { if( iel >= n_elem ) throw Exceptions::ELEMENTS; return elem_bnd[iel]; }

  //! @brief Upper bound of finite element <a>iel</a>
  double elem_up
    ( size_t const iel )
    const
    { if( iel >= n_elem ) throw Exceptions::ELEMENTS; return elem_bnd[iel+1]; }

  //! @brief Length of finite element <a>iel</a>
  double elem_width
    ( size_t const iel )
    const
    { if( iel >= n_elem ) throw Exceptions::ELEMENTS; return elem_len[iel]; }

protected:

  //! @brief Validate finite-element data
  void _validate
    ()
    const
    {
      if( lo_dom >= up_dom ) throw Exceptions::BOUNDS;
      if( !n_elem )          throw Exceptions::ELEMENTS;
      if( !n_node )          throw Exceptions::NODES;
      if( elem_bnd.size() != n_elem + 1 || elem_len.size() != n_elem )
        throw Exceptions::ELEMENTS;
      double const tol = 64. * DBL_EPSILON
                       * std::max( 1., std::max( std::fabs(lo_dom), std::fabs(up_dom) ) );
      if( std::fabs( elem_bnd.front() - lo_dom ) > tol
       || std::fabs( elem_bnd.back () - up_dom ) > tol )
        throw Exceptions::BOUNDS;
      for( size_t i=0; i<n_elem; ++i ){
        if( elem_bnd[i+1] <= elem_bnd[i] ) throw Exceptions::BOUNDS;
        if( elem_len[i] <= 0. )            throw Exceptions::ELEMENTS;
        if( std::fabs( elem_len[i] - ( elem_bnd[i+1] - elem_bnd[i] ) ) > tol )
          throw Exceptions::ELEMENTS;
      }
    }

  //! @brief Set uniform finite-element boundaries
  void _set_uniform
    ( double const& lo, double const& up, size_t const ne )
    {
      lo_dom = lo; up_dom = up; n_elem = ne;
      if( lo_dom >= up_dom ) throw Exceptions::BOUNDS;
      if( !n_elem )          throw Exceptions::ELEMENTS;
      w_elem = ( up_dom - lo_dom ) / static_cast<double>( n_elem );
      elem_bnd.resize( n_elem + 1 );
      elem_len.assign( n_elem, w_elem );
      for( size_t i=0; i<=n_elem; ++i )
        elem_bnd[i] = lo_dom + w_elem * static_cast<double>( i );
      elem_bnd.back() = up_dom;
    }

  //! @brief Set nonuniform finite-element boundaries
  void _set_boundaries
    ( std::vector<double> const& bnd )
    {
      if( bnd.size() < 2 ) throw Exceptions::ELEMENTS;
      elem_bnd = bnd;
      lo_dom = elem_bnd.front();
      up_dom = elem_bnd.back();
      n_elem = elem_bnd.size() - 1;
      elem_len.resize( n_elem );
      for( size_t i=0; i<n_elem; ++i )
        elem_len[i] = elem_bnd[i+1] - elem_bnd[i];
      w_elem = ( up_dom - lo_dom ) / static_cast<double>( n_elem );
      _validate();
    }

  //! @brief Set nonuniform finite-element lengths
  void _set_lengths
    ( double const& lo, std::vector<double> const& len )
    {
      if( len.empty() ) throw Exceptions::ELEMENTS;
      elem_len = len;
      elem_bnd.resize( elem_len.size() + 1 );
      elem_bnd[0] = lo;
      for( size_t i=0; i<elem_len.size(); ++i ){
        if( elem_len[i] <= 0. ) throw Exceptions::ELEMENTS;
        elem_bnd[i+1] = elem_bnd[i] + elem_len[i];
      }
      lo_dom = elem_bnd.front();
      up_dom = elem_bnd.back();
      n_elem = elem_len.size();
      w_elem = ( up_dom - lo_dom ) / static_cast<double>( n_elem );
      _validate();
    }

public:

  //! @brief Node construction
  bool set_nodes
    ( double const& TOLNODE = 1E-7, double const& TOLEVAL = DBL_EPSILON )
    {
      if( !set_lgnodes( type, n_node, -1.0, 1.0, TOLNODE ) // collocation nodes are in [-1,1], not in [lo_dom,up_dom]
       || !set_lagrange( TOLEVAL ) ) return false;
      return true;
    }
};

inline std::ostream&
operator<<
( std::ostream& os, OCDom const& opt )
{
  os << std::scientific << std::setprecision(4)
     << "[" << opt.lo_dom << " : " << opt.up_dom << "], TYP=";
  switch( opt.type ){
    default:
    case OCDom::TYPE::LG:  os << "LG";  break;
    case OCDom::TYPE::LGR: os << "LGR"; break;
    case OCDom::TYPE::LGL: os << "LGL"; break;
    case OCDom::TYPE::CGL: os << "CGL"; break;
  }
  os << ", ELE=" << opt.n_elem << ", NOD=" << opt.n_node;
  return os;
}

//! @brief Base class for propagation of collocation coefficients in distributed subexpressions through a DAG using MC++
////////////////////////////////////////////////////////////////////////
//! mc::OCBase is a C++ based class defining the environment for
//! propagation of collocation coefficients in distributed
//! subexpressions through a DAG using MC++. 
//! OCBase deliberately contains only the domain map and tensor-product
//! helper routines needed by the collocation arithmetic.  Full model
//! registration, setup, residual evaluation, derivative caching, and
//! interface handling live in OCEnv (see ocenv.hpph).
////////////////////////////////////////////////////////////////////////
class OCBase
{
  template <typename U> friend class OCVar;

public:
  //! @brief Enumeration variable type
  enum VarType
  {
    DOMAIN = 0,  //!< Domain/independent
    STATE,       //!< State
    INPUT,       //!< Input
    CONSTANT     //!< Constant
  };

  //! @brief Exceptions shared by OCBase and derived OCEnv.
  class Exceptions
  {
  public:
    enum TYPE{
      SETUP = 1,   //!< Incomplete setup before evaluation
      INDEX,       //!< Index mismatch in data access
      DOMAIN,      //!< Variable-domain misspecification in expression
      ENV,         //!< Environment mismatch in collocation operation
      CSTVAL,      //!< Undefined constant values
      UNDEF,       //!< Undefined collocation operation
      INTERNAL     //!< Internal error
    };
    Exceptions( TYPE ierr=UNDEF ) : _ierr( ierr ){}
    int ierr(){ return _ierr; }
    std::string what(){
      switch( _ierr ){
      case SETUP:
        return "OCBase::Exceptions  Incomplete setup before evaluation";
      case INDEX:
        return "OCBase::Exceptions  Index mismatch in data access";
      case DOMAIN:
        return "OCBase::Exceptions  Variable-domain misspecification in expression";
      case ENV:
        return "OCBase::Exceptions  Environment mismatch in collocation operation";
      case CSTVAL:
        return "OCBase::Exceptions  Undefined constant values";
      case UNDEF:
        return "OCBase::Exceptions  Undefined collocation operation";
      case INTERNAL:
      default:
        return "OCBase::Exceptions  Internal error";
      }
    }
  private:
    TYPE _ierr;
  };

  typedef std::map< FFVar, OCDom, lt_FFVar >                         t_Dom;
  //typedef std::map< FFVar, t_Dom, lt_FFVar >                         t_Coll;

  virtual ~OCBase() {}

  //! @brief Return the active local domain map.
  t_Dom const& var_domain() const { return _mDom; }

  //! @brief Compute the stride for evaluating collocated expressions within multidimensional stacked arrays.
  static std::pair<size_t, size_t> stride
    ( std::map<FFVar const*, OCDom const*, lt_FFVar> const& dom,
      FFVar const* index, bool const check=true );

  //! @brief Lift collocation array from one tensor-product domain map to another.
  template <typename T>
  static bool lift_colloc
    ( std::vector<T>& coeff_out, const std::vector<T>& coeff_in,
      const std::map<FFVar const*, OCDom const*, lt_FFVar>& dom_out,
      const std::map<FFVar const*, OCDom const*, lt_FFVar>& dom_in,
      const bool check = true );

protected:
  //! @brief map of active local domain variables.
  t_Dom _mDom;
};

//! @brief Arithmetic for propagation of collocation coefficients in distributed subexpressions through a DAG using MC++
////////////////////////////////////////////////////////////////////////
//! mc::OCVar is a C++ class implementing an arithmetic for
//! propagation of collocation coefficients in distributed
//! subexpressions through a DAG using MC++. The template parameter
//! is the evaluation type for collocation coefficients.
////////////////////////////////////////////////////////////////////////
template <typename T>
class OCVar
////////////////////////////////////////////////////////////////////////
{
  friend FFPartial;
  friend FFIntegral;

  template <typename U> friend std::ostream& operator<<
    ( std::ostream &, OCVar<U> const& );
  template <typename U> friend class SupModel;

  template <typename U> friend std::ostream& operator<<
    ( std::ostream &, OCVar<U> const& );

  template <typename U> friend OCVar<U> operator+
    ( OCVar<U> const& );
  template <typename U> friend OCVar<U> && operator+
    ( OCVar<U> && );

  template <typename U> friend OCVar<U> operator+
    ( OCVar<U> const&, OCVar<U> const& );
  template <typename U> friend OCVar<U> && operator+
    ( OCVar<U> const&, OCVar<U> && );
  template <typename U> friend OCVar<U> && operator+
    ( OCVar<U> &&, OCVar<U> const& );    
  template <typename U> friend OCVar<U> && operator+
    ( OCVar<U> &&, OCVar<U> && );    

  template <typename U> friend OCVar<U> operator+
    ( double const&, OCVar<U> const& );
  template <typename U> friend OCVar<U> && operator+
    ( double const&, OCVar<U> && );

  template <typename U> friend OCVar<U> operator+
    ( OCVar<U> const&, double const& );
  template <typename U> friend OCVar<U> && operator+
    ( OCVar<U> &&, double const& );

  template <typename U> friend OCVar<U> operator-
    ( OCVar<U> const& );
  template <typename U> friend OCVar<U> && operator-
    ( OCVar<U> && );

  template <typename U> friend OCVar<U> operator-
    ( OCVar<U> const&, OCVar<U> const& );
  template <typename U> friend OCVar<U> && operator-
    ( OCVar<U> const&, OCVar<U> && );
  template <typename U> friend OCVar<U> && operator-
    ( OCVar<U> &&, OCVar<U> const& );
  template <typename U> friend OCVar<U> && operator-
    ( OCVar<U> &&, OCVar<U> && );

  template <typename U> friend OCVar<U> operator-
    ( double const&, OCVar<U> const& );
  template <typename U> friend OCVar<U> && operator-
    ( double const&, OCVar<U> && );

  template <typename U> friend OCVar<U> operator-
    ( OCVar<U> const&, double const& );
  template <typename U> friend OCVar<U> && operator-
    ( OCVar<U> &&, double const& );

  template <typename U> friend OCVar<U> operator*
    ( OCVar<U> const&, OCVar<U> const& );
  template <typename U> friend OCVar<U> && operator*
    ( OCVar<U> const&, OCVar<U> && );
  template <typename U> friend OCVar<U> && operator*
    ( OCVar<U> &&, OCVar<U> const& );
  template <typename U> friend OCVar<U> && operator*
    ( OCVar<U> &&, OCVar<U> && );

  template <typename U> friend OCVar<U> operator*
    ( double const&, OCVar<U> const& );
  template <typename U> friend OCVar<U> && operator*
    ( double const&, OCVar<U> && );

  template <typename U> friend OCVar<U> operator*
    ( OCVar<U> const&, double const& );
  template <typename U> friend OCVar<U> && operator*
    ( OCVar<U> &&, double const& );

  template <typename U> friend OCVar<U> operator/
    ( OCVar<U> const&, OCVar<U> const& );
  template <typename U> friend OCVar<U> && operator/
    ( OCVar<U> const&, OCVar<U> && );
  template <typename U> friend OCVar<U> && operator/
    ( OCVar<U> &&, OCVar<U> const& );
  template <typename U> friend OCVar<U> && operator/
    ( OCVar<U> &&, OCVar<U> && );

  template <typename U> friend OCVar<U> operator/
    ( double const&, OCVar<U> const& );
  template <typename U> friend OCVar<U> && operator/
    ( double const&, OCVar<U> && );

  template <typename U> friend OCVar<U> operator/
    ( OCVar<U> const&, double const& );
  template <typename U> friend OCVar<U> && operator/
    ( OCVar<U> &&, double const& );

  template <typename U> friend OCVar<U> inv
    ( OCVar<U> const& );
  template <typename U> friend OCVar<U> && inv
    ( OCVar<U> && );
  template <typename U> friend OCVar<U> max
    ( OCVar<U> const&, OCVar<U> const& );
  template <typename U> friend OCVar<U> max
    ( OCVar<U> const&, double const& );
  template <typename U> friend OCVar<U> && max
    ( OCVar<U> &&, double const& );
  template <typename U> friend OCVar<U> min
    ( OCVar<U> const&, OCVar<U> const& );
  template <typename U> friend OCVar<U> min
    ( OCVar<U> const&, double const& );
  template <typename U> friend OCVar<U> && min
    ( OCVar<U> &&, double const& );
  template <typename U> friend OCVar<U> sqr
    ( OCVar<U> const& );
  template <typename U> friend OCVar<U> && sqr
    ( OCVar<U> && );
  template <typename U> friend OCVar<U> sqrt
    ( OCVar<U> const& );
  template <typename U> friend OCVar<U> && sqrt
    ( OCVar<U> && );
  template <typename U> friend OCVar<U> fabs
    ( OCVar<U> const& );
  template <typename U> friend OCVar<U> && fabs
    ( OCVar<U> && );  
  template <typename U> friend OCVar<U> exp
    ( OCVar<U> const& );
  template <typename U> friend OCVar<U> && exp
    ( OCVar<U> && );
  template <typename U> friend OCVar<U> log
    ( OCVar<U> const& );
  template <typename U> friend OCVar<U> && log
    ( OCVar<U> && );
  template <typename U> friend OCVar<U> xlog
    ( OCVar<U> const& );
  template <typename U> friend OCVar<U> && xlog
    ( OCVar<U> && );
  template <typename U> friend OCVar<U> cos
    ( OCVar<U> const& );
  template <typename U> friend OCVar<U> && cos
    ( OCVar<U> && );
  template <typename U> friend OCVar<U> sin
    ( OCVar<U> const& );
  template <typename U> friend OCVar<U> && sin
    ( OCVar<U> && );
  template <typename U> friend OCVar<U> tan
    ( OCVar<U> const& );
  template <typename U> friend OCVar<U> && tan
    ( OCVar<U> && );
  template <typename U> friend OCVar<U> asin
    ( OCVar<U> const& );
  template <typename U> friend OCVar<U> && asin
    ( OCVar<U> && );
  template <typename U> friend OCVar<U> acos
    ( OCVar<U> const& );
  template <typename U> friend OCVar<U> && acos
    ( OCVar<U> && );
  template <typename U> friend OCVar<U> atan
    ( OCVar<U> const& );
  template <typename U> friend OCVar<U> && atan
    ( OCVar<U> && );
  template <typename U> friend OCVar<U> sinh
    ( OCVar<U> const& );
  template <typename U> friend OCVar<U> && sinh
    ( OCVar<U> && );
  template <typename U> friend OCVar<U> cosh
    ( OCVar<U> const& );
  template <typename U> friend OCVar<U> && cosh
    ( OCVar<U> && );
  template <typename U> friend OCVar<U> tanh
    ( OCVar<U> const& );
  template <typename U> friend OCVar<U> && tanh
    ( OCVar<U> && );
  template <typename U> friend OCVar<U> erf
    ( OCVar<U> const& );
  template <typename U> friend OCVar<U> && erf
    ( OCVar<U> && );
  template <typename U> friend OCVar<U> erfc
    ( OCVar<U> const& );
  template <typename U> friend OCVar<U> && erfc
    ( OCVar<U> && );
  template <typename U> friend OCVar<U> fstep
    ( OCVar<U> const& );
  template <typename U> friend OCVar<U> && fstep
    ( OCVar<U> && );
  template <typename U> friend OCVar<U> bstep
    ( OCVar<U> const& );
  template <typename U> friend OCVar<U> && bstep
    ( OCVar<U> && );
  template <typename U> friend OCVar<U> pow
    ( OCVar<U> const&, int const& );
  template <typename U> friend OCVar<U> && pow
    ( OCVar<U> &&, int const& );
  template <typename U, typename V> friend OCVar<U> pow
    ( OCVar<U> const&, V const& );
  template <typename U, typename V> friend OCVar<U> && pow
    ( OCVar<U> &&, V const& );
  template <typename U> friend OCVar<U> cheb
    ( OCVar<U> const&, unsigned int const& n );
  template <typename U> friend OCVar<U> && cheb
    ( OCVar<U> &&, unsigned int const& n );

public:

  typedef std::map< FFVar const*, OCDom const*, lt_FFVar >   t_Dom;
 
private:

  //! @brief Pointer to collocation propagation environment
  OCBase const*        _env;

  //! @brief Pointer to variable domain
  t_Dom               _dom;

  //! @brief Current finite-element index for each distributed domain in this local coefficient block
  std::map< FFVar const*, size_t, lt_FFVar > _elem;

  //! @brief Array of collocation coefficients
  std::vector<T>      _coef;

  //! @brief Update coefficients in distributed variable
  template <typename U>
  OCVar& _update
    ( U const* coef );

  //! @brief Set current finite-element indices for local block scaling
  void _set_elem
    ( std::map<FFVar,size_t,lt_FFVar> const& ndx_el );

  //! @brief Element width for a domain in this local block
  double _elem_width
    ( FFVar const* pvar )
    const;

  //! @brief Set distributed variable
  OCVar& _set
    ( OCBase const& env, std::set< FFVar, lt_FFVar > const& dep, T const* coef );

  //! @brief Set distributed variable with coefficient move
  void _set
    ( OCBase const* env, t_Dom const& dep, std::vector<T> && coef );

public:

  //! @brief Constructor for a real constant
  OCVar
    ( double const& cst=0. )
    : _env( nullptr ),
      _coef( 1, T(cst) )
    {}

  //! @brief Constructor for a lumped coefficient of type T
  template <typename U=T, 
            typename std::enable_if<!std::is_same<typename std::remove_cv<U>::type, double>::value, int>::type=0>
  OCVar
    ( T const& coef )
    : _env( nullptr ),
      _coef( {coef} )
    {}

  //! @brief Constructor for a lumped variable
  OCVar
    ( OCBase const& env, T const& coef )
    : _env( &env ),
      _coef( {coef} )
    {}

  //! @brief Constructor for a distributed variable
  OCVar
    ( OCBase const& env, std::set< FFVar, lt_FFVar > const& dep, T const* coef )
    : _env( &env )
    {
      _set( env, dep, coef );
    }

  //! @brief Copy constructor
  OCVar
    ( OCVar const& other )
    : _env  ( other._env  ),
      _dom  ( other._dom  ),
      _elem ( other._elem ),
      _coef ( other._coef )
    {}

  //! @brief Move constructor
  OCVar
    ( OCVar && other )
    : _env  ( other._env ),
      _dom  ( std::move( other._dom ) ),
      _elem ( std::move( other._elem ) ),
      _coef ( std::move( other._coef ) )
    {}

  //! @brief Destructor
  virtual ~OCVar()
    {}

  //! @brief Set as constant
  OCVar<T>& set
    ( double const& cst )
    {
      _env  = nullptr;
      _dom.clear();
      _elem.clear();
      _coef.assign( 1, T(cst) );
      return *this;
    }

  //! @brief Set as lumped variable
  OCVar<T>& set
    ( OCBase const& env, T const& coef )
    {
      _env  = &env;
      _dom.clear();
      _elem.clear();
      _coef = {coef};
      return *this;
    }

  //! @brief Set as distributed variable
  OCVar<T>& set
    ( OCBase const& env, std::set< FFVar, lt_FFVar > const& dep, T const* coef=nullptr )
    {
      return _set( env, dep, coef );
    }

  //! @brief Set as distributed variable with an explicit collocation domain map
  OCVar<T>& set
    ( OCBase const& env, t_Dom const& dep, T const* coef=nullptr )
    {
      size_t ncoef = 1;
      for( auto const& [pvar,pdom] : dep ){
        if( !pvar || !pdom || !pdom->n_node )
          throw std::runtime_error( "OCVar::set: invalid collocation domain" );
        ncoef *= pdom->n_node;
      }

      _env = &env;
      _dom = dep;
      _elem.clear();
      std::vector<T>().swap( _coef );
      if( coef ) _coef.assign( coef, coef+ncoef );
      else       _coef.resize( ncoef );
      return *this;
    }

  //! @brief Update coefficients in distributed variable
  template <typename U>
  OCVar<T>& update
    ( U const* coef )
    {
      return _update( coef );
    }

  //! @brief Update coefficients in distributed variable and record local finite-element indices
  template <typename U>
  OCVar<T>& update
    ( U const* coef, std::map<FFVar,size_t,lt_FFVar> const& ndx_el )
    {
      _set_elem( ndx_el );
      return _update( coef );
    }

  //! @brief Update coefficients in distributed variable
  template <typename U>
  OCVar<T>& update
    ( std::vector<U> const& coef )
    {
      return _update( coef.data() );
    }

  //! @brief Update coefficients in distributed variable and record local finite-element indices
  template <typename U>
  OCVar<T>& update
    ( std::vector<U> const& coef, std::map<FFVar,size_t,lt_FFVar> const& ndx_el )
    {
      _set_elem( ndx_el );
      return _update( coef.data() );
    }

  //! @brief Assignment operator for constant
  OCVar<T>& operator=
    ( double const& cst )
    {
      return set( cst );
    }

  //! @brief Assigment operator for variable
  OCVar<T>& operator=
    ( OCVar<T> const& other )
    {
      _env  = other._env;
      _dom  = other._dom;
      _elem = other._elem;
      _coef = other._coef;
      return *this;
    }

  //! @brief Move operator for variable
  OCVar<T>& operator=
    ( OCVar<T> && other )
    {
      _env  = other._env;
      _dom  = std::move( other._dom );
      _elem = std::move( other._elem );
      _coef = std::move( other._coef );
      return *this;
    }

  template <typename U=T, 
            typename std::enable_if<!std::is_same<typename std::remove_cv<U>::type, double>::value, int>::type=0>
  OCVar<T>& operator+=
    ( T const& cst )
    {
      for( auto& val : _coef )
        val += cst;
      return *this;
    }
  template <typename U=T, 
            typename std::enable_if<!std::is_same<typename std::remove_cv<U>::type, double>::value, int>::type=0>
  OCVar<T>& operator-=
    ( T const& cst )
    {
      for( auto& val : _coef )
        val -= cst;
      return *this;
    }
  template <typename U=T, 
            typename std::enable_if<!std::is_same<typename std::remove_cv<U>::type, double>::value, int>::type=0>
  OCVar<T>& operator*=
    ( T const& cst )
    {
      for( auto& val : _coef )
        val *= cst;
      return *this;
    }
  template <typename U=T, 
            typename std::enable_if<!std::is_same<typename std::remove_cv<U>::type, double>::value, int>::type=0>
  OCVar<T>& operator/=
    ( T const& cst )
    {
      for( auto& val : _coef )
        val /= cst;
      return *this;
    }
  
  OCVar<T>& operator+=
    ( double const& );
  OCVar<T>& operator+=
    ( OCVar<T> const& );
  //OCVar<T>& operator+=
  //  ( OCVar<T> && );

  OCVar<T>& operator-=
    ( double const& );
  OCVar<T>& operator-=
    ( OCVar<T> const& );
  OCVar<T>& operator-=
    ( OCVar<T> && );
 
  OCVar<T>& operator*=
    ( double const& );
  OCVar<T>& operator*=
    ( OCVar<T> const& );
  //OCVar<T>& operator*=
  //  ( OCVar<T> && );

  OCVar<T>& operator/=
    ( double const& );
  OCVar<T>& operator/=
    ( OCVar<T> const& );
  OCVar<T>& operator/=
    ( OCVar<T> && );

  //! @brief Lift domain
  bool lift
    ( std::set< FFVar, lt_FFVar > const& dep, bool const check=true );

  // Environment
  OCBase const* env
    ()
    const
    { return _env; };

  // Coefficients
  std::vector<T> const& coef
    ()
    const
    { return _coef; };

  // Domain
  t_Dom const& dom
    ()
    const
    { return _dom; };
};

//! @brief C++ class defining partial differentiation as external DAG operations in MC++
////////////////////////////////////////////////////////////////////////
//! mc::FFPartial is a C++ class for defining partial differentiation
//! operations as external DAG operations in MC++.
////////////////////////////////////////////////////////////////////////
class FFPartial
: public FFOp
{
  typedef SMon<FFVar,lt_FFVar> t_SMon;

protected:

  //! @brief Vector of operand variables
  mutable std::vector<FFVar>  _Var;

  //! @brief Map of independent variables
  mutable t_SMon              _Indep;

  FFVar** _set
    ( size_t nVar, FFVar const* pVar, t_SMon const& Indep )
    const
    {
      _Var.assign( pVar, pVar+nVar );
      _Indep = Indep;
      data = nullptr;
      owndata = false;
      return insert_external_operation( *this, nVar, nVar, pVar );
    }

public:

  //! @brief Default constructor
  FFPartial
    ()
    : FFOp( EXTERN )
    {}

  //! @brief Destructor
  virtual ~FFPartial
    ()
    {}

  //! @brief Copy constructor
  FFPartial
    ( FFPartial const& other )
    : FFOp   ( other ),
      _Var   ( other._Var ),
      _Indep ( other._Indep )
    {}

  // Define operation
  std::vector<FFVar> operator()
    ( std::vector<FFVar> const& Var, t_SMon const& Indep )
    {
#ifdef MC__FFPARTIAL_CHECK
      assert( !Var.empty() && Indep.tord );
#endif
      FFVar** ppDer = _set( Var.size(), Var.data(), Indep );
      std::vector<FFVar> Der( Var.size() );
      for( size_t i=0; i<Var.size(); ++i ) Der[i] = *ppDer[i];
      return Der;//std::move( Der );
    }

  FFVar operator()
    ( FFVar const& Var, std::map<FFVar,unsigned,lt_FFVar> const& Indep )
    {
#ifdef MC__FFPARTIAL_CHECK
      assert( Indep.tord );
#endif
      return *(_set( 1, &Var, Indep )[0]);
    }

  FFVar operator()
    ( FFVar const& Var, t_SMon const& Indep )
    {
#ifdef MC__FFPARTIAL_CHECK
      assert( Indep.tord );
#endif
      return *(_set( 1, &Var, Indep )[0]);
    }

  std::vector<FFVar> operator()
    ( std::vector<FFVar> const& Var, FFVar const& Indep )
    {
#ifdef MC__FFPARTIAL_CHECK
      assert( !Var.empty() );
#endif
      FFVar** ppDer = _set( Var.size(), Var.data(), Indep );
      std::vector<FFVar> Der( Var.size() );
      for( size_t i=0; i<Var.size(); ++i ) Der[i] = *ppDer[i];
      return Der;//std::move( Der );
    }

  FFVar operator()
    ( FFVar const& Var, FFVar const& Indep )
    {
      return *(_set( 1, &Var, Indep )[0]);
    }

  t_SMon const& Indep
    ()
    const
    {
      return _Indep;
    }

  std::vector<FFVar> const& Var
    ()
    const
    {
      return _Var;
    }

  // Evaluation overloads
  virtual void feval
    ( std::type_info const& idU, unsigned const nRes, void* vRes, unsigned const nVar,
      void const* vVar, unsigned const* mVar )
    const;

  void eval
    ( size_t const nRes, FFVar* vRes, size_t const nVar, FFVar const* vVar, unsigned const* mVar )
    const;

  void eval
    ( size_t const nRes, FFDep* vRes, size_t const nVar, FFDep const* vVar, unsigned const* mVar )
    const;

  void eval
    ( size_t const nRes, FADType<FFVar>* vRes, size_t const nVar, FADType<FFVar> const* vVar,
      unsigned const* mVar )
    const;

  void eval
    ( size_t const nRes, SLiftVar* vRes, size_t const nVar, SLiftVar const* vVar, unsigned const* mVar )
    const;

  void eval
    ( size_t const nRes, FFExpr* vRes, size_t const nVar, FFExpr const* vVar, unsigned const* mVar )
    const;

  template <typename U>
  void eval
    ( size_t const nRes, OCVar<U>* vRes, size_t const nVar, OCVar<U> const* vVar, unsigned const* mVar )
    const;

  // Derivatives
  //void deriv
  //  ( size_t const nRes, FFVar const* vRes, size_t const nVar, FFVar const* vVar, FFVar** vDer )
  //  const;

  // Ordering
  bool lt
    ( FFOp const* other )
    const;

  // Properties
  std::string name
    ()
    const
    {
      return "Partial" + _Indep.display(0);
    }

  // Commutativity
  bool commutative
    ()
    const
    {
      return false;
    }
};

inline bool
FFPartial::lt
( FFOp const* other )
const
{
#ifdef MC__FFPARTIAL_TRACE
  std::cout << "FFPartial::lt\n";
#endif
  FFPartial const* op = dynamic_cast<FFPartial const*>(other);

  // Compare independent variables
  return lt_SMon<lt_FFVar>()( _Indep, op->_Indep );
}

inline void
FFPartial::feval
( std::type_info const& idU, unsigned const nRes, void* vRes, unsigned const nVar,
  void const* vVar, unsigned const* mVar )
const
{
  if( idU == typeid( FFVar ) )
    return eval( nRes, static_cast<FFVar*>(vRes), nVar, static_cast<FFVar const*>(vVar), mVar );
  else if( idU == typeid( FFDep ) )
    return eval( nRes, static_cast<FFDep*>(vRes), nVar, static_cast<FFDep const*>(vVar), mVar );
  else if( idU == typeid( SLiftVar ) )
    return eval( nRes, static_cast<SLiftVar*>(vRes), nVar, static_cast<SLiftVar const*>(vVar), mVar );
  else if( idU == typeid( FFExpr ) )
    return eval( nRes, static_cast<FFExpr*>(vRes), nVar, static_cast<FFExpr const*>(vVar), mVar );
  else if( idU == typeid( OCVar<double> ) )
    return eval( nRes, static_cast<OCVar<double>*>(vRes), nVar, static_cast<OCVar<double> const*>(vVar), mVar );
  else if( idU == typeid( OCVar<FFDep> ) )
    return eval( nRes, static_cast<OCVar<FFDep>*>(vRes), nVar, static_cast<OCVar<FFDep> const*>(vVar), mVar );
  else if( idU == typeid( OCVar< FADType<double> > ) )
    return eval( nRes, static_cast<OCVar< FADType<double> >*>(vRes),
                 nVar, static_cast<OCVar< FADType<double> > const*>(vVar), mVar );

  throw std::runtime_error( "FFPartial::feval ** No evaluation method for type"+std::string(idU.name())+"\n" );
}

inline void
FFPartial::eval
( size_t const nRes, FFExpr* vRes, size_t const nVar, FFExpr const* vVar,
  unsigned const* mVar )
const
{
#ifdef MC__FFPARTIAL_TRACE
  std::cout << "FFPartial::eval: FFExpr\n";
#endif

  switch( FFExpr::options.LANG ){
   case FFExpr::Options::DAG:
    for( unsigned j=0; j<nRes; ++j ){
      std::ostringstream os; os << name();
      if( nRes > 1 ) os << "[" << j << "]";
      vRes[j] = FFExpr::compose( os.str(), nVar, vVar );
    }
    break;
   case FFExpr::Options::GAMS:
   default:
    throw typename FFExpr::Exceptions( FFExpr::Exceptions::UNDEF );
  }
}

inline void
FFPartial::eval
( size_t const nRes, SLiftVar* vRes, size_t const nVar, SLiftVar const* vVar,
  unsigned const* mVar )
const
{
#ifdef MC__FFPARTIAL_TRACE
  std::cout << "FFPartial::eval: SLiftVar\n";
#endif

  // Lift partial differentiation operation
  vVar->env()->lift( nRes, vRes, nVar, vVar );
}

inline void
FFPartial::eval
( size_t const nRes, FFDep* vRes, size_t const nVar, FFDep const* vVar,
  unsigned const* mVar )
const
{
#ifdef MC__FFPARTIAL_TRACE
  std::cout << "FFPartial::eval: FFDep\n";
#endif
#ifdef MC__FFPARTIAL_CHECK
  assert( nRes == nVar );
#endif

  // Differentiation is a linear operator
  for( size_t i=0; i<nVar; ++i )
    vRes[i] = vVar[i];
}

inline void
FFPartial::eval
( size_t const nRes, FFVar* vRes, size_t const nVar, FFVar const* vVar,
  unsigned const* mVar )
const
{
#ifdef MC__FFPARTIAL_TRACE
  std::cout << "FFPartial::eval: FFVar\n";
#endif
#ifdef MC__FFPARTIAL_CHECK
  assert( nRes == nVar );
#endif

  FFVar** ppRes = _set( nVar, vVar, _Indep );
  for( unsigned j=0; j<nRes; ++j )
    vRes[j] = *(ppRes[j]);
}

template <typename U>
inline void
FFPartial::eval
( size_t const nRes, OCVar<U>* vRes, size_t const nVar, OCVar<U> const* vVar,
  unsigned const* mVar )
const
{
#ifdef MC__FFPARTIAL_TRACE
  std::cout << "FFPartial::eval: OCVar<U>\n";
#endif
#ifdef MC__FFPARTIAL_CHECK
  assert( nRes == nVar );
  bool const check = true;
#else
  bool const check = false;
#endif

  for( size_t i=0; i<nVar; ++i ){
    if( !_Indep.tord ){
      vRes[i] = vVar[i];
      continue;
    }
    
    auto env = vVar[i].env();
    if( !env )
      throw std::runtime_error( "FFPartial::eval<OCVar<U>> ** No collocation environment\n" );

    auto const& din = vVar[i].dom();   // std::map<FFVar const*, OCDom const*, lt_FFVar>
    std::vector<U> yin, yout;
    std::vector<U> const* pin = &vVar[i].coef();

    // Loop over independent variables and apply differentiation of desired order
    for( auto const& [var, ord] : _Indep.expr ){
      // Retrieve independent variable and domain
      auto itd = din.find( &var );
      if( itd == din.end() )
        throw std::runtime_error( "FFPartial::eval<OCVar<U>> ** Independent variable not in domain\n" );
      auto const& [pvar, pdom] = *itd;

      // Get the stride for current independent variable
      auto [stride,outer] = env->stride( din, pvar, check );
      if( !stride || !outer )
        throw std::runtime_error( "FFPartial::eval<OCVar<U>> ** Stride calculation was unsuccessful\n" );

      // Apply differentiation for current independent variable and order
      double const elem_w = vVar[i]._elem_width( pvar );
      double const rescale = ord? std::pow( 2./elem_w, ord ): 1.;
      if( !pdom->eval_lagrange( yout, ord, *pin, stride, outer, rescale ) )
        throw std::runtime_error( "FFPartial::eval<OCVar<U>> ** Differentiation was unsuccessful\n" );
        
      // Swap coefficient vectors for next round of differentiation
      yout.swap( yin );
      pin = &yin;
    }

    // Set resulting collocation variable; keep the local element indices so
    // subsequent nested partial/integral operations use the same block widths.
    vRes[i]._set( env, din, std::move( yin ) );
    vRes[i]._elem = vVar[i]._elem;
  }
}

//! @brief C++ class defining integration as external DAG operations in MC++
////////////////////////////////////////////////////////////////////////
//! mc::FFIntegral is a C++ class for defining integration operations
//! as external DAG operations in MC++.
////////////////////////////////////////////////////////////////////////
class FFIntegral
: public FFOp
{
  typedef SMon<FFVar,lt_FFVar> t_SMon;

protected:

  //! @brief Vector of operand variables
  mutable std::vector<FFVar>   _Var;

  //! @brief Map of independent variables
  mutable t_SMon               _Indep;

  FFVar** _set
    ( size_t nVar, FFVar const* pVar, t_SMon const& Indep )
    const
    {
      _Var.assign( pVar, pVar+nVar );
      _Indep = Indep;
      data = nullptr;
      owndata = false;
      return insert_external_operation( *this, nVar, nVar, pVar );
    }

public:

  //! @brief Default constructor
  FFIntegral
    ( bool const sparse=true )
    : FFOp( EXTERN )
    {}

  //! @brief Destructor
  virtual ~FFIntegral
    ()
    {
      this->sparse = sparse;
    }


  //! @brief Copy constructor
  FFIntegral
    ( FFIntegral const& other )
    : FFOp   ( other ),
      _Var   ( other._Var ),
      _Indep ( other._Indep )
    {}

  // Define operation
  std::vector<FFVar> operator()
    ( std::vector<FFVar> const& Var, t_SMon const& Indep )
    {
#ifdef MC__FFINTEGRAL_CHECK
      assert( !Var.empty() && !Indep.empty() );
#endif
      FFVar** ppDer = _set( Var.size(), Var.data(), Indep );
      std::vector<FFVar> Der( Var.size() );
      for( size_t i=0; i<Var.size(); ++i ) Der[i] = *ppDer[i];
      return Der;//std::move( Der );
    }

  FFVar operator()
    ( FFVar const& Var, t_SMon const& Indep )
    {
#ifdef MC__FFINTEGRAL_CHECK
      assert( Indep.tord );
#endif
      return *(_set( 1, &Var, Indep )[0]);
    }

  std::vector<FFVar> operator()
    ( std::vector<FFVar> const& Var, FFVar const& Indep )
    {
#ifdef MC__FFINTEGRAL_CHECK
      assert( !Var.empty() );
#endif
      FFVar** ppDer = _set( Var.size(), Var.data(), {Indep} );
      std::vector<FFVar> Der( Var.size() );
      for( size_t i=0; i<Var.size(); ++i ) Der[i] = *ppDer[i];
      return Der;//std::move( Der );
    }

  FFVar operator()
    ( FFVar const& Var, FFVar const& Indep )
    {
      return *(_set( 1, &Var, {Indep} )[0]);
    }

  t_SMon const& Indep
    ()
    const
    {
      return _Indep;
    }

  std::vector<FFVar> const& Var
    ()
    const
    {
      return _Var;
    }

  // Evaluation overloads
  virtual void feval
    ( std::type_info const& idU, unsigned const nRes, void* vRes, unsigned const nVar,
      void const* vVar, unsigned const* mVar )
    const;

  void eval
    ( size_t const nRes, FFVar* vRes, size_t const nVar, FFVar const* vVar, unsigned const* mVar )
    const;
    
  void eval
    ( size_t const nRes, FFDep* vRes, size_t const nVar, FFDep const* vVar, unsigned const* mVar )
    const;

  //void eval
  //  ( size_t const nRes, FADType<FFVar>* vRes, size_t const nVar, FADType<FFVar> const* vVar,
  //    unsigned const* mVar )
  //  const;

  void eval
    ( size_t const nRes, SLiftVar* vRes, size_t const nVar, SLiftVar const* vVar, unsigned const* mVar )
    const;

  void eval
    ( size_t const nRes, FFExpr* vRes, size_t const nVar, FFExpr const* vVar, unsigned const* mVar )
    const;

  template <typename U>
  void eval
    ( size_t const nRes, OCVar<U>* vRes, size_t const nVar, OCVar<U> const* vVar, unsigned const* mVar )
    const;

  // Derivatives
  //void deriv
  //  ( size_t const nRes, FFVar const* vRes, size_t const nVar, FFVar const* vVar, FFVar** vDer )
  //  const;

  // Ordering
  bool lt
    ( FFOp const* other )
    const;

  // Properties
  std::string name
    ()
    const
    {
      return "Integral" + _Indep.display(0);
    }

  // Commutativity
  bool commutative
    ()
    const
    {
      return false;
    }
};

inline bool
FFIntegral::lt
( FFOp const* other )
const
{
#ifdef MC__FFINTEGRAL_TRACE
  std::cout << "FFIntegral::lt\n";
#endif
  FFIntegral const* op = dynamic_cast<FFIntegral const*>(other);

  // Compare independent variables
  return lt_SMon<lt_FFVar>()( _Indep, op->_Indep );
}

inline void
FFIntegral::feval
( std::type_info const& idU, unsigned const nRes, void* vRes, unsigned const nVar,
  void const* vVar, unsigned const* mVar )
const
{
  if( idU == typeid( FFVar ) )
    return eval( nRes, static_cast<FFVar*>(vRes), nVar, static_cast<FFVar const*>(vVar), mVar );
  else if( idU == typeid( FFDep ) )
    return eval( nRes, static_cast<FFDep*>(vRes), nVar, static_cast<FFDep const*>(vVar), mVar );
  else if( idU == typeid( SLiftVar ) )
    return eval( nRes, static_cast<SLiftVar*>(vRes), nVar, static_cast<SLiftVar const*>(vVar), mVar );
  else if( idU == typeid( FFExpr ) )
    return eval( nRes, static_cast<FFExpr*>(vRes), nVar, static_cast<FFExpr const*>(vVar), mVar );
  else if( idU == typeid( OCVar<double> ) )
    return eval( nRes, static_cast<OCVar<double>*>(vRes), nVar, static_cast<OCVar<double> const*>(vVar), mVar );
  else if( idU == typeid( OCVar<FFDep> ) )
    return eval( nRes, static_cast<OCVar<FFDep>*>(vRes), nVar, static_cast<OCVar<FFDep> const*>(vVar), mVar );
  else if( idU == typeid( OCVar< FADType<double> > ) )
    return eval( nRes, static_cast<OCVar< FADType<double> >*>(vRes), nVar, static_cast<OCVar< FADType<double> > const*>(vVar), mVar );

  throw std::runtime_error( "FFIntegral::feval ** No evaluation method for type"+std::string(idU.name())+"\n" );
}

inline void
FFIntegral::eval
( size_t const nRes, FFExpr* vRes, size_t const nVar, FFExpr const* vVar,
  unsigned const* mVar )
const
{
#ifdef MC__FFINTEGRAL_TRACE
  std::cout << "FFIntegral::eval: FFExpr\n";
#endif

  switch( FFExpr::options.LANG ){
   case FFExpr::Options::DAG:
    for( unsigned j=0; j<nRes; ++j ){
      std::ostringstream os; os << name();
      if( nRes > 1 ) os << "[" << j << "]";
      vRes[j] = FFExpr::compose( os.str(), nVar, vVar );
    }
    break;
   case FFExpr::Options::GAMS:
   default:
    throw typename FFExpr::Exceptions( FFExpr::Exceptions::UNDEF );
  }
}

inline void
FFIntegral::eval
( size_t const nRes, SLiftVar* vRes, size_t const nVar, SLiftVar const* vVar,
  unsigned const* mVar )
const
{
#ifdef MC__FFINTEGRAL_TRACE
  std::cout << "FFIntegral::eval: SLiftVar\n";
#endif

  // Lift integration operation
  vVar->env()->lift( nRes, vRes, nVar, vVar );
}

inline void
FFIntegral::eval
( size_t const nRes, FFDep* vRes, size_t const nVar, FFDep const* vVar,
  unsigned const* mVar )
const
{
#ifdef MC__FFINTEGRAL_TRACE
  std::cout << "FFIntegral::eval: FFDep\n";
#endif
#ifdef MC__FFINTEGRAL_CHECK
  assert( nRes == nVar );
#endif

  // Integration is a linear operator
  for( size_t i=0; i<nVar; ++i )
    vRes[i] = vVar[i];
}

inline void
FFIntegral::eval
( size_t const nRes, FFVar* vRes, size_t const nVar, FFVar const* vVar,
  unsigned const* mVar )
const
{
#ifdef MC__FFINTEGRAL_TRACE
  std::cout << "FFIntegral::eval: FFVar\n";
#endif
#ifdef MC__FFINTEGRAL_CHECK
  assert( nRes == nVar );
#endif

  FFVar** ppRes = _set( nVar, vVar, _Indep );
  for( unsigned j=0; j<nRes; ++j )
    vRes[j] = *(ppRes[j]);
}

template <typename U>
inline void
FFIntegral::eval
( size_t const nRes, OCVar<U>* vRes, size_t const nVar, OCVar<U> const* vVar,
  unsigned const* mVar )
const
{
#ifdef MC__FFINTEGRAL_TRACE
  std::cout << "FFIntegral::eval: OCVar<U>\n";
#endif
#ifdef MC__FFINTEGRAL_CHECK
  assert( nRes == nVar );
  bool const check = true;
#else
  bool const check = false;
#endif

  for( size_t i=0; i<nVar; ++i ){
    if( !_Indep.tord ){
      vRes[i] = vVar[i];
      continue;
    }
    
    auto env = vVar[i].env();
    if( !env )
      throw std::runtime_error( "FFIntegral::eval<OCVar<U>> ** No collocation environment\n" );

    auto const& din = vVar[i].dom();   // std::map<FFVar const*, OCDom const*, lt_FFVar>
    auto dout = din;
    std::vector<U> yin, yout;
    std::vector<U> const* pin = &vVar[i].coef();

    // Loop over independent variables and apply differentiation of desired order
    for( auto const& [var, ord] : _Indep.expr ){
      if( ord > 1 )
        throw std::runtime_error( "FFIntegral::eval<OCVar<U>> ** Integration order larger than one\n" );

      // Retrieve independent variable and domain
      auto itd = dout.find( &var );
      if( itd == dout.end() )
        throw std::runtime_error( "FFIntegral::eval<OCVar<U>> ** Independent variable not in domain\n" );
      auto const& [pvar, pdom] = *itd;

      // Get the stride for current independent variable
      auto [stride,outer] = env->stride( dout, pvar, check );
      if( !stride || !outer )
        throw std::runtime_error( "FFIntegral::eval<OCVar<U>> ** Stride calculation was unsuccessful\n" );

      // Apply integration for current independent variable
      double const elem_w = vVar[i]._elem_width( pvar );
      double const rescale = elem_w / 2.;
      if( !pdom->quad_lagrange( yout, *pin, stride, outer, rescale ) )
        throw std::runtime_error( "FFIntegral::eval<OCVar<U>> ** Integration was unsuccessful\n" );

      // Drop integrated variable from dependency set
      dout.erase( itd ); 
      
      // Swap coefficient vectors for next round of differentiation
      yout.swap( yin );
      pin = &yin;
    }

    // Set resulting collocation variable; retain element indices for the
    // remaining, non-integrated domains.
    vRes[i]._set( env, dout, std::move( yin ) );
    vRes[i]._elem = vVar[i]._elem;
    for( auto const& [var,ord] : _Indep.expr )
      vRes[i]._elem.erase( &var );
  }
}

////////////////////////////////////////////////////////////////////////

inline void
BASE_OC::_lgnodes
( size_t const N, double const& TOL )
{
  // Trivial cases
  if( !N ){
    _x.clear();
    _w.clear();
    _P.clear();
    return;
  }

  _x.resize( N );
  _w.resize( N );

  if( N == 1 ){
    _x[0] = 0.0;
    _w[0] = 2.0;
    return;
  }

  // Chebyshev–Gauss initial guess for Legendre–Gauss nodes
  // x_i ≈ -cos( pi (2i+1) / (2N) ), i=0..N-1
  for( size_t i = 0; i < N; ++i )
    _x[i] = -std::cos( PI * ( (2.0*i + 1.0) / (2.0*N) ) );

  // Newton iteration on f(x) = P_N(x); f'(x) = P'_N(x)
  // Recurrences:
  //   P_{k+1}(x) = ((2k+1)x P_k - k P_{k-1})/(k+1)
  //   P'_N(x)    = N * ( P_{N-1}(x) - x P_N(x) ) / (1 - x^2)
  _P.resize( N+1 );
  for( double err = TOL + 1.0; err > TOL; ){
    err = 0.0;

    for( size_t i = 0; i < N; ++i ){

      // Build P_0..P_N at x=xi
      _P[0] = 1.0;
      _P[1] = _x[i];
      for( size_t k = 1; k < N; ++k )
        _P[k+1] = ( (2.0*k + 1.0)*_x[i]*_P[k] - k*_P[k-1] ) / (k + 1.0);

      // Derivative via stable identity
      auto dP = [&](size_t n) -> double
      {
        if( n == 0 ) return 0.0;
        const double denom = std::max( DBL_EPSILON, 1.0 - _x[i]*_x[i] );
        return n * ( _P[n-1] - _x[i]*_P[n] ) / denom;
      };

      // Newton step
      const double dx = _P[N] / dP(N);
      _x[i] -= dx;

      // Error monitor
      const double adx = std::abs(dx);
      if( adx > err ) err = adx;
    }
  }

  // Ensure ascending order [-1,1]
  std::sort( _x.begin(), _x.end() );
 
  // LG weights formula: 
  // w_i = 2 * (1 - x_i^2) / ( N^2 * [ P_{N-1}(x_i) ]^2 )
  for( size_t i = 0; i < N; ++i ){
    // Update Legendre polynomials up to P_N(xi)
    _P[0] = 1.0;
    _P[1] = _x[i];
    for( size_t k = 1; k < N-1; ++k )
      _P[k+1] = ( (2.0*k + 1.0)*_x[i]*_P[k] - k*_P[k-1] ) / (k + 1.0);

    _w[i] = 2.0 * std::max( 0.0, 1.0 - _x[i]*_x[i] ) / ( N*N * (_P[N-1]*_P[N-1]) );
  }
}


/*
  // -------------------------
  // Compute Gauss–Legendre weights
  // -------------------------
  const double N2 = static_cast<double>(N) * static_cast<double>(N);
  for( size_t i = 0; i < N; ++i ){
    const double xi = _x[i];

    // Evaluate P_{N-1}(xi) by recurrence
    double Pkm1 = 1.0;   // P_0
    double Pk   = xi;    // P_1
    if( N > 2 ){
      for( size_t k = 1; k < N-1; ++k ){
        const double Pkp1 = ( (2.0*k + 1.0)*xi*Pk - k*Pkm1 ) / (k + 1.0);
        Pkm1 = Pk;
        Pk   = Pkp1;
      }
    }
    // For N=2, loop skipped and Pk = P_1(xi) = xi, which is P_{N-1}
    const double PN1 = (N == 2) ? xi : Pk;

    const double one_m_x2 = std::max( 0.0, 1.0 - xi*xi ); // guard tiny negatives
    _w[i] = 2.0 * one_m_x2 / ( N2 * (PN1 * PN1) );
  }
*/
inline void
BASE_OC::_lgrnodes
( size_t const N, double const& TOL )
{
  // Trivial cases
  if( !N ){
    _x.clear();
    _w.clear();
    _P.clear();
    return;
  }

  _x.resize( N );
  _w.resize( N );

  if( N == 1 ){
    // Single-node left Radau: node at -1 with weight 2
    _x[0] = -1.0;
    _w[0] =  2.0;
    return;
  }

  // Chebyshev–Radau initial guess (-1 included)
  // x_i = -cos( pi (2i) / (2N-1) ), i=0..N-1
  for( size_t i = 0; i < N; ++i )
    _x[i] = -std::cos( (2.0*i) * PI / (2.0*N - 1.0) );

  // Enforce exact left endpoint
  _x[0] = -1.0;

  // We need P_k(x) up to k = N at each xi
  _P.resize( N+1 );
  for( double err = TOL + 1.0; err > TOL; ){
    err = 0.0;

    // Update only interior nodes (indices 1..N-1)
    for( size_t i = 1; i < N; ++i ){

      // Build Legendre polynomials up to P_N(xi)
      _P[0] = 1.0;
      _P[1] = _x[i];
      for( size_t k = 1; k < N; ++k )
        _P[k+1] = ( (2.0*k + 1.0)*_x[i]*_P[k] - k*_P[k-1] ) / (k + 1.0);

      // Derivatives via stable identity:
      // P'_n(x) = n * ( P_{n-1}(x) - x P_n(x) ) / (1 - x^2)
      auto dP = [&](size_t n) -> double
      {
        if( n == 0 ) return 0.0;
        const double denom = std::max( DBL_EPSILON, 1.0 - _x[i]*_x[i] );
        return n * ( _P[n-1] - _x[i]*_P[n] ) / denom;
      };

      // Newton step for left Radau nodes:
      // f(x)  = P_N(x) + P_{N-1}(x)
      // f'(x) = P'_N(x) + P'_{N-1}(x)
      const double f  = _P[N]  + _P[N-1];
      const double fp = dP(N)  + dP(N-1);
      const double dx = f / fp;
      _x[i] -= dx;

      // Keep node inside (-1,1) (simple safeguard)
      if( _x[i] <= -1.0 ) _x[i] = -1.0 + 1e2*DBL_EPSILON;
      if( _x[i] >=  1.0 ) _x[i] =  1.0 - 1e2*DBL_EPSILON;

      // Update error estimate
      const double adx = std::abs(dx);
      if( adx > err ) err = adx;
    }
  }

  // Ensure ascending order [-1,1]
  std::sort( _x.begin(), _x.end() );

  // LGR weights formula: 
  // w_0 = 2 / N^2
  // For i = 1..N-1: w_i = (1 - x_i) / ( N^2 * [ P_{N-1}(x_i) ]^2 )
  _w[0] = 2.0 / (N*N);
  for( size_t i = 1; i < N; ++i ){
    // Update Legendre polynomials up to P_N(xi)
    _P[0] = 1.0;
    _P[1] = _x[i];
    for( size_t k = 1; k < N-1; ++k )
      _P[k+1] = ( (2.0*k + 1.0)*_x[i]*_P[k] - k*_P[k-1] ) / (k + 1.0);
  
    _w[i] = (1.0 - _x[i]) / ( N*N * (_P[N-1]*_P[N-1]) );
  }
}

inline void
BASE_OC::_lglnodes
( size_t const N, double const& TOL )
{
  // Trivial cases
  if( !N ){
    _x.clear();
    _w.clear();
    _P.clear();
    return;
  }

  _x.resize( N );
  _w.resize( N );

  if( N == 1 ){
    _x[0] = 0.0;
    _w[0] = 2.0;                 // integral over [-1,1]
    return;
  }
  if( N == 2 ){
    _x[0] = -1.0; _x[1] =  1.0;
    _w[0] =  1.0; _w[1] =  1.0;  // 2 / (1*2) = 1
    return;
  }

  // Chebyshev–Lobatto initial guess (-1 and +1 included)
  for( size_t i = 0; i < N; ++i )
    _x[i] = -std::cos( i * PI / (N - 1.0) );

  // Enforce exact endpoints
  _x[0]   = -1.0;
  _x[N-1] =  1.0;

  // Iterate Newton correction on interior nodes: roots of P'_{N-1}(x)
  const size_t NN = N - 1; // degree used for LGL (P_NN with NN=N-1)
  _P.resize( N );
  for( double err = TOL + 1.0; err > TOL; ){
    err = 0.0;

    // Update only interior nodes (indices 1..N-2)
    for( size_t i = 1; i + 1 < N; ++i ){

      // Build Legendre polynomials up to P_NN(xi)
      _P[0] = 1.0;
      _P[1] = _x[i];
      for( size_t k = 1; k < NN; ++k )
        _P[k+1] = ( (2.0*k + 1.0)*_x[i]*_P[k] - k*_P[k-1] ) / (k + 1.0);

      // Derivatives via stable identities:
      // P'_n(x) = n * ( P_{n-1}(x) - x P_n(x) ) / (1 - x^2)
      // P''_n(x) = (2x P'_n(x) - n(n+1) P_n(x)) / (1 - x^2)
      auto dP_and_d2P = [&](size_t n) -> std::pair<double,double>
      {
        if( n == 0 ) return { 0.0, 0.0 };
        const double denom = std::max( DBL_EPSILON, 1.0 - _x[i]*_x[i] );
        const double dP  = n * ( _P[n-1] - _x[i]*_P[n] ) / denom;
        const double d2P = ( 2.0*_x[i]*dP - n*(n+1.0)*_P[n] ) / denom;
        return { dP, d2P };
      };

      // Newton step on f(x) = P'_{NN}(x)
      const auto [dP, d2P] = dP_and_d2P(NN);
      const double dx = dP / d2P;
      _x[i] -= dx;

      // Error monitor
      const double adx = std::abs(dx);
      if( adx > err ) err = adx;
    }
  }

  // Ensure ascending order [-1,1]
  std::sort( _x.begin(), _x.end() );

  // LGL weights formula: w_i = 2 / (n(n+1)) * 1 / [P_n(x_i)]^2, with n = N-1
  // (Endpoints included; P_n(±1)=±1 so it reduces to 2/(n(n+1)) there)
  double const c = 2.0 / (NN*NN + NN);
  _w[0] = _w[NN] = c;
  for( size_t i = 1; i < NN; ++i ){
    // Update Legendre polynomials up to P_NN(xi)
    _P[0] = 1.0;
    _P[1] = _x[i];
    for( size_t k = 1; k < NN; ++k )
      _P[k+1] = ( (2.0*k + 1.0)*_x[i]*_P[k] - k*_P[k-1] ) / (k + 1.0);

    _w[i] = c / (_P[NN] * _P[NN]);
  }
}


inline void
BASE_OC::_cglnodes
( size_t const N )
{
  // Chebyshev-Gauss-Lobatto nodes (Chebyshev extrema), ordered increasingly.
  // The quadrature weights below are the corresponding Clenshaw-Curtis
  // interpolatory weights for integrating the Lagrange interpolant on [-1,1].
  if( !N ){
    _x.clear();
    _w.clear();
    _P.clear();
    return;
  }

  _x.resize( N );
  _w.assign( N, 0.0 );
  _P.clear();

  if( N == 1 ){
    _x[0] = 0.0;
    _w[0] = 2.0;
    return;
  }

  size_t const n = N - 1;
  for( size_t j = 0; j <= n; ++j )
    _x[j] = -std::cos( PI * static_cast<double>(j) / static_cast<double>(n) );

  if( n == 1 ){
    _w[0] = 1.0;
    _w[1] = 1.0;
    return;
  }

  // Trefethen's clencurt formula, adapted to the ascending node order above.
  if( n % 2 == 0 ){
    _w[0] = _w[n] = 1.0 / ( static_cast<double>(n)*n - 1.0 );
    for( size_t j = 1; j < n; ++j ){
      double v = 1.0;
      double const theta = PI * static_cast<double>(j) / static_cast<double>(n);
      for( size_t k = 1; k <= n/2 - 1; ++k ){
        double const kk = static_cast<double>(k);
        v -= 2.0 * std::cos( 2.0 * kk * theta ) / ( 4.0 * kk * kk - 1.0 );
      }
      v -= std::cos( static_cast<double>(n) * theta )
         / ( static_cast<double>(n)*n - 1.0 );
      _w[j] = 2.0 * v / static_cast<double>(n);
    }
  }
  else{
    _w[0] = _w[n] = 1.0 / ( static_cast<double>(n)*n );
    for( size_t j = 1; j < n; ++j ){
      double v = 1.0;
      double const theta = PI * static_cast<double>(j) / static_cast<double>(n);
      for( size_t k = 1; k <= (n-1)/2; ++k ){
        double const kk = static_cast<double>(k);
        v -= 2.0 * std::cos( 2.0 * kk * theta ) / ( 4.0 * kk * kk - 1.0 );
      }
      _w[j] = 2.0 * v / static_cast<double>(n);
    }
  }
}

inline bool
BASE_OC::set_lgnodes
( int const TYPE, size_t const N, double const& a, double const& b, double const& TOL )
{
  if( !N || b <= a ) return false;

  switch( TYPE ){
    default:
    case  0: _lgnodes( N, TOL );                   break;
    case  1: _lgrnodes( N, TOL );                  break;
    case -1: _lgrnodes( N, TOL );
             std::for_each( _x.begin(), _x.end(), []( double& x ){ x *= -1.; } );
             std::reverse( _x.begin(), _x.end() );
             std::reverse( _w.begin(), _w.end() ); break;
    case  2: _lglnodes( N, TOL );                  break;
    case  3: _cglnodes( N );                       break;
  }
  
  // Rescale from [-1,1] -> [a,b]
  const double half = 0.5*(b-a);
  const double mid  = 0.5*(a+b);
  if( a != -1.0 || b != 1.0 ){
    std::for_each( _x.begin(), _x.end(), [&]( double& xi ){ xi = mid + half*xi; } );
    std::for_each( _w.begin(), _w.end(), [&]( double& wi ){ wi *= half; } );
  }
  return true;
}

inline bool
BASE_OC::set_lagrange
( double const& TOL )
{
  if( _x.empty() ) return false;
  size_t const n = _x.size();

  // Compute w and c (needed for second-derivative limits at nodes)
  _Lw.assign( n, 1.0 );
  _Lc.resize( n );
  for( size_t i = 0; i < n; ++i ){
    double s = 0.0;
    for( size_t j = 0; j < n; ++j ){
      if( i == j ) continue;
      double const d = _x[i] - _x[j];
      if( d <= TOL && d >= -TOL ) return false;
      _Lw[i] /= d;
      s += 1.0 / d;
    }
    _Lc[i] = s;
  }

  return true;
}

inline bool
BASE_OC::set_lagrange_ops
( int const MAXORD )
{
  if( _x.empty() ) return false;
  if( MAXORD < 0 ) return false;

  // Make sure barycentric weights are available
  if( _Lw.empty() || _Lw.size() != _x.size() ){
    if( !set_lagrange() ) return false;
  }

  const size_t n = _x.size();
  const int M = std::min( MAXORD, (int)n-1 ); // derivatives of order >= n are identically zero

  _D.clear();
  _D.resize( (size_t)M + 1 );

  // D^0 = I
  _D[0] = arma::eye( (arma::uword)n, (arma::uword)n );

  if( M >= 1 ){
    arma::mat D1( (arma::uword)n, (arma::uword)n, arma::fill::zeros );

    // D^1 from barycentric weights
    for( size_t i = 0; i < n; ++i ){
      double rowsum = 0.0;
      for( size_t j = 0; j < n; ++j ){
        if( i == j ) continue;
        const double v = _Lw[j] / ( _Lw[i] * (_x[i] - _x[j]) );
        D1( (arma::uword)i, (arma::uword)j ) = v;
        rowsum += v;
      }
      D1( (arma::uword)i, (arma::uword)i ) = -rowsum;
    }

    _D[1] = std::move(D1);

    // Higher orders by stable Welfert recursion
    for( int m = 2; m <= M; ++m ){
      arma::mat Dm( (arma::uword)n, (arma::uword)n, arma::fill::zeros );
      const arma::mat& Dm1 = _D[m-1];

      for( size_t i = 0; i < n; ++i ){
        double rowsum = 0.0;
        const double Dii = Dm1( (arma::uword)i, (arma::uword)i );

        for( size_t j = 0; j < n; ++j ){
          if( i == j ) continue;
          const double denom = _x[i] - _x[j];
          const double wij   = _Lw[j] / _Lw[i];
          const double v = ( (double)m / denom ) * ( wij * Dii - Dm1( (arma::uword)i, (arma::uword)j ) );
          Dm( (arma::uword)i, (arma::uword)j ) = v;
          rowsum += v;
        }
        Dm( (arma::uword)i, (arma::uword)i ) = -rowsum;
      }

      _D[m] = std::move(Dm);
    }
  }

  _Dmax = M;
  return true;
}

inline bool
BASE_OC::_ensure_lagrange_ops
( int const MAXORD )
const
{
  if( MAXORD <= _Dmax ) return true;
  return const_cast<BASE_OC*>(this)->set_lagrange_ops( MAXORD );
}

// ==========================
// Templated Lagrange utilities
// ==========================

template <typename U>
inline bool
BASE_OC::w_eval_lagrange
( std::vector<U>& w, double const& t, int const ORD, double const& TOL )
const
{
  if( _x.empty()
   || _Lw.empty()
   || _Lw.size() != _x.size()
   || ORD < 0 )
    return false;

  const size_t n = _x.size();
  w.assign( n, U(0) );

  // Trivial cases
  if( n == 1 ){
    if( ORD == 0 ) w[0] = U(1);
    return true;
  }

  // Degree is n-1, so derivatives of order >= n are identically zero
  if( ORD >= (int)n ){
    return true;
  }

  // Exact-node case
  for( size_t k = 0; k < n; ++k ){
    if( std::abs( t - _x[k] ) > TOL ) continue;

    if( ORD == 0 ){
      w[k] = U(1);
      return true;
    }

    if( !_ensure_lagrange_ops( ORD ) ) return false;
    for( size_t j = 0; j < n; ++j )
      w[j] = U( _D[ORD]( (arma::uword)k, (arma::uword)j ) );
    return true;
  }

  // General t (not a node): compute weight vector for p^(ORD)(t) = dot(w,y)
  // via stable ratio/Leibniz recursion on barycentric form.
  const int M = ORD;

  std::vector<long double> inv( n, 0.0L );
  std::vector<long double> alpha( n, 0.0L );
  for( size_t i = 0; i < n; ++i ){
    inv[i]   = 1.0L / ( (long double)t - (long double)_x[i] );
    alpha[i] = (long double)_Lw[i] * inv[i];
  }

  std::vector<long double> S( (size_t)M + 1, 0.0L );
  std::vector< std::vector<long double> > A( (size_t)M + 1, std::vector<long double>(n, 0.0L) );

  for( int m = 0; m <= M; ++m ){
    long double sm = 0.0L;
    for( size_t i = 0; i < n; ++i ){
      sm += alpha[i];
      A[(size_t)m][i] = alpha[i];
    }
    S[(size_t)m] = sm;
    // alpha^(m+1) from alpha^(m)
    for( size_t i = 0; i < n; ++i )
      alpha[i] = -(long double)(m+1) * alpha[i] * inv[i];
  }

  if( S[0] == 0.0L ) return false;

  std::vector< std::vector<long double> > W( (size_t)M + 1, std::vector<long double>(n, 0.0L) );
  for( size_t i = 0; i < n; ++i )
    W[0][i] = A[0][i] / S[0];

  for( int m = 1; m <= M; ++m ){
    std::vector<long double> tmp = A[(size_t)m];
    long double binom = 1.0L; // C(m,0)
    for( int k = 1; k <= m; ++k ){
      binom *= (long double)(m - k + 1) / (long double)k; // C(m,k)
      const long double coeff = binom * S[(size_t)k];
      const std::vector<long double>& Wmk = W[(size_t)(m-k)];
      for( size_t i = 0; i < n; ++i )
        tmp[i] -= coeff * Wmk[i];
    }
    for( size_t i = 0; i < n; ++i )
      W[(size_t)m][i] = tmp[i] / S[0];
  }

  for( size_t i = 0; i < n; ++i )
    w[i] = U( (double)W[(size_t)M][i] );

  return true;
}

template <typename U>
inline bool
BASE_OC::eval_lagrange
( U& yout, double const& t, int const ORD, std::vector<U> const& yin, double const& TOL, double const& RESCALE )
const
{
  if( _x.empty() || yin.size() != _x.size() || ORD < 0 )
    return false;

  // Degree is n-1, so derivatives of order >= n are identically zero
  if( ORD >= (int)_x.size() ){
    yout = U(0);
    return true;
  }

  if( !w_eval_lagrange( _weval, t, ORD, TOL ) )
    return false;

  yout = _weval[0] * yin[0];
  for( size_t i = 1; i < _weval.size(); ++i )
    yout += _weval[i] * yin[i];
  if( RESCALE != 1 )
    yout *= RESCALE;

  return true;
}

template <typename U>
inline bool
BASE_OC::eval_lagrange
( std::vector<U>& yout, int const ORD, std::vector<U> const& yin, double const& RESCALE )
const
{
  if( _x.empty() || yin.size() != _x.size() || ORD < 0 )
    return false;

  const size_t n = _x.size();
  yout.resize( n );

  // Trivial cases
  if( n == 1 ){
    if( yout.size() < 1 ) yout.resize(1);
    yout[0] = ( ORD == 0 ? yin[0] : 0. );
    return true;
  }

  if( ORD >= (int)n ){
    yout.assign( _x.size(), 0. );
    return true;
  }

  if( ORD == 0 ){
    yout = yin;
    return true;
  }

  if( _Lw.empty() || _Lw.size() != _x.size() || !_ensure_lagrange_ops( ORD ) )
    return false;

  // yout = D^(ORD) * yin
  for( size_t i = 0; i < n; ++i ){
    yout[i] = _D[ORD](i,0) * yin[0];
    for( size_t j = 1; j < n; ++j )
      yout[i] += _D[ORD](i,j) * yin[j];
      if( RESCALE != 1 )
        yout[i] *= RESCALE;
  }

  return true;
}

template <typename U>
inline bool
BASE_OC::eval_lagrange
( std::vector<U>& yout, double const& t, int const ORD, std::vector<U> const& yin,
  size_t const STRIDE, size_t const OUTER, double const& TOL, double const& RESCALE )
const
{
  if( _x.empty() || ORD < 0 || !STRIDE || !OUTER )
    return false;

  const size_t n    = _x.size();
  const size_t nin  = n * STRIDE * OUTER;
  const size_t nout = STRIDE * OUTER;

  if( yin.size() != nin )
    return false;

  yout.resize( nout );

  // Degree is n-1, so derivatives of order >= n are identically zero
  if( ORD >= (int)n ){
    yout.assign( nout, U(0) );
    return true;
  }

  if( !w_eval_lagrange( _weval, t, ORD, TOL ) || _weval.size() != n )
    return false;

  // Apply yout = w * yin on each 1D slice
  for( size_t o = 0; o < OUTER; ++o ){
    const size_t block_in  = o * n * STRIDE;
    const size_t block_out = o * STRIDE;

    for( size_t s = 0; s < STRIDE; ++s ){
      const size_t base = block_in + s;
      const size_t iout = block_out + s;

      yout[iout] = _weval[0] * yin[base];
      for( size_t i = 1; i < n; ++i )
        yout[iout] += _weval[i] * yin[base + i * STRIDE];
      if( RESCALE != 1 )
        yout[iout] *= RESCALE;
    }
  }

  return true;
}

template <typename U>
inline bool
BASE_OC::_eval_lagrange
( U* yout, double const& t, U const* yin, size_t const STRIDE, size_t const OUTER,
  bool const sub, double const& TOL )
const
{
  if( _x.empty() || !STRIDE || !OUTER || !yin || !yout )
    return false;

  const size_t n    = _x.size();
  //const size_t nin  = n * STRIDE * OUTER;
  //const size_t nout = STRIDE * OUTER;

  if( !w_eval_lagrange( _weval, t, 0, TOL ) || _weval.size() != n )
    return false;

  // Apply yout = w * yin on each 1D slice
  for( size_t o = 0; o < OUTER; ++o ){
    const size_t block_in  = o * n * STRIDE;
    const size_t block_out = o * STRIDE;

    for( size_t s = 0; s < STRIDE; ++s ){
      const size_t base = block_in + s;
      const size_t iout = block_out + s;

      if( sub ){  // subtract from existing content of yout
        for( size_t i = 0; i < n; ++i )
          yout[iout] -= _weval[i] * yin[base + i * STRIDE];     
      }

      else{       // overwrite existing content of yout
        yout[iout] = _weval[0] * yin[base];
        for( size_t i = 1; i < n; ++i )
          yout[iout] += _weval[i] * yin[base + i * STRIDE];
      }
    }
  }

  return true;
}

template <typename U>
inline bool
BASE_OC::eval_lagrange
( std::vector<U>& yout, int const ORD, std::vector<U> const& yin,
  size_t const STRIDE, size_t const OUTER, double const& RESCALE )
const
{
  if( _x.empty() || ORD < 0 || !STRIDE || !OUTER )
    return false;

  if( &yout == &yin )
    return false;

  const size_t n    = _x.size();
  const size_t ntot = n * STRIDE * OUTER;

  if( yin.size() != ntot )
    return false;

  yout.resize( ntot );

  // Trivial cases
  if( n == 1 ){
    if( ORD == 0 ) yout = yin;
    else           yout.assign( ntot, U(0) );
    return true;
  }

  if( ORD >= (int)n ){
    yout.assign( ntot, U(0) );
    return true;
  }

  if( ORD == 0 ){
    yout = yin;
    return true;
  }

  if( _Lw.empty() || _Lw.size() != n || !_ensure_lagrange_ops( ORD ) )
    return false;

  // Apply yout = D^(ORD) * yin on each 1D slice
  for( size_t o = 0; o < OUTER; ++o ){
    const size_t block = o * n * STRIDE;

    for( size_t s = 0; s < STRIDE; ++s ){
      const size_t base = block + s;

      for( size_t i = 0; i < n; ++i ){
        const size_t iout = base + i * STRIDE;
        
        yout[iout] = _D[ORD](i,0) * yin[base];
        for( size_t j = 1; j < n; ++j )
          yout[iout] += _D[ORD](i,j) * yin[base + j * STRIDE];
        if( RESCALE != 1 )
          yout[iout] *= RESCALE;
      }
    }
  }

  return true;
}

template <typename U>
inline bool
BASE_OC::eval_lagrange_grad
( U& val, std::vector<U>& w, double const& t, int const ORD, std::vector<U> const& y, double const& TOL )
const
{
  if( y.size() != _x.size() ) return false;
  if( !w_eval_lagrange( w, t, ORD, TOL ) ) return false;

  U acc = U(0);
  for( size_t i = 0; i < w.size(); ++i )
    acc += w[i] * y[i];
  val = acc;
  return true;
}

template <typename U>
inline bool
BASE_OC::quad_lagrange
( U& yout, std::vector<U> const& yin, double const& RESCALE )
const
{
  if( _x.empty()
   || _w.empty()
   || _x.size()  != _w.size()
   || yin.size() != _x.size() )
    return false;

  yout = _w[0] * yin[0];
  for( size_t i = 1; i < _x.size(); ++i )
    yout += _w[i] * yin[i];
  if( RESCALE != 1 )
    yout *= RESCALE;
  return true;
}

template <typename U>
inline bool
BASE_OC::quad_lagrange
( std::vector<U>& yout, std::vector<U> const& yin,
  size_t const STRIDE, size_t const OUTER, double const& RESCALE )
const
{
  if( _x.empty()
   || _w.empty()
   || _x.size() != _w.size()
   || !STRIDE
   || !OUTER )
    return false;

  const size_t n = _x.size();
  const size_t nin  = n * STRIDE * OUTER;
  const size_t nout = STRIDE * OUTER;

  if( yin.size() != nin )
    return false;

  yout.resize( nout );

  for( size_t o = 0; o < OUTER; ++o ){
    const size_t block_in  = o * n * STRIDE;
    const size_t block_out = o * STRIDE;

    for( size_t s = 0; s < STRIDE; ++s ){
      const size_t base = block_in + s;
      const size_t iout = block_out + s;

      yout[iout] = _w[0] * yin[base];
      for( size_t i = 1; i < n; ++i )
        yout[iout] += _w[i] * yin[base + i * STRIDE];
      if( RESCALE != 1 )
        yout[iout] *= RESCALE;
    }
  }

  return true;
}

template <typename U>
inline bool
BASE_OC::w_quad_lagrange
( std::vector<U>& wq )
const
{
  if( _x.empty() || _w.empty() || _x.size() != _w.size() )
    return false;
  const size_t n = _x.size();
  wq.resize(n);
  for( size_t i = 0; i < n; ++i )
    wq[i] = U(_w[i]);
  return true;
}


template <typename U>
inline bool
BASE_OC::quad_lagrange_grad
( U& val, std::vector<U>& wq, std::vector<U> const& y )
const
{
  if( y.size() != _x.size() ) return false;
  if( !w_quad_lagrange( wq ) ) return false;

  U acc = U(0);
  for( size_t i = 0; i < wq.size(); ++i )
    acc += wq[i] * y[i];
  val = acc;
  return true;
}

// New overload: Lagrange evaluator of ORD-th derivative at a single point t,
// applied to strided/blocked data yin, with optional RESCALE and sub flag.
// This is the counterpart to eval_lagrange(yout,t,ORD,yin,STRIDE,OUTER,...) but
// operates on raw pointers and supports in-place subtraction, making it compatible
// with the existing _eval_lagrange interface used in the continuity residuals.
template <typename U>
inline bool
BASE_OC::_eval_lagrange
( U* yout, double const& t, int const ORD, U const* yin, size_t const STRIDE,
  size_t const OUTER, bool const sub, double const& TOL, double const& RESCALE )
const
{
  if( _x.empty() || !STRIDE || !OUTER || !yin || !yout || ORD < 0 )
    return false;

  const size_t n = _x.size();

  // Derivatives of order >= n are identically zero
  if( ORD >= (int)n ){
    if( !sub )
      for( size_t o = 0; o < OUTER; ++o )
        for( size_t s = 0; s < STRIDE; ++s )
          yout[ o * STRIDE + s ] = U(0);
    // if sub: subtracting zero changes nothing — leave yout unchanged
    return true;
  }

  // Compute the weight vector w such that p^(ORD)(t) = dot(w, y)
  if( !w_eval_lagrange( _weval, t, ORD, TOL ) || _weval.size() != n )
    return false;

  // Apply w to each 1D slice of the strided/blocked data
  for( size_t o = 0; o < OUTER; ++o ){
    const size_t block_in  = o * n * STRIDE;
    const size_t block_out = o * STRIDE;

    for( size_t s = 0; s < STRIDE; ++s ){
      const size_t base = block_in + s;
      const size_t iout = block_out + s;

      // Accumulate dot product
      U val = _weval[0] * yin[base];
      for( size_t i = 1; i < n; ++i )
        val += _weval[i] * yin[base + i * STRIDE];
      if( RESCALE != 1. )
        val *= RESCALE;

      if( sub ) yout[iout] -= val;
      else      yout[iout]  = val;
    }
  }

  return true;
}

////////////////////////////////////////////////////////////////////////

template <typename T>
std::ostream& operator<<
( std::ostream& out, OCVar<T> const& var )
{
  out << "deps: {";
  for( auto const& [v,d] : var._dom )
    out << " " << *v;
  out << " }, "; 
  out << "coef: [" << std::scientific << std::setprecision(5);
  for( auto const& x : var._coef )
    out << std::setw(13) << x;
  out << " ]";
  return out;
}


template <class T>
inline void
OCVar<T>::_set_elem
( std::map<FFVar,size_t,lt_FFVar> const& ndx_el )
{
  _elem.clear();
  for( auto const& [pvar,pdom] : _dom ){
    if( !pvar || !pdom ) continue;
    auto it = ndx_el.find( *pvar );
    if( it != ndx_el.end() ) _elem[pvar] = it->second;
  }
}

template <class T>
inline double
OCVar<T>::_elem_width
( FFVar const* pvar )
const
{
  auto itd = _dom.find( pvar );
  if( itd == _dom.end() || !itd->second )
    throw OCBase::Exceptions( OCBase::Exceptions::INTERNAL );

  auto ite = _elem.find( pvar );
  if( ite != _elem.end() )
    return itd->second->elem_width( ite->second );

  return itd->second->w_elem;
}

template <class T>
template <class U>
inline OCVar<T>&
OCVar<T>::_update
( U const* coef )
{
  // Reuse the coefficient storage.  For OCVar<FADType<double>> this avoids
  // reallocating every local tensor on each derivative evaluation.  The FADType
  // objects themselves are rebuilt and seeded only by OCEnv::setup() (or by a
  // later cache rebuild), so in-place assignment is safe as long as the cached
  // color count has not changed.
  assert( coef && _coef.size() );
  size_t const ncoef = _coef.size();
  for( size_t i = 0; i < ncoef; ++i )
    _coef[i] = coef[i];

  return *this;
}

template <class T>
inline OCVar<T>&
OCVar<T>::_set
( OCBase const& env, std::set< FFVar, lt_FFVar > const& dep, T const* coef )
{
  // Set environment
  _env = &env;
       
  // Set domain
  _dom.clear();
  size_t ncoef = 1;
  for( auto const& var : dep ){
    auto itvar = _env->_mDom.find( var );
    assert( itvar != _env->_mDom.end() );
    _dom[&itvar->first] = &itvar->second; 
    ncoef *= itvar->second.n_node;
  }
      
  // Set coefficients.  Re-create the storage rather than resizing/assigning
  // in place; this avoids reusing cached FADType objects with stale
  // derivative-buffer lengths.
  std::vector<T>().swap( _coef );
  if( coef )
    _coef.assign( coef, coef+ncoef );
  else
    _coef.resize( ncoef );

  return *this;
}

template <class T>
inline void
OCVar<T>::_set
( OCBase const* env, OCVar<T>::t_Dom const& dep, std::vector<T> && coef )
{
  // Set environment
  _env = env;
       
  // Set domain
  _dom = dep;
  _elem.clear();
      
  // Set coefficients
  _coef = std::move( coef );
}

template <class T>
inline bool
OCVar<T>::lift
( std::set< FFVar, lt_FFVar > const& dep, bool const check )
{
  // Check new dependents form a superset
  for( auto const& [pvar,pdom] : _dom )
    if( dep.find( *pvar ) == dep.cend() ) return false;

  t_Dom dom_out = _dom;
  for( auto const& var : dep ){
    if( dom_out.find( &var ) != dom_out.end() ) continue;
    auto itvar = _env->_mDom.find( var );
    if( itvar == _env->_mDom.end() ) return false;
    dom_out[&itvar->first] = &itvar->second;
  }
      
  std::vector<T> coef_out;
  if( _env->lift_colloc( coef_out, _coef, dom_out, _dom, check ) ){
    _coef.swap( coef_out );
    _dom.swap( dom_out );
    return true;
  }
  return false;
}

template <typename T>
inline
OCVar<T> operator+
( OCVar<T> const& var )
{
  return var;
}

template <typename T>
inline
OCVar<T> && operator+
( OCVar<T> && var )
{
  return( std::move(var) );  
}

template <typename T>
inline
OCVar<T>& OCVar<T>::operator+=
( double const& cst )
{
  if( cst == 0. )
    return *this;

  for( auto& val : _coef )
    val += cst;
  return *this;
}

template <typename T>
inline
OCVar<T>& OCVar<T>::operator+=
( OCVar<T> const& var )
{
  if( !_env && !var._env ){
    if( _coef.size() != 1 || var._coef.size() != 1 )
      throw OCBase::Exceptions( OCBase::Exceptions::INTERNAL );
    _coef[0] += var._coef[0];
    return *this;
  }

  if( !var._env || var._dom.empty() ){
    if( var._coef.size() != 1 )
      throw OCBase::Exceptions( OCBase::Exceptions::INTERNAL );
    *this += var._coef[0];
    return *this;
  }

  if( !_env || _dom.empty() ){
    if( _coef.size() != 1 )
      throw OCBase::Exceptions( OCBase::Exceptions::INTERNAL );
    T cst = _coef[0];
    *this = var;
    *this += cst;
    return *this;
  }

  if( _env != var._env )
    throw OCBase::Exceptions( OCBase::Exceptions::ENV );
  // Fast path: identical domains and grids (the common case) -> coefficient-wise sum,
  // no domain copy and no Lagrange lift (lift_colloc has no identity short-circuit).
  if( _dom == var._dom ){
    if( _coef.size() != var._coef.size() )
      throw OCBase::Exceptions( OCBase::Exceptions::INTERNAL );
    for( size_t i=0; i<_coef.size(); ++i )
      _coef[i] += var._coef[i];
    _elem.insert( var._elem.cbegin(), var._elem.cend() );
    return *this;
  }

  bool check = true;
  // Union the operand domains; for a dimension shared by both operands keep the
  // FINER grid (larger n_node).  std::map::insert would keep this->_dom's grid for a
  // shared key, making the result depend on operand order and collapsing a coarse
  // coefficient input onto its own grid; taking the max lifts the coarse operand UP
  // to the finer (e.g. state) grid instead.  No-op when the grids already match.
  auto domlift = _dom;
  for( auto const& kv : var._dom ){
    auto it = domlift.find( kv.first );
    if( it == domlift.end() )
      domlift.insert( kv );
    else if( kv.second && it->second && kv.second->n_node > it->second->n_node )
      it->second = kv.second;
  }
  std::vector<T> coeflift;
  if( _env->lift_colloc( coeflift, _coef, domlift, _dom, check ) ){
    _coef.swap( coeflift );
    _dom.swap( domlift );
    _elem.insert( var._elem.cbegin(), var._elem.cend() );
  }
  if( _env->lift_colloc( coeflift, var._coef, _dom, var._dom, check ) ){
    if( _coef.size() != coeflift.size() )
      throw OCBase::Exceptions( OCBase::Exceptions::INTERNAL );
    for( size_t i=0; i<_coef.size(); ++i )
      _coef[i] += coeflift[i];
  }
  else{
    if( _coef.size() != var._coef.size() )
      throw OCBase::Exceptions( OCBase::Exceptions::INTERNAL );
    for( size_t i=0; i<_coef.size(); ++i )
      _coef[i] += var._coef[i];  
  }
  _elem.insert( var._elem.cbegin(), var._elem.cend() );
  return *this;
}

template <typename T>
inline
OCVar<T> operator+
( OCVar<T> const& var1, OCVar<T> const& var2 )
{
  if( !var1._env && !var2._env ){
    if( var1._coef.size() != 1 || var2._coef.size() != 1 )
      throw OCBase::Exceptions( OCBase::Exceptions::INTERNAL );
    return( var1._coef[0] + var2._coef[0] );
  }
  
  if( var1._env ){
    OCVar<T> var3( var1 );
    var3 += var2;
    return var3;
  }
  
  OCVar<T> var3( var2 );
  var3 += var1;
  return var3;
}

template <typename T>
inline
OCVar<T> && operator+
( OCVar<T> const& var1, OCVar<T> && var2 )
{
  if( !var1._env && !var2._env ){
    if( var1._coef.size() != 1 || var2._coef.size() != 1 )
      throw OCBase::Exceptions( OCBase::Exceptions::INTERNAL );
    var2._coef[0] += var1._coef[0]; 
    return std::move( var2 );
  }

  var2 += var1;
  return std::move( var2 );
}

template <typename T>
inline
OCVar<T> && operator+
( OCVar<T> && var1, OCVar<T> const& var2 )
{
  if( !var1._env && !var2._env ){
    if( var1._coef.size() != 1 || var2._coef.size() != 1 )
      throw OCBase::Exceptions( OCBase::Exceptions::INTERNAL );
    var1._coef[0] += var2._coef[0]; 
    return std::move( var1 );
  }
   
  var1 += var2;
  return std::move( var1 );
}

template <typename T>
inline
OCVar<T> && operator+
( OCVar<T> && var1, OCVar<T> && var2 )
{
  if( !var1._env && !var2._env ){
    if( var1._coef.size() != 1 || var2._coef.size() != 1 )
      throw OCBase::Exceptions( OCBase::Exceptions::INTERNAL );
    var1._coef[0] += var2._coef[0]; 
    return std::move( var1 );
  }

  if( var1._env ){
    var1 += var2;
    return std::move( var1 );
  }

  var2 += var1;
  return std::move( var2 );
}

template <typename T>
inline
OCVar<T> operator+
( OCVar<T> const& var1, double const& cst2 )
{
  if( !var1._env ){
    if( var1._coef.size() != 1 )
      throw OCBase::Exceptions( OCBase::Exceptions::INTERNAL );
    return( var1._coef[0] + cst2 );
  }
  
  OCVar<T> var3( var1 );
  var3 += cst2;
  return var3;
}

template <typename T>
inline
OCVar<T> && operator+
( OCVar<T> && var1, double const& cst2 )
{
  var1 += cst2;
  return std::move( var1 );
}

template <typename T>
inline
OCVar<T> operator+
( double const& cst1, OCVar<T> const& var2 )
{
  return var2 + cst1;
}

template <typename T>
inline
OCVar<T> && operator+
( double const& cst1, OCVar<T> && var2 )
{
  var2 += cst1;
  return std::move( var2 );
}

template <typename T> 
inline
OCVar<T> operator-
( OCVar<T> const& var )
{
  OCVar<T> var2( var );
  return operator-( std::move( var2 ) );
}


template <typename T> 
inline
OCVar<T> && operator-
( OCVar<T> && var )
{
  for( auto& val : var._coef )
    val *= -1;
  return std::move( var );
}

template <typename T>
inline
OCVar<T>& OCVar<T>::operator-=
( double const& cst )
{
  *this += -cst;
  return *this;
}

template <typename T>
inline
OCVar<T>& OCVar<T>::operator-=
( OCVar<T> const& var )
{
  OCVar<T> var2( var );
  return operator+=( operator-( std::move( var2 ) ) );
}

template <typename T>
inline
OCVar<T>& OCVar<T>::operator-=
( OCVar<T> && var )
{
  return operator+=( operator-( std::move( var ) ) );
}

template <typename T>
inline
OCVar<T> operator-
( OCVar<T> const& var1, OCVar<T> const& var2 )
{
  if( !var1._env && !var2._env ){
    if( var1._coef.size() != 1 || var2._coef.size() != 1 )
      throw OCBase::Exceptions( OCBase::Exceptions::INTERNAL );
    return( var1._coef[0] - var2._coef[0] );
  }

  if( var1._env ){
    OCVar<T> var3( var1 );
    var3 -= var2;    
    return var3;
  }

  OCVar<T> var3( -var2 );
  var3 += var1;
  return var3;
}

template <typename T>
inline
OCVar<T> && operator-
( OCVar<T> const& var1, OCVar<T> && var2 )
{
  if( !var1._env && !var2._env ){
    if( var1._coef.size() != 1 || var2._coef.size() != 1 )
      throw OCBase::Exceptions( OCBase::Exceptions::INTERNAL );
    var2._coef[0] = var1._coef[0] - var2._coef[0]; 
    return std::move( var2 );
  }
  
  var2 -= var1;
  return operator-( std::move( var2 ) );
}

template <typename T>
inline
OCVar<T> && operator-
( OCVar<T> && var1, OCVar<T> const& var2 )
{
  if( !var1._env && !var2._env ){
    if( var1._coef.size() != 1 || var2._coef.size() != 1 )
      throw OCBase::Exceptions( OCBase::Exceptions::INTERNAL );
    var1._coef[0] -= var2._coef[0]; 
    return std::move( var1 );
  }
   
  var1 -= var2;
  return std::move( var1 );
}

template <typename T>
inline
OCVar<T> && operator-
( OCVar<T> && var1, OCVar<T> && var2 )
{
  if( !var1._env && !var2._env ){
    if( var1._coef.size() != 1 || var2._coef.size() != 1 )
      throw OCBase::Exceptions( OCBase::Exceptions::INTERNAL );
    var1._coef[0] -= var2._coef[0]; 
    return std::move( var1 );
  }
  
  var1 -= std::move( var2 );
  return std::move( var1 );
}

template <typename T>
inline
OCVar<T> operator-
( OCVar<T> const& var1, double const& cst2 )
{
  if( !var1._env ){
    if( var1._coef.size() != 1 )
      throw OCBase::Exceptions( OCBase::Exceptions::INTERNAL );
    return( var1._coef[0] - cst2 );
  }
      
  OCVar<T> var3( var1 );
  var3 += -cst2;
  return var3;
}

template <typename T>
inline
OCVar<T> && operator-
( OCVar<T> && var1, double const& cst2 )
{
  var1 -= cst2;
  return std::move( var1 );
}

template <typename T>
inline
OCVar<T> operator-
( double const& cst1, OCVar<T> const& var2 )
{
  if( !var2._env ){
    if( var2._coef.size() != 1 )
      throw OCBase::Exceptions( OCBase::Exceptions::INTERNAL );
    return( cst1 - var2._coef[0] );
  }
  
  OCVar<T> var3( -var2 );
  var3 += cst1;
  return var3;
}

template <typename T>
inline
OCVar<T> && operator-
( double const& cst1, OCVar<T> && var2 )
{
  var2 -= cst1;
  return operator-( std::move( var2 ) );
}

template <typename T>
inline
OCVar<T>& OCVar<T>::operator*=
( double const& cst )
{
  if( cst == 1. )
    return *this;

  for( auto& val : _coef )
    val *= cst;
  return *this;
}

template <typename T>
inline
OCVar<T>& OCVar<T>::operator*=
( OCVar<T> const& var )
{
  if( !_env && !var._env ){
    if( _coef.size() != 1 || var._coef.size() != 1 )
      throw OCBase::Exceptions( OCBase::Exceptions::INTERNAL );
    _coef[0] *= var._coef[0];
    return *this;
  }

  if( !var._env || var._dom.empty() ){
    if( var._coef.size() != 1 )
      throw OCBase::Exceptions( OCBase::Exceptions::INTERNAL );
    *this *= var._coef[0];
    return *this;
  }

  if( !_env || _dom.empty() ){
    if( _coef.size() != 1 )
      throw OCBase::Exceptions( OCBase::Exceptions::INTERNAL );
    T cst = _coef[0];
    *this = var;
    *this *= cst;
    return *this;
  }

  if( _env != var._env )
    throw OCBase::Exceptions( OCBase::Exceptions::ENV );
  // Fast path: identical domains and grids (the common case) -> coefficient-wise
  // product, no domain copy and no Lagrange lift.  (lift_colloc has no identity
  // short-circuit, so without this it runs a full interpolation pass per operand even
  // when both operands already share the grid.)
  if( _dom == var._dom ){
    if( _coef.size() != var._coef.size() )
      throw OCBase::Exceptions( OCBase::Exceptions::INTERNAL );
    for( size_t i=0; i<_coef.size(); ++i )
      _coef[i] *= var._coef[i];
    return *this;
  }

  bool check = true;
  // Union the operand domains; for a dimension shared by both operands keep the
  // FINER grid (larger n_node) so the product is independent of operand order and a
  // coarse coefficient input (e.g. n_node=1) lifts UP to the finer (state) grid
  // rather than collapsing the product onto the coarse grid.  No-op when grids match.
  auto domlift = _dom;
  for( auto const& kv : var._dom ){
    auto it = domlift.find( kv.first );
    if( it == domlift.end() )
      domlift.insert( kv );
    else if( kv.second && it->second && kv.second->n_node > it->second->n_node )
      it->second = kv.second;
  }
  std::vector<T> coeflift;
  if( _env->lift_colloc( coeflift, _coef, domlift, _dom, check ) ){
    _coef.swap( coeflift );
    _dom.swap( domlift );
  }
  if( _env->lift_colloc( coeflift, var._coef, _dom, var._dom, check ) ){
    if( _coef.size() != coeflift.size() )
      throw OCBase::Exceptions( OCBase::Exceptions::INTERNAL );
    for( size_t i=0; i<_coef.size(); ++i )
      _coef[i] *= coeflift[i];
  }
  else{
    if( _coef.size() != var._coef.size() )
      throw OCBase::Exceptions( OCBase::Exceptions::INTERNAL );
    for( size_t i=0; i<_coef.size(); ++i )
      _coef[i] *= var._coef[i];  
  }
  return *this;
}

template <typename T>
inline
OCVar<T> operator*
( OCVar<T> const& var1, OCVar<T> const& var2 )
{
  if( !var1._env && !var2._env ){
    if( var1._coef.size() != 1 || var2._coef.size() != 1 )
      throw OCBase::Exceptions( OCBase::Exceptions::INTERNAL );
    return( var1._coef[0] * var2._coef[0] );
  }
  
  if( var1._env ){
    OCVar<T> var3( var1 );
    var3 *= var2;
    return var3;
  }
  
  OCVar<T> var3( var2 );
  var3 *= var1;
  return var3;
}

template <typename T>
inline
OCVar<T> && operator*
( OCVar<T> const& var1, OCVar<T> && var2 )
{
  if( !var1._env && !var2._env ){
    if( var1._coef.size() != 1 || var2._coef.size() != 1 )
      throw OCBase::Exceptions( OCBase::Exceptions::INTERNAL );
    var2._coef[0] *= var1._coef[0]; 
    return std::move( var2 );
  }

  var2 *= var1;
  return std::move( var2 );
}

template <typename T>
inline
OCVar<T> && operator*
( OCVar<T> && var1, OCVar<T> const& var2 )
{
  if( !var1._env && !var2._env ){
    if( var1._coef.size() != 1 || var2._coef.size() != 1 )
      throw OCBase::Exceptions( OCBase::Exceptions::INTERNAL );
    var1._coef[0] *= var2._coef[0]; 
    return std::move( var1 );
  }
   
  var1 *= var2;
  return std::move( var1 );
}

template <typename T>
inline
OCVar<T> && operator*
( OCVar<T> && var1, OCVar<T> && var2 )
{
  if( !var1._env && !var2._env ){
    if( var1._coef.size() != 1 || var2._coef.size() != 1 )
      throw OCBase::Exceptions( OCBase::Exceptions::INTERNAL );
    var1._coef[0] *= var2._coef[0]; 
    return std::move( var1 );
  }

  if( var1._env ){
    var1 *= var2;
    return std::move( var1 );
  }

  var2 *= var1;
  return std::move( var2 );
}

template <typename T>
inline
OCVar<T> operator*
( OCVar<T> const& var1, double const& cst2 )
{
  if( !var1._env ){
    if( var1._coef.size() != 1 )
      throw OCBase::Exceptions( OCBase::Exceptions::INTERNAL );
    return( var1._coef[0] * cst2 );
  }
  
  OCVar<T> var3( var1 );
  var3 *= cst2;
  return var3;
}

template <typename T>
inline
OCVar<T> && operator*
( OCVar<T> && var1, double const& cst2 )
{
  var1 *= cst2;
  return std::move( var1 );
}

template <typename T>
inline
OCVar<T> operator*
( double const& cst1, OCVar<T> const& var2 )
{
  return var2 * cst1;
}

template <typename T>
inline
OCVar<T> && operator*
( double const& cst1, OCVar<T> && var2 )
{
  var2 *= cst1;
  return std::move( var2 );
}

template <typename T>
inline
OCVar<T>& OCVar<T>::operator/=
( double const& cst )
{
  if( cst == 1. )
    return *this;

  return( operator*=( 1./cst ) );
}

template <typename T>
inline
OCVar<T>& OCVar<T>::operator/=
( OCVar<T> const& var )
{
  return( operator*=( inv( var ) ) );
}

template <typename T>
inline
OCVar<T>& OCVar<T>::operator/=
( OCVar<T> && var )
{
  return( operator*=( inv( std::move(var) ) ) );
}

template <typename T>
inline
OCVar<T> operator/
( OCVar<T> const& var1, OCVar<T> const& var2 )
{
  OCVar<T> var3( var1 );
  var3 /= var2;
  return var3;
}

template <typename T>
inline
OCVar<T> && operator/
( OCVar<T> const& var1, OCVar<T> && var2 )
{
  OCVar<T> var3( var1 );
  var3 /= std::move( var2 );
  return var3;
}

template <typename T>
inline
OCVar<T> && operator/
( OCVar<T> && var1, OCVar<T> const& var2 )
{
  var1 /= var2;
  return std::move( var1 );
}

template <typename T>
inline
OCVar<T> && operator/
( OCVar<T> && var1, OCVar<T> && var2 )
{
  var1 /= std::move( var2 );
  return std::move( var1 );
}

template <typename T>
inline
OCVar<T> operator/
( OCVar<T> const& var1, double const& cst2 )
{
  OCVar<T> var3( var1 );
  var3 /= cst2;
  return var3;
}

template <typename T>
inline
OCVar<T> && operator/
( OCVar<T> && var1, double const& cst2 )
{
  var1 /= cst2;
  return std::move( var1 );
}

template <typename T>
inline
OCVar<T> operator/
( double const& cst1, OCVar<T> const& var2 )
{
  OCVar<T> var3( inv( var2 ) );
  var3 *= cst1;
  return var3;
}

template <typename T>
inline
OCVar<T> && operator/
( double const& cst1, OCVar<T> && var2 )
{
  inv( std::move( var2 ) );
  var2 *= cst1;
  return std::move( var2 );
}

template <typename T>
inline
OCVar<T> inv
( OCVar<T> const& var )
{
  OCVar<T> var2( var );
  for( size_t i=0; i<var._coef.size(); ++i )
    var2._coef[i] = 1. / var._coef[i];  
  return var2;  
}

template <typename T>
inline
OCVar<T> && inv
( OCVar<T> && var )
{
  for( size_t i=0; i<var._coef.size(); ++i )
    var._coef[i] = 1. / var._coef[i];  
  return std::move( var );
}

template <typename T>
inline
OCVar<T> max
( OCVar<T> const& var1, double const& cst2 )
{
  OCVar<T> var3( var1 );
  return max( std::move( var3 ), cst2 );
}

template <typename T>
inline
OCVar<T> && max
( OCVar<T> && var1, double const& cst2 )
{
  for( size_t i=0; i<var1._coef.size(); ++i )
    var1._coef[i] = Op<T>::max( var1._coef[i], cst2 );
  return std::move( var1 );
}

template <typename T>
inline
OCVar<T> max
( OCVar<T> const& var1, OCVar<T> const& var2 )
{
  OCVar<T> var3( var1 - var2 );
  max( std::move( var3 ), 0. );
  var3 += var2;
  return var3;
}

template <typename T>
inline
OCVar<T> min
( OCVar<T> const& var1, double const& cst2 )
{
  OCVar<T> var3( var1 );
  return min( std::move( var3 ), cst2 );
}

template <typename T>
inline
OCVar<T> && min
( OCVar<T> && var1, double const& cst2 )
{
  for( size_t i=0; i<var1._coef.size(); ++i )
    var1._coef[i] = Op<T>::min( var1._coef[i], cst2 );
  return std::move( var1 );
}

template <typename T>
inline
OCVar<T> min
( OCVar<T> const& var1, OCVar<T> const& var2 )
{
  OCVar<T> var3( var1 - var2 );
  min( std::move( var3 ), 0. );
  var3 += var2;
  return var3;
}

template <typename T>
inline
OCVar<T> prod
( unsigned int const nvars, OCVar<T> const* pvars )
{
  if( !nvars ) return T(1.);
  OCVar<T> tmp = pvars[0];
  for( unsigned int i=1; i<nvars; ++i ) tmp *= pvars[i];
  return tmp;
}

template <typename T>
inline
OCVar<T> cheb
( OCVar<T> const& var, unsigned int const& n )
{
  OCVar<T> var2( var );
  for( size_t i=0; i<var._coef.size(); ++i )
    var2._coef[i] = Op<T>::cheb( var._coef[i], n );  
  return var2;  
}

template <typename T>
inline
OCVar<T> && cheb
( OCVar<T> && var, unsigned int const& n )
{
  for( size_t i=0; i<var._coef.size(); ++i )
    var._coef[i] = Op<T>::cheb( var._coef[i], n );  
  return std::move( var );
}

template <typename T>
inline
OCVar<T> pow
( OCVar<T> const& var, int const& n )
{
  OCVar<T> var2( var );
  for( size_t i=0; i<var._coef.size(); ++i )
    var2._coef[i] = Op<T>::pow( var._coef[i], n );  
  return var2;  
}

template <typename T>
inline
OCVar<T> && pow
( OCVar<T> && var, int const& n )
{
  for( size_t i=0; i<var._coef.size(); ++i )
    var._coef[i] = Op<T>::pow( var._coef[i], n );  
  return std::move( var );
}

template <typename T, typename E>
inline
OCVar<T> pow
( OCVar<T> const& var, E const& e )
{
  OCVar<T> var2( var );
  for( size_t i=0; i<var._coef.size(); ++i )
    var2._coef[i] = Op<T>::pow( var._coef[i], e );  
  return var2;  
}

template <typename T, typename E>
inline
OCVar<T> && pow
( OCVar<T> && var, E const& e )
{
  for( size_t i=0; i<var._coef.size(); ++i )
    var._coef[i] = Op<T>::pow( var._coef[i], e );  
  return std::move( var );
}

template <typename T>
inline
OCVar<T> exp
( OCVar<T> const& var )
{
  OCVar<T> var2( var );
  for( size_t i=0; i<var._coef.size(); ++i )
    var2._coef[i] = Op<T>::exp( var._coef[i] );  
  return var2;  
}

template <typename T>
inline
OCVar<T> && exp
( OCVar<T> && var )
{
  for( size_t i=0; i<var._coef.size(); ++i )
    var._coef[i] = Op<T>::exp( var._coef[i] );  
  return std::move( var );
}

template <typename T>
inline
OCVar<T> log
( OCVar<T> const& var )
{
  OCVar<T> var2( var );
  for( size_t i=0; i<var._coef.size(); ++i )
    var2._coef[i] = Op<T>::log( var._coef[i] );  
  return var2;  
}

template <typename T>
inline
OCVar<T> && log
( OCVar<T> && var )
{
  for( size_t i=0; i<var._coef.size(); ++i )
    var._coef[i] = Op<T>::log( var._coef[i] );  
  return std::move( var );
}

template <typename T>
inline
OCVar<T> sqr
( OCVar<T> const& var )
{
  OCVar<T> var2( var );
  for( size_t i=0; i<var._coef.size(); ++i )
    var2._coef[i] = Op<T>::sqr( var._coef[i] );  
  return var2;  
}

template <typename T>
inline
OCVar<T> && sqr
( OCVar<T> && var )
{
  for( size_t i=0; i<var._coef.size(); ++i )
    var._coef[i] = Op<T>::sqr( var._coef[i] );  
  return std::move( var );
}

template <typename T>
inline
OCVar<T> sqrt
( OCVar<T> const& var )
{
  OCVar<T> var2( var );
  for( size_t i=0; i<var._coef.size(); ++i )
    var2._coef[i] = Op<T>::sqrt( var._coef[i] );  
  return var2;  
}

template <typename T>
inline
OCVar<T> && sqrt
( OCVar<T> && var )
{
  for( size_t i=0; i<var._coef.size(); ++i )
    var._coef[i] = Op<T>::sqrt( var._coef[i] );  
  return std::move( var );
}

template <typename T>
inline
OCVar<T> sin
( OCVar<T> const& var )
{
  OCVar<T> var2( var );
  for( size_t i=0; i<var._coef.size(); ++i )
    var2._coef[i] = Op<T>::sin( var._coef[i] );  
  return var2;  
}

template <typename T>
inline
OCVar<T> && sin
( OCVar<T> && var )
{
  for( size_t i=0; i<var._coef.size(); ++i )
    var._coef[i] = Op<T>::sin( var._coef[i] );  
  return std::move( var );
}

template <typename T>
inline
OCVar<T> cos
( OCVar<T> const& var )
{
  OCVar<T> var2( var );
  for( size_t i=0; i<var._coef.size(); ++i )
    var2._coef[i] = Op<T>::cos( var._coef[i] );  
  return var2;  
}

template <typename T>
inline
OCVar<T> && cos
( OCVar<T> && var )
{
  for( size_t i=0; i<var._coef.size(); ++i )
    var._coef[i] = Op<T>::cos( var._coef[i] );  
  return std::move( var );
}

template <typename T>
inline
OCVar<T> tan
( OCVar<T> const& var )
{
  OCVar<T> var2( var );
  for( size_t i=0; i<var._coef.size(); ++i )
    var2._coef[i] = Op<T>::tan( var._coef[i] );  
  return var2;  
}

template <typename T>
inline
OCVar<T> && tan
( OCVar<T> && var )
{
  for( size_t i=0; i<var._coef.size(); ++i )
    var._coef[i] = Op<T>::tan( var._coef[i] );  
  return std::move( var );
}

template <typename T>
inline
OCVar<T> asin
( OCVar<T> const& var )
{
  OCVar<T> var2( var );
  for( size_t i=0; i<var._coef.size(); ++i )
    var2._coef[i] = Op<T>::asin( var._coef[i] );  
  return var2;  
}

template <typename T>
inline
OCVar<T> && asin
( OCVar<T> && var )
{
  for( size_t i=0; i<var._coef.size(); ++i )
    var._coef[i] = Op<T>::asin( var._coef[i] );  
  return std::move( var );
}

template <typename T>
inline
OCVar<T> acos
( OCVar<T> const& var )
{
  OCVar<T> var2( var );
  for( size_t i=0; i<var._coef.size(); ++i )
    var2._coef[i] = Op<T>::acos( var._coef[i] );  
  return var2;  
}

template <typename T>
inline
OCVar<T> && acos
( OCVar<T> && var )
{
  for( size_t i=0; i<var._coef.size(); ++i )
    var._coef[i] = Op<T>::acos( var._coef[i] );  
  return std::move( var );
}

template <typename T>
inline
OCVar<T> atan
( OCVar<T> const& var )
{
  OCVar<T> var2( var );
  for( size_t i=0; i<var._coef.size(); ++i )
    var2._coef[i] = Op<T>::atan( var._coef[i] );  
  return var2;  
}

template <typename T>
inline
OCVar<T> && atan
( OCVar<T> && var )
{
  for( size_t i=0; i<var._coef.size(); ++i )
    var._coef[i] = Op<T>::atan( var._coef[i] );  
  return std::move( var );
}

template <typename T>
inline
OCVar<T> sinh
( OCVar<T> const& var )
{
  OCVar<T> var2( var );
  for( size_t i=0; i<var._coef.size(); ++i )
    var2._coef[i] = Op<T>::sinh( var._coef[i] );  
  return var2;  
}

template <typename T>
inline
OCVar<T> && sinh
( OCVar<T> && var )
{
  for( size_t i=0; i<var._coef.size(); ++i )
    var._coef[i] = Op<T>::sinh( var._coef[i] );  
  return std::move( var );
}

template <typename T>
inline
OCVar<T> cosh
( OCVar<T> const& var )
{
  OCVar<T> var2( var );
  for( size_t i=0; i<var._coef.size(); ++i )
    var2._coef[i] = Op<T>::cosh( var._coef[i] );  
  return var2;  
}

template <typename T>
inline
OCVar<T> && cosh
( OCVar<T> && var )
{
  for( size_t i=0; i<var._coef.size(); ++i )
    var._coef[i] = Op<T>::cosh( var._coef[i] );  
  return std::move( var );
}

template <typename T>
inline
OCVar<T> tanh
( OCVar<T> const& var )
{
  OCVar<T> var2( var );
  for( size_t i=0; i<var._coef.size(); ++i )
    var2._coef[i] = Op<T>::tanh( var._coef[i] );  
  return var2;  
}

template <typename T>
inline
OCVar<T> && tanh
( OCVar<T> && var )
{
  for( size_t i=0; i<var._coef.size(); ++i )
    var._coef[i] = Op<T>::tanh( var._coef[i] );  
  return std::move( var );
}

template <typename T>
inline
OCVar<T> xlog
( OCVar<T> const& var )
{
  OCVar<T> var2( var );
  for( size_t i=0; i<var._coef.size(); ++i )
    var2._coef[i] = Op<T>::xlog( var._coef[i] );  
  return var2;  
}

template <typename T>
inline
OCVar<T> && xlog
( OCVar<T> && var )
{
  for( size_t i=0; i<var._coef.size(); ++i )
    var._coef[i] = Op<T>::xlog( var._coef[i] );  
  return std::move( var );
}

template <typename T>
inline
OCVar<T> fabs
( OCVar<T> const& var )
{
  OCVar<T> var2( var );
  for( size_t i=0; i<var._coef.size(); ++i )
    var2._coef[i] = Op<T>::fabs( var._coef[i] );  
  return var2;  
}

template <typename T>
inline
OCVar<T> && fabs
( OCVar<T> && var )
{
  for( size_t i=0; i<var._coef.size(); ++i )
    var._coef[i] = Op<T>::fabs( var._coef[i] );  
  return std::move( var );
}

template <typename T>
inline
OCVar<T> erf
( OCVar<T> const& var )
{
  OCVar<T> var2( var );
  for( size_t i=0; i<var._coef.size(); ++i )
    var2._coef[i] = Op<T>::erf( var._coef[i] );  
  return var2;  
}

template <typename T>
inline
OCVar<T> && erf
( OCVar<T> && var )
{
  for( size_t i=0; i<var._coef.size(); ++i )
    var._coef[i] = Op<T>::erf( var._coef[i] );  
  return std::move( var );
}

template <typename T>
inline
OCVar<T> erfc
( OCVar<T> const& var )
{
  OCVar<T> var2( var );
  for( size_t i=0; i<var._coef.size(); ++i )
    var2._coef[i] = Op<T>::erfc( var._coef[i] );  
  return var2;  
}

template <typename T>
inline
OCVar<T> && erfc
( OCVar<T> && var )
{
  for( size_t i=0; i<var._coef.size(); ++i )
    var._coef[i] = Op<T>::erfc( var._coef[i] );  
  return std::move( var );
}

template <typename T>
inline
OCVar<T> fstep
( OCVar<T> const& var )
{
  OCVar<T> var2( var );
  for( size_t i=0; i<var._coef.size(); ++i )
    var2._coef[i] = Op<T>::fstep( var._coef[i] );  
  return var2;  
}

template <typename T>
inline
OCVar<T> && fstep
( OCVar<T> && var )
{
  for( size_t i=0; i<var._coef.size(); ++i )
    var._coef[i] = Op<T>::fstep( var._coef[i] );  
  return std::move( var );
}

template <typename T>
inline
OCVar<T> bstep
( OCVar<T> const& var )
{
  OCVar<T> var2( var );
  for( size_t i=0; i<var._coef.size(); ++i )
    var2._coef[i] = Op<T>::bstep( var._coef[i] );  
  return var2;  
}

template <typename T>
inline
OCVar<T> && bstep
( OCVar<T> && var )
{
  for( size_t i=0; i<var._coef.size(); ++i )
    var._coef[i] = Op<T>::bstep( var._coef[i] );  
  return std::move( var );
}

template <typename T>
inline bool
OCBase::lift_colloc
( std::vector<T>& coeff_out, const std::vector<T>& coeff_in,
  const std::map<FFVar const*, OCDom const*, lt_FFVar>& dom_out,
  const std::map<FFVar const*, OCDom const*, lt_FFVar>& dom_in,
  const bool check )
{
  if( check && dom_in.size() > dom_out.size() ) return false;

  // Total output size
  size_t nout = 1;
  for( const auto& d : dom_out ){
    if( check && ( !d.first || !d.second || d.second->n_node == 0 ) ) return false;
    nout *= d.second->n_node;
  }

  // Scalar input: replicate over all lifted dimensions.
  if( dom_in.empty() ){
    if( check && coeff_in.size() != 1 ) return false;
    coeff_out.assign( nout, coeff_in.empty()? T(0.): coeff_in[0] );
    return true;
  }

  if( check && dom_out.empty() ) return false;

  std::vector<std::pair<FFVar const*,OCDom const*>> outdom;
  std::vector<std::pair<FFVar const*,OCDom const*>> indom;
  outdom.reserve( dom_out.size() );
  indom .reserve( dom_in .size() );
  for( auto const& d : dom_out ) outdom.push_back( d );
  for( auto const& d : dom_in  ) indom .push_back( d );

  size_t expected_in_size = 1;
  std::vector<size_t> match_out( indom.size(), static_cast<size_t>(-1) );
  for( size_t i = 0; i < indom.size(); ++i ){
    if( check && ( !indom[i].first || !indom[i].second || indom[i].second->n_node == 0 ) )
      return false;
    expected_in_size *= indom[i].second->n_node;

    for( size_t j = 0; j < outdom.size(); ++j ){
      if( outdom[j].first->id() == indom[i].first->id() ){
        match_out[i] = j;
        break;
      }
    }
    if( match_out[i] == static_cast<size_t>(-1) ) return false;
  }
  if( check && coeff_in.size() != expected_in_size ) return false;

  coeff_out.resize( nout );

  std::vector<size_t> idx_out( outdom.size(), 0 );
  std::vector<std::vector<double>> weights( indom.size() );

  for( size_t flat_out = 0; flat_out < nout; ++flat_out ){

    size_t tmp = flat_out;
    for( size_t j = 0; j < outdom.size(); ++j ){
      idx_out[j] = tmp % outdom[j].second->n_node;
      tmp /= outdom[j].second->n_node;
    }

    for( size_t i = 0; i < indom.size(); ++i ){
      auto const* dom_eval = outdom[ match_out[i] ].second;
      auto const* dom_coef = indom [ i ].second;
      double const xi = dom_eval->lgnodes().at( idx_out[ match_out[i] ] );
      if( !dom_coef->w_eval_lagrange( weights[i], xi, 0 ) ) return false;
      if( check && weights[i].size() != dom_coef->n_node ) return false;
    }

    T val( 0. );
    std::function<void(size_t,size_t,size_t,double)> accumulate =
      [&]( size_t idim, size_t flat_in, size_t stride_in, double wprod )
      {
        if( idim == indom.size() ){
          val += coeff_in[flat_in] * wprod;
          return;
        }
        size_t const n = indom[idim].second->n_node;
        for( size_t k = 0; k < n; ++k )
          accumulate( idim+1, flat_in + k*stride_in, stride_in*n,
                      wprod * weights[idim][k] );
      };

    accumulate( 0, 0, 1, 1. );
    coeff_out[flat_out] = val;
  }

  return true;
}

inline std::pair<size_t, size_t>
OCBase::stride
( std::map<FFVar const*, OCDom const*, lt_FFVar> const& dom,
  FFVar const* index, bool const check )
{
  auto itndx = dom.find( index );
  if( check && ( itndx == dom.end() || !itndx->second || itndx->second->n_node == 0 ) )
    return {0, 0};

  // stride = product of dimensions before index
  size_t stride = 1;
  for( auto it = dom.begin(); it != itndx; ++it ){
    if( check && ( !it->second || it->second->n_node == 0 ) )
      return {0, 0};
    stride *= it->second->n_node;
  }
  //size_t const n_index = itndx->second->n_node;

  // outer = product of dimensions after index
  size_t outer = 1;
  for( auto it = std::next(itndx); it != dom.end(); ++it ){
    if( check && ( !it->second || it->second->n_node == 0 ) ) 
      return {0, 0};
    outer *= it->second->n_node;
  }

  //if( check && coef.size() != stride * n_index * outer )
  //  return {0, 0};

  return {stride, outer};
}

} // namespace mc


#ifdef MC__USE_FADBAD

#include "mcfadbad.hpp"

namespace fadbad
{

//! @brief Specialization of the structure fadbad::Op for use of the type mc::OCVar of MC++ as a template parameter of the classes fadbad::F, fadbad::B and fadbad::T of FADBAD++
template <typename T> struct Op<mc::OCVar<T>>
{
  typedef mc::OCVar<T> OCV;
  typedef double Base;
  static Base myInteger( int const i ) { return Base(i); }
  static Base myZero() { return myInteger(0); }
  static Base myOne()  { return myInteger(1);}
  static Base myTwo()  { return myInteger(2); }
  static double myPI() { return mc::PI; }
  static OCV myPos( OCV const& x ) { return  x; }
  static OCV myNeg( OCV const& x ) { return -x; }
  template <typename U> static OCV& myCadd( OCV& x, U const& y ) { return x += y; }
  template <typename U> static OCV& myCsub( OCV& x, U const& y ) { return x -= y; }
  template <typename U> static OCV& myCmul( OCV& x, U const& y ) { return x *= y; }
  template <typename U> static OCV& myCdiv( OCV& x, U const& y ) { return x /= y; }
  static OCV myInv( OCV const& x ) { return mc::inv( x ); }
  static OCV mySqr( OCV const& x ) { return mc::sqr( x ); }
  template <typename X, typename Y> static OCV myPow( X const& x, Y const& y ) { return mc::pow( x, y ); }
  static OCV mySqrt( OCV const& x ){ return mc::sqrt(x); }
  static OCV myLog( OCV const& x ) { return mc::log( x ); }
  static OCV myExp( OCV const& x ) { return mc::exp( x ); }
  static OCV mySin( OCV const& x ) { return mc::sin( x ); }
  static OCV myCos( OCV const& x ) { return mc::cos( x ); }
  static OCV myTan( OCV const& x )  { return mc::tan( x ); }
  static OCV myAsin( OCV const& x ) { return mc::asin( x ); }
  static OCV myAcos( OCV const& x ) { return mc::acos( x ); }
  static OCV myAtan( OCV const& x ) { return mc::atan( x ); }
  static OCV mySinh( OCV const& x ) { return mc::sinh( x ); }
  static OCV myCosh( OCV const& x ) { return mc::cosh( x ); }
  static OCV myTanh( OCV const& x ) { return mc::tanh( x ); }
  static bool myEq( OCV const& x, OCV const& y ) { throw mc::OCBase::Exceptions( mc::OCBase::Exceptions::UNDEF ); }
  static bool myNe( OCV const& x, OCV const& y ) { throw mc::OCBase::Exceptions( mc::OCBase::Exceptions::UNDEF ); }
  static bool myLt( OCV const& x, OCV const& y ) { throw mc::OCBase::Exceptions( mc::OCBase::Exceptions::UNDEF ); }
  static bool myLe( OCV const& x, OCV const& y ) { throw mc::OCBase::Exceptions( mc::OCBase::Exceptions::UNDEF ); }
  static bool myGt( OCV const& x, OCV const& y ) { throw mc::OCBase::Exceptions( mc::OCBase::Exceptions::UNDEF ); }
  static bool myGe( OCV const& x, OCV const& y ) { throw mc::OCBase::Exceptions( mc::OCBase::Exceptions::UNDEF ); }
};

} // end namespace fadbad

#endif


namespace mc
{

//! @brief C++ structure for specialization of the mc::Op templated structure for use of mc::OCVar in DAG evaluation and as template parameter in other MC++ types
template <typename T> struct Op<mc::OCVar<T>>
{
  typedef mc::OCVar<T> OCV;
  static OCV point( const double c ) { return OCV(c); }
  static OCV zeroone() { throw mc::OCBase::Exceptions( mc::OCBase::Exceptions::UNDEF ); }
  static void I(OCV& x, OCV const&y) { x = y; }
  static double l(OCV const& x) { throw mc::OCBase::Exceptions( mc::OCBase::Exceptions::UNDEF ); }
  static double u(OCV const& x) { throw mc::OCBase::Exceptions( mc::OCBase::Exceptions::UNDEF ); }
  static double abs (OCV const& x) { throw mc::OCBase::Exceptions( mc::OCBase::Exceptions::UNDEF ); }
  static double mid (OCV const& x) { throw mc::OCBase::Exceptions( mc::OCBase::Exceptions::UNDEF ); }
  static double diam(OCV const& x) { throw mc::OCBase::Exceptions( mc::OCBase::Exceptions::UNDEF ); }
  static OCV inv (OCV const& x) { return mc::inv(x);  }
  static OCV sqr (OCV const& x) { return mc::sqr(x);  }
  static OCV sqrt(OCV const& x) { return mc::sqrt(x); }
  static OCV exp (OCV const& x) { return mc::exp(x);  }
  static OCV log (OCV const& x) { return mc::log(x);  }
  static OCV xlog(OCV const& x) { return mc::xlog(x); }
  static OCV lmtd(OCV const& x, OCV const& y) { return (x-y)/(mc::log(x)-mc::log(y)); }
  static OCV rlmtd(OCV const& x, OCV const& y) { return (mc::log(x)-mc::log(y))/(x-y); }
  static OCV fabs(OCV const& x) { return mc::fabs(x); }
  static OCV sin (OCV const& x) { return mc::sin(x); }
  static OCV cos (OCV const& x) { return mc::cos(x); }
  static OCV tan (OCV const& x) { return mc::tan(x); }
  static OCV asin(OCV const& x) { return mc::asin(x); }
  static OCV acos(OCV const& x) { return mc::acos(x); }
  static OCV atan(OCV const& x) { return mc::atan(x); }
  static OCV sinh(OCV const& x) { return mc::sinh(x); }
  static OCV cosh(OCV const& x) { return mc::cosh(x); }
  static OCV tanh(OCV const& x) { return mc::tanh(x); }
  static OCV erf (OCV const& x) { return mc::erf(x); }
  static OCV erfc(OCV const& x) { return mc::erfc(x); }
  static OCV fstep(OCV const& x) { throw mc::fstep(x); }
  static OCV bstep(OCV const& x) { throw mc::bstep(x); }
  template <typename Y> static OCV min (OCV const& x, Y const& y) { return mc::min(x,y); }
  template <typename Y> static OCV max (OCV const& x, Y const& y) { return mc::max(x,y); }
  template <typename X, typename Y> static OCV pow(const X& x, const Y& y) { return mc::pow(x,y); }
  static OCV cheb(OCV const& x, const unsigned n) { return mc::cheb(x,n); }
  static OCV prod (const unsigned n, const OCV* x) { return mc::prod(n,x); }
  static OCV monom (const unsigned n, const OCV* x, const unsigned* k) { throw mc::OCBase::Exceptions( mc::OCBase::Exceptions::UNDEF ); }
  static bool inter(OCV& xIy, OCV const& x, OCV const& y) { throw mc::OCBase::Exceptions( mc::OCBase::Exceptions::UNDEF ); }
  static OCV hull(OCV const& x, OCV const& y) { throw mc::OCBase::Exceptions( mc::OCBase::Exceptions::UNDEF ); }
  static bool eq(OCV const& x, OCV const& y) { throw mc::OCBase::Exceptions( mc::OCBase::Exceptions::UNDEF ); }
  static bool ne(OCV const& x, OCV const& y) { throw mc::OCBase::Exceptions( mc::OCBase::Exceptions::UNDEF ); }
  static bool lt(OCV const& x, OCV const& y) { throw mc::OCBase::Exceptions( mc::OCBase::Exceptions::UNDEF ); }
  static bool le(OCV const& x, OCV const& y) { throw mc::OCBase::Exceptions( mc::OCBase::Exceptions::UNDEF ); }
  static bool gt(OCV const& x, OCV const& y) { throw mc::OCBase::Exceptions( mc::OCBase::Exceptions::UNDEF ); }
  static bool ge(OCV const& x, OCV const& y) { throw mc::OCBase::Exceptions( mc::OCBase::Exceptions::UNDEF ); }
};

} // namespace mc

#endif
