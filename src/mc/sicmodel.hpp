// Copyright (C) Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the Eclipse Public License.

/*!
\page page_CHEBYSHEV_SPARSEINT Sparse Interval Chebyshev Model Arithmetic for Factorable Functions

The classes mc::SICModel and mc::SICVar provide an implementation of sparse interval Chebyshev model arithmetic [Zha & Chachuat, 2021]. By contrast to sparse Chebyshev model arithmetic, the distinction between the polynomial part and the remainder part is blurred in a sparse interval Chebyshev model by considering a multivariate polynomial with interval coefficients, \f$\mathcal P:\mathbb{R}^n\to\mathbb{IR}\f$, so that:
\f{align*}
  f({x}) \in \mathcal P({x}), \quad \forall {x}\in D.
\f}
These interval coefficients are propagated in the arithmetic of the templated type, typically a verified interval arithmetic. Sparse interval Chebyshev models yield a superset of sparse Chebyhev models. In general, they enable tighter inclusions by taking advantage of the distribution of the remainder term over the participating monomials. 

\section sec_CHEBYSHEV_SPARSEINT_use How do I compute a sparse interval Chebyshev model of a factorable function?

Suppose we want to compute a 4th-order sparse interval Chebyshev model for the real-valued function \f$f(x,y)=x\exp(y^2)-y^2\f$ with \f$(x,y)\in [1,2]\times[0,1]\f$. For simplicity, bounds on the remainder terms are computed using the default interval type mc::Interval here:

\code
      #include "interval.hpp"
      #include "sicmodel.hpp"
      typedef mc::Interval I;
      typedef mc::SICModel<I> SICM;
      typedef mc::SICVar<I> SICV;
\endcode


First, an mc::SICModel object is initialized by passing the order of the polynomial inclusion (4th order here):

\code
      SICM mod( 4 );
\endcode

Next, the variables \f$x\f$ and \f$y\f$ are added to the sparse interval Chebyshev model:

\code
      SICV X( &mod, 0, I(1.,2.) );
      SICV Y( &mod, 1, I(0.,1.) );
\endcode

Essentially, the first line means that <tt>X</tt> is a variable of class mc::SICVar, participating in the Chebyshev model <tt>mod</tt>, belonging to the interval \f$[1,2]\f$, and having index 0 (C-style indexing in use). The same holds for the variable <tt>Y</tt>, participating in the model <tt>mod</tt>, belonging to the interval \f$[0,1]\f$, and having index 1.

Having defined the variables, a sparse Chebyshev model of \f$f(x,y)=x\exp(x+y^2)-y^2\f$ on \f$[1,2]\times[0,1]\f$ at the mid-point \f$(\frac{3}{2},\frac{1}{2})\f$ is computed as:

\code
      SICV F = X*exp(pow(Y,2))-pow(Y,2);
\endcode

This model can be displayed to the standard output as:

\code
      std::cout << "f sparse interval Chebyshev model: " << F << std::endl;
\endcode

which produces the following output:

\verbatim
f sparse interval Chebyshev model: 
[ 1.9601057e+00, 1.9675908e+00]   0  1
[ 7.7530158e-01, 7.8393061e-01]   1  T1[0]
[ 7.0133649e-01, 7.0133649e-01]   1  T1[1]
[ 4.0044550e-01, 4.0044550e-01]   2  T1[0]·T1[1]
[ 3.0580072e-01, 3.0580072e-01]   2  T2[1]
[ 1.4360024e-01, 1.4360024e-01]   3  T1[0]·T2[1]
[ 8.4292715e-02, 8.4292715e-02]   3  T3[1]
[ 2.8097572e-02, 2.8097572e-02]   4  T1[0]·T3[1]
[ 1.5334992e-02, 2.1468989e-02]   4  T4[1]
   B     =  [-5.0886708e-01, 4.4365637e+00]
\endverbatim

The interval coefficients and corresponding monomials in the Chebyshev model are displayed first, with the middle integer being the total monomial order. The Chebyshev model range estimator is reported last.

Other operations involve retreiving the model range, or computing the range of the interval-valued polynomial at a given point:

\code
      I B = F.B();
      double x[2] = { 1.5, 0.5 };
      I Pval = F.IP( x );
\endcode

See the documentations of mc::SICModel and mc::SICVar for a complete list of member functions. 


\section sec_CHEBYSHEV_SPARSEINT_fct Which functions are overloaded for sparse interval Chebyshev model arithmetic?

mc::SICVar overloads the usual functions <tt>exp</tt>, <tt>log</tt>, <tt>sqr</tt>, <tt>sqrt</tt>, <tt>pow</tt>, <tt>inv</tt>, <tt>cos</tt>, <tt>sin</tt>, <tt>tan</tt>, <tt>acos</tt>, <tt>asin</tt>, <tt>atan</tt>, <tt>cosh</tt>, <tt>sinh</tt>, <tt>tanh</tt>, <tt>fabs</tt>. The functions <tt>erf</tt>, <tt>erfc</tt>, <tt>min</tt>, <tt>max</tt>, <tt>fstep</tt>, <tt>bstep</tt>, <tt>inter</tt>, and <tt>hull</tt> are not currently overloaded in mc::SICVar.

\section sec_CHEBYSHEV_SPARSEINT_opt How are the options set for the computation of a sparse interval Chebyshev model?

The class mc::SICModel has a public member called mc::SICModel::options that can be used to set/modify the options. For instance:

\code
      mod.options.BOUNDER_TYPE = SICM::Options::NAIVE;
      mod.options.SCALE_VARIABLES = true;
\endcode

The available options are the following:

<TABLE border="1">
<CAPTION><EM>Options in mc::SICModel::Options: name, type and description</EM></CAPTION>
     <TR><TH><b>Name</b>  <TD><b>Type</b><TD><b>Default</b>
         <TD><b>Description</b>
     <TR><TH><tt>BASIS</tt> <TD><tt>unsigned int</tt> <TD>1
         <TD>Basis representation of the monomials: 0-Monomial basis; 1-Chebyshev basis.
     <TR><TH><tt>HOT_SPLIT</tt> <TD><tt>mc::SICModel::Options::ALLOCATION</tt> <TD>mc::SICModel::Options::FULL
         <TD>Strategy for allocating higher-order terms to monomials in sparse product terms.
     <TR><TH><tt>PRODBND_SPLIT</tt> <TD><tt>bool</tt> <TD>false
         <TD>Whether to distribute uncertainty in product with interval bound.
     <TR><TH><tt>REMEZ_USE</tt> <TD><tt>bool</tt> <TD>false
         <TD>Whether to use the Remez algorithm for computing a minimax approximation for univariate terms.
     <TR><TH><tt>REMEZ_MAXIT</tt> <TD><tt>unsigned int</tt> <TD>10
         <TD>Maximal number of iterations in Remez algorithm for computing a minimax approximation for univariate terms.
     <TR><TH><tt>REMEZ_TOL</tt> <TD><tt>double</tt> <TD>1e-5
         <TD>Stopping tolerance in Remez algorithm for computing a minimax approximation for univariate terms.
     <TR><TH><tt>REMEZ_MIG</tt> <TD><tt>double</tt> <TD>1e-10
         <TD>Threshold for interval width below which Remez algorithm is not useds.
     <TR><TH><tt>INTERP_EXTRA</tt> <TD><tt>unsigned int</tt> <TD>0
         <TD>Extra terms in Chebyshev interpolation of univariates: 0-Chebyshev interpolation of order NORD; extra terms allow approximation of Chebyshev truncated series.
     <TR><TH><tt>INTERP_THRES</tt> <TD><tt>double</tt> <TD>1e2*machprec()
         <TD>Threshold for coefficient values in Chebyshev expansion for bounding of transcendental univariates.
     <TR><TH><tt>BOUNDER_TYPE</tt> <TD><tt>mc::SICModel::Options::BOUNDER</tt> <TD>mc::SICModel::Options::LSB
         <TD>Chebyshev model range bounder.
     <TR><TH><tt>MIXED_IA</tt> <TD><tt>bool</tt> <TD>false
         <TD>Whether to intersect internal bounds with underlying bounds in the templated arithmetics.
     <TR><TH><tt>MIN_FACTOR</tt> <TD><tt>double</tt> <TD>0e0
         <TD>Threshold for monomial coefficients below which the term is removed and appended to the remainder term.
     <TR><TH><tt>REF_POLY</tt> <TD><tt>double</tt> <TD>0.
         <TD>Scalar in \f$[0,1]\f$ related to the choice of the polynomial part in the overloaded functions mc::inter and mc::hull (see \ref sec_CHEBYSHEV_SPARSEINT_fct). A value of 0. amounts to selecting the polynomial part of the left operand, whereas a value of 1. selects the right operand.
     <TR><TH><tt>DISPLAY_DIGITS</tt> <TD><tt>unsigned int</tt> <TD>5
         <TD>Number of digits in output stream for Chebyshev model coefficients.
</TABLE>


\section sec_CHEBYSHEV_SPARSEINT_err Errors What errors can I encounter during computation of a sparse interval Chebyshev model?

Errors are managed based on the exception handling mechanism of the C++ language. Each time an error is encountered, a class object of type mc::SICModel::Exceptions is thrown, which contains the type of error. Possible errors encountered during the computation of a sparse interval Chebyshev model are:

<TABLE border="1">
<CAPTION><EM>Errors from class mc::CModel::Exceptions during the Computation of a Chebyshev Model</EM></CAPTION>
     <TR><TH><b>Number</b> <TD><b>Description</b>
     <TR><TH><tt>1</tt> <TD>Division by zero scalar
     <TR><TH><tt>2</tt> <TD>Inverse operation with zero in range
     <TR><TH><tt>3</tt> <TD>Log operation with non-positive numbers in range
     <TR><TH><tt>4</tt> <TD>Square-root operation with negative numbers in range
     <TR><TH><tt>5</tt> <TD>Real power operation with negative numbers in range
     <TR><TH><tt>6</tt> <TD>Tangent operation with (k+1/2)·PI in range
     <TR><TH><tt>7</tt> <TD>Cosine inverse operation with range outside [-1,1]
     <TR><TH><tt>8</tt> <TD>Sine inverse operation with range outside [-1,1]
     <TR><TH><tt>9</tt> <TD>Failed to compose sparse Chebyshev variable
     <TR><TH><tt>-1</tt> <TD>Failed to construct sparse Chebyshev variable
     <TR><TH><tt>-2</tt> <TD>Inconsistent bounds with template parameter arithmetic
     <TR><TH><tt>-3</tt> <TD>Operation between Chebyshev variables linked to different Chebyshev models
     <TR><TH><tt>-4</tt> <TD>Internal error
     <TR><TH><tt>-33</tt> <TD>Feature not yet implemented in mc::SICModel
</TABLE>

Further exceptions may be thrown by the template parameter class itself.

*/

#ifndef MC__SICMODEL_H
#define MC__SICMODEL_H

#include <iostream>
#include <iomanip>
#include <typeinfo>
#include <sstream>
#include <string>
#include <vector>
#include <set>
#include <map>
#include <stdarg.h>
#include <cassert>
#include <climits>
#include <limits>
#include <stdlib.h>
#include <complex>
#include <numeric>

#include "mcfunc.hpp"
#include "mclapack.hpp"
#include "mcop.hpp"
#include "spolymon.hpp"
#include "remez.hpp"

#undef  MC__SICMODEL_DEBUG
#undef  MC__SICMODEL_DEBUG_SCALE
#define MC__SICMODEL_CHECK
#undef  MC__SICMODEL_CHECK_PMODEL
#undef  MC__SICVAR_DEBUG_EXP
#undef  MC__SICVAR_DEBUG_BERNSTEIN
#define MC__SICMODEL_USE_PROD
#undef  MC__SICMODEL_DEBUG_MINIMAX

//#undef MC__SCVAR_FABS_SQRT
//#undef MC__SCVAR_FORCE_REM_DERIV
//#undef MC__SCVAR_SPARSE_PRODUCT_FULL

namespace mc
{

template <typename T> class SCModel;
template <typename T> class SICVar;

//! @brief C++ class for the computation of sparse interval Chebyshev models for factorable function: environment
////////////////////////////////////////////////////////////////////////
//! mc::SICModel is a C++ class for definition of interval Chebyshev
//! model environment in sparse format. Propagation of sparse interval
//! Chebyshev models for factorable functions is via the C++ class
//! mc::SICVar. The template parameter corresponds to the type used to
//! propagate the polynomial coefficients.

////////////////////////////////////////////////////////////////////////
template <typename T>
class SICModel
////////////////////////////////////////////////////////////////////////
{
  friend class SICVar<T>;
  template <typename U> friend class SICModel;
  template <typename U> friend SICVar<U> inv
    ( const SICVar<U>& );
  template <typename U> friend SICVar<U> sqr
    ( const SICVar<U>& );
  template <typename U> friend SICVar<U> sqrt
    ( const SICVar<U>& );
  template <typename U> friend SICVar<U> exp
    ( const SICVar<U>& );
  template <typename U> friend SICVar<U> log
    ( const SICVar<U>& );
  template <typename U> friend SICVar<U> xlog
    ( const SICVar<U>& );
  template <typename U> friend SICVar<U> pow
    ( const SICVar<U>&, const int );
  template <typename U> friend SICVar<U> pow
    ( const SICVar<U>&, double const );
  template <typename U> friend SICVar<U> pow
    ( double const, const SICVar<U>& );
  template <typename U> friend SICVar<U> pow
    ( const SICVar<U>&, const SICVar<U>& );
  template <typename U> friend SICVar<U> prod
    ( const unsigned int, const SICVar<U>* );
  template <typename U> friend SICVar<U> monom
    ( const unsigned int, const SICVar<U>*, const unsigned* );
  template <typename U> friend SICVar<U> cheb
    ( const SICVar<U>&, const unsigned );
  template <typename U> friend SICVar<U> cos
    ( const SICVar<U>& );
  template <typename U> friend SICVar<U> sin
    ( const SICVar<U>& );
  template <typename U> friend SICVar<U> tan
    ( const SICVar<U>& );
  template <typename U> friend SICVar<U> acos
    ( const SICVar<U>& );
  template <typename U> friend SICVar<U> asin
    ( const SICVar<U>& );
  template <typename U> friend SICVar<U> atan
    ( const SICVar<U>& );
  template <typename U> friend SICVar<U> cosh
    ( const SICVar<U>& );
  template <typename U> friend SICVar<U> sinh
    ( const SICVar<U>& );
  template <typename U> friend SICVar<U> tanh
    ( const SICVar<U>& );
  template <typename U> friend SICVar<U> fabs
    ( const SICVar<U>& );

  typedef std::map< SPolyMon, T, lt_SPolyMon > t_coefmon;
  typedef std::set< unsigned > t_var;

protected:

  //! @brief Max order of polynomial model
  unsigned _maxord;

  //! @brief Max number of variables in polynomial model
  unsigned _maxvar;

  //! @brief Set of particpating variable indices
  std::set<unsigned> _setvar;

  //! @brief Bounds on model variables
  std::vector<T> _bndvar;

  //! @brief Reference points for the model variables
  std::vector<double> _refvar;

  //! @brief Scaling factors of the model variables
  std::vector<double> _scalvar; 

  //! @brief Resize the model variable data containers
  void _resize
    ( const unsigned ivar )
    { if( ivar < _maxvar ) return; _maxvar = ivar+1;
      _bndvar.resize( _maxvar ); _refvar.resize( _maxvar ); _scalvar.resize( _maxvar ); }
      
  //! @brief Resize the model variable data containers
  void _set
    ( const unsigned i, T const& X )
    { _resize( i ); _setvar.insert( i );
      _bndvar[i] = X; _refvar[i] = Op<T>::mid(X); _scalvar[i] = 0.5*Op<T>::diam(X); }

  //! @brief Resize and return a pointer to the Chebyshev coefficient interpolant array
  std::vector<double>& _resize_coefuniv
    ()
    { _coefuniv.resize( _maxord + options.INTERP_EXTRA + 1 );
      return _coefuniv; }

public:
  /** @addtogroup SICHEBYSHEV Sparse Interval Chebyshev Model Arithmetic for Factorable Functions
   *  @{
   */
  //! @brief Constructor of Sparse Chebyshev model environment for maximal order <tt>maxord</tt>
  SICModel
    ( const unsigned maxord=3, const unsigned nvarres=0 )
    : _maxord( maxord ), _maxvar( 0 )
    { _coefuniv.reserve( _maxord );
      if( !nvarres ) return;
      _bndvar.reserve( nvarres );
      _refvar.reserve( nvarres );
      _scalvar.reserve( nvarres ); }

  //! @brief Copy constructor of Sparse Chebyshev model environment
  SICModel
    ( SICModel<T> const& mod )
    : _maxord( mod._maxord ), _maxvar( mod._maxvar ), _setvar( mod._setvar ),
      _bndvar( mod._bndvar ), _refvar( mod._refvar ), _scalvar( mod._scalvar )
    { _coefuniv.reserve( _maxord ); }

  //! @brief Copy constructor of Sparse Chebyshev model environment
  SICModel
    ( SCModel<T> const& mod )
    : _maxord( mod.maxord() ), _maxvar( mod.maxvar() ), _setvar( mod.setvar() ),
      _bndvar( mod.bndvar() ), _refvar( mod.refvar() ), _scalvar( mod.scalvar() )
    { _coefuniv.reserve( _maxord ); }

  //! @brief Destructor of Sparse Chebyshev model environment
  ~SICModel()
    {}

  //! @brief Maximal order of polynomial model
  unsigned maxord
    ()
    const
    { return _maxord; };

  //! @brief Number of variables in polynomial model
  unsigned maxvar
    ()
    const
    { return _maxvar; };

  //! @brief Const reference to reference points for the model variables
  std::set<unsigned> const& setvar() const
    { return _setvar; }

  //! @brief Const reference to the bounds of the model variables
  std::vector<T> const& bndvar() const
    { return _bndvar; }

  //! @brief Const reference to reference points for the model variables
  std::vector<double> const& refvar() const
    { return _refvar; }

  //! @brief Const reference to the scaling factors of the model variables
  std::vector<double> const& scalvar() const
    { return _scalvar; }

  //! @brief Get Chebyshev basis functions in U arithmetic for variable array <a>X</a>
  template <typename U> U** get_basis
    ( const unsigned maxord, const U*bndvar, const bool scaled=false ) const;

  //! @brief Get Chebyshev monomial bounds in U arithmetic for variable array <a>X</a> and monomial indexes in <a>ndxmon</a>
  template <typename U> void get_bndmon
    ( std::map<SPolyMon,U,lt_SPolyMon>& bndmon, U const* bndvar,
      bool const scaled=false ) const;

  //! @brief Exceptions of mc::SICModel
  class Exceptions
  {
  public:
    //! @brief Enumeration type for SICModel exception handling
    enum TYPE{
      DIV=1,	//!< Division by zero scalar
      INV,	//!< Inverse operation with zero in range
      LOG,	//!< Log operation with non-positive numbers in range
      SQRT,	//!< Square-root operation with negative numbers in range
      DPOW,    //!< Real power operation with negative numbers in range
      TAN,	//!< Tangent operation with (k+1/2)·PI in range
      ACOS,	//!< Cosine inverse operation with range outside [-1,1]
      ASIN,	//!< Sine inverse operation with range outside [-1,1]
      COMPOSE,	//!< Failed to compose Chebyshev variable
      INIT=-1,	//!< Failed to construct Chebyshev variable
      INCON=-2, //!< Inconsistent bounds with template parameter arithmetic
      SCMODEL=-3,//!< Operation between Chebyshev variables linked to different Chebyshev models
      INTERNAL = -4,//!< Internal error
      UNDEF=-33 //!< Feature not yet implemented in mc::SCModel
    };
    //! @brief Constructor for error <a>ierr</a>
    Exceptions( TYPE ierr ) : _ierr( ierr ){}
    //! @brief Error flag
    int ierr(){ return _ierr; }
    //! @brief Error description
    std::string what(){
      switch( _ierr ){
      case DIV:
        return "mc::SCModel\t Division by zero scalar";
      case INV:
        return "mc::SCModel\t Inverse operation with zero in range";
      case LOG:
        return "mc::SCModel\t Log operation with non-positive numbers in range";
      case SQRT:
        return "mc::SCModel\t Square-root operation with negative numbers in range";
      case DPOW:
        return "mc::SCModel\t Real power operation with negative numbers in range";
      case TAN:
        return "mc::SCModel\t Tangent operation with (k+1/2)·PI in range";
      case ACOS:
        return "mc::SCModel\t Cosine inverse operation with range outside [-1,1]";
      case ASIN:
        return "mc::SCModel\t Sine inverse operation with range outside [-1,1]";
      case COMPOSE:
        return "mc::SCModel\t Chebyshev composition failed";
      case INIT:
        return "mc::SCModel\t Chebyshev variable initialization failed";
      case INCON:
        return "mc::SCModel\t Inconsistent bounds with template parameter arithmetic";
      case SCMODEL:
        return "mc::SCModel\t Operation between Chebyshev variables in different Chebyshev model environment";
      case UNDEF:
        return "mc::SCModel\t Feature not yet implemented in mc::SCModel class";
      case INTERNAL:
      default:
        return "mc::SCModel\t Internal error";
      }
    }

  private:
    TYPE _ierr;
  };

  //! @brief Options of mc::SICModel
  struct Options
  {
    //! @brief Constructor of mc::SICModel::Options
    Options():
      BASIS(CHEB), HOT_SPLIT(FULL), PRODBND_SPLIT(false),
      REMEZ_USE(false), REMEZ_MAXIT(10), REMEZ_TOL(1e-5), REMEZ_MIG(1e-10),
      INTERP_EXTRA(0), INTERP_THRES(1e2*machprec()), BOUNDER_TYPE(LSB),
      MIXED_IA(false), MIN_FACTOR(0.), REF_POLY(0.),
      DISPLAY_DIGITS(7)
      {}
    //! @brief Copy constructor of mc::SICModel::Options
    template <typename U> Options
      ( U&options )
      : BASIS( options.BASIS ),
        HOT_SPLIT( options.HOT_SPLIT ),
        PRODBND_SPLIT( options.PRODBND_SPLIT ),
        REMEZ_USE( options.REMEZ_USE ),
        REMEZ_MAXIT( options.REMEZ_MAXIT ),
        REMEZ_TOL( options.REMEZ_TOL ),
        REMEZ_MIG( options.REMEZ_MIG ),
        INTERP_EXTRA( options.INTERP_EXTRA ),
        INTERP_THRES( options.INTERP_THRES ),
        BOUNDER_TYPE( options.BOUNDER_TYPE ),
        MIN_FACTOR( options.MIN_FACTOR ),
        MIXED_IA( options.MIXED_IA ),
        REF_POLY(options.REF_POLY),
        DISPLAY_DIGITS( options.DISPLAY_DIGITS )
      {}
    //! @brief Assignment of mc::SICModel::Options
    template <typename U> Options& operator =
      ( U&options ){
        BASIS            = options.BASIS;
        HOT_SPLIT        = options.HOT_SPLIT;
        PRODBND_SPLIT    = options.PRODBND_SPLIT;
        REMEZ_USE        = options.REMEZ_USE;
        REMEZ_MAXIT      = options.REMEZ_MAXIT;
        REMEZ_TOL        = options.REMEZ_TOL;
        REMEZ_MIG        = options.REMEZ_MIG;
        INTERP_EXTRA     = options.INTERP_EXTRA;
        INTERP_THRES     = options.INTERP_THRES;
        BOUNDER_TYPE     = (BOUNDER)options.BOUNDER_TYPE;
        MIN_FACTOR       = options.MIN_FACTOR;
        MIXED_IA         = options.MIXED_IA;
        REF_POLY         = options.REF_POLY;
        DISPLAY_DIGITS   = options.DISPLAY_DIGITS;
        return *this;
      }
    //! @brief Chebyshev model range bounder option
    enum BOUNDER{
      NAIVE=0,	//!< Naive polynomial range bounder
      LSB 	//!< Lin & Stadtherr range bounder
    };
    //! @brief Strategies for allocating higher-order terms to monomials in sparse product terms
    enum ALLOCATION{
      NONE=0,	//!< No split of HOT and allocation to constant coefficient
      SIMPLE,	//!< Split HOT between current term and new variable
      FULL 	//!< Split HOT among all variables
    };
    //! @brief Available basis representations
    enum MONBASIS{
      POW=0,	//!< Power basis
      CHEB	//!< Chebyshev basis
    };
    //! @brief Basis representation of the monomials
    unsigned BASIS;
    //! @brief Strategy for allocating higher-order terms to monomials in sparse product terms
    ALLOCATION HOT_SPLIT;
    //! @brief Whether to distribute uncertainty in product with interval bound
    bool PRODBND_SPLIT;
    //! @brief Whether to use the Remez algorithm for computing a minimax approximation for univariate terms
    bool REMEZ_USE;
    //! @brief Maximal number of iterations in Remez algorithm for computing a minimax approximation for univariate terms
    unsigned REMEZ_MAXIT;
    //! @brief Stopping tolerance in Remez algorithm for computing a minimax approximation for univariate terms
    double REMEZ_TOL;
    //! @brief Threshold for interval width below which Remez algorithm is not used
    double REMEZ_MIG;
    //! @brief Extra terms in chebyshev interpolation of univariates: 0-Chebyshev interpolation of order MAXORD; extra terms allow approximation of Chebyshev truncated series
    unsigned INTERP_EXTRA;
    //! @brief Threshold for coefficient values in Chebyshev interpolation of univariates
    double INTERP_THRES;
    //! @brief Chebyshev model range bounder
    BOUNDER BOUNDER_TYPE;
    //! @brief Array of Chebyshev model range bounder names (for display)
    static const std::string BOUNDER_NAME[2];
    //! @brief Whether to intersect internal bounds with underlying bounds in T arithmetics
    bool MIXED_IA;
    //! @brief Threshold for monomial coefficients below which the term is removed and appended to the remainder term.
    double MIN_FACTOR;
    //! @brief Scalar in \f$[0,1]\f$ related to the choice of the polynomial part in the overloaded functions mc::inter and mc::hull. A value of 0 amounts to selecting the polynomial part of the left operand, whereas a value of 1 selects the right operand.
    double REF_POLY;
    //! @brief Number of digits in output stream for Chebyshev model coefficients.
    unsigned DISPLAY_DIGITS;
  } options;
  
  //!brief Unit ball in T arithmetic
  static T TOne;
  
  //!brief Zero-one in T arithmetic
  static T TZerOne;

  /** @} */

private:

  //! @brief Interval coefficients in Chebyshev interpolant of univariate functions
  std::vector<T> _coefuniv;

  //! @brief Coefficients in Chebyshev interpolant of univariate functions
  std::vector<double> _ainterp;

  //! @brief Points in Chebyshev interpolant of univariate functions
  std::vector<double> _xinterp;
      
  //! @brief Function values in Chebyshev interpolant of univariate functions
  std::vector<double> _finterp;

  //! @brief Get Chebyshev basis functions in U arithmetic for variable <a>bndvar</a>
  template <typename U>
  U* _get_bndpow
    ( const unsigned maxord, U const& bndvar, double const ref, double const scal )
    const;

  //! @brief Get Chebyshev basis functions in U arithmetic for [-1,1] scaled variable <a>bndvar</a>
  template <typename U>
  U* _get_bndpow
    ( const unsigned maxord, U const& bndvar )
    const;

  //! @brief Conversion from power series to Chebyshev series
  template <typename VEC>
  VEC _to_chebyshev
    ( VEC const& veccoef )
    const;

  //! @brief Construct minimax polynomial coefficient for univariate <a>f</a> in <a>_coefuniv</a> and return corresponding error for univariate <a>f</a> using Remez algorithm
  template <typename PUNIV>
  double _minimax
    ( PUNIV const& f, unsigned const maxord );

  //! @brief Compose minimax polynomial for univariate <a>f</a> with Chebyshev variable <a>CV</a> and return resulting Chebyshev variable in <a>CV2</a>
  template <typename PUNIV>
  bool _minimax
    ( PUNIV const& f, SICVar<T> const& CV, SICVar<T>& CV2 );

  //! @brief Return error at bound between univariate <a>f</a> and corresponding interpolating coefficients in <a>_coefuniv</a>
  template <typename PUNIV>
  double _rematbound
    ( PUNIV const& f );

  //! @brief Construct interpolating polynomial coefficient for univariate <a>f</a> in <a>_coefuniv</a> up to order <a>maxord</a>
  template <typename PUNIV>
  void _chebinterp
    ( PUNIV const& f, unsigned const maxord );

  //! @brief Construct interpolating polynomial coefficient for univariate <a>f</a> in <a>_coefuniv</a> until coefficient magnitude is less than <a>tol</a>
  template <typename PUNIV>
  void _chebinterp
    ( PUNIV const& f, double const tol, unsigned& maxord );

  //! @brief Apply the Clenshaw algorithm
  SICVar<T> _clenshaw
    ( SICVar<T> const& CVinner, unsigned const maxord );

  //! @brief Apply the Horner algorithm
  SICVar<T> _horner
    ( SICVar<T> const& CVinner, unsigned const maxord );
    
  //! @brief Compose interpolating polynomial for univariate <a>f</a> with Chebyshev variable <a>CV</a> and return resulting Chebyshev variable in <a>CV2</a>
  template <typename PUNIV>
  bool _chebinterp
    ( PUNIV const& f, SICVar<T> const& CV, SICVar<T>& CV2, bool const rematbound );

  //! @brief Apply Chebyshev composition to variable <a>CVI</a> using the coefficients <a>coefmon</a> of the outer function
  template <typename U> static SICVar<T> _composition
    ( std::vector<U> const& coefouter, unsigned const maxord,
      SICVar<T> const& CVinner );

  //! @brief Recursive calculation of nonnegative integer powers
  SICVar<T> _intpow
    ( SICVar<T> const& CV, int const n )
    const;

  //! @brief Polynomial range bounder - Naive approach
  template <typename C, typename U> U _polybound_naive
    ( std::map<SPolyMon,C,lt_SPolyMon> const& coefmon, U const* const* bndbasis,
      unsigned const minord=0 )
    const;

  //! @brief Polynomial range bounder - Lin & Stadtherr approach
  template <typename C, typename U> U _polybound_LSB
    ( std::map<SPolyMon,C,lt_SPolyMon> const& coefmon, U const* const* bndbasis )
  const;

  //! @brief Polynomial range bounder using specified bounder <a>type</a> with basis functions <a>bndbasis</a> in U arithmetic and monomial coefficients <a>coefmon</a> in C arithmetic
  template <typename C, typename U> U _polybound
    ( std::map<SPolyMon,C,lt_SPolyMon> const& coefmon, U const *const* bndbasis,
      int const type )
    const;

//  //! @brief Polynomial product with range - Naive approach
//  std::map<SPolyMon,T,lt_SPolyMon>& _polyprod_naive
//    ( std::map<SPolyMon,T,lt_SPolyMon>& coefmon, T const& bndvar,
//      unsigned const minord=0 )
//    const;

//  //! @brief Polynomial product with range - Lin & Stadtherr approach
//  std::map<SPolyMon,T,lt_SPolyMon>& _polyprod_LSB
//    ( std::map<SPolyMon,T,lt_SPolyMon>& coefmon, T const& bndvar )
//    const;

//  //! @brief Polynomial product with range using specified bounder <a>type</a> with monomial coefficients <a>coefmon</a> in T arithmetic
//  std::map<SPolyMon,T,lt_SPolyMon>& _polyprod
//    ( std::map<SPolyMon,T,lt_SPolyMon>& coefmon, T const& bndvar,
//      int const type )
//    const;

  //! @brief Recursive product of univariate Chebyshev polynomials
  void _sprod1D
    ( std::map<unsigned,t_coefmon> const& sp1map,
      std::map<unsigned,t_coefmon> const& sp2map,
      t_coefmon& coefmon, std::set<unsigned> const& ndxvar,
      typename t_var::const_iterator itvar ) const;

  //! @brief Scaling of Chebyshev coefficient maps
  void _sscal1D
    ( const t_coefmon& coefmon0, T const&coefscal,
      t_coefmon& coefmon ) const;

  //! @brief Lifting of Chebyshev coefficient maps
  void _slift1D
    ( const t_coefmon& coefmon0, double const& dscal,
      t_coefmon& coefmon ) const;

  //! @brief Lifting of Chebyshev coefficient maps
  void _slift1D
    ( t_coefmon const& coefmon0, double const& dscal,
      t_coefmon& coefmon, typename t_var::const_iterator itvar,
      unsigned const ndxord ) const;

  //! @brief Build 1D vector of Chebyshev coefficients
  void _svec1D
    ( typename t_var::const_iterator itvar,
      std::pair<SPolyMon,T> const& mon,
      std::map<unsigned,t_coefmon>& mapspoly ) const;

  //! @brief Display of recursive univariate Chebyshev polynomials
  void _sdisp1D
    ( const std::vector<SICVar<T>>& coefmon, const unsigned ndxvar,
      const std::string&name="", std::ostream&os=std::cout ) const;

  //! @brief Display of recursive univariate Chebyshev polynomials
  void _sdisp1D
    ( std::map<unsigned,t_coefmon> const& coefmon, typename t_var::const_iterator itvar,
      const std::string&name="", std::ostream&os=std::cout ) const;

  //! @brief Display of recursive univariate Chebyshev polynomials
  void _sdisp1D
    ( const t_coefmon& coefmon, const std::string&name="",
      std::ostream&os=std::cout ) const;

  //! @brief Build 1D vector of Chebyshev coefficients
  void _svec1Dfull
    ( unsigned const ndxvar, std::pair<SPolyMon,T> const& coefmon,
      std::vector<SICVar<T>>& vec ) const;

  //! @brief Product of multivariate Chebyshev polynomials in sparse format
  void _sprod
    ( SICVar<T> const& CV1, SICVar<T> const& CV2, t_coefmon& coefmon,
      double&coefrem ) const;

  //! @brief Squaring of multivariate Chebyshev polynomials in sparse format
  void _ssqr
    ( SICVar<T> const& CV, t_coefmon& coefmon, double&coefrem ) const;
};

template <typename T>
inline
const std::string SICModel<T>::Options::BOUNDER_NAME[2]
  = { "NAIVE", "LSB" };

template <typename T>
inline
T SICModel<T>::TZerOne
  = Op<T>::zeroone();

template <typename T>
inline
T SICModel<T>::TOne
  = 2.*Op<T>::zeroone()-1.;

//! @brief C++ class for the computation of sparse interval Chebyshev models for factorable function: propagation
////////////////////////////////////////////////////////////////////////
//! mc::SICVar is a C++ class for propagation of sparse interval
//! Chebyshev models through factorable functions. Variables of type
//! mc::SICVar are registered in a Chebyshev model environment of class
//! mc::SCModel. The template parameter corresponds to the type used in
//! computing the interval coefficients.
////////////////////////////////////////////////////////////////////////
template <typename T>
class SICVar
////////////////////////////////////////////////////////////////////////
{
  template <typename U> friend class SICVar;
  template <typename U> friend class SICModel;

  template <typename U> friend SICVar<U> operator-
    ( const SICVar<U>& );
  template <typename U> friend SICVar<U> operator*
    ( const SICVar<U>&, const SICVar<U>& );
  template <typename U> friend std::ostream& operator<<
    ( std::ostream&, const SICVar<U>& );

  template <typename U> friend U funcptr
    ( const U, const unsigned  );
  template <typename U> friend void interpolation
    ( double*, const SICVar<U>&, const unsigned  );
  template <typename U> friend SICVar<U> composition
    ( double const*, const SICVar<U>& );

  template <typename U> friend SICVar<U> inv
    ( const SICVar<U>& );
  template <typename U> friend SICVar<U> sqr
    ( const SICVar<U>& );
  template <typename U> friend SICVar<U> sqrt
    ( const SICVar<U>& );
  template <typename U> friend SICVar<U> exp
    ( const SICVar<U>& );
  template <typename U> friend SICVar<U> log
    ( const SICVar<U>& );
  template <typename U> friend SICVar<U> xlog
    ( const SICVar<U>& );
  template <typename U> friend SICVar<U> pow
    ( const SICVar<U>&, const int );
  template <typename U> friend SICVar<U> pow
    ( const SICVar<U>&, double const );
  template <typename U> friend SICVar<U> pow
    ( double const, const SICVar<U>& );
  template <typename U> friend SICVar<U> pow
    ( const SICVar<U>&, const SICVar<U>& );
  template <typename U> friend SICVar<U> prod
    ( const unsigned int, const SICVar<U>* );
  template <typename U> friend SICVar<U> monom
    ( const unsigned int, const SICVar<U>*, const unsigned* );
  template <typename U> friend SICVar<U> cheb
    ( const SICVar<U>&, const unsigned );
  template <typename U> friend SICVar<U> cos
    ( const SICVar<U>& );
  template <typename U> friend SICVar<U> sin
    ( const SICVar<U>& );
  template <typename U> friend SICVar<U> tan
    ( const SICVar<U>& );
  template <typename U> friend SICVar<U> acos
    ( const SICVar<U>& );
  template <typename U> friend SICVar<U> asin
    ( const SICVar<U>& );
  template <typename U> friend SICVar<U> atan
    ( const SICVar<U>& );
  template <typename U> friend SICVar<U> cosh
    ( const SICVar<U>& );
  template <typename U> friend SICVar<U> sinh
    ( const SICVar<U>& );
  template <typename U> friend SICVar<U> tanh
    ( const SICVar<U>& );
  template <typename U> friend SICVar<U> fabs
    ( const SICVar<U>& );
  template <typename U> friend SICVar<U> hull
    ( const SICVar<U>&, const SICVar<U>& );
  template <typename U> friend bool inter
    ( SICVar<U>&, const SICVar<U>&, const SICVar<U>& );

public:
  // Monomial representation: <total order, <variable index, order>>
  typedef std::map< SPolyMon, T, lt_SPolyMon > t_coefmon;
  typedef std::set< unsigned > t_var;
  
private:
  //! @brief Pointer to Chebyshev model environment
  SICModel<T> *_CM;

  //! @brief Set of participating variables in polynomial expression
  t_var _ndxvar;

  //! @brief Array of size <tt>_nmon</tt> with monomial coefficients of variable
  t_coefmon _coefmon;

  //! @brief Pointer to variable bound in underlying T arithmetic (possibly NULL if not computed)
  mutable T* _bndT;

  //! @brief Pointer to polynomial bound of variable (possibly NULL if not available)
  mutable T* _bndpol;

  //! @brief Initialize private/protected members of model variable
  void _init
    ();

  //! @brief Reinitialize private/protected members of model variable
  void _reinit
    ();

  //! @brief Clean up private/protected members of variable
  void _cleanup
    ();

  //! @brief Set variable bound in unerlying T arithmetic
  void _set_bndT
    ( T const& bndT );

  //! @brief Set variable bound in unerlying T arithmetic
  void _set_bndT
    ( T const* bndT );

  //! @brief Unset variable bound in underlying T arithmetic
  void _unset_bndT
    ();

  //! @brief Set polynomial bound in variable as <tt>bndpol</tt>
  void _set_bndpol
    ( T const& bndpol );

  //! @brief Set polynomial bound in variable as <tt>bndpol</tt>
  void _set_bndpol
    ( T const* bndpol );

  //! @brief Unset polynomial bound in variable
  void _unset_bndpol
    ();

  //! @brief Original bounds on variable <tt>ivar</tt>
  T const&  _bndvar
    ( const unsigned ivar )
    const
    { return _CM->_bndvar[ivar]; };

  //! @brief Reference point for variable <tt>ivar</tt> in Chebyshev model
  double _refvar
    ( const unsigned ivar )
    const
    { return _CM->_refvar[ivar]; };

  //! @brief Scaling for variable <tt>ivar</tt> in Cheyshev model
  double _scalvar
    ( const unsigned ivar ) const
    { return _CM->_scalvar[ivar]; };

  //! @brief Set Chebyshev variable equal to <a>CV</a>
  SICVar<T>& _set
    ( SICVar<T> const& CV );

  //! @brief Set Chebyshev variable equal to <a>CV</a>
  SICVar<T>& _set
    ( SICVar<T> && CV );

  //! @brief Set Chebyshev variable with index <a>ix</a> (starting from 0) and bounded by <a>X</a>
  SICVar<T>& _set
    ( const unsigned ivar, T const& X, const bool updMod=true );

public:
  /** @addtogroup SICHEBYSHEV Chebyshev Model Arithmetic for Factorable Functions
   *  @{
   */
  //! @brief Display sparse polynomial
  std::string display
    ( typename SICVar::t_coefmon coefmon, int const& BASIS=1, int const& IDISP=5 )
    const;
    
  //! @brief Get pointer to linked Chebyshev model environment
  SICModel<T>* env() const
    { return _CM; }

  //! @brief Max order of model environment
  unsigned maxord
    ()
    const
    { return( _CM? _CM->_maxord: 0 ); };

  //! @brief Max number of variables in model environment
  unsigned maxvar
    ()
    const
    { return( _CM? _CM->_maxvar: 0 ); };

  //! @brief Max order of participating variables in polynomial variable
  unsigned nord
    ()
    const
    { //std::cout << "NORD = " << _coefmon.crbegin()->first.tord << std::endl; 
      return _coefmon.crbegin()->first.tord; }

  //! @brief Number of participating variables in polynomial variable
  unsigned nvar
    ()
    const
    { return _ndxvar.size(); }

  //! @brief Total number of monomial terms in polynomial variable
  unsigned nmon
    ()
    const
    { return _coefmon.size(); }

  //! @brief Get const map of monomial coefficients
  t_coefmon const& coefmon
    ()
    const
    { return _coefmon; }

  //! @brief Get map of monomial coefficients
  t_coefmon& coefmon
    ()
    { return _coefmon; }

  //! @brief Get const set of participating variables
  t_var const& ndxvar
    ()
    const
    { return _ndxvar; }

  //! @brief Get set of participating variables
  t_var& ndxvar
    ()
    { return _ndxvar; }

  //! @brief Constructor of Chebyshev variable for a real scalar
  SICVar
    ( SICModel<T>*CM );

  //! @brief Constructor of Chebyshev variable for a real scalar
  SICVar
    ( double const d=0., SICModel<T>*CM=0 );

  //! @brief Constructor of Chebyshev variable for a remainder bound
  SICVar
    ( T const& B, SICModel<T>*CM=0 );

  //! @brief Constructor of Chebyshev variable with index <a>ix</a> (starting from 0) and bounded by <a>X</a>
  SICVar
    ( SICModel<T>*CM, const unsigned ix, T const& X );

  //! @brief Copy constructor of Chebyshev variable
  SICVar
    ( SICVar<T> const& CV );

  //! @brief Copy constructor of Chebyshev variable
  SICVar
    ( SICVar<T> && CV );

  //! @brief Destructor of Chebyshev variable
  ~SICVar
    ()
    { _unset_bndpol(); _unset_bndT(); }

  //! @brief Set Chebyshev variable with index <a>ix</a> (starting from 0) and bounded by <a>X</a>
  SICVar<T>& set
    ( SICModel<T>*CM, const unsigned ix, T const& X )
    { set( CM ); _set( ix, X ); return *this; }

  //! @brief Set environment in Chebyshev model
  SICVar<T>& set
    ( SICModel<T>*CM )
    { _CM = CM; return *this; }
    
  //! @brief Set multivariate polynomial coefficients in variable as <tt>coefmon</tt>
  SICVar<T>& set
    ( t_coefmon const& coefmon, t_var const& ndxvar )
    { _coefmon = coefmon; _ndxvar = ndxvar;
      _unset_bndT(); _unset_bndpol();
      return *this; } // this is assuming the same order and number of variables
    
  //! @brief Set multivariate polynomial coefficients in variable as <tt>coefmon</tt>
  SICVar<T>& set
    ( std::map< SPolyMon, double, lt_SPolyMon > const& coefmon, t_var const& ndxvar,
      T const& bndrem, T const* bndT=nullptr )
    { _coefmon.clear();
      for( auto [mon,coef] : coefmon ) _coefmon.insert( std::make_pair(mon,coef) ); 
      _ndxvar = ndxvar;
      _unset_bndpol();
      *this += bndrem;
      if( bndT ) _set_bndT( bndT );
      else       _unset_bndT();
      return *this; } // this is assuming the same order and number of variables
#if 0
  //! @brief Set multivariate polynomial coefficients and remainder term equal to those in variable <tt>var</tt>, possibly defined in another polynomial model environment with fewer variables or with a different expansion order. Coefficients involving other variables or higher order are initialized to 0 if <tt>reset=true</tt> (default), otherwise they are left unmodified. Higher-order terms in TV are bounded and added to the remainder bound.
  SICVar<T>& set
    ( SICVar<T> const& var, bool const reset=true );

  //! @brief Copy multivariate polynomial coefficients from current variable into variable <tt>var</tt>, possibly defined in another polynomial model environment with less variables or with a lower expansion order. Copied coefficients are reset to 0 in current Taylor variable if <tt>reset=true</tt>, otherwise they are left unmodified (default).
  SICVar<T>& get
    ( SICVar<T>&var, const bool reset=false );
#endif
  //! @brief Compute bound on variable using bounder <a>type</a>
  T bound
    ( int const type )
    const
    { if( !_bndT ) return _polybound(type);
      else{ T bndT; return Op<T>::inter( bndT, _polybound(type), *_bndT )?
                           bndT: _polybound(type); } }

  //! @brief Retreive bound on variable using default bounder
  T bound
    ()
    const
    { if( !_bndT ) return bndpol();
      else{ T bndT; return Op<T>::inter( bndT, bndpol(), *_bndT )? bndT: bndpol(); } }

  //! @brief Retreive bound on variable using default bounder in U arithmetic
  template <typename U>
  U bound
    ( const U*const*bndbasis )
    const
    { return _polybound( bndbasis ); }

  //! @brief Retreive bound on variable using bounder <a>type</a> in U arithmetic
  template <typename U>
  U bound
    ( const U*const*bndbasis, const int type ) const
    { return _polybound( bndbasis, type ); }

  //! @brief Retreive bound on multivariate polynomial using default bounder
  T bndpol
    ()
    const
    { if( !_bndpol ) _bndpol = new T( _polybound() );
      return *_bndpol; }

  //! @brief Retreive bound on all terms with (total) order <tt>minord</tt> or higher in polynomial part
  T bndord
    ( const unsigned minord )
    const
    { if( !_CM ) return !minord && !_coefmon.empty() && !_coefmon.begin()->first.tord? _coefmon.begin()->second: 0.;
      return _CM->_polybound_naive( _coefmon, (T const* const*)0, minord ); }

  //! @brief Shortcut to mc::SICVar::bound
  T B
    ( const int type )
    const
    { return bound( type ); }

  //! @brief Shortcut to mc::SICVar::bound
  T B
    ()
    const
    { return bound(); }

  //! @brief Evaluate mid-point polynomial part at <tt>x</tt>
  double P
    ( double const*x )
    const;

  //! @brief Evaluate interval polynomial part at <tt>x</tt>
  T IP
    ( double const*x )
    const;
#if 0
  //! @brief Return new Chebyshev variable with same multivariate polynomial part but zero remainder
  SICVar<T> polynomial
    ()
    const
    { SICVar<T> var = *this; var._bndrem = 0.; return var; }

  //! @brief Shortcut to mc::SICVar::polynomial
  SICVar<T> P
    ()
    const
    { return polynomial(); }
#endif
  //! @brief Get coefficient of constant term in Chebyshev variable. The value of this coefficient is reset to 0 if <tt>reset=true</tt>, otherwise it is left unmodified (default).
  T constant
    ( bool const reset=false );

  //! @brief Get coefficients of linear term for variable <tt>ivar</tt> in Chebyshev variable. The value of this coefficient is reset to 0 if <tt>reset=true</tt>, otherwise it is left unmodified (default).
  T linear
    ( unsigned const ivar, bool const reset=false );

  //! @brief Simplify the model by appending coefficient less than TOL as well as terms with order greater than or equal to ORD to the remainder term
  SICVar<T>& simplify
    ( int const& BASIS=1, double const& TOL=0e0, int const& TORD=-1 );
    
  //! @brief Scale coefficients in Chebyshev variable for the modified range <a>X</a> of variable i
  SICVar<T>& scale
    ( unsigned const i, T const& X );
    
  //! @brief Scale coefficients in Chebyshev variable for the modified variable ranges <a>X</a>
  SICVar<T>& scale
    ( T const* X );
    
  //! @brief Rescale coefficients in Chebyshev variable for their original variable ranges
  SICVar<T>::t_coefmon unscale
    ()
    const;
        
  //! @brief Return new coefficient map in monomial basis representation
  t_coefmon to_monomial
    ( bool const scaled=false )
    const;
    
  //! @brief Return new coefficient map in monomial basis representation after removing terms with coefficient less than TOL or order greater than or equal to ORD, and also return a bound on the removed terms
  std::pair<t_coefmon,T> to_monomial
    ( bool const scaled, double const TOL, int const TORD=-1 )
    const;
 /** @} */

  SICVar<T>& operator=
    ( SICVar<T> const& );
  SICVar<T>& operator=
    ( SICVar<T> && );
  SICVar<T>& operator=
    ( double const );
  SICVar<T>& operator=
    ( T const&  );
  template <typename U> SICVar<T>& operator+=
    ( SICVar<U> const& );
  template <typename U> SICVar<T>& operator+=
    ( U const&  );
  SICVar<T>& operator+=
    ( double const );
  template <typename U> SICVar<T>& operator-=
    ( const SICVar<U>& );
  template <typename U> SICVar<T>& operator-=
    ( U const&  );
  SICVar<T>& operator-=
    ( double const );
  SICVar<T>& operator*=
    ( const SICVar<T>& );
  SICVar<T>& operator*=
    ( double const );
  SICVar<T>& operator*=
    ( T const&  );
  SICVar<T>& operator/=
    ( SICVar<T> const& );
  SICVar<T>& operator/=
    ( double const );

private:
  //! @brief Polynomial range bounder using specified bounder <a>type</a> with basis functions <a>bndbasis</a> in U arithmetic
  template <typename U>
  U _polybound
    ( U const* const* bndbasis, int const type ) const
    { if( !_CM ) return !_coefmon.empty() && !_coefmon.begin()->first.tord? _coefmon.begin()->second: 0.;
      return _CM->_polybound( _coefmon, bndbasis, type ); }

  //! @brief Polynomial range bounder using default bounder in U arithmetic
  template <typename U>
  U _polybound
    ( U const* const* bndbasis ) const
    { return _polybound( bndbasis, _CM? _CM->options.BOUNDER_TYPE: 0 ); }

  //! @brief Polynomial range bounder using specified bounder <a>type</a>
  T _polybound
    ( int const type ) const
    { return _polybound( (const T*const*)nullptr, type ); }

  //! @brief Polynomial range bounder using default bounder
  T _polybound
    () const
    { return _polybound( _CM? _CM->options.BOUNDER_TYPE: 0 ); }

  //! @brief Scale current variable in order for its range to be within [-1,1], with <a>c</a> and <a>w</a> respectively the center and width, respectively, of the orginal variable range
  SICVar<T> _rescale
    ( double const w, double const c ) const
    { return( !isequal(w,0.)? (*this-c)/w: c ); }

  //! @brief Scale variable <a>X</a>
  template <typename U> static
  U _rescale
    ( U& X, double const w, double const c )
    { return( !isequal(w,0.)? (X-c)/w: c ); }

  //! @brief Return an array of Chebyshev variables representing the coefficients in the univariate polynomial for variable <a>ivar</a> only
  SICVar<T>* _single
    ( unsigned const ivar );

  //! @brief Cancel zero entries in coefficient map <a>coefmon</a>
  void _simplify
    ( t_coefmon& coefmon )
    const;

  //! @brief Simplify the coefficient map <a>coefmon</a> by removing entries with coefficient less than TOL or order greater than or equal to ORD, and return a bound on the removed entries
  T _simplify_monomial
    ( typename SICVar<T>::t_coefmon& coefmon, bool const scaled,
      double const& TOL=0e0, int const& TORD=-1 )
    const;

  //! @brief Scale coefficient map <a>coefmon</a> for the modified range <a>X</a> of variable i
  void _scale
    ( unsigned const i, T const& X, typename SICVar<T>::t_coefmon& coefmon )
    const;

  //! @brief Convert Chebyshev basis into monomial basis for variable <a>ivar</a> in coefficient map <a>coefmon</a>
  void _to_monomial
    ( unsigned const ivar, typename SICVar<T>::t_coefmon& coefmon )
    const;
};

////////////////////////////////// SICModel //////////////////////////////////////

template <typename T>
template <typename U>
inline U*
SICModel<T>::_get_bndpow
( unsigned const maxord, U const& bndvar, double const ref, double const scal )
const
{
  U *bndmon = new U[maxord+1];
  bndmon[0] = 1;
  if( maxord ) bndmon[1] = ( bndvar - ref ) / scal;
  for( unsigned i=2; i<=maxord; i++ )
    switch( options.BASIS ){
      case Options::POW:  bndmon[i] = Op<U>::pow( bndmon[1], i );  break;
      case Options::CHEB: bndmon[i] = Op<U>::cheb( bndmon[1], i ); break;
    }
 return bndmon;
}

template <typename T>
template <typename U>
inline U*
SICModel<T>::_get_bndpow
( unsigned const maxord, U const& bndvar )
const
{
  U* bndmon = new U[maxord+1];
  bndmon[0] = 1;
  if( maxord ) bndmon[1] = bndvar;
  for( unsigned i=2; i<=maxord; i++ )
    switch( options.BASIS ){
      case Options::POW:  bndmon[i] = Op<U>::pow( bndmon[1], i );  break;
      case Options::CHEB: bndmon[i] = Op<U>::cheb( bndmon[1], i ); break;
    }
 return bndmon;
}

template <typename T>
template <typename U>
inline U**
SICModel<T>::get_basis
( unsigned const maxord, U const* bndvar, bool const scaled )
const
{
  assert ( bndvar );
  U** bndbasis = new U*[_maxvar];
  for( unsigned i=0; i<_maxvar; i++){
    if( _setvar.find(i) == _setvar.end() ){
      bndbasis[i] = nullptr;
      continue;
    }
    bndbasis[i] = scaled? _get_bndpow( maxord, bndvar[i] ):
                          _get_bndpow( maxord, bndvar[i], _refvar[i], _scalvar[i] );
  }
  return bndbasis;
}

template <typename T>
template <typename U>
inline void
SICModel<T>::get_bndmon
( std::map<SPolyMon,U,lt_SPolyMon>& bndmon, U const* bndvar,
  bool const scaled )
const
{
  if( bndmon.empty() ) return;

  unsigned const nord = bndmon.rbegin()->first.tord;
  U** bndbasis = get_basis( nord, bndvar, scaled );

#ifdef MC__SICMODEL_USE_PROD
  std::vector<U> Umon; Umon.reserve(_maxvar);
  auto it = bndmon.begin();
  if( !it->first.tord ){ it->second = 1; ++it; }
  for( ; it!=bndmon.end(); ++it ){
    Umon.clear();
    for( auto& [ivar,iord] : it->first.second )
      Umon.push_back( bndbasis[ivar][iord] );
    it->second = prod( Umon.size(), Umon.data() );
  }
#else
  auto it = bndmon.begin();
  if( !it->first.tord ){ it->second = 1; ++it; }
  for( ; it!=bndmon.end(); ++it ){
    it->second = 1;
    for( auto& [ivar,iord] : it->first.second )
      it->second *= bndbasis[ivar][iord];
  }
#endif

  for( unsigned j : _setvar ) delete[] bndbasis[j];
  delete[] bndbasis;
}

template <typename T>
template <typename PUNIV>
inline double
SICModel<T>::_minimax
( PUNIV const& f, unsigned const maxord ) // range of f assumed as [-1,1]
{
  boost::math::tools::remez_minimax<double> problem
    ( f, maxord, 0, -1., 1., false, false, 0, 64 );

  for( unsigned iter=0; iter<options.REMEZ_MAXIT; ++iter ){
    problem.iterate();
#ifdef MC__SICMODEL_DEBUG_MINIMAX
    const boost::math::tools::polynomial<double> a = problem.numerator();
    std::cout << iter << ": [ " << std::right << std::scientific << std::setprecision(15);
    for( unsigned k=0; k<a.size(); ++k ) std::cout << std::setw(23) << a[k];
    std::cout << " ] +/- " << std::setw(23) << problem.max_error() << std::endl;
#endif
    if( problem.max_change() < options.REMEZ_TOL ) break;
  }

  switch( options.BASIS ){
  case Options::CHEB:
    _ainterp = _to_chebyshev(problem.numerator().data());
    break;
  case Options::POW:
    _ainterp = problem.numerator().data();
    break;
  }
  return problem.max_error();
}

template <typename T>
template <typename VEC>
inline VEC
SICModel<T>::_to_chebyshev
( VEC const& veccoef )
const
{
  // Converts COEFF in monomial basis to the Chebyshev basis.
  VEC vecnew = veccoef;
  int const N = vecnew.size()-1;
  if( N <= 1 ) return vecnew;
  // TP = .5D0**(N-1)
  // COEFF(N) = TP * COEFF(N)
  // COEFF(N-1) = TP * COEFF(N-1)
  double TP = std::pow(0.5,N-1);
  vecnew[N]   *= TP;
  vecnew[N-1] *= TP;
  // do 20 J = N-2, 0, -1
  //   TP = 2.D0 * TP
  //   COEFF(J) = TP * COEFF(J)
  //   COEFF(J+1) = 2.D0 * COEFF(J+1)
  for( int J=N-2; J>=0; --J ){
    TP *= 2;
    vecnew[J]   *= TP;
    vecnew[J+1] *= 2;
    // do 10 I = J, N-2
    //   COEFF(I) = COEFF(I) + COEFF(I+2)
    for( int I=J; I<N-1; ++I )
      vecnew[I] += vecnew[I+2];
  }
  return vecnew;
}

template <typename T>
template <typename PUNIV>
inline bool
SICModel<T>::_minimax
( PUNIV const& f, SICVar<T> const& CV, SICVar<T>& CV2 )
{
  try{
    // Minimax interpolant
    double m( Op<SICVar<T>>::mid(CV) ), r( 0.5*Op<SICVar<T>>::diam(CV) );
    auto const& fscaled = [=]( const double& x ){ return f( r*x + m ); };
    unsigned const nord = _maxord + options.INTERP_EXTRA;
    double rem = _minimax( fscaled, nord );
    assert( _ainterp.size() == nord+1 );

    // Order reduction
    _coefuniv.resize( _maxord+1 );
    for( unsigned i=0; i<=_maxord; i++ )
      _coefuniv[i] = _ainterp[i];
    for( unsigned i=nord; i>_maxord; i-- ){
      unsigned k=0, d=1;
      for( ; (i-_maxord-1)/(d*2); ++k, d*=2 ){}
      if( k ){
        _coefuniv[k] += 2.*TOne*_ainterp[i];
        _ainterp[i-d] -= _ainterp[i];
      }
      else
        _coefuniv[0] += TOne*_ainterp[i];
    }

    // Composition
    SICVar<T> CVI( this );
    CVI = CV._rescale(r,m);
    switch( options.BASIS ){
    case Options::CHEB:
      CV2 = _clenshaw( CVI, _maxord ) + TOne * rem;
      break;
    case Options::POW:
      CV2 = _horner( CVI, _maxord ) + TOne * rem;
      break;
    }
    //CV2._ndxvar = CV._ndxvar;
    if( options.MIN_FACTOR >= 0. ) CV2.simplify( options.BASIS, options.MIN_FACTOR );
  }
  catch(...){
    return false;
  }
  return true;
}

template <typename T>
template <typename PUNIV>
inline double
SICModel<T>::_rematbound
( PUNIV const& f ) // range of f assumed as [-1,1]
{
  double ub(0), lb(0);
  for( unsigned i=0; i<=_maxord; i++ ){
    ub += _ainterp[i];
    lb += i%2? -_ainterp[i]: _ainterp[i];
  }
  return std::max( std::fabs( f(1.) - ub ), std::fabs( f(-1.) - lb ) );
}

template <typename T>
template <typename PUNIV>
inline void
SICModel<T>::_chebinterp
( PUNIV const& f, unsigned const maxord ) // range of f assumed as [-1,1]
{
  _xinterp.resize( maxord+1 );
  _finterp.resize( maxord+1 );
  double mulconst( PI/(2.*double(maxord+1)) );
  for( unsigned i=0; i<=maxord; i++ ){
    _xinterp[i]  = std::cos( mulconst*(2.*double(i)+1.) );
    _finterp[i] = f( _xinterp[i] );
  }

  _ainterp.resize( maxord+1 );
  switch( maxord ){
  case 0:
    _ainterp[0] = _finterp[0];
    return;
  case 1:
    _ainterp[0] = 0.5 * ( _finterp[0] + _finterp[1] );
    _ainterp[1] = ( _finterp[1] - _finterp[0] ) / ( _xinterp[1] - _xinterp[0] );
    return;
  default:
    for( unsigned i=0; i<=maxord; i++ ){
      double mulconst2( std::cos( mulconst*double(i) ) ),
             mulconst3( 4*mulconst2*mulconst2-2 ), 
             b0( 0 ), b1( 0 );
      b0 = _finterp[maxord];
      b1 = _finterp[maxord-1] + mulconst3*b0;
      for( unsigned j=maxord-2; j>1; j-=2 ){
        b0 = _finterp[j] + mulconst3*b1 - b0;
        b1 = _finterp[j-1] + mulconst3*b0 - b1;
      }
      if( !(maxord%2) )
        b0 = _finterp[0] + mulconst3*b1 - b0 - b1;
      else{
        b0 = _finterp[1] + mulconst3*b1 - b0;
        b0 = _finterp[0] + mulconst3*b0 - b1 - b0;
      }
      _ainterp[i] = 2./double(maxord+1)*mulconst2*b0;
    }
    _ainterp[0] *= 0.5;
    return;
  }
}

template <typename T>
template <typename PUNIV>
inline void
SICModel<T>::_chebinterp
( PUNIV const& f, double const tol, unsigned& maxord ) // range of f assumed as [-1,1]
{
  _chebinterp( f, maxord );
  for( ; std::fabs(_ainterp[maxord])>tol || (maxord && std::fabs(_ainterp[maxord-1])>tol); ){
    maxord *= 2;
    _chebinterp( f, maxord );
#ifdef MC__SICMODEL_DEBUG_CHEBINTERP
    for( unsigned i=0; i<=maxord; i++ )
      std::cout << "a[" << i << "] = " << _ainterp[i] << std::endl;
    { int dum; std::cout << "ENTER <1> TO CONTINUE"; std::cin >> dum; }
#endif
  }
}

template <typename T>
inline SICVar<T>
SICModel<T>::_horner
( SICVar<T> const& CVinner, unsigned const maxord )
{
  assert( maxord < _coefuniv.size() );
  
  // composition based on Horner's scheme
  switch( maxord ){
  case 0:
    return _coefuniv[0];
    break;
  default:
    SICVar<T> CV2 = _coefuniv[maxord] * CVinner;
    for( unsigned ord=maxord-1; ord>0; --ord ){
      CV2 += _coefuniv[ord];
      CV2 *= CVinner;
    }
    CV2 += _coefuniv[0];
    return CV2;
  }
}

template <typename T>
inline SICVar<T>
SICModel<T>::_clenshaw
( SICVar<T> const& CVinner, unsigned const maxord )
{
  assert( maxord < _coefuniv.size() );
  
  if( !maxord )
    return _coefuniv[0];

  else if( maxord == 1 )
    return CVinner * _coefuniv[1] + _coefuniv[0];

  //composition based on http://en.wikipedia.org/wiki/Clenshaw_algorithm#Special_case_for_Chebyshev_series
  SICVar<T> CVinnerx2 = 2. * CVinner;
#ifdef MC__SICMODEL_DEBUG_COMPOSITION
  std::cout << "CVinner:" << CVinner;
  std::cout << "CVinnerx2:" << CVinnerx2;
#endif
  SICVar<T> CV1 = _coefuniv[maxord];
#ifdef MC__SICMODEL_DEBUG_COMPOSITION
  std::cout << "CV1:" << CV1;
#endif
  SICVar<T> CV2 = _coefuniv[maxord-1] + CVinnerx2 * CV1;
#ifdef MC__SICMODEL_DEBUG_COMPOSITION
  std::cout << "CV2:" << CV2;
#endif

  for( unsigned i=maxord-2; i>1; i-=2 ){
#ifdef MC__SICMODEL_DEBUG_COMPOSITION
    std::cout << "CVinnerx2 * CV2" << CVinnerx2 * CV2;
    std::cout << "_coefuniv[i]:" << _coefuniv[i];
    std::cout << "_coefuniv[i] + CVinnerx2 * CV2" << _coefuniv[i] + CVinnerx2 * CV2;
#endif
    CV1 = _coefuniv[i]   + CVinnerx2 * CV2 - CV1;
#ifdef MC__SICMODEL_DEBUG_COMPOSITION
    std::cout << "CV1:" << CV1;
#endif
    CV2 = _coefuniv[i-1] + CVinnerx2 * CV1 - CV2;
#ifdef MC__SICMODEL_DEBUG_COMPOSITION
    std::cout << "CV2:" << CV2;
#endif
  }

  if( !(maxord%2) ){
#ifdef MC__SICMODEL_DEBUG_COMPOSITION
    CV1 = _coefuniv[0] + CVinner * CV2 - CV1;
    std::cout << "CV1:" << CV1;
    return CV1;
#else
    return _coefuniv[0] + CVinner * CV2 - CV1;
#endif
  }

  CV1 = _coefuniv[1] + CVinnerx2 * CV2 - CV1;
#ifdef MC__SICMODEL_DEBUG_COMPOSITION
  std::cout << "CV1:" << CV1;
  std::cout << "CVinner * CV1:" << CVinner * CV1;
  CV2 = _coefuniv[0] + CVinner * CV1 - CV2;
  std::cout << "CV2:" << CV2;
  return CV2;
#else
  return _coefuniv[0] + CVinner * CV1 - CV2;
#endif
}

template <typename T>
template <typename PUNIV>
inline bool
SICModel<T>::_chebinterp
( PUNIV const& f, SICVar<T> const& CV, SICVar<T>& CV2, bool const rematbound )
{
  try{
    double m( Op<SICVar<T>>::mid(CV) ), r( 0.5*Op<SICVar<T>>::diam(CV) ), rem( 0. );
    auto const& fscaled = [=]( const double& x ){ return f( r*x + m ); };
    _coefuniv.resize( _maxord+1 );
    if( rematbound ){
      unsigned const nord = _maxord + options.INTERP_EXTRA;
      _chebinterp( fscaled, nord );
      for( unsigned i=0; i<=_maxord; i++ )
        _coefuniv[i] = _ainterp[i];
      rem = _rematbound( fscaled );
    }
    else{
      unsigned nord = _maxord+2;
      double tol = options.INTERP_THRES;
      _chebinterp( fscaled, tol, nord );
      for( unsigned i=0; i<=_maxord; i++ )
        _coefuniv[i] = _ainterp[i];
      rem = 2*tol;
      for( unsigned i=_maxord+1; i<=nord; i++ )
        rem += std::fabs( _ainterp[i] );
    }

    //SICVar<T> CVI( this );
    //CVI = CV._rescale( r, m );
    CV2 = _clenshaw( CV._rescale( r, m ), _maxord );
    CV2 += TOne * rem;
    //CV2._ndxvar = CV._ndxvar;
    if( options.MIN_FACTOR >= 0. ) CV2.simplify( options.BASIS, options.MIN_FACTOR );
  }
  catch(...){
    return false;
  }
  return true;
}

template <typename T>
template <typename U>
inline SICVar<T>
SICModel<T>::_composition
( std::vector<U> const& coefouter, unsigned const maxord, SICVar<T> const& CVinner )
{
  //composition based on http://en.wikipedia.org/wiki/Clenshaw_algorithm#Special_case_for_Chebyshev_series
  if( !maxord )
    return coefouter[0];

  else if( maxord == 1 )
    return CVinner * coefouter[1] + coefouter[0];

  SICVar<T> CVinnerx2 = 2. * CVinner;
#ifdef MC__SICMODEL_DEBUG_COMPOSITION
  std::cout << "CVinner:" << CVinner;
  std::cout << "CVinnerx2:" << CVinnerx2;
#endif
  SICVar<T> CV1 = coefouter[maxord];
#ifdef MC__SICMODEL_DEBUG_COMPOSITION
  std::cout << "CV1:" << CV1;
#endif
  SICVar<T> CV2 = coefouter[maxord-1] + CVinnerx2 * CV1;
#ifdef MC__SICMODEL_DEBUG_COMPOSITION
  std::cout << "CV2:" << CV2;
#endif

  for( unsigned i=maxord-2; i>1; i-=2 ){
#ifdef MC__SICMODEL_DEBUG_COMPOSITION
    std::cout << "CVinnerx2 * CV2" << CVinnerx2 * CV2;
    std::cout << "coefouter[i]:" << coefouter[i];
    std::cout << "coefouter[i] + CVinnerx2 * CV2" << coefouter[i] + CVinnerx2 * CV2;
#endif
    CV1 = coefouter[i]   + CVinnerx2 * CV2 - CV1;
#ifdef MC__SICMODEL_DEBUG_COMPOSITION
    std::cout << "CV1:" << CV1;
#endif
    CV2 = coefouter[i-1] + CVinnerx2 * CV1 - CV2;
#ifdef MC__SICMODEL_DEBUG_COMPOSITION
    std::cout << "CV2:" << CV2;
#endif
  }

  if( !(maxord%2) ){
#ifdef MC__SICMODEL_DEBUG_COMPOSITION
    CV1 = coefouter[0] + CVinner * CV2 - CV1;
    std::cout << "CV1:" << CV1;
    return CV1;
#else
    return coefouter[0] + CVinner * CV2 - CV1;
#endif
  }

  CV1 = coefouter[1] + CVinnerx2 * CV2 - CV1;
#ifdef MC__SICMODEL_DEBUG_COMPOSITION
  CV2 = coefouter[0] + CVinner * CV1 - CV2;
  std::cout << "CV2:" << CV2;
  return CV2;
#else
  return coefouter[0] + CVinner * CV1 - CV2;
#endif
}

template <typename T>
inline SICVar<T>
SICModel<T>::_intpow
( SICVar<T> const& CV, int const n ) const
{
  if( n == 0 ) return 1.;
  else if( n == 1 ) return CV;
  return n%2 ? sqr( _intpow( CV, n/2 ) ) * CV : sqr( _intpow( CV, n/2 ) );
}

template <typename T>
inline void
SICModel<T>::_sscal1D
( t_coefmon const& coefmon0, T const& Iscal, t_coefmon& coefmon )
const
{
  if( Op<T>::abs(Iscal) == 0. ) return;
//  double mscal = Op<T>::mid(Iscal);
//  T rngsum = 0.;
  coefmon = coefmon0;
//  if( isequal(dscal,1.) ) return;
//  switch( options.HOT_SPLIT ){ 
//  case Options::FULL:
//  case Options::SIMPLE:
    for( auto& [mon,coef] : coefmon )
      coef *= Iscal;
//    break;
//  case Options::NONE:
//    for( auto& [mon,coef] : coefmon ){
//      rngsum += coef * ( Iscal - mscal );
//      coef *= mscal;
//    }
//    auto [itmon,ins] = coefmon.insert( std::make_pair( SPolyMon(), rngsum ) );
//    if( !ins ) itmon->second += rngsum;
//#ifdef MC__SICMODEL_DEBUG_SSACL
//    std::cout << "_sscal1D: " << itmon->first.display(options.BASIS) << "  " << itmon->second << std::endl;
//#endif
//    break;
//  }
  return;
}

template <typename T>
inline void
SICModel<T>::_slift1D
( t_coefmon const& coefmon0, double const& dscal, t_coefmon& coefmon )
const
{
  if( isequal(dscal,0.) ) return;
  for( auto const& [mon0,coef0] : coefmon0 ){
    auto [itmon,ins] = coefmon.insert( std::make_pair(mon0,coef0) );
    if( ins ){
      if( isequal(dscal,1.) ) continue;
      itmon->second *= dscal;
    }
    else{
      itmon->second += coef0*dscal;
    }
  }
}

template <typename T>
inline void
SICModel<T>::_slift1D
( t_coefmon const& coefmon0, double const& dscal, t_coefmon& coefmon,
  typename t_var::const_iterator itvar, unsigned const ndxord )
const
{
  for( auto const& [mon0,coef0] : coefmon0 ){
    // If total order too large then allocate among monomials
    if( mon0.tord + ndxord > _maxord ){
      switch( options.BASIS ){
      case Options::CHEB:
        switch( options.HOT_SPLIT ){ 
        case Options::FULL:
        {
          T add2mon = (dscal/double(mon0.expr.size()+1))*coef0*TOne;
          for( auto const& [ivar,iord] : mon0.expr ){
            auto [itmon1,ins1] = coefmon.insert( std::make_pair( SPolyMon(ivar,iord), add2mon ) );
            if( !ins1 ) itmon1->second += add2mon;
#ifdef MC__SICMODEL_DEBUG_SLIFT
            std::cout << "_slift1D: " << itmon1->first.display(options.BASIS) << "  " << add2mon << std::endl;
#endif
          }
          if( ndxord <= _maxord ){
            auto [itmon2,ins2] = coefmon.insert( std::make_pair( SPolyMon(*itvar,ndxord), add2mon ) );
            if( !ins2 ) itmon2->second += add2mon;
#ifdef MC__SICMODEL_DEBUG_SLIFT
            std::cout << "_slift1D: " << itmon2->first.display(options.BASIS) << "  " << add2mon << std::endl;
#endif
          }
          else{
            auto [itmon2,ins2] = coefmon.insert( std::make_pair( SPolyMon(), add2mon ) );
            if( !ins2 ) itmon2->second += add2mon;
#ifdef MC__SICMODEL_DEBUG_SLIFT
            std::cout << "_slift1D: " << itmon2->first.display(options.BASIS) << "  " << add2mon << std::endl;
#endif
          }
          continue;
        }
        case Options::SIMPLE:
        {
          T add2mon = (0.5*dscal)*coef0*TOne;
          auto [itmon1,ins1] = coefmon.insert( std::make_pair( mon0, add2mon ) );
          if( !ins1 ) itmon1->second += add2mon;
#ifdef MC__SICMODEL_DEBUG_SLIFT
          std::cout << "_slift1D: " << itmon1->first.display(options.BASIS) << "  " << itmon1->second << std::endl;
#endif
          if( ndxord <= _maxord ){
            auto [itmon2,ins2] = coefmon.insert( std::make_pair( SPolyMon(*itvar,ndxord), add2mon ) );
            if( !ins2 ) itmon2->second += add2mon;
#ifdef MC__SICMODEL_DEBUG_SLIFT
            std::cout << "_slift1D: " << itmon2->first.display(options.BASIS) << "  " << add2mon << std::endl;
#endif
          }
          else{
            auto [itmon2,ins2] = coefmon.insert( std::make_pair( SPolyMon(), add2mon ) );
            if( !ins2 ) itmon2->second += add2mon;
#ifdef MC__SICMODEL_DEBUG_SLIFT
            std::cout << "_slift1D: " << itmon2->first.display(options.BASIS) << "  " << add2mon << std::endl;
#endif
          }
          continue;
        }
        case Options::NONE:
        {
          T add2mon = dscal*coef0*TOne;
          auto [itmon,ins] = coefmon.insert( std::make_pair( SPolyMon(), add2mon ) );
          if( !ins ) itmon->second += add2mon;
#ifdef MC__SICMODEL_DEBUG_SLIFT
          std::cout << "_slift1D: " << itmon->first.display(options.BASIS) << "  " << itmon->second << std::endl;
#endif
          continue;
        }}

      case Options::POW:
        switch( options.HOT_SPLIT ){ 
        case Options::FULL:
        {
          T add2mon1 = (0.5*dscal)*coef0*(ndxord%2? TOne: TZerOne);
          for( auto const& [ivar,iord] : mon0.expr ){
            auto [itmon1,ins1] = coefmon.insert( std::make_pair( SPolyMon(ivar,iord), add2mon1 ) );
            if( !ins1 ) itmon1->second += add2mon1;
          }
          T add2mon2 = (0.5*dscal)*coef0*(mon0.gcexp()%2? TOne: TZerOne);         
          if( ndxord <= _maxord ){
            auto [itmon2,ins2] = coefmon.insert( std::make_pair( SPolyMon(*itvar,ndxord), add2mon2 ) );
            if( !ins2 ) itmon2->second += add2mon2;
          }
          else{
            auto [itmon2,ins2] = coefmon.insert( std::make_pair( SPolyMon(), add2mon2 ) );
            if( !ins2 ) itmon2->second += add2mon2; 
          }
          continue;
        }
        case Options::SIMPLE:
        {
          T add2mon1 = (0.5*dscal)*coef0*(ndxord%2? TOne: TZerOne);
          auto [itmon1,ins1] = coefmon.insert( std::make_pair( mon0, add2mon1 ) );
          if( !ins1 ) itmon1->second += add2mon1;
          T add2mon2 = (0.5*dscal)*coef0*(mon0.gcexp()%2? TOne: TZerOne);         
          if( ndxord <= _maxord ){
            auto [itmon2,ins2] = coefmon.insert( std::make_pair( SPolyMon(*itvar,ndxord), add2mon2 ) );
            if( !ins2 ) itmon2->second += add2mon2;
          }
          else{
            auto [itmon2,ins2] = coefmon.insert( std::make_pair( SPolyMon(), add2mon2 ) );
            if( !ins2 ) itmon2->second += add2mon2; 
          }
          continue;
        }
        case Options::NONE:
        {
          T add2mon = dscal*coef0*(ndxord%2 || mon0.gcexp()%2? TOne: TZerOne); continue;
          auto [itmon,ins] = coefmon.insert( std::make_pair( SPolyMon(), add2mon ) );
          if( !ins ) itmon->second += add2mon;           
          continue;
        }}
      }
    }

    // Else introduce new product monomial
    SPolyMon mon = mon0; // local copy for modification
#ifdef MC__SPOLYEXPR_DEBUG_SPROD
    std::cout << "mon: " << mon.display(options.BASIS) << "  (" << *itvar << "," << ndxord << ")" << std::endl;
#endif
    mon.tord += ndxord;
    assert( mon.expr.insert( std::make_pair( *itvar, ndxord ) ).second );
    auto [itmon,ins] = coefmon.insert( std::make_pair( mon, coef0*dscal ) );
    if( !ins ) itmon->second += coef0*dscal;
  }
}

template <typename T>
inline void
SICModel<T>::_sdisp1D
( t_coefmon const& coefmon, std::string const& name, std::ostream& os )
const
{
  os << name;
  for( auto it=coefmon.begin(); it!=coefmon.end(); ++it ){
    if( it != coefmon.begin() ) os << " + ";
    os << it->second;
      for( auto const& [ivar,iord] : it->first.expr )
        os << "·T" << iord << "[" << ivar << "]";
  }
  //os << std::endl;
}

template <typename T>
inline void
SICModel<T>::_sdisp1D
( std::map<unsigned,t_coefmon> const& coefmap, typename t_var::const_iterator itvar,  
  std::string const& name, std::ostream& os )
const
{
  os << name;
  bool first = true;
  for( auto const& [i,coefmon] : coefmap ){
    if( !first ) os << " + T" << i << "[" << *itvar << "] ·";
    os << " { ";
    _sdisp1D( coefmon, "", os );
    os << " }";
    first = false;
  }    
  os << std::endl;
}

template <typename T>
inline void
SICModel<T>::_sdisp1D
( std::vector<SICVar<T>> const& vec, unsigned const ivar, 
  std::string const& name, std::ostream& os )
const
{
  os << name;
  for( unsigned i=0; i<vec.size(); i++ ){
    if( i ) os << " + T" << i << "[" << ivar << "] ·";
    os << " { ";
    bool first = true;
    for( auto const& [mon,coef] : vec[i].coefmon() ){
      if( !first ) os << " + ";
      first = false;
      os << coef;
      for( auto const& [ivar,iord] : mon.expr )
        os << "·T" << iord << "[" << ivar << "]";
    }
    os << " }";
  }
  os << std::endl;
}

template <typename T>
inline void
SICModel<T>::_sprod1D
( std::map<unsigned,t_coefmon> const& sp1map,
  std::map<unsigned,t_coefmon> const& sp2map,
  t_coefmon& coefmon, std::set<unsigned> const& ndxvar,
  typename t_var::const_iterator itvar )
const
{
  // construct product matrix of polynomial coefficients
  auto itvarnext = itvar;
  if( !ndxvar.empty() ) std::advance( itvarnext, 1 );
  std::map<std::pair<unsigned,unsigned>,t_coefmon> sp12map;
  for( auto& [ndx1,coefmon1] : sp1map ){
    // empty monomial in sp1
    if( coefmon1.empty() ) continue; 
    // constant monomial in sp1
    if( coefmon1.size() == 1 && !coefmon1.begin()->first.tord ){
      for( auto& [ndx2,coefmon2] : sp2map ){
        auto [it12,ins] = sp12map.insert( std::make_pair( std::make_pair(ndx1,ndx2), t_coefmon() ) );
        assert( ins ); // map is initially empty
        _sscal1D( coefmon2, coefmon1.begin()->second, it12->second );
      }
      continue;
    }
    // general monomial in sp1
    for( auto& [ndx2,coefmon2] : sp2map ){
      // empty monomial in sp2
      if( coefmon2.empty() ) continue; // no term
      // constant monomial in sp2
      if( coefmon2.size() == 1 && !coefmon2.begin()->first.tord ){
        auto [it12,ins] = sp12map.insert( std::make_pair( std::make_pair(ndx1,ndx2), t_coefmon() ) );
        assert( ins ); // map is initially empty
        _sscal1D( coefmon1, coefmon2.begin()->second, it12->second );
        continue;
      }
#ifdef MC__SPOLYEXPR_DEBUG_SPROD
      std::cout << "Term (" << ndx1 << "," << ndx2 << "):\n";
#endif
      std::map<unsigned,t_coefmon> sp11map, sp22map;
      for( auto& mon1 : coefmon1 )
        _svec1D( itvarnext, mon1, sp11map );
#ifdef MC__SPOLYEXPR_DEBUG_SPROD
      _sdisp1D( sp11map, itvarnext, "Poly #1: " );
#endif
      for( auto& mon2 : coefmon2 )
        _svec1D( itvarnext, mon2, sp22map );
#ifdef MC__SPOLYEXPR_DEBUG_SPROD
      _sdisp1D( sp22map, itvarnext, "Poly #2: " );
#endif
      auto [it12,ins] = sp12map.insert( std::make_pair( std::make_pair(ndx1,ndx2), t_coefmon() ) );
      assert( ins ); // map is initially empty
      _sprod1D( sp11map, sp22map, it12->second, ndxvar, itvarnext );
    }
  }
  
  // construct 1D product result and allocate higher-order terms as appropriate
  coefmon.clear();
  for( auto const& [ndx12,coefmon12] : sp12map ){
    auto const& [ndx1,ndx2] = ndx12;
    // Product involving two constant terms
    if( !ndx1 && !ndx2 )
      _slift1D( coefmon12, 1., coefmon );
    // Product involving exactly one constant term
    else if( !ndx1 || !ndx2 )
      _slift1D( coefmon12, 1., coefmon, itvar, ndx1+ndx2 );
    // Product between non-constant terms
    else{
      switch( options.BASIS ){
        // Chebyshev basis functions
        case Options::CHEB:
          _slift1D( coefmon12, .5, coefmon, itvar, ndx1+ndx2 );
          if( ndx1 == ndx2 )
            _slift1D( coefmon12, .5, coefmon );
          else if( ndx1 > ndx2 )
            _slift1D( coefmon12, .5, coefmon, itvar, ndx1-ndx2 );
          else
            _slift1D( coefmon12, .5, coefmon, itvar, ndx2-ndx1 );
          break;
        // Power basis functions
        case Options::POW:
          _slift1D( coefmon12, 1., coefmon, itvar, ndx1+ndx2 );
          break;
      }
    }
  }
#ifdef MC__POLYMODEL_DEBUG_SPROD
  _sdisp1D( coefmon, "Prod: " );
  std::endl;
#endif
}

template <typename T>
inline void
SICModel<T>::_svec1D
( typename t_var::const_iterator itvar, std::pair<SPolyMon,T> const& mon,
  std::map<unsigned,t_coefmon>& mapspoly )
const
{
  auto const& [ivar,iord] = *mon.first.expr.begin();
  if( !mon.first.tord || ivar != *itvar ) // no dependence on variable *itvar 
    mapspoly[ 0 ].insert( mon );
  else     // dependence on variable *itvar of order iord
    mapspoly[ iord ].insert( std::make_pair( SPolyMon( mon.first.tord - iord,
      std::map<unsigned,unsigned>( ++mon.first.expr.begin(), mon.first.expr.end() ) ),
      mon.second ) );
}

template <typename T>
inline void
SICModel<T>::_svec1Dfull
( unsigned const ivar, std::pair<SPolyMon,T> const& coefmon,
  std::vector<SICVar<T>>& vec )
const
{
  auto& [mon,coef] = coefmon;
  auto ie = mon.expr.find( ivar );
  if( ie == mon.expr.end() ) // no dependence on variable ivar 
    vec[ 0 ].coefmon().insert( coefmon );
  else{
    auto& [ivar,iord] = *ie;
    SPolyMon monmod( mon.tord - iord, mon.expr );
    monmod.expr.erase( ivar ); // remove T[ivar] entry
    vec[ iord ].coefmon().insert( std::make_pair( monmod, coef ) );
  }
}

template <typename T>
template <typename C, typename U>
inline
U
SICModel<T>::_polybound_LSB
( std::map<SPolyMon,C,lt_SPolyMon> const& coefmon, U const* const* bndbasis )
const
{
  //std::cout << "bndpol_naive = " << _polybound_naive( coefmon, bndbasis ) << std::endl;
  
  // Constant or linear model
  if( coefmon.empty() || coefmon.rbegin()->first.tord < 2 )
    return _polybound_naive( coefmon, bndbasis );

  // Quadratic terms in combination with linear terms
  static double const TOL = 1e-8;
  U bndpol = coefmon.cbegin()->first.tord? 0.: coefmon.cbegin()->second;
  auto it1 = coefmon.lower_bound( SPolyMon( 1, std::map<unsigned,unsigned>() ) );
  auto it2 = coefmon.lower_bound( SPolyMon( 2, std::map<unsigned,unsigned>() ) );
  auto it3 = coefmon.lower_bound( SPolyMon( 3, std::map<unsigned,unsigned>() ) );
  std::map<SPolyMon,C,lt_SPolyMon> coeflin; coeflin.insert( it1, it2 );
  for( ; it2!=it3; ++it2 ){
    auto ie2 = it2->first.expr.begin();
    if( ie2->second == 1 ){ // off-diagonal quadratic terms
      bndpol += it2->second * ( bndbasis? bndbasis[ie2->first][1] * bndbasis[(++ie2)->first][1]: TOne );
      continue;
    }
    SPolyMon explin( 1, std::map<unsigned,unsigned>() );
    explin.expr.insert( std::make_pair( ie2->first, 1 ) ); 
    it1 = coeflin.find( explin );
    switch( options.BASIS ){
      case Options::CHEB:
        if( it1 != coeflin.end() && std::fabs( Op<C>::mid(it2->second) ) > TOL ){
          double ai = Op<C>::mid( it1->second ), aii = Op<C>::mid( it2->second );
          bndpol += ( it2->second - aii ) * ( bndbasis? bndbasis[ie2->first][2]: TOne )
                  + ( it1->second - ai  ) * ( bndbasis? bndbasis[ie2->first][1]: TOne )
//                  + aii * ( bndbasis? bndbasis[ie2->first][2]: TOne )
//                  + ai  * ( bndbasis? bndbasis[ie2->first][1]: TOne );
                  + (2.*aii) * Op<U>::sqr( (bndbasis? bndbasis[ie2->first][1]: TOne) + ai/(aii*4.) )
                  - aii - ai*ai/(8.*aii);
//          bndpol += (2.*aii)*Op<U>::sqr((bndbasis? bndbasis[ie2->first][1]: TOne)+ai/(aii*4.))-aii-ai*ai/8./aii;
          //std::cout << "ai[-1,1]+aii[-1,1] = " << aii * TOne + ai * TOne << std::endl;
          //std::cout << "2aii([-1,1]+ai/(4aii))^2-aii-ai^2/(8aii) = " << (2.*aii)*Op<T>::sqr(TOne+ai/(aii*4.))-aii-ai*ai/8./aii << std::endl;
          coeflin.erase( it1 );
        }
        else if( it1 != coeflin.end() ){
          bndpol += it2->second * ( bndbasis? bndbasis[ie2->first][2]: TOne )
                  + it1->second * ( bndbasis? bndbasis[ie2->first][1]: TOne );
          coeflin.erase( it1 );
        }
        else
          bndpol += it2->second * ( bndbasis? bndbasis[ie2->first][2]: TOne );
        break;
      case Options::POW:
        if( it1 != coeflin.end() && std::fabs( Op<C>::l(it2->second) ) > TOL ){
          double ai = Op<C>::mid( it1->second ), aii = Op<C>::mid( it2->second );
          bndpol += ( it2->second - aii ) * ( bndbasis? bndbasis[ie2->first][2]: TZerOne )
                  + ( it1->second - ai  ) * ( bndbasis? bndbasis[ie2->first][1]: TOne )
                  + aii * Op<U>::sqr( (bndbasis? bndbasis[ie2->first][1]: TOne) + ai/(aii*2.) )
                  - ai*ai/(aii*4.);
          coeflin.erase( it1 );
        }
        else if( it1 != coeflin.end() ){
          bndpol += it2->second * ( bndbasis? bndbasis[ie2->first][2]: TZerOne )
                  + it1->second * ( bndbasis? bndbasis[ie2->first][1]: TOne );
          coeflin.erase( it1 );
        }
        else
          bndpol += it2->second * ( bndbasis? bndbasis[ie2->first][2]: TZerOne );
        break;
    }
  }

  // Remaining linear terms
  for( it1=coeflin.begin(); it1!=coeflin.end(); ++it1 ){
    auto ie1 = it1->first.expr.begin();
    bndpol += it1->second * ( bndbasis? bndbasis[ie1->first][1]: TOne );
  }

  // Thrid and higher-order terms
  if( coefmon.rbegin()->first.tord > 2 )
    bndpol += _polybound_naive( coefmon, bndbasis, 3 );

  //std::cout << "bndpol_LSB = " << bndpol << std::endl;
  return bndpol;
}

template <typename T>
template <typename C, typename U>
inline
U
SICModel<T>::_polybound_naive
( std::map<SPolyMon,C,lt_SPolyMon> const& coefmon, U const* const* bndbasis,
  unsigned const minord )
const
{
  // Empty model
  if( coefmon.empty() || coefmon.rbegin()->first.tord < minord ) return 0.;
  auto it = coefmon.lower_bound( SPolyMon( minord, std::map<unsigned,unsigned>() ) );

  // Polynomial bounding in T arithmetic
  if( !bndbasis ){
    switch( options.BASIS ){
      case Options::CHEB:
      {
        T bndpol = 0.;
        if( !it->first.tord ){ bndpol = it->second; ++it; }
        for( ; it!=coefmon.end(); ++it )
          bndpol += it->second * TOne;
        return bndpol;
      }
      case Options::POW:
      {
        T bndpol = 0.;
        if( !it->first.tord ){ bndpol = it->second; ++it; }
        for( ; it!=coefmon.end(); ++it )
          bndpol += it->second * (it->first.gcexp()%2? TOne: TZerOne);
        return bndpol;
      }
    }
  }

  // Polynomial bounding in U arithmetic
  U bndpol( 0. );
  if( !it->first.tord ){ bndpol = it->second; ++it; }
  for( ; it!=coefmon.end(); ++it ){
    U bndmon( 1. );
    for( auto& [ivar,iord] : it->first.expr )
      bndmon *= bndbasis[ivar][iord];
    bndpol += it->second  * bndmon;
  }
  return bndpol;
}

template <typename T>
template <typename C, typename U>
inline
U
SICModel<T>::_polybound
( std::map<SPolyMon,C,lt_SPolyMon> const& coefmon, U const* const* bndbasis,
  int const type )
const
{
  switch( type ){
  case Options::LSB:
    return _polybound_LSB( coefmon, bndbasis );
  case Options::NAIVE: default:
    return _polybound_naive( coefmon, bndbasis );
  }
}

//template <typename T>
//inline
//std::map<SPolyMon,T,lt_SPolyMon>&
//SICModel<T>::_polyprod_LSB
//( std::map<SPolyMon,T,lt_SPolyMon>& coefmon, T const& bndvar )
//const
//{
//  // Constant or linear model
//  if( coefmon.empty() || coefmon.rbegin()->first.tord < 2 )
//    return _polyprod_naive( coefmon, bndvar );

//  // Quadratic terms in combination with linear terms
//  static double const TOL = 1e8;
//  T bndcst = 0.;
//  bool bndcstadd = false;
//  auto it1 = coefmon.lower_bound( SPolyMon( 1, std::map<unsigned,unsigned>() ) );
//  auto it2 = coefmon.lower_bound( SPolyMon( 2, std::map<unsigned,unsigned>() ) );
//  auto it3 = coefmon.lower_bound( SPolyMon( 3, std::map<unsigned,unsigned>() ) );
//  std::map<SPolyMon,T,lt_SPolyMon> coeflin; coeflin.insert( it1, it2 );
//  for( ; it2!=it3; ++it2 ){
//    auto ie2 = it2->first.expr.begin();
//    if( ie2->second == 1 ){ // off-diagonal quadratic terms
//      it2->second *= bndvar;
//      continue;
//    }
//    SPolyMon explin( 1, std::map<unsigned,unsigned>() );
//    explin.expr.insert( std::make_pair( ie2->first, 1 ) ); 
//    it1 = coeflin.find( explin );
//    if( it1 != coeflin.end() && std::fabs( Op<T>::mid(it2->second) ) > TOL ){
//      double ai = Op<T>::mid( it1->second ), aii = Op<T>::mid( it2->second );
//      it2->second -= aii;
//      it2->second *= bndvar;
//      it1->second -= ai;
//      it1->second *= bndvar;
//      switch( options.BASIS ){
//        case Options::CHEB:
//          bndcst += (2.*aii) * Op<T>::sqr( TOne + ai/(aii*4.) ) - aii - ai*ai/(8.*aii);
//          break;
//        case Options::POW:
//          bndcst += aii * Op<T>::sqr( TOne + ai/(aii*2.) ) - ai*ai/(aii*4.);
//          break;
//      }
//      bndcstadd = true;
//      coeflin.erase( it1 );
//    }
//    else if( it1 != coeflin.end() ){
//      it2->second *= bndvar;
//      it1->second *= bndvar;
//      coeflin.erase( it1 );
//    }
//    else
//      it2->second *= bndvar;
//  }
//  if( bndcstadd ){
//    auto [it0,ins] = coefmon.insert( std::make_pair( SPolyMon(), bndcst ) );
//    if( !ins ) it0->second += bndcst;
//  }
//  
//  // Remaining linear terms
//  for( it1=coeflin.begin(); it1!=coeflin.end(); ++it1 )
//    it1->second *= bndvar;

//  // Thrid and higher-order terms
//  if( coefmon.rbegin()->first.tord > 2 )
//    return _polyprod_naive( coefmon, bndvar, 3 );

//  return coefmon;
//}

//template <typename T>
//inline
//std::map<SPolyMon,T,lt_SPolyMon>&
//SICModel<T>::_polyprod_naive
//( std::map<SPolyMon,T,lt_SPolyMon>& coefmon, T const& bndvar,
//  unsigned const minord )
//const
//{
//  // Empty model
//  if( coefmon.empty() || coefmon.rbegin()->first.tord < minord )
//    return coefmon;

//  // Polynomial product with bound
//  auto it = coefmon.lower_bound( SPolyMon( minord, std::map<unsigned,unsigned>() ) );
//  for( ; it!=coefmon.end(); ++it )
//    it->second *= bndvar;
//  return coefmon;
//}

//template <typename T>
//inline
//std::map<SPolyMon,T,lt_SPolyMon>&
//SICModel<T>::_polyprod
//( std::map<SPolyMon,T,lt_SPolyMon>& coefmon, T const& bndvar,
//  int const type )
//const
//{
//  switch( type ){
//  case Options::LSB:
//    return _polyprod_LSB( coefmon, bndvar );
//  case Options::NAIVE: default:
//    return _polyprod_naive( coefmon, bndvar );
//  }
//}

////////////////////////////////// SICVar ///////////////////////////////////////

template <typename T>
inline
void
SICVar<T>::_init
()
{
  _bndpol = nullptr;
  _bndT   = nullptr;
  return;
}

template <typename T>
inline
void
SICVar<T>::_cleanup
()
{
  _unset_bndpol();
  _unset_bndT();
  _ndxvar.clear();
  _coefmon.clear();
}

template <typename T>
inline
void
SICVar<T>::_reinit
()
{
  _cleanup();
  _init();
}

template <typename T>
inline
void
SICVar<T>::_unset_bndpol
()
{
  delete _bndpol;
  _bndpol = nullptr;
}

template <typename T>
inline
void
SICVar<T>::_set_bndpol
( const T*bndpol )
{
  if( !bndpol ){
    if( _bndpol ) delete _bndpol;
    _bndpol = 0;
  }
  else if( !_bndpol )
    _bndpol = new T( *bndpol );
  else
    *_bndpol = *bndpol;
}

template <typename T>
inline
void
SICVar<T>::_set_bndpol
( const T&bndpol )
{
  if( !_bndpol )
    _bndpol = new T( bndpol );
  else
    *_bndpol = bndpol;
}

template <typename T>
inline
void
SICVar<T>::_unset_bndT
()
{
  delete _bndT;
  _bndT = nullptr;
}

template <typename T>
inline
void
SICVar<T>::_set_bndT
( const T*bndT )
{
  if( !bndT ){
    if( _bndT ) delete _bndT;
    _bndT = 0;
  }
  else if( !_bndT )
    _bndT = new T( *bndT );
  else
    *_bndT = *bndT;
}

template <typename T>
inline
void
SICVar<T>::_set_bndT
( const T&bndT )
{
  if( !_bndT )
    _bndT = new T( bndT );
  else
    *_bndT = bndT;
}

template <typename T>
inline
SICVar<T>::SICVar
( SICVar<T> const& CV )
{
  _init();
  _set( CV );
}

template <typename T>
inline
SICVar<T>&
SICVar<T>::operator=
( SICVar<T> const& CV )
{
  _set( CV );
  return *this;
}

template <typename T>
inline
SICVar<T>&
SICVar<T>::_set
( SICVar<T> const& CV )
{
  // Same SICVar?
  if( this == &CV ) return *this;

  // Set coefficients and remainder
  _CM      = CV._CM;
  _ndxvar  = CV._ndxvar;
  _coefmon = CV._coefmon;

  // Set bounds
  _set_bndpol( CV._bndpol );
  _set_bndT( CV._bndT );

  return *this;
}

template <typename T>
inline
SICVar<T>::SICVar
( SICVar<T> && CV )
{
#ifdef MC__SICMODEL_TRACE
    std::cerr << "-- SICVar<T>( SICVar<T> && )\n";
#endif
  _init();
  _set( CV );
}

template <typename T>
inline
SICVar<T>&
SICVar<T>::operator=
( SICVar<T> && CV )
{
#ifdef MC__SICMODEL_TRACE
    std::cerr << "-- SICVar<T>& operator= ( SICVar<T> && )\n";
#endif
  _set( CV );
  return *this;
}

template <typename T>
inline
SICVar<T>&
SICVar<T>::_set
( SICVar<T> && CV )
{
  // Same SICVar?
  if( this == &CV ) return *this;

  // Set coefficients and remainder
  std::swap( _CM      , CV._CM      );
  std::swap( _ndxvar  , CV._ndxvar  );
  std::swap( _coefmon , CV._coefmon );

  // Set bounds
  std::swap( _bndpol  , CV._bndpol  );
  std::swap( _bndT    , CV._bndT    );

  return *this;
}

template <typename T>
inline
SICVar<T>::SICVar
( SICModel<T>*CM )
: _CM( CM )
{
  _init();
  _set_bndpol( 0. );
  if( _CM->options.MIXED_IA ) _set_bndT( 0. );
}

template <typename T>
inline
SICVar<T>::SICVar
( double const d, SICModel<T>*CM )
: _CM( CM ), _bndT(nullptr), _bndpol(nullptr)
{
  _init();
  if( d == 0. ) return;
  _coefmon.insert( std::make_pair( SPolyMon(), d ) );
  _set_bndpol( d );
  if( !_CM || _CM->options.MIXED_IA ) _set_bndT( d );
}

template <typename T>
inline
SICVar<T>&
SICVar<T>::operator=
( double const d )
{
  _reinit();
  _CM = nullptr;
  if( d != 0. ) _coefmon.insert( std::make_pair( SPolyMon(), d ) );
  _set_bndpol( d );
  _set_bndT( d );
  return *this;
}

template <typename T>
inline
SICVar<T>::SICVar
( T const& B, SICModel<T>*CM )
: _CM( CM )
{
  _init();
  _coefmon.insert( std::make_pair( SPolyMon(), B ) );
  _set_bndpol( B );
  if( !_CM || _CM->options.MIXED_IA ) _set_bndT( B );
}

template <typename T>
inline
SICVar<T>&
SICVar<T>::operator=
( T const& B )
{
  _reinit();
  _CM = nullptr;
  _coefmon.insert( std::make_pair( SPolyMon(), B ) );
  _set_bndpol( B );
  _set_bndT( B );
  return *this;
}

template <typename T>
inline
SICVar<T>::SICVar
( SICModel<T>* CM, unsigned const ivar, T const& X )
: _CM( CM )
{
  _init();
  _set( ivar, X );
}

template <typename T>
inline
SICVar<T>&
SICVar<T>::_set
( unsigned const i, T const& X, bool const updMod )
{
  if( !_CM ) throw typename SICModel<T>::Exceptions( SICModel<T>::Exceptions::INIT );  

  // Keep data for variable #ivar in model environment
  if( updMod ) _CM->_set( i, X );

  // Populate model variable
  _ndxvar.clear();
  _ndxvar.insert(i);
  _coefmon.clear();
  auto [itmon,ins] = _coefmon.insert( std::make_pair( SPolyMon(), _refvar(i) ) );
  if( _CM->_maxord && !isequal(_scalvar(i),0.) ){
    _coefmon.insert( std::make_pair( SPolyMon(i), _scalvar(i) ) );
    _set_bndpol( _bndvar(i) );
  }
  else{
    _set_bndpol( _refvar(i) );
    itmon->second += _bndvar(i) - _refvar(i);
  }
  if( _CM->options.MIN_FACTOR >= 0. )
    simplify( _CM->options.BASIS, _CM->options.MIN_FACTOR );
#if 0
  std::cout << "coefmon size: " << _coefmon.size() << std::endl;
#endif

  // Interval bounds
  if( _CM->options.MIXED_IA ) _set_bndT( _bndvar(i) );
  else                        _unset_bndT();
  return *this;
}
#if 0
template <typename T>
inline
SICVar<T>&
SICVar<T>::set
( const SICVar<T>& CV, const bool reset )
{
  if( reset ){
    _ndxvar.clear();
    _coefmon.clear();
  }
  auto&& itord = CV._coefmon.upper_bound( SPolyMon( maxord()+1, std::map<unsigned,unsigned>() ) );
  for( auto&& it=CV._coefmon.begin(); it!=itord; ++it ){
    for( auto&& [var,order] : it->first.expr ) _ndxvar.insert( var );
    _coefmon.insert( *it );
  }
  _bndrem  = CV._bndrem + CV.bndord(maxord()+1);  // CV may have a higher order than *this
  _unset_bndT();
  _unset_bndpol();

  return *this;
}

template <typename T>
inline
SICVar<T>&
SICVar<T>::get
( SICVar<T>& var, const bool reset )
{
  if( !_PM ){
    var._ndxmon.clear();
    if( var._PM && var._PM->sparse() ) var._ndxmon.insert(0);
    var._coefmon[0] = _coefmon[0];
    // Reset monomial coefficients to 0
    if( reset ) _coefmon[0] = 0.;
    return *this;
  }
  if( !var._PM || nvar() < var.nvar() || nord() < var.nord() ) // looks fishy...
    return *this;

  // Copy monomial coefficients from *this into var
  unsigned*iexp = new unsigned[nvar()];
  // Case: sparse representation
  for( auto it=_ndxmon.begin(); it!=_ndxmon.end() && *it<_posord(nord()+1); ++it ){
    for( unsigned ivar=0; ivar<nvar(); ivar++ )
      iexp[ivar] = ( ivar<var.nvar()? var._expmon(*it)[ivar]: 0 );
    //unsigned imon = _loc_expmon(iexp);
    auto pmon = _ndxmon.insert( _loc_expmon(iexp) );
    var._coefmon[*it] = _coefmon[*(pmon.first)];
    // Reset monomial coefficients to 0
    if( reset ) _coefmon[*(pmon.first)] = 0.;
  }
  // Case: dense representation
  for( unsigned jmon=0; _ndxmon.empty() && jmon<var.nmon(); jmon++ ){
    for( unsigned ivar=0; ivar<nvar(); ivar++ )
      iexp[ivar] = ( ivar<var.nvar()? var._expmon(jmon)[ivar]: 0 );
    var._coefmon[jmon] = _coefmon[_loc_expmon(iexp)];
    // Reset monomial coefficients to 0
    if( reset ) _coefmon[_loc_expmon(iexp)] = 0.;
  }
  delete[] iexp;
  *var._bndrem = 0.;
  var._unset_bndT();
  var._unset_bndpol();
  var._bndord_uptd = false;

  if( reset ){
    _unset_bndT();
    _unset_bndpol();
    _bndord_uptd = false;
  }

  return *this;
}
#endif

template <typename T>
inline
T
SICVar<T>::IP
( double const* x ) const
{
  T Pval = 0.;
  for( auto& [mon,coef] : _coefmon ){
    double monval = 1.;
    for( auto& [ivar,iord] : mon.expr ){
      switch( _CM->options.BASIS ){
        case SICModel<T>::Options::CHEB:
          monval *= isequal(_scalvar(ivar),0.)? _refvar(ivar): 
                    mc::cheb((x[ivar]-_refvar(ivar))/_scalvar(ivar),iord);
          break;
        case SICModel<T>::Options::POW:
          monval *= isequal(_scalvar(ivar),0.)? _refvar(ivar): 
                    std::pow((x[ivar]-_refvar(ivar))/_scalvar(ivar),(int)iord);
          break;
      }
    }
    
    Pval += coef * monval;
  }
  return Pval;
}

template <typename T>
inline
double
SICVar<T>::P
( double const* x ) const
{
  double Pval = 0.;
  for( auto& [mon,coef] : _coefmon ){
    double monval = 1.;
    for( auto& [ivar,iord] : mon.expr ){
      switch( _CM->options.BASIS ){
        case SICModel<T>::Options::CHEB:
          monval *= isequal(_scalvar(ivar),0.)? _refvar(ivar): 
                    mc::cheb((x[ivar]-_refvar(ivar))/_scalvar(ivar),iord);
          break;
        case SICModel<T>::Options::POW:
          monval *= isequal(_scalvar(ivar),0.)? _refvar(ivar): 
                    std::pow((x[ivar]-_refvar(ivar))/_scalvar(ivar),(int)iord);
          break;
      }
    }
    
    Pval += Op<T>::mid(coef) * monval;
  }
  return Pval;
}

template <typename T>
inline
T
SICVar<T>::constant( const bool reset )
{
  auto it_0 = ( _coefmon.empty() || _coefmon.begin()->first.tord? _coefmon.end(): _coefmon.begin() );
  T coefcst = ( it_0 != _coefmon.end()? it_0->second: 0. );
  if( reset && it_0 != _coefmon.end() ){
    _coefmon.erase( it_0 );
    if( _bndpol ) *_bndpol -= coefcst;
    if( _bndT )   *_bndT -= coefcst;
  }
  return coefcst;
}

template <typename T>
inline
T
SICVar<T>::linear
( const unsigned i, const bool reset )
{
  if( i>=nvar() || !nord() ) return 0.;
  auto it_i = ( _coefmon.empty()? _coefmon.end(): _coefmon.find( SPolyMon( i ) ) );
  T coeflin = ( it_i == _coefmon.end() || isequal(_scalvar(i),0.)? 0.: it_i->second/_scalvar(i) );
  if( reset && it_i != _coefmon.end() ){
    _coefmon.erase( it_i );
    _unset_bndpol();
    _unset_bndT();
  }
  return coeflin;
}

template <typename T>
inline
void
SICVar<T>::_scale
( unsigned const ivar, T const& Xivar, typename SICVar<T>::t_coefmon& coefmon )
const
{
  // Nothing to do if model _CM is NULL, i is outside of variable range,
  // or variable range X did not change
  if( !_CM || _ndxvar.find( ivar ) == _ndxvar.end()
   || ( isequal( Op<T>::l(Xivar), Op<T>::l(_bndvar(ivar)) )
     && isequal( Op<T>::u(Xivar), Op<T>::u(_bndvar(ivar)) ) ) ) return;

  // Get coefficients in univariate polynomial representation w.r.t variable #ivar
  std::vector<SICVar<T>> veccoef( nord()+1, SICVar<T>(_CM) );
  for( auto const& [mon,coef] : coefmon ){
    auto ie = mon.expr.find( ivar );
    if( ie == mon.expr.end() ){ // no dependence on variable ivar
#ifdef MC__POLYMODEL_DEBUG_SCALE
      std::cout << "Inserting: " << coef << "  " << mon.display(1) << " (" << mon.tord << ")" << std::endl; 
      if( veccoef[ 0 ]._coefmon.find( mon ) != veccoef[ 0 ]._coefmon.end() ) std::cout << "Already present!" << std::endl;
#endif
      auto&& [it,ins] = veccoef[ 0 ]._coefmon.insert( std::make_pair( mon, coef ) );
      assert( ins );
      for( auto const& [jvar,jord] : mon.expr )
        veccoef[ 0 ]._ndxvar.insert( jvar );
    }
    else{
      auto const& [ivar,iord] = *ie;
      SPolyMon monmod( mon.tord - iord, mon.expr );
      monmod.expr.erase( ivar ); // remove T[ivar] entry
      veccoef[ iord ]._coefmon.insert( std::make_pair( monmod, coef ) );
      for( auto const& [jvar,jord] : monmod.expr )
        veccoef[ iord ]._ndxvar.insert( jvar );
    }
  }
  for( auto&& coef : veccoef ){
    coef._unset_bndpol();
    coef._unset_bndT();
  }
#ifdef MC__POLYMODEL_DEBUG_SCALE
  _CM->_sdisp1D( veccoef, ivar, "Var #i: " );
#endif
 
  // Nothing to scale if independent of current variable #i
  bool nodep = true;
  for( unsigned k=1; nodep && k<=nord(); k++ )
    if( !veccoef[k].coefmon().empty() ) nodep = false;
  if( nodep ) return;

  // Compose with rescaled inner variable
  SICVar<T> CVivar( _CM ); CVivar._set( ivar, Xivar, false );
#ifdef MC__POLYMODEL_DEBUG_SCALE
  std::cout << "CVivar[" << ivar << "]:" << CVivar;
#endif
  if( !isequal(_scalvar(ivar),0.) ){
    CVivar -= _refvar(ivar);
    CVivar *= Op<T>::diam(Xivar) / (2.*_scalvar(ivar) );
    CVivar += Op<T>::mid(Xivar);
  }
#ifdef MC__POLYMODEL_DEBUG_SCALE
  std::cout << "CVivar[" << ivar << "]:" << CVivar;
#endif
  CVivar = CVivar._rescale( _scalvar(ivar), _refvar(ivar) );
#ifdef MC__POLYMODEL_DEBUG_SCALE
  std::cout << "CVivar[" << ivar << "]:" << CVivar;
#endif
  coefmon = _CM->_composition( veccoef, nord(), CVivar )._coefmon;
}

template <typename T>
inline
SICVar<T>&
SICVar<T>::scale
( unsigned const i, T const& X )
{
  _scale( i, X, _coefmon );
  _unset_bndpol();
  _unset_bndT();
  return *this;
}

template <typename T>
inline
SICVar<T>&
SICVar<T>::scale
( T const* X )
{
  // Return *this if null pointer to model _CM or variable ranges X
  if( !X || !_CM ) return *this;
  for( unsigned i : _ndxvar ) scale( i, X[i] );
  if( _CM && _CM->options.MIN_FACTOR >= 0. )
    simplify( _CM->options.BASIS, _CM->options.MIN_FACTOR );
  return *this;
}

template <typename T>
inline
typename SICVar<T>::t_coefmon
SICVar<T>::unscale
()
const
{
  // Return *this if null pointer to model _CM or variable ranges X
  if( !_CM ) return _coefmon;
  t_coefmon coefmon = _coefmon;
  for( unsigned i : _ndxvar ) _scale( i, T(-1,1), coefmon );
  if( _CM && _CM->options.MIN_FACTOR >= 0. )
    simplify( _CM->options.BASIS, _CM->options.MIN_FACTOR );
  return coefmon;
}

template <typename T>
inline
SICVar<T>&
SICVar<T>::simplify
( int const& BASIS, double const& TOL, int const& TORD )
{
  if( _coefmon.empty() ) return *this;
  for( auto it=_coefmon.begin(); it!=_coefmon.end(); ){
    // Eliminate any terms with zero coefficient
    if( Op<T>::abs(it->second) == 0e0 ){
      it = _coefmon.erase( it );
      continue;
    }
    // Eliminate any non-constant terms with small enough coefficient or large enough total order
    if( ( TOL > 0e0 && it->first.tord && Op<T>::abs(it->second) <= TOL )
     || ( TORD >= 0 && (int)it->first.tord > TORD ) ){
      T bndmon = ( BASIS == SICModel<T>::Options::POW && !(it->first.gcexp()%2) ?
                   SICModel<T>::TZerOne : SICModel<T>::TOne );
      if( _coefmon.empty() || _coefmon.begin()->first.tord ) 
        _coefmon.insert( std::make_pair( SPolyMon(), bndmon ) );
      else
        _coefmon.begin()->second += bndmon;
      _unset_bndpol();
      it = _coefmon.erase( it );
      continue;
    }
    ++it; // only increment if current element was not erased
  }
  return *this;
}

template <typename T>
inline
T
SICVar<T>::_simplify_monomial
( typename SICVar<T>::t_coefmon& coefmon, bool const scaled,
  double const& TOL, int const& TORD )
const
{
  T bndrem( 0e0 );
  if( coefmon.empty() ) return bndrem;
  
  for( auto it=coefmon.begin(); it!=coefmon.end(); ){
    auto const& [mon,coef] = *it;

    // Eliminate any terms with zero coefficient
    if( Op<T>::abs(coef) == 0e0 ){
      it = coefmon.erase( it );
      continue;
    }

    // Eliminate any non-constant terms with small enough coefficient or large enough total order
    //if( ( TOL > 0e0 && mon.tord && std::fabs(coef) <= TOL )
    if( ( TOL > 0e0 && Op<T>::abs(coef) <= TOL )
     || ( TORD >= 0 && (int)mon.tord > TORD ) ){
      // compute monomial bound
      T bndmon( 1e0 );
      for( auto const& [ivar,iord] : mon.expr )
        bndmon *= Op<T>::pow( scaled? SICModel<T>::TOne: _bndvar(ivar), (int)iord );
      bndrem += bndmon * coef;
      it = coefmon.erase( it );
      continue;
    }
    ++it; // only increment if current element was not erased
  }
  return bndrem;
}

template <typename T>
inline
void
SICVar<T>::_simplify
( typename SICVar<T>::t_coefmon& coefmon )
const
{
  if( coefmon.empty() ) return;
  for( auto it=coefmon.begin(); it!=coefmon.end(); ){
    // Eliminate any terms with zero coefficient
    if( Op<T>::abs(it->second) == 0e0 ){
      it = coefmon.erase( it );
      continue;
    }
    ++it; // only increment if current element was not erased
  }
}

template <typename T>
inline
void
SICVar<T>::_to_monomial
( unsigned const ivar, typename SICVar<T>::t_coefmon& coefmon )
const
{
  // Get coefficients in univariate polynomial representation w.r.t variable #ivar
  std::vector<t_coefmon> veccoef( nord()+1, t_coefmon() );
  for( auto [mon,coef] : coefmon ){
    auto ie = mon.expr.find( ivar );
    if( ie == mon.expr.end() ) // no dependence on variable ivar 
      veccoef[ 0 ].insert( std::make_pair( mon, coef ) );
    else{
      auto const& [ivar,iord] = *ie;
      SPolyMon monmod( mon.tord - iord, mon.expr );
      monmod.expr.erase( ivar ); // remove T[ivar] entry
      veccoef[ iord ].insert( std::make_pair( monmod, coef ) );
    }
  }
 
  // Nothing to convert if independent of current variable #i
  bool nodep = true;
  for( unsigned k=1; nodep && k<=nord(); k++ )
    if( !veccoef[k].empty() ) nodep = false;
  if( nodep ) return;

  // Convert current variable #i to monomial
  // Based on SCONCM function (http://www.netlib.org/math/docpdf/ch11-03.pdf)
  double TP = 1e0;
  // do 20 J = 0, N-2
  for( int j=0; j<(int)nord()-1; j++ ){
    // do 10 I = N-2, J, -1
    for( int i=(int)nord()-2; i>=j; i-- ){
      // COEFF(I) = COEFF(I) - COEFF(I+2)
      for( auto& [mon,coef] : veccoef[i+2] ){
        auto [itmon, ins] = veccoef[i].insert( std::make_pair( mon, -coef ) );
        if( !ins ) itmon->second -= coef;
      }
    // 10    continue
    }
    // COEFF(J+1) = .5E0 * COEFF(J+1)
    for( auto& [mon,coef] : veccoef[j+1] ) coef /= 2e0;
    // COEFF(J) = TP * COEFF(J)
    for( auto& [mon,coef] : veccoef[j] ) coef *= TP;
    // TP = 2.E0 * TP
    TP *= 2e0;
  // 20 continue
  }
  // COEFF(N) = TP * COEFF(N)
  for( auto& [mon,coef] : veccoef[nord()] ) coef *= TP;
  // COEFF(N-1) = TP * COEFF(N-1)
  for( auto& [mon,coef] : veccoef[nord()-1] ) coef *= TP;

  // Collect monomials back together into <a>coefmon</a>
  coefmon = veccoef[0];
  for( unsigned iord=1; iord<=nord(); iord++ ){
    for( auto& [mon,coef] : veccoef[iord] ){
      SPolyMon monmod( mon.tord + iord, mon.expr );
      monmod.expr[ ivar ] = iord;
      coefmon[ monmod ] = coef;
    }
  }
}

template <typename T>
inline
typename SICVar<T>::t_coefmon
SICVar<T>::to_monomial
( bool const scaled )
const
{
  if( !_CM || _coefmon.empty() || !nord() ) return _coefmon;
  //std::cout << *this;

  // Convert to monomial form
  t_coefmon coefmon = _coefmon;
  if( !scaled ){
    for( unsigned i : _ndxvar ){
      _scale( i, SICModel<T>::TOne, coefmon );
    }
  }
  
  for( unsigned i : _ndxvar ){
    _to_monomial( i, coefmon );
    _simplify( coefmon );
  }
  
  return coefmon;
}

template <typename T>
inline
std::pair<typename SICVar<T>::t_coefmon,T>
SICVar<T>::to_monomial
( bool const scaled, double const TOL, int const TORD )
const
{
  auto&& coefmon = to_monomial( scaled );
  auto&& bndrem  = _simplify_monomial( coefmon, scaled, TOL, TORD );
  return std::make_pair( coefmon, bndrem );
}

template <typename T>
inline
std::string
SICVar<T>::display
( typename SICVar::t_coefmon coefmon, int const& BASIS, int const& IDISP )
const
{
  std::ostringstream out;
  out << std::endl << std::scientific << std::setprecision(IDISP) << std::right;

  // Sparse multivariate polynomial
  for( auto const& [mon,coef] : coefmon )
    out << "[" << std::setw(IDISP+7) << Op<T>::l(coef)
        << "," << std::setw(IDISP+7) << Op<T>::u(coef) << "]  "
        << std::setw(2) << mon.tord << "  " << mon.display( BASIS )
        << std::endl;

  return out.str();
}

template <typename T>
inline
std::ostream&
operator<<
( std::ostream& out, SICVar<T> const& CV )
{
  unsigned IDISP = CV._CM? CV._CM->options.DISPLAY_DIGITS: 7;
  int BASIS = CV._CM? CV._CM->options.BASIS: SICModel<T>::Options::CHEB;
  out << CV.display( CV._coefmon, BASIS, IDISP )
      << std::scientific << std::setprecision(IDISP) << std::right;

  // Range bounder
  out << std::right << "   B     =  "
      << "[" << std::setw(IDISP+7) << Op<T>::l(CV.B())
      << "," << std::setw(IDISP+7) << Op<T>::u(CV.B()) << "]"
      << std::endl;

  // Index set
#if 0
  out << std::right << "   I     =  {";
  for( auto const& ndx : CV._ndxvar )
    out << std::right << " " << ndx;
  out << std::right << " }" << std::endl;
#endif
  return out;
}

template <typename T>
inline
SICVar<T>
operator+
( SICVar<T> const& CV )
{
  return CV;
}

template <typename T>
template <typename U>
inline
SICVar<T>&
SICVar<T>::operator+=
( SICVar<U> const& CV )
{
  if( _CM && CV._CM && _CM != CV._CM )
    throw typename SICModel<T>::Exceptions( SICModel<T>::Exceptions::SICMODEL );
  if( CV._CM && !_CM ) _CM = CV._CM;

  _ndxvar.insert( CV._ndxvar.begin(), CV._ndxvar.end() );
  for( auto&& [mon,coef] : CV._coefmon ){
    // No warm-start with map::insert unfortunately...
    auto [itmon, ins] = _coefmon.insert( std::make_pair( mon, coef ) );
    if( !ins ) itmon->second += coef;
  }
  if( _CM && _CM->options.MIN_FACTOR >= 0. )
    simplify( _CM->options.BASIS, _CM->options.MIN_FACTOR );

  _unset_bndpol();
  if( _bndT && CV._bndT ) *_bndT += *CV._bndT;
  else _unset_bndT();

  return *this;
}

template <typename T, typename U>
inline SICVar<T>
operator+
( SICVar<T> const& CV1, SICVar<U> const& CV2 )
{
  if( CV1.nmon() >= CV2.nmon() ){
    SICVar<T> CV3( CV1 );
    CV3 += CV2;
    return CV3;
  }
  
  SICVar<T> CV3( CV2 );
  CV3 += CV1;
  return CV3;
}

template <typename T>
inline
SICVar<T>&
SICVar<T>::operator +=
( double const c )
{
  if( _coefmon.empty() || _coefmon.begin()->first.tord ) 
    _coefmon.insert( std::make_pair( SPolyMon(), c ) );
  else
    _coefmon.begin()->second += c;
  if( _bndpol ) *_bndpol += c;
  if( _bndT )   *_bndT += c;
  return *this;
}

template <typename T>
inline
SICVar<T>
operator+
( SICVar<T> const& CV1, double const c )
{
  SICVar<T> CV3( CV1 );
  CV3 += c;
  return CV3;
}

template <typename T>
inline
SICVar<T>
operator+
( double const c, SICVar<T> const& CV2 )
{
  SICVar<T> CV3( CV2 );
  CV3 += c;
  return CV3;
}

template <typename T>
template <typename U>
inline
SICVar<T>&
SICVar<T>::operator+=
( U const& B )
{
  if( Op<T>::abs(B) == 0. ) return *this;
  if( _coefmon.empty() || _coefmon.begin()->first.tord ) 
    _coefmon.insert( std::make_pair( SPolyMon(), B ) );
  else
    _coefmon.begin()->second += B;
  if( _bndpol ) *_bndpol += B;
  if( _bndT )   *_bndT += B;
  return *this;
}

template <typename T>
inline
SICVar<T>
operator+
( SICVar<T> const& CV1, T const& B )
{
  SICVar<T> CV3( CV1 );
  CV3 += B;
  return CV3;
}

template <typename T>
inline
SICVar<T>
operator+
( T const& B, SICVar<T> const& CV2 )
{
  SICVar<T> CV3( CV2 );
  CV3 += B;
  return CV3;
}

template <typename T>
inline
SICVar<T>
operator-
( SICVar<T> const& CV )
{
  //std::cout << "CV:" << CV;
  SICVar<T> CV2;
  CV2.set( CV._CM );
  CV2._ndxvar = CV._ndxvar;
  for( auto& [mon,coef] : CV._coefmon )
    CV2._coefmon.insert( std::make_pair( mon, -coef ) );
  if( CV._bndpol ) CV2._set_bndpol( - *CV._bndpol );
  if( CV._bndT )   CV2._set_bndT( - *CV._bndT );
  //std::cout << "-CV:" << CV2;
  return CV2;
}

template <typename T>
template <typename U>
inline
SICVar<T>&
SICVar<T>::operator-=
( SICVar<U> const& CV )
{
  if( _CM && CV._CM && _CM != CV._CM )
    throw typename SICModel<T>::Exceptions( SICModel<T>::Exceptions::SICMODEL );
  if( CV._CM && !_CM ) _CM = CV._CM;

  _ndxvar.insert( CV._ndxvar.begin(), CV._ndxvar.end() );
  for( auto& [mon,coef] : CV._coefmon ){
    // No warm-start with map::insert unfortunately...
    auto [itmon, ins] = _coefmon.insert( std::make_pair( mon, -coef ) );
    if( !ins ) itmon->second -= coef;
  }
  if( _CM && _CM->options.MIN_FACTOR >= 0. )
    simplify( _CM->options.BASIS, _CM->options.MIN_FACTOR );

  _unset_bndpol();
  if( _bndT && CV._bndT ) *_bndT -= *CV._bndT;
  else _unset_bndT();

  return *this;
}

template <typename T, typename U>
inline
SICVar<T>
operator-
( SICVar<T> const& CV1, SICVar<U> const& CV2 )
{
  if( CV1.nmon() >= CV2.nmon() ){
    SICVar<T> CV3( CV1 );
    CV3 -= CV2;
    return CV3;
  }

  SICVar<T> CV3( -CV2 );
  CV3 += CV1;
  return CV3;
}

template <typename T>
inline
SICVar<T>&
SICVar<T>::operator-=
( double const c )
{
  if( _coefmon.empty() || _coefmon.begin()->first.tord ) 
    _coefmon.insert( std::make_pair( SPolyMon(), -c ) );
  else{
    _coefmon.begin()->second -= c;
    if( Op<T>::abs(_coefmon.begin()->second) == 0. )
      _coefmon.erase( _coefmon.begin() );
  }
  if( _bndpol ) *_bndpol -= c;
  if( _bndT )   *_bndT -= c;
  return *this;
}

template <typename T>
inline
SICVar<T>
operator-
( SICVar<T> const& CV1, double const c )
{
  SICVar<T> CV3( CV1 );
  CV3 -= c;
  return CV3;
}

template <typename T>
inline
SICVar<T>
operator-
( double const c, SICVar<T> const& CV2 )
{
  SICVar<T> CV3( -CV2 );
  CV3 += c;
  return CV3;
}

template <typename T>
template <typename U>
inline
SICVar<T>&
SICVar<T>::operator-=
( U const& B )
{
  if( Op<T>::abs(B) == 0. ) return *this;
  if( _coefmon.empty() || _coefmon.begin()->first.tord ) 
    _coefmon.insert( std::make_pair( SPolyMon(), -B ) );
  else{
    _coefmon.begin()->second -= B;
    if( Op<T>::abs(_coefmon.begin()->second) == 0. )
      _coefmon.erase( _coefmon.begin() );
  }
  if( _bndpol ) *_bndpol -= B;
  if( _bndT )   *_bndT -= B;
  return *this;
}

template <typename T>
inline SICVar<T>
operator-
( SICVar<T> const& CV1, T const& B )
{
  SICVar<T> CV3( CV1 );
  CV3 -= B;
  return CV3;
}

template <typename T>
inline
SICVar<T>
operator-
( T const& B, SICVar<T> const& CV2 )
{
  SICVar<T> CV3( -CV2 );
  CV3 += B;
  return CV3;
}

template <typename T>
inline
SICVar<T>&
SICVar<T>::operator*=
( SICVar<T> const& CV )
{
  if( this == &CV ){
    *this = sqr( CV );
    return *this;
  }

  if( _CM && CV._CM && _CM != CV._CM )
    throw typename SICModel<T>::Exceptions( SICModel<T>::Exceptions::SICMODEL );

  if( CV._CM && !_CM ) _CM = CV._CM;
  else if( !CV._CM && !_CM ) return operator*=( CV.bound() );
  _unset_bndpol();
  
#ifdef MC__SICMODEL_DEBUG_SPROD
  std::cout << "Poly #1:" << *this;
  std::cout << "Poly #2:" << CV;
#endif

  // Extract uncertain parts
  SICVar<T> tmpCV0( *this );
  for( auto& [mon,coef] : tmpCV0._coefmon ) coef -= Op<T>::mid(coef);
  tmpCV0._unset_bndpol();
  T BCV0 = tmpCV0.bound();
#ifdef MC__SICMODEL_DEBUG_SPROD
  std::cout << "Uncertain part Poly #1:" << tmpCV0;
#endif

  SICVar<T> tmpCV( CV );
  for( auto& [mon,coef] : tmpCV._coefmon ) coef -= Op<T>::mid(coef);
  tmpCV._unset_bndpol();
  T BCV = tmpCV.bound();
#ifdef MC__SICMODEL_DEBUG_SPROD
  std::cout << "Uncertain part Poly #2:" << tmpCV;
#endif

  // Uncertainty propagation
  auto it0 = _coefmon.begin(), it1 = tmpCV0._coefmon.begin();
  for( ; it0!=_coefmon.end(); ++it0, ++it1 ) it1->second = Op<T>::mid(it0->second);
  tmpCV0._unset_bndpol();

  auto it2 = CV._coefmon.begin(), it3 = tmpCV._coefmon.begin();
  for( ; it2!=CV._coefmon.end(); ++it2, ++it3 ) it3->second = Op<T>::mid(it2->second);
  tmpCV._unset_bndpol();

  // Coefficient maps for first participating variable
  _ndxvar.insert( CV._ndxvar.begin(), CV._ndxvar.end() );
  auto itvar = _ndxvar.begin();
  std::map<unsigned,t_coefmon> sp1map, sp2map;
  for( auto& [mon,coef] : _coefmon )
    _CM->_svec1D( itvar, std::make_pair(mon,Op<T>::mid(coef)), sp1map );
#ifdef MC__SICMODEL_DEBUG_SPROD
  _CM->_sdisp1D( sp1map, itvar, "Mid part Poly #1: " );
#endif
  for( auto& [mon,coef] : CV._coefmon )
    _CM->_svec1D( itvar, std::make_pair(mon,Op<T>::mid(coef)), sp2map );
#ifdef MC__SICMODEL_DEBUG_SPROD
  _CM->_sdisp1D( sp2map, itvar, "Mid part Poly #2: " );
#endif

  // Recursive product of univariate Chebyshev polynomials
  _CM->_sprod1D( sp1map, sp2map, _coefmon, _ndxvar, itvar );
  if( _CM && _CM->options.MIN_FACTOR >= 0. )
    simplify( _CM->options.BASIS, _CM->options.MIN_FACTOR );
#ifdef MC__SICMODEL_DEBUG_SPROD
  _unset_bndpol();
  std::cout << "Product of mid parts:" << *this;
#endif
  
  // Uncertainty propagation
#if 0
  switch( _CM->options.HOT_SPLIT ){
  case SICModel<T>::Options::NONE:
#ifdef MC__SICMODEL_DEBUG_SPROD
    std::cout << "Uncertainty propagation 1:" << tmpCV0.bound() * BCV << std::endl;
    std::cout << "Uncertainty propagation 2:" << tmpCV.bound() * BCV0 << std::endl;
#endif
    operator+=( tmpCV0.bound() * BCV + tmpCV.bound() * BCV0 + BCV * BCV0 );
    break;
  case SICModel<T>::Options::SIMPLE:
  case SICModel<T>::Options::FULL:
#endif
    tmpCV0 *= BCV;
    tmpCV  *= BCV0;
#ifdef MC__SICMODEL_DEBUG_SPROD
    std::cout << "Uncertainty propagation 1:" << tmpCV0;
    std::cout << "Uncertainty propagation 2:" << tmpCV;
#endif
    operator+=( tmpCV0 + tmpCV + BCV * BCV0 );
#if 0
    break;
  }
#endif  
  // Update bounds
  if( _bndT && CV._bndT ) *_bndT *= *CV._bndT;
  else _unset_bndT();

#ifdef MC__SICMODEL_DEBUG_SPROD
  std::cout << "Product model:" << *this;
#endif
  return *this;
}

template <typename T>
inline
SICVar<T>
operator*
( SICVar<T> const& CV1, SICVar<T> const& CV2 )
{
  if( &CV1 == &CV2 ) return sqr( CV1 );
  if( CV1.nmon() >= CV2.nmon() ){
    SICVar<T> CV3( CV1 );
    CV3 *= CV2;
    return CV3;
  }

  SICVar<T> CV3( CV2 );
  CV3 *= CV1;
  return CV3;
}

template <typename T>
inline
SICVar<T>&
SICVar<T>::operator*=
( double const c )
{
  if( c == 0. ){ *this = 0.; return *this; }
  if( c == 1. ) return *this;
  for( auto& [mon,coef] : _coefmon ) coef *= c;
  if( _CM && _CM->options.MIN_FACTOR >= 0. )
    simplify( _CM->options.BASIS, _CM->options.MIN_FACTOR );
  if( _bndpol ) *_bndpol *= c;
  if( _bndT ) *_bndT *= c;
  return *this;
}

template <typename T>
inline
SICVar<T>
operator*
( SICVar<T> const& CV1, double const c )
{
  SICVar<T> CV3( CV1 );
  CV3 *= c;
  return CV3;
}

template <typename T>
inline
SICVar<T>
operator*
( double const c, SICVar<T> const& CV2 )
{
  SICVar<T> CV3( CV2 );
  CV3 *= c;
  return CV3;
}

template <typename T>
inline
SICVar<T>&
SICVar<T>::operator*=
( T const& B )
{
  if( Op<T>::abs(B) == 0. ){ *this  = 0.; return *this; }
  if( !_CM || _CM->options.PRODBND_SPLIT ){
    double const Bmid = Op<T>::mid(B);
    T const bndmod = bound();
    operator*=( Bmid );
    operator+=( ( B - Bmid ) * bndmod );
  }
  else{
    for( auto& [mon,coef] : _coefmon ) coef *= B;
    if( _CM && _CM->options.MIN_FACTOR >= 0. )
      simplify( _CM->options.BASIS, _CM->options.MIN_FACTOR );
    _unset_bndpol();
    if( _bndT ) *_bndT *= B;
  }
  return *this;
}

template <typename T> inline SICVar<T>
operator*
( SICVar<T> const& CV1, T const& I )
{
  SICVar<T> CV3( CV1 );
  CV3 *= I;
  return CV3;
}

template <typename T>
inline
SICVar<T>
operator*
( T const& I, SICVar<T> const& CV2 )
{
  SICVar<T> CV3( CV2 );
  CV3 *= I;
  return CV3;
}

template <typename T>
inline
SICVar<T>
sqr
( SICVar<T> const& CV )
{
  if( !CV._CM ) return Op<T>::sqr( CV.B() );
  
  // Extract uncertain parts
  SICVar<T> tmpCV( CV );
  for( auto& [mon,coef] : tmpCV._coefmon ) coef -= Op<T>::mid(coef);
  tmpCV._unset_bndpol();
  T BCV = tmpCV.bound();
#ifdef MC__SICMODEL_DEBUG_SPROD
  std::cout << "Uncertain part:" << tmpCV;
#endif

  // Uncertainty propagation
  auto it0 = CV._coefmon.begin(), it1 = tmpCV._coefmon.begin();
  for( ; it0!=CV._coefmon.end(); ++it0, ++it1 ) it1->second = Op<T>::mid(it0->second);
  tmpCV._unset_bndpol();

  // Coefficient maps for first participating variable
  auto itvar = CV._ndxvar.begin();
  std::map<unsigned,typename SICVar<T>::t_coefmon> sp1map, sp2map;
  for( auto& [mon,coef] : CV._coefmon )
    CV._CM->_svec1D( itvar, std::make_pair(mon,Op<T>::mid(coef)), sp1map );
#ifdef MC__SPOLYEXPR_DEBUG_SQR
  SICVar<T> tmp2CV( CV );
  for( auto& [mon,coef] : tmp2CV._coefmon ) coef = Op<T>::mid(coef);
  std::cout << "Mid part Var: " << tmp2CV;
  CV._CM->_sdisp1D( sp1map, itvar, "Mid part Var: " );
#endif

  // Recursive product of univariate Chebyshev polynomials
  SICVar<T> CVSQR;
  CVSQR.set( CV._CM );
  CVSQR._ndxvar = CV._ndxvar;
  CV._CM->_sprod1D( sp1map, sp1map, CVSQR._coefmon, CVSQR._ndxvar, itvar );
  if( CV._CM->options.MIN_FACTOR >= 0. )
    CVSQR.simplify( CV._CM->options.BASIS, CV._CM->options.MIN_FACTOR );
  CVSQR._unset_bndpol();
#ifdef MC__SICMODEL_DEBUG_SPROD
  CVSQR._unset_bndT();
  std::cout << "Product of mid parts:" << CVSQR;
#endif
  
  // Uncertainty propagation
  switch( CV._CM->options.HOT_SPLIT ){
  case SICModel<T>::Options::NONE:
#ifdef MC__SICMODEL_DEBUG_SPROD
    std::cout << "Uncertainty propagation:" << tmpCV.bound() * BCV;
#endif
    CVSQR += 2. * tmpCV.bound() * BCV + Op<T>::sqr( BCV );
    break;
  case SICModel<T>::Options::SIMPLE:
  case SICModel<T>::Options::FULL:
    tmpCV *= BCV;
#ifdef MC__SICMODEL_DEBUG_SPROD
    std::cout << "Uncertainty propagation:" << tmpCV;
#endif
    CVSQR += 2 * tmpCV + Op<T>::sqr( BCV );
    break;
  }
  
  // Bound propagation
  if( CV._CM->options.MIXED_IA ) CVSQR._set_bndT( Op<T>::sqr( CV.bound() ) );
  else CVSQR._unset_bndT();
  return CVSQR;
}

template <typename T>
inline
SICVar<T>&
SICVar<T>::operator/=
( SICVar<T> const& CV )
{
  *this *= inv(CV);
  return *this;
}

template <typename T> inline SICVar<T>
operator/
( SICVar<T> const& CV1, SICVar<T> const& CV2 )
{
  return CV1 * inv(CV2);
}

template <typename T> inline SICVar<T>&
SICVar<T>::operator/=
( double const c )
{
  if( isequal( c, 0. ) )
    throw typename SICModel<T>::Exceptions( SICModel<T>::Exceptions::DIV );
  if( c == 1. ) return *this;
  *this *= (1./c);
  return *this;
}

template <typename T> inline SICVar<T>
operator/
( SICVar<T> const& CV, double const c )
{
  if ( isequal( c, 0. ))
    throw typename SICModel<T>::Exceptions( SICModel<T>::Exceptions::DIV );
  if( c == 1. ) return CV;
  return CV * (1./c);
}

template <typename T> inline SICVar<T>
operator/
( double const c, SICVar<T> const& CV )
{
  if( c == 0. ) return 0.;
  if( c == 1. ) return inv(CV);
  return inv(CV) * c;
}

template <typename T> inline SICVar<T>
inv
( SICVar<T> const& CV )
{
  if( !CV._CM )
    return SICVar<T>( Op<T>::inv( CV.B() ) );
  if ( Op<T>::l(CV.B()) <= 0. && Op<T>::u(CV.B()) >= 0. )
    throw typename SICModel<T>::Exceptions( SICModel<T>::Exceptions::INV );

  SICVar<T> CV2( CV._CM );
  auto const& f = [=]( const double& x ){ return 1./x; };
  if( ( !CV._CM->options.REMEZ_USE
     || CV._CM->options.REMEZ_MIG > Op<T>::diam(CV.B())
     || !CV._CM->_minimax( f, CV, CV2 ) )
     && !CV._CM->_chebinterp( f, CV, CV2, true ) ) // with remainder at bound
    throw typename SICModel<T>::Exceptions( SICModel<T>::Exceptions::COMPOSE );

  if( CV._CM->options.MIXED_IA ) CV2._set_bndT( Op<T>::inv( CV.B() ) );
  return CV2;
}

template <typename T> inline SICVar<T>
sqrt
( SICVar<T> const& CV )
{
  if( !CV._CM )
    return SICVar<T>( Op<T>::sqrt( CV.B() ) );
  if ( Op<T>::l(CV.B()) < 0. )
    throw typename SICModel<T>::Exceptions( SICModel<T>::Exceptions::SQRT );

  SICVar<T> CV2( CV._CM );
  auto const& f = [=]( const double& x ){ return std::sqrt( x ); };
  if( ( !CV._CM->options.REMEZ_USE
     || CV._CM->options.REMEZ_MIG > Op<T>::diam(CV.B())
     || !CV._CM->_minimax( f, CV, CV2 ) )
     && !CV._CM->_chebinterp( f, CV, CV2, true ) ) // with remainder at bound
    throw typename SICModel<T>::Exceptions( SICModel<T>::Exceptions::COMPOSE );

  if( CV._CM->options.MIXED_IA ) CV2._set_bndT( Op<T>::sqrt( CV.B() ) );
  return CV2;
}

template <typename T> inline SICVar<T>
exp
( SICVar<T> const& CV )
{ 
  if( !CV._CM )
    return SICVar<T>( Op<T>::exp( CV.B() ) );

  SICVar<T> CV2( CV._CM );
  auto const& f = [=]( const double& x ){ return std::exp( x ); };
  if( ( !CV._CM->options.REMEZ_USE
     || CV._CM->options.REMEZ_MIG > Op<T>::diam(CV.B())
     || !CV._CM->_minimax( f, CV, CV2 ) )
     && !CV._CM->_chebinterp( f, CV, CV2, true ) ) // with remainder at bound
    throw typename SICModel<T>::Exceptions( SICModel<T>::Exceptions::COMPOSE );

  if( CV._CM->options.MIXED_IA ) CV2._set_bndT( Op<T>::exp( CV.B() ) );
  return CV2;
}

template <typename T> inline SICVar<T>
log
( SICVar<T> const& CV )
{
  if( !CV._CM )
    return SICVar<T>( Op<T>::log( CV.B() ) );
  if ( Op<T>::l(CV.B()) <= 0. )
    throw typename SICModel<T>::Exceptions( SICModel<T>::Exceptions::LOG );

  SICVar<T> CV2( CV._CM );
  auto const& f = [=]( const double& x ){ return std::log( x ); };
  if( ( !CV._CM->options.REMEZ_USE
     || CV._CM->options.REMEZ_MIG > Op<T>::diam(CV.B())
     || !CV._CM->_minimax( f, CV, CV2 ) )
     && !CV._CM->_chebinterp( f, CV, CV2, true ) ) // with remainder at bound
    throw typename SICModel<T>::Exceptions( SICModel<T>::Exceptions::COMPOSE );

  if( CV._CM->options.MIXED_IA ) CV2._set_bndT( Op<T>::log( CV.B() ) );
  return CV2;
}

template <typename T> inline SICVar<T>
xlog
( SICVar<T> const& CV )
{
#if 0
  return CV * log( CV );
#endif

  if( !CV._CM )
    return SICVar<T>( Op<T>::xlog( CV.B() ) );
  if ( Op<T>::l(CV.B()) <= 0. )
    throw typename SICModel<T>::Exceptions( SICModel<T>::Exceptions::LOG );

  SICVar<T> CV2( CV._CM );
  auto const& f = [=]( const double& x ){ return x*std::log( x ); };
  if( ( !CV._CM->options.REMEZ_USE
     || CV._CM->options.REMEZ_MIG > Op<T>::diam(CV.B())
     || !CV._CM->_minimax( f, CV, CV2 ) )
     && !CV._CM->_chebinterp( f, CV, CV2, false ) ) // remainder not at bound
    throw typename SICModel<T>::Exceptions( SICModel<T>::Exceptions::COMPOSE );

  if( CV._CM->options.MIXED_IA ) CV2._set_bndT( Op<T>::xlog( CV.B() ) );
  return CV2;
}

template <typename T> inline SICVar<T>
pow
( SICVar<T> const& CV, const int n )
{
  if( !CV._CM )
    return SICVar<T>( Op<T>::pow( CV.B(), n ) );

  if( n < 0 ) return pow( inv( CV ), -n );
  SICVar<T> CV2( CV._CM->_intpow( CV, n ) );
  
  if( CV._CM->options.MIXED_IA ) CV2._set_bndT( Op<T>::pow( CV.B(), n ) );
  return CV2;
}

template <typename T> inline SICVar<T>
pow
( const SICVar<T> &CV, double const a )
{
#if 0
  return exp( a * log( CV ) );
#endif

  if( !CV._CM )
    return SICVar<T>( Op<T>::pow( CV.B(), a ) );
  if ( Op<T>::l(CV.B()) <= 0. )
    throw typename SICModel<T>::Exceptions( SICModel<T>::Exceptions::DPOW );

  SICVar<T> CV2( CV._CM );
  auto const& f = [=]( const double& x ){ return std::pow( x, a ); };
  if( ( !CV._CM->options.REMEZ_USE
     || CV._CM->options.REMEZ_MIG > Op<T>::diam(CV.B())
     || !CV._CM->_minimax( f, CV, CV2 ) )
     && !CV._CM->_chebinterp( f, CV, CV2, true ) ) // with remainder at bound
    throw typename SICModel<T>::Exceptions( SICModel<T>::Exceptions::COMPOSE );

  if( CV._CM->options.MIXED_IA ) CV2._set_bndT( Op<T>::pow( CV.B(), a ) );
  return CV2;
}

template <typename T> inline SICVar<T>
pow
( const SICVar<T> &CV1, const SICVar<T> &CV2 )
{
  return exp( CV2 * log( CV1 ) );
}

template <typename T> inline SICVar<T>
pow
( double const a, const SICVar<T> &CV )
{
  return exp( CV * std::log( a ) );
}

template <typename T> inline SICVar<T>
prod
( const unsigned n, const SICVar<T>*CV )
{
  switch( n ){
   case 0:  return 1.;
   case 1:  return CV[0];
   default: return CV[0] * prod( n-1, CV+1 );
  }
}

template <typename T> inline SICVar<T>
monom
( const unsigned n, const SICVar<T>*CV, const unsigned*k )
{
  switch( n ){
   case 0:  return 1.;
   case 1:  return pow( CV[0], (int)k[0] );
   default: return pow( CV[0], (int)k[0] ) * monom( n-1, CV+1, k+1 );
  }
}

template <typename T> inline SICVar<T>
cheb
( const SICVar<T> &CV, const unsigned n )
{
  switch( n ){
    case 0:  return 1.;
    case 1:  return CV;
    default: break;
  }
  SICVar<T> CV2( 2.*(CV*cheb(CV,n-1))-cheb(CV,n-2) );
  
  if( CV._CM->options.MIXED_IA ) CV2._set_bndT( Op<T>::cheb( CV.B(), n ) );
  return CV2;
}

template <typename T> inline SICVar<T>
cos
( const SICVar<T> &CV )
{
  if( !CV._CM )
    return SICVar<T>( Op<T>::cos( CV.B() ) );

  SICVar<T> CV2( CV._CM );
  auto const& f = [=]( const double& x ){ return std::cos( x ); };
  if( ( !CV._CM->options.REMEZ_USE
     || CV._CM->options.REMEZ_MIG > Op<T>::diam(CV.B())
     || !CV._CM->_minimax( f, CV, CV2 ) )
     && !CV._CM->_chebinterp( f, CV, CV2, false ) ) // remainder not at bound
    throw typename SICModel<T>::Exceptions( SICModel<T>::Exceptions::COMPOSE );

  if( CV._CM->options.MIXED_IA ) CV2._set_bndT( Op<T>::cos( CV.B() ) );
  return CV2;
}

template <typename T> inline SICVar<T>
sin
( const SICVar<T> &CV )
{
#if 0
  return cos( CV - PI/2. );
#endif

  if( !CV._CM )
    return SICVar<T>( Op<T>::sin( CV.B() ) );

  SICVar<T> CV2( CV._CM );
  auto const& f = [=]( const double& x ){ return std::sin( x ); };
  if( ( !CV._CM->options.REMEZ_USE
     || CV._CM->options.REMEZ_MIG > Op<T>::diam(CV.B())
     || !CV._CM->_minimax( f, CV, CV2 ) )
     && !CV._CM->_chebinterp( f, CV, CV2, false ) ) // remainder not at bound
    throw typename SICModel<T>::Exceptions( SICModel<T>::Exceptions::COMPOSE );

  if( CV._CM->options.MIXED_IA ) CV2._set_bndT( Op<T>::sin( CV.B() ) );
  return CV2;
}

template <typename T> inline SICVar<T>
tan
( const SICVar<T> &CV )
{
#if 0
  return sin( CV ) / cos( CV );
#endif

  if( !CV._CM )
    return SICVar<T>( Op<T>::tan( CV.B() ) );
  if ( Op<T>::l(Op<T>::cos(CV.B())) <= 0. && Op<T>::u(Op<T>::cos(CV.B())) >= 0. )
    throw typename SICModel<T>::Exceptions( SICModel<T>::Exceptions::TAN );

  SICVar<T> CV2( CV._CM );
  auto const& f = [=]( const double& x ){ return std::tan( x ); };
  if( ( !CV._CM->options.REMEZ_USE
     || CV._CM->options.REMEZ_MIG > Op<T>::diam(CV.B())
     || !CV._CM->_minimax( f, CV, CV2 ) )
     && !CV._CM->_chebinterp( f, CV, CV2, false ) ) // remainder not at bound
    throw typename SICModel<T>::Exceptions( SICModel<T>::Exceptions::COMPOSE );

  if( CV._CM->options.MIXED_IA ) CV2._set_bndT( Op<T>::tan( CV.B() ) );
  return CV2;
}

template <typename T> inline SICVar<T>
acos
( const SICVar<T> &CV )
{
  if( !CV._CM )
    return SICVar<T>( Op<T>::acos( CV.B() ) );
  if ( Op<T>::l(CV.B()) < -1. || Op<T>::u(CV.B()) > 1. )
    throw typename SICModel<T>::Exceptions( SICModel<T>::Exceptions::ACOS );

  SICVar<T> CV2( CV._CM );
  auto const& f = [=]( const double& x ){ return std::acos( x ); };
  if( ( !CV._CM->options.REMEZ_USE
     || CV._CM->options.REMEZ_MIG > Op<T>::diam(CV.B())
     || !CV._CM->_minimax( f, CV, CV2 ) )
     && !CV._CM->_chebinterp( f, CV, CV2, false ) ) // remainder not at bound
    throw typename SICModel<T>::Exceptions( SICModel<T>::Exceptions::COMPOSE );

  if( CV._CM->options.MIXED_IA ) CV2._set_bndT( Op<T>::acos( CV.B() ) );
  return CV2;
}

template <typename T> inline SICVar<T>
asin
( const SICVar<T> &CV )
{
#if 0
  return PI/2. - acos( CV );
#endif

  if( !CV._CM )
    return SICVar<T>( Op<T>::asin( CV.B() ) );
  if ( Op<T>::l(CV.B()) < -1. || Op<T>::u(CV.B()) > 1. )
    throw typename SICModel<T>::Exceptions( SICModel<T>::Exceptions::ASIN );

  SICVar<T> CV2( CV._CM );
  auto const& f = [=]( const double& x ){ return std::asin( x ); };
  if( ( !CV._CM->options.REMEZ_USE
     || CV._CM->options.REMEZ_MIG > Op<T>::diam(CV.B())
     || !CV._CM->_minimax( f, CV, CV2 ) )
     && !CV._CM->_chebinterp( f, CV, CV2, false ) ) // remainder not at bound
    throw typename SICModel<T>::Exceptions( SICModel<T>::Exceptions::COMPOSE );

  if( CV._CM->options.MIXED_IA ) CV2._set_bndT( Op<T>::asin( CV.B() ) );
  return CV2;
}

template <typename T> inline SICVar<T>
atan
( const SICVar<T> &CV )
{
#if 0
  return asin( CV / sqrt( sqr( CV ) + 1. ) );
#endif

  if( !CV._CM )
    return SICVar<T>( Op<T>::atan( CV.B() ) );

  SICVar<T> CV2( CV._CM );
  auto const& f = [=]( const double& x ){ return std::atan( x ); };
  if( ( !CV._CM->options.REMEZ_USE
     || CV._CM->options.REMEZ_MIG > Op<T>::diam(CV.B())
     || !CV._CM->_minimax( f, CV, CV2 ) )
     && !CV._CM->_chebinterp( f, CV, CV2, false ) ) // remainder not at bound
    throw typename SICModel<T>::Exceptions( SICModel<T>::Exceptions::COMPOSE );

  if( CV._CM->options.MIXED_IA ) CV2._set_bndT( Op<T>::atan( CV.B() ) );
  return CV2;
}

template <typename T> inline SICVar<T>
cosh
( const SICVar<T> &CV )
{
#if 0
  return 0.5*(mc::exp(CV)+mc::exp(-CV));
#endif

  if( !CV._CM )
    return SICVar<T>( Op<T>::cosh( CV.B() ) );

  SICVar<T> CV2( CV._CM );
  auto const& f = [=]( const double& x ){ return std::cosh( x ); };
  if( ( !CV._CM->options.REMEZ_USE
     || CV._CM->options.REMEZ_MIG > Op<T>::diam(CV.B())
     || !CV._CM->_minimax( f, CV, CV2 ) )
     && !CV._CM->_chebinterp( f, CV, CV2, false ) ) // remainder not at bound
    throw typename SICModel<T>::Exceptions( SICModel<T>::Exceptions::COMPOSE );

  if( CV._CM->options.MIXED_IA ) CV2._set_bndT( Op<T>::cosh( CV.B() ) );
  return CV2;
}

template <typename T> inline SICVar<T>
sinh
( const SICVar<T> &CV )
{
#if 0
  return 0.5*(mc::exp(CV)-mc::exp(-CV));
#endif

  if( !CV._CM )
    return SICVar<T>( Op<T>::sinh( CV.B() ) );

  SICVar<T> CV2( CV._CM );
  auto const& f = [=]( const double& x ){ return std::sinh( x ); };
  if( ( !CV._CM->options.REMEZ_USE
     || CV._CM->options.REMEZ_MIG > Op<T>::diam(CV.B())
     || !CV._CM->_minimax( f, CV, CV2 ) )
     && !CV._CM->_chebinterp( f, CV, CV2, false ) ) // remainder not at bound
    throw typename SICModel<T>::Exceptions( SICModel<T>::Exceptions::COMPOSE );

  if( CV._CM->options.MIXED_IA ) CV2._set_bndT( Op<T>::sinh( CV.B() ) );
  return CV2;
}

template <typename T> inline SICVar<T>
tanh
( const SICVar<T> &CV )
{
#if 0
  return (mc::exp(2*CV)-1)/(mc::exp(2*CV)+1);
#endif

  if( !CV._CM )
    return SICVar<T>( Op<T>::tanh( CV.B() ) );

  SICVar<T> CV2( CV._CM );
  auto const& f = [=]( const double& x ){ return std::tanh( x ); };
  if( ( !CV._CM->options.REMEZ_USE
     || CV._CM->options.REMEZ_MIG > Op<T>::diam(CV.B())
     || !CV._CM->_minimax( f, CV, CV2 ) )
     && !CV._CM->_chebinterp( f, CV, CV2, false ) ) // remainder not at bound
    throw typename SICModel<T>::Exceptions( SICModel<T>::Exceptions::COMPOSE );

  if( CV._CM->options.MIXED_IA ) CV2._set_bndT( Op<T>::tanh( CV.B() ) );
  return CV2;
}

template <typename T>
inline
SICVar<T>
fabs
( const SICVar<T> &CV )
{
  if( !CV._CM )
    return SICVar<T>( Op<T>::fabs( CV.B() ) );
  if ( Op<T>::l(CV.B()) >= 0. )
    return CV;
  if ( Op<T>::u(CV.B()) <= 0. )
    return -CV;

#if 0
  SICVar<T> CV2( sqrt( sqr(CV) ) );
  if( CV._CM->options.MIXED_IA ) CV2._set_bndT( Op<T>::fabs( CV.B() ) );
  return CV2;
#endif

  SICVar<T> CV2( CV._CM );
  auto const& f = [=]( const double& x ){ return std::fabs( x ); };
  if( ( !CV._CM->options.REMEZ_USE
     || CV._CM->options.REMEZ_MIG > Op<T>::diam(CV.B())
     || !CV._CM->_minimax( f, CV, CV2 ) )
     && !CV._CM->_chebinterp( f, CV, CV2, false ) ) // remainder not at bound
    throw typename SICModel<T>::Exceptions( SICModel<T>::Exceptions::COMPOSE );

  if( CV._CM->options.MIXED_IA ) CV2._set_bndT( Op<T>::fabs( CV.B() ) );
  return CV2;
}

template <typename T>
inline
SICVar<T>
hull
( SICVar<T> const& CV1, SICVar<T> const& CV2 )
{
  throw typename SICModel<T>::Exceptions( SICModel<T>::Exceptions::UNDEF );
#if 0 
  // Neither operands associated to SICModel -- Make intersection in T type     
  if( !CV1._CM && !CV2._CM ){
    T R1 = CV1.B();
    T R2 = CV2.B();
    return Op<T>::hull(R1, R2);
  }

  // First operand not associated to SICModel
  else if( !CV1._CM )
    return hull( CV2, CV1 );

  // Second operand not associated to SICModel
  else if( !CV2._CM ){
    SICVar<T> CVR = CV1.P();
    return CVR + Op<T>::hull( CV1.R(), CV2._coefmon[0]+CV2._bndrem-CVR.B() );
  }

  // SICModel for first and second operands are inconsistent
  else if( CV1._CM != CV2._CM )
    throw typename SICModel<T>::Exceptions( SICModel<T>::Exceptions::SICMODEL );

  // Perform union
  SICVar<T> CV1C( CV1 ), CV2C( CV2 );
  double const eta = CV1._CM->options.REF_POLY;
  T R1C = CV1C.C().R(), R2C = CV2C.C().R(); 
  CV1C.set(T(0.));
  CV2C.set(T(0.));
  T BCVD = (CV1C-CV2C).B();
  return (1.-eta)*CV1C + eta*CV2C + Op<T>::hull( R1C+eta*BCVD, R2C+(eta-1.)*BCVD );
#endif
}

template <typename T>
inline
bool
inter
( SICVar<T>&CVR, SICVar<T> const& CV1, SICVar<T> const& CV2 )
{
    throw typename SICModel<T>::Exceptions( SICModel<T>::Exceptions::UNDEF );
#if 0 
  // Neither operands associated to SICModel -- Make intersection in T type     
  if( !CV1._CM && !CV2._CM ){
    T R1 = CV1.B();
    T R2 = CV2.B();
    T RR( 0. );
    bool flag = Op<T>::inter(RR, R1, R2);
    CVR = RR;
    return flag;
  }

  // First operand not associated to SICModel
  else if( !CV1._CM )
    return inter( CVR, CV2, CV1 );

  // Second operand not associated to SICModel
  else if( !CV2._CM ){
    // First intersect in T arithmetic
    T B2 = CV2.B(), BR;
    if( CV1._CM->options.MIXED_IA && !Op<T>::inter( BR, CV1.B(), B2 ) )
      return false;

    // Perform intersection in PM arithmetic
    T R1 = CV1.R();
    CVR = CV1.P();
    if( !Op<T>::inter(CVR._bndrem, R1, B2-CVR.B()) )
      return false;
//    CVR._center();

    if( CVR._CM->options.MIXED_IA ) CVR._set_bndT( BR );
    else CVR._unset_bndT();
    return true;
  }

  // SICModel for first and second operands are inconsistent
  else if( CV1._CM != CV2._CM )
    throw typename SICModel<T>::Exceptions( SICModel<T>::Exceptions::SICMODEL );

  // First intersect in T arithmetic
  T BR;
  if( CV1._CM->options.MIXED_IA && !Op<T>::inter( BR, CV1.B(), CV2.B() ) )
    return false;

  // Perform intersection in PM arithmetic
  SICVar<T> CV1C( CV1 ), CV2C( CV2 );
  double const eta = CV1._CM->options.REF_POLY;
  T R1C = CV1C.C().R(), R2C = CV2C.C().R(); 
  CV1C.set(T(0.));
  CV2C.set(T(0.));
  CVR = (1.-eta)*CV1C + eta*CV2C;
  CV1C -= CV2C;
  T BCVD = CV1C.B();
  if( !Op<T>::inter( CVR._bndrem, R1C+eta*BCVD, R2C+(eta-1.)*BCVD ) )
    return false;
  if( CVR._CM->options.MIN_FACTOR >= 0. )
    CVR.simplify( CVR._CM->options.BASIS, CVR._CM->options.MIN_FACTOR );
//  CVR._center();

  if( CVR._CM->options.MIXED_IA ) CVR._set_bndT( BR );
  else CVR._unset_bndT();
  return true;
#endif
}

} // namespace mc

#include "mcfadbad.hpp"

namespace fadbad
{

//! @brief Specialization of the structure fadbad::Op for use of the type mc::SICVar of MC++ as a template parameter of the classes fadbad::F, fadbad::B and fadbad::T of FADBAD++
template< typename T > struct Op< mc::SICVar<T> >
{ 
  typedef mc::SICVar<T> CV;
  typedef double Base;
  static Base myInteger( const int i ) { return Base(i); }
  static Base myZero() { return myInteger(0); }
  static Base myOne() { return myInteger(1);}
  static Base myTwo() { return myInteger(2); }
  static double myPI() { return mc::PI; }
  static CV myPos( const CV& x ) { return  x; }
  static CV myNeg( const CV& x ) { return -x; }
  template <typename U> static CV& myCadd( CV& x, U const&  y ) { return x+=y; }
  template <typename U> static CV& myCsub( CV& x, U const&  y ) { return x-=y; }
  template <typename U> static CV& myCmul( CV& x, U const&  y ) { return x*=y; }
  template <typename U> static CV& myCdiv( CV& x, U const&  y ) { return x/=y; }
  static CV myInv( const CV& x ) { return mc::inv( x ); }
  static CV mySqr( const CV& x ) { return mc::pow( x, 2 ); }
  template <typename X, typename Y> static CV myPow( const X& x, const Y& y ) { return mc::pow( x, y ); }
  //static CV myCheb( const CV& x, const unsigned n ) { return mc::cheb( x, n ); }
  static CV mySqrt( const CV& x ) { return mc::sqrt( x ); }
  static CV myLog( const CV& x ) { return mc::log( x ); }
  static CV myExp( const CV& x ) { return mc::exp( x ); }
  static CV mySin( const CV& x ) { return mc::sin( x ); }
  static CV myCos( const CV& x ) { return mc::cos( x ); }
  static CV myTan( const CV& x ) { return mc::tan( x ); }
  static CV myAsin( const CV& x ) { return mc::asin( x ); }
  static CV myAcos( const CV& x ) { return mc::acos( x ); }
  static CV myAtan( const CV& x ) { return mc::atan( x ); }
  static CV mySinh( const CV& x ) { return mc::sinh( x ); }
  static CV myCosh( const CV& x ) { return mc::cosh( x ); }
  static CV myTanh( const CV& x ) { return mc::tanh( x ); }
  static bool myEq( const CV& x, const CV& y ) { return mc::Op<T>::eq(const_cast<CV*>(&x)->bound(),const_cast<CV*>(&y)->bound()); } 
  static bool myNe( const CV& x, const CV& y ) { return mc::Op<T>::ne(const_cast<CV*>(&x)->bound(),const_cast<CV*>(&y)->bound()); }
  static bool myLt( const CV& x, const CV& y ) { return mc::Op<T>::lt(const_cast<CV*>(&x)->bound(),const_cast<CV*>(&y)->bound()); }
  static bool myLe( const CV& x, const CV& y ) { return mc::Op<T>::le(const_cast<CV*>(&x)->bound(),const_cast<CV*>(&y)->bound()); }
  static bool myGt( const CV& x, const CV& y ) { return mc::Op<T>::gt(const_cast<CV*>(&x)->bound(),const_cast<CV*>(&y)->bound()); }
  static bool myGe( const CV& x, const CV& y ) { return mc::Op<T>::ge(const_cast<CV*>(&x)->bound(),const_cast<CV*>(&y)->bound()); }
};

} // end namespace fadbad


//#include "mcop.hpp"

namespace mc
{

//! @brief C++ structure for specialization of the mc::Op templated structure for use of mc::SICVar in DAG evaluation and as template parameter in other MC++ types
template< typename T > struct Op< mc::SICVar<T> >
{
  typedef mc::SICVar<T> CV;
  static CV point( double const c ) { return CV(c); }
  static CV zeroone() { return CV( mc::Op<T>::zeroone() ); }
  static void I(CV& x, const CV&y) { x = y; }
  static double l(const CV& x) { return mc::Op<T>::l(x.B()); }
  static double u(const CV& x) { return mc::Op<T>::u(x.B()); }
  static double abs (const CV& x) { return mc::Op<T>::abs(x.B());  }
  static double mid (const CV& x) { return mc::Op<T>::mid(x.B());  }
  static double diam(const CV& x) { return mc::Op<T>::diam(x.B()); }
  static CV inv (const CV& x) { return mc::inv(x);  }
  static CV sqr (const CV& x) { return mc::sqr(x);  }
  static CV sqrt(const CV& x) { return mc::sqrt(x); }
  static CV log (const CV& x) { return mc::log(x);  }
  static CV xlog(const CV& x) { return x*mc::log(x); }
  static CV lmtd(const CV& x, const CV& y) { return (x-y)/(mc::log(x)-mc::log(y)); }
  static CV rlmtd(const CV& x, const CV& y) { return (mc::log(x)-mc::log(y))/(x-y); }
  static CV fabs(const CV& x) { return mc::fabs(x); }
  static CV exp (const CV& x) { return mc::exp(x);  }
  static CV sin (const CV& x) { return mc::sin(x);  }
  static CV cos (const CV& x) { return mc::cos(x);  }
  static CV tan (const CV& x) { return mc::tan(x);  }
  static CV asin(const CV& x) { return mc::asin(x); }
  static CV acos(const CV& x) { return mc::acos(x); }
  static CV atan(const CV& x) { return mc::atan(x); }
  static CV sinh(const CV& x) { return mc::sinh(x); }
  static CV cosh(const CV& x) { return mc::cosh(x); }
  static CV tanh(const CV& x) { return mc::tanh(x); }
  static CV erf (const CV& x) { throw typename mc::SICModel<T>::Exceptions( SICModel<T>::Exceptions::UNDEF ); }
  static CV erfc(const CV& x) { throw typename mc::SICModel<T>::Exceptions( SICModel<T>::Exceptions::UNDEF ); }
  static CV fstep(const CV& x) { return CV( mc::Op<T>::fstep(x.B()) ); }
  static CV bstep(const CV& x) { return CV( mc::Op<T>::bstep(x.B()) ); }
  static CV hull(const CV& x, const CV& y) { return mc::hull(x,y); }
  static CV min (const CV& x, const CV& y) { return mc::Op<T>::min(x.B(),y.B());  }
  static CV max (const CV& x, const CV& y) { return mc::Op<T>::max(x.B(),y.B());  }
  static CV arh (const CV& x, double const k) { return mc::exp(-k/x); }
  template <typename X, typename Y> static CV pow(const X& x, const Y& y) { return mc::pow(x,y); }
  static CV cheb(const CV& x, const unsigned n) { return mc::cheb(x,n); }
  static CV prod (const unsigned n, const CV* x) { return mc::prod(n,x); }
  static CV monom (const unsigned n, const CV* x, const unsigned* k) { return mc::monom(n,x,k); }
  static bool inter(CV& xIy, const CV& x, const CV& y) { return mc::inter(xIy,x,y); }
  static bool eq(const CV& x, const CV& y) { return mc::Op<T>::eq(x.B(),y.B()); }
  static bool ne(const CV& x, const CV& y) { return mc::Op<T>::ne(x.B(),y.B()); }
  static bool lt(const CV& x, const CV& y) { return mc::Op<T>::lt(x.B(),y.B()); }
  static bool le(const CV& x, const CV& y) { return mc::Op<T>::le(x.B(),y.B()); }
  static bool gt(const CV& x, const CV& y) { return mc::Op<T>::gt(x.B(),y.B()); }
  static bool ge(const CV& x, const CV& y) { return mc::Op<T>::ge(x.B(),y.B()); }
};

} // namespace mc

#endif

