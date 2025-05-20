// Copyright (C) Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the Eclipse Public License.

/*!
\page page_SUPREL Superposition Relaxation Arithmetic for Factorable Functions
\author Yanlin Zha, Beno&icirc;t Chachuat

A superposition relaxation is any estimator that brackets the range of an estimated function between an underestimator and an overestimator that are both separable functions. Formally, given the function \f$f:\mathcal{X}\to\mathbb{R}\f$ defined on the domain \f$\mathcal{X}\subseteq \mathbb{R}^n\f$, a <I>superposition relaxation</I> of \f$f\f$ on \f${\bf X}\subseteq \mathcal{X}\f$, \f${\bf X}\in\mathbb{IR}^n\f$ is any pair of separable functions \f$f^{\rm u},f^{\rm o}:{\bf X}\to \mathbb{R}\f$ that bracket the variations of \f$f\f$ on \f${\bf X}\f$,
\f{align*}
\forall {\bf x}\in{\bf X}, \quad f^{\rm u}({\bf x}) =: \sum_{i=1}^n f^{\rm u}_i(x_i) \leq f({\bf x}) \leq f^{\rm o}({\bf x}) =: \sum_{i=1}^n f^{\rm o}_i(x_i),
\f}
for some univariate real-valued functions \f$f^{\rm u}_i,f^{\rm o}_i:X_i\to\mathbb{R}\f$. We thus call the separable function \f$f^{\rm u}\f$ a </I>superposition underestimator</I> and \f$f^{\rm o}\f$, a <I>superposition overestimator</I>.

The classes mc::SupModel and mc::SupVar provide an implementation of superposition relaxation arithmetic based on the operator/function overloading mechanism of C++. This makes superposition relaxations both simple and intuitive to compute, similar to computing function values in real arithmetics or function bounds in interval arithmetic (see \ref page_INTERVAL). mc::SupModel and mc::SupVar are templated in the type used to propagate the univariate (under- and over-)estimators. Available types include mc::PWLU for piecewise linear univariate estimators (on adaptive partitions) and mc:PWCU for piecewise constant univariate estimators (on grid partitions). In addition, mc::SupVar can be used as the template parameter of other available types in MC++; as well as types in <A href="http://www.fadbad.com/fadbad.html">FADBAD++</A> for computing superposition relaxations for a factorable function's derivatives or Taylor coefficients.

- \subpage page_PWCU
- \subpage page_PWLU

\section sec_SUPREL_use How do I compute a superposition relaxation of a factorable function?

Suppose we want to compute a superposition relaxation for the real-valued function \f$f(x,y)=x\exp(x+y^2)-y^2\f$ with \f$(x,y)\in [1,2]\times[0,1]\f$, using piecewise linear univariate estimators:

\code
      #include "interval.hpp"
      #include "pwlu.hpp"
      #include "supmodel.hpp"
      typedef mc::SupModel<mc::PWLU> SM;
      typedef mc::SupVar<mc::PWLU> SV;
\endcode

First, the number of independent variables in the factorable function (\f$x\f$ and \f$y\f$ here) is specified by defining the superposition model mc::SupModel:

\code
      SM mod( 2 );
\endcode

Next, the variables \f$x\f$ and \f$y\f$ are defined as follows:

\code
      SV X( mod, 0, mc::Interval(1.,2.), 8 );
      SV Y( mod, 1, mc::Interval(0.,1.), 8 );
\endcode

Essentially, the first line means that <tt>X</tt> is a variable of class mc::SupVar, participating in the superposition model <tt>mod</tt>, with index 0 (using C-style index) and range \f$[1,2]\f$ with 8 equal-size partition. The same holds for the variable <tt>Y</tt>, participating in the superposition model <tt>mod</tt>, with index 1 and range \f$[0,1]\f$ with 8 equal-size partition.

Having defined both variables, a superposition relaxation of \f$f(x,y)=x\exp(x+y^2)-y^2\f$ on \f$[1,2]\times[0,1]\f$ is simply computed as:

\code
      SV F = X*exp(X+pow(Y,2))-pow(Y,2);
\endcode

This model can be displayed to the standard output as:

\code
      std::cout << "Superposition relaxation of f:\n" << F;
\endcode

which produces the following output:

\verbatim
Superposition relaxation of f:
   U[0]: {  1.00000e+00:  3.74398e+00,   1.03190e+00:  3.97165e+00,   1.06380e+00:  4.22296e+00,   1.09440e+00:  4.45864e+00,
            1.12500e+00:  4.71687e+00,   1.15690e+00:  4.97005e+00,   1.18880e+00:  5.24773e+00,   1.21940e+00:  5.51272e+00,
            1.25000e+00:  5.80117e+00,   1.28190e+00:  6.08677e+00,   1.31380e+00:  6.39786e+00,   1.34440e+00:  6.69967e+00,
            1.37500e+00:  7.02598e+00,   1.40690e+00:  7.35205e+00,   1.43880e+00:  7.70477e+00,   1.46940e+00:  8.05211e+00,
            1.50000e+00:  8.42520e+00,   1.53190e+00:  8.80111e+00,   1.56380e+00:  9.20500e+00,   1.59440e+00:  9.60806e+00, 
            1.62500e+00:  1.00383e+01,   1.65690e+00:  1.04749e+01,   1.68880e+00:  1.09411e+01,   1.71940e+00:  1.14117e+01,
            1.75000e+00:  1.19111e+01,   1.78190e+00:  1.24211e+01,   1.81380e+00:  1.29625e+01,   1.84440e+00:  1.35143e+01,
            1.87500e+00:  1.40969e+01,   1.90690e+00:  1.46950e+01,   1.93880e+00:  1.53266e+01,   1.96940e+00:  1.59755e+01,
            2.00000e+00:  1.66575e+01 }
   U[1]: {  0.00000e+00: -1.87935e+00,   6.25000e-02: -1.88716e+00,   7.81657e-02: -1.87460e+00,   9.38314e-02: -1.86202e+00,
            1.09416e-01: -1.84931e+00,   1.25000e-01: -1.83659e+00,   1.40666e-01: -1.82775e+00,   1.56331e-01: -1.81890e+00, 
            1.71916e-01: -1.80989e+00,   1.87500e-01: -1.80088e+00,   2.03206e-01: -1.77689e+00,   2.18913e-01: -1.75288e+00,
            2.34456e-01: -1.72828e+00,   2.50000e-01: -1.70366e+00,   2.65706e-01: -1.68281e+00,   2.81413e-01: -1.66193e+00,
            2.96956e-01: -1.64042e+00,   3.12500e-01: -1.61888e+00,   3.28247e-01: -1.58145e+00,   3.43994e-01: -1.54394e+00,
            3.59497e-01: -1.50506e+00,   3.75000e-01: -1.46610e+00,   3.90747e-01: -1.43075e+00,   4.06494e-01: -1.39533e+00,
            4.21997e-01: -1.35843e+00,   4.37500e-01: -1.32146e+00,   4.53288e-01: -1.26714e+00,   4.69075e-01: -1.21267e+00,
            4.84538e-01: -1.15562e+00,   5.00000e-01: -1.09840e+00,   5.15788e-01: -1.04457e+00,   5.31575e-01: -9.90566e-01,
            5.47038e-01: -9.33832e-01,   5.62500e-01: -8.76908e-01,   5.78328e-01: -8.00647e-01,   5.94157e-01: -7.24075e-01,
            6.09578e-01: -6.43251e-01,   6.25000e-01: -5.62081e-01,   6.40828e-01: -4.84063e-01,   6.56657e-01: -4.05680e-01,
            6.72078e-01: -3.22847e-01,   6.87500e-01: -2.39610e-01,   7.03369e-01: -1.34535e-01,   7.19238e-01: -2.88435e-02,
            7.34619e-01:  8.32278e-02,   7.50000e-01:  1.95997e-01,   7.65869e-01:  3.05659e-01,   7.81738e-01:  4.16064e-01,
            7.97119e-01:  5.33013e-01,   8.12500e-01:  6.50803e-01,   8.28410e-01:  7.93092e-01,   8.44320e-01:  9.36608e-01,
            8.59660e-01:  1.08877e+00,   8.75000e-01:  1.24236e+00,   8.90910e-01:  1.39210e+00,   9.06820e-01:  1.54336e+00,
            9.22160e-01:  1.70306e+00,   9.37500e-01:  1.86453e+00,   9.53450e-01:  2.05214e+00,   9.69401e-01:  2.24225e+00,
            9.84700e-01:  2.44219e+00,   1.00000e+00:  2.64508e+00 }
   O[0]: {  1.00000e+00:  4.09862e+00,   1.06250e+00:  4.74738e+00,   1.12500e+00:  5.35028e+00,   1.18750e+00:  6.13348e+00,
            1.25000e+00:  6.87627e+00,   1.31250e+00:  7.83686e+00,   1.37500e+00:  8.76356e+00,   1.43750e+00:  9.96022e+00,
            1.50000e+00:  1.11305e+01,   1.56250e+00:  1.26433e+01,   1.62500e+00:  1.41381e+01,   1.68750e+00:  1.60764e+01,
            1.75000e+00:  1.80049e+01,   1.81250e+00:  2.05179e+01,   1.87500e+00:  2.30281e+01,   1.93750e+00:  2.63198e+01,
            2.00000e+00:  2.96114e+01 }
   O[1]: {  0.00000e+00: -1.38034e+00,   6.25000e-02: -1.35872e+00,   6.25000e-02: -1.35872e+00,   1.25000e-01: -1.35274e+00,
            1.87500e-01: -1.29865e+00,   2.50000e-01: -1.26046e+00,   3.12500e-01: -1.15749e+00,   3.12500e-01: -1.15749e+00,
            3.75000e-01: -1.07112e+00,   4.37500e-01: -8.83886e-01,   5.00000e-01: -7.15065e-01,   5.62500e-01: -3.70708e-01,
            6.25000e-01: -4.96009e-02,   6.87500e-01:  6.08522e-01,   7.50000e-01:  1.22967e+00,   8.12500e-01:  2.56792e+00,
            8.75000e-01:  3.82717e+00,   9.37500e-01:  6.80287e+00,   1.00000e+00:  9.55965e+00 }
   RNG:  [  1.85682e+00 :  3.91711e+01 ]
\endverbatim

A plot of the resulting superposition relaxation for the function \f$f\f$ is shown on the left plot in the figure below. Similarly, a superposition relaxation could be computed using piecewise constant univariate estimators (class mc::PWCU available from the header file pwcu.hpp), now producing the following output and illustrated on the right plot below.

\verbatim
Superposition relaxation of f:
   U[0]: { 1.00000e+00 <  2.97039e+00 :  5.35192e+00 >  1.12500e+00 <  3.78751e+00 :  6.58146e+00 >  1.25000e+00 <  4.72313e+00 :  7.94976e+00 >
           1.37500e+00 <  5.80664e+00 :  9.48985e+00 >  1.50000e+00 <  7.07226e+00 :  1.12402e+01 >  1.62500e+00 <  8.55980e+00 :  1.32454e+01 >
           1.75000e+00 <  1.03156e+01 :  1.55576e+01 >  1.87500e+00 <  1.23933e+01 :  1.82370e+01 >  2.00000e+00 }
   U[1]: { 0.00000e+00 < -1.87788e+00 : -1.80281e+00 >  1.25000e-01 < -1.86531e+00 : -1.63542e+00 >  2.50000e-01 < -1.76042e+00 : -1.36126e+00 > 
           3.75000e-01 < -1.54876e+00 : -9.54894e-01 >  5.00000e-01 < -1.20489e+00 : -3.78288e-01 >  6.25000e-01 < -6.90788e-01 :  4.20299e-01 >
           7.50000e-01 <  4.52993e-02 :  1.50367e+00 >  8.75000e-01 <  1.06617e+00 :  2.93059e+00 >  1.00000e+00 }
   O[0]: { 1.00000e+00 <  3.59737e+00 :  5.70625e+00 >  1.12500e+00 <  4.73747e+00 :  7.33694e+00 >  1.25000e+00 <  6.15876e+00 :  9.32013e+00 >
           1.37500e+00 <  7.95014e+00 :  1.17717e+01 >  1.50000e+00 <  1.02325e+01 :  1.48493e+01 >  1.62500e+00 <  1.31700e+01 :  1.87677e+01 >
           1.75000e+00 <  1.69853e+01 :  2.38184e+01 >  1.87500e+00 <  2.19809e+01 :  3.03986e+01 >  2.00000e+00 }
   O[1]: { 0.00000e+00 < -1.23888e+00 : -1.18051e+00 >  1.25000e-01 < -1.24301e+00 : -1.05851e+00 >  2.50000e-01 < -1.18351e+00 : -8.40816e-01 >
           3.75000e-01 < -1.02832e+00 : -4.58402e-01 >  5.00000e-01 < -7.08402e-01 :  2.30224e-01 >  6.25000e-01 < -8.22759e-02 :  1.52721e+00 >
           7.50000e-01 <  1.15221e+00 :  4.13264e+00 >  8.75000e-01 <  3.69514e+00 :  9.85460e+00 >  1.00000e+00 }
   RNG:  [  1.09251e+00 :  4.02532e+01 ]
\endverbatim

<CENTER><TABLE BORDER=0>
<TR>
<TD>\image html PWLSM-2D.png</TD>
<TD>\image html PWCSM-2D.png</TD>
</TR>
</TABLE></CENTER>

Other operations involve retreiving the piecewise linear univariate under- and over-estimators for each participating variable:

\code
      std::vector<mc::PWLU> const& Fuest = F.uest();
      std::vector<mc::PWLU> const& Foest = F.oest();
\endcode

or retrieving lower and upper bounds from the superposition relaxation:

\code
      double const& Fl = F.l();
      double const& Fu = F.u();
\endcode

The piecewise linear univariate under- and over-estimators can also be evaluated at given points:

\code
      std::cout << "uest(1.5,0.5) = " << F.uval({{0,1.5},{1,0.5}}) << std::endl;
      std::cout << "oest(1.5,0.5) = " << F.oval({{0,1.5},{1,0.5}}) << std::endl;
\endcode

\verbatim
  Suposition relaxation of f:
  uest(1.5,0.5) = 7.32679e+00
  oest(1.5,0.5) = 1.04154e+01
\endverbatim


See the documentations of mc::SupModel and mc::SupVar for a complete list of member functions. 


\section sec_SUPREL_err Errors Errors encountered during computation of a superposition relaxation

Errors are managed based on the exception handling mechanism of the C++ language. Each time an error is encountered, a class object of type mc::SupModel::Exceptions is thrown, which contains the type of error. It is the user's responsibility to test whether an exception was thrown during the computation of a superposition relaxation, and then make the appropriate changes. Should an exception be thrown and not caught by the calling program, the execution will abort.

Further exceptions may be thrown by the template class itself.


\section sec_SUPREL_refs References

- Y. Zha, M.E. Villanueva, B. Houska. <A HREF="https://arxiv.org/abs/1610.05862">Interval superposition arithmetic</A>, <I>ArXiv</I>, 1610.05862v2, 2018.
.
*/

#ifndef MC__SUPMODEL_HPP
#define MC__SUPMODEL_HPP

#include <iostream>
#include <iomanip> 
#include <vector>
#include <set>
#include <map>
#include <algorithm>
#include <cmath>
#include <limits>
#include <bitset>
#include <cassert>

#include "mcop.hpp"
#include "mcfunc.hpp"

//#define MC__SUPMODEL_TRACE

namespace mc
{

template <typename Summand> class SupVar;

//! @brief C++ class for superposition models of factorable functions: environment
////////////////////////////////////////////////////////////////////////
//! mc::SupModel is a C++ class for definition of the environment of a
//! superposition model. Superposition model propagation in factorable
//! functions is implemented via the C++ class mc::SupVar. The template
//! parameters correspond to 1) the type of the range enclosures and
//! 2) the type of the summands.
////////////////////////////////////////////////////////////////////////
template <typename Summand>
class SupModel
////////////////////////////////////////////////////////////////////////
{
  template <typename S> friend std::ostream& operator<<
    ( std::ostream&, const SupModel<S>& );

  template <typename S> friend class SupVar;

 private:

  //! @brief Number of variables
  unsigned int _nvar;
  //! @brief Whether variables are defined or not
  std::vector<bool> _defvar;
  //! @brief Variable lower bounds
  std::vector<double> _lbdvar;
  //! @brief Variable upper bounds
  std::vector<double> _ubdvar;

 public:

  //! @brief Constructor of superposition model with <a>nvar</a> variables
  SupModel
  ( size_t const& nvar )
  : _nvar(nvar)
  {
    _lbdvar.resize( _nvar );
    _ubdvar.resize( _nvar );
    _defvar.resize( _nvar, false );
  }

  ~SupModel() 
  {}

  //! @brief Retrieve number of variables in superposition model
  unsigned nvar
  ()
  const
  { return _nvar; };

  //! @brief Retrieve variable lower bounds in superposition model
  std::vector<double> lbdvar
  ()
  const
  { return _lbdvar; };

  //! @brief Retrieve variable upper bounds in superposition model
  std::vector<double> ubdvar
  ()
  const
  { return _ubdvar; };

 //! @brief Options of mc::SupModel
  struct Options
  {
    //! @brief Constructor
    Options():
      PROD_METH        ( PARTIAL ),
      PROD_CUT         ( 0 ),
      SUM_TOL          ( 1e2*DBL_EPSILON ),
      REF_WEIGHT       ( 0.5 ),
      MAX_SUBDIV       ( 0 ),
      USE_SHADOW       ( false ),
      DISPLAY_SHADOW   ( true ),
      DISPLAY_DIGITS   ( 5 )
      {}
    //! @brief Assign options
    template <typename OPTIONS>
    Options& operator=
      ( OPTIONS const& opt )
      {
        PROD_METH         = opt.PROD_METH;
        PROD_CUT          = opt.PROD_CUT;
        SUM_TOL           = opt.SUM_TOL;
        REF_WEIGHT        = opt.REF_WEIGHT;
        MAX_SUBDIV        = opt.MAX_SUBDIV;
        USE_SHADOW        = opt.USE_SHADOW;
        DISPLAY_SHADOW    = opt.DISPLAY_SHADOW;
        DISPLAY_DIGITS    = opt.DISPLAY_DIGITS;
        return *this;
      }
    //! @brief Reset options
    void reset
      ()
      {
        PROD_METH         = PARTIAL;
        PROD_CUT          = 0;
        SUM_TOL           = 1e2*DBL_EPSILON;
        REF_WEIGHT        = 0.5,
        MAX_SUBDIV        = 0;
        USE_SHADOW        = false;
        DISPLAY_SHADOW    = true;
        DISPLAY_DIGITS    = 5;
      }

    //! @brief Reformulation method for product term
    enum PROD_REF{
      NONE=0,	//!< DC decomposition w/o rescaling
      PARTIAL,	//!< DC decomposition w/ range rescaling
      FULL,	//!< DC decomposition w/ range and midpoint rescaling
      LOG	//!< log-transform w/ range and midpoint rescaling
    };
    //! @brief Reformulation method used for product terms - Default: PARTIAL
    int PROD_METH; 
    //! @brief Whether to cut superposition relaxations for product terms - Default: false
    bool PROD_CUT; 
    //! @brief Tolerance on range for univariate estimator propagation - Default: 1e2*DBL_EPSILON
    double SUM_TOL;
    //! @brief Weight in \f$[0,1]\f$ used in the overloaded functions mc::min and mc::max, where a value of 0 amounts to choosing the first operand as reference, a value >=1 the second operand, or a weighted intermediate between the first and second operands - Default: 0.5
    double REF_WEIGHT;
    //! @brief Maximal number of subdivisions, for univariate estimators on adaptive grids - Default: 0 (no restriction)
    size_t MAX_SUBDIV; 
    //! @brief Whether to enable shadow estimators - Default: false
    bool USE_SHADOW; 
    //! @brief Whether to display shadow estimators - Default: true
    bool DISPLAY_SHADOW; 
    //! @brief Number of digits displayed with << operator - Default: 5   
    unsigned int DISPLAY_DIGITS;
  } options;

  //! @brief Exceptions of mc::SupModel
  class Exceptions
  {
   public:
    //! @brief Enumeration type for SCModel exception handling
    enum TYPE{
      DIV=1,	  //!< Division by zero scalar
      INV,	  //!< Inverse operation with zero in range
      LOG,	  //!< Log operation with non-positive numbers in range
      SQRT,	  //!< Square-root operation with negative numbers in range
      RPOW,	  //!< real power operation with negative numbers in range
      TAN,	  //!< Tangent operation with (k+1/2)·PI in range
      INTERN=-1,  //!< Internal error
      INDEX=-2,   //!< Variable index out of range
      MODEL=-3,	  //!< Operation between variables linked to different models
      UNDEF=-33   //!< Feature not yet implemented
    };
    //! @brief Constructor for error <a>ierr</a>
    Exceptions( TYPE ierr ) : _ierr( ierr ){}
    //! @brief Error flag
    int ierr(){ return _ierr; }
    //! @brief Error description
    std::string what(){
      switch( _ierr ){
      case DIV:
        return "mc::SupModel\t Division by zero scalar";
      case INV:
        return "mc::SupModel\t Inverse operation with zero in range";
      case LOG:
        return "mc::SupModel\t Log operation with non-positive numbers in range";
      case SQRT:
        return "mc::SupModel\t Square-root operation with negative numbers in range";
      case RPOW:
        return "mc::SupModel\t Real power operation with negative numbers in range";
      case TAN:
        return "mc::SupModel\t Tangent operation with (k+1/2)·PI in range";
      case INDEX:
        return "mc::SupModel\t Variable index out of range";
      case MODEL:
        return "mc::SupModel\t Operation between variables linked to different superposition models";
      case UNDEF:
        return "mc::SupModel\t Feature not yet implemented";
      case INTERN:
      default:
        return "mc::SupModel\t Internal error";
      }
    }
   private:
    TYPE _ierr;
  };

  //! @brief Compose variable <a>var</a> with a univariate outer function <a>f</a> with derivative <a>Df</a> and given convexity <a>cvx</a> and monotonicity <a>inc</a>
  template <typename UNIV, typename DUNIV>
  void compose
    ( SupVar<Summand>& var, UNIV const& f, DUNIV const& Df, bool const cvx, bool const inc )
    const
    {
      // No need to propagate shadow estimators in single dependency case
      if( options.USE_SHADOW && var._sdep.size() > 1 ){
        if( var._uest2.empty() && ( ( cvx &&  inc ) || ( !cvx && !inc ) ) )
          var._uest2 = var._uest; // duplicate current underestimator
        if( !var._uest2.empty() ) SupVar<Summand>::_ub( _ubd2, var._sdep, var._uest2 );
        if( var._oest2.empty() && ( ( cvx && !inc ) || ( !cvx &&  inc ) ) )
          var._oest2 = var._oest; // duplicate current overestimator
        if( !var._oest2.empty() ) SupVar<Summand>::_lb( _lbd2, var._sdep, var._oest2 );
        _compose( var._sdep, var._uest2, var._oest2, _ubd2, _lbd2, 1, f, Df, cvx, inc );
        var._lbd2.second = var._ubd2.second = false; // reset shadow estimator bounds
      }

      _compose( var._sdep, var._uest, var._oest, var.l(0), var.u(0), 0, f, Df, cvx, inc );
      var._lbd.second = var._ubd.second = false; // reset primal estimator bounds
    }

  //! @brief Intersect variable <a>var</a> with upper bound value <a>val</a>
  void max
    ( SupVar<Summand>& var, double const& val )
    const
    {
      // No need to propagate shadow estimators in single dependency case
      if( options.USE_SHADOW && var._sdep.size() > 1 ){
        if( var._uest2.empty() )
          var._uest2 = var._uest; // duplicate current underestimator
        SupVar<Summand>::_ub( _ubd2, var._sdep, var._uest2 );
        if( !var._oest2.empty() )
          SupVar<Summand>::_lb( _lbd2, var._sdep, var._oest2 );
        _max( var._sdep, var._uest2, var._oest2, _ubd2, _lbd2, 1, val );
        var._lbd2.second = var._ubd2.second = false; // reset shadow estimator bounds
      }

      _max( var._sdep, var._uest, var._oest, var.l(0), var.u(0), 0, val );
      var._lbd.second = var._ubd.second = false; // reset primal estimator bounds
    }

  //! @brief Intersect variable <a>var</a> with lower bound value <a>val</a>
  void min
    ( SupVar<Summand>& var, double const& val )
    const
    {
      // No need to propagate shadow estimators in single dependency case
      if( options.USE_SHADOW && var._sdep.size() > 1 ){
        if( !var._uest2.empty() )
          SupVar<Summand>::_ub( _ubd2, var._sdep, var._uest2 );
        if( var._oest2.empty() )
          var._oest2 = var._oest; // duplicate current overestimator
        SupVar<Summand>::_lb( _lbd2, var._sdep, var._oest2 );
        _min( var._sdep, var._uest2, var._oest2, _ubd2, _lbd2, 1, val );
        var._lbd2.second = var._ubd2.second = false; // reset shadow estimator bounds
      }

      _min( var._sdep, var._uest, var._oest, var.l(0), var.u(0), 0, val );
      var._lbd.second = var._ubd.second = false; // reset primal estimator bounds
    }

  //! @brief Intersect underestimator in variable <a>var</a> with the natural lower bound value <a>val</a>
  void uref
    ( SupVar<Summand>& var, double const& val )
    const
    {
      static std::vector<Summand> est0 = std::vector<Summand>();

      // No need to propagate shadow estimators in single dependency case
      if( options.USE_SHADOW && var._sdep.size() > 1 ){
        if( var._uest2.empty() )
          var._uest2 = var._uest; // duplicate current underestimator
        SupVar<Summand>::_ub( _ubd2, var._sdep, var._uest2 );
        _max( var._sdep, var._uest2, est0, _ubd2, _lbd2, 1, val );
        var._lbd2.second = false; // reset shadow underestimator bounds
      }

      _max( var._sdep, var._uest, est0, var.l(0), 0., 0, val );
      var._lbd.second = false; // reset primal underestimator bounds
    }

  //! @brief Intersect overestimator in variable <a>var</a> with the natural upper bound value <a>val</a>
  void oref
    ( SupVar<Summand>& var, double const& val )
    const
    {
      static std::vector<Summand> est0 = std::vector<Summand>();

      // No need to propagate shadow estimators in single dependency case
      if( options.USE_SHADOW && var._sdep.size() > 1 ){
        if( var._oest2.empty() )
          var._oest2 = var._oest; // duplicate current overestimator
        SupVar<Summand>::_lb( _lbd2, var._sdep, var._oest2 );
        _min( var._sdep, est0, var._oest2, _ubd2, _lbd2, 1, val );
        var._ubd2.second = false; // reset shadow overestimator bounds
      }

      _min( var._sdep, est0, var._oest, 0., var.u(0), 0, val );
      var._ubd.second = false; // reset primal overestimator bounds
    }

 private:

  //! @brief Intermediate storage
  mutable std::vector<Summand> _est3;
  mutable std::vector<Summand> _est4;
  mutable double _lbd2;
  mutable double _ubd2;
  mutable double _lbd3;
  mutable double _ubd3;
  mutable double _lbd4;
  mutable double _ubd4;

  void _max
    ( std::set<unsigned> const& sdep, std::vector<Summand>& uest, std::vector<Summand>& oest,
       double const& ubnd, double const& obnd, bool const shadow, double const& val )
    const
    {
      auto const& fmax = [=]( const double& x ){ return x > val? x: val; };

      // Single dependency case
      if( sdep.size() == 1 ){
        auto const& i = *sdep.cbegin();

        if( !uest.empty() && uest[i].w() > options.SUM_TOL )
          uest[i].max( val ).reduce( 1, options.MAX_SUBDIV );
        else if( !uest.empty() )
          uest[i].set( fmax( uest[i].l() ) );
          
        if( !oest.empty() && oest[i].w() > options.SUM_TOL )
          oest[i].max( val ).reduce( 0, options.MAX_SUBDIV );
        else if( !oest.empty() )
          oest[i].set( fmax( oest[i].u() ) );

        return;
      }

      // Multiple dependency case
      double wu( 0. ), wo( 0. );
      for( auto const& i: sdep ){
        if( !_defvar[i] )
          throw Exceptions( Exceptions::INTERN );

        if( !uest.empty() && uest[i].w() > options.SUM_TOL )
          wu += uest[i].w();

        if( !oest.empty() && oest[i].w() > options.SUM_TOL )
          wo += oest[i].w();
      }
        
      for( auto const& i: sdep ){

        if( !uest.empty() && uest[i].w() > options.SUM_TOL ){
          double const qu = uest[i].w() / wu;
          uest[i] += ubnd - (shadow? uest[i].u(): uest[i].l());
          uest[i].max( val ).reduce( 1, options.MAX_SUBDIV );
          uest[i] -= ( 1 - qu ) * fmax( ubnd );
        }
        else if( !uest.empty() )
          uest[i].set( fmax( uest[i].l() ) );
          
        if( !oest.empty() && oest[i].w() > options.SUM_TOL ){
          double const qo = oest[i].w() / wo;
          oest[i] -= (shadow? oest[i].l(): oest[i].u()); //oest[i].u() 
          oest[i] /= qo;
          oest[i] += obnd;
          oest[i].max( val ).reduce( 0, options.MAX_SUBDIV );
          oest[i] *= qo;
        }
        else if( !oest.empty() )
          oest[i].set( fmax( oest[i].u() ) );
      }
    }

  void _min
    ( std::set<unsigned> const& sdep, std::vector<Summand>& uest, std::vector<Summand>& oest,
       double const& ubnd, double const& obnd, bool const shadow, double const& val )
    const
    {
      auto const& fmin = [=]( const double& x ){ return x < val? x: val; };

      // Single dependency case
      if( sdep.size() == 1 ){
        auto const& i = *sdep.cbegin();

        if( !oest.empty() && oest[i].w() > options.SUM_TOL )
          oest[i].min( val ).reduce( 0, options.MAX_SUBDIV );
        else if( !oest.empty() )
          oest[i].set( fmin( oest[i].u() ) );

        if( !uest.empty() && uest[i].w() > options.SUM_TOL )
          uest[i].min( val ).reduce( 1, options.MAX_SUBDIV );
        else if( !uest.empty() )
          uest[i].set( fmin( uest[i].l() ) );

        return;
      }
      
      // Multiple dependency case
      double wu( 0. ), wo( 0. );
      for( auto const& i: sdep ){
        if( !_defvar[i] )
          throw Exceptions( Exceptions::INTERN );

        if( !uest.empty() && uest[i].w() > options.SUM_TOL )
          wu += uest[i].w();

        if( !oest.empty() && oest[i].w() > options.SUM_TOL )
          wo += oest[i].w();
      }
        
      for( auto const& i: sdep ){

        if( !oest.empty() && oest[i].w() > options.SUM_TOL ){
          double const qo = oest[i].w() / wo;
          oest[i] += obnd - (shadow? oest[i].l(): oest[i].u());
          oest[i].min( val ).reduce( 0, options.MAX_SUBDIV );
          oest[i] -= ( 1 - qo ) * fmin( obnd );
        }
        else if( !oest.empty() )
          oest[i].set( fmin( oest[i].u() ) );

        if( !uest.empty() && uest[i].w() > options.SUM_TOL ){
          double const qu = uest[i].w() / wu;
          uest[i] -= (shadow? uest[i].u(): uest[i].l()); //uest[i].l();
          uest[i] /= qu;
          uest[i] += ubnd;
          uest[i].min( val ).reduce( 1, options.MAX_SUBDIV );
          uest[i] *= qu;
        }
        else if( !uest.empty() )
          uest[i].set( fmin( uest[i].l() ) );
      }
    }
    
  template <typename UNIV, typename DUNIV>
  void _compose
    ( std::set<unsigned> const& sdep, std::vector<Summand>& uest, std::vector<Summand>& oest,
       double const& ubnd, double const& obnd, bool const shadow,
       UNIV const& f, DUNIV const& Df, bool const cvx, bool const inc )
    const
    {
      // Single dependency case
      if( sdep.size() == 1 ){
        auto const& i = *sdep.cbegin();

        if( ( cvx && inc ) || ( !cvx && !inc ) ){
          if( !uest.empty() && uest[i].w() > options.SUM_TOL )
            uest[i].compose( f, Df, cvx?1:0, cvx, inc ).reduce( cvx?1:0, options.MAX_SUBDIV );
          else if( !uest.empty() )
            uest[i].set( f( uest[i].l() ) );

          if( !oest.empty() && oest[i].w() > options.SUM_TOL )
            oest[i].compose( f, Df, cvx?0:1, cvx, inc ).reduce( cvx?0:1, options.MAX_SUBDIV );
          else if( !oest.empty() )
            oest[i].set( f( oest[i].u() ) );
        }
        
        else{
          if( !oest.empty() && oest[i].w() > options.SUM_TOL )
            oest[i].compose( f, Df, cvx?1:0, cvx, inc ).reduce( cvx?1:0, options.MAX_SUBDIV );
          else if( !oest.empty() )
            oest[i].set( f( oest[i].u() ) );

          if( !uest.empty() && uest[i].w() > options.SUM_TOL )
            uest[i].compose( f, Df, cvx?0:1, cvx, inc ).reduce( cvx?0:1, options.MAX_SUBDIV );
          else if( !uest.empty() )
            uest[i].set( f( uest[i].l() ) );
        }
      
        if( !inc ) std::swap( uest, oest );
        return;
      }
      
      // Multiple dependency case
      double wu( 0. ), wo( 0. );
      for( auto const& i: sdep ){
        if( !_defvar[i] )
          throw Exceptions( Exceptions::INTERN );

        if( !uest.empty() && uest[i].w() > options.SUM_TOL )
          wu += uest[i].w();

        if( !oest.empty() && oest[i].w() > options.SUM_TOL )
          wo += oest[i].w();
      }
        
      for( auto const& i: sdep ){

        if( ( cvx && inc ) || ( !cvx && !inc ) ){
          if( !uest.empty() && uest[i].w() > options.SUM_TOL ){
            double const qu = uest[i].w() / wu;
            uest[i] += ubnd - (shadow? uest[i].u(): uest[i].l());
            uest[i].compose( f, Df, cvx?1:0, cvx, inc ).reduce( cvx?1:0, options.MAX_SUBDIV );
            uest[i] -= ( 1 - qu ) * f( ubnd );
          }
          else if( !uest.empty() ){
            uest[i].set( f( uest[i].l() ) );
          }
          
          if( !oest.empty() && oest[i].w() > options.SUM_TOL ){
            double const qo = oest[i].w() / wo;
            oest[i] -= (shadow? oest[i].l(): oest[i].u()); //oest[i].u() 
            oest[i] /= qo;
            oest[i] += obnd;
            oest[i].compose( f, Df, cvx?0:1, cvx, inc ).reduce( cvx?0:1, options.MAX_SUBDIV );
            oest[i] *= qo;
          }
          else if( !oest.empty() ){
            oest[i].set( f( oest[i].u() ) );
          }
        }
        
        else{
          if( !oest.empty() && oest[i].w() > options.SUM_TOL ){
            double const qo = oest[i].w() / wo;
            oest[i] += obnd - (shadow? oest[i].l(): oest[i].u());
            oest[i].compose( f, Df, cvx?1:0, cvx, inc ).reduce( cvx?1:0, options.MAX_SUBDIV );
            oest[i] -= ( 1 - qo ) * f( obnd );
          }
          else if( !oest.empty() ){
            oest[i].set( f( oest[i].u() ) );
          }
          
          if( !uest.empty() && uest[i].w() > options.SUM_TOL ){
            double const qu = uest[i].w() / wu;
            uest[i] -= (shadow? uest[i].u(): uest[i].l()); //uest[i].l();
            uest[i] /= qu;
            uest[i] += ubnd;
            uest[i].compose( f, Df, cvx?0:1, cvx, inc ).reduce( cvx?0:1, options.MAX_SUBDIV );
            uest[i] *= qu;
          }
          else if( !uest.empty() ){
            uest[i].set( f( uest[i].l() ) );
          }
        }
      }
      
      if( !inc ) std::swap( uest, oest );
    }
 
  void _add
    ( std::set<unsigned> const& sdep,  std::vector<Summand>& est, 
      std::set<unsigned> const& sdep2, std::vector<Summand> const& est2,
      bool const under )
    const
    {
      for( auto& i : sdep2 ){
        if( !sdep.count(i) ) est[i]  = est2[i];
        else                 est[i] += est2[i];
        est[i].reduce( under, options.MAX_SUBDIV );
      }
    }
};

template <typename Summand>
std::ostream& operator<<
( std::ostream& out, const SupModel<Summand>& mod )
{
  out << std::endl;
  out << "Superposition model:" << std::endl;
  out << "   " << "No. variables:  " << mod._nvar << std::endl;
  out << "   " << "Variable bounds: " << std::endl;
  for( unsigned int i=0; i<mod._nvar; i++ ){
    out << "       " << i << ":";
    if( mod._defvar[i] ){
      out << std::right << std::scientific << std::setprecision(mod.options.DISPLAY_DIGITS);
      out << "[ "  << std::setw(mod.options.DISPLAY_DIGITS+7) << mod._lbdvar[i]
          << " : " << std::setw(mod.options.DISPLAY_DIGITS+7) << mod._ubdvar[i] << " ]";
    }
    else
      out << "-";
    out << std::endl;
  }
  out << std::endl;
  return out;
}

//! @brief C++ class for superposition models of factorable functions: variable
////////////////////////////////////////////////////////////////////////
//! mc::SupModel is a C++ class for propagation of a superposition
//! model. The template parameters correspond to 1) the type of the
//! range enclosures and 2) the type of the summands.
////////////////////////////////////////////////////////////////////////
template <typename Summand> 
class SupVar
{
  template <typename S> friend class SupModel;

  template <typename S> friend std::ostream& operator<<
    ( std::ostream &, SupVar<S> const& );

  template <typename S> friend SupVar<S> operator+
    ( SupVar<S> const& );
  template <typename S> friend SupVar<S> && operator+
    ( SupVar<S> && );

  template <typename S> friend SupVar<S> operator+
    ( SupVar<S> const&, SupVar<S> const& );
  template <typename S> friend SupVar<S> && operator+
    ( SupVar<S> const&, SupVar<S> && );
  template <typename S> friend SupVar<S> && operator+
    ( SupVar<S> &&, SupVar<S> const& );    
  template <typename S> friend SupVar<S> && operator+
    ( SupVar<S> &&, SupVar<S> && );    

  template <typename S> friend SupVar<S> operator+
    ( double const&, SupVar<S> const& );
  template <typename S> friend SupVar<S> && operator+
    ( double const&, SupVar<S> && );

  template <typename S> friend SupVar<S> operator+
    ( SupVar<S> const&, double const& );
  template <typename S> friend SupVar<S> && operator+
    ( SupVar<S> &&, double const& );

  template <typename S> friend SupVar<S> operator-
    ( SupVar<S> const& );
  template <typename S> friend SupVar<S> && operator-
    ( SupVar<S> && );

  template <typename S> friend SupVar<S> operator-
    ( SupVar<S> const&, SupVar<S> const& );
  template <typename S> friend SupVar<S> && operator-
    ( SupVar<S> const&, SupVar<S> && );
  template <typename S> friend SupVar<S> && operator-
    ( SupVar<S> &&, SupVar<S> const& );
  template <typename S> friend SupVar<S> && operator-
    ( SupVar<S> &&, SupVar<S> && );

  template <typename S> friend SupVar<S> operator-
    ( double const&, SupVar<S> const& );
  template <typename S> friend SupVar<S> && operator-
    ( double const&, SupVar<S> && );

  template <typename S> friend SupVar<S> operator-
    ( SupVar<S> const&, double const& );
  template <typename S> friend SupVar<S> && operator-
    ( SupVar<S> &&, double const& );

  template <typename S> friend SupVar<S> operator*
    ( SupVar<S> const&, SupVar<S> const& );
  template <typename S> friend SupVar<S> && operator*
    ( SupVar<S> const&, SupVar<S> && );
  template <typename S> friend SupVar<S> && operator*
    ( SupVar<S> &&, SupVar<S> const& );
  template <typename S> friend SupVar<S> && operator*
    ( SupVar<S> &&, SupVar<S> && );

  template <typename S> friend SupVar<S> operator*
    ( double const&, SupVar<S> const& );
  template <typename S> friend SupVar<S> && operator*
    ( double const&, SupVar<S> && );

  template <typename S> friend SupVar<S> operator*
    ( SupVar<S> const&, double const& );
  template <typename S> friend SupVar<S> && operator*
    ( SupVar<S> &&, double const& );

  template <typename S> friend SupVar<S> operator/
    ( SupVar<S> const&, SupVar<S> const& );
  template <typename S> friend SupVar<S> && operator/
    ( SupVar<S> const&, SupVar<S> && );
  template <typename S> friend SupVar<S> && operator/
    ( SupVar<S> &&, SupVar<S> const& );
  template <typename S> friend SupVar<S> && operator/
    ( SupVar<S> &&, SupVar<S> && );

  template <typename S> friend SupVar<S> operator/
    ( double const&, SupVar<S> const& );
  template <typename S> friend SupVar<S> && operator/
    ( double const&, SupVar<S> && );

  template <typename S> friend SupVar<S> operator/
    ( SupVar<S> const&, double const& );
  template <typename S> friend SupVar<S> && operator/
    ( SupVar<S> &&, double const& );

  template <typename S> friend SupVar<S> max
    ( SupVar<S> const&, SupVar<S> const& );
  template <typename S> friend SupVar<S> max
    ( SupVar<S> const&, double const& );
  template <typename S> friend SupVar<S> && max
    ( SupVar<S> &&, double const& );
  template <typename S> friend SupVar<S> min
    ( SupVar<S> const&, SupVar<S> const& );
  template <typename S> friend SupVar<S> min
    ( SupVar<S> const&, double const& );
  template <typename S> friend SupVar<S> && min
    ( SupVar<S> &&, double const& );

  template <typename S> friend SupVar<S> sqrt
    ( SupVar<S> const& );
  template <typename S> friend SupVar<S> && sqrt
    ( SupVar<S> && );
  template <typename S> friend SupVar<S> fabs
    ( SupVar<S> const& );
  template <typename S> friend SupVar<S> && fabs
    ( SupVar<S> && );  
  template <typename S> friend SupVar<S> exp
    ( SupVar<S> const& );
  template <typename S> friend SupVar<S> && exp
    ( SupVar<S> && );
  template <typename S> friend SupVar<S> log
    ( SupVar<S> const& );
  template <typename S> friend SupVar<S> && log
    ( SupVar<S> && );
  template <typename S> friend SupVar<S> xlog
    ( SupVar<S> const& );
  template <typename S> friend SupVar<S> && xlog
    ( SupVar<S> && );
  template <typename S> friend SupVar<S> sin
    ( SupVar<S> const& );
  template <typename S> friend SupVar<S> && sin
    ( SupVar<S> && );
  template <typename S> friend SupVar<S> cos
    ( SupVar<S> const& );
  template <typename S> friend SupVar<S> && cos
    ( SupVar<S> && );
  template <typename S> friend SupVar<S> tan
    ( SupVar<S> const& );
  template <typename S> friend SupVar<S> && tan
    ( SupVar<S> && );
  template <typename S> friend SupVar<S> tanh
    ( SupVar<S> const& );
  template <typename S> friend SupVar<S> && tanh
    ( SupVar<S> && );
  template <typename S> friend SupVar<S> pow
    ( SupVar<S> const&, int const& );
  template <typename S> friend SupVar<S> && pow
    ( SupVar<S> &&, int const& );
  template <typename S> friend SupVar<S> pow
    ( SupVar<S> const&, double const& );
  template <typename S> friend SupVar<S> && pow
    ( SupVar<S> &&, double const& );
  template <typename S> friend SupVar<S> intersect
    ( SupVar<S> const& , double const& );
  template <typename S> friend SupVar<S> && intersect
    ( SupVar<S> && , double const& );  
  template <typename S> friend SupVar<S> cheb
    ( SupVar<S> const&, unsigned int const& n );
  template <typename S> friend SupVar<S> && cheb
    ( SupVar<S> &&, unsigned int const& n );

 private:
 
  //! @brief Pointer to superposition model
  SupModel<Summand>* _mod;
  //! @brief Constant estimator (in which case both _uest and _oest are empty)
  double _cst;
  //! @brief Set of dependencies
  std::set<unsigned int> _sdep;

  //! @brief Vector of underestimator summands
  std::vector<Summand> _uest;
  //! @brief Vector of overestimator summands
  std::vector<Summand> _oest;
  
  //! @brief lower bound on primal underestimator
  mutable std::pair<double,bool> _lbd;
  //! @brief upper bound on primal overestimator
  mutable std::pair<double,bool> _ubd;

  //! @brief Vector of shadow underestimator summands
  std::vector<Summand> _uest2;
  //! @brief Vector of shadow overestimator summands
  std::vector<Summand> _oest2;
  
  //! @brief lower bound on shadow underestimator
  mutable std::pair<double,bool> _lbd2;
  //! @brief upper bound on shadow overestimator
  mutable std::pair<double,bool> _ubd2;

 public:

  SupVar<Summand>& operator+=
    ( SupVar<Summand> const& );
  //SupVar<Summand>& operator+=
  //  ( SupVar<Summand> && );
  SupVar<Summand>& operator+=
    ( double const& );

  SupVar<Summand>& operator-=
    ( SupVar<Summand> const& );
  SupVar<Summand>& operator-=
    ( double const& );
 
  SupVar<Summand>& operator*=
    ( SupVar<Summand> const& );
  SupVar<Summand>& operator*=
    ( SupVar<Summand> && );
  SupVar<Summand>& operator*=
    ( double const& );

  SupVar<Summand>& operator/=
    ( SupVar<Summand> const& );
  SupVar<Summand>& operator/=
    ( SupVar<Summand> && );
  SupVar<Summand>& operator/=
    ( double const& );

  SupVar<Summand>& operator=
  ( double const& cst )
  {
    _mod = nullptr;
    _cst = cst;
    _sdep.clear();

    _uest.clear();
    _oest.clear();

    _lbd = std::make_pair( 0., false );
    _ubd = std::make_pair( 0., false );

    _uest2.clear();
    _oest2.clear();

    _lbd2 = std::make_pair( 0., false );
    _ubd2 = std::make_pair( 0., false );

    return *this;
  }

  SupVar<Summand>& operator=
  ( SupVar<Summand> const& var )
  {
#ifdef MC__SUPMODEL_TRACE
    std::cerr << "-- SupVar<Summand>& operator= ( SupVar<Summand> const& )\n";
#endif
    if( this == &var )
      return *this;

    _mod = var._mod;
    if( !_mod || var._sdep.empty() )
      _cst = var._cst;
    _sdep = var._sdep;

    _uest = var._uest;
    _oest = var._oest;

    _lbd  = var._lbd;
    _ubd  = var._ubd;

    _uest2 = var._uest2;
    _oest2 = var._oest2;

    _lbd2  = var._lbd2;
    _ubd2  = var._ubd2;

    return *this;
  } 

  SupVar<Summand>& operator=
  ( SupVar<Summand> && var )
  {
#ifdef MC__SUPMODEL_TRACE
    std::cerr << "-- SupVar<Summand>& operator= ( SupVar<Summand> && )\n";
#endif
    if( this == &var )
      return *this;

    _mod = var._mod;
    if( !_mod || var._sdep.empty() )
      _cst = std::move( var._cst );
    _sdep = std::move( var._sdep );

    _uest = std::move( var._uest );
    _oest = std::move( var._oest );

    _lbd  = std::move( var._lbd );
    _ubd  = std::move( var._ubd );

    _uest2 = std::move( var._uest2 );
    _oest2 = std::move( var._oest2 );

    _lbd2  = std::move( var._lbd2 );
    _ubd2  = std::move( var._ubd2 );

    return *this;
  }

  SupVar
  ( SupModel<Summand>& mod )
  : _mod( &mod ),
    _cst( 0. ),
    _uest( mod._nvar ),
    _oest( mod._nvar ),
    _lbd( 0., false ),
    _ubd( 0., false ),
    // do not presize shadow estimators
    _lbd2( 0., false ),
    _ubd2( 0., false )
  {}

  template <typename Range>
  SupVar
  ( SupModel<Summand>& mod, unsigned int ndx, Range const& bnd, size_t const ndiv=1 )
  : _mod( &mod ),
    _cst( 0. ),
    _lbd( Op<Range>::l(bnd), true ),
    _ubd( Op<Range>::u(bnd), true ),
    // do not presize shadow estimators
    _lbd2( 0., false ),
    _ubd2( 0., false )
  {
    if( ndx >= mod._nvar )
      throw typename SupModel<Summand>::Exceptions( SupModel<Summand>::Exceptions::INDEX );

    _mod->_defvar[ndx] = true;
    _mod->_lbdvar[ndx] = _lbd.first;
    _mod->_ubdvar[ndx] = _ubd.first;

    _sdep = { ndx };
    _uest.resize( mod._nvar );
    _oest.resize( mod._nvar );
    _uest[ndx].set( _lbd.first, _ubd.first, ndiv );
    _oest[ndx].set( _lbd.first, _ubd.first, ndiv );
  }

  SupVar
  ( SupModel<Summand>& mod, unsigned int ndx, Summand const& est )
  : _mod( &mod ),
    _cst( 0. ),
    // do not presize shadow estimators
    _lbd2( 0., false ),
    _ubd2( 0., false )
  {
    if( ndx >= mod._nvar )
      throw typename SupModel<Summand>::Exceptions( SupModel<Summand>::Exceptions::INDEX );

    _mod->_defvar[ndx] = true;
    _mod->_lbdvar[ndx] = est.l();
    _mod->_ubdvar[ndx] = est.u();

    _sdep = { ndx };
    _uest.resize( mod._nvar );
    _oest.resize( mod._nvar );
    _uest[ndx] = _oest[ndx] = est;

    _lbd = std::make_pair( _mod->_lbdvar[ndx], true );
    _ubd = std::make_pair( _mod->_ubdvar[ndx], true );
  }

  SupVar
  ( double const& cst=0. )
  : _mod( nullptr ),
    _cst( cst ),
    _lbd( 0., false ),
    _ubd( 0., false ),
    _lbd2( 0., false ),
    _ubd2( 0., false )
  {}

  SupVar
  ( SupVar<Summand> const& var )
  : _mod( var._mod )
  {
#ifdef MC__SUPMODEL_TRACE
    std::cerr << "-- SupVar( SupVar<Summand> const& )\n";
#endif
    if( this == &var ) return;

    if( !_mod || var._sdep.empty() )
      _cst = var._cst;
    _sdep = var._sdep;

    _uest = var._uest;
    _oest = var._oest;

    _lbd  = var._lbd;
    _ubd  = var._ubd;

    _uest2 = var._uest2;
    _oest2 = var._oest2;

    _lbd2  = var._lbd2;
    _ubd2  = var._ubd2;
  }

  SupVar
  ( SupVar<Summand> && var )
  : _mod( var._mod )
  {
#ifdef MC__SUPMODEL_TRACE
    std::cerr << "-- SupVar( SupVar<Summand> && var )\n";
#endif
    if( this == &var ) return;

    _cst  = std::move( var._cst );
    _sdep = std::move( var._sdep );

    _uest = std::move( var._uest );
    _oest = std::move( var._oest );

    _lbd  = std::move( var._lbd );
    _ubd  = std::move( var._ubd );

    _uest2 = std::move( var._uest2 );
    _oest2 = std::move( var._oest2 );

    _lbd2  = std::move( var._lbd2 );
    _ubd2  = std::move( var._ubd2 );
  }

  ~SupVar
  () 
  {}

  SupVar<Summand>& set
  ( SupModel<Summand>& mod )
  {
    if( _mod == &mod )
      return *this;

    _mod = &mod;
    _sdep.clear();
    _cst = 0.;
    
    _uest.resize( mod._nvar );
    _oest.resize( mod._nvar );

    _lbd = std::make_pair( 0., false );
    _ubd = std::make_pair( 0., false );

    _uest2.clear();
    _oest2.clear();

    _lbd2 = std::make_pair( 0., false );
    _ubd2 = std::make_pair( 0., false );

    return *this;
  }

  template <typename Range>
  SupVar<Summand>& set
  ( SupModel<Summand>& mod, unsigned int ndx, Range const& bnd, size_t const ndiv=1 )
  {
    if( ndx >= mod._nvar )
      throw typename SupModel<Summand>::Exceptions( SupModel<Summand>::Exceptions::INDEX );

    _mod = &mod;
    _cst = 0.;

    _uest.resize( mod._nvar );
    _oest.resize( mod._nvar );

    _lbd = std::make_pair( Op<Range>::l(bnd), true ),
    _ubd = std::make_pair( Op<Range>::u(bnd), true ),

    _sdep = { ndx };
    _mod->_defvar[ndx] = true;
    _mod->_lbdvar[ndx] = _lbd.first;
    _mod->_ubdvar[ndx] = _ubd.first;

    _uest[ndx].set( _lbd.first, _ubd.first, ndiv );
    _oest[ndx].set( _lbd.first, _ubd.first, ndiv );

    _uest2.clear();
    _oest2.clear();

    _lbd2 = std::make_pair( 0., false );
    _ubd2 = std::make_pair( 0., false );
 
    return *this;
  }

  SupVar<Summand>& set
  ( SupModel<Summand>& mod, unsigned int ndx, Summand const& est )
  {
    if( ndx >= mod._nvar )
      throw typename SupModel<Summand>::Exceptions( SupModel<Summand>::Exceptions::INDEX );

    _mod = &mod;
    _cst = 0.;

    _sdep = { ndx };
    _mod->_defvar[ndx] = true;
    _mod->_lbdvar[ndx] = est.l();
    _mod->_ubdvar[ndx] = est.u();

    _uest.resize( mod._nvar );
    _oest.resize( mod._nvar );
    _uest[ndx] = _oest[ndx] = est;

    _lbd = std::make_pair( _mod->_lbdvar[ndx], true );
    _ubd = std::make_pair( _mod->_ubdvar[ndx], true );
 
    _uest2.clear();
    _oest2.clear();

    _lbd2 = std::make_pair( 0., false );
    _ubd2 = std::make_pair( 0., false );
   
    return *this;
  }

  SupVar<Summand>& set
  ( double const& cst=0. )
  {
    _mod = nullptr;
    _cst = cst;
    _sdep.clear();

    _uest.clear();
    _oest.clear();
 
    _lbd = std::make_pair( 0., false );
    _ubd = std::make_pair( 0., false );

    _uest2.clear();
    _oest2.clear();

    _lbd2 = std::make_pair( 0., false );
    _ubd2 = std::make_pair( 0., false );
 
    return *this;
  }

  std::vector<Summand>& uest
  ( unsigned int opt=0 )
  { return opt? _uest2: _uest; }

  std::vector<Summand>& oest
  ( unsigned int opt=0 )
  { return opt? _oest2: _oest; }

  std::vector<Summand> const& uest
  ( unsigned int opt=0 )
  const
  { return opt? _uest2: _uest; }

  std::vector<Summand> const& oest
  ( unsigned int opt=0 )
  const
  { return opt? _oest2: _oest; }

  //! @brief Retrieve constant field - relevant only when dependent set is empty
  double const& cst
  ()
  const
  { return _cst; }

  //! @brief Retrieve number of dependents
  unsigned int ndep
  ()
  const
  { return _sdep.size(); }

  //! @brief Retrieve set of dependents
  std::set<unsigned int> const& sdep
  ()
  const
  { return _sdep; }
 
   //! @brief Lower bound on superposition underestimator. 0: primal estimator only; 1: shadow estimator only; 2: best estimator (default); >=3: best estimator recomputed
  double const& l
  ( unsigned int opt=2 )
  const
  {
    if( !_mod || _sdep.empty() )
      return _cst;

    switch( opt ){
      default:
      // lower bound on primal underestimator
      case 0:
        if( !_lbd.second || opt > 2 ){
          _lbd.second = true;
          _lbd.first  = 0.;
          for( auto const& i : _sdep )
            _lbd.first += _uest[i].l();
        }
        if( opt == 0 || _uest2.empty() ) return _lbd.first;

      // lower bound on shadow underestimator
      case 1:
        if( _uest2.empty() ) return -DBL_MAX; // -infinity
        if( !_lbd2.second || opt > 2 ){
          _lbd2.second = true;
          _lbd2.first  = 0.;
          for( auto const& i : _sdep )
            _lbd2.first += _uest2[i].l();
        }
        if( opt == 1 ) return _lbd2.first;

        // best lower bound between primal and shadow underestimator
        return _lbd.first > _lbd2.first? _lbd.first: _lbd2.first; 
    }
  }

   //! @brief Upper bound on superposition overestimator. 0: primal estimator only; 1: shadow estimator only; 2: best estimator (default); >=3: best estimator recomputed
  double const& u
  ( unsigned int opt=2 )
  const
  {
    if( !_mod || _sdep.empty() )
      return _cst;

    switch( opt ){
      default:
      // lower bound on primal underestimator
      case 0:
        if( !_ubd.second || opt > 2 ){
          _ubd.second = true;
          _ubd.first  = 0.;
          for( auto const& i : _sdep )
            _ubd.first += _oest[i].u();
        }
        if( opt == 0 || _oest2.empty() ) return _ubd.first;

      // lower bound on shadow underestimator
      case 1:
        if( _oest2.empty() ) return DBL_MAX; // +infinity
        if( !_ubd2.second || opt > 2 ){
          _ubd2.second = true;
          _ubd2.first  = 0.;
          for( auto const& i : _sdep )
            _ubd2.first += _oest2[i].u();
        }
        if( opt == 1 ) return _ubd2.first;

        // best lower bound between primal and shadow underestimator
        return _ubd.first < _ubd2.first? _ubd.first: _ubd2.first; 
    }
  }
  
   //! @brief Evaluate superposition underestimator at point. 0: primal estimator only; 1: shadow estimator only; >=2: best estimator (default)
  double uval
  ( std::map<unsigned int,double> const& x, unsigned int opt=2 )
  const
  {
    if( !_mod || _sdep.empty() )
      return _cst;

    double y(0.), y2(0.);
    switch( opt ){
      default:
      // value of primal underestimator
      case 0:
        for( auto it=_sdep.cbegin(); it!=_sdep.cend(); ++it )
          if( it == _sdep.cbegin() ) y  = _uest.at(*it).l( x.at(*it) );
          else                       y += _uest.at(*it).l( x.at(*it) );
        if( opt == 0 || _uest2.empty() ) return y;

      // value of shadow underestimator
      case 1:
        if( _uest2.empty() ) return -DBL_MAX; // -infinity
        for( auto it=_sdep.cbegin(); it!=_sdep.cend(); ++it )
          if( it == _sdep.cbegin() ) y2  = _uest2.at(*it).l( x.at(*it) );
          else                       y2 += _uest2.at(*it).l( x.at(*it) );
        if( opt == 1 ) return y2;

        // best (max) value between primal and shadow underestimator
        return y > y2? y: y2;
    }
  }
  
   //! @brief Evaluate superposition overestimator at point. 0: primal estimator only; 1: shadow estimator only; >=2: best estimator (default)
  double oval
  ( std::map<unsigned int,double> const& x, unsigned int opt=2 )
  const
  {
    if( !_mod || _sdep.empty() )
      return _cst;

    double y(0.), y2(0.);
    switch( opt ){
      default:
      // value of primal underestimator
      case 0:
        for( auto it=_sdep.cbegin(); it!=_sdep.cend(); ++it )
          y += _oest.at(*it).u( x.at(*it) );
        if( opt == 0 || _oest2.empty() ) return y;

      // value of shadow underestimator
      case 1:
        if( _oest2.empty() ) return DBL_MAX; // +infinity
        for( auto it=_sdep.cbegin(); it!=_sdep.cend(); ++it )
          y2 += _oest2.at(*it).u( x.at(*it) );
        if( opt == 1 ) return y2;

        // best (min) value between primal and shadow underestimator
        return y < y2? y: y2;
    }
  }

 private:
/* 
  static void _add
    ( std::set<unsigned> const& sdep,  std::vector<Summand>& est, 
      std::set<unsigned> const& sdep2, std::vector<Summand> const& est2 )
    {
      for( auto& i : sdep2 )
        if( !sdep.count(i) ) est[i]  = est2[i];
        else                 est[i] += est2[i];
    }
*/
  static void _lb
  ( double& lb, std::set<unsigned> const& sdep,  std::vector<Summand> const& est )
  {
    lb = 0.;
    for( auto const& i : sdep )
      lb += est[i].l();
  }

  static void _ub
  ( double& ub, std::set<unsigned> const& sdep,  std::vector<Summand> const& est )
  {
    ub = 0.;
    for( auto const& i : sdep )
      ub += est[i].u();
  }
};

template <typename Summand> 
std::ostream& operator<<
( std::ostream& out, SupVar<Summand> const& var )
{
  if( !var._mod || var._sdep.empty() ){
    out << std::right << std::scientific << std::setprecision(5)
        << "[ "  << std::setw(12) << var._cst
        << std::endl;
    return out;
  }
  
  for( auto const& i : var._sdep )
    out << std::right << std::setw(5) << "U[" << i << "]: " << var._uest[i] << " (" << var._uest[i].size() << ")" << std::endl;
  for( auto const& i : var._sdep )
    out << std::right << std::setw(5) << "O[" << i << "]: " << var._oest[i] << " (" << var._oest[i].size() << ")" << std::endl;

  if( var._mod->options.DISPLAY_SHADOW && !var._uest2.empty() )
    for( auto const& i : var._sdep )
      out << std::right << std::setw(5) << "U2[" << i << "]: " << var._uest2[i] << " (" << var._uest2[i].size() << ")" << std::endl;
  if( var._mod->options.DISPLAY_SHADOW && !var._oest2.empty() )
    for( auto const& i : var._sdep )
      out << std::right << std::setw(5) << "O2[" << i << "]: " << var._oest2[i] << " (" << var._oest2[i].size() << ")" << std::endl;

  out << std::right << "   RNG:  ";
  out << std::right << std::scientific << std::setprecision(var._mod->options.DISPLAY_DIGITS);
  out << "[ "  << std::setw(var._mod->options.DISPLAY_DIGITS+7) << var.l()
      << " : " << std::setw(var._mod->options.DISPLAY_DIGITS+7) << var.u() << " ]"
      << std::endl;

  return out;
}

template <typename Summand>
inline
SupVar<Summand> operator+
( SupVar<Summand> const& var )
{
  return var;
}

template <typename Summand>
inline
SupVar<Summand> && operator+
( SupVar<Summand> && var )
{
  return( std::move(var) );  
}

template <typename Summand>
inline
SupVar<Summand>& SupVar<Summand>::operator+=
( double const& cst )
{
  if( cst == 0. ){
    return *this;
  }
  
  if( !_mod || _sdep.empty() ){
    _cst += cst;
    return *this;
  }

  double const& offset = cst / (double)_sdep.size();
  for( auto& i : _sdep ){
    _uest[i] += offset;
    _oest[i] += offset;
    if( !_uest2.empty() ) _uest2[i] += offset;
    if( !_oest2.empty() ) _oest2[i] += offset;
  }
  if( _lbd.second )  _lbd.first  += cst;
  if( _ubd.second )  _ubd.first  += cst;
  if( _lbd2.second ) _lbd2.first += cst;
  if( _ubd2.second ) _ubd2.first += cst;

  return *this;
}

template <typename Summand>
inline
SupVar<Summand>& SupVar<Summand>::operator+=
( SupVar<Summand> const& var )
{
  if( !_mod && !var._mod ){
    _cst += var._cst;
    return *this;
  }

  if( !var._mod || var._sdep.empty() ){
    *this += var._cst;
    return *this;
  }

  if( !_mod || _sdep.empty() ){
    double copy_cst = _cst;
    *this = var;
    *this += copy_cst;
    return *this;
  }

  if( _mod != var._mod )
    throw typename SupModel<Summand>::Exceptions( SupModel<Summand>::Exceptions::MODEL );

#ifdef MC__SUPMODEL_TRACE
    std::cerr << "-- SupVar<Summand>& operator+=( SupVar<Summand> const& )\n";
#endif

  // Updated dependency set
  auto sdep = _sdep;
  for( auto const& i : var._sdep ){
    if( sdep.count(i) ) continue;
    sdep.insert(i);
  }

  // No shadow underestimator
  if( _uest2.empty() && var._uest2.empty() ){
    _mod->_add( _sdep, _uest, var._sdep, var._uest, 1 );
    _lbd.second = false;
  }

  // Other operand has shadow underestimator
  else if( _uest2.empty() ){
    _uest2 = _uest; // Duplicate current primal underestimator
    _mod->_add( _sdep, _uest,  var._sdep, var._uest,  1 );
    _mod->_add( _sdep, _uest2, var._sdep, var._uest2, 1 );
    _lb( _lbd.first,  sdep, _uest  );
    _lb( _lbd2.first, sdep, _uest2 );
    _lbd.second = _lbd2.second = true;
    if( _lbd.first < _lbd2.first ){
      std::swap( _uest, _uest2 );
      std::swap( _lbd.first, _lbd2.first );
    }
  }

  // Current operand has shadow underestimator
  else if( var._uest2.empty() ){
    _mod->_add( _sdep, _uest,  var._sdep, var._uest, 1 );
    _mod->_add( _sdep, _uest2, var._sdep, var._uest, 1 );
    _lb( _lbd.first,  sdep, _uest  );
    _lb( _lbd2.first, sdep, _uest2 );
    _lbd.second = _lbd2.second = true;
    if( _lbd.first < _lbd2.first ){
      std::swap( _uest, _uest2 );
      std::swap( _lbd.first, _lbd2.first );
    }
  }

  // Both operands have shadow underestimators
  else{
    _mod->_est3 = _uest;  // Duplicate current primal underestimator
    _mod->_est4 = _uest2; // Duplicate current shadow underestimator
    _mod->_add( _sdep, _uest,  var._sdep, var._uest,  1 );
    _mod->_add( _sdep, _uest2, var._sdep, var._uest2, 1 );
    _mod->_add( _sdep, _mod->_est3, var._sdep, var._uest2, 1 );
    _mod->_add( _sdep, _mod->_est4, var._sdep, var._uest,  1 );

    _lb( _lbd.first,  sdep, _uest  );
    _lb( _lbd2.first, sdep, _uest2 );
    _lbd.second = _lbd2.second = true;
    if( _lbd.first < _lbd2.first ){
      std::swap( _uest, _uest2 );
      std::swap( _lbd.first, _lbd2.first );
    }

    // Check if better primal underestimator
    _lb( _mod->_lbd3, sdep, _mod->_est3 );
    _lb( _mod->_lbd4, sdep, _mod->_est4 );
    if( _mod->_lbd4 < _mod->_lbd3 && _lbd.first < _mod->_lbd3 ){
      std::swap( _uest, _mod->_est3 );
      std::swap( _lbd.first, _mod->_lbd3 );
    }
    else if( _lbd.first < _mod->_lbd4 ){
      std::swap( _uest, _mod->_est4 );
      std::swap( _lbd.first, _mod->_lbd4 );
    }

    // Check if better shadow underestimator
    _ub( _mod->_ubd3, sdep, _mod->_est3 );
    _ub( _mod->_ubd4, sdep, _mod->_est4 );
    _ub( _mod->_ubd2, sdep, _uest2 );
    if( _mod->_ubd4 < _mod->_ubd3 && _mod->_ubd2 < _mod->_ubd3 ){
      std::swap( _uest2, _mod->_est3 );
      std::swap( _lbd2.first, _mod->_lbd3 );
    }
    else if( _mod->_ubd2 < _mod->_ubd4 ){
      std::swap( _uest2, _mod->_est4 );
      std::swap( _lbd2.first, _mod->_lbd4 );
    }
  }

  // No shadow overestimator
  if( _oest2.empty() && var._oest2.empty() ){
    _mod->_add( _sdep, _oest, var._sdep, var._oest, 0 );
    _ubd.second = false;
  }

  // Other operand has shadow overestimator
  else if( _oest2.empty() ){
    _oest2 = _oest; // Duplicate current primal overestimator
    _mod->_add( _sdep, _oest,  var._sdep, var._oest,  0 );
    _mod->_add( _sdep, _oest2, var._sdep, var._oest2, 0 );
    _ub( _ubd.first,  sdep, _oest  );
    _ub( _ubd2.first, sdep, _oest2 );
    _ubd.second = _ubd2.second = true;
    if( _ubd.first > _ubd2.first ){
      std::swap( _oest, _oest2 );
      std::swap( _ubd.first, _ubd2.first );
    }
  }

  // Current operand has shadow overestimator
  else if( var._oest2.empty() ){
    _mod->_add( _sdep, _oest,  var._sdep, var._oest, 0 );
    _mod->_add( _sdep, _oest2, var._sdep, var._oest, 0 );
    _ub( _ubd.first,  sdep, _oest  );
    _ub( _ubd2.first, sdep, _oest2 );
    _ubd.second = _ubd2.second = true;
    if( _ubd.first > _ubd2.first ){
      std::swap( _oest, _oest2 );
      std::swap( _ubd.first, _ubd2.first );
    }
  }

  // Both operands have shadow overestimators
  else{
    _mod->_est3 = _oest;  // Duplicate current primal overestimator
    _mod->_est4 = _oest2; // Duplicate current shadow overestimator
    _mod->_add( _sdep, _oest,  var._sdep, var._oest,  0 );
    _mod->_add( _sdep, _oest2, var._sdep, var._oest2, 0 );
    _mod->_add( _sdep, _mod->_est3, var._sdep, var._oest2, 0 );
    _mod->_add( _sdep, _mod->_est4, var._sdep, var._oest,  0 );

    _ub( _ubd.first,  sdep, _oest  );
    _ub( _ubd2.first, sdep, _oest2 );
    _ubd.second = _ubd2.second = true;
    if( _ubd.first > _ubd2.first ){
      std::swap( _oest, _oest2 );
      std::swap( _ubd.first, _ubd2.first );
    }

    // Check if better primal overestimator
    _ub( _mod->_ubd3, sdep, _mod->_est3 );
    _ub( _mod->_ubd4, sdep, _mod->_est4 );
    if( _mod->_ubd4 > _mod->_ubd3 && _ubd.first > _mod->_ubd3 ){
      std::swap( _oest, _mod->_est3 );
      std::swap( _ubd.first, _mod->_ubd3 );
    }
    else if( _ubd.first > _mod->_ubd4 ){
      std::swap( _oest, _mod->_est4 );
      std::swap( _ubd.first, _mod->_ubd4 );
    }

    // Check if better shadow overestimator
    _lb( _mod->_lbd3, sdep, _mod->_est3 );
    _lb( _mod->_lbd4, sdep, _mod->_est4 );
    _lb( _mod->_lbd2, sdep, _oest2 );
    if( _mod->_lbd4 > _mod->_lbd3 && _mod->_lbd2 > _mod->_lbd3 ){
      std::swap( _oest2, _mod->_est3 );
      std::swap( _ubd2.first, _mod->_ubd3 );
    }
    else if( _mod->_lbd2 > _mod->_lbd4 ){
      std::swap( _oest2, _mod->_est4 );
      std::swap( _ubd2.first, _mod->_ubd4 );
    }
  }

  // Update dependencies in current operand
  std::swap( _sdep, sdep );

  return *this;
}

template <typename Summand>
inline
SupVar<Summand> operator+
( SupVar<Summand> const& var1, SupVar<Summand> const& var2 )
{
  if( !var1._mod && !var2._mod ) 
    return( var1._cst + var2._cst );
  if( var1._mod ){
    SupVar<Summand> var3( var1 );
    var3 += var2;
    return var3;
  }
  SupVar<Summand> var3( var2 );
    var3 += var1;
    return var3;
}

template <typename Summand>
inline
SupVar<Summand> && operator+
( SupVar<Summand> const& var1, SupVar<Summand> && var2 )
{
  if( !var1._mod && !var2._mod ){
    var2._cst += var1._cst; 
    return std::move( var2 );
  }

  var2 += var1;
  return std::move( var2 );
}

template <typename Summand>
inline
SupVar<Summand> && operator+
( SupVar<Summand> && var1, SupVar<Summand> const& var2 )
{
  if( !var1._mod && !var2._mod ){
    var1._cst += var2._cst; 
    return std::move( var1 );
  }
   
  var1 += var2;
  return std::move( var1 );
}

template <typename Summand>
inline
SupVar<Summand> && operator+
( SupVar<Summand> && var1, SupVar<Summand> && var2 )
{
  if( !var1._mod && !var2._mod ){
    var1._cst += var2._cst; 
    return std::move( var1 );
  }

  if( var1._mod || !var2._mod ){
    var1 += var2;
    return std::move( var1 );
  }

  var2 += var1;
  return std::move( var2 );
}

template <typename Summand>
inline
SupVar<Summand> operator+
( SupVar<Summand> const& var1, double const& cst2 )
{
  if( !var1._mod )
    return( var1._cst + cst2 );

  SupVar<Summand> var3( var1 );
  var3 += cst2;
  return var3;
}

template <typename Summand>
inline
SupVar<Summand> && operator+
( SupVar<Summand> && var1, double const& cst2 )
{
  var1 += cst2;
  return std::move( var1 );
}

template <typename Summand>
inline
SupVar<Summand> operator+
( double const& cst1, SupVar<Summand> const& var2 )
{
  if( !var2._mod ) 
    return( cst1 + var2._cst );

  SupVar<Summand> var3( var2 );
  var3 += cst1;
  return var3;
}

template <typename Summand>
inline
SupVar<Summand> && operator+
( double const& cst1, SupVar<Summand> && var2 )
{
  var2 += cst1;
  return std::move( var2 );
}

template <typename Summand> 
inline
SupVar<Summand> operator-
( SupVar<Summand> const& var )
{
  SupVar<Summand> var2( var );
  return operator-( std::move( var2 ) );
}


template <typename Summand> 
inline
SupVar<Summand> && operator-
( SupVar<Summand> && var )
{
  if( !var._mod || var._sdep.empty() ){
    var._cst *= -1;
    return std::move( var );
  }

  std::swap( var._uest, var._oest );
  for( auto& i : var._sdep ){
    var._uest[i] *= -1;
    var._oest[i] *= -1;
  }

  std::swap( var._lbd, var._ubd );
  if( var._lbd.second ) var._lbd.first *= -1;
  if( var._ubd.second ) var._ubd.first *= -1;

  if( !var._uest2.empty() || !var._oest2.empty() ){
    std::swap( var._uest2, var._oest2 );
    for( auto& i : var._sdep ){
      if( !var._uest2.empty() ) var._uest2[i] *= -1;
      if( !var._oest2.empty() ) var._oest2[i] *= -1;
    }

    std::swap( var._lbd2, var._ubd2 );
    if( var._lbd2.second ) var._lbd2.first *= -1;
    if( var._ubd2.second ) var._ubd2.first *= -1;
  }

  return std::move( var );
}

template <typename Summand>
inline
SupVar<Summand>& SupVar<Summand>::operator-=
( double const& cst )
{
  *this += -cst;
  return *this;
}

template <typename Summand>
inline
SupVar<Summand>& SupVar<Summand>::operator-=
( SupVar<Summand> const& var )
{
  SupVar<Summand> var2( var );
  return operator+=( operator-( std::move( var2 ) ) );
}

template <typename Summand>
inline
SupVar<Summand> operator-
( SupVar<Summand> const& var1, SupVar<Summand> const& var2 )
{
  if( !var1._mod && !var2._mod ) 
    return( var1._cst - var2._cst );

  if( var1._mod ){
    SupVar<Summand> var3( var1 );
    var3 -= var2;    
    return var3;
  }

  SupVar<Summand> var3( -var2 );
  var3 += var1;
  return var3;
}

template <typename Summand>
inline
SupVar<Summand> && operator-
( SupVar<Summand> const& var1, SupVar<Summand> && var2 )
{
  if( !var1._mod && !var2._mod ){
    var2._cst = var1._cst - var2._cst;
    return std::move( var2 );
  }
  
  var2 -= var1;
  return operator-( std::move( var2 ) );
}

template <typename Summand>
inline
SupVar<Summand> && operator-
( SupVar<Summand> && var1, SupVar<Summand> const& var2 )
{
  if( !var1._mod && !var2._mod ){
    var1._cst -= var2._cst; 
    return std::move( var1 );
  }
   
  var1 -= var2;
  return std::move( var1 );
}

template <typename Summand>
inline
SupVar<Summand> && operator-
( SupVar<Summand> && var1, SupVar<Summand> && var2 )
{
  if( !var1._mod && !var2._mod ){
    var1._cst -= var2._cst;
    return std::move( var1 );
  }
  
  var1 -= std::move( var2 );
  return std::move( var1 );
}

template <typename Summand>
inline
SupVar<Summand> operator-
( SupVar<Summand> const& var1, double const& cst2 )
{
  if( !var1._mod ) 
    return( var1._cst - cst2 );
    
  SupVar<Summand> var3( var1 );
  var3 += -cst2;
  return var3;
}

template <typename Summand>
inline
SupVar<Summand> && operator-
( SupVar<Summand> && var1, double const& cst2 )
{
  var1 -= cst2;
  return std::move( var1 );
}

template <typename Summand>
inline
SupVar<Summand> operator-
( double const& cst1, SupVar<Summand> const& var2 )
{
  if( !var2._mod ) 
    return( cst1 - var2._cst );

  SupVar<Summand> var3( -var2 );
  var3 += cst1;
  return var3;
}

template <typename Summand>
inline
SupVar<Summand> && operator-
( double const& cst1, SupVar<Summand> && var2 )
{
  var2 -= cst1;
  return operator-( std::move( var2 ) );
}

template <typename Summand>
inline
SupVar<Summand>& SupVar<Summand>::operator*=
( double const& cst )
{
  if( cst == 0. ){
    *this = 0.;
    return *this;
  }

  if( cst == 1. ){
    return *this;
  }

  if( !_mod || _sdep.empty() ){
    _cst *= cst;
    return *this;
  }

  if( _sdep.empty() )
    throw typename SupModel<Summand>::Exceptions( SupModel<Summand>::Exceptions::INTERN );
    
#ifdef MC__SUPMODEL_TRACE
    std::cerr << "-- SupVar<Summand>& operator*=( double const& )\n";
#endif

  if( cst < 0. ){
    std::swap( _uest, _oest ); 
    std::swap( _lbd, _ubd );
    if( !_uest2.empty() || !_oest2.empty() ){
      std::swap( _uest2, _oest2 ); 
      std::swap( _lbd2, _ubd2 ); 
    }
  }
  
  for( auto& i : _sdep ){
    _uest[i] *= cst;
    _oest[i] *= cst;
    if( !_uest2.empty() ) _uest2[i] *= cst;
    if( !_oest2.empty() ) _oest2[i] *= cst;
  }
  if( _lbd.second ) _lbd.first *= cst;
  if( _ubd.second ) _ubd.first *= cst;
  if( _lbd2.second ) _lbd2.first *= cst;
  if( _ubd2.second ) _ubd2.first *= cst;

  return *this;
}

template <typename Summand>
inline
SupVar<Summand>& SupVar<Summand>::operator*=
( SupVar<Summand> const& var )
{
  if( !_mod && !var._mod ){
    _cst *= var._cst;
    return *this;
  }

  if( !var._mod || var._sdep.empty() ){
    *this *= var._cst;
    return *this;
  }

  if( !_mod || _sdep.empty() ){
    double copy_cst = _cst;
    *this = var;
    *this *= copy_cst;
    return *this;
  }

  if( &var == this ){
    pow( std::move(*this), 2 );
    return *this;
  }
  
  if( _mod != var._mod )
    throw typename SupModel<Summand>::Exceptions( SupModel<Summand>::Exceptions::MODEL );

#ifdef MC__SUPMODEL_TRACE
    std::cerr << "-- SupVar<Summand>& operator*=( SupVar<Summand> const& )\n";
#endif

  double const u0 = u(), l0 = l(), u1 = var.u(), l1 = var.l();
#ifdef MC__SUPMODEL_DEBUG_PROD
  std::cout << "variable bounds: " << "[ " << l0 << " : " << u0 << " ] x [ " << l1 << " : " << u1 << " ]" << std::endl;
#endif

  switch( _mod->options.PROD_METH ){
    default:
    case SupModel<Summand>::Options::NONE:
    {
      SupVar<Summand> var1( var );
      var1 -= *this;
      pow( std::move( var1 ), 2 );
      *this += var;
      pow( std::move( *this ), 2 );
      (*this -= var1) *= 0.25;
      break;
    }

    case SupModel<Summand>::Options::PARTIAL:
    {
      double const wid0 = u()-l(), wid1 = var.u()-var.l();
      SupVar<Summand> var1( var );
      var1 /= wid1;
      SupVar<Summand> var2( var1 );
      *this /= wid0;
      var1 -= *this;
      pow( std::move( var1 ), 2 );
      *this += var2;
      pow( std::move( *this ), 2 );
      (*this -= var1) *= 0.25*wid0*wid1;
      break;
    }

    case SupModel<Summand>::Options::FULL:
    {
      double const wid0 = u()-l(),       wid1 = var.u()-var.l(),
                   mid0 = 0.5*(u()+l()), mid1 = 0.5*(var.u()+var.l());
      SupVar<Summand> sub( mid0 * var + mid1 * *this - mid0*mid1 );
      SupVar<Summand> var1( var );
      (var1 -= mid1) /= wid1;
      SupVar<Summand> var2( var1 );
      (*this -= mid0) /= wid0;
      var1 -= *this;
      pow( std::move( var1 ), 2 );
      *this += var2;
      pow( std::move( *this ), 2 );
      ((*this -= var1) *= 0.25*wid0*wid1) += sub;
      break;
    }

    case SupModel<Summand>::Options::LOG:
    {
      double const wid0 = u()-l(), wid1 = var.u()-var.l(),
                   lbd0 = l(),     lbd1 = var.l();
      static double const a = 0.5, b = 1.0;
      SupVar<Summand> sub( (lbd0-a/b*wid0) * var + (lbd1-a/b*wid1) * *this
                         + (a/b*(wid0*lbd1+wid1*lbd0-a/b*wid0*wid1)-lbd0*lbd1) );
      SupVar<Summand> var1( var );
      ((var1  -= lbd1) *= b/wid1) += a;
      ((*this -= lbd0) *= b/wid0) += a;
      log( std::move( var1 ) );
      log( std::move( *this ) );
      exp( std::move( *this += var1 ) );
      (*this *= (wid0/b*wid1/b)) += sub;
      break;
    }
  }
 
  _lbd.second = false;
  _ubd.second = false;
  if( !_mod->options.PROD_CUT )
    return *this;

  // Bounds refinement
  auto const& fmin2 = [=]( const double& x1, const double& x2 )
                      { return x1<x2? x1: x2; };
  auto const& fmin4 = [=]( const double& x1, const double& x2, const double& x3, const double& x4 )
                      { return fmin2(fmin2(x1,x2),fmin2(x3,x4)); };
  auto const& fmax2 = [=]( const double& x1, const double& x2 )
                      { return x1>x2? x1: x2; };
  auto const& fmax4 = [=]( const double& x1, const double& x2, const double& x3, const double& x4 )
                      { return fmax2(fmax2(x1,x2),fmax2(x3,x4)); };

  double const& lnat = fmin4( l0*l1, l0*u1, u0*l1, u0*u1 );
#ifdef MC__SUPMODEL_DEBUG_PROD
  std::cout << "Natural lower bound: " << "[ " << l0 << " : " << u0 << " ] x [ " << l1 << " : " << u1 << " ] -> "
            << lnat << std::endl;
#endif
  if( l() < lnat )
    _mod->uref( *this, lnat );

  double const& unat = fmax4( l0*l1, l0*u1, u0*l1, u0*u1 );
#ifdef MC__SUPMODEL_DEBUG_PROD
  std::cout << "Natural upper bound: " << "[ " << l0 << " : " << u0 << " ] x [ " << l1 << " : " << u1 << " ] -> "
            << unat << std::endl;
#endif
  if( u() > unat )
    _mod->oref( *this, unat );

  return *this;
}

template <typename Summand>
inline
SupVar<Summand>& SupVar<Summand>::operator*=
( SupVar<Summand> && var )
{
  if( !_mod && !var._mod ){
    _cst *= var._cst;
    return *this;
  }

  if( !var._mod || var._sdep.empty() ){
    *this *= var._cst;
    return *this;
  }

  if( !_mod || _sdep.empty() ){
    var *= _cst;
    std::swap( *this, var );
    return *this;
  }

  if( &var == this ){
    pow( std::move(*this), 2 );
    return *this;
  }
  
  if( _mod != var._mod )
    throw typename SupModel<Summand>::Exceptions( SupModel<Summand>::Exceptions::MODEL );

#ifdef MC__SUPMODEL_TRACE
    std::cerr << "-- SupVar<Summand>& operator*=( SupVar<Summand> && )\n";
#endif

  double const u0 = u(), l0 = l(), u1 = var.u(), l1 = var.l();
#ifdef MC__SUPMODEL_DEBUG_PROD
  std::cout << "variable bounds: " << "[ " << l0 << " : " << u0 << " ] x [ " << l1 << " : " << u1 << " ]" << std::endl;
#endif
  
  switch( _mod->options.PROD_METH ){
    default:
    case SupModel<Summand>::Options::NONE:
    {
      SupVar<Summand> var1( var );
      var1 -= *this;
      pow( std::move( var1 ), 2 );
      *this += var;
      pow( std::move( *this ), 2 );
      (*this -= var1) *= 0.25;
      break;
    }

    case SupModel<Summand>::Options::PARTIAL:
    {
      double const wid0 = u()-l(), wid1 = var.u()-var.l();
      var /= wid1;
      SupVar<Summand> var1( var );
      *this /= wid0;
      var -= *this;
      pow( std::move( var ), 2 );
      *this += var1;
      pow( std::move( *this ), 2 );
      (*this -= var) *= 0.25*wid0*wid1;
      break;
    }

    case SupModel<Summand>::Options::FULL:
    {
      double const wid0 = u()-l(),       wid1 = var.u()-var.l(),
                   mid0 = 0.5*(u()+l()), mid1 = 0.5*(var.u()+var.l());
      SupVar<Summand> sub( mid0 * var + mid1 * *this - mid0*mid1 );
      (var -= mid1) /= wid1;
      SupVar<Summand> var1( var );
      (*this -= mid0) /= wid0;
      var -= *this;
      pow( std::move( var ), 2 );
      *this += var1;
      pow( std::move( *this ), 2 );
      ((*this -= var) *= 0.25*wid0*wid1) += sub;
      break;
    }

    case SupModel<Summand>::Options::LOG:
    {
      double const wid0 = u()-l(), wid1 = var.u()-var.l(),
                   lbd0 = l(),     lbd1 = var.l();
      static double const a = 0.5, b = 1.0;
      SupVar<Summand> sub( (lbd0-a/b*wid0) * var + (lbd1-a/b*wid1) * *this
                         + (a/b*(wid0*lbd1+wid1*lbd0-a/b*wid0*wid1)-lbd0*lbd1) );
      ((var   -= lbd1) *= b/wid1) += a;
      ((*this -= lbd0) *= b/wid0) += a;
      log( std::move( var ) );
      log( std::move( *this ) );
      exp( std::move( *this += var ) );
      (*this *= (wid0/b*wid1/b)) += sub;
      break;
    }
  }
 
  _lbd.second = false;
  _ubd.second = false;
  if( !_mod->options.PROD_CUT )
    return *this;

  // Bounds refinement
  auto const& fmin2 = [=]( const double& x1, const double& x2 )
                      { return x1<x2? x1: x2; };
  auto const& fmin4 = [=]( const double& x1, const double& x2, const double& x3, const double& x4 )
                      { return fmin2(fmin2(x1,x2),fmin2(x3,x4)); };
  auto const& fmax2 = [=]( const double& x1, const double& x2 )
                      { return x1>x2? x1: x2; };
  auto const& fmax4 = [=]( const double& x1, const double& x2, const double& x3, const double& x4 )
                      { return fmax2(fmax2(x1,x2),fmax2(x3,x4)); };

  double const& lnat = fmin4( l0*l1, l0*u1, u0*l1, u0*u1 );
#ifdef MC__SUPMODEL_DEBUG_PROD
  std::cout << "Natural lower bound: " << "[ " << l0 << " : " << u0 << " ] x [ " << l1 << " : " << u1 << " ] -> "
            << lnat << std::endl;
#endif
  if( l() < lnat )
    _mod->uref( *this, lnat );

  double const& unat = fmax4( l0*l1, l0*u1, u0*l1, u0*u1 );
#ifdef MC__SUPMODEL_DEBUG_PROD
  std::cout << "Natural upper bound: " << "[ " << l0 << " : " << u0 << " ] x [ " << l1 << " : " << u1 << " ] -> "
            << unat << std::endl;
#endif
  if( u() > unat )
    _mod->oref( *this, unat );

  return *this;
}

template <typename Summand>
inline
SupVar<Summand> operator*
( SupVar<Summand> const& var1, SupVar<Summand> const& var2 )
{
  if( !var1._mod && !var2._mod ) 
    return( var1._cst * var2._cst );
  if( var1._mod ){
    SupVar<Summand> var3( var1 );
    var3 *= var2;
    return var3;
  }
  SupVar<Summand> var3( var2 );
    var3 *= var1;
    return var3;
}

template <typename Summand>
inline
SupVar<Summand> && operator*
( SupVar<Summand> && var1, SupVar<Summand> const& var2 )
{
  var1 *= var2;
  return std::move( var1 );
}

template <typename Summand>
inline
SupVar<Summand> && operator*
( SupVar<Summand> const& var1, SupVar<Summand> && var2 )
{
  var2 *= var1;
  return std::move( var2 );
}

template <typename Summand>
inline
SupVar<Summand> && operator*
( SupVar<Summand> && var1, SupVar<Summand> && var2 )
{
  if( !var1._mod && var2._mod ){
    var2 *= std::move( var1 );
    return std::move( var2 );
  }

  var1 *= std::move( var2 );
  return std::move( var1 );
}

template <typename Summand>
inline
SupVar<Summand> operator*
( SupVar<Summand> const& var1, double const& cst2 )
{
  if( !var1._mod ) 
    return( var1._cst * cst2 );

  SupVar<Summand> var3( var1 );
  var3 *= cst2;
  return var3;
}

template <typename Summand>
inline
SupVar<Summand> && operator*
( SupVar<Summand> && var1, double const& cst2 )
{
  var1 *= cst2;
  return std::move( var1 );
}

template <typename Summand>
inline
SupVar<Summand> operator*
( double const& cst1, SupVar<Summand> const& var2 )
{
  if( !var2._mod ) 
    return( cst1 * var2._cst );

  SupVar<Summand> var3( var2 );
  var3 *= cst1;
  return var3;
}

template <typename Summand>
inline
SupVar<Summand> && operator*
( double const& cst1, SupVar<Summand> && var2 )
{
  var2 *= cst1;
  return std::move( var2 );
}

template <typename Summand>
inline
SupVar<Summand>& SupVar<Summand>::operator/=
( double const& cst )
{
  if( cst == 0. )
    throw typename SupModel<Summand>::Exceptions( SupModel<Summand>::Exceptions::DIV );
  return( operator*=( 1/cst ) );
}

template <typename Summand>
inline
SupVar<Summand>& SupVar<Summand>::operator/=
( SupVar<Summand> const& var )
{
  return( operator*=( pow( var, -1 ) ) );
}

template <typename Summand>
inline
SupVar<Summand>& SupVar<Summand>::operator/=
( SupVar<Summand> && var )
{
  return( operator*=( pow( std::move(var), -1 ) ) );
}

template <typename Summand>
inline
SupVar<Summand> operator/
( SupVar<Summand> const& var1, SupVar<Summand> const& var2 )
{
  SupVar<Summand> var3( var1 );
  var3 /= var2;
  return var3;
}

template <typename Summand>
inline
SupVar<Summand> && operator/
( SupVar<Summand> const& var1, SupVar<Summand> && var2 )
{
  SupVar<Summand> var3( var1 );
  var3 /= std::move( var2 );
  return var3;
}

template <typename Summand>
inline
SupVar<Summand> && operator/
( SupVar<Summand> && var1, SupVar<Summand> const& var2 )
{
  var1 /= var2;
  return std::move( var1 );
}

template <typename Summand>
inline
SupVar<Summand> && operator/
( SupVar<Summand> && var1, SupVar<Summand> && var2 )
{
  var1 /= std::move( var2 );
  return std::move( var1 );
}

template <typename Summand>
inline
SupVar<Summand> operator/
( SupVar<Summand> const& var1, double const& cst2 )
{
  SupVar<Summand> var3( var1 );
  var3 /= cst2;
  return var3;
}

template <typename Summand>
inline
SupVar<Summand> && operator/
( SupVar<Summand> && var1, double const& cst2 )
{
  var1 /= cst2;
  return std::move( var1 );
}

template <typename Summand>
inline
SupVar<Summand> operator/
( double const& cst1, SupVar<Summand> const& var2 )
{
  SupVar<Summand> var3( pow( var2, -1 ) );
  var3 *= cst1;
  return var3;
}

template <typename Summand>
inline
SupVar<Summand> && operator/
( double const& cst1, SupVar<Summand> && var2 )
{
  pow( std::move(var2), -1 );
  var2 *= cst1;
  return std::move( var2 );
}

template <typename Summand>
inline
SupVar<Summand> inv
( SupVar<Summand> const& var )
{
  return pow( var, -1 );
}

template <typename Summand>
inline
SupVar<Summand> && inv
( SupVar<Summand> && var )
{
  return pow( std::move(var), -1 );
}

template <typename Summand>
inline
SupVar<Summand> sqr
( SupVar<Summand> const& var )
{
  return pow( var, 2 );
}

template <typename Summand>
inline
SupVar<Summand> && sqr
( SupVar<Summand> && var )
{
  return pow( std::move(var), 2 );
}

template <typename Summand>
inline
SupVar<Summand> sqrt
( SupVar<Summand> const& var )
{
  if ( var.l() < 0. )
    throw typename mc::SupModel<Summand>::Exceptions( mc::SupModel<Summand>::Exceptions::SQRT ); 

  if( !var._mod )
    return std::sqrt( var._cst );

#ifdef MC__SUPMODEL_TRACE
std::cout << "SupVar<Summand> sqrt: copy" << std::endl;
#endif

  auto const& f  = [=]( const double& x ){ return std::sqrt(x); };
  auto const& df = [=]( const double& x ){ return x>0.? 0.5/std::sqrt(x): 0.5/std::sqrt(DBL_MAX); };

  SupVar<Summand> var2( var );
  var2._mod->compose( var2, f, df, 0, 1 ); // concave & non-decreasing
  return var2;
}

template <typename Summand>
inline
SupVar<Summand> && sqrt
( SupVar<Summand> && var )
{
  if ( var.l() < -1e2*DBL_EPSILON ){//0. )
    std::cout << "var.l = " << var.l() << std::endl;
    throw typename mc::SupModel<Summand>::Exceptions( mc::SupModel<Summand>::Exceptions::SQRT ); 
  }
  
  if( !var._mod ){
    var._cst = std::sqrt( var._cst );
    return std::move(var);
  }

#ifdef MC__SUPMODEL_TRACE
std::cout << "SupVar<Summand> sqrt: move" << std::endl;
#endif

  auto const& f  = [=]( const double& x ){ return std::sqrt(x); };
  auto const& df = [=]( const double& x ){ return x>0.? 0.5/std::sqrt(x): 0.5/std::sqrt(DBL_MAX); };

  var._mod->compose( var, f, df, 0, 1 ); // concave & non-decreasing
  return std::move(var);
}

template <typename Summand>
inline
SupVar<Summand> exp
( SupVar<Summand> const& var )
{
  if( !var._mod )
    return std::exp( var._cst );

#ifdef MC__SUPMODEL_TRACE
std::cout << "SupVar<Summand> exp: copy" << std::endl;
#endif

  auto const& f = [=]( const double& x ){ return std::exp(x); };

  SupVar<Summand> var2( var );
  var2._mod->compose( var2, f, f, 1, 1 ); // convex & non-decreasing
  return var2;
}

template <typename Summand>
inline
SupVar<Summand> && exp
( SupVar<Summand> && var )
{
  if( !var._mod ){
    var._cst = std::exp( var._cst );
    return std::move(var);
  }

#ifdef MC__SUPMODEL_TRACE
std::cout << "SupVar<Summand> exp: move" << std::endl;
#endif

  auto const& f = [=]( const double& x ){ return std::exp(x); };

  var._mod->compose( var, f, f, 1, 1 ); // convex & non-decreasing
  return std::move(var);
}

template <typename Summand>
inline
SupVar<Summand> log
( SupVar<Summand> const& var )
{
  if ( var.l() <= 0. )
    throw typename mc::SupModel<Summand>::Exceptions( mc::SupModel<Summand>::Exceptions::LOG ); 

  if( !var._mod )
    return std::log( var._cst );

#ifdef MC__SUPMODEL_TRACE
std::cout << "SupVar<Summand> log: copy" << std::endl;
#endif

  auto const& f  = [=]( const double& x ){ return std::log(x); };
  auto const& df = [=]( const double& x ){ return 1./x; };

  SupVar<Summand> var2( var );
  var2._mod->compose( var2, f, df, 0, 1 ); // concave & non-decreasing
  return var2;
}

template <typename Summand>
inline
SupVar<Summand> && log
( SupVar<Summand> && var )
{
  if ( var.l() <= 0. )
    throw typename mc::SupModel<Summand>::Exceptions( mc::SupModel<Summand>::Exceptions::LOG ); 

  if( !var._mod ){
    var._cst = std::log( var._cst );
    return std::move(var);
  }

#ifdef MC__SUPMODEL_TRACE
std::cout << "SupVar<Summand> log: move" << std::endl;
#endif

  auto const& f  = [=]( const double& x ){ return std::log(x); };
  auto const& df = [=]( const double& x ){ return 1./x; };

  var._mod->compose( var, f, df, 0, 1 ); // concave & non-decreasing
  return std::move(var);
}

template <typename Summand>
inline
SupVar<Summand> xlog
( SupVar<Summand> const& var )
{
  if ( var.l() < 0. )
    throw typename mc::SupModel<Summand>::Exceptions( mc::SupModel<Summand>::Exceptions::LOG ); 

  if( !var._mod )
    return var._cst>0.? var._cst*std::log( var._cst ): 0.;

#ifdef MC__SUPMODEL_TRACE
std::cout << "SupVar<Summand> xlog: copy" << std::endl;
#endif

  if( var.l() >= std::exp(-1) ){
    auto const& f  = [=]( const double& x ){ return x>0.? x*std::log(x): 0.; };
    auto const& df = [=]( const double& x ){ return x>0.? std::log(x)+1.: std::log(DBL_MIN)+1.; };

    SupVar<Summand> var2( var );
    var2._mod->compose( var2, f, df, 1, 0 ); // convex & non-increasing
    return var;
  }

  if( var.u() <= std::exp(-1) ){
    auto const& f  = [=]( const double& x ){ return x>0.? x*std::log(x): 0.; };
    auto const& df = [=]( const double& x ){ return x>0.? std::log(x)+1.: std::log(DBL_MIN)+1.; };

    SupVar<Summand> var2( var );
    var2._mod->compose( var2, f, df, 1, 1 ); // convex & non-decreasing
    return var;
  }

  SupVar<Summand> varl( var ), varr( var );
  auto const& fl  = [=]( const double& x ){ return x>0. && x<std::exp(-1.)? x*std::log(x)+std::exp(-1.): 0.; };
  auto const& dfl = [=]( const double& x ){ return x>0.? (x<std::exp(-1.)? std::log(x)+1.: 0.): std::log(DBL_MIN)+1.; };
  auto const& fr  = [=]( const double& x ){ return x>std::exp(-1.)? x*std::log(x): -std::exp(-1.); };
  auto const& dfr = [=]( const double& x ){ return x>std::exp(-1.)? std::log(x)+1.: 0.; };

  var._mod->compose( varr, fr, dfr, 1, 1 ); // convex & non-decreasing
  var._mod->compose( varl, fl, dfl, 1, 0 ); // convex & non-increasing
  return varr += varl;
}

template <typename Summand>
inline
SupVar<Summand> && xlog
( SupVar<Summand> && var )
{
  if ( var.l() < 0. )
    throw typename mc::SupModel<Summand>::Exceptions( mc::SupModel<Summand>::Exceptions::LOG ); 

  if( !var._mod ){
    if( var._cst > 0. ){
      var._cst *= std::log( var._cst );
      return std::move( var );
    }
  }
  
#ifdef MC__SUPMODEL_TRACE
std::cout << "SupVar<Summand> xlog: move" << std::endl;
#endif

  if( var.l() >= std::exp(-1) ){
    auto const& f  = [=]( const double& x ){ return x>0.? x*std::log(x): 0.; };
    auto const& df = [=]( const double& x ){ return x>0.? std::log(x)+1.: std::log(DBL_MIN)+1.; };

    var._mod->compose( var, f, df, 1, 0 ); // convex & non-increasing
    return std::move( var );
  }

  if( var.u() <= std::exp(-1) ){
    auto const& f  = [=]( const double& x ){ return x>0.? x*std::log(x): 0.; };
    auto const& df = [=]( const double& x ){ return x>0.? std::log(x)+1.: std::log(DBL_MIN)+1.; };

    var._mod->compose( var, f, df, 1, 1 ); // convex & non-decreasing
    return std::move( var );
  }

  SupVar<Summand> varl( var );
  auto const& fl  = [=]( const double& x ){ return x>0. && x<std::exp(-1.)? x*std::log(x)+std::exp(-1.): 0.; };
  auto const& dfl = [=]( const double& x ){ return x>0.? (x<std::exp(-1.)? std::log(x)+1.: 0.): std::log(DBL_MIN)+1.; };
  auto const& fr  = [=]( const double& x ){ return x>std::exp(-1.)? x*std::log(x): -std::exp(-1.); };
  auto const& dfr = [=]( const double& x ){ return x>std::exp(-1.)? std::log(x)+1.: 0.; };

  var._mod->compose( var,  fr, dfr, 1, 1 ); // convex & non-decreasing
  var._mod->compose( varl, fl, dfl, 1, 0 ); // convex & non-increasing
  return std::move( var += std::move( varl ) );
}

template <typename Summand>
inline
SupVar<Summand> tanh
( SupVar<Summand> const& var )
{
  if( !var._mod )
    return std::tanh( var._cst );

#ifdef MC__SUPMODEL_TRACE
std::cout << "SupVar<Summand> tanh: copy" << std::endl;
#endif

  if( var.l() >= 0 ){
    SupVar<Summand> var2( var );
    auto const& f  = [=]( const double& x ){ return std::tanh(x); };
    auto const& df = [=]( const double& x ){ double z=std::tanh(x); return 1.-z*z; };

    var2._mod->compose( var2, f, df, 0, 1 ); // concave & non-decreasing  
    return var2;
  }

  if( var.u() <= 0 ){
    SupVar<Summand> var2( var );
    auto const& f  = [=]( const double& x ){ return std::tanh(x); };
    auto const& df = [=]( const double& x ){ double z=std::tanh(x); return 1.-z*z; };

    var2._mod->compose( var2, f, df, 1, 1 ); // convex & non-decreasing  
    return var2;
  }

  SupVar<Summand> varl( var ), varr( var );
  auto const& fl  = [=]( const double& x ){ double z=std::tanh(x); return z>x?z:x; };
  auto const& dfl = [=]( const double& x ){ double z=std::tanh(x); return z>x?1.-z*z:1.; };
  auto const& fr  = [=]( const double& x ){ double z=std::tanh(x); return z<x?z-x:0.; };
  auto const& dfr = [=]( const double& x ){ double z=std::tanh(x); return z<x?-z*z:0.; };

  varr._mod->compose( varr, fl, dfl, 1, 1 ); // convex & non-decreasing
  varl._mod->compose( varl, fr, dfr, 0, 0 ); // concave & non-increasing
  return varr += varl;
}

template <typename Summand>
inline
SupVar<Summand> && tanh
( SupVar<Summand> && var )
{
  if( !var._mod ){
    var._cst = std::tanh( var._cst );
    return std::move( var );
  }
  
#ifdef MC__SUPMODEL_TRACE
std::cout << "SupVar<Summand> tanh: move" << std::endl;
#endif

  if( var.l() >= 0 ){
    auto const& f  = [=]( const double& x ){ return std::tanh(x); };
    auto const& df = [=]( const double& x ){ double z=std::tanh(x); return 1.-z*z; };

    var._mod->compose( var, f, df, 0, 1 ); // concave & non-decreasing  
    return std::move( var );
  }

  if( var.u() <= 0 ){
    auto const& f  = [=]( const double& x ){ return std::tanh(x); };
    auto const& df = [=]( const double& x ){ double z=std::tanh(x); return 1.-z*z; };

    var._mod->compose( var, f, df, 1, 1 ); // convex & non-decreasing  
    return std::move( var );
  }

  SupVar<Summand> varr( var );
  auto const& fl  = [=]( const double& x ){ double z=std::tanh(x); return z>x?z:x; };
  auto const& dfl = [=]( const double& x ){ double z=std::tanh(x); return z>x?1.-z*z:1.; };
  auto const& fr  = [=]( const double& x ){ double z=std::tanh(x); return z<x?z-x:0.; };
  auto const& dfr = [=]( const double& x ){ double z=std::tanh(x); return z<x?-z*z:0.; };

  var._mod->compose( var,  fl, dfl, 1, 1 ); // convex & non-decreasing
  var._mod->compose( varr, fr, dfr, 0, 0 ); // concave & non-increasing
  return std::move( var += varr );
}

template <typename Summand>
inline
SupVar<Summand>
pow
( SupVar<Summand> const& var, int const& n )
{
  if( !n )
    return 1.; 

  if( n == 1 )
    return var; 

  if( !var._mod )
    return std::pow( var._cst, n );

#ifdef MC__SUPMODEL_TRACE
std::cout << "SupVar<Summand> pow: copy" << std::endl;
#endif

  if( n < 0 ){
    if( var.l() * var.u() <= 0. )
      throw typename mc::SupModel<Summand>::Exceptions( mc::SupModel<Summand>::Exceptions::INV ); 

    SupVar<Summand> var2( var );
    auto const& f  = [=]( const double& x ){ return std::pow(x,n); };
    auto const& df = [=]( const double& x ){ return (double)n*std::pow(x,n-1); };

    if( var.l() > 0 )
      var2._mod->compose( var2, f, df, 1, 0 ); // convex & non-increasing
    else
      var2._mod->compose( var2, f, df, n%2?0:1, n%2?0:1 ); // concave/convex & non-in/decreasing
    return var2;
  }

  if( var.l() >= 0 ){
    SupVar<Summand> var2( var );
    auto const& f  = [=]( const double& x ){ return std::pow(x,n); };
    auto const& df = [=]( const double& x ){ return (double)n*std::pow(x,n-1); };

    var2._mod->compose( var2, f, df, 1, 1 ); // convex & non-decreasing  
    return var2;
  }

  if( var.u() <= 0 ){
    SupVar<Summand> var2( var );
    auto const& f  = [=]( const double& x ){ return std::pow(x,n); };
    auto const& df = [=]( const double& x ){ return (double)n*std::pow(x,n-1); };

    var2._mod->compose( var2, f, df, n%2?0:1, n%2?1:0 ); // concave/convex & non-de/increasing  
    return var2;
  }
   
  SupVar<Summand> varl( var ), varr( var );
  auto const& fl  = [=]( const double& x ){ return std::pow(x<0.?x:0.,n); };
  auto const& dfl = [=]( const double& x ){ return (double)n*std::pow(x<0.?x:0.,n-1); };
  auto const& fr  = [=]( const double& x ){ return std::pow(x>0.?x:0.,n); };
  auto const& dfr = [=]( const double& x ){ return (double)n*std::pow(x>0.?x:0.,n-1); };

  varr._mod->compose( varr, fr, dfr, 1, 1 ); // convex & non-decreasing
  varl._mod->compose( varl, fl, dfl, n%2?0:1, n%2?1:0 ); // concave/convex & non-de/increasing
  return varr += varl;
}

template <typename Summand>
inline
SupVar<Summand> &&
pow
( SupVar<Summand> && var, int const& n )
{
  if( !n ){
    var._cst = 1.;
    return std::move( var );
  }
  
  if( n == 1 )
    return std::move( var );

  if( !var._mod ){
    var._cst = std::pow( var._cst, n );
    return std::move( var );
  }

#ifdef MC__SUPMODEL_TRACE
std::cout << "SupVar<Summand> pow: move" << std::endl;
#endif

  if( n < 0 ){
    if( var.l() * var.u() <= 0. )
      throw typename mc::SupModel<Summand>::Exceptions( mc::SupModel<Summand>::Exceptions::INV ); 

    auto const& f  = [=]( const double& x ){ return std::pow(x,n); };
    auto const& df = [=]( const double& x ){ return (double)n*std::pow(x,n-1); };

    if( var.l() > 0 )
      var._mod->compose( var, f, df, 1, 0 ); // convex & non-increasing
    else
      var._mod->compose( var, f, df, n%2?0:1, n%2?0:1 ); // concave/convex & non-in/decreasing
    return std::move( var );
  }
  
  if( var.l() >= 0 ){
    auto const& f  = [=]( const double& x ){ return std::pow(x,n); };
    auto const& df = [=]( const double& x ){ return (double)n*std::pow(x,n-1); };

    var._mod->compose( var, f, df, 1, 1 ); // convex & non-decreasing  
    return std::move( var );
  }

  if( var.u() <= 0 ){
    auto const& f  = [=]( const double& x ){ return std::pow(x,n); };
    auto const& df = [=]( const double& x ){ return (double)n*std::pow(x,n-1); };

    var._mod->compose( var, f, df, n%2?0:1, n%2?1:0 ); // concave/convex & non-de/increasing  
    return std::move( var );
  }
   
  SupVar<Summand> varl( var );
  auto const& fl  = [=]( const double& x ){ return std::pow(x<0.?x:0.,n); };
  auto const& dfl = [=]( const double& x ){ return (double)n*std::pow(x<0.?x:0.,n-1); };
  auto const& fr  = [=]( const double& x ){ return std::pow(x>0.?x:0.,n); };
  auto const& dfr = [=]( const double& x ){ return (double)n*std::pow(x>0.?x:0.,n-1); };

  var._mod->compose( var,  fr, dfr, 1, 1 ); // convex & non-decreasing
  var._mod->compose( varl, fl, dfl, n%2?0:1, n%2?1:0 ); // concave/convex & non-de/increasing
  return std::move( var += varl );
}

template <typename Summand>
inline
SupVar<Summand> pow
( SupVar<Summand> const& var, double const& d )
{
  if( d == 0. )
    return 1.; 

  if( d == 1. )
    return var;

  if( var.l() < 0. || ( var.l() == 0. && d < 0. ) )
    throw typename mc::SupModel<Summand>::Exceptions( mc::SupModel<Summand>::Exceptions::RPOW ); 

  if( !var._mod )
    return std::pow( var._cst, d );

#ifdef MC__SUPMODEL_TRACE
std::cout << "SupVar<Summand> pow( SupVar<Summand> const&, double const& )" << std::endl;
#endif

  auto const& f  = [=]( const double& x ){ return std::pow(x,(double)d); };
  auto const& df = [=]( const double& x ){ return x>0.? d*std::pow(x,d-1.): d*std::pow(DBL_MAX,d-1.); };

  SupVar<Summand> var2( var );
  if( d > 1. )
    var2._mod->compose( var2, f, df, 1, 1 ); // convex & non-decreasing
  else if( d > 0. )
    var2._mod->compose( var2, f, df, 0, 1 ); // concave & non-decreasing
  else
    var2._mod->compose( var2, f, df, 1, 0 ); // convex & non-increasing
  return var2;
}

template <typename Summand>
inline
SupVar<Summand> && pow
( SupVar<Summand> && var, double const& d )
{
  if( d == 0. ){
    var._cst = 1.;
    return std::move( var );
  }
  
  if( d == 1 )
    return std::move( var );

  if( var.l() < 0. || ( var.l() == 0. && d < 0. ) )
    throw typename mc::SupModel<Summand>::Exceptions( mc::SupModel<Summand>::Exceptions::RPOW ); 

  if( !var._mod ){
    var._cst = std::pow( var._cst, d );
    return std::move( var );
  }

#ifdef MC__SUPMODEL_TRACE
std::cout << "SupVar<Summand> && pow( SupVar<Summand> &&, double const& )" << std::endl;
#endif

  auto const& f  = [=]( const double& x ){ return std::pow(x,(double)d); };
  auto const& df = [=]( const double& x ){ return x>0.? d*std::pow(x,d-1.): d*std::pow(DBL_MAX,d-1.); };

   if( d > 1. )
    var._mod->compose( var, f, df, 1, 1 ); // convex & non-decreasing
  else if( d > 0. )
    var._mod->compose( var, f, df, 0, 1 ); // concave & non-decreasing
  else
    var._mod->compose( var, f, df, 1, 0 ); // convex & non-increasing
 return std::move(var);
}

template <typename Summand>
inline
SupVar<Summand>
pow
( double const& d, SupVar<Summand> const& var )
{
  return exp( var * std::log( d ) );
}

template <typename Summand>
inline
SupVar<Summand>
pow
( double const& d, SupVar<Summand> && var )
{
  return exp( var * std::log( d ) );
}

template <typename Summand>
inline
SupVar<Summand>
pow
( SupVar<Summand> const& var, SupVar<Summand> const& vexp )
{
  return exp( vexp * log( var ) );
}

template <typename Summand>
inline
SupVar<Summand>
pow
( SupVar<Summand> && var, SupVar<Summand> const& vexp )
{
  return exp( vexp * log( var ) );
}

template <typename Summand>
inline
SupVar<Summand>
prod
( unsigned int const& n, SupVar<Summand> const* var )
{
  switch( n ){
   case 0:  return 1.;
   case 1:  return var[0];
   default: return var[0] * prod( n-1, var+1 );
  }
}

template <typename Summand>
inline
SupVar<Summand>
monom
( unsigned int const& n, SupVar<Summand> const* var, unsigned int const* exp )
{
  switch( n ){
   case 0:  return 1.;
   case 1:  return pow( var[0], (int)exp[0] );
   default: return pow( var[0], (int)exp[0] ) * monom( n-1, var+1, exp+1 );
  }
}

template <typename Summand>
inline
SupVar<Summand>
cheb
( SupVar<Summand> const& var, unsigned int const& n )
{
  if( !var._mod )
    return cheb( var._cst, n );

  switch( n ){
    case 0:  return 1.;
    case 1:  return var;
    default: break;
  }

  SupVar<Summand> var2( 2.*var*cheb(var,n-1) - cheb(var,n-2) );
  return var2;
}

template <typename Summand>
inline
SupVar<Summand> fabs
( SupVar<Summand> const& var )
{
  if( var.l() >= 0 ){
    return var;
  }

  if( var.u() <= 0 ){
    return -var;
  }

  SupVar<Summand> varl( var ), varr( var );
  auto const& fl  = [=]( const double& x ){ return x<0.?-x:0.; };
  auto const& dfl = [=]( const double& x ){ return x<0.?-1.:0.; };
  auto const& fr  = [=]( const double& x ){ return x>0.?x:0.; };
  auto const& dfr = [=]( const double& x ){ return x>0.?1.:0.; };

  varr._mod->compose( varr, fr, dfr, 1, 1 ); // convex & non-decreasing
  varl._mod->compose( varl, fl, dfl, 1, 0 ); // convex & non-increasing
  return varr += varl;
}

template <typename Summand>
inline
SupVar<Summand> && fabs
( SupVar<Summand> && var )
{
  if( var.l() >= 0 ){
    return std::move( var );
  }

  if( var.u() <= 0 ){
    return operator-( std::move( var ) );
  }

  SupVar<Summand> varl( var );
  auto const& fl  = [=]( const double& x ){ return x<0.?-x:0.; };
  auto const& dfl = [=]( const double& x ){ return x<0.?-1.:0.; };
  auto const& fr  = [=]( const double& x ){ return x>0.?x:0.; };
  auto const& dfr = [=]( const double& x ){ return x>0.?1.:0.; };

  var._mod->compose( var, fr, dfr, 1, 1 ); // convex & non-decreasing
  varl._mod->compose( varl, fl, dfl, 1, 0 ); // convex & non-increasing
  return std::move( var += varl );
}

template <typename Summand>
inline
SupVar<Summand>
max
( SupVar<Summand> const& var1, SupVar<Summand> const& var2 )
{
  if( !var1._mod && !var2._mod )
    return var1._cst > var2._cst? var1._cst: var2._cst;

  if( !var2._mod )
    return max( var1, var2._cst );

  if( !var1._mod )
    return max( var2, var1._cst );

  double const& w = var1._mod->options.REF_WEIGHT;
  if( w <= 0. )
    return max( var2 - var1, 0. ) + var1;
  else if( w >= 1. )
    return max( var1 - var2, 0. ) + var2;
  else
    return (1 - w ) * ( max( var2 - var1, 0. ) + var1 ) + w * ( max( var1 - var2, 0. ) + var2 );
}

template <typename Summand>
inline
SupVar<Summand>
max
( SupVar<Summand> const& var1, double const& cst2 )
{
  SupVar<Summand> var2( var1 );
  return max( std::move( var2 ), cst2 ); 
}

template <typename Summand>
inline
SupVar<Summand> &&
max
( SupVar<Summand> && var1, double const& cst2 )
{
  if( !var1._mod ){
    if( var1._cst < cst2 ) var1._cst = cst2;
    return std::move( var1 );
  }

#ifdef MC__SUPMODEL_TRACE
std::cout << "SupVar<Summand> && max( SupVar<Summand> &&, double const& )" << std::endl;
#endif

  var1._mod->max( var1, cst2 );
  return std::move( var1 );
}

template <typename Summand>
inline
SupVar<Summand>
max
( double const& cst1, SupVar<Summand> const& var2 )
{
  return max( var2, cst1 );
}

template <typename Summand>
inline
SupVar<Summand> &&
max
( double const& cst1, SupVar<Summand> && var2 )
{
  return max( std::move( var2 ), cst1 );
}

template <typename Summand>
inline
SupVar<Summand>
relu
( SupVar<Summand> const& var )
{
  return max( var, 0. );
}

template <typename Summand>
inline
SupVar<Summand> &&
relu
( SupVar<Summand> && var )
{
  return max( std::move( var ), 0. );
}

template <typename Summand>
inline
SupVar<Summand>
min
( SupVar<Summand> const& var1, SupVar<Summand> const& var2 )
{
  if( !var1._mod && !var2._mod )
    return var1._cst < var2._cst? var1._cst: var2._cst;

  if( !var2._mod )
    return min( var1, var2._cst );

  if( !var1._mod )
    return min( var2, var1._cst );

  double const& w = var1._mod->options.REF_WEIGHT;
  if( w <= 0. )
    return min( var2 - var1, 0. ) + var1;
  else if( w >= 1. )
    return min( var1 - var2, 0. ) + var2;
  else
    return (1 - w ) * ( min( var2 - var1, 0. ) + var1 ) + w * ( min( var1 - var2, 0. ) + var2 );
}

template <typename Summand>
inline
SupVar<Summand>
min
( SupVar<Summand> const& var1, double const& cst2 )
{
  SupVar<Summand> var2( var1 );
  return min( std::move( var2 ), cst2 ); 
}

template <typename Summand>
inline
SupVar<Summand> &&
min
( SupVar<Summand> && var1, double const& cst2 )
{
  if( !var1._mod ){
    if( var1._cst > cst2 ) var1._cst = cst2;
    return std::move( var1 );
}

#ifdef MC__SUPMODEL_TRACE
std::cout << "SupVar<Summand> && min( SupVar<Summand> &&, double const& )" << std::endl;
#endif

  var1._mod->min( var1,  cst2 );
  return std::move( var1 );
}

template <typename Summand>
inline
SupVar<Summand>
min
( double const& cst1, SupVar<Summand> const& var2 )
{
  SupVar<Summand> var1( var2 );
  return min( std::move( var1 ), cst1 );
}

template <typename Summand>
inline
SupVar<Summand> &&
min
( double const& cst1, SupVar<Summand> && var2 )
{
  return min( std::move( var2 ), cst1 );
}
/*
template <typename Summand>
inline
bool inter
( SupVar<Summand>& var12, SupVar<Summand> const& var1, SupVar<Summand> const& var2 )
{

}

template <typename Summand>
inline
SupVar<Summand> hull
( SupVar<Summand> const& var1, SupVar<Summand> const& var2 )
{

}
*/
} // namespace mc


#include "mcfadbad.hpp"

namespace fadbad
{

//! @brief Specialization of the structure fadbad::Op for use of the type mc::SupVar of MC++ as a template parameter of the classes fadbad::F, fadbad::B and fadbad::T of FADBAD++
template <typename Summand> struct Op<mc::SupVar<Summand>>
{
  typedef mc::SupVar<Summand> SV;
  typedef double Base;
  static Base myInteger( int const i ) { return Base(i); }
  static Base myZero() { return myInteger(0); }
  static Base myOne()  { return myInteger(1);}
  static Base myTwo()  { return myInteger(2); }
  static double myPI() { return mc::PI; }
  static SV myPos( SV const& x ) { return  x; }
  static SV myNeg( SV const& x ) { return -x; }
  template <typename U> static SV& myCadd( SV& x, U const& y ) { return x += y; }
  template <typename U> static SV& myCsub( SV& x, U const& y ) { return x -= y; }
  template <typename U> static SV& myCmul( SV& x, U const& y ) { return x *= y; }
  template <typename U> static SV& myCdiv( SV& x, U const& y ) { return x /= y; }
  static SV myInv( SV const& x ) { return mc::inv( x ); }
  static SV mySqr( SV const& x ) { return mc::sqr( x ); }
  template <typename X, typename Y> static SV myPow( X const& x, Y const& y ) { return mc::pow( x, y ); }
  static SV mySqrt( SV const& x ){ return mc::sqrt(x); }
  static SV myLog( SV const& x ) { return mc::log( x ); }
  static SV myExp( SV const& x ) { return mc::exp( x ); }
  static SV mySin( SV const& x ) { throw typename mc::SupModel<Summand>::Exceptions( mc::SupModel<Summand>::Exceptions::UNDEF ); } //{ return mc::sin( x ); }
  static SV myCos( SV const& x ) { throw typename mc::SupModel<Summand>::Exceptions( mc::SupModel<Summand>::Exceptions::UNDEF ); } //{ return mc::cos( x ); }
  static SV myTan( SV const& x ) { throw typename mc::SupModel<Summand>::Exceptions( mc::SupModel<Summand>::Exceptions::UNDEF ); } //{ return mc::tan( x ); }
  static SV myAsin( SV const& x ) { throw typename mc::SupModel<Summand>::Exceptions( mc::SupModel<Summand>::Exceptions::UNDEF ); }
  static SV myAcos( SV const& x ) { throw typename mc::SupModel<Summand>::Exceptions( mc::SupModel<Summand>::Exceptions::UNDEF ); }
  static SV myAtan( SV const& x ) { throw typename mc::SupModel<Summand>::Exceptions( mc::SupModel<Summand>::Exceptions::UNDEF ); }
  static SV mySinh( SV const& x ) { throw typename mc::SupModel<Summand>::Exceptions( mc::SupModel<Summand>::Exceptions::UNDEF ); }
  static SV myCosh( SV const& x ) { throw typename mc::SupModel<Summand>::Exceptions( mc::SupModel<Summand>::Exceptions::UNDEF ); }
  static SV myTanh( SV const& x ) { return mc::tanh( x ); }
  static bool myEq( SV const& x, SV const& y ) { throw typename mc::SupModel<Summand>::Exceptions( mc::SupModel<Summand>::Exceptions::UNDEF ); }
  static bool myNe( SV const& x, SV const& y ) { throw typename mc::SupModel<Summand>::Exceptions( mc::SupModel<Summand>::Exceptions::UNDEF ); }
  static bool myLt( SV const& x, SV const& y ) { throw typename mc::SupModel<Summand>::Exceptions( mc::SupModel<Summand>::Exceptions::UNDEF ); }
  static bool myLe( SV const& x, SV const& y ) { throw typename mc::SupModel<Summand>::Exceptions( mc::SupModel<Summand>::Exceptions::UNDEF ); }
  static bool myGt( SV const& x, SV const& y ) { throw typename mc::SupModel<Summand>::Exceptions( mc::SupModel<Summand>::Exceptions::UNDEF ); }
  static bool myGe( SV const& x, SV const& y ) { throw typename mc::SupModel<Summand>::Exceptions( mc::SupModel<Summand>::Exceptions::UNDEF ); }
};

} // end namespace fadbad

namespace mc
{

//! @brief C++ structure for specialization of the mc::Op templated structure for use of mc::SupVar in DAG evaluation and as template parameter in other MC++ types
template <typename Summand> struct Op<mc::SupVar<Summand>>
{
  typedef mc::SupVar<Summand> SV;
  static SV point( const double c ) { throw typename SupModel<Summand>::Exceptions( SupModel<Summand>::Exceptions::UNDEF ); }
  static SV zeroone() { throw typename SupModel<Summand>::Exceptions( SupModel<Summand>::Exceptions::UNDEF ); }
  static void I(SV& x, SV const&y) { x = y; }
  static double l(SV const& x) { return x.l(); }
  static double u(SV const& x) { return x.u(); }
  static double abs (SV const& x) { return std::fabs(x.u())>std::fabs(x.l())? std::fabs(x.u()): std::fabs(x.l()); }
  static double mid (SV const& x) { return 0.5*(std::fabs(x.l())+std::fabs(x.u())); }
  static double diam(SV const& x) { return x.u()-x.l(); }
  static SV inv (SV const& x) { return mc::inv(x);  }
  static SV sqr (SV const& x) { return mc::sqr(x);  }
  static SV sqrt(SV const& x) { return mc::sqrt(x); }
  static SV exp (SV const& x) { return mc::exp(x);  }
  static SV log (SV const& x) { return mc::log(x);  }
  static SV xlog(SV const& x) { return mc::xlog(x); }
  static SV lmtd(SV const& x, SV const& y) { return (x-y)/(mc::log(x)-mc::log(y)); }
  static SV rlmtd(SV const& x, SV const& y) { return (mc::log(x)-mc::log(y))/(x-y); }
  static SV fabs(SV const& x) { return mc::fabs(x); }
  static SV sin (SV const& x) { throw typename SupModel<Summand>::Exceptions( SupModel<Summand>::Exceptions::UNDEF ); }
  static SV cos (SV const& x) { throw typename SupModel<Summand>::Exceptions( SupModel<Summand>::Exceptions::UNDEF ); }
  static SV tan (SV const& x) { throw typename SupModel<Summand>::Exceptions( SupModel<Summand>::Exceptions::UNDEF ); }
  static SV asin(SV const& x) { throw typename SupModel<Summand>::Exceptions( SupModel<Summand>::Exceptions::UNDEF ); }
  static SV acos(SV const& x) { throw typename SupModel<Summand>::Exceptions( SupModel<Summand>::Exceptions::UNDEF ); }
  static SV atan(SV const& x) { throw typename SupModel<Summand>::Exceptions( SupModel<Summand>::Exceptions::UNDEF ); }
  static SV sinh(SV const& x) { throw typename SupModel<Summand>::Exceptions( SupModel<Summand>::Exceptions::UNDEF ); }
  static SV cosh(SV const& x) { throw typename SupModel<Summand>::Exceptions( SupModel<Summand>::Exceptions::UNDEF ); }
  static SV tanh(SV const& x) { return mc::tanh(x); }
  static SV erf (SV const& x) { throw typename SupModel<Summand>::Exceptions( SupModel<Summand>::Exceptions::UNDEF ); }
  static SV erfc(SV const& x) { throw typename SupModel<Summand>::Exceptions( SupModel<Summand>::Exceptions::UNDEF ); }
  static SV fstep(SV const& x) { throw typename SupModel<Summand>::Exceptions( SupModel<Summand>::Exceptions::UNDEF ); }
  static SV bstep(SV const& x) { throw typename SupModel<Summand>::Exceptions( SupModel<Summand>::Exceptions::UNDEF ); }
  template <typename Y> static SV min (SV const& x, Y const& y) { return mc::min(x,y); }
  template <typename Y> static SV max (SV const& x, Y const& y) { return mc::max(x,y); }
  template <typename X, typename Y> static SV pow(const X& x, const Y& y) { return mc::pow(x,y); }
  static SV cheb(SV const& x, const unsigned n) { return mc::cheb(x,n); }
  static SV prod (const unsigned n, const SV* x) { return mc::prod(n,x); }
  static SV monom (const unsigned n, const SV* x, const unsigned* k) { return mc::monom(n,x,k); }
  static bool inter(SV& xIy, SV const& x, SV const& y) { throw typename SupModel<Summand>::Exceptions( SupModel<Summand>::Exceptions::UNDEF ); } //{ return mc::inter(xIy,x,y); }
  static SV hull(SV const& x, SV const& y) { throw typename SupModel<Summand>::Exceptions( SupModel<Summand>::Exceptions::UNDEF ); } //{ return mc::hull(x,y); }
  static bool eq(SV const& x, SV const& y) { throw typename mc::SupModel<Summand>::Exceptions( mc::SupModel<Summand>::Exceptions::UNDEF ); }
  static bool ne(SV const& x, SV const& y) { throw typename mc::SupModel<Summand>::Exceptions( mc::SupModel<Summand>::Exceptions::UNDEF ); }
  static bool lt(SV const& x, SV const& y) { throw typename mc::SupModel<Summand>::Exceptions( mc::SupModel<Summand>::Exceptions::UNDEF ); }
  static bool le(SV const& x, SV const& y) { throw typename mc::SupModel<Summand>::Exceptions( mc::SupModel<Summand>::Exceptions::UNDEF ); }
  static bool gt(SV const& x, SV const& y) { throw typename mc::SupModel<Summand>::Exceptions( mc::SupModel<Summand>::Exceptions::UNDEF ); }
  static bool ge(SV const& x, SV const& y) { throw typename mc::SupModel<Summand>::Exceptions( mc::SupModel<Summand>::Exceptions::UNDEF ); }
};

} // namespace mc

#endif
