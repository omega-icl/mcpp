// Copyright (C) Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the Eclipse Public License.

// To test the numerical-issue fixer
//#define MC__UPWLE_COMPUTATION_TOL  1e-3
//#define MC__UPWLE_SMOOTH_TOL 1e-10  
#define MC__UPWLE_COMPUTATION_TOL  1e-14
#define MC__UPWLE_SMOOTH_TOL 1e-13
//#define MC__UPWLE_DEBUG
//#define MC__UPWLE_DEBUG_TRACE
/*!
\page page_UModels Univariate Estimator Model Arithmetic for Univariate Factorable Functions
\author Yanlin Zha, Beno&icirc;t Chachuat

Univariate estimator model arithmetic is a computational method to construct non-convex enclosures in the form of continuous piecewise linear enclosures over the entire domain. (to be modified) Formally, 

The classes mc::UnivarPWL (and etc.) provide an implementation of univariate (over- or under-) estimator model arithmetic based on the operator/function overloading mechanism of C++. This makes univarate estimator arithmetic both simple and intuitive to compute, similar to computing function values in real arithmetics or function bounds in interval arithmetic (see \ref page_INTERVAL). mc::UnivarPWL (and etc.) are templated in the type used to propagate the double or long double coefficients. 
In addition, mc::UnivarPWL can be used as the template parameter of other available types in MC++; as well as types in <A href="http://www.fadbad.com/fadbad.html">FADBAD++</A> for computing  of the derivatives or Taylor coefficients of a univariate factorable function.


\section sec_UPWLE_use How do I compute an UnivarPWL Model of a univariate factorable function?

Suppose we want to compute an UnivarPWL for the real-valued function (to be added and modified)

\code
      #include "interval.hpp"
      #include "ismodel.hpp"
      typedef mc::Interval I;
      typedef mc::ISModel<I> ISM;
      typedef mc::ISVar<I> ISV;
\endcode

First, the number of independent variables in the factorable function (\f$x\f$ and \f$y\f$ here) as well as the partition size of the ISM (10 partitions here) are specified by defining an mc::ISModel object as:

\code
      ISM mod( 2, 8 );
\endcode

Next, the variables \f$x\f$ and \f$y\f$ are defined as follows:

\code
      ISV X( &mod, 0, I(1.,2.) );
      ISV Y( &mod, 1, I(0.,1.) );
\endcode

Essentially, the first line means that <tt>X</tt> is a variable of class mc::ISVar, participating in the ISM <tt>mod</tt>, with range \f$[1,2]\f$ and index 0 (using the convention that C/C++ indexing starts at 0). The same holds for the Chebyshev variable <tt>Y</tt>, participating in the model <tt>mod</tt>, with range \f$[0,1]\f$ and index 1.

Having defined both variables, an ISM of \f$f(x,y)=x\exp(x+y^2)-y^2\f$ on \f$[1,2]\times[0,1]\f$ is simply computed as:

\code
      auto F = X*exp(X+pow(Y,2))-pow(Y,2);
\endcode

This model can be displayed to the standard output as:

\code
      std::cout << "ISM of f:\n" << F;
\endcode

which produces the following output:

\verbatim
ISM of f:
    0: [ -8.09107e-01 :  5.38394e+00 ][  4.54199e-01 :  7.17827e+00 ][  1.99772e+00 :  9.30281e+00 ]
       [  3.87146e+00 :  1.18159e+01 ][  6.13375e+00 :  1.47856e+01 ][  8.85265e+00 :  1.82912e+01 ]
       [  1.21074e+01 :  2.24248e+01 ][  1.59903e+01 :  2.72940e+01 ]
    1: [ -5.84612e+00 : -2.89623e-02 ][ -5.77362e+00 :  3.24938e-01 ][ -5.48222e+00 :  9.33765e-01 ]
       [ -4.93589e+00 :  1.86409e+00 ][ -4.06806e+00 :  3.22442e+00 ][ -2.77024e+00 :  5.18473e+00 ]
       [ -8.72422e-01 :  8.00930e+00 ][  1.88965e+00 :  1.21115e+01 ]
    B: [ -6.65522e+00 :  3.94054e+01 ]
\endverbatim

A plot of the resulting ISM enclosure (interval boxes) is shown on the figure below for the function \f$f\f$ (color map).

<CENTER><TABLE BORDER=0>
<TR>
<TD>\image html ISM-2D.png</TD>
</TR>
</TABLE></CENTER>

Other operations involve retreiving the interval coefficient matrix:

\code
      auto&& mat = F.C();
\endcode

See the documentations of mc::ISModel and mc::ISVar for a complete list of member functions. 


\section sec_ISM_err Errors Errors encountered during computation of an ISM

Errors are managed based on the exception handling mechanism of the C++ language. Each time an error is encountered, a class object of type mc::SCModel::Exceptions is thrown, which contains the type of error. It is the user's responsibility to test whether an exception was thrown during the computation of a Chebyshev model, and then make the appropriate changes. Should an exception be thrown and not caught by the calling program, the execution will abort.

Possible errors encountered during the computation of a Chebyshev model are:

<TABLE border="1">
<CAPTION><EM>Errors during the computation of an ISM</EM></CAPTION>
     <TR><TH><b>Number</b> <TD><b>Description</b>
     <TR><TH><tt>1</tt> <TD>Division by zero scalar
     <TR><TH><tt>2</tt> <TD>Inverse operation with zero in range
     <TR><TH><tt>3</tt> <TD>Log operation with non-positive numbers in range
     <TR><TH><tt>4</tt> <TD>Square-root operation with negative numbers in range
     <TR><TH><tt>5</tt> <TD>Tangent operation with (k+1/2)·PI in range
     <TR><TH><tt>-1</tt> <TD>Variable index is out of range
     <TR><TH><tt>-2</tt> <TD>Operation between variables belonging to different models not permitted
     <TR><TH><tt>-3</tt> <TD>Feature not yet implemented
     <TR><TH><tt>-33</tt> <TD>Internal error
</TABLE>


Further exceptions may be thrown by the template class itself.


\section sec_ISM_refs References

- J. Su, Y. Zha, K. Wang, M.E. Villanueva, R. Paulen, B. Houska. <A HREF="https://doi.org/10.1016/j.ifacol.2019.06.124">Interval superposition arithmetic for guaranteed parameter estimation</A>, <I>IFAC-PapersOnLine</I> <B>52</B>(1):574-579, 2019.
- Y. Zha, M.E. Villanueva, B. Houska. <A HREF="https://arxiv.org/abs/1610.05862">Interval superposition arithmetic</A>, <I>ArXiv</I>, 1610.05862v2, 2018.
.
*/

#ifndef MC__UNIVARODELS_HPP
#define MC__UNIVARODELS_HPP
#include <fstream>
#include <iostream>
#include <iomanip> 
#include <vector>
#include <algorithm>
#include <cmath>
#include <limits>
#include <bitset>
#include <cassert>
//#include <cfloat>

#include "mcop.hpp"
#include "mcfunc.hpp"

//#define MC__UPWLE_DEBUG_TRACE
#undef  MC__UNIVARODELS_DEBUG_PROD

namespace mc
{
//! @brief C++ class for Univariate Piecewise Linear Over- and Under-estimators
////////////////////////////////////////////////////////////////////////
//! mc::UnivarPWLE is a C++ class for propagation of univariate piecewise linear 
//! estimators (UPWLE) through (univariate) factorable functions. The template  
//! parameter corresponds to the type used to propagate the breakpoints of UPWLEs.
////////////////////////////////////////////////////////////////////////

#ifdef MC__ASM_DEBUG_LOGGING
std::ofstream MC__ASM_DEBUG_LOGGER;
bool asm_write_log(const std::string& log_str) {
  MC__ASM_DEBUG_LOGGER.open("MC__ASM_DEBUG_LOGGER.txt", std::ios::app); 
  if (!MC__ASM_DEBUG_LOGGER.is_open()){
    std::cerr << "Open log file failed! in writing " << log_str << std::endl;
    return false;
  }
  MC__ASM_DEBUG_LOGGER << log_str;
  MC__ASM_DEBUG_LOGGER.close();
  return true;
}
#ifdef MC__ASM_DEBUG_EVAL_LOGGING
  std::ofstream MC__ASM_EVAL_LOGGER_N;  
#endif
#endif

template <typename T> 
class UnivarPWLE 
{ //typedef std::pair<std::vector<double>,std::vector<double>> PWL;

 public:
  //template <typename U> friend class ASModel;
  //template <typename U> friend class ASVar;

  template <typename U> friend std::ostream& operator<<
    ( std::ostream &, UnivarPWLE<U> const& );

  template <typename U> friend UnivarPWLE<U> operator+
    ( UnivarPWLE<U> const& );
  template <typename U> friend UnivarPWLE<U> operator+
    ( UnivarPWLE<U> && );

  template <typename U> friend UnivarPWLE<U> operator+
    ( UnivarPWLE<U> const&, UnivarPWLE<U> const& );
  template <typename U> friend UnivarPWLE<U> operator+
    ( UnivarPWLE<U> const&, UnivarPWLE<U> && );
  template <typename U> friend UnivarPWLE<U> operator+
    ( UnivarPWLE<U> &&, UnivarPWLE<U> const& );    
  template <typename U> friend UnivarPWLE<U> operator+
    ( UnivarPWLE<U> &&, UnivarPWLE<U> && );    

  template <typename U> friend UnivarPWLE<U> operator+
    ( double const&, UnivarPWLE<U> const& );
  template <typename U> friend UnivarPWLE<U> operator+
    ( double const&, UnivarPWLE<U> && );

  template <typename U> friend UnivarPWLE<U> operator+
    ( UnivarPWLE<U> const&, double const& );
  template <typename U> friend UnivarPWLE<U> operator+
    ( UnivarPWLE<U> &&, double const& );


  template <typename U> friend UnivarPWLE<U> operator-
    ( UnivarPWLE<U> const& );
  template <typename U> friend UnivarPWLE<U> operator-
    ( UnivarPWLE<U> && );

  template <typename U> friend UnivarPWLE<U> operator-
    ( UnivarPWLE<U> const&, UnivarPWLE<U> const& );
  template <typename U> friend UnivarPWLE<U> operator-
    ( UnivarPWLE<U> const&, UnivarPWLE<U> && );
  template <typename U> friend UnivarPWLE<U> operator-
    ( UnivarPWLE<U> &&, UnivarPWLE<U> const& );
  template <typename U> friend UnivarPWLE<U> operator-
    ( UnivarPWLE<U> &&, UnivarPWLE<U> && );

  template <typename U> friend UnivarPWLE<U> operator-
    ( double const&, UnivarPWLE<U> const& );
  template <typename U> friend UnivarPWLE<U> operator-
    ( double const&, UnivarPWLE<U> && );


  template <typename U> friend UnivarPWLE<U> operator-
    ( UnivarPWLE<U> const&, double const& );
  template <typename U> friend UnivarPWLE<U> operator-
    ( UnivarPWLE<U> &&, double const& );



  template <typename U> friend UnivarPWLE<U> operator*
    ( UnivarPWLE<U> const&, UnivarPWLE<U> const& );
  // template <typename U> friend UnivarPWLE<U> operator*
  //   ( UnivarPWLE<U> const&, UnivarPWLE<U> && );
  // template <typename U> friend UnivarPWLE<U> operator*
  //   ( UnivarPWLE<U> &&, UnivarPWLE<U> const& );
  // template <typename U> friend UnivarPWLE<U> operator*
  //   ( UnivarPWLE<U> &&, UnivarPWLE<U> && );


  template <typename U> friend UnivarPWLE<U> operator*
    ( double const&, UnivarPWLE<U> const& );
  template <typename U> friend UnivarPWLE<U> operator*
    ( double const&, UnivarPWLE<U> && );


  template <typename U> friend UnivarPWLE<U> operator*
    ( UnivarPWLE<U> const&, double const& );
  template <typename U> friend UnivarPWLE<U> operator*
    ( UnivarPWLE<U> &&, double const& );


  template <typename U> friend UnivarPWLE<U> operator/
    ( UnivarPWLE<U> const&, UnivarPWLE<U> const& );
  // template <typename U> friend UnivarPWLE<U> operator/
  //   ( UnivarPWLE<U> const&, UnivarPWLE<U> && );
  // template <typename U> friend UnivarPWLE<U> operator/
  //   ( UnivarPWLE<U> &&, UnivarPWLE<U> const& );
  // template <typename U> friend UnivarPWLE<U> operator/
  //   ( UnivarPWLE<U> &&, UnivarPWLE<U> && );


  template <typename U> friend UnivarPWLE<U> operator/
    ( double const&, UnivarPWLE<U> const& );
  template <typename U> friend UnivarPWLE<U> operator/
    ( double const&, UnivarPWLE<U> && );

  template <typename U> friend UnivarPWLE<U> operator/
    ( UnivarPWLE<U> const&, double const& );
  template <typename U> friend UnivarPWLE<U> operator/
    ( UnivarPWLE<U> &&, double const& );

  // template <typename U> friend UnivarPWLE<U> dotprod
  //   ( std::vector<UnivarPWLE<U>> const&, double const& );

  // template <typename U> friend UnivarPWLE<U> max
  //   ( UnivarPWLE<U> const&, UnivarPWLE<U> const& );
  // template <typename U> friend UnivarPWLE<U> max
  //   ( UnivarPWLE<U> const&, double const& );
  // template <typename U> friend UnivarPWLE<U> max
  //   ( UnivarPWLE<U> &&, double const& );
  // template <typename U> friend UnivarPWLE<U> min
  //   ( UnivarPWLE<U> const&, UnivarPWLE<U> const& );
  // template <typename U> friend UnivarPWLE<U> min
  //   ( UnivarPWLE<U> const&, double const& );
  // template <typename U> friend UnivarPWLE<U> min
  //   ( UnivarPWLE<U> &&, double const& );

  // template <typename U> friend UnivarPWLE<U> inv
  //   ( UnivarPWLE<U> const& );
  // template <typename U> friend UnivarPWLE<U> inv
  //   ( UnivarPWLE<U>&& );
  // template <typename U> friend UnivarPWLE<U> sqr
  //   ( UnivarPWLE<U> const& );
  // template <typename U> friend UnivarPWLE<U> sqr
  //   ( UnivarPWLE<U> && );
  // template <typename U> friend UnivarPWLE<U> sqrt
  //   ( UnivarPWLE<U> const& );
  // template <typename U> friend UnivarPWLE<U> sqrt
  //   ( UnivarPWLE<U> && );
  // template <typename U> friend UnivarPWLE<U> fabs
  //   ( UnivarPWLE<U> const& );
  // template <typename U> friend UnivarPWLE<U> fabs
  //   ( UnivarPWLE<U> && );  
  template <typename U> friend UnivarPWLE<U> relu
    ( UnivarPWLE<U> const& );
  template <typename U> friend UnivarPWLE<U> relu
    ( UnivarPWLE<U> && );    
  // template <typename U> friend UnivarPWLE<U> exp
  //   ( UnivarPWLE<U> const& );
  // template <typename U> friend UnivarPWLE<U> exp
  //   ( UnivarPWLE<U> && );
  // template <typename U> friend UnivarPWLE<U> log
  //   ( UnivarPWLE<U> const& );
  // template <typename U> friend UnivarPWLE<U> log
  //   ( UnivarPWLE<U> && );
  // template <typename U> friend UnivarPWLE<U> xlog
  //   ( UnivarPWLE<U> const& );
  // template <typename U> friend UnivarPWLE<U> xlog
  //   ( UnivarPWLE<U> && );
  // template <typename U> friend UnivarPWLE<U> sin
  //   ( UnivarPWLE<U> const& );
  // template <typename U> friend UnivarPWLE<U> sin
  //   ( UnivarPWLE<U> && );
  // template <typename U> friend UnivarPWLE<U> cos
  //   ( UnivarPWLE<U> const& );
  // template <typename U> friend UnivarPWLE<U> cos
  //   ( UnivarPWLE<U> && );
  // template <typename U> friend UnivarPWLE<U> tan
  //   ( UnivarPWLE<U> const& );
  // template <typename U> friend UnivarPWLE<U> tan
  //   ( UnivarPWLE<U> && );
  // template <typename U> friend UnivarPWLE<U> tanh
  //   ( UnivarPWLE<U> const& );
  // template <typename U> friend UnivarPWLE<U> tanh
  //   ( UnivarPWLE<U> && );
  // template <typename U> friend UnivarPWLE<U> pow
  //   ( UnivarPWLE<U> const&, int const& n );
  // template <typename U> friend UnivarPWLE<U> pow
  //   ( UnivarPWLE<U> &&, int const& n );
  // template <typename U> friend UnivarPWLE<U> pow
  //   ( UnivarPWLE<U> const&, double const&  );
  // template <typename U> friend UnivarPWLE<U> pow
  //   ( UnivarPWLE<U> &&, double const&  );
  // template <typename U> friend UnivarPWLE<U> intersect
  //   ( UnivarPWLE<U> const& , U);
  // template <typename U> friend UnivarPWLE<U> intersect
  //   ( UnivarPWLE<U> && , U);  
  // template <typename U> friend UnivarPWLE<U> affine_transform
  //   ( std::vector<UnivarPWLE<U> > const& , std::vector<double> const&, const double);    
  // template <typename U> friend UnivarPWLE<U> cheb
  //   ( UnivarPWLE<U> const&, unsigned int const& n );
  // template <typename U> friend UnivarPWLE<U> cheb
  //   ( UnivarPWLE<U> &&, unsigned int const& n );

 public:
  //! @brief the upper limit of the number of breakpoints;
  static unsigned int nbpsMax;
  //! @brief var vector, the vector of the coordinates in x-axis for breakpoints 
  std::vector<T> first;
  //! @brief val vector, the vector of the coordinates in y-axis for breakpoints 
  std::vector<T> second;
  //! @brief bounds: (((lb,isComputed),(ub,isComputed)),isBoundComputed) 
  // mutable std::pair<std::pair<std::pair<T,bool>,std::pair<T,bool>>,bool> _bnd;
  mutable std::pair<T,bool> _lbnd;
  mutable std::pair<T,bool> _ubnd;

 private:
  //! @brief constant value if the variable is constant
  T _cst;
  //! @brief flag if the PWL is underestimating the variable
  bool _isUnder;
  //! @brief Internal Intermediate containing the forward difference of the coordinates of breakpoints in x-axis
  mutable std::pair<std::vector<T>,bool> _xFwdDiff;
  //! @brief Internal Intermediate containing the forward difference of the coordinates of breakpoints in y-axis
  mutable std::pair<std::vector<T>,bool> _yFwdDiff;



 public:

  /*
  //! @brief Default Constructor of UnivarPWLE 
  UnivarPWLE
  ()
  {

  }
  */
  //! @brief Constructor of UnivarPWLE for a variable x_i on [a,b]
  UnivarPWLE
  ( T a, T b, bool isUnder)
  : first(2, a), second(2, b),_isUnder(isUnder)
  {
  // : first(3, a), second(3, a), _bnd(std::make_pair(std::make_pair( std::min(a,b), true ),std::make_pair( std::max(a,b), true )),true),_cst(0.),_isUnder(isUnder),
  // _xFwdDiff(0,false),_yFwdDiff(0,false)
    // first[0] = 0.5*(a + b); 
    // first[1] = std::min(a,b) - first[0];
    // first[2] = std::max(a,b) - first[0];
    // second[0] = first[0];
    // second[1] = first[1];
    // second[2] = first[2];

    first[1]  = b;
    second[0] = a;
#ifdef MC__UPWLE_DEBUG
    if(b - a <= T(0.)) std::cout << "b " << b << " a " << a << std::endl;   
    assert(b - a > T(0.));
#endif
  }

  //! @brief Constructor of UnivarPWLE for a unknown variable x_i on [a,b] with the constant bound c
  UnivarPWLE
   ( T a, T b, T cst, bool isUnder)
  : first(2, a), second(1, cst),_isUnder(isUnder)
  {
  // : first(3, a), second(1, cst), _bnd(std::make_pair(std::make_pair( cst, true ),std::make_pair( cst, true )),true),_cst(T(0.)),_isUnder(isUnder),
  // _xFwdDiff(0,false),_yFwdDiff(0,false)    
    // first[0] = 0.5*(a + b); 
    // first[1] = std::min(a,b) - first[0];
    // first[2] = std::max(a,b) - first[0];
     first[1] = b;
#ifdef MC__UPWLE_DEBUG    
    if(b - a <= T(0.)) std::cout << "b " << b << " a " << a << std::endl;   
    assert(b - a > T(0.));
#endif  
  }

  // //! @brief Constructor of UnivarPWLE for a constant c
  UnivarPWLE
  ( T cst = T(0.))
  : first(0), _cst(cst)
  {}

  // : first(0), second(0), _bnd(std::make_pair(std::make_pair( cst, true ),std::make_pair( cst, true )),true),_cst(cst),
  //   _isUnder(true),_xFwdDiff(0,false),_yFwdDiff(0,false)

  // UnivarPWLE
  // (const std::vector<T> & firstIn, bool isUnder)
  // : first(firstIn), second(1,T(0.)),_isUnder(isUnder)
  // {
  // // : first(firstIn), second(1,T(0.)), _lbnd(T(0.), true ),_ubnd(T(0.), true ),
  // //   _isUnder(isUnder),_xFwdDiff(0,false),_yFwdDiff(0,false)    
  //   // if(firstIn.empty){
  //   //   _cst = 0.;
  //   // }
  // }

  // UnivarPWLE
  // (std::vector<T> && firstIn, bool isUnder)
  // : first(std::move(firstIn)), second(1,T(0.)),_isUnder(isUnder)
  // {
  // // : first(std::move(firstIn)), second(1,T(0.)), _lbnd(T(0.), true ),_ubnd(T(0.), true ),
  // //   _isUnder(isUnder),_xFwdDiff(0,false),_yFwdDiff(0,false)    
  //   // if(firstIn.empty){
  //   //   _cst = 0.;
  //   // }    
  // }




  //! @brief Copy constructor of UnivarPWLE for a var (lvalue)
  UnivarPWLE
  ( UnivarPWLE<T> const& var )
  {
#ifdef MC__ASModel_TRACE
    std::cerr << "-- UnivarPWLE( UnivarPWLE<T> const& )\n";
#endif
#ifdef TEST_MOVE
    std::cout << "Copy Constructor" << std::endl;
#endif
    if( this == &var ) return;
    first  = var.first;
    if(first.empty()){
      _cst = var._cst;
      return;
    }    
    second = var.second;
    _isUnder = var._isUnder;
    if(second.size() > 2){
      _lbnd = var._lbnd;
      _ubnd = var._ubnd;
      _xFwdDiff = std::make_pair(std::vector<T>(0),false);
      _yFwdDiff = std::make_pair(std::vector<T>(0),false);   
    //if(var._xFwdDiff.second) _xFwdDiff = var._xFwdDiff;
    //if(var._yFwdDiff.second) _yFwdDiff = var._yFwdDiff;        
    }
    



  }


  //! @brief Copy constructor of UnivarPWLE for a var (rvalue)
  UnivarPWLE
  ( UnivarPWLE<T> && var )
  {
#ifdef MC__ASModel_TRACE
    std::cerr << "-- UnivarPWLE( UnivarPWLE<T> && var )\n";
#endif
#ifdef TEST_MOVE
    std::cout << "Move Constructor" << std::endl;
#endif
    if( this == &var ) return;
    first  = std::move(var.first);    
    if(first.empty()){
      _cst = var._cst;
      return;
    }        
    second = std::move(var.second);
    _isUnder = var._isUnder;
    if(second.size() > 2){
      _lbnd = std::move(var._lbnd);
      _ubnd = std::move(var._ubnd);
      //_xFwdDiff = std::make_pair(std::vector<T>(0),false);
      //_yFwdDiff = std::make_pair(std::vector<T>(0),false);   
      _xFwdDiff = var._xFwdDiff.second?std::move(var._xFwdDiff):std::make_pair(std::vector<T>(0),false);
      _yFwdDiff = var._yFwdDiff.second?std::move(var._yFwdDiff):std::make_pair(std::vector<T>(0),false);        
    }

  }
  

//   //! @brief Copy constructor (with scaling) of UnivarPWLE for a var (lvalue)
//   UnivarPWLE
//   ( UnivarPWLE<T> const& var, const double mtpr ) // multiplier
//   : _cst(var._cst),_isUnder(var._isUnder)
//   {
// #ifdef MC__ASModel_TRACE
//     std::cerr << "-- ASVar( ASVar<T> const&, const double )\n";
// #endif
// #ifdef TEST_MOVE
//     std::cout << "Copy and Scale Constructor" << std::endl;
// #endif
//     if( var.first.empty() ) {_cst *= mtpr; return;}

//     if( mtpr == 0. ) {_cst = 0.; 
//     first.resize(0); second.resize(0);
//     _xFwdDiff = std::make_pair(std::vector<T>(0),false);
//     _yFwdDiff = std::make_pair(std::vector<T>(0),false);
//     _bnd = std::make_pair(std::make_pair(std::make_pair( 0., false ),std::make_pair( 0., false )),false); 
//     return;
//     }

//             // if( !_ndep )
//             //   throw typename ASModel<T>::Exceptions(  ASModel<T>::Exceptions::DEPINDEX  );
    
//     first = var.first;
//     const unsigned int nbps = var.second.size();
//     second.resize(nbps); 

//     for (unsigned int i = 0; i < nbps; i++){
//         second[i] = var.second[i] * mtpr;
//     }

//     if(var._bnd.second){
//       _bnd = var._bnd;
//       if (mtpr > 0){
//         //if(_bnd.first.first.second)  
//         _bnd.first.first.first *= mtpr;
//         //if(_bnd.first.second.second) 
//         _bnd.first.second.first *= mtpr;        
//       }
//       else{
//         //if(_bnd.first.first.second)  
//         _bnd.first.first.first *= mtpr;
//         //if(_bnd.first.second.second) 
//         _bnd.first.second.first *= mtpr;
//         std::swap(_bnd.first.first,_bnd.first.second);   
//       }      
//     }
//     else 
//       _bnd = std::make_pair(std::make_pair(std::make_pair( 0., false ),std::make_pair( 0., false )),false); 

//     if(var._yFwdDiff.second){ 
//       _yFwdDiff = std::make_pair(std::vector<T>(nbps - 2),true);  
//       for (unsigned int i = 0; i < nbps - 2; i++){
//           _yFwdDiff.first[i] = mtpr * var._yFwdDiff.first[i];
//       }          
//     }
//     else{
//       _yFwdDiff = std::make_pair(std::vector<T>(0),false);  
//     }      

//     _xFwdDiff = var._xFwdDiff.second? var._xFwdDiff:std::make_pair(std::vector<T>(0),false);  

//   }


#ifdef ASM_LIFITIME_DEBUG     
  //! @brief Destructor of ASM (for checking the lifetime of ASMs)
  ~UnivarPWLE() 
  {
    std::cout<< "UnivarPWLE delated" <<std::endl;
  }
#endif    

  //! @brief Exceptions of mc::UnivarPWLE
  class Exceptions
  {
   public:
    //! @brief Enumeration type for UnivarPWLE exception handling
    enum TYPE{
      INDEX=1,	   //!< breakpoint index out of range
      ESTMATCH,    //!< Addition is defined for two under- or over-estimators
      UNDEF,	     //!< Feature not yet implemented
      DIV          //!< Division by zero scalar      
    };
    //! @brief Constructor for error <a>ierr</a>
    Exceptions( TYPE ierr ) : _ierr( ierr ){}
    //! @brief Error flag
    int ierr(){ return _ierr; }
    //! @brief Error description
    std::string what(){
      switch( _ierr ){
      case INDEX:
        return "mc::UnivarPWLE\t breakpoint index out of range";
      case ESTMATCH:
        return "mc::UnivarPWLE\t Addition cannot be performed between an over- and an under-estimator";
      case UNDEF:
        return "mc::UnivarPWLE\t Feature not yet implemented";
      case DIV:
        return "mc::UnivarPWLE\t Divided by zero";        
      default:
        return "mc::UnivarPWLE\t Internal error";
      }
    }
   private:
    TYPE _ierr;
  };


  UnivarPWLE<T>& operator=
  ( T const& cst )
  {
    
    if(first.empty()){
      _cst = cst;
    }
    else{
      _cst = cst;
      first.resize(0);
      second.resize(0);
      //_bnd = std::make_pair(std::make_pair(std::make_pair( cst, true ),std::make_pair( cst, true )),true);
      _xFwdDiff = std::make_pair(std::vector<T>(0),false);
      _yFwdDiff = std::make_pair(std::vector<T>(0),false);    

      // _bnd.first.first.first = cst;
      // _bnd.first.second.first = cst;    
      // _bnd.first.first.second = true;
      // _bnd.first.second.second = true;
      // _bnd.second = true;
  
      // _yFwdDiff.second = false;
      // _xFwdDiff.second = false;
        
      // first[2] = first.back(); 
      // first.resize(3);
      // second[0] = cst;
      // second.resize(1); 
    }    
        
    return *this;
  }


  UnivarPWLE<T>& operator=
  ( UnivarPWLE<T> const& var )
  {
#ifdef MC__ASModel_TRACE
    std::cerr << "-- UnivarPWLE<T>& operator= ( UnivarPWLE<T> const& )\n";
#endif
#ifdef TEST_MOVE
    std::cout << "Copy Assign" << std::endl;
#endif
    if( this == &var )
      return *this;
    first = var.first;
    if(first.empty()){
      _cst = var._cst;
      return *this;
    }
    second = var.second;
    _isUnder = var._isUnder;
    if (second.size() > 2){
      _lbnd = var._lbnd;
      _ubnd = var._ubnd;

      _xFwdDiff = std::make_pair(std::vector<T>(0),false);
      _yFwdDiff = std::make_pair(std::vector<T>(0),false);

      // if(var._xFwdDiff.second) _xFwdDiff = var._xFwdDiff;
      // if(var._yFwdDiff.second) _yFwdDiff = var._yFwdDiff; 
        
    }    

    // first  = var.first;
    // second = var.second;   
    // _cst = var._cst;
    // _bnd = var._bnd;
    // _isUnder = var._isUnder;  

    // _xFwdDiff = std::make_pair(std::vector<T>(0),false);
    // _yFwdDiff = std::make_pair(std::vector<T>(0),false);
    //_xFwdDiff = var._xFwdDiff.second? var._xFwdDiff:std::make_pair(std::vector<T>(0),false);  
    //_yFwdDiff = var._yFwdDiff.second? var._yFwdDiff:std::make_pair(std::vector<T>(0),false); 

    return *this;
  } 

  UnivarPWLE<T>& operator=
  ( UnivarPWLE<T> && var )
  {
#ifdef MC__ASModel_TRACE
    std::cerr << "-- UnivarPWLE<T>& operator= ( UnivarPWLE<T> && )\n";
#endif
#ifdef TEST_MOVE
    std::cout << "Move Assign" << std::endl;
#endif
    if( this == &var )
      return *this;
    first  = std::move(var.first);
    if(first.empty()){
      _cst = var._cst;
      return *this;
    }
    second = std::move(var.second);  
    _isUnder = var._isUnder;
    if (second.size() > 2){
      _lbnd = std::move(var._lbnd);
      _ubnd = std::move(var._ubnd);

      _xFwdDiff = var._xFwdDiff.second? std::move(var._xFwdDiff):std::make_pair(std::vector<T>(0),false);  
      _yFwdDiff = var._yFwdDiff.second? std::move(var._yFwdDiff):std::make_pair(std::vector<T>(0),false); 

      // if(var._xFwdDiff.second) _xFwdDiff = var._xFwdDiff;
      // if(var._yFwdDiff.second) _yFwdDiff = var._yFwdDiff; 
      
    }    

    return *this;
  }


  UnivarPWLE<T>& operator+=
    ( UnivarPWLE<T> const& );
  UnivarPWLE<T>& operator+=
    ( UnivarPWLE<T> && );       
  UnivarPWLE<T>& operator+=
    ( double const& );
   UnivarPWLE<T>& operator+=
    ( std::pair<T,T> const&  );     
  // UnivarPWLE<T>& operator+=
  //   ( T const& );


  UnivarPWLE<T>& operator-=
    ( UnivarPWLE<T> const& );
  UnivarPWLE<T>& operator-=
    ( double const& );
  // UnivarPWLE<T>& operator-=
  //   ( T const& );
 
  UnivarPWLE<T>& operator*=
    ( UnivarPWLE<T> const& );
  UnivarPWLE<T>& operator*=
    ( double const& );
  // UnivarPWLE<T>& operator*=
  //   ( T const& );

  UnivarPWLE<T>& operator/=
    ( UnivarPWLE<T> const& );
  UnivarPWLE<T>& operator/=
    ( double const& );


  bool empty () const {
    return first.empty();
  } 

  bool isCst () const {
    if (first.empty()) return true;
    else if (second.size() <= 1) return true; 
    else return false;
  } 
  
  bool isLinear () const {
    if (first.empty()) return false;
    else if (second.size() == 2) return true; 
    else return false;
  }   

  std::pair<T,bool> get_cst () const {
    if (first.empty()) return std::make_pair( _cst, true );
    else if (second.size() <= 1) return std::make_pair( second.front(), true ); 
    else return std::make_pair( 0., false );
  } 

  T lbVar () const {
    if (first.empty()) {std::cout << "lb of var is unknown in the constant function"<<std::endl; return 0.;}
    else if(first.size() == 2){
      return first[0];
    }
    else 
      return first[0] + first[1];
  }   

  T ubVar () const {
    if (first.empty()) {std::cout << "ub of var is unknown in the constant function"<<std::endl; return 0.;}
    else if(first.size() == 2){
      return first[1];
    }
    else 
      return first[0] + first.back();
  }   

  std::ostream& display
  ( std::ostream& out=std::cout )
  const
  {
    
    if(first.empty()){
      std::cout << " cannot display a constant function without knowing the range of the var" << std::endl;
    }
    else{
      const unsigned int sizeTmp = second.size();
      if (sizeTmp == 1){
        if (_isUnder){
          out << std::setw(14) << first[1] + first[0]  << std::setw(14) << second[0]  << std::endl;
          out << std::setw(14) << first[2] + first[0]  << std::setw(14) << second[0]  << std::endl;
        }
        else{
          out << std::setw(14) << first[2] + first[0]  << std::setw(14) << second[0]  << std::endl;
          out << std::setw(14) << first[1] + first[0]  << std::setw(14) << second[0]  << std::endl;        
        }  
      }
      else if(isLinear()){
        if (_isUnder){
          out << std::setw(14) << first[0]  << std::setw(14) << second[0]  << std::endl;
          out << std::setw(14) << first[1]  << std::setw(14) << second[1]  << std::endl;
        }
        else{
          out << std::setw(14) << first[1]   << std::setw(14) << second[1]  << std::endl;
          out << std::setw(14) << first[0]   << std::setw(14) << second[0]  << std::endl;        
        }          
      }
      else{
        if (_isUnder){
          for( unsigned int j=1; j < sizeTmp; j++){
            out << std::setw(14) << first[j] + first[0]  << std::setw(14) << second[j] + second[0]  << std::endl;
          }     
        }
        else{
          for( unsigned int j=1; j < sizeTmp; j++){
            out << std::setw(14) << first[sizeTmp-j] + first[0]  << std::setw(14) << second[sizeTmp-j] + second[0]  << std::endl;
          }          
        }     
      }
    }
    return out;
  }


  std::ostream& display
  ( std::ostream& out, UnivarPWLE<T> const & BEstmrIn )
  const
  {
    UnivarPWLE<T> BEstmr(BEstmrIn);
    std::vector<T> Afirst(first);
    std::vector<T> Asecond(second);    
    if(first.empty() || BEstmr.empty()){
      std::cout << " cannot display a constant function without knowing the range of one of the two vars" << std::endl;
      return out;
    }
    if(isLinear() && BEstmr.isLinear()){
      out << std::setw(14) <<  Afirst[0] << std::setw(14) << BEstmr.first[0]
          << std::setw(14) << Asecond[0] + BEstmr.second[0] << std::endl; 
      out << std::setw(14) << Afirst[1] << std::setw(14) << BEstmr.first[0] 
          << std::setw(14) << Asecond[1] + BEstmr.second[0] << std::endl;                 
      out << std::setw(14) << Afirst[1]  << std::setw(14) << BEstmr.first[1]  
          << std::setw(14) << Asecond[1] + BEstmr.second[1] << std::endl; 
      out << std::setw(14) <<  Afirst[0]  << std::setw(14) << BEstmr.first[1]
          << std::setw(14) << Asecond[0] + BEstmr.second[1] << std::endl; 
      out << std::setw(14) <<  Afirst[0]   << std::setw(14) << BEstmr.first[0]
          << std::setw(14) << Asecond[0] + BEstmr.second[0] << std::endl;
      out << std::endl << std::endl;         
      return out;
      
    }
    else if(!isLinear() && BEstmr.isLinear()){
      BEstmr.first.push_back(BEstmr.first[0]+BEstmr.first[1]);
      BEstmr.first[2] *= 0.5;
      BEstmr.first[0] -= BEstmr.first[2];
      BEstmr.first[1] -= BEstmr.first[2];
      std::swap(BEstmr.first[2],BEstmr.first[0]);
      std::swap(BEstmr.first[1],BEstmr.first[2]);

      BEstmr.second.push_back(BEstmr.second[0]+BEstmr.second[1]);
      BEstmr.second[2] *= 0.5;
      BEstmr.second[0] -= BEstmr.second[2];
      BEstmr.second[1] -= BEstmr.second[2];
      std::swap(BEstmr.second[2],BEstmr.second[0]);
      std::swap(BEstmr.second[1],BEstmr.second[2]);                
    }
    else if(isLinear() && !BEstmr.isLinear()){
      Afirst.push_back(Afirst[0]+Afirst[1]);
      Afirst[2] *= 0.5;
      Afirst[0] -= Afirst[2];
      Afirst[1] -= Afirst[2];
      std::swap(Afirst[2],Afirst[0]);
      std::swap(Afirst[1],Afirst[2]);
      
      Asecond.push_back(Asecond[0]+Asecond[1]);
      Asecond[2] *= 0.5;
      Asecond[0] -= Asecond[2];
      Asecond[1] -= Asecond[2];
      std::swap(Asecond[2],Asecond[0]);
      std::swap(Asecond[1],Asecond[2]);  
    }

    if(first.size() == 2){
      Afirst.push_back(Afirst[0]+Afirst[1]);
      Afirst[2] *= 0.5;
      Afirst[0] -= Afirst[2];
      Afirst[1] -= Afirst[2];
      std::swap(Afirst[2],Afirst[0]);
      std::swap(Afirst[1],Afirst[2]);      
    }   

    if(BEstmr.first.size() == 2){
      BEstmr.first.push_back(BEstmr.first[0]+BEstmr.first[1]);
      BEstmr.first[2] *= 0.5;
      BEstmr.first[0] -= BEstmr.first[2];
      BEstmr.first[1] -= BEstmr.first[2];
      std::swap(BEstmr.first[2],BEstmr.first[0]);
      std::swap(BEstmr.first[1],BEstmr.first[2]);      
    }   

    const unsigned int AsizeTmp = second.size();
    const unsigned int BsizeTmp = BEstmr.second.size();
    // if(AsizeTmp == 0 || BsizeTmp == 0){
    //   std::cout << " cannot display a constant function without knowing the range of one of the two vars" << std::endl;
    // }
    if(AsizeTmp == 1 && BsizeTmp == 1){

      out << std::setw(14) << Afirst[1] + Afirst[0] << std::setw(14) << BEstmr.first[1] + BEstmr.first[0]
          << std::setw(14) << Asecond[0] + BEstmr.second[0] << std::endl; 
      out << std::setw(14) << Afirst[2] + Afirst[0] << std::setw(14) << BEstmr.first[1] + BEstmr.first[0]
          << std::setw(14) << Asecond[0] + BEstmr.second[0] << std::endl;                 
      out << std::setw(14) << Afirst[2] + Afirst[0] << std::setw(14) << BEstmr.first[2] + BEstmr.first[0]
          << std::setw(14) << Asecond[0] + BEstmr.second[0] << std::endl; 
      out << std::setw(14) << Afirst[1] + Afirst[0] << std::setw(14) << BEstmr.first[2] + BEstmr.first[0]
          << std::setw(14) << Asecond[0] + BEstmr.second[0] << std::endl; 
      out << std::setw(14) << Afirst[1] + Afirst[0] << std::setw(14) << BEstmr.first[1] + BEstmr.first[0]
          << std::setw(14) << Asecond[0] + BEstmr.second[0] << std::endl;
      out << std::endl << std::endl;                                 
    }
    else if(AsizeTmp == 1){
      for( unsigned int j2=2; j2< BsizeTmp; j2++){  
        out << std::setw(14) << Afirst[1] + Afirst[0] << std::setw(14) << BEstmr.first[j2-1] + BEstmr.first[0]
            << std::setw(14) << Asecond[0] + BEstmr.second[0] + BEstmr.second[j2-1] << std::endl; 
        out << std::setw(14) << Afirst[2] + Afirst[0] << std::setw(14) << BEstmr.first[j2-1] + BEstmr.first[0]
            << std::setw(14) << Asecond[0] + BEstmr.second[0] + BEstmr.second[j2-1]  << std::endl;                 
        out << std::setw(14) << Afirst[2] + Afirst[0] << std::setw(14) << BEstmr.first[j2] + BEstmr.first[0]
            << std::setw(14) << Asecond[0] + BEstmr.second[0] + BEstmr.second[j2]  << std::endl; 
        out << std::setw(14) << Afirst[1] + Afirst[0] << std::setw(14) << BEstmr.first[j2] + BEstmr.first[0]
            << std::setw(14) << Asecond[0] + BEstmr.second[0] + BEstmr.second[j2]  << std::endl; 
        out << std::setw(14) << Afirst[1] + Afirst[0] << std::setw(14) << BEstmr.first[j2-1] + BEstmr.first[0]
            << std::setw(14) << Asecond[0] + BEstmr.second[0] + BEstmr.second[j2-1]  << std::endl;
        out << std::endl << std::endl;      
      }

    }
    else if(BsizeTmp == 1){
      for( unsigned int j1=2; j1< AsizeTmp; j1++){  
        out << std::setw(14) << Afirst[j1-1] + Afirst[0] << std::setw(14) << BEstmr.first[1] + BEstmr.first[0]
            << std::setw(14) << Asecond[0] + BEstmr.second[0] + Asecond[j1-1] << std::endl; 
        out << std::setw(14) << Afirst[j1] + Afirst[0] << std::setw(14) << BEstmr.first[1] + BEstmr.first[0]
            << std::setw(14) << Asecond[0] + BEstmr.second[0] + Asecond[j1]  << std::endl;                 
        out << std::setw(14) << Afirst[j1] + Afirst[0] << std::setw(14) << BEstmr.first[2] + BEstmr.first[0]
            << std::setw(14) << Asecond[0] + BEstmr.second[0] + Asecond[j1]  << std::endl; 
        out << std::setw(14) << Afirst[j1-1] + Afirst[0] << std::setw(14) << BEstmr.first[2] + BEstmr.first[0]
            << std::setw(14) << Asecond[0] + BEstmr.second[0] + Asecond[j1-1]  << std::endl; 
        out << std::setw(14) << Afirst[j1-1] + Afirst[0] << std::setw(14) << BEstmr.first[1] + BEstmr.first[0]
            << std::setw(14) << Asecond[0] + BEstmr.second[0] + Asecond[j1-1]  << std::endl;  
        out << std::endl << std::endl;    
      }
    }  
    else{
      for( unsigned int j1=2; j1< AsizeTmp; j1++){
        for( unsigned int j2=2; j2< BsizeTmp; j2++){        
          double l1 = Afirst[j1-1] + Afirst[0];
          double u1 = Afirst[j1]   + Afirst[0];          
          double l2 = BEstmr.first[j2-1] + BEstmr.first[0];
          double u2 = BEstmr.first[j2]   + BEstmr.first[0];
          double u11 = Asecond[0] + Asecond[j1-1] + BEstmr.second[0] + BEstmr.second[j2-1];
          double u21 = Asecond[0] + Asecond[j1] + BEstmr.second[0] + BEstmr.second[j2-1];
          double u12 = Asecond[0] + Asecond[j1-1] + BEstmr.second[0] + BEstmr.second[j2];
          double u22 = Asecond[0] + Asecond[j1] + BEstmr.second[0] + BEstmr.second[j2];                    
          out << std::setw(14) << l1 << std::setw(14) << l2 << std::setw(14) << u11 << std::endl;
          out << std::setw(14) << u1 << std::setw(14) << l2 << std::setw(14) << u21 << std::endl;
          out << std::setw(14) << u1 << std::setw(14) << u2 << std::setw(14) << u22 << std::endl;
          out << std::setw(14) << l1 << std::setw(14) << u2 << std::setw(14) << u12 << std::endl;
          out << std::setw(14) << l1 << std::setw(14) << l2 << std::setw(14) << u11 << std::endl;
          out << std::endl << std::endl; 
        }
      }
    }
    return out;
  }



  void flip (){
    if(!first.empty())
    _isUnder = !_isUnder;
  }


  void set_zero (){

    if(first.empty()){
      _cst = T(0.);
    }
    else{
      second[0] = T(0.);
      second.resize(1);       
      if(first.size() > 2){
        first[1] += first[0];  
        first[0] += first.back();
        std::swap(first[1],first[0]); 
        first.resize(2);
        //_yFwdDiff.second = false;
        _yFwdDiff.first.resize(0);
        //_xFwdDiff.second = false;
        _xFwdDiff.first.resize(0);        
      }    
    }
  }

  void condense_by_relax (const unsigned int k){
    _condense_by_relax(k);
  }

  void recal_B()const // for debug, note that we do not check the size of first and second
  {
    // _bnd.first.first.second = false;
    // _bnd.first.second.second = false;
    // _bnd.second = false;
    _lbnd.second = false;
    _ubnd.second = false;    
    _yFwdDiff.second = false;
    _xFwdDiff.second = false;
    _B(true,true);
    _eval_fwdDiff(); 
  }


  void B(bool isUnder, bool isOver) const
  {
    // note that this function should be called when second.size > 2
    _B(isUnder,isOver);
  }

  T get_ub()const
  {
    if(first.empty())
      return _cst;
    switch (second.size()){
      case 1:
      {
        return second[0];
      }
      case 2:
      {
        return std::max(second[0],second[1]);
      }
      //default:
    }
    if(_ubnd.second)
      return _ubnd.first;
    else{
      _B(false,true); return _ubnd.first;
    }
  }

  T get_lb()const
  {
    if(first.empty())
      return _cst;
    switch (second.size()){
      case 1:
      {
        return second[0];
      }
      case 2:
      {
        return std::min(second[0],second[1]);
      }
      //default:
    }
    if(_lbnd.second)
      return _lbnd.first;
    else{
      _B(true,false); return _lbnd.first;
    }
  }


  bool get_flag()const
  {
    return _isUnder;
  }


#ifdef MC__ASM_DEBUG_EVAL_LOGGING
  void write2log
  ( std::ostream& out)
  const 
  {
  
    out << std::scientific << std::setprecision(14);
    if(first.empty()){
      // out << std::setw(24) << _cst << std::endl; 
      // out << std::setw(24) << _cst << std::endl; 
      out << _cst << std::endl; 
      out << _cst << std::endl;           
      return ;
    }
    if(second.size() == 1){
      // out << std::setw(24) << first[0]  << std::setw(24) << first[1] << std::endl;   
      // out << std::setw(24) << second[0] << std::endl; 
      out << first[0]  << "," << first[1] << std::endl;   
      out << second[0] << std::endl; 
      return ;
    }
    if(second.size() == 2){
      // out << std::setw(24) << first[0]  << std::setw(24) << first[1] << std::endl;   
      // out << std::setw(24) << second[0] << std::setw(24) << second[1]<< std::endl;   
      out << first[0]  << "," << first[1] << std::endl;      
      out << second[0]  << "," << second[1] << std::endl;      
      return ;
    }  

    auto&& x = first;
    auto&& y = second;
  
  
    out << std::right;
    for( unsigned int i=1; i< first.size()-1; i++ ){
      // out << std::setw(24) << x[0] + x[i];
      out << x[0] + x[i] << ",";
    }
    out << x[0] + x.back() <<std::endl;

    out << std::right;
    for( unsigned int i=1; i< second.size()-1; i++ ){
      // out << std::setw(24) << y[0] + y[i];
      out << y[0] + y[i] << ",";
    }
    out << y[0] + y.back() <<std::endl; 
  
  // //  var._B(true,true);
  //   out << std::endl;
  //   out << std::right << std::setw(5) << "[LB, UB]: " << "[" <<
  //   var.get_lb() << ", " << var.get_ub() << "]" << std::endl;    
    return ;
  }
#endif

  void check_fwdDiff(){ // this function is for debug.
    if(!isLinear()){
      for(unsigned int i = 0; i < _xFwdDiff.first.size(); i++){
      if(std::fabs(first[i+2] - first[i+1] - _xFwdDiff.first[i]) > 1e-12)
        std::cout << i << " valx " << first[i+2] - first[i+1] << " valx " << _xFwdDiff.first[i] << std::endl;
      if(std::fabs(second[i+2] - second[i+1] - _yFwdDiff.first[i]) > 1e-12)
        std::cout << i << " valy " << second[i+2] - second[i+1] << " valy " << _yFwdDiff.first[i] << std::endl;        
    }
  }
  else{
      std::cout << "this PWLE is linear" << std::endl;
    }
  }


  T eval (T x) const {
    if(first.empty()){
      //std::cout << "x is out of the range" << std::endl;
      return _cst;
    }
    if(second.size()==1) return second[0];
    if(second.size()==2){
      if(x + 4e-15 < first[0] || x - 4e-15 > first[1]){
        std::cout << "x is out of the range" << std::endl;
        return T(first[0]);
      }
      return (x - first[0])*((second[1] - second[0]) / (first[1] - first[0]) ) + second[0];
    }
    if(x + 4e-15 < first[0] + first[1] || (x - 4e-15 > first[0] + first[first.size()-1])){
      std::cout << "x is out of the range"  << std::endl;
      std::cout << "x - (first[0] + first[1]) "  << std::scientific << std::setprecision(14) << x - (first[0] + first[1]) << std::endl;
      std::cout << "x - (first[0] + first[first.size()-1]) "  << std::scientific << std::setprecision(14) << x - (first[0] + first[first.size()-1]) << std::endl;      
      return T(first[0] + first[1]);
    }

    const T xOffset = x - first[0];
    unsigned int lbInd = 1;
    unsigned int ubInd = first.size() - 1;
    unsigned int mid = (lbInd + ubInd)/2;
    while(ubInd - lbInd != 1){
      
      if(xOffset < first[mid])
        ubInd = mid;
      else if(xOffset > first[mid])
        lbInd = mid;
      else 
        return second[0] + second[mid];
      mid = (lbInd + ubInd)/2;
    }
    //std::cout << ubInd << lbInd << std::endl;
    T y = second[lbInd];
    T step = first[ubInd] - first[lbInd];
    T ptVarDiff = first[ubInd] - xOffset;     
    _iterative_interpolate(ptVarDiff,step,y,second[ubInd]);
    y += second[0];    
    return y;
  }

  // void eval (T lb, T ub, T stepsize, std::vector<T> & vals){
  //   ;
  // }

  //void eval (T stepsize, std::vector<T> & vals){
  void eval (const unsigned int N, std::vector<T> & vals) const {  
    if(first.empty()){
      vals.resize(N,_cst);
      return;
    } 
    if(second.size()==1){
      vals.resize(N,second[0]);
      return;
    }
    if(second.size()==2){
      vals.resize(N);
      const T h = (second[1] - second[0])/(T) N;
      T lb = second[0];
      for (unsigned int i = 0; i < N ; i++){
        vals[i] = lb;
        lb += h; 
      }
      return;
    }    
    unsigned int lbInd = 1;
    unsigned int ubInd = first.size() - 1;    
    //const int N = std::ceil((first[ubInd]-first[lbInd])/stepsize);
    vals.resize(N);
    const T stepsize = (first[ubInd]-first[lbInd])/N;
    T element = first[1];
    vals[0] = second[0] + second[1];
    lbInd = 2;
    T ptVarDiff = first[lbInd] - first[1];
    T      step = first[lbInd] - first[1];
    for(unsigned int i = 1; i < N; i++){
      element += stepsize;
      if(element < first[lbInd])
        ptVarDiff -= stepsize;        
      else{
        while(element > first[lbInd]){
          lbInd ++;
        }
        step = first[lbInd] - first[lbInd-1];
        ptVarDiff = first[lbInd] - element; 
      }
      vals[i] = second[lbInd-1];
      _iterative_interpolate(ptVarDiff,step,vals[i],second[lbInd]);
      vals[i] += second[0];
    }
      
  }  


// (const T & ptVarDiff, const T & step, T & PtValLast, const T & PtVal)
// {
//   const T ratio  = ptVarDiff / step;   
//   const T offset = ratio * (PtVal - PtValLast);
//   PtValLast = PtVal - offset;   
// }

  // std::vector<std::vector<T>> align_then_relax_to_PWC
  // ( const std::vector<T> & xIn, const std::vector<T> & yIn)
  // const 
  // {
  //   std::vector<std::vector<T>> out(3);
  //   out[0].resize(10);
  //   out[1].resize(10);
  //   out[2].resize(10);


  //   return out;      
  // }

  std::vector<T> get_PWC
  ( unsigned int N)
  const
  {
    if (first.empty()){
      std::vector<T> output(N,_cst);
      return  output;
    }  
  
    if (second.size() == 1){
      std::vector<T> output(N,second[0]);
      return  output;
    }

    std::vector<T> output(N);
    
    if (second.size() == 2){
      const T h = (second[1] - second[0])/ (T) N;
      if(_isUnder){
        if (h > T(0)){
          T lb = second[0];
          for (unsigned int i = 0; i < N ; i++){
            output[i] = lb;
            lb += h; 
          }
        }
        else{
          T lb = second[0] + h;
          for (unsigned int i = 0; i < N-1; i++){
            output[i] = lb;
            lb += h; 
          }
          output[N-1] = second[1];
        }
      }
      else{
        if (h < T(0)){
          T lb = second[0];
          for (unsigned int i = 0; i < N ; i++){
            output[i] = lb;
            lb += h; 
          }
        }
        else{
          T lb = second[0] + h;
          for (unsigned int i = 0; i < N-1; i++){
            output[i] = lb;
            lb += h; 
          }
          output[N-1] = second[1];
        }        
      }

      return output;    
    }    

    const long double h = (first.back() - first[1])/N;

    long double right = first[1] + h; 
    unsigned int cnt = 1;    
    T minOrMax = second[cnt];
    T ptValLast = second[cnt];
    for (unsigned int i = 0; i < N-1; i++){
      //std::cout << "in loop" << std::endl; 
      while ( first[cnt] < right ){
        minOrMax = _isUnder? std::min(second[cnt],minOrMax):std::max(second[cnt],minOrMax); 
        cnt ++;
      }
      // std::cout << "out while" << std::endl; 
      if(first[cnt] == right){
        ptValLast = second[cnt];
        // std::cout << "cnt " << cnt << " at " << i << std::endl; 
      }
      else{
        ptValLast = second[cnt-1];
        _iterative_interpolate(first[cnt] - right, first[cnt] - first[cnt-1], ptValLast, second[cnt]);
        // std::cout << "cnt " << cnt << " at " << i << std::endl;         
      }  
      minOrMax = _isUnder? std::min(ptValLast,minOrMax):std::max(ptValLast,minOrMax);    
      // std::cout << "out" << std::endl; 
      output[i] = second[0] + minOrMax;
      // std::cout << "output[i] " << output[i] << std::endl;
      right += h;
      minOrMax = ptValLast;
    }
    for (unsigned int i = cnt; i < first.size(); i++){
        minOrMax = _isUnder? std::min(second[i],minOrMax):std::max(second[i],minOrMax); 
    }   
    // std::cout << "out" << std::endl; 
    output[N-1] = second[0] + minOrMax;    
    return output;

  }




 private:

  //! @brief Return the bound of the UnivarPWLE <a>var</a>
  //! @param[in]  isMax specifies the upper bound to be filled in _bnd. Note the difference between it and the rec in ISModel::_B().  
  //! @param[in]  isMin specifies the lower bound to be filled in _bnd. Note the difference between it and the rec in ISModel::_B().  
  //! @param[in]  second is the vector of values (i.e., the coordinates of breakpoints in the y-axis)
  //! @param[out] _bnd stores the bounds of the var  
  void _B
  (const bool isMin,const bool isMax)
  const
  {
    // note that this functions should only be called when var is not constant
#ifdef MC__UPWLE_DEBUG  
    assert( !first.empty() );
    assert( second.size() > 2 );
#endif
    const T _temp_mid = second[0];
    const unsigned int offset = 1; 
    const int nEstmtrPts = second.size() - offset;    // excluding the average
    T ub(0.),lb(0.);
    if(nEstmtrPts > 0){
      ub = second[1]; lb = second[1];  
      for (unsigned int i = 2; i < second.size(); i++) {
        if(isMax) ub = ub > second[i] ? ub:second[i]; 
        if(isMin) lb = lb < second[i] ? lb:second[i];          
      }    
    } 
  
    if(isMin){
      _lbnd.second = true;       
      _lbnd.first  = _temp_mid + lb;      
    }

    if(isMax){
      _ubnd.second = true;       
      _ubnd.first  = _temp_mid + ub;        
    }
  }

  //! @brief Evaluate the forward difference of the UnivarPWLE <a>var</a>
  //! @param[in]  first 
  //! @param[in]  second 
  //! @param[out] _xFwdDiff 
  //! @param[out] _yFwdDiff 
  void _eval_fwdDiff()
  const
  {
    // note that this functions should only be called when var is not constant
    assert( !first.empty() );
    assert( second.size() > 2 );    

    const unsigned int offset = 1; 
    const unsigned int nEstmtrPts = second.size() - offset;    // excluding the average

    if(!(_xFwdDiff.second)){
      _xFwdDiff.first.resize(nEstmtrPts-1);
      for (unsigned int i = 1; i < nEstmtrPts; i++) {
        _xFwdDiff.first[i-1] =  first[i+1] - first[i];  
      }
      _xFwdDiff.second = true;  
    }

    if(!(_yFwdDiff.second)){
        _yFwdDiff.first.resize(nEstmtrPts-1);
      for (unsigned int i = 1; i < nEstmtrPts; i++) {
        _yFwdDiff.first[i-1] = second[i+1] - second[i]; 
      }  
      _yFwdDiff.second = true; 
    }
    // std::cout << "difference evaled" << std::endl;
    // for(unsigned int i = 0; i < nEstmtrPts - 1; i++)
    //     std::cout << "y" << i << ": " << _yFwdDiff.first[i] << std::endl;
  }


/*
  void _inv
  ( std::vector<std::vector<T>>& mat, unsigned const& ndep )
  const;
  void _sqr
  ( std::vector<std::vector<T>>& mat, unsigned const& ndep )
  const;
  void _sqrt
  ( std::vector<std::vector<T>>& mat, unsigned const& ndep )
  const;
  void _exp
  ( std::vector<std::vector<T>>& mat, unsigned const& ndep )
  const;
  void _log
  ( std::vector<std::vector<T>>& mat, unsigned const& ndep )
  const;
  void _xlog
  ( std::vector<std::vector<T>>& mat, unsigned const& ndep )
  const;
  void _sin
  ( std::vector<std::vector<T>>& mat, unsigned const& ndep )
  const;
  void _cos
  ( std::vector<std::vector<T>>& mat, unsigned const& ndep )
  const;

  // void _sin_slope
  // ( std::vector<std::vector<T>>& mat,std::vector<std::vector<std::vector<long double>>>& slope, unsigned const& ndep )
  // const;
  void _cos_slope
  ( std::vector<std::vector<T>>& mat,std::vector<std::vector<std::vector<long double>>>& slope, const std::vector<long double>& partitionSize, unsigned const& ndep )
  const;
  void _exp
  ( std::vector<std::vector<T>>& mat,std::vector<std::vector<std::vector<long double>>>& slope, const std::vector<long double>& partitionSize, unsigned const& ndep )
  const;
  void _log
  ( std::vector<std::vector<T>>& mat,std::vector<std::vector<std::vector<long double>>>& slope, const std::vector<long double>& partitionSize, unsigned const& ndep )
  const;
  void _xlog
  ( std::vector<std::vector<T>>& mat,std::vector<std::vector<std::vector<long double>>>& slope, const std::vector<long double>& partitionSize, unsigned const& ndep )
  const;
  void _sqrt
  ( std::vector<std::vector<T>>& mat,std::vector<std::vector<std::vector<long double>>>& slope, const std::vector<long double>& partitionSize, unsigned const& ndep )
  const;
  void _inv
  ( std::vector<std::vector<T>>& mat,std::vector<std::vector<std::vector<long double>>>& slope, const std::vector<long double>& partitionSize, unsigned const& ndep )
  const;

  void _sqr
  ( std::vector<std::vector<T>>& mat,std::vector<std::vector<std::vector<long double>>>& slope, const std::vector<long double>& partitionSize, unsigned const& ndep )
  const;

  void _tanh
  ( std::vector<std::vector<T>>& mat,std::vector<std::vector<std::vector<long double>>>& slope, const std::vector<long double>& partitionSize, unsigned const& ndep )
  const;


  void _pow
  ( std::vector<std::vector<T>>& mat,std::vector<std::vector<std::vector<long double>>>& slope, const std::vector<long double>& partitionSize, int const&iexp, unsigned const& ndep )
  const;

  void _prod
  ( std::vector<std::vector<T>> const& mat1, std::vector<std::vector<T>> const& mat2,
    std::vector<std::vector<T>>& mat3, unsigned& ndep3 )
  const;
  void _intersect
  ( std::vector<std::vector<T>>& mat, unsigned const& ndep, T _rangebnd )
  const;
  void _intersect
  ( std::vector<std::vector<T>>& mat, std::vector<std::vector<std::vector<long double>>>& slope, const std::vector<long double>& partitionSize, unsigned const& ndep, T _rangebnd )
  const;

  void _tanh
  ( std::vector<std::vector<T>>& mat, unsigned const& ndep )
  const;
  void _pow
  ( std::vector<std::vector<T>>& mat, int const&iexp, unsigned const& ndep )
  const;

  // @unsure: the algorithm for computing asymmetric over-/under-estimators is applicable for sin and cos 
  // when they are concave or convex on the range of the input ASM  


  // @unsure: there is a similar algorithm which can refining the multiplicative result, 
  // by bounding x_iy_j+x_jy_i in 0.5*(x_i^2+y_i^2)*[-1,1] + 0.5*(x_j^2+y_j^2)*[-1,1]
*/

  // template <typename U> friend U _interpolate_at_ind
  // (UnivarPWLE<T> const &, T , T , unsigned int );

  static T _interpolate_at_ind
  (std::vector<T> const & x, const T width, const T rangeDiff, const unsigned int ind);


  static void _underestimate_in_window4
  (const std::vector<T> & x, const std::vector<T> & y,std::vector<T> & data);

  static void _underestimate_in_window4
  (const std::vector<T> & x, const std::vector<T> & y, const unsigned int pos, const T p1AreaLoss, const T p2AreaLoss, std::vector<T> & data, const bool isUnder);

  static bool _cmp_indexed_areaLoss
  (const std::pair<T, int>& a, const std::pair<T, int>& b);


  static std::pair<T, int> _heap_pop
  (std::vector<std::pair<T, int>> & heap, std::vector<int> & hind);

  static void _modify_heap
  (std::vector<std::pair<T, int>> & heap, const int pos, std::vector<int> & hind);

  static void _heap_bubble_up
  (std::vector<std::pair<T, int>> & heap, int pos, std::vector<int> & hind);

  static void _heap_bubble_down
  (std::vector<std::pair<T, int>> & heap, int pos, std::vector<int> & hind);


  // static void _iterative_align
  // (UnivarPWLE<T> const & estmtrA, UnivarPWLE<T> const & estmtrB, std::vector<unsigned int> & ind, std::vector<T> & data, char mod = 1 )
  // const; 

  static void _iterative_smoothening_yFilter
  (std::vector<T> & xOut, std::vector<T> & yOut, std::vector<T> & yDiff, const bool isUnder, T & varSumPtValLast, T & varSum2bSet, const T & ptVar, unsigned int & nbpsAdded, bool & toMerge);

  static void _iterative_smoothening_xFilter
  (const T & ptVarDiff, const T & stepLast, const T & stepNext, T & PtValLast, const T & xStepLast, const T & xStepNext, const bool isUnder);

  static void _iterative_interpolate
  (const T & ptVarDiff, const T & step, T & PtValLast, const T & PtVal);


  static bool _smoothening_endpointFixer
  (std::vector<T> & data, std::vector<T> & endpoint, const std::vector<T> & yIn,const T & xEndpoint, const T & yEndpoint, const T & xCrntPwl2,const unsigned int nbpsProcsed, const unsigned int nBptPwl2, bool isUnder );


  void _condense_by_relax
  (const unsigned int k);

  void _add_and_altMerge
  ( const std::vector<T> & xIn, const std::vector<T> & yIn,std::pair<std::vector<T>,bool> xInFwdDiff, std::pair<std::vector<T>,bool> yInFwdDiff);

  void _relu();


  // static std::ostream& _dispmat
  // ( std::vector<PWLMAT> const & matlist, unsigned const len=3, std::ostream& out=std::cout );
  
  // std::ostream& _dispvar
  // ( std::vector<PWLMAT> const & matlist, unsigned const& ndep, const int& opt=0,
  //   std::ostream& out=std::cout )
  // const;

  // static std::ostream& _dispmat
  // ( std::vector<std::vector<T>> const& mat, unsigned const len=3, std::ostream& out=std::cout );
  
  // std::ostream& _dispvar
  // ( std::vector<std::vector<T>> const& mat, unsigned const& ndep, const int& opt=0,
  //   std::ostream& out=std::cout )
  // const;


};


template <typename T> inline
unsigned int UnivarPWLE<T>::nbpsMax;


template <typename T> 
std::ostream& operator<<
( std::ostream& out, UnivarPWLE<T> const& var )
{

  out << std::scientific << std::setprecision(14);
  out << " the PWL" << (var._isUnder? " under-": " over-") << "estimator";
  if(var.first.empty()){
    out << " is the constant " << var._cst << std::endl; 
    return out;
  }
  if(var.second.size() == 1){
    out << " is the constant " << var.second[0] << " on [" 
        << var.first[0] << ", " 
        << var.first[1] << "]" << std::endl;   
    return out;
  }

   if(var.second.size() == 2){
    out << " is linear [ " << var.second[0] << ", " << var.second[1]  << " ] on [" 
        << var.first[0] << ", " 
        << var.first[1] << "]" << std::endl;  
    return out;
  }

  out << ":" << std::endl;

  auto&& x = var.first;
  auto&& y = var.second;

  out << std::right << std::setw(4) <<"No.: ";
  for( unsigned int i=1; (i<var.second.size()) && (i<11) ; i++ ){
    if( (i-1) && !((i-1)%10)) out << std::endl << "      ";
    out << std::setw(6) << i << ",  ";
  }
  out << std::endl;


  out << std::right << std::setw(5) <<"x: ";
  for( unsigned int i=1; !var.first.empty() && i<var.second.size(); i++ ){
    if( (i-1) && !((i-1)%10)) out << std::endl << "      ";
    out << std::setw(6) << x[0] + x[i] << ",  ";
  }
  out << std::endl;

  out << std::right << std::setw(5) <<"y: ";
  for( unsigned int i=1; !var.first.empty() && i<var.second.size(); i++ ){
    if( (i-1) && !((i-1)%10)) out << std::endl << "      ";
//    out << std::setw(6) << std::setprecision(4) << y[0] + y[i]<< ",  ";
    out << std::setw(6) << y[0] + y[i]<< ",  ";
  }
  out << std::endl;

//  var._B(true,true);
  out << std::endl;
  out << std::right << std::setw(5) << "[LB, UB]: " << "[" <<
  var.get_lb() << ", " << var.get_ub() << "]" << std::endl;    
  return out;
}


//! static member function 
//! @brief Linear interpolation for aligning a univariate linear estimator with a univariate piecewise linear estimator w.r.t. the same variable
//! @param[in] pos is
//! @param[in] p1AreaLoss is 
//! @param[in] p2AreaLoss is t
//! @param[out] data 
template <typename T>
inline
void UnivarPWLE<T>::_underestimate_in_window4
(const std::vector<T> & x, const std::vector<T> & y,std::vector<T> & data)
{
  const T p1AreaLoss = (x[2]-x[1])*(y[1]-y[0])-(x[1]-x[0])*(y[2]-y[1]);
  const T p2AreaLoss = (x[3]-x[2])*(y[2]-y[1])-(x[2]-x[1])*(y[3]-y[2]);
  
  _underestimate_in_window4(x,y,0,p1AreaLoss,p2AreaLoss,data,true);
}


//! static member function 
//! @brief Linear interpolation for aligning a univariate linear estimator with a univariate piecewise linear estimator w.r.t. the same variable
//! @param[in] pos is
//! @param[in] p1AreaLoss is 
//! @param[in] p2AreaLoss is t
//! @param[out] data 
template <typename T>
inline
void UnivarPWLE<T>::_underestimate_in_window4
(const std::vector<T> & x, const std::vector<T> & y, const unsigned int pos, const T p1AreaLoss, const T p2AreaLoss, std::vector<T> & data, const bool isUnder)
{

  T & vtxX  = data[0];
  T & vtxY  = data[1];
  T & aLoss = data[2];        

  //std::vector<T> & x = first;
  //std::vector<T> & y = second;
  const unsigned int id0(pos),id1(pos+1),id2(pos+2),id3(pos+3);

  if(p1AreaLoss >= 0. && p2AreaLoss >= 0.){ // ccv ccv
      if(p1AreaLoss < p2AreaLoss){
          vtxX = x[id2];
          vtxY = y[id2];
          aLoss = p1AreaLoss;
      }
      else{
          vtxX = x[id1];
          vtxY = y[id1];
          aLoss = p2AreaLoss;
      }
  }
  else if(p2AreaLoss >= 0.){ // cvx ccv
      const T p3AreaLoss = isUnder?((y[id1]-y[id0])*(x[id3]-x[id0]) - (y[id3]-y[id0])*(x[id1]-x[id0])):((y[id3]-y[id0])*(x[id1]-x[id0]) - (y[id1]-y[id0])*(x[id3]-x[id0]));
      if(p3AreaLoss >=0.){
          const T ratio = p1AreaLoss / (p1AreaLoss - p3AreaLoss);
          const T ratioCpl = 1-ratio;
          //output = cross_point([x[0],x[1],x[2],x[3]],[y[0],y[1],y[2],y[3]])
          vtxX = x[id2]*(ratioCpl) + x[id3]*ratio; // output[0];
          vtxY = y[id2]*(ratioCpl) + y[id3]*ratio;
          aLoss = p2AreaLoss * ratio; // (output[0]-x[2])*(y[2]-y[1])-(x[2]-x[1])*(output[1]-y[2]);
      }
      else{
          vtxX = x[id1];
          vtxY = y[id1];
          aLoss = p2AreaLoss;
      }
  }
  else if(p1AreaLoss >= 0.){ // ccv cvx
      const T p3AreaLoss = isUnder?((y[id3]-y[id0])*(x[id3]-x[id2]) - (y[id3]-y[id2])*(x[id3]-x[id0])):((y[id3]-y[id2])*(x[id3]-x[id0]) - (y[id3]-y[id0])*(x[id3]-x[id2]));
      if (p3AreaLoss>=0.){
          const T ratio = p2AreaLoss / (p2AreaLoss - p3AreaLoss);
          const T ratioCpl = 1-ratio;
          //output = cross_point([x[0],x[1],x[2],x[3]],[y[0],y[1],y[2],y[3]])
          vtxX = x[id1]*(ratioCpl) + x[id0]*ratio; //utput[0];
          vtxY = y[id1]*(ratioCpl) + y[id0]*ratio; //output[1];
          aLoss = p1AreaLoss * ratio; //(x[2]-x[1])*(y[1]-output[1])-(x[1]-output[0])*(y[2]-y[1]);
      }
      else{
          vtxX = x[id2];
          vtxY = y[id2];
          aLoss = p1AreaLoss;
      }
  }
  else{
      const T p3AreaLoss = isUnder?((y[id3]-y[id0])*(x[id3]-x[id2]) - (y[id3]-y[id2])*(x[id3]-x[id0])):((y[id3]-y[id2])*(x[id3]-x[id0]) - (y[id3]-y[id0])*(x[id3]-x[id2]));
      const T ratio = p2AreaLoss / (p2AreaLoss - p3AreaLoss );
      const T ration1 = 1-ratio;
      //output = cross_point([x[0],x[1],x[2],x[3]],[y[0],y[1],y[2],y[3]])
      vtxX = x[id1]*(ration1) + x[id0]*ratio; //output[0]
      vtxY = y[id1]*(ration1) + y[id0]*ratio; //output[1]
      aLoss = ratio * p1AreaLoss; //-(x[2]-output[0])*(output[1]-y[1])+(output[0]-x[1])*(y[2]-output[1]);
  }
}


//! @brief Linear interpolation for aligning a univariate linear estimator with a univariate piecewise linear estimator w.r.t. the same variable
//! @param[in] nbpsMax is
//! @param[in] isUnder is 
//! @param[out] first 
//! @param[out] second
template <typename T>
inline
void UnivarPWLE<T>::_condense_by_relax
(const unsigned int k)
{

  //std::cout <<  "        in condense" << std::endl;
  // k is nbpsMax 
  std::vector<T> & xori = first;
  std::vector<T> & yori = second;
  const unsigned int n = first.size() - 1;
  //_isUnder;

    T p1AreaLoss;
    T p2AreaLoss;
    std::vector<T> data(3);
    std::vector<T> vtx2addX(n-3);
    std::vector<T> vtx2addY(n-3);
    std::vector<std::pair<T, int>> areaLoss(n-3);
    
    _eval_fwdDiff(); 
    //std::cout << "diffsize " << _xFwdDiff.first.size() << std::endl;
    //std::cout << "diffsize " << _yFwdDiff.first.size() << std::endl;
   
    p1AreaLoss = _isUnder?(_xFwdDiff.first[1]*_yFwdDiff.first[0]-_xFwdDiff.first[0]*_yFwdDiff.first[1]):(_xFwdDiff.first[0]*_yFwdDiff.first[1]-_xFwdDiff.first[1]*_yFwdDiff.first[0]);
    for(unsigned int j = 0; j < n - 3; j++){       
      p2AreaLoss = _isUnder?(_xFwdDiff.first[j+2]*_yFwdDiff.first[j+1]-_xFwdDiff.first[j+1]*_yFwdDiff.first[j+2]):(_xFwdDiff.first[j+1]*_yFwdDiff.first[j+2]-_xFwdDiff.first[j+2]*_yFwdDiff.first[j+1]);
      _underestimate_in_window4(first,second,j+1,p1AreaLoss,p2AreaLoss,data,_isUnder);
      vtx2addX[j] = data[0];//vtxX
      vtx2addY[j] = data[1];//vtxY
      areaLoss[j] = std::make_pair(data[2],j);// AreaLoss(aLoss,j)
      //std::cout << "j: " << j << " area: " << data[2] << std::endl;
      p1AreaLoss = p2AreaLoss;
    }

    if(_isUnder){
      if(!_lbnd.second) _B(true,false);
      _ubnd.second = false;
      // else {std::cout << "not recomputed" << std::endl;
      // auto temp = _bnd.first.first.first;
      // _B(true,false);
      // if(temp != _bnd.first.first.first) {
      //   std::cout << "und bound mismatch" << std::endl;
      //   throw typename UnivarPWLE<T>::Exceptions( UnivarPWLE<T>::Exceptions::UNDEF );
      // }
      // }
    }
    else{
      if(!_ubnd.second) _B(false,true);
      _lbnd.second = false;
      // else {std::cout << "not recomputed" << std::endl;
      // auto temp = _bnd.first.second.first;
      // _B(false,true);
      // if(temp != _bnd.first.second.first) {
      //   std::cout << "ove bound mismatch" << std::endl;
      //   throw typename UnivarPWLE<T>::Exceptions( UnivarPWLE<T>::Exceptions::UNDEF );
      // }
      // }
    }
    
    const T yMinOrMax = (_isUnder?(_lbnd.first):(_ubnd.first))-yori[0];

    std::make_heap(areaLoss.begin(), areaLoss.end(), _cmp_indexed_areaLoss);

    std::vector<int> hind(n-3);
    for (unsigned int i = 0; i < n-3; i++) {
        hind[areaLoss[i].second] = (int) i;
    }        
    
    //std::pop_heap(v2.begin(), v2.end(), cmpa); 
    //v2.pop_back();

    std::vector<int> nxtInd(n); // = [i for i in range(1,n+1)]
    std::vector<int> prvInd(n); // = [i for i in range(-1,n-1,1)]
    for (int i = 0; i < ((int) n); i++) {
        prvInd[i] = i - 1;  
        nxtInd[i] = i + 1;        
    }                 

    auto & aMin = areaLoss[0]; // bind the reference to the top of the heap, heappop(areaLoss,hind) in python 
    // std::cout << "aMin " << aMin.first << " at " << aMin.second << std::endl;
    // std::cout << "areaLoss len" << areaLoss.size() << std::endl;
    // std::cout << "aMin1 " << areaLoss[1].first << " at " << areaLoss[1].second << std::endl;
    if (_isUnder){
      while (vtx2addY[aMin.second] + 1e-14 < yMinOrMax){
        areaLoss[0].first = DBL_MAX;
        _modify_heap(areaLoss,0,hind);
      
        // std::cout << "aMin " << aMin.first << " at " << aMin.second << std::endl;
        // if(areaLoss.size() <= 2) {
        //   std::cout << " here we only have 1-2 choices" << std::endl;
        //   throw typename UnivarPWLE<T>::Exceptions( UnivarPWLE<T>::Exceptions::UNDEF );
        //}        
        //aMin = areaLoss[0]; // as the reference binds to that, we do not need assignment
      }
    }
    else{
      while (vtx2addY[aMin.second] - 1e-14 > yMinOrMax){
        areaLoss[0].first = DBL_MAX;
        _modify_heap(areaLoss,0,hind);
        //aMin = areaLoss[0]; // as the reference binds to that, we do not need assignment
      }       
    }    
    //std::cout <<  "            aMin.second " << aMin.second << " with " << std::scientific << std::setprecision(14) << aMin.first << std::endl;    
    // updating the indexes so that we can skip the second element in the window, and 
    //   update the third one with the computed one
    // std::cout << "update the third one with the computed one" << std::endl;
    int pos = aMin.second;
    int id_1 = nxtInd[aMin.second];
    int id   = nxtInd[id_1];
    prvInd[id]   = aMin.second;
    nxtInd[aMin.second] = id;
    xori[id+1] = vtx2addX[pos];
    yori[id+1] = vtx2addY[pos];

   // update the forward difference
    _xFwdDiff.first[aMin.second] = xori[id+1] - xori[aMin.second+1];
    _xFwdDiff.first[id] = xori[nxtInd[id]+1] - xori[id+1];
    _yFwdDiff.first[aMin.second] = yori[id+1] - yori[aMin.second+1];
    _yFwdDiff.first[id] = yori[nxtInd[id]+1] - yori[id+1];    

    //T accAreaLoss = aMin.first;
    //std::cout << "n-k-1 " << n-k-1 << std::endl;    
    for(unsigned int i = 1; i < n-k+1; i ++ ){
      _heap_pop(areaLoss,hind);
      const int idx1 = pos;
      const int idx2 = id; 
      const int idx3 = nxtInd[id]; 

      if (pos > 0){
        const int idx0 = prvInd[idx1];
        const std::vector<T> xIn{xori[idx0+1],xori[idx1+1],xori[idx2+1],xori[idx3+1]}; 
        const std::vector<T> yIn{yori[idx0+1],yori[idx1+1],yori[idx2+1],yori[idx3+1]};
        //_underestimate_in_window4(xIn,yIn,data);
        p1AreaLoss = _isUnder?(_xFwdDiff.first[idx1])*(_yFwdDiff.first[idx0])-(_xFwdDiff.first[idx0])*(_yFwdDiff.first[idx1]):(_xFwdDiff.first[idx0])*(_yFwdDiff.first[idx1])-(_xFwdDiff.first[idx1])*(_yFwdDiff.first[idx0]);
        p2AreaLoss = _isUnder?(_xFwdDiff.first[idx2])*(_yFwdDiff.first[idx1])-(_xFwdDiff.first[idx1])*(_yFwdDiff.first[idx2]):(_xFwdDiff.first[idx1])*(_yFwdDiff.first[idx2])-(_xFwdDiff.first[idx2])*(_yFwdDiff.first[idx1]);
        _underestimate_in_window4(xIn,yIn,0,p1AreaLoss,p2AreaLoss,data,_isUnder);
        //print("idx0",idx0,"aLoss",aLoss)
        areaLoss[hind[idx0]].first = data[2]; // aLoss;
        //print("check id0: ",areaLoss[hind[idx0]].id == idx0)
        vtx2addX[idx0] = data[0]; // vtxX;
        vtx2addY[idx0] = data[1]; // vtxY;
        _modify_heap(areaLoss,hind[idx0],hind);
        // std::cout << "data[0] " << data[0] << std::endl;
        // std::cout << "data[1] " << data[1] << std::endl;     
        // std::cout << "data[2] " << data[2] << std::endl;            
        //print("check mod0 ",h_check(areaLoss))
        //print("AL0: ",areaLoss)
        if (idx0 > 0){
          const int idx_1 = prvInd[idx0];
          const std::vector<T> xIn{xori[idx_1+1],xori[idx0+1],xori[idx1+1],xori[idx2+1]}; 
          const std::vector<T> yIn{yori[idx_1+1],yori[idx0+1],yori[idx1+1],yori[idx2+1]};                
          //_underestimate_in_window4(xIn,yIn,data);
          p2AreaLoss = p1AreaLoss;// (_xFwdDiff.first[idx1])*(_yFwdDiff.first[idx0])-(_xFwdDiff.first[idx0])*(_yFwdDiff.first[idx1]);
          p1AreaLoss = _isUnder?(_xFwdDiff.first[idx0])*(_yFwdDiff.first[idx_1])-(_xFwdDiff.first[idx_1])*(_yFwdDiff.first[idx0]):(_xFwdDiff.first[idx_1])*(_yFwdDiff.first[idx0])-(_xFwdDiff.first[idx0])*(_yFwdDiff.first[idx_1]);
          _underestimate_in_window4(xIn,yIn,0,p1AreaLoss,p2AreaLoss,data,_isUnder);          
          //print("idx_1",idx_1,"aLoss",aLoss)
          areaLoss[hind[idx_1]].first = data[2];
          vtx2addX[idx_1] = data[0];
          vtx2addY[idx_1] = data[1];
          _modify_heap(areaLoss,hind[idx_1],hind);
            //print("check mod-1 ",h_check(areaLoss))
            //print("AL-1: ",areaLoss)
          // std::cout << "data[0] " << data[0] << std::endl;
          // std::cout << "data[1] " << data[1] << std::endl;     
          // std::cout << "data[2] " << data[2] << std::endl;                
        }
      } 
      if (idx3 < ((int )n - 1)){
          const int idx4 = nxtInd[idx3];
          const std::vector<T> xIn{xori[idx1+1],xori[idx2+1],xori[idx3+1],xori[idx4+1]}; 
          const std::vector<T> yIn{yori[idx1+1],yori[idx2+1],yori[idx3+1],yori[idx4+1]};     
          p1AreaLoss = _isUnder?(_xFwdDiff.first[idx2])*(_yFwdDiff.first[idx1])-(_xFwdDiff.first[idx1])*(_yFwdDiff.first[idx2]):(_xFwdDiff.first[idx1])*(_yFwdDiff.first[idx2])-(_xFwdDiff.first[idx2])*(_yFwdDiff.first[idx1]);
          p2AreaLoss = _isUnder?(_xFwdDiff.first[idx3])*(_yFwdDiff.first[idx2])-(_xFwdDiff.first[idx2])*(_yFwdDiff.first[idx3]):(_xFwdDiff.first[idx2])*(_yFwdDiff.first[idx3])-(_xFwdDiff.first[idx3])*(_yFwdDiff.first[idx2]);
          _underestimate_in_window4(xIn,yIn,0,p1AreaLoss,p2AreaLoss,data,_isUnder);                          
          //_underestimate_in_window4(xIn,yIn,data);
          //print("id_1",id_1,"aLoss",aLoss)
          areaLoss[hind[id_1]].first = data[2];
          areaLoss[hind[id_1]].second = idx1;
          //hind[idx1] = hind[id_1]
          _modify_heap(areaLoss,hind[id_1],hind);
          //print("check mod1 ",h_check(areaLoss))
          vtx2addX[idx1] = data[0]; // vtxX
          vtx2addY[idx1] = data[1]; // vtxY 
          //print("AL1: ",areaLoss)
          // std::cout << "data[0] " << data[0] << std::endl;
          // std::cout << "data[1] " << data[1] << std::endl;     
          // std::cout << "data[2] " << data[2] << std::endl;    
          if (idx4 < ((int )n - 1)){
            const int idx5 = nxtInd[idx4];
            const std::vector<T> xIn{xori[idx2+1],xori[idx3+1],xori[idx4+1],xori[idx5+1]}; 
            const std::vector<T> yIn{yori[idx2+1],yori[idx3+1],yori[idx4+1],yori[idx5+1]};                     
            p1AreaLoss = p2AreaLoss; //(_xFwdDiff.first[idx3])*(_yFwdDiff.first[idx2])-(_xFwdDiff.first[idx2])*(_yFwdDiff.first[idx3]);
            p2AreaLoss = _isUnder?(_xFwdDiff.first[idx4])*(_yFwdDiff.first[idx3])-(_xFwdDiff.first[idx3])*(_yFwdDiff.first[idx4]):(_xFwdDiff.first[idx3])*(_yFwdDiff.first[idx4])-(_xFwdDiff.first[idx4])*(_yFwdDiff.first[idx3]);
            _underestimate_in_window4(xIn,yIn,0,p1AreaLoss,p2AreaLoss,data,_isUnder);    
            //_underestimate_in_window4(xIn,yIn,data);
            //print("idx2",idx2,"aLoss",aLoss)
            areaLoss[hind[idx2]].first = data[2]; // aLoss
            vtx2addX[idx2] = data[0]; //vtxX
            vtx2addY[idx2] = data[1]; //vtxY
            _modify_heap(areaLoss,hind[idx2],hind);
            // std::cout << "data[0] " << data[0] << std::endl;
            // std::cout << "data[1] " << data[1] << std::endl;     
            // std::cout << "data[2] " << data[2] << std::endl;                      
          }
      } 

      // std::cout << "aMin " << areaLoss[0].first << " at " << areaLoss[0].second << std::endl;
      // std::cout << "aMin1 " << areaLoss[1].first << " at " << areaLoss[1].second << std::endl;


      // if(aMin.first != areaLoss[0].first) std::cout << "error" << std::endl; // test to make sure the reference binds to the heap top
      // aMin = areaLoss[0]; // as the reference binds to that, we do not need assignment, heappop(areaLoss,hind) in python

      if (_isUnder){
        while (vtx2addY[aMin.second] + 1e-14 < yMinOrMax){
            areaLoss[0].first = DBL_MAX;
            _modify_heap(areaLoss,0,hind);
            //aMin = areaLoss[0];  // as the reference binds to that, we do not need assignment
        }
      }
      else{
        while (vtx2addY[aMin.second] - 1e-14 > yMinOrMax){
            areaLoss[0].first = DBL_MAX;
            _modify_heap(areaLoss,0,hind);
            //aMin = areaLoss[0];  // as the reference binds to that, we do not need assignment
        }           
      }
      //std::cout <<  "              aMin.second " << aMin.second << " with " << std::scientific << std::setprecision(14) << aMin.first << std::endl;
      pos = aMin.second;
      id_1 = nxtInd[aMin.second];
      id   = nxtInd[id_1];
      prvInd[id]   = aMin.second;
      nxtInd[aMin.second] = id;
      xori[id+1] = vtx2addX[pos];
      //print("vtx2addX[pos]: ", vtx2addX[pos])
      yori[id+1] = vtx2addY[pos];
      _xFwdDiff.first[aMin.second] = xori[id+1] - xori[aMin.second+1];
      _xFwdDiff.first[id] = xori[nxtInd[id]+1] - xori[id+1];
      _yFwdDiff.first[aMin.second] = yori[id+1] - yori[aMin.second+1];
      _yFwdDiff.first[id] = yori[nxtInd[id]+1] - yori[id+1];
      //print("aMin.al: ", aMin.al)
      //accAreaLoss = accAreaLoss + aMin.first;
    }



    std::vector<T> xlst(k); 
    std::vector<T> ylst(k);

    xlst[0] = xori[0];
    xlst[1] = xori[1];    
    ylst[0] = yori[0];
    ylst[1] = yori[1];    
    unsigned int i = 0;
    for(unsigned int j = 2; j < k-1; j++){
      i = nxtInd[i];
      xlst[j] = xori[i+1];
      ylst[j] = yori[i+1];
      _xFwdDiff.first[j-1] = _xFwdDiff.first[i];
      _yFwdDiff.first[j-1] = _yFwdDiff.first[i];

      //_xFwdDiff.first[]
      //_yFwdDiff.first[]
    }
    i = nxtInd[i];
    xlst[k-1] = xori[i+1];
    ylst[k-1] = yori[i+1];    
    
    first.swap(xlst);
    second.swap(ylst);
    _xFwdDiff.first.resize(k-2);
    _yFwdDiff.first.resize(k-2);

    //std::cout <<  "leave condense" << std::endl;

}

template <typename T>
inline
bool UnivarPWLE<T>::_cmp_indexed_areaLoss
(const std::pair<T, int>& a, const std::pair<T, int>& b){
    return a.first > b.first;
}


template <typename T>
inline
std::pair<T, int> UnivarPWLE<T>::_heap_pop
(std::vector<std::pair<T, int>> & heap, std::vector<int> & hind)
{
  // Pop the smallest item off the heap, maintaining the heap invariant."""
  auto lastelt =  heap[heap.size()-1];
  heap.pop_back();    // raises appropriate IndexError if heap is empty
  if (!heap.empty()){
    auto returnitem = heap[0];
    hind[heap[0].second] = -1;
    heap[0] = lastelt;
    hind[lastelt.second] = 0;
    _heap_bubble_down(heap, 0, hind);
    return returnitem;
  }
  hind[lastelt.second] = -1;
  return lastelt;
}


template <typename T>
inline 
void UnivarPWLE<T>::_modify_heap
(std::vector<std::pair<T, int>> & heap, const int pos, std::vector<int> & hind)
{
  if (pos == 0){
    _heap_bubble_down(heap, pos, hind);
    //print("pos",pos,"heap[pos].al",heap[pos].al)
  }
  else{
    unsigned int parentpos = (pos - 1) >> 1;
    // print("pos",pos,"parentpos",parentpos)
    // print("heap[pos]",heap[pos].al,"heap[parentpos]",heap[parentpos].al)
    if (heap[pos].first < heap[parentpos].first){
      std::swap(heap[pos],heap[parentpos]);
      // print("heap[pos]",heap[pos],"heap[parentpos]",heap[parentpos])
      hind[heap[pos].second] = pos;
      hind[heap[parentpos].second] = parentpos;
      _heap_bubble_up(heap, parentpos, hind);
    }
    else
      _heap_bubble_down(heap, pos, hind);
  }

}

template <typename T>
inline 
void UnivarPWLE<T>::_heap_bubble_up
(std::vector<std::pair<T, int>> & heap, int pos, std::vector<int> & hind)
{
    auto newitem = heap[pos];
    // Follow the path to the root, moving parents down until finding a place
    // newitem fits.
    while (pos > 0){
        unsigned int parentpos = (pos - 1) >> 1;
        auto parent = heap[parentpos];
        if (newitem.first < parent.first){
            heap[pos] = parent;
            hind[parent.second] = pos;
            pos = parentpos;
            continue;
        }
        break;
    }
    heap[pos] = newitem;
    hind[newitem.second] = pos;
}


template <typename T>
inline 
void UnivarPWLE<T>::_heap_bubble_down
(std::vector<std::pair<T, int>> & heap, int pos, std::vector<int> & hind)
{
    // print("bubble_down")
    int endpos = heap.size();
    //int startpos = pos;
    auto newitem = heap[pos];
    // print("endpos ", endpos,"startpos ",startpos,"newitem ",newitem)
    //  Bubble up the smaller child until hitting a leaf.
    int childpos = 2*pos + 1;    // leftmost child position
    while (childpos < endpos){
        // Set childpos to index of smaller child.
        int rightpos = childpos + 1;
        // print("left: heap[childpos] ", heap[childpos], "heap[rightpos]", heap[rightpos])
        if (rightpos < endpos && (heap[childpos].first >= heap[rightpos].first)){
            childpos = rightpos;
        // Move the smaller child up.
        }
        if (newitem.first > heap[childpos].first){
            heap[pos] = heap[childpos];
            hind[heap[pos].second] = pos;
            pos = childpos;
            childpos = 2*pos + 1;
        }
        else break;
    }
    // The leaf at pos is empty now.  Put newitem there, and bubble it up
    // to its final resting place (by sifting its parents down).
    heap[pos] = newitem;
    hind[newitem.second] = pos;
}







//! @brief Linear interpolation for aligning a univariate linear estimator with a univariate piecewise linear estimator w.r.t. the same variable
//! @param[in] x 
//! @param[in] width is the width of the var x_i (i.e., w(x_i) = x_i^U - x_i^L)
//! @param[in] rangeDiff is the range of the linear estimator z_i(x_i) w.r.t. to the same variable x_i, 
//! @param[in] ind specifies the index of the breakpoint of the estimator y_i(x_i)
//! @return returns the value z_i(x_i[ind]) - z_i(x_i^L), which will be used with y_i(x_i[ind]) - y_i(x_i^L) (which is given by estmtr.second[ind])
template <typename T>
inline
T UnivarPWLE<T>::_interpolate_at_ind
(std::vector<T> const & x, const T width, const T rangeDiff, const unsigned int ind)
{
  T offsetVar = x[ind];
  offsetVar -= x[1]; // assuming that only the first element estmtr.first[0] is reserved for numerical regularisation, e.g.,
                                // here the estmtr.first[0] stores the mean value of x_i
  return (offsetVar / width)*rangeDiff;
}



//! @brief static, const
//! @param[in] ptVarDiff is 
//! @param[in] step 
//! @param[in,out] PtValLast
//! @param[in] PtVal 
template <typename T>
inline
void UnivarPWLE<T>::_iterative_interpolate
(const T & ptVarDiff, const T & xStep, T & PtValLast, const T & PtVal)
{
  const T ratio  = ptVarDiff / xStep;   
  const T offset = ratio * (PtVal - PtValLast);
  PtValLast = PtVal - offset;   
}


//! @brief static, const
//! @param[in] ptVarDiff 
//! @param[in] stepLast ,  
//! @param[in] stepNext 
//! @param[in,out] PtValLast 
//! @param[in] xStepLast ,  
//! @param[in] xStepNext 
//! @param[in] isUnder 
template <typename T>
inline
void UnivarPWLE<T>::_iterative_smoothening_xFilter
(const T & ptVarDiff, const T & stepLast, const T & stepNext, T & PtValLast,
 const T & xStepLast, const T & xStepNext, const bool isUnder)
{
  const T candidateL = (stepLast/xStepLast)*ptVarDiff;
  //const T candidateR = stepNext*( ptVarDiff/xStepNext + 1.0);
  const T candidateR = (stepNext/xStepNext)*ptVarDiff;
  if(isUnder)
    PtValLast  += (candidateL < candidateR ? candidateL:candidateR);
  else
    PtValLast  += (candidateL > candidateR ? candidateL:candidateR);
}


//! @brief static, const
//! @param[in]  
//! @param[in]    
//! @param[in,out] 
//! @param[in,out] 
template <typename T>
inline
bool UnivarPWLE<T>::_smoothening_endpointFixer
(std::vector<T> & data, std::vector<T> & endpoint, const std::vector<T> & yIn, 
 const T & xEndpoint, const T & yEndpoint, const T & xCrntPwl2,
 const unsigned int nbpsProcsed, const unsigned int nBptPwl2, bool isUnder )
{
  //std::cout << "in fixer" << std::endl;
  // Most lines of this functions are devoted to handle potential numerical issues 
  // by extracting (outer-estimation) one breakpoint from filtering several breakpoints which are too close
  // here the long-double is to avoid large rounding errors  

  T & xAddedLast  = data[0];  // value of the last breakpoint that has been processed 
  T & yLastPwl1   = data[1];  // value of the last estmtrA's breakpoint that has been processed 
  T & yLastPwl2   = data[2];  // value of the last estmtrB's breakpoint that has been processed  
  // T   xCrntPwl2      = estmtrB.first[nbpsProcsed];

  // Case 1 varA > varB 
  // we can only process the smaller one (i.e. varB), because we don't know if there is other points in B smaller than the larger one (i.e. varA)
  // it should be noted that, in this case varB cannot be the right endpoint. 
   
  // Case 1.1  varA = right endpoint in this case all later varB are all >= varA
  // hence that we need to either 
  // find the valB at right endpoint such that all points are below/over, or 
  // find the maximum/minimum valB and use this value at the current point 
  // then use endpointFixer to modify the endpoint. 
  T extremaPWL = yIn[nbpsProcsed]; // abuse the data-container to find and store the maximum/minimum valB
  if(!isUnder){  // in such cases B is overestimator so we need to find the max
    for ( unsigned int j = nbpsProcsed+1; j < nBptPwl2; j++ ){ 
      extremaPWL = extremaPWL > yIn[j] ? extremaPWL:yIn[j];
    }
    endpoint[0] = yEndpoint;
    endpoint[1] = extremaPWL;       
    if( yLastPwl2 >= extremaPWL){   // in this case we can directly fix the endpoint to be extremaPWL
      // yLastPwl1 = yEndpoint;
      // xAddedLast = xEndpoint;
      // yLastPwl2 = extremaPWL;
      return false;   
    }             
    else{
      yLastPwl1 += (yEndpoint - yLastPwl1) * ((xCrntPwl2 - xAddedLast)/(xEndpoint - xAddedLast));    //update PWL1 at xCrntPwl2
      xAddedLast = xCrntPwl2;   
      yLastPwl2 = extremaPWL;//yIn[nbpsProcsed];//
      //std::cout << "true" << std::endl;
      return true;      
    } 
  }
  else{         // in such cases B is underestimator so we need to find the min
    for ( unsigned int j = nbpsProcsed+1; j < nBptPwl2; j++ ){ 
      extremaPWL = extremaPWL < yIn[j] ? extremaPWL:yIn[j];  
    }
    endpoint[0] = yEndpoint;
    endpoint[1] = extremaPWL;  
    if( yLastPwl2 <= extremaPWL){  // in this case we can directly fix the endpoint to be extremaPWL
      //yLastPwl1 = yEndpoint;
      //xAddedLast  = xEndpoint;
      //yLastPwl2 = extremaPWL;    
      return false;       
    }
    else{
      yLastPwl1 += (yEndpoint - yLastPwl1) * ((xCrntPwl2 - xAddedLast)/(xEndpoint - xAddedLast));    //update PWL1 at xCrntPwl2
      xAddedLast = xCrntPwl2;      
      yLastPwl2 = extremaPWL;
      return true;          
    }                                
  } 

}




//! @brief static, const
//! @param[in,out] xOut is a univariate piecewise linear estimator y_i(x_i) w.r.t. to a variable x_i,  
//! @param[in,out] yOut is a univariate piecewise linear estimator z_i(x_i) w.r.t. to the variable x_i,  
//! @param[in] smthCond specifies the index of the next two breakpoints of y_i and z_i, to be updated after the current processing
//! @param[in] isUnder specifies the index of the next two breakpoints of y_i and z_i, to be updated after the current processing
//! @param[in,out] varSumPtValLast consists of the values at last (aligned) point, to be updated as vals at current (just aligned) point, indexing to the next vars
//! @param[in,out] varSum2bSet specifies the one out of three modes: =0 if both estmtrA and estmtrA are underestimators; =2 if both are overestimators; and = 1 if A is under and B is over
//! @param[in,out] nbpsAdded specifies the one out of three modes: =0 if both estmtrA and estmtrA are underestimators; =2 if both are overestimators; and = 1 if A is under and B is over
//! @param[in,out] toMerge specifies the one out of three modes: =0 if both estmtrA and estmtrA are underestimators; =2 if both are overestimators; and = 1 if A is under and B is over

template <typename T>
inline
void UnivarPWLE<T>::_iterative_smoothening_yFilter
(std::vector<T> & xOut, std::vector<T> & yOut, std::vector<T> & yDiff, const bool isUnder,
 T & varSumPtValLast, T & varSum2bSet, const T & ptVar, unsigned int & nbpsAdded, bool & toMerge)
{

  if( std::fabs( yDiff[nbpsAdded] )  >= T(MC__UPWLE_SMOOTH_TOL) ){
    // to add the current breakpoint, set the counter, update the lastpoint (for merging), and reset the merging flag              
    xOut[2+nbpsAdded] = ptVar;
    yOut[2+nbpsAdded] = varSum2bSet;
    varSumPtValLast   = varSum2bSet;
    nbpsAdded ++;       
    toMerge = false;    
  }  
  else {
    if(!toMerge){
      // to add the current breakpoint, set the counter, update the lastpoint (for merging), and reset the merging flag              
      xOut[2+nbpsAdded] = ptVar;
      yOut[2+nbpsAdded] = varSum2bSet;
      varSumPtValLast   = varSum2bSet;
      nbpsAdded ++;             
      toMerge = true;
    }
    else{
      // smoothening filter
      xOut[1+nbpsAdded] = ptVar;    // update the last coordinate in x-axis
      // the following is to update the last coordinate in y-axis
      if(isUnder){
        varSumPtValLast = varSumPtValLast < varSum2bSet? varSumPtValLast:varSum2bSet;
        varSumPtValLast = varSumPtValLast < yOut[nbpsAdded] ? varSumPtValLast:yOut[nbpsAdded];
      }          
      else{
        varSumPtValLast = varSumPtValLast > varSum2bSet? varSumPtValLast:varSum2bSet; 
        varSumPtValLast = varSumPtValLast > yOut[nbpsAdded] ? varSumPtValLast:yOut[nbpsAdded];          
      } 
      yOut[1+nbpsAdded] = varSumPtValLast;
      yDiff[nbpsAdded - 1] += varSumPtValLast - yOut[nbpsAdded];
      yOut[nbpsAdded]   = varSumPtValLast;
      yDiff[nbpsAdded]  = T(0.);
    }
    
  }  

}

//! @brief not static, not const
//! @param[in] xIn is a univariate piecewise linear estimator y_i(x_i) w.r.t. to a variable x_i,  
//! @param[in] yIn is a univariate piecewise linear estimator z_i(x_i) w.r.t. to the variable x_i,  
//! @param[in,out] first specifies the index of the next two breakpoints of y_i and z_i, to be updated after the current processing
//! @param[in,out] second consists of the values at last (aligned) point, to be updated as vals at current (just aligned) point, indexing to the next vars
//! @param[in] _isUnder specifies the one out of three modes: =0 if both estmtrA and estmtrA are underestimators; =2 if both are overestimators; and = 1 if A is under and B is over
template <typename T>
inline
void UnivarPWLE<T>::_add_and_altMerge
( const std::vector<T> & xIn, const std::vector<T> & yIn, 
  std::pair<std::vector<T>,bool> xInFwdDiff, std::pair<std::vector<T>,bool> yInFwdDiff)
{

#ifdef MC__UPWLE_DEBUG_TRACE
  std::cout <<"      in altMerge " << std::endl;
#endif  
  const unsigned int nBptPwl1 = second.size();
  const unsigned int nBptPwl2 = yIn.size();
//  const unsigned int offset = 1;      // the number of reserved elements, and the index of the first breakpoint


  const unsigned int outSizeMax = nBptPwl1 + nBptPwl2 - 3; // excluding the two left endpoint and one average
  std::vector<T> xOut(outSizeMax); 
  std::vector<T> yOut(outSizeMax); 
  //xOut.resize(nbps2BPrcsed);      // excluding two right endpoints and adding one left endpoint and one right endpoint 
  //yOut.resize(nbps2BPrcsed);      // excluding two right endpoints and adding one left endpoint and one right endpoint 
  
  // to process the average and the left endpoint 
  xOut[0] = first[0];
  xOut[1] = first[1];
  yOut[0] = second[0] + yIn[0];
  yOut[1] = second[1] + yIn[1];                  


  // initialization
  // const unsigned int nbps2BPrcsed = outSizeMax - 1;
  std::vector<T> yFwdDiff(outSizeMax - 2);    // yFwdDiff.size() == outSizeMax - 2 in worst case, 
  std::vector<T> xFwdDiff(outSizeMax - 2);    
  //yFwdDiff.resize(nbps2BPrcsed);


  // prepare data container to call _iterative_align for each component later
  std::vector<unsigned int> ptsIndex(2); //ptsIndex.resize(2);
  unsigned int & nbpsProcsedPwl1 = ptsIndex[0];
  unsigned int & nbpsProcsedPwl2 = ptsIndex[1];
  std::vector<T> ptsCurrent(3); //ptsCurrent.resize(3);      
  T  &   xBkPwl = ptsCurrent[0];
  T  &  yBkPwl1 = ptsCurrent[1];
  T  &  yBkPwl2 = ptsCurrent[2];

  nbpsProcsedPwl1 = 2;
  nbpsProcsedPwl2 = 2;            
   xBkPwl = first[1];     // underEstmtr[i].first[1] = overEstmtr[i].first[1]      
  yBkPwl1 = second[1];
  yBkPwl2 = yIn[1];
  T varSumPtValLast = yOut[1];
  //long double varSumPtVarLast = output.first[1];
  T varSum2bSet;
  bool toMerge = false; 
   
  unsigned int nbpsAdded = 0;   // to count the number of breakpoints that have been added
  T step(0.),step2(0.),step3(0.),step4(0.),xBkpDiff(0.);
  while(nbpsProcsedPwl1 <= nBptPwl1 - 1 && nbpsProcsedPwl2 <= nBptPwl2 - 1){
    //  Note that the counter nbpsProcsedPwl1 and nbpsProcsedPwl2 will be updated in _iterative_interpolate
    //  Ideally, nbpsProcsedPwl1 == nBptPwl1 - 1 at the same iteration as nbpsProcsedPwl2 == nBptPwl2 - 1,
    //  as in each time of iteration we look at both the next Pwl1 breakpoint and the next Pwl2 breakpoint,
    //  and we only add the one with less x coordinate.  
    //  However, there can be cases when there are more than two breakpoints very close to the right endpoint
    //  To address such cases, we implement the functions _smoothening_xFilter and _smoothening_endpointFixer

    xBkpDiff = first[nbpsProcsedPwl1] - xIn[nbpsProcsedPwl2];
    
    if(std::fabs(xBkpDiff) > T(MC__UPWLE_COMPUTATION_TOL)){
      if (xBkpDiff > T(0.)){ // next x coordinate is of Pwl2 
        //if (_xFwdDiff.second) step = _xFwdDiff.first[nbpsProcsedPwl1 - 2];  // note that _yFwdDiff does not store its average
        //else step =  first[nbpsProcsedPwl1] - first[nbpsProcsedPwl1 - 1];
        step =  first[nbpsProcsedPwl1] - xBkPwl;        
        _iterative_interpolate(xBkpDiff, step, yBkPwl1, second[nbpsProcsedPwl1]);
        yBkPwl2 = yIn[nbpsProcsedPwl2];
        xBkPwl  = xIn[nbpsProcsedPwl2];
        nbpsProcsedPwl2 ++ ;
        // std::cout << "yBkPwl1: " << yBkPwl1 << " , "
        //           << "yBkPwl2: " << yBkPwl2 << " , "
        //           << "xBkPwl: "  << xBkPwl +50  << std::endl;
      }
      else{  // next x coordinate is of Pwl1  
        //if (xInFwdDiff.second) step = - xInFwdDiff.first[nbpsProcsedPwl2 - 2];
        //else step =  xIn[nbpsProcsedPwl2 - 1] - xIn[nbpsProcsedPwl2];
        step =  xBkPwl - xIn[nbpsProcsedPwl2];              
        _iterative_interpolate(xBkpDiff, step, yBkPwl2,    yIn[nbpsProcsedPwl2]);
        yBkPwl1 = second[nbpsProcsedPwl1];
        xBkPwl  = first[nbpsProcsedPwl1];   
        nbpsProcsedPwl1 ++ ;   
        // std::cout << "yBkPwl1: " << yBkPwl1 << " , "
        //           << "yBkPwl2: " << yBkPwl2 << " , "
        //           << "xBkPwl: "  << xBkPwl +50 << std::endl;
      }
    }
    else{ 
      if (std::fabs(xBkpDiff) < T(1e-15)){
        // At this time we view the two breakpoints are overlapped
        xBkPwl  =  first[nbpsProcsedPwl1];
        yBkPwl1 = second[nbpsProcsedPwl1];      
        yBkPwl2 =    yIn[nbpsProcsedPwl2];           
      }
      else{  //  here we need to condense the current two points, as they are too close 
        //std::cout << "in the loop " << xBkpDiff << " " << nbpsProcsedPwl2 << std::endl;
        if (xBkpDiff > T(0.)){ // next x coordinate is of Pwl2
          // check if we have reached the right endpoint of PWL1
          if( nbpsProcsedPwl1 == (nBptPwl1 - 1) ){
            std::vector<T> BkpsEnd(2);

            if(_smoothening_endpointFixer(ptsCurrent,BkpsEnd,yIn,first[nbpsProcsedPwl1],second[nbpsProcsedPwl1], xIn[nbpsProcsedPwl2], nbpsProcsedPwl2, nBptPwl2, _isUnder )){
              //  (std::vector<T> & data, std::vector<T> & endpoint, const std::vector<T> & yIn, 
              //  const T & xEndpoint, const T & yEndpoint, const T & BPtVar,
              //  const unsigned int nbpsProcsed, const unsigned int nBptPwl2, bool isUnder )
              varSum2bSet = yBkPwl1 + yBkPwl2;
              yFwdDiff[nbpsAdded] = varSum2bSet - varSumPtValLast;   
              _iterative_smoothening_yFilter(xOut,yOut,yFwdDiff,_isUnder,varSumPtValLast,varSum2bSet,xBkPwl,nbpsAdded,toMerge);                
            }

            xBkPwl = first[nbpsProcsedPwl1];
            varSum2bSet = BkpsEnd[0] + BkpsEnd[1];
            yFwdDiff[nbpsAdded] = varSum2bSet - varSumPtValLast;            
            _iterative_smoothening_yFilter(xOut,yOut,yFwdDiff,_isUnder,varSumPtValLast,varSum2bSet,xBkPwl,nbpsAdded,toMerge);                
            break;
          }

          // as the right endpoint of PWL1 has not yet attained, we can perform the filltering
          if (_yFwdDiff.second){  // note that _yFwdDiff does not store its average
            step  = -_yFwdDiff.first[nbpsProcsedPwl1 - 2];
            step2 =  _yFwdDiff.first[nbpsProcsedPwl1 - 1];   
          }
          else{   
            //step  = yBkPwl1 - second[nbpsProcsedPwl1];      
            step  = second[nbpsProcsedPwl1 - 1] - second[nbpsProcsedPwl1];           
            step2 = second[nbpsProcsedPwl1 + 1] - second[nbpsProcsedPwl1]; 
          }
          if (_xFwdDiff.second){  
            step3 =  _xFwdDiff.first[nbpsProcsedPwl1 - 2]; 
            step4 = -_xFwdDiff.first[nbpsProcsedPwl1 - 1];           
          }           
          else{
            //step3 =  first[nbpsProcsedPwl1] - xBkPwl; 
            step3 =  first[nbpsProcsedPwl1] - first[nbpsProcsedPwl1 - 1];                    
            step4 =  first[nbpsProcsedPwl1] - first[nbpsProcsedPwl1 + 1];            
          }
          yBkPwl1 = second[nbpsProcsedPwl1];
          _iterative_smoothening_xFilter(xBkpDiff,step,step2,yBkPwl1,step3,step4,_isUnder); 
          xBkPwl  = xIn[nbpsProcsedPwl2];
          yBkPwl2 = yIn[nbpsProcsedPwl2];

          // Note that there can be breakpoints of Pwl2 lying between first[nbpsProcsedPwl1] - xIn[nbpsProcsedPwl2]
          unsigned int i = nbpsProcsedPwl2 + 1;
          while(xIn[i] < first[nbpsProcsedPwl1]){
            yBkPwl2 = _isUnder? std::min(yBkPwl2,yIn[i]):std::max(yBkPwl2,yIn[i]);
            i++;
          }
          i --;
          if(i != nbpsProcsedPwl2){
            varSum2bSet = yBkPwl1 + yBkPwl2;
            yFwdDiff[nbpsAdded] = varSum2bSet - varSumPtValLast;
            _iterative_smoothening_yFilter(xOut,yOut,yFwdDiff,_isUnder,varSumPtValLast,varSum2bSet,xBkPwl,nbpsAdded,toMerge);          

            xBkPwl  = xIn[i];
            xBkpDiff = first[nbpsProcsedPwl1] - xIn[i];
            yBkPwl1 = second[nbpsProcsedPwl1];
            _iterative_smoothening_xFilter(xBkpDiff,step,step2,yBkPwl1,step3,step4,_isUnder); 
            nbpsProcsedPwl2 = i;
          }
      

        }
        else{  // next x coordinate is of Pwl1 

          // check if we have reached the right endpoint of PWL2  
          if(nbpsProcsedPwl2 == (nBptPwl2 - 1)){

            std::vector<T> BkpsEnd(2);
            // Note that the function _smoothening_endpointFixer assumes that we reaches the endpoint of PWL2,
            // so a simple swap is needed
            std::swap(yBkPwl1,yBkPwl2);
            if(_smoothening_endpointFixer(ptsCurrent,BkpsEnd,second,xIn[nbpsProcsedPwl2],yIn[nbpsProcsedPwl2], first[nbpsProcsedPwl1], nbpsProcsedPwl1, nBptPwl1, _isUnder )){
              //  (std::vector<T> & data, std::vector<T> & endpoint, const std::vector<T> & yIn, 
              //  const T & xEndpoint, const T & yEndpoint, const T & BPtVar,
              //  const unsigned int nbpsProcsed, const unsigned int nBptPwl2, bool isUnder )
              varSum2bSet = yBkPwl1 + yBkPwl2; // because of the commutative law of addition, we do not need to swap back
              yFwdDiff[nbpsAdded] = varSum2bSet - varSumPtValLast;   
              _iterative_smoothening_yFilter(xOut,yOut,yFwdDiff,_isUnder,varSumPtValLast,varSum2bSet,xBkPwl,nbpsAdded,toMerge);                
            }

            xBkPwl = xIn[nbpsProcsedPwl2];
            varSum2bSet = BkpsEnd[0] + BkpsEnd[1];  // because of the commutative law of addition, we do not need to swap back
            yFwdDiff[nbpsAdded] = varSum2bSet - varSumPtValLast;            
            _iterative_smoothening_yFilter(xOut,yOut,yFwdDiff,_isUnder,varSumPtValLast,varSum2bSet,xBkPwl,nbpsAdded,toMerge);                
            break;
          
          }

          if (yInFwdDiff.second){  // note that _yFwdDiff does not store its average
            step  = yInFwdDiff.first[nbpsProcsedPwl2 - 2];
            step2 = yInFwdDiff.first[nbpsProcsedPwl2 - 1];   
          }
          else{   
            //step  = yBkPwl1 - second[nbpsProcsedPwl2];      
            step  = yIn[nbpsProcsedPwl2] - yIn[nbpsProcsedPwl2 - 1];           
            step2 = yIn[nbpsProcsedPwl2 + 1] - yIn[nbpsProcsedPwl2]; 
          }
          if (xInFwdDiff.second){  
            step3 =  xInFwdDiff.first[nbpsProcsedPwl2 - 2]; 
            step4 =  xInFwdDiff.first[nbpsProcsedPwl2 - 1];           
          }           
          else{
            //step3 =  first[nbpsProcsedPwl2] - xBkPwl; 
            step3 =  xIn[nbpsProcsedPwl2] - xIn[nbpsProcsedPwl2 - 1];                    
            step4 =  xIn[nbpsProcsedPwl2 + 1] - xIn[nbpsProcsedPwl2];            
          }

          // step  = yIn[nbpsProcsedPwl2] - yBkPwl2;
          // step3 = xIn[nbpsProcsedPwl2] - xBkPwl;
          // step2 = yIn[nbpsProcsedPwl2 + 1] - yIn[nbpsProcsedPwl2]; 
          // step4 = xIn[nbpsProcsedPwl2 + 1] - xIn[nbpsProcsedPwl2];

          yBkPwl2 = yIn[nbpsProcsedPwl2];
          _iterative_smoothening_xFilter(xBkpDiff,step,step2,yBkPwl2,step3,step4,_isUnder); 
          yBkPwl1 = second[nbpsProcsedPwl1];
          xBkPwl  =  first[nbpsProcsedPwl1];    

          // Note that there can be breakpoints of Pwl1 lying between first[nbpsProcsedPwl1] - xIn[nbpsProcsedPwl2]
          unsigned int i = nbpsProcsedPwl1 + 1;
          while(first[i] < xIn[nbpsProcsedPwl2]){
            yBkPwl1 = _isUnder? std::min(yBkPwl1,second[i]):std::max(yBkPwl1,second[i]);
            i++;
          }
          i --;
          if(i != nbpsProcsedPwl1){
            varSum2bSet = yBkPwl1 + yBkPwl2;
            yFwdDiff[nbpsAdded] = varSum2bSet - varSumPtValLast;
            _iterative_smoothening_yFilter(xOut,yOut,yFwdDiff,_isUnder,varSumPtValLast,varSum2bSet,xBkPwl,nbpsAdded,toMerge);          

            xBkPwl  = first[i];
            xBkpDiff = first[i] - xIn[nbpsProcsedPwl2];
            yBkPwl2 = yIn[nbpsProcsedPwl2];
            _iterative_smoothening_xFilter(xBkpDiff,step,step2,yBkPwl2,step3,step4,_isUnder);  
            nbpsProcsedPwl1 = i;
          }

        }
      }
      nbpsProcsedPwl1 ++;   
      nbpsProcsedPwl2 ++;        
    }

    varSum2bSet = yBkPwl1 + yBkPwl2;

  //  long double yFwdDiff = varSum2bSet - varSumPtValLast;
    yFwdDiff[nbpsAdded] = varSum2bSet - varSumPtValLast;
    _iterative_smoothening_yFilter(xOut,yOut,yFwdDiff,_isUnder,varSumPtValLast,varSum2bSet,xBkPwl,nbpsAdded,toMerge);      
  }



  if(nbpsAdded <= 1){ // in this case the resultant estimator is degenarate to constant 
                      // As both pwl1 and pwl2 have at least one inner breakpoint, 
                      // nbpsAdded == 1 if and only if the resultant estimator is constant whose
                      // inner points are filtered in _iterative_smoothening_yFilter
    //std::cout << "nbpsAdded: " << nbpsAdded << std::endl;
#ifdef MC__UPWLE_DEBUG_TRACE
  std::cout <<"        nbpsAdded = " << nbpsAdded << " >= 1" << std::endl;
  std::cout <<"        xOut.size() " << xOut.size() << std::endl;  
  std::cout <<"        yOut.size() " << yOut.size() << std::endl;    
#endif         
    yOut[0] += yOut[1];
    yOut.resize(1);
    second.swap(yOut);      
    // xOut[1] += xOut[0];    
    // xOut[0] += xOut[2];
    // std::swap(xOut)
    // xOut.resize(2);
    // first.swap(xOut); 
    first[0] =  xOut[0] + xOut[1];
    first[1] =  xOut[0] + xOut[2];   
    first.resize(2); 
#ifdef MC__UPWLE_DEBUG_TRACE
  std::cout <<"      out altMerge 1" << std::endl;
#endif     
    return ;   
  }
  //else if( nbpsAdded <= _nbps + 1){    // _nbps is the number of inner points inside the interval-bnd of the var, nbpsAdded includes the right endpoint
  xOut.resize(nbpsAdded+2);       
  yOut.resize(nbpsAdded+2);   
  first.swap(xOut);
  second.swap(yOut);   

  _yFwdDiff.first.swap(yFwdDiff);
  _yFwdDiff.second = true;
  _xFwdDiff.first.resize(0);
  _xFwdDiff.second = false;   

#ifdef MC__UPWLE_DEBUG_TRACE
  std::cout <<"      out altMerge 2" << std::endl;
#endif    
  // if(second.size() > nbpsMax){
  //   _yFwdDiff.swap(yFwdDiff);
  //   for(unsigned j = 1; j < xOut.size() - 1 ; j++){
  //     xFwdDiff[j-1] = xOut[j+1] - xOut[j]; 
  //   }
  //   _xFwdDiff.swap(xFwdDiff);
  //   _isFwdDiff = true;
  // }
   
  //}

  // _yFwdDiff.resize(nbpsAdded);
  // _condense_by_relax(output,_yFwdDiff);
  // output.first.resize(_nbps+3);
  // output.second.resize(_nbps+3);  
  // return std::move(output);
}

template <typename T>
inline
void UnivarPWLE<T>::_relu()
{
  //std::cout << "in _relu" << std::endl;
#ifdef MC__UPWLE_DEBUG  
  assert( !first.empty() );  // it shouldn't be a constant
  assert(  second.size() > 2 );  // it shouldn't be linear
#endif

  const unsigned int nBptPwl = second.size();
  const T negMid = (-second[0]);

  // _bnd.second = true;
  // _bnd.first.first.second = false;
  // _bnd.first.first.first = 0.;
  // _bnd.first.second.first = ub;  
  // _bnd.first.second.second = true;  


  //unsigned int allZeroCnt = 1;
  bool existsZero = false;
  bool startPositive = true;
  std::vector<T> xOut(2*nBptPwl-1);
  std::vector<T> yOut(2*nBptPwl-1);
  xOut[0] =  first[0];               // <- load the mean of x
  xOut[1] =  first[1];               
  yOut[0] = second[0];               // <- load the mean of y (it is not the mean any more after computation)
  if (negMid >= second[1]){
    yOut[1] = negMid;
    existsZero = true;
    startPositive = false;
    //allZeroCnt ++;
  }
  else{
    yOut[1] = second[1];
  }

  
  unsigned int state(10),state0(2);
  unsigned int state1 = 2*(negMid >= second[1]? 0:1);
  T yLast2BItplt(second[1]),xLast2BItplt(first[1]);  
  unsigned int  bptsAdded = 2;
   

  for(unsigned int i = 2; i < nBptPwl; i++){
    
    state0 = (negMid >= second[i]? 0:1);    // 0: negative, 1 positive
    state = state0 + state1; 
    //std::cout << "_in _relu _in main loop? " << state << std::endl;    
    switch(state){
      // four states 00, 11, 01, 10,
      case 0:
        existsZero = true;//allZeroCnt ++;
        break;
      case 2:
        //_iterative_smoothening_yFilter();
        existsZero = true;

        // xOut[bptsAdded] = first[i] - (second[i] - negMid)/(second[i] - yLast2BItplt)*(first[i] - xLast2BItplt);
        // if(std::fabs(xOut[bptsAdded] - xOut[bptsAdded-1]) < T(MC__UPWLE_COMPUTATION_TOL)){
        //   xOut[bptsAdded-1] = _isUnder?std::max(xOut[bptsAdded],xOut[bptsAdded-1]):std::min(xOut[bptsAdded],xOut[bptsAdded-1]);
        //   yOut[bptsAdded-1] = _isUnder?std::min(negMid,yOut[bptsAdded-1]):std::max(negMid,yOut[bptsAdded-1]);
        // }
        // else{
        //   yOut[bptsAdded] = negMid;
        //   bptsAdded ++;
        // }

        if(yLast2BItplt - negMid < T(MC__UPWLE_COMPUTATION_TOL)){
          if(_isUnder){
            yOut[bptsAdded - 1] = negMid; // this is fine for under estimator             
          }
          // Note that for over estimator there should not be a numerical error breaking the validity
          // this is because the state is recogised as a zero vertex has been added (but it is not)
        }
        else{
          xOut[bptsAdded] = first[i] - (second[i] - negMid)/(second[i] - yLast2BItplt)*(first[i] - xLast2BItplt);
          yOut[bptsAdded] = negMid;
          bptsAdded ++;
        }

        break;
      case 1:
        //_iterative_smoothening_yFilter();                
        if(yLast2BItplt - negMid == T(0.)){
        // in this case when last point is 0, there are two possible cases: (1) it has been added, and (2) it hasn't
        // if it is added, then we add the current point then break
          if(xLast2BItplt == xOut[bptsAdded-1]){
            yOut[bptsAdded] = second[i];
            xOut[bptsAdded] =  first[i];
            bptsAdded ++;
            break;    
          }        
          // otherwise, we add the previous point
          xOut[bptsAdded] = xLast2BItplt;
          yOut[bptsAdded] = negMid;
          bptsAdded ++;
        }
        else if( negMid - yLast2BItplt < T(MC__UPWLE_COMPUTATION_TOL)){
          if(!_isUnder){
            if(xLast2BItplt > xOut[bptsAdded-1]){
              yOut[bptsAdded] = negMid; // this is fine for over estimator             
              xOut[bptsAdded] = xLast2BItplt;  // since yLast2BItplt < 0, the last point is not yet added
              bptsAdded ++;
            }
          }
          else{
          // Note that for under estimator there introduces numerical errors
            //xOut[bptsAdded - 1] = first[i] - (second[i] - negMid)/(second[i] - yLast2BItplt)*(first[i] - xLast2BItplt);
            xOut[bptsAdded] = first[i] - (second[i] - negMid)/(second[i] - yLast2BItplt)*(first[i] - xLast2BItplt);
            if( std::fabs(xOut[bptsAdded] - xOut[bptsAdded-1]) < T(MC__UPWLE_COMPUTATION_TOL)){              
//              xOut[bptsAdded-1] = (bptsAdded==2)?first[1]:std::max(xOut[bptsAdded],xOut[bptsAdded-1]);
//              yOut[bptsAdded-1] = negMid; 
             ;
              // The following is to ensure the validity of underestimator
              // xOut[bptsAdded] = std::nextafter(xOut[bptsAdded-1],first.back());
              // yOut[bptsAdded] = negMid;
              // bptsAdded ++;
            }
            else{
              yOut[bptsAdded] = negMid;
              bptsAdded ++;
            }

          }  
        }
        else{
            xOut[bptsAdded] = first[i] - (second[i] - negMid)/(second[i] - yLast2BItplt)*(first[i] - xLast2BItplt);
            yOut[bptsAdded] = negMid;
            bptsAdded ++;
        }        
        
      case 3:
        //_iterative_smoothening_yFilter();
        yOut[bptsAdded] = second[i];
        xOut[bptsAdded] =  first[i];
        bptsAdded ++;
        break;  
      default: 
        break;
    }
    yLast2BItplt = second[i];
    xLast2BItplt =  first[i];
    state1 = state0*2;
  } 

  if (yLast2BItplt <= negMid){
    xOut[bptsAdded] = xLast2BItplt;
    if((!startPositive) && bptsAdded <= 2){
      second[0] = T(0.);
      second.resize(1);
      first.back() += first[0];
      first[0] += first[1];
      first[1] = first.back();
      first.resize(2);
      // _lbnd.second = true;
      // _lbnd.first = T(0.);
      // _ubnd.second = true;
      // _ubnd.first = T(0.);

      // _xFwdDiff.second = false;  // it should be noted that the difference can be updated in the main loop
      // _xFwdDiff.first.resize(0);
      // _yFwdDiff.second = false;   
      // _yFwdDiff.first.resize(0);     
      return ;
    }
    else{
      existsZero = true;
      yOut[bptsAdded] = negMid;
      bptsAdded ++ ;
    }
  }
 

  // if(allZeroCnt == nBptPwl){
  //   second[0] = 0;
  //   second.resize(1);
  //   first[2] = first[nBptPwl-1];
  //   first.resize(3);
  //   _bnd.first.second.second = true;
  //   _bnd.first.second.first = 0.;  
  //   return ;
  // }

  if(existsZero){
    //std::cout << "existsZero" << std::endl;
//    _bnd.second = false;
    _lbnd.second = true;
    _lbnd.first = T(0.);
    _xFwdDiff.second = false;  // it should be noted that the difference can be updated in the main loop
    _xFwdDiff.first.resize(0);
    _yFwdDiff.second = false;   
    _yFwdDiff.first.resize(0);     

    if(std::fabs(xOut[bptsAdded-1] - xOut[bptsAdded-2]) < T(MC__UPWLE_COMPUTATION_TOL)){
      xOut[bptsAdded-2] = first.back();
      yOut[bptsAdded-2] = _isUnder?std::min(yOut[bptsAdded-2],yOut[bptsAdded-1]):std::max(yOut[bptsAdded-2],yOut[bptsAdded-1]);
      bptsAdded --;
    }


    if(bptsAdded == 3){
      first[0] = xOut[0] + xOut[1];
      first[1] = xOut[0] + xOut[2];
      first.resize(2);
      second[0] = yOut[0] + yOut[1];
      second[1] = yOut[0] + yOut[2];
      second.resize(2);

    }
    else{
      xOut.resize(bptsAdded);
      first.swap(xOut);
      yOut.resize(bptsAdded);
      second.swap(yOut);

    }
    return ;



  }
  else{
    return ;  // it is usally not possible as we have already handled this case outside
              // only one case we may have this: we did not computed the lb of the PWLE in advance  
  }
  
  //std::cout << "in _relu assign:" << bptsAdded << std::endl;


  // if(bptsAdded > nbpsMax){
  //   _condense_by_relax(nbpsMax);
  // }

}



template <typename T>
inline
UnivarPWLE<T> operator+
( UnivarPWLE<T> const& var )
{
  return var;
}

template <typename T>
inline
UnivarPWLE<T> operator+
( UnivarPWLE<T> && varIn )
{
  UnivarPWLE<T> var( std::move(varIn) );  
  return var;
}


template <typename T>
inline
UnivarPWLE<T>& UnivarPWLE<T>::operator+=
( double const& cst )
{

#ifdef MC__UPWLE_DEBUG_TRACE  
  std::cout << "UPWLE operator += cst " << std::endl;
#endif
  if( cst == 0. ){
    return *this;
  }
  if( first.empty() ){
    _cst += cst;
    return *this;
  }
  
#ifdef MC__UPWLE_DEBUG_TRACE  
  if(first.size() == 2){
      std::cout << "first " << first[0] << " : " << first[1]<<  std::endl;
  }
  else if(first.size() > 2){
    std::cout << "first " << std::endl;
    for(unsigned int jj = 1; jj < first.size(); jj++ ){
      std::cout << first[0] + first[jj] << "    ";
    }
    std::cout << std::endl;
  }
#endif

  second[0] += cst;
  // std::cout << "second[0] += cst " << second[0] << std::endl;
  if( second.size() == 2){
    second[1] += cst;
    // std::cout << "second[1] += cst " << second[1] << std::endl;
    return *this;
  }  

  if( _lbnd.second ){
    _lbnd.first += cst;
  }

  if( _ubnd.second ){
    _ubnd.first += cst;
  }
  return *this;
}

// Reserved for adding constant intervals/other sets
       // template <typename T>
       // inline
       // UnivarPWLE<T>& UnivarPWLE<T>::operator+=
       // ( T const& cstInterval )
       // {
       
       //   if( !_mod ){
       //     _cst += cst;
       //     return *this;
       //   }
       //   if( !_ndep )
       //     throw typename ISModel<T>::Exceptions( ISModel<T>::Exceptions::INTERN );
       
       //   const T offset = cst / (double)_ndep;
       
       //   for( unsigned int i=0; i<_nvar; i++ ){
       //     if( _mat[i].empty() ) continue;
       //     for( unsigned int j=0; j<_ndiv; j++ )
       //       _mat[i][j] += offset;
       //   }
       //   if( _bnd.second ) _bnd.first += cst;
       
       //   return *this;
       // }

template <typename T>
inline
UnivarPWLE<T>& UnivarPWLE<T>::operator+=
( std::pair<T,T> const& endpoints )
{

#ifdef MC__UPWLE_DEBUG_TRACE 
//std::cout << "operator+= endpoints " << std::endl;
#endif

  const unsigned int nBptPwl = second.size();
  if( nBptPwl == 1){
    second.resize(2);
    second[1] = second[0] + endpoints.second;
    second[0] += endpoints.first;
    return *this;
    //std::cout << "nBptPwl == 2 " << std::endl;  
  }
  else if(nBptPwl == 2){
    second[0] += endpoints.first;
    second[1] += endpoints.second;
    if(std::fabs(second[1] - second[2]) < T(MC__UPWLE_COMPUTATION_TOL) ){
      _cst = _isUnder? std::min(second[1],second[2]):std::max(second[1],second[2]); 
      first.resize(0);
      second.resize(0);
    }   
    return *this;
  } 
#ifdef MC__UPWLE_DEBUG_TRACE     
  //std::cout << "nBptPwl > 2 " << std::endl;             

    // std::cout << "first "  << std::endl;    
    // for(unsigned int jj = 1; jj < first.size(); jj++){
    //   std::cout << first[jj] + first[0]<< "  ";     
    // }
    // std::cout << std::endl;
    // std::cout << "second "  << std::endl;    
    // for(unsigned int jj = 1; jj < second.size(); jj++){
    //   std::cout << second[jj] + second[0]<< "  ";     
    // }
    // std::cout << std::endl; 
#endif

#ifdef MC__UPWLE_DEBUG
    for(unsigned int j = 0; j < second.size(); j++){
      assert(second[j]< 1e7);;
    }     
#endif

  const T mid2 = 0.5*(endpoints.first + endpoints.second);
  second[0] += mid2;    
  const T varLftVal = endpoints.first - mid2;
  second[1] += varLftVal;
  second.back() += (endpoints.second - mid2);    
  const T ratio = (endpoints.second - endpoints.first)/(first.back() - first[1]);
  for(unsigned int j = 2; j < nBptPwl - 1; j++){
    second[j] += (varLftVal + (first[j]-first[1])*ratio);//_interpolate_at_ind(first, width, diff, j);
  }

#ifdef MC__UPWLE_DEBUG_TRACE     
    // std::cout << "output first "  << std::endl;    
    // for(unsigned int jj = 1; jj < first.size(); jj++){
    //   std::cout << first[jj] + first[0]<< "  ";     
    // }
    // std::cout << std::endl;
    // std::cout << "output second "  << std::endl;    
    // for(unsigned int jj = 1; jj < second.size(); jj++){
    //   std::cout << second[jj] + second[0]<< "  ";     
    // }
    // std::cout << std::endl;
#endif

#ifdef MC__UPWLE_DEBUG
    for(unsigned int j = 0; j < second.size(); j++){
      if(second[j]> 1e7){ 
        std::cout << endpoints.first << "  " << endpoints.second << std::endl;
        std::cout <<  first[1] << "  " << first.back() << std::endl;
      }
      assert(second[j]< 1e7);
    } 
#endif

  _lbnd.second = false;
  _ubnd.second = false; 

  //if(_xFwdDiff.second) {_xFwdDiff.first.resize(0); _xFwdDiff.second = false};
  if(_yFwdDiff.second) {_yFwdDiff.first.resize(0); _yFwdDiff.second = false;}  
    

  return *this;
}


template <typename T>
inline
UnivarPWLE<T>& UnivarPWLE<T>::operator+=
( UnivarPWLE<T> const& var )
{
#ifdef MC__UPWLE_DEBUG_TRACE  
  std::cout << "***********************" << std::endl;
  std::cout << "in PWLE += var "  << std::endl;
  std::cout << (_isUnder?"    under ":"    over ") << std::endl;
  std::cout << "***********************" << std::endl; 
#endif

  // process the constant case
  if( first.empty() && var.first.empty() ){
    _cst += var._cst;
    return *this;
  }
  else if( first.empty() ){
    T copy_cst = _cst;
    *this = var;
    *this += copy_cst;
    return *this;
  }
  else if( var.first.empty() ){
    *this += var._cst;
    return *this;
  }

  // if both !first.empty() and !var.first.empty(), then _isUnder should be initialized
  if( _isUnder != var._isUnder ){
     std::cout << " Addition cannot be performed between an over- and an under-estimator" << std::endl;
     throw typename UnivarPWLE<T>::Exceptions( UnivarPWLE<T>::Exceptions::ESTMATCH );
  }



#ifdef MC__UPWLE_DEBUG_TRACE  
  std::cout << "    non-constant case "  << std::endl;
  if(var.first.size()>2){    
    std::cout << "    var first size "  << var.first.size() << std::endl;    
    for(unsigned int jj = 1; jj < var.first.size(); jj++){
      std::cout << "    jj = " << jj << "  "<< var.first[jj] + var.first[0]<< "  ";     
    }
    std::cout << std::endl;
#ifdef MC__UPWLE_DEBUG
    assert(var.first.back() - var.first[1] >=T(0.));
#endif
  }
  else if(var.first.size()==2){
    std::cout << "    var first size "  << var.first.size() << std::endl;    
    for(unsigned int jj = 0; jj < var.first.size(); jj++){
      std::cout << "    jj = " << jj << "  "<< var.first[jj] << "  ";     
    }
    std::cout << std::endl;    

    std::cout << "    first size"  << first.size() << std::endl;    
    for(unsigned int jj = 0; jj < first.size(); jj++){
      std::cout << "    jj = " << jj << "  "<< first[jj] << "  ";     
    }
    std::cout << std::endl;    
#ifdef MC__UPWLE_DEBUG
    assert(var.first.back() - var.first.front() >=T(0.));
#endif
  }
  else{
#ifdef MC__UPWLE_DEBUG
    std::cout << var.first[0] << std::endl;
    assert(false);
#endif
  }

  if(first.size()>2){    
    
    // std::cout << "first "  << first.size() << std::endl;    
    // for(unsigned int jj = 1; jj < first.size(); jj++){
    //   std::cout << "jj = " << jj << "  "<< first[jj] + first[0]<< "  ";     
    // }
    // std::cout << std::endl;
#ifdef MC__UPWLE_DEBUG
    assert(first.back() - first[1] >=T(0.));
#endif
  }
  else if(first.size()==2){

    // std::cout << "first "  << first.size() << std::endl;    
    // for(unsigned int jj = 0; jj < first.size(); jj++){
    //   std::cout << "jj = " << jj << "  "<< first[jj] << "  ";     
    // }
    // std::cout << std::endl;    

    //  std::cout << "var first "  << var.first.size() << std::endl;    
    // for(unsigned int jj = 0; jj < var.first.size(); jj++){
    //   std::cout << "jj = " << jj << "  "<< var.first[jj] << "  ";     
    // }
    // std::cout << std::endl;    
#ifdef MC__UPWLE_DEBUG
    assert(first.back() - first.front() >=T(0.));
#endif
  }
  else{
    
#ifdef MC__UPWLE_DEBUG
    std::cout << first[0] << std::endl;
    assert(false);
#endif
  }
#endif



  if(var.second.size() == 1){
    *this += std::make_pair(var.second[0],var.second[0]);
    return *this;
  }
  else if(var.second.size() == 2){
    *this += std::make_pair(var.second[0],var.second[1]);
    return *this;
  }

  if( second.size() <= 2  ){
    if(second.size() == 1){
      auto copy_cst = second[0];
      *this = var;
      *this += copy_cst;
      return *this;
    }    
    auto copy_cst = std::make_pair(second[0],second[1]);
    *this = var;
    *this += copy_cst;
    return *this;
  } 

#ifdef MC__UPWLE_DEBUG
  assert(second.size() != 3);
  assert(var.second.size() != 3);
#endif

#ifdef MC__UPWLE_DEBUG_TRACE
  std::cout <<"    before altMerge " << std::endl;
  std::cout <<"    first.size() " << first.size() << std::endl;
  std::cout <<"    second.size() " << second.size() << std::endl;  
  std::cout <<"    var first.size() " << var.first.size() << std::endl;
  std::cout <<"    var second.size() " << var.second.size() << std::endl;  
#endif

  _add_and_altMerge(var.first,var.second,var._xFwdDiff,var._yFwdDiff);

#ifdef MC__UPWLE_DEBUG_TRACE
  std::cout <<"    after altMerge " << std::endl;
  std::cout <<"    first.size() " << first.size() << std::endl;
  std::cout <<"    second.size() " << second.size() << std::endl;  
#endif
      
  _lbnd.second = false;
  _ubnd.second = false; 
  
  if(second.size() > nbpsMax){
    _condense_by_relax(nbpsMax);
  }
  else if(second.size() == 3){

#ifdef MC__UPWLE_DEBUG_TRACE
  std::cout <<"    merge > 3 " << std::endl;
#endif

    second[1] += second[0];
    second[0] += second[2];
    std::swap(second[0],second[1]);
    second.resize(2);

    first[1] += first[0]; 
    first[0] += first[2];
    std::swap(first[0],first[1]);
    first.resize(2);
#ifdef MC__UPWLE_DEBUG    
    assert(first[1] - first[0] > -1e-13);
#endif
  }  

  if(!_xFwdDiff.second) _xFwdDiff.first.resize(0);
  if(!_yFwdDiff.second) _yFwdDiff.first.resize(0);  

  
  return *this;

}



// template <typename T>
// inline
// UnivarPWLE<T>& UnivarPWLE<T>::operator+=
// ( UnivarPWLE<T> && var )
// {
//   //std::cout << "under?2 " << _isUnder  << std::endl;
//   if( first.empty() && var.first.empty() ){
//     _cst += var._cst;
//     return *this;
//   }
//   else if( first.empty() ){
//     T copy_cst = _cst;
//     *this = var;
//     *this += copy_cst;
//     return *this;
//   }
//   else if( var.first.empty() ){
//     *this += var._cst;
//     return *this;
//   }


//   if( _isUnder != var._isUnder ){
//      std::cout << "Addition cannot be performed between an over- and an under-estimator" << std::endl;
//      throw typename UnivarPWLE<T>::Exceptions( UnivarPWLE<T>::Exceptions::ESTMATCH );
//   }

//   const unsigned int nBptPwl1 = second.size();
//   const unsigned int nBptPwl2 = var.second.size();
//   if( nBptPwl2 == 1 ){
     
//     second[0] += var.second[0];
//     if(nBptPwl1 == 2) second[1] += var.second[0];
//     else if(nBptPwl1 > 2){
//       if(_lbnd.second) _lbnd.first += var.second[0];
//       if(_ubnd.second) _ubnd.first += var.second[0];       
//     }
//     return *this;
//   } 

//   if( nBptPwl1 == 1){
//     auto const tmpMid = second[0];
//     second = std::move(var.second); 
//     second[0] += tmpMid; 
//     if(nBptPwl2 == 2){
//       second[1] += tmpMid;
//     }
//     else{ // note that nBptPwl2 != 1
//       first = std::move(var.first);
//       if(var._lbnd.second){ _lbnd = std::make_pair(tmpMid + _lbnd.first,true); }
//       else { _lbnd = std::make_pair(0.,false); }
//       if(var._ubnd.second){ _ubnd = std::make_pair(tmpMid + _ubnd.first,true); }
//       else { _ubnd = std::make_pair(0.,false); }
//       _xFwdDiff = var._xFwdDiff.second? std::move(var._xFwdDiff):std::make_pair(std::vector<T>(0),false);
//       _yFwdDiff = var._yFwdDiff.second? std::move(var._yFwdDiff):std::make_pair(std::vector<T>(0),false);        
//     }  
//     return *this;
//   } 
// //  std::cout << "both nBptPwls != 1 " << std::endl;

  
//   if( nBptPwl1 == 2 && nBptPwl2 == 2 ){
//     //std::cout << "both nBptPwls == 2 " << std::endl;
//     second[0] += var.second[0];
//     second[1] += var.second[1];
//     if(std::fabs(second[1] - second[2]) < T(MC__UPWLE_COMPUTATION_TOL) ){
//       _cst = _isUnder? std::min(second[1],second[2]):std::max(second[1],second[2]); 
//       first.resize(0);
//       second.resize(0);
//     }
    
//     // _lbnd.second = false;
//     // _ubnd.second = false; 

//     // _yFwdDiff.second = false;
//     // _yFwdDiff.first.resize(0);
//     // _xFwdDiff.second = false;
//     // _xFwdDiff.first.resize(0);

//     for(unsigned int j = 0; j < second.size(); j++){
//       assert(second[j]< 1e7);;
//     } 

//     return *this;
//   } 
// //  std::cout << "nBptPwl2 = 3 " << std::endl;
//   if( nBptPwl2 == 2){
//     //std::cout << "nBptPwl2 == 2 " << std::endl;       
//     const T mid2 = 0.5*(var.second[0] + var.second[1]);
//     second[0] += mid2;    
//     const T varLftVal = var.second[0] - mid2;
//     second[1] += varLftVal;
//     //second[nBptPwl1 - 1] += var.second[1] - mid2;
//     second.back() += (var.second[1] - mid2);    
//     const T width = var.first[1]  - var.first[0];
//     const T diff  = var.second[1] - var.second[0];
//     for(unsigned int j = 2; j < nBptPwl1 - 1; j++){
//       second[j] += varLftVal + _interpolate_at_ind(first, width, diff, j);
//     }

//     // _bnd.second = false;
//     // _bnd.first.first.second = false; 
//     // _bnd.first.second.second = false;   


//     for(unsigned int j = 0; j < second.size(); j++){
//       assert(second[j]< 1e7);;
//     } 


//     _lbnd.second = false;
//     _ubnd.second = false; 

//     //if(_xFwdDiff.second) {_xFwdDiff.first.resize(0); _xFwdDiff.second = false};
//     if(_yFwdDiff.second) {_yFwdDiff.first.resize(0); _yFwdDiff.second = false;}  
    

//     return *this;
//   } 

// //  std::cout << "nBptPwl1 = 3 " << std::endl;  
//   if( nBptPwl1 == 2){
//     //std::cout << "nBptPwl1 == 2 " << std::endl;       
//     second.resize(nBptPwl2);
//     const T mid1  = 0.5*(second[0] + second[1]);
//     const T diff  = second[1] - second[0];    
//     second.back() = second[1] - mid1 + var.second.back();  
//     second[1] = second[0] - mid1;// + var.second[1];
//     second[0] = mid1 + var.second[0];
       

//     // yOut[0] = second[0] + var.second[0];
//     // yOut[1] = second[1] + var.second[1];
//     // yOut[nBptPwl2 - 1] = second[2] + var.second[nBptPwl2 - 1];    
//     const T width = first[1]  - first[0];

//     first = std::move(var.first);

//     for(unsigned int j = 2; j < nBptPwl2 - 1; j++){
//       second[j] = second[1] + _interpolate_at_ind(first, width, diff, j) + var.second[j];
//     }      
//     second[1] += + var.second[1];//     
//     // second.swap(yOut);
//     // _bnd.second = false;
//     // _bnd.first.first.second = false; 
//     // _bnd.first.second.second = false;    

//     // _yFwdDiff.second = false;
//     // _yFwdDiff.first.resize(0);

//       _lbnd = std::make_pair(0.,false); 
//       _ubnd = std::make_pair(0.,false);
//       _yFwdDiff = std::make_pair(std::vector<T>(0),false);  
//       _xFwdDiff = var._xFwdDiff.second? std::move(var._xFwdDiff):std::make_pair(std::vector<T>(0),false);        
      
//     //if(_xFwdDiff.second) {_xFwdDiff.first.resize(0); _xFwdDiff.second = false};
//     //if(_yFwdDiff.second) {_yFwdDiff.first.resize(0); _yFwdDiff.second = false;}  

//     return *this;
//   } 


// //  std::cout << "both nBptPwls > 3  " << std::endl;  
//   _add_and_altMerge(var.first,var.second,var._xFwdDiff,var._yFwdDiff);


//   //  std::cout << "after _add_and_altMerge" << std::endl;
//   //  auto temp = _bnd.first.first.first;
//   //  _B(true,false);
//   //   if(temp != _bnd.first.first.first) {
//   //     std::cout << "und bound mismatch" << std::endl;
//   //     throw typename UnivarPWLE<T>::Exceptions( UnivarPWLE<T>::Exceptions::UNDEF );
//   //   }
      
//   _lbnd.second = false;
//   _ubnd.second = false; 
  

//   if(second.size() > nbpsMax){
//     _condense_by_relax(nbpsMax);
//   }
//   else if(second.size() == 3){
//     second[1] += second[0];
//     second[0] += second[2];
//     std::swap(second[0],second[1]);
//     second.resize(2);

//     first[1] += first[0];
//     first[0] += first[2];
//     std::swap(first[0],first[1]);
//     first.resize(2);

//   }   
//   if(!_xFwdDiff.second) _xFwdDiff.first.resize(0);
//   if(!_yFwdDiff.second) _yFwdDiff.first.resize(0);  

//   // _bnd.second = false;
//   // _bnd.first.first.second = false; 
//   // _bnd.first.second.second = false;  

//   return *this;

// }                


template <typename T>
inline
UnivarPWLE<T> operator+
( UnivarPWLE<T> const& var1, UnivarPWLE<T> const& var2 )
{
  if( (var1.first.empty()) && (var2.first.empty()) ) 
    return( var1._cst + var2._cst );
  if( (var1.first.empty()) ){
    UnivarPWLE<T> var3( var1 );
    var3 += var2;
    return var3;
  }
  UnivarPWLE<T> var3( var2 );
    var3 += var1;
    return var3;
}


template <typename T>
inline
UnivarPWLE<T> operator+
( UnivarPWLE<T> const& var1, UnivarPWLE<T> && var2 )
{
  //std::cout << "+ 1move" << std::endl;
  if( (var1.first.empty()) && (var2.first.empty()) ) 
    return( var1._cst + var2._cst );
   
  if( (var2.first.empty()) ){
    UnivarPWLE<T> var3( std::move(var2) );
    var3 += var1;
    return var3;
  }
  UnivarPWLE<T> var3( var1 );
  var3 += var2;
  return var3;
}


template <typename T>
inline
UnivarPWLE<T> operator+
( UnivarPWLE<T> && var1, UnivarPWLE<T> const& var2 )
{
  //std::cout << "+ 2move" << std::endl;
  if( (var1.first.empty()) && (var2.first.empty())  ) 
    return( var1._cst + var2._cst );
   
  if( (var1.first.empty()) ){
    UnivarPWLE<T> var3( std::move(var1) );
    var3 += var2;
    return var3;
  }
  UnivarPWLE<T> var3( var2 );
  var3 += var1;
  return var3;
}


template <typename T>
inline
UnivarPWLE<T> operator+
( UnivarPWLE<T> && var1, UnivarPWLE<T> && var2 )
{
  //std::cout << "+ 3move" << std::endl;
  if( (var1.first.empty()) && (var2.first.empty()) ) 
    return( var1._cst + var2._cst );
  if( (var1.first.empty()) ){
    UnivarPWLE<T> var3( std::move(var1) );
    var3 += var2;
    return var3;
  }
  UnivarPWLE<T> var3( std::move(var2) );
    var3 += var1;
    return var3;
}



template <typename T>
inline
UnivarPWLE<T> operator+
( UnivarPWLE<T> const& var1, double const& cst2 )
{
  if( (var1.first.empty()) ) 
    return( var1._cst + cst2 );
  UnivarPWLE<T> var3( var1 );
  var3 += cst2;
  return var3;
}


template <typename T>
inline
UnivarPWLE<T> operator+
( UnivarPWLE<T> && var1, double const& cst2 )
{
  if( (var1.first.empty()) ) 
    return( var1._cst + cst2 );
  UnivarPWLE<T> var3( std::move(var1) );
  var3 += cst2;
  return var3;
}


template <typename T>
inline
UnivarPWLE<T> operator+
( double const& cst1, UnivarPWLE<T> const& var2 )
{
  if( (var2.first.empty()) ) 
    return( cst1 + var2._cst );
  UnivarPWLE<T> var3( var2 );
  var3 += cst1;
  return var3;
}

template <typename T>
inline
UnivarPWLE<T> operator+
( double const& cst1, UnivarPWLE<T> && var2 )
{
  if( (var2.first.empty()) ) 
    return( cst1 + var2._cst );
  UnivarPWLE<T> var3( std::move(var2) );
  var3 += cst1;
  return var3;
}


template <typename T> 
inline
UnivarPWLE<T> operator-
( UnivarPWLE<T> const& var )
{
  UnivarPWLE<T> var2( var );
  if( (var2.first.empty()) ){
    var2._cst *= -1;
    return var2;
  }
    
  for(auto & yBkp:var2.second){
    yBkp *= -1;
  }
    
  var2._isUnder = !(var2._isUnder);

  if(var2.second.size() > 2){
    // if(var2._yFwdDiff.second){
    //   // for(auto & yBkp:var2._yFwdDiff.first){
    //   //   yBkp *= -1;
    //   // }
    //   var2._yFwdDiff.second = false;
    //   var2._yFwdDiff.first.resize(0);
    // }
    std::swap(var2._lbnd,var2._ubnd);
    if(var2._lbnd.second) var2._lbnd.first *= -1;
    if(var2._ubnd.second) var2._ubnd.first *= -1;
    
  }

  //
  //auto&& bnd = var2._bnd;
  //if( bnd.second ) bnd.first *= -1;

  return var2;
}





template <typename T> 
inline
UnivarPWLE<T> operator-
( UnivarPWLE<T> && var )
{
  UnivarPWLE<T> var2( std::move(var) );
  if( (var2.first.empty()) ){
    var2._cst *= -1;
    return var2;
  }
    
  for(auto & yBkp:var2.second){
    yBkp *= -1;
  }
    
  var2._isUnder = !(var2._isUnder);


  if(var2.second.size() > 2){
    if(var2._yFwdDiff.second){
      for(auto & yBkp:var2._yFwdDiff.first){
        yBkp *= -1;
      }
    }
    std::swap(var2._lbnd,var2._ubnd);
    if(var2._lbnd.second) var2._lbnd.first *= -1;
    if(var2._ubnd.second) var2._ubnd.first *= -1;

  }

  return var2;
}





template <typename T>
inline
UnivarPWLE<T>& UnivarPWLE<T>::operator-=
( double const& cst )
{
  *this += -cst;
  return *this;
}

template <typename T>
inline
UnivarPWLE<T>& UnivarPWLE<T>::operator-=
( UnivarPWLE<T> const& var )
{

  // if( first.empty() && var.first.empty() ){
  //   _cst -= var._cst;
  //   return *this;
  // }
  // if( first.empty() ){
  //   T copy_cst = _cst;
  //   *this = -var;
  //   *this += copy_cst;
  //   return *this;
  // }
  // if( var.first.empty() ){
  //   *this -= var._cst;
  //   return *this;
  // }

  // if( _isUnder == var._isUnder )
  //    throw typename UnivarPWLE<T>::Exceptions( UnivarPWLE<T>::Exceptions::ESTMATCH );

  // _bnd.second = false;
  // _bnd.first.first.second = false;
  // _bnd.first.second.second = false;    

  // use slope information to tighten interval bounds

  *this += (-var);
  return *this;
}


template <typename T>
inline
UnivarPWLE<T> operator-
( UnivarPWLE<T> const& var1, UnivarPWLE<T> const& var2 )
{
  if(  var1.first.empty() && var2.first.empty()  ) 
    return( var1._cst - var2._cst );
  if(  (var1.first.empty()) ){
    UnivarPWLE<T> var3( var1 );  
    var3 -= var2;
    return var3;
  }
  UnivarPWLE<T> var3( -var2 );
    var3 += var1;
    return var3;
}


template <typename T>
inline
UnivarPWLE<T> operator-
( UnivarPWLE<T> const& var1, UnivarPWLE<T> && var2 )
{
  if( var1.first.empty() && var2.first.empty() ) 
    return( var1._cst - var2._cst );
  if( (var2.first.empty()) ){
    UnivarPWLE<T> var3( std::move(-var2) );
    var3 += var1;  
    return var3;
  }
  UnivarPWLE<T> var3( var1 );
  var3 -= var2;    
  return var3;
}


template <typename T>
inline
UnivarPWLE<T> operator-
( UnivarPWLE<T> && var1, UnivarPWLE<T> const& var2 )
{
  if( var1.first.empty() && var2.first.empty() ) 
    return( var1._cst - var2._cst );
  if( (var1.first.empty()) ){
    UnivarPWLE<T> var3( std::move(var1) );
    var3 -= var2;    
    return var3;
  }
  UnivarPWLE<T> var3( -var2 );
  var3 += var1;
  return var3;
}


template <typename T>
inline
UnivarPWLE<T> operator-
( UnivarPWLE<T> && var1, UnivarPWLE<T> && var2 )
{
  if( var1.first.empty() && var2.first.empty() ) 
    return( var1._cst - var2._cst );
  if( (var1.first.empty()) ){
    UnivarPWLE<T> var3( std::move(var1) );
    var3 -= var2;    
    return var3;
  }
  UnivarPWLE<T> var3( std::move(-var2) );
  var3 += var1;
  return var3;
}


template <typename T>
inline
UnivarPWLE<T> operator-
( UnivarPWLE<T> && var1, double const& cst2 )
{
  if( (var1.first.empty()) ) 
    return( var1._cst - cst2 );
  UnivarPWLE<T> var3( std::move(var1) );
  var3 += -cst2;
  return var3;
}


template <typename T>
inline
UnivarPWLE<T> operator-
( double const& cst1, UnivarPWLE<T> && var2 )
{
  if( (var2.first.empty()) ) 
    return( cst1 - var2._cst );
  UnivarPWLE<T> var3( std::move(-var2) );
  var3 += cst1;
  return var3;
}



template <typename T>
inline
UnivarPWLE<T> operator-
( UnivarPWLE<T> const& var1, double const& cst2 )
{
  if( (var1.first.empty()) ) 
    return( var1._cst - cst2 );
  UnivarPWLE<T> var3( var1 );
  var3 += -cst2;
  return var3;
}

template <typename T>
inline
UnivarPWLE<T> operator-
( double const& cst1, UnivarPWLE<T> const& var2 )
{
  if( (var2.first.empty()) ) 
    return( cst1 - var2._cst );
  UnivarPWLE<T> var3( -var2 );
  var3 += cst1;
  return var3;
}

template <typename T>
inline
UnivarPWLE<T>& UnivarPWLE<T>::operator*=
( double const& cst )
{
  if( cst == 0. ){
    *this = 0.;
    return *this;
  }
  if( cst == 1. ){
    return *this;
  }
  if( first.empty() ){
    _cst *= cst;
    return *this;
  }

  for(auto & yBkp:second){
    yBkp *= cst;
  }

  if(second.size()<=2){
    if(cst < 0.){
      _isUnder = !(_isUnder);
    }
    return *this;
  }


  if(_lbnd.second) _lbnd.first *= cst;
  if(_ubnd.second) _ubnd.first *= cst;
  if(cst < 0.){
    _isUnder = !(_isUnder);
    std::swap(_lbnd,_ubnd);
  }

  if(_yFwdDiff.second){ _yFwdDiff.second = false; _yFwdDiff.first.resize(0);
    // for(auto & yBkp:_yFwdDiff.first){
    //     yBkp *= cst;
    // }
  }

  return *this;
}

template <typename T>
inline
UnivarPWLE<T>& UnivarPWLE<T>::operator*=
( UnivarPWLE<T> const& var )
{
  throw typename UnivarPWLE<T>::Exceptions( UnivarPWLE<T>::Exceptions::UNDEF );
  return *this;
}

template <typename T>
inline
UnivarPWLE<T> operator*
( UnivarPWLE<T> const& var1, UnivarPWLE<T> const& var2 )
{
  throw typename UnivarPWLE<T>::Exceptions( UnivarPWLE<T>::Exceptions::UNDEF );

  return var1;
  //return intersect(var3,bnd1*bnd2);
}

template <typename T>
inline
UnivarPWLE<T> operator*
( UnivarPWLE<T> const& var1, double const& cst2 )
{
  if( var1.first.empty() ) 
    return( var1._cst * cst2 );

  UnivarPWLE<T> var3( var1 );
  var3 *= cst2;
  return var3;
}

template <typename T>
inline
UnivarPWLE<T> operator*
( double const& cst1, UnivarPWLE<T> const& var2 )
{
  if( var2.first.empty() ) 
    return( cst1 * var2._cst );

  UnivarPWLE<T> var3( var2 );
  var3 *= cst1;
  return var3;
}


template <typename T>
inline
UnivarPWLE<T> operator*
( UnivarPWLE<T> && var1, double const& cst2 )
{
  if( var1.first.empty() ) 
    return( var1._cst * cst2 );
  UnivarPWLE<T> var3( std::move(var1) );
  var3 *= cst2;
  return var3;
}


template <typename T>
inline
UnivarPWLE<T> operator*
( double const& cst1, UnivarPWLE<T> && var2 )
{
  if( var2.first.empty() ) 
    return( cst1 * var2._cst );
  UnivarPWLE<T> var3( std::move(var2) );
  var3 *= cst1;
  return var3;
}

template <typename T>
inline
UnivarPWLE<T>& UnivarPWLE<T>::operator/=
( double const& cst )
{
  if( cst == 0. )
    throw typename UnivarPWLE<T>::Exceptions( UnivarPWLE<T>::Exceptions::DIV );
  *this *= 1/cst; // inv( cst );
  return *this;
}

template <typename T>
inline
UnivarPWLE<T>& UnivarPWLE<T>::operator/=
( UnivarPWLE<T> const& var )
{
  throw typename UnivarPWLE<T>::Exceptions( UnivarPWLE<T>::Exceptions::UNDEF );
  //*this *= inv( var );
  return *this;
}

template <typename T>
inline
UnivarPWLE<T> operator/
( UnivarPWLE<T> const& var1, UnivarPWLE<T> const& var2 )
{
//  UnivarPWLE<T> var3( var1 );
//  var3 /= var2;
  throw typename UnivarPWLE<T>::Exceptions( UnivarPWLE<T>::Exceptions::UNDEF );
  return var1;
}

template <typename T>
inline
UnivarPWLE<T> operator/
( UnivarPWLE<T> const& var1, double const& cst2 )
{
  UnivarPWLE<T> var3( var1 );
  var3 /= cst2;
  return var3;
}

template <typename T>
inline
UnivarPWLE<T> operator/
( double const& cst1, UnivarPWLE<T> const& var2 )
{
  throw typename UnivarPWLE<T>::Exceptions( UnivarPWLE<T>::Exceptions::UNDEF );
  // UnivarPWLE<T> var3( inv( var2 ) );
  // var3 *= cst1;
  return var2;
}

template <typename T>
inline
UnivarPWLE<T> relu
( UnivarPWLE<T> const& var )
{
#ifdef MC__UPWLE_DEBUG_TRACE   
  std::cout << "in PWLE relu &" << std::endl;
#endif   
  if( var.first.empty() )
    return std::max(var._cst,T(0.));

  if( var.second.size() == 1 ){
    UnivarPWLE<T> var2( var );
    var2.second[0] = std::max(var2.second[0],T(0.));
#ifdef MC__UPWLE_DEBUG_TRACE  
    //std::cout << "out PWLE relu" << std::endl;
#endif
    return var2;
  } 

#ifdef MC__UPWLE_DEBUG 
  assert(var.second.size() == var.first.size());
#endif

  if(var.second.size() == 2){

    T lb = std::min(var.second[0],var.second[1]);
    if(lb >= -T(MC__UPWLE_COMPUTATION_TOL)){
      UnivarPWLE<T> var2( var );
#ifdef MC__UPWLE_DEBUG_TRACE  
      std::cout << "out PWLE relu" << std::endl;
      assert(var2.first.back() - var2.first[0] > -1e-13);
#endif
      return var2;
    }
         
    T ub = std::max(var.second[0],var.second[1]);
    if(ub <= T(MC__UPWLE_COMPUTATION_TOL)){
      //const T lb = var.first[0] + var.first[1];
      //const T ub = var.first[0] + var.first.back();      
      //UnivarPWLE<T> var2(lb,ub,T(0.),var._isUnder);  
      UnivarPWLE<T> var2(var.first[0],var.first[1],T(0.),var._isUnder);    
      //UnivarPWLE<T> var2( var.first,var._isUnder );
#ifdef MC__UPWLE_DEBUG_TRACE  
      std::cout << "out PWLE relu" << std::endl;
      assert(var2.first.back() - var2.first[0] > -1e-13);
#endif
      return var2;
    }

  
    //std::cout << "in PWLE relu" << std::endl;
    UnivarPWLE<T> var2(T(0.));
    var2.first.resize(4);
    var2.second.resize(4);
    var2._isUnder = var._isUnder;    
    var2._lbnd = std::make_pair(T(0.),true);
    var2._ubnd = std::make_pair(ub,true);    

    //std::cout << "Test Relu" << std::endl;
    var2.first[0] = 0.5*(var.first[0] + var.first[1]);
    var2.first[1] = var.first[0] - var2.first[0];
    var2.first[3] = var.first[1] - var2.first[0];

    var2.second[0] = 0.5*(var.second[0] + var.second[1]);
    if(var.second[1] > var.second[0]){
      var2.second[3] = var.second[1] - var2.second[0];   
      var2.second[1] = - var2.second[0];
    }
    else{
      var2.second[3] = - var2.second[0];   
      var2.second[1] = var.second[0] - var2.second[0];
    }
    var2.second[2] = - var2.second[0];
    

    var2.first[2]  = var2.first[3] - var.second[1] * ((var.first[1] - var.first[0])/ (var.second[1] - var.second[0]));

#ifdef MC__UPWLE_DEBUG    
    assert(var2.second[0]< 1e7);
    assert(var2.second[1]< 1e7);
    assert(var2.second[2]< 1e7);
    assert(var2.second[3]< 1e7);
#endif              
    var2._xFwdDiff = std::make_pair(std::vector<T>(0),false);
    var2._yFwdDiff = std::make_pair(std::vector<T>(0),false);   
#ifdef MC__UPWLE_DEBUG_TRACE  
    std::cout << "out PWLE relu" << std::endl;
    assert(var2.first.back() - var2.first[1] > -1e-13);
#endif
    return var2; 
  
  }
  else{
    if(var.get_ub() <= T(0.)){ // T(MC__UPWLE_COMPUTATION_TOL)
      const T lb = var.first[0] + var.first[1];
      const T ub = var.first[0] + var.first.back();      
      UnivarPWLE<T> var2(lb,ub,T(0.),var._isUnder);
#ifdef MC__UPWLE_DEBUG_TRACE  
      assert(ub - lb > -1e-13);
      //std::cout << "out PWLE relu" << std::endl;
#endif
      return var2;
    }
    if(var._lbnd.second){   
    // here we do not always check whether lb > 0 and that will be processed in _relu() anyway
      if(var._lbnd.first >= T(0.)){
        UnivarPWLE<T> var2( var );
#ifdef MC__UPWLE_DEBUG_TRACE  
        std::cout << "out PWLE relu" << std::endl;
        assert(var2.first.back() - var2.first[1] > -1e-13);
#endif
        return var2;
      }
    }    
  }


  UnivarPWLE<T> var2( var );
  var2._relu();    
  //std::cout << var2 << std::endl;
  //  var2._bnd.second = false; // processed in relu
  //  var2._bnd.first.first.second = false; 
  //  var2._bnd.first.second.second = false;    
  //var2._xFwdDiff;
  //var2._yFwdDiff;
#ifdef MC__UPWLE_DEBUG_TRACE  
  std::cout << "out PWLE relu" << std::endl;
  assert(var2.first.back() - var2.first[1] > -1e-13);
#endif
  return var2;
}

template <typename T>
inline
UnivarPWLE<T> relu
( UnivarPWLE<T> && var )
{
#ifdef MC__UPWLE_DEBUG_TRACE  
  std::cout << "in PWLE relu &&" << std::endl;
  //std::cout << var << std::endl;
#endif
  if( var.first.empty() )
    return std::max(var._cst,T(0.));


  if( var.second.size() == 1 ){
    //std::cout << "in PWLE relu" << std::endl;
    UnivarPWLE<T> var2( std::move(var) );
    var2.second[0] = std::max(var2.second[0],T(0.));
    // var2._bnd.second = true;
    // var2._bnd.first.first = std::make_pair(var2.second[0],true);
    // var2._bnd.first.second = std::make_pair(var2.second[0],true);
    //std::cout << "out PWLE relu" << std::endl;
    return var2;
  } 
#ifdef MC__UPWLE_DEBUG   
  assert(var.second.size() == var.first.size());
#endif
  if(var.second.size() == 2){
    //std::cout << "in PWLE relu" << std::endl;
    UnivarPWLE<T> var2( std::move(var) );
    T lb = std::min(var2.second[0],var2.second[1]);
    if(lb >= -T(MC__UPWLE_COMPUTATION_TOL)){
#ifdef MC__UPWLE_DEBUG_TRACE        
      std::cout << "out PWLE relu 2-X" << std::endl;
      assert(var2.first.back() - var2.first[0] > -1e-13);
#endif
      return var2;
    }
         
    T ub = std::max(var2.second[0],var2.second[1]);
    if(ub <= T(MC__UPWLE_COMPUTATION_TOL)){
      var2.second[0] = T(0.);
      var2.second.resize(1);
      // var2.first[1] += var2.first[0];
      // var2.first[0] += var2.first.back(); 
      //var2.first.resize(2);          
      //var2.set_zero();
#ifdef MC__UPWLE_DEBUG_TRACE  
      std::cout << "out PWLE relu 2-0" << std::endl;
#endif
      return var2;
    }

    //std::cout << "in PWLE relu ???" << std::endl;
    var2._lbnd = std::make_pair(T(0.),true);
    var2._ubnd = std::make_pair(ub,true);    

    //std::cout << "Test Relu" << std::endl;

    //std::cout << var2.first[0] << std::endl;
    //std::cout << var2.first[1] << std::endl;
    //std::cout << var2.second[0] << std::endl;
    //std::cout << var2.second[1] << std::endl; 
    //std::cout << std::endl; 
    var2.first.resize(4);
    var2.second.resize(4);

    var2.first[2] = 0.5*(var2.first[0] + var2.first[1]);
    var2.first[3] = var2.first[1] - var2.first[2];
    var2.first[1] = var2.first[0] - var2.first[2];
    

    var2.first[0] = var2.first[3] - var2.second[1] * ((var2.first[3] - var2.first[1])/ (var2.second[1] - var2.second[0]));
#ifdef MC__UPWLE_DEBUG_TRACE  
    std::cout << "first before swap" << std::endl;
    for(unsigned int jj = 0; jj < var2.first.size(); jj++ ){
      std::cout << var2.first[jj] << "    ";
    }
    std::cout << std::endl;
#endif
    std::swap(var2.first[0],var2.first[2]);
#ifdef MC__UPWLE_DEBUG_TRACE      
    std::cout << "first after swap" << std::endl;
    for(unsigned int jj =1; jj < var2.first.size(); jj++ ){
      std::cout << var2.first[0] + var2.first[jj] << "    ";
    }
    std::cout << std::endl;
#endif
    var2.second[2] = 0.5*(var2.second[0] + var2.second[1]);
    if(var2.second[1] > var2.second[0]){
      var2.second[3] = var2.second[1] - var2.second[2];   
      var2.second[1] = - var2.second[2];
    }
    else{
      var2.second[3] = - var2.second[2];   
      var2.second[1] = var2.second[0] - var2.second[2];
    } 
    var2.second[0] = - var2.second[2];
    std::swap(var2.second[0],var2.second[2]);
#ifdef MC__UPWLE_DEBUG_TRACE  
    // std::cout << var2.first[0] << std::endl;
    // std::cout << var2.first[1] << std::endl;
    // std::cout << var2.first[2] << std::endl;
    // std::cout << var2.first[3] << std::endl;    
    // std::cout << std::endl; 
    // std::cout << var2.second[0] << std::endl;
    // std::cout << var2.second[1] << std::endl;
    // std::cout << var2.second[2] << std::endl;
    // std::cout << var2.second[3] << std::endl;    
    // std::cout << std::endl; 
#endif
#ifdef MC__UPWLE_DEBUG
    assert(var2.second[0]< 1e7);
    assert(var2.second[1]< 1e7);
    assert(var2.second[2]< 1e7);
    assert(var2.second[3]< 1e7);
#endif

    var2._xFwdDiff = std::make_pair(std::vector<T>(0),false);
    var2._yFwdDiff = std::make_pair(std::vector<T>(0),false);   
#ifdef MC__UPWLE_DEBUG_TRACE  
    std::cout << "out PWLE relu 2 - T" << std::endl;
    assert(var2.first.back() - var2.first[1] > -1e-13);
#endif
    return var2; 
  
  }
  else{
    if(var.get_ub() <= T(0.)){ // T(MC__UPWLE_COMPUTATION_TOL)
//      const T lb = var.first[0] + var.first[1];
//      const T ub = var.first[0] + var.first.back();      
//      UnivarPWLE<T> var2(lb,ub,T(0.),var._isUnder);
      //UnivarPWLE<T> var2( std::move(var.first),var._isUnder );
      var.second[0] = T(0.);
      var.second.resize(1);
      var.first[1] += var.first[0];
      var.first[0] += var.first.back(); 
      std::swap(var.first[0],var.first[1]);
      var.first.resize(2);
      var._xFwdDiff.first.resize(0);
      var._yFwdDiff.first.resize(0);        
//      var.set_zero();
#ifdef MC__UPWLE_DEBUG_TRACE  
      std::cout << "out PWLE relu" << std::endl;
      assert(var.first.back() - var.first[0] > -1e-13);
#endif
      return var;
    }
    if(var._lbnd.second){   
    // here we do not always check whether lb > 0 and that will be processed in _relu() anyway
      if(var._lbnd.first >= T(0.)){
#ifdef MC__UPWLE_DEBUG_TRACE  
        std::cout << "out PWLE relu" << std::endl;
        assert(var.first.back() - var.first[1] > -1e-13);
#endif
        return var;
      }
    }    
  }


  //_relu( var2);
  //std::cout << "Test Relu" << std::endl;
  //std::cout << "in PWLE relu" << std::endl;
  var._relu();        
  //std::cout << var << std::endl;  
  //  var._bnd.second = false;  // processed in relu
  //  var._bnd.first.first.second = false; 
  //  var._bnd.first.second.second = false;    
#ifdef MC__UPWLE_DEBUG_TRACE  
  //std::cout << "out PWLE relu" << std::endl;
#endif
  return var;
}



//! @brief C++ class for Univariate Piecewise Linear Over- and Under-estimators
////////////////////////////////////////////////////////////////////////
//! mc::UnivarPWL is a C++ class for propagation of pairs of univariate piecewise linear 
//! estimators (UPWLE) through (univariate) factorable functions. The template  
//! parameter corresponds to the type used to propagate the breakpoints of UPWLEs.
////////////////////////////////////////////////////////////////////////
template <typename T> 
class UnivarPWL 
{
 public:
  //template <typename U> friend class ASModel;
  //template <typename U> friend class ASVar;

  template <typename U> friend std::ostream& operator<<
    ( std::ostream &, UnivarPWL<U> const& );

  template <typename U> friend UnivarPWL<U> operator+
    ( UnivarPWL<U> const& );
  template <typename U> friend UnivarPWL<U> operator+
    ( UnivarPWL<U> && );

  template <typename U> friend UnivarPWL<U> operator+
    ( UnivarPWL<U> const&, UnivarPWL<U> const& );
  template <typename U> friend UnivarPWL<U> operator+
    ( UnivarPWL<U> const&, UnivarPWL<U> && );
  template <typename U> friend UnivarPWL<U> operator+
    ( UnivarPWL<U> &&, UnivarPWL<U> const& );    
  template <typename U> friend UnivarPWL<U> operator+
    ( UnivarPWL<U> &&, UnivarPWL<U> && );    

  template <typename U> friend UnivarPWL<U> operator+
    ( double const&, UnivarPWL<U> const& );
  template <typename U> friend UnivarPWL<U> operator+
    ( double const&, UnivarPWL<U> && );

  template <typename U> friend UnivarPWL<U> operator+
    ( UnivarPWL<U> const&, double const& );
  template <typename U> friend UnivarPWL<U> operator+
    ( UnivarPWL<U> &&, double const& );


  template <typename U> friend UnivarPWL<U> operator-
    ( UnivarPWL<U> const& );
  template <typename U> friend UnivarPWL<U> operator-
    ( UnivarPWL<U> && );

  template <typename U> friend UnivarPWL<U> operator-
    ( UnivarPWL<U> const&, UnivarPWL<U> const& );
  template <typename U> friend UnivarPWL<U> operator-
    ( UnivarPWL<U> const&, UnivarPWL<U> && );
  template <typename U> friend UnivarPWL<U> operator-
    ( UnivarPWL<U> &&, UnivarPWL<U> const& );
  template <typename U> friend UnivarPWL<U> operator-
    ( UnivarPWL<U> &&, UnivarPWL<U> && );

  template <typename U> friend UnivarPWL<U> operator-
    ( double const&, UnivarPWL<U> const& );
  template <typename U> friend UnivarPWL<U> operator-
    ( double const&, UnivarPWL<U> && );


  template <typename U> friend UnivarPWL<U> operator-
    ( UnivarPWL<U> const&, double const& );
  template <typename U> friend UnivarPWL<U> operator-
    ( UnivarPWL<U> &&, double const& );



  template <typename U> friend UnivarPWL<U> operator*
    ( UnivarPWL<U> const&, UnivarPWL<U> const& );
  // template <typename U> friend UnivarPWL<U> operator*
  //   ( UnivarPWL<U> const&, UnivarPWL<U> && );
  // template <typename U> friend UnivarPWL<U> operator*
  //   ( UnivarPWL<U> &&, UnivarPWL<U> const& );
  // template <typename U> friend UnivarPWL<U> operator*
  //   ( UnivarPWL<U> &&, UnivarPWL<U> && );


  template <typename U> friend UnivarPWL<U> operator*
    ( double const&, UnivarPWL<U> const& );
  template <typename U> friend UnivarPWL<U> operator*
    ( double const&, UnivarPWL<U> && );


  template <typename U> friend UnivarPWL<U> operator*
    ( UnivarPWL<U> const&, double const& );
  template <typename U> friend UnivarPWL<U> operator*
    ( UnivarPWL<U> &&, double const& );


  template <typename U> friend UnivarPWL<U> operator/
    ( UnivarPWL<U> const&, UnivarPWL<U> const& );
  // template <typename U> friend UnivarPWL<U> operator/
  //   ( UnivarPWL<U> const&, UnivarPWL<U> && );
  // template <typename U> friend UnivarPWL<U> operator/
  //   ( UnivarPWL<U> &&, UnivarPWL<U> const& );
  // template <typename U> friend UnivarPWL<U> operator/
  //   ( UnivarPWL<U> &&, UnivarPWL<U> && );


  template <typename U> friend UnivarPWL<U> operator/
    ( double const&, UnivarPWL<U> const& );
  template <typename U> friend UnivarPWL<U> operator/
    ( double const&, UnivarPWL<U> && );

  template <typename U> friend UnivarPWL<U> operator/
    ( UnivarPWL<U> const&, double const& );
  template <typename U> friend UnivarPWL<U> operator/
    ( UnivarPWL<U> &&, double const& );

  // template <typename U> friend UnivarPWL<U> dotprod
  //   ( std::vector<UnivarPWL<U>> const&, double const& );

  // template <typename U> friend UnivarPWL<U> max
  //   ( UnivarPWL<U> const&, UnivarPWL<U> const& );
  // template <typename U> friend UnivarPWL<U> max
  //   ( UnivarPWL<U> const&, double const& );
  // template <typename U> friend UnivarPWL<U> max
  //   ( UnivarPWL<U> &&, double const& );
  // template <typename U> friend UnivarPWL<U> min
  //   ( UnivarPWL<U> const&, UnivarPWL<U> const& );
  // template <typename U> friend UnivarPWL<U> min
  //   ( UnivarPWL<U> const&, double const& );
  // template <typename U> friend UnivarPWL<U> min
  //   ( UnivarPWL<U> &&, double const& );

  // template <typename U> friend UnivarPWL<U> inv
  //   ( UnivarPWL<U> const& );
  // template <typename U> friend UnivarPWL<U> inv
  //   ( UnivarPWL<U>&& );
  // template <typename U> friend UnivarPWL<U> sqr
  //   ( UnivarPWL<U> const& );
  // template <typename U> friend UnivarPWL<U> sqr
  //   ( UnivarPWL<U> && );
  // template <typename U> friend UnivarPWL<U> sqrt
  //   ( UnivarPWL<U> const& );
  // template <typename U> friend UnivarPWL<U> sqrt
  //   ( UnivarPWL<U> && );
  // template <typename U> friend UnivarPWL<U> fabs
  //   ( UnivarPWL<U> const& );
  // template <typename U> friend UnivarPWL<U> fabs
  //   ( UnivarPWL<U> && );  
  template <typename U> friend UnivarPWL<U> relu
    ( UnivarPWL<U> const& );
  template <typename U> friend UnivarPWL<U> relu
    ( UnivarPWL<U> && );    
  // template <typename U> friend UnivarPWL<U> exp
  //   ( UnivarPWL<U> const& );
  // template <typename U> friend UnivarPWL<U> exp
  //   ( UnivarPWL<U> && );
  // template <typename U> friend UnivarPWL<U> log
  //   ( UnivarPWL<U> const& );
  // template <typename U> friend UnivarPWL<U> log
  //   ( UnivarPWL<U> && );
  // template <typename U> friend UnivarPWL<U> xlog
  //   ( UnivarPWL<U> const& );
  // template <typename U> friend UnivarPWL<U> xlog
  //   ( UnivarPWL<U> && );
  // template <typename U> friend UnivarPWL<U> sin
  //   ( UnivarPWL<U> const& );
  // template <typename U> friend UnivarPWL<U> sin
  //   ( UnivarPWL<U> && );
  // template <typename U> friend UnivarPWL<U> cos
  //   ( UnivarPWL<U> const& );
  // template <typename U> friend UnivarPWL<U> cos
  //   ( UnivarPWL<U> && );
  // template <typename U> friend UnivarPWL<U> tan
  //   ( UnivarPWL<U> const& );
  // template <typename U> friend UnivarPWL<U> tan
  //   ( UnivarPWL<U> && );
  // template <typename U> friend UnivarPWL<U> tanh
  //   ( UnivarPWL<U> const& );
  // template <typename U> friend UnivarPWL<U> tanh
  //   ( UnivarPWL<U> && );
  // template <typename U> friend UnivarPWL<U> pow
  //   ( UnivarPWL<U> const&, int const& n );
  // template <typename U> friend UnivarPWL<U> pow
  //   ( UnivarPWL<U> &&, int const& n );
  // template <typename U> friend UnivarPWL<U> pow
  //   ( UnivarPWL<U> const&, double const&  );
  // template <typename U> friend UnivarPWL<U> pow
  //   ( UnivarPWL<U> &&, double const&  );
  // template <typename U> friend UnivarPWL<U> intersect
  //   ( UnivarPWL<U> const& , U);
  // template <typename U> friend UnivarPWL<U> intersect
  //   ( UnivarPWL<U> && , U);  
  // template <typename U> friend UnivarPWL<U> affine_transform
  //   ( std::vector<UnivarPWL<U> > const& , std::vector<double> const&, const double);    
  // template <typename U> friend UnivarPWL<U> cheb
  //   ( UnivarPWL<U> const&, unsigned int const& n );
  // template <typename U> friend UnivarPWL<U> cheb
  //   ( UnivarPWL<U> &&, unsigned int const& n );

 public:
  //! @brief the upper limit of the number of breakpoints;
  static unsigned int nbpsMax;
  //! @brief underestimator 
  UnivarPWLE<double> undEst;
  //! @brief overestimator
  UnivarPWLE<double> oveEst;


 private:
  //! @brief Variable bound
  mutable std::pair<std::pair<double,double>,bool> _bnd;
  bool _isEmpty;

 public:

  // UnivarPWL
  // ()
  // : _isEmpty(true)
  // {

  // }

  UnivarPWL
  ()
  : undEst(0),oveEst(0), _bnd(std::make_pair(0.,0.),false),_isEmpty(true)
  {}

  //! @brief Constructor of UnivarPWL for a variable x_i on [a,b]
  UnivarPWL
  (T const& bnd)
  : undEst(Op<T>::l(bnd),Op<T>::u(bnd),true),oveEst(Op<T>::l(bnd),Op<T>::u(bnd),false),
  _bnd(std::make_pair(Op<T>::l(bnd),Op<T>::u(bnd)),true),_isEmpty(false)
{  
#ifdef MC__UPWL_DEBUG_TRACE  
    std::cout << "bnd: " << bnd << std::endl;
#endif
}

  //! @brief Scaling Constructor of UnivarPWL for a variable x_i on [a,b]
  UnivarPWL
  (T const& bnd, double mtpr)
  : undEst(Op<T>::l(bnd),Op<T>::u(bnd),true),oveEst(Op<T>::l(bnd),Op<T>::u(bnd),false),
  _bnd(std::make_pair(Op<T>::l(bnd)*mtpr,Op<T>::u(bnd)*mtpr),true),_isEmpty(false)
  {
    undEst *= mtpr;
    oveEst *= mtpr;
    if(mtpr < 0.){
      std::swap(undEst,oveEst);
      std::swap(_bnd.first.first,_bnd.first.second);
    }
#ifdef MC__UPWL_DEBUG_TRACE  
    std::cout << "bnd: " << bnd << std::endl;
    std::cout << "mtpr: " << mtpr << std::endl;
#endif
  }






  //! @brief Constructor of UnivarPWL for a constant c
  UnivarPWL
  ( const double & cst)
  : undEst(cst),oveEst(cst),_bnd(std::make_pair(cst,cst),true),_isEmpty(false)
  {  }


  //! @brief Copy constructor of UnivarPWL for a var (lvalue)
  UnivarPWL
  ( UnivarPWL<T> const& var )
  : undEst(0),oveEst(0),_bnd(std::make_pair(0.,0.),true),_isEmpty(false)
  {
#ifdef MC__UPWL_DEBUG_TRACE
    std::cerr << "-- UnivarPWL( UnivarPWL<T> const& )\n";
#endif
#ifdef TEST_MOVE
    std::cout << "Copy Constructor" << std::endl;
#endif
    if( this == &var ) return;
    undEst = var.undEst;
    oveEst = var.oveEst;
    _bnd = var._bnd;
    _isEmpty = var._isEmpty;
  }


  //! @brief Copy constructor of UnivarPWL for a var (rvalue)
  UnivarPWL
  ( UnivarPWL<T> && var )
  : undEst(0),oveEst(0),_bnd(std::make_pair(0.,0.),false),_isEmpty(false)  
  {
#ifdef MC__UPWL_DEBUG_TRACE
    std::cerr << "-- UnivarPWL( UnivarPWL<T> && var )\n";
#endif
#ifdef TEST_MOVE
    std::cout << "Move Constructor" << std::endl;
#endif
    if( this == &var ) return;
    undEst = std::move(var.undEst);
    oveEst = std::move(var.oveEst);
    _bnd = std::move(var._bnd);     
    _isEmpty = var._isEmpty;
  }


  //! @brief Copy constructor of UnivarPWL for a var (lvalue)
  UnivarPWL
  ( UnivarPWLE<double> const& est1, UnivarPWLE<double> const& est2)
  : undEst(est1),oveEst(est2),_bnd(std::make_pair(0.,0.),false),_isEmpty(false)
  {
    // if(est1._bnd.second && est1._bnd.first.first.second){
    //   if(est2._bnd.second && est2._bnd.first.second.second)
    //   _bnd = std::make_pair(std::make_pair(est1._bnd.first.first.first,est2._bnd.first.second.first),true);
    // }
  }

  //! @brief Copy constructor (with scaling) of UnivarPWL for a var (lvalue)
  UnivarPWL
  ( UnivarPWL<T> const& var, const double mtpr ) // multiplier
  : undEst(var.undEst,mtpr),oveEst(var.oveEst,mtpr),_isEmpty(false)
  {
#ifdef MC__UPWL_DEBUG_TRACE
    std::cerr << "-- UnivarPWL( UnivarPWL<T> const&, const double )\n";
#endif
#ifdef TEST_MOVE
    std::cout << "Copy and Scale Constructor" << std::endl;
#endif
    if(var._bnd.second){
    _bnd = std::make_pair(std::make_pair(var._bnd.first.first * mtpr,var._bnd.first.second * mtpr),true);
    }  
    if (mtpr < 0.){
      std::swap(undEst,oveEst);
      std::swap(_bnd.first.first,_bnd.first.second);
    }
  }



#ifdef UnivarPWL_LIFITIME_DEBUG     
  //! @brief Destructor of UnivarPWL (for checking the lifetime of UnivarPWLs)
  ~UnivarPWL() 
  {
    std::cout<< "UnivarPWL delated" <<std::endl;
  }
#endif    

  //! @brief Exceptions of mc::UnivarPWL
  class Exceptions
  {
   public:
    //! @brief Enumeration type for UnivarPWL exception handling
    enum TYPE{
      INDEX=1,	   //!< breakpoint index out of range
      ESTMATCH,    //!< Addition is defined for two under- or over-estimators
      UNDEF,	     //!< Feature not yet implemented
      DIV          //!< Division by zero scalar      
    };
    //! @brief Constructor for error <a>ierr</a>
    Exceptions( TYPE ierr ) : _ierr( ierr ){}
    //! @brief Error flag
    int ierr(){ return _ierr; }
    //! @brief Error description
    std::string what(){
      switch( _ierr ){
      case INDEX:
        return "mc::UnivarPWL\t breakpoint index out of range";
      case ESTMATCH:
        return "mc::UnivarPWL\t Addition cannot be performed between an over- and an under-estimator";
      case UNDEF:
        return "mc::UnivarPWL\t Feature not yet implemented";
      case DIV:
        return "mc::UnivarPWL\t Divided by zero";        
      default:
        return "mc::UnivarPWL\t Internal error";
      }
    }
   private:
    TYPE _ierr;
  };


  UnivarPWL<T>& operator=
  ( T const& cst )
  {
    undEst = cst;
    oveEst = cst;
    _isEmpty = false;
    _bnd = std::make_pair(std::make_pair(cst,cst),true);
    return *this;
  }


  UnivarPWL<T>& operator=
  ( UnivarPWL<T> const& var )
  {
#ifdef MC__UPWL_DEBUG_TRACE
    std::cerr << "-- UnivarPWL<T>& operator= ( UnivarPWL<T> const& )\n";
#endif
#ifdef TEST_MOVE
    std::cout << "Copy Assign" << std::endl;
#endif
    if( this == &var )
      return *this;

    undEst = var.undEst;
    oveEst = var.oveEst;
    _bnd = var._bnd;
    _isEmpty = var._isEmpty;
    return *this;
  } 

  UnivarPWL<T>& operator=
  ( UnivarPWL<T> && var )
  {
#ifdef MC__UPWL_DEBUG_TRACE
    std::cerr << "-- UnivarPWL<T>& operator= ( UnivarPWL<T> && )\n";
#endif
#ifdef TEST_MOVE
    std::cout << "Move Assign" << std::endl;
#endif
    if( this == &var )
      return *this;

    undEst = std::move(var.undEst);
    oveEst = std::move(var.oveEst);
    _bnd = std::move(var._bnd);
    _isEmpty = var._isEmpty;
    return *this;
  }


          // ASVar<T>& set
          // ( ASModel<T>* const mod )
          // {
        
        
          //   return *this;
          // }
 
  UnivarPWL<T>& operator+=
    ( UnivarPWL<T> const& );
  UnivarPWL<T>& operator+=
    ( double const& );
  // UnivarPWL<T>& operator+=
  //   ( T const& );


  UnivarPWL<T>& operator-=
    ( UnivarPWL<T> const& );
  UnivarPWL<T>& operator-=
    ( double const& );
  // UnivarPWL<T>& operator-=
  //   ( T const& );
 
  UnivarPWL<T>& operator*=
    ( UnivarPWL<T> const& );
  UnivarPWL<T>& operator*=
    ( double const& );
  // UnivarPWL<T>& operator*=
  //   ( T const& );

  UnivarPWL<T>& operator/=
    ( UnivarPWL<T> const& );
  UnivarPWL<T>& operator/=
    ( double const& );


  bool empty() const {
    return _isEmpty;
  }

  void flag_nonEmpty() {
    _isEmpty = false;
    return ;
  }



  void B() const {
    _B();
  }

  std::pair<double,double> getBnd() const{
    _B();
    return _bnd.first;
  }  

  void debug_check_overNunder_flags()
  const
  {
    if(oveEst.get_flag()){
      std::cout << "ove flag is not correct " << oveEst.get_flag() << std::endl;
    }
    if(!undEst.get_flag()){
      std::cout << "und flag is not correct " << undEst.get_flag() << std::endl;
    }
  }


  void debug_check_overNunder_flags(std::vector<double> shadow_info)
  const
  {

    if(shadow_info[1] > 0 && oveEst.get_flag()){
      std::cout << "ove flag is not correct " << oveEst.get_flag() << std::endl;
    }
    if(shadow_info[0] > 0 && !undEst.get_flag()){
      std::cout << "und flag is not correct " << undEst.get_flag() << std::endl;
    }
  }


 private:

  //! @brief Return the bound of the UnivarPWL <a>var</a>
  void _B
  ()
  const
  {
    // undEst.B(true,false);
    // oveEst.B(false,true);
    // std::cout << "lb" << undEst._bnd.first.first.first << std::endl;
    // std::cout << "ub" << oveEst._bnd.first.second.first << std::endl;
//    _bnd = std::make_pair(std::make_pair(undEst._bnd.first.first.first,oveEst._bnd.first.second.first),true);
     _bnd = std::make_pair(std::make_pair(undEst.get_lb(),oveEst.get_ub()),true);
  }

};


template <typename T> inline
unsigned int UnivarPWL<T>::nbpsMax;


template <typename T>
inline
UnivarPWL<T> operator+
( UnivarPWL<T> const& var )
{
  return var;
}

template <typename T>
inline
UnivarPWL<T> operator+
( UnivarPWL<T> && varIn )
{
  UnivarPWL<T> var( std::move(varIn) );  
  return var;
}


template <typename T>
inline
UnivarPWL<T>& UnivarPWL<T>::operator+=
( double const& cst )
{
  if( cst == 0. ){
    return *this;
  }
#ifdef MC__UPWL_DEBUG_TRACE  
  std::cout << "UPWL += cst" << std::endl;
#endif
  undEst += cst;
  oveEst += cst;  
  _bnd.second = false;
  
  // if( _bnd.second ){
  //   _bnd.first.first += cst;
  //   _bnd.first.second += cst;
  // }

  return *this;
}

// Reserved for adding constant intervals/other sets


template <typename T>
inline
UnivarPWL<T>& UnivarPWL<T>::operator+=
( UnivarPWL<T> const& var )
{
  undEst += var.undEst;
  oveEst += var.oveEst;  
  _bnd.second = false;

  return *this;
}        



template <typename T>
inline
UnivarPWL<T> operator+
( UnivarPWL<T> const& var1, UnivarPWL<T> const& var2 )
{
    UnivarPWL<T> var3( var1 );
    var3 += var2;
    return var3;
}


template <typename T>
inline
UnivarPWL<T> operator+
( UnivarPWL<T> const& var1, UnivarPWL<T> && var2 )
{

    UnivarPWL<T> var3( std::move(var2) );
    var3 += var1;
    return var3;
 
}


template <typename T>
inline
UnivarPWL<T> operator+
( UnivarPWL<T> && var1, UnivarPWL<T> const& var2 )
{
    UnivarPWL<T> var3( std::move(var1) );
    var3 += var2;
    return var3;
}


template <typename T>
inline
UnivarPWL<T> operator+
( UnivarPWL<T> && var1, UnivarPWL<T> && var2 )
{

    UnivarPWL<T> var3( std::move(var1) );
    var3 += var2;
    return var3;
}



template <typename T>
inline
UnivarPWL<T> operator+
( UnivarPWL<T> const& var1, double const& cst2 )
{
  UnivarPWL<T> var3( var1 );
  var3 += cst2;
  return var3;
}


template <typename T>
inline
UnivarPWL<T> operator+
( UnivarPWL<T> && var1, double const& cst2 )
{
  UnivarPWL<T> var3( std::move(var1) );
  var3 += cst2;
  return var3;
}


template <typename T>
inline
UnivarPWL<T> operator+
( double const& cst1, UnivarPWL<T> const& var2 )
{
  UnivarPWL<T> var3( var2 );
  var3 += cst1;
  return var3;
}

template <typename T>
inline
UnivarPWL<T> operator+
( double const& cst1, UnivarPWL<T> && var2 )
{
  UnivarPWL<T> var3( std::move(var2) );
  var3 += cst1;
  return var3;
}


template <typename T> 
inline
UnivarPWL<T> operator-
( UnivarPWL<T> const& var )
{
  UnivarPWL<T> var2(-var.oveEst,-var.undEst);
  var2._bnd.second = false;
  return var2;
}



template <typename T> 
inline
UnivarPWL<T> operator-
( UnivarPWL<T> && var )
{
  std::swap(var.oveEst,var.undEst);
  var.oveEst = -var.oveEst;
  var.undEst = -var.undEst;
  var._bnd.second = false;

  // if(var._bnd.second){
  //   var._bnd.first.first *= -1.;
  //   var._bnd.first.second *= -1.;    
  //   std::swap(var._bnd.first.first,var._bnd.first.second);
  // }
  return var;
}



template <typename T>
inline
UnivarPWL<T>& UnivarPWL<T>::operator-=
( double const& cst )
{
  *this += -cst;
  return *this;
}

template <typename T>
inline
UnivarPWL<T>& UnivarPWL<T>::operator-=
( UnivarPWL<T> const& var )
{
  *this += (-var);
  _bnd.second = false;   
  return *this;
}


template <typename T>
inline
UnivarPWL<T> operator-
( UnivarPWL<T> const& var1, UnivarPWL<T> const& var2 )
{
    UnivarPWL<T> var3( var1 );  
    var3 -= var2;
    return var3;
}


template <typename T>
inline
UnivarPWL<T> operator-
( UnivarPWL<T> const& var1, UnivarPWL<T> && var2 )
{
    UnivarPWL<T> var3( std::move(-var2) );
    var3 += var1;  
    return var3;
}


template <typename T>
inline
UnivarPWL<T> operator-
( UnivarPWL<T> && var1, UnivarPWL<T> const& var2 )
{
  UnivarPWL<T> var3( std::move(var1) );
  var3 -= var2;    
  return var3;
}


template <typename T>
inline
UnivarPWL<T> operator-
( UnivarPWL<T> && var1, UnivarPWL<T> && var2 )
{
  UnivarPWL<T> var3( std::move(var1) );
  var3 -= var2;    
  return var3;
}


template <typename T>
inline
UnivarPWL<T> operator-
( UnivarPWL<T> && var1, double const& cst2 )
{
  UnivarPWL<T> var3( std::move(var1) );
  var3 += -cst2;
  return var3;
}


template <typename T>
inline
UnivarPWL<T> operator-
( double const& cst1, UnivarPWL<T> && var2 )
{
  UnivarPWL<T> var3( std::move(-var2) );
  var3 += cst1;
  return var3;
}



template <typename T>
inline
UnivarPWL<T> operator-
( UnivarPWL<T> const& var1, double const& cst2 )
{
  UnivarPWL<T> var3( var1 );
  var3 += -cst2;
  return var3;
}

template <typename T>
inline
UnivarPWL<T> operator-
( double const& cst1, UnivarPWL<T> const& var2 )
{
  UnivarPWL<T> var3( -var2 );
  var3 += cst1;
  return var3;
}

template <typename T>
inline
UnivarPWL<T>& UnivarPWL<T>::operator*=
( double const& cst )
{
  if( cst == 0. ){
    UnivarPWL<T> _tmp(0.);
    *this = _tmp;
    return *this;
  }
  if( cst == 1. ){
    return *this;
  }

  undEst *= cst;
  oveEst *= cst;
  _bnd.second = false;  
  
  // if(_bnd.second){
  //   _bnd.first.first *= cst;
  //   _bnd.first.second *= cst;
  // }   
  if(cst < 0.){
    std::swap(undEst,oveEst);  
    // std::swap(_bnd.first.first,_bnd.first.second);  
  }


  return *this;
}

template <typename T>
inline
UnivarPWL<T>& UnivarPWL<T>::operator*=
( UnivarPWL<T> const& var )
{
  throw typename UnivarPWL<T>::Exceptions( UnivarPWL<T>::Exceptions::UNDEF );
  return *this;
}

template <typename T>
inline
UnivarPWL<T> operator*
( UnivarPWL<T> const& var1, UnivarPWL<T> const& var2 )
{
  throw typename UnivarPWL<T>::Exceptions( UnivarPWL<T>::Exceptions::UNDEF );

  return var1;
  //return intersect(var3,bnd1*bnd2);
}

template <typename T>
inline
UnivarPWL<T> operator*
( UnivarPWL<T> const& var1, double const& cst2 )
{
  UnivarPWL<T> var3( var1 );
  var3 *= cst2;
  return var3;
}

template <typename T>
inline
UnivarPWL<T> operator*
( double const& cst1, UnivarPWL<T> const& var2 )
{
  UnivarPWL<T> var3( var2 );
  var3 *= cst1;
  return var3;
}


template <typename T>
inline
UnivarPWL<T> operator*
( UnivarPWL<T> && var1, double const& cst2 )
{
  UnivarPWL<T> var3( std::move(var1) );
  var3 *= cst2;
  return var3;
}


template <typename T>
inline
UnivarPWL<T> operator*
( double const& cst1, UnivarPWL<T> && var2 )
{
  UnivarPWL<T> var3( std::move(var2) );
  var3 *= cst1;
  return var3;
}

template <typename T>
inline
UnivarPWL<T>& UnivarPWL<T>::operator/=
( double const& cst )
{
  if( cst == 0. )
    throw typename UnivarPWL<T>::Exceptions( UnivarPWL<T>::Exceptions::DIV );
  *this *= 1.0/cst; // inv( cst );
  return *this;
}

template <typename T>
inline
UnivarPWL<T>& UnivarPWL<T>::operator/=
( UnivarPWL<T> const& var )
{
  throw typename UnivarPWL<T>::Exceptions( UnivarPWL<T>::Exceptions::UNDEF );
  //*this *= inv( var );
  return *this;
}

template <typename T>
inline
UnivarPWL<T> operator/
( UnivarPWL<T> const& var1, UnivarPWL<T> const& var2 )
{
//  UnivarPWL<T> var3( var1 );
//  var3 /= var2;
  throw typename UnivarPWL<T>::Exceptions( UnivarPWL<T>::Exceptions::UNDEF );
  return var1;
}

template <typename T>
inline
UnivarPWL<T> operator/
( UnivarPWL<T> const& var1, double const& cst2 )
{
  UnivarPWL<T> var3( var1 );
  var3 /= cst2;
  return var3;
}

template <typename T>
inline
UnivarPWL<T> operator/
( double const& cst1, UnivarPWL<T> const& var2 )
{
  throw typename UnivarPWL<T>::Exceptions( UnivarPWL<T>::Exceptions::UNDEF );
  // UnivarPWL<T> var3( inv( var2 ) );
  // var3 *= cst1;
  return var2;
}

template <typename T>
inline
UnivarPWL<T> relu
( UnivarPWL<T> const& var )
{
  UnivarPWL<T> var2( var );

  var2.undEst = relu(var2.undEst);  
  var2.oveEst = relu(var2.oveEst); 
  var2._bnd.second = false;
  return var2;
}

template <typename T>
inline
UnivarPWL<T> relu
( UnivarPWL<T> && var )
{
  UnivarPWL<T> var2( std::move(var) );
  var2.undEst = relu(var2.undEst);  
  var2.oveEst = relu(var2.oveEst); 
  var2._bnd.second = false;
  return var;
}




}//namespace mc

// #include "mcfadbad.hpp"

// namespace fadbad
// {

// //! @brief Specialization of the structure fadbad::Op for use of the type mc::ISVar of MC++ as a template parameter of the classes fadbad::F, fadbad::B and fadbad::T of FADBAD++
// template< typename T > struct Op< mc::ISVar<T> >
// { 
//   typedef mc::ISVar<T> ISV;
//   typedef double Base;
//   static Base myInteger( const int i ) { return Base(i); }
//   static Base myZero() { return myInteger(0); }
//   static Base myOne() { return myInteger(1);}
//   static Base myTwo() { return myInteger(2); }
//   static double myPI() { return mc::PI; }
//   static ISV myPos( const ISV& x ) { return  x; }
//   static ISV myNeg( const ISV& x ) { return -x; }
//   template <typename U> static ISV& myCadd( ISV& x, const U& y ) { return x+=y; }
//   template <typename U> static ISV& myCsub( ISV& x, const U& y ) { return x-=y; }
//   template <typename U> static ISV& myCmul( ISV& x, const U& y ) { return x*=y; }
//   template <typename U> static ISV& myCdiv( ISV& x, const U& y ) { return x/=y; }
//   static ISV myInv( const ISV& x ) { return mc::inv( x ); }
//   static ISV mySqr( const ISV& x ) { return mc::sqr( x ); }
//   template <typename X, typename Y> static ISV myPow( const X& x, const Y& y ) { return mc::pow( x, y ); }
//   //static ISV myCheb( const ISV& x, const unsigned n ) { return mc::cheb( x, n ); }
//   static ISV mySqrt( const ISV& x ) { return mc::sqrt(x); }
//   static ISV myLog( const ISV& x ) { return mc::log( x ); }
//   static ISV myExp( const ISV& x ) { return mc::exp( x ); }
//   static ISV mySin( const ISV& x ) { return mc::sin( x ); }
//   static ISV myCos( const ISV& x ) { return mc::cos( x ); }
//   static ISV myTan( const ISV& x ) { throw typename mc::ISModel<T>::Exceptions( mc::ISModel<T>::Exceptions::UNDEF ); } //{ return mc::tan( x ); }
//   static ISV myAsin( const ISV& x ) { throw typename mc::ISModel<T>::Exceptions( mc::ISModel<T>::Exceptions::UNDEF ); }
//   static ISV myAcos( const ISV& x ) { throw typename mc::ISModel<T>::Exceptions( mc::ISModel<T>::Exceptions::UNDEF ); }
//   static ISV myAtan( const ISV& x ) { throw typename mc::ISModel<T>::Exceptions( mc::ISModel<T>::Exceptions::UNDEF ); }
//   static ISV mySinh( const ISV& x ) { throw typename mc::ISModel<T>::Exceptions( mc::ISModel<T>::Exceptions::UNDEF ); }
//   static ISV myCosh( const ISV& x ) { throw typename mc::ISModel<T>::Exceptions( mc::ISModel<T>::Exceptions::UNDEF ); }
//   static ISV myTanh( const ISV& x ) { throw typename mc::ISModel<T>::Exceptions( mc::ISModel<T>::Exceptions::UNDEF ); }
//   static bool myEq( const ISV& x, const ISV& y ) { return mc::Op<T>::eq(const_cast<ISV*>(&x)->bound(),const_cast<ISV*>(&y)->bound()); } 
//   static bool myNe( const ISV& x, const ISV& y ) { return mc::Op<T>::ne(const_cast<ISV*>(&x)->bound(),const_cast<ISV*>(&y)->bound()); }
//   static bool myLt( const ISV& x, const ISV& y ) { return mc::Op<T>::lt(const_cast<ISV*>(&x)->bound(),const_cast<ISV*>(&y)->bound()); }
//   static bool myLe( const ISV& x, const ISV& y ) { return mc::Op<T>::le(const_cast<ISV*>(&x)->bound(),const_cast<ISV*>(&y)->bound()); }
//   static bool myGt( const ISV& x, const ISV& y ) { return mc::Op<T>::gt(const_cast<ISV*>(&x)->bound(),const_cast<ISV*>(&y)->bound()); }
//   static bool myGe( const ISV& x, const ISV& y ) { return mc::Op<T>::ge(const_cast<ISV*>(&x)->bound(),const_cast<ISV*>(&y)->bound()); }
// };

// } // end namespace fadbad

// namespace mc
// {

// //! @brief C++ structure for specialization of the mc::Op templated structure for use of mc::ISVar in DAG evaluation and as template parameter in other MC++ types
// template< typename T > struct Op< mc::ISVar<T> >
// {
//   typedef mc::ISVar<T> ISV;
//   static ISV point( const double c ) { return ISV(c); }
//   static ISV zeroone() { return ISV( mc::Op<T>::zeroone() ); }
//   static void I(ISV& x, const ISV&y) { x = y; }
//   static double l(const ISV& x) { return mc::Op<T>::l(x.B()); }
//   static double u(const ISV& x) { return mc::Op<T>::u(x.B()); }
//   static double abs (const ISV& x) { return mc::Op<T>::abs(x.B());  }
//   static double mid (const ISV& x) { return mc::Op<T>::mid(x.B());  }
//   static double diam(const ISV& x) { return mc::Op<T>::diam(x.B()); }
//   static ISV inv (const ISV& x) { return mc::inv(x);  }
//   static ISV sqr (const ISV& x) { return mc::sqr(x);  }
//   static ISV sqrt(const ISV& x) { return mc::sqrt(x); }
//   static ISV exp (const ISV& x) { return mc::exp(x);  }
//   static ISV&& exp (ISV&& x) { return mc::exp(std::move(x));  }
//   static ISV log (const ISV& x) { return mc::log(x);  }
//   static ISV xlog(const ISV& x) { return mc::xlog(x); }
//   static ISV lmtd(const ISV& x, const ISV& y) { return (x-y)/(mc::log(x)-mc::log(y)); }
//   static ISV rlmtd(const ISV& x, const ISV& y) { return (mc::log(x)-mc::log(y))/(x-y); }
//   static ISV fabs(const ISV& x) { return mc::fabs(x); }
//   static ISV sin (const ISV& x) { return mc::sin(x);  }
//   static ISV cos (const ISV& x) { return mc::cos(x);  }
//   static ISV tan (const ISV& x) { throw typename ISModel<T>::Exceptions( ISModel<T>::Exceptions::UNDEF ); } //{ return mc::tan(x);  }
//   static ISV asin(const ISV& x) { throw typename ISModel<T>::Exceptions( ISModel<T>::Exceptions::UNDEF ); }
//   static ISV acos(const ISV& x) { throw typename ISModel<T>::Exceptions( ISModel<T>::Exceptions::UNDEF ); }
//   static ISV atan(const ISV& x) { throw typename ISModel<T>::Exceptions( ISModel<T>::Exceptions::UNDEF ); }
//   static ISV sinh(const ISV& x) { throw typename ISModel<T>::Exceptions( ISModel<T>::Exceptions::UNDEF ); }
//   static ISV cosh(const ISV& x) { throw typename ISModel<T>::Exceptions( ISModel<T>::Exceptions::UNDEF ); }
//   static ISV tanh(const ISV& x) { return mc::tanh(x); } // {ISV y(std::move(x)); return mc::tanh(std::move(y));} //
//   static ISV&& tanh(ISV&& x) { return mc::tanh(std::forward<ISV>(x)); }// { return mc::tanh(std::move(x)); }
//   static ISV erf (const ISV& x) { throw typename ISModel<T>::Exceptions( ISModel<T>::Exceptions::UNDEF ); }
//   static ISV erfc(const ISV& x) { throw typename ISModel<T>::Exceptions( ISModel<T>::Exceptions::UNDEF ); }
//   static ISV fstep(const ISV& x) { throw typename ISModel<T>::Exceptions( ISModel<T>::Exceptions::UNDEF ); }
//   static ISV bstep(const ISV& x) { throw typename ISModel<T>::Exceptions( ISModel<T>::Exceptions::UNDEF ); }
//   static ISV hull(const ISV& x, const ISV& y) { throw typename ISModel<T>::Exceptions( ISModel<T>::Exceptions::UNDEF ); } //{ return mc::hull(x,y); }
//   static ISV min (const ISV& x, const ISV& y) { return mc::min(x,y); }
//   static ISV max (const ISV& x, const ISV& y) { return mc::max(x,y); }
//   static ISV arh (const ISV& x, const double k) { return mc::exp(-k/x); }
//   template <typename X, typename Y> static ISV pow(const X& x, const Y& y) { return mc::pow(x,y); }
//   static ISV cheb(const ISV& x, const unsigned n) { return mc::cheb(x,n); }
//   static ISV prod (const unsigned n, const ISV* x) { return mc::prod(n,x); }
//   static ISV monom (const unsigned n, const ISV* x, const unsigned* k) { return mc::monom(n,x,k); }
//   static bool inter(ISV& xIy, const ISV& x, const ISV& y) { throw typename ISModel<T>::Exceptions( ISModel<T>::Exceptions::UNDEF ); } //{ return mc::inter(xIy,x,y); }
//   static bool eq(const ISV& x, const ISV& y) { return mc::Op<T>::eq(x.B(),y.B()); }
//   static bool ne(const ISV& x, const ISV& y) { return mc::Op<T>::ne(x.B(),y.B()); }
//   static bool lt(const ISV& x, const ISV& y) { return mc::Op<T>::lt(x.B(),y.B()); }
//   static bool le(const ISV& x, const ISV& y) { return mc::Op<T>::le(x.B(),y.B()); }
//   static bool gt(const ISV& x, const ISV& y) { return mc::Op<T>::gt(x.B(),y.B()); }
//   static bool ge(const ISV& x, const ISV& y) { return mc::Op<T>::ge(x.B(),y.B()); }
// };

// } // namespace mc

#endif
