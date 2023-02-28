// Copyright (C) Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
#define NOTTOTRACKSHADOW
#define MC__ISM_COMPUTATION_TOL  1e-13

/*!
\page page_ISM Interval Superposition Model Arithmetic for Factorable Functions
\author Yanlin Zha, Beno&icirc;t Chachuat

Interval superposition is a set-arithmetic for computing non-convex enclosures in the form of piecewise interval enclosures over a sub-paving of the domain. Formally, given an interval domain \f${\bf X} = X_1\times\cdots\times X_n = [\underline{x}_1,\overline{x}_1]\times\cdots\times[\underline{x}_n,\overline{x}_n]\f$ with \f$n\geq 1\f$, and equi-partitions of size \f$p\geq 1\f$ for each \f$X_i\f$ as
\f{align*}
X_{i,j} = [\underline{x}_i+(j-1)h_i,\overline{x}_i+jh_i], \quad \text{with}\ h_i=\frac{\overline{x}_i-\underline{x}_i}{p}, \quad i=1\ldots n, \quad j=1\ldots p,
\f}
an interval superposition model (ISM) of the real-valued function \f$f:\mathbb{R}^n\to\mathbb{R}\f$ on \f${\bf X}\f$ is any interval-valued function \f$\mathcal{I}[f]:{\bf X}\to\mathbb{IR}\f$ given by 
\f{align*}
\mathcal{I}[f](x) = \sum_{i=1}^{n} \sum_{j=1}^{q} A[f]_{i,j}\varphi_{i,j}(x_i), \quad \text{with}\ {\bf A}[f]\in\mathbb{IR}^{n\times p}, \quad \text{and}\ \varphi_{i,j}(x_i) = \left\{\begin{aligned} 1 & \ \text{if $x_i\in X_{i,j}$},\\ 0 & \ \text{otherwise},\end{aligned}\right.
\f}
such that 
\f{align*}
\forall {\bf x}\in {\bf X}, \quad f({\bf x})\in \mathcal{I}[f]({\bf x})\,.
\f}
The interval matrix \f${\bf A}[f]\f$ completely determines the enclosure \f$\mathcal{I}[f]\f$, and arithmetic rules have been developed that enable propagation of ISM through expression trees that comprise standard unary and binary operations, such as \f$\exp\f$, \f$\sin\f$, \f$+\f$, \f$\times\f$, \f$\ldots\f$ [Zha <I>et al.</I>, 2018]. As well as enabling tighter enclosures for factorable functions compared with polynomial models on wide parameter domains, ISM can be beneficial in set-inversion algorithms, e.g. with applications to guaranteed parameter estimation [Su <I>et al.</I>, 2018]. Nevertheless, ISM arithmetic comes with a larger computational burden than classical interval arithmetic. 

The classes mc::ISModel and mc::ISVar provide an implementation of ISM arithmetic based on the operator/function overloading mechanism of C++. This makes ISM both simple and intuitive to compute, similar to computing function values in real arithmetics or function bounds in interval arithmetic (see \ref page_INTERVAL). mc::ISModel and mc::ISVar are templated in the type used to propagate the interval coefficients. In addition, mc::ISVar can be used as the template parameter of other available types in MC++; as well as types in <A href="http://www.fadbad.com/fadbad.html">FADBAD++</A> for computing ISM of the partial derivatives or Taylor coefficients of a factorable function.


\section sec_ISM_use How do I compute an ISM of a factorable function?

Suppose we want to compute an ISM for the real-valued function \f$f(x,y)=x\exp(x+y^2)-y^2\f$ with \f$(x,y)\in [1,2]\times[0,1]\f$. For simplicity, bounds on the remainder terms are computed using the default interval type mc::Interval here:

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

#ifndef MC__ISMODEL_HPP
#define MC__ISMODEL_HPP

#include <iostream>
#include <iomanip> 
#include <vector>
#include <algorithm>
#include <cmath>
#include <limits>
#include <bitset>
#include <cassert>
#include "mcop.hpp"
#include "mcfunc.hpp"

//#define MC__ISMODEL_TRACE
#undef  MC__ISMODEL_DEBUG_PROD

namespace mc
{

template <typename T> class ISVar;

//! @brief C++ class for interval superposition models of factorable functions: environment
////////////////////////////////////////////////////////////////////////
//! mc::ISModel is a C++ class for definition of interval superposition
//! model (ISM) environment. ISM propagation of factorable functions is
//! implemented via the C++ class mc::ISVar. The template parameter
//! corresponds to the type used to propagate the interval coefficients.
////////////////////////////////////////////////////////////////////////
template <typename T> 
class ISModel
////////////////////////////////////////////////////////////////////////
{
  template <typename U> friend std::ostream& operator<<
    ( std::ostream&, const ISModel<U>& );

  template <typename U> friend class ISVar;

  template <typename U> friend ISVar<U> max
    ( ISVar<U> const&, ISVar<U> const& );
  template <typename U> friend ISVar<U> max
    ( ISVar<U> const&, double const& );
  template <typename U> friend ISVar<U> max
    ( ISVar<U> &&, double const& );
  template <typename U> friend ISVar<U> min
    ( ISVar<U> const&, ISVar<U> const& );
  template <typename U> friend ISVar<U> min
    ( ISVar<U> const&, double const& );
  template <typename U> friend ISVar<U> min
    ( ISVar<U> &&, double const& );

  template <typename U> friend ISVar<U> inv
    ( ISVar<U> const& );
  template <typename U> friend ISVar<U> inv
    ( ISVar<U> && );
  template <typename U> friend ISVar<U> sqr
    ( ISVar<U> const& );
  template <typename U> friend ISVar<U> sqr
    ( ISVar<U> && );
  template <typename U> friend ISVar<U> sqrt
    ( ISVar<U> const& );
  template <typename U> friend ISVar<U> sqrt
    ( ISVar<U> && );
  template <typename U> friend ISVar<U> fabs
    ( ISVar<U> const& );
  template <typename U> friend ISVar<U> fabs
    ( ISVar<U> && );  
  template <typename U> friend ISVar<U> relu
    ( ISVar<U> const& );
  template <typename U> friend ISVar<U> relu
    ( ISVar<U> && );    
  template <typename U> friend ISVar<U> exp
    ( ISVar<U> const& );
  template <typename U> friend ISVar<U> exp
    ( ISVar<U> && );
  template <typename U> friend ISVar<U> log
    ( ISVar<U> const& );
  template <typename U> friend ISVar<U> log
    ( ISVar<U> && );
  template <typename U> friend ISVar<U> xlog
    ( ISVar<U> const& );
  template <typename U> friend ISVar<U> xlog
    ( ISVar<U> && );
  template <typename U> friend ISVar<U> cos
    ( ISVar<U> const& );
  template <typename U> friend ISVar<U> cos
    ( ISVar<U> && );
  template <typename U> friend ISVar<U> sin
    ( ISVar<U> const& );
  template <typename U> friend ISVar<U> sin
    ( ISVar<U> && );
  template <typename U> friend ISVar<U> pow
    ( ISVar<U> const&, int const& n );
  template <typename U> friend ISVar<U> pow
    ( ISVar<U> &&, int const& n );
  template <typename U> friend ISVar<U> tanh
    ( ISVar<U> const& );
  template <typename U> friend ISVar<U> tanh
    ( ISVar<U> && );
  template <typename U> friend ISVar<U> intersect
    ( ISVar<U> const& , U);
  template <typename U> friend ISVar<U> intersect
    ( ISVar<U> && , U);
 private:

  //! @brief Number of partitions (to support up to 2^64 partitions)
  unsigned long long _ndiv;
  //! @brief Number of variables
  unsigned int _nvar;
  //! @brief Whether variables are defined or not
  std::vector<bool> _defvar;
  //! @brief Variable bounds
  std::vector<T> _bndvar;
   // Internal containing the parition sizes (length) of each variable
  std::vector<long double> _psize;

  // Internal Intermediate containing lower bounds of coefficient matrix rows
  mutable std::vector<double> _L1;
  // Internal Intermediate containing lower bounds of coefficient matrix rows
  mutable std::vector<double> _L2;
  // Internal Intermediate containing upper bounds of coefficient matrix rows
  mutable std::vector<double> _U1;
  // Internal Intermediate containing upper bounds of coefficient matrix rows
  mutable std::vector<double> _U2;
  // Internal Intermediate containing reference of coefficient matrix rows
  mutable std::vector<double> _c1;
  // Internal Intermediate containing references of coefficient matrix rows
  mutable std::vector<double> _c2;
  // Internal Intermediate containing radius of coefficient matrix rows
  mutable std::vector<double> _r1;
  // Internal Intermediate containing radius of coefficient matrix rows
  mutable std::vector<double> _r2;

 public:

  //! @brief Constructor of ISM with <a>nvar</a> variables and <a>ndiv</a> partitions
  ISModel
  ( const unsigned int& nvar, const unsigned int& ndiv )
  : _ndiv(ndiv), _nvar(nvar)
  {
    _bndvar.resize( _nvar );
    _defvar.resize( _nvar, false );
    
    _L1.resize( _nvar );
    _L2.resize( _nvar );
    _U1.resize( _nvar );
    _U2.resize( _nvar );
    _c1.resize( _nvar );
    _c2.resize( _nvar );
    _r1.resize( _nvar );
    _r2.resize( _nvar );

    _psize.resize( _nvar , 0. );

  }

  ~ISModel() 
  {
#ifdef FISM_LIFITIME_DEBUG     
    std::cout<< "ISM delated, nvar = " <<_nvar <<std::endl;
#endif    
  }

  //! @brief Retrieve number of variables in ISM
  unsigned nvar
  ()
  const
  { return _nvar; };

  //! @brief Retrieve number of partitions in ISM
  unsigned ndiv
  ()
  const
  { return _ndiv; };

  //! @brief Retrieve partition sizes in ISM
  std::vector<long double> psize
  ()
  const
  { return _psize; };

 //! @brief Options of mc::ISModel
  static struct Options
  {
    //! @brief Constructor
    Options():
      ASYREM_USE(true), DCDEC_USE(true), SCALING_TYPE(FULL), INTERSECTION_USE(true), ENVEL_USE(true), ROOT_MAXIT(100), ROOT_TOL(1e-10), SLOPE_USE(false),SHADOW_USE(false)
      {}

    //! @brief FISModel re-scaling option in binary product
    enum SCALING{
      NONE=0,	  //!< without using scaling
      PARTIAL,	//!< only re-scaling the radius of two mutiplicants(mutipliers) 
      FULL,	    //!< re-scaling the range of two mutiplicants(mutipliers) to [-1,1] 
      ADAPT     //!< adapted re-scaling
    };
    //! @brief Whether to use asymmetric inclusions for convex/concave terms as available
    bool ASYREM_USE;
    //! @brief Whether to use DC decomposition in product rule and composition rule with non-monotonic univariates
    bool DCDEC_USE;
    //! @brief Whether to use re-scaling in binary product, and which type is used 
    SCALING SCALING_TYPE; 
    //! @brief Whether to use intersection
    bool INTERSECTION_USE; 
    //! @brief Whether to use convex/concave envelopes of nonconvex terms as available  
    bool ENVEL_USE;
    //! @brief Maximal number of iterations in root search - Default: 100
    unsigned ROOT_MAXIT;
    //! @brief Termination tolerance in root search - Default: 1e-10
    double ROOT_TOL;
    //! @brief Whether to use slope-based enhancement - Default: false
    bool SLOPE_USE;
    //! @brief Whether to use shadow enhancement - Default: false
    bool SHADOW_USE;    
  } options;

  //! @brief Exceptions of mc::ISModel
  class Exceptions
  {
   public:
    //! @brief Enumeration type for SCModel exception handling
    enum TYPE{
      DIV=1,	 //!< Division by zero scalar
      INV,	     //!< Inverse operation with zero in range
      LOG,	     //!< Log operation with non-positive numbers in range
      SQRT,	     //!< Square-root operation with negative numbers in range
      TAN,	     //!< Tangent operation with (k+1/2)·PI in range
      ROOT,      //!< Error during root search for obtaining the convex/concave envelope of a univariate term
      INTERN=-1, //!< Internal error
      INDEX=-2,  //!< Variable index out of range
      MODEL=-3,	 //!< Operation between variables belonging to different models
      UNDEF=-33  //!< Feature not yet implemented
    };
    //! @brief Constructor for error <a>ierr</a>
    Exceptions( TYPE ierr ) : _ierr( ierr ){}
    //! @brief Error flag
    int ierr(){ return _ierr; }
    //! @brief Error description
    std::string what(){
      switch( _ierr ){
      case DIV:
        return "mc::ISModel\t Division by zero scalar";
      case INV:
        return "mc::ISModel\t Inverse operation with zero in range";
      case LOG:
        return "mc::ISModel\t Log operation with non-positive numbers in range";
      case SQRT:
        return "mc::ISModel\t Square-root operation with negative numbers in range";
      case TAN:
        return "mc::ISModel\t Tangent operation with (k+1/2)·PI in range";
      case INDEX:
        return "mc::ISModel\t Variable index out of range";
      case MODEL:
        return "mc::ISModel\t Operation between variables belonging to different models not permitted";
      case UNDEF:
        return "mc::ISModel\t Feature not yet implemented";
      case INTERN:
      default:
        return "mc::ISModel\t Internal error";
      }
    }
   private:
    TYPE _ierr;
  };

 private:

  //! @brief Return the bound of the ISM variable <a>var</a>
  T _B
  ( std::vector<std::vector<T>> const& mat, const unsigned int rec=0 )
  const
  {
    assert( !mat.empty() );
    T B( 0. );
    for( unsigned int i=0; i<_nvar; i++ ){
      if( mat[i].empty() ) continue;
      if( rec == 3 ){                     //  used for performaing intersection only;
        _L2[i] = _B_innerL( mat[i] );
        _U2[i] = _B_innerU( mat[i] );
      }
      else{                               //  basic function
        T Brow = _B( mat[i] );
        if( rec == 1 ){
          _L1[i] = Op<T>::l( Brow );
          _U1[i] = Op<T>::u( Brow );
        }
        else if( rec == 2 ){
          _L2[i] = Op<T>::l( Brow );
          _U2[i] = Op<T>::u( Brow );
        }
        B += Brow;
      }
    }   
    return B;
  }

  T _B
  ( std::vector<T> const& row )
  const
  {
    assert( !row.empty() );
    T iBnd( row[0] );
    for( unsigned int j=1; j<_ndiv; j++ )
      iBnd = Op<T>::hull( iBnd, row[j] );
    return iBnd;
  }

  double _B_innerL
  ( std::vector<T> const& row )
  const
  {
    assert( !row.empty() );
    double iBndL (Op<T>::l(row[0]));
    for( unsigned int j=1; j<_ndiv; j++ ){
      iBndL = std::max(iBndL,Op<T>::l(row[j]));
    }
    return iBndL;
  }

  double _B_innerU
  ( std::vector<T> const& row )
  const
  {
    assert( !row.empty() );
    double iBndU (Op<T>::u(row[0]));
    for( unsigned int j=1; j<_ndiv; j++ ){
      iBndU = std::min(iBndU,Op<T>::u(row[j]));
    }
    return iBndU;
  }


  template <typename PUNIV>
  void _asym
  ( std::vector<std::vector<T>>& mat, unsigned const& ndep, PUNIV const& f,
    double const& zopt, bool const cvx )
  const;
  template <typename PUNIV>
  void _asym
  ( std::vector<std::vector<T>>& mat, unsigned const& ndep, PUNIV const& f,
    double const& zopt, bool const cvx, T const& bnd )
  const;


  template <typename PUNIV, typename BNUIA>
  void _asym_slope  ( std::vector<std::vector<T>>& mat,std::vector<std::vector<std::vector<long double>>>& slope, const std::vector<long double>& partitionSize, unsigned const& ndep, PUNIV const& f, BNUIA const& fDerv, double const& zopt, bool const cvx )
  const;

  template <typename PUNIV , typename BNUIA>
  void _asym_slope  ( std::vector<std::vector<T>>& mat,std::vector<std::vector<std::vector<long double>>>& slope, const std::vector<long double>& partitionSize, unsigned const& ndep, PUNIV const& f, BNUIA const& fDerv, double const& zopt, bool const cvx, T const& bnd )
  const;

  //  template <typename PUNIV>
  // void _asym_slope_sqr ( std::vector<std::vector<T>>& mat,std::vector<std::vector<std::vector<long double>>>& slope, const std::vector<long double>& partitionSize, unsigned const& ndep, PUNIV const& f, double const& zopt, bool const cvx, T const& bnd )
  // const;

  //  template <typename PUNIV>
  // void _asym_slope_relu ( std::vector<std::vector<T>>& mat,std::vector<std::vector<std::vector<long double>>>& slope, const std::vector<long double>& partitionSize, unsigned const& ndep, PUNIV const& f, double const& zopt, bool const cvx )
  // const;

  void _asym_slope_relu ( std::vector<std::vector<T>>& mat,std::vector<std::vector<std::vector<long double>>>& slope, const std::vector<long double>& partitionSize, unsigned const& ndep)
  const;

  void _asym_slope_relu_ws ( std::vector<std::vector<T>>& mat,std::vector<std::vector<std::vector<long double>>>& slope, const std::vector<long double>& partitionSize, unsigned const& ndep,
                          std::vector<std::vector<std::vector<long double>>>& shadow, std::vector<std::vector<std::vector<long double>>>& shadow_slope)
  const;

  void _asym_slope_relu_shadow ( std::vector<std::vector<T>>& mat,std::vector<std::vector<std::vector<long double>>>& slope, const std::vector<long double>& partitionSize, unsigned const& ndep, std::vector<std::vector<std::vector<long double>>>& shadow, std::vector<std::vector<std::vector<long double>>>& shadow_slope, std::vector<double> & shadow_info )
  const;

  void _add_aggregate_shadow ( std::vector<std::vector<T>>& Amat, const std::vector<std::vector<T>>& Bmat, std::vector<std::vector<std::vector<long double>>>& Aslope, const std::vector<std::vector<std::vector<long double>>>& Bslope,
  std::vector<std::vector<std::vector<long double>>>& Ashadow, const std::vector<std::vector<std::vector<long double>>>& Bshadow, 
  std::vector<std::vector<std::vector<long double>>>& Ashadow_slope, const std::vector<std::vector<std::vector<long double>>>& Bshadow_slope, 
  std::vector<double>& Ashadow_info, const std::vector<double>& Bshadow_info,const std::vector<long double>& partitionSize, unsigned const& ndep,T const& bndB)
  const;

  //  template <typename PUNIV>
  // void _asym_slope_max ( std::vector<std::vector<T>>& mat,std::vector<std::vector<std::vector<long double>>>& slope, const std::vector<long double>& partitionSize, unsigned const& ndep, PUNIV const& f, double const& zopt, bool const cvx )
  // const;


  template <typename PUNIV>
  void _asymDCdecNTC
  ( std::vector<std::vector<T>>& mat, unsigned const& ndep, PUNIV const& f,
    double const& zopt, bool const cvx_ccv )
  const;
  template <typename PUNIV>
  void _asymDCdecNTC
  ( std::vector<std::vector<T>>& mat, unsigned const& ndep, PUNIV const& f,
    double const& zopt, bool const cvx_ccv, T const& bnd )
  const;


  template <typename PUNIV, typename BNUIA>
  void _asymNTCsym_slope  ( std::vector<std::vector<T>>& mat,std::vector<std::vector<std::vector<long double>>>& slope, const std::vector<long double>& partitionSize, unsigned const& ndep, PUNIV const& f, BNUIA const& fDerv, double const& zopt, bool const cvx )
  const;

  template <typename PUNIV , typename BNUIA>
  void _asymNTCsym_slope  ( std::vector<std::vector<T>>& mat,std::vector<std::vector<std::vector<long double>>>& slope, const std::vector<long double>& partitionSize, unsigned const& ndep, PUNIV const& f, BNUIA const& fDerv, double const& zopt, bool const cvx, T const& bnd )
  const;



 template <typename PUNIV>
  void _asymNTCsym
  ( std::vector<std::vector<T>>& mat, unsigned const& ndep, PUNIV const& f)
  const;
  template <typename PUNIV>
  void _asymNTCsym
  ( std::vector<std::vector<T>>& mat, unsigned const& ndep, PUNIV const& f, T const& bnd )
  const;

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
  // when they are concave or convex on the range of the input ISM  


  // @unsure: there is a similar algorithm which can refining the multiplicative result, 
  // by bounding x_iy_j+x_jy_i in 0.5*(x_i^2+y_i^2)*[-1,1] + 0.5*(x_j^2+y_j^2)*[-1,1]

  static std::ostream& _dispmat
  ( std::vector<std::vector<T>> const& mat, unsigned const len=3, std::ostream& out=std::cout );
  
  std::ostream& _dispvar
  ( std::vector<std::vector<T>> const& mat, unsigned const& ndep, const int& opt=0,
    std::ostream& out=std::cout )
  const;

  template <typename PUNIV>
  double _goldsect
  ( const double xL, const double xU, PUNIV const& f, const double TOL, const unsigned MAXIT )
  const;
  
  template <typename PUNIV>
  double _goldsect_iter
  ( unsigned& iter, const double a, const double fa, const double b,
    const double fb, const double c, const double fc, PUNIV const& f,
    const double TOL, const unsigned MAXIT )
  const;
};

template <typename T> inline
typename ISModel<T>::Options ISModel<T>::options;

template <typename T>
std::ostream& operator<<
( std::ostream& out, const ISModel<T>& mod)
{
  out << std::endl;
  out << "ISM settings:" << std::endl;
  out << "   " << "no. variables:  " << mod._nvar << std::endl;
  out << "   " << "no. partitions: " << mod._ndiv << std::endl;
  out << "   " << "variable bounds: " << std::endl;
  for( unsigned int i=0; i<mod._nvar; i++ )
    out << "       " << i << ": " << (mod._defvar[i]? mod._bndvar[i]: "-") << std::endl;
  out << std::endl;
  return out;
}

template <typename T>
inline
std::ostream& ISModel<T>::_dispmat
( std::vector<std::vector<T>> const& mat, unsigned const len, std::ostream& out )
{
  if( mat.empty() ) return out;
  unsigned irow = 0;
  for( auto&& row : mat ){
    if( row.empty() ) continue;
    out << std::right << std::setw(5) << irow++ <<": ";
    unsigned icol = 0;
    for( auto&&el : row ){
      if( icol && !(icol%len) ) out << std::endl << "       ";
      out << std::setw(0) << el;
      icol++;
    }
    out << std::endl;
  }
  return out;
}

template <typename T>
inline
std::ostream& ISModel<T>::_dispvar
( std::vector<std::vector<T>> const& mat, unsigned const& ndep, const int& opt,
  std::ostream& out )
const
{
  if( ndep > 2 ) return out;
  assert( !mat.empty() );
  out << std::scientific << std::setprecision(5) << std::right;
  
  if( ndep == 1 ){
    for( unsigned int i=0; i<_nvar; i++ ){
      if( mat[i].empty() ) continue;
      assert( _defvar[i] );
      double l = Op<T>::l(_bndvar[i]);
      double h = Op<T>::diam(_bndvar[i]) / (double)_ndiv;
      for( unsigned int j=0; j<_ndiv; j++, l+=h ){
        if( opt == 0 ){
          out << std::setw(14) << l   << std::setw(14) << Op<T>::l(mat[i][j]) << std::endl;
          out << std::setw(14) << l   << std::setw(14) << Op<T>::u(mat[i][j]) << std::endl;
          out << std::setw(14) << l+h << std::setw(14) << Op<T>::u(mat[i][j]) << std::endl;
          out << std::setw(14) << l+h << std::setw(14) << Op<T>::l(mat[i][j]) << std::endl;
          out << std::setw(14) << l   << std::setw(14) << Op<T>::l(mat[i][j]) << std::endl;
          out << std::endl;
        }
        else if( opt > 0 ){
          out << std::setw(14) << l   << std::setw(14) << Op<T>::u(mat[i][j]) << std::endl;
          out << std::setw(14) << l+h << std::setw(14) << Op<T>::u(mat[i][j]) << std::endl;
          out << std::endl;
        }
        else{
          out << std::setw(14) << l   << std::setw(14) << Op<T>::l(mat[i][j]) << std::endl;
          out << std::setw(14) << l+h << std::setw(14) << Op<T>::l(mat[i][j]) << std::endl;
          out << std::endl;
        }
      }
      break;
    }
  }
  
  else if( ndep == 2 ){
    for( unsigned int i1=0; i1<_nvar; i1++ ){
      if( mat[i1].empty() ) continue;
      assert( _defvar[i1] );
      //std::cout << "var1: " << _bndvar[i1] << std::endl;
      double l1 = Op<T>::l(_bndvar[i1]);
      double h1 = Op<T>::diam(_bndvar[i1]) / (double)_ndiv;
      for( unsigned int j1=0; j1<_ndiv; j1++, l1+=h1 ){
        for( unsigned int i2=i1+1; i2<_nvar; i2++ ){
          if( mat[i2].empty() ) continue;
          assert( _defvar[i2] );
          //std::cout << "var2: " << _bndvar[i2] << std::endl;
          double l2 = Op<T>::l(_bndvar[i2]);
          double h2 = Op<T>::diam(_bndvar[i2]) / (double)_ndiv;
          for( unsigned int j2=0; j2<_ndiv; j2++, l2+=h2 ){
            if( opt == 0 ){
              out << std::setw(14) << l1    << std::setw(14) << l2    << std::setw(14) << Op<T>::l(mat[i1][j1]+mat[i2][j2]) << std::endl;
              out << std::setw(14) << l1    << std::setw(14) << l2    << std::setw(14) << Op<T>::u(mat[i1][j1]+mat[i2][j2]) << std::endl;
              out << std::setw(14) << l1+h1 << std::setw(14) << l2    << std::setw(14) << Op<T>::u(mat[i1][j1]+mat[i2][j2]) << std::endl;
              out << std::setw(14) << l1+h1 << std::setw(14) << l2    << std::setw(14) << Op<T>::l(mat[i1][j1]+mat[i2][j2]) << std::endl;
              out << std::setw(14) << l1    << std::setw(14) << l2    << std::setw(14) << Op<T>::l(mat[i1][j1]+mat[i2][j2]) << std::endl;
              out << std::endl;
              out << std::setw(14) << l1    << std::setw(14) << l2+h2 << std::setw(14) << Op<T>::l(mat[i1][j1]+mat[i2][j2]) << std::endl;
              out << std::setw(14) << l1    << std::setw(14) << l2+h2 << std::setw(14) << Op<T>::u(mat[i1][j1]+mat[i2][j2]) << std::endl;
              out << std::setw(14) << l1+h1 << std::setw(14) << l2+h2 << std::setw(14) << Op<T>::u(mat[i1][j1]+mat[i2][j2]) << std::endl;
              out << std::setw(14) << l1+h1 << std::setw(14) << l2+h2 << std::setw(14) << Op<T>::l(mat[i1][j1]+mat[i2][j2]) << std::endl;
              out << std::setw(14) << l1    << std::setw(14) << l2+h2 << std::setw(14) << Op<T>::l(mat[i1][j1]+mat[i2][j2]) << std::endl;
              out << std::endl;
              out << std::setw(14) << l1    << std::setw(14) << l2    << std::setw(14) << Op<T>::l(mat[i1][j1]+mat[i2][j2]) << std::endl;
              out << std::setw(14) << l1    << std::setw(14) << l2    << std::setw(14) << Op<T>::u(mat[i1][j1]+mat[i2][j2]) << std::endl;
              out << std::setw(14) << l1    << std::setw(14) << l2+h2 << std::setw(14) << Op<T>::u(mat[i1][j1]+mat[i2][j2]) << std::endl;
              out << std::setw(14) << l1    << std::setw(14) << l2+h2 << std::setw(14) << Op<T>::l(mat[i1][j1]+mat[i2][j2]) << std::endl;
              out << std::setw(14) << l1    << std::setw(14) << l2    << std::setw(14) << Op<T>::l(mat[i1][j1]+mat[i2][j2]) << std::endl;
              out << std::endl;
              out << std::setw(14) << l1+h1 << std::setw(14) << l2    << std::setw(14) << Op<T>::l(mat[i1][j1]+mat[i2][j2]) << std::endl;
              out << std::setw(14) << l1+h1 << std::setw(14) << l2    << std::setw(14) << Op<T>::u(mat[i1][j1]+mat[i2][j2]) << std::endl;
              out << std::setw(14) << l1+h1 << std::setw(14) << l2+h2 << std::setw(14) << Op<T>::u(mat[i1][j1]+mat[i2][j2]) << std::endl;
              out << std::setw(14) << l1+h1 << std::setw(14) << l2+h2 << std::setw(14) << Op<T>::l(mat[i1][j1]+mat[i2][j2]) << std::endl;
              out << std::setw(14) << l1+h1 << std::setw(14) << l2    << std::setw(14) << Op<T>::l(mat[i1][j1]+mat[i2][j2]) << std::endl;
              out << std::endl << std::endl;
            }
            else if( opt > 0 ){
              out << std::setw(14) << l1    << std::setw(14) << l2    << std::setw(14) << Op<T>::u(mat[i1][j1]+mat[i2][j2]) << std::endl;
              out << std::setw(14) << l1+h1 << std::setw(14) << l2    << std::setw(14) << Op<T>::u(mat[i1][j1]+mat[i2][j2]) << std::endl;
              out << std::setw(14) << l1+h1 << std::setw(14) << l2+h2 << std::setw(14) << Op<T>::u(mat[i1][j1]+mat[i2][j2]) << std::endl;
              out << std::setw(14) << l1    << std::setw(14) << l2+h2 << std::setw(14) << Op<T>::u(mat[i1][j1]+mat[i2][j2]) << std::endl;
              out << std::setw(14) << l1    << std::setw(14) << l2    << std::setw(14) << Op<T>::u(mat[i1][j1]+mat[i2][j2]) << std::endl;
              out << std::endl << std::endl;            
            }
            else{
              out << std::setw(14) << l1    << std::setw(14) << l2    << std::setw(14) << Op<T>::l(mat[i1][j1]+mat[i2][j2]) << std::endl;
              out << std::setw(14) << l1+h1 << std::setw(14) << l2    << std::setw(14) << Op<T>::l(mat[i1][j1]+mat[i2][j2]) << std::endl;
              out << std::setw(14) << l1+h1 << std::setw(14) << l2+h2 << std::setw(14) << Op<T>::l(mat[i1][j1]+mat[i2][j2]) << std::endl;
              out << std::setw(14) << l1    << std::setw(14) << l2+h2 << std::setw(14) << Op<T>::l(mat[i1][j1]+mat[i2][j2]) << std::endl;
              out << std::setw(14) << l1    << std::setw(14) << l2    << std::setw(14) << Op<T>::l(mat[i1][j1]+mat[i2][j2]) << std::endl;
              out << std::endl << std::endl;            
            }
          }
          break;
        }
      }
      break;
    }
  }

  return out;
}






template <typename T>
template <typename PUNIV, typename BNUIA>
inline
void
ISModel<T>::_asym_slope
( std::vector<std::vector<T>>& mat,std::vector<std::vector<std::vector<long double>>>& slope, const std::vector<long double>& partitionSize,
    unsigned const& ndep, PUNIV const& f, BNUIA const& fDerv, double const& zopt, bool const cvx  )
const
{
  assert( !mat.empty() );
  T bnd = _B( mat, 1 );
  return _asym_slope( mat, slope, partitionSize, ndep, f, fDerv, zopt, cvx, bnd );
}




template <typename T>
template <typename PUNIV>
inline
void
ISModel<T>::_asym
( std::vector<std::vector<T>>& mat, unsigned const& ndep, PUNIV const& f,
  double const& zopt, bool const cvx )
const
{
  assert( !mat.empty() );
  T bnd = _B( mat, 1 );
  return _asym( mat, ndep, f, zopt, cvx, bnd );
}

template <typename T>
template <typename PUNIV>
inline
void
ISModel<T>::_asym
( std::vector<std::vector<T>>& mat, unsigned const& ndep, PUNIV const& f,
  double const& zopt, bool const cvx, T const& bnd )
const
{
  // anchor points
  int imid( -1 );
  mid( Op<T>::l(bnd), Op<T>::u(bnd), zopt, imid );
  double sum_r1( 0. );
  double C1( 0. ), C2(0.);
  for( unsigned int i=0; i<_nvar; i++ ){
    if( mat[i].empty() ) continue;
    switch( imid ){
      case ICONV: _c1[i] = _L1[i], C1 += _c1[i]; break;
      case ICONC: _c1[i] = _U1[i], C1 += _c1[i]; break;
      case ICUT:  _c1[i] = _L1[i]; C1 += _c1[i];
                  _c2[i] = _U1[i]; C2 += _c2[i]; break;
    }
    _r1[i] = ( _U1[i] - _L1[i] );
    sum_r1 += _r1[i];
  }

  double fopt = f( imid == ICUT? zopt: C1 );
  double fopt_over_ndep = fopt / ndep;

  for( unsigned int i=0; i<_nvar; i++ ){
    if( mat[i].empty() ) continue;   
    if( isequal( _r1[i], 0. ) ){
      for( unsigned int j=0; j<_ndiv; j++ )
        mat[i][j] = fopt_over_ndep;
      continue;
    }
    else{
      double scal_r1 = _r1[i] / sum_r1; 
      for( unsigned int j=0; j<_ndiv; j++ ){
        if( cvx ){
          T Zu = ( mat[i][j] - _c1[i] ) / _r1[i] * sum_r1 + C1;
          double Du = scal_r1 * ( std::max( f( Op<T>::l(Zu) ), f( Op<T>::u(Zu) ) ) - fopt ) + fopt_over_ndep;
          if( imid != ICUT ){
            T Zl = mat[i][j] - _c1[i] + C1;
            double El =  f( mid( Op<T>::l(Zl), Op<T>::u(Zl), zopt ) ) - fopt_over_ndep * (ndep-1.);//std::min( f( Op<T>::l(Zl) ), f( Op<T>::u(Zl) ) ) - fopt_over_ndep * (ndep-1.);
            mat[i][j] = T( El, Du );
#ifdef MC__USE_FILIB
            if(mat[i][j].isEmpty()){
                if(std::fabs(El-Du) <= MC__ISM_COMPUTATION_TOL)
                    mat[i][j] = T( El );}
#ifdef FILIB__COMPUTATION_DEBUG   
            if(mat[i][j].isEmpty())
              std::cout<<"Case 1: " <<std::setprecision(18)<<El<<" > "<<Du<<std::endl;
#endif
#endif              
          }
          else if( !options.DCDEC_USE ){
            mat[i][j] = T( fopt_over_ndep, Du );
#ifdef MC__USE_FILIB
            if(mat[i][j].isEmpty()){
                if(std::fabs(fopt_over_ndep - Du) <= MC__ISM_COMPUTATION_TOL)
                    mat[i][j] = T( fopt_over_ndep );}
#ifdef FILIB__COMPUTATION_DEBUG   
            if(mat[i][j].isEmpty())
              std::cout<<"Case 2: " <<std::setprecision(18)<<fopt_over_ndep<<" > "<<Du<<std::endl;             
#endif 
#endif            
          }
          else{
            auto const& flinc = [=]( const T& x ){ return Op<T>::l(x) < zopt? 0.: f( Op<T>::l(x) ) - fopt; };
            auto const& fldec = [=]( const T& x ){ return Op<T>::u(x) > zopt? 0.: f( Op<T>::u(x) ) - fopt; };
            double El = flinc( mat[i][j] - _c1[i] + C1 ) + fldec( mat[i][j] - _c2[i] + C2 ) + fopt_over_ndep;
            mat[i][j] = T( El, Du );
#ifdef MC__USE_FILIB
            if(mat[i][j].isEmpty()){
                if(std::fabs(El - Du) <= MC__ISM_COMPUTATION_TOL)
                    mat[i][j] = T( El );}
#ifdef FILIB__COMPUTATION_DEBUG   
            if(mat[i][j].isEmpty())
              std::cout<<"Case 3: " <<std::setprecision(18)<<El<<" > "<<Du<<std::endl;             
#endif
#endif   
          }
        }
        else{
          T Zl = ( mat[i][j] - _c1[i] ) / _r1[i] * sum_r1 + C1;
          double Dl = scal_r1 * ( std::min( f( Op<T>::l(Zl) ), f( Op<T>::u(Zl) ) ) - fopt ) + fopt_over_ndep;
          if( imid != ICUT ){
            T Zu = mat[i][j] - _c1[i] + C1;
            double Eu = f( mid( Op<T>::l(Zu), Op<T>::u(Zu), zopt ) ) - fopt_over_ndep * (ndep-1.);
            mat[i][j] = T( Dl, Eu );
#ifdef MC__USE_FILIB
            if(mat[i][j].isEmpty()){
                if(std::fabs(Dl - Eu) <= MC__ISM_COMPUTATION_TOL)
                    mat[i][j] = T( Eu );}
#ifdef FILIB__COMPUTATION_DEBUG   
            if(mat[i][j].isEmpty())
              std::cout<<"Case 4: " <<std::setprecision(18)<<Dl<<" > "<<Eu<<std::endl;              
#endif
#endif   
          }
          else if( !options.DCDEC_USE ){
            mat[i][j] = T( Dl, fopt_over_ndep );
#ifdef MC__USE_FILIB
            if(mat[i][j].isEmpty()){
                if(std::fabs(Dl - fopt_over_ndep) <= MC__ISM_COMPUTATION_TOL)
                    mat[i][j] = T( Dl );}
#ifdef FILIB__COMPUTATION_DEBUG   
            if(mat[i][j].isEmpty())
              std::cout<<"Case 5: " <<std::setprecision(18)<<Dl<<" > "<<fopt_over_ndep<<std::endl;              
#endif   
#endif
          }
          else{
            auto const& fuinc = [=]( const T& x ){ return Op<T>::u(x) > zopt? 0.: f( Op<T>::u(x) ) - fopt; };
            auto const& fudec = [=]( const T& x ){ return Op<T>::l(x) < zopt? 0.: f( Op<T>::l(x) ) - fopt; };
            double Eu = fuinc( mat[i][j] - _c2[i] + C2 ) + fudec( mat[i][j] - _c1[i] + C1 ) + fopt_over_ndep;
            mat[i][j] = T( Dl, Eu );
#ifdef MC__USE_FILIB
            if(mat[i][j].isEmpty()){
                if(std::fabs(Dl - Eu) <= MC__ISM_COMPUTATION_TOL)
                    mat[i][j] = T( Eu );}
#ifdef FILIB__COMPUTATION_DEBUG   
            if(mat[i][j].isEmpty())
              std::cout<<"Case 6: " <<std::setprecision(18)<<Dl<<" > "<<Eu<<std::endl;  
#endif
#endif                
          }
        }
#ifdef MC__USE_FILIB        
#ifdef FILIB__COMPUTATION_DEBUG        
        if(mat[i][j].isEmpty())
          std::cout<<"_asym:" <<i<<" , "<<j<<std::endl;  
#endif   
#endif
      }
    }
  }
}


template <typename T>
template <typename PUNIV, typename BNUIA>
inline
void
ISModel<T>::_asym_slope
( std::vector<std::vector<T>>& mat, std::vector<std::vector<std::vector<long double>>>& slope, const std::vector<long double>& partitionSize,
    unsigned const& ndep, PUNIV const& f, BNUIA const& fDerv, double const& zopt, bool const cvx, T const& bnd )
const
{
  // anchor points
  int imid( -1 );
  mid( Op<T>::l(bnd), Op<T>::u(bnd), zopt, imid );
  long double sum_r1( Op<T>::u(bnd) - Op<T>::l(bnd) );
  double C1( 0. ), C2(0.);
  for( unsigned int i=0; i<_nvar; i++ ){
    if( mat[i].empty() ) continue;
    switch( imid ){
      case ICONV: _c1[i] = _L1[i], C1 += _c1[i]; break;
      case ICONC: _c1[i] = _U1[i], C1 += _c1[i]; break;
      case ICUT:  _c1[i] = _L1[i]; C1 += _c1[i];
                  _c2[i] = _U1[i]; C2 += _c2[i]; break;
    }
    _r1[i] = ( _U1[i] - _L1[i] );
  }


  double fopt = f( imid == ICUT? zopt: C1 );
  double fopt_over_ndep = fopt / ndep;

  for( unsigned int i=0; i<_nvar; i++ ){
    if( mat[i].empty() ) continue;   
    if( isequal( _r1[i], 0. ) ){
      for( unsigned int j=0; j<_ndiv; j++ ){
        mat[i][j] = T(fopt_over_ndep);
        slope[i][j][0] = 0.;
        slope[i][j][1] = 0.;
      }
      continue;
    }
    else{
      long double scal_r1 = _r1[i] / sum_r1; 
      for( unsigned int j=0; j<_ndiv; j++ ){
        if( cvx ){
          switch( imid ){
            case ICONV: {
              long double zU = Op<T>::u(mat[i][j]);
              long double delta_u = std::fabs(slope[i][j][1]*partitionSize[i]);
              long double Du = f(( zU - _c1[i] ) / _r1[i] * sum_r1 + C1);
              long double slope1_to_be_multiplied = fDerv(( zU - delta_u - _c1[i] ) / _r1[i] * sum_r1 + C1);
              if (delta_u > 1e2 * MC__ISM_COMPUTATION_TOL)
                slope1_to_be_multiplied = (Du - f(( zU - delta_u - _c1[i] ) / _r1[i] * sum_r1 + C1))*(scal_r1/delta_u);
              slope[i][j][1] = slope1_to_be_multiplied*slope[i][j][1];
              Du = scal_r1 * (Du  - fopt ) + fopt_over_ndep;

              long double zL = Op<T>::l(mat[i][j]);
              long double El = f( zL - _c1[i] + C1 ) - fopt_over_ndep * (ndep-1.);
              long double slope0_to_be_multiplied = fDerv( zL- _c1[i] + C1 );
              slope[i][j][0] = slope0_to_be_multiplied*slope[i][j][0];

              if(Du - El <= MC__ISM_COMPUTATION_TOL){
                mat[i][j] = T(std::min(Du,El),std::max(Du,El)); 
                slope[i][j][0] = 0.;
                slope[i][j][1] = 0.;                
              }      
              else
                mat[i][j] = T( El, Du );

            }
            break;
            case ICONC: {
              long double zL = Op<T>::l(mat[i][j]);
              long double delta_l = std::fabs(slope[i][j][0]*partitionSize[i]);
              long double Du = f(( zL - _c1[i] ) / _r1[i] * sum_r1 + C1);
              long double slope1_to_be_set = fDerv(( zL + delta_l - _c1[i] ) / _r1[i] * sum_r1 + C1);
              if (delta_l > 1e2 * MC__ISM_COMPUTATION_TOL)
                slope1_to_be_set = (f(( zL + delta_l - _c1[i] ) / _r1[i] * sum_r1 + C1) - Du)*(scal_r1/delta_l);
              slope1_to_be_set  = slope1_to_be_set*slope[i][j][0];
              Du = scal_r1 * (Du  - fopt ) + fopt_over_ndep;

              long double zU = Op<T>::u(mat[i][j]);
              long double El = f( zU - _c1[i] + C1 ) - fopt_over_ndep * (ndep-1.);
              long double slope0_to_be_multiplied = fDerv( zU - _c1[i] + C1 );
              slope[i][j][0] = slope0_to_be_multiplied*slope[i][j][1];
              slope[i][j][1] = slope1_to_be_set;

              if(Du - El <= MC__ISM_COMPUTATION_TOL){
                slope[i][j][0] = 0.;
                slope[i][j][1] = 0.; 
                mat[i][j] = T(std::min(Du,El),std::max(Du,El));     
              }
              else
                mat[i][j] = T( El, Du );

            }
            break;
            case ICUT: {
              // f(x) = fInc(x) + fDec(x) + fopt = fInc(x) + fDec(x)  - sum scal_r1 * fopt + fopt_over_ndep
              auto const& fInc     = [=]( const double& x ){ return x < zopt? 0.: f( x ) - fopt; };
              auto const& fDec     = [=]( const double& x ){ return x > zopt? 0.: f( x ) - fopt; };      
              auto const& fIncDerv = [=]( const double& x ){ return x < zopt? 0.: fDerv( x ); };
              auto const& fDecDerv = [=]( const double& x ){ return x > zopt? 0.: fDerv( x ); };  

              long double zU = Op<T>::u(mat[i][j]);
              long double delta_u = std::fabs(slope[i][j][1]*partitionSize[i]);
              long double DuInc = fInc(( zU - _c2[i] ) / _r1[i] * sum_r1 + C2);
              long double slope1Inc_to_be_set = fIncDerv(( zU - delta_u - _c2[i] ) / _r1[i] * sum_r1 + C2);
              if (delta_u > 1e2 * MC__ISM_COMPUTATION_TOL)
                slope1Inc_to_be_set = (DuInc - fInc(( zU - delta_u - _c2[i] ) / _r1[i] * sum_r1 + C2))*(scal_r1/delta_u);
              slope1Inc_to_be_set = slope1Inc_to_be_set*slope[i][j][1];
              DuInc = scal_r1 * DuInc;    // For sum_i DuInc > fInc(x), its fIncOpt = 0;

              long double zL = Op<T>::l(mat[i][j]);
              long double delta_l = std::fabs(slope[i][j][0]*partitionSize[i]);
              long double DuDec = fDec(( zL - _c1[i] ) / _r1[i] * sum_r1 + C1);
              long double slope1Dec_to_be_set = fDecDerv(( zL + delta_l - _c1[i] ) / _r1[i] * sum_r1 + C1);
              if (delta_l > 1e2 * MC__ISM_COMPUTATION_TOL)
                slope1Dec_to_be_set = (fDec(( zL + delta_l - _c1[i] ) / _r1[i] * sum_r1 + C1) - DuDec)*(scal_r1/delta_l);
              slope1Dec_to_be_set  = slope1Dec_to_be_set*slope[i][j][0];
              DuDec = scal_r1 * DuDec;  // For sum_i DuDec > fDec(x), its fDecOpt = 0;
                            
              // The original one is: 
              // T Zu = ( mat[i][j] - _c1[i] ) / _r1[i] * sum_r1 + C1; 
              // double Du = scal_r1 * ( std::max( f( Op<T>::l(Zu) ), f( Op<T>::u(Zu) ) ) - fopt ) + fopt_over_ndep;
              // which is the same as Du = std::max(DuInc, DuDec) - scal_r1 * fopt + fopt_over_ndep
              long double slope1_to_be_set = slope1Inc_to_be_set + slope1Dec_to_be_set;    // Note that we should not directly assign values to slope[i][j][1] = 0. as we will use them later
              long double upperbound_to_be_set = DuDec + DuInc;
              if((slope1Inc_to_be_set < 0.) != (slope1Dec_to_be_set < 0.)){   
                upperbound_to_be_set = upperbound_to_be_set - std::min(std::fabs(slope1Inc_to_be_set),std::fabs(slope1Dec_to_be_set))*partitionSize[i];
              }

              long double upperboundCandidate = std::max(DuInc, DuDec);
              if ( upperbound_to_be_set > upperboundCandidate + MC__ISM_COMPUTATION_TOL){
                  upperbound_to_be_set = upperboundCandidate;
                  slope1_to_be_set = 0.;
              }
                     
              // if (upperbound_to_be_set > Du){
              //   upperbound_to_be_set = Du;
              //   slope[i][j][1] = 0.;
              // }
              
              // The original one is 
              // double El = flinc( mat[i][j] - _c1[i] + C1 ) + fldec( mat[i][j] - _c2[i] + C2 ) + fopt_over_ndep;          
              long double ElInc = fInc( zL - _c1[i] + C1 ); // For sum_i ElInc < fInc(x), its fIncOpt = 0;
              long double slope0Inc_to_be_set = fIncDerv( zL- _c1[i] + C1 );
              slope0Inc_to_be_set = slope0Inc_to_be_set*slope[i][j][0];

              long double ElDec = fDec( zU - _c2[i] + C2 ); // For sum_i ElDec < fDec(x), its fDecOpt = 0;
              long double slope0Dec_to_be_set = fDecDerv( zU - _c2[i] + C2 );
              slope0Dec_to_be_set = slope0Dec_to_be_set*slope[i][j][1];

              long double slope0_to_be_set = slope0Inc_to_be_set + slope0Dec_to_be_set;
              
              long double lowerbound_to_be_set = ElDec + ElInc;
              if((slope0Inc_to_be_set < 0.) != (slope0Dec_to_be_set < 0.)){            
               lowerbound_to_be_set = lowerbound_to_be_set + std::min(std::fabs(slope0Inc_to_be_set),std::fabs(slope0Dec_to_be_set))*partitionSize[i];
              }

              if ( lowerbound_to_be_set <  MC__ISM_COMPUTATION_TOL){
                  lowerbound_to_be_set = 0.;
                  slope0_to_be_set = 0.;                   
              }

               
              // if (lowerbound_to_be_set < El){
              //   lowerbound_to_be_set = Eul;
              //   slope[i][j][0] = 0.;
              // }


              if(upperbound_to_be_set - lowerbound_to_be_set <= MC__ISM_COMPUTATION_TOL){
                mat[i][j] = T(std::min(lowerbound_to_be_set,upperbound_to_be_set),std::max(lowerbound_to_be_set,upperbound_to_be_set))  + fopt_over_ndep;                               
                slope[i][j][0] = 0.;
                slope[i][j][1] = 0.;
              }
              else
                mat[i][j] = T(lowerbound_to_be_set , upperbound_to_be_set ) + fopt_over_ndep;     
                slope[i][j][0] = slope0_to_be_set;
                slope[i][j][1] = slope1_to_be_set;                           
            }
            break;
          }          

        }
        else{
          switch( imid ){
            case ICONV: {   // decreasing concave function
              long double zU = Op<T>::u(mat[i][j]);
              long double delta_u = std::fabs(slope[i][j][1]*partitionSize[i]);
              long double Dl = f(( zU - _c1[i] ) / _r1[i] * sum_r1 + C1);
              long double slope0_to_be_set = fDerv(( zU - delta_u - _c1[i] ) / _r1[i] * sum_r1 + C1);
              if (delta_u > 1e2 * MC__ISM_COMPUTATION_TOL)
                slope0_to_be_set = (Dl - f(( zU - delta_u - _c1[i] ) / _r1[i] * sum_r1 + C1))*(scal_r1/delta_u);
              slope0_to_be_set = slope0_to_be_set*slope[i][j][1];
              Dl = scal_r1 * (Dl  - fopt ) + fopt_over_ndep;

              long double zL = Op<T>::l(mat[i][j]);
              long double Eu = f( zL - _c1[i] + C1 ) - fopt_over_ndep * (ndep-1.);
              long double slope1_to_be_set = fDerv( zL- _c1[i] + C1 );
              slope1_to_be_set = slope1_to_be_set*slope[i][j][0];

              if(Dl - Eu <= MC__ISM_COMPUTATION_TOL){
                mat[i][j] = T(std::min(Dl,Eu),std::max(Dl,Eu)); 
                slope[i][j][0] = 0.;
                slope[i][j][1] = 0.;                
              }      
              else
                mat[i][j] = T( Dl, Eu );
                slope[i][j][0] = slope0_to_be_set;
                slope[i][j][1] = slope1_to_be_set;    
            }
            break;
            case ICONC: {   // increasing concave function
              long double zL = Op<T>::l(mat[i][j]);
              long double delta_l = std::fabs(slope[i][j][0]*partitionSize[i]);
              long double Dl = f(( zL - _c1[i] ) / _r1[i] * sum_r1 + C1);
              long double slope0_to_be_set = fDerv(( zL + delta_l - _c1[i] ) / _r1[i] * sum_r1 + C1);
              if (delta_l > 1e2 * MC__ISM_COMPUTATION_TOL)
                slope0_to_be_set = (f(( zL + delta_l - _c1[i] ) / _r1[i] * sum_r1 + C1) - Dl)*(scal_r1/delta_l);
              slope0_to_be_set  = slope0_to_be_set*slope[i][j][0];
              Dl = scal_r1 * (Dl  - fopt ) + fopt_over_ndep;

              long double zU = Op<T>::u(mat[i][j]);
              long double Eu = f( zU - _c1[i] + C1 ) - fopt_over_ndep * (ndep-1.);
              long double slope1_to_be_set = fDerv( zU - _c1[i] + C1 );
              slope1_to_be_set = slope1_to_be_set*slope[i][j][1];

              if(Dl - Eu <= MC__ISM_COMPUTATION_TOL){
                slope[i][j][0] = 0.;
                slope[i][j][1] = 0.; 
                mat[i][j] = T(std::min(Dl,Eu),std::max(Dl,Eu));     
              }
              else
                mat[i][j] = T( Dl, Eu );
                slope[i][j][0] = slope0_to_be_set;
                slope[i][j][1] = slope1_to_be_set;    
            }
            break;
            case ICUT: {
              std::cout << "   ERROR: slope enhancer for the scenario: concave-cut has not yet implemented " << std::endl;
              throw typename ISModel<T>::Exceptions( ISModel<T>::Exceptions::UNDEF ); 
              // // f(x) = fInc(x) + fDec(x) + fopt = fInc(x) + fDec(x)  - sum scal_r1 * fopt + fopt_over_ndep
              // auto const& fInc     = [=]( const double& x ){ return x < zopt? 0.: f( x ) - fopt; };
              // auto const& fDec     = [=]( const double& x ){ return x > zopt? 0.: f( x ) - fopt; };      
              // auto const& fIncDerv = [=]( const double& x ){ return x < zopt? 0.: fDerv( x ); };
              // auto const& fDecDerv = [=]( const double& x ){ return x > zopt? 0.: fDerv( x ); };  

              // long double zU = Op<T>::u(mat[i][j]);
              // long double delta_u = std::fabs(slope[i][j][1]*partitionSize[i]);
              // long double DuInc = fInc(( zU - _c1[i] ) / _r1[i] * sum_r1 + C1);
              // long double slope1Inc_to_be_set = fIncDerv(( zU - delta_u - _c1[i] ) / _r1[i] * sum_r1 + C1);
              // if (delta_u > 1e2 * MC__ISM_COMPUTATION_TOL)
              //   slope1Inc_to_be_set = (DuInc - fInc(( zU - delta_u - _c1[i] ) / _r1[i] * sum_r1 + C1))*(scal_r1/delta_u);
              // slope1Inc_to_be_set = slope1Inc_to_be_set*slope[i][j][1];
              // DuInc = scal_r1 * DuInc;    // For sum_i DuInc > fInc(x), its fIncOpt = 0;

              // long double zL = Op<T>::l(mat[i][j]);
              // long double delta_l = std::fabs(slope[i][j][0]*partitionSize[i]);
              // long double DuDec = fDec(( zL - _c1[i] ) / _r1[i] * sum_r1 + C1);
              // long double slope1Dec_to_be_set = fDecDerv(( zL + delta_l - _c1[i] ) / _r1[i] * sum_r1 + C1);
              // if (delta_l > 1e2 * MC__ISM_COMPUTATION_TOL)
              //   slope1Dec_to_be_set = (fDec(( zL + delta_l - _c1[i] ) / _r1[i] * sum_r1 + C1) - DuDec)*(scal_r1/delta_l);
              // slope1Dec_to_be_set  = slope1Dec_to_be_set*slope[i][j][0];
              // DuDec = scal_r1 * DuDec;  // For sum_i DuDec > fDec(x), its fDecOpt = 0;
                            
              // // The original one is: 
              // // T Zu = ( mat[i][j] - _c1[i] ) / _r1[i] * sum_r1 + C1; 
              // // double Du = scal_r1 * ( std::max( f( Op<T>::l(Zu) ), f( Op<T>::u(Zu) ) ) - fopt ) + fopt_over_ndep;
              // // which is the same as Du = std::max(DuInc, DuDec) - scal_r1 * fopt + fopt_over_ndep

              // slope[i][j][1] = slope1Inc_to_be_set + slope1Dec_to_be_set;
              // long double upperbound_to_be_set = DuDec + DuInc - std::min(std::fabs(slope1Inc_to_be_set),std::fabs(slope1Dec_to_be_set))*partitionSize[i];
              // // if (upperbound_to_be_set > Du){
              // //   upperbound_to_be_set = Du;
              // //   slope[i][j][1] = 0.;
              // // }
              
              // // The original one is 
              // // double El = flinc( mat[i][j] - _c1[i] + C1 ) + fldec( mat[i][j] - _c2[i] + C2 ) + fopt_over_ndep;

              
              // long double ElInc = fInc( zL - _c1[i] + C1 ); // For sum_i ElInc < fInc(x), its fIncOpt = 0;
              // long double slope0Inc_to_be_set = fIncDerv( zL- _c1[i] + C1 );
              // slope0Inc_to_be_set = slope0Inc_to_be_set*slope[i][j][0];

              // long double ElDec = fDec( zU - _c1[i] + C1 ); // For sum_i ElDec < fDec(x), its fDecOpt = 0;
              // long double slope0Dec_to_be_set = fDecDerv( zU - _c1[i] + C1 );
              // slope0Dec_to_be_set = slope0Dec_to_be_set*slope[i][j][1];

              // slope[i][j][0] = slope0Inc_to_be_set + slope0Dec_to_be_set;
              // long double lowerbound_to_be_set = ElDec + ElInc + std::min(std::fabs(slope0Inc_to_be_set),std::fabs(slope0Dec_to_be_set))*partitionSize[i];
              // // if (lowerbound_to_be_set < El){
              // //   lowerbound_to_be_set = Eul;
              // //   slope[i][j][0] = 0.;
              // // }

              // if(upperbound_to_be_set - lowerbound_to_be_set <= MC__ISM_COMPUTATION_TOL){
              //   mat[i][j] = T(std::min(lowerbound_to_be_set,upperbound_to_be_set),std::max(lowerbound_to_be_set,upperbound_to_be_set))  + fopt_over_ndep;                
              //   slope[i][j][0] = 0.;
              //   slope[i][j][1] = 0.;
              // }
              // else
              //   mat[i][j] = T(lowerbound_to_be_set , upperbound_to_be_set ) + fopt_over_ndep;                
            }
            break;
          }          
          //throw typename ISModel<T>::Exceptions( ISModel<T>::Exceptions::UNDEF ); 
          
        }
      }
    }
  }
}




template <typename T>
inline
void
ISModel<T>::_asym_slope_relu
( std::vector<std::vector<T>>& mat, std::vector<std::vector<std::vector<long double>>>& slope, const std::vector<long double>& partitionSize,
    unsigned const& ndep)
const
{

  assert( !mat.empty() );

  // // Get the maxima of all components of the input underestimator, stored in _L2[i]
  // long double sigma_o = 0.;       // <- the maximum of the input underestimator
  // for( unsigned int i=0; i<_nvar; i++ ){
  //   if( mat[i].empty() ) continue;
  //   long double _tmp_row = ((long double)(Op<T>::l( mat[i][0]) + std::fabs(slope[i][0][0]*partitionSize[i])));
  //    for( unsigned int j=1; j<_ndiv; j++ ){
  //      _tmp_row = std::max(_tmp_row, ((long double)(Op<T>::l( mat[i][j]) + std::fabs(slope[i][j][0]*partitionSize[i]))));
  //    }
  //    _L2[i]= _tmp_row;
  //    sigma_o += _L2[i];
  // }

  // // Get the minima of all components of the input overestimator, stored in _U2[i] 
  // long double sigma_u = 0.;       // <- the minimum of the input overestimator
  // for( unsigned int i=0; i<_nvar; i++ ){
  //   if( mat[i].empty() ) continue;
  //   long double _tmp_row = ((long double)(Op<T>::u( mat[i][0]) - std::fabs(slope[i][0][1]*partitionSize[i])));
  //    for( unsigned int j=1; j<_ndiv; j++ ){
  //      _tmp_row = std::min(_tmp_row, ((long double)(Op<T>::u( mat[i][j]) - std::fabs(slope[i][j][1]*partitionSize[i]) )));
  //    }
  //    _U2[i]= _tmp_row;
  //    sigma_u += _U2[i];
  // }

  // bool updateUnderEstimator = (sigma_o > MC__ISM_COMPUTATION_TOL);   // <- only compute the underestimator when some elements are greater than 0, o.w. we set all elements 0.
  // bool updateOverEstimator = (sigma_u <  -MC__ISM_COMPUTATION_TOL);   // <- only update the overestimator when some elements are less than 0
  
  // if(false){//(!updateOverEstimator) && (!updateUnderEstimator)){       // only need to set underestimator to be zero.
  //   for( unsigned int i=0; i<_nvar; i++ ){
  //     if( mat[i].empty() ) continue;   
  //     for( unsigned int j=0; j<_ndiv; j++ ){
  //       mat[i][j] = T(0.,1.)*Op<T>::u(mat[i][j]);//T(0.,Op<T>::u(mat[i][j]));   // To avoid NaN in Filib
  //       slope[i][j][0] = 0.;
  //     }
  //   }
  // }
  // else{

  T bnd = _B( mat, 1 );
  auto const& f = [=]( const double& x ){ return std::max( x, 0. ); };
  auto const& fDerv = [=]( const double& x ){ return x > 0? 1.: 0.; };  

  // anchor points
  long double lambda ( Op<T>::l(bnd) );
  long double mu ( Op<T>::u(bnd) );
  long double sum_r1( mu - lambda );
  long double C1( 0. ), C2(0.);
  for( unsigned int i=0; i<_nvar; i++ ){
    if( mat[i].empty() ) continue;
    _c1[i] = _L1[i]; C1 += _c1[i];
    _c2[i] = _U1[i]; C2 += _c2[i];     
    // switch( imid ){
    //   case ICONV: _c1[i] = _L1[i], C1 += _c1[i]; break;
    //   case ICONC: _c1[i] = _U1[i], C1 += _c1[i]; break;
    //   case ICUT:  _c1[i] = _L1[i]; C1 += _c1[i];
    //               _c2[i] = _U1[i]; C2 += _c2[i]; break;
    // }
    _r1[i] = ( _U1[i] - _L1[i] );
  }
  if(std::fabs(C1 - lambda) > 1e2*MC__ISM_COMPUTATION_TOL || std::fabs(C2 - mu) > 1e2*MC__ISM_COMPUTATION_TOL){
    std::cout << "numerical error in asym_slope_relu" << std::endl;
    std::cout << std::setprecision(18) << std::fabs(C1 - lambda) - MC__ISM_COMPUTATION_TOL << std::endl;
    std::cout << std::setprecision(18) << std::fabs(C2 - mu) - MC__ISM_COMPUTATION_TOL << std::endl;
    throw typename ISModel<T>::Exceptions( ISModel<T>::Exceptions::UNDEF ); 
  }
      

  for( unsigned int i=0; i<_nvar; i++ ){
    if( mat[i].empty() ) continue;   
    if( isequal( _r1[i], 0. ) ){
      for( unsigned int j=0; j<_ndiv; j++ ){
        mat[i][j] = T(0.);
        slope[i][j][0] = 0.;
        slope[i][j][1] = 0.;
      }
      continue;
    }
    else{
      const long double rowOffsetUnder = - _c1[i] + C1;
      const long double theta_i = _r1[i] / sum_r1;//(_c2[i] - _U2[i])/(C2 - sigma_u);
      const long double theta_i_times_mu = theta_i*C2;
      const long double rowOffsetOver =  - _c2[i] + theta_i_times_mu;

      for( unsigned int j=0; j<_ndiv; j++ ){
        const long double zU = Op<T>::u(mat[i][j]);
        const long double upPtOverEstmt = zU + rowOffsetOver;
        const long double Du = f(upPtOverEstmt);
        
        const long double delta_u = std::fabs(slope[i][j][1]*partitionSize[i]);        
        long double slope1_to_be_multiplied = fDerv(upPtOverEstmt - delta_u);
        if (delta_u > 1e2 * MC__ISM_COMPUTATION_TOL)
          slope1_to_be_multiplied = (Du - f(upPtOverEstmt - delta_u))/delta_u; 
        slope[i][j][1] = slope1_to_be_multiplied*slope[i][j][1];

        const long double zL = Op<T>::l(mat[i][j]);            
        const long double loPtUnderEstmt = zL + rowOffsetUnder;          
        const long double El = f( loPtUnderEstmt );
        const long double slope0_to_be_multiplied = fDerv( loPtUnderEstmt );
        slope[i][j][0] = slope0_to_be_multiplied*slope[i][j][0];


        if(Du - El <= MC__ISM_COMPUTATION_TOL){
          mat[i][j] = T(std::min(Du,El),std::max(Du,El)); 
          slope[i][j][0] = 0.;
          slope[i][j][1] = 0.;
          // slope[i][j][0] = std::min(slope[i][j][0],slope[i][j][1]);
          // if (slope[i][j][0]*partitionSize[i] >1e2*MC__ISM_COMPUTATION_TOL){
          //   slope[i][j][0] = 0.;
          //   slope[i][j][1] = 0.;
          // }
          // else{
          //   slope[i][j][1] = slope[i][j][0];
          // }                          
        }      
        else
          mat[i][j] = T( El, Du );
      }
    }
  }



//         long double theta_i = _r1[i] / sum_r1;//(_c2[i] - _U2[i])/(C2-sigma_u);//_r1[i] / sum_r1; 
//         long double theta_i_times_mu = theta_i*C2;
//         long double rowOffset = - _c2[i] + theta_i_times_mu;
// //        if(_U2[i] + rowOffset < -MC__ISM_COMPUTATION_TOL){          // Only update the ith component of the overestimator when there are some elements less than zero
//           for( unsigned int j=0; j<_ndiv; j++ ){
//             long double zU = Op<T>::u(mat[i][j]);
//             long double upPtOverEstmt = zU + rowOffset;
//             long double Du = f(upPtOverEstmt);

//             if(Du <= MC__ISM_COMPUTATION_TOL){
//               mat[i][j] = T(std::min((double)Du,0.),std::max((double)Du,0.)); 
//               slope[i][j][1] = 0.;           
//               continue;  
//             }      
            
//             mat[i][j] = T(Op<T>::l(mat[i][j]),Du);

//             long double delta_u = std::fabs(slope[i][j][1]*partitionSize[i]);
//             long double loPtOverEstmt = upPtOverEstmt - delta_u;
//             long double slope1_to_be_multiplied = fDerv(loPtOverEstmt);
//             if (delta_u > 1e2 * MC__ISM_COMPUTATION_TOL)
//               slope1_to_be_multiplied = (Du - f(loPtOverEstmt))/delta_u;
//             slope[i][j][1] = slope1_to_be_multiplied*slope[i][j][1];
// //          }
// //        }       

//         rowOffset = - _c1[i] + C1; // reuse the containter
// //        if(_L2[i] + rowOffset > MC__ISM_COMPUTATION_TOL){   // Only update the ith component of the underestimator when there are some elements greater than zero
// //          for( unsigned int j=0; j<_ndiv; j++ ){
//             long double zL = Op<T>::l(mat[i][j]);
//             long double loPtUnderEstmt = zL + rowOffset;          
//             //long double upPtUnderEstmt = loPtUnderEstmt + ;
//             long double El = f( loPtUnderEstmt );

//             if(Op<T>::u(mat[i][j]) - El <= MC__ISM_COMPUTATION_TOL){
//               mat[i][j] = T(std::min(Op<T>::u(mat[i][j]),(double)El),std::max(Op<T>::u(mat[i][j]),(double)El)); 
//               slope[i][j][0] = 0.;     
//               continue;       
//             }      
            
//             mat[i][j] = T( El, Op<T>::u(mat[i][j]) );

//             long double slope0_to_be_multiplied = fDerv( loPtUnderEstmt );
//             slope[i][j][0] = slope0_to_be_multiplied*slope[i][j][0];

//           }
//         // }
//         // else{
//         //   for( unsigned int j=0; j<_ndiv; j++ ){
//         //     mat[i][j] = T(std::min(Op<T>::u(mat[i][j]),0.),std::max(Op<T>::u(mat[i][j]),0.));
//         //     slope[i][j][0] = 0.;
//         //   }          
//         // }


}



template <typename T>
inline
void
ISModel<T>::_asym_slope_relu_ws
( std::vector<std::vector<T>>& mat, std::vector<std::vector<std::vector<long double>>>& slope, const std::vector<long double>& partitionSize,
    unsigned const& ndep,
    std::vector<std::vector<std::vector<long double>>>& shadow, std::vector<std::vector<std::vector<long double>>>& shadow_slope)
const
{
#ifndef NOTTOTRACKSHADOW
  std::cout << "relu_ws" << std::endl;
  assert( !mat.empty() );
#endif

  // Get the maxima of all components of the input underestimator, stored in _L2[i]
  long double sigma_o = 0.;       // <- the maximum of the input underestimator
  for( unsigned int i=0; i<_nvar; i++ ){
    if( mat[i].empty() ) continue;
    long double _tmp_row = ((long double)(Op<T>::l( mat[i][0]) + std::fabs(slope[i][0][0]*partitionSize[i])));
     for( unsigned int j=1; j<_ndiv; j++ ){
       _tmp_row = std::max(_tmp_row, ((long double)(Op<T>::l( mat[i][j]) + std::fabs(slope[i][j][0]*partitionSize[i]))));
     }
     _L2[i]= _tmp_row;
     sigma_o += _L2[i];
  }

  // // Get the minima of all components of the input overestimator, stored in _U2[i] 
  // long double sigma_u = 0.;       // <- the minimum of the input overestimator
  // for( unsigned int i=0; i<_nvar; i++ ){
  //   if( mat[i].empty() ) continue;
  //   long double _tmp_row = ((long double)(Op<T>::u( mat[i][0]) - std::fabs(slope[i][0][1]*partitionSize[i])));
  //    for( unsigned int j=1; j<_ndiv; j++ ){
  //      _tmp_row = std::min(_tmp_row, ((long double)(Op<T>::u( mat[i][j]) - std::fabs(slope[i][j][1]*partitionSize[i]) )));
  //    }
  //    _U2[i]= _tmp_row;
  //    sigma_u += _U2[i];
  // }

  // bool updateUnderEstimator = (sigma_o > MC__ISM_COMPUTATION_TOL);   // <- only compute the underestimator when some elements are greater than 0, o.w. we set all elements 0.
  // bool updateOverEstimator = (sigma_u <  -MC__ISM_COMPUTATION_TOL);   // <- only update the overestimator when some elements are less than 0
   
  // if(false){//(!updateOverEstimator) && (!updateUnderEstimator)){       // only need to set underestimator to be zero.
  //   for( unsigned int i=0; i<_nvar; i++ ){
  //     if( mat[i].empty() ) continue;   
  //     for( unsigned int j=0; j<_ndiv; j++ ){
  //       mat[i][j] = T(0.,1.)*Op<T>::u(mat[i][j]);//T(0.,Op<T>::u(mat[i][j]));   // To avoid NaN in Filib
  //       slope[i][j][0] = 0.;
  //     }
  //   }
  // }
  // else{

  T bnd = _B( mat, 1 );
  auto const& f = [=]( const double& x ){ return std::max( x, 0. ); };
  //auto const& fDerv = [=]( const double& x ){ return x > 1e2*MC__ISM_COMPUTATION_TOL? 1.: 0.; };  
  auto const& fDerv = [=]( const double& x ){ return x > 0.? 1.: 0.; };  
  // anchor points
  long double lambda ( Op<T>::l(bnd) );
  long double mu ( Op<T>::u(bnd) );
  long double sum_r1( mu - lambda );
  long double C1( 0. ), C2(0.);
  for( unsigned int i=0; i<_nvar; i++ ){
    if( mat[i].empty() ) continue;
    _c1[i] = _L1[i]; C1 += _c1[i];
    _c2[i] = _U1[i]; C2 += _c2[i];     
    // switch( imid ){
    //   case ICONV: _c1[i] = _L1[i], C1 += _c1[i]; break;
    //   case ICONC: _c1[i] = _U1[i], C1 += _c1[i]; break;
    //   case ICUT:  _c1[i] = _L1[i]; C1 += _c1[i];
    //               _c2[i] = _U1[i]; C2 += _c2[i]; break;
    // }
    _r1[i] = ( _U1[i] - _L1[i] );
  }
  if(std::fabs(C1 - lambda) > 1e2*MC__ISM_COMPUTATION_TOL || std::fabs(C2 - mu) > 1e2*MC__ISM_COMPUTATION_TOL){
    std::cout << "numerical error in asym_slope_relu_ws" << std::endl;
    std::cout << std::setprecision(18) << std::fabs(C1 - lambda) - MC__ISM_COMPUTATION_TOL << std::endl;
    std::cout << std::setprecision(18) << std::fabs(C2 - mu) - MC__ISM_COMPUTATION_TOL << std::endl;
    throw typename ISModel<T>::Exceptions( ISModel<T>::Exceptions::UNDEF ); 
  }
      
  const long double shadow_global_offset = (1.0 - 1.0/((double) ndep))*sigma_o;
#ifndef NOTTOTRACKSHADOW
  std::cout << "    shadow_global_offset = " << shadow_global_offset << std::endl;
#endif
  for( unsigned int i=0; i<_nvar; i++ ){
    if( mat[i].empty() ) continue;   
    if( isequal( _r1[i], 0. ) ){
      shadow[0][i].clear();
      shadow[0][i].resize(_ndiv);
      shadow_slope[0][i].clear();
      shadow_slope[0][i].resize(_ndiv);

      for( unsigned int j=0; j<_ndiv; j++ ){
        mat[i][j] = T(0.);
        slope[i][j][0] = 0.;
        slope[i][j][1] = 0.;
        shadow[0][i][j] = 0.;
        shadow_slope[0][i][j] = 0.; 
      }
      continue;
    }
    else{
      const long double rowOffsetUnder = - _c1[i] + C1;
      const long double theta_i = _r1[i] / sum_r1;//(_c2[i] - _U2[i])/(C2 - sigma_u);
      const long double theta_i_times_mu = theta_i*C2;
      const long double rowOffsetOver =  - _c2[i] + theta_i_times_mu;
      const long double rowOffsetShadow = - _L2[i] + sigma_o;
      // shadow[0][i].clear();
      // shadow[0][i].resize(_ndiv);      
      // shadow_slope[0][i].clear();
      // shadow_slope[0][i].resize(_ndiv);
      for( unsigned int j=0; j<_ndiv; j++ ){
        const long double zU = Op<T>::u(mat[i][j]);
        const long double upPtOverEstmt = zU + rowOffsetOver;
        const long double Du = f(upPtOverEstmt);
        
        const long double delta_u = std::fabs(slope[i][j][1]*partitionSize[i]);        
        long double slope1_to_be_multiplied = fDerv(upPtOverEstmt - delta_u);
        if (delta_u > 1e2 * MC__ISM_COMPUTATION_TOL)
          slope1_to_be_multiplied = (Du - f(upPtOverEstmt - delta_u))/delta_u; 
        slope[i][j][1] = slope1_to_be_multiplied*slope[i][j][1];

        const long double zL = Op<T>::l(mat[i][j]);            
        const long double loPtUnderEstmt = zL + rowOffsetUnder;          
        const long double El = f( loPtUnderEstmt );
        const long double slope0_to_be_set = fDerv( loPtUnderEstmt )*slope[i][j][0];
        
        //if(shadow_global_offset > 0){

          // compute shadow
          // Note that the function max(zi - zi^U + mu) - max(zi - zi^L + lambda) is no longer convex but still nondecreasing
          const long double delta_l = std::fabs(slope[i][j][0]*partitionSize[i]); 
          const long double loPtShadow = zL + rowOffsetShadow;//zL - _c2[i] + C2;
          const long double upPtShadow = loPtShadow + delta_l;//zL - _c2[i] + C2;        
          const long double shadowLoBnd = f( loPtShadow ) - El;
          //const long double shadowUpBnd = f( upPtShadow ) - f( loPtUnderEstmt + delta_l);
          shadow[0][i][j] = shadowLoBnd - f(shadow_global_offset);
          //shadow[2][i][j] = shadowUpBnd - f(shadow_global_offset);
            
          //if(shadowUpBnd > shadowLoBnd + 1e2*MC__ISM_COMPUTATION_TOL ){
            // The update of the derv needs special consideration
            const long double upBndUnderEstmt = f(loPtUnderEstmt + delta_l);
            long double slope0shadow_to_be_multiplied =  0.;//std::min(shadowFuncDervAtRightPt, shadowFuncDervAtLeftPt);
            if (shadowLoBnd == 0.)//f(loPtUnderEstmtShadow) < 1e2*MC__ISM_COMPUTATION_TOL)
              slope0shadow_to_be_multiplied = 0.;
            else if (f(loPtUnderEstmt+delta_l + 1e5*MC__ISM_COMPUTATION_TOL) == 0. ){
              slope0shadow_to_be_multiplied = std::min(fDerv(loPtShadow) - fDerv(loPtUnderEstmt),fDerv(upPtShadow) - fDerv(loPtUnderEstmt + delta_l));
              //std::cout << f(loPtUnderEstmt+delta_l + 1e2*MC__ISM_COMPUTATION_TOL) << "," << shadowFuncValueAtLeftPt << std::endl;
            }
            else if (f(loPtUnderEstmt-1e2*MC__ISM_COMPUTATION_TOL) == 0.){
              if ( delta_l > 1e2 * MC__ISM_COMPUTATION_TOL)
                slope0shadow_to_be_multiplied = (f(upPtShadow) - upBndUnderEstmt - (f(loPtShadow) -f(loPtUnderEstmt)))/(delta_l);
              else 
                ;//slope0shadow_to_be_multiplied = fDerv(upPtShadow) - fDerv(loPtUnderEstmt + delta_l);
            }
            shadow_slope[0][i][j] = slope0shadow_to_be_multiplied * slope[i][j][0];
    
            // long double slope1shadow_to_be_multiplied =  1.;
            // if (upBndUnderEstmt > 0.)
            //   slope1shadow_to_be_multiplied = 0.;
            // else if (f( loPtShadow ) > 1e2*MC__ISM_COMPUTATION_TOL ){
            //   slope1shadow_to_be_multiplied = 1.;
            //   //std::cout << f(loPtUnderEstmt+delta_l + 1e2*MC__ISM_COMPUTATION_TOL) << "," << shadowFuncValueAtLeftPt << std::endl;
            // }
            // else if (f(upPtShadow-1e2*MC__ISM_COMPUTATION_TOL) == 0.){
            //   if ( delta_l > 1e2 * MC__ISM_COMPUTATION_TOL)
            //     slope1shadow_to_be_multiplied = (f(upPtShadow) - f(loPtUnderEstmt + delta_l) - (f(loPtShadow) -f(loPtUnderEstmt)))/(delta_l);
            //   else 
            //     slope1shadow_to_be_multiplied = fDerv(loPtShadow) - fDerv(loPtUnderEstmt);
            // }
            // shadow_slope[2][i][j] = slope1shadow_to_be_multiplied * slope[i][j][0];
    
            slope[i][j][0] = slope0_to_be_set;
            
            //const long double shadowFuncDervAtRightPt = fDerv( zU - _L2[i] + sigma_o ) - fDerv(zU + rowOffsetUnder);
            //const long double shadowFuncDervAtRightPt = fDerv( zL + delta_l - _L2[i] + sigma_o ) - fDerv(zL + delta_l + rowOffsetUnder);
            //const long double shadowFuncDervAtLeftPt  = std::min(fDerv(zL - _L2[i] + sigma_o) -f(zL + rowOffsetUnder)); 
            //const long double shadowFuncAtLeftPt = f( zL - _L2[i] + sigma_o - 1e2*MC__ISM_COMPUTATION_TOL) - f(zL + rowOffsetUnder - 1e2*MC__ISM_COMPUTATION_TOL);
            //if(shadowFuncDervAtRightPt == 1. && shadowFuncDervAtLeftPt == 1.)
            //  slope0shadow_to_be_multiplied = 1.;
            //if (shadowFuncAtLeftPt > 1e2*MC__ISM_COMPUTATION_TOL )
            //  slope0shadow_to_be_multiplied = 1.;
            //if(shadowFuncDervAtRightPt == 0. && shadowFuncAtLeftPt > 1e2*MC__ISM_COMPUTATION_TOL){
            //  if( delta_l > 1e2 * MC__ISM_COMPUTATION_TOL){
            //    slope0shadow_to_be_multiplied = (f(zL + delta_l - _L2[i] + sigma_o) - f(zL + delta_l + rowOffsetUnder) - (f(zL - _L2[i] + sigma_o) -f(zL + rowOffsetUnder)))/(delta_l);
            //  }
            //  else
            //    slope0shadow_to_be_multiplied = 0.;
            //}
            //shadow_slope[0][i][j] = slope0shadow_to_be_multiplied * original_slope0;  
          //}
          //else if(std::fabs(std::min(shadow[2][i][j],shadow[0][i][j])) < 1e2*MC__ISM_COMPUTATION_TOL){
          //  shadow[0][i][j] = 0;
          //  shadow[2][i][j] = 0;
          //  shadow_slope[0][i][j] = 0;
          //  shadow_slope[2][i][j] = 0;
          //}
          //else{
          //}
  
          // if(std::fabs(shadow_slope[2][i][j]) * partitionSize[i] > shadow[2][i][j]-shadow[0][i][j]){
          //   if (partitionSize[i] > 1e5 * MC__ISM_COMPUTATION_TOL)
          //     shadow_slope[2][i][j] = (shadow[2][i][j]-shadow[0][i][j])/partitionSize[i];
          //   else{
          //     shadow[2][i][j] = shadow[2][i][j]+0.5*(std::fabs(shadow_slope[2][i][j]) * partitionSize[i]); 
          //     shadow[0][i][j] = shadow[0][i][j]-0.5*(std::fabs(shadow_slope[2][i][j]) * partitionSize[i]); 
          //   }   
          // }
            
          // if(std::fabs(shadow_slope[0][i][j]) * partitionSize[i] > shadow[2][i][j]-shadow[0][i][j]){
          //   if (partitionSize[i] > 1e5 * MC__ISM_COMPUTATION_TOL)
          //     shadow_slope[0][i][j] = (shadow[2][i][j]-shadow[0][i][j])/partitionSize[i];
          //   else{
          //     shadow[2][i][j] = shadow[2][i][j]+0.5*(std::fabs(shadow_slope[0][i][j]) * partitionSize[i]); 
          //     shadow[0][i][j] = shadow[0][i][j]-0.5*(std::fabs(shadow_slope[0][i][j]) * partitionSize[i]); 
          //   }         
          // } 
        //}

 
        if(Du - El <= -1e2*MC__ISM_COMPUTATION_TOL){
          std::cout << "    error in relu_xs " << Du<<"," <<El<<std::endl; 
          throw typename ISModel<T>::Exceptions( ISModel<T>::Exceptions::UNDEF );   
        }
        if(Du - El <= MC__ISM_COMPUTATION_TOL){  
          mat[i][j] = T(std::min(Du,El),std::max(Du,El));
          slope[i][j][0] = 0.;
          slope[i][j][1] = 0.;   
//          shadow_slope[0][i][j] = 0.;       
          //slope[i][j][0] = std::min(slope[i][j][0],slope[i][j][1]);
          //if (slope[i][j][0]*partitionSize[i] >1e2*MC__ISM_COMPUTATION_TOL){
          //  slope[i][j][0] = 0.;
          //  slope[i][j][1] = 0.;
          //}
          //else{
          //  slope[i][j][1] = slope[i][j][0];
          //}
               
          //if(std::max(Du,El)<MC__ISM_COMPUTATION_TOL){
          //slope[i][j][0] = 0.;
          //slope[i][j][1] = 0.;            
          //} 
          //slope[i][j][0] = 0.;
          //slope[i][j][1] = 0.;
          //shadow[0][i][j] = 0.;   
          //shadow_slope[0][i][j] = 0.;
        }      
        else
          mat[i][j] = T( El, Du );
      }
    }
  }



//         long double theta_i = _r1[i] / sum_r1;//(_c2[i] - _U2[i])/(C2-sigma_u);//_r1[i] / sum_r1; 
//         long double theta_i_times_mu = theta_i*C2;
//         long double rowOffset = - _c2[i] + theta_i_times_mu;
// //        if(_U2[i] + rowOffset < -MC__ISM_COMPUTATION_TOL){          // Only update the ith component of the overestimator when there are some elements less than zero
//           for( unsigned int j=0; j<_ndiv; j++ ){
//             long double zU = Op<T>::u(mat[i][j]);
//             long double upPtOverEstmt = zU + rowOffset;
//             long double Du = f(upPtOverEstmt);

//             if(Du <= MC__ISM_COMPUTATION_TOL){
//               mat[i][j] = T(std::min((double)Du,0.),std::max((double)Du,0.)); 
//               slope[i][j][1] = 0.;           
//               continue;  
//             }      
            
//             mat[i][j] = T(Op<T>::l(mat[i][j]),Du);

//             long double delta_u = std::fabs(slope[i][j][1]*partitionSize[i]);
//             long double loPtOverEstmt = upPtOverEstmt - delta_u;
//             long double slope1_to_be_multiplied = fDerv(loPtOverEstmt);
//             if (delta_u > 1e2 * MC__ISM_COMPUTATION_TOL)
//               slope1_to_be_multiplied = (Du - f(loPtOverEstmt))/delta_u;
//             slope[i][j][1] = slope1_to_be_multiplied*slope[i][j][1];
// //          }
// //        }       

//         rowOffset = - _c1[i] + C1; // reuse the containter
// //        if(_L2[i] + rowOffset > MC__ISM_COMPUTATION_TOL){   // Only update the ith component of the underestimator when there are some elements greater than zero
// //          for( unsigned int j=0; j<_ndiv; j++ ){
//             long double zL = Op<T>::l(mat[i][j]);
//             long double loPtUnderEstmt = zL + rowOffset;          
//             //long double upPtUnderEstmt = loPtUnderEstmt + ;
//             long double El = f( loPtUnderEstmt );

//             if(Op<T>::u(mat[i][j]) - El <= MC__ISM_COMPUTATION_TOL){
//               mat[i][j] = T(std::min(Op<T>::u(mat[i][j]),(double)El),std::max(Op<T>::u(mat[i][j]),(double)El)); 
//               slope[i][j][0] = 0.;     
//               continue;       
//             }      
            
//             mat[i][j] = T( El, Op<T>::u(mat[i][j]) );

//             long double slope0_to_be_multiplied = fDerv( loPtUnderEstmt );
//             slope[i][j][0] = slope0_to_be_multiplied*slope[i][j][0];

//           }
//         // }
//         // else{
//         //   for( unsigned int j=0; j<_ndiv; j++ ){
//         //     mat[i][j] = T(std::min(Op<T>::u(mat[i][j]),0.),std::max(Op<T>::u(mat[i][j]),0.));
//         //     slope[i][j][0] = 0.;
//         //   }          
//         // }


}

template <typename T>
inline
void
ISModel<T>::_add_aggregate_shadow
( std::vector<std::vector<T>>& Amat, const std::vector<std::vector<T>>& Bmat, std::vector<std::vector<std::vector<long double>>>& Aslope, const std::vector<std::vector<std::vector<long double>>>& Bslope,
  std::vector<std::vector<std::vector<long double>>>& Ashadow, const std::vector<std::vector<std::vector<long double>>>& Bshadow, 
  std::vector<std::vector<std::vector<long double>>>& Ashadow_slope, const std::vector<std::vector<std::vector<long double>>>& Bshadow_slope, 
  std::vector<double>& Ashadow_info, const std::vector<double>& Bshadow_info, const std::vector<long double>& partitionSize, unsigned const& ndep, T const& bndB)
const
{
#ifndef NOTTOTRACKSHADOW  
  std::cout << "start agrt" << std::endl;
#endif
  const bool A_under_shadow = Ashadow_info[1];
  const bool B_under_shadow = Bshadow_info[1];
  const bool A_over_shadow = Ashadow_info[2];
  const bool B_over_shadow = Bshadow_info[2];

   // process the aggreagation mode
  const unsigned char agrgatMode = A_under_shadow + 2*B_under_shadow + 4*A_over_shadow + 8*B_over_shadow;
  if (agrgatMode == 0 ){     
    std::cout << "    error in aggregating in preprocessing addition" << std::endl; 
    throw typename ISModel<T>::Exceptions( ISModel<T>::Exceptions::UNDEF );   
  }
   
   bool toUpdateUnderEstmt = false;
   bool toUpdateOverEstmt = false;
   
  switch (agrgatMode){
    case 0x01: 
    case 0x02: 
    case 0x03: toUpdateUnderEstmt = true; break;

    case 0x04: 
    case 0x08: 
    case 0x0C: toUpdateOverEstmt = true; break;

    case 0x05: //directionOver = false; directionUnder = false; toUpdateUnderEstmt = true; toUpdateOverEstmt = true; break;
    case 0x0A: //directionOver = false; directionUnder = false; toUpdateUnderEstmt = true; toUpdateOverEstmt = true; break;
    case 0x06: //; break;    
    case 0x09: toUpdateUnderEstmt = true; toUpdateOverEstmt = true; break;

    case 0x07: 
    case 0x0B: toUpdateUnderEstmt = true; toUpdateOverEstmt = true; break;

    case 0x0D: 
    case 0x0E: 
    case 0x0F: toUpdateUnderEstmt = true; toUpdateOverEstmt = true; break;   
  }
#ifndef NOTTOTRACKSHADOW  
  std::cout << "    agrt mode setted, Mode: " << (int)agrgatMode << std::endl;
  std::cout << "    Au    Bu    Ao    Bo    updO    updU    " << "\n"
            << "     " << A_under_shadow << "     " << B_under_shadow 
            << "     " << A_over_shadow  << "     " << B_over_shadow 
            << "     " << toUpdateOverEstmt << "        " << toUpdateUnderEstmt  << std::endl;
#endif
  // Step 1: cross addition
  T bndA = _B(Amat,1);
  const long double muA = Op<T>::u(bndA);
  const long double muB = Op<T>::u(bndB);  
  const long double lambdaA = Op<T>::l(bndA);
  const long double lambdaB = Op<T>::l(bndB);    
#ifndef NOTTOTRACKSHADOW   
  std::cout << "    A in [ " << lambdaA << " , " <<  muA <<" ]" << std::endl;
  std::cout << "    B in [ " << lambdaB << " , " <<  muB <<" ]" << std::endl;
#endif
  bool underEstmtUpdated = false;
  bool overEstmtUpdated  = false;

  std::vector<std::vector<long double>> NewOverEstmt; 
  std::vector<std::vector<long double>> NewOverEstmt_slope;   
  std::vector<std::vector<long double>> NewOverEstmt_shadow; 
  std::vector<std::vector<long double>> NewOverEstmt_shadow_slope;

  std::vector<std::vector<long double>> NewUnderEstmt; 
  std::vector<std::vector<long double>> NewUnderEstmt_slope;   
  std::vector<std::vector<long double>> NewUnderEstmt_shadow; 
  std::vector<std::vector<long double>> NewUnderEstmt_shadow_slope;  
  //bool shadowPointingToB = false;
  std::vector<std::vector<std::vector<long double>>> LocalBshadow = Bshadow;
  std::vector<std::vector<std::vector<long double>>> LocalBshadow_slope = Bshadow_slope;


  double Ashadow_sign = 0;
  double Bshadow_sign = 0;
  if(Ashadow_info[0] > 0 )
    Ashadow_sign = 1;
  else if(Ashadow_info[0] < 0 ) 
    Ashadow_sign = -1;
  else{
    std::cout << "            Ashadow_sign error at shadow addition: Ashadow_info[0] = " << Ashadow_info[0] <<std::endl;
    throw typename ISModel<T>::Exceptions( ISModel<T>::Exceptions::UNDEF ); 
  }
  if(Bshadow_info[0] > 0 )
    Bshadow_sign = 1;
  else if(Bshadow_info[0] < 0 ) 
    Bshadow_sign = -1;
  else{
    std::cout << "            Bshadow_sign error at shadow addition: Bshadow_info[0] = " << Bshadow_info[0] << std::endl;
    throw typename ISModel<T>::Exceptions( ISModel<T>::Exceptions::UNDEF ); 
  }

  NewOverEstmt.resize(_nvar);
  NewOverEstmt_slope.resize(_nvar);
  NewOverEstmt_shadow.resize(_nvar);    
  NewOverEstmt_shadow_slope.resize(_nvar);    


  if (toUpdateOverEstmt){
#ifndef NOTTOTRACKSHADOW  
    std::cout << "        toUpdateOverEstmt" << std::endl;
#endif
    std::vector<std::vector<long double>> NewAOverEstmt; 
    std::vector<std::vector<long double>> NewAOverEstmt_slope;     
    std::vector<std::vector<long double>> NewBOverEstmt; 
    std::vector<std::vector<long double>> NewBOverEstmt_slope;   

    NewAOverEstmt.resize(_nvar);
    NewAOverEstmt_slope.resize(_nvar);
    NewBOverEstmt.resize(_nvar);
    NewBOverEstmt_slope.resize(_nvar);


    bool Aupdated = false;
    bool Bupdated = false;

    if(A_over_shadow){
    // B + A_shadow -> New B 
#ifndef NOTTOTRACKSHADOW      
      std::cout << "            A_over_shadow" << std::endl;
#endif      
      long double shadowEnhancedMax = 0.;
      const double shadow_sign = Ashadow_sign;

      if(shadow_sign > 0){
#ifndef NOTTOTRACKSHADOW         
        std::cout << "                shadow_sign > 0" << std::endl;       
#endif          
        for( unsigned int i=0; i<_nvar; i++ ){
          bool updateShadowEnhancedMax = false;
          if( !Amat[i].empty() && !Bmat[i].empty() ){
#ifndef NOTTOTRACKSHADOW                  
            std::cout << "                    !Amat[i].empty() && !Bmat[i].empty(), processing row No. " << i << std::endl;
#endif
            NewBOverEstmt[i].resize(_ndiv); 
            NewBOverEstmt_slope[i].resize(_ndiv);
            for( unsigned long long j=0; j<_ndiv; j++ ){
              NewBOverEstmt[i][j] = Op<T>::u(Bmat[i][j]) + Ashadow[1][i][j];
              if((Ashadow_slope[1][i][j] < 0.) != (Bslope[i][j][1] < 0.)){   
                long double tightener = std::min(std::fabs(Ashadow_slope[1][i][j]),std::fabs(Bslope[i][j][1]))*(partitionSize[i]);
                NewBOverEstmt[i][j] -= tightener;
              }
              NewBOverEstmt_slope[i][j] = Ashadow_slope[1][i][j] + Bslope[i][j][1];
            }
            updateShadowEnhancedMax = true;
          }
          else if( !Amat[i].empty() ){
#ifndef NOTTOTRACKSHADOW               
            std::cout << "                    !Amat[i].empty() && Bmat[i].empty(), processing row No. " << i << std::endl;
#endif            
            NewBOverEstmt[i].resize(_ndiv); 
            NewBOverEstmt_slope[i].resize(_ndiv);            
            for( unsigned long long j=0; j<_ndiv; j++ ){
              NewBOverEstmt[i][j] = Ashadow[1][i][j];
              NewBOverEstmt_slope[i][j] =  Ashadow_slope[1][i][j];
            }
            updateShadowEnhancedMax = true;
          }
          else if( !Bmat[i].empty() ){
#ifndef NOTTOTRACKSHADOW                
            std::cout << "                    Amat[i].empty() && !Bmat[i].empty(), processing row No. " << i << std::endl;
#endif                
            NewBOverEstmt[i].resize(_ndiv); 
            NewBOverEstmt_slope[i].resize(_ndiv);            
            for( unsigned long long j=0; j<_ndiv; j++ ){
              NewBOverEstmt[i][j] = Op<T>::u(Bmat[i][j]);
              NewBOverEstmt_slope[i][j] =  Bslope[i][j][1];
            }
            updateShadowEnhancedMax = true;
          }          
          if(updateShadowEnhancedMax){
            const auto tmp_max = std::max_element(NewBOverEstmt[i].begin(),NewBOverEstmt[i].end());     
            shadowEnhancedMax += (*tmp_max);
          }
        }
      }
      else if(shadow_sign < 0){
#ifndef NOTTOTRACKSHADOW           
        std::cout << "                shadow_sign < 0" << std::endl;  
#endif 
        for( unsigned int i=0; i<_nvar; i++ ){
          bool updateShadowEnhancedMax = false;
          if( !Amat[i].empty() && !Bmat[i].empty() ){
#ifndef NOTTOTRACKSHADOW   
            std::cout << "                    !Amat[i].empty() && !Bmat[i].empty(), processing row No. " << i << std::endl;
#endif             
            NewBOverEstmt[i].resize(_ndiv); 
            NewBOverEstmt_slope[i].resize(_ndiv);
            for( unsigned long long j=0; j<_ndiv; j++ ){
              NewBOverEstmt[i][j] = Op<T>::u(Bmat[i][j]) - Ashadow[1][i][j];
              if((Ashadow_slope[1][i][j] < 0.) == (Bslope[i][j][1] < 0.)){   
                long double tightener = std::min(std::fabs(Ashadow_slope[1][i][j]),std::fabs(Bslope[i][j][1]))*(partitionSize[i]);
                NewBOverEstmt[i][j] -= tightener;
              }
              NewBOverEstmt_slope[i][j] = - Ashadow_slope[1][i][j] + Bslope[i][j][1];
            }
            updateShadowEnhancedMax = true;
          }
          else if( !Amat[i].empty() ){
#ifndef NOTTOTRACKSHADOW                
            std::cout << "                    !Amat[i].empty() && Bmat[i].empty(), processing row No. " << i << std::endl;
#endif             
            NewBOverEstmt[i].resize(_ndiv); 
            NewBOverEstmt_slope[i].resize(_ndiv);            
            for( unsigned long long j=0; j<_ndiv; j++ ){
              NewBOverEstmt[i][j] = -Ashadow[1][i][j];
              NewBOverEstmt_slope[i][j] =  -Ashadow_slope[1][i][j];
            }
            updateShadowEnhancedMax = true;
          }
          else if( !Bmat[i].empty() ){
#ifndef NOTTOTRACKSHADOW                
            std::cout << "                    Amat[i].empty() && !Bmat[i].empty(), processing row No. " << i << std::endl;
#endif                
            NewBOverEstmt[i].resize(_ndiv); 
            NewBOverEstmt_slope[i].resize(_ndiv);            
            for( unsigned long long j=0; j<_ndiv; j++ ){
              NewBOverEstmt[i][j] = Op<T>::u(Bmat[i][j]);
              NewBOverEstmt_slope[i][j] =  Bslope[i][j][1];
            }
            updateShadowEnhancedMax = true;
          }    

          if(updateShadowEnhancedMax){
            const auto tmp_max = std::max_element(NewBOverEstmt[i].begin(),NewBOverEstmt[i].end());     
            shadowEnhancedMax += (*tmp_max);
          }
        }        
      }
      else {
        std::cout << "error in binary addition aggr overestimator Ashadow when processing shadow sign for overestimator" <<std::endl;
        throw typename ISModel<T>::Exceptions( ISModel<T>::Exceptions::UNDEF ); 
      }
  
      if(shadowEnhancedMax < muB - MC__ISM_COMPUTATION_TOL){
#ifndef NOTTOTRACKSHADOW                
        std::cout << "                Bupdated = true: " << std::setprecision(18) << shadowEnhancedMax << " , " << muB <<std::endl;
#endif
        Bupdated = true;
      }
    }

    if(B_over_shadow){
    // A + B_shadow -> New A 
      long double shadowEnhancedMax = 0.;
      const double shadow_sign = Bshadow_sign;
#ifndef NOTTOTRACKSHADOW        
      std::cout << "            B_under_shadow" << std::endl;
#endif      
      if(shadow_sign > 0){
#ifndef NOTTOTRACKSHADOW
        std::cout << "                shadow_sign > 0" << std::endl;
#endif
        for( unsigned int i=0; i<_nvar; i++ ){
          bool updateShadowEnhancedMax = false;
          if( !Amat[i].empty() && !Bmat[i].empty() ){
#ifndef NOTTOTRACKSHADOW                    
            std::cout << "                    !Amat[i].empty() && !Bmat[i].empty(), processing row No. " << i << std::endl;
#endif            
            NewAOverEstmt[i].resize(_ndiv); 
            NewAOverEstmt_slope[i].resize(_ndiv);
            for( unsigned long long j=0; j<_ndiv; j++ ){
              NewAOverEstmt[i][j] = Op<T>::u(Amat[i][j]) + Bshadow[1][i][j];
              if((Bshadow_slope[1][i][j] < 0.) != (Aslope[i][j][1] < 0.)){   
                long double tightener = std::min(std::fabs(Bshadow_slope[1][i][j]),std::fabs(Aslope[i][j][1]))*(partitionSize[i]);
                NewAOverEstmt[i][j] -= tightener;
              }
              NewAOverEstmt_slope[i][j] = Bshadow_slope[1][i][j] + Aslope[i][j][1];
            }
            updateShadowEnhancedMax = true;
          }
          else if( !Amat[i].empty() ){
#ifndef NOTTOTRACKSHADOW                    
            std::cout << "                    !Amat[i].empty() && Bmat[i].empty(), processing row No. " << i << std::endl;
#endif            
            NewAOverEstmt[i].resize(_ndiv); 
            NewAOverEstmt_slope[i].resize(_ndiv);            
            for( unsigned long long j=0; j<_ndiv; j++ ){
              NewAOverEstmt[i][j] = Op<T>::u(Amat[i][j]);
              NewAOverEstmt_slope[i][j] =  Aslope[i][j][1];
            }
            updateShadowEnhancedMax = true;
          }
          else if( !Bmat[i].empty() ){
#ifndef NOTTOTRACKSHADOW                    
            std::cout << "                    Amat[i].empty() && !Bmat[i].empty(), processing row No. " << i << std::endl;
#endif                
            NewAOverEstmt[i].resize(_ndiv); 
            NewAOverEstmt_slope[i].resize(_ndiv);            
            for( unsigned long long j=0; j<_ndiv; j++ ){
              NewAOverEstmt[i][j] = Bshadow[1][i][j];
              NewAOverEstmt_slope[i][j] =  Bshadow_slope[1][i][j];
            }
            updateShadowEnhancedMax = true;
          }          
        
          if(updateShadowEnhancedMax){
            const auto tmp_max = std::max_element(NewAOverEstmt[i].begin(),NewAOverEstmt[i].end());     
            shadowEnhancedMax += (*tmp_max);
          }
        }
      }
      else if(shadow_sign < 0){
#ifndef NOTTOTRACKSHADOW                  
        std::cout << "                shadow_sign < 0" << std::endl;
#endif           
        for( unsigned int i=0; i<_nvar; i++ ){
          bool updateShadowEnhancedMax = false;
          if( !Amat[i].empty() && !Bmat[i].empty() ){
#ifndef NOTTOTRACKSHADOW                    
            std::cout << "                    !Amat[i].empty() && !Bmat[i].empty(), processing row No. " << i << std::endl;
#endif              
            NewAOverEstmt[i].resize(_ndiv); 
            NewAOverEstmt_slope[i].resize(_ndiv);
            for( unsigned long long j=0; j<_ndiv; j++ ){
              NewAOverEstmt[i][j] = Op<T>::u(Amat[i][j]) - Bshadow[1][i][j];
              if((Bshadow_slope[1][i][j] < 0.) == (Aslope[i][j][1] < 0.)){   
                long double tightener = std::min(std::fabs(Bshadow_slope[1][i][j]),std::fabs(Aslope[i][j][1]))*(partitionSize[i]);
                NewAOverEstmt[i][j] -= tightener;
              }
              NewAOverEstmt_slope[i][j] = - Bshadow_slope[1][i][j] + Aslope[i][j][1];
            }
            updateShadowEnhancedMax = true;
          }
          else if( !Amat[i].empty() ){
#ifndef NOTTOTRACKSHADOW                 
            std::cout << "                    !Amat[i].empty() && Bmat[i].empty(), processing row No. " << i << std::endl;
#endif              
            NewAOverEstmt[i].resize(_ndiv); 
            NewAOverEstmt_slope[i].resize(_ndiv);            
            for( unsigned long long j=0; j<_ndiv; j++ ){
              NewAOverEstmt[i][j] = Op<T>::u(Amat[i][j]);
              NewAOverEstmt_slope[i][j] =  Aslope[i][j][1];
            }
            updateShadowEnhancedMax = true;
          }
          else if( !Bmat[i].empty() ){
#ifndef NOTTOTRACKSHADOW              
            std::cout << "                    Amat[i].empty() && !Bmat[i].empty(), processing row No. " << i << std::endl;
#endif             
            NewAOverEstmt[i].resize(_ndiv); 
            NewAOverEstmt_slope[i].resize(_ndiv);            
            for( unsigned long long j=0; j<_ndiv; j++ ){
              NewAOverEstmt[i][j] = -Bshadow[1][i][j];
              NewAOverEstmt_slope[i][j] =  -Bshadow_slope[1][i][j];
            }
            updateShadowEnhancedMax = true;
          }          
        
          if(updateShadowEnhancedMax){
            const auto tmp_max = std::max_element(NewAOverEstmt[i].begin(),NewAOverEstmt[i].end());     
            shadowEnhancedMax += (*tmp_max);
          }  
        }        
      }
      else {
        std::cout << "error in binary addition aggr overestimator Bshadow when processing shadow sign for overestimator " << shadow_sign << std::endl;
        throw typename ISModel<T>::Exceptions( ISModel<T>::Exceptions::UNDEF ); 
      }
  
      if(shadowEnhancedMax < muA - MC__ISM_COMPUTATION_TOL){
#ifndef NOTTOTRACKSHADOW          
        std::cout << "                Aupdated = true: " << std::setprecision(18) << shadowEnhancedMax << " , " << muA <<std::endl;
#endif               
        Aupdated = true;
      }
        
    }


    // Step 2.1: addition (overestimator)
    long double shadowAMin(0.),shadowBMin(0.),shadowCMin(0.);
#ifndef NOTTOTRACKSHADOW        
    std::cout << "        Addition" << std::endl;
#endif      
    for( unsigned int i=0; i<_nvar; i++ ){
      if( !Amat[i].empty() && !Bmat[i].empty() ){
#ifndef NOTTOTRACKSHADOW             
        std::cout << "                    !Amat[i].empty() && !Bmat[i].empty(), processing row No. " << i << std::endl;
#endif              
        NewOverEstmt[i].resize(_ndiv); 
        NewOverEstmt_slope[i].resize(_ndiv);
        NewOverEstmt_shadow[i].resize(_ndiv);            
        NewOverEstmt_shadow_slope[i].resize(_ndiv);                
        for( unsigned long long j=0; j<_ndiv; j++ ){
          long double const tmp_AoverEstmt_at_ij = Aupdated? NewAOverEstmt[i][j]:Op<T>::u(Amat[i][j]);
          long double const tmp_BoverEstmt_at_ij = Bupdated? NewBOverEstmt[i][j]:Op<T>::u(Bmat[i][j]);
          long double const tmp_AoverEstmtSlp_at_ij = Aupdated? NewAOverEstmt_slope[i][j]:Aslope[i][j][1];
          long double const tmp_BoverEstmtSlp_at_ij = Bupdated? NewBOverEstmt_slope[i][j]:Bslope[i][j][1];
          
          NewOverEstmt[i][j] = tmp_AoverEstmt_at_ij + tmp_BoverEstmt_at_ij;
          if((tmp_AoverEstmtSlp_at_ij < 0.) != (tmp_BoverEstmtSlp_at_ij < 0.)){   
            long double tightener = std::min(std::fabs(tmp_AoverEstmtSlp_at_ij),std::fabs(tmp_BoverEstmtSlp_at_ij))*(partitionSize[i]);
            NewOverEstmt[i][j] -= tightener;
          }
          NewOverEstmt_slope[i][j] = tmp_AoverEstmtSlp_at_ij + tmp_BoverEstmtSlp_at_ij;
          //std::cout << "                        estmator and slope updated at elm No. " << j << std::endl;


          long double tmp_BoverEstmtShad_at_ij = 0.;
          long double tmp_BoverEstmtShadSlp_at_ij =0.;
          if(Bshadow_info[2]!=0){
            if(!Aupdated){
              tmp_BoverEstmtShad_at_ij = Bshadow_sign*Bshadow[1][i][j];
              tmp_BoverEstmtShadSlp_at_ij = Bshadow_sign*Bshadow_slope[1][i][j];             
            }
            else
            {
              tmp_BoverEstmtShadSlp_at_ij = -Bshadow_sign*Bshadow_slope[1][i][j];
              tmp_BoverEstmtShad_at_ij = -Bshadow_sign*Bshadow[1][i][j] + std::fabs(tmp_BoverEstmtShadSlp_at_ij)*partitionSize[i];       
            }
            LocalBshadow[1][i][j] = tmp_BoverEstmtShad_at_ij;
            LocalBshadow_slope[1][i][j] = tmp_BoverEstmtShadSlp_at_ij;                   
          }

          long double tmp_AoverEstmtShad_at_ij = 0.;
          long double tmp_AoverEstmtShadSlp_at_ij = 0.;
          if(Ashadow_info[2]!=0){
            if(!Bupdated){
              tmp_AoverEstmtShad_at_ij = Ashadow_sign*Ashadow[1][i][j];
              tmp_AoverEstmtShadSlp_at_ij = Ashadow_sign*Ashadow_slope[1][i][j];  
            }
            else{
              tmp_AoverEstmtShadSlp_at_ij = -Ashadow_sign*Ashadow_slope[1][i][j];
              tmp_AoverEstmtShad_at_ij = -Ashadow_sign*Ashadow[1][i][j] + std::fabs(tmp_AoverEstmtShadSlp_at_ij)*partitionSize[i];
            }
            Ashadow[1][i][j] = tmp_AoverEstmtShad_at_ij;
            Ashadow_slope[1][i][j] = tmp_AoverEstmtShadSlp_at_ij;
          }



          // long double tmp_BoverEstmtShad_at_ij = 0.;
          // long double tmp_BoverEstmtShadSlp_at_ij =0.;
          // if(Bshadow_info[2]!=0){
          //   tmp_BoverEstmtShad_at_ij = Bshadow_sign*(Aupdated?(- Bshadow[1][i][j] ):Bshadow[1][i][j]); 
          //   tmp_BoverEstmtShadSlp_at_ij = Bshadow_sign*(Aupdated?(- Bshadow_slope[1][i][j] ):Bshadow_slope[1][i][j]);
          //   LocalBshadow[1][i][j] = tmp_BoverEstmtShad_at_ij;
          //   LocalBshadow_slope[1][i][j] = tmp_BoverEstmtShadSlp_at_ij;            
          // }

          // long double tmp_AoverEstmtShad_at_ij = 0.;
          // long double tmp_AoverEstmtShadSlp_at_ij = 0.;
          // if(Ashadow_info[2]!=0){
          //   tmp_AoverEstmtShad_at_ij = Ashadow_sign*(Bupdated?(- Ashadow[1][i][j] ):Ashadow[1][i][j]);
          //   tmp_AoverEstmtShadSlp_at_ij = Ashadow_sign*(Bupdated? (- Ashadow_slope[1][i][j]  ):Ashadow_slope[1][i][j]);
          //   Ashadow[1][i][j] = tmp_AoverEstmtShad_at_ij;
          //   Ashadow_slope[1][i][j] = tmp_AoverEstmtShadSlp_at_ij;
          // }
          
          //std::cout << "                        ABshadow and slope updated at elm No. " << j << std::endl;

          NewOverEstmt_shadow[i][j] = tmp_AoverEstmtShad_at_ij + tmp_BoverEstmtShad_at_ij;
          if((tmp_AoverEstmtShadSlp_at_ij < 0.) != (tmp_BoverEstmtShadSlp_at_ij < 0.)){   
            long double tightener = std::min(std::fabs(tmp_AoverEstmtShadSlp_at_ij),std::fabs(tmp_BoverEstmtShadSlp_at_ij))*(partitionSize[i]);
            NewOverEstmt_shadow[i][j] -= tightener;
          }
          NewOverEstmt_shadow_slope[i][j] = tmp_AoverEstmtShadSlp_at_ij + tmp_BoverEstmtShadSlp_at_ij;                      
        }
      }
      else if( !Amat[i].empty() ){
#ifndef NOTTOTRACKSHADOW           
        std::cout << "                    !Amat[i].empty() && Bmat[i].empty(), processing row No. " << i << std::endl;
#endif
        NewOverEstmt[i].resize(_ndiv); 
        NewOverEstmt_slope[i].resize(_ndiv);
        NewOverEstmt_shadow[i].resize(_ndiv);            
        NewOverEstmt_shadow_slope[i].resize(_ndiv);           
        for( unsigned long long j=0; j<_ndiv; j++ ){
          long double const tmp_AoverEstmt_at_ij = Aupdated? NewAOverEstmt[i][j]:Op<T>::u(Amat[i][j]);
          long double const tmp_BoverEstmt_at_ij = Bupdated? NewBOverEstmt[i][j]:0.;
          long double const tmp_AoverEstmtSlp_at_ij = Aupdated? NewAOverEstmt_slope[i][j]:Aslope[i][j][1];
          long double const tmp_BoverEstmtSlp_at_ij = Bupdated? NewBOverEstmt_slope[i][j]:0.;          
          NewOverEstmt[i][j] = tmp_AoverEstmt_at_ij + tmp_BoverEstmt_at_ij;
          if((tmp_AoverEstmtSlp_at_ij < 0.) != (tmp_BoverEstmtSlp_at_ij < 0.)){
            long double tightener = std::min(std::fabs(tmp_AoverEstmtSlp_at_ij),std::fabs(tmp_BoverEstmtSlp_at_ij))*(partitionSize[i]);
            NewOverEstmt[i][j] -= tightener;            
          }
          NewOverEstmt_slope[i][j] = tmp_AoverEstmtSlp_at_ij + tmp_BoverEstmtSlp_at_ij;

          long double tmp_AoverEstmtShad_at_ij = 0.;
          long double tmp_AoverEstmtShadSlp_at_ij = 0.;
          if(Ashadow_info[2]!=0){
            if(!Bupdated){
              tmp_AoverEstmtShad_at_ij = Ashadow_sign*Ashadow[1][i][j];
              tmp_AoverEstmtShadSlp_at_ij = Ashadow_sign*Ashadow_slope[1][i][j];  
            }
            else{
              tmp_AoverEstmtShadSlp_at_ij = -Ashadow_sign*Ashadow_slope[1][i][j];
              tmp_AoverEstmtShad_at_ij = -Ashadow_sign*Ashadow[1][i][j] + std::fabs(tmp_AoverEstmtShadSlp_at_ij)*partitionSize[i];
            }
            Ashadow[1][i][j] = tmp_AoverEstmtShad_at_ij;
            Ashadow_slope[1][i][j] = tmp_AoverEstmtShadSlp_at_ij;
          }
          // if(Ashadow_info[2]!=0){
          //   tmp_AoverEstmtShad_at_ij = Ashadow_sign*(Bupdated?(- Ashadow[1][i][j] ):Ashadow[1][i][j]);
          //   tmp_AoverEstmtShadSlp_at_ij = Ashadow_sign*(Bupdated? (- Ashadow_slope[1][i][j]  ):Ashadow_slope[1][i][j]);
          //   Ashadow[1][i][j] = tmp_AoverEstmtShad_at_ij;
          //   Ashadow_slope[1][i][j] = tmp_AoverEstmtShadSlp_at_ij;
          // }

          NewOverEstmt_shadow[i][j] = tmp_AoverEstmtShad_at_ij;
          NewOverEstmt_shadow_slope[i][j] = tmp_AoverEstmtShadSlp_at_ij;
        }
      }
      else if( !Bmat[i].empty() ){
#ifndef NOTTOTRACKSHADOW             
        std::cout << "                    Amat[i].empty() && !Bmat[i].empty(), processing row No. " << i << std::endl;
#endif
        NewOverEstmt[i].resize(_ndiv); 
        NewOverEstmt_slope[i].resize(_ndiv);   
        NewOverEstmt_shadow[i].resize(_ndiv);            
        NewOverEstmt_shadow_slope[i].resize(_ndiv);                   
        for( unsigned long long j=0; j<_ndiv; j++ ){
          long double const tmp_AoverEstmt_at_ij = Aupdated? NewAOverEstmt[i][j]:0.;
          long double const tmp_BoverEstmt_at_ij = Bupdated? NewBOverEstmt[i][j]:Op<T>::u(Bmat[i][j]);
          long double const tmp_AoverEstmtSlp_at_ij = Aupdated? NewAOverEstmt_slope[i][j]:0.;
          long double const tmp_BoverEstmtSlp_at_ij = Bupdated? NewBOverEstmt_slope[i][j]:Bslope[i][j][1];          
          NewOverEstmt[i][j] = tmp_AoverEstmt_at_ij + tmp_BoverEstmt_at_ij;
          if((tmp_AoverEstmtSlp_at_ij < 0.) != (tmp_BoverEstmtSlp_at_ij < 0.)){
            long double tightener = std::min(std::fabs(tmp_AoverEstmtSlp_at_ij),std::fabs(tmp_BoverEstmtSlp_at_ij))*(partitionSize[i]);
            NewOverEstmt[i][j] -= tightener;            
          }          
          NewOverEstmt_slope[i][j] =  tmp_AoverEstmtSlp_at_ij + tmp_BoverEstmtSlp_at_ij;

          long double tmp_BoverEstmtShad_at_ij = 0.;
          long double tmp_BoverEstmtShadSlp_at_ij =0.;
          if(Bshadow_info[2]!=0){
            if(!Aupdated){
              tmp_BoverEstmtShad_at_ij = Bshadow_sign*Bshadow[1][i][j];
              tmp_BoverEstmtShadSlp_at_ij = Bshadow_sign*Bshadow_slope[1][i][j];             
            }
            else
            {
              tmp_BoverEstmtShadSlp_at_ij = -Bshadow_sign*Bshadow_slope[1][i][j];
              tmp_BoverEstmtShad_at_ij = -Bshadow_sign*Bshadow[1][i][j] + std::fabs(tmp_BoverEstmtShadSlp_at_ij)*partitionSize[i];       
            }
            LocalBshadow[1][i][j] = tmp_BoverEstmtShad_at_ij;
            LocalBshadow_slope[1][i][j] = tmp_BoverEstmtShadSlp_at_ij;                   
          }
          // if(Bshadow_info[2]!=0){
          //   tmp_BoverEstmtShad_at_ij = Bshadow_sign*(Aupdated?(- Bshadow[1][i][j] ):Bshadow[1][i][j]); 
          //   tmp_BoverEstmtShadSlp_at_ij = Bshadow_sign*(Aupdated?(- Bshadow_slope[1][i][j] ):Bshadow_slope[1][i][j]);
          //   LocalBshadow[1][i][j] = tmp_BoverEstmtShad_at_ij;
          //   LocalBshadow_slope[1][i][j] = tmp_BoverEstmtShadSlp_at_ij;            
          // }

          NewOverEstmt_shadow[i][j] = tmp_BoverEstmtShad_at_ij;
          NewOverEstmt_shadow_slope[i][j] = tmp_BoverEstmtShadSlp_at_ij;          
        }
      }          

      if(A_over_shadow && B_over_shadow ){     
        //const auto [tmp_min,tmp_max] = std::minmax_element(shadow[1][i].begin(),shadow[1][i].end());
        bool updateShadowCMin = false;
        if( !Ashadow[1][i].empty() ){
          const auto tmp_minA = std::min_element(Ashadow[1][i].begin(),Ashadow[1][i].end()); 
          shadowAMin += (*tmp_minA); 
          updateShadowCMin = true;
  
        }    
        if( !LocalBshadow[1][i].empty() ){   
          const auto tmp_minB = std::min_element(LocalBshadow[1][i].begin(),LocalBshadow[1][i].end());    
          shadowBMin += (*tmp_minB);   
          updateShadowCMin = true;
   
        }
        if(updateShadowCMin){
          const auto tmp_minC = std::min_element(NewOverEstmt_shadow[i].begin(),NewOverEstmt_shadow[i].end());     
          shadowCMin += (*tmp_minC);
        }
      }
    }


    // Step 2.2: aggregation
    if(A_over_shadow && B_over_shadow ){
      if (shadowCMin > shadowBMin){
        if(shadowBMin > shadowAMin){
          NewOverEstmt_shadow.swap(Ashadow[1]);
          NewOverEstmt_shadow_slope.swap(Ashadow_slope[1]);
        }
        else{
          //shadowPointingToB = true;
          NewOverEstmt_shadow.swap(LocalBshadow[1]);
          NewOverEstmt_shadow_slope.swap(LocalBshadow_slope[1]);            
        }
       }
       else {
        if(shadowCMin > shadowAMin){ 
          NewOverEstmt_shadow.swap(Ashadow[1]);
          NewOverEstmt_shadow_slope.swap(Ashadow_slope[1]);
        }
      }  
    }
    // to keep simple we do not further aggregate the new shadow over-estimator
    if (Aupdated || Bupdated)
      overEstmtUpdated = true;
  
  } 
  else{
#ifndef NOTTOTRACKSHADOW    
    std::cout << "        nottoUpdateOverEstmt" << std::endl;
#endif
    NewOverEstmt.resize(_nvar);
    NewOverEstmt_slope.resize(_nvar);
#ifndef NOTTOTRACKSHADOW      
    std::cout << "        Addition Over no shadow" << std::endl;
#endif


    for( unsigned int i=0; i<_nvar; i++ ){
      if( !Amat[i].empty() && !Bmat[i].empty() ){
#ifndef NOTTOTRACKSHADOW              
        std::cout << "                    !Amat[i].empty() && !Bmat[i].empty(), processing row No. " << i << std::endl;
#endif        
        NewOverEstmt[i].resize(_ndiv); 
        NewOverEstmt_slope[i].resize(_ndiv);
        for( unsigned long long j=0; j<_ndiv; j++ ){
          long double const tmp_AoverEstmt_at_ij = Op<T>::u(Amat[i][j]);
          long double const tmp_BoverEstmt_at_ij = Op<T>::u(Bmat[i][j]);
          long double const tmp_AoverEstmtSlp_at_ij = Aslope[i][j][1];
          long double const tmp_BoverEstmtSlp_at_ij = Bslope[i][j][1];
          
          NewOverEstmt[i][j] = tmp_AoverEstmt_at_ij + tmp_BoverEstmt_at_ij;
          if((tmp_AoverEstmtSlp_at_ij < 0.) != (tmp_BoverEstmtSlp_at_ij < 0.)){   
            long double tightener = std::min(std::fabs(tmp_AoverEstmtSlp_at_ij),std::fabs(tmp_BoverEstmtSlp_at_ij))*(partitionSize[i]);
            NewOverEstmt[i][j] -= tightener;
          }
          NewOverEstmt_slope[i][j] = tmp_AoverEstmtSlp_at_ij + tmp_BoverEstmtSlp_at_ij;
          //std::cout << std::setprecision(18) << i << "," << j << ": [ "  <<NewOverEstmt[i][j] << " ]   ";
        }
        //std::cout << std::endl;          
      }
      else if( !Amat[i].empty() ){
#ifndef NOTTOTRACKSHADOW         
        std::cout << "                    !Amat[i].empty() && Bmat[i].empty(), processing row No. " << i << std::endl;
#endif             
        NewOverEstmt[i].resize(_ndiv); 
        NewOverEstmt_slope[i].resize(_ndiv);
        for( unsigned long long j=0; j<_ndiv; j++ ){
          NewOverEstmt[i][j] = Op<T>::u(Amat[i][j]);
          NewOverEstmt_slope[i][j] =  Aslope[i][j][1];
        }
      }
      else if( !Bmat[i].empty() ){
#ifndef NOTTOTRACKSHADOW           
        std::cout << "                    Amat[i].empty() && !Bmat[i].empty(), processing row No. " << i << std::endl;
#endif
        NewOverEstmt[i].resize(_ndiv); 
        NewOverEstmt_slope[i].resize(_ndiv);   
        for( unsigned long long j=0; j<_ndiv; j++ ){
          NewOverEstmt[i][j] = Op<T>::u(Bmat[i][j]);
          NewOverEstmt_slope[i][j] =  Bslope[i][j][1];
        }
      }          
    }


  }







  NewUnderEstmt.resize(_nvar);
  NewUnderEstmt_slope.resize(_nvar);
  NewUnderEstmt_shadow.resize(_nvar);    
  NewUnderEstmt_shadow_slope.resize(_nvar);   

  if (toUpdateUnderEstmt){
#ifndef NOTTOTRACKSHADOW   
    std::cout << "        toUpdateUnderEstmt" << std::endl;
#endif
    std::vector<std::vector<long double>> NewAUnderEstmt; 
    std::vector<std::vector<long double>> NewAUnderEstmt_slope;     
    std::vector<std::vector<long double>> NewBUnderEstmt; 
    std::vector<std::vector<long double>> NewBUnderEstmt_slope;   

    NewAUnderEstmt.resize(_nvar);
    NewAUnderEstmt_slope.resize(_nvar);
    NewBUnderEstmt.resize(_nvar);
    NewBUnderEstmt_slope.resize(_nvar);
    


    bool Aupdated = false;
    bool Bupdated = false;

    if(A_under_shadow){
    // B + A_shadow -> New B 
#ifndef NOTTOTRACKSHADOW   
      std::cout << "            A_under_shadow" << std::endl;
#endif
      long double shadowEnhancedMin = 0.;
      const double shadow_sign = Ashadow_sign;
      if(shadow_sign > 0){
#ifndef NOTTOTRACKSHADOW   
        std::cout << "                shadow_sign > 0" << std::endl;
#endif
        for( unsigned int i=0; i<_nvar; i++ ){
          bool updateShadowEnhancedMin = false;  
          if( !Amat[i].empty() && !Bmat[i].empty() ){
#ifndef NOTTOTRACKSHADOW               
            std::cout << "                    !Amat[i].empty() && !Bmat[i].empty(), processing row No. " << i << std::endl;
#endif
            NewBUnderEstmt[i].resize(_ndiv); 
            NewBUnderEstmt_slope[i].resize(_ndiv);
            for( unsigned long long j=0; j<_ndiv; j++ ){
              NewBUnderEstmt[i][j] = Op<T>::l(Bmat[i][j]) + Ashadow[0][i][j];
              if((Ashadow_slope[0][i][j] < 0.) != (Bslope[i][j][0] < 0.)){   
                long double tightener = std::min(std::fabs(Ashadow_slope[0][i][j]),std::fabs(Bslope[i][j][0]))*(partitionSize[i]);
                NewBUnderEstmt[i][j] += tightener;
              }
              NewBUnderEstmt_slope[i][j] = Ashadow_slope[0][i][j] + Bslope[i][j][0];
            }
            updateShadowEnhancedMin = true;
          }
          else if( !Amat[i].empty() ){
#ifndef NOTTOTRACKSHADOW               
            std::cout << "                    !Amat[i].empty() && Bmat[i].empty(), processing row No. " << i << std::endl;
#endif
            NewBUnderEstmt[i].resize(_ndiv); 
            NewBUnderEstmt_slope[i].resize(_ndiv);              
            for( unsigned long long j=0; j<_ndiv; j++ ){
              NewBUnderEstmt[i][j] = Ashadow[0][i][j];
              NewBUnderEstmt_slope[i][j] = Ashadow_slope[0][i][j];// Bslope[i][j][0];
            }
            updateShadowEnhancedMin = true;
          }
          else if( !Bmat[i].empty() ){
#ifndef NOTTOTRACKSHADOW   
            std::cout << "                    Amat[i].empty() && !Bmat[i].empty(), processing row No. " << i << std::endl;
#endif
            NewBUnderEstmt[i].resize(_ndiv); 
            NewBUnderEstmt_slope[i].resize(_ndiv);            
            for( unsigned long long j=0; j<_ndiv; j++ ){
              NewBUnderEstmt[i][j] = Op<T>::l(Bmat[i][j]);
              NewBUnderEstmt_slope[i][j] = Bslope[i][j][0];
            }
            updateShadowEnhancedMin = true;
          }          
          if(updateShadowEnhancedMin){
            const auto tmp_min = std::min_element(NewBUnderEstmt[i].begin(),NewBUnderEstmt[i].end());     
            shadowEnhancedMin += (*tmp_min);
          }
        }
      }
      else if(shadow_sign < 0){
#ifndef NOTTOTRACKSHADOW   
        std::cout << "                shadow_sign < 0" << std::endl;
#endif
        for( unsigned int i=0; i<_nvar; i++ ){
          bool updateShadowEnhancedMin = false;
          if( !Amat[i].empty() && !Bmat[i].empty() ){
#ifndef NOTTOTRACKSHADOW               
            std::cout << "                    !Amat[i].empty() && !Bmat[i].empty(), processing row No. " << i << std::endl;
#endif
            NewBUnderEstmt[i].resize(_ndiv); 
            NewBUnderEstmt_slope[i].resize(_ndiv);
            for( unsigned long long j=0; j<_ndiv; j++ ){
              NewBUnderEstmt[i][j] = Op<T>::l(Bmat[i][j]) - Ashadow[0][i][j];
              if((Ashadow_slope[0][i][j] < 0.) == (Bslope[i][j][0] < 0.)){   
                long double tightener = std::min(std::fabs(Ashadow_slope[0][i][j]),std::fabs(Bslope[i][j][0]))*(partitionSize[i]);
                NewBUnderEstmt[i][j] += tightener;
              }
              NewBUnderEstmt_slope[i][j] = - Ashadow_slope[0][i][j] + Bslope[i][j][0];
            }
            updateShadowEnhancedMin = true;
          }
          else if( !Amat[i].empty() ){
#ifndef NOTTOTRACKSHADOW   
            std::cout << "                    !Amat[i].empty() && Bmat[i].empty(), processing row No. " << i << std::endl;
#endif
            NewBUnderEstmt[i].resize(_ndiv); 
            NewBUnderEstmt_slope[i].resize(_ndiv);            
            for( unsigned long long j=0; j<_ndiv; j++ ){
              NewBUnderEstmt[i][j] = -Ashadow[0][i][j] ;
              NewBUnderEstmt_slope[i][j] = -Ashadow_slope[0][i][j] ;
            }
            updateShadowEnhancedMin = true;
          }
          else if( !Bmat[i].empty() ){
#ifndef NOTTOTRACKSHADOW   
            std::cout << "                    Amat[i].empty() && !Bmat[i].empty(), processing row No. " << i << std::endl;
#endif
            NewBUnderEstmt[i].resize(_ndiv); 
            NewBUnderEstmt_slope[i].resize(_ndiv);            
            for( unsigned long long j=0; j<_ndiv; j++ ){
              NewBUnderEstmt[i][j] = Op<T>::l(Bmat[i][j]);
              NewBUnderEstmt_slope[i][j] = Bslope[i][j][0];
            }
            updateShadowEnhancedMin = true;
          }    

          if(updateShadowEnhancedMin){
            const auto tmp_min = std::min_element(NewBUnderEstmt[i].begin(),NewBUnderEstmt[i].end());     
            shadowEnhancedMin += (*tmp_min);
          }
        }        
      }
      else {
        std::cout << "error in binary addition aggr overestimator Ashadow when processing shadow sign for underestimator" <<std::endl;
        throw typename ISModel<T>::Exceptions( ISModel<T>::Exceptions::UNDEF ); 
      }
  
      if(shadowEnhancedMin > lambdaB + MC__ISM_COMPUTATION_TOL){
#ifndef NOTTOTRACKSHADOW   
        std::cout << "                Bupdated = true: " << std::setprecision(18) << shadowEnhancedMin << " , " << lambdaB <<std::endl;
#endif
        Bupdated = true;
      }  
    }

    if(B_under_shadow){
    // A + B_shadow -> New A 
#ifndef NOTTOTRACKSHADOW   
      std::cout << "            B_under_shadow" << std::endl;
#endif
      long double shadowEnhancedMin = 0.;
      const double shadow_sign = Bshadow_sign;
      if(shadow_sign > 0){
#ifndef NOTTOTRACKSHADOW   
        std::cout << "                shadow_sign > 0" << std::endl;
#endif
        for( unsigned int i=0; i<_nvar; i++ ){
          bool updateShadowEnhancedMin = false;
          if( !Amat[i].empty() && !Bmat[i].empty() ){
#ifndef NOTTOTRACKSHADOW   
            std::cout << "                    !Amat[i].empty() && !Bmat[i].empty(), processing row No. " << i << std::endl;
#endif            
            NewAUnderEstmt[i].resize(_ndiv); 
            NewAUnderEstmt_slope[i].resize(_ndiv);
            for( unsigned long long j=0; j<_ndiv; j++ ){
              NewAUnderEstmt[i][j] = Op<T>::l(Amat[i][j]) + Bshadow[0][i][j];
              if((Bshadow_slope[0][i][j] < 0.) != (Aslope[i][j][0] < 0.)){   
                long double tightener = std::min(std::fabs(Bshadow_slope[0][i][j]),std::fabs(Aslope[i][j][0]))*(partitionSize[i]);
                NewAUnderEstmt[i][j] += tightener;
              }
              NewAUnderEstmt_slope[i][j] = Bshadow_slope[0][i][j] + Aslope[i][j][0];
            }
            updateShadowEnhancedMin = true;
          }
          else if( !Amat[i].empty() ){
#ifndef NOTTOTRACKSHADOW   
            std::cout << "                    !Amat[i].empty() && Bmat[i].empty(), processing row No. " << i << std::endl;
#endif               
            NewAUnderEstmt[i].resize(_ndiv); 
            NewAUnderEstmt_slope[i].resize(_ndiv);            
            for( unsigned long long j=0; j<_ndiv; j++ ){
              NewAUnderEstmt[i][j] = Op<T>::l(Amat[i][j]);
              NewAUnderEstmt_slope[i][j] =  Aslope[i][j][0];
            }
            updateShadowEnhancedMin = true;
          }
          else if( !Bmat[i].empty() ){
#ifndef NOTTOTRACKSHADOW               
            std::cout << "                    Amat[i].empty() && !Bmat[i].empty(), processing row No. " << i << std::endl;
#endif                
            NewAUnderEstmt[i].resize(_ndiv); 
            NewAUnderEstmt_slope[i].resize(_ndiv);            
            for( unsigned long long j=0; j<_ndiv; j++ ){
              NewAUnderEstmt[i][j] = Bshadow[0][i][j];
              NewAUnderEstmt_slope[i][j] =  Bshadow_slope[0][i][j];
            }
            updateShadowEnhancedMin = true;
          }          
        
          if(updateShadowEnhancedMin){
            const auto tmp_min = std::min_element(NewAUnderEstmt[i].begin(),NewAUnderEstmt[i].end());     
            shadowEnhancedMin += (*tmp_min);
          }
        }
      }
      else if(shadow_sign < 0){
#ifndef NOTTOTRACKSHADOW               
        std::cout << "                shadow_sign < 0" << std::endl;
#endif         
        for( unsigned int i=0; i<_nvar; i++ ){
          bool updateShadowEnhancedMin = false;
          if( !Amat[i].empty() && !Bmat[i].empty() ){
#ifndef NOTTOTRACKSHADOW                 
            std::cout << "                    !Amat[i].empty() && !Bmat[i].empty(), processing row No. " << i << std::endl;
#endif                
            NewAUnderEstmt[i].resize(_ndiv); 
            NewAUnderEstmt_slope[i].resize(_ndiv);
            for( unsigned long long j=0; j<_ndiv; j++ ){
              NewAUnderEstmt[i][j] = Op<T>::l(Amat[i][j]) - Bshadow[0][i][j];
              if((Bshadow_slope[0][i][j] < 0.) == (Aslope[i][j][0] < 0.)){   
                long double tightener = std::min(std::fabs(Bshadow_slope[0][i][j]),std::fabs(Aslope[i][j][0]))*(partitionSize[i]);
                NewAUnderEstmt[i][j] += tightener;
              }
              NewAUnderEstmt_slope[i][j] = - Bshadow_slope[0][i][j] + Aslope[i][j][0];
            }
            updateShadowEnhancedMin = true;
          }
          else if( !Amat[i].empty() ){
#ifndef NOTTOTRACKSHADOW    
            std::cout << "                    !Amat[i].empty() && Bmat[i].empty(), processing row No. " << i << std::endl;
#endif                   
            NewAUnderEstmt[i].resize(_ndiv); 
            NewAUnderEstmt_slope[i].resize(_ndiv);            
            for( unsigned long long j=0; j<_ndiv; j++ ){
              NewAUnderEstmt[i][j] = Op<T>::l(Amat[i][j]);
              NewAUnderEstmt_slope[i][j] =  Aslope[i][j][0];
            }
            updateShadowEnhancedMin = true;
          }
          else if( !Bmat[i].empty() ){
#ifndef NOTTOTRACKSHADOW              
            std::cout << "                    Amat[i].empty() && !Bmat[i].empty(), processing row No. " << i << std::endl;
#endif
            NewAUnderEstmt[i].resize(_ndiv); 
            NewAUnderEstmt_slope[i].resize(_ndiv);            
            for( unsigned long long j=0; j<_ndiv; j++ ){
              NewAUnderEstmt[i][j] = -Bshadow[0][i][j];
              NewAUnderEstmt_slope[i][j] =  -Bshadow_slope[0][i][j];
            }
            updateShadowEnhancedMin = true;
          }          
        
          if(updateShadowEnhancedMin){
            const auto tmp_min = std::min_element(NewAUnderEstmt[i].begin(),NewAUnderEstmt[i].end());     
            shadowEnhancedMin += (*tmp_min);
          }  
        }        
      }
      else {
        std::cout << "error in binary addition aggr overestimator Bshadow when processing shadow sign for underestimator " << shadow_sign <<std::endl;
        throw typename ISModel<T>::Exceptions( ISModel<T>::Exceptions::UNDEF ); 
      }
  
      if(shadowEnhancedMin > lambdaA + MC__ISM_COMPUTATION_TOL){
#ifndef NOTTOTRACKSHADOW          
        std::cout << "                Aupdated = true: " << std::setprecision(18) << shadowEnhancedMin << " , " << lambdaA <<std::endl;
#endif        
        Aupdated = true;
      }
    }


    // Step 2.1: addition (underestimator) 
    long double shadowAMax(0.),shadowBMax(0.),shadowCMax(0.);
#ifndef NOTTOTRACKSHADOW      
    std::cout << "        Addition" << std::endl;
#endif
    for( unsigned int i=0; i<_nvar; i++ ){
      if( !Amat[i].empty() && !Bmat[i].empty() ){
#ifndef NOTTOTRACKSHADOW          
        std::cout << "                    !Amat[i].empty() && !Bmat[i].empty(), processing row No. " << i << std::endl;
#endif        
        NewUnderEstmt[i].resize(_ndiv); 
        NewUnderEstmt_slope[i].resize(_ndiv);
        NewUnderEstmt_shadow[i].resize(_ndiv);
        NewUnderEstmt_shadow_slope[i].resize(_ndiv);
        for( unsigned long long j=0; j<_ndiv; j++ ){
          long double const tmp_AoverEstmt_at_ij = Aupdated? NewAUnderEstmt[i][j]:Op<T>::l(Amat[i][j]);
          long double const tmp_BoverEstmt_at_ij = Bupdated? NewBUnderEstmt[i][j]:Op<T>::l(Bmat[i][j]);
          long double const tmp_AoverEstmtSlp_at_ij = Aupdated? NewAUnderEstmt_slope[i][j]:Aslope[i][j][0];
          long double const tmp_BoverEstmtSlp_at_ij = Bupdated? NewBUnderEstmt_slope[i][j]:Bslope[i][j][0];
          
          NewUnderEstmt[i][j] = tmp_AoverEstmt_at_ij + tmp_BoverEstmt_at_ij;
          if((tmp_AoverEstmtSlp_at_ij < 0.) != (tmp_BoverEstmtSlp_at_ij < 0.)){   
            long double tightener = std::min(std::fabs(tmp_AoverEstmtSlp_at_ij),std::fabs(tmp_BoverEstmtSlp_at_ij))*(partitionSize[i]);
            NewUnderEstmt[i][j] += tightener;
          }
          NewUnderEstmt_slope[i][j] = tmp_AoverEstmtSlp_at_ij + tmp_BoverEstmtSlp_at_ij;


          long double tmp_BoverEstmtShad_at_ij = 0.;
          long double tmp_BoverEstmtShadSlp_at_ij =0.;
          if(Bshadow_info[1]!=0){
            if(!Aupdated){
              tmp_BoverEstmtShad_at_ij = Bshadow_sign*Bshadow[0][i][j];
              tmp_BoverEstmtShadSlp_at_ij = Bshadow_sign*Bshadow_slope[0][i][j];             
            }
            else
            {
              tmp_BoverEstmtShadSlp_at_ij = -Bshadow_sign*Bshadow_slope[0][i][j];
              tmp_BoverEstmtShad_at_ij = -Bshadow_sign*Bshadow[0][i][j] - std::fabs(tmp_BoverEstmtShadSlp_at_ij)*partitionSize[i];       
            }
            LocalBshadow[0][i][j] = tmp_BoverEstmtShad_at_ij;
            LocalBshadow_slope[0][i][j] = tmp_BoverEstmtShadSlp_at_ij;                   
          }
          // std::cout << "    B shadow updated " << std::endl;
          long double tmp_AoverEstmtShad_at_ij = 0.;
          long double tmp_AoverEstmtShadSlp_at_ij = 0.;
          if(Ashadow_info[1]!=0){
            if(!Bupdated){
              tmp_AoverEstmtShad_at_ij = Ashadow_sign*Ashadow[0][i][j];
              tmp_AoverEstmtShadSlp_at_ij = Ashadow_sign*Ashadow_slope[0][i][j];  
            }
            else{
              tmp_AoverEstmtShadSlp_at_ij = -Ashadow_sign*Ashadow_slope[0][i][j];
              tmp_AoverEstmtShad_at_ij = -Ashadow_sign*Ashadow[0][i][j] - std::fabs(tmp_AoverEstmtShadSlp_at_ij)*partitionSize[i];
            }
            Ashadow[0][i][j] = tmp_AoverEstmtShad_at_ij;
            Ashadow_slope[0][i][j] = tmp_AoverEstmtShadSlp_at_ij;
          }
          // std::cout << "    A shadow updated " << std::endl;
          // long double const tmp_BoverEstmtShad_at_ij = Bshadow_sign*(Aupdated?(- Bshadow[0][i][j] ):Bshadow[0][i][j]);  
          // long double const tmp_AoverEstmtShad_at_ij = Ashadow_sign*(Bupdated?(- Ashadow[0][i][j] ):Ashadow[0][i][j]);
          // LocalBshadow[0][i][j] = tmp_BoverEstmtShad_at_ij;
          // Ashadow[0][i][j] = tmp_AoverEstmtShad_at_ij;
          // long double const tmp_BoverEstmtShadSlp_at_ij = Bshadow_sign*(Aupdated?(- Bshadow_slope[0][i][j] ):Bshadow_slope[0][i][j]);
          // long double const tmp_AoverEstmtShadSlp_at_ij = Ashadow_sign*(Bupdated? (- Ashadow_slope[0][i][j]  ):Ashadow_slope[0][i][j]);
          // LocalBshadow_slope[0][i][j] = tmp_BoverEstmtShadSlp_at_ij;
          // Ashadow_slope[0][i][j] = tmp_AoverEstmtShadSlp_at_ij;

          NewUnderEstmt_shadow[i][j] = tmp_AoverEstmtShad_at_ij + tmp_BoverEstmtShad_at_ij;
          if((tmp_AoverEstmtShadSlp_at_ij < 0.) != (tmp_BoverEstmtShadSlp_at_ij < 0.)){   
            long double tightener = std::min(std::fabs(tmp_AoverEstmtShadSlp_at_ij),std::fabs(tmp_BoverEstmtShadSlp_at_ij))*(partitionSize[i]);
            NewUnderEstmt_shadow[i][j] += tightener;
          }
          NewUnderEstmt_shadow_slope[i][j] = tmp_AoverEstmtShadSlp_at_ij + tmp_BoverEstmtShadSlp_at_ij;                      
        }
      }
      else if( !Amat[i].empty() ){
#ifndef NOTTOTRACKSHADOW          
        std::cout << "                    !Amat[i].empty() && Bmat[i].empty(), processing row No. " << i << std::endl;
#endif        
        NewUnderEstmt[i].resize(_ndiv); 
        NewUnderEstmt_slope[i].resize(_ndiv);
        NewUnderEstmt_shadow[i].resize(_ndiv);            
        NewUnderEstmt_shadow_slope[i].resize(_ndiv);           
        for( unsigned long long j=0; j<_ndiv; j++ ){
          long double const tmp_AoverEstmt_at_ij = Aupdated? NewAUnderEstmt[i][j]:Op<T>::l(Amat[i][j]);
          long double const tmp_BoverEstmt_at_ij = Bupdated? NewBUnderEstmt[i][j]:0.;
          long double const tmp_AoverEstmtSlp_at_ij = Aupdated? NewAUnderEstmt_slope[i][j]:Aslope[i][j][0];
          long double const tmp_BoverEstmtSlp_at_ij = Bupdated? NewBUnderEstmt_slope[i][j]:0.;          
          NewUnderEstmt[i][j] = tmp_AoverEstmt_at_ij + tmp_BoverEstmt_at_ij;
          if((tmp_AoverEstmtSlp_at_ij < 0.) != (tmp_BoverEstmtSlp_at_ij < 0.)){
            long double tightener = std::min(std::fabs(tmp_AoverEstmtSlp_at_ij),std::fabs(tmp_BoverEstmtSlp_at_ij))*(partitionSize[i]);
            NewUnderEstmt[i][j] += tightener;            
          }
          NewUnderEstmt_slope[i][j] = tmp_AoverEstmtSlp_at_ij + tmp_BoverEstmtSlp_at_ij;

          //if( A_under_shadow ){
            long double tmp_AoverEstmtShad_at_ij = 0.;
            long double tmp_AoverEstmtShadSlp_at_ij = 0.;
            if(Ashadow_info[1]!=0){
              if(!Bupdated){
                tmp_AoverEstmtShad_at_ij = Ashadow_sign*Ashadow[0][i][j];
                tmp_AoverEstmtShadSlp_at_ij = Ashadow_sign*Ashadow_slope[0][i][j];  
              }
              else{
                tmp_AoverEstmtShadSlp_at_ij = -Ashadow_sign*Ashadow_slope[0][i][j];
                tmp_AoverEstmtShad_at_ij = -Ashadow_sign*Ashadow[0][i][j] - std::fabs(tmp_AoverEstmtShadSlp_at_ij)*partitionSize[i];
              }
              Ashadow[0][i][j] = tmp_AoverEstmtShad_at_ij;
              Ashadow_slope[0][i][j] = tmp_AoverEstmtShadSlp_at_ij;
            }

            // if(Ashadow_info[1]!=0){
            //   tmp_AoverEstmtShad_at_ij = Ashadow_sign*(Bupdated?(- Ashadow[0][i][j] ):Ashadow[0][i][j]);
            //   tmp_AoverEstmtShadSlp_at_ij = Ashadow_sign*(Bupdated? (- Ashadow_slope[0][i][j]  ):Ashadow_slope[0][i][j]);
            //   Ashadow[0][i][j] = tmp_AoverEstmtShad_at_ij;
            //   Ashadow_slope[0][i][j] = tmp_AoverEstmtShadSlp_at_ij;
            // }
  
            // long double const tmp_AoverEstmtShad_at_ij = Ashadow_sign*(Bupdated?(- Ashadow[0][i][j] ):Ashadow[0][i][j]);
            // Ashadow[0][i][j] = tmp_AoverEstmtShad_at_ij;
            // long double const tmp_AoverEstmtShadSlp_at_ij = Ashadow_sign*(Bupdated? (- Ashadow_slope[0][i][j]  ):Ashadow_slope[0][i][j]);
            // Ashadow_slope[0][i][j] = tmp_AoverEstmtShadSlp_at_ij;
  
            NewUnderEstmt_shadow[i][j] = tmp_AoverEstmtShad_at_ij;
            NewUnderEstmt_shadow_slope[i][j] = tmp_AoverEstmtShadSlp_at_ij;
          //}
        }
      }
      else if( !Bmat[i].empty() ){
#ifndef NOTTOTRACKSHADOW            
        std::cout << "                    Amat[i].empty() && !Bmat[i].empty(), processing row No. " << i << std::endl;
#endif      
        NewUnderEstmt[i].resize(_ndiv); 
        NewUnderEstmt_slope[i].resize(_ndiv);   
        NewUnderEstmt_shadow[i].resize(_ndiv);            
        NewUnderEstmt_shadow_slope[i].resize(_ndiv);                   
        for( unsigned long long j=0; j<_ndiv; j++ ){
          long double const tmp_AoverEstmt_at_ij = Aupdated? NewAUnderEstmt[i][j]:0.;
          long double const tmp_BoverEstmt_at_ij = Bupdated? NewBUnderEstmt[i][j]:Op<T>::l(Bmat[i][j]);
          long double const tmp_AoverEstmtSlp_at_ij = Aupdated? NewAUnderEstmt_slope[i][j]:0.;
          long double const tmp_BoverEstmtSlp_at_ij = Bupdated? NewBUnderEstmt_slope[i][j]:Bslope[i][j][0];          
          NewUnderEstmt[i][j] = tmp_AoverEstmt_at_ij + tmp_BoverEstmt_at_ij;
          if((tmp_AoverEstmtSlp_at_ij < 0.) != (tmp_BoverEstmtSlp_at_ij < 0.)){
            long double tightener = std::min(std::fabs(tmp_AoverEstmtSlp_at_ij),std::fabs(tmp_BoverEstmtSlp_at_ij))*(partitionSize[i]);
            NewUnderEstmt[i][j] += tightener;            
          }          
          NewUnderEstmt_slope[i][j] =  tmp_AoverEstmtSlp_at_ij + tmp_BoverEstmtSlp_at_ij;

          // long double const tmp_BoverEstmtShad_at_ij = Bshadow_sign*(Aupdated?(- Bshadow[0][i][j] ):Bshadow[0][i][j]);  
          // LocalBshadow[0][i][j] = tmp_BoverEstmtShad_at_ij;
          // long double const tmp_BoverEstmtShadSlp_at_ij = Bshadow_sign*(Aupdated?(- Bshadow_slope[0][i][j] ):Bshadow_slope[0][i][j]);
          // LocalBshadow_slope[0][i][j] = tmp_BoverEstmtShadSlp_at_ij;
          //if( B_under_shadow ){
            long double tmp_BoverEstmtShad_at_ij = 0.;
            long double tmp_BoverEstmtShadSlp_at_ij =0.;
            if(Bshadow_info[1]!=0){
              if(!Aupdated){
                tmp_BoverEstmtShad_at_ij = Bshadow_sign*Bshadow[0][i][j];
                tmp_BoverEstmtShadSlp_at_ij = Bshadow_sign*Bshadow_slope[0][i][j];             
              }
              else
              {
                tmp_BoverEstmtShadSlp_at_ij = -Bshadow_sign*Bshadow_slope[0][i][j];
                tmp_BoverEstmtShad_at_ij = -Bshadow_sign*Bshadow[0][i][j] - std::fabs(tmp_BoverEstmtShadSlp_at_ij)*partitionSize[i];       
              }
              LocalBshadow[0][i][j] = tmp_BoverEstmtShad_at_ij;
              LocalBshadow_slope[0][i][j] = tmp_BoverEstmtShadSlp_at_ij;                   
            }

            // if(Bshadow_info[1]!=0){
            //   tmp_BoverEstmtShad_at_ij = Bshadow_sign*(Aupdated?(- Bshadow[0][i][j] ):Bshadow[0][i][j]); 
            //   tmp_BoverEstmtShadSlp_at_ij = Bshadow_sign*(Aupdated?(- Bshadow_slope[0][i][j] ):Bshadow_slope[0][i][j]);
            //   LocalBshadow[0][i][j] = tmp_BoverEstmtShad_at_ij;
            //   LocalBshadow_slope[0][i][j] = tmp_BoverEstmtShadSlp_at_ij;            
            // }
            // NewUnderEstmt_shadow[i][j] = tmp_BoverEstmtShad_at_ij;
            // NewUnderEstmt_shadow_slope[i][j] = tmp_BoverEstmtShadSlp_at_ij;                  
          //}

    
        }
      }          
    
      //const auto [tmp_min,tmp_max] = std::minmax_element(shadow[1][i].begin(),shadow[1][i].end());
      if( A_under_shadow && B_under_shadow){
        bool updateShadowCMax = false;
        if( !Ashadow[0][i].empty() ){
          const auto tmp_maxA = std::max_element(Ashadow[0][i].begin(),Ashadow[0][i].end()); 
          shadowAMax += (*tmp_maxA); 
          updateShadowCMax = true;
  
        }    
        if( !LocalBshadow[0][i].empty() ){   
          const auto tmp_maxB = std::max_element(LocalBshadow[0][i].begin(),LocalBshadow[0][i].end());    
          shadowBMax += (*tmp_maxB);   
          updateShadowCMax = true;
   
        }
        if(updateShadowCMax){
          const auto tmp_maxC = std::max_element(NewUnderEstmt_shadow[i].begin(),NewUnderEstmt_shadow[i].end());     
          shadowCMax += (*tmp_maxC);
        }
      }
    }
    
    // Step 2.2: aggregation
    if( A_under_shadow && B_under_shadow){
      if (shadowCMax < shadowBMax){
        if(shadowBMax < shadowAMax){
          NewUnderEstmt_shadow.swap(Ashadow[0]);
          NewUnderEstmt_shadow_slope.swap(Ashadow_slope[0]);
        }
        else{
          //shadowPointingToB = true;     
          NewUnderEstmt_shadow.swap(LocalBshadow[0]);
          NewUnderEstmt_shadow_slope.swap(LocalBshadow_slope[0]);                 
        }
       }
       else {
        if(shadowCMax < shadowAMax){
          NewUnderEstmt_shadow.swap(Ashadow[0]);
          NewUnderEstmt_shadow_slope.swap(Ashadow_slope[0]);
        }
      }  
    }
    // to keep simple we do not further aggregate the new shadow over-estimator
    if (Aupdated || Bupdated)
      underEstmtUpdated = true;
  
  }
  else{
#ifndef NOTTOTRACKSHADOW       
    std::cout << "        nottoUpdateUnderEstmt" << std::endl;
#endif    
    NewUnderEstmt.resize(_nvar);
    NewUnderEstmt_slope.resize(_nvar);
#ifndef NOTTOTRACKSHADOW       
    std::cout << "        Addition Under no shadow" << std::endl;
#endif    


    for( unsigned int i=0; i<_nvar; i++ ){
      
      if( !Amat[i].empty() && !Bmat[i].empty() ){
#ifndef NOTTOTRACKSHADOW           
        std::cout << "                    !Amat[i].empty() && !Bmat[i].empty(), processing row No. " << i << std::endl;
#endif        
        NewUnderEstmt[i].resize(_ndiv); 
        NewUnderEstmt_slope[i].resize(_ndiv);
        for( unsigned long long j=0; j<_ndiv; j++ ){
          long double const tmp_AoverEstmt_at_ij = Op<T>::l(Amat[i][j]);
          long double const tmp_BoverEstmt_at_ij = Op<T>::l(Bmat[i][j]);
          long double const tmp_AoverEstmtSlp_at_ij = Aslope[i][j][0];
          long double const tmp_BoverEstmtSlp_at_ij = Bslope[i][j][0];
          
          NewUnderEstmt[i][j] = tmp_AoverEstmt_at_ij + tmp_BoverEstmt_at_ij;
          if((tmp_AoverEstmtSlp_at_ij < 0.) != (tmp_BoverEstmtSlp_at_ij < 0.)){   
            long double tightener = std::min(std::fabs(tmp_AoverEstmtSlp_at_ij),std::fabs(tmp_BoverEstmtSlp_at_ij))*(partitionSize[i]);
            NewUnderEstmt[i][j] += tightener;
          }
          NewUnderEstmt_slope[i][j] = tmp_AoverEstmtSlp_at_ij + tmp_BoverEstmtSlp_at_ij;
        }
      }
      else if( !Amat[i].empty() ){
#ifndef NOTTOTRACKSHADOW   
        std::cout << "                    !Amat[i].empty() && Bmat[i].empty(), processing row No. " << i << std::endl;
#endif        
        NewUnderEstmt[i].resize(_ndiv); 
        NewUnderEstmt_slope[i].resize(_ndiv);
        for( unsigned long long j=0; j<_ndiv; j++ ){
          NewUnderEstmt[i][j] = Op<T>::l(Amat[i][j]);
          NewUnderEstmt_slope[i][j] =  Aslope[i][j][0];
        }
      }
      else if( !Bmat[i].empty() ){
#ifndef NOTTOTRACKSHADOW           
        std::cout << "                    Amat[i].empty() && !Bmat[i].empty(), processing row No. " << i << std::endl;
#endif
        NewUnderEstmt[i].resize(_ndiv); 
        NewUnderEstmt_slope[i].resize(_ndiv);   
        for( unsigned long long j=0; j<_ndiv; j++ ){
          NewUnderEstmt[i][j] = Op<T>::l(Bmat[i][j]);
          NewUnderEstmt_slope[i][j] =  Bslope[i][j][0];
        }
      }          
    }


  }
 

#ifndef NOTTOTRACKSHADOW   
  std::cout << "    to align" << std::endl;
#endif 
 // Step 3: assemble Amat, Aslope, Ashadow, Ashadow_slope
  if((!underEstmtUpdated) && (!overEstmtUpdated)){
    //In this case we do not need to align, but the shadow/shadow slopes should be prepared
#ifndef NOTTOTRACKSHADOW      
    std::cout << "        (!underEstmtUpdated) && (!overEstmtUpdated)" << std::endl;
#endif    
    Ashadow_info.resize(2*_nvar+5);
    Ashadow_info[0]=1.;
    Ashadow_info[1]=(double)toUpdateUnderEstmt;
    Ashadow_info[2]=(double)toUpdateOverEstmt;
    Ashadow_info[3]=0;
    Ashadow_info[4]=0;      

    for( unsigned int i=0; i<_nvar; i++ ){
      if (Amat[i].empty() && Bmat[i].empty()) continue;
#ifndef NOTTOTRACKSHADOW        
      std::cout << "            Amat[i].empty = "<< Amat[i].empty() << " Bmat[i].empty() = " << Bmat[i].empty() << ", processing row No. " << i << std::endl;
#endif
      Ashadow_info[2*i+5]=0;
      Ashadow_info[2*i+6]=0;
      bool makeSpaceForA = false;
      if(Amat[i].empty()){
        makeSpaceForA = true;
        Amat[i].resize(_ndiv);
        Aslope[i].resize(_ndiv);
      }   
      for( unsigned long long j=0; j<_ndiv; j++ ){
        if(NewOverEstmt[i][j] - NewUnderEstmt[i][j] > 0.)
          Amat[i][j] = T(NewUnderEstmt[i][j],NewOverEstmt[i][j]);
        else if(NewOverEstmt[i][j] - NewUnderEstmt[i][j] > -MC__ISM_COMPUTATION_TOL)
          Amat[i][j] = T(std::min(NewUnderEstmt[i][j],NewOverEstmt[i][j]),std::max(NewUnderEstmt[i][j],NewOverEstmt[i][j]));
        else {
          std::cout << "            Numerical error in add agrt alinging: ub- lb = " << std::setprecision(18) << NewOverEstmt[i][j] - NewUnderEstmt[i][j] << std::endl;
          throw typename ISModel<T>::Exceptions( ISModel<T>::Exceptions::UNDEF ); 
        }
        if(makeSpaceForA){
          Amat[i][j] = T(0.,0.);
          Aslope[i][j].resize(2);
        }
        Aslope[i][j][0] = NewUnderEstmt_slope[i][j];
        Aslope[i][j][1] = NewOverEstmt_slope[i][j];  
      }

    }
  
    Ashadow[0].swap(NewUnderEstmt_shadow);
    Ashadow[1].swap(NewOverEstmt_shadow);
    Ashadow_slope[0].swap(NewUnderEstmt_shadow_slope);
    Ashadow_slope[1].swap(NewOverEstmt_shadow_slope);

    //std::cout << "    End processing " << std::endl;



 }
 else{
#ifndef NOTTOTRACKSHADOW     
    std::cout << "        underEstmtUpdated || overEstmtUpdated" << std::endl;
#endif


  // We extract superposition components from the shadow estimators before we compute addimissible ranges for all rows


  //Get the maxima of all components of the input underestimator, stored in _L2[i]
  if(toUpdateUnderEstmt && underEstmtUpdated){
    long double sigma_o = 0.;       // <- the maximum of the input underestimator
    long double underEstmtShadow_lambda = 0.;
    //std::vector<long double>
    for( unsigned int i=0; i<_nvar; i++ ){
      if( NewUnderEstmt_shadow[i].empty() ) continue;
      long double _tmp_rowMin = NewUnderEstmt_shadow[i][0];
      long double _tmp_rowMax = _tmp_rowMin + std::fabs(NewUnderEstmt_shadow_slope[i][0])*partitionSize[i];
       for( unsigned int j=1; j<_ndiv; j++ ){
         _tmp_rowMax = std::max(_tmp_rowMax, NewUnderEstmt_shadow[i][j] + std::fabs(NewUnderEstmt_shadow_slope[i][j])*partitionSize[i]);
         _tmp_rowMin = std::min(_tmp_rowMin, NewUnderEstmt_shadow[i][j]);
       }
       _L2[i] = _tmp_rowMax;
       _L1[i] = _tmp_rowMin;
       sigma_o += _L2[i];
       underEstmtShadow_lambda += _L1[i];
    }
    // Only do the extraction if there exists nonzero explicit superposition components
    // bool toExtractUnderEstmtShadow = false;
    // for( unsigned int i=0; i<_nvar; i++ ){
    //   if(_L2[i]-_L1[i]+underEstmtShadow_lambda)
    //   {
    //     toExtractUnderEstmtShadow = true;
    //     break;
    //   }
      
    // } 
    auto const& f = [=]( const double& x ){ return std::max( x, 0. ); };
    auto const& fDerv = [=]( const double& x ){ return x > 0.? 1.: 0.; };  

    const long double shadow_global_offset = (1.0 - 1.0/((double) ndep))*sigma_o;
    for( unsigned int i=0; i<_nvar; i++ ){
      if( NewUnderEstmt_shadow[i].empty() ) continue;   
      
      const long double rowOffsetUnder  = - _L1[i] + underEstmtShadow_lambda;
      const long double rowOffsetShadow = - _L2[i] + sigma_o;
      for( unsigned int j=0; j<_ndiv; j++ ){  
        const long double zL = NewUnderEstmt_shadow[i][j];            
        const long double loPtUnderEstmt = zL + rowOffsetUnder;          
        const long double El = f( loPtUnderEstmt );
        const long double slope0_to_be_set = fDerv( loPtUnderEstmt )*NewUnderEstmt_shadow_slope[i][j];

        // compute shadow
        // Note that the function max(zi - zi^U + mu) - max(zi - zi^L + lambda) is no longer convex but still nondecreasing
        const long double delta_l = std::fabs(NewUnderEstmt_shadow_slope[i][j]*partitionSize[i]); 
        const long double loPtShadow = zL + rowOffsetShadow;//zL - _c2[i] + C2;
        const long double upPtShadow = loPtShadow + delta_l;//zL - _c2[i] + C2;        
        const long double shadowLoBnd = f( loPtShadow ) - El;
        //const long double shadowUpBnd = f( upPtShadow ) - f( loPtUnderEstmt + delta_l);
        NewUnderEstmt_shadow[i][j] = shadowLoBnd - f(shadow_global_offset);
        //shadow[2][i][j] = shadowUpBnd - shadow_global_offset;
        
        // The update of the derv needs special consideration
        const long double upBndUnderEstmt = f(loPtUnderEstmt + delta_l);
        long double slope0shadow_to_be_multiplied =  0.;
        if (shadowLoBnd == 0.)
          slope0shadow_to_be_multiplied = 0.;
        else if (f(loPtUnderEstmt+delta_l + 1e5*MC__ISM_COMPUTATION_TOL) == 0. ){
          slope0shadow_to_be_multiplied = std::min(fDerv(loPtShadow) - fDerv(loPtUnderEstmt),fDerv(upPtShadow) - fDerv(loPtUnderEstmt + delta_l));
          
        }
        else if (f(loPtUnderEstmt-1e2*MC__ISM_COMPUTATION_TOL) == 0.){
          if ( delta_l > 1e2 * MC__ISM_COMPUTATION_TOL)
            slope0shadow_to_be_multiplied = (f(upPtShadow) - upBndUnderEstmt - (f(loPtShadow) -f(loPtUnderEstmt)))/(delta_l);
          else 
            ;//slope0shadow_to_be_multiplied = fDerv(upPtShadow) - fDerv(loPtUnderEstmt + delta_l);
        }
        NewUnderEstmt_shadow_slope[i][j] = slope0shadow_to_be_multiplied * NewUnderEstmt_shadow_slope[i][j];
        NewUnderEstmt[i][j] = NewUnderEstmt[i][j] + El;  
        
        if(NewUnderEstmt_slope[i][j] < 0 != slope0_to_be_set < 0 ){
          long double tightener = std::min(std::fabs(NewUnderEstmt_slope[i][j]),std::fabs(slope0_to_be_set))*(partitionSize[i]);
          NewUnderEstmt[i][j] += tightener;
        }
        NewUnderEstmt_slope[i][j] = NewUnderEstmt_slope[i][j] + slope0_to_be_set;           
      }
    }
  }



  // if(toUpdateOverEstmt && overEstmtUpdated){
  //   long double sigma_u = 0.;       // <- the maximum of the input underestimator
  //   long double overEstmtShadow_mu = 0.;

  //   for( unsigned int i=0; i<_nvar; i++ ){
  //     if( NewOverEstmt_shadow[i].empty() ) continue;
  //     long double _tmp_rowMax = NewOverEstmt_shadow[i][0];
  //     long double _tmp_rowMin = _tmp_rowMin - std::fabs(NewOverEstmt_shadow_slope[i][0])*partitionSize[i];
  //      for( unsigned int j=1; j<_ndiv; j++ ){
  //        _tmp_rowMin = std::min(_tmp_rowMin, NewOverEstmt_shadow[i][j] - std::fabs(NewOverEstmt_shadow_slope[i][j])*partitionSize[i]);
  //        _tmp_rowMax = std::max(_tmp_rowMax, NewOverEstmt_shadow[i][j]);
  //      }
  //      _U2[i] = _tmp_rowMin;
  //      _U1[i] = _tmp_rowMax;
  //      sigma_u += _U2[i];
  //      overEstmtShadow_mu += _U1[i];
  //   }

  //   auto const& f = [=]( const double& x ){ return std::min( x, 0. ); };
  //   auto const& fDerv = [=]( const double& x ){ return x < 0.? 1.: 0.; };  

  //   const long double shadow_global_offset = (1.0 - 1.0/((double) ndep))*sigma_u;
  //   for( unsigned int i=0; i<_nvar; i++ ){
  //     if( NewOverEstmt_shadow[i].empty() ) continue;   
      
  //     const long double rowOffsetOver   = - _U1[i] + overEstmtShadow_mu;
  //     const long double rowOffsetShadow = - _U2[i] + sigma_u;
  //     // shadow[0][i].clear();
  //     // shadow[0][i].resize(_ndiv);      
  //     // shadow_slope[0][i].clear();
  //     // shadow_slope[0][i].resize(_ndiv);
  //     for( unsigned int j=0; j<_ndiv; j++ ){  
  //       const long double zU = NewOverEstmt_shadow[i][j];            
  //       const long double loPtUnderEstmt = zU + rowOffsetOver;           
  //       const long double Eu = f( loPtUnderEstmt );
  //       const long double slope1_to_be_set = fDerv( loPtUnderEstmt )*NewOverEstmt_shadow_slope[i][j];

  //       // compute shadow
  //       // Note that the function max(zi - zi^U + mu) - max(zi - zi^L + lambda) is no longer convex but still nondecreasing
  //       const long double delta_u = std::fabs(NewOverEstmt_shadow_slope[i][j]*partitionSize[i]); 
  //       const long double loPtShadow = zU + rowOffsetShadow;//zL - _c2[i] + C2;
  //       const long double upPtShadow = loPtShadow - delta_u;//zL - _c2[i] + C2;        
  //       const long double shadowLoBnd = f( loPtShadow ) - Eu;
  //       //const long double shadowUpBnd = f( upPtShadow ) - f( loPtUnderEstmt - delta_u);
  //       NewOverEstmt_shadow[i][j] = shadowLoBnd - f(shadow_global_offset);
  //       //shadow[2][i][j] = shadowUpBnd - shadow_global_offset;
        
  //       // The update of the derv needs special consideration
  //       const long double upBndUnderEstmt = f(loPtUnderEstmt - delta_u);
  //       long double slope1shadow_to_be_multiplied =  0.;
  //       if (shadowLoBnd == 0.)
  //         slope1shadow_to_be_multiplied = 0.;
  //       else if (f(loPtUnderEstmt-delta_u + 1e5*MC__ISM_COMPUTATION_TOL) == 0. ){
  //         slope1shadow_to_be_multiplied = 1.;
          
  //       }
  //       else if (f(loPtUnderEstmt+1e2*MC__ISM_COMPUTATION_TOL) == 0.){
  //         if ( delta_u > 1e2 * MC__ISM_COMPUTATION_TOL)
  //           slope1shadow_to_be_multiplied = ( f(loPtShadow) -f(loPtUnderEstmt) - ( f(upPtShadow) - upBndUnderEstmt ))/(delta_u);
  //         else 
  //           slope1shadow_to_be_multiplied = std::fabs(fDerv(upPtShadow) - fDerv(loPtUnderEstmt - delta_u));
  //       }
  //       NewOverEstmt_shadow_slope[i][j] = slope1shadow_to_be_multiplied * NewOverEstmt_shadow_slope[i][j];
  //       NewOverEstmt[i][j] = NewOverEstmt[i][j] + Eu;  
        
  //       if(NewOverEstmt_slope[i][j] < 0 != slope1_to_be_set < 0 ){
  //         long double tightener = std::min(std::fabs(NewOverEstmt_slope[i][j]),std::fabs(slope1_to_be_set))*(partitionSize[i]);
  //         NewOverEstmt[i][j] -= tightener;
  //       }
  //       NewOverEstmt_slope[i][j] = NewOverEstmt_slope[i][j] + slope1_to_be_set;           
  //     }
  //   }
  // }




    // compute the addimissible range for all rows
    long double sum_adms_range = 0.;       // <- the minimum of the input overestimator
    std::vector<long double> adms_range;
    adms_range.resize(_nvar);
    for( unsigned int i=0; i<_nvar; i++ ){
      if( Amat[i].empty() && Bmat[i].empty()) continue;
#ifndef NOTTOTRACKSHADOW      
      std::cout << "            !Amat[i].empty() && !Bmat[i].empty(), processing row No. " << i << std::endl;
#endif      
      long double tmp_overEstimator_at_ij     = NewOverEstmt[i][0];
      long double tmp_overEstimatorSlp_at_ij  = NewOverEstmt_slope[i][0];
      long double tmp_underEstimator_at_ij    = NewUnderEstmt[i][0];
      long double tmp_underEstimatorSlp_at_ij = NewUnderEstmt_slope[i][0];
      long double adms_range_row = 0.;    
      if (tmp_overEstimatorSlp_at_ij < 0 != tmp_underEstimatorSlp_at_ij < 0)       
        adms_range_row = tmp_overEstimator_at_ij - tmp_underEstimator_at_ij - partitionSize[i]*(std::fabs(tmp_overEstimatorSlp_at_ij) - std::fabs(tmp_underEstimatorSlp_at_ij)*partitionSize[i]);
      else 
        adms_range_row = std::min(tmp_overEstimator_at_ij-std::fabs(tmp_overEstimatorSlp_at_ij)*partitionSize[i]-tmp_underEstimator_at_ij,
                                  tmp_overEstimator_at_ij-std::fabs(tmp_underEstimatorSlp_at_ij)*partitionSize[i]-tmp_underEstimator_at_ij);

      for( unsigned int j=1; j<_ndiv; j++ ){
        //std::cout << "align 2 values in j = "<< j << std::endl;
        long double tmp_overEstimator_at_ij     = NewOverEstmt[i][j];
        long double tmp_overEstimatorSlp_at_ij  = NewOverEstmt_slope[i][j];
        long double tmp_underEstimator_at_ij    = NewUnderEstmt[i][j];
        long double tmp_underEstimatorSlp_at_ij = NewUnderEstmt_slope[i][j];
        if (tmp_overEstimatorSlp_at_ij < 0 != tmp_underEstimatorSlp_at_ij < 0)       
          adms_range_row = std::min(adms_range_row, tmp_overEstimator_at_ij - tmp_underEstimator_at_ij - partitionSize[i]*(std::fabs(tmp_overEstimatorSlp_at_ij) + std::fabs(tmp_underEstimatorSlp_at_ij)));
        else
          adms_range_row = std::min(adms_range_row,std::min(tmp_overEstimator_at_ij-std::fabs(tmp_overEstimatorSlp_at_ij)*partitionSize[i]-tmp_underEstimator_at_ij,
                                    tmp_overEstimator_at_ij-std::fabs(tmp_underEstimatorSlp_at_ij)*partitionSize[i]-tmp_underEstimator_at_ij));
      }
      adms_range[i]= adms_range_row;
      sum_adms_range += adms_range[i];
    }
    if (sum_adms_range < - 1e2*(_nvar)*MC__ISM_COMPUTATION_TOL){
      //if(std::fabs(sum_adms_range) < 1e-2*(muA+muB - lambdaA - lambdaB)){
      //  sum_adms_range = 1e-8*(muA+muB - lambdaA - lambdaB);
      //}
      //else{
        for( unsigned int i=0; i<_nvar; i++ ){
          if( NewOverEstmt[i].empty()) continue;
          for( unsigned int j=1; j<_ndiv; j++ ){
            std::cout << std::setprecision(18) << i << "," << j << ": [ " <<NewUnderEstmt[i][j] << " , " <<NewOverEstmt[i][j] << " ]   ";
          }
          std::cout << std::endl;
        }  
        std::cout << "error in aggregating in preprocessing relu(x) "<< std::setprecision(18) << sum_adms_range << std::endl; 
        throw typename ISModel<T>::Exceptions( ISModel<T>::Exceptions::UNDEF );   
      //}
    }
#ifndef NOTTOTRACKSHADOW    
    std::cout << "                addimissible range finished: spare range = " << sum_adms_range <<  std::endl;
#endif
    Ashadow_info.resize(2*_nvar+5);
    Ashadow_info[0]=1.;
    Ashadow_info[1]=(double)toUpdateUnderEstmt;
    Ashadow_info[2]=(double)toUpdateOverEstmt;
    Ashadow_info[3]=0;
    Ashadow_info[4]=0;      
#ifndef NOTTOTRACKSHADOW
    std::cout << "                assembly the mat and slope" << std::endl;
#endif    
    // assembly the mat and slope
    const long double sumAdmsRange_over_ndep = sum_adms_range/ndep;
    for( unsigned int i=0; i<_nvar; i++ ){
      if (Amat[i].empty() && Bmat[i].empty()) continue;
      Ashadow_info[2*i+5]=0;
      Ashadow_info[2*i+6]=0;
      const long double rowOffset_over = - adms_range[i] + sumAdmsRange_over_ndep;
      for( unsigned long long j=0; j<_ndiv; j++ ){
        long double tmp_overEstimator_at_ij     = NewOverEstmt[i][j] + rowOffset_over;
        long double tmp_underEstimator_at_ij    = NewUnderEstmt[i][j];
        tmp_underEstimator_at_ij    = std::min(tmp_underEstimator_at_ij,tmp_overEstimator_at_ij);         //<- to avoid NaN in filib

        Amat[i][j] = T(tmp_underEstimator_at_ij,tmp_overEstimator_at_ij);
        Aslope[i][j][0] = NewUnderEstmt_slope[i][j];
        Aslope[i][j][1] = NewOverEstmt_slope[i][j];  
      }

    }
#ifndef NOTTOTRACKSHADOW
    std::cout << "                assemblyed" << std::endl;
#endif
    Ashadow[0].swap(NewUnderEstmt_shadow);
    Ashadow[1].swap(NewOverEstmt_shadow);
    Ashadow_slope[0].swap(NewUnderEstmt_shadow_slope);
    Ashadow_slope[1].swap(NewOverEstmt_shadow_slope);


 }

  //std::cout << "    End aggrt " << std::endl;
}
 

template <typename T>
inline
void
ISModel<T>::_asym_slope_relu_shadow
( std::vector<std::vector<T>>& mat, std::vector<std::vector<std::vector<long double>>>& slope, const std::vector<long double>& partitionSize, unsigned const& ndep, 
  std::vector<std::vector<std::vector<long double>>>& shadow, std::vector<std::vector<std::vector<long double>>>& shadow_slope, std::vector<double>& shadow_info)
const
{
#ifndef NOTTOTRACKSHADOW  
  std::cout << "relu_shadow" << std::endl;
#endif

  assert( !mat.empty() );
  if(shadow_info.size()<2){
#ifndef NOTTOTRACKSHADOW 
    std::cout << "    shadow_info.size()<2" << std::endl;
#endif
    shadow_info.resize(2*_nvar+5);
    shadow_info[0]=1.;
    shadow_info[1]=1;
    shadow_info[2]=0;
    shadow_info[3]=0;
    shadow_info[4]=0; 
    for(unsigned int i=0;i<_nvar;i++){
      if( mat[i].empty() ) continue;
      shadow_info[5+2*i]=0;
      shadow_info[6+2*i]=0;
      shadow[0][i].resize(_ndiv);    
      shadow_slope[0][i].resize(_ndiv);
      shadow[2][i].resize(_ndiv);    
      shadow_slope[2][i].resize(_ndiv);
    }   
    return _asym_slope_relu_ws(mat,slope,partitionSize,ndep,shadow,shadow_slope);
  }

#ifndef NOTTOTRACKSHADOW     
  std::cout << "    shadow_info.size()>1" << std::endl;
#endif
  const bool overShandowEnhancer  = shadow_info[2];
  const bool underShandowEnhancer = shadow_info[1];
  const double shadow_sign = shadow_info[0];
  bool overEstmtUpdated  = false;
  bool underEstmtUpdated = false; 
  bool overEstmtElimit     = false;
  bool underEstmtElimit    = false;

  T bnd = _B( mat, 1 );
  const long double mu = Op<T>::u(bnd);
  const long double lambda = Op<T>::l(bnd);
#ifndef NOTTOTRACKSHADOW 
  std::cout << "    overShandowEnhancer: " << overShandowEnhancer << std::endl;
  std::cout << "    underShandowEnhancer: " << underShandowEnhancer << std::endl;
  std::cout << "    shadow_sign: " << shadow_sign << std::endl;
#endif


  if( overShandowEnhancer ){  // Process the shadow-enhancement for the overestimator
#ifndef NOTTOTRACKSHADOW 
    std::cout << "        overShandowEnhancer: " << std::endl;
#endif
    long double shadowEnhancedMax = 0.;
    if(shadow_sign > 0){
#ifndef NOTTOTRACKSHADOW       
      std::cout << "            shadow_sign > 0: " << std::endl;
#endif
      for( unsigned int i=0; i<_nvar; i++ ){
        if( mat[i].empty() ) continue;
#ifndef NOTTOTRACKSHADOW  
        std::cout << "                processing the row NO. " << i << std::endl;      
#endif
        for( unsigned long long j=0; j<_ndiv; j++ ){
          shadow[1][i][j] += Op<T>::u(mat[i][j]);
          if((shadow_slope[1][i][j] < 0.) != (slope[i][j][1] < 0.)){   
            long double tightener = std::min(std::fabs(shadow_slope[1][i][j]),std::fabs(slope[i][j][1]))*(partitionSize[i]);
            shadow[1][i][j] -= tightener;
          }
          shadow_slope[1][i][j] += slope[i][j][1];
  
        }
        //const auto [tmp_min,tmp_max] = std::minmax_element(shadow[1][i].begin(),shadow[1][i].end());   
        const auto tmp_max = std::max_element(shadow[1][i].begin(),shadow[1][i].end());  
        shadow_info[2*i+6] = *tmp_max;      
        shadowEnhancedMax += (*tmp_max);
      }
    }
    else if(shadow_sign < 0){
#ifndef NOTTOTRACKSHADOW 
      std::cout << "            shadow_sign < 0: " << std::endl;
#endif
      for( unsigned int i=0; i<_nvar; i++ ){
        if( mat[i].empty() ) continue;
#ifndef NOTTOTRACKSHADOW 
        std::cout << "                processing the row NO. " << i << std::endl;      
#endif
        for( unsigned long long j=0; j<_ndiv; j++ ){
          shadow[1][i][j] =  Op<T>::u(mat[i][j]) - shadow[1][i][j];
          if((shadow_slope[1][i][j] < 0.) == (slope[i][j][1] < 0.)){   
            long double tightener = std::min(std::fabs(shadow_slope[1][i][j]),std::fabs(slope[i][j][1]))*(partitionSize[i]);
            shadow[1][i][j] -= tightener;
          }
          shadow_slope[1][i][j] = slope[i][j][1] - shadow_slope[1][i][j];
  
        }
        //const auto [tmp_min,tmp_max] = std::minmax_element(shadow[1][i].begin(),shadow[1][i].end());   
        const auto tmp_max = std::max_element(shadow[1][i].begin(),shadow[1][i].end());  
        shadow_info[2*i+6] = *tmp_max;      
        shadowEnhancedMax += (*tmp_max);
      } 
      shadow_info[0] = 1.; 
    }
    else {
      std::cout << "error in asym_slope_relu when processing shadow sign for overestimator" << shadow_sign << std::endl;
      throw typename ISModel<T>::Exceptions( ISModel<T>::Exceptions::UNDEF ); 
    }

    shadow_info[4] = shadowEnhancedMax;
    if(shadowEnhancedMax < 0.){
#ifndef NOTTOTRACKSHADOW 
      std::cout << "            overEstmtElimit = true " << std::endl;
#endif
      overEstmtElimit = true;
    }
    else if(shadowEnhancedMax < mu - MC__ISM_COMPUTATION_TOL){
#ifndef NOTTOTRACKSHADOW 
      std::cout << "            overEstmtUpdated = true " << std::endl;
#endif
      overEstmtUpdated = true;
    }
      

  }

  if(overEstmtElimit){
    shadow_info[0] = 200.;
#ifndef NOTTOTRACKSHADOW 
    std::cout << "        shadow_info[0] = 200? : " << shadow_info[0] << std::endl;
#endif
    return;
  }
    


  if( underShandowEnhancer ){  // Process the shadow-enhancement for the underestimator
#ifndef NOTTOTRACKSHADOW    
    std::cout << "        underShandowEnhancer: " << std::endl;
#endif
    long double shadowEnhancedMin = 0.;
    if(shadow_sign > 0){
#ifndef NOTTOTRACKSHADOW    
      std::cout << "            shadow_sign > 0: " << std::endl;
#endif
      for( unsigned int i=0; i<_nvar; i++ ){
        if( mat[i].empty() ) continue;
#ifndef NOTTOTRACKSHADOW    
        std::cout << "                processing the row NO. " << i << std::endl;   
#endif
        for( unsigned long long j=0; j<_ndiv; j++ ){
          shadow[0][i][j] += Op<T>::l(mat[i][j]);
          if((shadow_slope[0][i][j] < 0.) != (slope[i][j][0] < 0.)){   
            long double tightener = std::min(std::fabs(shadow_slope[0][i][j]),std::fabs(slope[i][j][0]))*(partitionSize[i]);
            shadow[0][i][j] += tightener;
          }
          shadow_slope[0][i][j] += slope[i][j][0];
        }
        //const auto [tmp_min,tmp_max] = std::minmax_element(shadow[0][i].begin(),shadow[0][i].end());     
        const auto tmp_min = std::min_element(shadow[0][i].begin(),shadow[0][i].end());  
        shadow_info[2*i+5] = *tmp_min;      
        shadowEnhancedMin += (*tmp_min);
      }
    }
    else if(shadow_sign < 0){
#ifndef NOTTOTRACKSHADOW    
      std::cout << "            shadow_sign < 0: " << std::endl;
#endif
      for( unsigned int i=0; i<_nvar; i++ ){
        if( mat[i].empty() ) continue;
#ifndef NOTTOTRACKSHADOW    
        std::cout << "                processing the row NO. " << i << std::endl;    
#endif
        for( unsigned long long j=0; j<_ndiv; j++ ){
          shadow[0][i][j] =  Op<T>::l(mat[i][j]) - shadow[0][i][j];
          if((shadow_slope[0][i][j] < 0.) == (slope[i][j][0] < 0.)){   
            long double tightener = std::min(std::fabs(shadow_slope[0][i][j]),std::fabs(slope[i][j][0]))*(partitionSize[i]);
            shadow[0][i][j] -= tightener;
          }
          shadow_slope[0][i][j] = slope[i][j][0] - shadow_slope[0][i][j];
  
        }
        //const auto [tmp_min,tmp_max] = std::minmax_element(shadow[0][i].begin(),shadow[0][i].end());   
        const auto tmp_min = std::min_element(shadow[0][i].begin(),shadow[0][i].end());  
        shadow_info[2*i+5] = *tmp_min;      
        shadowEnhancedMin += (*tmp_min);
      }
      shadow_info[0] = 1.; 
    }
    else {
      std::cout << "error in asym_slope_relu when processing shadow sign for underestimator" << shadow_sign <<std::endl;
      throw typename ISModel<T>::Exceptions( ISModel<T>::Exceptions::UNDEF ); 
    }


    shadow_info[3] = shadowEnhancedMin;
    if(shadowEnhancedMin > 0.){
#ifndef NOTTOTRACKSHADOW    
      std::cout << "            underEstmtElimit = true " << std::endl;
#endif
      underEstmtElimit = true;
    }  
    // if (!underEstmtElimit){
    //   // Get the maxima of all components of the input underestimator, stored in _L2[i]
    //   long double sigma_o = 0.;       // <- the maximum of the input underestimator
    //   for( unsigned int i=0; i<_nvar; i++ ){
    //     if( mat[i].empty() ) continue;
    //     long double _tmp_row = ((long double)(Op<T>::l( mat[i][0]) + std::fabs(slope[i][0][0]*partitionSize[i])));
    //      for( unsigned int j=1; j<_ndiv; j++ ){
    //        _tmp_row = std::max(_tmp_row, ((long double)(Op<T>::l( mat[i][j]) + std::fabs(slope[i][j][0]*partitionSize[i]))));
    //      }
    //      _L2[i]= _tmp_row;
    //      std::max(_L2[i] - _L1[i] + lambda,0.);
    //      sigma_o += _L2[i];
    //   }

    //   std::max(_L2[i] - _L1[i] + lambda,0.);
    //   underEstmtUpdated = true;
    // }
    if(shadowEnhancedMin > lambda + MC__ISM_COMPUTATION_TOL){
#ifndef NOTTOTRACKSHADOW    
      std::cout << "            underEstmtUpdated = true " << std::endl;
#endif
      underEstmtUpdated = true;      
    }  

  }


  if(underEstmtElimit){ 
    // this means the underestimator can be updated to allow positive result, so we should only aggregate estimators and finish the processing
#ifndef NOTTOTRACKSHADOW    
    std::cout << "        underEstmtElimit " << std::endl;
#endif
    // compute the addimissible range for all rows
    long double sum_adms_range = 0.;       // <- the minimum of the input overestimator
    std::vector<long double> adms_range;
    adms_range.resize(_nvar);
    for( unsigned int i=0; i<_nvar; i++ ){
      if( mat[i].empty() ) continue;
#ifndef NOTTOTRACKSHADOW    
      std::cout << "                processing the row NO. " << i << std::endl;    
#endif
      long double tmp_overEstimator_at_ij     = overEstmtUpdated?shadow[1][i][0]:Op<T>::u(mat[i][0]);
      long double tmp_overEstimatorSlp_at_ij  = overEstmtUpdated?shadow_slope[1][i][0]:slope[i][0][1];
      long double tmp_underEstimator_at_ij    = shadow[0][i][0];
      long double tmp_underEstimatorSlp_at_ij = shadow_slope[0][i][0];
      long double adms_range_row = 0.;
      if (tmp_overEstimatorSlp_at_ij < 0 != tmp_underEstimatorSlp_at_ij < 0)       
        adms_range_row = tmp_overEstimator_at_ij - tmp_underEstimator_at_ij - partitionSize[i]*(std::fabs(tmp_overEstimatorSlp_at_ij) + std::fabs(tmp_underEstimatorSlp_at_ij));
      else 
        adms_range_row = std::min(tmp_overEstimator_at_ij-std::fabs(tmp_overEstimatorSlp_at_ij)*partitionSize[i]-tmp_underEstimator_at_ij,
                                  tmp_overEstimator_at_ij-std::fabs(tmp_underEstimatorSlp_at_ij)*partitionSize[i]-tmp_underEstimator_at_ij);

       for( unsigned int j=1; j<_ndiv; j++ ){
        long double tmp_overEstimator_at_ij     = overEstmtUpdated?shadow[1][i][j]:Op<T>::u(mat[i][j]);
        long double tmp_overEstimatorSlp_at_ij  = overEstmtUpdated?shadow_slope[1][i][j]:slope[i][j][1];
        long double tmp_underEstimator_at_ij    = shadow[0][i][j];
        long double tmp_underEstimatorSlp_at_ij = shadow_slope[0][i][j];        
        if (tmp_overEstimatorSlp_at_ij < 0 != tmp_underEstimatorSlp_at_ij < 0)       
          adms_range_row = std::min(adms_range_row, tmp_overEstimator_at_ij - tmp_underEstimator_at_ij - partitionSize[i]*(std::fabs(tmp_overEstimatorSlp_at_ij) + std::fabs(tmp_underEstimatorSlp_at_ij)));
        else
          adms_range_row = std::min(adms_range_row,std::min(tmp_overEstimator_at_ij-std::fabs(tmp_overEstimatorSlp_at_ij)*partitionSize[i]-tmp_underEstimator_at_ij,
                                    tmp_overEstimator_at_ij-std::fabs(tmp_underEstimatorSlp_at_ij)*partitionSize[i]-tmp_underEstimator_at_ij));
       }
       adms_range[i]= adms_range_row;
       sum_adms_range += adms_range[i];
    }
    if (sum_adms_range < -  1e2*(_nvar)*MC__ISM_COMPUTATION_TOL){
#ifndef NOTTOTRACKSHADOW    
      std::cout << "error in aggregating in preprocessing relu_shadow(x): case: underEstmtElimit" << std::setprecision(18) << sum_adms_range << std::endl; 
#endif
      throw typename ISModel<T>::Exceptions( ISModel<T>::Exceptions::UNDEF );   
    }
#ifndef NOTTOTRACKSHADOW    
    std::cout << "        addimissible range set " << std::endl;
#endif


    // assembly the mat and slope
    const long double sumAdmsRange_over_ndep = sum_adms_range/ndep;
    for( unsigned int i=0; i<_nvar; i++ ){ 
      if( mat[i].empty() ) continue;
#ifndef NOTTOTRACKSHADOW    
      std::cout << "                processing the row NO. " << i << std::endl;    
#endif
      const long double rowOffset_over = - adms_range[i] + sumAdmsRange_over_ndep;

      for( unsigned int j=0; j<_ndiv; j++ ){

        long double tmp_overEstimator_at_ij     = (overEstmtUpdated?shadow[1][i][j]:Op<T>::u(mat[i][j])) + rowOffset_over;
        long double tmp_overEstimatorSlp_at_ij  = overEstmtUpdated?shadow_slope[1][i][j]:slope[i][j][1];
        long double tmp_underEstimator_at_ij    = std::min(shadow[0][i][j],tmp_overEstimator_at_ij);    //<- to avoid NaN in filib
        long double tmp_underEstimatorSlp_at_ij = shadow_slope[0][i][j];  
        slope[i][j][0] = tmp_underEstimatorSlp_at_ij;
        slope[i][j][1] = tmp_overEstimatorSlp_at_ij;
        mat[i][j] = T(tmp_underEstimator_at_ij,tmp_overEstimator_at_ij);
      }

    }    

    // process shadow and shadowInfo
    // Although that information can be directly passed without outer-approximation,
    // because of the complexity, we employ a naive way that is to discard them


    shadow_info.resize(1);
    shadow_info[0]=1.;   
    for( unsigned int i=0; i<_nvar; i++ ){
#ifndef NOTTOTRACKSHADOW
      std::cout << "                reset shadow and shadowslopes, the row NO. " << i << std::endl;    
#endif
      for (unsigned int k = 0; k < 4; k++){
        shadow[k][i].clear();
        shadow[k][i].resize(0);
        shadow_slope[k][i].clear();
        shadow_slope[k][i].resize(0);           
      }
      shadow_info[5+2*i]=0;
      shadow_info[6+2*i]=0;          
    }

//     for( unsigned int i=0; i<_nvar; i++ ){
//       if( mat[i].empty() ) continue;
// #ifndef NOTTOTRACKSHADOW    
//       std::cout << "                processing shadowInfo the row NO. " << i << std::endl;    
// #endif
//       for( unsigned int j=0; j<_ndiv; j++ ){
//         shadow[0][i].clear();
//         shadow[0][i].resize(0);
//         shadow[1][i].clear();
//         shadow[1][i].resize(0);
//         shadow_slope[0][i].clear();
//         shadow_slope[0][i].resize(0);
//         shadow_slope[1][i].clear();
//         shadow_slope[1][i].resize(0);        
//       }
//     }
//     shadow_info.resize(1);
//     shadow_info[0] = 1.;


    return;
  }

  // Not let us aggregate 
  if(overEstmtUpdated || underEstmtUpdated){
#ifndef NOTTOTRACKSHADOW    
    std::cout << "        overEstmtUpdated || underEstmtUpdated " << std::endl;
#endif
    // compute the addimissible range for all rows
    long double sum_adms_range = 0.;       // <- the minimum of the input overestimator
    std::vector<long double> adms_range;
    adms_range.resize(_nvar);
    for( unsigned int i=0; i<_nvar; i++ ){
      if( mat[i].empty() ) continue;
#ifndef NOTTOTRACKSHADOW    
      std::cout << "                processing the row NO. " << i << std::endl;    
#endif
      long double tmp_overEstimator_at_ij     = overEstmtUpdated?shadow[1][i][0]:Op<T>::u(mat[i][0]);
      long double tmp_overEstimatorSlp_at_ij  = overEstmtUpdated?shadow_slope[1][i][0]:slope[i][0][1];
      long double tmp_underEstimator_at_ij    = underEstmtUpdated?shadow[0][i][0]:Op<T>::l(mat[i][0]);
      long double tmp_underEstimatorSlp_at_ij = underEstmtUpdated?shadow_slope[0][i][0]:slope[i][0][0];
      long double adms_range_row = 0.;
      if (tmp_overEstimatorSlp_at_ij < 0 != tmp_underEstimatorSlp_at_ij < 0)       
        adms_range_row = tmp_overEstimator_at_ij - tmp_underEstimator_at_ij - partitionSize[i]*(std::fabs(tmp_overEstimatorSlp_at_ij) + std::fabs(tmp_underEstimatorSlp_at_ij));
      else 
        adms_range_row = std::min(tmp_overEstimator_at_ij-std::fabs(tmp_overEstimatorSlp_at_ij)*partitionSize[i]-tmp_underEstimator_at_ij,
                                  tmp_overEstimator_at_ij-std::fabs(tmp_underEstimatorSlp_at_ij)*partitionSize[i]-tmp_underEstimator_at_ij);

      for( unsigned int j=1; j<_ndiv; j++ ){
        long double tmp_overEstimator_at_ij     = overEstmtUpdated?shadow[1][i][j]:Op<T>::u(mat[i][j]);
        long double tmp_overEstimatorSlp_at_ij  = overEstmtUpdated?shadow_slope[1][i][j]:slope[i][j][1];
        long double tmp_underEstimator_at_ij    = underEstmtUpdated?shadow[0][i][j]:Op<T>::l(mat[i][j]);
        long double tmp_underEstimatorSlp_at_ij = underEstmtUpdated?shadow_slope[0][i][j]:slope[i][j][0];        
        if (tmp_overEstimatorSlp_at_ij < 0 != tmp_underEstimatorSlp_at_ij < 0)       
          adms_range_row = std::min(adms_range_row, tmp_overEstimator_at_ij - tmp_underEstimator_at_ij - partitionSize[i]*(std::fabs(tmp_overEstimatorSlp_at_ij) + std::fabs(tmp_underEstimatorSlp_at_ij)));
        else
          adms_range_row = std::min(adms_range_row,std::min(tmp_overEstimator_at_ij-std::fabs(tmp_overEstimatorSlp_at_ij)*partitionSize[i]-tmp_underEstimator_at_ij,
                                    tmp_overEstimator_at_ij-std::fabs(tmp_underEstimatorSlp_at_ij)*partitionSize[i]-tmp_underEstimator_at_ij));

      }
      adms_range[i]= adms_range_row;
      sum_adms_range += adms_range[i];
    }
    if (sum_adms_range < - 1e2*(_nvar)*MC__ISM_COMPUTATION_TOL){
#ifndef NOTTOTRACKSHADOW    
      std::cout << "error in aggregating in preprocessing relu_shadow(x):  case: overEstmtUpdated || underEstmtUpdated" << std::setprecision(18) << sum_adms_range << std::endl; 
#endif 
      throw typename ISModel<T>::Exceptions( ISModel<T>::Exceptions::UNDEF );   
    }
#ifndef NOTTOTRACKSHADOW       
    std::cout << "        addimissible range set " << std::endl;
#endif    
    // assembly the mat and slope
    const long double sumAdmsRange_over_ndep = sum_adms_range/ndep;
    for( unsigned int i=0; i<_nvar; i++ ){
      if( mat[i].empty() ) continue;
#ifndef NOTTOTRACKSHADOW          
      std::cout << "                processing the row NO. " << i << std::endl;    
#endif
      const long double rowOffset_over = - adms_range[i] + sumAdmsRange_over_ndep;
      for( unsigned int j=0; j<_ndiv; j++ ){
        long double tmp_overEstimator_at_ij     = rowOffset_over + (overEstmtUpdated?shadow[1][i][j]:Op<T>::u(mat[i][j]));
        long double tmp_overEstimatorSlp_at_ij  = overEstmtUpdated?shadow_slope[1][i][j]:slope[i][j][1];
        long double tmp_underEstimator_at_ij    = underEstmtUpdated?shadow[0][i][j]:Op<T>::l(mat[i][j]);
        tmp_underEstimator_at_ij    = std::min(tmp_underEstimator_at_ij,tmp_overEstimator_at_ij);         //<- to avoid NaN in filib
        long double tmp_underEstimatorSlp_at_ij = underEstmtUpdated?shadow_slope[0][i][j]:slope[i][j][0];  
        slope[i][j][0] = tmp_underEstimatorSlp_at_ij;
        slope[i][j][1] = tmp_overEstimatorSlp_at_ij;
        mat[i][j] = T(tmp_underEstimator_at_ij,tmp_overEstimator_at_ij);
#ifdef MC__USE_FILIB         
        if(mat[i][j].isEmpty())
          std::cout << i << "," << j << ": " << std::setprecision(18)
                    << tmp_underEstimator_at_ij << " , " << tmp_overEstimator_at_ij;
#endif
      }

    }

    // process shadow and shadowInfo
    // Note that the information can be directly passed without outer-approximation,
    // because of the complexity, we employ a naive way that is to discard them
    
    shadow_info.resize(2*_nvar+5);
    shadow_info[0]=1.;
    shadow_info[1]=1;
    shadow_info[2]=0;
    shadow_info[3]=0;
    shadow_info[4]=0;      
    for( unsigned int i=0; i<_nvar; i++ ){
#ifndef NOTTOTRACKSHADOW
      std::cout << "                reset shadow and shadowslopes, the row NO. " << i << std::endl;    
#endif
      for (unsigned int k = 0; k < 4; k++){
        shadow[k][i].clear();
        shadow[k][i].resize(0);
        shadow_slope[k][i].clear();
        shadow_slope[k][i].resize(0);           
      }
      if( mat[i].empty() ) continue;
      shadow[0][i].resize(_ndiv);
      shadow[2][i].resize(_ndiv);      
      shadow_slope[0][i].resize(_ndiv);
      shadow_slope[2][i].resize(_ndiv);      
      shadow_info[5+2*i]=0;
      shadow_info[6+2*i]=0;          
    }
    return _asym_slope_relu_ws(mat,slope,partitionSize,ndep,shadow,shadow_slope);


 

  // // Get the maxima of all components of the input underestimator, stored in _L2[i]
  // long double sigma_o = 0.;       // <- the maximum of the input underestimator
  // for( unsigned int i=0; i<_nvar; i++ ){
  //   if( mat[i].empty() ) continue;
  //   long double _tmp_row = ((long double)(Op<T>::l( mat[i][0]) + std::fabs(slope[i][0][0]*partitionSize[i])));
  //    for( unsigned int j=1; j<_ndiv; j++ ){
  //      _tmp_row = std::max(_tmp_row, ((long double)(Op<T>::l( mat[i][j]) + std::fabs(slope[i][j][0]*partitionSize[i]))));
  //    }
  //    _L2[i]= _tmp_row;
  //    sigma_o += _L2[i];
  // }

  // // Get the minima of all components of the input overestimator, stored in _U2[i] 
  // long double sigma_u = 0.;       // <- the minimum of the input overestimator
  // for( unsigned int i=0; i<_nvar; i++ ){
  //   if( mat[i].empty() ) continue;
  //   long double _tmp_row = ((long double)(Op<T>::u( mat[i][0]) - std::fabs(slope[i][0][1]*partitionSize[i])));
  //    for( unsigned int j=1; j<_ndiv; j++ ){
  //      _tmp_row = std::min(_tmp_row, ((long double)(Op<T>::u( mat[i][j]) - std::fabs(slope[i][j][1]*partitionSize[i]) )));
  //    }
  //    _U2[i]= _tmp_row;
  //    sigma_u += _U2[i];
  // }

  }
  if( (!overEstmtUpdated) && (!underEstmtUpdated) ){
#ifndef NOTTOTRACKSHADOW
      std::cout << "    No enhancement processed "  << std::endl;    
#endif
    // process shadow and shadowInfo
    // Note that the information can be directly passed without outer-approximation,
    // because of the complexity, we employ a naive way that is to discard them
    
    shadow_info.resize(2*_nvar+5);
    shadow_info[0]=1.;
    shadow_info[1]=1;
    shadow_info[2]=0;
    shadow_info[3]=0;
    shadow_info[4]=0;      
    for( unsigned int i=0; i<_nvar; i++ ){
#ifndef NOTTOTRACKSHADOW
      std::cout << "                reset shadow and shadowslopes, the row NO. " << i << std::endl;    
#endif
      for (unsigned int k = 0; k < 4; k++){
        shadow[k][i].clear();
        shadow[k][i].resize(0);
        shadow_slope[k][i].clear();
        shadow_slope[k][i].resize(0);           
      }
      if( mat[i].empty() ) continue;
      shadow[0][i].resize(_ndiv);
      shadow[2][i].resize(_ndiv);      
      shadow_slope[0][i].resize(_ndiv);
      shadow_slope[2][i].resize(_ndiv);      
      shadow_info[5+2*i]=0;
      shadow_info[6+2*i]=0;          
    }
    return _asym_slope_relu_ws(mat,slope,partitionSize,ndep,shadow,shadow_slope);

  }


  std::cout << "this part of function relu_shadow should not be processed" << std::endl;
  throw typename ISModel<T>::Exceptions( ISModel<T>::Exceptions::UNDEF ); 
  // process relu
 
}





template <typename T>
template <typename PUNIV>
inline
void
ISModel<T>::_asymDCdecNTC
( std::vector<std::vector<T>>& mat, unsigned const& ndep, PUNIV const& f,
  double const& z_iflec, bool const cvx_ccv )
const
{
  assert( !mat.empty() );
  T bnd = _B( mat, 1 );
  return _asymDCdecNTC( mat, ndep, f, z_iflec, cvx_ccv, bnd );
}

template <typename T>
template <typename PUNIV>
inline
void
ISModel<T>::_asymDCdecNTC
( std::vector<std::vector<T>>& mat, unsigned const& ndep, PUNIV const& f,
  double const& z_iflec, bool const cvx_ccv, T const& bnd )
const
{
  // anchor points
  
  double _lambda( Op<T>::l(bnd)), _mu(Op<T>::u(bnd));
  double _range(_mu-_lambda);
  

  double f_iflec = f( z_iflec );
  double f_iflec_over_ndep = f_iflec / ndep;
  
  if(isequal( _range, 0. ) ){
      for( unsigned int i=0; i<_nvar; i++ ){
        if( mat[i].empty() ) continue;
        for( unsigned int j=0; j<_ndiv; j++ ){ 
          mat[i][j] = f_iflec_over_ndep * T(1.,1.);
        }
      }
  }
  else{

    auto const& fright = [=]( const double & x ){ return f(-x) - f_iflec; };  // For constructing estimators of the right part

  
    double f_lambda = f( _lambda );
    double fright_mu = fright( -_mu );
    //double sum_c1( 0. );
    //double f_iflec_lambda = f_lambda*(ndep-1) + f_iflec;
    for( unsigned int i=0; i<_nvar; i++ ){
      if( mat[i].empty() ) continue;
      //double _temp(0.);
      _r1[i] = ( _U1[i] - _L1[i] );                       // Delta_i
      _r2[i] = _r1[i]/_range;                             // theta_i
     
      //_c1[i] = _r2[i]*f_iflec + (1.0 - _r2[i])*f_lambda;  // a_i 
      //_temp = f( _lambda + _r1[i]);
      // sum_c1 += _temp;

      //_c1[i] = f( _lambda + _r1[i]);
      //sum_c1 += _c1[i];
      _c2[i] = _r2[i]*f_iflec + (1.0 - _r2[i])*f_lambda; //_c1[i] +  _r2[i] * f_iflec_lambda; // v_i
      _c1[i] = (1.0 - _r2[i])*fright_mu;

    }
    

    /*
    for( unsigned int i=0; i<_nvar; i++ ){
      if( mat[i].empty() ) continue;
      _c2[i] -= _r2[i]*sum_c1;                               // v_i
    }
    */

    //std::cout<< (_c1[0] - _c2[0])<<std::endl;
    //std::cout<< (_c1[1] - _c2[1])<<std::endl;
    for( unsigned int i=0; i<_nvar; i++ ){
      if( mat[i].empty() ) continue;   
      if( isequal( _r1[i], 0. ) ){
        for( unsigned int j=0; j<_ndiv; j++ )
          mat[i][j] = T(0.);
        continue;
      }
      else{
        double _theta_fiflec = _r2[i]*f_iflec;
        //std::cout<<_theta_fiflec<<std::endl;
        double _term_f_lambda = (1. -_r2[i])*f_lambda;
        double _term_fright_mu = _c1[i];//(1. -_r2[i])*fright_mu;
        if (cvx_ccv){
          for( unsigned int j=0; j<_ndiv; j++ ){ 
            double El = std::min(_c2[i], f( _lambda - _L1[i] + Op<T>::l(mat[i][j])) ) - _term_f_lambda ;
            double Du1 = _r2[i]*f( _mu + (Op<T>::u(mat[i][j]) - _U1[i])/_r2[i] );
            double Du2 = f( _mu - _U1[i] + Op<T>::u(mat[i][j])) + _theta_fiflec - f_iflec;
            double Du = std::min(_theta_fiflec , std::max(Du1,Du2));// + _c1[i] - _c2[i];
            //std::cout<<Du1<<" and "<<Du2<<" and "<<_theta_fiflec<<" and "<<Du<<" and "<<El<<std::endl;
            double Eu = std::max(_c1[i], fright( -_mu + _U1[i] - Op<T>::u(mat[i][j])) ) - _term_fright_mu ;
            double Dl1 = _r2[i]*fright( -_lambda + (-Op<T>::l(mat[i][j]) + _L1[i])/_r2[i] );
            double Dl2 = fright( -_lambda + _L1[i] - Op<T>::l(mat[i][j]));
            double Dl = std::max(0. , std::min(Dl1,Dl2)); 
            mat[i][j] = T( El+Dl, Du+Eu );
            //mat[i][j] = T( El, Du);   Used for evaluating the performance of the estimator for tanh(min(g,c))
          }
        }
        else{
          for( unsigned int j=0; j<_ndiv; j++ ){
            double Eu = std::max(_c2[i], f( _lambda - _L1[i] + Op<T>::l(mat[i][j])) ) - _term_f_lambda ;
            double Dl1 = _r2[i]*f( _mu + (Op<T>::u(mat[i][j]) - _U1[i])/_r2[i] );
            double Dl2 = f( _mu - _U1[i] + Op<T>::u(mat[i][j])) + _theta_fiflec - f_iflec;
            double Dl = std::max(_theta_fiflec , std::min(Dl1,Dl2));// + _c1[i] - _c2[i];
            //std::cout<<Du1<<" and "<<Du2<<" and "<<_theta_fiflec<<" and "<<Du<<" and "<<El<<std::endl;
            double El = std::min(_c1[i], fright( -_mu + _U1[i] - Op<T>::u(mat[i][j])) ) - _term_fright_mu;
            double Du1 = _r2[i]*fright( - _lambda + (-Op<T>::l(mat[i][j]) + _L1[i])/_r2[i] );
            double Du2 = fright( -_lambda + _L1[i] - Op<T>::l(mat[i][j]));
            double Du = std::min(0. , std::max(Du1,Du2));
            mat[i][j] = T( El+Dl, Du+Eu );
          }          
        }
        /*
        if( cvx ){
          T Zu = ( mat[i][j] - _c1[i] ) / _r1[i] * sum_r1 + C1;
          double Du = scal_r1 * ( std::max( f( Op<T>::l(Zu) ), f( Op<T>::u(Zu) ) ) - fopt ) + fopt_over_ndep;
          if( imid != ICUT ){
            T Zl = mat[i][j] - _c1[i] + C1;
            double El =  f( mid( Op<T>::l(Zl), Op<T>::u(Zl), zopt ) ) - fopt_over_ndep * (ndep-1.);//std::min( f( Op<T>::l(Zl) ), f( Op<T>::u(Zl) ) ) - fopt_over_ndep * (ndep-1.);
            mat[i][j] = T( El, Du );
          }
          else if( !options.DCDEC_USE ){
            mat[i][j] = T( fopt_over_ndep, Du );
          }
          else{
            auto const& flinc = [=]( const T& x ){ return Op<T>::l(x) < zopt? 0.: f( Op<T>::l(x) ) - fopt; };
            auto const& fldec = [=]( const T& x ){ return Op<T>::u(x) > zopt? 0.: f( Op<T>::u(x) ) - fopt; };
            double El = flinc( mat[i][j] - _c1[i] + C1 ) + fldec( mat[i][j] - _c2[i] + C2 ) + fopt_over_ndep;
            mat[i][j] = T( El, Du );
          }
        }
        else{
          T Zl = ( mat[i][j] - _c1[i] ) / _r1[i] * sum_r1 + C1;
          double Dl = scal_r1 * ( std::min( f( Op<T>::l(Zl) ), f( Op<T>::u(Zl) ) ) - fopt ) + fopt_over_ndep;
          if( imid != ICUT ){
            T Zu = mat[i][j] - _c1[i] + C1;
            double Eu = f( mid( Op<T>::l(Zu), Op<T>::u(Zu), zopt ) ) - fopt_over_ndep * (ndep-1.);
            mat[i][j] = T( Dl, Eu );
          }
          else if( !options.DCDEC_USE ){
            mat[i][j] = T( Dl, fopt_over_ndep );
          }
          else{
            auto const& fuinc = [=]( const T& x ){ return Op<T>::u(x) > zopt? 0.: f( Op<T>::u(x) ) - fopt; };
            auto const& fudec = [=]( const T& x ){ return Op<T>::l(x) < zopt? 0.: f( Op<T>::l(x) ) - fopt; };
            double Eu = fuinc( mat[i][j] - _c2[i] + C2 ) + fudec( mat[i][j] - _c1[i] + C1 ) + fopt_over_ndep;
            mat[i][j] = T( Dl, Eu );
          }
        }
        */
      }  
    }
  }

}

template <typename T>
template <typename PUNIV>
inline
void
ISModel<T>::_asymNTCsym
( std::vector<std::vector<T>>& mat, unsigned const& ndep, PUNIV const& f)
const
{
  assert( !mat.empty() );
  T bnd = _B( mat, 1 );
  return _asymNTCsym( mat, ndep, f, bnd );
}

template <typename T>
template <typename PUNIV>
inline
void
ISModel<T>::_asymNTCsym
( std::vector<std::vector<T>>& mat, unsigned const& ndep, PUNIV const& f, T const& bnd )
const
{
  // anchor points
  
  double _lambda( Op<T>::l(bnd)), _mu(Op<T>::u(bnd));
  double _range(_mu-_lambda);
  

  double f_lambda = f( _lambda );
  double f_lambda_over_ndep = f_lambda / ndep;
  
  if(isequal( _range, 0. ) ){
      for( unsigned int i=0; i<_nvar; i++ ){
        if( mat[i].empty() ) continue;
        for( unsigned int j=0; j<_ndiv; j++ ){ 
          mat[i][j] = f_lambda_over_ndep * T(1.,1.);
        }
      }
  }
  else{

  
    double f_mu = f( _mu );
    //double sum_c1( 0. );
    //double f_iflec_lambda = f_lambda*(ndep-1) + f_iflec;
    for( unsigned int i=0; i<_nvar; i++ ){
      if( mat[i].empty() ) continue;
      //double _temp(0.);
      _r1[i] = ( _U1[i] - _L1[i] );                       // Delta_i
      _r2[i] = _r1[i]/_range;                             // theta_i
     
      _c1[i] = (1.0 - _r2[i])*f_lambda; //_c1[i] +  _r2[i] * f_iflec_lambda; // v_i
      _c2[i] = (1.0 - _r2[i])*f_mu;

    }
    


    for( unsigned int i=0; i<_nvar; i++ ){
      if( mat[i].empty() ) continue;   
      if( isequal( _r1[i], 0. ) ){
        for( unsigned int j=0; j<_ndiv; j++ )
          mat[i][j] = T(0.);
        continue;
      }
      else{
        //std::cout<<_theta_fiflec<<std::endl;
        for( unsigned int j=0; j<_ndiv; j++ ){ 
          double El = f(Op<T>::l(mat[i][j]) - _L1[i] + _lambda)  - _c1[i] ;
//          double D = _r2[i]*f( _mu + (Op<T>::u(mat[i][j]) - _U1[i])/_r2[i] );
          double Du = _r2[i]*f( _mu + (Op<T>::u(mat[i][j]) - _U1[i])/_r2[i] );
          double Dl = _r2[i]*f( _lambda + (Op<T>::l(mat[i][j]) - _L1[i])/_r2[i] );
          double Eu = f(Op<T>::u(mat[i][j]) - _U1[i] +  _mu  )  - _c2[i];
          //std::cout<<Du1<<" and "<<Du2<<" and "<<_theta_fiflec<<" and "<<Du<<" and "<<El<<std::endl;
          
          // this makes use of the constructor of T, which always compares the value of the two imputs
//          mat[i][j] = T( std::min(El , D) , std::max(Eu , D) );   
          mat[i][j] = T( std::min(El , Dl) , std::max(Eu , Du) );
#ifdef MC__USE_FILIB 
          if(mat[i][j].isEmpty()){
            if(std::fabs(std::min(El , Dl) - std::max(Eu , Du)) <= MC__ISM_COMPUTATION_TOL)
                mat[i][j] = T( El );}
#ifdef FILIB__COMPUTATION_DEBUG             
          if(std::min(El , Dl) > std::max(Eu , Du))
             std::cout<<std::setprecision(18) << "warning: "<< std::min(El , Dl) << " > " << std::max(Eu , Du) <<std::endl;  
          if(mat[i][j].isEmpty())
             std::cout<<"_asymNTCsym:" <<i<<" , "<<j<<": "<<std::setprecision(18) << "warning: "<< std::min(El , Dl) << " > " << std::max(Eu , Du) <<std::endl;     
#endif
#endif
        }
      }  
    }
  }

}

template <typename T>
template <typename PUNIV, typename BNUIA>
inline
void
ISModel<T>::_asymNTCsym_slope
( std::vector<std::vector<T>>& mat,std::vector<std::vector<std::vector<long double>>>& slope, const std::vector<long double>& partitionSize,
    unsigned const& ndep, PUNIV const& f, BNUIA const& fDerv, double const& zopt, bool const cvx  )
const
{
  assert( !mat.empty() );
  T bnd = _B( mat, 1 );
  return _asymNTCsym_slope( mat, slope, partitionSize, ndep, f, fDerv, zopt, cvx, bnd );
}


template <typename T>
template <typename PUNIV, typename BNUIA>
inline
void
ISModel<T>::_asymNTCsym_slope
( std::vector<std::vector<T>>& mat, std::vector<std::vector<std::vector<long double>>>& slope, const std::vector<long double>& partitionSize,
    unsigned const& ndep, PUNIV const& f, BNUIA const& fDerv, double const& zopt, bool const cvx, T const& bnd )
const
{
  // anchor points
  
  double _lambda( Op<T>::l(bnd)), _mu(Op<T>::u(bnd));
  double _range(_mu-_lambda);
  

  double f_lambda = f( _lambda );
  
  if(isequal( _range, 0. ) ){
      double f_lambda_over_ndep = f_lambda/ndep;
      for( unsigned int i=0; i<_nvar; i++ ){
        if( mat[i].empty() ) continue;
        for( unsigned int j=0; j<_ndiv; j++ ){ 
          mat[i][j] = f_lambda_over_ndep * T(1.,1.);
        }
      }
  }
  else{

    double f_mu = f( _mu );

    for( unsigned int i=0; i<_nvar; i++ ){
      if( mat[i].empty() ) continue;
      //double _temp(0.);
      _r1[i] = ( _U1[i] - _L1[i] );                       // Delta_i
      _r2[i] = _r1[i]/_range;                             // theta_i
     
      _c1[i] = _r2[i]*f_lambda; // theta_i f_lambda
      _c2[i] = _r2[i]*f_mu;     // theta_i f_mu
    }
    

    //double lambda_minus_iflec = _lambda - zopt;
    //double mu_minus_iflec = _mu - zopt;


    for( unsigned int i=0; i<_nvar; i++ ){
      if( mat[i].empty() ) continue;   
      if( isequal( _r1[i], 0. ) ){
        for( unsigned int j=0; j<_ndiv; j++ ){
          mat[i][j] = T(0.);
          slope[i][j][0] = (long double)0.;
          slope[i][j][1] = (long double)0.;
          }
        continue;
      }
      else{
        //std::cout<<_theta_fiflec<<std::endl;
        for( unsigned int j=0; j<_ndiv; j++ ){ 
          long double zL = Op<T>::l(mat[i][j]),zU = Op<T>::u(mat[i][j]);
          zL = zL - _L1[i];
          zU = zU - _U1[i];
          long double delta_zL = std::fabs(slope[i][j][0]*partitionSize[i]);
          long double delta_zU = std::fabs(slope[i][j][1]*partitionSize[i]);

          long double linUndEstm_upPt = zL + delta_zL;
          long double linOveEstm_loPt = zU - delta_zU;

          long double slope0_to_be_multiplied = 0.;  
          long double slope1_to_be_multiplied = 0.;
          double lowerbound = 0.; 
          double upperbound = 0.;

          long double fastVar_UndEstm_upPt = linUndEstm_upPt / _r2[i] + _lambda;
          if (fastVar_UndEstm_upPt < ((long double) zopt) - 1e1*MC__ISM_COMPUTATION_TOL){
          //  both two candidate under-estimators are convex, and therefore the one that varies faster is always greater than the other
            slope0_to_be_multiplied = fDerv( zL + _lambda );
            lowerbound = f(zL + _lambda) - f_lambda + _c1[i];
          }
          else{
            long double f_fastVar_UndEstm_loPt = _r2[i] * f(zL / _r2[i] + _lambda);
            long double slowVar_UndEstm_loPt = zL + _lambda;
            long double f_slowVar_UndEstm_loPt = f(slowVar_UndEstm_loPt) - f_lambda + _c1[i];
            if (f_fastVar_UndEstm_loPt < f_slowVar_UndEstm_loPt - 1e1*MC__ISM_COMPUTATION_TOL){
              lowerbound = f_fastVar_UndEstm_loPt;
              if (delta_zL > 1e2*MC__ISM_COMPUTATION_TOL)
                slope0_to_be_multiplied = (_r2[i] * f(fastVar_UndEstm_upPt)  - f_fastVar_UndEstm_loPt)/delta_zL;
              //else 
                //slope0_to_be_multiplied = (long double) 0.;
            } 
            else{
              long double slowVar_UndEstm_upPt = linUndEstm_upPt + _lambda;
              long double f_linUndEstm_upPt = std::min( _r2[i] * f(fastVar_UndEstm_upPt), f(slowVar_UndEstm_upPt) - f_lambda + _c1[i]);
              long double slope0_candidate_1 = fDerv( slowVar_UndEstm_loPt );
              long double slope0_candidate_2 = slope0_candidate_1;

              if (delta_zL > 1e2*MC__ISM_COMPUTATION_TOL)
                slope0_candidate_2 = (f_linUndEstm_upPt - f_slowVar_UndEstm_loPt)/delta_zL;

              slope0_to_be_multiplied = std::min(slope0_candidate_1, slope0_candidate_2);
              lowerbound = std::min(f_fastVar_UndEstm_loPt,f_slowVar_UndEstm_loPt);
            }  
          }


          long double fastVar_OveEstm_loPt = linOveEstm_loPt / _r2[i] + _mu;
          if (fastVar_OveEstm_loPt > ((long double) zopt) + 1e2*MC__ISM_COMPUTATION_TOL){
          //  both two candidate under-estimators are concave
            slope1_to_be_multiplied = fDerv( zU + _mu );
            upperbound = f(zU + _mu) - f_mu + _c2[i];
          }
          else{
            long double f_fastVar_OveEstm_upPt = _r2[i] * f(zU / _r2[i] + _mu);
            long double slowVar_OveEstm_upPt = zU + _mu;
            long double f_slowVar_OveEstm_upPt = f(slowVar_OveEstm_upPt) - f_mu + _c2[i];
            if (f_fastVar_OveEstm_upPt > f_slowVar_OveEstm_upPt + 1e1*MC__ISM_COMPUTATION_TOL){
              upperbound = f_fastVar_OveEstm_upPt;
              if (delta_zU > 1e2*MC__ISM_COMPUTATION_TOL)
                slope1_to_be_multiplied = ( f_fastVar_OveEstm_upPt - _r2[i] * f(fastVar_OveEstm_loPt)  )/delta_zU;
              //else 
                //slope0_to_be_multiplied = (long double) 0.;
            } 
            else{
              long double slowVar_OveEstm_loPt = linOveEstm_loPt + _mu;
              long double f_linOveEstm_loPt = std::max( _r2[i] * f(fastVar_OveEstm_loPt), f(slowVar_OveEstm_loPt) - f_mu + _c2[i]);
              long double slope1_candidate_1 = fDerv( slowVar_OveEstm_upPt );
              long double slope1_candidate_2 = slope1_candidate_1;

              if (delta_zU > 1e2*MC__ISM_COMPUTATION_TOL)
                slope1_candidate_2 = (f_slowVar_OveEstm_upPt - f_linOveEstm_loPt)/delta_zU;

              slope1_to_be_multiplied = std::min(slope1_candidate_1, slope1_candidate_2);
              upperbound = std::max(f_fastVar_OveEstm_upPt,f_slowVar_OveEstm_upPt);
            }  
          }


//           double El = f(Op<T>::l(mat[i][j]) - _L1[i] + _lambda)  - _c1[i] ;
// //          double D = _r2[i]*f( _mu + (Op<T>::u(mat[i][j]) - _U1[i])/_r2[i] );
//           double Du = _r2[i]*f( _mu + (Op<T>::u(mat[i][j]) - _U1[i])/_r2[i] );
//           double Dl = _r2[i]*f( _lambda + (Op<T>::l(mat[i][j]) - _L1[i])/_r2[i] );
//           double Eu = f(Op<T>::u(mat[i][j]) - _U1[i] +  _mu  )  - _c2[i];
//           //std::cout<<Du1<<" and "<<Du2<<" and "<<_theta_fiflec<<" and "<<Du<<" and "<<El<<std::endl;
          
          // this makes use of the constructor of T, which always compares the value of the two imputs
//          mat[i][j] = T( std::min(El , D) , std::max(Eu , D) );   
          mat[i][j] = T( lowerbound , upperbound);
          slope[i][j][0] = slope0_to_be_multiplied * slope[i][j][0];
          slope[i][j][1] = slope1_to_be_multiplied * slope[i][j][1];
#ifdef MC__USE_FILIB 
          if(mat[i][j].isEmpty()){
            if(std::fabs(lowerbound - upperbound) <= MC__ISM_COMPUTATION_TOL)
                mat[i][j] = T( lowerbound );
                slope[i][j][0] = (long double)0.;
                slope[i][j][1] = (long double)0.;  
          }
#ifdef FILIB__COMPUTATION_DEBUG             
          if(lowerbound > upperbound)
             std::cout<<std::setprecision(18) << "warning: "<< lowerbound << " > " << upperbound <<std::endl;  
          if(mat[i][j].isEmpty())
             std::cout<<"_asymNTCsym_slope:" <<i<<" , "<<j<<": "<<std::setprecision(18) << "warning: "<< lowerbound << " > " << upperbound <<std::endl;     
#endif
#endif
        }
      }  
    }
  }

}


template <typename T>
template <typename PUNIV>
inline
double
ISModel<T>::_goldsect
( const double xL, const double xU, PUNIV const& f, const double TOL, const unsigned MAXIT )
const
{
  const double phi = 2.-(1.+std::sqrt(5.))/2.;
  const double fL = f(xL), fU = f(xU);
  if( fL*fU > 0 ) throw Exceptions( Exceptions::ROOT );
  const double xm = xU-phi*(xU-xL), fm = f(xm);
  unsigned iter = 0;
  return _goldsect_iter( iter, xL, fL, xm, fm, xU, fU, f, TOL, MAXIT );
}

template <typename T>
template <typename PUNIV>
inline
double
ISModel<T>::_goldsect_iter
( unsigned& iter, const double a, const double fa, const double b,
  const double fb, const double c, const double fc, PUNIV const& f,
  const double TOL, const unsigned MAXIT )
const
// a and c are the current bounds; the minimum is between them.
// b is a center point
{
  std::cout << "Iter #" << iter << std::endl;
  const double phi = 2.-(1.+std::sqrt(5.))/2.;
  bool b_then_x = ( c-b > b-a );
  double x = ( b_then_x? b+phi*(c-b): b-phi*(b-a) );
  if( std::fabs(c-a) < TOL*(std::fabs(b)+std::fabs(x)) || iter >= MAXIT ) return (c+a)/2.;
  double fx = f(x);
  if( b_then_x )
    return( fa*fx<0? _goldsect_iter( ++iter, a, fa, b, fb, x, fx, f, TOL, MAXIT ):
                     _goldsect_iter( ++iter, b, fb, x, fx, c, fc, f, TOL, MAXIT ) );
  return( fa*fb<0? _goldsect_iter( ++iter, a, fa, x, fx, b, fb, f, TOL, MAXIT ):
                   _goldsect_iter( ++iter, x, fx, b, fb, c, fc, f, TOL, MAXIT ) );
}

template <typename T>
inline
void ISModel<T>::_tanh
( std::vector<std::vector<T>>& mat, unsigned const& ndep )
const
{
  assert( !mat.empty() );

  // Bounds
  T bnd = _B( mat, 1 );
  double L( Op<T>::l(bnd) ), U( Op<T>::u(bnd) );

  // Asymetric inclusion
  auto const& f = [=]( const double& x ){ return std::tanh( x ); };
  if( U <= 0. )
    return _asym( mat, ndep, f, -DBL_MAX, true, bnd ); // convex part, min -INF
  if( L >= 0. )
    return _asym( mat, ndep, f, DBL_MAX, false, bnd ); // concave part, max INF
  else{
     
    //auto const& fcv = [=]( const double& x ){ return std::tanh( std::min( x, 0. ) ) + std::max( x, 0. ); };
    //std::cout<<"fcv"<<std::endl;
    //return _asym( mat, ndep, fcv, -DBL_MAX, true, bnd );
    
    //auto const& fmin = [=]( const double& x ){ return std::min(0.,std::tanh( x )); /*return std::tanh( -x )+2.0;*/ }; // for debugging concave-convex function with nonzero function value at the inflection point
    //return _asymDCdecNTC( mat, ndep, frev, 0., false, bnd );
    return _asymNTCsym( mat, ndep, f, bnd );
    //return _asymDCdecNTC( mat, ndep, f, 0., true, bnd );
  }
    
  /*
  if( ndep == 1 || !options.ENVEL_USE ){
    auto const& fcv = [=]( const double& x ){ return std::tanh( std::min( x, 0. ) ) + std::max( x, 0. ); };
    auto const& fcc = [=]( const double& x ){ return std::tanh( std::max( x, 0. ) ) + std::min( x, 0. ); };
    auto matcp1 = mat;
    auto matcp2 = mat;
    _dispmat( mat );
    _asym( matcp1, ndep, fcv, -DBL_MAX, true, bnd ); // convex term, min -INF
    _dispmat( matcp1 );
    _asym( matcp2, ndep, fcc,  DBL_MAX, false, bnd ); // concave term, max INF
    _dispmat( matcp2 );
    for( unsigned int i=0; i<_nvar; i++ ){
      if( mat[i].empty() ) continue;   
      for( unsigned int j=0; j<_ndiv; j++ )
        mat[i][j] = matcp1[i][j] + matcp2[i][j] - mat[i][j];
    }
  */
#if 0
    auto const& fcv = [=]( const double& x ){ return std::tanh( std::min( x, 0. ) ) - std::min( x, 0. ) + .5 * x; };
    auto const& fcc = [=]( const double& x ){ return std::tanh( std::max( x, 0. ) ) - std::max( x, 0. ) + .5 * x; };
    double zopt = 0.5 * std::log( 4. / ( 2. - std::sqrt(2.) ) - 1. );
    auto matcp = mat;
    _asym( mat, ndep, fcv, -zopt, true, bnd ); // convex term, min -zopt
    _asym( matcp, ndep, fcc, zopt, false, bnd ); // concave term, max zopt
    for( unsigned int i=0; i<_nvar; i++ ){
      if( matcp[i].empty() ) continue;   
      for( unsigned int j=0; j<_ndiv; j++ )
    //    mat[i][j] += matcp[i][j];
    //}
#endif
    //return;
  }
/*  
  auto const& fl = [=]( const double& x )
    { double tanhx = std::tanh(x); return (x-U)*(1-tanhx*tanhx)-(tanhx-std::tanh(U)); };
  double zcv = _goldsect( L, 0., fl, options.ROOT_TOL, options.ROOT_MAXIT );
  auto const& fcv = [=]( const double& x ) // convex envelope
    { double tanhcv = std::tanh(zcv); return x<zcv? std::tanh(x): tanhcv+(x-zcv)*(1-tanhcv*tanhcv); };
  auto const& fu = [=]( const double& x )
    { double tanhx = std::tanh(x); return (x-L)*(1-tanhx*tanhx)-(tanhx-std::tanh(L)); };
  double zcc = _goldsect( 0., U, fu, options.ROOT_TOL, options.ROOT_MAXIT );
  auto const& fcc = [=]( const double& x ) // concave envelope
    { double tanhcc = std::tanh(zcc); return x>zcc? std::tanh(x): tanhcc+(x-zcc)*(1-tanhcc*tanhcc); };
  auto matcp = mat;
  _asym( mat, ndep, fcv, -DBL_MAX, true, bnd ); // convex term, min -INF
  _dispmat( mat );
  _asym( matcp, ndep, fcc, DBL_MAX, false, bnd ); // concave term, max INF
  _dispmat( matcp );
  for( unsigned int i=0; i<_nvar; i++ ){
    if( mat[i].empty() ) continue;   
    for( unsigned int j=0; j<_ndiv; j++ )
      mat[i][j] = Op<T>::hull( mat[i][j], matcp[i][j] );
  }
}
*/


template <typename T>
inline
void ISModel<T>::_tanh
( std::vector<std::vector<T>>& mat,std::vector<std::vector<std::vector<long double>>>& slope, const std::vector<long double>& partitionSize, unsigned const& ndep )
const
{
  assert( !mat.empty() );

  // Bounds
  T bnd = _B( mat, 1 );
  double L( Op<T>::l(bnd) ), U( Op<T>::u(bnd) );

  // Asymetric inclusion
  auto const& f = [=]( const double& x ){ return std::tanh( x ); };
  auto const& fDerv = [=]( const double& x ){ return std::pow(1./(std::cosh( x )),2); };  
  if( U <= 0. )
    return _asym_slope( mat, slope, partitionSize, ndep, f, fDerv, -DBL_MAX, true, bnd ); // _asym( mat, ndep, f, -DBL_MAX, true, bnd ); // convex part, min -INF
  if( L >= 0. )
    return _asym_slope( mat, slope, partitionSize, ndep, f, fDerv, DBL_MAX, false, bnd ); // _asym( mat, ndep, f, DBL_MAX, false, bnd ); // concave part, max INF
  else{
    return _asymNTCsym_slope( mat, slope, partitionSize, ndep, f, fDerv, 0.0, false, bnd );
  }
    
}



template <typename T>
inline
void ISModel<T>::_pow
( std::vector<std::vector<T>>& mat, std::vector<std::vector<std::vector<long double>>>& slope, const std::vector<long double>& partitionSize, int const&iexp, unsigned const& ndep )
const
{
  assert( !mat.empty() && iexp && iexp!=1 );
  
  // Bounds
  T bnd = _B( mat, 1 );
  double L( Op<T>::l(bnd) ), U( Op<T>::u(bnd) );
  if ( iexp < 0 && L*U <= 0. )
    throw typename ISModel<T>::Exceptions( ISModel<T>::Exceptions::INV );
  //std::cout << "in pow asym: " << iexp << std::endl;
  // Asymetric inclusion
  auto const& f = [=]( const double& x ){ return std::pow( x, iexp ); };
  auto const& fDerv = [=]( const double& x ){ return iexp * std::pow( x, iexp - 1 ); };
  if( L > 0. && iexp < 0 )
    return _asym_slope( mat, slope, partitionSize, ndep, f, fDerv, DBL_MAX, true, bnd ); // convex part, min +INF 
  if( U < 0. && iexp < 0 && (-iexp)%2 )
    return _asym_slope( mat, slope, partitionSize, ndep, f, fDerv, -DBL_MAX, false, bnd ); // concave part, max -INF
  if( U < 0. && iexp < 0 )
    return _asym_slope( mat, slope, partitionSize, ndep, f, fDerv, -DBL_MAX, true, bnd ); // convex part, min -INF
  //std::cout << "in pow asym positive" << std::endl;
  if( iexp > 0 && !(iexp%2) ){//std::cout << "in pow asym even" << std::endl;
    return _asym_slope( mat, slope, partitionSize, ndep, f, fDerv, 0., true, bnd );} // convex part, min 0
  if( L >= 0. && iexp > 0 && iexp%2 ){//std::cout << "in pow asym 1" << std::endl;
    return _asym_slope( mat, slope, partitionSize, ndep, f, fDerv, -DBL_MIN, true, bnd );} // convex part, min 0
  if( U <= 0. && iexp > 0 && iexp%2 ){//std::cout << "in pow asym 2" << std::endl;
    return _asym_slope( mat, slope, partitionSize, ndep, f, fDerv, DBL_MIN, false, bnd );} // concave part, max 0
//std::cout << "in pow asym outer" << std::endl;
  auto const& fcv = [=]( const double& x ){ return std::pow( std::max( x, 0. ), iexp ); };
  auto const& fcc = [=]( const double& x ){ return std::pow( std::min( x, 0. ), iexp ); };
  auto const& fcvDerv = [=]( const double& x ){ return x > 0. ? ( (double) iexp)*std::pow(x,iexp-1): 0.; };
  auto const& fccDerv = [=]( const double& x ){ return x < 0. ? ( (double) iexp)*std::pow(x,iexp-1): 0.; };  
  auto matcp = mat;
  _asym_slope( mat, slope, partitionSize, ndep, fcv, fcvDerv, -DBL_MIN, true, bnd ); // convex part, min 0
  _asym_slope( matcp, slope, partitionSize, ndep, fcc, fccDerv, DBL_MIN, false, bnd ); // concave part, max 0
  for( unsigned int i=0; i<_nvar; i++ ){
    if( matcp[i].empty() ) continue;   
    for( unsigned int j=0; j<_ndiv; j++ )
      mat[i][j] += matcp[i][j];
  }
  return;
}

template <typename T>
inline
void ISModel<T>::_pow
( std::vector<std::vector<T>>& mat, int const&iexp, unsigned const& ndep )
const
{
  assert( !mat.empty() && iexp && iexp!=1 );
  
  // Bounds
  T bnd = _B( mat, 1 );
  double L( Op<T>::l(bnd) ), U( Op<T>::u(bnd) );
  if ( iexp < 0 && L*U <= 0. )
    throw typename ISModel<T>::Exceptions( ISModel<T>::Exceptions::INV );
  //std::cout << "in pow asym: " << iexp << std::endl;
  // Asymetric inclusion
  auto const& f = [=]( const double& x ){ return std::pow( x, iexp ); };
  if( L > 0. && iexp < 0 )
    return _asym( mat, ndep, f, DBL_MAX, true, bnd ); // convex part, min +INF
  if( U < 0. && iexp < 0 && (-iexp)%2 )
    return _asym( mat, ndep, f, -DBL_MAX, false, bnd ); // concave part, max -INF
  if( U < 0. && iexp < 0 )
    return _asym( mat, ndep, f, -DBL_MAX, true, bnd ); // convex part, min -INF
  //std::cout << "in pow asym positive" << std::endl;
  if( iexp > 0 && !(iexp%2) ){//std::cout << "in pow asym even" << std::endl;
    return _asym( mat, ndep, f, 0, true, bnd );} // convex part, min 0
  if( L >= 0. && iexp > 0 && iexp%2 ){//std::cout << "in pow asym 1" << std::endl;
    return _asym( mat, ndep, f, -DBL_MIN, true, bnd );} // convex part, min 0
  if( U <= 0. && iexp > 0 && iexp%2 ){//std::cout << "in pow asym 2" << std::endl;
    return _asym( mat, ndep, f, DBL_MIN, false, bnd );} // concave part, max 0
//std::cout << "in pow asym outer" << std::endl;
  auto const& fcv = [=]( const double& x ){ return std::pow( std::max( x, 0. ), iexp ); };
  auto const& fcc = [=]( const double& x ){ return std::pow( std::min( x, 0. ), iexp ); };
  auto matcp = mat;
  _asym( mat, ndep, fcv, -DBL_MIN, true, bnd ); // convex part, min 0
  _asym( matcp, ndep, fcc, DBL_MIN, false, bnd ); // concave part, max 0
  for( unsigned int i=0; i<_nvar; i++ ){
    if( matcp[i].empty() ) continue;   
    for( unsigned int j=0; j<_ndiv; j++ )
      mat[i][j] += matcp[i][j];
  }
  return;
}

template <typename T>
inline
void ISModel<T>::_inv
( std::vector<std::vector<T>>& mat, std::vector<std::vector<std::vector<long double>>>& slope, const std::vector<long double>& partitionSize, unsigned const& ndep )
const
{
  assert( !mat.empty() );

  // Bounds
  T bnd = _B( mat, 1 );
  double L( Op<T>::l(bnd) ), U( Op<T>::u(bnd) );
  if ( L*U <= 0. )
    throw Exceptions( Exceptions::INV );

  auto const& f = [=]( const double& x ){ return 1.0 / x; };
  auto const& fDerv = [=]( const double& x ){ return -1.0 / (x*x); };
  if( L > 0. )
    return _asym_slope( mat, slope, partitionSize, ndep, f, fDerv, DBL_MAX, true, bnd ); // concave term, max +INF
  else if( U < 0. )
    return _asym_slope( mat, slope, partitionSize, ndep, f, fDerv, -DBL_MAX, false, bnd ); // concave term, max +INF  
}

template <typename T>
inline
void ISModel<T>::_inv
( std::vector<std::vector<T>>& mat, unsigned const& ndep )
const
{
  assert( !mat.empty() );

  // Bounds
  T bnd = _B( mat, 1 );
  double L( Op<T>::l(bnd) ), U( Op<T>::u(bnd) );
  if ( L*U <= 0. )
    throw Exceptions( Exceptions::INV );
  
  // Asymetric inclusion
  if( options.ASYREM_USE ){
    auto const& f = [=]( const double& x ){ return  1 / x; };
    if( L > 0. )
      return _asym( mat, ndep, f, DBL_MAX, true, bnd ); // convex part, min +INF
    else if( U < 0. )
      return _asym( mat, ndep, f, -DBL_MAX, false, bnd ); // concave part, max -INF
  }

  // Central points
  double w( 0. ), den( L+U );
  for( unsigned int i=0; i<_nvar; i++ ){
    if( mat[i].empty() ) continue;
    _c1[i] = ( _L1[i] * U + _U1[i] * L ) / den;
    w += _c1[i];
  }
  
  // Remainder
  double rem( 0. ), s; //, t;
  for( unsigned int i=0; i<_nvar; i++ ){
    if( mat[i].empty() ) continue;
    if( L > 0. ){
      s = std::max( 1./std::fabs(w/(_c1[i]-_L1[i])-1.), 1./std::fabs(w/(_U1[i]-_c1[i])+1.) ); // <- NEEDS CHECKING 
      rem += s*std::fabs(U-w-_U1[i]+_c1[i]);
    }
    else if( U < 0. ){
      s = std::max( 1./std::fabs(w/(_c1[i]-_L1[i])-1.), 1./std::fabs(w/(_U1[i]-_c1[i])+1.) ); // <- NEEDS CHECKING 
      rem += s*std::fabs(L-w-_L1[i]+_c1[i]);
    }
  }
  rem /= std::fabs(w) * std::min( std::fabs(L), std::fabs(U) );

  // Interval matrix coefficients
  bnd = T(-1.,1.)*(rem/double(ndep)) - (ndep-1.)/double(ndep)*inv(w);
    for( unsigned int i=0; i<_nvar; i++ ){
    if( mat[i].empty() ) continue;
    for( unsigned int j=0; j<_ndiv; j++ )
      mat[i][j] = Op<T>::inv( (w-_c1[i])+mat[i][j] ) + bnd;
  }
}

// _sqr enhanced by using slopes
template <typename T>
inline
void ISModel<T>::_sqr
( std::vector<std::vector<T>>& mat, std::vector<std::vector<std::vector<long double>>>& slope, const std::vector<long double>& partitionSize, unsigned const& ndep )
const
{
  assert( !mat.empty() );

  // Bounds
  T bnd = _B( mat, 1 );
  
  auto const& f = [=]( const double& x ){ return x*x; };
  auto const& fDerv = [=]( const double& x ){ return 2.0*x; };
  //std::cout << "SQR_ASYM" << std::endl;
  return _asym_slope( mat, slope, partitionSize, ndep, f,fDerv, 0., true, bnd ); // convex term, min zero
}

template <typename T>
inline
void ISModel<T>::_sqr
( std::vector<std::vector<T>>& mat, unsigned const& ndep )
const
{
  assert( !mat.empty() );

  // Bounds
  T bnd = _B( mat, 1 );
  
  // Asymetric inclusion
  if( options.ASYREM_USE ){
    auto const& f = [=]( const double& x ){ return x*x; };
    //std::cout << "SQR_ASYM" << std::endl;
    return _asym( mat, ndep, f, 0., true, bnd ); // convex term, min zero
  }
  
  // Central points
  double w( 0. ), s( 0. );
  for( unsigned int i=0; i<_nvar; i++ ){
    if( mat[i].empty() ) continue;
    _c1[i] = 0.5 * ( _L1[i] + _U1[i] );
    w += _c1[i];
    _r1[i] = 0.5 * ( _U1[i] - _L1[i] );
    s += _r1[i];
  }

  // Remainder
  double rem( 0. );
  for( unsigned int i=0; i<_nvar; i++ ){
    if( mat[i].empty() ) continue;
    rem += _r1[i]*(s-_r1[i]);
  }

  // Interval matrix coefficients
  bnd = T(-1.,1.)*(rem/double(ndep)) - (ndep-1.)/double(ndep)*sqr(w);
  for( unsigned int i=0; i<_nvar; i++ ){
    if( mat[i].empty() ) continue;
    for( unsigned int j=0; j<_ndiv; j++ )
      mat[i][j] = Op<T>::sqr((w-_c1[i])+mat[i][j]) + bnd;
  }
}


template <typename T>
inline
void ISModel<T>::_sqrt
( std::vector<std::vector<T>>& mat, std::vector<std::vector<std::vector<long double>>>& slope, const std::vector<long double>& partitionSize, unsigned const& ndep )
const
{
  assert( !mat.empty() );

  // Bounds
  T bnd = _B( mat, 1 );
  double L( Op<T>::l(bnd) );
  if ( L < 0. )
    throw Exceptions( Exceptions::SQRT );

  //T bnd1(Op<T>::l(bnd),Op<T>::u(bnd));

  auto const& f = [=]( const double& x ){ return std::sqrt( x ); };
  auto const& fDerv = [=]( const double& x ){ return 0.5/std::sqrt( x ); };
  return _asym_slope( mat, slope, partitionSize, ndep, f, fDerv, DBL_MAX, false, bnd ); // concave term, max +INF
}

template <typename T>
inline
void ISModel<T>::_sqrt
( std::vector<std::vector<T>>& mat, unsigned const& ndep )
const
{
  assert( !mat.empty() );

  // Bounds
  T bnd = _B( mat, 1 );
  double L( Op<T>::l(bnd) );
  if ( L < 0. )
    throw Exceptions( Exceptions::SQRT );

  // Asymetric inclusion
  if( options.ASYREM_USE ){
    auto const& f = [=]( const double& x ){ return std::sqrt( x ); };
    return _asym( mat, ndep, f, DBL_MAX, false, bnd ); // concave term, max +INF
  }

  // Central points
  double w( 0. ), s1( 0. ), s2( 0. );
  for( unsigned int i=0; i<_nvar; i++ ){
    if( mat[i].empty() ) continue;
    _c1[i] = 0.5 * ( _L1[i] + _U1[i] );
    w += _c1[i];
    _r1[i] = 0.5 * ( _U1[i] - _L1[i] );
    s1 += _r1[i];
    s2 += std::sqrt(_r1[i]);
  }
  
  // Remainder
  double rem( s2 - std::sqrt(s1) );

  // Interval matrix coefficients
  bnd = T(-1.,1.)*(rem/double(ndep)) - (ndep-1.)/double(ndep)*std::sqrt(w);
  for( unsigned int i=0; i<_nvar; i++ ){
    if( mat[i].empty() ) continue;
    for( unsigned int j=0; j<_ndiv; j++ )
      mat[i][j] = Op<T>::sqrt((w-_c1[i])+mat[i][j]) + bnd;
  }
}


template <typename T>
inline
void ISModel<T>::_exp
( std::vector<std::vector<T>>& mat, std::vector<std::vector<std::vector<long double>>>& slope, const std::vector<long double>& partitionSize, unsigned const& ndep )
const
{
  assert( !mat.empty() );

  // Bounds
  T bnd = _B( mat, 1 );

  //T bnd1(Op<T>::l(bnd),Op<T>::u(bnd));

  auto const& f = [=]( const double& x ){ return std::exp( x ); };
  auto const& fDerv = [=]( const double& x ){ return std::exp( x ); };
  return _asym_slope( mat, slope, partitionSize, ndep, f, fDerv, -DBL_MAX, true, bnd ); // concave term, max +INF
}

template <typename T>
inline
void ISModel<T>::_exp
( std::vector<std::vector<T>>& mat, unsigned const& ndep )
const
{
  assert( !mat.empty() );

  // Bounds
  T bnd = _B( mat, 1 );

  //T bnd1(Op<T>::l(bnd),Op<T>::u(bnd));

  // Asymetric inclusion
  if( options.ASYREM_USE ){
    auto const& f = [=]( const double& x ){ return std::exp( x ); };
    return _asym( mat, ndep, f, -DBL_MAX, true, bnd ); // concave term, max +INF
  }
  
  // Central points
  double w( 0. ), s( 0. ), p( 1. );
  for( unsigned int i=0; i<_nvar; i++ ){
    if( mat[i].empty() ) continue;
    _c1[i] = std::log( 0.5 * ( std::exp(_L1[i]) + std::exp(_U1[i]) ) );
    w += _c1[i];
    _r1[i] = ( std::exp(_U1[i]) - std::exp(_L1[i]) ) / ( std::exp(_U1[i]) + std::exp(_L1[i]) );
    s += _r1[i];
    p *= _r1[i] + 1.;
  }
  
  // Remainder
  double rem( std::exp(w) * ( p - s - 1. ) );

  // Interval matrix coefficients
  bnd = T(-1.,1.)*(rem/double(ndep)) - (ndep-1.)/double(ndep)*std::exp(w);
  for( unsigned int i=0; i<_nvar; i++ ){
    if( mat[i].empty() ) continue;
    for( unsigned int j=0; j<_ndiv; j++ )
      mat[i][j] = Op<T>::exp((w-_c1[i])+mat[i][j]) + bnd;
  }
  //_intersect(mat,ndep,Op<T>::exp(bnd1));
}


template <typename T>
inline
void ISModel<T>::_log
( std::vector<std::vector<T>>& mat, std::vector<std::vector<std::vector<long double>>>& slope, const std::vector<long double>& partitionSize, unsigned const& ndep )
const
{
  assert( !mat.empty() );

  // Bounds
  T bnd = _B( mat, 1 );
  double L( Op<T>::l(bnd) );
  if ( L <= 0. )
    throw Exceptions( Exceptions::LOG );


  auto const& f = [=]( const double& x ){ return std::log( x ); };
  auto const& fDerv = [=]( const double& x ){ return 1.0/x; };
  return _asym_slope( mat, slope, partitionSize, ndep, f, fDerv, DBL_MAX, false, bnd ); // concave term, max +INF
}

template <typename T>
inline
void ISModel<T>::_log
( std::vector<std::vector<T>>& mat, unsigned const& ndep )
const
{
  assert( !mat.empty() );

  // Bounds
  T bnd = _B( mat, 1 );
  double L( Op<T>::l(bnd) );
  if ( L <= 0. )
    throw Exceptions( Exceptions::LOG );

  // Asymetric inclusion
  if( options.ASYREM_USE ){
    auto const& f = [=]( const double& x ){ return std::log( x ); };
    return _asym( mat, ndep, f, DBL_MAX, false, bnd ); // concave term, max +INF
  }

  // Central points
  double w( 0. ), s( 0. );
  for( unsigned int i=0; i<_nvar; i++ ){
    if( mat[i].empty() ) continue;
    _c1[i] = 0.5 * ( _L1[i] + _U1[i] );
    w += _c1[i];
    _r1[i] = 0.5 * ( _U1[i] - _L1[i] );
    s += _r1[i];
  }
  
  // Remainder
  double p( 1. );
  for( unsigned int i=0; i<_nvar; i++ ){
    if( mat[i].empty() ) continue;
    p *= w + _r1[i];
  }
  double wn( std::pow(w,ndep-1) );
  double rem( -std::log( 1.- ( p - wn*(w+s) ) / ( wn * L ) ) );

  // Interval matrix coefficients
  bnd = T(-1.,1.)*(rem/double(ndep)) - (ndep-1.)/double(ndep)*std::log(w);
  for( unsigned int i=0; i<_nvar; i++ ){
    if( mat[i].empty() ) continue;
    for( unsigned int j=0; j<_ndiv; j++ )
      mat[i][j] = Op<T>::log((w-_c1[i])+mat[i][j]) + bnd;
  }
}


template <typename T>
inline
void ISModel<T>::_xlog
( std::vector<std::vector<T>>& mat, std::vector<std::vector<std::vector<long double>>>& slope, const std::vector<long double>& partitionSize, unsigned const& ndep )
const
{
  assert( !mat.empty() );

  // Bounds
  T bnd = _B( mat, 1 );
  double L( Op<T>::l(bnd) );
  if ( L <= 0. )
    throw Exceptions( Exceptions::LOG );


  auto const& f = [=]( const double& x ){  return x * std::log( x ); };
  auto const& fDerv = [=]( const double& x ){ return std::log( x ) + 1.0; };
  std::cout << "xlog slope" << std::endl;
  return _asym_slope( mat, slope, partitionSize, ndep, f, fDerv, std::exp(-1.), true, bnd ); // convex term, min 1/e
}


template <typename T>
inline
void ISModel<T>::_xlog
( std::vector<std::vector<T>>& mat, unsigned const& ndep )
const
{
  assert( !mat.empty() );

  // Bounds
  T bnd = _B( mat, 1 );
  double L( Op<T>::l(bnd) );
  if ( L <= 0. )
    throw Exceptions( Exceptions::LOG );

  // Asymetric inclusion
  auto const& f = [=]( const double& x ){ return x * std::log( x ); };
  return _asym( mat, ndep, f, std::exp(-1.), true, bnd ); // convex term, min 1/e
}

template <typename T>
inline
void ISModel<T>::_sin
( std::vector<std::vector<T>>& mat, unsigned const& ndep )
const
{
  assert( !mat.empty() );

  // Central points
  T bnd = _B( mat, 1 );
  double w( 0. ), r( 0. ), s( 0. ), p( 1. );
  for( unsigned int i=0; i<_nvar; i++ ){
    if( mat[i].empty() ) continue;
    _c1[i] = 0.5 * ( _L1[i] + _U1[i] );
    w += _c1[i];
    r = 0.25 * ( _U1[i] - _L1[i] );
    if( r <= 0.5*PI )
      _r1[i] = 2. * std::sin( r );
    else
      _r1[i] = 2.;
    s += _r1[i];
    p *= _r1[i] + 1.;
  }

  // Remainder
  double rem( (std::fabs(std::sin(w))+std::fabs(std::cos(w))) * ( p - s - 1. ) );

  // Interval matrix coefficients
  bnd = T(-1.,1.)*(rem/double(ndep)) - (ndep-1.)/double(ndep)*std::sin(w);
  for( unsigned int i=0; i<_nvar; i++ ){
    if( mat[i].empty() ) continue;
    for( unsigned int j=0; j<_ndiv; j++ )
      mat[i][j] = Op<T>::sin((w-_c1[i])+mat[i][j]) + bnd;
  }
}

template <typename T>
inline
void ISModel<T>::_cos
( std::vector<std::vector<T>>& mat, unsigned const& ndep )
const
{
  assert( !mat.empty() );

  // Central points
  T bnd = _B( mat, 1 );
  double w( 0. ), r( 0. ), s( 0. ), p( 1. );
  for( unsigned int i=0; i<_nvar; i++ ){
    if( mat[i].empty() ) continue;
    _c1[i] = 0.5 * ( _L1[i] + _U1[i] );
    w += _c1[i];
    r = 0.25 * ( _U1[i] - _L1[i] );
    if( r <= 0.5*PI )
      _r1[i] = 2. * std::sin( r );
    else
      _r1[i] = 2.;
    s += _r1[i];
    p *= _r1[i] + 1.;
  }

  // Remainder
  double rem( (std::fabs(std::sin(w))+std::fabs(std::cos(w))) * ( p - s - 1. ) );

  // Interval matrix coefficients
  bnd = T(-1.,1.)*(rem/double(ndep)) - (ndep-1.)/double(ndep)*std::cos(w);
  for( unsigned int i=0; i<_nvar; i++ ){
    if( mat[i].empty() ) continue;
    for( unsigned int j=0; j<_ndiv; j++ )
      mat[i][j] = Op<T>::cos((w-_c1[i])+mat[i][j]) + bnd;
  }
}



template <typename T>
inline
void ISModel<T>::_cos_slope
( std::vector<std::vector<T>>& mat, std::vector<std::vector<std::vector<long double>>>& slope, const std::vector<long double>& partitionSize, unsigned const& ndep )
const
{
  assert( !mat.empty() );

  // this function is implemented only for the cases when ndep = 1
  // T bnd = _B( mat, 1 );
  // long double lambda = Op<T>::l(bnd);
  // long double mu = Op<T>::u(bnd);
  for( unsigned int i=0; i<_nvar; i++ ){
    if( mat[i].empty() ) continue;
    for( unsigned int j=0; j<_ndiv; j++ ){
      long double zL = Op<T>::l(mat[i][j]), zU = Op<T>::u(mat[i][j]);
      const int k = std::ceil(-(1.+zL/PI)/2.); // -pi <= xL+2*k*pi < pi
      const long double l = zL+2.*PI*k, u = zU+2.*PI*k;
      if( l <= 0 ){
        if( u <= 0 ){
          long double delta_zL = std::fabs(slope[i][j][0]*partitionSize[i]);
          long double slope0_candidate_1 = -std::sin(l);
          long double slope0_candidate_2 = slope0_candidate_1;
          if (delta_zL > 1e2*MC__ISM_COMPUTATION_TOL)
            slope0_candidate_2 = ( std::cos(l + delta_zL) - std::cos(l) )/delta_zL;
          else
            slope0_candidate_2 = std::min(-std::sin(l),-std::sin(l+delta_zL));
          long double slope0_to_be_multiplied = std::min(slope0_candidate_1,slope0_candidate_2);      // positive slope
          slope[i][j][0] = slope[i][j][0]*slope0_to_be_multiplied;

          long double delta_zU = std::fabs(slope[i][j][1]*partitionSize[i]);
          long double slope1_candidate_1 = -std::sin(u);
          long double slope1_candidate_2 = slope1_candidate_1;
          if (delta_zU > 1e2*MC__ISM_COMPUTATION_TOL)
            slope1_candidate_2 = ( std::cos(u) - std::cos(u - delta_zU) )/delta_zU;
          else 
            slope1_candidate_2 = std::min(-std::sin(u),-std::sin(u-delta_zU));
          long double slope1_to_be_multiplied = std::min(slope1_candidate_1,slope1_candidate_2);      // positive slope
          slope[i][j][1] = slope[i][j][1]*slope1_to_be_multiplied;

          mat[i][j] = T(std::cos(l), std::cos(u));
          //return Interval( std::cos(l), std::cos(u) );
        }   
        else if( u >= PI ){
          // note that there are better ways to do so, c.f. envelope
          mat[i][j] = T(-1., 1.);
          slope[i][j][0] = 0.;
          slope[i][j][1] = 0.;
          //return Interval( -1., 1. );
        } 
        else { 
          const long double delta_zL = std::fabs(slope[i][j][0]*partitionSize[i]);
          long double slope0_left_candidate_1 = -std::sin(l);
          // // Since we have known that l in [-pi,0] u in [0,pi] 
          long double slope0_left_candidate_2 = std::max((long double)0.,-std::sin(l+delta_zL));
          if (delta_zL > 1e2*MC__ISM_COMPUTATION_TOL){
            if (l + delta_zL < 0.)
              slope0_left_candidate_2 = std::max(slope0_left_candidate_2,( std::cos(l + delta_zL) - std::cos(l) )/delta_zL);
            else
              slope0_left_candidate_2 = std::max(slope0_left_candidate_2,( 1.0 - std::cos(l) )/delta_zL);
          }
          long double slope0_left_to_be_set = slope[i][j][0]*std::min(slope0_left_candidate_1,slope0_left_candidate_2);     // positive slope
          long double lowerbound_left = std::cos(l);
          // std::cout << "slope0_left_candidate_1 = " << slope0_left_candidate_1 << std::endl;
          // std::cout << "slope0_left_candidate_2 = " << slope0_left_candidate_2 << std::endl;
          // std::cout << "slope0_left_to_be_set = " << slope0_left_to_be_set << std::endl;

          //long double slope1_right = 0.;
          //long double upperbound_right = 1.;

          long double delta_zU = std::fabs(slope[i][j][1]*partitionSize[i]);
          long double slope0_right_candidate_1 = -std::sin(u);
          long double slope0_right_candidate_2 = std::min((long double)0.,-std::sin(u-delta_zU));
          if (delta_zU > 1e2*MC__ISM_COMPUTATION_TOL){
            if (u - delta_zU > 0)
              slope0_right_candidate_2 = std::min(slope0_right_candidate_2,(std::cos(u) - std::cos(u - delta_zU))/delta_zU);
            else
              slope0_right_candidate_2 = std::min(slope0_right_candidate_2,(std::cos(u) - 1.)/delta_zU); 
          }       
          long double slope0_right_to_be_set = slope[i][j][1]*std::max(slope0_right_candidate_1,slope0_right_candidate_2);  // negative slope
          long double lowerbound_right = std::cos(u);
          // std::cout << "slope0_right_candidate_1 = " << slope0_right_candidate_1 << std::endl;
          // std::cout << "slope0_right_candidate_2 = " << slope0_right_candidate_2 << std::endl;
          // std::cout << "slope0_right_to_be_set = " << slope0_right_to_be_set << std::endl;


          //long double slope1_left = 0.;
          //long double upperbound_left = 1.;
          slope[i][j][1] = 0.;

          long double lowerbound_to_be_set = -1.0 + lowerbound_left + lowerbound_right;
          slope[i][j][0] = slope0_left_to_be_set + slope0_right_to_be_set;
          // if(slope[i][j][0]>0)
          //   lowerbound_to_be_set += slope0_right_to_be_set*partitionSize[i];
          // else 
          //   lowerbound_to_be_set += slope0_left_to_be_set*partitionSize[i];
            
          if((slope0_left_to_be_set < 0.) != (slope0_right_to_be_set < 0.)){   
            lowerbound_to_be_set += std::min(std::fabs(slope0_left_to_be_set),std::fabs(slope0_right_to_be_set))*partitionSize[i];
          }
          if ( 1.0 - lowerbound_to_be_set  < MC__ISM_COMPUTATION_TOL)
            mat[i][j] = T(std::min((long double)1.,lowerbound_to_be_set),1.);
          else if (lowerbound_to_be_set < std::min(std::cos(l), std::cos(u))){
            mat[i][j] = T(std::min(std::cos(l), std::cos(u)),1.);
            slope[i][j][0] = 0.;
          }
          else 
            mat[i][j] = T(lowerbound_to_be_set, 1.);

        } 
        
      }
      else if( u <= PI ){
        long double delta_zL = std::fabs(slope[i][j][0]*partitionSize[i]);
        long double slope1_candidate_1 = -std::sin(l);
        long double slope1_candidate_2 = slope1_candidate_1;
        if (delta_zL > 1e2*MC__ISM_COMPUTATION_TOL)
          slope1_candidate_2 = ( std::cos(l + delta_zL) - std::cos(l) )/delta_zL;
        else 
          slope1_candidate_2 = std::max(-std::sin(l),-std::sin(l+delta_zL));
        long double slope1_to_be_set   = slope[i][j][0]*std::max(slope1_candidate_1,slope1_candidate_2);     // negative slope, still anchoring at the same point w.r.t. to x

        long double delta_zU = std::fabs(slope[i][j][1]*partitionSize[i]);
        long double slope0_candidate_1 = -std::sin(u);
        long double slope0_candidate_2 = slope0_candidate_1;
        if (delta_zU > 1e2*MC__ISM_COMPUTATION_TOL)
          slope0_candidate_2 = ( std::cos(u) - std::cos(u - delta_zU) )/delta_zU;
        else 
          slope0_candidate_2 = std::max(-std::sin(u),-std::sin(u-delta_zU));
        long double slope0_to_be_set  = slope[i][j][1]*std::max(slope0_candidate_1,slope0_candidate_2);     // negative slope, still anchoring at the same point w.r.t. to x

        slope[i][j][0] = slope0_to_be_set;
        slope[i][j][1] = slope1_to_be_set;
        
        mat[i][j] = T(std::cos(u), std::cos(l));
        //return Interval( std::cos(u), std::cos(l) );
      }    
      else if( u >= 2.*PI ){
        // note that there are better ways to do so c.f. envelope
        mat[i][j] = T(-1., 1.);
        slope[i][j][0] = 0.;
        slope[i][j][1] = 0.;
        //return Interval( -1., 1. );
      }
      else{ 
          const long double delta_zL = std::fabs(slope[i][j][0]*partitionSize[i]);
          long double slope1_left_candidate_1 = -std::sin(l);
          // // Since we have known that l in [0, pi] u in [pi,2pi] 
          long double slope1_left_candidate_2 = std::min((long double)0.,-std::sin(l+delta_zL));
          if (delta_zL > 1e2*MC__ISM_COMPUTATION_TOL){
            if (l + delta_zL < PI)
              slope1_left_candidate_2 = std::min(slope1_left_candidate_2,( std::cos(l + delta_zL) - std::cos(l) )/delta_zL);
            else
              slope1_left_candidate_2 = std::min(slope1_left_candidate_2,( -1.0 - std::cos(l))/delta_zL);
          }
          long double slope1_left_to_be_set = slope[i][j][0]*std::max(slope1_left_candidate_1,slope1_left_candidate_2);     // negative slope
          long double upperbound_left = std::cos(l);

          //long double slope0_right = 0.;
          //long double lowerbound_right = -1.;
          
          // std::cout << "slope1_left_candidate_1 = " << slope1_left_candidate_1 << std::endl;
          // std::cout << "slope1_left_candidate_2 = " << slope1_left_candidate_2 << std::endl;
          // std::cout << "slope1_left_to_be_set = " << slope1_left_to_be_set << std::endl;

          long double delta_zU = std::fabs(slope[i][j][1]*partitionSize[i]);
          long double slope1_right_candidate_1 = -std::sin(u);
          long double slope1_right_candidate_2 = std::max((long double)0.,-std::sin(u-delta_zU));
          if (delta_zU > 1e2*MC__ISM_COMPUTATION_TOL){
            if (u - delta_zU > PI)
              slope1_right_candidate_2 = std::max(slope1_right_candidate_2,(std::cos(u) - std::cos(u - delta_zU))/delta_zU);
            else
              slope1_right_candidate_2 = std::max(slope1_right_candidate_2,(std::cos(u) + 1.)/delta_zU); 
          }       
          long double slope1_right_to_be_set = slope[i][j][1]*std::min(slope1_right_candidate_1,slope1_right_candidate_2);  // postive slope
          long double upperbound_right = std::cos(u);

          //long double slope0_left = 0.;
          //long double lowerbound_left = -1.;
          slope[i][j][0] = 0.;
          // std::cout << "slope1_right_candidate_1 = " << slope1_right_candidate_1 << std::endl;
          // std::cout << "slope1_right_candidate_2 = " << slope1_right_candidate_2 << std::endl;
          // std::cout << "slope1_right_to_be_set = " << slope1_right_to_be_set << std::endl;


          long double upperbound_to_be_set = 1.0 + upperbound_left + upperbound_right;
          slope[i][j][1] = slope1_left_to_be_set + slope1_right_to_be_set;
          // if(slope[i][j][0]>0)
          //   lowerbound_to_be_set += slope0_right_to_be_set*partitionSize[i];
          // else 
          //   lowerbound_to_be_set += slope0_left_to_be_set*partitionSize[i];
            
          if((slope1_left_to_be_set < 0.) != (slope1_right_to_be_set < 0.)){   
            upperbound_to_be_set -= std::min(std::fabs(slope1_left_to_be_set),std::fabs(slope1_right_to_be_set))*partitionSize[i];
          }
          if (upperbound_to_be_set + 1.0 < MC__ISM_COMPUTATION_TOL)
            mat[i][j] = T(-1,std::max((long double)-1.,upperbound_to_be_set));
          else if (upperbound_to_be_set > std::max(std::cos(l), std::cos(u))){
            mat[i][j] = T(-1.,std::max(std::cos(l), std::cos(u)));
            slope[i][j][1] = 0.;
          }
          else 
            mat[i][j] = T(-1.,upperbound_to_be_set);
        }
      //slope[i][j][1] = 0.;
      //slope[i][j][0] = 0.;
    }
  }
}


template <typename T>
inline
void ISModel<T>::_prod
( std::vector<std::vector<T>> const& mat1, std::vector<std::vector<T>> const& mat2,
  std::vector<std::vector<T>>& mat3, unsigned& ndep3 )
const
{
  assert( !mat1.empty() && !mat2.empty() );
  // Bounds
  _B( mat1, 1 );
  _B( mat2, 2 );
  //T bnd1 = _B( mat1, 1 ), bnd2 = _B( mat2, 2 );

  // Central points and radii
  double C1( 0. ), R1( 0. );
  for( unsigned int i=0; i<_nvar; i++ ){
    if( mat1[i].empty() ) continue;
    _c1[i] = 0.5 * ( _L1[i] + _U1[i] );
    C1 += _c1[i];
    _r1[i] = 0.5 * ( _U1[i] - _L1[i] );
    R1 += _r1[i];
  }
#ifdef MC__ISMODEL_DEBUG_PROD
  std::cerr << "C1 = " << C1 << std::endl;
  std::cerr << "R1 = " << R1 << std::endl;
#endif
  double C2( 0. ), R2( 0. );
  for( unsigned int i=0; i<_nvar; i++ ){
    if( mat2[i].empty() ) continue;
    _c2[i] = 0.5 * ( _L2[i] + _U2[i] );
    C2 += _c2[i];
    _r2[i] = 0.5 * ( _U2[i] - _L2[i] );
    R2 += _r2[i];
  }
#ifdef MC__ISMODEL_DEBUG_PROD
  std::cerr << "C2 = " << C2 << std::endl;
  std::cerr << "R2 = " << R2 << std::endl;
#endif
  double C12( 0. ), R12( 0. );
  ndep3 = 0;
  for( unsigned int i=0; i<_nvar; i++ ){
    if( !mat1[i].empty() || !mat2[i].empty() ) ndep3++;
    if(  mat1[i].empty() ||  mat2[i].empty() ) continue;
    C12 +=  _c1[i] * _c2[i];
    R12 +=  _r1[i] * _r2[i];
  }
  double w( C1*C2 - C12 );
#ifdef MC__ISMODEL_DEBUG_PROD
  std::cerr << "C12 = " << C12 << std::endl;
  std::cerr << "R12 = " << R12 << std::endl;
  std::cerr << "ndep3 = " << ndep3 << std::endl;
#endif
  
  // Remainder
  double rem( R1*R2 - R12 );
  // Interval matrix coefficients
  double _asybnd (-w/double(ndep3));
  T bnd = T(-1.,1.)*(rem/double(ndep3)) + _asybnd;
#ifdef MC__ISMODEL_DEBUG_PROD
  std::cerr << "bnd = " << bnd << std::endl;
#endif
  for( unsigned int i=0; i<_nvar; i++ ){
    if( mat1[i].empty() && mat2[i].empty() ) continue;
    if( !mat1[i].empty() && !mat2[i].empty() ){
      if( mat3[i].empty() ) mat3[i].resize(_ndiv);
#ifdef MC__ISMODEL_DEBUG_PROD
      std::cerr << "1 & 2" << std::endl;
#endif
      for( unsigned int j=0; j<_ndiv; j++ ){
        // mat3[i][j] = (mat1[i][j]+(C1-_c1[i]))*(mat2[i][j]+(C2-_c2[i])) - (C1-_c1[i])*(C2-_c2[i])+bnd;
        T _bnd(0.);
        _bnd = bnd; // a new internal variable which does not modify the value of the variable bnd
        if(options.ASYREM_USE){
          T _temp(0.);
          T _asyremU1 (_r1[i]==0? _temp: Op<T>::sqr(mat1[i][j]-_c1[i])/(_r1[i]*2.)*(R2-_r2[i])*T(-1.,1.) );
          T _asyremU2 (_r2[i]==0? _temp: Op<T>::sqr(mat2[i][j]-_c2[i])/(_r2[i]*2.)*(R1-_r1[i])*T(-1.,1.) );
          _bnd = _asyremU1+_asyremU2+_asybnd;
        }
        mat3[i][j] = (mat1[i][j]+(C1-_c1[i]))*(mat2[i][j]+(C2-_c2[i])) - (C1-_c1[i])*(C2-_c2[i]) + _bnd;
      }
    }
    else if( !mat1[i].empty() ){
      if( mat3[i].empty() ) mat3[i].resize(_ndiv);
#ifdef MC__ISMODEL_DEBUG_PROD
      std::cerr << "1 ONLY" << std::endl;
#endif
      for( unsigned int j=0; j<_ndiv; j++ ){
      // mat3[i][j] = mat1[i][j]*C2 + bnd;
      // mat3[i][j] = mat1[i][j]+0.;
      // mat3[i][j] = mat1[i][j]*C2; 
      // Here is a strange bug: if we execute mat3[i][j] = mat1[i][j]*C2, then the value of mat1[i][j] is changed
      // This might be caused by the convensional use of storing coefficients, in ISVar<T>& ISVar<T>::operator*= ( ISVar<T> const& var )
      // Therefore, we add the value after computing _asyremU1  
        T _bnd(0.);
        _bnd = bnd; // a new internal variable which does not modify the value of the variable bnd
        if(options.ASYREM_USE){
          T _temp(0.);
          T _asyremU1 (_r1[i]==0? _temp: Op<T>::sqr(mat1[i][j]-_c1[i])/(_r1[i]*2.)*(R2)*T(-1.,1.) );
          _bnd = _asyremU1+_asybnd;
        }
        mat3[i][j] = mat1[i][j]*C2 + _bnd;
      }
      //std::cout<<"i"<<i<<","<<R2<<","<<_c1[i]<<","<<_r1[i]<<std::endl;       
    }
    else if( !mat2[i].empty() ){
      if( mat3[i].empty() ) mat3[i].resize(_ndiv);
      for( unsigned int j=0; j<_ndiv; j++ ){
        //mat3[i][j] = mat2[i][j]*C1 + bnd;
        T _bnd(0.);
        _bnd = bnd; // a new internal variable which does not modify the value of the variable bnd
        if(options.ASYREM_USE){
          T _temp(0.);
          T _asyremU2 (_r2[i]==0? _temp: Op<T>::sqr(mat2[i][j]-_c2[i])/(_r2[i]*2.)*(R1)*T(-1.,1.) );
          _bnd = _asyremU2+_asybnd;
        }
        mat3[i][j] = mat2[i][j]*C1 + _bnd; 
      }
      //std::cout<<"i"<<i<<","<<R1<<","<<_c2[i]<<","<<_r2[i]<<std::endl;
    }
  }
}





template <typename T>
inline
void ISModel<T>::_intersect
( std::vector<std::vector<T>>& mat, std::vector<std::vector<std::vector<long double>>>& slope, const std::vector<long double>& partitionSize, unsigned const& ndep, T _rangebnd )
const  
{
  assert( !mat.empty() );
  //std::cout << "Intersection for slope is processing" << std::endl;
  // Bounds
  T bnd = _B( mat, 1 );   // Note that _U1 and _L1 have been updated here 

  // std::vector<long double> rowlb,rowub;
  // rowlb.resize(_nvar);
  // rowub.resize(_nvar);  
  // for( unsigned int i=0; i<_nvar; i++ ){
  //   if( mat[i].empty() ) continue;
  //   rowlb[i] = _L1[i];
  //   rowub[i] = _U1[i];
  // }
  long double lambda = Op<T>::l(bnd);
  long double mu =  Op<T>::u(bnd);
  long double z_u = Op<T>::l(_rangebnd);
  long double z_o = Op<T>::u(_rangebnd);

  
  

  if(z_o < mu - MC__ISM_COMPUTATION_TOL){

    // Update all maximum of all components of the input overestimator, stored in _L2[i]
    long double sigma_o = 0.;
    for( unsigned int i=0; i<_nvar; i++ ){
      if( mat[i].empty() ) continue;
      long double _tmp_row = ((long double)(Op<T>::l( mat[i][0]) + std::fabs(slope[i][0][0]*partitionSize[i])));
       for( unsigned int j=1; j<_ndiv; j++ ){
         _tmp_row = std::max(_tmp_row, ((long double)(Op<T>::l( mat[i][j]) + std::fabs(slope[i][j][0]*partitionSize[i]))));
       }
       _L2[i]= _tmp_row;
       sigma_o += _L2[i];
    }

    if (z_o < sigma_o - MC__ISM_COMPUTATION_TOL)  // when there is something wrong, as the input _rangebnd must be valid over the whole domain
      throw typename ISModel<T>::Exceptions( ISModel<T>::Exceptions::UNDEF ); 
    else if (z_o < sigma_o + MC__ISM_COMPUTATION_TOL) // when z_o does not make improvement  
      ;
    else {
      // Step 2: anchor points
      long double w_o (0.);
      w_o = std::min(mu,z_o);

      // Step 3: update the upper bounds
      for( unsigned int i=0; i<_nvar; i++ ){
        if( mat[i].empty() ) continue;
        long double in_first_term( -_U1[i] + mu );
        long double second_term( w_o + ( _L2[i] - _U1[i] )/( sigma_o - mu ) * ( sigma_o - w_o) - _L2[i] );
        //temp += second_term;
        for( unsigned int j=0; j<_ndiv; j++ ){
          long double _upbd(0.);
          long double _tmp_u = Op<T>::u( mat[i][j]);
          _upbd = std::min(in_first_term + _tmp_u, w_o) - second_term;

          if (in_first_term >= w_o - _tmp_u) slope[i][j][1] = 0.;      

          //std::cout<<"i"<<i<<"j"<<j<<std::endl;
          //std::cout<<"prev"<<Op<T>::u( mat[i][j])<<std::endl;
          //std::cout<<"uppe"<<_upbd <<std::endl;
          //T _update(_lb, _upbd);
          long double _tmp_l = Op<T>::l( mat[i][j]);
          mat[i][j] = T(_tmp_l,_upbd);
#ifdef MC__USE_FILIB           
          if(mat[i][j].isEmpty()){
            if(std::fabs(_tmp_l - _upbd) <= 1e3*MC__ISM_COMPUTATION_TOL){
              mat[i][j] = T( _upbd,_tmp_l );
              slope[i][j][1] = 0.;
            }
                
            else
               std::cout <<std::setprecision(18)<< "over " << _upbd << "    " << _tmp_l
                                                << "    " << i << "," << j << std::endl;
          }           
#endif                    
          //std::cout<<"late"<<Op<T>::l(mat[i][j])<<std::endl;
        }
      }
    }
  } 

  if(z_u  >  lambda + MC__ISM_COMPUTATION_TOL){
    // Update all minimum of all components of the input overestimator, stored in _U2[i]
    long double sigma_u = 0.;
    for( unsigned int i=0; i<_nvar; i++ ){
      if( mat[i].empty() ) continue;
      long double _tmp_row = ((long double)(Op<T>::u( mat[i][0]) - std::fabs(slope[i][0][1]*partitionSize[i])));
       for( unsigned int j=1; j<_ndiv; j++ ){
         _tmp_row = std::min(_tmp_row, ((long double)(Op<T>::u( mat[i][j]) - std::fabs(slope[i][j][1]*partitionSize[i]) )));
       }
       _U2[i]= _tmp_row;
       sigma_u += _U2[i];
    }
    
    if (z_u > sigma_u + MC__ISM_COMPUTATION_TOL)  // when there is something wrong, as the input _rangebnd must be valid over the whole domain
      throw typename ISModel<T>::Exceptions( ISModel<T>::Exceptions::UNDEF ); 
    else if (z_u > sigma_u - MC__ISM_COMPUTATION_TOL) // when z_u does not make improvement  
      ;
    else {
      // Step 2: anchor points
      long double w_u (0.);
      w_u = std::max(lambda,z_u);
      
      // Step 3: update the lower bounds
      for( unsigned int i=0; i<_nvar; i++ ){
        if( mat[i].empty() ) continue;
        long double in_first_term( -_L1[i] + lambda );
        long double second_term( w_u + ( _U2[i] - _L1[i] )/( sigma_u - lambda ) * ( sigma_u - w_u) - _U2[i] );
        for( unsigned int j=0; j<_ndiv; j++ ){
          long double _tmp_l = Op<T>::l( mat[i][j]);

          if (in_first_term <= w_u - _tmp_l) slope[i][j][0] = 0.;

          long double _lowerbound(std::max(in_first_term + _tmp_l, w_u) - second_term);
          long double _tmp_u = Op<T>::u( mat[i][j]);
          mat[i][j] = T(_lowerbound,_tmp_u );
#ifdef MC__USE_FILIB           
          if(mat[i][j].isEmpty()){
            if(std::fabs(_lowerbound - _tmp_u) <= 1e3*MC__ISM_COMPUTATION_TOL){
              mat[i][j] = T( _tmp_u, _lowerbound );
              slope[i][j][0] = 0.;
            }
            else
               std::cout <<std::setprecision(18)<< "under  " << _tmp_u << "    " << _lowerbound 
                                                << "    " << i << "," << j << std::endl;
          }          
#endif              
        }
      }
    }
  }

}



template <typename T>
inline
void ISModel<T>::_intersect
( std::vector<std::vector<T>>& mat, unsigned const& ndep, T _rangebnd)
const
{
  assert( !mat.empty() );

  // Bounds
  T bnd = _B( mat, 1 );   // Note that _U1 and _L1 have been updated here 

  // std::vector<long double> rowlb,rowub;
  // rowlb.resize(_nvar);
  // rowub.resize(_nvar);  
  // for( unsigned int i=0; i<_nvar; i++ ){
  //   if( mat[i].empty() ) continue;
  //   rowlb[i] = _L1[i];
  //   rowub[i] = _U1[i];
  // }
  long double lambda = Op<T>::l(bnd);
  long double mu =  Op<T>::u(bnd);
  long double z_u = Op<T>::l(_rangebnd);
  long double z_o = Op<T>::u(_rangebnd);

  
  

  if(z_o < mu - MC__ISM_COMPUTATION_TOL){

    // Update all maximum of all components of the input overestimator, stored in _L2[i]
    long double sigma_o = 0.;
    for( unsigned int i=0; i<_nvar; i++ ){
      if( mat[i].empty() ) continue;
      long double _tmp_row = ((long double)(Op<T>::l( mat[i][0])));
       for( unsigned int j=1; j<_ndiv; j++ ){
         _tmp_row = std::max(_tmp_row, ((long double)(Op<T>::l( mat[i][j]))));
       }
       _L2[i]= _tmp_row;
       sigma_o += _L2[i];
    }

    if (z_o < sigma_o - MC__ISM_COMPUTATION_TOL)  // when there is something wrong, as the input _rangebnd must be valid over the whole domain
      throw typename ISModel<T>::Exceptions( ISModel<T>::Exceptions::UNDEF ); 
    else if (z_o < sigma_o + MC__ISM_COMPUTATION_TOL) // when z_o does not make improvement  
      ;
    else {
      // Step 2: anchor points
      long double w_o (0.);
      w_o = std::min(mu,z_o);

      // Step 3: update the upper bounds
      for( unsigned int i=0; i<_nvar; i++ ){
        if( mat[i].empty() ) continue;
        long double in_first_term( -_U1[i] + mu );
        long double second_term( w_o + ( _L2[i] - _U1[i] )/( sigma_o - mu ) * ( sigma_o - w_o) - _L2[i] );
        //temp += second_term;
        for( unsigned int j=0; j<_ndiv; j++ ){
          long double _upbd(0.);
          long double _tmp_u = Op<T>::u( mat[i][j]);
          _upbd = std::min(in_first_term + _tmp_u, w_o) - second_term;
          //std::cout<<"i"<<i<<"j"<<j<<std::endl;
          //std::cout<<"prev"<<Op<T>::u( mat[i][j])<<std::endl;
          //std::cout<<"uppe"<<_upbd <<std::endl;
          //T _update(_lb, _upbd);
          long double _tmp_l = Op<T>::l( mat[i][j]);
          mat[i][j] = T(_tmp_l,_upbd);
#ifdef MC__USE_FILIB           
          if(mat[i][j].isEmpty()){
            if(std::fabs(_tmp_l - _upbd) <= 1e3*MC__ISM_COMPUTATION_TOL)
                mat[i][j] = T( _upbd,_tmp_l );
            else
               std::cout <<std::setprecision(18)<< "over " << _upbd << "    " << _tmp_l
                                                << "    " << i << "," << j << std::endl;
          }          
#endif                     
          //std::cout<<"late"<<Op<T>::l(mat[i][j])<<std::endl;
        }
      }
    }
  } 

  if(z_u  >  lambda + MC__ISM_COMPUTATION_TOL){
    // Update all minimum of all components of the input overestimator, stored in _U2[i]
    long double sigma_u = 0.;
    for( unsigned int i=0; i<_nvar; i++ ){
      if( mat[i].empty() ) continue;
      long double _tmp_row = ((long double)(Op<T>::u( mat[i][0])));
       for( unsigned int j=1; j<_ndiv; j++ ){
         _tmp_row = std::min(_tmp_row, ((long double)(Op<T>::u( mat[i][j]))));
       }
       _U2[i]= _tmp_row;
       sigma_u += _U2[i];
    }
    
    if (z_u > sigma_u + MC__ISM_COMPUTATION_TOL)  // when there is something wrong, as the input _rangebnd must be valid over the whole domain
      throw typename ISModel<T>::Exceptions( ISModel<T>::Exceptions::UNDEF ); 
    else if (z_u > sigma_u - MC__ISM_COMPUTATION_TOL) // when z_u does not make improvement  
      ;
    else {
      // Step 2: anchor points
      long double w_u (0.);
      w_u = std::max(lambda,z_u);
      
      // Step 3: update the lower bounds
      for( unsigned int i=0; i<_nvar; i++ ){
        if( mat[i].empty() ) continue;
        long double in_first_term( -_L1[i] + lambda );
        long double second_term( w_u + ( _U2[i] - _L1[i] )/( sigma_u - lambda ) * ( sigma_u - w_u) - _U2[i] );
        for( unsigned int j=0; j<_ndiv; j++ ){
          long double _tmp_l = Op<T>::l( mat[i][j]);
          long double _lowerbound(std::max(in_first_term + _tmp_l, w_u) - second_term);
          long double _tmp_u = Op<T>::u( mat[i][j]);
          mat[i][j] = T(_lowerbound,_tmp_u );
#ifdef MC__USE_FILIB 
          if(mat[i][j].isEmpty()){
            if(std::fabs(_lowerbound - _tmp_u) <= 1e3*MC__ISM_COMPUTATION_TOL)
                mat[i][j] = T( _tmp_u, _lowerbound );
            else
               std::cout <<std::setprecision(18)<< "under  " << _tmp_u << "    " << _lowerbound 
                                                << "    " << i << "," << j << std::endl;
          }          
#endif              
        }
      }
    }
  }
}


//! @brief C++ class for interval superposition models of factorable functions: variable
////////////////////////////////////////////////////////////////////////
//! mc::ISVar is a C++ class for propagation of interval superposition
//! models (ISM) through factorable functions. The template parameter
//! corresponds to the type used to propagate the interval coefficients.
////////////////////////////////////////////////////////////////////////
template <typename T> 
class ISVar
{
  template <typename U> friend class ISModel;

  template <typename U> friend std::ostream& operator<<
    ( std::ostream &, ISVar<U> const& );

  template <typename U> friend ISVar<U> operator+
    ( ISVar<U> const& );
  template <typename U> friend ISVar<U> operator+
    ( ISVar<U> && );

  template <typename U> friend ISVar<U> operator+
    ( ISVar<U> const&, ISVar<U> const& );
  template <typename U> friend ISVar<U> operator+
    ( ISVar<U> const&, ISVar<U> && );
  template <typename U> friend ISVar<U> operator+
    ( ISVar<U> &&, ISVar<U> const& );    
  template <typename U> friend ISVar<U> operator+
    ( ISVar<U> &&, ISVar<U> && );    

  template <typename U> friend ISVar<U> operator+
    ( double const&, ISVar<U> const& );
  template <typename U> friend ISVar<U> operator+
    ( double const&, ISVar<U> && );

  template <typename U> friend ISVar<U> operator+
    ( ISVar<U> const&, double const& );
  template <typename U> friend ISVar<U> operator+
    ( ISVar<U> &&, double const& );


  template <typename U> friend ISVar<U> operator-
    ( ISVar<U> const& );
  template <typename U> friend ISVar<U> operator-
    ( ISVar<U> && );

  template <typename U> friend ISVar<U> operator-
    ( ISVar<U> const&, ISVar<U> const& );
  template <typename U> friend ISVar<U> operator-
    ( ISVar<U> const&, ISVar<U> && );
  template <typename U> friend ISVar<U> operator-
    ( ISVar<U> &&, ISVar<U> const& );
  template <typename U> friend ISVar<U> operator-
    ( ISVar<U> &&, ISVar<U> && );

  template <typename U> friend ISVar<U> operator-
    ( double const&, ISVar<U> const& );
  template <typename U> friend ISVar<U> operator-
    ( double const&, ISVar<U> && );


  template <typename U> friend ISVar<U> operator-
    ( ISVar<U> const&, double const& );
  template <typename U> friend ISVar<U> operator-
    ( ISVar<U> &&, double const& );



  template <typename U> friend ISVar<U> operator*
    ( ISVar<U> const&, ISVar<U> const& );
  // template <typename U> friend ISVar<U> operator*
  //   ( ISVar<U> const&, ISVar<U> && );
  // template <typename U> friend ISVar<U> operator*
  //   ( ISVar<U> &&, ISVar<U> const& );
  // template <typename U> friend ISVar<U> operator*
  //   ( ISVar<U> &&, ISVar<U> && );


  template <typename U> friend ISVar<U> operator*
    ( double const&, ISVar<U> const& );
  template <typename U> friend ISVar<U> operator*
    ( double const&, ISVar<U> && );


  template <typename U> friend ISVar<U> operator*
    ( ISVar<U> const&, double const& );
  template <typename U> friend ISVar<U> operator*
    ( ISVar<U> &&, double const& );


  template <typename U> friend ISVar<U> operator/
    ( ISVar<U> const&, ISVar<U> const& );
  // template <typename U> friend ISVar<U> operator/
  //   ( ISVar<U> const&, ISVar<U> && );
  // template <typename U> friend ISVar<U> operator/
  //   ( ISVar<U> &&, ISVar<U> const& );
  // template <typename U> friend ISVar<U> operator/
  //   ( ISVar<U> &&, ISVar<U> && );


  template <typename U> friend ISVar<U> operator/
    ( double const&, ISVar<U> const& );
  template <typename U> friend ISVar<U> operator/
    ( double const&, ISVar<U> && );

  template <typename U> friend ISVar<U> operator/
    ( ISVar<U> const&, double const& );
  template <typename U> friend ISVar<U> operator/
    ( ISVar<U> &&, double const& );

  template <typename U> friend ISVar<U> max
    ( ISVar<U> const&, ISVar<U> const& );
  template <typename U> friend ISVar<U> max
    ( ISVar<U> const&, double const& );
  template <typename U> friend ISVar<U> max
    ( ISVar<U> &&, double const& );
  template <typename U> friend ISVar<U> min
    ( ISVar<U> const&, ISVar<U> const& );
  template <typename U> friend ISVar<U> min
    ( ISVar<U> const&, double const& );
  template <typename U> friend ISVar<U> min
    ( ISVar<U> &&, double const& );

  template <typename U> friend ISVar<U> inv
    ( ISVar<U> const& );
  template <typename U> friend ISVar<U> inv
    ( ISVar<U>&& );
  template <typename U> friend ISVar<U> sqr
    ( ISVar<U> const& );
  template <typename U> friend ISVar<U> sqr
    ( ISVar<U> && );
  template <typename U> friend ISVar<U> sqrt
    ( ISVar<U> const& );
  template <typename U> friend ISVar<U> sqrt
    ( ISVar<U> && );
  template <typename U> friend ISVar<U> fabs
    ( ISVar<U> const& );
  template <typename U> friend ISVar<U> fabs
    ( ISVar<U> && );  
  template <typename U> friend ISVar<U> relu
    ( ISVar<U> const& );
  template <typename U> friend ISVar<U> relu
    ( ISVar<U> && );    
  template <typename U> friend ISVar<U> exp
    ( ISVar<U> const& );
  template <typename U> friend ISVar<U> exp
    ( ISVar<U> && );
  template <typename U> friend ISVar<U> log
    ( ISVar<U> const& );
  template <typename U> friend ISVar<U> log
    ( ISVar<U> && );
  template <typename U> friend ISVar<U> xlog
    ( ISVar<U> const& );
  template <typename U> friend ISVar<U> xlog
    ( ISVar<U> && );
  template <typename U> friend ISVar<U> sin
    ( ISVar<U> const& );
  template <typename U> friend ISVar<U> sin
    ( ISVar<U> && );
  template <typename U> friend ISVar<U> cos
    ( ISVar<U> const& );
  template <typename U> friend ISVar<U> cos
    ( ISVar<U> && );
  template <typename U> friend ISVar<U> tan
    ( ISVar<U> const& );
  template <typename U> friend ISVar<U> tan
    ( ISVar<U> && );
  template <typename U> friend ISVar<U> tanh
    ( ISVar<U> const& );
  template <typename U> friend ISVar<U> tanh
    ( ISVar<U> && );
  template <typename U> friend ISVar<U> pow
    ( ISVar<U> const&, int const& n );
  template <typename U> friend ISVar<U> pow
    ( ISVar<U> &&, int const& n );
  template <typename U> friend ISVar<U> intersect
    ( ISVar<U> const& , U);
  template <typename U> friend ISVar<U> intersect
    ( ISVar<U> && , U);  
  template <typename U> friend ISVar<U> cheb
    ( ISVar<U> const&, unsigned int const& n );
  template <typename U> friend ISVar<U> cheb
    ( ISVar<U> &&, unsigned int const& n );

 private:
 
  //! @brief Pointer to interval superposition model
  ISModel<T>* _mod;
  //! @brief Number of partitions
  unsigned int _ndiv;
  //! @brief Number of variables
  unsigned int _nvar;
  //! @brief Number of dependencies
  unsigned int _ndep;
  //! @brief Interval coefficient matrix
  std::vector<std::vector<T>> _mat;
  //! @brief Constant value (in case _mat is empty)
  T _cst;
  //! @brief Variable bound
  mutable std::pair<T,bool> _bnd;

  //! @brief Slope matrix
  std::vector<std::vector<std::vector<long double>>> _slope;
  //! @brief Shadow matrix
  std::vector<std::vector<std::vector<long double>>> _shadow;
  //std::vector<std::vector<std::vector<T>>> _shadow_mat;
  //! @brief Shadow slop matrix
  std::vector<std::vector<std::vector<long double>>> _shadow_slope;  

  mutable std::vector<double> _shadow_info;

  // std::vector<std::vector<T>> _matForAgrt;  
  // std::vector<std::vector<std::vector<long double>>> _slopeForAgrt;
  // std::vector<std::vector<std::vector<long double>>> _shadowForAgrt;
  // std::vector<std::vector<std::vector<long double>>> _shadow_slopeForAgrt;     

 public:
 
  ISVar<T>& operator+=
    ( ISVar<T> const& );
  ISVar<T>& operator+=
    ( double const& );
  ISVar<T>& operator+=
    ( T const& );


  ISVar<T>& operator-=
    ( ISVar<T> const& );
  ISVar<T>& operator-=
    ( double const& );
  ISVar<T>& operator-=
    ( T const& );
 
  ISVar<T>& operator*=
    ( ISVar<T> const& );
  ISVar<T>& operator*=
    ( double const& );
  ISVar<T>& operator*=
    ( T const& );

  ISVar<T>& operator/=
    ( ISVar<T> const& );
  ISVar<T>& operator/=
    ( double const& );

  ISVar<T>& operator=
  ( double const& cst )
  {
    _mod = nullptr;
    _nvar = 0;
    _ndiv = 0;
    _ndep = 0;
    _mat.clear(); // may not resize _mat
    _cst = T(cst);
    _bnd = std::make_pair( 0., false );
    _slope.clear();  // added for allowing slope-based enhancement
    _shadow.clear(); // added for allowing shadow enhancement
    _shadow_slope.clear(); // added for allowing shadow enhancement
    _shadow_info.clear();
    _shadow_info.resize(1);
    _shadow_info[0] = 1.;    
    return *this;
  }

ISVar<T>& operator=
  ( T const& cst )
  {
    _mod = nullptr;
    _nvar = 0;
    _ndiv = 0;
    _ndep = 0;
    _mat.clear(); // may not resize _mat
    _cst = cst;
    _bnd = std::make_pair( 0., false );
    _slope.clear();  // added for allowing slope-based enhancement
    _shadow.clear(); // added for allowing shadow enhancement
    _shadow_slope.clear(); // added for allowing shadow enhancement
    _shadow_info.clear();
    _shadow_info.resize(1);
    _shadow_info[0] = 1.;    
    return *this;
  }


  ISVar<T>& operator=
  ( ISVar<T> const& var )
  {
#ifdef MC__ISMODEL_TRACE
    std::cerr << "-- ISVar<T>& operator= ( ISVar<T> const& )\n";
#endif
#ifdef TEST_MOVE
    std::cout << "Copy Assign" << std::endl;
#endif
    if( this == &var )
      return *this;
    _mod = var._mod;
    if( !_mod ) _cst = var._cst;
    _nvar = var._nvar;
    _ndiv = var._ndiv;
    _ndep = var._ndep;
    _mat = var._mat;
    _bnd = var._bnd;
    _slope = var._slope;
    _shadow = var._shadow;
    _shadow_slope = var._shadow_slope;   
    _shadow_info = var._shadow_info;    
    return *this;
  } 

  ISVar<T>& operator=
  ( ISVar<T> && var )
  {
#ifdef MC__ISMODEL_TRACE
    std::cerr << "-- ISVar<T>& operator= ( ISVar<T> && )\n";
#endif
#ifdef TEST_MOVE
    std::cout << "Move Assign" << std::endl;
#endif
    if( this == &var )
      return *this;
    _mod = var._mod;
    if( !_mod ) _cst = std::move(var._cst);
    _nvar = var._nvar;
    _ndiv = var._ndiv;
    _ndep = var._ndep;
    //_mat = var._mat;
    //std::swap( _mat, var._mat );
/*
    _mat.swap(var._mat);
    _bnd = var._bnd;
    //std::swap( _slope, var._slope );
    _slope.swap(var._slope);
    _shadow.swap(var._shadow);
    _shadow_slope.swap(var._shadow_slope);
*/
    _mat = std::move(var._mat);
    _bnd = var._bnd;
    _slope = std::move(var._slope);
    _shadow = std::move(var._shadow);
    _shadow_slope = std::move(var._shadow_slope);       
    _shadow_info = std::move(var._shadow_info);    
    return *this;
  }



  ISVar
  ( ISModel<T>* const mod )
  : _mod(mod), _ndiv(mod->_ndiv), _nvar(mod->_nvar), _ndep(0), _mat(_nvar), _bnd(0.,false), _slope(_nvar),_shadow(4), _shadow_slope(4),_shadow_info(1)
  {
    for(unsigned i = 0; i< 4; i++){
      _shadow[i].resize(_nvar);
      _shadow_slope[i].resize(_nvar);
    }
    _shadow_info[0] = 1.;
  }

 

  ISVar
  ( ISModel<T>* const mod, unsigned int ndx, T const& bnd )
  : _mod(mod), _ndiv(mod->_ndiv), _nvar(mod->_nvar), _ndep(1), _mat(_nvar), _bnd(bnd,true), _slope(_nvar), _shadow(4), _shadow_slope(4),_shadow_info(1)
  {
    if( ndx >= _nvar )
      throw typename ISModel<T>::Exceptions( ISModel<T>::Exceptions::INDEX );

    long double l = Op<T>::l(bnd);
    long double h = Op<T>::diam(bnd) / (long double)_ndiv;
    _mat[ndx].resize( _ndiv );
    for( unsigned int j=0; j<_ndiv-1; j++, l+=h )
      _mat[ndx][j] = T( l, l+h );
    _mat[ndx][_ndiv-1] = T( Op<T>::u(bnd)-h, Op<T>::u(bnd));  
    _mod->_defvar[ndx] = true;
    _mod->_bndvar[ndx] = bnd;

    if (_mod->options.SLOPE_USE){
      _mod->_psize[ndx] = h;
      _slope[ndx].resize( _ndiv );
      for( unsigned int j=0; j<_ndiv; j++){
        _slope[ndx][j].resize(2);
        _slope[ndx][j][0] = 1.;
        _slope[ndx][j][1] = 1.;
      }
    }

    if (_mod->options.SHADOW_USE){
      //for( unsigned int j=0; j<_ndiv; j++){
      for(unsigned i = 0; i< 4; i++){
        _shadow[i].resize(_nvar);
        _shadow_slope[i].resize(_nvar);
      }
      // _shadow[0].resize(_nvar);
      // _shadow[1].resize(_nvar);
      // _shadow_slope[0].resize(_nvar);
      // _shadow_slope[1].resize(_nvar);
        _shadow_info[0] = 1.;
        // initialze containers and info
      //}
    }    

  }


//   ISVar
//   ( double const& cst=0. )
//   : _mod(nullptr), _ndiv(0), _nvar(0), _ndep(0), _cst(cst), _bnd(0.,false)
//   {}

//   ISVar
//   ( ISVar<T> const& var )
//   : _mod(var._mod), _ndiv(var._ndiv), _nvar(var._nvar), _ndep(var._ndep)
//   {
// #ifdef MC__ISMODEL_TRACE
//     std::cerr << "-- ISVar( ISVar<T> const& )\n";
// #endif
//     if( this == &var ) return;
//     if( !_mod ) _cst = var._cst;
//     _mat = var._mat;
//     _bnd = var._bnd;
//   }

  ISVar
  ( double const& cst=0. )
  : _mod(nullptr), _ndiv(0), _nvar(0), _ndep(0), _cst(T(cst)), _bnd(0.,false), _slope(0), _shadow(0), _shadow_slope(0),_shadow_info(0)
  {}

  ISVar
  ( T const& cst )
  : _mod(nullptr), _ndiv(0), _nvar(0), _ndep(0), _cst(cst), _bnd(0.,false), _slope(0), _shadow(0), _shadow_slope(0),_shadow_info(0)
  {}

  ISVar
  ( ISVar<T> const& var )
  : _mod(var._mod), _ndiv(var._ndiv), _nvar(var._nvar), _ndep(var._ndep)
  {
#ifdef MC__ISMODEL_TRACE
    std::cerr << "-- ISVar( ISVar<T> const& )\n";
#endif
#ifdef TEST_MOVE
    std::cout << "Copy Constructor" << std::endl;
#endif

    if( this == &var ) return;
    if( !_mod ) _cst = var._cst;
    _mat = var._mat;
    _bnd = var._bnd;
    _slope = var._slope;
    _shadow = var._shadow;
    _shadow_slope = var._shadow_slope;
    _shadow_info = var._shadow_info;

  }

  ISVar
  ( ISVar<T> && var )
  : _mod(var._mod), _ndiv(var._ndiv), _nvar(var._nvar), _ndep(var._ndep)
  {
#ifdef MC__ISMODEL_TRACE
    std::cerr << "-- ISVar( ISVar<T> && var )\n";
#endif
#ifdef TEST_MOVE
    std::cout << "Move Constructor" << std::endl;
#endif
    if( this == &var ) return;
    if( !_mod ) _cst = std::move(var._cst);
    //_mat = var._mat;
    //std::swap( _mat, var._mat );
/*
    _mat.swap( var._mat );
    _bnd = var._bnd;
    //std::swap( _slope, var._slope );
    _slope.swap( var._slope );
    _shadow.swap( var._shadow );
    _shadow_slope.swap( var._shadow_slope);
*/
    _mat = std::move(var._mat);
    _bnd = var._bnd;
    _slope = std::move(var._slope);
    _shadow = std::move(var._shadow);
    _shadow_slope = std::move(var._shadow_slope);     
    _shadow_info = std::move(var._shadow_info);    
}

  ~ISVar
  () 
  {
#ifdef FISM_LIFITIME_DEBUG 
    std::cout<< "ISV delated, nvar = " <<_nvar <<std::endl;
#endif    
  }

  ISVar<T>& set
  ( ISModel<T>* const mod )
  {
    _mod = mod;
    _nvar = mod->_nvar;
    _ndiv = mod->_ndiv;
    _ndep = 0;
    _mat.clear();
    _mat.resize( _nvar );
    _bnd = std::make_pair( 0., false );
    _slope.clear();
    _slope.resize( _nvar );
    _shadow.clear();
    _shadow.resize(4);
    _shadow_slope.clear();
    _shadow_slope.resize(4);
    for(unsigned i = 0; i< 4; i++){
      _shadow[i].resize(_nvar);
      _shadow_slope[i].resize(_nvar);
    }
    _shadow_info[0] = 1.;

    return *this;
  }

  ISVar<T>& set
  ( ISModel<T>* const mod, unsigned int ndx, T const& bnd )
  {
    _mod = mod;
    _nvar = mod->_nvar;
    if( ndx >= _nvar )
      throw typename ISModel<T>::Exceptions( ISModel<T>::Exceptions::INDEX );

    _ndiv = mod->_ndiv;
    _ndep = 1;
    _mat.clear();
    _mat.resize( _nvar );
    long double l = Op<T>::l(bnd);
    long double h = Op<T>::diam(bnd) / (long double)_ndiv;
    _mat[ndx].resize( _ndiv );
    for( unsigned int j=0; j<_ndiv - 1; j++, l+=h )
      _mat[ndx][j] = T( l, l+h );
      _mat[ndx][_ndiv - 1] = T( Op<T>::u(bnd) - h, Op<T>::u(bnd) );
    _mod->_defvar[ndx] = true;
    _mod->_bndvar[ndx] = bnd;
    _bnd = std::make_pair( bnd, true );

    _slope.clear();
    _slope.resize( _nvar );
    if (_mod->options.SLOPE_USE){
      _mod->_psize[ndx] = h;
      _slope[ndx].resize( _ndiv );
      for( unsigned int j=0; j<_ndiv; j++){
        _slope[ndx][j].resize(2);
        _slope[ndx][j][0] = 1.;
        _slope[ndx][j][1] = 1.;
      }
    }

    if (_mod->options.SHADOW_USE){
      //for( unsigned int j=0; j<_ndiv; j++){
      _shadow.clear();
      _shadow.resize(4);
      _shadow_slope.clear();
      _shadow_slope.resize(4);
      for(unsigned i = 0; i< 4; i++){
        _shadow[i].resize(_nvar);
        _shadow_slope[i].resize(_nvar);
      }
      _shadow_info[0] = 1.;        
        ;//initialize containers and info
      //}
    }
 

    return *this;
  }

  std::vector<std::vector<T>> const& C
  ()
  const
  { return _mat; }

  std::vector<std::vector<std::vector<long double>>> const& get_slopes
  ()
  const
  { return _slope; }

  std::vector<std::vector<std::vector<long double>>> const& get_shadow
  ()
  const
  { return _shadow; }

  std::vector<std::vector<std::vector<long double>>> const& get_shadow_slopes
  ()
  const
  { return _shadow_slope; }

  std::vector<double> const& get_shadow_info
  ()
  const
  { return _shadow_info; }


  std::vector<std::vector<T>>& C
  ()
  { return _mat; }

  T const& cst
  ()
  const
  { return _cst; }

  unsigned int const& ndep
  ()
  const
  { return _ndep; }

  unsigned int& ndep
  ()
  { return _ndep; }

  double ub
  ()
  const
  { return Op<T>::u(bound()); }

  double lb
  ()
  const
  { return Op<T>::l(bound()); }

  T B
  ()
  const
  { return bound(); }

  T bound
  ()
  const
  {
    if( _bnd.second ) return _bnd.first;
    _bnd.first = ( _mod? _mod->_B( _mat ): _cst );
    _bnd.second = true;
    return _bnd.first;
  }

  T eval
  ( double const* const point )
  const
  {
    assert( point );
    if( !_mod ) return _cst;
    T val( 0. );

    long double sum_underEnhancer = 0.0;
    long double sum_overEnhancer = 0.0;
    const long double longdoublezero = 0.;
    long double shadow_sign = 0.;
    //std::cout << "eval at point" <<std::endl;    
    for( unsigned int i=0; i<_nvar; i++ ){
      if( _mat[i].empty() ) { 
        //std::cout << "skip variable i = " << i << std::endl;   
        continue;
      }
      double l = Op<T>::l( _mod->_bndvar[i] );
      double u = Op<T>::u( _mod->_bndvar[i] );
      //std::cout << "the lower bound is " << l << std::endl;
      //std::cout << "the upper bound is " << u << std::endl;      
      if(point[i] < l || point[i] > u){
        //std::cout<<"l: "<<l<<" p: "<<point[i]<<" u: "<<u<<std::endl;
        const double _eps (0.00000000000001);
        if (point[i] - u >= _eps || point[i] - l <= -_eps)
          assert( point[i] >= l && point[i] <= u );
      }
    //   double h = Op<T>::diam( _mod->_bndvar[i] ) / (double)_ndiv;
    //   int ndx = std::floor( ( point[i] - l ) /  h );
    //   if( ndx < 0 ) ndx = 0;
    //   if( ndx >= (int)_ndiv ) ndx = _ndiv-1;      
    //   val += _mat[i][ndx];
      long double h = Op<T>::diam( _mod->_bndvar[i] ) / (long double)_ndiv;
      int ndx = std::floor( ( point[i] - l ) /  h ); 
      if( ndx < 0 ) ndx = 0;
      if( ndx >= (int)_ndiv ) ndx = _ndiv-1;      
      val += _mat[i][ndx]; 
       
      //std::cout << "val computed at " << i <<std::endl;
      
      if ( _mod->options.SLOPE_USE){
        //std::cout << " in slope " << i <<std::endl;
        long double valL = Op<T>::l(val);
        long double valU = Op<T>::u(val);
        if (_slope[i][ndx][0]>MC__ISM_COMPUTATION_TOL){
          valL = valL + _slope[i][ndx][0]*(point[i] - h*((double) ndx) - l); 
          //std::cout << "pl " << std::setprecision(18) << _slope[i][ndx][0]*(point[i] - h*((double) ndx) - l) << "," << h << std::endl;
        }
        else if (_slope[i][ndx][0] < -MC__ISM_COMPUTATION_TOL){
          valL = valL + _slope[i][ndx][0]*(point[i] - h*((double)(ndx+1)) - l);
          //std::cout << "pu " << std::setprecision(18) << _slope[i][ndx][0]*(point[i] - h*((double)(ndx+1)) - l) << "," << h << std::endl;
        }
          
        if (_slope[i][ndx][1]>MC__ISM_COMPUTATION_TOL)
          valU = valU + _slope[i][ndx][1]*(point[i] - h*((double)(ndx+1)) - l);
        else if (_slope[i][ndx][1] < -MC__ISM_COMPUTATION_TOL)
          valU = valU + _slope[i][ndx][1]*(point[i] - h*((double) ndx) - l);

        if (valU - valL > MC__ISM_COMPUTATION_TOL)
          val = T(valL,valU);
        else if (valU - valL < -1e2*MC__ISM_COMPUTATION_TOL){
          std::cout << "Error in slope: at " << i << "," << ndx << std::endl; 
          std::cout << std::setprecision(18) << valU - valL << "," << valU << "," << valL << std::endl; 
        }
        else 
          val = T(std::min(valL,valU),std::max(valL,valU));

        if(_mod->options.SHADOW_USE && _shadow_info.size()>1){
          //std::cout << "shadow" << std::endl;
          long double underEnhancer = 0.;
          long double overEnhancer = 0.;
          if(_shadow_info[0]>0)
            shadow_sign = 1.;
          else if(_shadow_info[0]<0)
            shadow_sign = -1.;
          else
            throw typename ISModel<T>::Exceptions( ISModel<T>::Exceptions::UNDEF );
 

          if(_shadow_info[1] == 1){
            underEnhancer = _shadow[0][i][ndx];
            if (_shadow_slope[0][i][ndx]>MC__ISM_COMPUTATION_TOL){
              underEnhancer = underEnhancer + _shadow_slope[0][i][ndx]*(point[i] - h*((double) ndx) - l); 
              //std::cout << "pl " << std::setprecision(18) << _slope[i][ndx][0]*(point[i] - h*((double) ndx) - l) << "," << h << std::endl;
            }
            else if (_shadow_slope[0][i][ndx] < -MC__ISM_COMPUTATION_TOL){
              underEnhancer = underEnhancer + _shadow_slope[0][i][ndx]*(point[i] - h*((double)(ndx+1)) - l);
              //std::cout << "pu " << std::setprecision(18) << _slope[i][ndx][0]*(point[i] - h*((double)(ndx+1)) - l) << "," << h << std::endl;
            }
            //underEnhancer = std::max(shadow_sign*underEnhancer,longdoublezero);
          }
          if(_shadow_info[2] == 1){
            overEnhancer = _shadow[1][i][ndx];    
            if (_shadow_slope[1][i][ndx]>MC__ISM_COMPUTATION_TOL)
              overEnhancer = overEnhancer + _shadow_slope[1][i][ndx]*(point[i] - h*((double)(ndx+1)) - l);
            else if (_shadow_slope[1][i][ndx] < -MC__ISM_COMPUTATION_TOL)
              overEnhancer = overEnhancer + _shadow_slope[1][i][ndx]*(point[i] - h*((double) ndx) - l);
            //overEnhancer = std::min(shadow_sign*overEnhancer,longdoublezero);
          }
          //std::cout << "proccessed" <<std::endl;
          sum_underEnhancer += underEnhancer;
          sum_overEnhancer += overEnhancer;
        }
      }
        
    }

    //std::cout << "after for loop" <<std::endl;   

    if(_mod->options.SHADOW_USE ){
      if(_shadow_info.size()>1){
        long double valL = Op<T>::l(val)+std::max(shadow_sign*sum_underEnhancer,longdoublezero);

        // if(std::max(shadow_sign*sum_underEnhancer,longdoublezero)>0.)
        //   std::cout << std::max(shadow_sign*sum_underEnhancer,longdoublezero) << std::endl;

        long double valU = Op<T>::u(val)+std::min(shadow_sign*sum_overEnhancer,longdoublezero);       
        if (valU - valL > MC__ISM_COMPUTATION_TOL)
          val = T(valL,valU);
        else if (valU - valL < -1e2*MC__ISM_COMPUTATION_TOL){
          std::cout << "Error in shadowslope: " << std::endl; 
          std::cout << std::setprecision(18) << "valU - valL = " <<valU - valL << "," 
                    << "valU = " << valU << ", valL = " << valL << std::endl; 
        }
        else 
          val = T(std::min(valL,valU),std::max(valL,valU));    
      }

    }

    return val;
  }

  std::ostream& display
  ( const int& opt=0, std::ostream& out=std::cout )
  const
  {
    if( _mod ) _mod->_dispvar( _mat, _ndep, opt, out );
    return out;
  }
  
 private:


  ISVar
  ( ISModel<T> const* const mod )
  : _mod(mod), _ndiv(mod->_ndiv), _nvar(mod->_nvar ), _mat(_nvar), _bnd(0.,false), _slope(_nvar),_shadow(4), _shadow_slope(4),_shadow_info(1)
  {
    for(unsigned i = 0; i< 4; i++){
      _shadow[i].resize(_nvar);
      _shadow_slope[i].resize(_nvar);
    }  
  }

};

template <typename T> 
std::ostream& operator<<
( std::ostream& out, ISVar<T> const& var )
{
  auto&& mat = var._mat;
  for( unsigned int i=0; !var._mat.empty() && i<var._nvar; i++ ){
    if( var._mat[i].empty() ) continue;
    out << std::right << std::setw(5) << i <<": ";
    for( unsigned int j=0; j<var._ndiv; j++ ){
      if( j && !(j%3) ) out << std::endl << "       ";
      out << std::setw(0) << mat[i][j];
    }
    out << std::endl;
  }
  out << std::right << std::setw(7) << "B: " << var.B() << std::endl;    
  return out;
}

template <typename T>
inline
ISVar<T> operator+
( ISVar<T> const& var )
{
  return var;
}

template <typename T>
inline
ISVar<T> operator+
( ISVar<T> && varIn )
{
  ISVar<T> var( std::move(varIn) );  
  return var;
}


template <typename T>
inline
ISVar<T>& ISVar<T>::operator+=
( double const& cst )
{
  if( cst == 0. ){
    return *this;
  }
  if( !_mod ){
    _cst += cst;
    return *this;
  }
  if( !_ndep )
    throw typename ISModel<T>::Exceptions( ISModel<T>::Exceptions::INTERN );

  const double offset = cst / (double)_ndep;

  for( unsigned int i=0; i<_nvar; i++ ){
    if( _mat[i].empty() ) continue;
    for( unsigned int j=0; j<_ndiv; j++ )
      _mat[i][j] += offset;
  }
  if( _bnd.second ) _bnd.first += cst;

  return *this;
}

template <typename T>
inline
ISVar<T>& ISVar<T>::operator+=
( T const& cst )
{

  if( !_mod ){
    _cst += cst;
    return *this;
  }
  if( !_ndep )
    throw typename ISModel<T>::Exceptions( ISModel<T>::Exceptions::INTERN );

  const T offset = cst / (double)_ndep;

  for( unsigned int i=0; i<_nvar; i++ ){
    if( _mat[i].empty() ) continue;
    for( unsigned int j=0; j<_ndiv; j++ )
      _mat[i][j] += offset;
  }
  if( _bnd.second ) _bnd.first += cst;

  return *this;
}


template <typename T>
inline
ISVar<T>& ISVar<T>::operator+=
( ISVar<T> const& var )
{
  if( !_mod && !var._mod ){
    _cst += var._cst;
    return *this;
  }
  if( !_mod ){
    T copy_cst = _cst;
    *this = var;
    *this += copy_cst;
    return *this;
  }
  if( !var._mod ){
    *this += var._cst;
    return *this;
  }
  if( _mod != var._mod )
    throw typename ISModel<T>::Exceptions( ISModel<T>::Exceptions::MODEL );

  if (_mod->options.SHADOW_USE){
#ifndef NOTTOTRACKSHADOW    
    std::cout << "SHADOW ADDITION" << std::endl;
#endif
    bool toAggregate = false;
    if (_shadow_info.size() > 1 && var._shadow_info.size() > 1){
      if (_shadow_info[1] > 0 || _shadow_info[2] > 0 || var._shadow_info[1]>0 || var._shadow_info[2] > 0)
        toAggregate = true;
      else{
        std::cout << "error in shadow info in processing binary addition" << std::endl;
        throw typename ISModel<T>::Exceptions( ISModel<T>::Exceptions::UNDEF );
      }  
    }
    else if(_shadow_info.size() > 1){
      var._shadow_info.push_back(0.);
      var._shadow_info.push_back(0.);
      toAggregate = true;
    }
    else if(var._shadow_info.size() > 1){
      _shadow_info.push_back(0.);
      _shadow_info.push_back(0.);
      toAggregate = true;
    }  
    if (toAggregate){

#ifndef NOTTOTRACKSHADOW
      for (size_t k = 0;k<_shadow_info.size();k++){
        std::cout << "_shadow_info[" << k << "] " <<_shadow_info[k] << "\n";
      }
      std::cout << std::endl;
#endif
#ifndef NOTTOTRACKSHADOW           
      for (size_t k = 0;k<var._shadow_info.size();k++){
        std::cout << "var_shadow_info[" << k << "] " << var._shadow_info[k] << "\n";
      }
      std::cout << std::endl;
#endif

       //update ndep
      for( unsigned int i=0; i<_nvar; i++ ){
        if( !_mat[i].empty() && !var._mat[i].empty() )
          ;
        else if( !var._mat[i].empty() ){
          _ndep++;                // add dependency
        }
      }     

      T bndB = var.bound();
      _mod->_add_aggregate_shadow( _mat, var._mat, _slope, var._slope, _shadow, var._shadow, _shadow_slope, var._shadow_slope, _shadow_info, var._shadow_info, _mod->_psize, _ndep,bndB);
      return *this; 
    } 
  }

  for( unsigned int i=0; i<_nvar; i++ ){
    if( !_mat[i].empty() && !var._mat[i].empty() )
      for( unsigned int j=0; j<_ndiv; j++ )
        _mat[i][j] += var._mat[i][j];
    else if( !var._mat[i].empty() ){
      _ndep++;                // add dependency
      _mat[i] = var._mat[i];  // copy entire row
    }
  }
  _bnd.second = false;

  //use slope information to tighten interval bounds
  if (_mod->options.SLOPE_USE){
    for( unsigned int i=0; i<_nvar; i++ ){
      //std::cout << "slope cases " << i << ":" <<_slope[i].empty() <<" , "<< var._slope[i].empty() << std::endl;
      if( (!_slope[i].empty()) && (!var._slope[i].empty()) ){
        //std::cout << "    Slope addition" << std::endl;
        for( unsigned int j=0; j<_ndiv; j++ ){
          if((_slope[i][j][0] < 0.) != (var._slope[i][j][0] < 0.)){   
            long double tightener = std::min(std::fabs(_slope[i][j][0]),std::fabs(var._slope[i][j][0]))*(_mod->_psize[i]);
            if ( Op<T>::u(_mat[i][j]) - tightener - Op<T>::l(_mat[i][j]) < MC__ISM_COMPUTATION_TOL)
              _mat[i][j] = T(Op<T>::l(_mat[i][j]) + tightener);
            else
              _mat[i][j] = T(Op<T>::l(_mat[i][j]) + tightener , Op<T>::u(_mat[i][j]));
          }
 
          if((_slope[i][j][1] < 0.) != (var._slope[i][j][1] < 0.)){   
            long double tightener = std::min(std::fabs(_slope[i][j][1]),std::fabs(var._slope[i][j][1]))*(_mod->_psize[i]);
            if ( Op<T>::u(_mat[i][j]) - tightener - Op<T>::l(_mat[i][j]) < MC__ISM_COMPUTATION_TOL)
              _mat[i][j] = T(Op<T>::u(_mat[i][j]) - tightener);
            else 
              _mat[i][j] = T(Op<T>::l(_mat[i][j]) , Op<T>::u(_mat[i][j]) - tightener);
         }
          _slope[i][j][0] += var._slope[i][j][0];
          _slope[i][j][1] += var._slope[i][j][1];
        }
      }
      else if( !var._slope[i].empty() ){
        //std:: cout << "copy row " << i <<std::endl; 
        _slope[i] = var._slope[i];  // copy entire row
      }
    }
    //std::cout << "slope cases ended" << std::endl;
  }

  return *this;
}





template <typename T>
inline
ISVar<T> operator+
( ISVar<T> const& var1, ISVar<T> const& var2 )
{
  if( !var1._mod && !var2._mod ) 
    return( var1._cst + var2._cst );
  if( var1._mod ){
    ISVar<T> var3( var1 );
    var3 += var2;
    return var3;
  }
  ISVar<T> var3( var2 );
    var3 += var1;
    return var3;
}


template <typename T>
inline
ISVar<T> operator+
( ISVar<T> const& var1, ISVar<T> && var2 )
{
  //std::cout << "+ 1move" << std::endl;
  if( !var1._mod && !var2._mod ) 
    return( var1._cst + var2._cst );
   
  if( var2._mod ){
    ISVar<T> var3( std::move(var2) );
    var3 += var1;
    return var3;
  }
  ISVar<T> var3( var1 );
  var3 += var2;
  return var3;
}


template <typename T>
inline
ISVar<T> operator+
( ISVar<T> && var1, ISVar<T> const& var2 )
{
  //std::cout << "+ 2move" << std::endl;
  if( !var1._mod && !var2._mod ) 
    return( var1._cst + var2._cst );
   
  if( var1._mod ){
    ISVar<T> var3( std::move(var1) );
    var3 += var2;
    return var3;
  }
  ISVar<T> var3( var2 );
  var3 += var1;
  return var3;
}


template <typename T>
inline
ISVar<T> operator+
( ISVar<T> && var1, ISVar<T> && var2 )
{
  //std::cout << "+ 3move" << std::endl;
  if( !var1._mod && !var2._mod ) 
    return( var1._cst + var2._cst );
  if( var1._mod ){
    ISVar<T> var3( std::move(var1) );
    var3 += var2;
    return var3;
  }
  ISVar<T> var3( std::move(var2) );
    var3 += var1;
    return var3;
}



template <typename T>
inline
ISVar<T> operator+
( ISVar<T> const& var1, double const& cst2 )
{
  if( !var1._mod ) 
    return( var1._cst + cst2 );
  ISVar<T> var3( var1 );
  var3 += cst2;
  return var3;
}


template <typename T>
inline
ISVar<T> operator+
( ISVar<T> && var1, double const& cst2 )
{
  if( !var1._mod ) 
    return( var1._cst + cst2 );
  ISVar<T> var3( std::move(var1) );
  var3 += cst2;
  return var3;
}


template <typename T>
inline
ISVar<T> operator+
( double const& cst1, ISVar<T> const& var2 )
{
  if( !var2._mod ) 
    return( cst1 + var2._cst );
  ISVar<T> var3( var2 );
  var3 += cst1;
  return var3;
}

template <typename T>
inline
ISVar<T> operator+
( double const& cst1, ISVar<T> && var2 )
{
  if( !var2._mod ) 
    return( cst1 + var2._cst );
  ISVar<T> var3( std::move(var2) );
  var3 += cst1;
  return var3;
}



template <typename T> 
inline
ISVar<T> operator-
( ISVar<T> const& var )
{
  ISVar<T> var2( var );
  if( !var2._mod ){
    var2._cst *= -1;
  }
  else{ 
    auto&& mat2 = var2._mat;
    for( unsigned int i=0; i<var2._nvar; i++ ){
      if( mat2[i].empty() ) continue;
      for( unsigned int j=0; j<var2._ndiv; j++ )
        mat2[i][j] *= -1;
    }

    if(var2._mod->options.SLOPE_USE){
      auto&& slope2 = var2._slope;
      for( unsigned int i=0; i<var2._nvar; i++ ){
        if( mat2[i].empty() ) continue;
        for( unsigned int j=0; j<var2._ndiv; j++ ){
          slope2[i][j][0] = -1*slope2[i][j][0];
          slope2[i][j][1] = -1*slope2[i][j][1];
          std::swap(slope2[i][j][0],slope2[i][j][1]);
        }
      }      
    }

    if(var2._mod->options.SHADOW_USE){
      if(var2._shadow_info.size()>1){
        var2._shadow[0].swap(var2._shadow[1]);
//        var2._shadow[2].swap(var2._shadow[3]);
        var2._shadow_slope[0].swap(var2._shadow_slope[1]);
//        var2._shadow_slope[2].swap(var2._shadow_slope[3]);      
        var2._shadow_info[0] *= -1.;

        std::swap(var2._shadow_info[1],var2._shadow_info[2]);
        std::swap(var2._shadow_info[3],var2._shadow_info[4]);
        for ( unsigned int i=0; i< var2._nvar; i++ ){
          if( mat2[i].empty() ) continue;
          std::swap(var2._shadow_info[2*i+5],var2._shadow_info[2*i+6]);       
        }

      }
    }
  }
  auto&& bnd = var2._bnd;
  if( bnd.second ) bnd.first *= -1;

  return var2;
}


template <typename T> 
inline
ISVar<T> operator-
( ISVar<T> && var )
{
  ISVar<T> var2( std::move(var) );
  if( !var2._mod ){
    var2._cst *= -1;
  }
  else{ 
    auto&& mat2 = var2._mat;
    for( unsigned int i=0; i<var2._nvar; i++ ){
      if( mat2[i].empty() ) continue;
      for( unsigned int j=0; j<var2._ndiv; j++ )
        mat2[i][j] *= -1;
    }

    if(var2._mod->options.SLOPE_USE){
      auto&& slope2 = var2._slope;
      for( unsigned int i=0; i<var2._nvar; i++ ){
        if( mat2[i].empty() ) continue;
        for( unsigned int j=0; j<var2._ndiv; j++ ){
          slope2[i][j][0] = -1*slope2[i][j][0];
          slope2[i][j][1] = -1*slope2[i][j][1];
          std::swap(slope2[i][j][0],slope2[i][j][1]);
        }
      }      
    }

    if(var2._mod->options.SHADOW_USE){
      if(var2._shadow_info.size()>1){
        var2._shadow[0].swap(var2._shadow[1]);
//        var2._shadow[2].swap(var2._shadow[3]);
        var2._shadow_slope[0].swap(var2._shadow_slope[1]);
//        var2._shadow_slope[2].swap(var2._shadow_slope[3]);      
        var2._shadow_info[0] *= -1.;

        std::swap(var2._shadow_info[1],var2._shadow_info[2]);
        std::swap(var2._shadow_info[3],var2._shadow_info[4]);
        for ( unsigned int i=0; i< var2._nvar; i++ ){
          if( mat2[i].empty() ) continue;
          std::swap(var2._shadow_info[2*i+5],var2._shadow_info[2*i+6]);       
        }

      }
    }
  }
  auto&& bnd = var2._bnd;
  if( bnd.second ) bnd.first *= -1;

  return var2;
}



template <typename T>
inline
ISVar<T>& ISVar<T>::operator-=
( double const& cst )
{
  *this += -cst;
  return *this;
}

template <typename T>
inline
ISVar<T>& ISVar<T>::operator-=
( T const& cst )
{
  *this += -cst;
  return *this;
}

template <typename T>
inline
ISVar<T>& ISVar<T>::operator-=
( ISVar<T> const& var )
{
  if( !_mod && !var._mod ){
    _cst -= var._cst;
    return *this;
  }
  if( !_mod ){
    T copy_cst = _cst;
    *this = -var;
    *this += copy_cst;
    return *this;
  }
  if( !var._mod ){
    *this -= var._cst;
    return *this;
  }
  if( _mod != var._mod )
    throw typename ISModel<T>::Exceptions( ISModel<T>::Exceptions::MODEL );

  if(_mod->options.SHADOW_USE){
    std::cout << "operator -= error" << std::endl;
    throw typename ISModel<T>::Exceptions( ISModel<T>::Exceptions::UNDEF );
  }

  for( unsigned int i=0; i<_nvar; i++ ){
    if( !_mat[i].empty() && !var._mat[i].empty() )
      for( unsigned int j=0; j<_ndiv; j++ )
        _mat[i][j] -= var._mat[i][j];
    else if( !var._mat[i].empty() ){
      _ndep++;            // add dependency
      _mat[i].resize(_ndiv);
      for( unsigned int j=0; j<_ndiv; j++ )
        _mat[i][j] = -var._mat[i][j];
    }
  }
  _bnd.second = false;

  // use slope information to tighten interval bounds
  if (_mod->options.SLOPE_USE){
    for( unsigned int i=0; i<_nvar; i++ ){
      if( !_slope[i].empty() && !var._slope[i].empty() )
        for( unsigned int j=0; j<_ndiv; j++ ){

          if((_slope[i][j][0] < 0.) == (var._slope[i][j][1] < 0.)){   
            long double tightener = std::min(std::fabs(_slope[i][j][0]),std::fabs(var._slope[i][j][1]))*(_mod->_psize[i]);
            if ( Op<T>::u(_mat[i][j]) - tightener - Op<T>::l(_mat[i][j]) < MC__ISM_COMPUTATION_TOL)
              _mat[i][j] = T(Op<T>::l(_mat[i][j]) + tightener);
            else
              _mat[i][j] = T(Op<T>::l(_mat[i][j]) + tightener , Op<T>::u(_mat[i][j]));
          }

          if((_slope[i][j][1] < 0.) == (var._slope[i][j][0] < 0.)){   
            long double tightener = std::min(std::fabs(_slope[i][j][1]),std::fabs(var._slope[i][j][0]))*(_mod->_psize[i]);
            if ( Op<T>::u(_mat[i][j]) - tightener - Op<T>::l(_mat[i][j]) < MC__ISM_COMPUTATION_TOL)
              _mat[i][j] = T(Op<T>::u(_mat[i][j]) - tightener);
            else
              _mat[i][j] = T(Op<T>::l(_mat[i][j]) , Op<T>::u(_mat[i][j]) - tightener);
         }

          _slope[i][j][0] -= var._slope[i][j][1];
          _slope[i][j][1] -= var._slope[i][j][0];
        }
      else if( !var._slope[i].empty() ){
        _slope[i].resize(_ndiv);
        for( unsigned int j=0; j<_ndiv; j++ ){
          _slope[i][j].resize(2);
          _slope[i][j][0] = -var._slope[i][j][1];
          _slope[i][j][1] = -var._slope[i][j][0];
        }
      }
    }
  }


  return *this;
}

template <typename T>
inline
ISVar<T> operator-
( ISVar<T> const& var1, ISVar<T> const& var2 )
{
  if( !var1._mod && !var2._mod ) 
    return( var1._cst - var2._cst );
  if( var1._mod ){
    ISVar<T> var3( var1 );
    if(var1._mod->options.SHADOW_USE){
      var3 += (-var2);
    }
    else
      var3 -= var2;    
    return var3;
  }
  ISVar<T> var3( -var2 );
    var3 += var1;
    return var3;
}



template <typename T>
inline
ISVar<T> operator-
( ISVar<T> const& var1, ISVar<T> && var2 )
{
  if( !var1._mod && !var2._mod ) 
    return( var1._cst - var2._cst );
  if( var2._mod ){
    ISVar<T> var3( std::move(-var2) );
    var3 += var1;  
    return var3;
  }
  ISVar<T> var3( var1 );
  if(var1._mod->options.SHADOW_USE){
    var3 += (-var2);
  }
  else
    var3 -= var2;    
  return var3;
}


template <typename T>
inline
ISVar<T> operator-
( ISVar<T> && var1, ISVar<T> const& var2 )
{
  if( !var1._mod && !var2._mod ) 
    return( var1._cst - var2._cst );
  if( var1._mod ){
    ISVar<T> var3( std::move(var1) );
    if(var1._mod->options.SHADOW_USE){
      var3 += (-var2);
    }
    else
      var3 -= var2;    
    return var3;
  }
  ISVar<T> var3( -var2 );
    var3 += var1;
    return var3;
}


template <typename T>
inline
ISVar<T> operator-
( ISVar<T> && var1, ISVar<T> && var2 )
{
  if( !var1._mod && !var2._mod ) 
    return( var1._cst - var2._cst );
  if( var1._mod ){
    ISVar<T> var3( std::move(var1) );
    if(var1._mod->options.SHADOW_USE){
      var3 += (-var2);
    }
    else
      var3 -= var2;    
    return var3;
  }
  ISVar<T> var3( std::move(-var2) );
    var3 += var1;
    return var3;
}



template <typename T>
inline
ISVar<T> operator-
( ISVar<T> const& var1, double const& cst2 )
{
  if( !var1._mod ) 
    return( var1._cst - cst2 );
  ISVar<T> var3( var1 );
  var3 += -cst2;
  return var3;
}

template <typename T>
inline
ISVar<T> operator-
( ISVar<T> && var1, double const& cst2 )
{
  if( !var1._mod ) 
    return( var1._cst - cst2 );
  ISVar<T> var3( std::move(var1) );
  var3 += -cst2;
  return var3;
}

template <typename T>
inline
ISVar<T> operator-
( double const& cst1, ISVar<T> const& var2 )
{
  if( !var2._mod ) 
    return( cst1 - var2._cst );
  ISVar<T> var3( -var2 );
  var3 += cst1;
  return var3;
}

template <typename T>
inline
ISVar<T> operator-
( double const& cst1, ISVar<T> && var2 )
{
  if( !var2._mod ) 
    return( cst1 - var2._cst );
  ISVar<T> var3( std::move(-var2) );
  var3 += cst1;
  return var3;
}


template <typename T>
inline
ISVar<T>& ISVar<T>::operator*=
( double const& cst )
{
  if( cst == 0. ){
    *this = 0.;
    return *this;
  }
  if( cst == 1. ){
    return *this;
  }
  if( !_mod ){
    _cst *= cst;
    return *this;
  }
  if( !_ndep )
    throw typename ISModel<T>::Exceptions( ISModel<T>::Exceptions::INTERN );

  for( unsigned int i=0; i<_nvar; i++ ){
    if( _mat[i].empty() ) continue;
    for( unsigned int j=0; j<_ndiv; j++ )
      _mat[i][j] *= cst;
  }
  if( _bnd.second ) _bnd.first *= cst;

  if (_mod->options.SLOPE_USE){
    if (cst >=0.){
      for( unsigned int i=0; i<_nvar; i++ ){
        if( _slope[i].empty() ) continue;
        for( unsigned int j=0; j<_ndiv; j++ ){
          _slope[i][j][0] *= cst;
          _slope[i][j][1] *= cst;
        }
      } 
    }
    else {
      for( unsigned int i=0; i<_nvar; i++ ){
        if( _mat[i].empty() ) continue;
        for( unsigned int j=0; j<_ndiv; j++ ){
          _slope[i][j][0] *= cst;
          _slope[i][j][1] *= cst;
          std::swap(_slope[i][j][0],_slope[i][j][1]);
        }
      } 
    } 
  }

  if (_mod->options.SHADOW_USE && _ndep > 1){
#ifndef NOTTOTRACKSHADOW
    for (unsigned i = 0; i < _shadow_info.size();i++){
      std::cout << "_shadow_info[" << i << "] = " <<_shadow_info[i] <<"\n";    
    }
    std::cout << std::endl;
#endif
    if (_shadow_info.size() > 1){
      if(cst < 0.){
        std::swap(_shadow_info[1],_shadow_info[2]); 
        _shadow[0].swap(_shadow[1]);
//        _shadow[2].swap(_shadow[3]);
        _shadow_slope[0].swap(_shadow_slope[1]);
//        _shadow_slope[2].swap(_shadow_slope[3]);        
      }
#ifndef NOTTOTRACKSHADOW
      std::cout << "    shadow info swaped" << std::endl;    
#endif      
      for ( unsigned int i=0; i<_nvar+1; i++ ){
        if( _mat[i].empty() ) continue;
        _shadow_info[2*i+3] *= cst;
        _shadow_info[2*i+4] *= cst; 
        if(cst < 0.){
          std::swap(_shadow_info[2*i+3],_shadow_info[2*i+4]);   
        }    
      }

      for( unsigned int i=0; i<_nvar; i++ ){
        if( _shadow[0][i].empty() ) continue;
        for( unsigned int j=0; j<_ndiv; j++ ){
            _shadow[0][i][j] *= cst;
            _shadow_slope[0][i][j] *= cst;
            // _shadow[2][i][j] *= cst;
            // _shadow_slope[2][i][j] *= cst;          
        }
      }      
#ifndef NOTTOTRACKSHADOW
      std::cout << "    shadow under swaped" << std::endl;    
#endif      
      for( unsigned int i=0; i<_nvar; i++ ){
        if( _shadow[1][i].empty() ) continue;
        for( unsigned int j=0; j<_ndiv; j++ ){
            _shadow[1][i][j] *= cst;
            _shadow_slope[1][i][j] *= cst;
            // _shadow[3][i][j] *= cst;
            // _shadow_slope[3][i][j] *= cst;
        }
      } 
#ifndef NOTTOTRACKSHADOW
      std::cout << "    shadow over swaped" << std::endl;    
#endif            
    }       


  }
  return *this;
}



template <typename T>
inline
ISVar<T>& ISVar<T>::operator*=
( T const& cst )
{

  double cstl = Op<T>::l( cst );
  double cstu = Op<T>::u( cst );

  if( cstl == cstu){
    return ((*this)*=cstl);
  }
  else if( cstu-cstl <= MC__ISM_COMPUTATION_TOL ){
    (*this)*=(0.5*(cstl+cstu));
    return *this;
  }
  else {
    throw typename ISModel<T>::Exceptions( ISModel<T>::Exceptions::UNDEF );
    return *this;
  }
    
}


template <typename T>
inline
ISVar<T>& ISVar<T>::operator*=
( ISVar<T> const& var )
{
  if( !_mod && !var._mod ){
    _cst *= var._cst;
    return *this;
  }
  if( !_mod ){
    T copy_cst = _cst;
    *this = var;
    *this *= copy_cst;
    return *this;
  }
  if( !var._mod ){
    *this *= var._cst;
    return *this;
  }
  if( _mod != var._mod )
    throw typename ISModel<T>::Exceptions( ISModel<T>::Exceptions::MODEL );

  if( var._mod->options.SLOPE_USE ){
    throw typename ISModel<T>::Exceptions( ISModel<T>::Exceptions::UNDEF );
  }


  // To be developed: currently =* does not allow rescalling
  _mod->_prod( _mat, var._mat, _mat, _ndep );
  _bnd.second = false;

  return *this;
}

template <typename T>
inline
ISVar<T> operator*
( ISVar<T> const& var1, ISVar<T> const& var2 )
{
  if( !var1._mod && !var2._mod ) 
    return( var1._cst * var2._cst );
  if( !var2._mod ){
    ISVar<T> var3( var1 );
    var3 *= var2._cst;
    return var3;
  }
  if( !var1._mod ){
    ISVar<T> var3( var2 );
    var3 *= var1._cst;
    return var3;
  }

  // if( var1._mod->options.SLOPE_USE ){
  //   return sqr( 0.5*var1 + 0.5*var2 ) - sqr( 0.5*var1 - 0.5*var2 );
  // }

  auto bnd1 = var1.bound();
  auto bnd2 = var2.bound();
  if( var1._mod->options.DCDEC_USE ){
    //double coef1(0.),coef2(0.),mid1(0.),mid2(0.);
    
    double coef1 = 0.5 * (Op<T>::u(bnd1)-Op<T>::l(bnd1)); //0.5*(bnd1.first.u()-bnd1.first.l());
    double coef2 = 0.5 * (Op<T>::u(bnd2)-Op<T>::l(bnd2)); //0.5*(bnd2.first.u()-bnd2.first.l());
    double mid1 = 0.5 * (Op<T>::u(bnd1)+Op<T>::l(bnd1));
    double mid2 = 0.5 * (Op<T>::u(bnd2)+Op<T>::l(bnd2));
    
    //std::cout<<Op<T>::u(bnd1*bnd2)<<std::endl;
    //return coef1 * coef2 * 0.25 * ( sqr( var1/coef1 + var2/coef2 ) - sqr( var1/coef1 - var2/coef2 ) );      
    //return intersect( coef1 * coef2 * 0.25 * ( sqr( var1/coef1 + var2/coef2 ) - sqr( var1/coef1 - var2/coef2 ) ), bnd1*bnd2);
    //return 0.25 * ( sqr( var1 + var2 ) - sqr( var1 - var2 ) );
    //return intersect( 0.25 * ( sqr( var1 + var2 ) - sqr( var1 - var2 ) ), bnd1*bnd2);
    //return coef1 * coef2 * 0.25 * ( sqr( (var1-mid1)/coef1 + (var2-mid2)/coef2 ) - sqr( (var1-mid1)/coef1 - (var2-mid2)/coef2 ) ) + mid1*var2+mid2*var1-mid1*mid2;
    auto _type = var1._mod->options.SCALING_TYPE;
    //auto _type_temp = _type;
    if(_type == ISModel<T>::Options::ADAPT){
      T rst1(0.25 * ( Op<T>::sqr( bnd1 + bnd2 ) - Op<T>::sqr( bnd1 - bnd2 ) ));
      T rst2(coef1 * coef2 * 0.25 * ( Op<T>::sqr( bnd1/coef1 + bnd2/coef2 ) - Op<T>::sqr( bnd1/coef1 - bnd2/coef2 ) ));
      T rst3(coef1 * coef2 * 0.25 * ( Op<T>::sqr( (bnd1-mid1)/coef1 + (bnd2-mid2)/coef2 ) - Op<T>::sqr( (bnd1-mid1)/coef1 - (bnd2-mid2)/coef2 ) ) + mid1*bnd2+mid2*bnd1-mid1*mid2);
      //std::cout<<Op<T>::u(rst1)-Op<T>::l(rst1)<<"  "<<Op<T>::u(rst2)-Op<T>::l(rst2)<<"  "<<Op<T>::u(rst3)-Op<T>::l(rst3)<<std::endl; 
      if ((Op<T>::u(rst1)-Op<T>::l(rst1)) <= (Op<T>::u(rst2)-Op<T>::l(rst2))){
        if((Op<T>::u(rst1)-Op<T>::l(rst1)) <= (Op<T>::u(rst3)-Op<T>::l(rst3)))
          _type = ISModel<T>::Options::NONE;
        else
          _type = ISModel<T>::Options::FULL;  
      }
      else{
        if((Op<T>::u(rst2)-Op<T>::l(rst2)) <= (Op<T>::u(rst3)-Op<T>::l(rst3)))
          _type = ISModel<T>::Options::PARTIAL;
        else
          _type = ISModel<T>::Options::FULL;  
      }            
    }
    switch( _type ){          
      case ISModel<T>::Options::PARTIAL:{
//        return (var1._mod->options.INTERSECTION_USE?intersect(coef1 * coef2 * 0.25 * ( sqr( var1/coef1 + var2/coef2 ) - sqr( var1/coef1 - var2/coef2 ) ),bnd1*bnd2):coef1 * coef2 * 0.25 * ( sqr( var1/coef1 + var2/coef2 ) - sqr( var1/coef1 - var2/coef2 ) ) );
        return (var1._mod->options.INTERSECTION_USE?intersect((coef1 * coef2) * ( sqr( var1/(2.0*coef1) + var2/(2.0*coef2) ) - sqr( var1/(2.0*coef1) - var2/(2.0*coef2) ) ),bnd1*bnd2)
                                                             :(coef1 * coef2) * ( sqr( var1/(2.0*coef1) + var2/(2.0*coef2) ) - sqr( var1/(2.0*coef1) - var2/(2.0*coef2) ) ));
        break;
      }
      case ISModel<T>::Options::FULL:{
//        return (var1._mod->options.INTERSECTION_USE?intersect(coef1 * coef2 * 0.25 * ( sqr( (var1-mid1)/coef1 + (var2-mid2)/coef2 ) - sqr( (var1-mid1)/coef1 - (var2-mid2)/coef2 ) ) + mid1*var2+mid2*var1-mid1*mid2,bnd1*bnd2):coef1 * coef2 * 0.25 * ( sqr( (var1-mid1)/coef1 + (var2-mid2)/coef2 ) - sqr( (var1-mid1)/coef1 - (var2-mid2)/coef2 ) ) + mid1*var2+mid2*var1-mid1*mid2);
        return (var1._mod->options.INTERSECTION_USE? intersect((coef1 * coef2) * ( sqr( (var1-mid1)/(2.0*coef1) + (var2-mid2)/(2.0*coef2) ) - sqr( (var1-mid1)/(2.0*coef1) - (var2-mid2)/(2.0*coef2) ) ) + mid1*var2+mid2*var1-mid1*mid2,bnd1*bnd2)
                                                              :(coef1 * coef2) * ( sqr( (var1-mid1)/(2.0*coef1) + (var2-mid2)/(2.0*coef2) ) - sqr( (var1-mid1)/(2.0*coef1) - (var2-mid2)/(2.0*coef2) ) ) + mid1*var2+mid2*var1-mid1*mid2);
        break;
      }
      case ISModel<T>::Options::NONE: default:{
//        return (var1._mod->options.INTERSECTION_USE?intersect(0.25 * ( sqr( var1 + var2 ) - sqr( var1 - var2 ) ),bnd1*bnd2):0.25 * ( sqr( var1 + var2 ) - sqr( var1 - var2 ) ));
        return (var1._mod->options.INTERSECTION_USE?intersect(  sqr( (var1 + var2)*0.5 ) - sqr( (var1 - var2)*0.5 ) ,bnd1*bnd2)
                                                               :sqr( (var1 + var2)*0.5 ) - sqr( (var1 - var2)*0.5 ));
        break;
      }
    }
  }

  ISVar<T> var3( var2 );
  var3 *= var1;
  return (var1._mod->options.INTERSECTION_USE?intersect(var3,bnd1*bnd2):var3);
  //return intersect(var3,bnd1*bnd2);
}

template <typename T>
inline
ISVar<T> operator*
( ISVar<T> const& var1, double const& cst2 )
{
  if( !var1._mod ) 
    return( var1._cst * cst2 );
  ISVar<T> var3( var1 );
  var3 *= cst2;
  return var3;
}

template <typename T>
inline
ISVar<T> operator*
( ISVar<T> && var1, double const& cst2 )
{
  if( !var1._mod ) 
    return( var1._cst * cst2 );
  ISVar<T> var3( std::move(var1) );
  var3 *= cst2;
  return var3;
}



template <typename T>
inline
ISVar<T> operator*
( double const& cst1, ISVar<T> const& var2 )
{
  if( !var2._mod ) 
    return( cst1 * var2._cst );
  ISVar<T> var3( var2 );
  var3 *= cst1;
  return var3;
}


template <typename T>
inline
ISVar<T> operator*
( double const& cst1, ISVar<T> && var2 )
{
  if( !var2._mod ) 
    return( cst1 * var2._cst );
  ISVar<T> var3( std::move(var2) );
  var3 *= cst1;
  return var3;
}


template <typename T>
inline
ISVar<T>& ISVar<T>::operator/=
( double const& cst )
{
  if( cst == 0. )
    throw typename ISModel<T>::Exceptions( ISModel<T>::Exceptions::DIV );
  *this *= inv( cst );
  return *this;
}

template <typename T>
inline
ISVar<T>& ISVar<T>::operator/=
( ISVar<T> const& var )
{
  *this *= inv( var );
  return *this;
}

template <typename T>
inline
ISVar<T> operator/
( ISVar<T> const& var1, ISVar<T> const& var2 )
{
  ISVar<T> var3( var1 );
  var3 /= var2;
  return var3;
}

template <typename T>
inline
ISVar<T> operator/
( ISVar<T> const& var1, double const& cst2 )
{
  ISVar<T> var3( var1 );
  var3 /= cst2;
  return var3;
}

template <typename T>
inline
ISVar<T> operator/
( ISVar<T> && var1, double const& cst2 )
{
  ISVar<T> var3( std::move(var1) );
  var3 /= cst2;
  return var3;
}


template <typename T>
inline
ISVar<T> operator/
( double const& cst1, ISVar<T> const& var2 )
{
  ISVar<T> var3( inv( var2 ) );
  var3 *= cst1;
  return var3;
}

template <typename T>
inline
ISVar<T> operator/
( double const& cst1, ISVar<T> && var2 )
{
  ISVar<T> var3( inv( std::move(var2) ) );
  var3 *= cst1;
  return var3;
}


template <typename T>
inline
ISVar<T> inv
( ISVar<T> const& var )
{
  if( !var._mod )
    return 1./var._cst;

  ISVar<T> var2( var );
  if (var2._mod->options.SLOPE_USE ){
    var2._mod->_inv( var2._mat, var2._slope, var2._mod->_psize, var2._ndep );
    var2._bnd.second = false;
    return var2;
  }    
  var2._mod->_inv( var2._mat, var2._ndep );
  var2._bnd.second = false;
  return var2;
}

template <typename T>
inline
ISVar<T> inv
( ISVar<T> && varIn )
{
  if( !varIn._mod ){
    varIn._cst = 1./varIn._cst;
    return varIn;  
  }

  ISVar<T> var(std::move(varIn));

  if (var._mod->options.SLOPE_USE ){
    var._mod->_inv( var._mat, var._slope, var._mod->_psize, var._ndep );
    var._bnd.second = false;
    return var;
  }

  var._mod->_inv( var._mat, var._ndep );
  var._bnd.second = false;
  return var;
}

template <typename T>
inline
ISVar<T> sqr
( ISVar<T> const& var )
{
  if( !var._mod )
    return sqr(var._cst);

  ISVar<T> var2( var );
  if (var2._mod->options.SLOPE_USE ){
    var2._mod->_sqr( var2._mat, var2._slope, var2._mod->_psize, var2._ndep );
    var2._bnd.second = false;
    return var2;
  }  
  var2._mod->_sqr( var2._mat, var2._ndep );
  var2._bnd.second = false;
  return var2;
}

template <typename T>
inline
ISVar<T> sqr
( ISVar<T> && varIn )
{
  if( !varIn._mod ){
    varIn._cst = sqr(varIn._cst);
    return varIn;  
  }  

  ISVar<T> var(std::move(varIn));

  if (var._mod->options.SLOPE_USE ){
    var._mod->_sqr( var._mat, var._slope, var._mod->_psize, var._ndep );
    var._bnd.second = false;
    return var;
  }

  var._mod->_sqr( var._mat, var._ndep );
  var._bnd.second = false;
  return var;
}

template <typename T>
inline
ISVar<T> fabs
( ISVar<T> const& var )
{
  if( !var._mod ){
    double cstl = Op<T>::l( var._cst );
    double cstu = Op<T>::u( var._cst );
    if (cstl < 0){
      if (cstu > 0.){
        return T(0.,std::max(-cstl,cstu));
      }
      return -var._cst;
    }
    return var._cst;
  }

  ISVar<T> var2( var );
  auto const& f = [=]( const double& x ){ return std::fabs( x ); };
  if (var2._mod->options.SLOPE_USE ){
    auto const& fDerv = [=]( const double& x ){ 
      if (std::fabs(x) > MC__ISM_COMPUTATION_TOL) {
        if (x>0) return 1.;
        else return -1.;
      }
      else{
        return 0.; 
      }
    };
    var2._mod->_asym_slope( var2._mat, var2._slope, var2._mod->_psize, var2._ndep, f, fDerv, 0., true);
    var2._bnd.second = false;
    return var2;
  }  
  var2._mod->_asym( var2._mat, var2._ndep, f, 0., true ); // convex term, min 0
  var2._bnd.second = false;
  return var2;
}

template <typename T>
inline
ISVar<T> fabs
( ISVar<T> && varIn )
{

  ISVar<T> var(std::move(varIn));

  if( !var._mod ){
    double cstl = Op<T>::l( var._cst );
    double cstu = Op<T>::u( var._cst );
    if (cstl < 0){
      if (cstu > 0.){
        return T(0.,std::max(-cstl,cstu));
      }
      return -var;
    }
    return var;
  }

  // if( !varIn._mod ){
  //   varIn._cst = std::fabs(varIn._cst);
  //   return varIn;
  // }


  auto const& f = [=]( const double& x ){ return std::fabs( x ); };
  if (var._mod->options.SLOPE_USE ){
    auto const& fDerv = [=]( const double& x ){ 
      if (std::fabs(x) > MC__ISM_COMPUTATION_TOL) {
        if (x>0) return 1.;
        else return -1.;
      }
      else{
        return 0.; 
      }
    };               
    var._mod->_asym_slope( var._mat, var._slope, var._mod->_psize, var._ndep, f, fDerv, 0., true);    
    var._bnd.second = false;
    return var;
  }
  var._mod->_asym( var._mat, var._ndep, f, 0., true ); // convex term, min 0
  var._bnd.second = false;
  return var;
}

template <typename T>
inline
ISVar<T> relu
( ISVar<T> const& var )
{
  // if( !var._mod )
  //   //return relu(var._cst);
  //   return std::max(var._cst,0.);

  if( !var._mod ){
    double cstl = Op<T>::l( var._cst );
    double cstu = Op<T>::u( var._cst );
    return T(std::max(cstl,0.),std::max(cstu,0.));
  }


  T lazyBnd = var.bound();
  
  if(Op<T>::l(lazyBnd) > -MC__ISM_COMPUTATION_TOL){
    ISVar<T> var2( var );
    var2._bnd.second = false;
    return var2;
  } 
  
  if(Op<T>::u(lazyBnd) < MC__ISM_COMPUTATION_TOL){
    ISVar<T> var2( 0. );
    return var2;
  }  


  ISVar<T> var2( var );
  //auto const& f = [=]( const double& x ){ return std::max( x, 0. ); };
  if (var2._mod->options.SLOPE_USE ){ 
     //auto const& fDerv = [=]( const double& x ){ return x > 0? 1.: 0.; };  
     //var2._mod->_asym_slope( var2._mat, var2._slope, var2._mod->_psize, var2._ndep, f,fDerv, -DBL_MAX, true );
     //var2._mod->_asym_slope_relu( var2._mat, var2._slope, var2._mod->_psize, var2._ndep);
    if(var2._mod->options.SHADOW_USE && var2._ndep>1){
#ifndef NOTTOTRACKSHADOW     
      std::cout << "_asym_slope_relu_shadow" << std::endl;
#endif
      var2._mod->_asym_slope_relu_shadow( var2._mat, var2._slope, var2._mod->_psize, var2._ndep,var2._shadow,var2._shadow_slope,var2._shadow_info);
      if(std::fabs(var2._shadow_info[0]) > 100.){
        return var2*=0;
      }
    }
    else{
      var2._mod->_asym_slope_relu( var2._mat, var2._slope, var2._mod->_psize, var2._ndep);
    }
    var2._bnd.second = false;
    return var2;
  }
  auto const& f = [=]( const double& x ){ return std::max( x, 0. ); };
  var2._mod->_asym( var2._mat, var2._ndep, f, 0., true ); // convex term, min 0
  var2._bnd.second = false;
  return var2;
}

template <typename T>
inline
ISVar<T> relu
( ISVar<T> && varIn )
{

  ISVar<T> var(std::move(varIn));

  if( !var._mod ){
    double cstl = Op<T>::l( var._cst );
    double cstu = Op<T>::u( var._cst );
    var._cst = std::move(T(std::max(cstl,0.),std::max(cstu,0.)));
    return var;
  }


  // if( !varIn._mod ){
  //   //return relu(var._cst);
  //   varIn._cst = std::max(varIn._cst,0.);
  //   return varIn;
  // }


  T lazyBnd = var.bound();
  
  if(Op<T>::l(lazyBnd)> -MC__ISM_COMPUTATION_TOL){
    var._bnd.second = false;
    return var;
  } 
  
  if(Op<T>::u(lazyBnd) < MC__ISM_COMPUTATION_TOL){
    var = 0.;
    return var;
  }  

  //auto const& f = [=]( const double& x ){ return std::max( x, 0. ); };
  if (var._mod->options.SLOPE_USE){
    //auto const& fDerv = [=]( const double& x ){ return x > 0? 1.: 0.; }; 
    //var._mod->_asym_slope( var._mat, var._slope, var._mod->_psize, var._ndep, f,fDerv, -DBL_MAX, true );
    if(var._mod->options.SHADOW_USE){
#ifndef NOTTOTRACKSHADOW     
      std::cout << "_asym_slope_relu_shadow" << std::endl;
#endif
      var._mod->_asym_slope_relu_shadow( var._mat, var._slope, var._mod->_psize, var._ndep,var._shadow,var._shadow_slope,var._shadow_info);
      if(std::fabs(var._shadow_info[0]) > 100.){
        return var*=0;
      }
    }
    else{
      var._mod->_asym_slope_relu( var._mat, var._slope, var._mod->_psize, var._ndep);
    }
    var._bnd.second = false;
    return var;
  }
  auto const& f = [=]( const double& x ){ return std::max( x, 0. ); };
  var._mod->_asym( var._mat, var._ndep, f, 0., true ); // convex term, min 0
  var._bnd.second = false;
  return var;
}



template <typename T>
inline
ISVar<T> sqrt
( ISVar<T> const& var )
{
  if( !var._mod )
    return sqrt(var._cst);

  ISVar<T> var2( var );
  if (var2._mod->options.SLOPE_USE ){
    var2._mod->_sqrt( var2._mat, var2._slope, var2._mod->_psize, var2._ndep );
    var2._bnd.second = false;
    return var2;
  }    
  var2._mod->_sqrt( var2._mat, var2._ndep );
  var2._bnd.second = false;
  return var2;
}

template <typename T>
inline
ISVar<T> sqrt
( ISVar<T> && varIn )
{
  if( !varIn._mod ){
    varIn._cst = sqrt(varIn._cst);
    return varIn;
  }

  ISVar<T> var(std::move(varIn));

  if (var._mod->options.SLOPE_USE ){
    var._mod->_sqrt( var._mat, var._slope, var._mod->_psize, var._ndep );
    var._bnd.second = false;
    return var;
  }

  var._mod->_sqrt( var._mat, var._ndep );
  var._bnd.second = false;
  return var;
}

template <typename T>
inline
ISVar<T> exp
( ISVar<T> const& var )
{
  if( !var._mod )
    return exp(var._cst);

  ISVar<T> var2( var );
  if (var2._mod->options.SLOPE_USE ){
    var2._mod->_exp( var2._mat, var2._slope, var2._mod->_psize, var2._ndep );
    var2._bnd.second = false;
    return var2;
  }
  var2._mod->_exp( var2._mat, var2._ndep );
  var2._bnd.second = false;
  return var2;
}

template <typename T>
inline
ISVar<T> exp
( ISVar<T> && varIn )
{
  if( !varIn._mod ){
    varIn._cst = exp(varIn._cst);
    return varIn;
  }

  ISVar<T> var(std::move(varIn));

  if (var._mod->options.SLOPE_USE ){
    var._mod->_exp( var._mat, var._slope, var._mod->_psize, var._ndep );
    var._bnd.second = false;
    return var;
  }

  var._mod->_exp( var._mat, var._ndep );
  var._bnd.second = false;
  return var;
}

template <typename T>
inline
ISVar<T> log
( ISVar<T> const& var )
{
  if( !var._mod )
    return log(var._cst);

  ISVar<T> var2( var );
  if (var2._mod->options.SLOPE_USE ){
    var2._mod->_log( var2._mat, var2._slope, var2._mod->_psize, var2._ndep );
    var2._bnd.second = false;
    return var2;
  }    
  var2._mod->_log( var2._mat, var2._ndep );
  var2._bnd.second = false;
  return var2;
}

template <typename T>
inline
ISVar<T> log
( ISVar<T> && varIn )
{
  if( !varIn._mod ){
    varIn._cst = log(varIn._cst);
    return varIn;
  }

  ISVar<T> var(std::move(varIn));

  if (var._mod->options.SLOPE_USE ){
    var._mod->_log( var._mat, var._slope, var._mod->_psize, var._ndep );
    var._bnd.second = false;
    return var;
  }

  var._mod->_log( var._mat, var._ndep );
  var._bnd.second = false;
  return var;
}

template <typename T>
inline
ISVar<T> xlog
( ISVar<T> const& var )
{
  if( !var._mod ){
    double cstl = Op<T>::l( var._cst );
    double cstu = Op<T>::u( var._cst );
    if(cstl <= 0.)
      typename mc::ISModel<T>::Exceptions( mc::ISModel<T>::Exceptions::LOG ); 
    else{
      if (cstl < std::exp(-1.)){
        if (cstu > std::exp(-1.))
          return T(-std::exp(-1.),std::max(xlog(cstl),xlog(cstl)));
        else
          return T(xlog(cstu),xlog(cstl));
      }
      else 
        return T(xlog(cstl),xlog(cstu));
    }
  }

  ISVar<T> var2( var );

  if (var2._mod->options.SLOPE_USE ){
    var2._mod->_log( var2._mat, var2._slope, var2._mod->_psize, var2._ndep );
    var2._bnd.second = false;
    return var2;
  }    

  var2._mod->_xlog( var2._mat, var2._ndep );
  var2._bnd.second = false;
  return var2;
}

template <typename T>
inline
ISVar<T> xlog
( ISVar<T> && varIn )
{
 
  if( !varIn._mod ){
    double cstl = Op<T>::l( varIn._cst );
    double cstu = Op<T>::u( varIn._cst );
    if(cstl <= 0.)
      throw typename mc::ISModel<T>::Exceptions( mc::ISModel<T>::Exceptions::LOG ); 
    else{
      if (cstl < std::exp(-1.)){
        if (cstu > std::exp(-1.))
          varIn._cst = T(-std::exp(-1.),std::max(xlog(cstl),xlog(cstl)));
        else
          varIn._cst = T(xlog(cstu),xlog(cstl));
      }
      else 
        varIn._cst = T(xlog(cstl),xlog(cstu));
    }
    return varIn;
  }



  ISVar<T> var(std::move(varIn));

  if (var._mod->options.SLOPE_USE ){
    var._mod->_xlog( var._mat, var._slope, var._mod->_psize, var._ndep );
    var._bnd.second = false;
    return var;
  }

  var._mod->_xlog( var._mat, var._ndep );
  var._bnd.second = false;
  return var;
}

template <typename T>
inline
ISVar<T> sin
( ISVar<T> const& var )
{
  if( !var._mod )
    return sin(var._cst);

  if (var._mod->options.SLOPE_USE && var._ndep == 1){
    return cos(var - PI/2.);
  }

  ISVar<T> var2( var );
  var2._mod->_sin( var2._mat, var2._ndep );
  var2._bnd.second = false;
  return var2;
}

template <typename T>
inline
ISVar<T> sin
( ISVar<T> && varIn )
{
  if( !varIn._mod ){
    varIn._cst = sin(varIn._cst);
    return varIn;
  }

  ISVar<T> var(std::move(varIn));

  if (var._mod->options.SLOPE_USE && var._ndep == 1){
    return cos(var - PI/2.);
  }

  var._mod->_sin( var._mat, var._ndep );
  var._bnd.second = false;
  return var;
}

template <typename T>
inline
ISVar<T> cos
( ISVar<T> const& var )
{
  if( !var._mod )
    return cos(var._cst);

  ISVar<T> var2( var );
  if (var2._mod->options.SLOPE_USE && var2._ndep == 1){
    var2._mod->_cos_slope( var2._mat, var2._slope, var2._mod->_psize, var2._ndep );
    var2._bnd.second = false;
    return var2;
  }  

  var2._mod->_cos( var2._mat, var2._ndep );
  var2._bnd.second = false;
  return var2;
}

template <typename T>
inline
ISVar<T> cos
( ISVar<T> && varIn )
{
  if( !varIn._mod ){
    varIn._cst = cos(varIn._cst);
    return varIn;
  }

  ISVar<T> var(std::move(varIn));

  if (var._mod->options.SLOPE_USE && var._ndep == 1){
    var._mod->_cos_slope( var._mat, var._slope, var._mod->_psize, var._ndep );
    var._bnd.second = false;
    return var;
  }

  var._mod->_cos( var._mat, var._ndep );
  var._bnd.second = false;
  return var;
}

template <typename T>
inline
ISVar<T> tanh
( ISVar<T> const& var )
{
  if( !var._mod )
    return tanh(var._cst);
#ifdef TEST_MOVE
std::cout << "Tanh Copy" << std::endl;
#endif
  if( !var._mod->options.DCDEC_USE )
    return 1-2/(exp(2*var)+1);

  ISVar<T> var2( var );
  if (var2._mod->options.SLOPE_USE ){
    var2._mod->_tanh( var2._mat, var2._slope, var2._mod->_psize, var2._ndep );
    var2._bnd.second = false;
    return var2;
  }    
  var2._mod->_tanh( var2._mat, var2._ndep );
  var2._bnd.second = false;
  return var2;

}

template <typename T>
inline
ISVar<T> tanh
( ISVar<T> && varIn )
{
  if( !varIn._mod ){
    varIn._cst = tanh(varIn._cst);
    return varIn;//std::tanh(var._cst);
  }
#ifdef TEST_MOVE
std::cout << "Tanh Move" << std::endl;
#endif
  ISVar<T> var(std::move(varIn));

  if( !var._mod->options.DCDEC_USE )
    return 1-2/(exp(2*var)+1);

  if (var._mod->options.SLOPE_USE ){
    var._mod->_tanh( var._mat, var._slope, var._mod->_psize, var._ndep );
    var._bnd.second = false;
    return var;
  }

  var._mod->_tanh( var._mat, var._ndep );
  var._bnd.second = false;
  return var;

}

template <typename T>
inline
ISVar<T>
pow
( ISVar<T> const& var, int const& n )
{
 
  if( !var._mod )
    return Op<T>::pow( var._cst, n );

  std::cout << "Pow Copy" << std::endl;  

  if( n == 0 )
    return 1.;
  if( n == 1 )
    return var;


  if( !var._mod->options.DCDEC_USE ){
    if( n < 0 )
      return pow( inv( var ), -n );
    return // recursive call
      !(n%2) ?
      sqr( pow( var, n/2 ) ) :
      sqr( pow( var, n/2 ) ) * var;
  }

  ISVar<T> var2( var );
  if(var2._mod->options.SLOPE_USE){
    var2._mod->_pow( var2._mat, var2._slope, var2._mod->_psize, n, var2._ndep );
  }
  else{  
    var2._mod->_pow( var2._mat, n, var2._ndep );
  }  
  var2._bnd.second = false;
  return var2;
}


template <typename T>
inline
ISVar<T>
pow
( ISVar<T> && varIn, int const& n )
{

  if( !varIn._mod ){
    varIn._cst = Op<T>::pow( varIn._cst, n );
    return varIn;
  }
    
  if( n == 0 )
    return 1.;

  //std::cout << "Pow Move2" << std::endl;  
  ISVar<T> var(std::move(varIn));

  if( n == 1 )
    return var;

  if( !var._mod->options.DCDEC_USE ){
    if( n < 0 )
      return pow( inv( var ), -n );
    return // recursive call
      !(n%2) ?
      sqr( pow( var, n/2 ) ) :
      sqr( pow( var, n/2 ) ) * var;
  }
  
  if(var._mod->options.SLOPE_USE){
    var._mod->_pow( var._mat,var._slope, var._mod->_psize, n, var._ndep );
  }
  else{
    var._mod->_pow( var._mat, n, var._ndep );
  }  
  var._bnd.second = false;
  return var;
}

template <typename T>
inline
ISVar<T>
pow
( ISVar<T> const& var, double const& d )
{
  return exp( d * log( var ) );
}

template <typename T>
inline
ISVar<T>
pow
( ISVar<T> && var, double const& d )
{
  return exp( d * log( var ) );
}

template <typename T>
inline
ISVar<T>
pow
( ISVar<T> const& var, ISVar<T> const& exp )
{
  return exp( exp * log( var ) );
}

template <typename T>
inline
ISVar<T>
pow
( ISVar<T> && var, ISVar<T> const& exp )
{
  return exp( exp * log( var ) );
}

template <typename T>
inline
ISVar<T>
pow
( double const& d, ISVar<T> const& var )
{
  return exp( var * std::log( d ) );
}

template <typename T>
inline
ISVar<T>
pow
( double const& d, ISVar<T> && var )
{
  return exp( var * std::log( d ) );
}

template <typename T>
inline
ISVar<T>
prod
( unsigned int const& n, ISVar<T> const* var )
{
  switch( n ){
   case 0:  return 1.;
   case 1:  return var[0];
   default: return var[0] * prod( n-1, var+1 );
  }
}


template <typename T>
inline
ISVar<T> intersect
( ISVar<T> const& var , T _rangebnd )
{
  if( !var._mod ){
    double cstl = Op<T>::l( var._cst );
    double cstu = Op<T>::u( var._cst );
    double rgbndl = Op<T>::l( _rangebnd );
    double rgbndu = Op<T>::u( _rangebnd );
    return T(std::max(cstl,rgbndl),std::min(cstu,rgbndu));
  }
  

  ISVar<T> var2( var );
  if (var2._mod->options.SLOPE_USE ){
    var2._mod->_intersect( var2._mat, var2._slope, var2._mod->_psize, var2._ndep , _rangebnd);
    var2._bnd.second = false;
    return var2;
  }        
  var2._mod->_intersect( var2._mat, var2._ndep , _rangebnd);
  var2._bnd.second = false;
  return var2;
}

template <typename T>
inline
ISVar<T> intersect
( ISVar<T> && varIn , T _rangebnd )
{

  ISVar<T> var(std::move(varIn));

  if( !var._mod ){
    double cstl = Op<T>::l( var._cst );
    double cstu = Op<T>::u( var._cst );
    double rgbndl = Op<T>::l( _rangebnd );
    double rgbndu = Op<T>::u( _rangebnd );
    var._cst = std::move(T(std::max(cstl,rgbndl),std::min(cstu,rgbndu)));
    return var;
  }

  if (var._mod->options.SLOPE_USE ){
    var._mod->_intersect( var._mat, var._slope, var._mod->_psize, var._ndep , _rangebnd);
    var._bnd.second = false;
    return var;
  }

  var._mod->_intersect( var._mat, var._ndep , _rangebnd);
  var._bnd.second = false;
  return var;
}

template <typename T>
inline
ISVar<T>
monom
( unsigned int const& n, ISVar<T> const* var, unsigned int const* exp )
{
  switch( n ){
   case 0:  return 1.;
   case 1:  return pow( var[0], (int)exp[0] );
   default: return pow( var[0], (int)exp[0] ) * monom( n-1, var+1, exp+1 );
  }
}

template <typename T>
inline
ISVar<T>
cheb
( ISVar<T> const& var, unsigned int const& n )
{
  if( !var._mod )
    return T(-1.,1.);

  switch( n ){
    case 0:  return 1.;
    case 1:  return var;
    default: break;
  }
  ISVar<T> var2( 2.*var*cheb(var,n-1) - cheb(var,n-2) );
  return var2;
}

template <typename T>
inline
ISVar<T>
cheb
( ISVar<T> && varIn, unsigned int const& n )
{
  if( !varIn._mod )
    return T(-1.,1.);

  ISVar<T> var(std::move(varIn));

  switch( n ){
    case 0:  return 1.;
    case 1:  return var;
    default: break;
  }
  ISVar<T> var2( 2.*var*cheb(var,n-1) - cheb(var,n-2) );
  return var2;
}

template <typename T>
inline
ISVar<T>
max
( ISVar<T> const& var1, ISVar<T> const& var2 )
{
  if( !var1._mod && !var2._mod ){ 
    double cstl = Op<T>::l( var1._cst );
    double cstu = Op<T>::u( var1._cst );
    double rgbndl = Op<T>::l( var2._cst );
    double rgbndu = Op<T>::u( var2._cst );
    return T(std::max(cstl,rgbndl),std::max(cstu,rgbndu));
  }
  if( !var2._mod ){
    double cstl = Op<T>::l( var2._cst );
    double cstu = Op<T>::u( var2._cst );
    if (cstu - cstl < MC__ISM_COMPUTATION_TOL )    
      return max( var1, cstl );
    else
      throw typename mc::ISModel<T>::Exceptions( mc::ISModel<T>::Exceptions::UNDEF ); 
  } 
  if( !var1._mod ){
    double cstl = Op<T>::l( var1._cst );
    double cstu = Op<T>::u( var1._cst );
    if (cstu - cstl < MC__ISM_COMPUTATION_TOL )    
      return max( var2, cstl );
    else
      throw typename mc::ISModel<T>::Exceptions( mc::ISModel<T>::Exceptions::UNDEF ); 
  }      
  return max( var1 - var2, 0. ) + var2;
}

template <typename T>
inline
ISVar<T>
max
( ISVar<T> const& var1, double const& cst2 )
{

  if(isequal(cst2,0.)){
    return relu(var1);
  }

  if( !var1._mod ){
    double cstl = Op<T>::l( var1._cst );
    double cstu = Op<T>::u( var1._cst );
    return T(std::max(cstl,cst2),std::max(cstu,cst2));
  } 
  T lazyBnd = var1.bound();
  
  if(Op<T>::l(lazyBnd)> cst2){
    ISVar<T> var2( var1 );
    var2._bnd.second = false;
    return var2;
  } 

  if(Op<T>::u(lazyBnd) < cst2){
    ISVar<T> var2( cst2 );
    return var2;
  }  

  ISVar<T> var2( var1 );
  auto const& f = [=]( const double& x ){return std::max( x, cst2 ); };//Op<T>::max( x, cst2 ); }; // why not  return std::max( x, cst2 ); };
  if (var2._mod->options.SLOPE_USE ){ 
    auto const& fDerv = [=]( const double& x ){ return x > cst2? 1.: 0.; };  
    var2._mod->_asym_slope( var2._mat, var2._slope, var2._mod->_psize, var2._ndep, f,fDerv, -DBL_MAX, true );
    var2._bnd.second = false;
    return var2;
  }        
  var2._mod->_asym( var2._mat, var2._ndep, f, cst2, true ); // convex term, min cst2
  var2._bnd.second = false;
  return var2;
}

template <typename T>
inline
ISVar<T>
max
( ISVar<T> && varIn, double const& cst2 )
{

  if( !varIn._mod ){
    double cstl = Op<T>::l( varIn._cst );
    double cstu = Op<T>::u( varIn._cst );
    varIn._cst = T(std::max(cstl,cst2),std::max(cstu,cst2));
    return varIn;
  } 

  ISVar<T> var1(std::move(varIn));


  if(isequal(cst2,0.)){
    return relu(var1);
  }

  T lazyBnd = var1.bound();
  
  if(Op<T>::l(lazyBnd)> cst2){
    var1._bnd.second = false;
    return var1;
  } 
  
  if(Op<T>::u(lazyBnd) < cst2){
    var1 = cst2;
    return var1;
  }  

  auto const& f = [=]( const double& x ){ return std::max( x, cst2 ); };
  if (var1._mod->options.SLOPE_USE){
    auto const& fDerv = [=]( const double& x ){ return x > cst2? 1.: 0.; }; 
    var1._mod->_asym_slope( var1._mat, var1._slope, var1._mod->_psize, var1._ndep, f,fDerv, -DBL_MAX, true );
    var1._bnd.second = false;
    return var1;
  }  
  var1._mod->_asym( var1._mat, var1._ndep, f, cst2, true ); // convex term, min cst2
  var1._bnd.second = false;
  return var1;
}

template <typename T>
inline
ISVar<T>
max
( double const& cst1, ISVar<T> const& var2 )
{
  return max( var2, cst1 );
}

template <typename T>
inline
ISVar<T>
max
( double const& cst1, ISVar<T> && var2 )
{
  return max( var2, cst1 );
}

template <typename T>
inline
ISVar<T>
min
( ISVar<T> const& var1, ISVar<T> const& var2 )
{
  if( !var1._mod && !var2._mod ){ 
    double cstl = Op<T>::l( var1._cst );
    double cstu = Op<T>::u( var1._cst );
    double rgbndl = Op<T>::l( var2._cst );
    double rgbndu = Op<T>::u( var2._cst );
    return T(std::min(cstl,rgbndl),std::min(cstu,rgbndu));
  }
  if( !var2._mod ){
    double cstl = Op<T>::l( var2._cst );
    double cstu = Op<T>::u( var2._cst );
    if (cstu - cstl < MC__ISM_COMPUTATION_TOL )    
      return min( var1, cstu );
    else
      throw typename mc::ISModel<T>::Exceptions( mc::ISModel<T>::Exceptions::UNDEF ); 
  } 
  if( !var1._mod ){
    double cstl = Op<T>::l( var1._cst );
    double cstu = Op<T>::u( var1._cst );
    if (cstu - cstl < MC__ISM_COMPUTATION_TOL )    
      return min( var2, cstu );
    else
      throw typename mc::ISModel<T>::Exceptions( mc::ISModel<T>::Exceptions::UNDEF ); 
  }      
  return min( var1 - var2, 0. ) + var2;

}

template <typename T>
inline
ISVar<T>
min
( ISVar<T> const& var1, double const& cst2 )
{

  if( !var1._mod ){
    double cstl = Op<T>::l( var1._cst );
    double cstu = Op<T>::u( var1._cst );
    return T(std::min(cstl,cst2),std::min(cstu,cst2));
  } 

  ISVar<T> var2( var1 );
  auto const& f = [=]( const double& x ){ return std::min( x, cst2 ); };
  if (var2._mod->options.SLOPE_USE ){ 
    auto const& fDerv = [=]( const double& x ){ return x < cst2? 1.: 0.; };  
    var2._mod->_asym_slope( var2._mat, var2._slope, var2._mod->_psize, var2._ndep, f,fDerv, DBL_MAX, false );
    var2._bnd.second = false;
    return var2;
  }      
  var2._mod->_asym( var2._mat, var2._ndep, f, cst2, false ); // convex term, max cst2
  var2._bnd.second = false;
  return var2;
}

template <typename T>
inline
ISVar<T>
min
( ISVar<T> && varIn, double const& cst2 )
{
  if( !varIn._mod ){
    double cstl = Op<T>::l( varIn._cst );
    double cstu = Op<T>::u( varIn._cst );
    varIn._cst = T(std::min(cstl,cst2),std::min(cstu,cst2));
    return varIn;
  } 

  ISVar<T> var1(std::move(varIn));

  auto const& f = [=]( const double& x ){return std::min( x, cst2 ); };//Op<T>::min( x, cst2 ); }; // why not return std::min( x, cst2 ); };
  if (var1._mod->options.SLOPE_USE){
    auto const& fDerv = [=]( const double& x ){ return x < cst2? 1.: 0.; }; 
    var1._mod->_asym_slope( var1._mat, var1._slope, var1._mod->_psize, var1._ndep, f,fDerv, DBL_MAX, false );
    var1._bnd.second = false;
    return var1;
  }  
  var1._mod->_asym( var1._mat, var1._ndep, f, cst2, false ); // concave term, max cst2
  var1._bnd.second = false;
  return var1;
}

template <typename T>
inline
ISVar<T>
min
( double const& cst1, ISVar<T> const& var2 )
{
  return min( var2, cst1 );
}

template <typename T>
inline
ISVar<T>
min
( double const& cst1, ISVar<T> && var2 )
{
  return min( var2, cst1 );
}

}//namespace mc

#include "mcfadbad.hpp"

namespace fadbad
{

//! @brief Specialization of the structure fadbad::Op for use of the type mc::ISVar of MC++ as a template parameter of the classes fadbad::F, fadbad::B and fadbad::T of FADBAD++
template< typename T > struct Op< mc::ISVar<T> >
{ 
  typedef mc::ISVar<T> ISV;
  typedef double Base;
  static Base myInteger( const int i ) { return Base(i); }
  static Base myZero() { return myInteger(0); }
  static Base myOne() { return myInteger(1);}
  static Base myTwo() { return myInteger(2); }
  static double myPI() { return mc::PI; }
  static ISV myPos( const ISV& x ) { return  x; }
  static ISV myNeg( const ISV& x ) { return -x; }
  template <typename U> static ISV& myCadd( ISV& x, const U& y ) { return x+=y; }
  template <typename U> static ISV& myCsub( ISV& x, const U& y ) { return x-=y; }
  template <typename U> static ISV& myCmul( ISV& x, const U& y ) { return x*=y; }
  template <typename U> static ISV& myCdiv( ISV& x, const U& y ) { return x/=y; }
  static ISV myInv( const ISV& x ) { return mc::inv( x ); }
  static ISV mySqr( const ISV& x ) { return mc::sqr( x ); }
  template <typename X, typename Y> static ISV myPow( const X& x, const Y& y ) { return mc::pow( x, y ); }
  //static ISV myCheb( const ISV& x, const unsigned n ) { return mc::cheb( x, n ); }
  static ISV mySqrt( const ISV& x ) { return mc::sqrt(x); }
  static ISV myLog( const ISV& x ) { return mc::log( x ); }
  static ISV myExp( const ISV& x ) { return mc::exp( x ); }
  static ISV mySin( const ISV& x ) { return mc::sin( x ); }
  static ISV myCos( const ISV& x ) { return mc::cos( x ); }
  static ISV myTan( const ISV& x ) { throw typename mc::ISModel<T>::Exceptions( mc::ISModel<T>::Exceptions::UNDEF ); } //{ return mc::tan( x ); }
  static ISV myAsin( const ISV& x ) { throw typename mc::ISModel<T>::Exceptions( mc::ISModel<T>::Exceptions::UNDEF ); }
  static ISV myAcos( const ISV& x ) { throw typename mc::ISModel<T>::Exceptions( mc::ISModel<T>::Exceptions::UNDEF ); }
  static ISV myAtan( const ISV& x ) { throw typename mc::ISModel<T>::Exceptions( mc::ISModel<T>::Exceptions::UNDEF ); }
  static ISV mySinh( const ISV& x ) { throw typename mc::ISModel<T>::Exceptions( mc::ISModel<T>::Exceptions::UNDEF ); }
  static ISV myCosh( const ISV& x ) { throw typename mc::ISModel<T>::Exceptions( mc::ISModel<T>::Exceptions::UNDEF ); }
  static ISV myTanh( const ISV& x ) { throw typename mc::ISModel<T>::Exceptions( mc::ISModel<T>::Exceptions::UNDEF ); }
  static bool myEq( const ISV& x, const ISV& y ) { return mc::Op<T>::eq(const_cast<ISV*>(&x)->bound(),const_cast<ISV*>(&y)->bound()); } 
  static bool myNe( const ISV& x, const ISV& y ) { return mc::Op<T>::ne(const_cast<ISV*>(&x)->bound(),const_cast<ISV*>(&y)->bound()); }
  static bool myLt( const ISV& x, const ISV& y ) { return mc::Op<T>::lt(const_cast<ISV*>(&x)->bound(),const_cast<ISV*>(&y)->bound()); }
  static bool myLe( const ISV& x, const ISV& y ) { return mc::Op<T>::le(const_cast<ISV*>(&x)->bound(),const_cast<ISV*>(&y)->bound()); }
  static bool myGt( const ISV& x, const ISV& y ) { return mc::Op<T>::gt(const_cast<ISV*>(&x)->bound(),const_cast<ISV*>(&y)->bound()); }
  static bool myGe( const ISV& x, const ISV& y ) { return mc::Op<T>::ge(const_cast<ISV*>(&x)->bound(),const_cast<ISV*>(&y)->bound()); }
};

} // end namespace fadbad

namespace mc
{

//! @brief C++ structure for specialization of the mc::Op templated structure for use of mc::ISVar in DAG evaluation and as template parameter in other MC++ types
template< typename T > struct Op< mc::ISVar<T> >
{
  typedef mc::ISVar<T> ISV;
  static ISV point( const double c ) { return ISV(c); }
  static ISV zeroone() { return ISV( mc::Op<T>::zeroone() ); }
  static void I(ISV& x, const ISV&y) { x = y; }
  static double l(const ISV& x) { return mc::Op<T>::l(x.B()); }
  static double u(const ISV& x) { return mc::Op<T>::u(x.B()); }
  static double abs (const ISV& x) { return mc::Op<T>::abs(x.B());  }
  static double mid (const ISV& x) { return mc::Op<T>::mid(x.B());  }
  static double diam(const ISV& x) { return mc::Op<T>::diam(x.B()); }
  static ISV inv (const ISV& x) { return mc::inv(x);  }
  static ISV sqr (const ISV& x) { return mc::sqr(x);  }
  static ISV sqrt(const ISV& x) { return mc::sqrt(x); }
  static ISV exp (const ISV& x) { return mc::exp(x);  }
  static ISV log (const ISV& x) { return mc::log(x);  }
  static ISV xlog(const ISV& x) { return mc::xlog(x); }
  static ISV lmtd(const ISV& x, const ISV& y) { return (x-y)/(mc::log(x)-mc::log(y)); }
  static ISV rlmtd(const ISV& x, const ISV& y) { return (mc::log(x)-mc::log(y))/(x-y); }
  static ISV fabs(const ISV& x) { return mc::fabs(x); }
  static ISV sin (const ISV& x) { return mc::sin(x);  }
  static ISV cos (const ISV& x) { return mc::cos(x);  }
  static ISV tan (const ISV& x) { throw typename ISModel<T>::Exceptions( ISModel<T>::Exceptions::UNDEF ); } //{ return mc::tan(x);  }
  static ISV asin(const ISV& x) { throw typename ISModel<T>::Exceptions( ISModel<T>::Exceptions::UNDEF ); }
  static ISV acos(const ISV& x) { throw typename ISModel<T>::Exceptions( ISModel<T>::Exceptions::UNDEF ); }
  static ISV atan(const ISV& x) { throw typename ISModel<T>::Exceptions( ISModel<T>::Exceptions::UNDEF ); }
  static ISV sinh(const ISV& x) { throw typename ISModel<T>::Exceptions( ISModel<T>::Exceptions::UNDEF ); }
  static ISV cosh(const ISV& x) { throw typename ISModel<T>::Exceptions( ISModel<T>::Exceptions::UNDEF ); }
  static ISV tanh(const ISV& x) { return mc::tanh(x); } // {ISV y(std::move(x)); return mc::tanh(std::move(y));} //
  //static ISV tanh(ISV && x) { ISV y(std::move(x)); return mc::tanh(std::move(y)); }  
  static ISV erf (const ISV& x) { throw typename ISModel<T>::Exceptions( ISModel<T>::Exceptions::UNDEF ); }
  static ISV erfc(const ISV& x) { throw typename ISModel<T>::Exceptions( ISModel<T>::Exceptions::UNDEF ); }
  static ISV fstep(const ISV& x) { throw typename ISModel<T>::Exceptions( ISModel<T>::Exceptions::UNDEF ); }
  static ISV bstep(const ISV& x) { throw typename ISModel<T>::Exceptions( ISModel<T>::Exceptions::UNDEF ); }
  static ISV hull(const ISV& x, const ISV& y) { throw typename ISModel<T>::Exceptions( ISModel<T>::Exceptions::UNDEF ); } //{ return mc::hull(x,y); }
  static ISV min (const ISV& x, const ISV& y) { return mc::min(x,y); }
  static ISV max (const ISV& x, const ISV& y) { return mc::max(x,y); }
  static ISV arh (const ISV& x, const double k) { return mc::exp(-k/x); }
  template <typename X, typename Y> static ISV pow(const X& x, const Y& y) { return mc::pow(x,y); }
  static ISV cheb(const ISV& x, const unsigned n) { return mc::cheb(x,n); }
  static ISV prod (const unsigned n, const ISV* x) { return mc::prod(n,x); }
  static ISV monom (const unsigned n, const ISV* x, const unsigned* k) { return mc::monom(n,x,k); }
  static bool inter(ISV& xIy, const ISV& x, const ISV& y) { throw typename ISModel<T>::Exceptions( ISModel<T>::Exceptions::UNDEF ); } //{ return mc::inter(xIy,x,y); }
  static bool eq(const ISV& x, const ISV& y) { return mc::Op<T>::eq(x.B(),y.B()); }
  static bool ne(const ISV& x, const ISV& y) { return mc::Op<T>::ne(x.B(),y.B()); }
  static bool lt(const ISV& x, const ISV& y) { return mc::Op<T>::lt(x.B(),y.B()); }
  static bool le(const ISV& x, const ISV& y) { return mc::Op<T>::le(x.B(),y.B()); }
  static bool gt(const ISV& x, const ISV& y) { return mc::Op<T>::gt(x.B(),y.B()); }
  static bool ge(const ISV& x, const ISV& y) { return mc::Op<T>::ge(x.B(),y.B()); }
};

} // namespace mc

#endif
