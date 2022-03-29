// Copyright (C) Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the Eclipse Public License.

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

 private:

  //! @brief Number of partitions
  unsigned int _ndiv;
  //! @brief Number of variables
  unsigned int _nvar;
  //! @brief Whether variables are defined or not
  std::vector<bool> _defvar;
  //! @brief Variable bounds
  std::vector<T> _bndvar;
  
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
  }

  ~ISModel() 
  {}

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


 //! @brief Options of mc::ISModel
  static struct Options
  {
    //! @brief Constructor
    Options():
      ASYREM_USE(true)
      {}
    //! @brief Whether to use asymmetric inclusions for convex/concave terms where possible
    bool ASYREM_USE;
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
  void _prod
  ( std::vector<std::vector<T>> const& mat1, std::vector<std::vector<T>> const& mat2,
    std::vector<std::vector<T>>& mat3, unsigned& ndep3 )
  const;
  // @unsure: the algorithm for computing asymmetric over-/under-estimators is applicable for sin and cos 
  // when they are concave or convex on the range of the input ISM  


  // @unsure: there is a similar algorithm which can refining the multiplicative result, 
  // by bounding x_iy_j+x_jy_i in 0.5*(x_i^2+y_i^2)*[-1,1] + 0.5*(x_j^2+y_j^2)*[-1,1]

  std::ostream& _dispvar
  ( std::vector<std::vector<T>> const& mat, unsigned const& ndep, const int& opt=0,
    std::ostream& out=std::cout )
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
  double C1( 0. );
  for( unsigned int i=0; i<_nvar; i++ ){
    if( mat[i].empty() ) continue;
    switch( imid ){
      case ICONV: _c1[i] = _L1[i], C1 += _c1[i]; break;
      case ICONC: _c1[i] = _U1[i], C1 += _c1[i]; break;
      case ICUT:  _c1[i] = 0.5*(_L1[i]+_U1[i]); C1 += _c1[i]; break;
    }
    _r1[i] = ( _U1[i] - _L1[i] );
    sum_r1 += _r1[i];
  }

  double inv_ndep( 1. / ndep );
  T fopt = f( imid == ICUT? zopt: C1 );
  T fopt_over_ndep = fopt * inv_ndep;
  
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
        T D = scal_r1 * ( f( ( mat[i][j] - _c1[i] ) / scal_r1 + C1 ) - fopt ) + fopt_over_ndep;
        T E = ( imid == ICUT? fopt_over_ndep: f( mat[i][j] - _c1[i] + C1 ) - fopt_over_ndep * (ndep-1.) );
        mat[i][j] = ( cvx? T( Op<T>::l(E), Op<T>::u(D) ): T( Op<T>::l(D), Op<T>::u(E) ) );
      }
    }
  }
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
    throw typename ISModel<T>::Exceptions( ISModel<T>::Exceptions::INV );
  
  // Asymetric inclusion
  if( options.ASYREM_USE ){
    auto const& f = [=]( const T& x ){ return Op<T>::inv( x ); };
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
    auto const& f = [=]( const T& x ){ return Op<T>::sqr( x ); };
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
( std::vector<std::vector<T>>& mat, unsigned const& ndep )
const
{
  assert( !mat.empty() );

  // Bounds
  T bnd = _B( mat, 1 );
  double L( Op<T>::l(bnd) );
  if ( L < 0. )
    throw typename ISModel<T>::Exceptions( ISModel<T>::Exceptions::SQRT );

  // Asymetric inclusion
  if( options.ASYREM_USE ){
    auto const& f = [=]( const T& x ){ return Op<T>::sqrt( x ); };
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
( std::vector<std::vector<T>>& mat, unsigned const& ndep )
const
{
  assert( !mat.empty() );

  // Bounds
  T bnd = _B( mat, 1 );

  // Asymetric inclusion
  if( options.ASYREM_USE ){
    auto const& f = [=]( const T& x ){ return Op<T>::exp( x ); };
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
    throw typename ISModel<T>::Exceptions( ISModel<T>::Exceptions::LOG );

  // Asymetric inclusion
  if( options.ASYREM_USE ){
    auto const& f = [=]( const T& x ){ return Op<T>::log( x ); };
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
( std::vector<std::vector<T>>& mat, unsigned const& ndep )
const
{
  assert( !mat.empty() );

  // Bounds
  T bnd = _B( mat, 1 );
  double L( Op<T>::l(bnd) );
  if ( L < 0. )
    throw typename ISModel<T>::Exceptions( ISModel<T>::Exceptions::LOG );

  // Asymetric inclusion
  auto const& f = [=]( const T& x ){ return Op<T>::xlog( x ); };
  return _asym( mat, ndep, f, std::exp(-1.), true, bnd ); // concave term, min 1/e
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
  T bnd = T(-1.,1.)*(rem/double(ndep3)) - w/double(ndep3);
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
      for( unsigned int j=0; j<_ndiv; j++ )
        mat3[i][j] = (mat1[i][j]+(C1-_c1[i]))*(mat2[i][j]+(C2-_c2[i])) - (C1-_c1[i])*(C2-_c2[i]) + bnd;
    }
    else if( !mat1[i].empty() ){
      if( mat3[i].empty() ) mat3[i].resize(_ndiv);
#ifdef MC__ISMODEL_DEBUG_PROD
      std::cerr << "1 ONLY" << std::endl;
#endif
      for( unsigned int j=0; j<_ndiv; j++ )
        mat3[i][j] = mat1[i][j]*C2 + bnd;
    }
    else if( !mat2[i].empty() ){
      if( mat3[i].empty() ) mat3[i].resize(_ndiv);
      for( unsigned int j=0; j<_ndiv; j++ )
        mat3[i][j] = mat2[i][j]*C1 + bnd;
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
    ( ISVar<U> const&, ISVar<U> const& );
  template <typename U> friend ISVar<U> operator+
    ( double const&, ISVar<U> const& );
  template <typename U> friend ISVar<U> operator+
    ( ISVar<U> const&, double const& );

  template <typename U> friend ISVar<U> operator-
    ( ISVar<U> const& );
  template <typename U> friend ISVar<U> operator-
    ( ISVar<U> const&, ISVar<U> const& );
  template <typename U> friend ISVar<U> operator-
    ( double const&, ISVar<U> const& );
  template <typename U> friend ISVar<U> operator-
    ( ISVar<U> const&, double const& );

  template <typename U> friend ISVar<U> operator*
    ( ISVar<U> const&, ISVar<U> const& );
  template <typename U> friend ISVar<U> operator*
    ( double const&, ISVar<U> const& );
  template <typename U> friend ISVar<U> operator*
    ( ISVar<U> const&, double const& );

  template <typename U> friend ISVar<U> operator/
    ( ISVar<U> const&, ISVar<U> const& );
  template <typename U> friend ISVar<U> operator/
    ( double const&, ISVar<U> const& );
  template <typename U> friend ISVar<U> operator/
    ( ISVar<U> const&, double const& );

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
  template <typename U> friend ISVar<U> pow
    ( ISVar<U> const& , int const& n );
  template <typename U> friend ISVar<U> pow
    ( ISVar<U> && , int const& n );
  template <typename U> friend ISVar<U> cheb
    ( ISVar<U> const& , unsigned int const& n );
  template <typename U> friend ISVar<U> cheb
    ( ISVar<U> && , unsigned int const& n );

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
  double _cst;
  //! @brief Variable bound
  mutable std::pair<T,bool> _bnd;

 public:
 
  ISVar<T>& operator+=
    ( ISVar<T> const& );
  ISVar<T>& operator+=
    ( double const& );
 
  ISVar<T>& operator-=
    ( ISVar<T> const& );
  ISVar<T>& operator-=
    ( double const& );
 
  ISVar<T>& operator*=
    ( ISVar<T> const& );
  ISVar<T>& operator*=
    ( double const& );
 
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
    _cst = cst;
    _bnd = std::make_pair( 0., false );
    return *this;
  }

  ISVar<T>& operator=
  ( ISVar<T> const& var )
  {
#ifdef MC__ISMODEL_TRACE
    std::cerr << "-- ISVar<T>& operator= ( ISVar<T> const& )\n";
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
    return *this;
  }

  ISVar<T>& operator=
  ( ISVar<T> && var )
  {
#ifdef MC__ISMODEL_TRACE
    std::cerr << "-- ISVar<T>& operator= ( ISVar<T> && )\n";
#endif
    if( this == &var )
      return *this;
    _mod = var._mod;
    if( !_mod ) _cst = var._cst;
    _nvar = var._nvar;
    _ndiv = var._ndiv;
    _ndep = var._ndep;
    //_mat = var._mat;
    std::swap( _mat, var._mat );
    _bnd = var._bnd;
    return *this;
  }

  ISVar
  ( ISModel<T>* const mod, unsigned int ndx, T const& bnd )
  : _mod(mod), _ndiv(mod->_ndiv), _nvar(mod->_nvar), _ndep(1), _mat(_nvar), _bnd(bnd,true)
  {
    if( ndx >= _nvar )
      throw typename ISModel<T>::Exceptions( ISModel<T>::Exceptions::INDEX );

    double l = Op<T>::l(bnd);
    double h = Op<T>::diam(bnd) / (double)_ndiv;
    _mat[ndx].resize( _ndiv );
    for( unsigned int j=0; j<_ndiv; j++, l+=h )
      _mat[ndx][j] = T( l, l+h );
    _mod->_defvar[ndx] = true;
    _mod->_bndvar[ndx] = bnd;
  }

  ISVar
  ( double const& cst=0. )
  : _mod(nullptr), _ndiv(0), _nvar(0), _ndep(0), _cst(cst), _bnd(0.,false)
  {}

  ISVar
  ( ISVar<T> const& var )
  : _mod(var._mod), _ndiv(var._ndiv), _nvar(var._nvar), _ndep(var._ndep)
  {
#ifdef MC__ISMODEL_TRACE
    std::cerr << "-- ISVar( ISVar<T> const& )\n";
#endif
    if( this == &var ) return;
    if( !_mod ) _cst = var._cst;
    _mat = var._mat;
    _bnd = var._bnd;
  }

  ISVar
  ( ISVar<T> && var )
  : _mod(var._mod), _ndiv(var._ndiv), _nvar(var._nvar), _ndep(var._ndep)
  {
#ifdef MC__ISMODEL_TRACE
    std::cerr << "-- ISVar( ISVar<T> && var )\n";
#endif
    if( this == &var ) return;
    if( !_mod ) _cst = var._cst;
    //_mat = var._mat;
    std::swap( _mat, var._mat );
    _bnd = var._bnd;
  }

  ~ISVar
  () 
  {}

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
    double l = Op<T>::l(bnd);
    double h = Op<T>::diam(bnd) / (double)_ndiv;
    _mat[ndx].resize( _ndiv );
    for( unsigned int j=0; j<_ndiv; j++, l+=h )
      _mat[ndx][j] = T( l, l+h );

    _mod->_defvar[ndx] = true;
    _mod->_bndvar[ndx] = bnd;
    _bnd = std::make_pair( bnd, true );

    return *this;
  }

  std::vector<std::vector<T>> const& C
  ()
  const
  { return _mat; }

  double const& cst
  ()
  const
  { return _cst; }

  unsigned int const& ndep
  ()
  const
  { return _ndep; }

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
    for( unsigned int i=0; i<_nvar; i++ ){
      if( _mat[i].empty() ) continue;
      double l = Op<T>::l( _mod->_bndvar[i] );
      double u = Op<T>::u( _mod->_bndvar[i] );
      assert( point[i] >= l && point[i] <= u );
      double h = Op<T>::diam( _mod->_bndvar[i] ) / (double)_ndiv;
      int ndx = std::floor( ( point[i] - l ) /  h );
      if( ndx < 0 ) ndx = 0;
      if( ndx >= (int)_ndiv ) ndx = _ndiv-1;      
      //std::cout << "ndx = " << ndx << std::endl;
      //{ int dum; std::cin >> dum; }
      val += _mat[i][ndx];
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
  : _mod(mod), _ndiv(mod->_ndiv), _nvar(mod->_nvar ), _mat(_nvar), _bnd(0.,false)
  {}

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

  for( unsigned int i=0; i<_nvar; i++ ){
    if( _mat[i].empty() ) continue;
    for( unsigned int j=0; j<_ndiv; j++ )
      _mat[i][j] += cst / (double)_ndep;
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
    const double copy_cst = _cst;
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

  for( unsigned int i=0; i<_nvar; i++ ){
    if( !_mat[i].empty() && !var._mat[i].empty() )
      for( unsigned int j=0; j<_ndiv; j++ )
        _mat[i][j] += var._mat[i][j];
    else if( !var._mat[i].empty() ){
      _ndep++;                // add dependency
      _mat[i] = var._mat[i];  // copy entire row
    }
  }
  if( _bnd.second ) _bnd.first += var.B();

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
( ISVar<T> const& var )
{
  if( !_mod && !var._mod ){
    _cst -= var._cst;
    return *this;
  }
  if( !_mod ){
    const double copy_cst = _cst;
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
  if( _bnd.second ) _bnd.first -= var.B();

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

  return *this;
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
    const double copy_cst = _cst;
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

  //if( var1._mod->options.ASYREM_USE )
  //  return 0.25 * ( sqr( var1 + var2 ) - sqr( var1 - var2 ) ); 
  
  ISVar<T> var3( var2 );
  var3 *= var1;
  return var3;
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
( double const& cst1, ISVar<T> const& var2 )
{
  ISVar<T> var3( inv( var2 ) );
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
  var2._mod->_inv( var2._mat, var2._ndep );
  var2._bnd.second = false;
  return var2;
}

template <typename T>
inline
ISVar<T> inv
( ISVar<T> && var )
{
  if( !var._mod )
    return 1./var._cst;

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

//  ISVar<T> var2( var );
//  auto const& finc = [=]( const T& x ){ return Op<T>::sqr( Op<T>::max( x, 0. ) ); };
//  var2._mod->_asym( var2._mat, var2._ndep, finc, -DBL_MAX, true ); // convex term, min 0
//  var2._bnd.second = false;

//  ISVar<T> var3( var );
//  auto const& fdec = [=]( const T& x ){ return Op<T>::sqr( Op<T>::min( x, 0. ) ); };
//  var3._mod->_asym( var3._mat, var3._ndep, fdec, DBL_MAX, true ); // convex term, min 0
//  var3._bnd.second = false;
//  return var2 + var3;

  ISVar<T> var2( var );
  var2._mod->_sqr( var2._mat, var2._ndep );
  var2._bnd.second = false;
  return var2;
}

template <typename T>
inline
ISVar<T> sqr
( ISVar<T> && var )
{
  if( !var._mod )
    return sqr(var._cst);

  var._mod->_sqr( var._mat, var._ndep );
  var._bnd.second = false;
  return var;
}

template <typename T>
inline
ISVar<T> fabs
( ISVar<T> const& var )
{
  if( !var._mod )
    return std::fabs(var._cst);

  ISVar<T> var2( var );
  auto const& f = [=]( const T& x ){ return Op<T>::fabs( x ); };
  var2._mod->_asym( var2._mat, var2._ndep, f, 0., true ); // convex term, min 0
  var2._bnd.second = false;
  return var2;
}

template <typename T>
inline
ISVar<T> fabs
( ISVar<T> && var )
{
  if( !var._mod )
    return std::fabs(var._cst);

  auto const& f = [=]( const T& x ){ return Op<T>::fabs( x ); };
  var._mod->_asym( var._mat, var._ndep, f, 0., true ); // convex term, min 0
  var._bnd.second = false;
  return var;
}

template <typename T>
inline
ISVar<T> relu
( ISVar<T> const& var )
{
  if( !var._mod )
    return relu(var._cst);

  ISVar<T> var2( var );
  auto const& f = [=]( const T& x ){ return Op<T>::max( x, 0. ); };
  var2._mod->_asym( var2._mat, var2._ndep, f, 0., true ); // convex term, min 0
  var2._bnd.second = false;
  return var2;
}

template <typename T>
inline
ISVar<T> relu
( ISVar<T> && var )
{
  if( !var._mod )
    return relu(var._cst);

  auto const& f = [=]( const T& x ){ return Op<T>::max( x, 0. ); };
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
    return std::sqrt(var._cst);

  ISVar<T> var2( var );
  var2._mod->_sqrt( var2._mat, var2._ndep );
  var2._bnd.second = false;
  return var2;
}

template <typename T>
inline
ISVar<T> sqrt
( ISVar<T> && var )
{
  if( !var._mod )
    return std::sqrt(var._cst);

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
    return std::exp(var._cst);

  ISVar<T> var2( var );
  var2._mod->_exp( var2._mat, var2._ndep );
  var2._bnd.second = false;
  return var2;
}

template <typename T>
inline
ISVar<T> exp
( ISVar<T> && var )
{
  if( !var._mod )
    return std::exp(var._cst);

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
    return std::log(var._cst);

  ISVar<T> var2( var );
  var2._mod->_log( var2._mat, var2._ndep );
  var2._bnd.second = false;
  return var2;
}

template <typename T>
inline
ISVar<T> log
( ISVar<T> && var )
{
  if( !var._mod )
    return std::log(var._cst);

  var._mod->_log( var._mat, var._ndep );
  var._bnd.second = false;
  return var;
}

template <typename T>
inline
ISVar<T> xlog
( ISVar<T> const& var )
{
  if( !var._mod )
    return xlog(var._cst);

  ISVar<T> var2( var );
  var2._mod->_xlog( var2._mat, var2._ndep );
  var2._bnd.second = false;
  return var2;
}

template <typename T>
inline
ISVar<T> xlog
( ISVar<T> && var )
{
  if( !var._mod )
    return xlog(var._cst);

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
    return std::sin(var._cst);

  ISVar<T> var2( var );
  var2._mod->_sin( var2._mat, var2._ndep );
  var2._bnd.second = false;
  return var2;
}

template <typename T>
inline
ISVar<T> sin
( ISVar<T> && var )
{
  if( !var._mod )
    return std::sin(var._cst);

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
    return std::cos(var._cst);

  ISVar<T> var2( var );
  var2._mod->_cos( var2._mat, var2._ndep );
  var2._bnd.second = false;
  return var2;
}

template <typename T>
inline
ISVar<T> cos
( ISVar<T> && var )
{
  if( !var._mod )
    return std::cos(var._cst);

  var._mod->_cos( var._mat, var._ndep );
  var._bnd.second = false;
  return var;
}

template <typename T>
inline
ISVar<T> tanh
( ISVar<T> const& var )
{
  return 1-2/(exp(2*var)+1);
}


// @TODO: since power functions are convex/concave on certain domain, the Alg 2 is applicable here as well. 
template <typename T>
inline
ISVar<T>
pow
( ISVar<T> const& var, int const& n )
{
  if( !var._mod )
    return std::pow( var._cst, n );

  if( n < 0 )
    return pow( inv( var ), -n );
  if( n == 0 )
    return 1.;
  if( n == 1 )
    return var;

  return // recursive call
    !(n%2) ?
    sqr( pow( var, n/2 ) ) :
    sqr( pow( var, n/2 ) ) * var;
}


template <typename T>
inline
ISVar<T>
pow
( ISVar<T> && var, int const& n )
{
  if( !var._mod )
    return std::pow( var._cst, n );

  if( n < 0 )
    return pow( inv( var ), -n );
  if( n == 0 )
    return 1.;
  if( n == 1 )
    return var;

  return // recursive call
    !(n%2) ?
    sqr( pow( var, n/2 ) ) :
    sqr( pow( var, n/2 ) ) * var;
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
    return cheb( var._cst, n );

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
( ISVar<T> && var, unsigned int const& n )
{
  if( !var._mod )
    return cheb( var._cst, n );

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
  if( !var1._mod && !var2._mod ) 
    return std::max( var1._cst, var2._cst );
  if( !var2._mod ) 
    return max( var1, var2._cst );
  if( !var1._mod ) 
    return max( var2, var1._cst );
  return max( var1 - var2, 0. ) + var2;
}

template <typename T>
inline
ISVar<T>
max
( ISVar<T> const& var1, double const& cst2 )
{
  if( !var1._mod ) 
    return std::max( var1._cst, cst2 );

  ISVar<T> var2( var1 );
  auto const& f = [=]( const T& x ){ return Op<T>::max( x, cst2 ); };
  var2._mod->_asym( var2._mat, var2._ndep, f, cst2, true ); // convex term, min cst2
  var2._bnd.second = false;
  return var2;
}

template <typename T>
inline
ISVar<T>
max
( ISVar<T> && var1, double const& cst2 )
{
  if( !var1._mod ) 
    return std::max( var1._cst, cst2 );

  auto const& f = [=]( const T& x ){ return Op<T>::max( x, cst2 ); };
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
  if( !var1._mod && !var2._mod ) 
    return std::min( var1._cst, var2._cst );
  if( !var2._mod ) 
    return min( var1, var2._cst );
  if( !var1._mod ) 
    return min( var2, var1._cst );
  return min( var1 - var2, 0. ) + var2;
}

template <typename T>
inline
ISVar<T>
min
( ISVar<T> const& var1, double const& cst2 )
{
  if( !var1._mod ) 
    return std::min( var1._cst, cst2 );

  ISVar<T> var2( var1 );
  auto const& f = [=]( const T& x ){ return Op<T>::min( x, cst2 ); };
  var2._mod->_asym( var2._mat, var2._ndep, f, cst2, false ); // convex term, max cst2
  var2._bnd.second = false;
  return var2;
}

template <typename T>
inline
ISVar<T>
min
( ISVar<T> && var1, double const& cst2 )
{
  if( !var1._mod ) 
    return std::min( var1._cst, cst2 );

  auto const& f = [=]( const T& x ){ return Op<T>::min( x, cst2 ); };
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
  static ISV tanh(const ISV& x) { return mc::tanh(x); }
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
