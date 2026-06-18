// Copyright (C) Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the Eclipse Public License.

/*!
\page page_ELLIPSOID Ellipsoidal Calculus and Ellipsoidal Arithmetic for
Factorable Functions \author Mario E. Villanueva, Jai Rajyaguru, Boris Houska,
Beno&icirc;t Chachuat

An ellipsoid with center \f$c \in \mathbb R^n\f$ and shape matrix \f$Q \in
\mathbb S_{+}^n\f$ is defined as \f{align*} \mathcal E(c,Q) := & \left\{ \left.
c + Q^\frac{1}{2} v \ \right| \ \exists v \in \mathbb R^{n}: \ v^T v \leq 1
\right\} \subseteq \mathbb R^{n} . \f}

As well as constructors and data access/manipulations functions, the class
mc::Ellipsoid provides a set of functions for ellipsoidal calculus. Data
manipulation functions include:
- checking positive semi-definiteness of the shape matrix (mc::Ellipsoid::psdQ)
- determining the rank of or regularizing the shape matrix
(mc::Ellipsoid::rankQ, mc::Ellipsoid::regQ)
- computing the eigenvalues or singular values of the shape matrix
(mc::Ellipsoid::eigQ, mc::Ellipsoid::svdQ)
- taking the trace, square root or inverse of the shape matrix
(mc::Ellipsoid::trQ, mc::Ellipsoid::sqrtQ, mc::Ellipsoid::invQ)
.
Ellipsoidal calculus includes:
- applying a linear transformation to an ellipsoid (mc::mtimes)
- taking the (exact) intersection of an ellipsoid with a hyperplane
(mc::hpintersection)
- computing a (minimum volume) external ellipsoidal approximation of the
intersection between an ellipsoid and a halfspace, or between two ellipsoids
(mc::intersection_ea)
- computing an external ellipsoidal approximation of the geometric (Minkowski)
sum of several ellipsoids along a given direction, or a (minimum trace) external
ellipsoidal approximation of the geometric (Minkowski) sum of an ellipsoid with
an interval box (mc::minksum_ea)
.

Besides ellipsoidal calculus, the classes mc::EllImg and mc::EllVar provide an
implementation of ellipsoidal arithmetic in order to enclose the image
\f$\mathcal E(c_f,Q_f)\f$ of an \f$n_x\f$-dimensional ellispoid \f$\mathcal
E(c_x,Q_x)\f$ under a vector-valued function \f$
f:\mathbb{R}^{n_x}\to\mathbb{R}^{n_f} \f$: \f{align*} \mathcal E(c_f,Q_f)
\supseteq & \left\{ f(x) \,\mid\, x\in \mathcal{E}(c_{x},Q_x) \right\}. \f}
Notice that the exact image \f$\left\{ f(x) \,\mid\, x\in \mathcal{E}(c_{x},Q_x)
\right\}\f$ is not an ellipsoid in general.

The class mc::EllImg is derived from mc::Ellipsoid. The implementation of
mc::EllImg and mc::EllVar relies on the operator/function overloading mechanism
of C++. This makes the computation of the ellipsoidal enclosure for the image of
an ellipsoid under a factorable function both simple and intuitive, similar to
computing function values in real arithmetic or bounds based on interval, Taylor
or Chebyshev model arithmetics (see \ref page_INTERVAL, \ref page_TAYLOR, \ref
page_CHEBYSHEV). mc::EllImg stores the center of a <i>lifted</i> ellipsoid as a
column vector (arma::vec) and its shape as a symmetric matrix kept in a hashed,
lower-triangular dictionary-of-keys (a <tt>std::unordered_map</tt> from packed
(row,col) index to value), which supports O(1)-average incremental get/set as
the lifting proceeds; an Armadillo sparse matrix (arma::sp_mat) is materialised
on demand, e.g. via mc::EllImg::Q_lift or when the image is streamed to an
output. The lifted ellipsoid is the result of adding an extra dimension for each
elementary operation participating in the factorable function. The nonlinearity
of each univariate elementary operation is enclosed by a linear term plus an
interval remainder: by default a degree-1 Remez minimax line is used (option
<tt>REMEZ_USE</tt>), with a secant-line fallback, and the remainder radius is
obtained from a sampled (non-verified) residual bound. Note that the
implementation in mc::EllImg and mc::EllVar is <a>not verified</a> in the sense
that rounding errors are not accounted for during the propagation.

The classes mc::EllImg and mc::EllVar are templated in the interval type used to
bound the nonlinearity of the function. By default, mc::EllImg and mc::EllVar
can be used with the non-verified interval type mc::Interval of MC++. For
reliability, however, it is recommended to use verified interval arithmetic such
as <A
href="http://www.ti3.tu-harburg.de/Software/PROFILEnglisch.html">PROFIL</A>
(header file <tt>mcprofil.hpp</tt>), <A
href="http://www.math.uni-wuppertal.de/~xsc/software/filib.html">FILIB++</A>
(header file <tt>mcfilib.hpp</tt>), or the Boost interval library (header file
<tt>mcboost.hpp</tt>).

\section sec_ELLCALC How do I define an ellipsoid and apply ellipsoidal
calculus?

In order to define the ellipsoid \f$\mathcal E(c_x,Q_x)\f$ with
\f{align*}
  c_x = & \left(\begin{array}{c} 3\\4\end{array}\right),\ \text{and} &
  Q_x = & \left(\begin{array}{cc} 5 & 4 \\ 4 & 5 \end{array}\right).
\f}

we proceed as follows:

\code
    const unsigned int n = 2;
    arma::vec cx(n);
    arma::mat Qx(n,n,arma::fill::zeros);
    cx(0) = 3.;  Qx(0,0) = 5.;
    cx(1) = 4.;  Qx(1,0) = 4.;  Qx(0,1) = 4.;  Qx(1,1) = 5.;
    mc::Ellipsoid Ex( Qx, cx );
\endcode

This ellipsoid is simply displayed as:

\code
    std::cout << Ex << std::endl;
\endcode

In the present case, the following information is displayed:

\verbatim
center:
   3.0000
   4.0000
shape:
   5.0000   4.0000
   4.0000   5.0000
\endverbatim

In order to illustrate ellipsoidal calculus, suppose that we want to add the
interval box \f$[0,0.1]^2\f$ to the foregoing ellispoid and determine an
external ellispoidal approximation of this geometric sum, using mc::minksum_ea.

\section sec_ELLIMG How do I compute an ellipsoidal enclosure for the image set
of an ellipsoid under a factorable function?

Suppose that we want to compute an ellipsoidal enclosure for the image set of
function: \f[ f(x) = \left(\begin{array}{c}
\log(x_{1})+x^{2}_{2}        \\
\sin(x_{1})-\cos(x_{2})       \\
\end{array} \right) \qquad \text{with} \qquad
x \in \left\{
\left(\begin{array}{c}
3 \\
4 \\
\end{array} \right)+
\left(\begin{array}{cc}
5 & 4 \\
4 & 5 \\
\end{array} \right)^{1/2}v
:
v^{\rm T}v\leq 1 \right\}.
\f]

For simplicity, the underlying interval bounds are propagated using the default
interval type mc::Interval. The required header files are:

\code
#include "ellimage.hpp"
#include "interval.hpp"

typedef mc::Interval I;
typedef mc::EllImg<I> EI;
typedef mc::EllVar<I> EV;
typedef arma::vec dcv;
typedef arma::mat dsm;
\endcode

First, the host ellipsoidal set for the independent variables \f$x_{1}\f$ and
\f$x_{2}\f$ is specified as: \code dcv cx(2); dsm Qx(2,2,arma::fill::zeros);
cx(0) = 3.;  Qx(0,0) = 5.;
cx(1) = 4.;  Qx(1,0) = 4.;  Qx(0,1) = 4.;  Qx(1,1) = 5.;
EI Ex( Qx, cx );
\endcode

Then, the independent variables themselves are specified as:
\code
EV X1( Ex, 0 );
EV X2( Ex, 1 );
\endcode
If independent interval bounds are known for a variable they can be passed as an
optional third argument to the mc::EllVar constructor.

The dependent variables \f$f_{1}(x)\f$ and \f$f_{2}(x)\f$ are propagated as:
\code
EV F[2] = { log( X1 ) + mc::sqr( X2 ),
            sin( X1 ) - cos( X2 ) };
\endcode

and the lifted ellipsoid can be displayed as:

\code
std::cout << "lifted ellipsoidal image of f =" << Ex << std::endl;
\endcode

The display reports the lifted center followed by the (sparse) lifted shape
matrix.  Its center and shape can also be retrieved directly with
mc::EllImg::c_lift and mc::EllImg::Q_lift.  For the function above the lifting
introduces 6 auxiliary dimensions (one per elementary operation), giving an
8-dimensional lifted ellipsoid with center:

\verbatim
   3.00000e+00
   4.00000e+00
   1.85000e+01
   9.13697e-01
   1.94137e+01
  -1.93826e-01
   6.44906e-02
   2.58154e-01
\endverbatim

and (dense) lifted shape matrix Ex.Q_lift():

\verbatim
   5.9494e+00   4.7595e+00   3.8076e+01   2.5607e+00   4.0637e+01   1.2945e+00
-2.7081e+00  -4.0026e+00
   4.7595e+00   5.9494e+00   4.7595e+01   2.0485e+00   4.9644e+01   1.6181e+00
-2.1664e+00  -3.7846e+00
   3.8076e+01   4.7595e+01   4.3522e+02   1.6388e+01   4.5161e+02   1.2945e+01
-1.7332e+01  -3.0277e+01
   2.5607e+00   2.0485e+00   1.6388e+01   5.9579e+00   2.2346e+01   5.5716e-01
-1.1656e+00  -1.7227e+00
   4.0637e+01   4.9644e+01   4.5161e+02   2.2346e+01   4.7395e+02   1.3502e+01
-1.8497e+01  -3.1999e+01
   1.2945e+00   1.6181e+00   1.2945e+01   5.5716e-01   1.3502e+01   1.9497e+01
-5.8924e-01  -2.0087e+01 -2.7081e+00  -2.1664e+00  -1.7332e+01  -1.1656e+00
-1.8497e+01  -5.8924e-01   1.3491e+01   1.4081e+01 -4.0026e+00  -3.7846e+00
-3.0277e+01  -1.7227e+00  -3.1999e+01  -2.0087e+01   1.4081e+01   3.4167e+01
\endverbatim

Finally, the lifted ellipsoid is projected onto the space of the dependent
variables to obtain the desired image enclosure (note that mc::EllImg::get
returns a plain mc::Ellipsoid): \code mc::Ellipsoid Ef = Ex.get( 2, F );
std::cout << "Ellipsoidal enclosure Ef = " << Ef << std::endl;
\endcode

\verbatim
Ellipsoidal enclosure Ef =
center:
   1.94137e+01
   2.58154e-01
shape:
   4.7395e+02  -3.1999e+01
  -3.1999e+01   3.4167e+01
\endverbatim

After this, the ellipsoid can be manipulated according to the rules of
ellipsoidal calculus (see \ref page_ELLIPSOID). A comparison between the exact
image and the ellipsoidal enclosure is presented in the following figure, in
orange the exact image and in blue the boundary of the ellipsoidal enclosure.
<CENTER><TABLE BORDER=0>
<TR>
<TD>\image html ELL-2D.png</TD>
</TR>
</TABLE></CENTER>


\section sec_ELL_opt How are the options set for ellipsoidal calculus and
ellipsoidal arithmetic?

The classes mc::EllImg and mc::Ellipsoid have public members called
mc::EllImg::options and mc::Ellipsoid::options (both static) that can be used to
set/modify a number of options. Note that mc::EllImg::options and
mc::Ellipsoid::options are distinct option sets (mc::EllImg is derived from
mc::Ellipsoid but defines its own Options structure). For instance, options can
be set as follows: \code mc::EllImg<I>::options.PREALLOC  = 5;
        mc::EllImg<I>::options.REMEZ_USE = true;
        mc::Ellipsoid::options.PSDCHK    = false;
\endcode

The full set of available options is reported in the following tables.

<TABLE border="1">
<CAPTION><EM>Options in mc::EllImg::Options: name, type and
description</EM></CAPTION> <TR><TH><b>Name</b> <TD><b>Type</b><TD><b>Default</b>
         <TD><b>Description</b>
     <TR><TH><tt>PREALLOC</tt> <TD><tt>long</tt> <TD>0
         <TD> Number of rows to preallocate in the center vector when lifting;
<tt>0</tt> uses geometric (doubling) growth <TR><TH><tt>MINK_TOL</tt>
<TD><tt>double</tt> <TD>1e-10 <TD> Tolerance used in the Minkowski-sum step that
adds the linearisation remainder <TR><TH><tt>REMEZ_USE</tt> <TD><tt>bool</tt>
<TD>true <TD> Whether to use a degree-1 Remez minimax line for univariate terms
(otherwise a secant-line approximation is used) <TR><TH><tt>REMEZ_MAXIT</tt>
<TD><tt>unsigned</tt> <TD>5 <TD> Maximum number of iterations in the Remez
exchange algorithm <TR><TH><tt>REMEZ_TOL</tt> <TD><tt>double</tt> <TD>1e-5 <TD>
Stopping tolerance in the Remez exchange algorithm <TR><TH><tt>REMEZ_MIG</tt>
<TD><tt>double</tt> <TD>1e-10 <TD> Minimum domain diameter for invoking Remez;
smaller domains fall back to the secant-line approximation
     <TR><TH><tt>DCPROD_USE</tt> <TD><tt>bool</tt> <TD>false
         <TD> Whether to lift a bilinear product via a difference-of-convex
decomposition rather than the midpoint linearisation
</TABLE>

<TABLE border="1">
<CAPTION><EM>Options in mc::Ellipsoid::Options: name, type and
description</EM></CAPTION> <TR><TH><b>Name</b> <TD><b>Type</b><TD><b>Default</b>
         <TD><b>Description</b>
     <TR><TH><tt>PSDCHK</tt> <TD><tt>bool</tt> <TD>false
         <TD>Whether or not to check positive semi-definiteness of shape
matrices <TR><TH><tt>PSDTOL</tt> <TD><tt>double</tt> <TD>1e2*MACHPREC
         <TD>Absolute tolerance for positive semi-definiteness check of shape
matrices <TR><TH><tt>RKTOLA</tt> <TD><tt>double</tt> <TD>MACHPREC <TD>Absolute
tolerance for rank and regularization of shape matrices <TR><TH><tt>RKTOLR</tt>
<TD><tt>double</tt> <TD>MACHPREC <TD>Relative tolerance for rank and
regularization of shape matrices <TR><TH><tt>ROOTTOL</tt> <TD><tt>double</tt>
<TD>1e-10 <TD>Absolute stopping tolerance for root-finding method (objective
function value less than ROOTTOL) <TR><TH><tt>ROOTSECANT</tt> <TD><tt>bool</tt>
<TD>false <TD>Whether to use the secant method for root finding
     <TR><TH><tt>ROOTMAXIT</tt> <TD><tt>unsigned</tt> <TD>0
         <TD>Maximum number of iterations for root-finding method (no maximum
when ROOTMAXIT=0)
</TABLE>


\section sec_ELL_err What are the errors encountered during ellipsoidal calculus
and ellipsoidal arithmetic?

Errors are managed based on the exception handling mechanism of the C++
language. Each time an error is encountered, a class object of type
mc::EllImg::Exceptions or mc::Ellipsoid::Exceptions is thrown, which contains
the type of error. It is the user's responsibility to test whether an exception
was thrown during the computation of the lifted ellipsoid, and then make the
appropriate changes. Should an exception be thrown and not caught by the calling
program, the execution will abort.

Possible errors encountered during application of ellipsoidal calculus and
ellispoidal arithmetic are reported in the following tables.

<TABLE border="1">
<CAPTION><EM>Errors during the Computation of an Ellipsoidal Image
(mc::EllImg::Exceptions)</EM></CAPTION> <TR><TH><b>Number</b>
<TD><b>Description</b> <TR><TH> <tt> 1 </tt>  <TD> Division by zero scalar
     <TR><TH> <tt> 2 </tt>  <TD> Inverse operation with zero in domain
     <TR><TH> <tt> 3 </tt>  <TD> Log operation with non-positive numbers in
domain <TR><TH> <tt> 4 </tt>  <TD> Square-root operation with negative numbers
in domain <TR><TH> <tt> 5 </tt>  <TD> Tangent operation with zero in domain of
cosine, tan(x) = sin(x)/cos(x) <TR><TH> <tt> 6 </tt>  <TD> Inverse cosine (acos)
operation with domain outside [-1,1] <TR><TH> <tt> 7 </tt>  <TD> Inverse sine
(asin) operation with domain outside [-1,1] <TR><TH> <tt>-1 </tt>  <TD> Failed
to construct ellipsoidal variable <TR><TH> <tt>-3 </tt>  <TD> Operation between
variables mc::EllVar linked to different images mc::EllImg <TR><TH> <tt>-33</tt>
<TD> Feature not yet implemented in mc::EllImg
</TABLE>

<TABLE border="1">
<CAPTION><EM>Errors with mc::Ellipsoid
(mc::Ellipsoid::Exceptions)</EM></CAPTION> <TR><TH><b>Number</b>
<TD><b>Description</b> <TR><TH> <tt> 1 </tt> <TD> Non positive-semi definite
shape matrix <TR><TH> <tt> 2 </tt> <TD> Failure in a LAPACK linear algebra
routine <TR><TH> <tt> 3 </tt> <TD> Failure in a root-finding routine
</TABLE>

Moreover, exceptions may be thrown by the template parameter class itself.

\section sec_ELL_refs References

- Kurzhanskiy, A., and P. Varaiya, <A
href="http://code.google.com/p/ellipsoids">"Ellipsoidal Toolbox"</A>, Technical
Report UCB/EECS-2006-46, EECS Department, University of California, Berkeley,
May 2006.
- Villanueva, M.E., J. Rajyagurua, B. Houskab, B. Chachuat, <A
href="https://doi.org/10.1016/B978-0-444-63578-5.50123-7">"Ellipsoidal
Arithmetic for Multivariate Systems"</A>, <i>Computer Aided Chemical
Engineering</i>, <b>37</b>, 767-772, 2015.
*/

#ifndef MC__ELLIMAGE_H
#define MC__ELLIMAGE_H

#include <assert.h>

#include <algorithm>
#include <cmath>
#include <exception>
#include <fstream>
#include <iomanip>
#include <unordered_map>
#include <vector>

#include "ellipsoid.hpp"
#include "mcfunc.hpp"
#include "mcop.hpp"

#undef MC__ELLIMAGE_DEBUG

namespace mc
{

template <class T>
class EllVar;

//! @brief C++ class for ellipsoidal arithmetic - Ellipsoidal image environment
////////////////////////////////////////////////////////////////////////
//! mc::EllImg is a C++ class for definition of ellipsoidal image
//! environment, derived from mc::Ellipsoid. Computation of ellipsoidal
//! image for factorable functions is via the C++ class mc::EllVar. The
//! template parameter corresponds to the type used to propagate
//! variable range. Round-off errors are not accounted for in the
//! computations (non-verified implementation).
////////////////////////////////////////////////////////////////////////
template <class T>
class EllImg : public Ellipsoid
////////////////////////////////////////////////////////////////////////
{
  friend class EllVar<T>;

  template <class U>
  friend std::ostream& operator<<(std::ostream&, const EllImg<U>&);
  template <class U>
  friend EllVar<U> inv(const EllVar<U>&);
  template <class U>
  friend EllVar<U> exp(const EllVar<U>&);
  template <class U>
  friend EllVar<U> log(const EllVar<U>&);
  template <class U>
  friend EllVar<U> xlog(const EllVar<U>&);
  template <class U>
  friend EllVar<U> sqrt(const EllVar<U>&);
  template <class U>
  friend EllVar<U> sqr(const EllVar<U>&);
  template <class U>
  friend EllVar<U> pow(const EllVar<U>&, const int);
  template <class U>
  friend EllVar<U> cheb(const EllVar<U>&, const unsigned);
  template <class U>
  friend EllVar<U> cos(const EllVar<U>&);
  template <class U>
  friend EllVar<U> sin(const EllVar<U>&);
  template <class U>
  friend EllVar<U> tan(const EllVar<U>&);
  template <class U>
  friend EllVar<U> acos(const EllVar<U>&);
  template <class U>
  friend EllVar<U> asin(const EllVar<U>&);
  template <class U>
  friend EllVar<U> atan(const EllVar<U>&);
  template <class U>
  friend EllVar<U> cosh(const EllVar<U>&);
  template <class U>
  friend EllVar<U> sinh(const EllVar<U>&);
  template <class U>
  friend EllVar<U> tanh(const EllVar<U>&);
  template <class U>
  friend EllVar<U> erf(const EllVar<U>&);
  template <class U>
  friend EllVar<U> erfc(const EllVar<U>&);

 public:
  /** @ingroup ELLIPSOID
   *  @{
   */
  //! @brief Exceptions of mc::EllImg
  class Exceptions
  {
   public:
    //! @brief Enumeration type for EllImg exception handling
    enum TYPE
    {
      DIV = 1,    //!< Division by zero scalar
      INV,        //!< Inverse operation with zero in domain
      LOG,        //!< Log operation with non-positive numbers in domain
      SQRT,       //!< Square-root operation with negative numbers in domain
      TAN,        //!< Tangent operation with zero in domain of cosine, tan(x) =
                  //!< sin(x)/cos(x)
      ACOS,       //!< Cosine inverse operation with domain outside [-1,1]
      ASIN,       //!< Sine inverse operation with domain outside [-1,1]
      INIT = -1,  //!< Failed to construct ellipsoidal variable
      EIMG = -3,  //!< Operation between ellipsoidal variables linked to
                  //!< different ellipsoidal images
      UNDEF = -33  //!< Feature not yet implemented in mc::EllImg
    };
    //! @brief Constructor for error <a>ierr</a>
    Exceptions(TYPE ierr) : _ierr(ierr) {}
    //! @brief Error flag
    int
    ierr()
    {
      return _ierr;
    }
    //! @brief Error description
    std::string
    what()
    {
      switch (_ierr)
      {
        case DIV:
          return "mc::EllImg\t Division by zero scalar";
        case INV:
          return "mc::EllImg\t Inverse operation with zero in domain";
        case LOG:
          return "mc::EllImg\t Log operation with non-positive numbers in "
                 "domain";
        case SQRT:
          return "mc::EllImg\t Square-root operation with negative numbers in "
                 "domain";
        case TAN:
          return "mc::EllImg\t Tangent operation with zero in domain of cosine";
        case ACOS:
          return "mc::EllImg\t Inverse cosine operation with domain outside "
                 "[-1,1]";
        case ASIN:
          return "mc::EllImg\t Inverse sine operation with domain outside "
                 "[-1,1]";
        case EIMG:
          return "mc::EllImg\t Variables belong to different ellipsoidal "
                 "images";
        case UNDEF:
          return "mc::EllImg\t Feature not yet implemented in mc::EllImg class";
        default:
          return "mc::EllImg\t Undocumented error";
      }
    }

   private:
    //! @brief Type of error
    TYPE _ierr;
  };

  //! @brief Structure containing the options for EllImg
  static struct Options
  {
    //! @brief Constructor of mc::EllImg::Options
    Options()
        : PREALLOC(0),
          MINK_TOL(1e-10),
          REMEZ_USE(true),
          REMEZ_MAXIT(5),
          REMEZ_TOL(1e-5),
          REMEZ_MIG(1e-10),
          DCPROD_USE(false)
    {
    }
    //! @brief Copy constructor of mc::EllImg::Options
    template <typename U>
    Options(U& options)
        : PREALLOC(options.PREALLOC),
          MINK_TOL(options.MINK_TOL),
          REMEZ_USE(options.REMEZ_USE),
          REMEZ_MAXIT(options.REMEZ_MAXIT),
          REMEZ_TOL(options.REMEZ_TOL),
          REMEZ_MIG(options.REMEZ_MIG),
          DCPROD_USE(options.DCPROD_USE)
    {
    }
    //! @brief Assignment of mc::EllImg::Options
    template <typename U>
    Options&
    operator=(U& options)
    {
      PREALLOC    = options.PREALLOC;
      MINK_TOL    = options.MINK_TOL;
      REMEZ_USE   = options.REMEZ_USE;
      REMEZ_MAXIT = options.REMEZ_MAXIT;
      REMEZ_TOL   = options.REMEZ_TOL;
      REMEZ_MIG   = options.REMEZ_MIG;
      DCPROD_USE  = options.DCPROD_USE;
      return *this;
    }

    //! @brief Sets number of rows to preallocate in the shape matrix and center
    //! vector (Default: 0)
    long PREALLOC;
    //! @brief Tolerance in minkowski sum (Default: 1e-10)
    double MINK_TOL;
    //! @brief Whether to use a degree-1 Remez linear approximation for
    //! univariate terms (Default: true). If false, use the secant-line
    //! sampled-residual fallback.
    bool REMEZ_USE;
    //! @brief Maximal number of iterations in Remez algorithm for computing a
    //! minimax approximation for univariate terms
    unsigned REMEZ_MAXIT;
    //! @brief Stopping tolerance in Remez algorithm for computing a minimax
    //! approximation for univariate terms
    double REMEZ_TOL;
    //! @brief Minimum interval diameter for invoking Remez; smaller intervals
    //! use the secant-line sampled-residual fallback
    double REMEZ_MIG;
    //! @brief Whether to use DC decomposition to lift bilinear terms (Default:
    //! false)
    bool DCPROD_USE;
  } options;
  /** @} */

 private:
  //! @brief Lifted shape matrix, stored as a hashed dictionary-of-keys over
  //! the lower triangle (row >= col).  The lifting algorithm reads/writes a
  //! handful of entries per elementary operation and depends on them
  //! immediately, so an incremental store with O(1) average get/set is
  //! required.  arma::sp_mat (compressed-column) was previously used with
  //! single-element operator()(r,c) writes -- Armadillo's documented slow
  //! path, which reshuffles the CSC arrays on every insertion -- making
  //! assembly quadratic in the number of stored entries.  We keep the data
  //! in a std::unordered_map and materialise an arma::sp_mat only on demand
  //! (Q_lift / streaming) via efficient batch construction.
  std::unordered_map<unsigned long long, double> _Qmap;
  //! @brief Allocated lifted dimension (row/column count of the shape matrix)
  long _Qn;
  //! @brief Centre of the lifted Ellipsoid
  arma::vec _q;
  //! @brief Dimension of the dependent variables
  long _nx;
  //! @brief Map between pointers to EllVars (key) and row number
  long _curRow;

  //! @brief Pack a lower-triangular (r,c) index pair into a 64-bit hash key
  static unsigned long long
  _qkey(arma::uword const r, arma::uword const c)
  {
    return (static_cast<unsigned long long>(r) << 32) |
           static_cast<unsigned long long>(c);
  }

  //! @brief Canonical lower-triangular sparse index
  static void
  _lower_index(arma::uword& r, arma::uword& c, long const i, long const j)
  {
    if (i >= j)
    {
      r = static_cast<arma::uword>(i);
      c = static_cast<arma::uword>(j);
    }
    else
    {
      r = static_cast<arma::uword>(j);
      c = static_cast<arma::uword>(i);
    }
  }

  //! @brief Symmetric getter from lower-triangular storage (0 if absent)
  double
  qget(long const i, long const j) const
  {
    arma::uword r, c;
    _lower_index(r, c, i, j);
    auto const it = _Qmap.find(_qkey(r, c));
    return it == _Qmap.end() ? 0. : it->second;
  }

  //! @brief Whether a symmetric entry is structurally nonzero
  bool
  qlisted(long const i, long const j) const
  {
    arma::uword r, c;
    _lower_index(r, c, i, j);
    return _Qmap.find(_qkey(r, c)) != _Qmap.end();
  }

  //! @brief Symmetric setter into lower-triangular storage.  Writing an exact
  //! zero erases the entry (so qlisted() stays meaningful) without the
  //! O(nnz) full-matrix clean() that the sp_mat version performed per write.
  void
  qput(long const i, long const j, double const& v)
  {
    arma::uword r, c;
    _lower_index(r, c, i, j);
    unsigned long long const k = _qkey(r, c);
    if (v == 0.)
      _Qmap.erase(k);
    else
      _Qmap[k] = v;
  }

  //! @brief Scale every stored shape entry by a scalar (O(nnz))
  void
  _scaleQ(double const s)
  {
    if (s == 0.)
    {
      _Qmap.clear();
      return;
    }
    for (auto& kv : _Qmap) kv.second *= s;
  }

  //! @brief Materialise the lower-triangular shape as a sparse matrix
  arma::sp_mat
  _spQ() const
  {
    arma::uword const nnz = static_cast<arma::uword>(_Qmap.size());
    arma::uword const N   = static_cast<arma::uword>(_Qn > 0 ? _Qn : 0);
    if (!nnz || !N) return arma::sp_mat(N, N);
    arma::umat loc(2, nnz);
    arma::vec val(nnz);
    arma::uword p = 0;
    for (auto const& kv : _Qmap)
    {
      loc(0, p) = static_cast<arma::uword>(kv.first >> 32);
      loc(1, p) = static_cast<arma::uword>(kv.first & 0xffffffffULL);
      val(p)    = kv.second;
      ++p;
    }
    return arma::sp_mat(loc, val, N, N);  // batch construction
  }

  //! @brief Current lifted dimension
  long
  qdim() const
  {
    return _Qn;
  }

#ifdef MC__ELLIMAGE_DEBUG
  //! @brief Output for debugging
  std::ofstream _dbugout;
#endif

 public:
  /** @ingroup ELLIPSOID
   *  @{
   */
  //! @brief Returns the shape matrix of the lifted Ellipsoid
  arma::sp_mat
  Q_lift()
  {
    return arma::symmatl(_spQ());
  }

  //! @brief Returns the centre of the lifted Ellipsoid
  arma::vec
  c_lift()
  {
    return _q;
  }

  //! @brief Default constructor
  EllImg();

  //! @brief Constructor for ellipsoid with shape matrix \f$Q\f$ and center
  //! \f$c\f$
  EllImg(arma::mat const& Q, arma::vec const& c = arma::vec());

  //! @brief Constructor for ellipsoid of dimension \f$n\f$ with shape matrix
  //! \f$Q\f$ (lower triangular part stored contiguously and columnwise) and
  //! center \f$c\f$
  EllImg(unsigned const n, double const* Q, double const* c = nullptr);

  //! @brief Constructor for ellipsoid enclosing interval vector of radius
  //! \f$r\f$ centered at \f$c\f$
  EllImg(arma::vec const& r, arma::vec const& c = arma::vec());

  //! @brief Copy constructor
  EllImg(EllImg<T> const& E);

  //! @brief Destructor
  virtual ~EllImg();

  //! @brief Set an ellipsoid identical to <a>E</a>
  EllImg<T>&
  set(EllImg<T> const& E)
  {
    Ellipsoid::set(E.Q(), E.c());
    return _reset();
  }

  //! @brief Set an ellipsoid with shape matrix \f$Q\f$ and center \f$c\f$
  EllImg<T>&
  set(arma::mat const& Q, arma::vec const& c = arma::vec())
  {
    Ellipsoid::set(Q, c);
    return _reset();
  }

  //! @brief Set an ellipsoid of dimension \f$n\f$ with shape matrix \f$Q\f$
  //! (lower triangular part stored contiguously and columnwise) and center
  //! \f$c\f$
  EllImg<T>&
  set(unsigned const n, double const* Q, double const* c = 0)
  {
    Ellipsoid::set(n, Q, c);
    return _reset();
  }

  //! @brief Set an ellipsoidal enclosing interval vector of radius \f$r\f$
  //! centered at \f$c\f$
  EllImg<T>&
  set(arma::vec const& r, arma::vec const& c = arma::vec())
  {
    Ellipsoid::set(r, c);
    return _reset();
  }

  //! @brief Reset ellipsoidal image to underlying defining ellipsoid
  EllImg<T>&
  reset()
  {
    return _reset();
  }

  //! @brief Get projection of lifted ellipsoid on variables <a>var</a>
  Ellipsoid
  get(unsigned const nvar, EllVar<T> const* var)
  {
    return _get(nvar, var);
  }
  /** @} */

 private:
  //! @brief Compute the Minkowski sum of the ellipsoid and an interval
  void _minksum(long const, double const&);

  //! @brief Compute the Trace of the Shape matrix for the lifted ellipsoidal
  //! image
  double _trQ();

  //! @brief Set the dependency map for a lifted ellipsoid
  EllImg<T>& _set();

  //! @brief Reset the lifted ellipsoid
  EllImg<T>&
  _reset()
  {
    return _set();
  }

  //! @brief project the lifted ellipsoid
  Ellipsoid _get(unsigned const, EllVar<T> const*);

  //! @brief stretch dimensions of the lifted ellipsoid
  void _stretch(long int);

  //! @brief univariate outer composition
  template <typename PUNIV>
  EllVar<T> _univcompose(long const, T const&, PUNIV const&, T const&);

  //! @brief affine outer composition
  void _affcompose(long const, double const&, long const, double const&,
                   long const, double const&);

  //! @brief univariate minimax linearization
  template <typename PUNIV>
  std::tuple<double, double, double> _minimax(PUNIV const&, T const&);
};

template <typename T>
inline typename EllImg<T>::Options EllImg<T>::options;

//! @brief C++ class for ellipsoidal arithmetic - Ellipsoidal image propagation
////////////////////////////////////////////////////////////////////////
//! mc::EllVar is a C++ class for propagation of ellipsoidal image
//! through a factorable function. The template parameter corresponds
//! to the type used to propagate variable range. Round-off errors are
//! not accounted for in the computations (non-verified implementation).
////////////////////////////////////////////////////////////////////////
template <class T>
class EllVar
////////////////////////////////////////////////////////////////////////
{
  friend class EllImg<T>;

  template <class U>
  friend std::ostream& operator<<(std::ostream&, const EllVar<U>&);
  template <class U>
  friend EllVar<U> inv(const EllVar<U>&);
  template <class U>
  friend EllVar<U> exp(const EllVar<U>&);
  template <class U>
  friend EllVar<U> log(const EllVar<U>&);
  template <class U>
  friend EllVar<U> xlog(const EllVar<U>&);
  template <class U>
  friend EllVar<U> sqrt(const EllVar<U>&);
  template <class U>
  friend EllVar<U> sqr(const EllVar<U>&);
  template <class U>
  friend EllVar<U> pow(const EllVar<U>&, const int);
  template <class U>
  friend EllVar<U> cheb(const EllVar<U>&, const unsigned);
  template <class U>
  friend EllVar<U> cos(const EllVar<U>&);
  template <class U>
  friend EllVar<U> sin(const EllVar<U>&);
  template <class U>
  friend EllVar<U> tan(const EllVar<U>&);
  template <class U>
  friend EllVar<U> acos(const EllVar<U>&);
  template <class U>
  friend EllVar<U> asin(const EllVar<U>&);
  template <class U>
  friend EllVar<U> atan(const EllVar<U>&);
  template <class U>
  friend EllVar<U> cosh(const EllVar<U>&);
  template <class U>
  friend EllVar<U> sinh(const EllVar<U>&);
  template <class U>
  friend EllVar<U> tanh(const EllVar<U>&);
  template <class U>
  friend EllVar<U> erf(const EllVar<U>&);
  template <class U>
  friend EllVar<U> erfc(const EllVar<U>&);

 private:
  //! @brief pointer to underlying lifted ellipsoid
  EllImg<T>* _img;
  //! @brief row index in ellipsoid
  long _ndxRow;
  //! @brief range of the variable
  T _range;

 public:
  /** @ingroup ELLIPSOID
   *  @{
   */
  //! @brief Default constructor
  EllVar();
  //! @brief Copy constructor
  EllVar(EllVar<T> const&);
  //! @brief Constructor for constants
  EllVar(double const& d);
  //! @brief Constructor for intervals
  EllVar(T const& B);
  //! @brief Constructor for intervals
  EllVar(double const& l, double const& u);
  //! @brief Constructor for variable in ellipsoidal image
  EllVar(EllImg<T>&, const unsigned);
  //! @brief Constructor for variable in ellipsoidal image with tailored range
  EllVar(EllImg<T>&, const unsigned, const T&);

  //! @brief Destructor
  virtual ~EllVar() {};

  //! @brief set variable the ellipsoidal image environment
  EllVar<T>&
  set(EllImg<T>& img, const unsigned i)
  {
    return _set(img, i);
  }
  //! @brief set variable in ellipsoidal image environment and tailored range
  EllVar<T>&
  set(EllImg<T>& img, const unsigned i, const T& Irange)
  {
    return _set(img, i, Irange);
  }

  //! @brief get variable range
  T
  range() const
  {
    return _range;
  }
  //! @brief get pointer to ellipsoidal image
  EllImg<T>*
  image() const
  {
    return _img;
  }
  //! @brief get pointer to row index
  long
  index() const
  {
    return _ndxRow;
  }

  //! brief Unit ball in T arithmetic
  static T TOne;
  /** @} */

 private:
  //! @brief constructor for internal variables
  EllVar(EllImg<T>&);

  //! @brief set variable in the lifted ellipsoid
  EllVar<T>& _set(EllImg<T>&, const unsigned);
  //! @brief set variable in the lifted ellipsoid
  EllVar<T>& _set(EllImg<T>&, const unsigned, const T&);

  //! @brief lift ellipsoid and set index
  void _lift(long const ndxRow);

 public:
  // Public overloads
  EllVar<T>& operator=(EllVar<T>&&);
  EllVar<T>& operator=(EllVar<T> const&);
  EllVar<T>& operator=(double const&);
  EllVar<T>& operator=(T const&);

  EllVar<T>& operator+=(EllVar<T> const&);
  EllVar<T>& operator+=(T const&);
  EllVar<T>& operator+=(double const&);
  EllVar<T>& operator-=(EllVar<T> const&);
  EllVar<T>& operator-=(T const&);
  EllVar<T>& operator-=(double const&);
  EllVar<T>& operator*=(EllVar<T> const&);
  EllVar<T>& operator*=(T const&);
  EllVar<T>& operator*=(double const&);
  EllVar<T>& operator/=(EllVar<T> const&);
  EllVar<T>& operator/=(T const&);
  EllVar<T>& operator/=(double const&);
};

template <typename T>
inline T EllVar<T>::TOne = 2. * Op<T>::zeroone() - 1.;

///////////////////////////////// EllImg ///////////////////////////////////

template <class T>
inline EllImg<T>::EllImg() : _nx(0), _curRow(0)
{
#ifdef MC__ELLIMAGE_DEBUG
  _dbugout.open("debug.log", std::ios_base::out);
#endif
}

template <class T>
inline EllImg<T>::EllImg(arma::mat const& Q, arma::vec const& c)
    : Ellipsoid(Q, c)
{
#ifdef MC__ELLIMAGE_DEBUG
  _dbugout.open("debug.log", std::ios_base::out);
#endif
  _set();
}

template <class T>
inline EllImg<T>::EllImg(EllImg<T> const& E)
    : Ellipsoid(E)  //, options( E.options )
{
#ifdef MC__ELLIMAGE_DEBUG
  _dbugout.open("debug.log", std::ios_base::out);
#endif
  _set();
}

template <class T>
inline EllImg<T>::EllImg(unsigned const n, double const* Q, double const* c)
    : Ellipsoid(n, Q, c)
{
#ifdef MC__ELLIMAGE_DEBUG
  _dbugout.open("debug.log", std::ios_base::out);
#endif
  _set();
}

template <class T>
inline EllImg<T>::EllImg(arma::vec const& r, arma::vec const& c)
    : Ellipsoid(r, c)
{
#ifdef MC__ELLIMAGE_DEBUG
  _dbugout.open("debug.log", std::ios_base::out);
#endif
  _set();
}

template <class T>
inline EllImg<T>::~EllImg()
{
#ifdef MC__ELLIMAGE_DEBUG
  _dbugout.close();
#endif
}

template <class T>
inline EllImg<T>&
EllImg<T>::_set()
{
  // set Ellipsoid to  E0(Q0,q0)
  _nx = static_cast<long>(Q().n_cols);
  // Initialise the hashed lower-triangular store from the base ellipsoid.
  _Qmap.clear();
  _Qn = _nx;
  {
    arma::mat const& Q0 = Q();
    for (long j = 0; j < _nx; ++j)
      for (long i = j; i < _nx; ++i)
      {
        double const v = Q0(i, j);
        if (v != 0.)
          _Qmap[_qkey(static_cast<arma::uword>(i),
                      static_cast<arma::uword>(j))] = v;
      }
  }
  _q      = c();
  _curRow = _nx;  // resets _curRow to the number of dependent variables
  return *this;
}

template <class T>
inline Ellipsoid
EllImg<T>::_get(unsigned const nvar, EllVar<T> const* var)
{
  arma::mat Q0(nvar, nvar, arma::fill::zeros);
  arma::vec q0(nvar, arma::fill::zeros);
  for (long j = 0; j < static_cast<long>(nvar); ++j)
  {
    long prev = var[j]._ndxRow;
    q0(j)     = _q(prev);
    for (long i = j; i < static_cast<long>(nvar); ++i)
      Q0(i, j) = qget(var[i]._ndxRow, prev);
  }
  Q0 = arma::symmatl(Q0);
  return Ellipsoid(Q0, q0);
}

template <class T>
inline void
EllImg<T>::_stretch(long const ndxRow)
{
  // Checks if reallocation is needed
  if (_Qn > ndxRow) return;
  // The shape matrix is dictionary-of-keys, so "growing" it is just bumping
  // the logical dimension -- no storage reshuffle.  Only the dense centre
  // vector needs reallocation, which we do geometrically (doubling) to keep
  // the total cost over a sequence of lifts amortised O(1) rather than the
  // previous +1-per-lift behaviour (linear reallocations, each O(nnz)).
  long const inc =
      (options.PREALLOC > 0 ? options.PREALLOC : std::max<long>(_Qn, 1));
  long const nnew = std::max<long>(ndxRow + 1, _Qn + inc);
  _Qn             = nnew;
  _q.resize(static_cast<arma::uword>(nnew));
}

template <typename T>
inline void
EllImg<T>::_affcompose(long const i, double const& coefk, long const k,
                       double const& coefl, long const l, double const& shift)
{
#ifdef MC__ELLIMAGE_DEBUG
  _dbugout << std::scientific << std::setprecision(3) << std::right;
  _dbugout << "affine starts: i= " << i << " k= " << k << " l= " << l
           << std::endl;
  _dbugout << "q \n";
  _dbugout << _q << std::endl;
  _dbugout << "Q \n";
  _dbugout << _spQ() << std::endl;
#endif

  // Update centre
  _q(i) = shift;
  if (k >= 0) _q(i) += coefk * _q(k);
  if (l >= 0) _q(i) += coefl * _q(l);

  // Update shape matrix
  for (long j = 0; j < i; ++j)
  {
    if ((k < 0 || !qlisted(k, j)) && (l < 0 || !qlisted(l, j))) continue;
    double cov = 0.;
    if (k >= 0) cov += coefk * qget(k, j);
    if (l >= 0) cov += coefl * qget(l, j);
    qput(i, j, cov);
  }
  double cov = 0.;
  if (k >= 0) cov += coefk * coefk * qget(k, k);
  if (l >= 0) cov += coefl * coefl * qget(l, l);
  if (k >= 0 && l >= 0) cov += 2 * coefk * coefl * qget(k, l);
  qput(i, i, cov);

#ifdef MC__ELLIMAGE_DEBUG
  _dbugout << std::scientific << std::setprecision(3) << std::right;
  _dbugout << "affine ends: i= " << i << " k= " << k << " l= " << l
           << std::endl;
  _dbugout << "q \n";
  _dbugout << _q << std::endl;
  _dbugout << "Q \n";
  _dbugout << _spQ() << std::endl;
#endif
}

template <class T>
template <typename PUNIV>
inline EllVar<T>
EllImg<T>::_univcompose(long const ndxRow, T const& domain, PUNIV const& f,
                        T const& range)
{
  long i = _curRow;
#ifdef MC__ELLIMAGE_DEBUG
  _dbugout << std::scientific << std::setprecision(3) << std::right;
  _dbugout << "compose starts: curRow= " << i << " ndxRow= " << ndxRow
           << std::endl;
  _dbugout << "q \n";
  _dbugout << _q << std::endl;
  _dbugout << "Q \n";
  _dbugout << _spQ() << std::endl;
#endif
  // Construct new variable
  EllVar<T> var(*this);
  var._lift(i);
  // Minimax linear approximation c0 + c1*x of f over the domain, with
  // linearisation-error radius eta
  auto const& [c0, c1, eta] = _minimax(f, domain);
  // Update centre and shape per linear transformation
  _affcompose(i, c1, ndxRow, 0., -1, c0);
  // Minkowski sum with an interval centered at 0 with radius eta
  _minksum(i, eta);
  // Save range of Variable: intersect the ellipsoidal projection with the
  // analytic range.  If the two enclosures are inconsistent (empty
  // intersection), fall back to the analytic range rather than leaving
  // var._range indeterminate (the inter() result was previously ignored here,
  // unlike in EllVar::_set).
  T const proj = _q(i) + EllVar<T>::TOne * std::sqrt(qget(i, i));
  if (!Op<T>::inter(var._range, proj, range)) var._range = range;
#ifdef MC__ELLIMAGE_DEBUG
  _dbugout << std::scientific << std::setprecision(3) << std::right;
  _dbugout << "compose ends: curRow= " << i << " ndxRow= " << ndxRow
           << std::endl;
  _dbugout << "q \n";
  _dbugout << _q << std::endl;
  _dbugout << "Q \n";
  _dbugout << _spQ() << std::endl;
#endif

  return var;
}

template <class T>
template <typename PUNIV>
inline std::tuple<double, double, double>
EllImg<T>::_minimax(PUNIV const& f, T const& domain)
//( PUNIV const& f ) // range of f assumed as [-1,1]
{
  double const xL = Op<T>::l(domain);
  double const xU = Op<T>::u(domain);
  if (xU <= xL)
  {
    return std::make_tuple(f(xL), 0.0, 0.0);
  }

  // Given a candidate line c0 + c1*x, return a recentred intercept and a
  // linearisation-error radius eta with f(x) in [line(x)-eta, line(x)+eta] for
  // x in [xL,xU].  The residual is sampled on a dense grid; the line is
  // recentred on the midpoint of the residual range (tighter than a one-sided
  // chord band) and eta is inflated by the largest residual change between
  // adjacent nodes.  That inter-sample term is a heuristic Lipschitz-style
  // safeguard: for the smooth univariate terms handled here the residual is
  // monotone between nodes away from its single interior extremum, so the
  // node-to-node variation over-estimates any overshoot in between.  This is
  // not a verified bound (the module is explicitly non-verified), but is
  // materially safer than the previous near-zero (1+10*machprec) inflation,
  // which ignored inter-sample error entirely.
  auto bound_residual = [&](double const c0,
                            double const c1) -> std::tuple<double, double>
  {
    unsigned const NG = 1024;
    double lo = arma::datum::inf, hi = -arma::datum::inf, maxstep = 0.0,
           prev = 0.0;
    for (unsigned k = 0; k <= NG; ++k)
    {
      double const x   = xL + (xU - xL) * double(k) / double(NG);
      double const res = f(x) - (c0 + c1 * x);
      lo               = std::min(lo, res);
      hi               = std::max(hi, res);
      if (k) maxstep = std::max(maxstep, std::fabs(res - prev));
      prev = res;
    }
    double const centre = 0.5 * (lo + hi);
    double eta          = 0.5 * (hi - lo) + maxstep;
    eta *= 1.0 + 10.0 * machprec();
    return std::make_tuple(c0 + centre, eta);
  };

  // Secant (chord) line plus the residual bound above.
  auto secant_fallback = [&]() -> std::tuple<double, double, double>
  {
    double const fL = f(xL), fU = f(xU);
    double const c1       = (fU - fL) / (xU - xL);
    double const c0       = fL - c1 * xL;
    auto const [c0b, eta] = bound_residual(c0, c1);
    return std::make_tuple(c0b, c1, eta);
  };

  if (!options.REMEZ_USE || options.REMEZ_MIG > xU - xL)
    return secant_fallback();

  // Degree-1 Remez exchange adapted from SCModel::_minimax to the present
  // EllImg context.  The approximation is constructed on t in [-1,1] for
  // g(t)=f(m+r*t), then converted back to c0+c1*x in the original variable x.
  double const m = 0.5 * (xL + xU);
  double const r = 0.5 * (xU - xL);
  auto const g   = [&](double const& t) -> double { return f(m + r * t); };

  unsigned const n_points = 3;  // linear polynomial + alternating error
  arma::vec ref(n_points);
  arma::mat mat(n_points, n_points);
  arma::vec rhs(n_points);
  arma::vec sol(n_points);

  for (unsigned i = 0; i < n_points; ++i)
    ref(i) = -std::cos(PI * double(i) / double(n_points - 1));

  double a0 = 0.0, a1 = 0.0;
  double cur_err     = 0.0;
  double max_abs_err = 0.0;
  bool have_solution = false;

  for (unsigned iter = 0; iter < options.REMEZ_MAXIT; ++iter)
  {
    for (unsigned i = 0; i < n_points; ++i)
    {
      double const t = ref(i);
      rhs(i)         = g(t);
      mat(i, 0)      = 1.0;
      mat(i, 1)      = t;
      mat(i, 2)      = (i % 2 == 0) ? 1.0 : -1.0;
    }

    try
    {
      sol           = arma::solve(mat, rhs);
      a0            = sol(0);
      a1            = sol(1);
      cur_err       = sol(2);
      have_solution = true;
    }
    catch (...)
    {
      return secant_fallback();
    }

    unsigned const n_grid = 256;
    double t_at_max       = -1.0;
    max_abs_err           = 0.0;
    for (unsigned k = 0; k <= n_grid; ++k)
    {
      double const t   = -1.0 + 2.0 * double(k) / double(n_grid);
      double const err = g(t) - (a0 + a1 * t);
      if (std::fabs(err) > max_abs_err)
      {
        max_abs_err = std::fabs(err);
        t_at_max    = t;
      }
    }

#ifdef MC__ELLIMAGE_DEBUG_MINIMAX
    std::ostream& _dbugout = std::cout;
    _dbugout << iter << ": [ " << std::right << std::scientific
             << std::setprecision(15) << std::setw(23) << a0 << std::setw(23)
             << a1 << " ] +/- " << std::setw(23) << max_abs_err << std::endl;
#endif

    if (std::fabs(max_abs_err - std::fabs(cur_err)) < options.REMEZ_TOL) break;

    std::vector<double> points;
    points.reserve(n_points + 1);
    for (unsigned i = 0; i < n_points; ++i) points.push_back(ref(i));

    bool duplicate = false;
    for (double const p : points)
    {
      if (std::fabs(p - t_at_max) < 1e-14)
      {
        duplicate = true;
        break;
      }
    }
    if (duplicate) break;

    points.push_back(t_at_max);
    std::sort(points.begin(), points.end());

    std::vector<double> errs(points.size());
    for (unsigned i = 0; i < points.size(); ++i)
      errs[i] = g(points[i]) - (a0 + a1 * points[i]);

    int remove_idx = -1;
    for (unsigned i = 0; i < n_points; ++i)
    {
      double const s1 = (errs[i] >= 0.0) ? 1.0 : -1.0;
      double const s2 = (errs[i + 1] >= 0.0) ? 1.0 : -1.0;
      if (s1 == s2)
      {
        remove_idx = (std::fabs(errs[i]) < std::fabs(errs[i + 1])) ? i : i + 1;
        break;
      }
    }
    if (remove_idx == -1)
      remove_idx =
          (std::fabs(errs[0]) < std::fabs(errs[n_points])) ? 0 : n_points;

    unsigned cnt = 0;
    for (unsigned i = 0; i < n_points + 1; ++i)
    {
      if (int(i) == remove_idx) continue;
      ref(cnt++) = points[i];
    }
  }

  if (!have_solution) return secant_fallback();

  // Convert the minimax line from t-space back to x-space, then bound its
  // residual with the shared (recentred, inter-sample-inflated) estimator.
  double const c1x      = a1 / r;
  double const c0x      = a0 - a1 * m / r;
  auto const [c0b, eta] = bound_residual(c0x, c1x);
  return std::make_tuple(c0b, c1x, eta);
}

template <class T>
inline void
EllImg<T>::_minksum(long const i, double const& rad)
{
#ifdef MC__ELLIMAGE_DEBUG
  std::ostream& _dbugout = std::cout;
  _dbugout << std::scientific << std::setprecision(3) << std::right;
  _dbugout << "Minkowski sum starts: i= " << i << std::endl;
  _dbugout << "q \n";
  _dbugout << _q << std::endl;
  _dbugout << "Q \n";
  _dbugout << _spQ() << std::endl;
#endif

  double EPS = DBL_EPSILON, strQ = 0.0;
  for (long j = 0; j <= i; ++j)
    if (qlisted(j, j)) strQ += qget(j, j);
  strQ        = std::sqrt(strQ) + EPS;
  double strD = rad + EPS;
  _scaleQ(1. + strD / strQ);
  qput(i, i, qget(i, i) + std::pow(rad, 2) * (1. + strQ / strD));

  //  double TOL = options.MINK_TOL, EPS = DBL_EPSILON, strQ = 0.0;
  //  for( long j=0; j<=i ; ++j ){
  //    if( qlisted(j,j) ) strQ += qget(j,j)/(qget(j,j)+TOL); // for some reason
  //    this loops modifies the diagonal elements of the product block if the
  //    isListed is not used ...
  //  }
  //  strQ = std::sqrt(strQ);
  //  double sqrR  =  rad / std::sqrt(qget(i,i)+TOL);
  //  double kappa = strQ + sqrR + EPS;
  //  _Q          *= kappa / (strQ+EPS);
  //  qget(i,i)     += std::pow(rad,2) * kappa / (sqrR+EPS);

#ifdef MC__ELLIMAGE_DEBUG
  _dbugout << std::scientific << std::setprecision(3) << std::right;
  _dbugout << "Minkowski sum ends: i= " << i << std::endl;
  _dbugout << "q \n";
  _dbugout << _q << std::endl;
  _dbugout << "Q \n";
  _dbugout << _spQ() << std::endl;
#endif
}

template <class T>
inline double
EllImg<T>::_trQ()
{
  if (!_Qn) return 0.;
  double tr(qget(0, 0));
  for (long i = 1; i < _Qn; i++)
    tr += qget(static_cast<long>(i), static_cast<long>(i));
  return tr;
}

/////////////////////////////////  EllVar  ///////////////////////////////////

template <class T>
inline EllVar<T>::EllVar() : _img(nullptr), _ndxRow(-1), _range(0.)
{
}

template <class T>
inline EllVar<T>::EllVar(EllVar<T> const& var)
    : _img(var._img), _ndxRow(var._ndxRow), _range(var._range)
{
}

template <class T>
inline EllVar<T>::EllVar(double const& scalar)
    : _img(nullptr), _ndxRow(-1), _range(scalar)
{
}

template <class T>
inline EllVar<T>::EllVar(T const& range)
    : _img(nullptr), _ndxRow(-1), _range(range)
{
}

template <class T>
inline EllVar<T>::EllVar(double const& l, double const& u)
    : _img(nullptr), _ndxRow(-1), _range(l, u)
{
}

template <class T>
inline EllVar<T>::EllVar(EllImg<T>& img, unsigned const i)
{
  _set(img, i);
}

template <class T>
inline EllVar<T>::EllVar(EllImg<T>& img, unsigned const i, T const& range)
{
  _set(img, i, range);
}

template <class T>
inline EllVar<T>&
EllVar<T>::_set(EllImg<T>& img, unsigned const i)
{
  if (i >= img._nx)
    throw typename EllImg<T>::Exceptions(EllImg<T>::Exceptions::INIT);
  _img    = &img;
  _ndxRow = i;
  _range  = img._q(i) + TOne * std::sqrt(img.qget(i, i));
  return *this;
}

template <class T>
inline EllVar<T>&
EllVar<T>::_set(EllImg<T>& img, unsigned const i, T const& range)
{
  if (i >= img._nx)
    throw typename EllImg<T>::Exceptions(EllImg<T>::Exceptions::INIT);
  _img    = &img;
  _ndxRow = i;
  if (!Op<T>::inter(_range, img._q(i) + TOne * std::sqrt(img.qget(i, i)),
                    range))
    throw typename EllImg<T>::Exceptions(EllImg<T>::Exceptions::INIT);
  return *this;
}

template <class T>
inline EllVar<T>::EllVar(EllImg<T>& img) : _img(&img), _ndxRow(-1), _range(0.)
{
}

template <class T>
inline void
EllVar<T>::_lift(long const ndxRow)
{
  if (!_img) throw typename EllImg<T>::Exceptions(EllImg<T>::Exceptions::INIT);
  _ndxRow = ndxRow;
  _img->_stretch(ndxRow);
  _img->_curRow++;
}

// template <class T>
// inline
// void
// EllVar<T>::_lift
//( long const ndxRow, long const curRow )
//{
//   if( !_img ) throw typename EllImg<T>::Exceptions(
//   EllImg<T>::Exceptions::INIT ); _ndxRow = ndxRow; _img->_stretch( ndxRow );
//   _img->_curRow = curRow;
// }

template <class T>
inline EllVar<T>&
EllVar<T>::operator=(EllVar<T> const& var)
{
  _img    = var._img;
  _ndxRow = var._ndxRow;
  _range  = var._range;
  return *this;
}

template <class T>
inline EllVar<T>&
EllVar<T>::operator=(EllVar<T>&& var)
{
  std::swap(_img, var._img);
  std::swap(_ndxRow, var._ndxRow);
  std::swap(_range, var._range);
  return *this;
}

template <class T>
inline EllVar<T>&
EllVar<T>::operator=(double const& scalar)
{
  _img    = nullptr;
  _ndxRow = -1;
  _range  = scalar;
  return *this;
}

template <class T>
inline EllVar<T>&
EllVar<T>::operator=(T const& range)
{
  _img    = nullptr;
  _ndxRow = -1;
  _range  = range;
  return *this;
}

template <typename T>
inline EllVar<T>&
EllVar<T>::operator+=(double const& scalar)
{
  // Current variable is a range
  if (!_img)
  {
    _range += scalar;
    return *this;
  }

  // Lift ellipsoid
  long i = _img->_curRow;
  long k = _ndxRow;
  _lift(i);
  _img->_affcompose(i, 1., k, 0., -1, scalar);

  // Set variable range
  _range = _img->_q(i) + TOne * std::sqrt(_img->qget(i, i));

  return *this;
}

template <typename T>
inline EllVar<T>&
EllVar<T>::operator-=(double const& scalar)
{
  return operator+=(-scalar);
}

template <typename T>
inline EllVar<T>&
EllVar<T>::operator+=(T const& range)
{
  // Current variable is a range
  if (!_img)
  {
    _range += range;
    return *this;
  }

  // Lift ellipsoid
  long i = _img->_curRow;
  operator+=(Op<T>::mid(range));
  double rad = 0.5 * Op<T>::diam(range);
  if (rad < DBL_EPSILON) return *this;
  _img->_minksum(i, rad);
  _range = _img->_q(i) + TOne * std::sqrt(_img->qget(i, i));

  return *this;
}

template <typename T>
inline EllVar<T>&
EllVar<T>::operator-=(T const& range)
{
  return operator+=(-range);
}

template <typename T>
inline EllVar<T>&
EllVar<T>::operator+=(EllVar<T> const& var)
{
  // The operand isn't an ellipsoidal variables
  if (!var._img) return operator+=(var._range);

  // *this isn't an ellipsoidal variables
  if (!_img)
  {
    T range = _range;
    *this   = var;
    return operator+=(range);
  }

  // Both operands correspond to different ellipsoids
  if (_img != var._img)
    throw typename EllImg<T>::Exceptions(EllImg<T>::Exceptions::EIMG);

  // Lift ellipsoid
  long i = _img->_curRow;
  long k = _ndxRow;
  long l = var._ndxRow;
  _lift(i);
  _img->_affcompose(i, 1., k, 1., l, 0.);

  // Set variable range
  _range = _img->_q(i) + TOne * std::sqrt(_img->qget(i, i));

  return *this;
}

template <typename T>
inline EllVar<T>&
EllVar<T>::operator-=(EllVar<T> const& var)
{
  // The operand isn't an ellipsoidal variables
  if (!var._img) return operator-=(var._range);

  // *this isn't an ellipsoidal variables
  if (!_img)
  {
    T range = _range;
    *this   = var;

    // Lift ellipsoid
    long i = _img->_curRow;
    long k = _ndxRow;
    _lift(i);
    _img->_affcompose(i, -1., k, 0., -1, Op<T>::mid(range));
    double rad = 0.5 * Op<T>::diam(range);
    if (rad > DBL_EPSILON) _img->_minksum(i, rad);
    _range = _img->_q(i) + TOne * std::sqrt(_img->qget(i, i));
    return *this;
  }

  // Both operands correspond to different ellipsoids
  if (_img != var._img)
    throw typename EllImg<T>::Exceptions(EllImg<T>::Exceptions::EIMG);

  // Lift ellipsoid
  long i = _img->_curRow;
  long k = _ndxRow;
  long l = var._ndxRow;
  _lift(i);
  _img->_affcompose(i, 1., k, -1., l, 0.);

  // Set variable range
  _range = _img->_q(i) + TOne * std::sqrt(_img->qget(i, i));

  return *this;
}

template <typename T>
inline EllVar<T>&
EllVar<T>::operator*=(double const& scalar)
{
  // Current variable is a range
  if (!_img)
  {
    _range *= scalar;
    return *this;
  }

  // Lift ellipsoid
  long i = _img->_curRow;
  long k = _ndxRow;
  _lift(i);
  _img->_affcompose(i, scalar, k, 0., -1, 0.);

  // Set variable range
  _range = _img->_q(i) + TOne * std::sqrt(_img->qget(i, i));

  return *this;
}

template <typename T>
inline EllVar<T>&
EllVar<T>::operator*=(T const& range)
{
  double rad = 0.5 * Op<T>::diam(range);
  if (rad > DBL_EPSILON)
    throw typename EllImg<T>::Exceptions(EllImg<T>::Exceptions::UNDEF);

  return operator*=(Op<T>::mid(range));
}

template <typename T>
inline EllVar<T>&
EllVar<T>::operator*=(EllVar<T> const& var)
{
  // The operand isn't an ellipsoidal variables
  if (!var._img) return operator*=(var._range);

  // *this isn't an ellipsoidal variables
  if (!_img)
  {
    T range = _range;
    *this   = var;
    return operator*=(range);
  }

  // Both operands correspond to different ellipsoids
  if (_img != var._img)
    throw typename EllImg<T>::Exceptions(EllImg<T>::Exceptions::EIMG);

  // Lift ellipsoid
  long i = _img->_curRow;
  long k = _ndxRow;
  long l = var._ndxRow;
  _lift(i);
  double coefk = Op<T>::mid(var._range);
  double coefl = Op<T>::mid(_range);
  _img->_affcompose(i, coefk, k, coefl, l, -coefk * coefl);
  _img->_minksum(i, Op<T>::diam(_range) * Op<T>::diam(var._range) / 4.);

  // Set variable range
  _range = _img->_q(i) + TOne * std::sqrt(_img->qget(i, i));

  return *this;
}

template <typename T>
inline EllVar<T>&
EllVar<T>::operator/=(double const& scalar)
{
  if (scalar == 0.)
    throw typename EllImg<T>::Exceptions(EllImg<T>::Exceptions::DIV);
  return operator*=(1. / scalar);
}

template <typename T>
inline EllVar<T>&
EllVar<T>::operator/=(T const& range)
{
  return operator*=(Op<T>::inv(range));
}

template <typename T>
inline EllVar<T>&
EllVar<T>::operator/=(EllVar<T> const& var)
{
  return operator*=(inv(var));
}

/////////////////////////////////  Operators ///////////////////////////////////

template <class T>
inline EllVar<T>
operator+(EllVar<T> const& var)
{
  return var;
}

template <class T>
inline EllVar<T>
operator+(EllVar<T> const& var1, EllVar<T> const& var2)
{
  EllVar<T> var3(var1);
  var3 += var2;
  return var3;
}

template <class T>
inline EllVar<T>
operator+(EllVar<T> const& var1, double const& shift2)
{
  EllVar<T> var3(var1);
  var3 += shift2;
  return var3;
}

template <class T>
inline EllVar<T>
operator+(double const& shift1, EllVar<T> const& var2)
{
  EllVar<T> var3(var2);
  var3 += shift1;
  return var3;
}

template <class T>
inline EllVar<T>
operator+(EllVar<T> const& var1, T const& range2)
{
  EllVar<T> var3(var1);
  var3 += range2;
  return var3;
}

template <class T>
inline EllVar<T>
operator+(T const& range1, EllVar<T> const& var2)
{
  EllVar<T> var3(var2);
  var3 += range1;
  return var3;
}

template <class T>
inline EllVar<T>
operator-(EllVar<T> const& var1)
{
  EllVar<T> var3(var1);
  var3 *= -1.;
  return var3;
}

template <class T>
inline EllVar<T>
operator-(EllVar<T> const& var1, EllVar<T> const& var2)
{
  EllVar<T> var3(var1);
  var3 -= var2;
  return var3;
}

template <class T>
inline EllVar<T>
operator-(EllVar<T> const& var1, double const& shift2)
{
  EllVar<T> var3(var1);
  var3 -= shift2;
  return var3;
}

template <class T>
inline EllVar<T>
operator-(double const& shift1, EllVar<T> const& var2)
{
  EllVar<T> var3(shift1);
  var3 -= var2;
  return var3;
}

template <class T>
inline EllVar<T>
operator-(EllVar<T> const& var1, T const& range2)
{
  EllVar<T> var3(var1);
  var3 -= range2;
  return var3;
}

template <class T>
inline EllVar<T>
operator-(T const& range1, EllVar<T> const& var2)
{
  EllVar<T> var3(range1);
  var3 -= var2;
  return var3;
}

template <class T>
inline EllVar<T>
operator*(EllVar<T> const& var1, EllVar<T> const& var2)
{
  if (EllImg<T>::options.DCPROD_USE)
    return (sqr(var1 + var2) - sqr(var1 - var2)) / 4.;

  EllVar<T> var3(var1);
  var3 *= var2;
  return var3;
}

template <class T>
inline EllVar<T>
operator*(EllVar<T> const& var1, double const& shift2)
{
  EllVar<T> var3(var1);
  var3 *= shift2;
  return var3;
}

template <class T>
inline EllVar<T>
operator*(double const& shift1, EllVar<T> const& var2)
{
  EllVar<T> var3(var2);
  var3 *= shift1;
  return var3;
}

template <class T>
inline EllVar<T>
operator*(EllVar<T> const& var1, T const& range2)
{
  EllVar<T> var3(var1);
  var3 *= range2;
  return var3;
}

template <class T>
inline EllVar<T>
operator*(T const& range1, EllVar<T> const& var2)
{
  EllVar<T> var3(var2);
  var3 *= range1;
  return var3;
}

template <class T>
inline EllVar<T>
operator/(EllVar<T> const& var1, EllVar<T> const& var2)
{
  return var1 * inv(var2);
}

template <class T>
inline EllVar<T>
operator/(EllVar<T> const& var1, double const& shift2)
{
  EllVar<T> var3(var1);
  var3 /= shift2;
  return var3;
}

template <class T>
inline EllVar<T>
operator/(double const& shift1, EllVar<T> const& var2)
{
  EllVar<T> var3(shift1);
  var3 /= var2;
  return var3;
}

template <class T>
inline EllVar<T>
operator/(EllVar<T> const& var1, T const& range2)
{
  EllVar<T> var3(var1);
  var3 /= range2;
  return var3;
}

template <class T>
inline EllVar<T>
operator/(T const& range1, EllVar<T> const& var2)
{
  EllVar<T> var3(range1);
  var3 /= var2;
  return var3;
}

template <class T>
inline EllVar<T>
inv(EllVar<T> const& var1)
{
  // var1 is a range
  if (var1._img == nullptr) return Op<T>::inv(var1._range);

  // var1 is a variable
  T const& domain = var1._range;
  if (Op<T>::l(domain) <= 0. && Op<T>::u(domain) >= 0.)
    throw typename EllImg<T>::Exceptions(EllImg<T>::Exceptions::INV);
  return var1._img->_univcompose(
      var1._ndxRow, domain, [=](const double& x) { return 1. / x; },
      Op<T>::inv(domain));
}

template <class T>
inline EllVar<T>
exp(EllVar<T> const& var1)
{
  // var1 is a range
  if (var1._img == nullptr) return Op<T>::exp(var1._range);

  // var1 is a variable
  T const& domain = var1._range;
  return var1._img->_univcompose(
      var1._ndxRow, domain, [=](const double& x) { return std::exp(x); },
      Op<T>::exp(domain));
}

template <class T>
inline EllVar<T>
log(EllVar<T> const& var1)
{
  // var1 is a range
  if (var1._img == nullptr) return Op<T>::log(var1._range);

  // var1 is a variable
  T const& domain = var1._range;
  if (Op<T>::l(domain) <= 0.)
    throw typename EllImg<T>::Exceptions(EllImg<T>::Exceptions::LOG);
  return var1._img->_univcompose(
      var1._ndxRow, domain, [=](const double& x) { return std::log(x); },
      Op<T>::log(domain));
}

template <class T>
inline EllVar<T>
xlog(EllVar<T> const& var1)
{
  // var1 is a range
  if (var1._img == nullptr) return Op<T>::xlog(var1._range);

  // var1 is a variable
  T const& domain = var1._range;
  if (Op<T>::l(domain) < 0.)
    throw typename EllImg<T>::Exceptions(EllImg<T>::Exceptions::LOG);
  return var1._img->_univcompose(
      var1._ndxRow, domain, [=](const double& x) { return mc::xlog(x); },
      Op<T>::xlog(domain));
}

template <class T>
inline EllVar<T>
sqrt(EllVar<T> const& var1)
{
  // var1 is a range
  if (var1._img == nullptr) return Op<T>::sqrt(var1._range);

  // var1 is a variable
  T const& domain = var1._range;
  if (Op<T>::l(domain) < 0.)
    throw typename EllImg<T>::Exceptions(EllImg<T>::Exceptions::SQRT);
  return var1._img->_univcompose(
      var1._ndxRow, domain, [=](const double& x) { return std::sqrt(x); },
      Op<T>::sqrt(domain));
}

template <class T>
inline EllVar<T>
sqr(EllVar<T> const& var1)
{
  // var1 is a range
  if (var1._img == nullptr) return Op<T>::sqr(var1._range);

  // var1 is a variable
  T const& domain = var1._range;
  return var1._img->_univcompose(
      var1._ndxRow, domain, [=](const double& x) { return mc::sqr(x); },
      Op<T>::sqr(domain));
}

template <class T>
inline EllVar<T>
pow(EllVar<T> const& var1, int const n)
{
  // special cases
  if (!n) return 1.;
  if (n == 1) return var1;

  // var1 is a range
  if (var1._img == nullptr) return Op<T>::pow(var1._range, n);

  // var1 is a variable
  T const& domain = var1._range;
  // x^n with n<0 is singular iff 0 lies in the domain.  Test that directly
  // (as inv() above does) rather than via 0 in inv(domain): the latter relies
  // on Op<T>::inv of a zero-straddling interval being defined and itself
  // containing 0, which is interval-library dependent and may throw or misfire.
  if (n < 0 && Op<T>::l(domain) <= 0. && Op<T>::u(domain) >= 0.)
    throw typename EllImg<T>::Exceptions(EllImg<T>::Exceptions::INV);
  return var1._img->_univcompose(
      var1._ndxRow, domain, [=](const double& x) { return std::pow(x, n); },
      Op<T>::pow(domain, n));
}

template <class T>
inline EllVar<T>
cheb(EllVar<T> const& var1, unsigned const n)
{
  // special cases
  if (!n) return 1.;
  if (n == 1) return var1;

  // var1 is a range
  if (var1._img == nullptr) return Op<T>::cheb(var1._range, n);

  // var1 is a variable
  T const& domain = var1._range;
  return var1._img->_univcompose(
      var1._ndxRow, domain, [=](const double& x) { return mc::cheb(x, n); },
      Op<T>::cheb(domain, n));
}

template <class T>
inline EllVar<T>
cos(EllVar<T> const& var1)
{
  // var1 is a range
  if (var1._img == nullptr) return Op<T>::cos(var1._range);

  // var1 is a variable
  T const& domain = var1._range;
  return var1._img->_univcompose(
      var1._ndxRow, domain, [=](const double& x) { return std::cos(x); },
      Op<T>::cos(domain));
}

template <class T>
inline EllVar<T>
sin(EllVar<T> const& var1)
{
  // var1 is a range
  if (var1._img == nullptr) return Op<T>::sin(var1._range);

  // var1 is a variable
  T const& domain = var1._range;
  return var1._img->_univcompose(
      var1._ndxRow, domain, [=](const double& x) { return std::sin(x); },
      Op<T>::sin(domain));
}

template <class T>
inline EllVar<T>
tan(EllVar<T> const& var1)
{
  // var1 is a range
  if (var1._img == nullptr) return Op<T>::tan(var1._range);

  // var1 is a variable
  T const& domain = var1._range;
  if (Op<T>::l(cos(domain)) <= 0. && Op<T>::u(cos(domain)) >= 0.)
    throw typename EllImg<T>::Exceptions(EllImg<T>::Exceptions::TAN);
  return var1._img->_univcompose(
      var1._ndxRow, domain, [=](const double& x) { return std::tan(x); },
      Op<T>::tan(domain));
}

template <class T>
inline EllVar<T>
acos(EllVar<T> const& var1)
{
  // var1 is a range
  if (var1._img == nullptr) return Op<T>::acos(var1._range);

  // var1 is a variable
  T const& domain = var1._range;
  if (Op<T>::l(domain) < -1. || Op<T>::u(domain) > 1.)
    throw typename EllImg<T>::Exceptions(EllImg<T>::Exceptions::ACOS);
  return var1._img->_univcompose(
      var1._ndxRow, domain, [=](const double& x) { return std::acos(x); },
      Op<T>::acos(domain));
}

template <class T>
inline EllVar<T>
asin(EllVar<T> const& var1)
{
  // var1 is a range
  if (var1._img == nullptr) return Op<T>::asin(var1._range);

  // var1 is a variable
  T const& domain = var1._range;
  if (Op<T>::l(domain) < -1. || Op<T>::u(domain) > 1.)
    throw typename EllImg<T>::Exceptions(EllImg<T>::Exceptions::ASIN);
  return var1._img->_univcompose(
      var1._ndxRow, domain, [=](const double& x) { return std::asin(x); },
      Op<T>::asin(domain));
}

template <class T>
inline EllVar<T>
atan(EllVar<T> const& var1)
{
  // var1 is a range
  if (var1._img == nullptr) return Op<T>::atan(var1._range);

  // var1 is a variable
  T const& domain = var1._range;
  return var1._img->_univcompose(
      var1._ndxRow, domain, [=](const double& x) { return std::atan(x); },
      Op<T>::atan(domain));
}

template <class T>
inline EllVar<T>
cosh(EllVar<T> const& var1)
{
  // var1 is a range
  if (var1._img == nullptr) return Op<T>::cosh(var1._range);

  // var1 is a variable
  T const& domain = var1._range;
  return var1._img->_univcompose(
      var1._ndxRow, domain, [=](const double& x) { return std::cosh(x); },
      Op<T>::cosh(domain));
}

template <class T>
inline EllVar<T>
sinh(EllVar<T> const& var1)
{
  // var1 is a range
  if (var1._img == nullptr) return Op<T>::sinh(var1._range);

  // var1 is a variable
  T const& domain = var1._range;
  return var1._img->_univcompose(
      var1._ndxRow, domain, [=](const double& x) { return std::sinh(x); },
      Op<T>::sinh(domain));
}

template <class T>
inline EllVar<T>
tanh(EllVar<T> const& var1)
{
  // var1 is a range
  if (var1._img == nullptr) return Op<T>::tanh(var1._range);

  // var1 is a variable
  T const& domain = var1._range;
  return var1._img->_univcompose(
      var1._ndxRow, domain, [=](const double& x) { return std::tanh(x); },
      Op<T>::tanh(domain));
}

template <class T>
inline EllVar<T>
erf(EllVar<T> const& var1)
{
  // var1 is a range
  if (var1._img == nullptr) return Op<T>::erf(var1._range);

  // var1 is a variable
  T const& domain = var1._range;
  return var1._img->_univcompose(
      var1._ndxRow, domain, [=](const double& x) { return std::erf(x); },
      Op<T>::erf(domain));
}

template <class T>
inline EllVar<T>
erfc(EllVar<T> const& var1)
{
  // var1 is a range
  if (var1._img == nullptr) return Op<T>::erfc(var1._range);

  // var1 is a variable
  T const& domain = var1._range;
  return var1._img->_univcompose(
      var1._ndxRow, domain, [=](const double& x) { return std::erfc(x); },
      Op<T>::erfc(domain));
}

template <typename T>
inline std::ostream&
operator<<(std::ostream& os, EllVar<T> const& var)
{
  return os << std::scientific << var._range << "(index: " << var._ndxRow
            << ")";
}

template <class T>
inline std::ostream&
operator<<(std::ostream& os, EllImg<T> const& img)
{
  const int iprec = 5;
  return os << std::scientific << std::setprecision(iprec) << "\ncenter:\n"
            << img._q << "shape:\n"
            << img._spQ();
  // return os << static_cast<const Ellipsoid&>( img );
}

}  // namespace mc

namespace mc
{

//! @brief Specialization of the structure mc::Op to allow usage of the type
//! mc::Interval for DAG evaluation or as a template parameter in other MC++
//! classes
template <typename T>
struct Op<mc::EllVar<T> >
{
  typedef mc::EllVar<T> EV;
  static EV
  point(const double c)
  {
    return EV(c);
  }
  static EV
  zeroone()
  {
    return EV(mc::Op<T>::zeroone());
  }  // ??
  static void
  I(EV& x, const EV& y)
  {
    x = y;
  }  // ??
  static double
  l(const EV& x)
  {
    return mc::Op<T>::l(x.range());
  }
  static double
  u(const EV& x)
  {
    return mc::Op<T>::u(x.range());
  }
  static double
  abs(const EV& x)
  {
    return mc::Op<T>::abs(x.range());
  }
  static double
  mid(const EV& x)
  {
    return mc::Op<T>::mid(x.range());
  }
  static double
  diam(const EV& x)
  {
    return mc::Op<T>::diam(x.range());
  }
  static EV
  inv(const EV& x)
  {
    return mc::inv(x);
  }
  static EV
  sqr(const EV& x)
  {
    return mc::sqr(x);
  }
  static EV
  sqrt(const EV& x)
  {
    return mc::sqrt(x);
  }
  static EV
  log(const EV& x)
  {
    return mc::log(x);
  }
  static EV
  xlog(const EV& x)
  {
    return mc::xlog(x);
  }
  static EV
  lmtd(const EV& x, const EV& y)
  {
    return (x - y) / (mc::log(x) - mc::log(y));
  }
  static EV
  rlmtd(const EV& x, const EV& y)
  {
    return (mc::log(x) - mc::log(y)) / (x - y);
  }
  static EV
  fabs(const EV& x)
  {
    throw typename mc::EllImg<T>::Exceptions(EllImg<T>::Exceptions::UNDEF);
  }
  static EV
  exp(const EV& x)
  {
    return mc::exp(x);
  }
  static EV
  sin(const EV& x)
  {
    return mc::sin(x);
  }
  static EV
  cos(const EV& x)
  {
    return mc::cos(x);
  }
  static EV
  tan(const EV& x)
  {
    return mc::tan(x);
  }
  static EV
  asin(const EV& x)
  {
    return mc::asin(x);
  }
  static EV
  acos(const EV& x)
  {
    return mc::acos(x);
  }
  static EV
  atan(const EV& x)
  {
    return mc::atan(x);
  }
  static EV
  sinh(const EV& x)
  {
    return mc::sinh(x);
  }
  static EV
  cosh(const EV& x)
  {
    return mc::cosh(x);
  }
  static EV
  tanh(const EV& x)
  {
    return mc::tanh(x);
  }
  static EV
  erf(const EV& x)
  {
    return mc::erf(x);
  }
  static EV
  erfc(const EV& x)
  {
    return mc::erfc(x);
  }
  static EV
  fstep(const EV& x)
  {
    throw typename mc::EllImg<T>::Exceptions(EllImg<T>::Exceptions::UNDEF);
  }
  static EV
  bstep(const EV& x)
  {
    throw typename mc::EllImg<T>::Exceptions(EllImg<T>::Exceptions::UNDEF);
  }
  static EV
  min(const EV& x, const EV& y)
  {
    throw typename mc::EllImg<T>::Exceptions(EllImg<T>::Exceptions::UNDEF);
  }
  static EV
  max(const EV& x, const EV& y)
  {
    throw typename mc::EllImg<T>::Exceptions(EllImg<T>::Exceptions::UNDEF);
  }
  static EV
  arh(const EV& x, const double k)
  {
    return mc::exp(-k / x);
  }
  static EV
  pow(const EV& x, const unsigned n)
  {
    return mc::pow(x, n);
  }
  template <typename X, typename Y>
  static EV
  pow(const X& x, const Y& y)
  {
    return mc::exp(y * mc::log(x));
  }
  static EV
  cheb(const EV& x, const unsigned n)
  {
    return mc::cheb(x, n);
  }
  static EV
  prod(const unsigned int n, const EV* x)
  {
    throw typename mc::EllImg<T>::Exceptions(EllImg<T>::Exceptions::UNDEF);
  }
  static EV
  monom(const unsigned int n, const EV* x, const int* k)
  {
    throw typename mc::EllImg<T>::Exceptions(EllImg<T>::Exceptions::UNDEF);
  }
  static EV
  hull(const EV& x, const EV& y)
  {
    throw typename mc::EllImg<T>::Exceptions(EllImg<T>::Exceptions::UNDEF);
  }
  static bool
  inter(EV& xIy, const EV& x, const EV& y)
  {
    throw typename mc::EllImg<T>::Exceptions(EllImg<T>::Exceptions::UNDEF);
  }
  static bool
  eq(const EV& x, const EV& y)
  {
    return x.env() == y.env() && x.index() == y.index() &&
           x.range() == y.range();
  }
  static bool
  ne(const EV& x, const EV& y)
  {
    return !eq(x, y);
  }
  static bool
  lt(const EV& x, const EV& y)
  {
    throw typename mc::EllImg<T>::Exceptions(EllImg<T>::Exceptions::UNDEF);
  }
  static bool
  le(const EV& x, const EV& y)
  {
    throw typename mc::EllImg<T>::Exceptions(EllImg<T>::Exceptions::UNDEF);
  }
  static bool
  gt(const EV& x, const EV& y)
  {
    throw typename mc::EllImg<T>::Exceptions(EllImg<T>::Exceptions::UNDEF);
  }
  static bool
  ge(const EV& x, const EV& y)
  {
    throw typename mc::EllImg<T>::Exceptions(EllImg<T>::Exceptions::UNDEF);
  }
};

}  // namespace mc

#endif
