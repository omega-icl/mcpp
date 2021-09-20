// Copyright (C) 2009-2018 Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the Eclipse Public License.

/*!
\page page_POLYHEDRAL Polyhedral Arithmetic for Factorable Functions
\author Beno&icirc;t Chachuat

Consider the image set \f$\Gamma:\{{\bf g}({\bf x}) \,\mid\, {\bf x}^{\rm L}\leq{\bf x}\leq{\bf x}^{\rm U}\} \subseteq \mathbb R^m\f$, where \f$g_j:\mathbb{R}^n\to\mathbb R, j=1,\ldots,m\f$, are factorable, potentially nonconvex functions. The classes mc::PolImg and mc::PolVar implement an arithmetic for the construction of a polyedral enclosure \f$\Gamma\f$ of the set \f$\overline{\Gamma}\f$ in the form:
\f{align*}
\overline{\Gamma} := \left\{ {\bf G} {\bf v}\; \middle|\; \begin{array}{l} {\bf A} {\bf v} \;=\; {\bf 0}\\ {\bf B} {\bf v} \;\leq\; {\bf 0}\\ {\bf v}^{\rm L}\leq{\bf v}\leq{\bf v}^{\rm U} \end{array} \right\} \supseteq \Gamma\, ,
\f}
where the variable vector \f${\bf v}\f$ contains the original variable vector \f${\bf x}\f$, as recovered by using the projection matrix \f${\bf G}\f$.

These classes build upon the DAG classes mc::FFGraph and mc::FFVar. Both are templated in the interval type used to bound the nonlinearity of the function, By default, mc::PolImg and mc::PolVar can be used with the non-verified interval type mc::Interval of MC++. For reliability, however, it is recommended to use verified interval arithmetic such as <A href="http://www.ti3.tu-harburg.de/Software/PROFILEnglisch.html">PROFIL</A> (header file <tt>mcprofil.hpp</tt>) or <A href="http://www.math.uni-wuppertal.de/~xsc/software/filib.html">FILIB++</A> (header file <tt>mcfilib.hpp</tt>). Note that the implementation in mc::PolImg and mc::PolVar is <a>not verified</a> in the sense that rounding errors are not accounted for in the polyhedral cuts.


\section sec_POLIMG_THEOR What is the theory behind polyhedral relaxations?

The basic procedure follows three steps: (i) decompose the factorable functions \f$g_j\f$ into atom operations, both unary and binary operations, by introducing auxiliary variables and extra constraints (lifting); (ii) generate convex enclosures for the nonlinear operations, including products, divisions and outer-compositions with nonlinear univariates such as exp, log, sin, cos, <I>etc</I>; and (iii) outer-approximate the nonlinear parts of the convex enclosures using supporting hyperplanes at finitely many points. Each step is detailed further below.

- <B>Step 1. Decomposition/Lifting</B>
\f{align*}
\Gamma:\{{\bf g}({\bf x}) \,\mid\, {\bf x}^{\rm L}\leq{\bf x}\leq {\bf x}^{\rm U}\}
& \xrightarrow{\displaystyle\text{decomp.}} && 
\overline{\Gamma} := \left\{ {\bf G} {\bf v}\; \middle|\; \begin{array}{l} {\bf A} {\bf v} \;=\; {\bf 0}\\ v_k = v_iv_j,\ \forall (i,j,k)\in\mathcal{B}\\ v_k = \frac{v_i}{v_j},\ \forall (i,j,k)\in\mathcal{F}\\ v_k = \varphi(v_i),\ \forall (i,k)\in\mathcal{U}\\ {\bf v}^{\rm L}\leq{\bf v}\leq {\bf v}^{\rm U} \end{array} \right\}
\f}
The advantage of this decomposition, which may be applied to any factorable function, is that it accounts for common subexpressions, thus enables tighter relaxations. Its main drawback, on the other hand, is that it may introduce a large number of auxiliary variables and constraints.

- <B>Step 2. Relaxation</B>
  - The bilinear terms \f$v_k = v_iv_j\f$, \f$v_i^{\rm L}\leq v_i\leq v_i^{\rm U}\f$, \f$v_j^{\rm L}\leq v_j\leq v_j^{\rm U}\f$, can be replaced by their polyhedral envelopes as
\f{align*}
v_k=v_iv_j \quad \xrightarrow{\displaystyle\text{relax.}}\quad \left\{\begin{array}{l}
v_k \geq v_i^{\rm L}v_j+v_j^{\rm L}v_i-v_i^{\rm L}v_j^{\rm L}\\
v_k \geq v_i^{\rm U}v_j+v_j^{\rm U}v_i-v_i^{\rm U}v_j^{\rm U}\\
v_k \leq v_i^{\rm U}v_j+v_j^{\rm L}v_i-v_i^{\rm U}v_j^{\rm L}\\
v_k \leq v_i^{\rm L}v_j+v_j^{\rm U}v_i-v_i^{\rm L}v_j^{\rm U}
\end{array}\right.
\f}
  - The fractional terms \f$v_k = \frac{v_i}{v_j}\f$, \f$v_i^{\rm L}\leq v_i\leq v_i^{\rm U}\f$, \f$v_j^{\rm L}\leq v_j\leq v_j^{\rm U}\f$, \f$0\notin[v_j^{\rm L},v_j^{\rm U}]\f$, can be rewritten as bilinear terms \f$v_i = v_jv_k\f$ and relaxed as indicated above, with the following bounds for the variables \f$v_k\f$:
\f{align*}
\min\left\{\frac{v_j^{\rm L}}{v_k^{\rm L}}, \frac{v_j^{\rm L}}{v_k^{\rm U}}, \frac{v_j^{\rm U}}{v_k^{\rm L}}, \frac{v_j^{\rm U}}{v_k^{\rm U}}\right\} =: v_k^{\rm L} \leq v_k \leq v_k^{\rm U} := \max\left\{\frac{v_j^{\rm L}}{v_k^{\rm L}}, \frac{v_j^{\rm L}}{v_k^{\rm U}}, \frac{v_j^{\rm U}}{v_k^{\rm L}}, \frac{v_j^{\rm U}}{v_k^{\rm U}}\right\}
\f}
Note that, although straightforward, this approach of dealing with fractional terms does not generally yield the convex/concave envelopes (which turn out to be quite complicated nonlinear expressions -- see \ref sec_POLIMG_REFS).
\n
  - The univariate terms \f$v_k = \varphi(v_i)\f$, \f$v_i^{\rm L}\leq v_i\leq v_i^{\rm U}\f$, are relaxed differently depending on whether the function \f$\varphi\f$ is convex, concave, convexo-concave, etc., on \f$[v_i^{\rm L},v_i^{\rm U}]\f$.
\f{align*}
\text{convex case:} & \quad v_k=\varphi(v_i) \quad \xrightarrow{\displaystyle\text{relax}} \quad \left\{\begin{array}{l} v_k \geq \varphi(v_i)\\ v_k \leq \varphi(v_i^{\rm L}) + \frac{\varphi(v_i^{\rm U})-\varphi(v_i^{\rm L})}{v_i^{\rm U}-v_i^{\rm L}}(v_i-v_i^{\rm L}) \end{array}\right.\\
\text{concave case:} & \quad v_k=\varphi(v_i) \quad \xrightarrow{\displaystyle\text{relax}} \quad \left\{\begin{array}{l} v_k \geq \varphi(v_i^{\rm L}) + \frac{\varphi(v_i^{\rm U})-\varphi(v_i^{\rm L})}{v_i^{\rm U}-v_i^{\rm L}}(v_i-v_i^{\rm L})\\ v_k \leq \varphi(v_i) \end{array}\right.\\
\text{convexo-concave case:} & \quad v_k=\varphi(v_i) \quad \xrightarrow{\displaystyle\text{relax}} \quad \left\{\begin{array}{l}
v_k \geq \left\{\begin{array}{ll}
\varphi(v_i), & \text{if $v_i\leq v_{\rm m}^{\rm cv}$}\\
\varphi(v_i^{\rm U}) + \frac{\varphi(v_i^{\rm U})-\varphi(v_{\rm m}^{\rm cv})}{v_i^{\rm U}-v_{\rm m}^{\rm cv}}(v_i-v_i^{\rm U}), & \text{otherwise}
\end{array}\right.\\
v_k \leq \left\{\begin{array}{ll}
\varphi(v_i), & \text{if $v_i\geq v_{\rm m}^{\rm cc}$}\\
\varphi(v_i^{\rm L}) + \frac{\varphi(v_i^{\rm L})-\varphi(v_{\rm m}^{\rm cc})}{v_i^{\rm L}-v_{\rm m}^{\rm cc}}(v_i-v_i^{\rm L}), & \text{otherwise}
\end{array}\right.
\end{array}\right.\\
& \quad \text{with:}\ v_{\rm m}^{\rm cv},v_{\rm m}^{\rm cc}\in[v_i^{\rm L},v_i^{\rm U}]:\ \varphi'(v_{\rm m}^{\rm cv}) = \textstyle\frac{\varphi(v_i^{\rm U})-\varphi(v_{\rm m}^{\rm cv})}{v_i^{\rm U}-v_{\rm m}^{\rm cv}} \text{ and } \varphi'(v_{\rm m}^{\rm cc}) = \textstyle\frac{\varphi(v_i^{\rm L})-\varphi(v_{\rm m}^{\rm cc})}{v_i^{\rm L}-v_{\rm m}^{\rm cc}}
\f}
\image html exm_uni.png
\n
  .

- <B>Step 3. Polyhedral Outer-Approximation</B>\n
Every convex, nonlinear univariate constraint generated in Step 2 is outer-approximated by constructing supporting cuts at a number of well-chosen points. Although the resulting polyhedral relaxations are inherently weaker than the nonlinear relaxations, LP solvers are currently more robust and faster than NLP solvers.\n
An iterative scheme (a.k.a. sandwich algorithm) can be applied that adds linearization points in such a way that the maximum distance \f$\delta^{\rm max}\f$ between the nonlinear constraint and its polyhedral approximation decreases as the inverse of the square of the number \f$n\f$ of linearization points; that is, \f$\delta^{\rm max}\propto \frac{1}{n^2}\f$. This algorithm proceeds as follows:
  -# Construct cuts at both interval end-points \f$v_i^{\rm L}\f$ and \f$v_i^{\rm U}\f$
  -# <B>REPEAT</B>
    - Identify an interval \f$[v_i^\ell,v_i^{\ell+1}]\f$ with maximum outer-approximation error \f$\delta_{\rm max}\f$
    - Subdivide \f$[v_i^\ell,v_i^{\ell+1}]\f$ at a suitably chosen point \f$v_i^{\rm new}\f$
    .
    <B>UNTIL</B> \f$\delta_{\rm max} < \varepsilon^{\rm tol}\f$
  .
  In particular, several strategies have been proposed for selecting a new linearization point \f$v_i^{\rm new}\f$. The <I>interval bisection rule</I> and the <I>maximum error rule</I> are depicted below.

\image html  OAcvx_strategy.png


\section sec_POLIMG_COMP How to generate a polyhedral relaxation of a factorable function?

For illustration, suppose that we want to compute a polyhedral enclosure for the image set of function:
\f[
{\bf g}(x) := \left(\begin{array}{c}
\log(x_1)+x^2_2 \\
\sin(x_1)-\cos(x_2) \\
\end{array} \right) \qquad \text{with} \qquad 
x \in [1,5]\times [2,6].
\f]

For simplicity, the underlying interval bounds are propagated using the default interval type mc::Interval, the required header files are:
  
\code
#include "interval.hpp"
typedef mc::Interval I;

#include "polimage.hpp"
typedef mc::PolImg<I> PI;
typedef mc::PolVar<I> PV;
\endcode

First, the DAG of the vector function \f${\bf g}\f$ is defined:
\code 
  mc::FFGraph DAG;
  mc::FFVar X[2]; X[0].set( &DAG ); X[1].set( &DAG );
  mc::FFVar F[2]; F[0] = log(X[0])+pow(X[1],2); F[1] = sin(X[0])-cos(X[1]); 
  std::cout << DAG;
\endcode

The following output is displayed, including the auxiliary variables created during the DAG decomposition of the factorable functions:
\verbatim
  DAG VARIABLES:
    X0    => { Z1 Z4 }
    X1    => { Z0 Z3 }

  DAG INTERMEDIATES:
    Z0    <=  SQR( X1 )       => { Z2 }
    Z1    <=  LOG( X0 )       => { Z2 }
    Z2    <=  Z0 + Z1         => { }
    Z3    <=  COS( X1 )       => { Z5 }
    Z4    <=  SIN( X0 )       => { Z5 }
    Z5    <=  Z4 - Z3         => { }
\endverbatim

Next, the polyhedral relaxation environment is created and the main variables are initialized:
\code 
  mc::PolImg<I> Env;
  I IX[2] = { I(1,5), I(2,6) };
  mc::PolVar<I> PX[2]; PX[0].set( &Env, X[0], IX[0] ); PX[1].set( &Env, X[1], IX[1] );
\endcode

Then, the polyhedral relaxation of the image set is initialized by propagating bounds for the lifted variables through the DAG:
\code 
  mc::PolVar<I> PF[2]; DAG.eval( 2, F, PF, 2, X, PX );
  std::cout << Env;
\endcode

This produces the following output, where no cuts have been generated to define the polyhedral enclosure yet:
\verbatim
VARIABLES:
  X0	 in [  1.00000e+00 :  5.00000e+00 ]	 (DAG: X0)
  X1	 in [  2.00000e+00 :  6.00000e+00 ]	 (DAG: X1)
  X2	 in [  4.00000e+00 :  3.60000e+01 ]	 (DAG: Z0)
  X3	 in [  0.00000e+00 :  1.60944e+00 ]	 (DAG: Z1)
  X4	 in [  4.00000e+00 :  3.76094e+01 ]	 (DAG: Z2)
  X6	 in [ -1.00000e+00 :  9.60170e-01 ]	 (DAG: Z3)
  X5	 in [ -1.00000e+00 :  1.00000e+00 ]	 (DAG: Z4)
  X7	 in [ -1.96017e+00 :  2.00000e+00 ]	 (DAG: Z5)

NO AUXILIARY

NO BILINEAR OR FRACTIONAL TERM

NO CUT
\endverbatim

Finally, polyhedral cuts are generated and displayed as follows:
\code 
  Env.generate_cuts( 2, PF );
  std::cout << Env;
\endcode

\verbatim
VARIABLES:
  X0	 in [  1.00000e+00 :  5.00000e+00 ]	 (DAG: X0)
  X1	 in [  2.00000e+00 :  6.00000e+00 ]	 (DAG: X1)
  X2	 in [  4.00000e+00 :  3.60000e+01 ]	 (DAG: Z0)
  X3	 in [  0.00000e+00 :  1.60944e+00 ]	 (DAG: Z1)
  X4	 in [  4.00000e+00 :  3.76094e+01 ]	 (DAG: Z2)
  X6	 in [ -1.00000e+00 :  9.60170e-01 ]	 (DAG: Z3)
  X5	 in [ -1.00000e+00 :  1.00000e+00 ]	 (DAG: Z4)
  X7	 in [ -1.96017e+00 :  2.00000e+00 ]	 (DAG: Z5)

NO AUXILIARY

NO BILINEAR OR FRACTIONAL TERM

CUTS:
  + 1.00000e+00X2 + 1.00000e+00X3 - 1.00000e+00X4 = -0.00000e+00
  + 1.00000e+00X5 - 1.00000e+00X6 - 1.00000e+00X7 = -0.00000e+00
  + 4.00000e+00X3 - 1.60944e+00X0 >= -1.60944e+00
  + 1.00000e+00X0 - 1.00000e+00X3 >= 1.00000e+00
  + 1.00000e+00X0 - 5.00000e+00X3 >= -3.04719e+00
  + 1.00000e+00X0 - 2.75054e+00X3 >= -3.24492e-02
  + 1.00000e+00X0 - 1.80361e+00X3 >= 7.39860e-01
  + 1.00000e+00X0 - 3.81983e+00X3 >= -1.29953e+00
  + 4.00000e+00X2 - 3.20000e+01X1 <= -4.80000e+01
  + 1.00000e+00X5 - 5.40302e-01X0 <= 3.01169e-01
  + 1.00000e+00X5 + 0.00000e+00X0 <= 1.00000e+00
  + 1.00000e+00X5 - 2.73845e-01X0 <= 6.07581e-01
  + 1.00000e+00X5 - 4.08535e-01X0 <= 4.42948e-01
  + 1.00000e+00X5 - 1.37362e-01X0 <= 7.93681e-01
  + 1.00000e+00X5 + 0.00000e+00X0 <= 1.00000e+00
  + 1.00000e+00X5 + 6.31638e-01X0 <= 2.19927e+00
  + 1.00000e+00X5 + 3.22022e-01X0 <= 1.55814e+00
  + 1.00000e+00X5 + 4.79346e-01X0 <= 1.87021e+00
  + 1.00000e+00X5 + 1.61734e-01X0 <= 1.26716e+00
  + 1.00000e+00X6 - 3.44637e-01X1 <= -1.10542e+00
  + 1.00000e+00X6 - 2.79415e-01X1 <= -7.16323e-01
  + 1.00000e+00X2 - 4.00000e+00X1 >= -4.00000e+00
  + 1.00000e+00X2 - 1.20000e+01X1 >= -3.60000e+01
  + 1.00000e+00X2 - 8.00000e+00X1 >= -1.60000e+01
  + 1.00000e+00X2 - 6.00000e+00X1 >= -9.00000e+00
  + 1.00000e+00X2 - 1.00000e+01X1 >= -2.50000e+01
  + 1.00000e+00X5 + 5.35701e-01X0 >= 1.37717e+00
  + 1.00000e+00X5 + 0.00000e+00X0 >= -1.00000e+00
  + 1.00000e+00X5 + 2.71443e-01X0 >= 2.42071e-01
  + 1.00000e+00X5 + 4.04992e-01X0 >= 8.25290e-01
  + 1.00000e+00X5 + 1.36150e-01X0 >= -3.67693e-01
  + 1.00000e+00X5 + 0.00000e+00X0 >= -1.00000e+00
  + 1.00000e+00X5 - 2.83662e-01X0 >= -2.37724e+00
  + 1.00000e+00X5 - 1.42321e-01X0 >= -1.68082e+00
  + 1.00000e+00X5 - 2.13178e-01X0 >= -2.02739e+00
  + 1.00000e+00X5 - 7.12210e-02X0 >= -1.33816e+00
  + 1.00000e+00X6 + 9.09297e-01X1 >= 1.40245e+00
  + 1.00000e+00X6 + 0.00000e+00X1 >= -1.00000e+00
  + 1.00000e+00X6 + 4.78987e-01X1 >= 3.87705e-01
  + 1.00000e+00X6 + 7.05713e-01X1 >= 9.55690e-01
  + 1.00000e+00X6 + 2.41998e-01X1 >= -2.69168e-01
  + 1.00000e+00X6 + 0.00000e+00X1 >= -1.00000e+00
  + 1.00000e+00X6 - 8.07940e-01X1 >= -3.88747e+00
  + 1.00000e+00X6 - 4.18937e-01X1 >= -2.40524e+00
  + 1.00000e+00X6 - 6.19996e-01X1 >= -3.14699e+00
  + 1.00000e+00X6 - 2.11107e-01X1 >= -1.68558e+00
\endverbatim

By default, the cuts in the polyhedral relaxation of a convex univariate terms are generated according to the <I>maximum error rule</I> (see above), and a maximum of 5 cuts are generated for each nonlinear constraint. The cut generation is also controlled by the absolute and relative tolerances on the maximum outer-approximation error, which are both set to \f$10^{-3}\f$ by default. All these default values can be altered as explained below in the section \ref sec_POLIMG_OPT below.

All the cuts can be retreived by using the function mc::PolImg::Cuts, which returns a set of cuts as defined in the class mc::PolCuts.


\section sec_POLIMG_OPT What are the options in computing a polyhedral relaxation?

All the options are defined in the structure mc::PolImg::Options, and the default values can be altered via the public static member mc::PolImg::options; for example:

\code
  Env.options.AGGREG_LQ = false;
  Env.options.SANDWICH_MAXCUT = 7;
\endcode


\section sec_POLIMG_ERR What errors may be encoutered in computing a polyhedral relaxation?

Errors are managed based on the exception handling mechanism of the C++ language. Each time an error is encountered, an instance of the class mc::PolImg::Exceptions is thrown, which contains the type of error. Additional exceptions may be sent by the template argument class in propagating the bounds through the DAG, as well as by the DAG class itself. It is the responsibility of a user to test whether an exception was thrown, and then make the appropriate changes. Should an exception be thrown and not caught by the calling program, execution will abort.

\section sec_POLIMG_REFS References

- Neumaier, A., <A href="https://doi.org/10.1017/S0962492904000194">Complete search in continuous global optimization and constraint satisfaction</A>, <I>Acta Numerica</I>, <B>13</B>(3):271-369, 2004.
- Tawarmalani, M., and N.V. Sahinidis, <A href="http://dx.doi.org/10.1007/s10107-003-0467-6">Global optimization of mixed-integer nonlinear programs: A theoretical and computational study</A>, <I>Mathematical Programming</I>, <B>99</B>(3):563-591, 2004.
- Smith, E.M.B, and C.C. Pantelides, <A href="http://dx.doi.org/10.1016/S0098-1354(98)00286-5">A symbolic reformulation/spatial branch-and-bound algorithm for the global optimisation of nonconvex MINLPs</A>, <I>Computers & Chemical Engineering</I>, <B>23</B>(4-5):457-478, 1999.
.
*/

/*
TODO: 
- Semilinear relaxation of terms that are neither convex nor concave? (e.g., pow(x,3)) ==> OK
- Account for multiple occurence of variables ==> OK
- Enable DC-decomposition of bilinear terms (e.g. lists of product/division terms) ==> OK
- Fix bug in DC scaling ==> OK
- Implement alternative piecewise relaxation approaches ==> OK
- Add all remaining univariate terms (trigo, min/max)
- Distinguish between binary and integer variables?
*/

#ifndef MC__POLIMAGE_H
#define MC__POLIMAGE_H

#include <assert.h>
#include <exception>
#include <fstream>
#include <iomanip>
#include <queue>
#include <map>

#include "ffunc.hpp"
#include "mcop.hpp"

#undef  MC__POLIMG_DEBUG
#undef  MC__POLIMG_DEBUG_CUTS
//#undef  MC__POLIMG_PWMCCORMICK_1D

namespace mc
{

template< class T > class PolImg;
template< class T > class PolLQExpr;
template< class T > class PolCut;
template< class T > class lt_PolVar;

//! @brief C++ class for defining polytopic image variables
////////////////////////////////////////////////////////////////////////
//! mc::PolVar is a C++ class for defining polytopic image variables.
//! The template parameter corresponds to the type used to propagate
//! variable range.
////////////////////////////////////////////////////////////////////////
template< class T >
class PolVar
////////////////////////////////////////////////////////////////////////
{
  friend class PolImg<T>;
  friend class lt_PolVar<T>;

  template< class U > friend  PolVar<U> operator+( const PolVar<U>&, const PolVar<U>& );
  template< class U > friend  PolVar<U> operator+( const PolVar<U>&, const double );
  template< class U > friend  PolVar<U> operator-( const PolVar<U>& );
  template< class U > friend  PolVar<U> operator-( const PolVar<U>&, const PolVar<U>& );
  template< class U > friend  PolVar<U> operator*( const PolVar<U>&, const PolVar<U>& );
  template< class U > friend  PolVar<U> operator*( const PolVar<U>&, const double );			
  template< class U > friend  PolVar<U> operator/( const PolVar<U>&, const PolVar<U>& );
  template< class U > friend  PolVar<U> operator/( const double, const PolVar<U>& );
  template< class U > friend  PolVar<U> operator^( const PolVar<U>&, const PolVar<U>& );
    
  template< class U > friend  PolVar<U> inv ( const PolVar<U>& );
  template< class U > friend  PolVar<U> exp ( const PolVar<U>& );
  template< class U > friend  PolVar<U> log ( const PolVar<U>& );
  template< class U > friend  PolVar<U> xlog( const PolVar<U>& );
  template< class U > friend  PolVar<U> sqrt( const PolVar<U>& );
  template< class U > friend  PolVar<U> sqr ( const PolVar<U>& );
  template< class U > friend  PolVar<U> pow ( const PolVar<U>&, const int );  
  template< class U > friend  PolVar<U> pow ( const PolVar<U>&, const double );  
  template< class U > friend  PolVar<U> cheb( const PolVar<U>&, const unsigned );  
  template< class U > friend  PolVar<U> prod( const unsigned, const PolVar<U>* );  
  template< class U > friend  PolVar<U> cos ( const PolVar<U>& );
  template< class U > friend  PolVar<U> sin ( const PolVar<U>& );
  template< class U > friend  PolVar<U> tan ( const PolVar<U>& );
  template< class U > friend  PolVar<U> acos( const PolVar<U>& );
  template< class U > friend  PolVar<U> asin( const PolVar<U>& );
  template< class U > friend  PolVar<U> atan( const PolVar<U>& );
  template< class U > friend  PolVar<U> cosh( const PolVar<U>& );
  template< class U > friend  PolVar<U> sinh( const PolVar<U>& );
  template< class U > friend  PolVar<U> tanh( const PolVar<U>& );
  template< class U > friend  PolVar<U> fabs( const PolVar<U>& );
  template< class U > friend  PolVar<U> erf( const PolVar<U>& );
  template< class U > friend  PolVar<U> fstep( const PolVar<U>& );
  template< class U > friend  PolVar<U> max( const PolVar<U>&, const PolVar<U>& );  
  template< class U > friend  PolVar<U> min( const PolVar<U>&, const PolVar<U>& );  
  template< class U > friend  PolVar<U> lmtd( const PolVar<U>&, const PolVar<U>& );  
  template< class U > friend  PolVar<U> rlmtd( const PolVar<U>&, const PolVar<U>& );  

  public:
    /** @ingroup POLYTOPE
     *  @{
     */
    //! @brief Enumeration type for variables in factorable program
    enum TYPE{
      VARCONT=0,//!< DAG continuous variable
      VARINT,	//!< DAG integer variable
      AUXCONT,	//!< Auxiliary continuous variable
      AUXINT,	//!< Auxiliary integer variable   
      AUXCST	//!< Auxiliary constant
    };
    //! @brief Typedef for variable identifier in factorable program
    typedef std::pair< TYPE, unsigned > t_idVar;

    //! @brief Return string with variable name for identifier <a>id</a>
    std::string name
      () const
      {
        std::ostringstream ovar;
        switch( _id.first ){
          case VARCONT: ovar << "X"; break;
          case VARINT:  ovar << "Y"; break;
          case AUXCONT: ovar << "W"; break;
          case AUXINT:  ovar << "Z"; break;
          case AUXCST:  ovar << "C"; break;
        }
        ovar << _id.second;
        return ovar.str();
     }
  /** @} */ 
							  
  private:
    //! @brief pointer to underlying polytope
    PolImg<T>* _img;
    //! @brief underlying variable in DAG
    FFVar _var;
    //! @brief variable range
    T _range;
    //! @brief variable identifier (type and index)
    t_idVar _id;
    //! @brief flag indicating whether cuts have already been generated for the associated operation
    mutable bool _hascuts;
    //! @brief variable break-points
    std::set<double> _breakpts;
    //! @brief variable subdivision w.r.t. break-points
    mutable std::pair< std::vector<double>, std::vector<PolVar<T> > > _subdiv;

    //! @brief set as DAG variable <a>X</a> in polytope image <a>P</a> with range <a>B</a>
    void _set
      ( PolImg<T>*P, FFVar const& X, T const& B, const bool cont, const unsigned index )
      { _img = P; _var = X; _range = (!_var.cst()? B: _var.num().val());
        _id = std::make_pair( cont? VARCONT: VARINT, index );
        _breakpts.clear(); reset_subdiv(); reset_cuts(); }
    //! @brief set as auxiliary variable in polytope image <a>P</a> with range <a>B</a>
    void _set
      ( PolImg<T>*P, const T&B, const bool cont, const unsigned index )
      { _img = P; _var = 0; _range = B;
        _id = std::make_pair( cont? AUXCONT: AUXINT, index );
         _breakpts.clear(); reset_subdiv(); reset_cuts(); }
    //! @brief update variable bounds and type
    void _update
      ( const T&B )
      { _range = B; }
    //! @brief update variable bounds and type
    void _update
      ( const bool cont )
      { _id.first = (cont? VARCONT: VARINT); }

    //! @brief push variable subdivision
    void _push_subdiv
      ( double pt ) const
      { if( !_img ) return;
        auto itVar = _img->_Vars.find( const_cast<FFVar*>(&_var) );
        if( itVar == _img->_Vars.end() ) return;
        itVar->second->_subdiv.first.push_back( pt ); }
        //_subdiv.first.push_back( pt ); }
    //! @brief push variable subdivision
    void _push_subdiv
      ( PolVar<T> var ) const
      { if( !_img ) return;
        auto itVar = _img->_Vars.find( const_cast<FFVar*>(&_var) );
        if( itVar == _img->_Vars.end() ) return;
        itVar->second->_subdiv.second.push_back( var ); }
        //_subdiv.second.push_back( var ); }

  public:
    /** @ingroup POLYTOPE
     *  @{
     */
    //! @brief Constructor for a constant value <a>d</a> (default)
    PolVar
      ( const double d=0. )
      : _img(0), _var(FFVar(d)), _range(d), _id( AUXCST, 0 ),
        _hascuts(false) 
      {} 
    //! @brief Constructor for a constant value <a>n</a>
    PolVar
      ( const int n )
      : _img(0), _var(FFVar(n)), _range(n), _id( AUXCST, 0 ),
        _hascuts(false) 
      {}
    //! @brief Constructor for DAG variable <a>X</a> in polytope image <a>P</a> with range <a>B</a>
    PolVar
      ( PolImg<T>*P, FFVar const& X, T const& B=0., const bool cont=true )
      { set( P, X, B, cont ); }
    //! @brief Constructor for auxiliary variable in polytope image <a>P</a> with range <a>B</a>
    PolVar
      ( PolImg<T>*P, T const& B=0., const bool cont=true )
      { set( P, B, cont ); }
    //! @brief Copy constructor for a polytope image <a>P</a>
    PolVar
      ( PolVar<T> const& P )
      : _img(P._img), _var(P._var), _range(P._range), _id( P._id ),
        _hascuts( P._hascuts ), _breakpts( P._breakpts ), _subdiv( P._subdiv )
      {}

    //! @brief Destructor
    virtual ~PolVar
      ()
      {};

    //! @brief Update range <a>B</a> and type <a>cont</a> of variable in polytope image
    PolVar<T>& update
      ( T const& B )
      { if( !_img ) return *this;
        PolVar<T>* polvar = _img->_update_var( &_var, B );
        if( polvar != nullptr ) *this = *polvar;
        reset_subdiv(); reset_cuts(); return *this; }
    //! @brief Update range <a>B</a> and type <a>cont</a> of variable in polytope image
    PolVar<T>& update
      ( bool const cont )
      { if( !_img ) return *this;
        PolVar<T>* polvar = _img->_update_var( &_var, cont );
        if( polvar != nullptr ) *this = *polvar;
        reset_subdiv(); reset_cuts(); return *this; }
    //! @brief set as DAG variable <a>X</a> in polytope image <a>P</a> with range <a>B</a>
    PolVar<T>& set
      ( PolImg<T>*P, FFVar const& X, T const& B=0., const bool cont=true )
      { *this = *P->_append_var( &X, B, cont );
        _breakpts.clear(); reset_subdiv(); reset_cuts(); return *this; }
    //! @brief set as auxiliary variable in polytope image <a>P</a> with range <a>B</a>
    PolVar<T>& set
      ( PolImg<T>*P, T const& B=0., const bool cont=true )
      { *this = *P->_append_aux( B, cont );
        _breakpts.clear(); reset_subdiv(); reset_cuts(); return *this; }

    //! @brief get variable integrality
    bool discr
      () const
      { return (_id.first == VARINT || _id.first == AUXINT? true: false); }
    //! @brief get variable constness
    bool cst
      () const
      { return (_var.id().first == FFVar::CINT
             || _var.id().first == FFVar::CREAL? true: false); }
    //! @brief get variable range
    T range
      () const
      { return _range; }//_var.cst()? _var.num().val(): _range; }
      //{ return _range; }
    //! @brief get pointer to polytopic image
    PolImg<T>* image
      () const
      { return _img; }
    //! @brief get pointer to variable identifier
    t_idVar id
      () const
      { return _id; }
    //! @brief get reference to DAG variable
    FFVar& var
      ()
      { return _var; }
    //! @brief get const reference to DAG variable
    const FFVar& var
      () const
      { return _var; }

    //! @brief add variable break-points
    void add_breakpt
      ( const double bkpt );
    //! @brief get variable break-points
    const std::set<double>& breakpts
      () const
      { return _breakpts; }

    //! @brief get variable subdivision
    const std::pair< std::vector<double>, std::vector< PolVar<T> > >& subdiv
      () const
      { return _subdiv; }
    //! @brief reset variable subdivision
    void reset_subdiv
      () const
      { _subdiv.first.clear(); _subdiv.second.clear(); }
    //! @brief create variable subdivision
    const std::vector<double>& create_subdiv
      ( const double XL, const double XU, const bool reset=false ) const;
    //! @brief set SOS2 variable subdivision
    const std::vector< PolVar<T> >& SOS2_subdiv
      ( FFOp*pOp=0, const bool reset=false ) const;
    //! @brief set linear binary variable subdivision
    const std::vector< PolVar<T> >& BIN_subdiv
      ( FFOp*pOp=0, const bool reset=false ) const;

    //! @brief get cuts flag 
    bool has_cuts
      () const
      { return _hascuts; }
    //! @brief reset cuts flag to false
    void reset_cuts
      () const
      { _hascuts = false; }
    //! @brief set cuts flag to true
    void set_cuts
      () const
      { _hascuts = true; }
    //! @brief propagate cuts backwards through polyhedral image
    void generate_cuts
      () const;
    //! @brief propagate linear-quadratic expressions backwards through polyhedral image
    bool generate_cuts_aggreg
      ( PolLQExpr<T>*&pLQ ) const;
    //! @brief propagate default cuts backwards through polyhedral image
    void generate_cuts_default
      () const;

    //! @brief Public overloads
    PolVar<T>& operator= ( const PolVar<T>& );
    PolVar<T>& operator= ( const double );
    PolVar<T>& operator= ( const int );

    PolVar<T>& operator+=( const PolVar<T>& );
    PolVar<T>& operator-=( const PolVar<T>& );
    PolVar<T>& operator*=( const PolVar<T>& );
    PolVar<T>& operator/=( const PolVar<T>& );
  /** @} */ 
};

//! @brief C++ structure for ordering of polytopic variables
template <class T> 
struct lt_PolVar
{
  bool operator()
    ( const PolVar<T>*Var1, const PolVar<T>*Var2 ) const
    {
      // Order variables/constants w.r.t. their types first
      if( Var1->_id.first < Var2->_id.first ) return true;
      if( Var1->_id.first > Var2->_id.first ) return false;
      // If variables, order w.r.t. their index next
      switch( Var1->_id.first ){
        case PolVar<T>::VARCONT: case PolVar<T>::VARINT:
        case PolVar<T>::AUXCONT: case PolVar<T>::AUXINT:
          if( Var1->_id.second < Var2->_id.second ) return true;
          if( Var1->_id.second > Var2->_id.second ) return false;
          break;
        case PolVar<T>::AUXCST:
          lt_FFNum ltNum;
          return ltNum( &Var1->_var.num(), &Var2->_var.num() );
          break;
      }
      return false;
    }
};

template <class T> 
inline
std::ostream&
operator <<
( std::ostream&out, const PolVar<T>& P )
{
  return out << P.name();
}

template <class T> 
inline 
PolVar<T>&
PolVar<T>::operator=
( const PolVar<T>& P ) 
{
  _img       = P._img;
  _var       = P._var;
  _range     = P._range;
  _id        = P._id;
  _hascuts   = P._hascuts;
  _breakpts  = P._breakpts;
  _subdiv    = P._subdiv;
  return *this; 
}

template <class T> 
inline 
PolVar<T>&
PolVar<T>::operator=
( const double d )
{
  _img       = 0;
  _var       = FFVar(d);
  _range     = d;
  _id        = std::make_pair( AUXCST, 0 );
  _hascuts   = false;
  _breakpts.clear();
  reset_subdiv();
  return *this; 
}

template <class T> 
inline 
PolVar<T>&
PolVar<T>::operator=
( const int n ) 
{
  _img       = 0;
  _var       = FFVar(n);
  _range     = n;
  _id        = std::make_pair( AUXCST, 0 );
  _hascuts   = false;
  _breakpts.clear();
  reset_subdiv();
  return *this; 
}

template <typename T> inline PolVar<T>&
PolVar<T>::operator +=
( const PolVar<T>&P1 )
{
   PolVar<T> P2( *this );
   *this = P2 + P1;
   return *this;
}

template <typename T> inline PolVar<T>&
PolVar<T>::operator -=
( const PolVar<T>&P1 )
{
   PolVar<T> P2( *this );
   *this = P2 - P1;
   return *this;
}

template <typename T> inline PolVar<T>&
PolVar<T>::operator *=
( const PolVar<T>&P1 )
{
   PolVar<T> P2( *this );
   *this = P2 * P1;
   return *this;
}

template <typename T> inline PolVar<T>&
PolVar<T>::operator /=
( const PolVar<T>&P1 )
{
   PolVar<T> P2( *this );
   *this = P2 / P1;
   return *this;
}

template <typename T>
inline void
PolVar<T>::add_breakpt
( const double bkpt )
{
   if( !_img ) return;
   auto itVar = _img->_Vars.find( &_var );
   if( itVar == _img->_Vars.end() ) return;
   double atol = _img->options.BREAKPOINT_ATOL,
          rtol = _img->options.BREAKPOINT_RTOL;
   // Reject breakpoint if at lower/upper bound of current interval
   if( isequal( Op<T>::l(itVar->second->range()), bkpt, atol, rtol ) 
    || isequal( Op<T>::u(itVar->second->range()), bkpt, atol, rtol ) ) return;
   // Reject breakpoint if at lower/upper bound of current interval
   auto itL = _breakpts.lower_bound( bkpt );
   if( itL!=_breakpts.end() && isequal( *itL, bkpt, atol, rtol ) ) return;
   if( itL!=_breakpts.begin() && isequal( *(--itL), bkpt, atol, rtol ) ) return;
   //auto itU = _breakpts.upper_bound( bkpt );
   //if( itU!=_breakpts.end() && isequal( *itU, bkpt, atol, rtol ) ) return;
   itVar->second->_breakpts.insert( bkpt );
   itVar->second->reset_subdiv();
   _breakpts.insert( bkpt );
   reset_subdiv();
}

template <typename T>
inline const std::vector<double>&
PolVar<T>::create_subdiv
( const double XL, const double XU, const bool reset ) const
{
  if( !reset && !_subdiv.first.empty() ) return _subdiv.first;
  if(  reset && !_subdiv.first.empty() ) reset_subdiv();
  _push_subdiv( XL ); 
  if( !_breakpts.empty() ){
    for( auto it = _breakpts.upper_bound( XL );
      it != _breakpts.end() && *it < XU; ++it )
      _push_subdiv( *it ); 
  }
  _push_subdiv( XU ); 
  return _subdiv.first;
}

template <typename T>
inline const std::vector< PolVar<T> >&
PolVar<T>::BIN_subdiv
( FFOp*pOp, const bool reset ) const
{
  if( !reset && !_subdiv.second.empty() ) return _subdiv.second;
  if(  reset && !_subdiv.second.empty() ) _subdiv.second.clear();
  if( !_img ) return _subdiv.second;

  const unsigned NINTS = _subdiv.first.size()-1;
  double coef[NINTS];
  for( unsigned isub=0; isub<NINTS; isub++ ){
    coef[isub] = _subdiv.first[isub] - _subdiv.first[isub+1];
    _push_subdiv( PolVar<T>( _img, Op<T>::zeroone(), true ) ); 
  }
  _img->_add_cut( pOp, PolCut<T>::EQ, _subdiv.first[0], NINTS, _subdiv.second.data(), coef, *this, 1. );

  for( unsigned isub=0; isub<NINTS-1; isub++ ){
    _push_subdiv( PolVar<T>( _img, Op<T>::zeroone(), false ) );
    _img->_add_cut( pOp, PolCut<T>::LE, 0., _subdiv.second[NINTS+isub], 1., _subdiv.second[isub], -1. );
    _img->_add_cut( pOp, PolCut<T>::GE, 0., _subdiv.second[NINTS+isub], 1., _subdiv.second[isub+1], -1. );
  }

  return _subdiv.second;
}

template <typename T>
inline const std::vector< PolVar<T> >&
PolVar<T>::SOS2_subdiv
( FFOp*pOp, const bool reset ) const
{
  if( !reset && !_subdiv.second.empty() ) return _subdiv.second;
  if(  reset && !_subdiv.second.empty() ) _subdiv.second.clear();
  if( !_img ) return _subdiv.second;

  const unsigned NINTS = _subdiv.first.size();
  double coef[NINTS];
  for( unsigned isub=0; isub<NINTS; isub++ ){
    coef[isub] = 1.;
    _push_subdiv( PolVar<T>( _img, Op<T>::zeroone(), true ) ); 
  }
  _img->_add_cut( pOp, PolCut<T>::EQ, 1., NINTS, _subdiv.second.data(), coef );
  _img->_add_cut( pOp, PolCut<T>::EQ, 0., NINTS, _subdiv.second.data(), _subdiv.first.data(), *this, -1. );
  _img->_add_cut( pOp, PolCut<T>::SOS2, 1., NINTS, _subdiv.second.data(), coef );

  return _subdiv.second;
}

template <typename T> inline void
PolVar<T>::generate_cuts
() const
{
  if( has_cuts() ) return;
  set_cuts();
#ifdef MC__POLIMG_DEBUG_CUTS
  std::cout << "CUTS FOR " << _var << ": " << *_var.ops().first << std::endl;
#endif
  
  // Nothing to do if underlying dag variable is constant
  if( !_var.ops().first ){
    if( !_var.cst() )
      throw typename PolImg<T>::Exceptions( PolImg<T>::Exceptions::INTERN );
    //assert( _var->cst() );
    return;
  }    

  PolLQExpr<T>*pLQ = 0;
  if( _img->options.AGGREG_LQ && generate_cuts_aggreg( pLQ ) ){
    auto itCut = _img->_add_cut( _var.ops().first, PolCut<T>::EQ, -pLQ->cst(), *this, -1. );
    (*itCut)->append( pLQ->size(),  pLQ->var(),   pLQ->coef() ); 
    (*itCut)->append( pLQ->qsize(), pLQ->qvar1(), pLQ->qvar2(), pLQ->qcoef() );

    const PolVar<T>*pVar = pLQ->var();
    for( unsigned i=0; i<pLQ->size(); i++ ){
      auto itVar = _img->_Vars.find( const_cast<FFVar*>(&pVar[i]._var) );
      if( itVar != _img->_Vars.end() ) itVar->second->generate_cuts();
    }
    
    const PolVar<T>*pqVar1 = pLQ->qvar1();
    const PolVar<T>*pqVar2 = pLQ->qvar2();
    for( unsigned i=0; i<pLQ->qsize(); i++ ){
      auto itqVar1 = _img->_Vars.find( const_cast<FFVar*>(&pqVar1[i]._var) );
      if( itqVar1 != _img->_Vars.end() ) itqVar1->second->generate_cuts();
      auto itqVar2 = _img->_Vars.find( const_cast<FFVar*>(&pqVar2[i]._var) );
      if( itqVar2 != _img->_Vars.end() ) itqVar2->second->generate_cuts();
    }
  }

  else{
    generate_cuts_default();
    for( auto && op : _var.ops().first->pops ){ 
      auto itVar = _img->_Vars.find( op );
      if( itVar != _img->_Vars.end() ) itVar->second->generate_cuts();
    }
  }
}

template <typename T> inline void
PolVar<T>::generate_cuts_default
() const
{
  // CNST, VAR,
  // PLUS, SHIFT, NEG, MINUS, TIMES, SCALE, DIV, INV,
  // PROD, IPOW, DPOW, CHEB, SQR, SQRT, EXP, LOG, XLOG,
  // SIN, COS, TAN, ASIN, ACOS, ATAN, COSH, SINH, TANH,
  // FABS, ERF, FSTEP, MINF, MAXF, INTER, LMTD, RLMTD
  
  switch( _var.ops().first->type ){
   case FFOp::CNST:
   case FFOp::VAR:
    break;

   case FFOp::SHIFT:
   case FFOp::PLUS:
    _img->_add_cuts_PLUS( this, _var.ops().first->pops[0], _var.ops().first->pops[1] );
    break;

   case FFOp::MINUS:
    _img->_add_cuts_MINUS( this, _var.ops().first->pops[0], _var.ops().first->pops[1] );
    break;

   case FFOp::NEG:
    _img->_add_cuts_NEG( this, _var.ops().first->pops[0] );
    break;

   case FFOp::SCALE:
   case FFOp::TIMES:
    _img->_add_cuts_TIMES( this, _var.ops().first->pops[0], _var.ops().first->pops[1] );
    break;

   case FFOp::INV:
   case FFOp::DIV:
    _img->_add_cuts_DIV( this, _var.ops().first->pops[0], _var.ops().first->pops[1] );
    break;

   case FFOp::IPOW:
    _img->_add_cuts_IPOW( this, _var.ops().first->pops[0], _var.ops().first->pops[1]->num().n );
    break;

   case FFOp::DPOW:
    _img->_add_cuts_DPOW( this, _var.ops().first->pops[0], _var.ops().first->pops[1]->num().x );
    break;

   case FFOp::CHEB:
    _img->_add_cuts_CHEB( this, _var.ops().first->pops[0], _var.ops().first->pops[1]->num().n );
    break;

   case FFOp::PROD:
    _img->_add_cuts_PROD( this, _var.ops().first->pops );
    break;

   case FFOp::SQR:
    _img->_add_cuts_SQR( this, _var.ops().first->pops[0] );
    break;

   case FFOp::SQRT:
    _img->_add_cuts_SQRT( this, _var.ops().first->pops[0] );
    break;

   case FFOp::EXP:
    _img->_add_cuts_EXP( this, _var.ops().first->pops[0] );
    break;

   case FFOp::LOG:
    _img->_add_cuts_LOG( this, _var.ops().first->pops[0] );
    break;

   case FFOp::XLOG:
    _img->_add_cuts_XLOG( this, _var.ops().first->pops[0] );
    break;

   case FFOp::FABS:
    _img->_add_cuts_FABS( this, _var.ops().first->pops[0] );
    break;

   case FFOp::ERF:
    _img->_add_cuts_ERF( this, _var.ops().first->pops[0] );
    break;

   case FFOp::FSTEP:
    _img->_add_cuts_FSTEP( this, _var.ops().first->pops[0] );
    break;

   case FFOp::COS:
    _img->_add_cuts_COS( this, _var.ops().first->pops[0] );
    break;

   case FFOp::SIN:
    _img->_add_cuts_SIN( this, _var.ops().first->pops[0] );
    break;

   case FFOp::TAN:
    _img->_add_cuts_TAN( this, _var.ops().first->pops[0] );
    break;

   case FFOp::ACOS:
    _img->_add_cuts_ACOS( this, _var.ops().first->pops[0] );
    break;

   case FFOp::ASIN:
    _img->_add_cuts_ASIN( this, _var.ops().first->pops[0] );
    break;

   case FFOp::ATAN:
    _img->_add_cuts_ATAN( this, _var.ops().first->pops[0] );
    break;

   case FFOp::COSH:
    _img->_add_cuts_COSH( this, _var.ops().first->pops[0] );
    break;

   case FFOp::SINH:
    _img->_add_cuts_SINH( this, _var.ops().first->pops[0] );
    break;

   case FFOp::TANH:
    _img->_add_cuts_TANH( this, _var.ops().first->pops[0] );
    break;

   case FFOp::INTER:
    _img->_add_cuts_INTER( this, _var.ops().first->pops[0], _var.ops().first->pops[1] );
    break;

   case FFOp::MINF:
    _img->_add_cuts_MINF( this, _var.ops().first->pops[0], _var.ops().first->pops[1] );
    break;

   case FFOp::MAXF:
    _img->_add_cuts_MAXF( this, _var.ops().first->pops[0], _var.ops().first->pops[1] );
    break;

//   case FFOp::LMTD:
//    _img->_add_cuts_LMTD( this, _var.ops().first->pops[0], _var.ops().first->pops[1] );
//    break;

//   case FFOp::RLMTD:
//    _img->_add_cuts_RLMTD( this, _var.ops().first->pops[0], _var.ops().first->pops[1] );
//    break;

   default:
    std::cout << "  -> EXCEPTION FOR " << _var << ": " << *_var.ops().first << std::endl;
    throw typename PolImg<T>::Exceptions( PolImg<T>::Exceptions::INTERN );
  }

  for( auto itOp=_var.ops().first->pops.begin(); itOp!=_var.ops().first->pops.end(); ++itOp ){ 
    auto itVar = _img->_Vars.find( *itOp );
    if( itVar != _img->_Vars.end() ) itVar->second->generate_cuts();
  }

  return;
}

template <typename T> inline bool
PolVar<T>::generate_cuts_aggreg
( PolLQExpr<T>*&pLQ ) const
{
#ifdef MC__POLIMG_DEBUG_CUTS
  std::cout << "  -> LINEAR CUTS FOR " << _var << ": " << *_var.ops().first << std::endl;
#endif

  switch( _var.ops().first->type ){
   case FFOp::SHIFT:
   case FFOp::PLUS:
    if( !_img->_add_LQ_PLUS( pLQ, this, _var.ops().first->pops[0], _var.ops().first->pops[1] ) )
      return false;
    break;

   case FFOp::MINUS:
    if( !_img->_add_LQ_MINUS( pLQ, this, _var.ops().first->pops[0], _var.ops().first->pops[1] ) )
      return false;
    break;

   case FFOp::NEG:
    if( !_img->_add_LQ_NEG( pLQ, this, _var.ops().first->pops[0] ) ) 
      return false;
    break;

   case FFOp::SCALE:
   case FFOp::TIMES:
    if( !_img->_add_LQ_TIMES( pLQ, this, _var.ops().first->pops[0], _var.ops().first->pops[1] ) )
      return false;
    break;

   case FFOp::INV:
   case FFOp::DIV:
    if( !_img->_add_LQ_DIV( pLQ, this, _var.ops().first->pops[0], _var.ops().first->pops[1] ) )
      return false;
    break;

   case FFOp::SQR:
    if( !_img->_add_LQ_SQR( pLQ, this, _var.ops().first->pops[0] ) )
      return false;
    break;

   case FFOp::FABS:
    if( !_img->_add_LQ_FABS( pLQ, this, _var.ops().first->pops[0] ) )
      return false;
    break;

   case FFOp::FSTEP:
    if( !_img->_add_LQ_FSTEP( pLQ, this, _var.ops().first->pops[0] ) )
      return false;
    break;

   case FFOp::MINF:
    if( !_img->_add_LQ_MINF( pLQ, this, _var.ops().first->pops[0], _var.ops().first->pops[1] ) )
      return false;
    break;

   case FFOp::MAXF:
    if( !_img->_add_LQ_MAXF( pLQ, this, _var.ops().first->pops[0], _var.ops().first->pops[1] ) )
      return false;
    break;

   default:
    return false;

   // INTER?
  }

  for( auto&& op : _var.ops().first->pops ){ 
    auto itVar = _img->_Vars.find( op );
    if( itVar != _img->_Vars.end() ) itVar->second->generate_cuts_aggreg( pLQ );
  }
  return true;
}

//! @brief C++ structure for ordering of monomial pairs
template< class T >
struct lt_pqPolVar
{
  // Comparison operator
  bool operator
    ()
    ( std::pair<const PolVar<T>*, const PolVar<T>*> const& pVar1,
      std::pair<const PolVar<T>*, const PolVar<T>*> const& pVar2 )
    const
    {
      // Order based on first monomial first
      if( lt_PolVar<T>()( pVar1.first, pVar2.first ) ) return true;
      if( lt_PolVar<T>()( pVar2.first, pVar1.first ) ) return false;
      // Order based on second monomial next
      if( lt_PolVar<T>()( pVar1.second, pVar2.second ) ) return true;
      if( lt_PolVar<T>()( pVar2.second, pVar1.second ) ) return false;
      // Pairs are identical on reaching this point
      return false;
    }
};

//! @brief C++ structure for holding information about linear-quadratic subexpressions in a polytopic image
////////////////////////////////////////////////////////////////////////
//! mc::PolLQExpr is a C++ structure for holding information about
//! linear or quadratic subexpressions in a polytopic image.
////////////////////////////////////////////////////////////////////////
template< class T >
class PolLQExpr
////////////////////////////////////////////////////////////////////////
{
 private:
  //! @brief Map of variables and coefficients in linear terms 
  std::map< const PolVar<T>*, double, lt_PolVar<T> > _terms;
  //! @brief Vector of coefficients in linear terms
  std::vector< double > _coef;
  //! @brief Vector of variables in linear terms 
  std::vector< PolVar<T> > _var;
  //! @brief Map of variables and coefficients in linear terms 
  std::map< std::pair<const PolVar<T>*, const PolVar<T>*>, double, lt_pqPolVar<T> > _qterms;
  //! @brief Vector of coefficients in quadratic terms
  std::vector< double > _qcoef;
  //! @brief Vector of variables in quadratic terms 
  std::vector< PolVar<T> > _qvar1;
  //! @brief Vector of variables in quadratic terms 
  std::vector< PolVar<T> > _qvar2;
  //! @brief Constant term
  double _const;

 public:
  //! @brief Constructor for linear-quadratic expression
  PolLQExpr<T>
    ( const PolVar<T>*Var )
    : _const(0.)
    { _terms.insert( std::make_pair( Var, 1. ) ); }
  //! @brief Add constant term <a>cst</a>
  void add
    ( const double cst )
    { _const += cst; }
  //! @brief Add linear term <a>coef * Var</a>
  void add
    ( const double coef, const PolVar<T>*Var )
    { auto itVar = _terms.find( Var );
      if( itVar != _terms.end() )
        itVar->second += coef;
      else
        _terms.insert( std::make_pair( Var, coef ) ); }
  //! @brief Add quadratic term <a>coef * Var</a>
  void add
    ( const double coef, const PolVar<T>*Var1, const PolVar<T>*Var2 )
    { auto itVar = _qterms.find( std::make_pair( Var1, Var2 ) );
      if( itVar != _qterms.end() )
        itVar->second += coef;
      else
        _qterms.insert( std::make_pair( std::make_pair( Var1, Var2 ), coef ) ); }
  //! @brief Substitute variable <a>Res</a> with term <a>coef * Var1 * Var2</a>
  bool substitute
    ( const PolVar<T>*Res, const double coef, const PolVar<T>*Var1,
      const PolVar<T>*Var2 )
    { auto itRes = _terms.find( Res );
      if( itRes == _terms.end() )
        return false;
      add( itRes->second * coef, Var1, Var2 );
      _terms.erase( itRes );
      return true; }
  //! @brief Substitute variable <a>Res</a> with term <a>coef * Var</a>
  bool substitute
    ( const PolVar<T>*Res, const double coef, const PolVar<T>*Var )
    { auto itRes = _terms.find( Res );
      if( itRes == _terms.end() )
        return false;
      add( itRes->second * coef, Var );
      _terms.erase( itRes );
      return true; }
  //! @brief Substitute variable <a>Res</a> with term <a>cst</a>
  bool substitute
    ( const PolVar<T>*Res, const double cst )
    { auto itRes = _terms.find( Res );
      if( itRes == _terms.end() )
        return false;
      add( itRes->second * cst );
      _terms.erase( itRes );
      return true; }
  //! @brief Substitute variable <a>Res</a> with term <a>coef * Var + cst</a>
  bool substitute
    ( const PolVar<T>*Res, const double coef, const PolVar<T>*Var, const double cst )
    { auto itRes = _terms.find( Res );
      if( itRes == _terms.end() )
        return false;
      add( itRes->second * coef, Var );
      add( itRes->second * cst );
      _terms.erase( itRes );
      return true; }
  //! @brief Substitute variable <a>Res</a> with term <a>coef1 * Var1 + coef2 * Var2</a>
  bool substitute
    ( const PolVar<T>*Res, const double coef1, const PolVar<T>*Var1,
      const double coef2, const PolVar<T>*Var2 )
    { auto itRes = _terms.find( Res );
      if( itRes == _terms.end() )
        return false;
      add( itRes->second * coef1, Var1 );
      add( itRes->second * coef2, Var2 );
      _terms.erase( itRes );
      return true; }
  //! @brief Return constant term
  const double cst
    () const
    { return _const; }
  //! @brief Return coefficient array for linear terms
  const double* coef
    ()
    { _coef.resize( _terms.size() );
      auto it = _terms.cbegin();
      for( unsigned i=0; it != _terms.cend(); ++it, i++ )
        _coef[i] = it->second;
      return _coef.data(); }
  //! @brief Return variable array for linear terms
  const PolVar<T>* var
    ()
    { _var.resize( _terms.size() );
      auto it = _terms.cbegin();
      for( unsigned i=0; it != _terms.cend(); ++it, i++ )
        _var[i] = *it->first;
      return _var.data(); }
  //! @brief Return number of linear terms
  const unsigned size
    () const
    { return _terms.size(); }
  //! @brief Return coefficient array for quadratic terms
  const double* qcoef
    ()
    { _qcoef.resize( _qterms.size() );
      auto it = _qterms.cbegin();
      for( unsigned i=0; it != _qterms.cend(); ++it, i++ )
        _qcoef[i] = it->second;
      return _qcoef.data(); }
  //! @brief Return left variable array for quadratic terms
  const PolVar<T>* qvar1
    ()
    { _qvar1.resize( _qterms.size() );
      auto it = _qterms.cbegin();
      for( unsigned i=0; it != _qterms.cend(); ++it, i++ )
        _qvar1[i] = *it->first.first;
      return _qvar1.data(); }
  //! @brief Return right variable array for quadratic terms
  const PolVar<T>* qvar2
    ()
    { _qvar2.resize( _qterms.size() );
      auto it = _qterms.cbegin();
      for( unsigned i=0; it != _qterms.cend(); ++it, i++ )
        _qvar2[i] = *it->first.second;
      return _qvar2.data(); }
  //! @brief Return number of quadratic terms
  const unsigned qsize
    () const
    { return _qterms.size(); }
};

//! @brief C++ template class for defining cuts in a polytopic image
////////////////////////////////////////////////////////////////////////
//! mc::PolCut is a C++ template class defining cuts in a polytopic
//! image.
////////////////////////////////////////////////////////////////////////
template< class T >
class PolCut
////////////////////////////////////////////////////////////////////////
{
  // friends of class PolCut for operator overloading
  template< class U > friend std::ostream& operator<<( std::ostream&, const PolCut<U>& );

public:
  /** @ingroup LPRELAX
   *  @{
   */
  //! @brief Enumeration type for cuts and special sets
  enum TYPE{
    EQ=0,	//!< Equality constraint Ax=b
    LE,		//!< Inequality constraint Ax<=b
    GE,		//!< Inequality constraint Ax>=b
    SOS1,	//!< SOS1-type constraint 
    SOS2,	//!< SOS2-type constraint 
    NLIN    //!< Nonlinear constraint y=f(x)
  };
  /** @} */

private:
  //! @brief Pointer to defining operation
  FFOp* _op;
  //! @brief Type of cut
  TYPE _type;
  //! @brief Right-hand side
  double _rhs;
  //! @brief Participating variables in linear terms
  std::vector<PolVar<T>> _var;
  //! @brief Coefficients in linear terms
  std::vector<double> _coef;
  //! @brief Participating variables in quadratic terms
  std::vector<PolVar<T>> _qvar1;
  //! @brief Participating variables in quadratic terms
  std::vector<PolVar<T>> _qvar2;
  //! @brief Coefficients
  std::vector<double> _qcoef;

public:
  /** @ingroup LPRELAX
   *  @{
   */
  //! @brief Retreive type of cut
  TYPE type
    () const
    { return _type; }
  //! @brief Retreive right-hand side
  double rhs
    () const
    { return _rhs; }
  //! @brief Retreive reference to right-hand side
  double& rhs
    ()
    { return _rhs; }
  //! @brief Retreive number of linear terms
 unsigned nvar
    () const
    { return _var.size(); }
  //! @brief Retreive coefficients of linear terms
  const double* coef
    () const 
    { return _coef.data(); }
  //! @brief Retreive variables of linear terms
  const PolVar<T>* var
    () const 
    { return _var.data(); }
  //! @brief Retreive number of quadratic terms
 unsigned nqvar
    () const
    { return _qvar1.size(); }
  //! @brief Retreive coefficients of quadratic terms
  const double* qcoef
    () const 
    { return _qcoef.data(); }
  //! @brief Retreive variables of quadratic terms
  const PolVar<T>* qvar1
    () const 
    { return _qvar1.data(); }
  //! @brief Retreive variables of quadratic terms
  const PolVar<T>* qvar2
    () const 
    { return _qvar2.data(); }
  //! @brief Retreive pointer to corresponding operation in DAG
  FFOp*& op
    ()
    { return _op; }
  //! @brief Retreive pointer to corresponding operation in DAG
  const FFOp* op
    () const
    { return _op; }

  //! @brief Constructor for nonlinear cut w/ 1 variable
  PolCut
    ( FFOp*op, const PolVar<T>&X1, const PolVar<T>&X2 )
    : _op(op), _type(NLIN), _rhs(0.), _var(2)
    {
      _var[0] = X1;
      _var[1] = X2;
    }
  //! @brief Constructor for nonlinear cut w/ 2 variables
  PolCut
    ( FFOp*op, const PolVar<T>&X1, const PolVar<T>&X2, const PolVar<T>&X3 )
    : _op(op), _type(NLIN), _rhs(0.), _var(3)
    {
      _var[0] = X1;
      _var[1] = X2;
      _var[2] = X3;
    }
  //! @brief Constructor for nonlinear cut w/ 1 variable + 1 constant
  PolCut
    ( FFOp*op, const PolVar<T>&X1, const PolVar<T>&X2, const double X3 )
    : _op(op), _type(NLIN), _rhs(X3), _var(2)
    {
      _var[0] = X1;
      _var[1] = X2;
    }
  //! @brief Constructor for cut w/o participating variable
  PolCut
    ( FFOp*op, TYPE type, const double b )
    : _op(op), _type(type), _rhs(b)
    {}
  //! @brief Constructor for cut w/ 1 linear variable
  PolCut
    ( FFOp*op, TYPE type, const double b, const PolVar<T>&X1, const double a1 )
    : _op(op), _type(type), _rhs(b), _var(1), _coef(1)
    {
      _coef[0] = a1;
      _var[0]  = X1;
    }
  //! @brief Constructor for cut w/ 2 linear variables
  PolCut
    ( FFOp*op, TYPE type, const double b, const PolVar<T>&X1, const double a1,
      const PolVar<T>&X2, const double a2 )
    : _op(op), _type(type), _rhs(b), _var(2), _coef(2)
    {
      _coef[0] = a1;
      _coef[1] = a2;
      _var[0]  = X1;
      _var[1]  = X2;
    }
  //! @brief Constructor for cut w/ 3 linear variables
  PolCut
    ( FFOp*op, TYPE type, const double b, const PolVar<T>&X1, const double a1,
      const PolVar<T>&X2, const double a2, const PolVar<T>&X3, const double a3 )
    : _op(op), _type(type), _rhs(b), _var(3), _coef(3)
    {
      _coef[0] = a1;
      _coef[1] = a2;
      _coef[2] = a3;
      _var[0]  = X1;
      _var[1]  = X2;
      _var[2]  = X3;
    }
  //! @brief Constructor for cut w/ 4 linear variables
  PolCut
    ( FFOp*op, TYPE type, const double b, const PolVar<T>&X1, const double a1,
      const PolVar<T>&X2, const double a2, const PolVar<T>&X3, const double a3,
      const PolVar<T>&X4, const double a4 )
    : _op(op), _type(type), _rhs(b), _var(4), _coef(4)
    {
      _coef[0] = a1;
      _coef[1] = a2;
      _coef[2] = a3;
      _coef[3] = a4;
      _var[0]  = X1;
      _var[1]  = X2;
      _var[2]  = X3;
      _var[3]  = X4;
    }
  //! @brief Constructor for cut w/ <a>n</a> linear variables
  PolCut
    ( FFOp*op, TYPE type, const double b, const unsigned n, const PolVar<T>*X,
      const double*a )
    : _op(op), _type(type), _rhs(b), _var(n), _coef(n)
    {
      for( unsigned ivar=0; ivar<n; ivar++ ){
        _coef[ivar] = a[ivar];
        _var[ivar]  = X[ivar];
      }
    }
  //! @brief Constructor for cut w/ selection among <a>n</a> linear variables
  PolCut
    ( FFOp*op, TYPE type, const double b, const std::set<unsigned>&ndx,
      const PolVar<T>*X, const double*a )
    : _op(op), _type(type), _rhs(b), _var(ndx.size()), _coef(ndx.size())
    {
      std::set<unsigned>::const_iterator it = ndx.begin();
      for( unsigned ivar=0; it!=ndx.end(); ++it, ivar++ ){
        _coef[ivar] = a[*it];
        _var[ivar]  = X[*it];
      }
    }
  //! @brief Constructor for cut w/ variable and coefficient maps
  template <typename U> PolCut
    ( FFOp*op, TYPE type, const double b, const std::map<U,PolVar<T>>&X,
      const std::map<U,double>&a )
    : _op(op), _type(type), _rhs(b), _var(a.size()), _coef(a.size())
    {
      auto ita = a.cbegin();
      for( unsigned ivar=0; ita!=a.cend(); ++ita, ivar++ ){
        _coef[ivar] = ita->second;
        auto itX = X.find(ita->first);
        if( itX == X.cend() )
          throw typename PolImg<T>::Exceptions( PolImg<T>::Exceptions::BADCUT );
        _var[ivar] = itX->second;
      }
    }
  //! @brief Constructor for cut w/ <a>n+1</a> linear variables
  PolCut
    ( FFOp*op, TYPE type, const double b, const unsigned n, const PolVar<T>*X,
      const double*a, const PolVar<T>&X1, const double a1 )
    : _op(op), _type(type), _rhs(b), _var(n+1), _coef(n+1)
    {
      for( unsigned ivar=0; ivar<n; ivar++ ){
        _coef[ivar] = a[ivar];
        _var[ivar] = X[ivar];
      }
      _coef[n] = a1;
      _var[n]  = X1;
    }
  //! @brief Constructor for cut w/ <a>n+2</a> linear variables
  PolCut
    ( FFOp*op, TYPE type, const double b, const unsigned n, const PolVar<T>*X,
      const double*a, const PolVar<T>&X1, const double a1, const PolVar<T>&X2,
      const double a2 )
    : _op(op), _type(type), _rhs(b), _var(n+2), _coef(n+2)
    {
      for( unsigned ivar=0; ivar<n; ivar++ ){
        _coef[ivar] = a[ivar];
        _var[ivar] = X[ivar];
      }
      _coef[n]   = a1;
      _var[n]    = X1;
      _coef[n+1] = a2;
      _var[n+1]  = X2;
    }
  //! @brief Constructor for cut w/ <a>n+3</a> linear variables
  PolCut
    ( FFOp*op, TYPE type, const double b, const unsigned n, const PolVar<T>*X,
      const double*a, const PolVar<T>&X1, const double a1, const PolVar<T>&X2,
      const double a2, const PolVar<T>&X3, const double a3 )
    : _op(op), _type(type), _rhs(b), _var(n+3), _coef(n+3)
    {
      for( unsigned ivar=0; ivar<n; ivar++ ){
        _coef[ivar] = a[ivar];
        _var[ivar] = X[ivar];
      }
      _coef[n]   = a1;
      _var[n]    = X1;
      _coef[n+1] = a2;
      _var[n+1]  = X2;
      _coef[n+2] = a3;
      _var[n+2]  = X3;
    }
  //! @brief Constructor for cut w/ selection among <a>n</a> linear variables plus 1 linear variable
  PolCut
    ( FFOp*op, TYPE type, const double b, const std::set<unsigned>&ndx,
      const PolVar<T>*X, const double*a, const PolVar<T>&X1, const double a1 )
    : _op(op), _type(type), _rhs(b), _var(ndx.size()+1), _coef(ndx.size()+1)
    {
      std::set<unsigned>::const_iterator it = ndx.begin();
      for( unsigned ivar=0; it!=ndx.end(); ++it, ivar++ ){
        _coef[ivar] = a[*it];
        _var[ivar] = X[*it];
      }
      _coef[ndx.size()] = a1;
      _var[ndx.size()]  = X1;
    }
  //! @brief Constructor for cut w/ variable and coefficient maps plus 1 linear variable
  template <typename U> PolCut
    ( FFOp*op, TYPE type, const double b, const std::map<U,PolVar<T>>&X,
      const std::map<U,double>&a, const PolVar<T>&X1, const double a1 )
    : _op(op), _type(type), _rhs(b), _var(a.size()+1), _coef(a.size()+1)
    {
      auto ita = a.cbegin();
      for( unsigned ivar=0; ita!=a.cend(); ++ita, ivar++ ){
        _coef[ivar] = ita->second;
        auto itX = X.find(ita->first);
        if( itX == X.cend() )
          throw typename PolImg<T>::Exceptions( PolImg<T>::Exceptions::BADCUT );
        _var[ivar] = itX->second;
      }
      _coef[a.size()] = a1;
      _var[a.size()]  = X1;
    }
  //! @brief Constructor for cut w/ variable and coefficient maps plus 2 linear variables
  template <typename U> PolCut
    ( FFOp*op, TYPE type, const double b, const std::map<U,PolVar<T>>&X,
      const std::map<U,double>&a, const PolVar<T>&X1, const double a1,
      const PolVar<T>&X2, const double a2 )
    : _op(op), _type(type), _rhs(b), _var(a.size()+2), _coef(a.size()+2)
    {
      auto ita = a.cbegin();
      for( unsigned ivar=0; ita!=a.cend(); ++ita, ivar++ ){
        _coef[ivar] = ita->second;
        auto itX = X.find(ita->first);
        if( itX == X.cend() )
          throw typename PolImg<T>::Exceptions( PolImg<T>::Exceptions::BADCUT );
        _var[ivar] = itX->second;
      }
      _coef[a.size()]   = a1;
      _var[a.size()]    = X1;
      _coef[a.size()+1] = a2;
      _var[a.size()+1]  = X2;
    }
  //! @brief Destructor
  ~PolCut
    ()
    {}

  //! @brief Append linear term to cut
  PolCut<T>& append
    ( const PolVar<T>&X, const double a )
    {
      _coef.push_back( a );
      _var.push_back( X );
      return *this;
    }
  //! @brief Append n linear terms to cut
  PolCut<T>& append
    ( const unsigned n, const PolVar<T>*X, const double*a )
    {
      _coef.insert( _coef.end(), a, a+n );
      _var.insert( _var.end(), X, X+n );
      return *this;
    }
  //! @brief Append quadratic term to cut
  PolCut<T>& append
    ( const PolVar<T>&X1, const PolVar<T>&X2, const double a )
    {
      _qcoef.push_back( a );
      _qvar1.push_back( X1 );
      _qvar2.push_back( X2 );
      return *this;
    }
  //! @brief Append n quadratic terms to cut
  PolCut<T>& append
    ( const unsigned n, const PolVar<T>*X1, const PolVar<T>*X2, const double*a )
    {
      _qcoef.insert( _qcoef.end(), a, a+n );
      _qvar1.insert( _qvar1.end(), X1, X1+n );
      _qvar2.insert( _qvar2.end(), X2, X2+n );
      return *this;
    }
  /** @} */

private:
  //! @brief Private methods to block default compiler methods
  PolCut
    ();
};

template <class T> 
inline std::ostream&
operator <<
( std::ostream&out, const PolCut<T>&cut )
{
  const int iprec = 5;
  out << std::right << std::scientific << std::setprecision(iprec);
  
  switch( cut._type ){
    case PolCut<T>::EQ: case PolCut<T>::LE: case PolCut<T>::GE:
      for(unsigned k=0; k<cut.nvar(); k++ ){
        if( isequal( cut._coef[k], 0. ) )
          out << " + " << std::setw(iprec+6) << 0.;
        else if( cut._coef[k] > 0. )
          out << " + " << std::setw(iprec+6) << cut._coef[k];
        else
          out << " - " << std::setw(iprec+6) << -cut._coef[k];
        out << cut._var[k].name();  
      }
      for(unsigned k=0; k<cut.nqvar(); k++ ){
        if( isequal( cut._qcoef[k], 0. ) )
          out << " + " << std::setw(iprec+6) << 0.;
        else if( cut._qcoef[k] > 0. )
          out << " + " << std::setw(iprec+6) << cut._qcoef[k];
        else
          out << " - " << std::setw(iprec+6) << -cut._qcoef[k];
        out << cut._qvar1[k].name() << "." << cut._qvar2[k].name();  
      }
      if( cut._type == PolCut<T>::EQ )
        out << " = ";
      else if( cut._type == PolCut<T>::LE )
        out << " <= ";
      else
        out << " >= ";
      out << std::setw(iprec+6) << cut._rhs;      
      break;

    case PolCut<T>::SOS1: case PolCut<T>::SOS2:
      out << " {";
      for(unsigned k=0; k<cut.nvar(); k++ )
        out << " " << cut._var[k].name();
      out << " }";
      if( cut._type == PolCut<T>::SOS1 )
        out << " SOS1";
      if( cut._type == PolCut<T>::SOS2 )
        out << " SOS2";
      break;

    case PolCut<T>::NLIN:
      out << " ";
      switch( cut._op->type ){
        case FFOp::IPOW:
          out << cut._var[0].name() << " = " << "IPOW(" << cut._var[1].name()
              << ", " << cut._op->pops[1]->num().n << ")"; break;
        case FFOp::DPOW:
          out << cut._var[0].name() << " = " << "DPOW(" << cut._var[1].name()
              << ", " << cut._op->pops[1]->num().x << ")"; break;
        case FFOp::CHEB:
          out << cut._var[0].name() << " = " << "CHEB(" << cut._var[1].name()
              << ", " << cut._op->pops[1]->num().n << ")"; break;
        case FFOp::SQR:
          out << cut._var[0].name() << " = " << "SQR(" << cut._var[1].name() << ")"; break;
        case FFOp::SQRT:
          out << cut._var[0].name() << " = " << "SQRT(" << cut._var[1].name() << ")"; break;
        case FFOp::EXP:
          out << cut._var[0].name() << " = " << "EXP(" << cut._var[1].name() << ")"; break;
        case FFOp::LOG:
          out << cut._var[0].name() << " = " << "LOG(" << cut._var[1].name() << ")"; break;
        case FFOp::XLOG:
          out << cut._var[0].name() << " = " << "XLOG(" << cut._var[1].name() << ")"; break;
        case FFOp::FABS:
          out << cut._var[0].name() << " = " << "FABS(" << cut._var[1].name() << ")"; break;
        case FFOp::FSTEP:
          out << cut._var[0].name() << " = " << "FSTEP(" << cut._var[1].name() << ")"; break;
        case FFOp::COS:
          out << cut._var[0].name() << " = " << "COS(" << cut._var[1].name() << ")"; break;
        case FFOp::SIN:
          out << cut._var[0].name() << " = " << "SIN(" << cut._var[1].name() << ")"; break;
        case FFOp::TAN:
          out << cut._var[0].name() << " = " << "TAN(" << cut._var[1].name() << ")"; break;
        case FFOp::ACOS:
          out << cut._var[0].name() << " = " << "ACOS(" << cut._var[1].name() << ")"; break;
        case FFOp::ASIN:
          out << cut._var[0].name() << " = " << "ASIN(" << cut._var[1].name() << ")"; break;
        case FFOp::ATAN:
          out << cut._var[0].name() << " = " << "ATAN(" << cut._var[1].name() << ")"; break;
        case FFOp::COSH:
          out << cut._var[0].name() << " = " << "COSH(" << cut._var[1].name() << ")"; break;
        case FFOp::SINH:
          out << cut._var[0].name() << " = " << "SINH(" << cut._var[1].name() << ")"; break;
        case FFOp::TANH:
          out << cut._var[0].name() << " = " << "TANH(" << cut._var[1].name() << ")"; break;
        case FFOp::ERF:
          out << cut._var[0].name() << " = " << "ERF(" << cut._var[1].name() << ")"; break;
        case FFOp::MINF:
          out << cut._var[0].name() << " = " << "MIN(" << cut._var[1].name() << ", ";
          if( cut._var.size() > 2 ) out << cut._var[2].name();
          else                      out << cut._rhs;
          out << ")"; break;
        case FFOp::MAXF:
          out << cut._var[0].name() << " = " << "MAX(" << cut._var[1].name() << ", ";
          if( cut._var.size() > 2 ) out << cut._var[2].name();
          else                      out << cut._rhs;
          out << ")"; break;
        case FFOp::LMTD:
          out << cut._var[0].name() << " = " << "LMTD(" << cut._var[1].name() << ", "
              << cut._var[2].name() << ")"; break;
        case FFOp::RLMTD:
          out << cut._var[0].name() << " = " << "RMLTD(" << cut._var[1].name() << ", "
              << cut._var[2].name() << ")"; break;
        default:
          throw typename PolImg<T>::Exceptions( PolImg<T>::Exceptions::BADCUT );
      } 
    }
  return out;
}

//! @brief C++ structure for ordering of polytopic cuts
template <class T> 
struct lt_PolCut
{
  bool operator()
    ( const PolCut<T>*Cut1, const PolCut<T>*Cut2 ) const
    {
      // Order cuts w.r.t. their types first
      if( Cut1->type() < Cut2->type() ) return true;
      // then w.r.t. their defining operation next
      return ( Cut1->op() && Cut2->op() ?
        lt_FFOp()( Cut1->op(), Cut2->op() ):
        false );
    }
};

//! @brief C++ structure for storing subintervals in the outer approximation of a univariate term
////////////////////////////////////////////////////////////////////////
//! mc::OAsub is a C++ structure for storing subintervals in the
//! outer approximation of a convex/concave portion of a univariate
//! function.
////////////////////////////////////////////////////////////////////////
class OAsub
////////////////////////////////////////////////////////////////////////
{

public:
  //! @brief Constructor
  OAsub( const double LB, const double UB, const double MID, 
    const double GAP ):
    _xL( LB ), _xU( UB ), _xM( MID ), _gap( GAP )
    {}
  //! @brief Retreive interval lower bound
  const double xL() const
    { return _xL; }
  //! @brief Retreive interval upper bound
  const double xU() const
    { return _xU; }
  //! @brief Retreive bisection point
  const double xM() const
    { return _xM; }
  //! @brief Retreive maximum gap
  const double gap() const
    { return _gap; }

private:
  //! @brief Interval lower bound
  double _xL;
  //! @brief Interval upper bound
  double _xU;
  //! @brief Bisection point
  double _xM;
  //! @brief Maximum gap
  double _gap;
};

//! @brief C++ structure for comparison of subintervals in the outer approximation of a univariate term
struct lt_OAsub
{
  bool operator()
    ( const OAsub&Dom1, const OAsub&Dom2 ) const
    {
      return( Dom1.gap() < Dom2.gap() );
    }
};

//! @brief C++ class for polytopic image evaluation
////////////////////////////////////////////////////////////////////////
//! mc::PolImg is a C++ class for evaluation of the polytoptic image of
//! a factorable function. Propagation of the image is via polytopic 
//! arithmetic, as implemented in mc::PolVar. The template parameter
//! corresponds to the type used to propagate variable range. Round-off
//! errors are not accounted for in the computations (non-verified
//! implementation).
////////////////////////////////////////////////////////////////////////
template< class T >
class PolImg
////////////////////////////////////////////////////////////////////////
{
  friend class PolVar<T>;
 
  template <typename U> friend std::ostream& operator<<( std::ostream&, const PolImg<U>& );

  template< class U > friend  PolVar<U> operator+( const PolVar<U>&, const PolVar<U>& );
  template< class U > friend  PolVar<U> operator+( const PolVar<U>&, const double );
  template< class U > friend  PolVar<U> operator-( const PolVar<U>& );
  template< class U > friend  PolVar<U> operator-( const PolVar<U>&, const PolVar<U>& );
  template< class U > friend  PolVar<U> operator*( const PolVar<U>&, const PolVar<U>& );
  template< class U > friend  PolVar<U> operator*( const PolVar<U>&, const double );
  template< class U > friend  PolVar<U> operator/( const PolVar<U>&, const PolVar<U>& );
  template< class U > friend  PolVar<U> operator/( const double, const PolVar<U>& );
  template< class U > friend  PolVar<U> operator^( const PolVar<U>&, const PolVar<U>& );

  template< class U > friend  PolVar<U> inv ( const PolVar<U>& );
  template< class U > friend  PolVar<U> exp ( const PolVar<U>& );
  template< class U > friend  PolVar<U> log ( const PolVar<U>& );
  template< class U > friend  PolVar<U> xlog( const PolVar<U>& );
  template< class U > friend  PolVar<U> sqrt( const PolVar<U>& );
  template< class U > friend  PolVar<U> sqr ( const PolVar<U>& );
  template< class U > friend  PolVar<U> pow ( const PolVar<U>&, const int );  
  template< class U > friend  PolVar<U> pow ( const PolVar<U>&, const double );  
  template< class U > friend  PolVar<U> cheb( const PolVar<U>&, const unsigned );  
  template< class U > friend  PolVar<U> prod( const unsigned, const PolVar<U>* );  
  template< class U > friend  PolVar<U> cos ( const PolVar<U>& );
  template< class U > friend  PolVar<U> sin ( const PolVar<U>& );
  template< class U > friend  PolVar<U> tan ( const PolVar<U>& );
  template< class U > friend  PolVar<U> acos( const PolVar<U>& );
  template< class U > friend  PolVar<U> asin( const PolVar<U>& );
  template< class U > friend  PolVar<U> atan( const PolVar<U>& );
  template< class U > friend  PolVar<U> cosh( const PolVar<U>& );
  template< class U > friend  PolVar<U> sinh( const PolVar<U>& );
  template< class U > friend  PolVar<U> tanh( const PolVar<U>& );
  template< class U > friend  PolVar<U> fabs( const PolVar<U>& );
  template< class U > friend  PolVar<U> erf( const PolVar<U>& );
  template< class U > friend  PolVar<U> fstep( const PolVar<U>& );
  template< class U > friend  PolVar<U> max( const PolVar<U>&, const PolVar<U>& );  
  template< class U > friend  PolVar<U> min( const PolVar<U>&, const PolVar<U>& );  
  template< class U > friend  PolVar<U> lmtd( const PolVar<U>&, const PolVar<U>& );  
  template< class U > friend  PolVar<U> rlmtd( const PolVar<U>&, const PolVar<U>& );  

public:
  typedef std::map< FFVar const*, PolVar<T>*, lt_FFVar > t_Vars;
  typedef typename t_Vars::iterator it_Vars;
  typedef std::list< PolVar<T>* > t_Aux;
  typedef std::list< PolLQExpr<T>* > t_LQExpr;
  typedef std::multiset< PolCut<T>*, lt_PolCut<T> > t_Cuts;
  typedef double (*p_Univ)( const double, const double*, const int* );
  typedef std::pair<double,double> (*p_dUniv)( const double, const double*, const int* );
  typedef void (PolImg<T>::*p_Cut)( FFOp*, const double, const PolVar<T>&, const double,
    const double, const PolVar<T>&, const double, const double, const double*, const int* );
  typedef void (PolImg<T>::*p_Cut2)( FFOp*, const PolVar<T>&, const double,
    const double, const PolVar<T>&, const double, const double, const double*, const int* );
  typedef std::priority_queue< OAsub, std::vector<OAsub>, lt_OAsub > t_OA;

protected:
  //! @brief Map of DAG variables in polytopic image
  t_Vars _Vars;
  //! @brief Appends new pointer to DAG variable map in polytopic image
  PolVar<T>* _append_var
    ( FFVar const* var, T const& range, bool const cont );
  //! @brief Update range of variable in polytopic image
  PolVar<T>* _update_var
    ( FFVar const* var, T const& range );
  //! @brief Update type of variable in polytopic image
  PolVar<T>* _update_var
    ( FFVar const* var, bool const cont );
  //! @brief Erase all entries in _Vars
  void _erase_vars
    ();
  //! @brief Reset cut-related field _subdiv in _Vars
  void _reset_vars
    ();

  //! @brief List of auxiliary variables in polytopic image
  t_Aux _Aux;
  //! @brief Appends new pointer to auxiliary variable set in polytopic image
  PolVar<T>* _append_aux
    ( const T&range, const bool cont );
  //! @brief Erase all entries in _Aux
  void _erase_aux
    ();

  //! @brief List of linear-quadratic expressions in polytopic image
  t_LQExpr _LQExpr;
  //! @brief Appends new term to linear-quadratic expression in polytopic image
  PolLQExpr<T>* _append_LQ
    ( const PolVar<T>*varR );
  //! @brief Erase all linear-quadratic expressions in _LQExpr
  void _erase_LQ
    ();

  //! @brief Set of cuts in polytopic image
  t_Cuts _Cuts;
  //! @brief Erase all entries in _Cuts
  void _erase_cuts
    ();
  //! @brief Erase all entries corresponding to operation <a>op</a> in _Cuts
  void _erase_cuts
    ( FFOp* op );

  //! @brief Appends new nonlinear cut in _Cuts w/ 1 variable
  typename t_Cuts::iterator _add_cut
    ( FFOp*op, const PolVar<T>&X1, const PolVar<T>&X2 );
  //! @brief Appends new nonlinear cut in _Cuts w/ 2 variables
  typename t_Cuts::iterator _add_cut
    ( FFOp*op, const PolVar<T>&X1, const PolVar<T>&X2, const PolVar<T>&X3 );
  //! @brief Appends new nonlinear cut in _Cuts w/ 1 variable + 1 constant
  typename t_Cuts::iterator _add_cut
    ( FFOp*op, const PolVar<T>&X1, const PolVar<T>&X2, const double X3 );
  //! @brief Appends new relaxation cut in _Cuts w/o variable
  typename t_Cuts::iterator _add_cut
    ( FFOp*op, const typename PolCut<T>::TYPE type, const double b );
  //! @brief Appends new relaxation cut in _Cuts w/ 1 variable
  typename t_Cuts::iterator _add_cut
    ( FFOp*op, const typename PolCut<T>::TYPE type,
      const double b, const PolVar<T>&X1, const double a1 );
  //! @brief Appends new relaxation cut in _Cuts w/ 2 variables
  typename t_Cuts::iterator _add_cut
    ( FFOp*op, const typename PolCut<T>::TYPE type,
      const double b, const PolVar<T>&X1, const double a1,
      const PolVar<T>&X2, const double a2 );
  //! @brief Appends new relaxation cut in _Cuts w/ 3 variables
  typename t_Cuts::iterator _add_cut
    ( FFOp*op, const typename PolCut<T>::TYPE type,
      const double b, const PolVar<T>&X1, const double a1,
      const PolVar<T>&X2, const double a2, const PolVar<T>&X3,
      const double a3 );
  //! @brief Appends new relaxation cut in _Cuts w/ 4 variables
  typename t_Cuts::iterator _add_cut
    ( FFOp*op, const typename PolCut<T>::TYPE type,
      const double b, const PolVar<T>&X1, const double a1,
      const PolVar<T>&X2, const double a2, const PolVar<T>&X3,
      const double a3, const PolVar<T>&X4, const double a4 );
  //! @brief Appends new relaxation cut in _Cuts w/ <a>n</a> variables
  typename t_Cuts::iterator _add_cut
    ( FFOp*op, const typename PolCut<T>::TYPE type, const double b,
      const unsigned n, const PolVar<T>*X, const double*a );
  //! @brief Appends new relaxation cut in _Cuts w/ <a>n</a> variables
  typename t_Cuts::iterator _add_cut
    ( FFOp*op, const typename PolCut<T>::TYPE type, const double b,
      const std::set<unsigned>&ndx, const PolVar<T>*X, const double*a );
  //! @brief Appends new relaxation cut in _Cuts w/ variable and coefficient maps
  template <typename U> typename t_Cuts::iterator _add_cut
    ( FFOp*op, const typename PolCut<T>::TYPE type, const double b,
      const std::map<U,PolVar<T>>&X, const std::map<U,double>&a );
  //! @brief Appends new relaxation cut in _Cuts w/ <a>n+1</a> variables
  typename t_Cuts::iterator _add_cut
    ( FFOp*op, const typename PolCut<T>::TYPE type, const double b,
      const unsigned n, const PolVar<T>*X, const double*a,
      const PolVar<T>&X1, const double a1 );
  //! @brief Appends new relaxation cut in _Cuts w/ <a>n+1</a> variables
  typename t_Cuts::iterator _add_cut
    ( FFOp*op, const typename PolCut<T>::TYPE type, const double b,
      const unsigned n, const PolVar<T>*X, const double*a,
      const PolVar<T>&X1, const double a1,
      const PolVar<T>&X2, const double a2 );
  //! @brief Appends new relaxation cut in _Cuts w/ <a>n+1</a> variables
  typename t_Cuts::iterator _add_cut
    ( FFOp*op, const typename PolCut<T>::TYPE type, const double b,
      const unsigned n, const PolVar<T>*X, const double*a,
      const PolVar<T>&X1, const double a1,
      const PolVar<T>&X2, const double a2,
      const PolVar<T>&X3, const double a3 );
  //! @brief Appends new relaxation cut in _Cuts w/ <a>n+1</a> variables
  typename t_Cuts::iterator _add_cut
    ( FFOp*op, const typename PolCut<T>::TYPE type, const double b,
      const std::set<unsigned>&ndx, const PolVar<T>*X, const double*a,
      const PolVar<T>&X1, const double a1 );
  //! @brief Appends new relaxation cut in _Cuts w/ 1 variable and coefficient maps
  template <typename U> typename t_Cuts::iterator _add_cut
    ( FFOp*op, const typename PolCut<T>::TYPE type, const double b,
      const std::map<U,PolVar<T>>&X, const std::map<U,double>&a,
      const PolVar<T>&X1, const double a1 );
  //! @brief Appends new relaxation cut in _Cuts w/ 2 variables and coefficient maps
  template <typename U> typename t_Cuts::iterator _add_cut
    ( FFOp*op, const typename PolCut<T>::TYPE type, const double b,
      const std::map<U,PolVar<T>>&X, const std::map<U,double>&a,
      const PolVar<T>&X1, const double a1,
      const PolVar<T>&X2, const double a2 );
        //! @brief Appends new relaxation cut in _Cuts w/ 3 variables and coefficient maps

  //! @brief Computes max distance between function and outer-approximation
  std::pair< double, double > _distmax
    ( p_dUniv f, const double xL, const double xU, const double*rpar=0,
      const int*ipar=0 ) const;
  //! @brief Computes solution of scalar nonlinear equation using the Newton method
  double _newton
    ( const double x0, const double xL, const double xU, p_dUniv f,
      const double TOL, const unsigned MAXIT, const double*rusr=0, const int*iusr=0 ) const;
  //! @brief Computes solution of scalar nonlinear equation using the secant method
  double _secant
    ( const double x0, const double x1, const double xL, const double xU, p_Univ f,
      const double TOL, const unsigned MAXIT, const double*rusr=0, const int*iusr=0 ) const;
  //! @brief Computes solution of scalar nonlinear equation using the Golden-section method
  double _goldsect
    ( const double xL, const double xU, p_Univ f, const double TOL, const unsigned MAXIT,
      const double*rusr=0, const int*iusr=0 ) const;
  //! @brief Iteration of the Golden-section method
  double _goldsect_iter
    ( const bool init, const double a, const double fa, const double b,
      const double fb, const double c, const double fc, p_Univ f,
      const double TOL, const unsigned MAXIT, const double*rusr=0, const int*iusr=0 ) const;
  //! @brief Append cuts for nonlinear operation using outer-approximation (sandwich algorithm)
  void _sandwich_cuts
    ( FFOp*pOp, const PolVar<T>&X, const double XL, const double XU, const PolVar<T>&Y,
      const double YL, const double YU, const typename PolCut<T>::TYPE sense, p_dUniv f,
      const double*rpar=0, const int*ipar=0 );
  //! @brief Append linearization cut for nonlinear operation using outer-approximation
  void _linearization_cut
    ( FFOp*pOp, const double Xref, const PolVar<T>&X, const double XL, const double XU,
      const PolVar<T>&Y, const double YL, const double YU, const typename PolCut<T>::TYPE sense,
      p_dUniv f, const double*rpar, const int*ipar );

  //! @brief Append cuts for nonlinear univariate using piecewise-linear approximation
  void _semilinear_cuts
    ( FFOp*pOp, const PolVar<T>&X, const double XL, const double XU, const PolVar<T>&Y,
      const typename PolCut<T>::TYPE sense, p_dUniv f, const double*rpar=0, const int*ipar=0 );
  //! @brief Form subintervals for semilinear cuts
  std::vector<double>* _semilinear_sub
    ( const PolVar<T>&X, const double XL, const double XU );

  //! @brief Append cuts for bilinear term using piecewise-linear approximation
  bool _pwmccormick_cuts
    ( FFOp*pOp, const PolVar<T>&X1, const double X1L, const double X1U,
      const PolVar<T>&X2, const double X2L, const double X2U, const PolVar<T>&Y );

  //! @brief Propagate linear cut for binary + operation
  bool _add_LQ_PLUS
    ( PolLQExpr<T>*&pLQ, const PolVar<T>*VarR, FFVar*pVar1, FFVar*pVar2 );
  //! @brief Propagate linear cut for unary - operation
  bool _add_LQ_NEG
    ( PolLQExpr<T>*&pLQ, const PolVar<T>*VarR, FFVar*pVar );
  //! @brief Propagate linear cut for binary - operation
  bool _add_LQ_MINUS
    ( PolLQExpr<T>*&pLQ, const PolVar<T>*VarR, FFVar*pVar1, FFVar*pVar2 );
  //! @brief Propagate linear cut for binary * operation
  bool _add_LQ_TIMES
    ( PolLQExpr<T>*&pLQ, const PolVar<T>*VarR, FFVar*pVar1, FFVar*pVar2 );
  //! @brief Propagate linear cut for binary / operation
  bool _add_LQ_DIV
    ( PolLQExpr<T>*&pLQ, const PolVar<T>*VarR, FFVar*pVar1, FFVar*pVar2 );
  //! @brief Propagate linear cut for fabs function
  bool _add_LQ_FABS
    ( PolLQExpr<T>*&pLQ, const PolVar<T>*VarR, FFVar*pVar );
  //! @brief Propagate linear cut for fstep function
  bool _add_LQ_FSTEP
    ( PolLQExpr<T>*&pLQ, const PolVar<T>*VarR, FFVar*pVar );
  //! @brief Propagate linear cut for binary min function
  bool _add_LQ_MINF
    ( PolLQExpr<T>*&pLQ, const PolVar<T>*VarR, FFVar*pVar1, FFVar*pVar2 );
  //! @brief Propagate linear cut for binary max function
  bool _add_LQ_MAXF
    ( PolLQExpr<T>*&pLQ, const PolVar<T>*VarR, FFVar*pVar1, FFVar*pVar2 );
  //! @brief Propagate linear cut for unary - operation
  bool _add_LQ_SQR
    ( PolLQExpr<T>*&pLQ, const PolVar<T>*VarR, FFVar*pVar );

  //! @brief Append linear cuts for binary intersection
  void _add_cuts_INTER
    ( const PolVar<T>*VarR, FFVar*pVar1, FFVar*pVar2 );
  //! @brief Append linear cuts for binary + operation
  void _add_cuts_PLUS
    ( const PolVar<T>*VarR, FFVar*pVar1, FFVar*pVar2 );
  //! @brief Append linear cuts for unary - operation
  void _add_cuts_NEG
    ( const PolVar<T>*VarR, FFVar*pVar );
  //! @brief Append linear cuts for binary - operation
  void _add_cuts_MINUS
    ( const PolVar<T>*VarR, FFVar*pVar1, FFVar*pVar2 );
  //! @brief Append linear cuts for binary * operation
  void _add_cuts_TIMES
    ( const PolVar<T>*VarR, FFVar*pVar1, FFVar*pVar2, FFOp*pOp );
  //! @brief Append linear cuts for binary * operation
  void _add_cuts_TIMES
    ( const PolVar<T>*VarR, FFVar*pVar1, FFVar*pVar2 );
  //! @brief Append linear cuts for binary / operation
  void _add_cuts_DIV
    ( const PolVar<T>*VarR, FFVar*pVar1, FFVar*pVar2 );
  //! @brief Append linear cuts for ipow function
  void _add_cuts_IPOW
    ( const PolVar<T>*VarR, FFVar*pVar, const int iExp );
  //! @brief Append linear cuts for dpow function
  void _add_cuts_DPOW
    ( const PolVar<T>*VarR, FFVar*pVar, const double dExp );
  //! @brief Append linear cuts for cheb function
  void _add_cuts_CHEB
    ( const PolVar<T>*VarR, FFVar*pVar, const unsigned iOrd );
  //! @brief Append linear cuts for n-ary product
  void _add_cuts_PROD
    ( const PolVar<T>*VarR, const std::vector<FFVar*>vVar );
  //! @brief Append linear cuts for exp function
  void _add_cuts_EXP
    ( const PolVar<T>*VarR, FFVar*pVar );
  //! @brief Append linear cuts for log function
  void _add_cuts_LOG
    ( const PolVar<T>*VarR, FFVar*pVar );
  //! @brief Append linear cuts for xlog function
  void _add_cuts_XLOG
    ( const PolVar<T>*VarR, FFVar*pVar );
  //! @brief Append linear cuts for erf function
  void _add_cuts_ERF
    ( const PolVar<T>*VarR, FFVar*pVar );
  //! @brief Append linear cuts for sqr function
  void _add_cuts_SQR
    ( const PolVar<T>*VarR, FFVar*pVar, FFOp*pOp );
  //! @brief Append linear cuts for sqr function
  void _add_cuts_SQR
    ( const PolVar<T>*VarR, FFVar*pVar );
  //! @brief Append linear cuts for sqrt function
  void _add_cuts_SQRT
    ( const PolVar<T>*VarR, FFVar*pVar );
  //! @brief Append linear cuts for fabs function
  void _add_cuts_FABS
    ( const PolVar<T>*VarR, FFVar*pVar );
  //! @brief Append linear cuts for fstep function
  void _add_cuts_FSTEP
    ( const PolVar<T>*VarR, FFVar*pVar );
  //! @brief Append linear cuts for cos function
  void _add_cuts_COS
    ( const PolVar<T>*VarR, FFVar*pVar );
  //! @brief Append linear cuts for sin function
  void _add_cuts_SIN
    ( const PolVar<T>*VarR, FFVar*pVar );
  //! @brief Append linear cuts for tan function
  void _add_cuts_TAN
    ( const PolVar<T>*VarR, FFVar*pVar );
  //! @brief Append linear cuts for acos function
  void _add_cuts_ACOS
    ( const PolVar<T>*VarR, FFVar*pVar );
  //! @brief Append linear cuts for asin function
  void _add_cuts_ASIN
    ( const PolVar<T>*VarR, FFVar*pVar );
  //! @brief Append linear cuts for atan function
  void _add_cuts_ATAN
    ( const PolVar<T>*VarR, FFVar*pVar );
  //! @brief Append linear cuts for cosh function
  void _add_cuts_COSH
    ( const PolVar<T>*VarR, FFVar*pVar );
  //! @brief Append linear cuts for sinh function
  void _add_cuts_SINH
    ( const PolVar<T>*VarR, FFVar*pVar );
  //! @brief Append linear cuts for tanh function
  void _add_cuts_TANH
    ( const PolVar<T>*VarR, FFVar*pVar );

  //! @brief Append semi-linear cuts for binary min function
  void _add_cuts_MINF
    ( const PolVar<T>*VarR, FFVar*pVar1, FFVar*pVar2 );
  //! @brief Append semi-linear cuts for binary max function
  void _add_cuts_MAXF
    ( const PolVar<T>*VarR, FFVar*pVar1, FFVar*pVar2 );

  //! @brief Increment index even subset
  bool _subset_incr
    ( const unsigned ntot, std::vector<unsigned>&ndx );

public:
  /** @ingroup POLYTOPE
   *  @{
   */
  //! @brief Default Constructor
  PolImg
    ()
    {}

  //! @brief Destructor
  virtual ~PolImg
    ()
    { reset(); }

  //! @brief PolImg exceptions
  class Exceptions
  {
  public:
    //! @brief Enumeration type for exception handling
    enum TYPE{
      ROOT=1,        //!< Error during root search for obtaining the convex/concave envelope of a univariate term
      INTER,         //!< Error during intersection of two terms (terms do not intersect)
      DIV,           //!< Error during division operation (division by 0)
      ENVMIS=-1,     //!< Error due to an operation between variables participating in different polytopic images
      BADCUT=-2,     //!< Error due to an error during cut generation
      UNAVAIL=-3,    //!< Error due to calling a function/feature not yet implemented in MC++
      NOTALLOWED=-4, //!< Error due to calling a function/feature not yet implemented in MC++
      INTERN=-5	     //!< Internal error
    };
    //! @brief Constructor for error <a>ierr</a>
    Exceptions( TYPE ierr ) : _ierr( ierr ){}
    //! @brief Inline function returning the error flag
    int ierr(){ return _ierr; }

  private:
    TYPE _ierr;
  };

  //! @brief PolImg options class
  struct Options
  {
    //! @brief Constructor
    Options():
      AGGREG_LQ(false), ROOT_USE(true), ROOT_MAXIT(100), ROOT_TOL(1e-10),
      SANDWICH_ATOL(1e-10), SANDWICH_RTOL(1e-3), SANDWICH_MAXCUT(5),
      SANDWICH_RULE(MAXERR), FRACTIONAL_ATOL(machprec()),
      FRACTIONAL_RTOL(machprec()), BREAKPOINT_TYPE(NONE),
      BREAKPOINT_ATOL(1e-5), BREAKPOINT_RTOL(1e-3),
      RELAX_QUAD(true), RELAX_MONOM(1), RELAX_NLIN(true), RELAX_DISC(2)
      {}
    //! @brief Assignment operator
    Options& operator= ( const Options&options ){
        AGGREG_LQ      = options.AGGREG_LQ;
        ROOT_USE        = options.ROOT_USE;
        ROOT_MAXIT      = options.ROOT_MAXIT;
        ROOT_TOL        = options.ROOT_TOL;
        SANDWICH_ATOL   = options.SANDWICH_ATOL;
        SANDWICH_RTOL   = options.SANDWICH_RTOL;
        SANDWICH_MAXCUT = options.SANDWICH_MAXCUT;
        SANDWICH_RULE   = options.SANDWICH_RULE;
        FRACTIONAL_ATOL = options.FRACTIONAL_ATOL;
        FRACTIONAL_RTOL = options.FRACTIONAL_RTOL;
        BREAKPOINT_TYPE = options.BREAKPOINT_TYPE;
        BREAKPOINT_ATOL = options.BREAKPOINT_ATOL;
        BREAKPOINT_RTOL = options.BREAKPOINT_RTOL;
        RELAX_QUAD      = options.RELAX_QUAD;
        RELAX_MONOM     = options.RELAX_MONOM;
        RELAX_NLIN      = options.RELAX_NLIN;
        RELAX_DISC      = options.RELAX_DISC;
        return *this;
      }
    //! @brief Enumeration type for sandwich strategy
    enum SANDWICH{
      BISECT=0,	//!< Range bisection
      MAXERR	//!< Maximum error rule
    };
    //! @brief Enumeration type for bilinear term relaxation strategy
    enum REFINE{
      NONE=0,	//!< No semi-linear cuts (use secant approximation)
      BIN,	//!< Semilinear cuts as linear binary reformulation
      SOS2	//!< Semilinear cuts as SOS2 reformulation
    };
    //! @brief Whether or not to aggregate linear expressions in cuts - Default: true
    bool AGGREG_LQ;
    //! @brief Whether or not to use root search to construct envelopes of univariate terms - Default: true
    bool ROOT_USE;
    //! @brief Maximal number of iterations in root search - Default: 100
    unsigned ROOT_MAXIT;
    //! @brief Termination tolerance in root search - Default: 1e-10
    double ROOT_TOL;
    //! @brief Absolute tolerance in outer-approximation of univariate terms - Default: 1e-3
    double SANDWICH_ATOL;
    //! @brief Relative tolerance in outer-approximation of univariate terms - Default: 1e-3
    double SANDWICH_RTOL;
    //! @brief Maximal number of cuts in outer approximation of univariate terms - Default: 5
    unsigned SANDWICH_MAXCUT;
    //! @brief Rule for outer-approximation of nonlinear convex/concave terms - Default: MAXERR
    SANDWICH SANDWICH_RULE;
    //! @brief Absolute tolerance to prevent division by zero - Default: EPS
    double FRACTIONAL_ATOL;
    //! @brief Relative tolerance to prevent division by zero - Default: EPS
    double FRACTIONAL_RTOL;
    //! @brief Rule for piecewise linear cuts of nonlinear convex/concave terms - Default: NONE
    REFINE BREAKPOINT_TYPE;
    //! @brief Absolute tolerance in adding breakpoints in piecewise linear cuts - Default: 1e-8
    double BREAKPOINT_ATOL;
    //! @brief Relative tolerance in adding breakpoints in piecewise linear cuts - Default: 1e-5
    double BREAKPOINT_RTOL;
    //! @brief Whether or not to linearize quadratic terms (i.e. square, bilinear) - Default: true
    bool RELAX_QUAD;
    //! @brief Whether or not to linearize monomial terms (i.e. pow, cheb) - Default: 1 - monomials of order >2 linearized & quadratic term linearization controlled by RELAX_QUAD (2 - all monomials linearized ; 0 - no monomials linearized)
    int RELAX_MONOM;
    //! @brief Whether or not to linearize nonlinear terms (i.e. exp, log, sin, cos, tan, dpow) terms - Default: true
    bool RELAX_NLIN;
    //! @brief Whether or not to linearize discontinuous terms (i.e. min, max, abs, fstep) terms - Default: 2 - mixed-integer linearization (1 - continuous linearization; 0 - no linearization )
    int RELAX_DISC;
  };

  //! @brief PolImg options handle
  Options options;

  //! @brief Retreive reference to set of DAG variables in polytopic image
  t_Vars& Vars()
    { return _Vars; }

  //! @brief Retreive const reference to set of DAG variables in polytopic image
  const t_Vars& Vars() const
    { return _Vars; }
  
  //! @brief Retreive reference to set of auxiliary variables in polytopic image
  t_Aux& Aux()
    { return _Aux; }

  //! @brief Retreive reference to set of cuts in polytopic image
  t_Cuts& Cuts()
    { return _Cuts; }

  //! @brief Reset polytopic image (all participating variables and cuts)
  void reset()
    { _erase_cuts(); _erase_vars(); _erase_aux(); _erase_LQ(); }

  //! @brief Append relaxation cuts for the <a>ndep</a> dependents in <a>pdep</a>
  void generate_cuts
    ( const unsigned ndep, const PolVar<T>*pdep, const bool reset=false );
  //! @brief Append relaxation cuts for the dependents in <a>pdep</a> indexed by <a>ndxdep</a>
  void generate_cuts
    ( const std::set<unsigned>&ndxdep, const PolVar<T>*pdep, const bool reset=false );
  //! @brief Append relaxation cuts for the dependents in the map <a>mdep</a>
  template <typename U> void generate_cuts
    ( const std::map<U,PolVar<T>>&mdep, const bool reset=false );

  //! @brief Append new nonlinear cut w/ 1 variable
  typename t_Cuts::iterator add_cut
    ( const PolVar<T>&X1, const PolVar<T>&X2 )
    { return _add_cut( 0, X1, X2 ); }
  //! @brief Append new nonlinear cut w/ 2 variables
  typename t_Cuts::iterator add_cut
    ( const PolVar<T>&X1, const PolVar<T>&X2, const PolVar<T>&X3 )
    { return _add_cut( 0, X1, X2, X3 ); }
  //! @brief Append new nonlinear cut w/ 1 variable + 1 constant
  typename t_Cuts::iterator add_cut
    ( const PolVar<T>&X1, const PolVar<T>&X2, const double X3 )
    { return _add_cut( 0, X1, X2, X3 ); }
  //! @brief Append new relaxation cut w/o variable
  typename t_Cuts::iterator add_cut
    ( const typename PolCut<T>::TYPE type, const double b=0. )
    { return _add_cut( 0, type, b ); }
  //! @brief Append new relaxation cut w/ 1 variable
  typename t_Cuts::iterator add_cut
    ( const typename PolCut<T>::TYPE type, const double b,
      const PolVar<T>&X1, const double a1 )
    { return _add_cut( 0, type, b, X1, a1 ); }
  //! @brief Append new relaxation cut w/ 2 variables
  typename t_Cuts::iterator add_cut
    ( const typename PolCut<T>::TYPE type, const double b,
      const PolVar<T>&X1, const double a1,
      const PolVar<T>&X2, const double a2 )
    { return _add_cut( 0, type, b, X1, a1, X2, a2 ); }
  //! @brief Append new relaxation cut w/ 3 variables
  typename t_Cuts::iterator add_cut
    ( const typename PolCut<T>::TYPE type, const double b,
      const PolVar<T>&X1, const double a1,
      const PolVar<T>&X2, const double a2,
      const PolVar<T>&X3, const double a3 )
    { return _add_cut( 0, type, b, X1, a1, X2, a2, X3, a3 ); }
  //! @brief Append new relaxation cut w/ 4 variables
  typename t_Cuts::iterator add_cut
    ( const typename PolCut<T>::TYPE type, const double b,
      const PolVar<T>&X1, const double a1,
      const PolVar<T>&X2, const double a2,
      const PolVar<T>&X3, const double a3,
      const PolVar<T>&X4, const double a4 )
    { return _add_cut( 0, type, b, X1, a1, X2, a2, X3, a3, X4, a4 ); }
  //! @brief Append new relaxation cut w/ <a>n</a> variables
  typename t_Cuts::iterator add_cut
    ( const typename PolCut<T>::TYPE type, const double b,
      const unsigned n, const PolVar<T>*X, const double*a )
    { return _add_cut( 0, type, b, n, X, a ); }
  //! @brief Append new relaxation cut w/ selection amongst <a>n</a> variables
  typename t_Cuts::iterator add_cut
    ( const typename PolCut<T>::TYPE type, const double b,
      const std::set<unsigned>&ndx, const PolVar<T>*X, const double*a )
    { return _add_cut( 0, type, b, ndx, X, a ); }
  //! @brief Append new relaxation cut w/ variable and coefficient maps
  template <typename U> typename t_Cuts::iterator add_cut
    ( const typename PolCut<T>::TYPE type, const double b,
      const std::map<U,PolVar<T>>&X, const std::map<U,double>&a )
    { return _add_cut( 0, type, b, X, a ); }
  //! @brief Append new relaxation cut w/ <a>n+1</a> variables
  typename t_Cuts::iterator add_cut
    ( const typename PolCut<T>::TYPE type, const double b,
      const unsigned n, const PolVar<T>*X, const double*a,
      const PolVar<T>&X1, const double a1 )
    { return _add_cut( 0, type, b, n, X, a, X1, a1 ); }
  //! @brief Append new relaxation cut w/ <a>n+2</a> variables
  typename t_Cuts::iterator add_cut
    ( const typename PolCut<T>::TYPE type, const double b,
      const unsigned n, const PolVar<T>*X, const double*a,
      const PolVar<T>&X1, const double a1,
      const PolVar<T>&X2, const double a2 )
    { return _add_cut( 0, type, b, n, X, a, X1, a1, X2, a2 ); }
  //! @brief Append new relaxation cut w/ <a>n+3</a> variables
  typename t_Cuts::iterator add_cut
    ( const typename PolCut<T>::TYPE type, const double b,
      const unsigned n, const PolVar<T>*X, const double*a,
      const PolVar<T>&X1, const double a1,
      const PolVar<T>&X2, const double a2,
      const PolVar<T>&X3, const double a3 )
    { return _add_cut( 0, type, b, n, X, a, X1, a1, X2, a2, X3, a3 ); }
  //! @brief Append new relaxation cut w/ selection amongst <a>n+1</a> variables
  typename t_Cuts::iterator add_cut
    ( const typename PolCut<T>::TYPE type, const double b,
      const std::set<unsigned>&ndx, const PolVar<T>*X, const double*a,
      const PolVar<T>&X1, const double a1 )
    { return _add_cut( 0, type, b, ndx, X, a, X1, a1 ); }
  //! @brief Append new relaxation cut w/ variable and coefficient maps
  template <typename U> typename t_Cuts::iterator add_cut
    ( const typename PolCut<T>::TYPE type, const double b,
      const std::map<U,PolVar<T>>&X, const std::map<U,double>&a,
      const PolVar<T>&X1, const double a1 )
    { return _add_cut( 0, type, b, X, a, X1, a1 ); }
  //! @brief Append new relaxation cut w/ variable and coefficient maps
  template <typename U> typename t_Cuts::iterator add_cut
    ( const typename PolCut<T>::TYPE type, const double b,
      const std::map<U,PolVar<T>>&X, const std::map<U,double>&a,
      const PolVar<T>&X1, const double a1, const PolVar<T>&X2, const double a2 )
    { return _add_cut( 0, type, b, X, a, X1, a1, X2, a2 ); }

  //! @brief Append new relaxation cuts for product term
  void append_cuts_TIMES
    ( PolVar<T> const& VarR, PolVar<T> const& Var1, PolVar<T> const& Var2,
      FFOp* pOp=nullptr );
  //! @brief Append new relaxation cuts for square term
  void append_cuts_SQR
    ( PolVar<T> const& VarR, PolVar<T> const& Var1, FFOp* pOp=nullptr );

  //! @brief Erase cut with iterator <a>itcut</a> from set of cuts
  void erase_cut
    ( typename t_Cuts::iterator itcut )
    { return _Cuts.erase( itcut ); }
  //! @brief Erase all cuts and auxiliary variables as well as main variable subdivision
  void reset_cuts
    ()
    { _reset_vars(); _erase_aux(); _erase_cuts(); }
  //! @brief Erase all cuts
  void erase_cuts
    ()
    { _erase_cuts(); }
};

template <typename T>
inline PolVar<T>*
PolImg<T>::_append_var
( FFVar const* var, T const& range, const bool cont )
{
  auto itVar = _Vars.find( var );
  if( itVar != _Vars.end() ){
    itVar->second->_update( range );
    itVar->second->_update( cont );
    return itVar->second;
  }
  PolVar<T>* pVar = new PolVar<T>;
  pVar->_set( this, *var, range, cont, _Vars.size() );
  _Vars.insert( std::make_pair( var, pVar ) );
  return pVar;
}

template <typename T>
inline PolVar<T>*
PolImg<T>::_update_var
( FFVar const* var, const bool cont )
{
  auto itVar = _Vars.find( var );
  if( itVar != _Vars.end() ){
    itVar->second->_update( cont );
    return itVar->second;
  }
  return nullptr;
}

template <typename T>
inline PolVar<T>*
PolImg<T>::_update_var
( FFVar const* var, T const& range )
{
  auto itVar = _Vars.find( var );
  if( itVar != _Vars.end() ){
    itVar->second->_update( range );
    return itVar->second;
  }
  return nullptr;
}

template <typename T> inline void
PolImg<T>::_erase_vars
()
{
  for( auto itv = _Vars.begin(); itv != _Vars.end(); ++itv )
    delete itv->second;
  _Vars.clear();
}

template <typename T> inline void
PolImg<T>::_reset_vars
()
{
  for( auto itv = _Vars.begin(); itv != _Vars.end(); ++itv ){
    itv->second->reset_subdiv();
    itv->second->reset_cuts();
  }
}

template <typename T>
inline PolVar<T>*
PolImg<T>::_append_aux
( const T&range, const bool cont )
{
  PolVar<T>* pAux = new PolVar<T>;
  pAux->_set( this, range, cont, _Aux.size() );
  _Aux.push_back( pAux );
  return pAux;
}

template <typename T> inline void
PolImg<T>::_erase_aux
()
{
  for( auto itv = _Aux.begin(); itv != _Aux.end(); ++itv )
    delete *itv;
  _Aux.clear();
}

template <typename T>
inline PolLQExpr<T>*
PolImg<T>::_append_LQ
( const PolVar<T>*varR )
{
  PolLQExpr<T>* pLQ = new PolLQExpr<T>( varR );
  _LQExpr.push_back( pLQ );
  return pLQ;
}

template <typename T>
inline void
PolImg<T>::_erase_LQ
()
{
  for( auto itl = _LQExpr.begin(); itl != _LQExpr.end(); ++itl )
    delete *itl;
  _LQExpr.clear();
}

template <typename T> inline void
PolImg<T>::_erase_cuts
()
{
  for( auto itc = _Cuts.begin(); itc != _Cuts.end(); ++itc )
    delete *itc;
  _Cuts.clear();
}
 
template <typename T>
inline void
PolImg<T>::_erase_cuts
( FFOp* op )
{
  auto itc = _Cuts.begin();
  while( itc != _Cuts.end() ){
    auto itp = itc; ++itc;
    if( (*itp)->op() == op ){ delete *itp; _Cuts.erase( itp ); }
  } 
}

template <typename T>
inline typename PolImg<T>::t_Cuts::iterator
PolImg<T>::_add_cut
( FFOp*op, const PolVar<T>&X1, const PolVar<T>&X2 )
{
  PolCut<T>* pCut = new PolCut<T>( op, X1, X2 );
  return _Cuts.insert( pCut );
}

template <typename T>
inline typename PolImg<T>::t_Cuts::iterator
PolImg<T>::_add_cut
( FFOp*op, const PolVar<T>&X1, const PolVar<T>&X2, const PolVar<T>&X3 )
{
  PolCut<T>* pCut = new PolCut<T>( op, X1, X2, X3 );
  return _Cuts.insert( pCut );
}

template <typename T>
inline typename PolImg<T>::t_Cuts::iterator
PolImg<T>::_add_cut
( FFOp*op, const PolVar<T>&X1, const PolVar<T>&X2, const double X3 )
{
  PolCut<T>* pCut = new PolCut<T>( op, X1, X2, X3 );
  return _Cuts.insert( pCut );
}

template <typename T>
inline typename PolImg<T>::t_Cuts::iterator
PolImg<T>::_add_cut
( FFOp*op, const typename PolCut<T>::TYPE type,
  const double b )
{
  PolCut<T>* pCut = new PolCut<T>( op, type, b );
  return _Cuts.insert( pCut );
}

template <typename T>
inline typename PolImg<T>::t_Cuts::iterator
PolImg<T>::_add_cut
( FFOp*op, const typename PolCut<T>::TYPE type,
  const double b, const PolVar<T>&X1, const double a1 )
{
  PolCut<T>* pCut = new PolCut<T>( op, type, b, X1, a1 );
  return _Cuts.insert( pCut );
}

template <typename T>
inline typename PolImg<T>::t_Cuts::iterator
PolImg<T>::_add_cut
( FFOp*op, const typename PolCut<T>::TYPE type,
  const double b, const PolVar<T>&X1, const double a1,
  const PolVar<T>&X2, const double a2 )
{
  PolCut<T>* pCut = new PolCut<T>( op, type, b, X1, a1, X2, a2 );
  return _Cuts.insert( pCut );
}

template <typename T>
inline typename PolImg<T>::t_Cuts::iterator
PolImg<T>::_add_cut
( FFOp*op, const typename PolCut<T>::TYPE type,
  const double b, const PolVar<T>&X1, const double a1,
  const PolVar<T>&X2, const double a2, const PolVar<T>&X3,
  const double a3 )
{
  PolCut<T>* pCut = new PolCut<T>( op, type, b, X1, a1, X2, a2, X3, a3 );
  return _Cuts.insert( pCut );
}

template <typename T>
inline typename PolImg<T>::t_Cuts::iterator
PolImg<T>::_add_cut
( FFOp*op, const typename PolCut<T>::TYPE type,
  const double b, const PolVar<T>&X1, const double a1,
  const PolVar<T>&X2, const double a2, const PolVar<T>&X3,
  const double a3, const PolVar<T>&X4, const double a4 )
{
  PolCut<T>* pCut = new PolCut<T>( op, type, b, X1, a1, X2, a2, X3, a3, X4, a4 );
  return _Cuts.insert( pCut );
}

template <typename T>
inline typename PolImg<T>::t_Cuts::iterator
PolImg<T>::_add_cut
( FFOp*op, const typename PolCut<T>::TYPE type,
  const double b, const unsigned n,
  const PolVar<T>*X, const double*a )
{
  //if( !n ) throw Exceptions( Exceptions::INTERNAL );
  PolCut<T>* pCut = new PolCut<T>( op, type, b, n, X, a );
  return _Cuts.insert( pCut );
}

template <typename T>
inline typename PolImg<T>::t_Cuts::iterator
PolImg<T>::_add_cut
( FFOp*op, const typename PolCut<T>::TYPE type,
  const double b, const std::set<unsigned>&ndx,
  const PolVar<T>*X, const double*a )
{
  PolCut<T>* pCut = new PolCut<T>( op, type, b, ndx, X, a );
  return _Cuts.insert( pCut );
}

template <typename T> template <typename U>
inline typename PolImg<T>::t_Cuts::iterator
PolImg<T>::_add_cut
( FFOp*op, const typename PolCut<T>::TYPE type,
  const double b, const std::map<U,PolVar<T>>&X,
  const std::map<U,double>&a )
{
  PolCut<T>* pCut = new PolCut<T>( op, type, b, X, a );
  return _Cuts.insert( pCut );
}

template <typename T>
inline typename PolImg<T>::t_Cuts::iterator
PolImg<T>::_add_cut
( FFOp*op, const typename PolCut<T>::TYPE type,
  const double b, const unsigned n,
  const PolVar<T>*X, const double*a,
  const PolVar<T>&X1, const double a1 )
{
  //if( !n ) throw Exceptions( Exceptions::INTERNAL );
  PolCut<T>* pCut = new PolCut<T>( op, type, b, n, X, a, X1, a1 );
  return _Cuts.insert( pCut );
}

template <typename T>
inline typename PolImg<T>::t_Cuts::iterator
PolImg<T>::_add_cut
( FFOp*op, const typename PolCut<T>::TYPE type,
  const double b, const unsigned n,
  const PolVar<T>*X, const double*a,
  const PolVar<T>&X1, const double a1,
  const PolVar<T>&X2, const double a2 )
{
  //if( !n ) throw Exceptions( Exceptions::INTERNAL );
  PolCut<T>* pCut = new PolCut<T>( op, type, b, n, X, a, X1, a1, X2, a2 );
  return _Cuts.insert( pCut );
}

template <typename T>
inline typename PolImg<T>::t_Cuts::iterator
PolImg<T>::_add_cut
( FFOp*op, const typename PolCut<T>::TYPE type,
  const double b, const unsigned n,
  const PolVar<T>*X, const double*a,
  const PolVar<T>&X1, const double a1,
  const PolVar<T>&X2, const double a2,
  const PolVar<T>&X3, const double a3 )
{
  //if( !n ) throw Exceptions( Exceptions::INTERNAL );
  PolCut<T>* pCut = new PolCut<T>( op, type, b, n, X, a, X1, a1, X2, a2, X3, a3 );
  return _Cuts.insert( pCut );
}

template <typename T>
inline typename PolImg<T>::t_Cuts::iterator
PolImg<T>::_add_cut
( FFOp*op, const typename PolCut<T>::TYPE type,
  const double b, const std::set<unsigned>&ndx,
  const PolVar<T>*X, const double*a,
  const PolVar<T>&X1, const double a1 )
{
  PolCut<T>* pCut = new PolCut<T>( op, type, b, ndx, X, a, X1, a1 );
  return _Cuts.insert( pCut );
}

template <typename T> template <typename U>
inline typename PolImg<T>::t_Cuts::iterator
PolImg<T>::_add_cut
( FFOp*op, const typename PolCut<T>::TYPE type,
  const double b, const std::map<U,PolVar<T>>&X,
  const std::map<U,double>&a, const PolVar<T>&X1,
  const double a1 )
{
  PolCut<T>* pCut = new PolCut<T>( op, type, b, X, a, X1, a1 );
  return _Cuts.insert( pCut );
}

template <typename T> template <typename U>
inline typename PolImg<T>::t_Cuts::iterator
PolImg<T>::_add_cut
( FFOp*op, const typename PolCut<T>::TYPE type,
  const double b, const std::map<U,PolVar<T>>&X,
  const std::map<U,double>&a, const PolVar<T>&X1,
  const double a1, const PolVar<T>&X2, const double a2 )
{
  PolCut<T>* pCut = new PolCut<T>( op, type, b, X, a, X1, a1, X2, a2 );
  return _Cuts.insert( pCut );
}

template <typename T> inline void
PolImg<T>::generate_cuts
( const unsigned ndep, const PolVar<T>*pdep, const bool reset )
{
  // Reset cuts in polyhedral image?
  if( reset ) reset_cuts();

  // Propagate cuts through all dependent subtrees
  for( unsigned i=0; i<ndep; i++ ){
    auto itDep = _Vars.find( const_cast<FFVar*>(&pdep[i]._var) );
    if( itDep != _Vars.end() ) itDep->second->generate_cuts();
  }
}

template <typename T> inline void
PolImg<T>::generate_cuts
( const std::set<unsigned>&ndxdep, const PolVar<T>*pdep, const bool reset )
{
  // Reset cuts in polyhedral image?
  if( reset ) reset_cuts();

  // Propagate cuts through all dependent subtrees
  auto it = ndxdep.cbegin();
  for( ; it != ndxdep.cend(); ++it ){
    auto itDep = _Vars.find( const_cast<FFVar*>(&pdep[*it]._var) );
    if( itDep != _Vars.end() ) itDep->second->generate_cuts();
  }
}

template <typename T> template <typename U>
inline void
PolImg<T>::generate_cuts
( const std::map<U,PolVar<T>>&mdep, const bool reset )
{
  // Reset cuts in polyhedral image?
  if( reset ) reset_cuts();

  // Propagate cuts through all dependent subtrees
  auto it = mdep.cbegin();
  for( ; it != mdep.cend(); ++it ){
    auto itDep = _Vars.find( const_cast<FFVar*>(&it->second._var) );
    if( itDep != _Vars.end() ) itDep->second->generate_cuts();
  }
}

template <typename T>
inline std::ostream&
operator <<
( std::ostream&out, const PolImg<T>&P )
{
  out << ( P._Vars.empty()? "\nNO VARIABLE\n": "\nVARIABLES:\n" );
  for( auto itv=P._Vars.begin(); itv!=P._Vars.end(); ++itv ){
    out << "  " << itv->second->name() << "\t in " << itv->second->range()
        << "\t (DAG: " << itv->second->var() << ")";
    //if( (*itv)->_Op ) out << "\t:= " << *(*itv)->_Op;
    out << std::endl;
  }
  out << ( P._Aux.empty()? "\nNO AUXILIARY\n": "\nAUXILIARIES:\n" );
  for( auto itv=P._Aux.begin(); itv!=P._Aux.end(); ++itv ){
    out << "  " << (*itv)->name() << "\t in " << (*itv)->range();
    //if( (*itv)->_Op ) out << "\t" << *(*itv)->_Op;
    out << std::endl;
  }

  out << ( P._Cuts.empty()? "\nNO CUT\n": "\nCUTS:\n" );
  for( auto itc=P._Cuts.begin(); itc!=P._Cuts.end(); ++itc )
    out << " " << **itc << std::endl;
  return out;
}

template <typename T>
inline PolVar<T>
operator^
( const PolVar<T>&Var1, const PolVar<T>&Var2 )
{
  if( Var1._img && Var2._img && Var1._img != Var2._img )
    throw typename PolImg<T>::Exceptions( PolImg<T>::Exceptions::ENVMIS );
  PolImg<T>* img = Var1._img? Var1._img: Var2._img;
  FFGraph* dag = Var1._var.cst()? Var2._var.dag(): Var1._var.dag();
  if( !dag || !dag->curOp() )
    throw typename PolImg<T>::Exceptions( PolImg<T>::Exceptions::NOTALLOWED );
  FFVar* pFFVarR = dag->curOp()->pres;
  T IVarR;
  if( !Op<T>::inter( IVarR, Var1._range, Var2._range ) )
    throw typename PolImg<T>::Exceptions( PolImg<T>::Exceptions::INTER );
  PolVar<T>* pVarR = img->_append_var( pFFVarR, IVarR, true );
  return *pVarR;
}

template <typename T> inline void
PolImg<T>::_add_cuts_INTER
( const PolVar<T>*VarR, FFVar*pVar1, FFVar*pVar2 )
{
  if( pVar1->cst() )
    _add_cut( VarR->_var.ops().first, PolCut<T>::EQ, pVar1->num().val(), *VarR, 1. );
  else{
    auto itVar1 = _Vars.find( pVar1 );
    _add_cut( VarR->_var.ops().first, PolCut<T>::EQ, 0., *VarR, 1., *itVar1->second, -1. );
  }
  if( pVar2->cst() )
    _add_cut( VarR->_var.ops().first, PolCut<T>::EQ, pVar2->num().val(), *VarR, 1. );
  else{
    auto itVar2 = _Vars.find( pVar2 );
    _add_cut( VarR->_var.ops().first, PolCut<T>::EQ, 0., *VarR, 1., *itVar2->second, -1. );
  }
}

template <typename T>
inline PolVar<T>
operator+
( const PolVar<T>&Var1, const double Cst2 )
{
  PolImg<T>* img = Var1._img;
  FFGraph* dag = Var1._var.dag();
  if( !dag || !dag->curOp() )
    throw typename PolImg<T>::Exceptions( PolImg<T>::Exceptions::NOTALLOWED );
  FFVar* pFFVarR = dag->curOp()->pres;
  PolVar<T>* pVarR = img->_append_var( pFFVarR, Var1.range() + Cst2, true );
  return *pVarR;
}

template <typename T>
inline PolVar<T>
operator+
( const PolVar<T>&Var1, const PolVar<T>&Var2 )
{
  if( Var1._img && Var2._img && Var1._img != Var2._img )
    throw typename PolImg<T>::Exceptions( PolImg<T>::Exceptions::ENVMIS );
  PolImg<T>* img = Var1._img? Var1._img: Var2._img;
  FFGraph* dag = Var1._var.cst()? Var2._var.dag(): Var1._var.dag();
  if( !dag || !dag->curOp() )
    throw typename PolImg<T>::Exceptions( PolImg<T>::Exceptions::NOTALLOWED );
  FFVar* pFFVarR = dag->curOp()->pres;
  PolVar<T>* pVarR = img->_append_var( pFFVarR, Var1.range() + Var2.range(), true );
  return *pVarR;
}

template <typename T> inline void
PolImg<T>::_add_cuts_PLUS
( const PolVar<T>*VarR, FFVar*pVar1, FFVar*pVar2 )
{
  if( pVar1->cst() && pVar2->cst() )
    _add_cut( VarR->_var.ops().first, PolCut<T>::EQ, pVar1->num().val()+pVar2->num().val(), *VarR, 1. );
  else if( pVar2->cst() ){
    auto itVar1 = _Vars.find( pVar1 );
    _add_cut( VarR->_var.ops().first, PolCut<T>::EQ, pVar2->num().val(), *VarR, 1., *itVar1->second, -1. );
  }
  else if( pVar1->cst() ){
    auto itVar2 = _Vars.find( pVar2 );
    _add_cut( VarR->_var.ops().first, PolCut<T>::EQ, pVar1->num().val(), *VarR, 1., *itVar2->second, -1. );
  }
  else{
    auto itVar1 = _Vars.find( pVar1 );
    auto itVar2 = _Vars.find( pVar2 );
    _add_cut( VarR->_var.ops().first, PolCut<T>::EQ, 0., *VarR, 1., *itVar1->second, -1., *itVar2->second, -1. );
  }
}

template <typename T> inline bool
PolImg<T>::_add_LQ_PLUS
( PolLQExpr<T>*&pLQ, const PolVar<T>*VarR, FFVar*pVar1, FFVar*pVar2 )
{
  if( !pLQ ) pLQ = _append_LQ( VarR );
  if( pVar1->cst() && pVar2->cst() )
    pLQ->substitute( VarR, pVar1->num().val()+pVar2->num().val() );
  else if( pVar2->cst() ){
    auto itVar1 = _Vars.find( pVar1 );
    pLQ->substitute( VarR, 1., itVar1->second, pVar2->num().val() );
  }
  else if( pVar1->cst() ){
    auto itVar2 = _Vars.find( pVar2 );
    pLQ->substitute( VarR, 1., itVar2->second, pVar1->num().val() );
  }
  else{
    auto itVar1 = _Vars.find( pVar1 );
    auto itVar2 = _Vars.find( pVar2 );
    pLQ->substitute( VarR, 1., itVar1->second, 1., itVar2->second );
  }
  return true;
}

template <typename T>
inline PolVar<T>
operator-
( const PolVar<T>&Var1 )
{
  FFGraph* dag = Var1._var.dag();
  if( !dag || !dag->curOp() )
    throw typename PolImg<T>::Exceptions( PolImg<T>::Exceptions::NOTALLOWED );
  PolImg<T>* img = Var1._img;
  FFVar* pFFVarR = dag->curOp()->pres;
  PolVar<T>* pVarR = img->_append_var( pFFVarR, -Var1.range(), true );
  return *pVarR;
}

template <typename T> inline void
PolImg<T>::_add_cuts_NEG
( const PolVar<T>*VarR, FFVar*pVar1 )
{
  if( pVar1->cst() )
    _add_cut( VarR->_var.ops().first, PolCut<T>::EQ, -pVar1->num().val(), *VarR, 1. );
  else{
    auto itVar1 = _Vars.find( pVar1 );
    _add_cut( VarR->_var.ops().first, PolCut<T>::EQ, 0., *VarR, 1., *itVar1->second, 1. );
  }
}

template <typename T> inline bool
PolImg<T>::_add_LQ_NEG
( PolLQExpr<T>*&pLQ, const PolVar<T>*VarR, FFVar*pVar1 )
{
  if( !pLQ ) pLQ = _append_LQ( VarR );
  if( pVar1->cst() )
    pLQ->substitute( VarR, -pVar1->num().val() );
  else{
    auto itVar1 = _Vars.find( pVar1 );
    pLQ->substitute( VarR, -1., itVar1->second );
  }
  return true;
}

template <typename T>
inline PolVar<T>
operator-
( const PolVar<T>&Var1, const PolVar<T>&Var2 )
{
  if( Var1._img && Var2._img && Var1._img != Var2._img )
    throw typename PolImg<T>::Exceptions( PolImg<T>::Exceptions::ENVMIS );
  PolImg<T>* img = Var1._img? Var1._img: Var2._img;
  FFGraph* dag = Var1._var.cst()? Var2._var.dag(): Var1._var.dag();
  if( !dag || !dag->curOp() )
    throw typename PolImg<T>::Exceptions( PolImg<T>::Exceptions::NOTALLOWED );
  FFVar* pFFVarR = dag->curOp()->pres;
  PolVar<T>* pVarR = img->_append_var( pFFVarR, Var1.range() - Var2.range(), true );
  return *pVarR;
}

template <typename T> inline void
PolImg<T>::_add_cuts_MINUS
( const PolVar<T>*VarR, FFVar*pVar1, FFVar*pVar2 )
{
  if( pVar1->cst() && pVar2->cst() )
    _add_cut( VarR->_var.ops().first, PolCut<T>::EQ, pVar1->num().val()-pVar2->num().val(), *VarR, 1. );
  else if( pVar2->cst() ){
    auto itVar1 = _Vars.find( pVar1 );
    _add_cut( VarR->_var.ops().first, PolCut<T>::EQ, -pVar2->num().val(), *VarR, 1., *itVar1->second, -1. );
  }
  else if( pVar1->cst() ){
    auto itVar2 = _Vars.find( pVar2 );
    _add_cut( VarR->_var.ops().first, PolCut<T>::EQ, pVar1->num().val(), *VarR, 1., *itVar2->second, 1. );
  }
  else{
    auto itVar1 = _Vars.find( pVar1 );
    auto itVar2 = _Vars.find( pVar2 );
    _add_cut( VarR->_var.ops().first, PolCut<T>::EQ, 0., *VarR, 1., *itVar1->second, -1., *itVar2->second, 1. );
  }
}

template <typename T> inline bool
PolImg<T>::_add_LQ_MINUS
( PolLQExpr<T>*&pLQ, const PolVar<T>*VarR, FFVar*pVar1, FFVar*pVar2 )
{
  if( !pLQ ) pLQ = _append_LQ( VarR );
  if( pVar1->cst() && pVar2->cst() )
    pLQ->substitute( VarR, pVar1->num().val()-pVar2->num().val() );
  else if( pVar2->cst() ){
    auto itVar1 = _Vars.find( pVar1 );
    pLQ->substitute( VarR, 1., itVar1->second, -pVar2->num().val() );
  }
  else if( pVar1->cst() ){
    auto itVar2 = _Vars.find( pVar2 );
    pLQ->substitute( VarR, -1., itVar2->second, pVar1->num().val() );
  }
  else{
    auto itVar1 = _Vars.find( pVar1 );
    auto itVar2 = _Vars.find( pVar2 );
    pLQ->substitute( VarR, 1., itVar1->second, -1., itVar2->second );
  }
  return true;
}

template <typename T>
inline PolVar<T>
operator*
( const PolVar<T>&Var1, const double Cst2 )
{
  PolImg<T>* img = Var1._img;
  FFGraph* dag = Var1._var.dag();
  if( !dag || !dag->curOp() )
    throw typename PolImg<T>::Exceptions( PolImg<T>::Exceptions::NOTALLOWED );
  FFVar* pFFVarR = dag->curOp()->pres;
  PolVar<T>* pVarR = img->_append_var( pFFVarR, Var1.range() * Cst2, true );
  return *pVarR;
}

template <typename T>
inline PolVar<T>
operator*
( const PolVar<T>&Var1, const PolVar<T>&Var2 )
{
  if( Var1._img && Var2._img && Var1._img != Var2._img )
    throw typename PolImg<T>::Exceptions( PolImg<T>::Exceptions::ENVMIS );
  PolImg<T>* img = Var1._img? Var1._img: Var2._img;
  FFGraph* dag = Var1._var.cst()? Var2._var.dag(): Var1._var.dag();
  if( !dag || !dag->curOp() )
    throw typename PolImg<T>::Exceptions( PolImg<T>::Exceptions::NOTALLOWED );
  FFVar* pFFVarR = dag->curOp()->pres;
  PolVar<T>* pVarR = img->_append_var( pFFVarR, Var1.range() * Var2.range(), Var1.discr() && Var2.discr()? false: true );
  return *pVarR;
}

template <typename T> inline void
PolImg<T>::_add_cuts_TIMES
( const PolVar<T>*VarR, FFVar*pVar1, FFVar*pVar2, FFOp*pOp )
{
  assert( VarR && pVar1 && pVar2 );
  if( pVar1->cst() && pVar2->cst() )
    _add_cut( pOp, PolCut<T>::EQ, pVar1->num().val() * pVar2->num().val(), *VarR, 1. );

  else if( pVar2->cst() ){
    auto itVar1 = _Vars.find( pVar1 );
    assert( itVar1 != _Vars.end() );
    _add_cut( pOp, PolCut<T>::EQ, 0., *VarR, 1., *itVar1->second, -pVar2->num().val() );
  }

  else if( pVar1->cst() ){
    auto itVar2 = _Vars.find( pVar2 );
    assert( itVar2 != _Vars.end() );
    _add_cut( pOp, PolCut<T>::EQ, 0., *VarR, 1., *itVar2->second, -pVar1->num().val() );
  }

  else if( !options.RELAX_QUAD ){
    auto itVar1 = _Vars.find( pVar1 );
    auto itVar2 = _Vars.find( pVar2 );
    assert( itVar1 != _Vars.end() && itVar2 != _Vars.end() );
    auto itCut = _add_cut( pOp, PolCut<T>::EQ, 0., *VarR, -1. );
    (*itCut)->append( *itVar1->second, *itVar2->second, 1. );
  }
  
  else{
    auto itVar1 = _Vars.find( pVar1 );
    auto itVar2 = _Vars.find( pVar2 );
    assert( itVar1 != _Vars.end() && itVar2 != _Vars.end() );
#ifndef MC__POLIMG_PWMCCORMICK_1D
    if( options.BREAKPOINT_TYPE == Options::NONE 
     || !_pwmccormick_cuts( pOp, *itVar1->second, Op<T>::l(itVar1->second->_range),
                            Op<T>::u(itVar1->second->_range), *itVar2->second, Op<T>::l(itVar2->second->_range),
                            Op<T>::u(itVar2->second->_range), *VarR ) ){
      _add_cut( pOp, PolCut<T>::GE, -Op<T>::u(itVar1->second->_range)*Op<T>::u(itVar2->second->_range),
        *VarR, 1., *itVar1->second, -Op<T>::u(itVar2->second->_range), *itVar2->second, -Op<T>::u(itVar1->second->_range) );
      _add_cut( pOp, PolCut<T>::GE, -Op<T>::l(itVar1->second->_range)*Op<T>::l(itVar2->second->_range),
        *VarR, 1., *itVar1->second, -Op<T>::l(itVar2->second->_range), *itVar2->second, -Op<T>::l(itVar1->second->_range) );
      _add_cut( pOp, PolCut<T>::LE, -Op<T>::u(itVar1->second->_range)*Op<T>::l(itVar2->second->_range),
        *VarR, 1., *itVar1->second, -Op<T>::l(itVar2->second->_range), *itVar2->second, -Op<T>::u(itVar1->second->_range) );
      _add_cut( pOp, PolCut<T>::LE, -Op<T>::l(itVar1->second->_range)*Op<T>::u(itVar2->second->_range),
        *VarR, 1., *itVar1->second, -Op<T>::u(itVar2->second->_range), *itVar2->second, -Op<T>::l(itVar1->second->_range) );
    }
#else
    switch( options.BREAKPOINT_TYPE ){
      case PolImg<T>::Options::BIN:
      case PolImg<T>::Options::SOS2:{
        const unsigned NKNOTS1 = itVar1->second->create_subdiv( Op<T>::l(itVar1->second->_range), Op<T>::u(itVar1->second->_range) ).size();
        if( NKNOTS1 > 2 )
          _pwmccormick_cuts( pOp, *itVar1->second, Op<T>::l(itVar1->second->_range),
            Op<T>::u(itVar1->second->_range), *itVar2->second, Op<T>::l(itVar2->second->_range),
            Op<T>::u(itVar2->second->_range), *VarR );
        const unsigned NKNOTS2 = itVar2->second->create_subdiv( Op<T>::l(itVar2->second->_range), Op<T>::u(itVar2->second->_range) ).size();
        if( NKNOTS2 > 2 )
          _pwmccormick_cuts( pOp, *itVar2->second, Op<T>::l(itVar2->second->_range),
            Op<T>::u(itVar2->second->_range), *itVar1->second, Op<T>::l(itVar1->second->_range),
            Op<T>::u(itVar1->second->_range), *VarR );
        if( NKNOTS1 > 2 || NKNOTS2 > 2 ) break; // The standard McCormick cuts are implied
      }
      case PolImg<T>::Options::NONE: default:
        _add_cut( pOp, PolCut<T>::GE, -Op<T>::u(itVar1->second->_range)*Op<T>::u(itVar2->second->_range),
          *VarR, 1., *itVar1->second, -Op<T>::u(itVar2->second->_range), *itVar2->second, -Op<T>::u(itVar1->second->_range) );
        _add_cut( pOp, PolCut<T>::GE, -Op<T>::l(itVar1->second->_range)*Op<T>::l(itVar2->second->_range),
          *VarR, 1., *itVar1->second, -Op<T>::l(itVar2->second->_range), *itVar2->second, -Op<T>::l(itVar1->second->_range) );
        _add_cut( pOp, PolCut<T>::LE, -Op<T>::u(itVar1->second->_range)*Op<T>::l(itVar2->second->_range),
          *VarR, 1., *itVar1->second, -Op<T>::l(itVar2->second->_range), *itVar2->second, -Op<T>::u(itVar1->second->_range) );
        _add_cut( pOp, PolCut<T>::LE, -Op<T>::l(itVar1->second->_range)*Op<T>::u(itVar2->second->_range),
          *VarR, 1., *itVar1->second, -Op<T>::u(itVar2->second->_range), *itVar2->second, -Op<T>::l(itVar1->second->_range) );
    }
#endif
  }
}

template <typename T> inline void
PolImg<T>::_add_cuts_TIMES
( const PolVar<T>*VarR, FFVar*pVar1, FFVar*pVar2 )
{
  return _add_cuts_TIMES( VarR, pVar1, pVar2, VarR->_var.ops().first );
}

template <typename T> inline void
PolImg<T>::append_cuts_TIMES
( PolVar<T> const& VarR, PolVar<T> const& Var1, PolVar<T> const& Var2, FFOp* pOp )
{
  assert( &Var1._var && &Var2._var );
  return _add_cuts_TIMES( &VarR, const_cast<FFVar*>(&Var1._var), const_cast<FFVar*>(&Var2._var), pOp );

//  if( !options.RELAX_QUAD ){
//    auto itCut = _add_cut( pOp, PolCut<T>::EQ, 0., VarR, -1. );
//    (*itCut)->append( Var1, Var2, 1. );
//  }
//  
//  else{
//  
//#ifndef MC__POLIMG_PWMCCORMICK_1D
//    if( options.BREAKPOINT_TYPE == Options::NONE 
//     || !_pwmccormick_cuts( pOp, Var1, Op<T>::l(Var1._range), Op<T>::u(Var1._range),
//                            Var2, Op<T>::l(Var2._range), Op<T>::u(Var2._range), VarR ) ){
//      _add_cut( pOp, PolCut<T>::GE, -Op<T>::u(Var1._range)*Op<T>::u(Var2._range),
//                VarR, 1., Var1, -Op<T>::u(Var2._range), Var2, -Op<T>::u(Var1._range) );
//      _add_cut( pOp, PolCut<T>::GE, -Op<T>::l(Var1._range)*Op<T>::l(Var2._range),
//                VarR, 1., Var1, -Op<T>::l(Var2._range), Var2, -Op<T>::l(Var1._range) );
//      _add_cut( pOp, PolCut<T>::LE, -Op<T>::u(Var1._range)*Op<T>::l(Var2._range),
//                VarR, 1., Var1, -Op<T>::l(Var2._range), Var2, -Op<T>::u(Var1._range) );
//      _add_cut( pOp, PolCut<T>::LE, -Op<T>::l(Var1._range)*Op<T>::u(Var2._range),
//                VarR, 1., Var1, -Op<T>::u(Var2._range), Var2, -Op<T>::l(Var1._range) );
//    }
//#else
//    switch( options.BREAKPOINT_TYPE ){
//      case PolImg<T>::Options::BIN:
//      case PolImg<T>::Options::SOS2:{
//        const unsigned NKNOTS1 = Var1.create_subdiv( Op<T>::l(Var1._range), Op<T>::u(Var1._range) ).size();
//        if( NKNOTS1 > 2 )
//          _pwmccormick_cuts( pOp, Var1, Op<T>::l(Var1._range), Op<T>::u(Var1._range),
//                             Var2, Op<T>::l(Var2._range), Op<T>::u(Var2._range), VarR );
//        const unsigned NKNOTS2 = Var2.create_subdiv( Op<T>::l(Var2._range), Op<T>::u(Var2._range) ).size();
//        if( NKNOTS2 > 2 )
//          _pwmccormick_cuts( pOp, *Var2, Op<T>::l(Var2._range), Op<T>::u(Var2._range),
//                             Var1, Op<T>::l(Var1._range), Op<T>::u(Var1._range), VarR );
//        if( NKNOTS1 > 2 || NKNOTS2 > 2 ) break; // The standard McCormick cuts are implied
//      }
//      case PolImg<T>::Options::NONE: default:
//        _add_cut( pOp, PolCut<T>::GE, -Op<T>::u(Var1._range)*Op<T>::u(Var2._range),
//                  VarR, 1., Var1, -Op<T>::u(Var2._range), Var2, -Op<T>::u(Var1._range) );
//        _add_cut( pOp, PolCut<T>::GE, -Op<T>::l(Var1._range)*Op<T>::l(Var2._range),
//                  VarR, 1., Var1, -Op<T>::l(Var2._range), Var2, -Op<T>::l(Var1._range) );
//        _add_cut( pOp, PolCut<T>::LE, -Op<T>::u(Var1._range)*Op<T>::l(Var2._range),
//                  VarR, 1., Var1, -Op<T>::l(Var2._range), Var2, -Op<T>::u(Var1._range) );
//        _add_cut( pOp, PolCut<T>::LE, -Op<T>::l(Var1._range)*Op<T>::u(Var2._range),
//                  VarR, 1., Var1, -Op<T>::u(Var2._range), Var2, -Op<T>::l(Var1._range) );
//    }
//#endif
//  }
}

template <typename T> inline bool
PolImg<T>::_add_LQ_TIMES
( PolLQExpr<T>*&pLQ, const PolVar<T>*VarR, FFVar*pVar1, FFVar*pVar2 )
{
  if( options.RELAX_QUAD && !pVar1->cst() && !pVar2->cst() ) return false;
  if( !pLQ ) pLQ = _append_LQ( VarR );
  if( pVar1->cst() && pVar2->cst() )
    pLQ->substitute( VarR, pVar1->num().val()*pVar2->num().val() );
  else if( pVar2->cst() ){
    auto itVar1 = _Vars.find( pVar1 );
    pLQ->substitute( VarR, pVar2->num().val(), itVar1->second );
  }
  else if( pVar1->cst() ){
    auto itVar2 = _Vars.find( pVar2 );
    pLQ->substitute( VarR, pVar1->num().val(), itVar2->second );
  }
  else{
    auto itVar1 = _Vars.find( pVar1 );
    auto itVar2 = _Vars.find( pVar2 );
    pLQ->substitute( VarR, 1., itVar1->second, itVar2->second );
  }
  return true;
}

template <typename T>
inline PolVar<T>
operator/
( const double Cst1, const PolVar<T>&Var2 )
{
  PolImg<T>* img = Var2._img;
  FFGraph* dag = Var2._var.dag();
  if( !dag || !dag->curOp() )
    throw typename PolImg<T>::Exceptions( PolImg<T>::Exceptions::NOTALLOWED );
  FFVar* pFFVarR = dag->curOp()->pres;
  PolVar<T>* pVarR = img->_append_var( pFFVarR, Cst1 / Var2.range(), true );
  return *pVarR;
}

template <typename T>
inline PolVar<T>
operator/
( const PolVar<T>&Var1, const PolVar<T>&Var2 )
{
  if( Var1._img && Var2._img && Var1._img != Var2._img )
    throw typename PolImg<T>::Exceptions( PolImg<T>::Exceptions::ENVMIS );
  PolImg<T>* img = Var1._img? Var1._img: Var2._img;
  FFGraph* dag = Var1._var.cst()? Var2._var.dag(): Var1._var.dag();
  if( !dag || !dag->curOp() )
    throw typename PolImg<T>::Exceptions( PolImg<T>::Exceptions::NOTALLOWED );
  FFVar* pFFVarR = dag->curOp()->pres;
  PolVar<T>* pVarR = img->_append_var( pFFVarR, Var1.range() / Var2.range(), true );
  return *pVarR;
}

template <typename T> inline void
PolImg<T>::_add_cuts_DIV
( const PolVar<T>*VarR, FFVar*pVar1, FFVar*pVar2 )
{
  if( pVar1->cst() && pVar2->cst() )
    _add_cut( VarR->_var.ops().first, PolCut<T>::EQ, pVar1->num().val(), *VarR, pVar2->num().val() );

  else if( pVar2->cst() ){
    auto itVar1 = _Vars.find( pVar1 );
    _add_cut( VarR->_var.ops().first, PolCut<T>::EQ, 0., *VarR, pVar2->num().val(), *itVar1->second, -1. );
  }

  else if( pVar1->cst() ){
    auto itVar2 = _Vars.find( pVar2 );
    double Cst1 = pVar1->num().val();
    struct loc{
      static std::pair<double,double> scalinv
        ( const double x, const double*rusr, const int*iusr )
        { return std::make_pair( *rusr/x, -*rusr/(x*x) ); }
    };
    // -- Convex Case
    if( (pVar1->num().val() >= 0. && Op<T>::l(itVar2->second->_range) > 0.)
     || (pVar1->num().val() <= 0. && Op<T>::u(itVar2->second->_range) < 0.) ){
      _semilinear_cuts( VarR->_var.ops().first, *itVar2->second, Op<T>::l(itVar2->second->_range),
        Op<T>::u(itVar2->second->_range), *VarR, PolCut<T>::LE, loc::scalinv, &Cst1, 0 );
      _sandwich_cuts( VarR->_var.ops().first, *itVar2->second, Op<T>::l(itVar2->second->_range),
        Op<T>::u(itVar2->second->_range), *VarR, Op<T>::l(VarR->_range), Op<T>::u(VarR->_range),
        PolCut<T>::GE, loc::scalinv, &Cst1, 0 );
    }
    // -- Concave Case
    else{
      _semilinear_cuts( VarR->_var.ops().first, *itVar2->second, Op<T>::l(itVar2->second->_range),
        Op<T>::u(itVar2->second->_range), *VarR, PolCut<T>::GE, loc::scalinv, &Cst1, 0 );
      _sandwich_cuts( VarR->_var.ops().first, *itVar2->second, Op<T>::l(itVar2->second->_range),
        Op<T>::u(itVar2->second->_range), *VarR, Op<T>::l(VarR->_range), Op<T>::u(VarR->_range),
        PolCut<T>::LE, loc::scalinv, &Cst1, 0 );
    }
  }

  else{
    auto itVar1 = _Vars.find( pVar1 );
    auto itVar2 = _Vars.find( pVar2 );
#ifndef MC__POLIMG_PWMCCORMICK_1D
    if( options.BREAKPOINT_TYPE == Options::NONE 
     || !_pwmccormick_cuts( VarR->_var.ops().first, *VarR, Op<T>::l(VarR->_range), Op<T>::u(VarR->_range),
                            *itVar2->second, Op<T>::l(itVar2->second->_range), Op<T>::u(itVar2->second->_range),
			    *itVar1->second ) ){
      _add_cut( VarR->_var.ops().first, PolCut<T>::GE, -Op<T>::u(VarR->_range)*Op<T>::u(itVar2->second->_range),
        *itVar1->second, 1., *VarR, -Op<T>::u(itVar2->second->_range), *itVar2->second, -Op<T>::u(VarR->_range) );
      _add_cut( VarR->_var.ops().first, PolCut<T>::GE, -Op<T>::l(VarR->_range)*Op<T>::l(itVar2->second->_range),
        *itVar1->second, 1., *VarR, -Op<T>::l(itVar2->second->_range), *itVar2->second, -Op<T>::l(VarR->_range) );
      _add_cut( VarR->_var.ops().first, PolCut<T>::LE, -Op<T>::u(VarR->_range)*Op<T>::l(itVar2->second->_range),
        *itVar1->second, 1., *VarR, -Op<T>::l(itVar2->second->_range), *itVar2->second, -Op<T>::u(VarR->_range) );
      _add_cut( VarR->_var.ops().first, PolCut<T>::LE, -Op<T>::l(VarR->_range)*Op<T>::u(itVar2->second->_range),
        *itVar1->second, 1., *VarR, -Op<T>::u(itVar2->second->_range), *itVar2->second, -Op<T>::l(VarR->_range) );
    }
#else
    switch( options.BREAKPOINT_TYPE ){
      case PolImg<T>::Options::BIN:
      case PolImg<T>::Options::SOS2:{
        const unsigned NKNOTSR = VarR->create_subdiv( Op<T>::l(VarR->_range), Op<T>::u(VarR->_range) ).size();
        if( NKNOTSR > 2 )
          _pwmccormick_cuts( VarR->_var.ops().first, *VarR, Op<T>::l(VarR->_range), Op<T>::u(VarR->_range),
            *itVar2->second, Op<T>::l(itVar2->second->_range), Op<T>::u(itVar2->second->_range), *itVar1->second );
        const unsigned NKNOTS2 = itVar2->second->create_subdiv( Op<T>::l(itVar2->second->_range), Op<T>::u(itVar2->second->_range) ).size();
        if( NKNOTS2 > 2 )
          _pwmccormick_cuts( VarR->_var.ops().first, *itVar2->second, Op<T>::l(itVar2->second->_range),
            Op<T>::u(itVar2->second->_range), *VarR, Op<T>::l(VarR->_range), Op<T>::u(VarR->_range), *itVar1->second );
        if( NKNOTSR > 2 || NKNOTS2 > 2 ) break; // The standard McCormick cuts are implied
      }
      case PolImg<T>::Options::NONE: default:
        _add_cut( VarR->_var.ops().first, PolCut<T>::GE, -Op<T>::u(VarR->_range)*Op<T>::u(itVar2->second->_range),
          *itVar1->second, 1., *VarR, -Op<T>::u(itVar2->second->_range), *itVar2->second, -Op<T>::u(VarR->_range) );
        _add_cut( VarR->_var.ops().first, PolCut<T>::GE, -Op<T>::l(VarR->_range)*Op<T>::l(itVar2->second->_range),
          *itVar1->second, 1., *VarR, -Op<T>::l(itVar2->second->_range), *itVar2->second, -Op<T>::l(VarR->_range) );
        _add_cut( VarR->_var.ops().first, PolCut<T>::LE, -Op<T>::u(VarR->_range)*Op<T>::l(itVar2->second->_range),
          *itVar1->second, 1., *VarR, -Op<T>::l(itVar2->second->_range), *itVar2->second, -Op<T>::u(VarR->_range) );
        _add_cut( VarR->_var.ops().first, PolCut<T>::LE, -Op<T>::l(VarR->_range)*Op<T>::u(itVar2->second->_range),
          *itVar1->second, 1., *VarR, -Op<T>::u(itVar2->second->_range), *itVar2->second, -Op<T>::l(VarR->_range) );
    }
#endif
  }
}

template <typename T> inline bool
PolImg<T>::_add_LQ_DIV
( PolLQExpr<T>*&pLQ, const PolVar<T>*VarR, FFVar*pVar1, FFVar*pVar2 )
{
  if( !pVar2->cst() ) return false;
  if( pVar1->cst() )
    pLQ->substitute( VarR, pVar1->num().val()/pVar2->num().val() );
  else{
    auto itVar1 = _Vars.find( pVar1 );
    pLQ->substitute( VarR, 1./pVar2->num().val(), itVar1->second );
  }
  return true;
}

template <typename T> inline std::pair< double, double >
PolImg<T>::_distmax
( p_dUniv f, const double xL, const double xU, const double*rpar,
  const int*ipar ) const
{
  std::pair<double,double> fL = f(xL,rpar,ipar), fU = f(xU,rpar,ipar);
  double Ddf = fU.second-fL.second,
         Adf = std::max(std::fabs(fL.second),std::fabs(fU.second));

  double xmid = 0.5 * ( xL + xU );
  double xmax = ( std::fabs(Ddf) - Adf*options.FRACTIONAL_RTOL > options.FRACTIONAL_ATOL?
    ( fU.second * xU - fL.second * xL - fU.first + fL.first ) / Ddf: xmid );
  std::pair<double,double> fmax = f(xmax,rpar,ipar);
  double dmax = std::fabs( fmax.first - fL.second * ( xmax - xL ) - fL.first );

  switch( options.SANDWICH_RULE ){
  case Options::BISECT:
    return std::make_pair( xmid, dmax );
  case Options::MAXERR: default:
    return std::make_pair( xmax, dmax );
  }
}

template <typename T>
inline void
PolImg<T>::_sandwich_cuts
( FFOp*pOp, const PolVar<T>&X, const double XL, const double XU, const PolVar<T>&Y,
  const double YL, const double YU, const typename PolCut<T>::TYPE sense, p_dUniv f,
  const double*rpar, const int*ipar )
{
  t_OA OA;
  unsigned NBCUTS = 0;

  switch( options.BREAKPOINT_TYPE ){
    // OA cuts @xL,xU + @breakpoints
    case Options::BIN:
    case Options::SOS2:{
      const std::vector<double>& XKNOT = X.create_subdiv( XL, XU );
      const unsigned NKNOTS = XKNOT.size();
      if( NKNOTS > 2 ){
        for( ; NBCUTS<NKNOTS; NBCUTS++ ){
          _linearization_cut( pOp, XKNOT[NBCUTS], X, XL, XU, Y, YL, YU, sense, f, rpar, ipar );
          if( !NBCUTS ) continue;
          std::pair<double,double> worst = _distmax( f, XKNOT[NBCUTS-1], XKNOT[NBCUTS], rpar, ipar );
          OA.push( OAsub( XKNOT[NBCUTS-1], XKNOT[NBCUTS], worst.first, worst.second ) );
#ifdef MC__POLIMG_DEBUG_SANDWICH
          std::cerr << OA.top() << std::endl;
#endif
        }
        break;
      }
    }

    // OA cuts @xL,xU only
    case Options::NONE: default:{
      _linearization_cut( pOp, XL, X, XL, XU, Y, YL, YU, sense, f, rpar, ipar );
      if( mc::isequal( XL, XU ) ) return;
      _linearization_cut( pOp, XU, X, XL, XU, Y, YL, YU, sense, f, rpar, ipar );
      NBCUTS += 2;
      std::pair<double,double> worst = _distmax( f, XL, XU, rpar, ipar );
      OA.push( OAsub( XL, XU, worst.first, worst.second ) );
#ifdef MC__POLIMG_DEBUG_SANDWICH
      std::cerr << OA.top() << std::endl;
#endif
      break;
    }
  }

  const double dtol = options.SANDWICH_ATOL
    + options.SANDWICH_RTOL*std::max(std::fabs(YL),std::fabs(YU));
#ifdef MC__POLIMG_DEBUG_SANDWICH
  std::cerr << "DTOL: " << dtol << std::endl;
  std::cerr << OA.top() << std::endl;
#endif

  // OA cut @xM
  while( OA.top().gap() > dtol &&
         ( !options.SANDWICH_MAXCUT || NBCUTS < options.SANDWICH_MAXCUT ) ){
    // x - y/exp(xref) <= xref - 1, @xref=xmax
    _linearization_cut( pOp, OA.top().xM(), X, OA.top().xL(), OA.top().xU(),
      Y, YL, YU, sense, f, rpar, ipar );
    NBCUTS++;
    std::pair<double,double> worst = _distmax( f, OA.top().xL(), OA.top().xM(), rpar, ipar );
    OA.push( OAsub( OA.top().xL(), OA.top().xM(), worst.first, worst.second ) );
#ifdef MC__POLIMG_DEBUG_SANDWICH
    std::cerr << "  Pushed: "  << FPOuter( OA.top().xL(), OA.top().xM(),
      worst.first, worst.second ) << std::endl;
#endif
    worst = _distmax( f, OA.top().xM(), OA.top().xU(), rpar, ipar );
    OA.push( OAsub( OA.top().xM(), OA.top().xU(), worst.first, worst.second ) );
#ifdef MC__POLIMG_DEBUG_SANDWICH
    std::cerr << "  Pushed: " << FPOuter( OA.top().xM(), OA.top().xU(),
      worst.first, worst.second ) << std::endl;
#endif
    OA.pop();
#ifdef MC__POLIMG_DEBUG_SANDWICH
    std::cerr << OA.top() << std::endl;
#endif
  }
#ifdef MC__POLIMG_DEBUG_SANDWICH
  std::cerr << "NBCUTS: " << NBCUTS << std::endl;
#endif
}

template <typename T>
inline void
PolImg<T>::_linearization_cut
( FFOp*pOp, const double Xref, const PolVar<T>&X, const double XL, const double XU,
  const PolVar<T>&Y, const double YL, const double YU, const typename PolCut<T>::TYPE sense,
  p_dUniv f, const double*rpar, const int*ipar )
{
  // sense: f convex -> PolCut<T>::GE - f concave -> PolCut<T>::LE
  // Y - f'(Xref)·X >= f(Xref) - f'(Xref)·Xref
  const std::pair<double,double> Fref = f(Xref,rpar,ipar);
  _add_cut( pOp, sense, Fref.first-Xref*Fref.second, Y, 1., X, -Fref.second );
}

template <typename T> inline double
PolImg<T>::_newton
( const double x0, const double xL, const double xU, p_dUniv f,
  const double TOL, const unsigned MAXIT, const double*rusr, const int*iusr ) const
{
  double xk = std::max(xL,std::min(xU,x0)), dk;
  std::pair<double,double> fk = f(xk,rusr,iusr);
  
  for( unsigned int it=0; it<MAXIT; it++ ){
    if( std::fabs(fk.first) < TOL ) return xk;
    if( fk.second == 0 ) throw Exceptions( Exceptions::ROOT );
    dk = fk.first/fk.second;
    if( mc::isequal(xk,xL) && dk>0 ) return xk;
    if( mc::isequal(xk,xU) && dk<0 ) return xk;
    xk = std::max(xL,std::min(xU,xk-dk));
    fk = f(xk,rusr,iusr);
  }

  throw Exceptions( Exceptions::ROOT );
}

template <typename T> inline double
PolImg<T>::_secant
( const double x0, const double x1, const double xL, const double xU, p_Univ f,
  const double TOL, const unsigned MAXIT, const double*rusr, const int*iusr ) const
{
  double xkm = std::max(xL,std::min(xU,x0));
  double fkm = f(xkm,rusr,iusr);
  double xk = std::max(xL,std::min(xU,x1));
  
  for( unsigned int it=0; it<MAXIT; it++ ){
    double fk = f(xk,rusr,iusr);
    if( std::fabs(fk) < TOL ) return xk;
    double Bk = (fk-fkm)/(xk-xkm);
    if( Bk == 0 ) throw Exceptions( Exceptions::ROOT );
    if( isequal(xk,xL) && fk/Bk>0 ) return xk;
    if( isequal(xk,xU) && fk/Bk<0 ) return xk;
    xkm = xk;
    fkm = fk;
    xk = std::max(xL,std::min(xU,xk-fk/Bk));
  }

  throw Exceptions( Exceptions::ROOT );
}

template <typename T> inline double
PolImg<T>::_goldsect
( const double xL, const double xU, p_Univ f, const double TOL, const unsigned MAXIT,
  const double*rusr, const int*iusr ) const
{
  const double phi = 2.-(1.+std::sqrt(5.))/2.;
  const double fL = f(xL,rusr,iusr), fU = f(xU,rusr,iusr);
  if( fL*fU > 0 ) throw Exceptions( Exceptions::ROOT );
  const double xm = xU-phi*(xU-xL), fm = f(xm,rusr,iusr);
  return _goldsect_iter( true, xL, fL, xm, fm, xU, fU, f, TOL, MAXIT, rusr, iusr );
}

template <typename T> inline double
PolImg<T>::_goldsect_iter
( const bool init, const double a, const double fa, const double b,
  const double fb, const double c, const double fc, p_Univ f,
  const double TOL, const unsigned MAXIT, const double*rusr, const int*iusr ) const
// a and c are the current bounds; the minimum is between them.
// b is a center point
{
  static unsigned int iter;
  iter = ( init? 1: iter+1 );
  const double phi = 2.-(1.+std::sqrt(5.))/2.;
  bool b_then_x = ( c-b > b-a );
  double x = ( b_then_x? b+phi*(c-b): b-phi*(b-a) );
  if( std::fabs(c-a) < TOL*(std::fabs(b)+std::fabs(x)) 
   || iter > MAXIT ) return (c+a)/2.;
  double fx = f(x,rusr,iusr);
  if( b_then_x )
    return( fa*fx<0? _goldsect_iter( false, a, fa, b, fb, x, fx, f, TOL, MAXIT, rusr, iusr ):
                     _goldsect_iter( false, b, fb, x, fx, c, fc, f, TOL, MAXIT, rusr, iusr ) );
  return( fa*fb<0? _goldsect_iter( false, a, fa, x, fx, b, fb, f, TOL, MAXIT, rusr, iusr ):
                   _goldsect_iter( false, x, fx, b, fb, c, fc, f, TOL, MAXIT, rusr, iusr ) );
}

template <typename T>
inline void
PolImg<T>::_semilinear_cuts
( FFOp*pOp, const PolVar<T>&X, const double XL, const double XU, const PolVar<T>&Y,
  const typename PolCut<T>::TYPE sense, p_dUniv f, const double*rpar, const int*ipar )
{
  switch( options.BREAKPOINT_TYPE ){
   case Options::BIN:{
    const std::vector<double>& XKNOT = X.create_subdiv( XL, XU );
    const unsigned NKNOTS = XKNOT.size()-1;
    if( NKNOTS > 1 ){
      // Represent variable range using linear binary transformation
      const std::vector< PolVar<T> >& subvar = X.BIN_subdiv( pOp );
      // Append semilinear cuts
      double coef[NKNOTS];
      double rhs = f( XKNOT[0], rpar, ipar ).first;
      for( unsigned isub=0; isub<NKNOTS; isub++ )
        coef[isub] = f( XKNOT[isub], rpar, ipar ).first - f( XKNOT[isub+1], rpar, ipar ).first;
      _add_cut( pOp, sense, rhs, NKNOTS, subvar.data(), coef, Y, 1. );
    }
    break;
   }
   case Options::SOS2:{
    const std::vector<double>& XKNOT = X.create_subdiv( XL, XU );
    const unsigned NKNOTS = XKNOT.size();
    if( NKNOTS > 2 ){
      // Represent variable range using SOS2 transformation
      const std::vector< PolVar<T> >& subvar = X.SOS2_subdiv( pOp );
      // Append semilinear cuts
      double coef[NKNOTS];
      for( unsigned isub=0; isub<NKNOTS; isub++ )
        coef[isub] = -f( XKNOT[isub], rpar, ipar ).first;
      _add_cut( pOp, sense, 0., NKNOTS, subvar.data(), coef, Y, 1. );
    }
    break;
   }
   default:
    break;
  }
  double dX = XU-XL, YL = f(XL,rpar,ipar).first, dY = f(XU,rpar,ipar).first-YL;
  _add_cut( pOp, sense, dX*YL-dY*XL, Y, dX, X, -dY );
}

#ifdef MC__POLIMG_PWMCCORMICK_1D
template <typename T>
inline bool
PolImg<T>::_pwmccormick_cuts
( FFOp*pOp, const PolVar<T>&X1, const double X1L, const double X1U,
  const PolVar<T>&X2, const double X2L, const double X2U, const PolVar<T>&Y )
{
  switch( options.BREAKPOINT_TYPE ){
   case Options::BIN:{
    const std::vector<double>& XKNOT = X1.create_subdiv( X1L, X1U );
    const unsigned NINTS = XKNOT.size()-1;
    if( NINTS > 1 ){
      double coef[NINTS];
      PolVar<T> auxr[NINTS];
      for( unsigned k=0; k<NINTS; k++ ){
        coef[k] = XKNOT[k+1] - XKNOT[k];
	auxr[k]  = *_append_aux( Op<T>::zeroone()*(X2U-X2L), true );
      }
      // Represent variable subdivision using linear binary transformation
      const std::vector< PolVar<T> >& aux1 = X1.BIN_subdiv( pOp );
      // Cut (46) in Gounaris et al. (2009), DOI: 10.1021/ie8016048
      _add_cut( pOp, PolCut<T>::EQ, X1L*X2L, NINTS, auxr, coef, X1, X2L, X2, X1L, Y, -1. );
      // Cuts (50) in Gounaris et al. (2009), DOI: 10.1021/ie8016048     
      for( unsigned k=0; k<NINTS; k++ ){
        _add_cut( pOp, PolCut<T>::GE, 0.,  aux1[k], X2U-X2L, auxr[k], -1. );
        _add_cut( pOp, PolCut<T>::LE, X2U, aux1[k], X2U-X2L, auxr[k], -1., X2, 1. );
        if( k ) _add_cut( pOp, PolCut<T>::LE, 0., auxr[k], 1., auxr[k-1], -1. );
	else    _add_cut( pOp, PolCut<T>::GE, X2L, X2, 1., auxr[0], -1. );
      }
    }
    break;
   }
   case Options::SOS2:
   default:
    break;
  }
  return true;
}
#else
template <typename T>
inline bool
PolImg<T>::_pwmccormick_cuts
( FFOp*pOp, const PolVar<T>&X1, const double X1L, const double X1U,
  const PolVar<T>&X2, const double X2L, const double X2U, const PolVar<T>&Y )
{
  switch( options.BREAKPOINT_TYPE ){
   case Options::BIN:{
    const std::vector<double>& X1KNOT = X1.create_subdiv( X1L, X1U );
    const unsigned NX1INTS = X1KNOT.size()-1;
    const std::vector<double>& X2KNOT = X2.create_subdiv( X2L, X2U );
    const unsigned NX2INTS = X2KNOT.size()-1;

    if( NX1INTS > 1 && NX2INTS == 1 ){
      double coef[NX1INTS];
      PolVar<T> auxr[NX1INTS];
      for( unsigned k=0; k<NX1INTS; k++ ){
        coef[k] = X1KNOT[k+1] - X1KNOT[k];
	auxr[k]  = *_append_aux( Op<T>::zeroone()*(X2U-X2L), true );
      }
      // Represent variable subdivision using linear binary transformation
      const std::vector< PolVar<T> >& aux1 = X1.BIN_subdiv( pOp );
      // Cut (46) in Gounaris et al. (2009), DOI: 10.1021/ie8016048
      _add_cut( pOp, PolCut<T>::EQ, X1L*X2L, NX1INTS, auxr, coef, X1, X2L, X2, X1L, Y, -1. );
      // Cuts (50) in Gounaris et al. (2009), DOI: 10.1021/ie8016048     
      for( unsigned k=0; k<NX1INTS; k++ ){
        _add_cut( pOp, PolCut<T>::GE, 0.,  aux1[k], X2U-X2L, auxr[k], -1. );
        _add_cut( pOp, PolCut<T>::LE, X2U, aux1[k], X2U-X2L, auxr[k], -1., X2, 1. );
        if( k ) _add_cut( pOp, PolCut<T>::LE, 0., auxr[k], 1., auxr[k-1], -1. );
	else    _add_cut( pOp, PolCut<T>::GE, X2L, X2, 1., auxr[0], -1. );
      }
      break;
    }

    else if( NX1INTS == 1 && NX2INTS > 1 ){
      double coef[NX2INTS];
      PolVar<T> auxr[NX2INTS];
      for( unsigned k=0; k<NX2INTS; k++ ){
        coef[k] = X2KNOT[k+1] - X2KNOT[k];
	auxr[k]  = *_append_aux( Op<T>::zeroone()*(X1U-X1L), true );
      }
      // Represent variable subdivision using linear binary transformation
      const std::vector< PolVar<T> >& aux2 = X2.BIN_subdiv( pOp );
      // Cut (46) in Gounaris et al. (2009), DOI: 10.1021/ie8016048
      _add_cut( pOp, PolCut<T>::EQ, X1L*X2L, NX2INTS, auxr, coef, X1, X2L, X2, X1L, Y, -1. );
      // Cuts (50) in Gounaris et al. (2009), DOI: 10.1021/ie8016048     
      for( unsigned k=0; k<NX2INTS; k++ ){
        _add_cut( pOp, PolCut<T>::GE, 0.,  aux2[k], X1U-X1L, auxr[k], -1. );
        _add_cut( pOp, PolCut<T>::LE, X1U, aux2[k], X1U-X1L, auxr[k], -1., X1, 1. );
        if( k ) _add_cut( pOp, PolCut<T>::LE, 0., auxr[k], 1., auxr[k-1], -1. );
	else    _add_cut( pOp, PolCut<T>::GE, X1L, X1, 1., auxr[0], -1. );
      }
      break;
    }
    
    else if( NX1INTS > 1 && NX2INTS > 1 ){
      double coef[NX1INTS*NX2INTS];
      PolVar<T> auxr[NX1INTS*NX2INTS];
      for( unsigned i=0, k=0; i<NX1INTS; i++ ){
        for( unsigned j=0; j<NX2INTS; j++, k++ ){
          coef[k] = ( X1KNOT[i+1] - X1KNOT[i] ) * ( X2KNOT[j+1] - X2KNOT[j] );
	  auxr[k]  = *_append_aux( Op<T>::zeroone(), true );
	}
      }
      // Represent variable subdivision using linear binary transformation
      const std::vector< PolVar<T> >& aux1 = X1.BIN_subdiv( pOp );
      const std::vector< PolVar<T> >& aux2 = X2.BIN_subdiv( pOp );
      // Cut (46) in Gounaris et al. (2009), DOI: 10.1021/ie8016048
      _add_cut( pOp, PolCut<T>::EQ, X1L*X2L, NX1INTS*NX2INTS, auxr, coef, X1, X2L, X2, X1L, Y, -1. );
      // Cuts (50) in Gounaris et al. (2009), DOI: 10.1021/ie8016048
      for( unsigned i=0, k=0; i<NX1INTS; i++ ){
        for( unsigned j=0; j<NX2INTS; j++, k++ ){
          _add_cut( pOp, PolCut<T>::LE, 1.,  aux1[i], 1., aux2[j], 1., auxr[k], -1. );
          if( j ) _add_cut( pOp, PolCut<T>::LE, 0., auxr[k], 1., auxr[k-1], -1. );
          else    _add_cut( pOp, PolCut<T>::GE, 0., aux1[i], 1., auxr[i*NX2INTS], -1. );
          if( i ) _add_cut( pOp, PolCut<T>::LE, 0., auxr[k], 1., auxr[k-NX2INTS], -1. );
          else    _add_cut( pOp, PolCut<T>::GE, 0., aux2[j], 1., auxr[j],   -1. );
	}
      }
      break;
    }
    
    return false;
   }

   case Options::SOS2:
   default:
    return false;
  }
  
  return true;
}
#endif
template <typename T>
inline PolVar<T>
inv
( const PolVar<T>&Var1 )
{
  return pow( Var1, -1 );
}

template <typename T>
inline PolVar<T>
sqr
( const PolVar<T>&Var1 )
{
  FFGraph* dag = Var1._var.dag();
#ifdef MC__POLIMG_CHECK
  if( !dag || !dag->curOp() )
    throw typename PolImg<T>::Exceptions( PolImg<T>::Exceptions::NOTALLOWED );
#endif

  FFVar* pFFVarR = dag->curOp()->pres;
  PolVar<T>* pVarR = Var1._img->_append_var( pFFVarR, Op<T>::sqr( Var1._range ), true );
  return *pVarR;
}

template <typename T> inline void
PolImg<T>::_add_cuts_SQR
( const PolVar<T>*VarR, FFVar*pVar1, FFOp* pOp )
{
  assert( VarR && pVar1 );
  if( pVar1->cst() )
    _add_cut( pOp, PolCut<T>::EQ, mc::sqr( pVar1->num().val() ), *VarR, 1. );

  else if( !options.RELAX_MONOM ){
    auto itVar1 = _Vars.find( pVar1 );
    assert( itVar1 != _Vars.end() );
    _add_cut( pOp, *VarR, *itVar1->second );
  }

  else if( !options.RELAX_QUAD ){
    auto itVar1 = _Vars.find( pVar1 );
    assert( itVar1 != _Vars.end() );
    auto itCut = _add_cut( pOp, PolCut<T>::EQ, 0., *VarR, -1. );
    (*itCut)->append( *itVar1->second, *itVar1->second, 1. );
  }

  else{
    auto itVar1 = _Vars.find( pVar1 );
    assert( itVar1 != _Vars.end() );
    struct loc{ static std::pair<double,double> sqr
      ( const double x, const double*rusr, const int*iusr )
      { return std::make_pair( x*x, 2.*x ); }
    };
    _semilinear_cuts( pOp, *itVar1->second, Op<T>::l(itVar1->second->_range),
      Op<T>::u(itVar1->second->_range), *VarR, PolCut<T>::LE, loc::sqr );
    _sandwich_cuts( pOp, *itVar1->second, Op<T>::l(itVar1->second->_range),
      Op<T>::u(itVar1->second->_range), *VarR, Op<T>::l(VarR->_range), Op<T>::u(VarR->_range),
      PolCut<T>::GE, loc::sqr );
  }
}

template <typename T> inline void
PolImg<T>::_add_cuts_SQR
( const PolVar<T>*VarR, FFVar*pVar1 )
{
  return _add_cuts_SQR( VarR, pVar1, VarR->_var.ops().first );
}

template <typename T> inline void
PolImg<T>::append_cuts_SQR
( PolVar<T> const& VarR, PolVar<T> const& Var1, FFOp* pOp )
{
  assert( &Var1._var );
  return _add_cuts_SQR( &VarR, const_cast<FFVar*>(&Var1._var), pOp );

//  if( !options.RELAX_QUAD ){
//    auto itCut = _add_cut( pOp, PolCut<T>::EQ, 0., VarR, -1. );
//    (*itCut)->append( Var1, Var1, 1. );
//  }

//  else{
//    struct loc{ static std::pair<double,double> sqr
//      ( const double x, const double*rusr, const int*iusr )
//      { return std::make_pair( x*x, 2.*x ); }
//    };
//    _semilinear_cuts( pOp, Var1, Op<T>::l(Var1._range), Op<T>::u(Var1._range), VarR,
//      PolCut<T>::LE, loc::sqr );
//    _sandwich_cuts( pOp, Var1, Op<T>::l(Var1._range), Op<T>::u(Var1._range), VarR,
//      Op<T>::l(VarR._range), Op<T>::u(VarR._range), PolCut<T>::GE, loc::sqr );
//  }
}

template <typename T> inline bool
PolImg<T>::_add_LQ_SQR
( PolLQExpr<T>*&pLQ, const PolVar<T>*VarR, FFVar*pVar1 )
{
  if( options.RELAX_MONOM == 2 || options.RELAX_QUAD ) return false;
  if( !pLQ ) pLQ = _append_LQ( VarR );
  if( pVar1->cst() )
    pLQ->substitute( VarR, sqr( pVar1->num().val() ) );
  else{
    auto itVar1 = _Vars.find( pVar1 );
    pLQ->substitute( VarR, 1., itVar1->second, itVar1->second );
  }
  return true;
}

template <typename T>
inline PolVar<T>
sqrt
( const PolVar<T>&Var1 )
{
  FFGraph* dag = Var1._var.dag();
#ifdef MC__POLIMG_CHECK
  if( !dag || !dag->curOp() )
    throw typename PolImg<T>::Exceptions( PolImg<T>::Exceptions::NOTALLOWED );
#endif

  FFVar* pFFVarR = dag->curOp()->pres;
  PolVar<T>* pVarR = Var1._img->_append_var( pFFVarR, Op<T>::sqrt( Var1._range ), true );
  return *pVarR;
}

template <typename T> inline void
PolImg<T>::_add_cuts_SQRT
( const PolVar<T>*VarR, FFVar*pVar1 )
{
  if( pVar1->cst() ){
    _add_cut( VarR->_var.ops().first, PolCut<T>::EQ, std::sqrt( pVar1->num().val() ), *VarR, 1. );
    return;
  }

  else if( !options.RELAX_NLIN ){
    auto itVar1 = _Vars.find( pVar1 );
    _add_cut( VarR->_var.ops().first, *VarR, *itVar1->second );
  }

  else{
    auto itVar1 = _Vars.find( pVar1 );
    struct loc{ 
      static std::pair<double,double> sqr
        ( const double x, const double*rusr, const int*iusr )
        { return std::make_pair( x*x, 2.*x ); }
      static std::pair<double,double> sqrt
        ( const double x, const double*rusr, const int*iusr )
        { return std::make_pair( std::sqrt(x), 1/(2*std::sqrt(x)) ); }
    };
    _semilinear_cuts( VarR->_var.ops().first, *itVar1->second, Op<T>::l(itVar1->second->_range),
      Op<T>::u(itVar1->second->_range), *VarR, PolCut<T>::GE, loc::sqrt );
    _sandwich_cuts( VarR->_var.ops().first, *VarR, Op<T>::l(VarR->_range), Op<T>::u(VarR->_range),
      *itVar1->second, Op<T>::l(itVar1->second->_range), Op<T>::u(itVar1->second->_range), PolCut<T>::GE, loc::sqr );
  }
}

template <typename T>
inline PolVar<T>
pow
( const PolVar<T>&Var1, const double dExp )
{
  FFGraph* dag = Var1._var.dag();
#ifdef MC__POLIMG_CHECK
  if( !dag || !dag->curOp() )
    throw typename PolImg<T>::Exceptions( PolImg<T>::Exceptions::NOTALLOWED );
#endif

  FFVar* pFFVarR = dag->curOp()->pres;
  PolVar<T>* pVarR = Var1._img->_append_var( pFFVarR, Op<T>::pow( Var1._range, dExp ), true );
  return *pVarR;
}

template <typename T> inline void
PolImg<T>::_add_cuts_DPOW
( const PolVar<T>*VarR, FFVar*pVar1, const double dExp )
{
  if( pVar1->cst() ){
    _add_cut( VarR->_var.ops().first, PolCut<T>::EQ, std::pow( pVar1->num().val(), dExp ), *VarR, 1. );
    return;
  }

  if( !options.RELAX_NLIN ){
    auto itVar1 = _Vars.find( pVar1 );
    _add_cut( VarR->_var.ops().first, *VarR, *itVar1->second );
    return;
  }

  struct loc{ static std::pair<double,double> dpow
    ( const double x, const double*rusr, const int*iusr )
    { return std::make_pair( std::pow(x,*rusr), *rusr*std::pow(x,*rusr-1) ); }
  };
  auto itVar1 = _Vars.find( pVar1 );

  // Convex case
  if( dExp >= 1 || dExp <= 0 ){
    _semilinear_cuts( VarR->_var.ops().first, *itVar1->second, Op<T>::l(itVar1->second->_range),
      Op<T>::u(itVar1->second->_range), *VarR, PolCut<T>::LE, loc::dpow, &dExp );
    _sandwich_cuts( VarR->_var.ops().first, *itVar1->second, Op<T>::l(itVar1->second->_range),
      Op<T>::u(itVar1->second->_range), *VarR, Op<T>::l(VarR->_range), Op<T>::u(VarR->_range),
      PolCut<T>::GE, loc::dpow, &dExp );
  }
  
  // Concave case
  else{
    _semilinear_cuts( VarR->_var.ops().first, *itVar1->second, Op<T>::l(itVar1->second->_range),
      Op<T>::u(itVar1->second->_range), *VarR, PolCut<T>::GE, loc::dpow, &dExp );
    _sandwich_cuts( VarR->_var.ops().first, *itVar1->second, Op<T>::l(itVar1->second->_range),
      Op<T>::u(itVar1->second->_range), *VarR, Op<T>::l(VarR->_range), Op<T>::u(VarR->_range),
      PolCut<T>::LE, loc::dpow, &dExp );
  }
}

template <typename T>
inline PolVar<T>
pow
( const PolVar<T>&Var1, const int iExp )
{
#ifdef MC__POLIMG_CHECK
  if( iExp == 0 || iExp == 1 || iExp == 2 )
    throw typename PolImg<T>::Exceptions( PolImg<T>::Exceptions::NOTALLOWED );
#endif
  FFGraph* dag = Var1._var.dag();
#ifdef MC__POLIMG_CHECK
  if( !dag || !dag->curOp() )
    throw typename PolImg<T>::Exceptions( PolImg<T>::Exceptions::NOTALLOWED );
#endif

  FFVar* pFFVarR = dag->curOp()->pres;
  PolVar<T>* pVarR = Var1._img->_append_var( pFFVarR, Op<T>::pow( Var1._range, iExp ), true );
  return *pVarR;
}

template <typename T> inline void
PolImg<T>::_add_cuts_IPOW
( const PolVar<T>*VarR, FFVar*pVar1, const int iExp )
{
  assert( iExp > 2 || iExp < 0 );
  if( pVar1->cst() ){
    _add_cut( VarR->_var.ops().first, PolCut<T>::EQ, std::pow( pVar1->num().val(), iExp ), *VarR, 1. );
    return;
  }
  
  // No linearization
  if( options.RELAX_MONOM == 0 ){ //|| options.RELAX_MONOM == 1 ){
    auto itVar1 = _Vars.find( pVar1 );
    _add_cut( VarR->_var.ops().first, *VarR, *itVar1->second );
    return;
  }

  auto itVar1 = _Vars.find( pVar1 );
  struct loc{
    static std::pair<double,double> pow
      ( const double x, const double*rusr, const int*iusr )
      { return std::make_pair( std::pow(x,*iusr), *iusr*std::pow(x,*iusr-1) ); }
    static std::pair<double,double> powcvu
      ( const double x, const double*rusr, const int*iusr )
      { return std::make_pair( std::pow(x,*iusr), *iusr*std::pow(x,*iusr-1) ); }
  };
  // Positive even exponent term
  if( iExp > 0 && !(iExp%2) ){
    _semilinear_cuts( VarR->_var.ops().first, *itVar1->second, Op<T>::l(itVar1->second->_range),
      Op<T>::u(itVar1->second->_range), *VarR, PolCut<T>::LE, loc::pow, 0, &iExp );
    _sandwich_cuts( VarR->_var.ops().first, *itVar1->second, Op<T>::l(itVar1->second->_range),
      Op<T>::u(itVar1->second->_range), *VarR, Op<T>::l(VarR->_range), Op<T>::u(VarR->_range),
      PolCut<T>::GE, loc::pow, 0, &iExp );
  }

  // Positive odd exponent term
  else if( iExp > 0 && options.ROOT_USE ){
    // -- Convex Portion
    if( Op<T>::l(itVar1->second->_range) >= 0. ){
      _semilinear_cuts( VarR->_var.ops().first, *itVar1->second, Op<T>::l(itVar1->second->_range),
        Op<T>::u(itVar1->second->_range), *VarR, PolCut<T>::LE, loc::pow, 0, &iExp );
      _sandwich_cuts( VarR->_var.ops().first, *itVar1->second, Op<T>::l(itVar1->second->_range),
        Op<T>::u(itVar1->second->_range), *VarR, Op<T>::l(VarR->_range), Op<T>::u(VarR->_range),
        PolCut<T>::GE, loc::pow, 0, &iExp );
    }
    // -- Concave Portion
    else if( Op<T>::u(itVar1->second->_range) <= 0. ){
      _semilinear_cuts( VarR->_var.ops().first, *itVar1->second, Op<T>::l(itVar1->second->_range),
        Op<T>::u(itVar1->second->_range), *VarR, PolCut<T>::GE, loc::pow, 0, &iExp );
      _sandwich_cuts( VarR->_var.ops().first, *itVar1->second, Op<T>::l(itVar1->second->_range),
        Op<T>::u(itVar1->second->_range), *VarR, Op<T>::l(VarR->_range), Op<T>::u(VarR->_range),
        PolCut<T>::LE, loc::pow, 0, &iExp );
    }
    // -- Nonconvex/Nonconcave Portion
    else{
      switch( options.BREAKPOINT_TYPE ){
        case PolImg<T>::Options::BIN:
        case PolImg<T>::Options::SOS2:{
          const unsigned NKNOTS = itVar1->second->create_subdiv( Op<T>::l(itVar1->second->_range),
            Op<T>::u(itVar1->second->_range) ).size();
          if( NKNOTS > 2 ){
            struct dc{
              static std::pair<double,double> pow1
                ( const double x, const double*rusr, const int*iusr )
                { return std::make_pair( x<0?std::pow(x,*iusr):0., x<0?*iusr*std::pow(x,*iusr-1):0. ); }
              static std::pair<double,double> pow2
                ( const double x, const double*rusr, const int*iusr )
                { return std::make_pair( x>0?std::pow(x,*iusr):0., x>0?*iusr*std::pow(x,*iusr-1):0. ); }
            };
            PolVar<T>* Var3 = _append_aux( Op<T>::pow( itVar1->second->_range, iExp ), true );
            PolVar<T>* Var4 = _append_aux( Op<T>::pow( itVar1->second->_range, iExp ), true );
            _semilinear_cuts( VarR->_var.ops().first, *itVar1->second, Op<T>::l(itVar1->second->_range),
              Op<T>::u(itVar1->second->_range), *Var3, PolCut<T>::GE, dc::pow1, 0, &iExp );
            _sandwich_cuts( VarR->_var.ops().first, *itVar1->second, Op<T>::l(itVar1->second->_range),
              Op<T>::u(itVar1->second->_range), *Var3, Op<T>::l(VarR->_range), 0., PolCut<T>::LE, dc::pow1, 0, &iExp );
            _semilinear_cuts( VarR->_var.ops().first, *itVar1->second, Op<T>::l(itVar1->second->_range),
              Op<T>::u(itVar1->second->_range), *Var4, PolCut<T>::LE, dc::pow2, 0, &iExp );
            _sandwich_cuts( VarR->_var.ops().first, *itVar1->second, Op<T>::l(itVar1->second->_range),
              Op<T>::u(itVar1->second->_range), *Var4, 0., Op<T>::u(VarR->_range), PolCut<T>::GE, dc::pow2, 0, &iExp );
            _add_cut( VarR->_var.ops().first, PolCut<T>::EQ, 0., *VarR, 1., *Var3, -1., *Var4, -1. );
            break;
          }
          // No break in order to append other "normal" cuts if no breakpoints
        }
        case PolImg<T>::Options::NONE: default:{
          struct fct{ static std::pair<double,double> powoddfunc
            ( const double x, const double*rusr, const int*iusr )
            { return std::make_pair(
                ((*iusr-1)*x-(*iusr)*(*rusr))*std::pow(x,*iusr-1) + std::pow(*rusr,*iusr),
                (*iusr)*(*iusr-1)*(x-(*rusr))*std::pow(x,*iusr-2) ); }
          };
          double xJcc = Op<T>::u(itVar1->second->_range);
          xJcc = _newton( Op<T>::l(itVar1->second->_range), Op<T>::l(itVar1->second->_range), 0.,
            fct::powoddfunc, options.ROOT_TOL, options.ROOT_MAXIT, &xJcc, &iExp );
          if( mc::isequal( xJcc, Op<T>::l(itVar1->second->_range) ) )
            _semilinear_cuts( VarR->_var.ops().first, *itVar1->second, Op<T>::l(itVar1->second->_range),
              Op<T>::u(itVar1->second->_range), *VarR, PolCut<T>::LE, loc::pow, 0, &iExp );
          else
            _sandwich_cuts( VarR->_var.ops().first, *itVar1->second, Op<T>::l(itVar1->second->_range), xJcc,
              *VarR, Op<T>::l(VarR->_range), Op<T>::u(VarR->_range), PolCut<T>::LE, loc::pow, 0, &iExp );

          double xJcv = Op<T>::l(itVar1->second->_range);
          xJcv = _newton( Op<T>::u(itVar1->second->_range), 0., Op<T>::u(itVar1->second->_range),
            fct::powoddfunc, options.ROOT_TOL, options.ROOT_MAXIT, &xJcv, &iExp );
          if( mc::isequal( xJcv, Op<T>::u(itVar1->second->_range) ) )
            _semilinear_cuts( VarR->_var.ops().first, *itVar1->second, Op<T>::l(itVar1->second->_range),
              Op<T>::u(itVar1->second->_range), *VarR, PolCut<T>::GE, loc::pow, 0, &iExp );
          else
            _sandwich_cuts( VarR->_var.ops().first, *itVar1->second, xJcv, Op<T>::u(itVar1->second->_range),
              *VarR, Op<T>::l(VarR->_range), Op<T>::u(VarR->_range), PolCut<T>::GE, loc::pow, 0, &iExp );
          break;
        }
      }
    }
  }

  // Negative exponent term
  else if( iExp < 0 ){
    // -- Convex Case
    if( !(iExp%2) || Op<T>::l(itVar1->second->_range) > 0. ){
      _semilinear_cuts( VarR->_var.ops().first, *itVar1->second, Op<T>::l(itVar1->second->_range),
        Op<T>::u(itVar1->second->_range), *VarR, PolCut<T>::LE, loc::pow, 0, &iExp );
      _sandwich_cuts( VarR->_var.ops().first, *itVar1->second, Op<T>::l(itVar1->second->_range),
        Op<T>::u(itVar1->second->_range), *VarR, Op<T>::l(VarR->_range), Op<T>::u(VarR->_range),
        PolCut<T>::GE, loc::pow, 0, &iExp );
    }
    // -- Concave Case
    else{
      _semilinear_cuts( VarR->_var.ops().first, *itVar1->second, Op<T>::l(itVar1->second->_range),
        Op<T>::u(itVar1->second->_range), *VarR, PolCut<T>::GE, loc::pow, 0, &iExp );
      _sandwich_cuts( VarR->_var.ops().first, *itVar1->second, Op<T>::l(itVar1->second->_range),
        Op<T>::u(itVar1->second->_range), *VarR, Op<T>::l(VarR->_range), Op<T>::u(VarR->_range),
        PolCut<T>::LE, loc::pow, 0, &iExp );
    }
  }
}

template <typename T>
inline PolVar<T>
cheb
( const PolVar<T>&Var1, const unsigned iOrd )
{
#ifdef MC__POLIMG_CHECK
  if( iOrd <= 2 )
    throw typename PolImg<T>::Exceptions( PolImg<T>::Exceptions::NOTALLOWED );
#endif
  FFGraph* dag = Var1._var.dag();
#ifdef MC__POLIMG_CHECK
  if( !dag || !dag->curOp() )
    throw typename PolImg<T>::Exceptions( PolImg<T>::Exceptions::NOTALLOWED );
#endif

  FFVar* pFFVarR = dag->curOp()->pres;
#ifdef MC__POLIMG_DEBUG_CHEB
  std::cout << "X = " << Var1._range << std::endl;
  std::cout << "Cheb(" << Var1._range << "," << iOrd << ") = " << Op<T>::cheb( Var1._range, iOrd ) << std::endl;
#endif
  PolVar<T>* pVarR = Var1._img->_append_var( pFFVarR, Op<T>::cheb( Var1._range, iOrd ), true );
  return *pVarR;
}

template <typename T> inline void
PolImg<T>::_add_cuts_CHEB
( const PolVar<T>*VarR, FFVar*pVar1, const unsigned iOrd )
{
  if( pVar1->cst() )
    _add_cut( VarR->_var.ops().first, PolCut<T>::EQ, mc::cheb( pVar1->num().val(), iOrd ), *VarR, 1. );
  
  // No linearization
  else if( options.RELAX_MONOM == 0 ){ //|| options.RELAX_MONOM == 1 ){
    auto itVar1 = _Vars.find( pVar1 );
    _add_cut( VarR->_var.ops().first, *VarR, *itVar1->second );
    return;
  }

  auto itVar1 = _Vars.find( pVar1 );
  struct loc{
    static std::pair<double,double> cheb
      ( const double x, const double*rusr, const int*iusr )
      { return std::make_pair( mc::cheb(x,*iusr), *iusr*mc::cheb2(x,*iusr-1) ); }
  };      
  const int iusr = iOrd;

  // Positive even order
  if( iOrd > 0 && !(iOrd%2) ){
    _semilinear_cuts( VarR->_var.ops().first, *itVar1->second, Op<T>::l(itVar1->second->_range),
      Op<T>::u(itVar1->second->_range), *VarR, PolCut<T>::LE, loc::cheb, 0, &iusr );
    double xJL = std::cos(mc::PI*(1.+1./(double)iOrd));
    _sandwich_cuts( VarR->_var.ops().first, *itVar1->second, Op<T>::l(itVar1->second->_range), xJL,
      *VarR, Op<T>::l(VarR->_range), Op<T>::u(VarR->_range), PolCut<T>::GE, loc::cheb, 0, &iusr );
    double xJU = std::cos(mc::PI/(double)iOrd);
    _sandwich_cuts( VarR->_var.ops().first, *itVar1->second, xJU, Op<T>::u(itVar1->second->_range),
      *VarR, Op<T>::l(VarR->_range), Op<T>::u(VarR->_range), PolCut<T>::GE, loc::cheb, 0, &iusr );
  }

  // Positive odd order
  else{
    double xJL = std::cos(mc::PI*(1.+1./(double)iOrd));
    _sandwich_cuts( VarR->_var.ops().first, *itVar1->second, Op<T>::l(itVar1->second->_range), xJL,
      *VarR, Op<T>::l(VarR->_range), Op<T>::u(VarR->_range), PolCut<T>::LE, loc::cheb, 0, &iusr );
    double xJU = std::cos(mc::PI/(double)iOrd);
    _sandwich_cuts( VarR->_var.ops().first, *itVar1->second, xJU, Op<T>::u(itVar1->second->_range),
      *VarR, Op<T>::l(VarR->_range), Op<T>::u(VarR->_range), PolCut<T>::GE, loc::cheb, 0, &iusr );
  }
}

template <typename T>
inline PolVar<T>
prod
( const unsigned nVar, const PolVar<T>*pVar )
{
#ifdef MC__POLIMG_CHECK
  if( nVar <= 2 )
    throw typename PolImg<T>::Exceptions( PolImg<T>::Exceptions::NOTALLOWED );
#endif
  FFGraph* dag = pVar->_var.dag();
#ifdef MC__POLIMG_CHECK
  if( !dag || !dag->curOp() )
    throw typename PolImg<T>::Exceptions( PolImg<T>::Exceptions::NOTALLOWED );
#endif

  FFVar* pFFVarR = dag->curOp()->pres;
  std::vector<T> IVar; IVar.reserve( nVar );
  for( unsigned i=0; i<nVar; i++ ) IVar.push_back( pVar[i]._range );
  PolVar<T>* pVarR = pVar->_img->_append_var( pFFVarR, Op<T>::prod( IVar.size(), IVar.data() ), true );
  return *pVarR;
}

template <typename T> inline bool
PolImg<T>::_subset_incr
( const unsigned ntot, std::vector<unsigned>&ndx )
{
  auto rit=ndx.rbegin();
  for( unsigned i=0; rit!=ndx.rend(); ++rit, i++ )
    if( *rit < ntot-i-1 ){ (*rit)++; return true; }
  if( ntot < ndx.size()+2 ) return false;
  for( unsigned i=0; i<ndx.size(); i++ ) ndx.at(i) = i;
  ndx.push_back( ndx.size() ); ndx.push_back( ndx.size() );
  return true;
}

template <typename T> inline void
PolImg<T>::_add_cuts_PROD
( const PolVar<T>*VarR, std::vector<FFVar*>vVar )
{
  std::vector<PolVar<T>> vPolVar;
  double scalR = 1., radR = 1.;
  double isCen = true;
  auto it=vVar.begin(); 
  for( unsigned i=0; it!=vVar.end(); ++it, i++ ){
    // Detect constant operands
    if( (*it)->cst() ){
      scalR *= (*it)->num().val();
      continue;
    }
    // Look for corresponding variable in polyhedral image
    auto itVar = _Vars.find( *it );
    assert( itVar->first );
    if( isequal( Op<T>::diam(itVar->second->_range), 0. ) ){
      scalR *= Op<T>::mid(itVar->second->_range);
      continue;
    }
    vPolVar.push_back( *(itVar->second) );
    if( !isequal( Op<T>::mid(itVar->second->_range), 0. ) )
      isCen = false;
    radR *= Op<T>::u(itVar->second->_range);
  }
#ifdef MC__POLIMG_DEBUG_PROD
  std::cout << "CENTERED:" << (isCen?"Y":"N") << std::endl;
#endif 
  
  // Case: Zero scaling factor in multilinear term
  if( isequal( scalR, 0. ) ){
    _add_cut( VarR->_var.ops().first, PolCut<T>::EQ, 0., *VarR, 1. );
    return;
  }

  // Case: Non-centered multilinear term
  if( !isCen ){
    throw typename PolImg<T>::Exceptions( PolImg<T>::Exceptions::NOTALLOWED );
//    for( unsigned i=0; i<vPolVar.size(); i++ ){
//    PolVar<T>* auxVar = _append_aux( Op<T>::pow( itVar1->second->_range, iExp ), true );
//    _add_cuts_TIMES( const PolVar<T>*VarR, FFVar*pVar1, FFVar*pVar2 )
  }

  // Case: Centered multilinear term
  std::vector<double> cPolVar( vPolVar.size() );
  for( std::vector<unsigned> ndxCut; ; ){
#ifdef MC__POLIMG_DEBUG_PROD
    std::cout << "{";
    for( auto&& idx : ndxCut ) std::cout << " " << idx;
    std::cout << " }\n";
#endif   
    // Subscase: Odd-sized multilinear term
    if( vPolVar.size()%2 ){
      for( unsigned i=0, k=0; i<vPolVar.size(); i++ ){
        if( ndxCut.size() > k && ndxCut[k] == i ){
          cPolVar[i] = 1. / Op<T>::u( vPolVar[i]._range );
          k++;
        }
        else
          cPolVar[i] = -1. / Op<T>::u( vPolVar[i]._range );
      }
      _add_cut( VarR->_var.ops().first, PolCut<T>::LE, vPolVar.size()-1.,
        vPolVar.size(), vPolVar.data(), cPolVar.data(), *VarR, 1./(scalR*radR) );
      _add_cut( VarR->_var.ops().first, PolCut<T>::GE, 1.-vPolVar.size(),
        vPolVar.size(), vPolVar.data(), cPolVar.data(), *VarR, 1./(scalR*radR) );
    }
    // Subscase: Even-sized multilinear term
    else{
      for( unsigned i=0, k=0; i<vPolVar.size(); i++ ){
        if( ndxCut.size() > k && ndxCut[k] == i ){
          cPolVar[i] = 1. / Op<T>::u( vPolVar[i]._range );
          k++;
        }
        else
          cPolVar[i] = -1. / Op<T>::u( vPolVar[i]._range );
      }
      _add_cut( VarR->_var.ops().first, PolCut<T>::GE, 1.-vPolVar.size(),
        vPolVar.size(), vPolVar.data(), cPolVar.data(), *VarR, -1./(scalR*radR) );
      _add_cut( VarR->_var.ops().first, PolCut<T>::GE, 1.-vPolVar.size(),
        vPolVar.size(), vPolVar.data(), cPolVar.data(), *VarR, 1./(scalR*radR) );
    }
    // Generate next index set
    if( !_subset_incr( vPolVar.size(), ndxCut ) ) break;
  }
#ifdef MC__POLIMG_DEBUG_PROD
  std::cout << "vVar: " << vVar.size() << "  vPolVar: " << vPolVar.size() << std::endl;
#endif
}

template <typename T>
inline PolVar<T>
exp
( const PolVar<T>&Var1 )
{
  FFGraph* dag = Var1._var.dag();
#ifdef MC__POLIMG_CHECK
  if( !dag || !dag->curOp() )
    throw typename PolImg<T>::Exceptions( PolImg<T>::Exceptions::NOTALLOWED );
#endif

  FFVar* pFFVarR = dag->curOp()->pres;
  PolVar<T>* pVarR = Var1._img->_append_var( pFFVarR, Op<T>::exp( Var1._range ), true );
  return *pVarR;
}

template <typename T> inline void
PolImg<T>::_add_cuts_EXP
( const PolVar<T>*VarR, FFVar*pVar1 )
{
  if( pVar1->cst() ){
    _add_cut( VarR->_var.ops().first, PolCut<T>::EQ, std::exp( pVar1->num().val() ), *VarR, 1. );
    return;
  }
  
  if( !options.RELAX_NLIN ){
    auto itVar1 = _Vars.find( pVar1 );
    _add_cut( VarR->_var.ops().first, *VarR, *itVar1->second );
    return;
  }

  auto itVar1 = _Vars.find( pVar1 );
  struct loc{ static std::pair<double,double> exp
    ( const double x, const double*rusr, const int*iusr )
    { return std::make_pair( std::exp(x), std::exp(x) ); }
  };
  _semilinear_cuts( VarR->_var.ops().first, *itVar1->second, Op<T>::l(itVar1->second->_range),
    Op<T>::u(itVar1->second->_range), *VarR, PolCut<T>::LE, loc::exp );
  _sandwich_cuts( VarR->_var.ops().first, *itVar1->second, Op<T>::l(itVar1->second->_range),
    Op<T>::u(itVar1->second->_range), *VarR, Op<T>::l(VarR->_range), Op<T>::u(VarR->_range),
    PolCut<T>::GE, loc::exp );
}

template <typename T>
inline PolVar<T>
log
( const PolVar<T>&Var1 )
{
  FFGraph* dag = Var1._var.dag();
#ifdef MC__POLIMG_CHECK
  if( !dag || !dag->curOp() )
    throw typename PolImg<T>::Exceptions( PolImg<T>::Exceptions::NOTALLOWED );
#endif

  FFVar* pFFVarR = dag->curOp()->pres;
  PolVar<T>* pVarR = Var1._img->_append_var( pFFVarR, Op<T>::log( Var1._range ), true );
  return *pVarR;
}

template <typename T> inline void
PolImg<T>::_add_cuts_LOG
( const PolVar<T>*VarR, FFVar*pVar1 )
{
  if( pVar1->cst() ){
    _add_cut( VarR->_var.ops().first, PolCut<T>::EQ, std::log( pVar1->num().val() ), *VarR, 1. );
    return;
  }

  if( !options.RELAX_NLIN ){
    auto itVar1 = _Vars.find( pVar1 );
    _add_cut( VarR->_var.ops().first, *VarR, *itVar1->second );
    return;
  }

  auto itVar1 = _Vars.find( pVar1 );
  struct loc{
    static std::pair<double,double> exp
      ( const double x, const double*rusr, const int*iusr )
      { return std::make_pair( std::exp(x), std::exp(x) ); }
    static std::pair<double,double> log
      ( const double x, const double*rusr, const int*iusr )
      { return std::make_pair( std::log(x), 1/x ); }
  };
  _semilinear_cuts( VarR->_var.ops().first, *itVar1->second, Op<T>::l(itVar1->second->_range),
    Op<T>::u(itVar1->second->_range), *VarR, PolCut<T>::GE, loc::log );
  _sandwich_cuts( VarR->_var.ops().first, *VarR, Op<T>::l(VarR->_range),
    Op<T>::u(VarR->_range), *itVar1->second, Op<T>::l(itVar1->second->_range), Op<T>::u(itVar1->second->_range),
    PolCut<T>::GE, loc::exp );
}

template <typename T>
inline PolVar<T>
xlog
( const PolVar<T>&Var1 )
{
  FFGraph* dag = Var1._var.dag();
#ifdef MC__POLIMG_CHECK
  if( !dag || !dag->curOp() )
    throw typename PolImg<T>::Exceptions( PolImg<T>::Exceptions::NOTALLOWED );
#endif

  FFVar* pFFVarR = dag->curOp()->pres;
  PolVar<T>* pVarR = Var1._img->_append_var( pFFVarR, Op<T>::xlog( Var1._range ), true );
  return *pVarR;
}

template <typename T> inline void
PolImg<T>::_add_cuts_XLOG
( const PolVar<T>*VarR, FFVar*pVar1 )
{
  if( pVar1->cst() ){
    _add_cut( VarR->_var.ops().first, PolCut<T>::EQ, mc::xlog( pVar1->num().val() ), *VarR, 1. );
    return;
  }
  
  auto itVar1 = _Vars.find( pVar1 );
  struct loc{ static std::pair<double,double> xlog
    ( const double x, const double*rusr, const int*iusr )
    { return std::make_pair( mc::xlog(x), std::log(x)+1. ); }
  };
  _semilinear_cuts( VarR->_var.ops().first, *itVar1->second, Op<T>::l(itVar1->second->_range),
    Op<T>::u(itVar1->second->_range), *VarR, PolCut<T>::LE, loc::xlog );
  _sandwich_cuts( VarR->_var.ops().first, *itVar1->second, Op<T>::l(itVar1->second->_range),
    Op<T>::u(itVar1->second->_range), *VarR, Op<T>::l(VarR->_range), Op<T>::u(VarR->_range),
    PolCut<T>::GE, loc::xlog );
}

template <typename T>
inline PolVar<T>
cos
( const PolVar<T>&Var1 )
{
  FFGraph* dag = Var1._var.dag();
#ifdef MC__POLIMG_CHECK
  if( !dag || !dag->curOp() )
    throw typename PolImg<T>::Exceptions( PolImg<T>::Exceptions::NOTALLOWED );
#endif

  FFVar* pFFVarR = dag->curOp()->pres;
  PolVar<T>* pVarR = Var1._img->_append_var( pFFVarR, Op<T>::cos( Var1._range ), true );
  return *pVarR;
}

template <typename T> inline void
PolImg<T>::_add_cuts_COS
( const PolVar<T>*VarR, FFVar*pVar1 )
{
  if( pVar1->cst() ){
    _add_cut( VarR->_var.ops().first, PolCut<T>::EQ, std::cos( pVar1->num().val() ), *VarR, 1. );
    return;
  }

  if( !options.RELAX_NLIN ){
    auto itVar1 = _Vars.find( pVar1 );
    _add_cut( VarR->_var.ops().first, *VarR, *itVar1->second );
    return;
  }

  auto itVar1 = _Vars.find( pVar1 );
  struct loc{
    static std::pair<double,double> cos
      ( const double x, const double*rusr, const int*iusr )
      { return std::make_pair( std::cos(x), -std::sin(x) ); }
  };
  struct fct{ static std::pair<double,double> cosfunc
    ( const double x, const double*rusr, const int*iusr )
    { return std::make_pair(
        (x-*rusr)*std::sin(x)+std::cos(x)-std::cos(*rusr),
        (x-*rusr)*std::cos(x) ); }
  };

  // Convex relaxation
  T IVar1 = itVar1->second->_range;
  int kL = std::ceil( -0.5*(1.+Op<T>::l(IVar1)/PI) ); // offset for xL to be in [-pi,pi]
  double dxL = 2.*PI*kL; 
  double xL1 = Op<T>::l(IVar1)+dxL, xU1 = Op<T>::u(IVar1)+dxL;
  assert( xL1 >= -PI && xL1 <= PI );
 
  if( xU1 >= PI ){
    const int kU = std::ceil( -0.5*(1.+Op<T>::u(IVar1)/PI) );
    const double dxU = 2.*PI*kU; 
    const double xU2 = Op<T>::u(IVar1)+dxU;
    assert( xU2 >= -PI && xU2 <= PI );

    double xJcv1 = Op<T>::l(IVar1);
    if( xL1 <= PI/2. ) xJcv1 = _newton( PI-dxL, Op<T>::l(IVar1), PI-dxL, fct::cosfunc,
      options.ROOT_TOL, options.ROOT_MAXIT, &xJcv1, 0 );
    double xJcv2 = Op<T>::u(IVar1);
    if( xU2 >= -PI/2. ) xJcv2 = _newton( -PI-dxU, -PI-dxU, Op<T>::u(IVar1), fct::cosfunc,
      options.ROOT_TOL, options.ROOT_MAXIT, &xJcv2, 0 );
    _sandwich_cuts( VarR->_var.ops().first, *itVar1->second, xJcv1, PI-dxL,
      *VarR, Op<T>::l(VarR->_range), Op<T>::u(VarR->_range), PolCut<T>::GE, loc::cos );
    _sandwich_cuts( VarR->_var.ops().first, *itVar1->second, -PI-dxU, xJcv2,
      *VarR, Op<T>::l(VarR->_range), Op<T>::u(VarR->_range), PolCut<T>::GE, loc::cos );
  }
  
  else if( xL1 >= PI/2. ){
    _sandwich_cuts( VarR->_var.ops().first, *itVar1->second, Op<T>::l(IVar1), Op<T>::u(IVar1),
      *VarR, Op<T>::l(VarR->_range), Op<T>::u(VarR->_range), PolCut<T>::GE, loc::cos );
  }
  
  else if( xL1 >= -PI/2. && xU1 <= PI/2. ){
    _semilinear_cuts( VarR->_var.ops().first, *itVar1->second, Op<T>::l(IVar1), Op<T>::u(IVar1),
      *VarR, PolCut<T>::GE, loc::cos );
  }

  else if( xL1 >= -PI/2. ){
    double xJcv1 = Op<T>::l(IVar1);
    xJcv1 = _newton( Op<T>::u(IVar1), Op<T>::l(IVar1), Op<T>::u(IVar1), fct::cosfunc,
      options.ROOT_TOL, options.ROOT_MAXIT, &xJcv1, 0 );
    if( mc::isequal( xJcv1, Op<T>::u(IVar1) ) ){
      _semilinear_cuts( VarR->_var.ops().first, *itVar1->second, Op<T>::l(IVar1), Op<T>::u(IVar1),
        *VarR, PolCut<T>::GE, loc::cos );
    }
    else{
      _sandwich_cuts( VarR->_var.ops().first, *itVar1->second, xJcv1, Op<T>::u(IVar1), 
        *VarR, Op<T>::l(VarR->_range), Op<T>::u(VarR->_range), PolCut<T>::GE, loc::cos );
    }
  }

  else{
    double xJcv1 = Op<T>::u(IVar1);
    xJcv1 = _newton( Op<T>::l(IVar1), Op<T>::l(IVar1), Op<T>::u(IVar1), fct::cosfunc,
      options.ROOT_TOL, options.ROOT_MAXIT, &xJcv1, 0 );
    if( mc::isequal( xJcv1, Op<T>::l(IVar1) ) ){
      _semilinear_cuts( VarR->_var.ops().first, *itVar1->second, Op<T>::l(IVar1), Op<T>::u(IVar1),
        *VarR, PolCut<T>::GE, loc::cos );
    }
    else{
      _sandwich_cuts( VarR->_var.ops().first, *itVar1->second, Op<T>::l(IVar1), xJcv1,
        *VarR, Op<T>::l(VarR->_range), Op<T>::u(VarR->_range), PolCut<T>::GE, loc::cos );
    }
  }

  // Concave relaxation
  kL = std::ceil( -0.5*(2.+Op<T>::l(IVar1)/PI) );
  dxL = 2.*PI*kL;
  xL1 = Op<T>::l(IVar1)+dxL, xU1 = Op<T>::u(IVar1)+dxL;
  assert( xL1 >= -2.*PI && xL1 <= 0. );

  if( xU1 >= 0. ){
    const int kU = std::ceil( -0.5*(2.+Op<T>::u(IVar1)/PI) );
    const double dxU = 2.*PI*kU; 
    const double xU2 = Op<T>::u(IVar1)+dxU;
    assert( xU2 >= -2.*PI && xU2 <= 0. );

    double xJcc1 = Op<T>::l(IVar1);
    if( xL1 <= -PI/2. ) xJcc1 = _newton( -dxL, Op<T>::l(IVar1), -dxL, fct::cosfunc,
      options.ROOT_TOL, options.ROOT_MAXIT, &xJcc1, 0 );
    double xJcc2 = Op<T>::u(IVar1);
    if( xU2 >= -3.*PI/2. ) xJcc2 = _newton( -2.*PI-dxU, -2.*PI-dxU, Op<T>::u(IVar1),
      fct::cosfunc, options.ROOT_TOL, options.ROOT_MAXIT, &xJcc2, 0 );
    _sandwich_cuts( VarR->_var.ops().first, *itVar1->second, xJcc1, -dxL,
      *VarR, Op<T>::l(VarR->_range), Op<T>::u(VarR->_range), PolCut<T>::LE, loc::cos );
    _sandwich_cuts( VarR->_var.ops().first, *itVar1->second, -2.*PI-dxU, xJcc2,
      *VarR, Op<T>::l(VarR->_range), Op<T>::u(VarR->_range), PolCut<T>::LE, loc::cos );
  }

  else if( xL1 >= -PI/2. ){
    _sandwich_cuts( VarR->_var.ops().first, *itVar1->second, Op<T>::l(IVar1), Op<T>::u(IVar1),
      *VarR, Op<T>::l(VarR->_range), Op<T>::u(VarR->_range), PolCut<T>::LE, loc::cos );
  }
  
  else if( xL1 >= -3.*PI/2. && xU1 <= -PI/2. ){
    _semilinear_cuts( VarR->_var.ops().first, *itVar1->second, Op<T>::l(IVar1), Op<T>::u(IVar1),
      *VarR, PolCut<T>::LE, loc::cos );
  }

  else if( xL1 >= -3.*PI/2. ){
    double xJcc1 = Op<T>::l(IVar1);
    xJcc1 = _newton( Op<T>::u(IVar1), Op<T>::l(IVar1), Op<T>::u(IVar1), fct::cosfunc,
      options.ROOT_TOL, options.ROOT_MAXIT, &xJcc1, 0 );
    if( mc::isequal( xJcc1, Op<T>::u(IVar1) ) ){
      _semilinear_cuts( VarR->_var.ops().first, *itVar1->second, Op<T>::l(IVar1), Op<T>::u(IVar1),
        *VarR, PolCut<T>::LE, loc::cos );
    }
    else{
      _sandwich_cuts( VarR->_var.ops().first, *itVar1->second, xJcc1, Op<T>::u(IVar1), 
        *VarR, Op<T>::l(VarR->_range), Op<T>::u(VarR->_range), PolCut<T>::LE, loc::cos );
    }
  }

  else{
    double xJcc1 = Op<T>::u(IVar1);
    xJcc1 = _newton( Op<T>::l(IVar1), Op<T>::l(IVar1), Op<T>::u(IVar1), fct::cosfunc,
      options.ROOT_TOL, options.ROOT_MAXIT, &xJcc1, 0 );
    if( mc::isequal( xJcc1, Op<T>::l(IVar1) ) ){
      _semilinear_cuts( VarR->_var.ops().first, *itVar1->second, Op<T>::l(IVar1), Op<T>::u(IVar1),
        *VarR, PolCut<T>::LE, loc::cos );
    }
    else{
      _sandwich_cuts( VarR->_var.ops().first, *itVar1->second, Op<T>::l(IVar1), xJcc1,
        *VarR, Op<T>::l(VarR->_range), Op<T>::u(VarR->_range), PolCut<T>::LE, loc::cos );
    }
  }
}

template <typename T>
inline PolVar<T>
sin
( const PolVar<T>&Var1 )
{
  FFGraph* dag = Var1._var.dag();
#ifdef MC__POLIMG_CHECK
  if( !dag || !dag->curOp() )
    throw typename PolImg<T>::Exceptions( PolImg<T>::Exceptions::NOTALLOWED );
#endif

  FFVar* pFFVarR = dag->curOp()->pres;
  PolVar<T>* pVarR = Var1._img->_append_var( pFFVarR, Op<T>::sin( Var1._range ), true );
  return *pVarR;
}

template <typename T> inline void
PolImg<T>::_add_cuts_SIN
( const PolVar<T>*VarR, FFVar*pVar1 )
{
  if( pVar1->cst() ){
    _add_cut( VarR->_var.ops().first, PolCut<T>::EQ, std::sin( pVar1->num().val() ), *VarR, 1. );
    return;
  }

  if( !options.RELAX_NLIN ){
    auto itVar1 = _Vars.find( pVar1 );
    _add_cut( VarR->_var.ops().first, *VarR, *itVar1->second );
    return;
  }

  auto itVar1 = _Vars.find( pVar1 );
  struct loc{
    static std::pair<double,double> sin
      ( const double x, const double*rusr, const int*iusr )
      { return std::make_pair( std::sin(x), std::cos(x) ); }
  };
  struct fct{ static std::pair<double,double> sinfunc
    ( const double x, const double*rusr, const int*iusr )
    { return std::make_pair(
        (*rusr-x)*std::cos(x)+std::sin(x)-std::sin(*rusr),
        (x-*rusr)*std::sin(x) ); }
  };

  // Convex relaxation
  T IVar1 = itVar1->second->_range;
  int kL = std::ceil( -0.5*(0.5+Op<T>::l(IVar1)/PI) ); // offset for xL to be in [-pi/2,3pi/2]
  double dxL = 2.*PI*kL; 
  double xL1 = Op<T>::l(IVar1)+dxL, xU1 = Op<T>::u(IVar1)+dxL;
  assert( xL1 >= -PI/2. && xL1 <= 3.*PI/2. );

  if( xU1 >= 3.*PI/2. ){
    const int kU = std::ceil( -0.5*(0.5+Op<T>::u(IVar1)/PI) );
    const double dxU = 2.*PI*kU; 
    const double xU2 = Op<T>::u(IVar1)+dxU;
    assert( xU2 >= -PI/2. && xU2 <= 3.*PI/2. );

    double xJcv1 = Op<T>::l(IVar1);
    if( xL1 <= PI ) xJcv1 = _newton( 3.*PI/2.-dxL, Op<T>::l(IVar1), 3.*PI/2.-dxL, fct::sinfunc,
      options.ROOT_TOL, options.ROOT_MAXIT, &xJcv1, 0 );
    double xJcv2 = Op<T>::u(IVar1);
    if( xU2 >= 0. ) xJcv2 = _newton( -PI/2.-dxU, -PI/2.-dxU, Op<T>::u(IVar1), fct::sinfunc,
      options.ROOT_TOL, options.ROOT_MAXIT, &xJcv2, 0 );
    _sandwich_cuts( VarR->_var.ops().first, *itVar1->second, xJcv1, 3.*PI/2.-dxL,
      *VarR, Op<T>::l(VarR->_range), Op<T>::u(VarR->_range), PolCut<T>::GE, loc::sin );
    _sandwich_cuts( VarR->_var.ops().first, *itVar1->second, -PI/2.-dxU, xJcv2,
      *VarR, Op<T>::l(VarR->_range), Op<T>::u(VarR->_range), PolCut<T>::GE, loc::sin );
  }
 
  else if( xL1 >= PI ){
    _sandwich_cuts( VarR->_var.ops().first, *itVar1->second, Op<T>::l(IVar1), Op<T>::u(IVar1),
      *VarR, Op<T>::l(VarR->_range), Op<T>::u(VarR->_range), PolCut<T>::GE, loc::sin );
  }
  
  else if( xL1 >= 0. && xU1 <= PI ){
    _semilinear_cuts( VarR->_var.ops().first, *itVar1->second, Op<T>::l(IVar1), Op<T>::u(IVar1),
      *VarR, PolCut<T>::GE, loc::sin );
  }

  else if( xL1 >= 0. ){
    double xJcv1 = Op<T>::l(IVar1);
    xJcv1 = _newton( Op<T>::u(IVar1), Op<T>::l(IVar1), Op<T>::u(IVar1), fct::sinfunc,
      options.ROOT_TOL, options.ROOT_MAXIT, &xJcv1, 0 );
    if( mc::isequal( xJcv1, Op<T>::u(IVar1) ) ){
      _semilinear_cuts( VarR->_var.ops().first, *itVar1->second, Op<T>::l(IVar1), Op<T>::u(IVar1),
        *VarR, PolCut<T>::GE, loc::sin );
    }
    else{
      _sandwich_cuts( VarR->_var.ops().first, *itVar1->second, xJcv1, Op<T>::u(IVar1), 
        *VarR, Op<T>::l(VarR->_range), Op<T>::u(VarR->_range), PolCut<T>::GE, loc::sin );
    }
  }

  else{
    double xJcv1 = Op<T>::u(IVar1);
    xJcv1 = _newton( Op<T>::l(IVar1), Op<T>::l(IVar1), Op<T>::u(IVar1), fct::sinfunc,
      options.ROOT_TOL, options.ROOT_MAXIT, &xJcv1, 0 );
    if( mc::isequal( xJcv1, Op<T>::l(IVar1) ) ){
      _semilinear_cuts( VarR->_var.ops().first, *itVar1->second, Op<T>::l(IVar1), Op<T>::u(IVar1),
        *VarR, PolCut<T>::GE, loc::sin );
    }
    else{
      _sandwich_cuts( VarR->_var.ops().first, *itVar1->second, Op<T>::l(IVar1), xJcv1,
        *VarR, Op<T>::l(VarR->_range), Op<T>::u(VarR->_range), PolCut<T>::GE, loc::sin );
    }
  }

  // Concave relaxation
  kL = std::ceil( -0.5*(1.5+Op<T>::l(IVar1)/PI) );
  dxL = 2.*PI*kL;
  xL1 = Op<T>::l(IVar1)+dxL, xU1 = Op<T>::u(IVar1)+dxL;
  assert( xL1 >= -3.*PI/2. && xL1 <= PI/2. );

  if( xU1 >= PI/2. ){
    const int kU = std::ceil( -0.5*(1.5+Op<T>::u(IVar1)/PI) );
    const double dxU = 2.*PI*kU; 
    const double xU2 = Op<T>::u(IVar1)+dxU;
    assert( xU2 >= -3.*PI/2. && xU2 <= PI/2. );

    double xJcc1 = Op<T>::l(IVar1);
    if( xL1 <= 0. ) xJcc1 = _newton( PI/2.-dxL, Op<T>::l(IVar1), PI/2.-dxL, fct::sinfunc,
      options.ROOT_TOL, options.ROOT_MAXIT, &xJcc1, 0 );
    double xJcc2 = Op<T>::u(IVar1);
    if( xU2 >= -PI ) xJcc2 = _newton( -3.*PI/2.-dxU, -3.*PI/2.-dxU, Op<T>::u(IVar1),
      fct::sinfunc, options.ROOT_TOL, options.ROOT_MAXIT, &xJcc2, 0 );
    _sandwich_cuts( VarR->_var.ops().first, *itVar1->second, xJcc1, PI/2.-dxL,
      *VarR, Op<T>::l(VarR->_range), Op<T>::u(VarR->_range), PolCut<T>::LE, loc::sin );
    _sandwich_cuts( VarR->_var.ops().first, *itVar1->second, -3.*PI/2.-dxU, xJcc2,
      *VarR, Op<T>::l(VarR->_range), Op<T>::u(VarR->_range), PolCut<T>::LE, loc::sin );
  }

  else if( xL1 >= 0. ){
    _sandwich_cuts( VarR->_var.ops().first, *itVar1->second, Op<T>::l(IVar1), Op<T>::u(IVar1),
      *VarR, Op<T>::l(VarR->_range), Op<T>::u(VarR->_range), PolCut<T>::LE, loc::sin );
  }
  
  else if( xL1 >= -PI && xU1 <= 0. ){
    _semilinear_cuts( VarR->_var.ops().first, *itVar1->second, Op<T>::l(IVar1), Op<T>::u(IVar1),
      *VarR, PolCut<T>::LE, loc::sin );
  }

  else if( xL1 >= -PI ){
    double xJcc1 = Op<T>::l(IVar1);
    xJcc1 = _newton( Op<T>::u(IVar1), Op<T>::l(IVar1), Op<T>::u(IVar1), fct::sinfunc,
      options.ROOT_TOL, options.ROOT_MAXIT, &xJcc1, 0 );
    if( mc::isequal( xJcc1, Op<T>::u(IVar1) ) ){
      _semilinear_cuts( VarR->_var.ops().first, *itVar1->second, Op<T>::l(IVar1), Op<T>::u(IVar1),
        *VarR, PolCut<T>::LE, loc::sin );
    }
    else{
      _sandwich_cuts( VarR->_var.ops().first, *itVar1->second, xJcc1, Op<T>::u(IVar1), 
        *VarR, Op<T>::l(VarR->_range), Op<T>::u(VarR->_range), PolCut<T>::LE, loc::sin );
    }
  }

  else{
    double xJcc1 = Op<T>::u(IVar1);
    xJcc1 = _newton( Op<T>::l(IVar1), Op<T>::l(IVar1), Op<T>::u(IVar1), fct::sinfunc,
      options.ROOT_TOL, options.ROOT_MAXIT, &xJcc1, 0 );
    if( mc::isequal( xJcc1, Op<T>::l(IVar1) ) ){
      _semilinear_cuts( VarR->_var.ops().first, *itVar1->second, Op<T>::l(IVar1), Op<T>::u(IVar1),
        *VarR, PolCut<T>::LE, loc::sin );
    }
    else{
      _sandwich_cuts( VarR->_var.ops().first, *itVar1->second, Op<T>::l(IVar1), xJcc1,
        *VarR, Op<T>::l(VarR->_range), Op<T>::u(VarR->_range), PolCut<T>::LE, loc::sin );
    }
  }
}

template <typename T>
inline PolVar<T>
tan
( const PolVar<T>&Var1 )
{
  FFGraph* dag = Var1._var.dag();
#ifdef MC__POLIMG_CHECK
  if( !dag || !dag->curOp() )
    throw typename PolImg<T>::Exceptions( PolImg<T>::Exceptions::NOTALLOWED );
#endif

  FFVar* pFFVarR = dag->curOp()->pres;
  PolVar<T>* pVarR = Var1._img->_append_var( pFFVarR, Op<T>::tan( Var1._range ), true );
  return *pVarR;
}

template <typename T> inline void
PolImg<T>::_add_cuts_TAN
( const PolVar<T>*VarR, FFVar*pVar1 )
{
  if( pVar1->cst() ){
    _add_cut( VarR->_var.ops().first, PolCut<T>::EQ, std::tan( pVar1->num().val() ), *VarR, 1. );
    return;
  }  

  if( !options.RELAX_NLIN ){
    auto itVar1 = _Vars.find( pVar1 );
    _add_cut( VarR->_var.ops().first, *VarR, *itVar1->second );
    return;
  }

  if( !options.ROOT_USE ) return;

  auto itVar1 = _Vars.find( pVar1 );
  struct loc{
    static std::pair<double,double> tan
      ( const double x, const double*rusr, const int*iusr )
      { return std::make_pair( std::tan(x), 1+mc::sqr(std::tan(x)) ); }
  };

  // -- Convex Portion
  if( Op<T>::l(itVar1->second->_range) >= 0. ){
    _semilinear_cuts( VarR->_var.ops().first, *itVar1->second, Op<T>::l(itVar1->second->_range),
      Op<T>::u(itVar1->second->_range), *VarR, PolCut<T>::LE, loc::tan, 0, 0 );
    _sandwich_cuts( VarR->_var.ops().first, *itVar1->second, Op<T>::l(itVar1->second->_range),
      Op<T>::u(itVar1->second->_range), *VarR, Op<T>::l(VarR->_range), Op<T>::u(VarR->_range),
      PolCut<T>::GE, loc::tan, 0, 0 );
  }
  // -- Concave Portion
  else if( Op<T>::u(itVar1->second->_range) <= 0. ){
    _semilinear_cuts( VarR->_var.ops().first, *itVar1->second, Op<T>::l(itVar1->second->_range),
      Op<T>::u(itVar1->second->_range), *VarR, PolCut<T>::GE, loc::tan, 0, 0 );
    _sandwich_cuts( VarR->_var.ops().first, *itVar1->second, Op<T>::l(itVar1->second->_range),
      Op<T>::u(itVar1->second->_range), *VarR, Op<T>::l(VarR->_range), Op<T>::u(VarR->_range),
      PolCut<T>::LE, loc::tan, 0, 0 );
  }
  // -- Nonconvex/Nonconcave Portion
  else{
    switch( options.BREAKPOINT_TYPE ){
      case PolImg<T>::Options::BIN:
      case PolImg<T>::Options::SOS2:{
        const unsigned NKNOTS = itVar1->second->create_subdiv( Op<T>::l(itVar1->second->_range),
          Op<T>::u(itVar1->second->_range) ).size();
        if( NKNOTS > 2 ){
          struct dc{
            static std::pair<double,double> tan1
              ( const double x, const double*rusr, const int*iusr )
              { return std::make_pair( x<0?std::tan(x)-x:0., x<0?mc::sqr(std::tan(x)):0. ); }
            static std::pair<double,double> tan2
              ( const double x, const double*rusr, const int*iusr )
              { return std::make_pair( x>0?std::tan(x)-x:0., x>0?mc::sqr(std::tan(x)):0. ); }
          };
          PolVar<T>* Var3 = _append_aux( Op<T>::tan( itVar1->second->_range )-itVar1->second->_range, true );
          PolVar<T>* Var4 = _append_aux( Op<T>::tan( itVar1->second->_range )-itVar1->second->_range, true );
          _semilinear_cuts( VarR->_var.ops().first, *itVar1->second, Op<T>::l(itVar1->second->_range),
            Op<T>::u(itVar1->second->_range), *Var3, PolCut<T>::GE, dc::tan1, 0, 0 );
          _sandwich_cuts( VarR->_var.ops().first, *itVar1->second, Op<T>::l(itVar1->second->_range),
            Op<T>::u(itVar1->second->_range), *Var3, Op<T>::l(VarR->_range), 0., PolCut<T>::LE, dc::tan1, 0, 0 );
          _semilinear_cuts( VarR->_var.ops().first, *itVar1->second, Op<T>::l(itVar1->second->_range),
            Op<T>::u(itVar1->second->_range), *Var4, PolCut<T>::LE, dc::tan2, 0, 0 );
          _sandwich_cuts( VarR->_var.ops().first, *itVar1->second, Op<T>::l(itVar1->second->_range),
            Op<T>::u(itVar1->second->_range), *Var4, 0., Op<T>::u(VarR->_range), PolCut<T>::GE, dc::tan2, 0, 0 );
          _add_cut( VarR->_var.ops().first, PolCut<T>::EQ, 0., *VarR, 1., *Var3, -1., *Var4, -1.,
              *itVar1->second, -1. );
          break;
        }
        // No break in order to append other "normal" cuts
      }
      case PolImg<T>::Options::NONE: default:{
        struct fct{ static std::pair<double,double> tanfunc
          ( const double x, const double*rusr, const int*iusr )
          { return std::make_pair(
              // f(z) = (z-a)-(tan(z)-tan(a))/(1+tan(z)^2) = 0
              (x-(*rusr))-(std::tan(x)-std::tan(*rusr))/(1.+sqr(std::tan(x))),
              // f'(z) = (tan(z)-tan(a))/(1+tan(z)^2)*2*tan(z)
              2.*std::tan(x)/(1.+sqr(std::tan(x)))*(std::tan(x)-std::tan(*rusr)) ); }
        };
        // Cuts form above
        double xJcc = Op<T>::u(itVar1->second->_range);
        xJcc = _newton( Op<T>::l(itVar1->second->_range), Op<T>::l(itVar1->second->_range), 0.,
          fct::tanfunc, options.ROOT_TOL, options.ROOT_MAXIT, &xJcc, 0 );
        if( mc::isequal( xJcc, Op<T>::l(itVar1->second->_range) ) )
          _semilinear_cuts( VarR->_var.ops().first, *itVar1->second, Op<T>::l(itVar1->second->_range),
            Op<T>::u(itVar1->second->_range), *VarR, PolCut<T>::LE, loc::tan, 0, 0 );
        else
          _sandwich_cuts( VarR->_var.ops().first, *itVar1->second, Op<T>::l(itVar1->second->_range), xJcc,
            *VarR, Op<T>::l(VarR->_range), Op<T>::u(VarR->_range), PolCut<T>::LE, loc::tan, 0, 0 );
        // Cuts from below
        double xJcv = Op<T>::l(itVar1->second->_range);
        xJcv = _newton( Op<T>::u(itVar1->second->_range), 0., Op<T>::u(itVar1->second->_range),
          fct::tanfunc, options.ROOT_TOL, options.ROOT_MAXIT, &xJcv, 0 );
        if( mc::isequal( xJcv, Op<T>::u(itVar1->second->_range) ) )
          _semilinear_cuts( VarR->_var.ops().first, *itVar1->second, Op<T>::l(itVar1->second->_range),
            Op<T>::u(itVar1->second->_range), *VarR, PolCut<T>::GE, loc::tan, 0, 0 );
        else
          _sandwich_cuts( VarR->_var.ops().first, *itVar1->second, xJcv, Op<T>::u(itVar1->second->_range),
            *VarR, Op<T>::l(VarR->_range), Op<T>::u(VarR->_range), PolCut<T>::GE, loc::tan, 0, 0 );
        break;
      }
    }
  }
}

template <typename T>
inline PolVar<T>
acos
( const PolVar<T>&Var1 )
{
  FFGraph* dag = Var1._var.dag();
#ifdef MC__POLIMG_CHECK
  if( !dag || !dag->curOp() )
    throw typename PolImg<T>::Exceptions( PolImg<T>::Exceptions::NOTALLOWED );
#endif

  FFVar* pFFVarR = dag->curOp()->pres;
  PolVar<T>* pVarR = Var1._img->_append_var( pFFVarR, Op<T>::acos( Var1._range ), true );
  return *pVarR;
}

template <typename T> inline void
PolImg<T>::_add_cuts_ACOS
( const PolVar<T>*VarR, FFVar*pVar1 )
{
  if( pVar1->cst() ){
    _add_cut( VarR->_var.ops().first, PolCut<T>::EQ, std::acos( pVar1->num().val() ), *VarR, 1. );
    return;
  }
    
  if( !options.RELAX_NLIN ){
    auto itVar1 = _Vars.find( pVar1 );
    _add_cut( VarR->_var.ops().first, *VarR, *itVar1->second );
    return;
  }

  if( !options.ROOT_USE ) return;

  auto itVar1 = _Vars.find( pVar1 );
  struct loc{
    static std::pair<double,double> acos
      ( const double x, const double*rusr, const int*iusr )
      { return std::make_pair( std::acos(x), -1./std::sqrt(1.-x*x) ); }
  };

  // -- Convex Portion
  if( Op<T>::u(itVar1->second->_range) <= 0. ){
    _semilinear_cuts( VarR->_var.ops().first, *itVar1->second, Op<T>::l(itVar1->second->_range),
      Op<T>::u(itVar1->second->_range), *VarR, PolCut<T>::LE, loc::acos, 0, 0 );
    _sandwich_cuts( VarR->_var.ops().first, *itVar1->second, Op<T>::l(itVar1->second->_range),
      Op<T>::u(itVar1->second->_range), *VarR, Op<T>::l(VarR->_range), Op<T>::u(VarR->_range),
      PolCut<T>::GE, loc::acos, 0, 0 );
  }
  // -- Concave Portion
  else if( Op<T>::l(itVar1->second->_range) >= 0. ){
    _semilinear_cuts( VarR->_var.ops().first, *itVar1->second, Op<T>::l(itVar1->second->_range),
      Op<T>::u(itVar1->second->_range), *VarR, PolCut<T>::GE, loc::acos, 0, 0 );
    _sandwich_cuts( VarR->_var.ops().first, *itVar1->second, Op<T>::l(itVar1->second->_range),
      Op<T>::u(itVar1->second->_range), *VarR, Op<T>::l(VarR->_range), Op<T>::u(VarR->_range),
      PolCut<T>::LE, loc::acos, 0, 0 );
  }
  // -- Nonconvex/Nonconcave Portion
  else{
    switch( options.BREAKPOINT_TYPE ){
      case PolImg<T>::Options::BIN:
      case PolImg<T>::Options::SOS2:{
        const unsigned NKNOTS = itVar1->second->create_subdiv( Op<T>::l(itVar1->second->_range),
          Op<T>::u(itVar1->second->_range) ).size();
        if( NKNOTS > 2 ){
          struct dc{
            static std::pair<double,double> acos1
              ( const double x, const double*rusr, const int*iusr )
              { return std::make_pair( x<0?std::acos(x)+x-PI/4.:PI/4., x<0?1.-1./std::sqrt(1.-x*x):0. ); }
            static std::pair<double,double> acos2
              ( const double x, const double*rusr, const int*iusr )
              { return std::make_pair( x>0?std::acos(x)+x-PI/4.:PI/4., x>0?1.-1./std::sqrt(1.-x*x):0. ); }
          };
          PolVar<T>* Var3 = _append_aux( Op<T>::acos(itVar1->second->_range)+itVar1->second->_range-PI/4., true );
          PolVar<T>* Var4 = _append_aux( Op<T>::acos(itVar1->second->_range)+itVar1->second->_range-PI/4., true );
          _semilinear_cuts( VarR->_var.ops().first, *itVar1->second, Op<T>::l(itVar1->second->_range),
            Op<T>::u(itVar1->second->_range), *Var3, PolCut<T>::LE, dc::acos1, 0, 0 );
          _sandwich_cuts( VarR->_var.ops().first, *itVar1->second, Op<T>::l(itVar1->second->_range),
            Op<T>::u(itVar1->second->_range), *Var3, PI/4., Op<T>::u(Var3->_range), PolCut<T>::GE,
            dc::acos1, 0, 0 );
          _semilinear_cuts( VarR->_var.ops().first, *itVar1->second, Op<T>::l(itVar1->second->_range),
            Op<T>::u(itVar1->second->_range), *Var4, PolCut<T>::GE, dc::acos2, 0, 0 );
          _sandwich_cuts( VarR->_var.ops().first, *itVar1->second, Op<T>::l(itVar1->second->_range),
            Op<T>::u(itVar1->second->_range), *Var4, Op<T>::l(Var4->_range), PI/4., PolCut<T>::LE,
            dc::acos2, 0, 0 );
          _add_cut( VarR->_var.ops().first, PolCut<T>::EQ, 0., *VarR, 1., *Var3, -1., *Var4, -1.,
              *itVar1->second, 1. );
          break;
        }
        // No break in order to append other "normal" cuts
      }
      case PolImg<T>::Options::NONE: default:{
        struct fct{ static double acosfunc
          ( const double x, const double*rusr, const int*iusr )
          { // f(z) = z-a+sqrt(1-z^2)*(acos(z)-acos(a)) = 0
            return x-(*rusr)+std::sqrt(1.-x*x)*(std::acos(x)-std::acos(*rusr)); }
        };
        // Cuts form below
        double xJcv = Op<T>::u(itVar1->second->_range);
        try{
          xJcv = _secant( 0., Op<T>::l(itVar1->second->_range), Op<T>::l(itVar1->second->_range), 0.,
            fct::acosfunc, options.ROOT_TOL, options.ROOT_MAXIT, &xJcv, 0 );
        }
        catch(...){
          xJcv = _goldsect( Op<T>::l(itVar1->second->_range), 0., fct::acosfunc, options.ROOT_TOL,
            options.ROOT_MAXIT, &xJcv, 0 );
        }
        if( mc::isequal( xJcv, Op<T>::l(itVar1->second->_range) ) )
          _semilinear_cuts( VarR->_var.ops().first, *itVar1->second, Op<T>::l(itVar1->second->_range),
            Op<T>::u(itVar1->second->_range), *VarR, PolCut<T>::GE, loc::acos, 0, 0 );
        else
          _sandwich_cuts( VarR->_var.ops().first, *itVar1->second, Op<T>::l(itVar1->second->_range), xJcv,
            *VarR, Op<T>::l(VarR->_range), Op<T>::u(VarR->_range), PolCut<T>::GE, loc::acos, 0, 0 );
        // Cuts from above
        double xJcc = Op<T>::l(itVar1->second->_range);
        try{
          xJcc = _secant( 0., Op<T>::u(itVar1->second->_range), 0., Op<T>::u(itVar1->second->_range),
            fct::acosfunc, options.ROOT_TOL, options.ROOT_MAXIT, &xJcc, 0 );
        }
        catch(...){
          xJcc = _goldsect( 0., Op<T>::u(itVar1->second->_range), fct::acosfunc, options.ROOT_TOL,
            options.ROOT_MAXIT, &xJcc, 0 );
        }
        if( mc::isequal( xJcc, Op<T>::u(itVar1->second->_range) ) )
          _semilinear_cuts( VarR->_var.ops().first, *itVar1->second, Op<T>::l(itVar1->second->_range),
            Op<T>::u(itVar1->second->_range), *VarR, PolCut<T>::LE, loc::acos, 0, 0 );
        else
          _sandwich_cuts( VarR->_var.ops().first, *itVar1->second, xJcc, Op<T>::u(itVar1->second->_range),
            *VarR, Op<T>::l(VarR->_range), Op<T>::u(VarR->_range), PolCut<T>::LE, loc::acos, 0, 0 );
        break;
      }
    }
  }
}

template <typename T>
inline PolVar<T>
asin
( const PolVar<T>&Var1 )
{
  FFGraph* dag = Var1._var.dag();
#ifdef MC__POLIMG_CHECK
  if( !dag || !dag->curOp() )
    throw typename PolImg<T>::Exceptions( PolImg<T>::Exceptions::NOTALLOWED );
#endif

  FFVar* pFFVarR = dag->curOp()->pres;
  PolVar<T>* pVarR = Var1._img->_append_var( pFFVarR, Op<T>::asin( Var1._range ), true );
  return *pVarR;
}

template <typename T> inline void
PolImg<T>::_add_cuts_ASIN
( const PolVar<T>*VarR, FFVar*pVar1 )
{
  if( pVar1->cst() ){
    _add_cut( VarR->_var.ops().first, PolCut<T>::EQ, std::asin( pVar1->num().val() ), *VarR, 1. );
    return;
  }
  
  if( !options.RELAX_NLIN ){
    auto itVar1 = _Vars.find( pVar1 );
    _add_cut( VarR->_var.ops().first, *VarR, *itVar1->second );
    return;
  }

  if( !options.ROOT_USE ) return;

  auto itVar1 = _Vars.find( pVar1 );
  struct loc{
    static std::pair<double,double> asin
      ( const double x, const double*rusr, const int*iusr )
      { return std::make_pair( std::asin(x), 1/std::sqrt(1-x*x) ); }
  };

  // -- Convex Portion
  if( Op<T>::l(itVar1->second->_range) >= 0. ){
    _semilinear_cuts( VarR->_var.ops().first, *itVar1->second, Op<T>::l(itVar1->second->_range),
      Op<T>::u(itVar1->second->_range), *VarR, PolCut<T>::LE, loc::asin, 0, 0 );
    _sandwich_cuts( VarR->_var.ops().first, *itVar1->second, Op<T>::l(itVar1->second->_range),
      Op<T>::u(itVar1->second->_range), *VarR, Op<T>::l(VarR->_range), Op<T>::u(VarR->_range),
      PolCut<T>::GE, loc::asin, 0, 0 );
  }
  // -- Concave Portion
  else if( Op<T>::u(itVar1->second->_range) <= 0. ){
    _semilinear_cuts( VarR->_var.ops().first, *itVar1->second, Op<T>::l(itVar1->second->_range),
      Op<T>::u(itVar1->second->_range), *VarR, PolCut<T>::GE, loc::asin, 0, 0 );
    _sandwich_cuts( VarR->_var.ops().first, *itVar1->second, Op<T>::l(itVar1->second->_range),
      Op<T>::u(itVar1->second->_range), *VarR, Op<T>::l(VarR->_range), Op<T>::u(VarR->_range),
      PolCut<T>::LE, loc::asin, 0, 0 );
  }
  // -- Nonconvex/Nonconcave Portion
  else{
    switch( options.BREAKPOINT_TYPE ){
      case PolImg<T>::Options::BIN:
      case PolImg<T>::Options::SOS2:{
        const unsigned NKNOTS = itVar1->second->create_subdiv( Op<T>::l(itVar1->second->_range),
          Op<T>::u(itVar1->second->_range) ).size();
        if( NKNOTS > 2 ){
          struct dc{
            static std::pair<double,double> asin1
              ( const double x, const double*rusr, const int*iusr )
              { return std::make_pair( x<0?std::asin(x)-x:0., x<0?1./std::sqrt(1.-x*x)-1.:0. ); }
            static std::pair<double,double> asin2
              ( const double x, const double*rusr, const int*iusr )
              { return std::make_pair( x>0?std::asin(x)-x:0., x>0?1./std::sqrt(1.-x*x)-1.:0. ); }
          };
          PolVar<T>* Var3 = _append_aux( Op<T>::asin(itVar1->second->_range)-itVar1->second->_range, true );
          PolVar<T>* Var4 = _append_aux( Op<T>::asin(itVar1->second->_range)-itVar1->second->_range, true );
          _semilinear_cuts( VarR->_var.ops().first, *itVar1->second, Op<T>::l(itVar1->second->_range),
            Op<T>::u(itVar1->second->_range), *Var3, PolCut<T>::GE, dc::asin1, 0, 0 );
          _sandwich_cuts( VarR->_var.ops().first, *itVar1->second, Op<T>::l(itVar1->second->_range),
            Op<T>::u(itVar1->second->_range), *Var3, Op<T>::l(Var3->_range), 0., PolCut<T>::LE,
            dc::asin1, 0, 0 );
          _semilinear_cuts( VarR->_var.ops().first, *itVar1->second, Op<T>::l(itVar1->second->_range),
            Op<T>::u(itVar1->second->_range), *Var4, PolCut<T>::LE, dc::asin2, 0, 0 );
          _sandwich_cuts( VarR->_var.ops().first, *itVar1->second, Op<T>::l(itVar1->second->_range),
            Op<T>::u(itVar1->second->_range), *Var4, 0., Op<T>::u(Var4->_range), PolCut<T>::GE,
            dc::asin2, 0, 0 );
          _add_cut( VarR->_var.ops().first, PolCut<T>::EQ, 0., *VarR, 1., *Var3, -1., *Var4, -1.,
              *itVar1->second, -1. );
          break;
        }
        // No break in order to append other "normal" cuts
      }
      case PolImg<T>::Options::NONE: default:{
        struct fct{ static double asinfunc
          ( const double x, const double*rusr, const int*iusr )
          { // f(z) = z-a-sqrt(1-z^2)*(asin(z)-asin(a)) = 0
            return x-(*rusr)-std::sqrt(1.-x*x)*(std::asin(x)-std::asin(*rusr)); }
        };
        // Cuts form above
        double xJcc = Op<T>::u(itVar1->second->_range);
        try{
          xJcc = _secant( 0., Op<T>::l(itVar1->second->_range), Op<T>::l(itVar1->second->_range), 0.,
            fct::asinfunc, options.ROOT_TOL, options.ROOT_MAXIT, &xJcc, 0 );
        }
        catch(...){
          xJcc = _goldsect( Op<T>::l(itVar1->second->_range), 0., fct::asinfunc, options.ROOT_TOL,
            options.ROOT_MAXIT, &xJcc, 0 );
        }
        if( mc::isequal( xJcc, Op<T>::l(itVar1->second->_range) ) )
          _semilinear_cuts( VarR->_var.ops().first, *itVar1->second, Op<T>::l(itVar1->second->_range),
            Op<T>::u(itVar1->second->_range), *VarR, PolCut<T>::LE, loc::asin, 0, 0 );
        else
          _sandwich_cuts( VarR->_var.ops().first, *itVar1->second, Op<T>::l(itVar1->second->_range), xJcc,
            *VarR, Op<T>::l(VarR->_range), Op<T>::u(VarR->_range), PolCut<T>::LE, loc::asin, 0, 0 );
        // Cuts from below
        double xJcv = Op<T>::l(itVar1->second->_range);
        try{
          xJcv = _secant( 0., Op<T>::u(itVar1->second->_range), 0., Op<T>::u(itVar1->second->_range),
            fct::asinfunc, options.ROOT_TOL, options.ROOT_MAXIT, &xJcv, 0 );
        }
        catch(...){
          xJcv = _goldsect( 0., Op<T>::u(itVar1->second->_range), fct::asinfunc, options.ROOT_TOL,
            options.ROOT_MAXIT, &xJcv, 0 );
        }
        if( mc::isequal( xJcv, Op<T>::u(itVar1->second->_range) ) )
          _semilinear_cuts( VarR->_var.ops().first, *itVar1->second, Op<T>::l(itVar1->second->_range),
            Op<T>::u(itVar1->second->_range), *VarR, PolCut<T>::GE, loc::asin, 0, 0 );
        else
          _sandwich_cuts( VarR->_var.ops().first, *itVar1->second, xJcv, Op<T>::u(itVar1->second->_range),
            *VarR, Op<T>::l(VarR->_range), Op<T>::u(VarR->_range), PolCut<T>::GE, loc::asin, 0, 0 );
        break;
      }
    }
  }
}

template <typename T>
inline PolVar<T>
atan
( const PolVar<T>&Var1 )
{
  FFGraph* dag = Var1._var.dag();
#ifdef MC__POLIMG_CHECK
  if( !dag || !dag->curOp() )
    throw typename PolImg<T>::Exceptions( PolImg<T>::Exceptions::NOTALLOWED );
#endif

  FFVar* pFFVarR = dag->curOp()->pres;
  PolVar<T>* pVarR = Var1._img->_append_var( pFFVarR, Op<T>::atan( Var1._range ), true );
  return *pVarR;
}

template <typename T> inline void
PolImg<T>::_add_cuts_ATAN
( const PolVar<T>*VarR, FFVar*pVar1 )
{
  if( pVar1->cst() ){
    _add_cut( VarR->_var.ops().first, PolCut<T>::EQ, std::atan( pVar1->num().val() ), *VarR, 1. );
    return;
  }
  
  if( !options.RELAX_NLIN ){
    auto itVar1 = _Vars.find( pVar1 );
    _add_cut( VarR->_var.ops().first, *VarR, *itVar1->second );
    return;
  }

  if( !options.ROOT_USE ) return;

  auto itVar1 = _Vars.find( pVar1 );
  struct loc{
    static std::pair<double,double> atan
      ( const double x, const double*rusr, const int*iusr )
      { return std::make_pair( std::atan(x), 1./(1.+x*x) ); }
  };

  // -- Convex Portion
  if( Op<T>::u(itVar1->second->_range) <= 0. ){
    _semilinear_cuts( VarR->_var.ops().first, *itVar1->second, Op<T>::l(itVar1->second->_range),
      Op<T>::u(itVar1->second->_range), *VarR, PolCut<T>::LE, loc::atan, 0, 0 );
    _sandwich_cuts( VarR->_var.ops().first, *itVar1->second, Op<T>::l(itVar1->second->_range),
      Op<T>::u(itVar1->second->_range), *VarR, Op<T>::l(VarR->_range), Op<T>::u(VarR->_range),
      PolCut<T>::GE, loc::atan, 0, 0 );
  }
  // -- Concave Portion
  else if( Op<T>::l(itVar1->second->_range) >= 0. ){
    _semilinear_cuts( VarR->_var.ops().first, *itVar1->second, Op<T>::l(itVar1->second->_range),
      Op<T>::u(itVar1->second->_range), *VarR, PolCut<T>::GE, loc::atan, 0, 0 );
    _sandwich_cuts( VarR->_var.ops().first, *itVar1->second, Op<T>::l(itVar1->second->_range),
      Op<T>::u(itVar1->second->_range), *VarR, Op<T>::l(VarR->_range), Op<T>::u(VarR->_range),
      PolCut<T>::LE, loc::atan, 0, 0 );
  }
  // -- Nonconvex/Nonconcave Portion
  else{
    switch( options.BREAKPOINT_TYPE ){
      case PolImg<T>::Options::BIN:
      case PolImg<T>::Options::SOS2:{
        const unsigned NKNOTS = itVar1->second->create_subdiv( Op<T>::l(itVar1->second->_range),
          Op<T>::u(itVar1->second->_range) ).size();
        if( NKNOTS > 2 ){
          struct dc{
            static std::pair<double,double> atan1
              ( const double x, const double*rusr, const int*iusr )
              { return std::make_pair( x<0?std::atan(x)-x:0., x<0?1./(1.+x*x)-1.:0. ); }
            static std::pair<double,double> atan2
              ( const double x, const double*rusr, const int*iusr )
              { return std::make_pair( x>0?std::atan(x)-x:0., x>0?1./(1.+x*x)-1.:0. ); }
          };
          PolVar<T>* Var3 = _append_aux( Op<T>::atan(itVar1->second->_range)-itVar1->second->_range, true );
          PolVar<T>* Var4 = _append_aux( Op<T>::atan(itVar1->second->_range)-itVar1->second->_range, true );
          _semilinear_cuts( VarR->_var.ops().first, *itVar1->second, Op<T>::l(itVar1->second->_range),
            Op<T>::u(itVar1->second->_range), *Var3, PolCut<T>::LE, dc::atan1, 0, 0 );
          _sandwich_cuts( VarR->_var.ops().first, *itVar1->second, Op<T>::l(itVar1->second->_range),
            Op<T>::u(itVar1->second->_range), *Var3, 0., Op<T>::u(Var3->_range), PolCut<T>::GE,
            dc::atan1, 0, 0 );
          _semilinear_cuts( VarR->_var.ops().first, *itVar1->second, Op<T>::l(itVar1->second->_range),
            Op<T>::u(itVar1->second->_range), *Var4, PolCut<T>::GE, dc::atan2, 0, 0 );
          _sandwich_cuts( VarR->_var.ops().first, *itVar1->second, Op<T>::l(itVar1->second->_range),
            Op<T>::u(itVar1->second->_range), *Var4, Op<T>::l(Var4->_range), 0., PolCut<T>::LE,
            dc::atan2, 0, 0 );
          _add_cut( VarR->_var.ops().first, PolCut<T>::EQ, 0., *VarR, 1., *Var3, -1., *Var4, -1.,
              *itVar1->second, -1. );
          break;
        }
        // No break in order to append other "normal" cuts
      }
      case PolImg<T>::Options::NONE: default:{
        struct fct{ static double atanfunc
          ( const double x, const double*rusr, const int*iusr )
          { // f(z) = z-a-(1+z^2)*(atan(z)-atan(a)) = 0
            return x-(*rusr)-(1.+x*x)*(std::atan(x)-std::atan(*rusr)); }
        };
        // Cuts form below
        double xJcv = Op<T>::u(itVar1->second->_range);
        try{
          xJcv = _secant( 0., Op<T>::l(itVar1->second->_range), Op<T>::l(itVar1->second->_range), 0.,
            fct::atanfunc, options.ROOT_TOL, options.ROOT_MAXIT, &xJcv, 0 );
        }
        catch(...){
          xJcv = _goldsect( Op<T>::l(itVar1->second->_range), 0., fct::atanfunc, options.ROOT_TOL,
            options.ROOT_MAXIT, &xJcv, 0 );
        }
        if( mc::isequal( xJcv, Op<T>::l(itVar1->second->_range) ) )
          _semilinear_cuts( VarR->_var.ops().first, *itVar1->second, Op<T>::l(itVar1->second->_range),
            Op<T>::u(itVar1->second->_range), *VarR, PolCut<T>::GE, loc::atan, 0, 0 );
        else
          _sandwich_cuts( VarR->_var.ops().first, *itVar1->second, Op<T>::l(itVar1->second->_range), xJcv,
            *VarR, Op<T>::l(VarR->_range), Op<T>::u(VarR->_range), PolCut<T>::GE, loc::atan, 0, 0 );
        // Cuts from above
        double xJcc = Op<T>::l(itVar1->second->_range);
        try{
          xJcc = _secant( 0., Op<T>::u(itVar1->second->_range), 0., Op<T>::u(itVar1->second->_range),
            fct::atanfunc, options.ROOT_TOL, options.ROOT_MAXIT, &xJcc, 0 );
        }
        catch(...){
          xJcc = _goldsect( 0., Op<T>::u(itVar1->second->_range), fct::atanfunc, options.ROOT_TOL,
            options.ROOT_MAXIT, &xJcc, 0 );
        }
        if( mc::isequal( xJcc, Op<T>::u(itVar1->second->_range) ) )
          _semilinear_cuts( VarR->_var.ops().first, *itVar1->second, Op<T>::l(itVar1->second->_range),
            Op<T>::u(itVar1->second->_range), *VarR, PolCut<T>::LE, loc::atan, 0, 0 );
        else
          _sandwich_cuts( VarR->_var.ops().first, *itVar1->second, xJcc, Op<T>::u(itVar1->second->_range),
            *VarR, Op<T>::l(VarR->_range), Op<T>::u(VarR->_range), PolCut<T>::LE, loc::atan, 0, 0 );
        break;
      }
    }
  }
}

template <typename T>
inline PolVar<T>
cosh
( const PolVar<T>&Var1 )
{
  FFGraph* dag = Var1._var.dag();
#ifdef MC__POLIMG_CHECK
  if( !dag || !dag->curOp() )
    throw typename PolImg<T>::Exceptions( PolImg<T>::Exceptions::NOTALLOWED );
#endif

  FFVar* pFFVarR = dag->curOp()->pres;
  PolVar<T>* pVarR = Var1._img->_append_var( pFFVarR, Op<T>::cosh( Var1._range ), true );
  return *pVarR;
}

template <typename T> inline void
PolImg<T>::_add_cuts_COSH
( const PolVar<T>*VarR, FFVar*pVar1 )
{
  if( pVar1->cst() ){
    _add_cut( VarR->_var.ops().first, PolCut<T>::EQ, std::cosh( pVar1->num().val() ), *VarR, 1. );
    return;
  }

  auto itVar1 = _Vars.find( pVar1 );
  struct loc{ static std::pair<double,double> cosh
    ( const double x, const double*rusr, const int*iusr )
    { return std::make_pair( std::cosh(x), std::sinh(x) ); }
  };
  _semilinear_cuts( VarR->_var.ops().first, *itVar1->second, Op<T>::l(itVar1->second->_range),
    Op<T>::u(itVar1->second->_range), *VarR, PolCut<T>::LE, loc::cosh );
  _sandwich_cuts( VarR->_var.ops().first, *itVar1->second, Op<T>::l(itVar1->second->_range),
    Op<T>::u(itVar1->second->_range), *VarR, Op<T>::l(VarR->_range), Op<T>::u(VarR->_range),
    PolCut<T>::GE, loc::cosh );
}

template <typename T>
inline PolVar<T>
sinh
( const PolVar<T>&Var1 )
{
  FFGraph* dag = Var1._var.dag();
#ifdef MC__POLIMG_CHECK
  if( !dag || !dag->curOp() )
    throw typename PolImg<T>::Exceptions( PolImg<T>::Exceptions::NOTALLOWED );
#endif

  FFVar* pFFVarR = dag->curOp()->pres;
  PolVar<T>* pVarR = Var1._img->_append_var( pFFVarR, Op<T>::sinh( Var1._range ), true );
  return *pVarR;
}

template <typename T> inline void
PolImg<T>::_add_cuts_SINH
( const PolVar<T>*VarR, FFVar*pVar1 )
{
  if( pVar1->cst() ){
    _add_cut( VarR->_var.ops().first, PolCut<T>::EQ, std::sinh( pVar1->num().val() ), *VarR, 1. );
    return;
  }

  if( !options.ROOT_USE ) return;

  auto itVar1 = _Vars.find( pVar1 );
  struct loc{
    static std::pair<double,double> sinh
      ( const double x, const double*rusr, const int*iusr )
      { return std::make_pair( std::sinh(x), std::cosh(x)); }
  };

  // -- Convex Portion
  if( Op<T>::l(itVar1->second->_range) >= 0. ){
    _semilinear_cuts( VarR->_var.ops().first, *itVar1->second, Op<T>::l(itVar1->second->_range),
      Op<T>::u(itVar1->second->_range), *VarR, PolCut<T>::LE, loc::sinh, 0, 0 );
    _sandwich_cuts( VarR->_var.ops().first, *itVar1->second, Op<T>::l(itVar1->second->_range),
      Op<T>::u(itVar1->second->_range), *VarR, Op<T>::l(VarR->_range), Op<T>::u(VarR->_range),
      PolCut<T>::GE, loc::sinh, 0, 0 );
  }
  // -- Concave Portion
  else if( Op<T>::u(itVar1->second->_range) <= 0. ){
    _semilinear_cuts( VarR->_var.ops().first, *itVar1->second, Op<T>::l(itVar1->second->_range),
      Op<T>::u(itVar1->second->_range), *VarR, PolCut<T>::GE, loc::sinh, 0, 0 );
    _sandwich_cuts( VarR->_var.ops().first, *itVar1->second, Op<T>::l(itVar1->second->_range),
      Op<T>::u(itVar1->second->_range), *VarR, Op<T>::l(VarR->_range), Op<T>::u(VarR->_range),
      PolCut<T>::LE, loc::sinh, 0, 0 );
  }
  // -- Nonconvex/Nonconcave Portion
  else{
    switch( options.BREAKPOINT_TYPE ){
      case PolImg<T>::Options::BIN:
      case PolImg<T>::Options::SOS2:{
        const unsigned NKNOTS = itVar1->second->create_subdiv( Op<T>::l(itVar1->second->_range),
          Op<T>::u(itVar1->second->_range) ).size();
        if( NKNOTS > 2 ){
          struct dc{
            static std::pair<double,double> sinh1
              ( const double x, const double*rusr, const int*iusr )
              { return std::make_pair( x<0?std::sinh(x)-x:0., x<0?std::cosh(x)-1.:0. ); }
            static std::pair<double,double> sinh2
              ( const double x, const double*rusr, const int*iusr )
              { return std::make_pair( x>0?std::sinh(x)-x:0., x>0?std::cosh(x)-1.:0. ); }
          };
          PolVar<T>* Var3 = _append_aux( Op<T>::sinh( itVar1->second->_range ), true );
          PolVar<T>* Var4 = _append_aux( Op<T>::sinh( itVar1->second->_range ), true );
          _semilinear_cuts( VarR->_var.ops().first, *itVar1->second, Op<T>::l(itVar1->second->_range),
            Op<T>::u(itVar1->second->_range), *Var3, PolCut<T>::GE, dc::sinh1, 0, 0 );
          _sandwich_cuts( VarR->_var.ops().first, *itVar1->second, Op<T>::l(itVar1->second->_range),
            Op<T>::u(itVar1->second->_range), *Var3, Op<T>::l(VarR->_range), 0., PolCut<T>::LE, dc::sinh1, 0, 0 );
          _semilinear_cuts( VarR->_var.ops().first, *itVar1->second, Op<T>::l(itVar1->second->_range),
            Op<T>::u(itVar1->second->_range), *Var4, PolCut<T>::LE, dc::sinh2, 0, 0 );
          _sandwich_cuts( VarR->_var.ops().first, *itVar1->second, Op<T>::l(itVar1->second->_range),
            Op<T>::u(itVar1->second->_range), *Var4, 0., Op<T>::u(VarR->_range), PolCut<T>::GE, dc::sinh2, 0, 0 );
          _add_cut( VarR->_var.ops().first, PolCut<T>::EQ, 0., *VarR, 1., *Var3, -1., *Var4, -1.,
              *itVar1->second, -1. );
          break;
        }
        // No break in order to append other "normal" cuts
      }
      case PolImg<T>::Options::NONE: default:{
        struct fct{ static std::pair<double,double> sinhfunc
          ( const double x, const double*rusr, const int*iusr )
          { return std::make_pair(
              // f(z) = (z-a)*cosh(z)-(sinh(z)-sinh(a)) = 0
              (x-(*rusr))*std::cosh(x)-(std::sinh(x)-std::sinh(*rusr)),
              // f'(z) = (z-a)*sinh(z)
              (x-(*rusr))*std::sinh(x) ); }
        };
        // Cuts form above
        double xJcc = Op<T>::u(itVar1->second->_range);
        xJcc = _newton( Op<T>::l(itVar1->second->_range), Op<T>::l(itVar1->second->_range), 0.,
          fct::sinhfunc, options.ROOT_TOL, options.ROOT_MAXIT, &xJcc, 0 );
        if( mc::isequal( xJcc, Op<T>::l(itVar1->second->_range) ) )
          _semilinear_cuts( VarR->_var.ops().first, *itVar1->second, Op<T>::l(itVar1->second->_range),
            Op<T>::u(itVar1->second->_range), *VarR, PolCut<T>::LE, loc::sinh, 0, 0 );
        else
          _sandwich_cuts( VarR->_var.ops().first, *itVar1->second, Op<T>::l(itVar1->second->_range), xJcc,
            *VarR, Op<T>::l(VarR->_range), Op<T>::u(VarR->_range), PolCut<T>::LE, loc::sinh, 0, 0 );
        // Cuts from below
        double xJcv = Op<T>::l(itVar1->second->_range);
        xJcv = _newton( Op<T>::u(itVar1->second->_range), 0., Op<T>::u(itVar1->second->_range),
          fct::sinhfunc, options.ROOT_TOL, options.ROOT_MAXIT, &xJcv, 0 );
        if( mc::isequal( xJcv, Op<T>::u(itVar1->second->_range) ) )
          _semilinear_cuts( VarR->_var.ops().first, *itVar1->second, Op<T>::l(itVar1->second->_range),
            Op<T>::u(itVar1->second->_range), *VarR, PolCut<T>::GE, loc::sinh, 0, 0 );
        else
          _sandwich_cuts( VarR->_var.ops().first, *itVar1->second, xJcv, Op<T>::u(itVar1->second->_range),
            *VarR, Op<T>::l(VarR->_range), Op<T>::u(VarR->_range), PolCut<T>::GE, loc::sinh, 0, 0 );
        break;
      }
    }
  }
}

template <typename T>
inline PolVar<T>
tanh
( const PolVar<T>&Var1 )
{
  FFGraph* dag = Var1._var.dag();
#ifdef MC__POLIMG_CHECK
  if( !dag || !dag->curOp() )
    throw typename PolImg<T>::Exceptions( PolImg<T>::Exceptions::NOTALLOWED );
#endif

  FFVar* pFFVarR = dag->curOp()->pres;
  PolVar<T>* pVarR = Var1._img->_append_var( pFFVarR, Op<T>::tanh( Var1._range ), true );
  return *pVarR;
}

template <typename T> inline void
PolImg<T>::_add_cuts_TANH
( const PolVar<T>*VarR, FFVar*pVar1 )
{
  if( pVar1->cst() ){
    _add_cut( VarR->_var.ops().first, PolCut<T>::EQ, std::tanh( pVar1->num().val() ), *VarR, 1. );
    return;
  }
  
  if( !options.ROOT_USE ) return;

  auto itVar1 = _Vars.find( pVar1 );
  struct loc{
    static std::pair<double,double> tanh
      ( const double x, const double*rusr, const int*iusr )
      { return std::make_pair( std::tanh(x), 1.-sqr(std::tanh(x)) ); }
  };

  // -- Convex Portion
  if( Op<T>::u(itVar1->second->_range) <= 0. ){
    _semilinear_cuts( VarR->_var.ops().first, *itVar1->second, Op<T>::l(itVar1->second->_range),
      Op<T>::u(itVar1->second->_range), *VarR, PolCut<T>::LE, loc::tanh, 0, 0 );
    _sandwich_cuts( VarR->_var.ops().first, *itVar1->second, Op<T>::l(itVar1->second->_range),
      Op<T>::u(itVar1->second->_range), *VarR, Op<T>::l(VarR->_range), Op<T>::u(VarR->_range),
      PolCut<T>::GE, loc::tanh, 0, 0 );
  }
  // -- Concave Portion
  else if( Op<T>::l(itVar1->second->_range) >= 0. ){
    _semilinear_cuts( VarR->_var.ops().first, *itVar1->second, Op<T>::l(itVar1->second->_range),
      Op<T>::u(itVar1->second->_range), *VarR, PolCut<T>::GE, loc::tanh, 0, 0 );
    _sandwich_cuts( VarR->_var.ops().first, *itVar1->second, Op<T>::l(itVar1->second->_range),
      Op<T>::u(itVar1->second->_range), *VarR, Op<T>::l(VarR->_range), Op<T>::u(VarR->_range),
      PolCut<T>::LE, loc::tanh, 0, 0 );
  }
  // -- Nonconvex/Nonconcave Portion
  else{
    switch( options.BREAKPOINT_TYPE ){
      case PolImg<T>::Options::BIN:
      case PolImg<T>::Options::SOS2:{
        const unsigned NKNOTS = itVar1->second->create_subdiv( Op<T>::l(itVar1->second->_range),
          Op<T>::u(itVar1->second->_range) ).size();
        if( NKNOTS > 2 ){
          struct dc{
            static std::pair<double,double> tanh1
              ( const double x, const double*rusr, const int*iusr )
              { return std::make_pair( x<0?std::tanh(x)-x:0., x<0?-sqr(std::tanh(x)):0. ); }
            static std::pair<double,double> tanh2
              ( const double x, const double*rusr, const int*iusr )
              { return std::make_pair( x>0?std::tanh(x)-x:0., x>0?-sqr(std::tanh(x)):0. ); }
          };
          PolVar<T>* Var3 = _append_aux( Op<T>::tanh(itVar1->second->_range)-itVar1->second->_range, true );
          PolVar<T>* Var4 = _append_aux( Op<T>::tanh(itVar1->second->_range)-itVar1->second->_range, true );
          _semilinear_cuts( VarR->_var.ops().first, *itVar1->second, Op<T>::l(itVar1->second->_range),
            Op<T>::u(itVar1->second->_range), *Var3, PolCut<T>::LE, dc::tanh1, 0, 0 );
          _sandwich_cuts( VarR->_var.ops().first, *itVar1->second, Op<T>::l(itVar1->second->_range),
            Op<T>::u(itVar1->second->_range), *Var3, 0., Op<T>::u(Var3->_range), PolCut<T>::GE,
            dc::tanh1, 0, 0 );
          _semilinear_cuts( VarR->_var.ops().first, *itVar1->second, Op<T>::l(itVar1->second->_range),
            Op<T>::u(itVar1->second->_range), *Var4, PolCut<T>::GE, dc::tanh2, 0, 0 );
          _sandwich_cuts( VarR->_var.ops().first, *itVar1->second, Op<T>::l(itVar1->second->_range),
            Op<T>::u(itVar1->second->_range), *Var4, Op<T>::l(Var4->_range), 0., PolCut<T>::LE,
            dc::tanh2, 0, 0 );
          _add_cut( VarR->_var.ops().first, PolCut<T>::EQ, 0., *VarR, 1., *Var3, -1., *Var4, -1.,
              *itVar1->second, -1. );
          break;
        }
        // No break in order to append other "normal" cuts
      }
      case PolImg<T>::Options::NONE: default:{
        struct fct{ static double tanhfunc
          ( const double x, const double*rusr, const int*iusr )
          { // f(z) = (z-a)*(1-tanh(z)^2)-(tanh(z)-tanh(a)) = 0
            return (x-*rusr)*(1-sqr(std::tanh(x)))-(std::tanh(x)-std::tanh(*rusr)); }
        };
        // Cuts form below
        double xJcv = Op<T>::u(itVar1->second->_range);
        try{
          xJcv = _secant( 0., Op<T>::l(itVar1->second->_range), Op<T>::l(itVar1->second->_range), 0.,
            fct::tanhfunc, options.ROOT_TOL, options.ROOT_MAXIT, &xJcv, 0 );
        }
        catch(...){
          xJcv = _goldsect( Op<T>::l(itVar1->second->_range), 0., fct::tanhfunc, options.ROOT_TOL,
            options.ROOT_MAXIT, &xJcv, 0 );
        }
        if( mc::isequal( xJcv, Op<T>::l(itVar1->second->_range) ) )
          _semilinear_cuts( VarR->_var.ops().first, *itVar1->second, Op<T>::l(itVar1->second->_range),
            Op<T>::u(itVar1->second->_range), *VarR, PolCut<T>::GE, loc::tanh, 0, 0 );
        else
          _sandwich_cuts( VarR->_var.ops().first, *itVar1->second, Op<T>::l(itVar1->second->_range), xJcv,
            *VarR, Op<T>::l(VarR->_range), Op<T>::u(VarR->_range), PolCut<T>::GE, loc::tanh, 0, 0 );
        // Cuts from above
        double xJcc = Op<T>::l(itVar1->second->_range);
        try{
          xJcc = _secant( 0., Op<T>::u(itVar1->second->_range), 0., Op<T>::u(itVar1->second->_range),
            fct::tanhfunc, options.ROOT_TOL, options.ROOT_MAXIT, &xJcc, 0 );
        }
        catch(...){
          xJcc = _goldsect( 0., Op<T>::u(itVar1->second->_range), fct::tanhfunc, options.ROOT_TOL,
            options.ROOT_MAXIT, &xJcc, 0 );
        }
        if( mc::isequal( xJcc, Op<T>::u(itVar1->second->_range) ) )
          _semilinear_cuts( VarR->_var.ops().first, *itVar1->second, Op<T>::l(itVar1->second->_range),
            Op<T>::u(itVar1->second->_range), *VarR, PolCut<T>::LE, loc::tanh, 0, 0 );
        else
          _sandwich_cuts( VarR->_var.ops().first, *itVar1->second, xJcc, Op<T>::u(itVar1->second->_range),
            *VarR, Op<T>::l(VarR->_range), Op<T>::u(VarR->_range), PolCut<T>::LE, loc::tanh, 0, 0 );
        break;
      }
    }
  }
}

template <typename T>
inline PolVar<T>
erf
( const PolVar<T>&Var1 )
{
  FFGraph* dag = Var1._var.dag();
#ifdef MC__POLIMG_CHECK
  if( !dag || !dag->curOp() )
    throw typename PolImg<T>::Exceptions( PolImg<T>::Exceptions::NOTALLOWED );
#endif

  FFVar* pFFVarR = dag->curOp()->pres;
  PolVar<T>* pVarR = Var1._img->_append_var( pFFVarR, Op<T>::erf( Var1._range ), true );
  return *pVarR;
}

template <typename T> inline void
PolImg<T>::_add_cuts_ERF
( const PolVar<T>*VarR, FFVar*pVar1 )
{
  if( pVar1->cst() ){
    _add_cut( VarR->_var.ops().first, PolCut<T>::EQ, std::erf( pVar1->num().val() ), *VarR, 1. );
    return;
  }
  
  if( !options.ROOT_USE ) return;

  auto itVar1 = _Vars.find( pVar1 );
  struct loc{
    static std::pair<double,double> erf
      ( const double x, const double*rusr, const int*iusr )
      { return std::make_pair( std::erf(x), 2./std::sqrt(PI)*std::exp(-sqr(x)) ); }
  };

  // -- Convex Portion
  if( Op<T>::u(itVar1->second->_range) <= 0. ){
    _semilinear_cuts( VarR->_var.ops().first, *itVar1->second, Op<T>::l(itVar1->second->_range),
      Op<T>::u(itVar1->second->_range), *VarR, PolCut<T>::LE, loc::erf, 0, 0 );
    _sandwich_cuts( VarR->_var.ops().first, *itVar1->second, Op<T>::l(itVar1->second->_range),
      Op<T>::u(itVar1->second->_range), *VarR, Op<T>::l(VarR->_range), Op<T>::u(VarR->_range),
      PolCut<T>::GE, loc::erf, 0, 0 );
  }
  // -- Concave Portion
  else if( Op<T>::l(itVar1->second->_range) >= 0. ){
    _semilinear_cuts( VarR->_var.ops().first, *itVar1->second, Op<T>::l(itVar1->second->_range),
      Op<T>::u(itVar1->second->_range), *VarR, PolCut<T>::GE, loc::erf, 0, 0 );
    _sandwich_cuts( VarR->_var.ops().first, *itVar1->second, Op<T>::l(itVar1->second->_range),
      Op<T>::u(itVar1->second->_range), *VarR, Op<T>::l(VarR->_range), Op<T>::u(VarR->_range),
      PolCut<T>::LE, loc::erf, 0, 0 );
  }
  // -- Nonconvex/Nonconcave Portion
  else{
    switch( options.BREAKPOINT_TYPE ){
      case PolImg<T>::Options::BIN:
      case PolImg<T>::Options::SOS2:{
        const unsigned NKNOTS = itVar1->second->create_subdiv( Op<T>::l(itVar1->second->_range),
          Op<T>::u(itVar1->second->_range) ).size();
        if( NKNOTS > 2 ){
          struct dc{
            static std::pair<double,double> erf1
              ( const double x, const double*rusr, const int*iusr )
              { return std::make_pair( x<0?std::erf(x)-2./std::sqrt(PI)*x:0., x<0?2./std::sqrt(PI)*(std::exp(-sqr(x))-1.):0. ); }
            static std::pair<double,double> erf2
              ( const double x, const double*rusr, const int*iusr )
              { return std::make_pair( x>0?std::erf(x)-2./std::sqrt(PI)*x:0., x>0?2./std::sqrt(PI)*(std::exp(-sqr(x))-1.):0. ); }
          };
          PolVar<T>* Var3 = _append_aux( Op<T>::erf(itVar1->second->_range)-2./std::sqrt(PI)*itVar1->second->_range, true );
          PolVar<T>* Var4 = _append_aux( Op<T>::erf(itVar1->second->_range)-2./std::sqrt(PI)*itVar1->second->_range, true );
          _semilinear_cuts( VarR->_var.ops().first, *itVar1->second, Op<T>::l(itVar1->second->_range),
            Op<T>::u(itVar1->second->_range), *Var3, PolCut<T>::LE, dc::erf1, 0, 0 );
          _sandwich_cuts( VarR->_var.ops().first, *itVar1->second, Op<T>::l(itVar1->second->_range),
            Op<T>::u(itVar1->second->_range), *Var3, 0., Op<T>::u(Var3->_range), PolCut<T>::GE,
            dc::erf1, 0, 0 );
          _semilinear_cuts( VarR->_var.ops().first, *itVar1->second, Op<T>::l(itVar1->second->_range),
            Op<T>::u(itVar1->second->_range), *Var4, PolCut<T>::GE, dc::erf2, 0, 0 );
          _sandwich_cuts( VarR->_var.ops().first, *itVar1->second, Op<T>::l(itVar1->second->_range),
            Op<T>::u(itVar1->second->_range), *Var4, Op<T>::l(Var4->_range), 0., PolCut<T>::LE,
            dc::erf2, 0, 0 );
          _add_cut( VarR->_var.ops().first, PolCut<T>::EQ, 0., *VarR, 1., *Var3, -1., *Var4, -1.,
              *itVar1->second, -2./std::sqrt(PI) );
          break;
        }
        // No break in order to append other "normal" cuts
      }
      case PolImg<T>::Options::NONE: default:{
        struct fct{ static double erffunc
          ( const double x, const double*rusr, const int*iusr )
          { // f(z) = (z-a)*exp(-z^2)-sqrt(pi)/2.*(erf(z)-erf(a)) = 0
            return (x-*rusr)*std::exp(-x*x)-std::sqrt(PI)/2.*(std::erf(x)-std::erf(*rusr)); }
        };
        // Cuts form below
        double xJcv = Op<T>::u(itVar1->second->_range);
        try{
          xJcv = _secant( 0., Op<T>::l(itVar1->second->_range), Op<T>::l(itVar1->second->_range), 0.,
            fct::erffunc, options.ROOT_TOL, options.ROOT_MAXIT, &xJcv, 0 );
        }
        catch(...){
          xJcv = _goldsect( Op<T>::l(itVar1->second->_range), 0., fct::erffunc, options.ROOT_TOL,
            options.ROOT_MAXIT, &xJcv, 0 );
        }
        if( mc::isequal( xJcv, Op<T>::l(itVar1->second->_range) ) )
          _semilinear_cuts( VarR->_var.ops().first, *itVar1->second, Op<T>::l(itVar1->second->_range),
            Op<T>::u(itVar1->second->_range), *VarR, PolCut<T>::GE, loc::erf, 0, 0 );
        else
          _sandwich_cuts( VarR->_var.ops().first, *itVar1->second, Op<T>::l(itVar1->second->_range), xJcv,
            *VarR, Op<T>::l(VarR->_range), Op<T>::u(VarR->_range), PolCut<T>::GE, loc::erf, 0, 0 );
        // Cuts from above
        double xJcc = Op<T>::l(itVar1->second->_range);
        try{
          xJcc = _secant( 0., Op<T>::u(itVar1->second->_range), 0., Op<T>::u(itVar1->second->_range),
            fct::erffunc, options.ROOT_TOL, options.ROOT_MAXIT, &xJcc, 0 );
        }
        catch(...){
          xJcc = _goldsect( 0., Op<T>::u(itVar1->second->_range), fct::erffunc, options.ROOT_TOL,
            options.ROOT_MAXIT, &xJcc, 0 );
        }
        if( mc::isequal( xJcc, Op<T>::u(itVar1->second->_range) ) )
          _semilinear_cuts( VarR->_var.ops().first, *itVar1->second, Op<T>::l(itVar1->second->_range),
            Op<T>::u(itVar1->second->_range), *VarR, PolCut<T>::LE, loc::erf, 0, 0 );
        else
          _sandwich_cuts( VarR->_var.ops().first, *itVar1->second, xJcc, Op<T>::u(itVar1->second->_range),
            *VarR, Op<T>::l(VarR->_range), Op<T>::u(VarR->_range), PolCut<T>::LE, loc::erf, 0, 0 );
        break;
      }
    }
  }
}

template <typename T>
inline PolVar<T>
fabs
( const PolVar<T>&Var1 )
{
  FFGraph* dag = Var1._var.dag();
#ifdef MC__POLIMG_CHECK
  if( !dag || !dag->curOp() )
    throw typename PolImg<T>::Exceptions( PolImg<T>::Exceptions::NOTALLOWED );
#endif

  FFVar* pFFVarR = dag->curOp()->pres;
  PolVar<T>* pVarR = Var1._img->_append_var( pFFVarR, Op<T>::fabs( Var1._range ), true );
  return *pVarR;
}

template <typename T> inline void
PolImg<T>::_add_cuts_FABS
( const PolVar<T>*VarR, FFVar*pVar1 )
{
  // Constant operand
  if( pVar1->cst() ){
    _add_cut( VarR->_var.ops().first, PolCut<T>::EQ, std::fabs( pVar1->num().val() ), *VarR, 1. );
    return;
  }

  // No relaxation
  if( !options.RELAX_DISC ){
    auto itVar1 = _Vars.find( pVar1 );
    _add_cut( VarR->_var.ops().first, *VarR, *itVar1->second );
    return;
  }

  auto itVar1 = _Vars.find( pVar1 );

  // Positive branch only
  if( Op<T>::l(itVar1->second->_range) >= 0. ){
    _add_cut( VarR->_var.ops().first, PolCut<T>::EQ, 0., *VarR, 1., *itVar1->second, -1. );
    return;
  }

  // Negative branch only
  if( Op<T>::u(itVar1->second->_range) <= 0. ){
    _add_cut( VarR->_var.ops().first, PolCut<T>::EQ, 0., *VarR, 1., *itVar1->second,  1. );
    return;
  }

  // Linear overestimator
  if( options.RELAX_DISC == 1 ){
    double XL = Op<T>::l(itVar1->second->_range),  dX = Op<T>::diam(itVar1->second->_range),
           YL = std::fabs(Op<T>::l(itVar1->second->_range)), dY = std::fabs(Op<T>::u(itVar1->second->_range))-YL;
    _add_cut( VarR->_var.ops().first, PolCut<T>::GE, dY*XL-dX*YL, *VarR, -dX, *itVar1->second, dY );
  }

  // Piecewise-linear overestimator
  else{
    const double M1 =  2*Op<T>::u(itVar1->second->_range);
    const double M2 = -2*Op<T>::l(itVar1->second->_range);
    PolVar<T>* VarB = _append_aux( Op<T>::zeroone(), false );
    _add_cut( VarR->_var.ops().first, PolCut<T>::LE, M2, *VarR, 1., *itVar1->second, -1., *VarB,  M2 );
    _add_cut( VarR->_var.ops().first, PolCut<T>::LE, 0., *VarR, 1., *itVar1->second,  1., *VarB, -M1 );
  }
  
  // Linear underestimators
  _add_cut( VarR->_var.ops().first, PolCut<T>::GE, 0., *VarR, 1., *itVar1->second,  -1. );
  _add_cut( VarR->_var.ops().first, PolCut<T>::GE, 0., *VarR, 1., *itVar1->second,   1. );
}

template <typename T> inline bool
PolImg<T>::_add_LQ_FABS
( PolLQExpr<T>*&pLQ, const PolVar<T>*VarR, FFVar*pVar1 )
{
  auto itVar1 = _Vars.find( pVar1 );
  if( !pVar1->cst() || Op<T>::l(itVar1->second->_range)*Op<T>::u(itVar1->second->_range) < 0. ) return false;
  if( !pLQ ) pLQ = _append_LQ( VarR );
  if( pVar1->cst() )
    pLQ->substitute( VarR, std::fabs(pVar1->num().val()) );
  else{
    if( Op<T>::l(itVar1->second->_range) >= 0 )
      pLQ->substitute( VarR, 1., itVar1->second );
    else
      pLQ->substitute( VarR, -1., itVar1->second );
  }
  return true;
}

template <typename T>
inline PolVar<T>
fstep
( const PolVar<T>&Var1 )
{
  FFGraph* dag = Var1._var.dag();
#ifdef MC__POLIMG_CHECK
  if( !dag || !dag->curOp() )
    throw typename PolImg<T>::Exceptions( PolImg<T>::Exceptions::NOTALLOWED );
#endif

  FFVar* pFFVarR = dag->curOp()->pres;
  PolVar<T>* pVarR = Var1._img->_append_var( pFFVarR, Op<T>::fstep( Var1._range ), true );
  return *pVarR;
}

template <typename T> inline void
PolImg<T>::_add_cuts_FSTEP
( const PolVar<T>*VarR, FFVar*pVar1 )
{
  // Constant operand
  if( pVar1->cst() ){
    _add_cut( VarR->_var.ops().first, PolCut<T>::EQ, mc::fstep( pVar1->num().val() ), *VarR, 1. );
    return;
  }

  // No relaxation
  if( !options.RELAX_DISC ){
    auto itVar1 = _Vars.find( pVar1 );
    _add_cut( VarR->_var.ops().first, *VarR, *itVar1->second );
    return;
  }

  auto itVar1 = _Vars.find( pVar1 );
  
  // Positive branch only
  if( Op<T>::l(itVar1->second->_range) >= 0. ){
    _add_cut( VarR->_var.ops().first, PolCut<T>::EQ, mc::fstep( Op<T>::l(itVar1->second->_range) ), *VarR, 1. );
  }

  // Negative branch only
  else if( Op<T>::u(itVar1->second->_range) < 0. ){
    _add_cut( VarR->_var.ops().first, PolCut<T>::EQ, mc::fstep( Op<T>::u(itVar1->second->_range) ), *VarR, 1. );
  }
  
  // Linear overestimator
  else if( options.RELAX_DISC == 1 ){
    _add_cut( VarR->_var.ops().first, PolCut<T>::GE, 0., *VarR, Op<T>::u(itVar1->second->_range),
                 *itVar1->second,  -1. );
    _add_cut( VarR->_var.ops().first, PolCut<T>::GE, Op<T>::l(itVar1->second->_range),
                 *VarR, Op<T>::l(itVar1->second->_range), *itVar1->second,  1. );
  }

  // Piecewise-linear overestimator
  else{
    PolVar<T>* VarB = _append_aux( Op<T>::zeroone(), false );
    _add_cut( VarR->_var.ops().first, PolCut<T>::EQ, 0., *VarR, 1., *VarB, -1. );
    _add_cut( VarR->_var.ops().first, PolCut<T>::GE, 0., *VarR, Op<T>::u(itVar1->second->_range),
                 *itVar1->second,  -1., *VarB, Op<T>::u(itVar1->second->_range) );
    _add_cut( VarR->_var.ops().first, PolCut<T>::GE, Op<T>::l(itVar1->second->_range), 
                 *itVar1->second,  1., *VarB, Op<T>::l(itVar1->second->_range) );
  }
}

template <typename T> inline bool
PolImg<T>::_add_LQ_FSTEP
( PolLQExpr<T>*&pLQ, const PolVar<T>*VarR, FFVar*pVar1 )
{
  auto itVar1 = _Vars.find( pVar1 );
  if( !pVar1->cst() || (Op<T>::l(itVar1->second->_range) < 0. && Op<T>::u(itVar1->second->_range) >= 0.) ) return false;
  if( !pLQ ) pLQ = _append_LQ( VarR );
  if( pVar1->cst() )
    pLQ->substitute( VarR, mc::fstep(pVar1->num().val()) );
  else
    pLQ->substitute( VarR, mc::fstep(Op<T>::l(itVar1->second->_range)) );
  return true;
}

template <typename T>
inline PolVar<T>
min
( const PolVar<T>&Var1, const PolVar<T>&Var2 )
{
  if( Var1._img && Var2._img && Var1._img != Var2._img )
    throw typename PolImg<T>::Exceptions( PolImg<T>::Exceptions::ENVMIS );
  PolImg<T>* img = Var1._img? Var1._img: Var2._img;
  FFGraph* dag = Var1._var.cst()? Var2._var.dag(): Var1._var.dag();
  if( !dag || !dag->curOp() )
    throw typename PolImg<T>::Exceptions( PolImg<T>::Exceptions::NOTALLOWED );
  FFVar* pFFVarR = dag->curOp()->pres;
  PolVar<T>* pVarR = img->_append_var( pFFVarR, Op<T>::min( Var1.range(), Var2.range() ), true );
  return *pVarR;
}

template <typename T> inline void
PolImg<T>::_add_cuts_MINF
( const PolVar<T>*VarR, FFVar*pVar1, FFVar*pVar2 )
{
  // Both operands constant
  if( pVar1->cst() && pVar2->cst() ){
    _add_cut( VarR->_var.ops().first, PolCut<T>::EQ, std::min( pVar1->num().val(), pVar2->num().val() ), *VarR, 1. );
    return;
  }

  // No relaxation
  it_Vars itVar1, itVar2;
  if( !options.RELAX_DISC ){
    if( pVar1->cst() ){
      itVar2 = _Vars.find( pVar2 );
      _add_cut( VarR->_var.ops().first, *VarR, *itVar2->second, pVar1->num().val() );
    }
    else if( pVar2->cst() ){
      itVar1 = _Vars.find( pVar1 );
      _add_cut( VarR->_var.ops().first, *VarR, *itVar1->second, pVar2->num().val() );
    }
    else{
      itVar1 = _Vars.find( pVar1 );
      itVar2 = _Vars.find( pVar2 );
      _add_cut( VarR->_var.ops().first, *VarR, *itVar1->second, *itVar2->second );
    }
    return;
  }

  // Set big-M values
  double M1 = 0., M2 = 0.;
  if( pVar1->cst() ){
    itVar2 = _Vars.find( pVar2 );
    M1 = pVar1->num().val() - Op<T>::l(itVar2->second->_range);
    M2 = Op<T>::u(itVar2->second->_range) - pVar1->num().val();
  }
  else if( pVar2->cst() ){
    itVar1 = _Vars.find( pVar1 );
    M1 = Op<T>::u(itVar1->second->_range) - pVar2->num().val();
    M2 = pVar2->num().val() - Op<T>::l(itVar1->second->_range);
  }
  else{
    itVar1 = _Vars.find( pVar1 );
    itVar2 = _Vars.find( pVar2 );
    M1 = Op<T>::u(itVar1->second->_range) - Op<T>::l(itVar2->second->_range);
    M2 = Op<T>::u(itVar2->second->_range) - Op<T>::l(itVar1->second->_range);
  }

  // Special cases
  if( M1 <= 0 ){
    if( pVar1->cst() )
      _add_cut( VarR->_var.ops().first, PolCut<T>::EQ, pVar1->num().val(), *VarR, 1. );    
    else
      _add_cut( VarR->_var.ops().first, PolCut<T>::EQ, 0., *VarR, 1., *itVar1->second, -1. );
    return;
  }
  if( M2 <= 0 ){
    if( pVar2->cst() )
      _add_cut( VarR->_var.ops().first, PolCut<T>::EQ, pVar2->num().val(), *VarR, 1. );    
    else
      _add_cut( VarR->_var.ops().first, PolCut<T>::EQ, 0., *VarR, 1., *itVar2->second, -1. );
    return;
  }

  // Linear underestimator and piecewise-linear overestimator
  PolVar<T>* VarB = _append_aux( Op<T>::zeroone(), false );
  if( pVar1->cst() ){
    _add_cut( VarR->_var.ops().first, PolCut<T>::LE, pVar1->num().val(),    *VarR, 1. );
    _add_cut( VarR->_var.ops().first, PolCut<T>::GE, pVar1->num().val()-M1, *VarR, 1., *VarB, -M1 );
  }
  else{
    _add_cut( VarR->_var.ops().first, PolCut<T>::LE, 0.,  *VarR, 1., *itVar1->second, -1. );
    _add_cut( VarR->_var.ops().first, PolCut<T>::GE, -M1, *VarR, 1., *itVar1->second, -1., *VarB, -M1 );
  }
  if( pVar2->cst() ){
    _add_cut( VarR->_var.ops().first, PolCut<T>::LE, pVar2->num().val(), *VarR, 1. );
    _add_cut( VarR->_var.ops().first, PolCut<T>::GE, pVar2->num().val(), *VarR, 1., *VarB, M2 );
  }
  else{
    _add_cut( VarR->_var.ops().first, PolCut<T>::LE, 0., *VarR, 1., *itVar2->second, -1. );
    _add_cut( VarR->_var.ops().first, PolCut<T>::GE, 0., *VarR, 1., *itVar2->second, -1., *VarB,  M2 );
  }
}

template <typename T> inline bool
PolImg<T>::_add_LQ_MINF
( PolLQExpr<T>*&pLQ, const PolVar<T>*VarR, FFVar*pVar1, FFVar*pVar2 )
{
    if( pVar1->cst() && pVar2->cst() ){
    pLQ->substitute( VarR, std::min( pVar1->num().val(), pVar2->num().val() ) );    
    return true;
  }
  
  if( pVar1->cst() ){
    auto itVar2 = _Vars.find( pVar2 );
    if( pVar1->num().val() < Op<T>::l(itVar2->second->_range) ){
      pLQ->substitute( VarR, pVar1->num().val() );
      return true;
    }
    else if( pVar1->num().val() > Op<T>::u(itVar2->second->_range) ){
      pLQ->substitute( VarR, 1., itVar2->second );
      return true;
    }
    return false;
  }
  
  if( pVar2->cst() ){
    auto itVar1 = _Vars.find( pVar1 );
    if( pVar2->num().val() < Op<T>::l(itVar1->second->_range) ){
      pLQ->substitute( VarR, pVar2->num().val() );
      return true;
    }
    else if( pVar2->num().val() > Op<T>::u(itVar1->second->_range) ){
      pLQ->substitute( VarR, 1., itVar1->second );
      return true;
    }
    return false;
  }

  auto itVar1 = _Vars.find( pVar1 );
  auto itVar2 = _Vars.find( pVar2 );
  if( Op<T>::u(itVar1->second->_range) < Op<T>::l(itVar2->second->_range) ){
    pLQ->substitute( VarR, 1., itVar1->second );
    return true;
  }
  else if( Op<T>::l(itVar1->second->_range) > Op<T>::u(itVar2->second->_range) ){
    pLQ->substitute( VarR, 1., itVar2->second );
    return true;
  }
  return false;
}

template <typename T>
inline PolVar<T>
max
( const PolVar<T>&Var1, const PolVar<T>&Var2 )
{
  if( Var1._img && Var2._img && Var1._img != Var2._img )
    throw typename PolImg<T>::Exceptions( PolImg<T>::Exceptions::ENVMIS );
  PolImg<T>* img = Var1._img? Var1._img: Var2._img;
  FFGraph* dag = Var1._var.cst()? Var2._var.dag(): Var1._var.dag();
  if( !dag || !dag->curOp() )
    throw typename PolImg<T>::Exceptions( PolImg<T>::Exceptions::NOTALLOWED );
  FFVar* pFFVarR = dag->curOp()->pres;
  PolVar<T>* pVarR = img->_append_var( pFFVarR, Op<T>::max( Var1.range(), Var2.range() ), true );
  return *pVarR;
}

template <typename T> inline void
PolImg<T>::_add_cuts_MAXF
( const PolVar<T>*VarR, FFVar*pVar1, FFVar*pVar2 )
{
  // Both operands constant
  if( pVar1->cst() && pVar2->cst() ){
    _add_cut( VarR->_var.ops().first, PolCut<T>::EQ, std::max( pVar1->num().val(), pVar2->num().val() ), *VarR, 1. );
    return;
  }

  // No relaxation
  it_Vars itVar1, itVar2;
  if( !options.RELAX_DISC ){
    if( pVar1->cst() ){
      itVar2 = _Vars.find( pVar2 );
      _add_cut( VarR->_var.ops().first, *VarR, *itVar2->second, pVar1->num().val() );
    }
    else if( pVar2->cst() ){
      itVar1 = _Vars.find( pVar1 );
      _add_cut( VarR->_var.ops().first, *VarR, *itVar1->second, pVar2->num().val() );
    }
    else{
      itVar1 = _Vars.find( pVar1 );
      itVar2 = _Vars.find( pVar2 );
      _add_cut( VarR->_var.ops().first, *VarR, *itVar1->second, *itVar2->second );
    }
    return;
  }

  // Set big-M values
  double M1 = 0., M2 = 0.;
  if( pVar1->cst() ){
    itVar2 = _Vars.find( pVar2 );
    M1 = pVar1->num().val() - Op<T>::l(itVar2->second->_range);
    M2 = Op<T>::u(itVar2->second->_range) - pVar1->num().val();
  }
  else if( pVar2->cst() ){
    itVar1 = _Vars.find( pVar1 );
    M1 = Op<T>::u(itVar1->second->_range) - pVar2->num().val();
    M2 = pVar2->num().val() - Op<T>::l(itVar1->second->_range);
  }
  else{
    itVar1 = _Vars.find( pVar1 );
    itVar2 = _Vars.find( pVar2 );
    M1 = Op<T>::u(itVar1->second->_range) - Op<T>::l(itVar2->second->_range);
    M2 = Op<T>::u(itVar2->second->_range) - Op<T>::l(itVar1->second->_range);
  }

  // Special cases
  if( M1 <= 0 ){
    if( pVar2->cst() )
      _add_cut( VarR->_var.ops().first, PolCut<T>::EQ, pVar2->num().val(), *VarR, 1. );    
    else
      _add_cut( VarR->_var.ops().first, PolCut<T>::EQ, 0., *VarR, 1., *itVar2->second, -1. );
    return;
  }
  if( M2 <= 0 ){
    if( pVar1->cst() )
      _add_cut( VarR->_var.ops().first, PolCut<T>::EQ, pVar1->num().val(), *VarR, 1. );    
    else
      _add_cut( VarR->_var.ops().first, PolCut<T>::EQ, 0., *VarR, 1., *itVar1->second, -1. );
    return;
  }

  // Linear underestimator and piecewise-linear overestimator
  PolVar<T>* VarB = _append_aux( Op<T>::zeroone(), false );
  if( pVar1->cst() ){
    _add_cut( VarR->_var.ops().first, PolCut<T>::GE, pVar1->num().val(),    *VarR, 1. );
    _add_cut( VarR->_var.ops().first, PolCut<T>::LE, pVar1->num().val()-M2, *VarR, 1., *VarB, M2 );
  }
  else{
    _add_cut( VarR->_var.ops().first, PolCut<T>::GE, 0.,  *VarR, 1., *itVar1->second, -1. );
    _add_cut( VarR->_var.ops().first, PolCut<T>::LE, M2, *VarR, 1., *itVar1->second, -1., *VarB, M2 );
  }
  if( pVar2->cst() ){
    _add_cut( VarR->_var.ops().first, PolCut<T>::GE, pVar2->num().val(), *VarR, 1. );
    _add_cut( VarR->_var.ops().first, PolCut<T>::LE, pVar2->num().val(), *VarR, 1., *VarB, -M1 );
  }
  else{
    _add_cut( VarR->_var.ops().first, PolCut<T>::GE, 0., *VarR, 1., *itVar2->second, -1. );
    _add_cut( VarR->_var.ops().first, PolCut<T>::LE, 0., *VarR, 1., *itVar2->second, -1., *VarB,  -M1 );
  }
}

template <typename T> inline bool
PolImg<T>::_add_LQ_MAXF
( PolLQExpr<T>*&pLQ, const PolVar<T>*VarR, FFVar*pVar1, FFVar*pVar2 )
{
  if( pVar1->cst() && pVar2->cst() ){
    pLQ->substitute( VarR, std::max( pVar1->num().val(), pVar2->num().val() ) );    
    return true;
  }
  
  if( pVar1->cst() ){
    auto itVar2 = _Vars.find( pVar2 );
    if( pVar1->num().val() > Op<T>::u(itVar2->second->_range) ){
      pLQ->substitute( VarR, pVar1->num().val() );
      return true;
    }
    else if( pVar1->num().val() < Op<T>::l(itVar2->second->_range) ){
      pLQ->substitute( VarR, 1., itVar2->second );
      return true;
    }
    return false;
  }
  
  if( pVar2->cst() ){
    auto itVar1 = _Vars.find( pVar1 );
    if( pVar2->num().val() > Op<T>::u(itVar1->second->_range) ){
      pLQ->substitute( VarR, pVar2->num().val() );
      return true;
    }
    else if( pVar2->num().val() < Op<T>::l(itVar1->second->_range) ){
      pLQ->substitute( VarR, 1., itVar1->second );
      return true;
    }
    return false;
  }

  auto itVar1 = _Vars.find( pVar1 );
  auto itVar2 = _Vars.find( pVar2 );
  if( Op<T>::l(itVar1->second->_range) > Op<T>::u(itVar2->second->_range) ){
    pLQ->substitute( VarR, 1., itVar1->second );
    return true;
  }
  else if( Op<T>::u(itVar1->second->_range) < Op<T>::l(itVar2->second->_range) ){
    pLQ->substitute( VarR, 1., itVar2->second );
    return true;
  }  
  return false;
}

} // namespace mc

//#include "mcop.hpp"

namespace mc
{

//! @brief Specialization of the mc::Op templated structure for use of type mc::PolVar as template type in other MC++ classes
template< typename T > struct Op< mc::PolVar<T> >
{
  typedef mc::PolVar<T> PV;
  static PV point( const double c ) { return PV(c); }
  static PV zeroone(){ throw std::runtime_error("operation not permitted"); }// { return PV( mc::Op<T>::zeroone() ); }
  static void I(PV& x, const PV&y) { x = y; }
  static double l(const PV& x) { return mc::Op<T>::l(x.range()); }
  static double u(const PV& x) { return mc::Op<T>::u(x.range()); }
  static double abs (const PV& x) { return mc::Op<T>::abs(x.range());  }
  static double mid (const PV& x) { return mc::Op<T>::mid(x.range());  }
  static double diam(const PV& x) { return mc::Op<T>::diam(x.range()); }
  static PV inv (const PV& x){ return mc::inv(x);  }
  static PV sqr (const PV& x){ return mc::sqr(x);  }
  static PV sqrt(const PV& x){ return mc::sqrt(x); }
  static PV log (const PV& x){ return mc::log(x);  }
  static PV xlog(const PV& x){ return mc::xlog(x); }
  static PV lmtd(const PV& x, const PV& y)
    { throw std::runtime_error("operation not permitted"); }
//  { return mc::lmtd(x,y); }
  static PV rlmtd(const PV& x, const PV& y)
    { throw std::runtime_error("operation not permitted"); }
//  { return mc::rlmtd(x,y); }
  static PV fabs(const PV& x){ return mc::fabs(x); }
  static PV exp (const PV& x){ return mc::exp(x);  }
  static PV cos (const PV& x){ return mc::cos(x);  }
  static PV sin (const PV& x){ return mc::sin(x);  }
  static PV tan (const PV& x){ return mc::tan(x);  }
  static PV acos(const PV& x){ return mc::acos(x); }
  static PV asin(const PV& x){ return mc::asin(x); }
  static PV atan(const PV& x){ return mc::atan(x); }
  static PV cosh(const PV& x){ return mc::cosh(x); }
  static PV sinh(const PV& x){ return mc::sinh(x); }
  static PV tanh(const PV& x){ return mc::tanh(x); }
  static PV erf (const PV& x){ return mc::erf(x);  }
  static PV erfc(const PV& x){ return 1-mc::erf(x); }
  static PV fstep(const PV& x){ return mc::fstep(x); }
  static PV bstep(const PV& x){ return mc::fstep(-x); }
  static PV min (const PV& x, const PV& y)
    { return mc::min(x,y);  }
  static PV max (const PV& x, const PV& y)
    { return mc::max(x,y);  }
  static PV arh (const PV& x, const double k){ return mc::exp(-k/x); }
  template <typename EXP> static PV pow(const PV& x, const EXP& y) { return mc::pow(x,y); }
  static PV cheb (const PV& x, const unsigned n) { return mc::cheb(x,n); }
  static PV prod (const unsigned int n, const PV* x) { return mc::prod(n,x); }
  static PV monom (const unsigned int n, const PV* x, const unsigned* k)
    { throw std::runtime_error("operation not permitted"); }
//  { return mc::monom(n,x,k); }
  static PV hull(const PV& x, const PV& y){ return mc::Op<T>::hull(x.range(),y.range()); }
  static bool inter(PV& xIy, const PV& x, const PV& y){
    try{ xIy = x^y; return true; }
    catch(typename mc::PolImg<T>::Exceptions &eObj){ return false; }
  }
  static bool eq(const PV& x, const PV& y)
    { return mc::Op<T>::eq(x.range(),y.range()); }
  static bool ne(const PV& x, const PV& y)
    { return mc::Op<T>::ne(x.range(),y.range()); }
  static bool lt(const PV& x, const PV& y)
    { return mc::Op<T>::lt(x.range(),y.range()); }
  static bool le(const PV& x, const PV& y)
    { return mc::Op<T>::le(x.range(),y.range()); }
  static bool gt(const PV& x, const PV& y)
    { return mc::Op<T>::gt(x.range(),y.range()); }
  static bool ge(const PV& x, const PV& y)
    { return mc::Op<T>::ge(x.range(),y.range()); }
};

} // namespace mc

#endif
