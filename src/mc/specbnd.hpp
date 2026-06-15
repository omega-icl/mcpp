// Copyright (C) 2021 Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the Eclipse Public License.

/*!
\page page_SPECBND Eigenvalue Arithmetic for Factorable Functions
\author Nikola Peric, Akshay Shah, Jai Rajyaguru & Beno&icirc;t Chachuat

Given a factorable, multivariate function \f$f:\mathbb{R}^n\to\mathbb{R}\f$, that is twice-continuously differentiable on a box \f$X:=[x^{\rm L},x^{\rm U}]\f$, M&ouml;nnigmann's technique [M&ouml;nnigmann, 2008; 2011] provides a means for computing spectral bounds of its Hessian matrix \f$\nabla^2f(x)\f$ at any point \f$x\in X\f$&mdash;without actually computing \f$\nabla^2f(x)\f$. Applications of this technique are in determining whether a function is convex or concave on a particular domain, as well as in constructing convex/concave relaxations for complete search approaches in global optimization. Alternative techniques for determining spectral bounds include the interval variant of Gershgorin's circle criterion [Adjiman <I>et al.</I>, 1998] as well as Hertz & Rohn's method [Hertz, 1992], which rely on an interval enclosure of the set of all possible Hessian matrices \f$[\nabla^2f] \supseteq \{\nabla^2f(x) \mid x\in X\}\f$. See also \ref page_MCCORMICK for an alternative way of constructing convex/concave relaxations. 

The class mc::Specbnd provides an implementation of eigenvalue arithmetic. It relies on the operator/function overloading mechanism of C++, which makes the computation of spectral bounds both simple and intuitive, similar to computing function values in real arithmetics. Moreover, mc::Specbnd can be used as the template parameter of other available types in MC++; for instance, mc::Specbnd can be used in order to propagate spectral bounds on the remainder term of polynomial models in mc::TVar, mc::CVar or mc::SCVar. Likewise, mc::Specbnd can be used as the template parameter of the types fadbad::F, fadbad::B and fadbad::T of <A href="http://www.fadbad.com/fadbad.html">FADBAD++</A> for computing spectral bounds of either the partial derivatives or the Taylor coefficients of a factorable function too (see \ref sec_SPECBND_fadbad).

mc:Specbnd is templated in the type used to propagate the necessary bounds. By default, mc::Specbnd can be used with the non-verified interval type mc::Interval of MC++. For reliability, however, it is strongly recommended to use verified interval arithmetic such as <A href="http://www.ti3.tu-harburg.de/Software/PROFILEnglisch.html">PROFIL</A> (header file <tt>mcprofil.hpp</tt>), <A href="https://www.boost.org/doc/libs/1_68_0/libs/numeric/interval/doc/interval.htm">Boost Interval Arithmetic Library</A> (header file <tt>mcboost.hpp</tt>) or <A href="http://www2.math.uni-wuppertal.de/wrswt/software/filib.html">FILIB++</A> (header file <tt>mcfilib.hpp</tt>). Other types, such as mc::McCormick, mc::TVar, mc::CVar or mc::SCVar, can also be used as template parameters of mc::Specbnd, thereby making it possible to compute, respectively, convex/concave bounds and polynomial inclusions of the spectrum of \f$\nabla^2f\f$.

As well as propagating spectral bounds for factorable functions using M&ouml;nnigmann's technique, mc::Specbnd provides support for computing spectral bounds based on the interval Hessian matrix of a twice continuously-differentiable, factorable function using either Gershgorin's circle criterion or Hertz & Rohn's method. As established by Darup <I>et al.</I> [2012], the spectral bound arithmetic may or may not produce tighter spectral bounds than these latter techniques, depending on the factorable function \f$f\f$ at hand and its variable range \f$X\f$.

Results obtained for the factorable function \f$f(x,y)=1+x-\sin(2x+3y)-\cos(3x-5y)\f$  for \f$x\in [-0.5,0.5]^2\f$ are shown in the figure below. This function (red line) is shown in the left plot and the minimal and maximal eigenvalues (red and green lines) as well as the computed spectral bounds (blue line) of the Hessian matrix \f$\nabla^2f\f$ are shown on the right plot.

<CENTER><TABLE BORDER=0>
<TR>
<TD>\image html SB-2D_function.png</TD>
<TD>\image html SB-2D_spectral.png</TD>
</TR>
</TABLE></CENTER>


\section sec_SPECBND_I How do I compute (interval) bounds on the spectrum of the Hessian matrix of a factorable function?

Suppose we want to compute a spectral (interval) bound for the Hessian matrix of the real-valued function \f$f(x_1,x_2,x_3)=\exp(x_1-2x_2^2+3x_3^3)\f$ for \f$(x_1,x_2,x_3)\in [-0.3,0.2]\times[-0.1,0.6]\times[-0.4,0.5]\f$. For simplicity, this bound is calculated using the default interval type, mc::Interval:

\code
      #include "interval.hpp"
      #include "specbnd.hpp"
      typedef mc::Interval I;
      typedef mc::Specbnd<I> SB;
\endcode

First, the variables \f$x_1\f$, \f$x_2\f$ and \f$x_3\f$ are defined as follows:

\code
      SB X1( I(-0.3,0.2), 0, 3 );
      SB X2( I(-0.1,0.6), 1, 3 );
      SB X3( I(-0.4,0.5), 2, 3 );
\endcode

Essentially, the first line means that <tt>X1</tt> is a variable of class mc::Specbnd, belongs to the interval \f$[-0.3,0.2]\f$, and having index 0 out of 3 independent variables (recall that indexing in C/C++ starts at 0 by convention!). The same holds for the mc::Specbnd variable <tt>X2</tt> and <tt>X3</tt>.

Having defined the variables, spectral bounds on the Hessian matrix \f$\nabla^2f\f$ of \f$f\f$ on \f$[-0.3,0.2]\times[-0.1,0.6]\times[-0.4,0.5]\f$ are simply calculated as:

\code
      SB F = exp( X1 - 2*sqr(X2) + 3*pow(X3,3) );
\endcode

The computed spectral bounds can be retrieved as:

\code
      I specF = F.SI();
\endcode

The function and first-derivative bounds can also be retrieved using mc::Specbnd::I and mc::Specbnd::FI. The results of the spectral bound propagation can be displayed to the standard output as:

\code
      std::cout<< "Spectral bound (eigenvalue arithmetic): " << F << std::endl;
\endcode

which produces the following display:
\verbatim
Spectral bound (eigenvalue arithmetic): [ -1.99039e+01 :  3.70043e+01 ]
\endverbatim

As noted earlier, other methods are available for bounding the spectrum of interval Hessian matrices, such as Gershgorin's circle criterion and Hertz & Rohn's method. mc::Specbnd also provides a means for computing these bounds, e.g. for comparison with the M&ouml;nnigmann's eigenvalue arithmetic. Interval Hessian matrices can be computed using the forward and/or reverse AD types of <A href="http://www.fadbad.com/fadbad.html">FADBAD++</A>:

\code
      #include "mcfadbad.hpp"
      typedef fadbad::F<I> FI;
      typedef fadbad::F<FI> FFI;
      typedef fadbad::B<FI> BFI;
\endcode

Then, in order to compute spectral bounds for \f$f(x_1,x_2,x_3)=\exp(x_1-2x_2^2+3x_3^3)\f$ on \f$(x_1,x_2,x_3)\in [-0.3,0.2]\times[-0.1,0.6]\times[-0.4,0.5]\f$, we proceed as follows:

\code
      FI FX1 = I(-0.3,0.2); FX1.diff(0,3);
      FI FX2 = I(-0.1,0.6); FX2.diff(1,3);
      FI FX3 = I(-0.4,0.5); FX3.diff(2,3);

      FFI FFX1 = FX1; FFX1.diff(0,3);
      FFI FFX2 = FX2; FFX2.diff(1,3);
      FFI FFX3 = FX3; FFX3.diff(2,3);
      FFI FFZ = exp( FFX1 - 2*pow(FFX2,2) + 3*pow(FFX3,3) );
      std::pair<double,double> spbndF;
      
      SB::options.HESSBND = SB::Options::GERSHGORIN;
      spbndF = SB::spectral_bound( FFZ );
      std::cout << "Spectral bound (Gershgorin, forward-forward): " << I(spbndF.first,spbndF.second) << std::endl;

      SB::options.HESSBND = SB::Options::HERTZROHN;
      spbndF = SB::spectral_bound( FFZ );
      std::cout << "Spectral bound (Hertz&Rohn, forward-forward): " << I(spbndF.first,spbndF.second) << std::endl;
\endcode

producing:
\verbatim
Spectral bound (Gershgorin, forward-forward): [ -2.63904e+01 :  3.85859e+01 ]
Spectral bound (Hertz&Rohn, forward-forward): [ -2.11973e+01 :  3.05272e+01 ]
\endverbatim

In this case, the eigenvalue arithmetic provides tighter bounds than with Gershgorin's circle criterion. Moreover, the lower bound is also tighter than with Herz & Rohn's method, yet the upper bound is weaker.

Spectral bounds for interval Hessian matrices generated with the forward-reverse mode of AD can also be computed:

\code
      BFI BFX[3] = { FX1, FX2, FX3 };
      BFI BFZ = exp( BFX[0] - 2*pow(BFX[1],2) + 3*pow(BFX[2],3) );
      BFZ.diff(0,1);
      std::pair<double,double> spbndB;

      SB::options.HESSBND = SB::Options::GERSHGORIN;
      spbndB = SB::spectral_bound( BFX );
      std::cout << "Spectral bound (Gershgorin, forward-reverse): " << I(spbndB.first,spbndB.second) << std::endl;

      SB::options.HESSBND = SB::Options::HERTZROHN;
      spbndB = SB::spectral_bound( BFX );
      std::cout << "Spectral bound (Hertz&Rohn, forward-reverse):  " << I(spbndB.first,spbndB.second) << std::endl;
\endcode

producing:
\verbatim
Spectral bound (Gershgorin, forward-reverse): [ -2.63904e+01 :  3.85859e+01 ]
Spectral bound (Hertz&Rohn, forward-reverse): [ -2.11973e+01 :  3.05272e+01 ]
\endverbatim

In this case, the spectral bounds computed from the Hessian matrix obtained with the forward-forward and forward-reverse modes of AD are indeed identical, although such may not always be the case; see [Darup <I>et al.</I>, 2012]

\section sec_SPECBND_MC How do I compute convex/concave relaxations on the spectrum of the Hessian matrix of a factorable function?

Instead of using standard interval arithmetic for propagating spectral bounds, we can as well make use of the McCormick relaxation technique to propagate convex/concave bounds. Using MC++, this is done simply by selecting mc::McCormick as the template parameter in mc::Specbnd:

\code
      #include "mccormick.hpp"
      typedef mc::McCormick<I> MC;
      typedef mc::Specbnd<MC> SBMC;
\endcode

Then, the procedure for computing a spectral bound remains essentially the same as described in the previous section. In order to compute convex/concave spectral relaxations and subgradients at point \f$(0,0,0)\f$ of \f$f(x_1,x_2,x_3)=\exp(x_1-2x_2^2+3x_3^3)\f$ for \f$(x_1,x_2,x_3)\in [-0.3,0.2]\times[-0.1,0.6]\times[-0.4,0.5]\f$, we proceed as follows:

\code
      SBMC X1( MC(I(-0.3,0.2),0.).sub(3,0), 0, 3 );
      SBMC X2( MC(I(-0.1,0.6),0.).sub(3,1), 1, 3 );
      SBMC X3( MC(I(-0.4,0.5),0.).sub(3,2), 2, 3 );
      SBMC F = exp( X1 - 2*sqr(X2) + 3*pow(X3,3) );
      MC specF = F.SI();
      std::cout << F << std::endl;
\endcode

The only difference here concerns the initialization of the McCormick variables inside the variables X and Y, which passes the bounds for the variable as well as the point at which the convex/concave bounds and their subgradients are computed&mdash;See How do I compute McCormick relaxations of a factorable function?

More information on the spectral relaxations can again be obtained as:

\code
      std::cout<< "Spectral bound (eigenvalue arithmetic): " << F << std::endl;
\endcode

which displays:
\verbatim
Spectral bound (eigenvalue arithmetic): [ -1.99039e+01 :  3.70043e+01 ] [ -1.54413e+01 :  2.81720e+01 ] [ (-9.27293e+00, 0.00000e+00,-5.21602e+00) : ( 1.72398e+01, 0.00000e+00, 1.50542e+01) ]
\endverbatim

By construction, the interval bounds supporting the McCormick relaxations are identical to their spectral interval bound counterparts. The convex/concave bounds at \f$(0,0,0)\f$ come as additional information and are seen to be tighter than the interval bounds. The corresponding subgradients can also be used to construct affine relaxations.

A polynomial inclusion of the spectrum may be computed similarly by selecting mc::TVar, mc::CVar or mc::SCVar as the template parameter in mc::Specbnd.


\section sec_SPECBND_fct Which functions are overloaded in mc::Specbnd eigenvalue arithmetic?

mc::Specbnd overloads the usual functions <tt>exp</tt>, <tt>log</tt>, <tt>sqr</tt>, <tt>sqrt</tt>, <tt>pow</tt>, <tt>inv</tt>, <tt>cos</tt>, <tt>sin</tt>, <tt>tan</tt>, <tt>acos</tt>, <tt>asin</tt>, <tt>atan</tt>, <tt>cosh</tt>, <tt>sinh</tt>, <tt>tanh</tt>. The functions <tt>min</tt>, <tt>max</tt>, <tt>fabs</tt>, <tt>fstep</tt> and <tt>bstep</tt>  are not overloaded in mc::Specbnd since they are not (twice) continuously differentiable.


\section sec_SPECBND_fadbad How do I compute spectral bounds of the partial derivatives or the Taylor coefficients of a factorable function using FADBAD++?

The combination of mc::Specbnd with the classes fadbad::F, fadbad::B and fadbad::T of <A href="http://www.fadbad.com/fadbad.html">FADBAD++</A> to compute a spectral bound for the Hessian matrix of either the partial derivatives or the Taylor coefficients of a factorable function is essentially the same as with mc::McCormick (see \ref sec_MCCORMICK_fadbad).

Next, we present the case of fadbad::F only. Continuing the previous example, spectral bounds of the partial derivatvies of \f$f(x_1,x_2,x_3)=\exp(x_1-2x_2^2+3x_3^3)\f$ for \f$(x_1,x_2,x_3)\in [-0.3,0.2]\times[-0.1,0.6]\times[-0.4,0.5]\f$ can be computed as follows:

\code
      typedef fadbad::F<SB> FSB;
\endcode

\code
      FSB FSBX1 = X1; FSBX1.diff(0,3);
      FSB FSBX2 = X2; FSBX2.diff(1,3);
      FSB FSBX3 = X3; FSBX3.diff(2,3);
      FSB FSBF = exp( FSBX1 - 2*pow(FSBX2,2) + 3*pow(FSBX3,3) );
      std::cout << "Spectral bounds of df/dx1: " << FSBF.d(0) << std::endl;
      std::cout << "Spectral bounds of df/dx2: " << FSBF.d(1) << std::endl;
      std::cout << "Spectral bounds of df/dx3: " << FSBF.d(2) << std::endl;
\endcode

producing the output:

\verbatim
Spectral bounds of df/dx1: [ -1.99039e+01 :  3.70043e+01 ]
Spectral bounds of df/dx2: [ -1.16096e+02 :  8.92716e+01 ]
Spectral bounds of df/dx3: [ -1.28567e+02 :  2.06229e+02 ]
\endverbatim



\section sec_SPECBND_opt How are the options set for the computation of a spectral bound?

The class mc::Specbnd has a public static member called mc::Specbnd::options that can be used to set/modify the options; e.g.,

\code
      mc::Specbnd<I>::options.HESSBND = mc::Specbnd<I>::Options::HERTZROHN;
\endcode

The available options are the following:

<TABLE border="1">
<CAPTION><EM>Options in mc::Specbnd::Options: name, type and description</EM></CAPTION>
     <TR><TH><b>Name</b>  <TD><b>Type</b><TD><b>Default</b>
         <TD><b>Description</b>
     <TR><TH><tt>HESSBND</tt> <TD><tt>mc::Specbnd::Options::HESSBND_STRATEGY</tt> <TD>mc::Specbnd::Options::GERSHGORIN
         <TD>Strategy for computing spectral bounds in interval Hessian matrix using mc::Specbnd::spectral_bound
</TABLE>


\section sec_SPECBND_err Errors What errors can I encounter during computation of a spectral bound?

Errors are managed based on the exception handling mechanism of the C++ language. Each time an error is encountered, a class object of type mc::Specbnd::Exceptions is thrown, which contains the type of error. It is the user's responsibility to test whether an exception was thrown during the computation of a spectral bound, and then make the appropriate changes. Should an exception be thrown and not caught by the calling program, the execution will abort.

Possible errors encountered during the computation of a spectral bound are:

<TABLE border="1">
<CAPTION><EM>Errors during the Computation of a Spectral Bound</EM></CAPTION>
     <TR><TH><b>Number</b> <TD><b>Description</b>
     <TR><TH><tt>1</tt> <TD>Failed to compute spectrum in Specbnd::spectrum
     <TR><TH><tt>2</tt> <TD>Failed to compute spectral bound in Specbnd::spectral_bound
     <TR><TH><tt>-1</tt> <TD>Operation between variables with different numbers of dependents
     <TR><TH><tt>-33</tt> <TD>Feature not yet implemented in mc::Specbnd
</TABLE>

Moreover, exceptions may be thrown by the template parameter class itself.


\section sec_SPECBND_refs References

- Adjiman, C.S., S. Dallwig, C.A. Floudas, and A. Neumaier, <A href="http://dx.doi.org/10.1016/S0098-1354(98)00027-1">A global optimization method, \f$\rm\alpha BB\f$, for general twice-differentiable constrained NLPs-I. Theoretical advances</A>, <I>Computers & Chemical Engineering</I> <B>22</B>(9):1137-1158, 1998.
- Hertz, D., <A href="http://dx.doi.org/10.1109/9.126593">The extreme eigenvalues and stability of real symmetric interval matrices</A>, <I>IEEE Transactions on Automatic Control</I> <B>37</B>:532-535, 1992.
- M&ouml;nnigmann, M., <A href="http://dx.doi.org/10.1137/070704186">Efficient calculation of bounds on spectra of Hessian matrices</A>, <i>SIAM Journal on Scientific Computing</i>, <b>30</b>:2340-2357, 2008.
- M&ouml;nnigmann, M., <A href="http://dx.doi.org/10.1137/10078760X">Fast Calculation of Spectral Bounds for Hessian Matrices on Hyperrectangles</A>, <i>SIAM Journal on Matrix Analysis and Applications</i>, <b>32</b>:4, 1351-1366, 2011.
- Darup, M. S., M. Kastsian, S. Mross, and M. M&ouml;nnigmann, <A href="http://arxiv.org/pdf/1206.0196.pdf">Efficient Computation of Spectral Bounds for Hessian Matrices on Hyperrectangles for Global Optimization</A>, arXiv:1206.0196v1, 1 June 2012.
.

*/

#ifndef MC__SPECBND_H
#define MC__SPECBND_H

#include <iostream>
#include <cmath>
#include <set>
#include <algorithm>

#if defined( MC__USE_ARMADILLO )
 #include <armadillo>
#else
 #include "mclapack.hpp"
#endif
#include "mcfadbad.hpp"

#undef  MC__SPECBND_DEBUG_SPECTRUM
#undef  MC__SPECBND_DEBUG_HESSBND

namespace mc
{
//! @brief C++ template class computing spectral bounds for the Hessian matrix of a factorable function on a box
////////////////////////////////////////////////////////////////////////
//! mc::Specbnd<T> is a C++ template class computing spectral bounds
//! for the Hessian matrix of a factorable function on a box. The
//! template parameter T corresponds to the type used in the underlying
//! interval arithmetic calculations.
////////////////////////////////////////////////////////////////////////
template <typename T>
class Specbnd
////////////////////////////////////////////////////////////////////////
{
  //template <typename U> friend class Specbnd;
  template <class U> friend std::ostream& operator<<( std::ostream&, Specbnd<U> const& );

private:
  //! @brief Number of independent variables
  unsigned           _n;
  //! @brief Set of variable dependencies
  std::set<unsigned> _D;
  //! @brief Set of nonlinear variable dependencies
  std::set<unsigned> _N;
  //! @brief Gradient bounds
  fadbad::F<T>      _FI;
  //! @brief Spectral bound
  T                 _SI;

  //! @brief Interval operator for spectral bound propagation in univariate terms
  static T _LambdaS( T const* a, std::set<unsigned> Da );
  //! @brief Interval operator for spectral bound propagation in product terms
  static T _LambdaT( T const* a, T const* b, std::set<unsigned> Dab );
  //! @brief Interval operator for spectral bound propagation in product terms
  static T _LambdaO( T const& a, T const& b, T const& c );

  //! @brief Computing spectral bound of interval Hessian matrix using Gershgorin's circle criterion
  static std::pair<double,double> _gershgorin_bound
    ( const unsigned N, const T*D2F );
  //! @brief Computing spectral bound of interval Hessian matrix using Hertz & Rohn's method
  static std::pair<double,double> _hertzrohn_bound
    ( const unsigned N, const T*D2F );

public: 
  // other operator overloadings
  Specbnd<T>& operator+=
    ( Specbnd<T> const& );
  Specbnd<T>& operator+=
    ( double const& );
  Specbnd<T>& operator-=
    ( Specbnd<T> const& );
  Specbnd<T>& operator-=
    ( double const& );
  Specbnd<T>& operator*=
    ( Specbnd<T> const& );
  Specbnd<T>& operator*=
    ( double const& );
  Specbnd<T>& operator/=
    ( Specbnd<T> const& );
  Specbnd<T>& operator/=
    ( Specbnd<T> && );
  Specbnd<T>& operator/=
    ( double const& );
  Specbnd<T> & operator=
    ( Specbnd<T> const& );
  Specbnd<T> & operator=
    ( Specbnd<T> && );
  Specbnd<T> & operator=
    ( double const& c );
  Specbnd<T> & operator=
    ( T const& c );

//  /** @defgroup SPECBND Eigenvalue Arithmetic for Factorable Functions
//   *  @{
//   */
//  //! @brief Options of mc::Specbnd
//  static struct Options
//  {
//    //! @brief Constructor
//    Options():
//      HESSBND(GERSHGORIN)
//      {}
//    //! @brief Method to bound eignevalues in interval Hessian matrix using mc::Specbnd::spectral_bound
//    HESSBND_STRATEGY HESSBND;
//  } options;

  //! @brief Exceptions of mc::Specbnd
  class Exceptions
  {
  public:
    //! @brief Enumeration type for Specbnd exception handling
    enum TYPE{
      SPECTR=1, //!< Failed to compute spectrum in Specbnd::spectrum
      HESSBND,  //!< Failed to compute spectral bound in Specbnd::spectral_bound
      SIZE=-1,  //!< Operation between variables with different numbers of dependents
      UNDEF=-33 //!< Feature not yet implemented in mc::Specbnd
    };
    //! @brief Constructor for error <a>ierr</a>
    Exceptions( TYPE ierr ) : _ierr( ierr ){}
    //! @brief Inline function returning the error flag
    int ierr(){ return _ierr; }
    //! @brief Error description
    std::string what(){
      switch( _ierr ){
      case SPECTR:
        return "mc::Specbnd\t Computation of Hessian spectrum failed";
      case HESSBND:
        return "mc::Specbnd\t Computation of Hessian spectral bound failed";
      case SIZE:
        return "mc::Specbnd\t Operation between variables with different numbers of dependents";
      case UNDEF:
        return "mc::Specbnd\t Feature not yet implemented in mc::Specbnd class";
      default:
        return "mc::Specbnd\t Undocumented error";
      }
    }

  private:
    TYPE _ierr;
  };

  //! @brief Default constructor (needed to declare arrays of Specbnd)
  Specbnd
    ():
    _n (0),
    _FI(0.),
    _SI(0.)
    /* _D and _N empty */
    {}

  //! @brief Constructor for real scalar <tt>c</tt>
  Specbnd
    ( double const& c ):
    _n (0),
    _FI(c),
    _SI(0.)
    /* _D and _N empty */
    {}

  //! @brief Constructor for an interval <tt>B</tt>
  Specbnd
    ( T const& B ):
    _n (0),
    _FI(B),
    _SI(0.)
    /* _D and _N empty */
    {}

  //! @brief Constructor for a variable with range <tt>B</tt> and index <a>i</a> of <a>n</a> independent variables
  Specbnd
    ( T const& B, unsigned const i, unsigned const n ):
    _n  (n),
    _FI (B),
    _SI (0.)
    /* _N empty */
    {
      _D.insert(i);
      _FI.diff(i,_n);
    }

  //! @brief Copy constructor
  Specbnd
    ( Specbnd<T> const& x ):
    _n ( x._n ),
    _D ( x._D ),
    _N ( x._N ),
    _FI( x._FI ),
    _SI( x._SI )
    {}

  //! @brief Move constructor
  Specbnd
    ( Specbnd<T> && x ):
    _n ( x._n ),
    _D ( std::move(x._D) ),
    _N ( std::move(x._N) ),
    _FI( std::move(x._FI) ),
    _SI( std::move(x._SI) )
    {}

  //! @brief Destructor
  ~Specbnd
    ()
    {}

  //! @brief Set variable with range <tt>B</tt> and index <a>i</a> of <a>n</a> independent variables
  Specbnd<T>& set
    ( T const& B, unsigned const i, unsigned const n )
    {
      if( i >= n )
        throw typename Specbnd<T>::Exceptions( Specbnd<T>::Exceptions::SIZE );
      _n  = n;
      _D  = {i};
      _N.clear();
      _FI = B;
      _FI.diff(i,_n);
      _SI = 0.;
      return *this;
    }

  //! @brief Set function and first derivative bounds as well as spectral bounds to, respectively, <tt>FB</tt> and <tt>SB</tt>
  Specbnd<T>& set
    ( std::set<unsigned> const& D, fadbad::F<T> const& FB, T const& SB );

  //! @brief Set the index of a variable (and total number of variables)
  Specbnd<T>& dep
    ( unsigned const i, unsigned const n )
    {
      if( i >= n )
        throw typename Specbnd<T>::Exceptions( Specbnd<T>::Exceptions::SIZE );
      _n  = n;
      _D  = {i};
      _N.clear();
      _FI.diff( i, n );
      _SI = 0.;
      return *this;
    }

  //! @brief Check proper subsets: -1 if Sx < Sy; +1 if Sy < Sx; 0 if Sx = Sy; +2 other
  static int include
    ( std::set<unsigned> const& Sx, std::set<unsigned> const& Sy )
    {
      if( Sx.size() < Sy.size() ){
        if( std::includes( Sy.begin(), Sy.end(), Sx.begin(), Sx.end() ) ) return -1;
      }
      else if( Sx.size() > Sy.size() ){
        if( std::includes( Sx.begin(), Sx.end(), Sy.begin(), Sy.end() ) ) return 1;
      }
      else/* Sx.size() == Sy.size() */{
        if( std::includes( Sx.cbegin(), Sx.cend(), Sy.cbegin(), Sy.cend() ) ) return 0;
      }
      return 2;
    }

  //! @brief Check intersection
  static int inter
    ( std::set<unsigned> const& Sx, std::set<unsigned> const& Sy )
    {
      // Edge Case Optimization: If either set is empty, the intersection must be empty.
      if( Sx.empty() || Sy.empty() ) return false;

      auto itx = Sx.cbegin(), ity = Sy.cbegin(), ex = Sx.cend(), ey = Sy.cend();
      // Perform linear merge-like scan
      while( itx != ex && ity != ey ){
        // Found a common element, so the intersection is not empty.
        if( *itx == *ity )     return true;
        // Element in set1 is smaller, so advance Sx iterator to look for a match
        else if( *itx < *ity ) ++itx;
        // Element in set2 is smaller, so advance Sy iterator to look for a match
        else                   ++ity;
      }
      // If the loop finishes without finding a match, the intersection is empty
      return false;
    }

  //! @brief Check cover
  static bool cover
    ( std::set<unsigned> const& Sx, std::set<unsigned> const& Sy, unsigned N )
    {
      if( Sx.size() + Sy.size() < N ) return false;

      auto itx = Sx.cbegin(), ity = Sy.cbegin(), ex = Sx.cend(), ey = Sy.cend();
      // Iterate through the expected numbers from 0 to N
      for( unsigned k = 0; k < N; ++k ){
        bool found = false;
        // Check if the expected k value is at the head of Sx
        // Since sets are sorted, we only need to look at the current iterator
        if( itx != ex && *itx == k ){
          found = true;
          ++itx; // Advance iterator only if we matched the number
        }
        // Check if the expected k value is at the head of Sy
        if( ity != ey && *ity == k ){
          found = true;
          ++ity; // Advance iterator
        }
        // If neither had the expected k value, we found a gap
        if( !found ) return false;
      }

      // If we finished the loop, we successfully found all 0..N-1
      return true;
    }

  //! @brief Compose variable with negative function
  static void neg
    ( Specbnd<T>& x )
    {
      x._FI = operator-( std::move( x._FI ) );
      x._SI = operator-( std::move( x._SI ) );
    }

  //! @brief Add two variables
  static void add
    ( Specbnd& x, Specbnd const& y )
    {
      if( x._n && y._n && x._n != y._n )
        throw Exceptions( Exceptions::SIZE );
        
      // Cases 1-3 in Schultze Darup & Monnigmann (2015), Table 4.1
      if( x._N.empty() || y._N.empty() )
        x._SI += y._SI;
      // Case 4 in Schultze Darup & Monnigmann (2015), Table 4.1
      else if( !inter( x._N, y._N ) )
        x._SI = Op<T>::hull( x._SI, y._SI );
      else{
        switch( include( x._N, y._N ) ){
          // Case 5 in Schultze Darup & Monnigmann (2015), Table 4.1
          case  0: x._SI += y._SI;
                   break;
          // Case 6 in Schultze Darup & Monnigmann (2015), Table 4.1
          case  1: x._SI += Op<T>::hull( y._SI, 0. );
                   break;
          // Case 7 in Schultze Darup & Monnigmann (2015), Table 4.1
          case -1: x._SI  = Op<T>::hull( x._SI, 0. ) + y._SI;
                   break;
          // Case 8 in Schultze Darup & Monnigmann (2015), Table 4.1
          default: x._SI  = Op<T>::hull( x._SI, 0. ) + Op<T>::hull( y._SI, 0. );
                   break;
        }
      }
      x._FI += y._FI;
      x._n   = x._FI.size();
      x._D.insert( y._D.cbegin(), y._D.cend() );
      x._N.insert( y._N.cbegin(), y._N.cend() );
    }

  //! @brief Multiply two variables
  static void multiply
    ( Specbnd& x, Specbnd const& y )
    {
      if( x._n && y._n && x._n != y._n )
        throw Exceptions( Exceptions::SIZE );

      std::set<unsigned> Dxy = x._D; Dxy.insert( y._D.cbegin(), y._D.cend() );

      // Case 1 in Schultze Darup & Monnigmann (2015), Table 4.4
      if( x._N.empty() && y._N.empty() ){
        x._SI = Specbnd<T>::_LambdaT( (x._FI.size()? &x._FI[0]: nullptr),
                                      (y._FI.size()? &y._FI[0]: nullptr), Dxy );
      }
      else if( y._N.empty() ){
        // Case 2 in Schultze Darup & Monnigmann (2015), Table 4.4
        if( Dxy.size() == x._N.size() ){
          x._SI *= y._FI.val();
          x._SI += Specbnd<T>::_LambdaT( (x._FI.size()? &x._FI[0]: nullptr),
                                         (y._FI.size()? &y._FI[0]: nullptr), Dxy );
        }
        else{
          // Case 3 in Schultze Darup & Monnigmann (2015), Table 4.4
          if( x._D.size() != 1 || y._D.size() != 1 || *x._D.cbegin() == *y._D.cbegin() ){
            x._SI  = y._FI.val() * Op<T>::hull( x._SI, 0. );
            x._SI += Specbnd<T>::_LambdaT( (x._FI.size()? &x._FI[0]: nullptr),
                                           (y._FI.size()? &y._FI[0]: nullptr), Dxy );          
          }
          // Case 4 in Schultze Darup & Monnigmann (2015), Table 4.4
          else{
            x._SI = Specbnd<T>::_LambdaO( x._SI *= y._FI.val(), 0., 
                                          x._FI[*x._D.cbegin()]*y._FI[*y._D.cbegin()] );
          }
        }
      }
      else if( x._N.empty() ){
        // Case 5 in Schultze Darup & Monnigmann (2015), Table 4.4
        if( Dxy.size() == y._N.size() ){
          x._SI  = x._FI.val() * y._SI;
          x._SI += Specbnd<T>::_LambdaT( (x._FI.size()? &x._FI[0]: nullptr),
                                         (y._FI.size()? &y._FI[0]: nullptr), Dxy );
        }
        else{
          // Case 6 in Schultze Darup & Monnigmann (2015), Table 4.4
          if( x._D.size() != 1 || y._D.size() != 1 || *x._D.cbegin() == *y._D.cbegin() ){
            x._SI  = y._FI.val() * Op<T>::hull( x._SI, 0. );
            x._SI += Specbnd<T>::_LambdaT( (x._FI.size()? &x._FI[0]: nullptr),
                                           (y._FI.size()? &y._FI[0]: nullptr), Dxy );          
          }
          // Case 7 in Schultze Darup & Monnigmann (2015), Table 4.4
          else{
            x._SI = Specbnd<T>::_LambdaO( 0., y._SI * x._FI.val(), 
                                          x._FI[*x._D.cbegin()]*y._FI[*y._D.cbegin()] );
          }
        }
      }
      else if( !inter( x._N, y._N ) ){
        std::set<unsigned> Nxy = x._N; Nxy.insert( y._N.cbegin(), y._N.cend() );
        // Case 8 in Schultze Darup & Monnigmann (2015), Table 4.4
        if( Dxy.size() > Nxy.size() ){
          x._SI  = Op<T>::hull( Op<T>::hull( x._FI.val() * y._SI,
                                             y._FI.val() * x._SI ), 0. );
          x._SI += Specbnd<T>::_LambdaT( (x._FI.size()? &x._FI[0]: nullptr),
                                         (y._FI.size()? &y._FI[0]: nullptr), Dxy );
        }
        else{
          // Case 9 in Schultze Darup & Monnigmann (2015), Table 4.4
          if( x._D.size() != 1 || y._D.size() != 1 || *x._D.cbegin() == *y._D.cbegin() ){
            x._SI  = Op<T>::hull( x._FI.val() * y._SI,
                                  y._FI.val() * x._SI );
            x._SI += Specbnd<T>::_LambdaT( (x._FI.size()? &x._FI[0]: nullptr),
                                           (y._FI.size()? &y._FI[0]: nullptr), Dxy );
          }
          // Case 10 in Schultze Darup & Monnigmann (2015), Table 4.4
          else/* if( x._D.size() == 1 && y._D.size() == 1 && *x._D.cbegin() != *y._D.cbegin() )*/{
            x._SI = Specbnd<T>::_LambdaO( y._FI.val() * x._SI, x._FI.val() * y._SI, 
                                          x._FI[*x._D.cbegin()]*y._FI[*y._D.cbegin()] );
          }
        }
      }
      else{
        switch( include( x._N, y._N ) ){
          case  0: 
            // Case 11 in Schultze Darup & Monnigmann (2015), Table 4.4
            if( Dxy.size() == x._N.size() ){
              x._SI *= y._FI.val();
              x._SI += x._FI.val() * y._SI;
              x._SI += Specbnd<T>::_LambdaT( (x._FI.size()? &x._FI[0]: nullptr),
                                             (y._FI.size()? &y._FI[0]: nullptr), Dxy );
            }
            // Case 14 in Schultze Darup & Monnigmann (2015), Table 4.4
            else{
              x._SI  = Op<T>::hull( x._FI.val() * y._SI
                                  + y._FI.val() * x._SI, 0. );
              x._SI += Specbnd<T>::_LambdaT( (x._FI.size()? &x._FI[0]: nullptr),
                                             (y._FI.size()? &y._FI[0]: nullptr), Dxy );
            }
            break;
          case  1:
            // Case 12 in Schultze Darup & Monnigmann (2015), Table 4.4
            if( Dxy.size() == x._N.size() ){
              x._SI *= y._FI.val();
              x._SI += x._FI.val() * Op<T>::hull( y._SI, 0. );
              x._SI += Specbnd<T>::_LambdaT( (x._FI.size()? &x._FI[0]: nullptr),
                                             (y._FI.size()? &y._FI[0]: nullptr), Dxy );
            }
            // Case 15 in Schultze Darup & Monnigmann (2015), Table 4.4
            else{
              x._SI  = Op<T>::hull( x._FI.val() * Op<T>::hull( y._SI, 0. )
                                  + y._FI.val() * x._SI, 0. );
              x._SI += Specbnd<T>::_LambdaT( (x._FI.size()? &x._FI[0]: nullptr),
                                             (y._FI.size()? &y._FI[0]: nullptr), Dxy );
            }
            break;
          case -1:
            // Case 13 in Schultze Darup & Monnigmann (2015), Table 4.4
            if( Dxy.size() == x._N.size() ){
              x._SI  = Op<T>::hull( x._FI.val() * y._SI
                                  + y._FI.val() * Op<T>::hull( y._SI, 0. ), 0. );
              x._SI += Specbnd<T>::_LambdaT( (x._FI.size()? &x._FI[0]: nullptr),
                                             (y._FI.size()? &y._FI[0]: nullptr), Dxy );
            }
            // Case 16 in Schultze Darup & Monnigmann (2015), Table 4.4
            else{
              x._SI  = Op<T>::hull( x._FI.val() * Op<T>::hull( y._SI, 0. )
                                  + y._FI.val() * x._SI, 0. );
              x._SI += Specbnd<T>::_LambdaT( (x._FI.size()? &x._FI[0]: nullptr),
                                             (y._FI.size()? &y._FI[0]: nullptr), Dxy );
            }
            break;
          // Case 17 in Schultze Darup & Monnigmann (2015), Table 4.4
          default:
              x._SI  = x._FI.val() * Op<T>::hull( y._SI, 0. )
                     + y._FI.val() * Op<T>::hull( x._SI, 0. );
              x._SI += Specbnd<T>::_LambdaT( (x._FI.size()? &x._FI[0]: nullptr),
                                             (y._FI.size()? &y._FI[0]: nullptr), Dxy );
        }
      }
      x._FI *= y._FI;
      x._n = x._FI.size();
      x._D = std::move( Dxy );
      x._N = x._D;
    }

  //! @brief Compose variable with a univariate outer function
  template <typename UNIV, typename DUNIV, typename D2UNIV>
  static void compose
    ( Specbnd& x, UNIV const& f, DUNIV const& Df, D2UNIV const& D2f, bool const linear=false )
    {
      // Case 1 in Schultze Darup & Monnigmann (2015), Table 4.3
      if( x._N.empty() ){
        x._SI = D2f( x._FI.val() ) * _LambdaS( (x._FI.size()? &x._FI[0]: nullptr), x._D );
      }
      else{
        // Case 3 in Schultze Darup & Monnigmann (2015), Table 4.3
        if( !linear && x._N.size() < x._D.size() ) x._SI = Op<T>::hull( x._SI, T(0.) );
        // Case 2 in Schultze Darup & Monnigmann (2015), Table 4.3
        /* Use x._SI as is */
        x._SI *= Df( x._FI.val() );
        x._SI += _LambdaS( (x._FI.size()? &x._FI[0]: nullptr), linear? x._N: x._D )
                 * D2f( x._FI.val() );
      }
      x._FI  = f( x._FI );
      // Linear dependencies become nonlinear after nonlinear composition
      if( !linear && x._N.size() < x._D.size() ) x._N.insert( x._D.cbegin(), x._D.cend() );
    }

  //! @brief Compose variable with acos function
  static void acos
    ( Specbnd<T>& x )
    {
      // Case 1 in Schultze Darup & Monnigmann (2015), Table 4.3
      if( x._N.empty() ){
        x._SI = - _LambdaS( (x._FI.size()? &x._FI[0]: nullptr), x._D )
                * x._FI.val() * Op<T>::pow( Op<T>::sqrt( 1. - Op<T>::sqr( x._FI.val() ) ), -3 );
      }
      else{
        // Case 3 in Schultze Darup & Monnigmann (2015), Table 4.3
        if( x._N.size() < x._D.size() ) x._SI = Op<T>::hull( x._SI, T(0.) );
        // Case 2 in Schultze Darup & Monnigmann (2015), Table 4.3
        /* Use x._SI as is */
        x._SI += _LambdaS( (x._FI.size()? &x._FI[0]: nullptr), x._D )
                 * x._FI.val() * Op<T>::inv( 1. - Op<T>::sqr( x._FI.val() ) );
        x._SI *= - Op<T>::inv( Op<T>::sqrt( 1. - Op<T>::sqr( x._FI.val() ) ) );
      }
      x._FI = fadbad::acos( x._FI );
      // Linear dependencies become nonlinear after nonlinear composition
      if( x._N.size() < x._D.size() ) x._N.insert( x._D.cbegin(), x._D.cend() );
    }

  //! @brief Compose variable with asin function
  static void asin
    ( Specbnd<T>& x )
    {
      // Case 1 in Schultze Darup & Monnigmann (2015), Table 4.3
      if( x._N.empty() ){
        x._SI = _LambdaS( (x._FI.size()? &x._FI[0]: nullptr), x._D )
                * x._FI.val() * Op<T>::pow( Op<T>::sqrt( 1. - Op<T>::sqr( x._FI.val() ) ), -3 );
      }
      else{
        // Case 3 in Schultze Darup & Monnigmann (2015), Table 4.3
        if( x._N.size() < x._D.size() ) x._SI = Op<T>::hull( x._SI, T(0.) );
        // Case 2 in Schultze Darup & Monnigmann (2015), Table 4.3
        /* Use x._SI as is */
        x._SI += _LambdaS( (x._FI.size()? &x._FI[0]: nullptr), x._D )
                 * x._FI.val() * Op<T>::inv( 1. - Op<T>::sqr( x._FI.val() ) );
        x._SI *= Op<T>::inv( Op<T>::sqrt( 1. - Op<T>::sqr( x._FI.val() ) ) );
      }
      x._FI = fadbad::asin( x._FI );
      // Linear dependencies become nonlinear after nonlinear composition
      if( x._N.size() < x._D.size() ) x._N.insert( x._D.cbegin(), x._D.cend() );
    }

  //! @brief Compose variable with asin function
  static void atan
    ( Specbnd<T>& x )
    {
      // Case 1 in Schultze Darup & Monnigmann (2015), Table 4.3
      if( x._N.empty() ){
        x._SI = _LambdaS( (x._FI.size()? &x._FI[0]: nullptr), x._D )
                * (-2.) * x._FI.val() * Op<T>::pow( 1. + Op<T>::sqr( x._FI.val() ), -2 );
      }
      else{
        // Case 3 in Schultze Darup & Monnigmann (2015), Table 4.3
        if( x._N.size() < x._D.size() ) x._SI = Op<T>::hull( x._SI, T(0.) );
        // Case 2 in Schultze Darup & Monnigmann (2015), Table 4.3
        /* Use x._SI as is */
        x._SI += _LambdaS( (x._FI.size()? &x._FI[0]: nullptr), x._D )
                 * (-2.) * x._FI.val() * Op<T>::inv( 1. + Op<T>::sqr( x._FI.val() ) );
        x._SI *= Op<T>::inv( 1. + Op<T>::sqr( x._FI.val() ) );
      }
      x._FI = fadbad::atan( x._FI );
      // Linear dependencies become nonlinear after nonlinear composition
      if( x._N.size() < x._D.size() ) x._N.insert( x._D.cbegin(), x._D.cend() );
    }

  //! @brief Compose variable with tan function
  static void tan
    ( Specbnd<T>& x )
    {
      // Case 1 in Schultze Darup & Monnigmann (2015), Table 4.3
      if( x._N.empty() ){
        fadbad::F<T> z_FI = fadbad::tan( x._FI );
        x._SI = _LambdaS( (x._FI.size()? &x._FI[0]: nullptr), x._D )
                * 2. * z_FI.val() * ( 1. + Op<T>::sqr( z_FI.val() ) );
        std::swap( z_FI, x._FI );
      }
      else{
        // Case 3 in Schultze Darup & Monnigmann (2015), Table 4.3
        if( x._N.size() < x._D.size() ) x._SI = Op<T>::hull( x._SI, T(0.) );
        // Case 2 in Schultze Darup & Monnigmann (2015), Table 4.3
        /* Use x._SI as is */
        fadbad::F<T> z_FI = fadbad::tan( x._FI );
        x._SI += 2. * _LambdaS( (x._FI.size()? &x._FI[0]: nullptr), x._D ) * z_FI.val();
        x._SI *= ( 1. + Op<T>::sqr( z_FI.val() ) );
        std::swap( z_FI, x._FI );
      }
      // Linear dependencies become nonlinear after nonlinear composition
      if( x._N.size() < x._D.size() ) x._N.insert( x._D.cbegin(), x._D.cend() );
    }

  //! @brief Compose variable with tanh function
  static void tanh
    ( Specbnd<T>& x )
    {
      // Case 1 in Schultze Darup & Monnigmann (2015), Table 4.3
      if( x._N.empty() ){
        fadbad::F<T> z_FI = fadbad::tanh( x._FI );
        x._SI = _LambdaS( (x._FI.size()? &x._FI[0]: nullptr), x._D )
                * (-2.) * z_FI.val() * ( 1. - Op<T>::sqr( z_FI.val() ) );
        std::swap( z_FI, x._FI );
      }
      else{
        // Case 3 in Schultze Darup & Monnigmann (2015), Table 4.3
        if( x._N.size() < x._D.size() ) x._SI = Op<T>::hull( x._SI, T(0.) );
        // Case 2 in Schultze Darup & Monnigmann (2015), Table 4.3
        /* Use x._SI as is */
        fadbad::F<T> z_FI = fadbad::tanh( x._FI );
        x._SI += (-2.) * _LambdaS( (x._FI.size()? &x._FI[0]: nullptr), x._D ) * z_FI.val();
        x._SI *= ( 1. - Op<T>::sqr( z_FI.val() ) );
        std::swap( z_FI, x._FI );
      }
      // Linear dependencies become nonlinear after nonlinear composition
      if( x._N.size() < x._D.size() ) x._N.insert( x._D.cbegin(), x._D.cend() );
    }

  //! @brief Compose variable with exp function
  static void exp
    ( Specbnd<T>& x )
    {
      // Case 1 in Schultze Darup & Monnigmann (2015), Table 4.3
      if( x._N.empty() ){
        x._SI = Op<T>::exp( x._FI.val() ) * _LambdaS( (x._FI.size()? &x._FI[0]: nullptr), x._D );
      }
      else{
        // Case 3 in Schultze Darup & Monnigmann (2015), Table 4.3
        if( x._N.size() < x._D.size() ) x._SI = Op<T>::hull( x._SI, T(0.) );
        // Case 2 in Schultze Darup & Monnigmann (2015), Table 4.3
        /* Use x._SI as is */
        x._SI += _LambdaS( (x._FI.size()? &x._FI[0]: nullptr), x._D );
        x._SI *= Op<T>::exp( x._FI.val() );
      }
      x._FI  = fadbad::exp( x._FI );
      // Linear dependencies become nonlinear after nonlinear composition
      if( x._N.size() < x._D.size() ) x._N.insert( x._D.cbegin(), x._D.cend() );
    }

  //! @brief Compose variable with log function
  static void log
    ( Specbnd<T>& x )
    {
      // Case 1 in Schultze Darup & Monnigmann (2015), Table 4.3
      if( x._N.empty() ){
        x._SI = - _LambdaS( (x._FI.size()? &x._FI[0]: nullptr), x._D ) / Op<T>::sqr( x._FI.val() );
      }
      else{
        // Case 3 in Schultze Darup & Monnigmann (2015), Table 4.3
        if( x._N.size() < x._D.size() ) x._SI = Op<T>::hull( x._SI, T(0.) );
        // Case 2 in Schultze Darup & Monnigmann (2015), Table 4.3
        /* Use x._SI as is */
        x._SI -= _LambdaS( (x._FI.size()? &x._FI[0]: nullptr), x._D ) / x._FI.val();
        x._SI /= x._FI.val();
      }
      x._FI  = fadbad::log( x._FI );
      // Linear dependencies become nonlinear after nonlinear composition
      if( x._N.size() < x._D.size() ) x._N.insert( x._D.cbegin(), x._D.cend() );
    }

  //! @brief Compose variable with log function
  static void sqrt
    ( Specbnd<T>& x )
    {
      // Case 1 in Schultze Darup & Monnigmann (2015), Table 4.3
      if( x._N.empty() ){
        x._SI = _LambdaS( (x._FI.size()? &x._FI[0]: nullptr), x._D )
                / ( -4. *Op<T>::pow( x._FI.val(), 1.5 ) );
      }
      else{
        // Case 3 in Schultze Darup & Monnigmann (2015), Table 4.3
        if( x._N.size() < x._D.size() ) x._SI = Op<T>::hull( x._SI, T(0.) );
        // Case 2 in Schultze Darup & Monnigmann (2015), Table 4.3
        /* Use x._SI as is */
        x._SI += _LambdaS( (x._FI.size()? &x._FI[0]: nullptr), x._D ) / ( -2. * x._FI.val() );
        x._SI /= 2. * Op<T>::sqrt( x._FI.val() );
      }
      x._FI  = fadbad::sqrt( x._FI );
      // Linear dependencies become nonlinear after nonlinear composition
      if( x._N.size() < x._D.size() ) x._N.insert( x._D.cbegin(), x._D.cend() );
    }

  //! @brief Compose variable with inverse function
  static void inv
    ( Specbnd<T>& x )
    {
      // Case 1 in Schultze Darup & Monnigmann (2015), Table 4.3
      if( x._N.empty() ){
        x._SI = 2. * _LambdaS( (x._FI.size()? &x._FI[0]: nullptr), x._D )
                / Op<T>::pow( x._FI.val(), 3 );
      }
      else{
        x._SI = operator-( std::move( x._SI ) );
        // Case 3 in Schultze Darup & Monnigmann (2015), Table 4.3
        if( x._N.size() < x._D.size() ) x._SI = Op<T>::hull( x._SI, T(0.) );
        // Case 2 in Schultze Darup & Monnigmann (2015), Table 4.3
        /* Use x._SI as is */
        x._SI += 2. * _LambdaS( (x._FI.size()? &x._FI[0]: nullptr), x._D ) / x._FI.val();
        x._SI /= Op<T>::sqr( x._FI.val() );
      }
      x._FI  = 1./x._FI;
      // Linear dependencies become nonlinear after nonlinear composition
      if( x._N.size() < x._D.size() ) x._N.insert( x._D.cbegin(), x._D.cend() );
    }

  //! @brief Compose variable with integral power function
  static void ipow
    ( Specbnd<T>& x, int const m )
    {
      switch( m ){
        case 0:  x = 0.; return;
        case 1:  return;
        case 2:  return sqr(x);
        case -1: return inv(x);
        default: break;
      }

      // Case 1 in Schultze Darup & Monnigmann (2015), Table 4.3
      if( x._N.empty() ){
        x._SI = _LambdaS( (x._FI.size()? &x._FI[0]: nullptr), x._D )
                * Op<T>::pow( x._FI.val(), m-2 ) * ((m-1.) * m);
      }
      else{
        // Case 3 in Schultze Darup & Monnigmann (2015), Table 4.3
        if( x._N.size() < x._D.size() ) x._SI = Op<T>::hull( x._SI, T(0.) );
        // Case 2 in Schultze Darup & Monnigmann (2015), Table 4.3
        /* Use x._SI as is */
        x._SI *= x._FI.val();
        x._SI += _LambdaS( (x._FI.size()? &x._FI[0]: nullptr), x._D ) * (m-1.);
        x._SI *= Op<T>::pow( x._FI.val(), m-2 ) * (double)m;
      }
      x._FI  = fadbad::pow2( x._FI, m );
      // Linear dependencies become nonlinear after nonlinear composition
      if( x._N.size() < x._D.size() ) x._N.insert( x._D.cbegin(), x._D.cend() );
    }

  //! @brief Compose variable with square function
  static void sqr
    ( Specbnd<T>& x )
    {
      // Case 1 in Schultze Darup & Monnigmann (2015), Table 4.3
      if( x._N.empty() ){
        x._SI = 2. * _LambdaS( (x._FI.size()? &x._FI[0]: nullptr), x._D );
      }
      else{
        // Case 3 in Schultze Darup & Monnigmann (2015), Table 4.3
        if( x._N.size() < x._D.size() ) x._SI = Op<T>::hull( x._SI, T(0.) );
        // Case 2 in Schultze Darup & Monnigmann (2015), Table 4.3
        /* Use x._SI as is */
        x._SI *= x._FI.val();
        x._SI += _LambdaS( (x._FI.size()? &x._FI[0]: nullptr), x._D );
        x._SI *= 2.;
      }
      x._FI  = fadbad::sqr( x._FI );
      // Linear dependencies become nonlinear after nonlinear composition
      if( x._N.size() < x._D.size() ) x._N.insert( x._D.cbegin(), x._D.cend() );
    }

  //! @brief Number of independent variables
  unsigned n
    ()
    const
    {
      return _n;
    }

  //! @brief Function dependencies
  std::set<unsigned> const& D
    ()
    const
    {
      return _D;
    }

  //! @brief Nonlinear function dependencies
  std::set<unsigned> const& N
    ()
    const
    {
      return _N;
    }

  //! @brief Function range
  T I
    ()
    const
    {
      return _FI.val();
    }

  //! @brief Gradient range
  T const* FI
    ()
    const
    {
      return _FI.size()? &_FI[0]: nullptr;
    }

  //! @brief Return spectral range
  T SI
    ()
    const
    {
      return( _N.size()==_n? _SI: Op<T>::hull(_SI,0.) );
    }

  //! @brief Compute spectrum of symmetric real matrix <tt>D2F</tt> of size <tt>N</tt>
  static std::pair<double,double> spectrum
    ( unsigned const N, double const* D2F, unsigned NNZ=0, unsigned const* irow=nullptr,
      unsigned const* jcol=nullptr );

  //! @brief Compute spectral bound of symmetric interval matrix <tt>D2F</tt> of size <tt>N</tt> using Gershgorin's method
  static std::pair<double,double> spectral_bound_gershgorin
    ( unsigned const N, T const* D2F, unsigned NNZ=0, unsigned const* irow=nullptr,
      unsigned const* jcol=nullptr, double const* d=nullptr );

  //! @brief Compute spectral bound of symmetric interval matrix <tt>D2F</tt> of size <tt>N</tt> using Rohn's method
  static std::pair<double,double> spectral_bound_rohn
    ( unsigned const N, T const* D2F, unsigned NNZ=0, unsigned const* irow=nullptr,
      unsigned const* jcol=nullptr );

  //! @brief Compute spectral bound of symmetric interval matrix <tt>D2F</tt> of size <tt>N</tt> using Hertz's method
  static std::pair<double,double> spectral_bound_hertz
    ( unsigned const N, T const* D2F, unsigned NNZ=0, unsigned const* irow=nullptr,
      unsigned const* jcol=nullptr );

private:
  //! @brief Compute mid and radius of symmetric interval matrix
  static void _mid_rad
    ( unsigned const N, T const* D2F, unsigned NNZ, unsigned const* irow, unsigned const* jcol,
      arma::mat& H_mid, arma::mat& H_rad );
/*
    //! @brief Strategy for computing spectral bounds in interval Hessian matrix
    enum HESSBND_STRATEGY{
      GERSHGORIN=0,	//!< Gershgorin circle's method, O(N^2) complexity
      ROHN,		//!< Rohn's method, O(N^3) complexity
      HERTZ		//!< Hertz's (exact) method, O(2^(N-1)*N^3) complexity
    };

  //! @brief Compute bound on the real part of the spectrum of square (non-symmetric) interval matrix <tt>A</tt> of size <tt>N</tt>. The bounding method is selected via mc::Specbnd::Options::HESSBND.
  static std::pair<double,double> spectral_bound_re
    ( const unsigned N, const T*A );

  //! @brief Compute bound on the imaginary part of the spectrum of square (non-symmetric) interval matrix <tt>A</tt> of size <tt>N</tt>. The bounding method is selected via mc::Specbnd::Options::HESSBND.
  static std::pair<double,double> spectral_bound_im
    ( const unsigned N, const T*A );
*/
  /** @} */
};

////////////////////////////////////////////////////////////////////////

//template <typename T> inline typename Specbnd<T>::Options Specbnd<T>::options;

template <class T> inline T
Specbnd<T>::_LambdaS
( T const* a, std::set<unsigned> Da )
{
  if( Da.empty() || !a ) return 0.;

  unsigned i0 = *Da.cbegin();
  if( Da.size() == 1 )   return Op<T>::sqr( a[i0] );

  double upbnd( 0. );
  for( auto i : Da )
    upbnd += Op<T>::u( Op<T>::sqr( a[i] ) );
  return Op<T>::zeroone() * upbnd;
}

template <class T> inline T
Specbnd<T>::_LambdaT
( T const* a, T const* b, std::set<unsigned> Dab )
{
  if( Dab.empty() || !a || !b ) return 0.;

  unsigned i0 = *Dab.cbegin();
  if( Dab.size() == 1 ) return 2.*a[i0]*b[i0];

  double upbnda( 0. ), upbndb( 0. );
  for( auto i : Dab ){
    upbnda += Op<T>::u( Op<T>::sqr( a[i] ) );
    upbndb += Op<T>::u( Op<T>::sqr( b[i] ) );
  }
  T lamb = (2.*Op<T>::zeroone()-1.) * std::sqrt( upbnda * upbndb );
  for( auto i : Dab )
    lamb += a[i] * b[i];
  return lamb;
}

template <class T> inline T
Specbnd<T>::_LambdaO
( T const& a, T const& b, T const& c )
{
  double const d = 4 * Op<T>::u( Op<T>::sqr( c ) );
  double const ambl = Op<T>::l(a) - Op<T>::l(b);
  double const ambu = Op<T>::u(a) - Op<T>::u(b);
  T D2 = Op<T>::hull( - std::sqrt( ambl*ambl + d ), std::sqrt( ambu*ambu + d ) ); 
  return( ( D2 += ( a + b ) ) *= 0.5 );
}

template <typename T> inline std::pair<double,double>
Specbnd<T>::spectrum
( unsigned const N, double const* D2F, unsigned NNZ, unsigned const* irow, unsigned const* jcol )
{
  if( !N || !D2F || (NNZ && (!irow || !jcol) ) )
    throw Exceptions( Exceptions::SPECTR );

  arma::mat H;
  if( irow ){
    H.zeros( N, N );
    for( unsigned k=0; k<NNZ; ++k )
      H(irow[k],jcol[k]) = D2F[k];
  }
  else{
    H = arma::mat( D2F, N, N );
  }
  H = symmatu(H);
#ifdef MC__SPECBND_DEBUG_SPECTRUM
  std::cout << H << std::endl;
#endif
  arma::vec D;
  try{
    D = arma::eig_sym( H );
  }
  catch(...){
    throw Exceptions( Exceptions::SPECTR );
  }
  return std::make_pair( D(0), D(N-1) );
}

template <typename T>
inline std::pair<double,double>
Specbnd<T>::spectral_bound_gershgorin
( unsigned const N, T const* D2F, unsigned NNZ, unsigned const* irow, unsigned const* jcol,
  double const* d )
{
  if( !N || !D2F || (NNZ && (!irow || !jcol) ) )
    throw Exceptions( Exceptions::SPECTR );

  auto const& rw = []( unsigned const i, unsigned const j, unsigned const n )
                     { return i*n+j; };

  // Initialize bounds to infinity
  double min_eig =  std::numeric_limits<double>::infinity();
  double max_eig = -std::numeric_limits<double>::infinity();

  for( unsigned i=0; i<N; i++ ){
    double ri = 0.;
    T D2Fij;

    if( irow ){
      T D2Fii(0.);
      for( unsigned ij=0; ij<NNZ; ++ij ){
        if( irow[ij] != i ) continue;
        if( jcol[ij] == i ) D2Fii = D2F[ij];
        else ri += d? Op<T>::abs( D2F[ij] )*d[jcol[ij]]/d[i]: Op<T>::abs( D2F[ij] );
      }

      min_eig = std::min( min_eig, Op<T>::l( D2Fii ) - ri );
      max_eig = std::max( max_eig, Op<T>::u( D2Fii ) + ri );
      continue;
    }
    
    for( unsigned j=0; j<N; j++ ){
      if( j == i ) continue;
      if( !Op<T>::inter( D2Fij, D2F[rw(i,j,N)], D2F[rw(j,i,N)] ) ){
        D2Fij = D2F[rw(i,j,N)];
      }
      ri += d? Op<T>::abs( D2Fij )*d[j]/d[i]: Op<T>::abs( D2Fij );
    }

    min_eig = std::min( min_eig, Op<T>::l( D2F[rw(i,i,N)] ) - ri );
    max_eig = std::max( max_eig, Op<T>::u( D2F[rw(i,i,N)] ) + ri );
  }

  return std::make_pair( min_eig, max_eig );
}

template <typename T>
inline void
Specbnd<T>::_mid_rad
( unsigned const N, T const* D2F, unsigned NNZ, unsigned const* irow, unsigned const* jcol,
  arma::mat& H_mid, arma::mat& H_rad )
{
  auto const& rw = []( unsigned const i, unsigned const j, unsigned const n )
                     { return i*n+j; };

  if( irow ){
    H_mid.zeros( N, N );
    H_rad.zeros( N, N );
    for( unsigned ij=0; ij<NNZ; ++ij ){
      if( irow[ij] == jcol[ij] ){
        H_mid(irow[ij],irow[ij]) = Op<T>::mid( D2F[ij] );
        H_rad(irow[ij],irow[ij]) = 0.5 * Op<T>::diam( D2F[ij] );
        continue;
      }
      unsigned ji = ij;
      for( ++ji; ji<NNZ; ++ji ){
        if( irow[ji] == jcol[ij] && jcol[ji] == irow[ij] ) break;
      }
      if( ji >= NNZ ) continue;
      T D2Fij;
      if( !Op<T>::inter( D2Fij, D2F[ij], D2F[ji] ) ) D2Fij = D2F[ij];
      H_mid(irow[ij],jcol[ij]) = H_mid(jcol[ij],irow[ij]) = Op<T>::mid( D2Fij );
      H_rad(irow[ij],jcol[ij]) = H_rad(jcol[ij],irow[ij]) = 0.5 * Op<T>::diam( D2Fij );
    }
    return;
  }

  H_mid.set_size( N, N );
  H_rad.set_size( N, N );
  for( unsigned j=0; j<N; j++ ){
    for( unsigned i=j; i<N; i++ ){
      if( i == j ){
        H_mid(i,i) = Op<T>::mid( D2F[rw(i,i,N)] );
        H_rad(i,i) = 0.5 * Op<T>::diam( D2F[rw(i,i,N)] );
       continue;
       }
      T D2Fij;
      if( !Op<T>::inter( D2Fij, D2F[rw(i,j,N)], D2F[rw(j,i,N)] ) ) D2Fij = D2F[rw(i,j,N)];
      H_mid(i,j) = H_mid(j,i) = Op<T>::mid( D2Fij );
      H_rad(i,j) = H_rad(j,i) = 0.5 * Op<T>::diam( D2Fij );
    }
  }
}

template <typename T>
inline std::pair<double,double>
Specbnd<T>::spectral_bound_rohn
( unsigned const N, T const* D2F, unsigned NNZ, unsigned const* irow, unsigned const* jcol )
{
  if( !N || !D2F || (NNZ && (!irow || !jcol)) )
    throw Exceptions( Exceptions::SPECTR );

  // 1. Compute Midpoint (H_mid) and Radius (H_rad) matrices
  arma::mat H_mid, H_rad;
  _mid_rad( N, D2F, NNZ, irow, jcol, H_mid, H_rad );

  // 2. Calculate Rohn Bounds: [ min(H_mid) - rho(H_rad),  max(H_mid) + rho(H_rad) ]
  try{
    arma::vec eig_H_mid = arma::eig_sym( H_mid );
    arma::vec eig_H_rad = arma::eig_sym( H_rad );
    return std::make_pair( eig_H_mid[0] - eig_H_rad[N-1], eig_H_mid[N-1] + eig_H_rad[N-1] );
  }
  catch(...){
    throw Exceptions( Specbnd<T>::Exceptions::HESSBND );
  }
}

template <typename T>
inline std::pair<double,double>
Specbnd<T>::spectral_bound_hertz
( unsigned const N, T const* D2F, unsigned NNZ, unsigned const* irow, unsigned const* jcol )
{
  if( !N || !D2F || (NNZ && (!irow || !jcol)) )
    throw Exceptions( Exceptions::SPECTR );

  // 1. Compute Midpoint (H_mid) and Radius (H_rad) matrices
  arma::mat H_mid, H_rad;
  _mid_rad( N, D2F, NNZ, irow, jcol, H_mid, H_rad );
    
  // Initialize bounds to infinity
  double global_min_eig =  std::numeric_limits<double>::infinity();
  double global_max_eig = -std::numeric_limits<double>::infinity();

  // 2. Iterate through 2^(n-1) signature vectors.
  // Optimization: Fix z[0] = 1 and vary the remaining n-1 elements.
  unsigned long long num_combinations = 1ULL << (N-1);

  for( unsigned long long i=0; i<num_combinations; ++i ){
        
    // Construct signature vector z
    arma::vec z(N);
    z(0) = 1.0; 

    for( unsigned bit=0; bit+1<N; ++bit ){
      z(bit + 1) = ((i >> bit) & 1) ? -1.0 : 1.0;
    }

    // 3. Construct the Interaction Term: D_z * Delta * D_z
    // Equivalent to element-wise multiplication: Delta % (z * z^T)
    arma::mat SignPattern = z * z.t(); 
    arma::mat Term = H_rad % SignPattern; 
        
    // 4. Construct Vertices
    // For Min Bound: H_minus = H_c - Term
    // For Max Bound: H_plus  = H_c + Term
    arma::mat H_minus = H_mid - Term;
    arma::mat H_plus  = H_mid + Term;

    // 5. Compute Eigenvalues
    try{
      // Compute min eigenvalue of H_minus
      arma::vec eig_minus = arma::eig_sym( H_minus );
      double local_min = eig_minus(0); // 0 is smallest
      if( local_min < global_min_eig ){
        global_min_eig = local_min;
      }

      // Compute max eigenvalue of H_plus
      arma::vec eig_plus = arma::eig_sym( H_plus );
      double local_max = eig_plus(N-1); // N-1 is largest
      if( local_max > global_max_eig ){
        global_max_eig = local_max;
      }
    }
    catch(...){
      throw Exceptions( Specbnd<T>::Exceptions::HESSBND );
    }
  }

  return std::make_pair( global_min_eig, global_max_eig );
}
/*
template <typename T>
inline std::pair<double,double>
Specbnd<T>::spectral_bound_re
( unsigned const N, T const* A )
{
  if( !N || !A ) throw Exceptions( Exceptions::HESSBND );
#ifdef MC__SPECBND_DEBUG_HESSBND
  mc::display( N, N, A, N, "\nMatrix A", std::cout );
#endif
  std::vector<T> ARe(N*N);
  for( unsigned int i=0; i<N; i++ )
    for( unsigned int j=0; j<N; j++ )
      ARe[i+N*j] = ( A[i+N*j] + A[i*N+j] ) / 2.;
#ifdef MC__SPECBND_DEBUG_HESSBND
  mc::display( N, N, ARe.data(), N, "\nMatrix ARe", std::cout );
#endif
  return spectral_bound( N, ARe.data() );
}

template <typename T> inline std::pair<double,double>
Specbnd<T>::spectral_bound_im
( const unsigned N, const T*A )
{
  if( !N || !A ) throw Exceptions( Exceptions::HESSBND );
  const unsigned N2 = 2*N;
#ifdef MC__SPECBND_DEBUG_HESSBND
  mc::display( N, N, A, N, "\nMatrix A", std::cout );
#endif
  std::vector<T> AIm(N2*N2);
  for( unsigned int i=0; i<N; i++ )
    for( unsigned int j=0; j<N; j++ ){
      AIm[i+N2*j] = AIm[N2*N+N+i+N2*j] = 0.;
      AIm[N+i+N2*j] = ( A[i+N*j] - A[i*N+j] ) / 2.;
      AIm[N2*N+i+N2*j] = - AIm[N+i+N2*j];
    }
#ifdef MC__SPECBND_DEBUG_HESSBND
  mc::display( N2, N2, AIm.data(), N2, "\nMatrix AIm", std::cout );
#endif
  return spectral_bound( N2, AIm.data() );
}
*/
template <class T> inline std::ostream&
operator<<
( std::ostream &out, const Specbnd<T> &y )
{
  out << y.SI();
  //out << "  " << y._SI << std::endl
  //    << "  " << y._FI.val() << std::endl;
  //for( unsigned int i=0; i<y._n; i++ )
  //  out << "  ( " << y._FI.deriv(i) << " )" << std::endl;
  return out;
}

template <class T> inline Specbnd<T>&
Specbnd<T>::operator=
( double const& c )
{ 
  _n  = 0;
  _D.clear();
  _N.clear();
  _FI = c;
  _SI = 0.;
  return *this;
}

template <class T> inline Specbnd<T>&
Specbnd<T>::operator=
( T const& I )
{ 
  _n  = 0;
  _D.clear();
  _N.clear();
  _FI = I;
  _SI = 0.;
  return *this;
}
/*
template <class T> inline Specbnd<T>&
Specbnd<T>::set
( std::set<unsigned> const& D, fadbad::F<T> const& FB, T const& SB )
{ 
  _n  = FB.size();
  _D  = D;
  _FI = FB;
  _SI = SB;
  return *this;
}
*/
template <class T> inline Specbnd<T>&
Specbnd<T>::operator=
( Specbnd<T> const& x )
{ 
  _n  = x._n;
  _D  = x._D;
  _N  = x._N;
  _FI = x._FI;
  _SI = x._SI;
  return *this;
}

template <class T> inline Specbnd<T>&
Specbnd<T>::operator=
( Specbnd<T> && x )
{
  _n  = x._n;
  _D  = std::move( x._D );
  _N  = std::move( x._N );
  _FI = std::move( x._FI );
  _SI = std::move( x._SI );
  return *this;
}

template <class T> inline Specbnd<T>&
Specbnd<T>::operator+=
( double const& c )
{ 
  _FI += c;
  return *this;
}

template <class T> inline Specbnd<T>&
Specbnd<T>::operator+=
( Specbnd<T> const& y )
{ 
  Specbnd::add( *this, y );
  return *this;
}

template <class T> inline Specbnd<T>
operator+
( Specbnd<T> const& y )
{
  return y;
}

template <class T> Specbnd<T>
operator+
( Specbnd<T> const& x, Specbnd<T> const& y )
{
  Specbnd<T> z( x );
  return( z += y );
}

template <class T> Specbnd<T> &&
operator+
( Specbnd<T> && x, Specbnd<T> const& y )
{
  return( std::move( x += y ) );
}

template <class T> Specbnd<T> &&
operator+
( Specbnd<T> const& x, Specbnd<T> && y )
{
  return( std::move( y += x ) );
}

template <class T> Specbnd<T> &&
operator+
( Specbnd<T> && x, Specbnd<T> && y )
{
  return( std::move( x += y ) );
}

template <class T> Specbnd<T>
operator+
( Specbnd<T> const& y, double const& c )
{
  Specbnd<T> z( y );
  return( z += c );
}

template <class T> Specbnd<T> &&
operator+
( Specbnd<T> && y, double const& c )
{
  return( std::move( y += c ) );
}

template <class T> Specbnd<T>
operator+
( const double& c, const Specbnd<T> &y )
{
  Specbnd<T> z( y );
  return( z += c );
}

template <class T> Specbnd<T> &&
operator+
( double const& c, Specbnd<T> && y )
{
  return( std::move( y += c ) );
}

template <class T> inline Specbnd<T>
operator-
( Specbnd<T> const& y )
{
  Specbnd<T> z( y );
  Specbnd<T>::neg( z );
  return z;
}

template <class T> inline Specbnd<T> &&
operator-
( Specbnd<T> && y )
{
  Specbnd<T>::neg( y );
  return std::move( y );
}

template <class T> inline Specbnd<T>&
Specbnd<T>::operator-=
( double const& c )
{ 
  _FI -= c;
  return *this;
}

template <class T> inline Specbnd<T>&
Specbnd<T>::operator-=
( Specbnd<T> const& y )
{
  Specbnd::add( *this, std::move( operator-( y ) ) );
  return *this;
}

template <class T> inline Specbnd<T>
operator-
( Specbnd<T> const& x, Specbnd<T> const& y )
{
  Specbnd<T> z( x );
  return( z -= y );
}

template <class T> inline Specbnd<T> &&
operator-
( Specbnd<T> && x, Specbnd<T> const& y )
{
  return( std::move( x -= y ) );
}

template <class T> inline Specbnd<T> &&
operator-
( Specbnd<T> const& x, Specbnd<T> && y )
{
  Specbnd<T>::neg( y -= x );
  return( std::move( y ) );
}

template <class T> inline Specbnd<T> &&
operator-
( Specbnd<T> && x, Specbnd<T> && y )
{
  return( std::move( x -= y ) );
}

template <class T> inline Specbnd<T>
operator-
( Specbnd<T> const& y, double const& c )
{
  Specbnd<T> z( y );
  return z -= c;
}

template <class T> inline Specbnd<T> &&
operator-
( Specbnd<T> && y, double const& c )
{
  return( std::move( y -= c ) );
}

template <class T> inline Specbnd<T>
operator-
( double const& c, Specbnd<T> const& y )
{
  Specbnd<T> z( y );
  Specbnd<T>::neg( z );
  return( z += c );
}

template <class T> inline Specbnd<T> &&
operator-
( double const& c, Specbnd<T> && y )
{
  Specbnd<T>::neg( y );
  return( std::move( y += c ) );
}

template <typename T> inline Specbnd<T>&
Specbnd<T>::operator*=
( double const& c )
{
  _FI *= c;
  _SI *= c;
  return *this;
}

template <class T> inline Specbnd<T>&
Specbnd<T>::operator*=
( Specbnd<T> const& y )
{ 
  Specbnd::multiply( *this, y );
  return *this;
}

template <class T> inline Specbnd<T>
operator*
( Specbnd<T> const& x, Specbnd<T> const& y )
{
  Specbnd<T> z( x );
  return( z *= y );
}

template <class T> inline Specbnd<T> &&
operator*
( Specbnd<T> && x, Specbnd<T> const& y )
{
  return( std::move( x *= y ) );
}

template <class T> inline Specbnd<T> &&
operator*
( Specbnd<T> && x, Specbnd<T> && y )
{
  return( std::move( x *= y ) );
}

template <class T> inline Specbnd<T> &&
operator*
( Specbnd<T> const& x, Specbnd<T> && y )
{
  return( std::move( y *= x ) );
}

template <class T> inline Specbnd<T>
operator*
( Specbnd<T> const& y, double const& c )
{
  Specbnd<T> z( y );
  return( z *= c );
}

template <class T> inline Specbnd<T> &&
operator*
( Specbnd<T> && y, double const& c )
{
  return( std::move( y *= c ) );
}

template <class T> inline Specbnd<T>
operator*
( double const& c, Specbnd<T> const& y )
{
  Specbnd<T> z( y );
  return( z *= c );
}

template <class T> inline Specbnd<T> &&
operator*
( double const& c, Specbnd<T> && y )
{
  return( std::move( y *= c ) );
}

template <class T> inline Specbnd<T> &&
sqr
( Specbnd<T> && y )
{
  Specbnd<T>::sqr( y );
  return std::move( y );
}

template <class T> inline Specbnd<T>
sqr
( Specbnd<T> const& y )
{
  Specbnd<T> z( y );
  Specbnd<T>::sqr( z );
  return z;
}

template <class T> inline Specbnd<T> &&
pow
( Specbnd<T> && y, const int m )
{
  Specbnd<T>::ipow( y, m );
  return std::move( y );
}

template <class T> inline Specbnd<T>
pow
( const Specbnd<T> &y, const int m )
{
  Specbnd<T> z( y );
  Specbnd<T>::ipow( z, m );
  return z;
}

template <class T> inline Specbnd<T> &&
inv
( Specbnd<T> && y )
{
  Specbnd<T>::inv( y );
  return std::move( y );
}

template <class T> inline Specbnd<T>
inv
( Specbnd<T> const& y )
{
  Specbnd<T> z( y );
  Specbnd<T>::inv( z );
  return z;
}

template <typename T> inline Specbnd<T>&
Specbnd<T>::operator/=
( double const& c )
{
  _FI /= c;
  _SI /= c;
  return *this;
}

template <typename T> inline Specbnd<T>&
Specbnd<T>::operator/=
( Specbnd<T> const& x )
{
  Specbnd<T> y( x );
  Specbnd<T>::inv( y );
  return( operator*=( y ) );
}

template <typename T> inline Specbnd<T>&
Specbnd<T>::operator/=
( Specbnd<T> && x )
{
  Specbnd<T>::inv( x );
  return( operator*=( std::move( x ) ) );
}

template <typename T>
inline
Specbnd<T> operator/
( Specbnd<T> const& x, Specbnd<T> const& y )
{
  Specbnd<T> z( x );
  z /= y;
  return z;
}

template <typename T>
inline
Specbnd<T> && operator/
( Specbnd<T> const& x, Specbnd<T> && y )
{
  Specbnd<T> z( x );
  z /= std::move( y );
  return z;
}

template <typename T>
inline
Specbnd<T> && operator/
( Specbnd<T> && x, Specbnd<T> const& y )
{
  x /= y;
  return std::move( x );
}

template <typename T>
inline
Specbnd<T> && operator/
( Specbnd<T> && x, Specbnd<T> && y )
{
  x /= std::move( y );
  return std::move( x );
}

template <typename T>
inline
Specbnd<T> operator/
( Specbnd<T> const& x, double const& c )
{
  Specbnd<T> z( x );
  z /= c;
  return z;
}

template <typename T>
inline
Specbnd<T> && operator/
( Specbnd<T> && x, double const& c )
{
  x /= c;
  return std::move( x );
}

template <typename T>
inline
Specbnd<T> operator/
( double const& c, Specbnd<T> const& y )
{
  Specbnd<T> z( inv( y ) );
  z *= c;
  return z;
}

template <typename T>
inline
Specbnd<T> && operator/
( double const& c, Specbnd<T> && y )
{
  inv( std::move( y ) );
  y *= c;
  return std::move( y );
}

template <class T> inline Specbnd<T> &&
sqrt
( Specbnd<T> && y )
{
  Specbnd<T>::sqrt( y );
  return std::move( y );
}

template <class T> inline Specbnd<T>
sqrt
( Specbnd<T> const& y )
{
  Specbnd<T> z( y );
  Specbnd<T>::sqrt( z );
  return z;
}

template <class T> inline Specbnd<T> &&
exp
( Specbnd<T> && y )
{
  Specbnd<T>::exp( y );
  return std::move( y );
}

template <class T> inline Specbnd<T>
exp
( Specbnd<T> const& y )
{
  Specbnd<T> z( y );
  Specbnd<T>::exp( z );
  return z;
}

template <class T> inline Specbnd<T> &&
log
( Specbnd<T> && y )
{
  Specbnd<T>::log( y );
  return std::move( y );
}

template <class T> inline Specbnd<T>
log
( Specbnd<T> const& y )
{
  Specbnd<T> z( y );
  Specbnd<T>::log( z );
  return z;
}

template <class T> inline Specbnd<T> &&
xlog
( Specbnd<T> && y )
{
  auto const& f   = [&]( fadbad::F<T> const& x ){ return x*mc::Op<fadbad::F<T>>::log(x); };
  auto const& Df  = [&]( T const& x ){ return 1.+Op<T>::log(x); };
  auto const& D2f = [&]( T const& x ){ return Op<T>::inv(x); };
  Specbnd<T>::compose( y, f, Df, D2f );
  return std::move( y );
}

template <class T> inline Specbnd<T>
xlog
( const Specbnd<T> &y )
{
  Specbnd<T> z( y );
  return xlog( std::move( z ) );
}

template <class T> inline Specbnd<T>
pow
( Specbnd<T> const& x, Specbnd<T> const& y )
{
  return exp( y * log( x ) );
}

template <class T> inline Specbnd<T> &&
pow
( Specbnd<T> && x, Specbnd<T> const& y )
{
  return std::move( exp( std::move( log( std::move( x ) ) ) * y ) );
}

template <class T> inline Specbnd<T>
pow
( double const& c, Specbnd<T> const& y )
{
  return exp( y * std::log( c ) );
}

template <class T> inline Specbnd<T> &&
pow
( double const& c, Specbnd<T> &&y )
{
  return exp( std::move( y ) * std::log( c ) );
}

template <class T> inline Specbnd<T>
pow
( Specbnd<T> const& x, double const& c )
{
  return exp( c * log( x ) );
}

template <class T> inline Specbnd<T> &&
pow
( Specbnd<T> && x, double const& c )
{
  return exp( c * log( std::move( x ) ) );
}

template <class T> inline Specbnd<T>
prod
( const unsigned int n, const Specbnd<T>*x )
{
  switch( n ){
   case 0:  return 1.;
   case 1:  return x[0];
   default: return x[0] * prod( n-1, x+1 );
  }
}

template <class T> inline Specbnd<T>
monom
( const unsigned int n, const Specbnd<T>*x, const unsigned*k )
{
  switch( n ){
   case 0:  return 1.;
   case 1:  return pow( x[0], (int)k[0] );
   default: return pow( x[0], (int)k[0] ) * monom( n-1, x+1, k+1 );
  }
}

template <typename T> inline Specbnd<T>
cheb
( const Specbnd<T> &x, const unsigned n )
{
  switch( n ){
    case 0:  return 1.;
    case 1:  return x;
    default: break;
  }
  return 2.*(x*cheb(x,n-1))-cheb(x,n-2);
}
/*
template <class T> inline Specbnd<T>
fabs
( const Specbnd<T> &y )
{
  throw typename Specbnd<T>::Exceptions( Specbnd<T>::Exceptions::UNDEF );
}
*/
template <class T> inline Specbnd<T> &&
cos
( Specbnd<T> && y )
{
  auto const& f   = [&]( fadbad::F<T> const& x ){ return Op<fadbad::F<T>>::cos(x); };
  auto const& Df  = [&]( T const& x ) -> T { return -Op<T>::sin(x); };
  auto const& D2f = [&]( T const& x ) -> T { return -Op<T>::cos(x); };
  Specbnd<T>::compose( y, f, Df, D2f );
  return std::move( y );
}

template <class T> inline Specbnd<T>
cos
( const Specbnd<T> &y )
{
  Specbnd<T> z( y );
  return cos( std::move( z ) );
}

template <class T> inline Specbnd<T> &&
sin
( Specbnd<T> && y )
{
  auto const& f   = [&]( fadbad::F<T> const& x ){ return Op<fadbad::F<T>>::sin(x); };
  auto const& Df  = [&]( T const& x ) -> T { return  Op<T>::cos(x); };
  auto const& D2f = [&]( T const& x ) -> T { return -Op<T>::sin(x); };
  Specbnd<T>::compose( y, f, Df, D2f );
  return std::move( y );
}

template <class T> inline Specbnd<T>
sin
( const Specbnd<T> &y )
{
  Specbnd<T> z( y );
  return sin( std::move( z ) );
}

template <class T> inline Specbnd<T> &&
tan
( Specbnd<T> && y )
{
  Specbnd<T>::tan( y );
  return std::move( y );
}

template <class T> inline Specbnd<T>
tan
( Specbnd<T> const& y )
{
  Specbnd<T> z( y );
  Specbnd<T>::tan( z );
  return z;
}

template <class T> inline Specbnd<T> &&
acos
( Specbnd<T> && y )
{
  Specbnd<T>::acos( y );
  return std::move( y );
}

template <class T> inline Specbnd<T>
acos
( Specbnd<T> const& y )
{
  Specbnd<T> z( y );
  Specbnd<T>::acos( z );
  return z;
}

template <class T> inline Specbnd<T> &&
asin
( Specbnd<T> && y )
{
  Specbnd<T>::asin( y );
  return std::move( y );
}

template <class T> inline Specbnd<T>
asin
( Specbnd<T> const& y )
{
  Specbnd<T> z( y );
  Specbnd<T>::asin( z );
  return z;
}

template <class T> inline Specbnd<T> &&
atan
( Specbnd<T> && y )
{
  Specbnd<T>::atan( y );
  return std::move( y );
}

template <class T> inline Specbnd<T>
atan
( Specbnd<T> const& y )
{
  Specbnd<T> z( y );
  Specbnd<T>::atan( z );
  return z;
}

template <class T> inline Specbnd<T> &&
cosh
( Specbnd<T> && y )
{
  auto const& f   = [&]( fadbad::F<T> const& x ){ return Op<fadbad::F<T>>::cosh(x); };
  auto const& Df  = [&]( T const& x ){ return Op<T>::sinh(x); };
  auto const& D2f = [&]( T const& x ){ return Op<T>::cosh(x); };
  Specbnd<T>::compose( y, f, Df, D2f );
  return std::move( y );
}

template <class T> inline Specbnd<T>
cosh
( const Specbnd<T> &y )
{
  Specbnd<T> z( y );
  return cosh( std::move( z ) );
}

template <class T> inline Specbnd<T> &&
sinh
( Specbnd<T> && y )
{
  auto const& f   = [&]( fadbad::F<T> const& x ){ return Op<fadbad::F<T>>::sinh(x); };
  auto const& Df  = [&]( T const& x ){ return Op<T>::cosh(x); };
  auto const& D2f = [&]( T const& x ){ return Op<T>::sinh(x); };
  Specbnd<T>::compose( y, f, Df, D2f );
  return std::move( y );
}

template <class T> inline Specbnd<T>
sinh
( const Specbnd<T> &y )
{
  Specbnd<T> z( y );
  return sinh( std::move( z ) );
}

template <class T> inline Specbnd<T> &&
tanh
( Specbnd<T> && y )
{
  Specbnd<T>::tanh( y );
  return std::move( y );
}

template <class T> inline Specbnd<T>
tanh
( Specbnd<T> const& y )
{
  Specbnd<T> z( y );
  Specbnd<T>::tanh( z );
  return z;
}

template <class T> inline Specbnd<T> &&
erf
( Specbnd<T> && y )
{
  auto const& f   = [&]( fadbad::F<T> const& x ){ return mc::Op<fadbad::F<T>>::erf(x); };
  auto const& Df  = [&]( T const& x ){ return (2./std::sqrt(PI))*Op<T>::exp(-Op<T>::sqr(x)); };
  auto const& D2f = [&]( T const& x ){ return (-4./std::sqrt(PI))*x*Op<T>::exp(-Op<T>::sqr(x)); };
  Specbnd<T>::compose( y, f, Df, D2f );
  return std::move( y );
}

template <class T> inline Specbnd<T>
erf
( Specbnd<T> const& y )
{
  Specbnd<T> z( y );
  return erf( std::move( z ) );
}

template <class T> inline Specbnd<T> &&
erfc
( Specbnd<T> && y )
{
  return 1.-erf(y);
}

template <class T> inline Specbnd<T>
erfc
( Specbnd<T> const& y )
{
  return 1.-erf(y);
}

} // namespace mc

namespace fadbad
{

//! @brief Specialization of the structure fadbad::Op for use of the type mc::Specbnd of MC++ as a template parameter of the classes fadbad::F, fadbad::B and fadbad::T of FADBAD++
template< typename T > struct Op< mc::Specbnd<T> >
{ 
  typedef mc::Specbnd<T> SB;
  typedef double Base;
  static Base myInteger( const int i ) { return Base(i); }
  static Base myZero() { return myInteger(0); }
  static Base myOne() { return myInteger(1);}
  static Base myTwo() { return myInteger(2); }
  static double myPI() { return mc::PI; }
  static SB myPos( const SB& x ) { return  x; }
  static SB myNeg( const SB& x ) { return -x; }
  template <typename U> static SB& myCadd( SB& x, const U& y ) { return x+=y; }
  template <typename U> static SB& myCsub( SB& x, const U& y ) { return x-=y; }
  template <typename U> static SB& myCmul( SB& x, const U& y ) { return x*=y; }
  template <typename U> static SB& myCdiv( SB& x, const U& y ) { return x/=y; }
  static SB myInv( const SB& x ) { return mc::inv( x ); }
  static SB mySqr( const SB& x ) { return mc::sqr( x ); }
  template <typename X, typename Y> static SB myPow( const X& x, const Y& y ) { return mc::pow( x, y ); }
  static SB myCheb( const SB& x, const unsigned n ) { return mc::cheb( x, n ); }
  static SB mySqrt( const SB& x ) { return mc::sqrt( x ); }
  static SB myLog( const SB& x ) { return mc::log( x ); }
  static SB myExp( const SB& x ) { return mc::exp( x ); }
  static SB mySin( const SB& x ) { return mc::sin( x ); }
  static SB myCos( const SB& x ) { return mc::cos( x ); }
  static SB myTan( const SB& x ) { return mc::tan( x ); }
  static SB myAsin( const SB& x ) { return mc::asin( x ); }
  static SB myAcos( const SB& x ) { return mc::acos( x ); }
  static SB myAtan( const SB& x ) { return mc::atan( x ); }
  static SB mySinh( const SB& x ) { return mc::sinh( x ); }
  static SB myCosh( const SB& x ) { return mc::cosh( x ); }
  static SB myTanh( const SB& x ) { return mc::tanh( x ); }
  static bool myEq( const SB& x, const SB& y ) { return mc::Op<T>::eq(x.SI(),y.SI()); } 
  static bool myNe( const SB& x, const SB& y ) { return mc::Op<T>::ne(x.SI(),y.SI()); }
  static bool myLt( const SB& x, const SB& y ) { return mc::Op<T>::lt(x.SI(),y.SI()); }
  static bool myLe( const SB& x, const SB& y ) { return mc::Op<T>::le(x.SI(),y.SI()); }
  static bool myGt( const SB& x, const SB& y ) { return mc::Op<T>::gt(x.SI(),y.SI()); }
  static bool myGe( const SB& x, const SB& y ) { return mc::Op<T>::ge(x.SI(),y.SI()); }
};

} // end namespace fadbad

//#include "mcop.hpp"

namespace mc
{

//! @brief C++ structure for specialization of the mc::Op templated structure for use of mc::Specbnd in other MC++ classes
template< typename T > struct Op< mc::Specbnd<T> >
{
  typedef mc::Specbnd<T> SB;
  static SB point( const double c ) { return SB(c); }
  static SB zeroone() { return SB( mc::Op<T>::zeroone() ); }
  static void I(SB& x, const SB&y) { x = y; }
  static double l(const SB& x) { return mc::Op<T>::l(x.I()); }
  static double u(const SB& x) { return mc::Op<T>::u(x.I()); }
  static double abs (const SB& x) { return mc::Op<T>::abs(x.I());  }
  static double mid (const SB& x) { return mc::Op<T>::mid(x.I());  }
  static double diam(const SB& x) { return mc::Op<T>::diam(x.I()); }
  static SB inv (const SB& x) { return mc::inv(x);  }
  static SB sqr (const SB& x) { return mc::sqr(x);  }
  static SB sqrt(const SB& x) { return mc::sqrt(x); }
  static SB exp (const SB& x) { return mc::exp(x);  }
  static SB log (const SB& x) { return mc::log(x);  }
  static SB xlog(const SB& x) { return x*mc::log(x); }
  static SB lmtd(const SB& x, const SB& y) { return (x-y)/(mc::log(x)-mc::log(y)); }
  static SB rlmtd(const SB& x, const SB& y) { return (mc::log(x)-mc::log(y))/(x-y); }
  static SB fabs(const SB& x) { throw typename mc::Specbnd<T>::Exceptions( Specbnd<T>::Exceptions::UNDEF ); }
  static SB sin (const SB& x) { return mc::sin(x);  }
  static SB cos (const SB& x) { return mc::cos(x);  }
  static SB tan (const SB& x) { return mc::tan(x);  }
  static SB asin(const SB& x) { return mc::asin(x); }
  static SB acos(const SB& x) { return mc::acos(x); }
  static SB atan(const SB& x) { return mc::atan(x); }
  static SB sinh(const SB& x) { return mc::sinh(x); }
  static SB cosh(const SB& x) { return mc::cosh(x); }
  static SB tanh(const SB& x) { return mc::tanh(x); }
  static SB erf (const SB& x) { return mc::erf(x); }
  static SB erfc(const SB& x) { return mc::erfc(x); }
  static SB fstep(const SB& x) { throw typename mc::Specbnd<T>::Exceptions( Specbnd<T>::Exceptions::UNDEF ); }
  static SB bstep(const SB& x) { throw typename mc::Specbnd<T>::Exceptions( Specbnd<T>::Exceptions::UNDEF ); }
  static SB hull(const SB& x, const SB& y) { throw typename mc::Specbnd<T>::Exceptions( Specbnd<T>::Exceptions::UNDEF ); }
  static SB min (const SB& x, const SB& y) { throw typename mc::Specbnd<T>::Exceptions( Specbnd<T>::Exceptions::UNDEF ); }
  static SB max (const SB& x, const SB& y) { throw typename mc::Specbnd<T>::Exceptions( Specbnd<T>::Exceptions::UNDEF ); }
  template <typename X, typename Y> static SB pow(const X& x, const Y& y) { return mc::pow(x,y); }
  static SB cheb (const SB& x, const unsigned n) { return mc::cheb(x,n); }
  static SB prod (const unsigned int n, const SB* x) { return mc::prod(n,x); }
  static SB monom (const unsigned int n, const SB* x, const unsigned* k) { return mc::monom(n,x,k); }
  static bool inter(SB& xIy, const SB& x, const SB& y) { throw typename mc::Specbnd<T>::Exceptions( Specbnd<T>::Exceptions::UNDEF ); }
  static bool eq(const SB& x, const SB& y) { return x.SI()==y.SI() && x.FI()==y.FI(); }
  static bool ne(const SB& x, const SB& y) { return x.SI()!=y.SI() || x.FI()==y.FI(); }
  static bool lt(const SB& x, const SB& y) { return x.SI()<y.SI()  && x.FI()<y.FI();  }
  static bool le(const SB& x, const SB& y) { return x.SI()<=y.SI() && x.FI()<=y.FI(); }
  static bool gt(const SB& x, const SB& y) { return x.SI()>y.SI()  && x.FI()>y.FI();  }
  static bool ge(const SB& x, const SB& y) { return x.SI()>=y.SI() && x.FI()>=y.FI(); }
};

} // namespace mc

#endif

