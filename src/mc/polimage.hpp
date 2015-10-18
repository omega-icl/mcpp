// Copyright (C) 2009-2014 Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the Eclipse Public License.

/*!
\page page_POLYTOPE Polyhedra and Polytope Arithmetic for Factorable Functions
\author Beno&icirc;t Chachuat

An ellipsoid with center \f$c \in \mathbb R^n\f$ and shape matrix \f$Q \in \mathbb S_{+}^n\f$ is defined as 
\f{align*}
  \mathcal E(c,Q) := & \left\{ \left. c + Q^\frac{1}{2} v \ \right| \ \exists v \in \mathbb R^{n}: \ v^T v \leq 1 \right\} \subseteq \mathbb R^{n} .
\f}

As well as constructors and data access/manipulations functions, the class mc::Ellipsoid provides a set of functions for ellipsoidal calculus. Data manipulation functions include:
- checking positive semi-definiteness of the shape matrix (mc::Ellipsoid::psdQ)
- determining the rank of or regularizing the shape matrix (mc::Ellipsoid::rankQ, mc::Ellipsoid::regQ)
- computing the eigenvalues or singular values of the shape matrix (mc::Ellipsoid::eigQ, mc::Ellipsoid::svdQ)
- taking the trace, square root or inverse of the shape matrix (mc::Ellipsoid::trQ, mc::Ellipsoid::sqrtQ, mc::Ellipsoid::invQ)
.
Ellipsoidal calculus includes:
- applying a linear transformation to an ellipsoid (mc::mtimes)
- taking the (exact) intersection of an ellipsoid with a hyperplane (mc::hpintersection)
- computing a (minimum volume) external ellipsoidal approximation of the intersection between an ellipsoid and a halfspace (mc::intersection_ea)
- computing an external ellipsoidal approximation of the geometric (Minkowski) sum of several ellipsoids along a given direction, or a (minimum trace) external ellipsoidal approximation of the geometric (Minkowski) sum of an ellipsoid with an interval box (mc::minksum_ea)
.

Besides ellipsoidal calculus, the classes mc::PolImg and mc::PolVar provide an implementation of ellipsoidal arithmetic in order to enclose the image \f$\mathcal E(c_f,Q_f)\f$ of an \f$n_x\f$-dimensional ellispoid \f$\mathcal E(c_x,Q_x)\f$ under a vector-valued function \f$ f:\mathbb{R}^{n_x}\to\mathbb{R}^{n_f} \f$:
\f{align*}
  \mathcal E(c_f,Q_f) \supseteq & \left\{ f(x) \,\mid\, x\in \mathcal{E}(c_{x},Q_x) \right\}.
\f}
Notice that the exact image \f$\left\{ f(x) \,\mid\, x\in \mathcal{E}(c_{x},Q_x) \right\}\f$ is not an ellipsoid in general.

The class mc::PolImg is derived from mc::Ellipsoid. The implementation of mc::PolImg and mc::PolVar relies on the operator/function overloading mechanism of C++. This makes the computation of the ellipsoidal enclosure for the image of an ellipsoid under a factorable function both simple and intuitive, similar to computing function values in real arithmetic or bounds based on interval, Taylor or Chebyshev model arithmetics (see \ref page_INTERVAL, \ref page_TAYLOR, \ref page_CHEBYSHEV). mc::PolImg stores a column vector CPPL::dcovector and a sparse symmetric matrix CPPL::dssmatrix provided by the LAPACK wrapper <A href="http://cpplapack.sourceforge.net/">CPPLAPACK</A>. The column vector stores the center and the sparse symmetric matrix the shape of a <i>lifted</i> ellipsoid, which is the result of adding extra dimension for each operation participating in the factorable function. Note that the implementation in mc::PolImg and mc::PolVar is <a>not verified</a> in the sense that rounding errors are not accounted for during the propagation.

The classes mc::PolImg and mc::PolVar are templated in the interval type used to bound the nonlinearity of the function, By default, mc::PolImg and mc::PolVar can be used with the non-verified interval type mc::Interval of MC++. For reliability, however, it is recommended to use verified interval arithmetic such as <A href="http://www.ti3.tu-harburg.de/Software/PROFILEnglisch.html">PROFIL</A> (header file <tt>mcprofil.hpp</tt>) or <A href="http://www.math.uni-wuppertal.de/~xsc/software/filib.html">FILIB++</A> (header file <tt>mcfilib.hpp</tt>). 

\section sec_ELLCALC How do I define an ellipsoid and apply ellipsoidal calculus?

In order to define the ellipsoid \f$\mathcal E(c_x,Q_x)\f$ with
\f{align*}
  c_x = & \left(\begin{array}{c} 3\\4\end{array}\right),\ \text{and} & 
  Q_x = & \left(\begin{array}{cc} 5 & 4 \\ 4 & 5 \end{array}\right).
\f}

we proceed as follows:

\code
    const unsigned n = 2;
    CPPL::dcovector cx(n); CPPL::csymatrix Qx(n);
    cx(0) = 3.;  Qx(0,0) = 5.;
    cx(1) = 4.;  Qx(1,0) = 4.;  Qx(1,1) = 5. ;
    mc::Ellipsoid Ex( Qx, cx ); 
\endcode

This ellipsoid is simply displayed as:

\code
    std::cout << Ex << std::endl;
\endcode

In the present case, the following information is displayed:

\verbatim
center:
 3.00000e+00
 4.00000e+00

shape:
 5.00000e+00 {4.00000e+00}
 4.00000e+00  5.00000e+00 
\endverbatim

In order to illustrate ellipsoidal calculus, suppose that we want to add the interval box \f$[0,0.1]^2\f$ to the foregoing ellispoid and determine an external ellispoidal approximation of this geometric sum. 

\code
   TBC
\endcode

\section sec_ELLIMG How do I compute an ellipsoidal enclosure for the image set of an ellipsoid under a factorable function?

Suppose that we want to compute an ellipsoidal enclosure for the image set of function:
\f[
f(x) = \left(\begin{array}{c}
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
v^{\rm T}v\leq 1 \right\}
\f].

For simplicity, the underlying interval bounds are propagated using the default interval type mc::Interval, the required header files 
are:
 
\code
#include "ellimage.hpp"
#include "interval.hpp"

typedef mc::Interval I		;
typedef mc::PolImg<I> EI	;
typedef mc::PolVar<I> EV	;
typedef CPPL::dcovector dcv	;
typedef CPPL::dsymatrix dsm	;
\endcode

First, the host ellipsoidal set for the independent variables \f$x_{1}\f$ and \f$x_{2}\f$ is specified as:
\code 
dcv cx(2)	;     dsm Qx(2)		;
cx(0) = 3.	;     Qx(0,0) = 5.	;
cx(1) = 4.	;     Qx(1,0) = 4.	;	Qx(1,1) = 5.	;
EI Ex( Qx, qx )	;
\endcode

Then, the independent variables themselves are specified as:
\code
EV X1( Ex, 0 ) ;
EV X2( Ex, 1 ) ;
\endcode
If independent interval bounds are known for each variable they can be passed as an optional third argument to the set function.

The dependent variables \f$f_{1}(x)\f$ and \f$f_{2}(c)\f$ are propagated as:
\code
EV F[2] = { log( X[0] ) + mc::sqr( X[1] ),
            sin( X[0] ) - cos( X[1] ) } ;
\endcode

and the lifted ellipsoid can be displayed as:

\code
std::cout << "lifted ellipsoidal image of f =";
Ex.output();
\endcode

\verbatim
lifted ellipsoidal image of f = 
center: 
 3
 4
 18.5
 0.913697
 19.4137
 -0.0590929
 0.012758
 0.0718509
Shape: 
 9.46133 {7.56906}{60.5525}{4.07224}{64.6247}{2.82094}{-4.61268}{-7.43362}
 7.56906  9.46133 {75.6906}{3.25779}{78.9484}{3.52618}{-3.69015}{-7.21632}
 60.5525  75.6906  752.099 {26.0623}{778.161}{28.2094}{-29.5212}{-57.7306}
 4.07224  3.25779  26.0623  2.52547 {28.5878}{1.21416}{-1.98534}{-3.1995}
 64.6247  78.9484  778.161  28.5878  806.749 {29.4236}{-31.5065}{-60.9301}
 2.82094  3.52618  28.2094  1.21416  29.4236  4.40317 {-1.37529}{-5.77847}
 -4.61268  -3.69015  -29.5212  -1.98534  -31.5065  -1.37529  4.34584 {5.72113}
 -7.43362  -7.21632  -57.7306  -3.1995  -60.9301  -5.77847  5.72113  11.4996 
\endverbatim

Finally, the lifted ellipsoid can be projected in the space of dependent variables to obtain the desired image enclosure:
\code
EI Ef = Eflift.get( 2, F );
std::cout << " Ellipsoidal enclosure Ef = " << Ef << std::endl;
\endcode

\verbatim
Ellipsoidal enclosure Ef =
center:
 1.94137e+01
 7.18509e-02
shape:
 8.06749e+02 {-6.09301e+01}
 -6.09301e+01  1.14996e+01 
\endverbatim

After this, the ellipsoid can be manipulated according to the rules of ellipsoidal calculus (see \ref page_POLYTOPE).
A comparison between the exact image and the ellipsoidal enclosure is presented in the following figure, in orange the exact image and in blue the boundary of the ellipsoidal enclosure.
<CENTER><TABLE BORDER=0>
<TR>
<TD>\image html ELL-2D.png</TD>
</TR>
</TABLE></CENTER>


\section sec_ELL_opt How are the options set for ellipsoidal calculus and ellipsoidal arithmetic?

The class mc::PolImg and mc::Ellipsoid have public members called mc::PolImg::options and mc::Ellipsoid::options (static), respectively, that can be used to set/modify a number of options. Note that mc::PolImg::options is a superset of mc::Ellipsoid::options since mc::PolImg is derived from mc::Ellispoid. For instance, options can be set as follows:
\code
	Ex.options.PREALLOC = 5; 
	Ex.options.CHEBUSE  = false;
\endcode

The full set of available options is reported in the following tables.

<TABLE border="1">
<CAPTION><EM>Options in mc::PolImg::Options: name, type and description</EM></CAPTION>
     <TR><TH><b>Name</b>  <TD><b>Type</b><TD><b>Default</b>
         <TD><b>Description</b>
     <TR><TH><tt>PREALLOC</tt> <TD><tt>unsigned long</tt> <TD>0
         <TD> Number of rows to preallocate in the shape matrix and center vector
     <TR><TH><tt>CHEBUSE</tt> <TD><tt>bool</tt> <TD>false
         <TD> Whether to use Chebyshev expansion to compute a linear approximation and bound the nonlinear dependencies of univariate terms
     <TR><TH><tt>CHEBORDER</tt> <TD><tt>unsigned</tt> <TD>5
         <TD>Order of the Chebyshev expansion (only of <tt>CHEBUSE = true</tt>)
</TABLE>

<TABLE border="1">
<CAPTION><EM>Options in mc::Ellipsoid::Options: name, type and description</EM></CAPTION>
     <TR><TH><b>Name</b>  <TD><b>Type</b><TD><b>Default</b>
         <TD><b>Description</b>
     <TR><TH><tt>PSDCHK</tt> <TD><tt>bool</tt> <TD>false
         <TD>Whether or not to check positive semi-definiteness of shape matrices
     <TR><TH><tt>PSDTOL</tt> <TD><tt>double</tt> <TD>1e2*MACHPREC
         <TD>Absolute tolerance for positive semi-definiteness check of shape matrices
     <TR><TH><tt>RKTOLA</tt> <TD><tt>double</tt> <TD>MACHPREC
         <TD>Absolute tolerance for rank and regularization of shape matrices
     <TR><TH><tt>RKTOLR</tt> <TD><tt>double</tt> <TD>MACHPREC*1e6
         <TD>Relative tolerance for rank and regularization of shape matrices
     <TR><TH><tt>ROOTTOL</tt> <TD><tt>double</tt> <TD>1e-10
         <TD>Absolute stopping tolerance for root-finding method (objective function value less than ROOTTOL)
     <TR><TH><tt>ROOTSECANT</tt> <TD><tt>bool</tt> <TD>false
         <TD>Whether to use the secant method for root finding
     <TR><TH><tt>ROOTMAXIT</tt> <TD><tt>bool</tt> <TD>0
         <TD>Maximum number of iteration for root-finding method (no maximum when ROOTMAXIT=0)
</TABLE>


\section sec_ELL_err What are the errors encountered during ellipsoidal calculus and ellipsoidal arithmetic?

Errors are managed based on the exception handling mechanism of the C++ language. Each time an error is encountered, a class object of type mc::PolImg::Exceptions or mc::Ellipsoid::Exceptions is thrown, which contains the type of error. It is the user's responsibility to test whether an exception was thrown during the computation of the lifted ellipsoid, and then make the appropriate changes. Should an exception be thrown and not caught by the calling program, the execution will abort.

Possible errors encountered during application of ellipsoidal calculus and ellispoidal arithmetic are reported in the following tables.

<TABLE border="1">
<CAPTION><EM>Errors during the Computation of an Ellipsoidal Image</EM></CAPTION>
     <TR><TH><b>Number</b> <TD><b>Description</b>
     <TR><TH> <tt> 1 </tt>  <TD> Division by zero
     <TR><TH> <tt> 2 </tt>  <TD> Inverse operation with zero in domain
     <TR><TH> <tt> 3 </tt>  <TD> Log operation with non-positive numbers in domain
     <TR><TH> <tt> 4 </tt>  <TD> Square-root operation with negative numbers in domain
     <TR><TH> <tt> 5 </tt>  <TD> Tangent operation with zero in domain of cosine, tan(x) = sin(x)/cos(x)
     <TR><TH> <tt> 6 </tt>  <TD> Sine/Cosine inverse operation with domain outside [-1,1]
     <TR><TH> <tt>-1 </tt>  <TD> Failed to construct ellipsoidal variable
     <TR><TH> <tt>-3 </tt>  <TD> Operation between variables mc::PolVar linked to different images mc::PolImg
     <TR><TH> <tt>-33</tt>  <TD> Feature not yet implemented in mc::PolImg
</TABLE>

<TABLE border="1">
<CAPTION><EM>Errors with mc::Ellipsoid</EM></CAPTION>
     <TR><TH><b>Number</b> <TD><b>Description</b>
     <TR><TH> <tt> 1 </tt> <TD> Non positive-semi definite shape matrix
     <TR><TH> <tt> 2 </tt> <TD> Failure in a LAPACK linear algebra routine
     <TR><TH> <tt> 3 </tt> <TD> Failure in a root-finding routine
</TABLE>

Moreover, exceptions may be thrown by the template parameter class itself.

\section sec_ELL_refs References

Kurzhanskiy, A., A.. and P. Varaiya, <A href="http://code.google.com/p/ellipsoids">"Ellipsoidal Toolbox"</A>, Technical Report UCB/EECS-2006-46, EECS Department, University of California, Berkeley, May 2006.

TODO: 
- Semilinear relaxation of terms that are neither convex nor concave? (e.g., pow(x,3)) ==> OK
- Account for multiple occurence of variables ==> OK
- Add all other univariate terms
- Enable DC-decomposition of bilinear terms (e.g. lists of product/division terms) ==> OK
- Fix bug in DC scaling ==> OK
- Implement alternative piecewise relaxation approach ==> OK
- Implement RLT and reduction constraints
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

#define MC__POLIMG_DEBUG

namespace mc
{

template< class T > class PolImg;
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

  template< class U > friend  PolVar<U> operator+( const PolVar<U>& );
  template< class U > friend  PolVar<U> operator+( const PolVar<U>&, const PolVar<U>& );
  template< class U > friend  PolVar<U> operator+( const PolVar<U>&, const double );
  template< class U > friend  PolVar<U> operator+( const double, const PolVar<U>& );
  template< class U > friend  PolVar<U> operator+( const PolVar<U>&, const U& );
  template< class U > friend  PolVar<U> operator+( const U&, const PolVar<U>& );
  template< class U > friend  PolVar<U> operator-( const PolVar<U>& );
  template< class U > friend  PolVar<U> operator-( const PolVar<U>&, const PolVar<U>& );
  template< class U > friend  PolVar<U> operator-( const PolVar<U>&, const double );
  template< class U > friend  PolVar<U> operator-( const double, const PolVar<U>& );
  template< class U > friend  PolVar<U> operator-( const PolVar<U>&, const U& );
  template< class U > friend  PolVar<U> operator-( const U&, const PolVar<U>& );
  template< class U > friend  PolVar<U> operator*( const PolVar<U>&, const PolVar<U>& );
  template< class U > friend  PolVar<U> operator*( const PolVar<U>&, const double );			
  template< class U > friend  PolVar<U> operator*( const double, const PolVar<U>& );
  template< class U > friend  PolVar<U> operator*( const PolVar<U>&, const U& );
  template< class U > friend  PolVar<U> operator*( const U&, const PolVar<U>& );
  template< class U > friend  PolVar<U> operator/( const PolVar<U>&, const PolVar<U>& );
  template< class U > friend  PolVar<U> operator/( const PolVar<U>&, const double );			
  template< class U > friend  PolVar<U> operator/( const double, const PolVar<U>& );
  template< class U > friend  PolVar<U> operator/( const PolVar<U>&, const U& );
  template< class U > friend  PolVar<U> operator/( const U&, const PolVar<U>& );
  template< class U > friend  PolVar<U> operator^( const PolVar<U>&, const PolVar<U>& );
  template< class U > friend  PolVar<U> operator^( const PolVar<U>&, const double );
  template< class U > friend  PolVar<U> operator^( const double, const PolVar<U>& );
    
  template< class U > friend  PolVar<U> inv ( const PolVar<U>& );
  template< class U > friend  PolVar<U> exp ( const PolVar<U>& );
  template< class U > friend  PolVar<U> log ( const PolVar<U>& );
  template< class U > friend  PolVar<U> sqrt( const PolVar<U>& );
  template< class U > friend  PolVar<U> sqr ( const PolVar<U>& );
  template< class U > friend  PolVar<U> pow ( const PolVar<U>&, const PolVar<U>& );  
  template< class U > friend  PolVar<U> cheb( const PolVar<U>&, const PolVar<U>& );  
  template< class U > friend  PolVar<U> cos ( const PolVar<U>& );
  template< class U > friend  PolVar<U> sin ( const PolVar<U>& );
  template< class U > friend  PolVar<U> tan ( const PolVar<U>& );
  template< class U > friend  PolVar<U> acos( const PolVar<U>& );
  template< class U > friend  PolVar<U> asin( const PolVar<U>& );
  template< class U > friend  PolVar<U> atan( const PolVar<U>& );

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
    //! @brief variable break-points
    std::set<double> _breakpts;
    //! @brief variable subdivision w.r.t. break-points
    mutable std::pair< std::vector<double>, std::vector<PolVar<T> > > _subdiv;
    //! @brief flag indicating whether cuts have already been generated
    mutable bool _hascuts;

    //! @brief set as DAG variable <a>X</a> in polytope image <a>P</a> with range <a>B</a>
    void _set
      ( PolImg<T>*P, FFVar&X, const T&B, const bool cont, const unsigned index )
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
      ( const T&B, const bool cont )
      { _range = B; _id.first = (cont? VARCONT: VARINT); }

    //! @brief push variable subdivision
    void _push_subdiv
      ( double pt ) const
      { if( !_img ) return;
        auto itVar = _img->_Vars.find( const_cast<FFVar*>(&_var) );
        if( itVar == _img->_Vars.end() ) return;
        itVar->second->_subdiv.first.push_back( pt );
        _subdiv.first.push_back( pt ); }
    //! @brief push variable subdivision
    void _push_subdiv
      ( PolVar<T> var ) const
      { if( !_img ) return;
        auto itVar = _img->_Vars.find( const_cast<FFVar*>(&_var) );
        if( itVar == _img->_Vars.end() ) return;
        itVar->second->_subdiv.second.push_back( var );
        _subdiv.second.push_back( var ); }

  public:
    /** @ingroup POLYTOPE
     *  @{
     */
    //! @brief Constructor for a constant value <a>d</a> (default)
    PolVar( const double d=0. )
      : _img(0), _var(FFVar(d)), _range(d), _id( AUXCST, 0 ), _hascuts(false) 
      {} 
    //! @brief Constructor for a constant value <a>n</a>
    PolVar( const int n )
      : _img(0), _var(FFVar(n)), _range(n), _id( AUXCST, 0 ), _hascuts(false)
      {}
    //! @brief Constructor for DAG variable <a>X</a> in polytope image <a>P</a> with range <a>B</a>
    PolVar( PolImg<T>*P, FFVar&X, const T&B=0., const bool cont=true )
      { set( P, X, B, cont ); }
    //! @brief Constructor for auxiliary variable in polytope image <a>P</a> with range <a>B</a>
    PolVar( PolImg<T>*P, const T&B=0., const bool cont=true )
      { set( P, B, cont ); }
    //! @brief Copy constructor for a polytope image <a>P</a>
    PolVar( const PolVar<T>&P )
      : _img(P._img), _var(P._var), _range(P._range), _id( P._id ),
        _breakpts( P._breakpts ), _subdiv( P._subdiv )
      {}

    //! @brief Destructor
    virtual ~PolVar
      ()
      {};

    //! @brief Update range <a>B</a> and type <a>cont</a> of variable in polytope image
    PolVar<T>& update
      ( const T&B=0., const bool cont=true )
      { if( !_img ) return *this;
        *this = *_img->_append_var( &_var, B, cont );
        reset_subdiv(); reset_cuts(); return *this; }
    //! @brief set as DAG variable <a>X</a> in polytope image <a>P</a> with range <a>B</a>
    PolVar<T>& set
      ( PolImg<T>*P, FFVar&X, const T&B=0., const bool cont=true )
      { *this = *P->_append_var( &X, B, cont );
        _breakpts.clear(); reset_subdiv(); reset_cuts(); return *this; }
    //! @brief set as auxiliary variable in polytope image <a>P</a> with range <a>B</a>
    PolVar<T>& set
      ( PolImg<T>*P, const T&B=0., const bool cont=true )
      { *this = *P->_append_aux( B, cont );
        _breakpts.clear(); reset_subdiv(); reset_cuts(); return *this; }

    //! @brief get variable range
    T range
      () const
      { return _range; }
    //! @brief get pointer to polytopic image
    PolImg<T>* image() const
      { return _img; }
    //! @brief get pointer to variable identifier
    t_idVar id() const
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
    const std::pair< std::vector<double>, std::vector<PolImg<T> > >& subdiv
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
    bool cuts() const
      { return _hascuts; }
    //! @brief reset cuts flag to false
    void reset_cuts() const
      { _hascuts = false; }
    //! @brief set cuts flag to true
    void set_cuts() const
      { _hascuts = true; }

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
PolVar<T>&
PolVar<T>::operator=
( const PolVar<T>& P ) 
{
  _img = P._img;
  _var = P._var;
  _range = P._range;
  _id = P._id;
  _breakpts = P._breakpts;
  _subdiv = P._subdiv;
  return *this; 
}

template <class T> 
inline 
PolVar<T>&
PolVar<T>::operator=
( const double d )
{
  _img = 0;
  _var = FFVar(d);
  _range = d;
  _id = std::make_pair( AUXCST, 0 );
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
  _img = 0;
  _var = FFVar(n);
  _range = n;
  _id = std::make_pair( AUXCST, 0 );
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
   double atol = _img->options.BREAKPOINT_ATOL, rtol = _img->options.BREAKPOINT_RTOL;
   if( isequal( Op<T>::l(itVar->second->range()), bkpt, atol, rtol ) 
    || isequal( Op<T>::u(itVar->second->range()), bkpt, atol, rtol ) ) return;
   auto itL = _breakpts.lower_bound( bkpt );
   if( itL!=_breakpts.end() && isequal( *itL, bkpt, atol, rtol ) ) return;
   auto itU = _breakpts.upper_bound( bkpt );
   if( itU!=_breakpts.end() && isequal( *itU, bkpt, atol, rtol ) ) return;
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

  const unsigned nsubint = _subdiv.first.size()-1;
  double coef[nsubint];
  for( unsigned isub=0; isub<nsubint; isub++ ){
    coef[isub] = _subdiv.first[isub] - _subdiv.first[isub+1];
    _push_subdiv( PolVar<T>( _img, Op<T>::zeroone(), true ) ); 
  }
  _img->_append_cut( pOp, PolCut<T>::EQ, _subdiv.first[0], nsubint, _subdiv.second.data(), coef, *this, 1. );

  for( unsigned isub=0; isub<nsubint-1; isub++ ){
    _push_subdiv( PolVar<T>( _img, Op<T>::zeroone(), false ) );
    _img->_append_cut( pOp, PolCut<T>::LE, 0., _subdiv.second[nsubint+isub], 1., _subdiv.second[isub], -1. );
    _img->_append_cut( pOp, PolCut<T>::GE, 0., _subdiv.second[nsubint+isub], 1., _subdiv.second[isub+1], -1. );
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

  const unsigned nsubint = _subdiv.first.size();
  double coef[nsubint];
  for( unsigned isub=0; isub<nsubint; isub++ ){
    coef[isub] = 1.;
    _push_subdiv( PolVar<T>( _img, Op<T>::zeroone(), true ) ); 
  }
  _img->_append_cut( pOp, PolCut<T>::EQ, 1., nsubint, _subdiv.second.data(), coef );

  for( unsigned isub=0; isub<nsubint; isub++ )
    coef[isub] = _subdiv.first[isub];
  _img->_append_cut( pOp, PolCut<T>::EQ, 0., nsubint, _subdiv.second.data(), coef, *this, -1. );
  _img->_append_cut( pOp, PolCut<T>::SOS2, 1., nsubint, _subdiv.second.data(), coef );

  return _subdiv.second;
}

//! @brief C++ class for defining polytopic image variables
////////////////////////////////////////////////////////////////////////
//! mc::PolBilin is a C++ structure for holding information about
//! bilinear terms in a polytopic image. This is useful in connection to
//! DC decomposition in order to refine piecewise-linear cuts.
////////////////////////////////////////////////////////////////////////
template< class T >
struct PolBilin
////////////////////////////////////////////////////////////////////////
{
  //! @brief Constructor for bilinear term
  PolBilin<T>
    ( const PolVar<T>*Var )
    : var(Var), scal1(0), scal2(0)
    {}
  //! @brief corresponding variable in polytopic image (product result)
  const PolVar<T>* var;
  //! @brief scaling constant for sum term in DC decomposition
  PolVar<T>* scal1;
  //! @brief scaling constant for difference term in DC decomposition
  PolVar<T>* scal2;
};

//! @brief C++ template class for defining cuts in the relaxation of factorable programs
////////////////////////////////////////////////////////////////////////
//! mc::PolCut is a C++ template class defining cuts in the relaxation
//! of factorable programs.
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
    SOS2	//!< SOS2-type constraint 
  };
  /** @} */

private:
  //! @brief Pointer to defining operation
  FFOp* _op;
  //! @brief Type of cut
  TYPE _type;
  //! @brief Right-hand side
  double _rhs;
  //! @brief Number of participating variables
  unsigned _nvar;
  //! @brief Participating variables
  PolVar<T>* _var;
  //! @brief Coefficients
  double* _coef;

public:
  /** @ingroup LPRELAX
   *  @{
   */
  //! @brief Retreive type of cut
  TYPE type() const
    { return _type; }
  //! @brief Retreive right-hand side
  double rhs() const
    { return _rhs; }
  //! @brief Retreive number of participating variables
 unsigned nvar() const
    { return _nvar; }
  //! @brief Retreive variable coefficients
  const double* coef() const 
    { return _coef; }
  //! @brief Retreive variable pointers
  const PolVar<T>* var() const 
    { return _var; }
  //! @brief Retreive pointer to corresponding operation in DAG
  FFOp*& op()
    { return _op; }
  //! @brief Retreive pointer to corresponding operation in DAG
  const FFOp* op() const
    { return _op; }

  //! @brief Default constructor for cut
  //PolCut
  //  ( FFOp*op, TYPE type=EQ ):
  //  _op(op), _type(type), _nvar(0), _coef(0), _idvar(0)
  //  {}

  //! @brief Constructor for cut w/ 1 participating variable
  PolCut
    ( FFOp*op, TYPE type, const double b, const PolVar<T>&X1, const double a1 )
    : _op(op), _type(type), _rhs(b)
    {
      _nvar = 1;
      _coef = new double[_nvar];
      _var  = new PolVar<T>[_nvar];
      _coef[0] = a1;
      _var[0] = X1;
    }
  //! @brief Constructor for cut w/ 2 participating variables
  PolCut
    ( FFOp*op, TYPE type, const double b, const PolVar<T>&X1, const double a1,
      const PolVar<T>&X2, const double a2 )
    : _op(op), _type(type), _rhs(b)
    {
      _nvar = 2;
      _coef = new double[_nvar];
      _var = new PolVar<T>[_nvar];
      _coef[0] = a1;
      _coef[1] = a2;
      _var[0] = X1;
      _var[1] = X2;
    }
  //! @brief Constructor for cut w/ 3 participating variables
  PolCut
    ( FFOp*op, TYPE type, const double b, const PolVar<T>&X1, const double a1,
      const PolVar<T>&X2, const double a2, const PolVar<T>&X3, const double a3 )
    : _op(op), _type(type), _rhs(b)
    {
      _nvar = 3;
      _coef = new double[_nvar];
      _var  = new PolVar<T>[_nvar];
      _coef[0] = a1;
      _coef[1] = a2;
      _coef[2] = a3;
      _var[0] = X1;
      _var[1] = X2;
      _var[2] = X3;
    }
  //! @brief Constructor for cut w/ <a>n</a> participating variables
  PolCut
    ( FFOp*op, TYPE type, const double b, const unsigned n, const PolVar<T>*X,
      const double*a ):
    _op(op), _type(type), _rhs(b)
    {
      _nvar = n;
      _coef = new double[_nvar];
      _var  = new PolVar<T>[_nvar];
      for( unsigned ivar=0; ivar<_nvar; ivar++ ){
        _coef[ivar] = a[ivar]; _var[ivar] = X[ivar];
      }
    }
  //! @brief Constructor for cut w/ <a>n</a> participating variables
  PolCut
    ( FFOp*op, TYPE type, const double b, const unsigned n, const PolVar<T>*X,
      const double*a, const PolVar<T>&X1, const double a1 ):
    _op(op), _type(type), _rhs(b)
    {
      _nvar = n+1;
      _coef = new double[_nvar+1];
      _var  = new PolVar<T>[_nvar+1];
      for( unsigned ivar=0; ivar<n; ivar++ ){
        _coef[ivar] = a[ivar]; _var[ivar] = X[ivar];
      }
      _coef[n] = a1;
      _var[n] = X1;
    }
  //! @brief Destructor
  ~PolCut()
    {
      delete[] _coef;
      delete[] _var;
    }
  /** @} */

private:
  //! @brief Private methods to block default compiler methods
  PolCut();
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
      break;

    case PolCut<T>::SOS1: case PolCut<T>::SOS2:
      out << " {";
      for(unsigned k=0; k<cut.nvar(); k++ )
        out << " " << cut._var[k].name();
      out << " }";
  }
  
  switch( cut._type ){
    case PolCut<T>::EQ: out << " = "; break;
    case PolCut<T>::LE: out << " <= "; break;
    case PolCut<T>::GE: out << " >= "; break;
    case PolCut<T>::SOS1: out << " SOS1"; return out;
    case PolCut<T>::SOS2: out << " SOS2"; return out;
  }
  
  out << std::setw(iprec+6) << cut._rhs;
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

//! @brief C++ structure for storing subintervals in an outer approximation
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

//! @brief C++ structure for comparison of subintervals in the outer approximation of a univariate convex/concave function
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

  template< class U > friend  PolVar<U> operator+( const PolVar<U>& );
  template< class U > friend  PolVar<U> operator+( const PolVar<U>&, const PolVar<U>& );
  template< class U > friend  PolVar<U> operator+( const PolVar<U>&, const double );
  template< class U > friend  PolVar<U> operator+( const double, const PolVar<U>& );
  template< class U > friend  PolVar<U> operator+( const PolVar<U>&, const U& );
  template< class U > friend  PolVar<U> operator+( const U&, const PolVar<U>& );
  template< class U > friend  PolVar<U> operator-( const PolVar<U>& );
  template< class U > friend  PolVar<U> operator-( const PolVar<U>&, const PolVar<U>& );
  template< class U > friend  PolVar<U> operator-( const PolVar<U>&, const double );
  template< class U > friend  PolVar<U> operator-( const double, const PolVar<U>& );
  template< class U > friend  PolVar<U> operator-( const PolVar<U>&, const U& );
  template< class U > friend  PolVar<U> operator-( const U&, const PolVar<U>& );
  template< class U > friend  PolVar<U> operator*( const PolVar<U>&, const PolVar<U>& );
  template< class U > friend  PolVar<U> operator*( const PolVar<U>&, const double );			
  template< class U > friend  PolVar<U> operator*( const double, const PolVar<U>& );
  template< class U > friend  PolVar<U> operator*( const PolVar<U>&, const U& );
  template< class U > friend  PolVar<U> operator*( const U&, const PolVar<U>& );
  template< class U > friend  PolVar<U> operator/( const PolVar<U>&, const PolVar<U>& );
  template< class U > friend  PolVar<U> operator/( const PolVar<U>&, const double );			
  template< class U > friend  PolVar<U> operator/( const double, const PolVar<U>& );
  template< class U > friend  PolVar<U> operator/( const PolVar<U>&, const U& );
  template< class U > friend  PolVar<U> operator/( const U&, const PolVar<U>& );
  template< class U > friend  PolVar<U> operator^( const PolVar<U>&, const PolVar<U>& );
  template< class U > friend  PolVar<U> operator^( const PolVar<U>&, const double );
  template< class U > friend  PolVar<U> operator^( const double, const PolVar<U>& );

  template< class U > friend  PolVar<U> inv ( const PolVar<U>& );
  template< class U > friend  PolVar<U> exp ( const PolVar<U>& );
  template< class U > friend  PolVar<U> log ( const PolVar<U>& );
  template< class U > friend  PolVar<U> sqrt( const PolVar<U>& );
  template< class U > friend  PolVar<U> sqr ( const PolVar<U>& );
  template< class U > friend  PolVar<U> pow ( const PolVar<U>&, const PolVar<U>& );  
  template< class U > friend  PolVar<U> cheb( const PolVar<U>&, const PolVar<U>& );  
  template< class U > friend  PolVar<U> cos ( const PolVar<U>& );
  template< class U > friend  PolVar<U> sin ( const PolVar<U>& );
  template< class U > friend  PolVar<U> tan ( const PolVar<U>& );
  template< class U > friend  PolVar<U> acos( const PolVar<U>& );
  template< class U > friend  PolVar<U> asin( const PolVar<U>& );
  template< class U > friend  PolVar<U> atan( const PolVar<U>& );

public:
  typedef std::map< FFVar*, PolVar<T>*, lt_FFVar > t_Vars;
  typedef std::list< PolVar<T>* > t_Aux;
  typedef std::map< const FFVar*, PolBilin<T>*, lt_FFVar > t_Bilin;
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
    ( FFVar*var, const T&range, const bool cont );
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

  //! @brief Map of bilinear terms in polytopic image
  t_Bilin _Bilin;
  //! @brief Appends new pointer to bilinear term map in polytopic image
  PolBilin<T>* _append_bilin
    ( const PolVar<T>*varR );
  //! @brief DC-decompose bilinear term in polytopic image
  void _decompose_bilin
    ( FFOp*op, PolBilin<T>* pBilin, const PolVar<T>&varR, const PolVar<T>&var1,
      const PolVar<T>&var2 );
  //! @brief Erase all entries in _Bilin
  void _erase_bilin
    ();

  //! @brief Set of cuts in polytopic image
  t_Cuts _Cuts;
  //! @brief Erase all entries in _Cuts
  void _erase_cuts
    ();
  //! @brief Erase all entries corresponding to operation <a>op</a> in _Cuts
  void _erase_cuts
    ( FFOp* op );
  //! @brief Appends new relaxation cut in _Cuts w/ 1 variable
  typename t_Cuts::iterator _append_cut
    ( FFOp*op, const typename PolCut<T>::TYPE type,
      const double b, const PolVar<T>&X1, const double a1 );
  //! @brief Appends new relaxation cut in _Cuts w/ 2 variables
  typename t_Cuts::iterator _append_cut
    ( FFOp*op, const typename PolCut<T>::TYPE type,
      const double b, const PolVar<T>&X1, const double a1,
      const PolVar<T>&X2, const double a2 );
  //! @brief Appends new relaxation cut in _Cuts w/ 3 variables
  typename t_Cuts::iterator _append_cut
    ( FFOp*op, const typename PolCut<T>::TYPE type,
      const double b, const PolVar<T>&X1, const double a1,
      const PolVar<T>&X2, const double a2, const PolVar<T>&X3,
      const double a3 );
  //! @brief Appends new relaxation cut in _Cuts w/ <a>n</a> variables
  typename t_Cuts::iterator _append_cut
    ( FFOp*op, const typename PolCut<T>::TYPE type, const double b,
      const unsigned n, const PolVar<T>*X, const double*a );
  //! @brief Appends new relaxation cut in _Cuts w/ <a>n+1</a> variables
  typename t_Cuts::iterator _append_cut
    ( FFOp*op, const typename PolCut<T>::TYPE type, const double b,
      const unsigned n, const PolVar<T>*X, const double*a,
      const PolVar<T>&X1, const double a1 );

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

  //! @brief Append cuts for nonlinear operation using piecewise-linear approximation
  void _semilinear_cuts
    ( FFOp*pOp, const PolVar<T>&X, const double XL, const double XU, const PolVar<T>&Y,
      const typename PolCut<T>::TYPE sense, p_dUniv f, const double*rpar=0, const int*ipar=0 );
  //! @brief Form subintervals for semilinear cuts
  std::vector<double>* _semilinear_sub
    ( const PolVar<T>&X, const double XL, const double XU );

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
      ROOT=1, //!< Error during root search for obtaining the convex/concave envelope of a univariate term
      INTER, //!< Error during intersection of two terms (terms do not intersect)
      DIV, //!< Error during division operation (division by 0)
      ENVMIS=-1, //!< Error due to an operation between variables participating in different polytopic images
      UNAVAIL=-2, //!< Error due to calling a function/feature not yet implemented in MC++
      NOTALLOWED=-3, //!< Error due to calling a function/feature not yet implemented in MC++
      INTERN=-4	//!< Internal error
    };
    //! @brief Constructor for error <a>ierr</a>
    Exceptions( TYPE ierr ) : _ierr( ierr ){}
    //! @brief Inline function returning the error flag
    int ierr(){ return _ierr; }

  private:
    TYPE _ierr;
  };

  //! @brief PolImg options
  struct Options
  {
    //! @brief Constructor
    Options():
      ROOT_USE(true), ROOT_MAXIT(100), ROOT_TOL(1e-10),
      SANDWICH_ATOL(1e-3), SANDWICH_RTOL(1e-3), SANDWICH_MAXCUT(5),
      SANDWICH_RULE(MAXERR), FRACTIONAL_ATOL(machprec()),
      FRACTIONAL_RTOL(machprec()), BREAKPOINT_TYPE(NONE),
      BREAKPOINT_ATOL(1e-8), BREAKPOINT_RTOL(1e-5), DCDECOMP_SCALE(false)
      {}
    //! @brief Assignment operator
    Options& operator= ( Options&options ){
        ROOT_USE        = options.ROOT_USE;
        ROOT_MAXIT      = options.ROOT_MAXIT;
        ROOT_TOL        = options.ROOT_TOL;
        SANDWICH_ATOL   = options.SANDWICH_ATOL;
        SANDWICH_RTOL   = options.SANDWICH_RTOL;
        SANDWICH_MAXCUT = options.SANDWICH_MAXCUT;
        SANDWICH_RULE   = options.SANDWICH_RULE;
        FRACTIONAL_ATOL = options.FRACTIONAL_ATOL;
        FRACTIONAL_RTOL = options.FRACTIONAL_RTOL;
        BREAKPOINT_TYPE     = options.BREAKPOINT_TYPE;
        BREAKPOINT_ATOL     = options.BREAKPOINT_ATOL;
        BREAKPOINT_RTOL     = options.BREAKPOINT_RTOL;
        DCDECOMP_SCALE      = options.DCDECOMP_SCALE;
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
    //! @brief Whether or not to use root search to construct envelopes of univariate terms
    bool ROOT_USE;
    //! @brief Maximal number of iterations in root search
    unsigned ROOT_MAXIT;
    //! @brief Termination tolerance in root search
    double ROOT_TOL;
    //! @brief Absolute tolerance in outer-approximation of univariate terms
    double SANDWICH_ATOL;
    //! @brief Relative tolerance in outer-approximation of univariate terms
    double SANDWICH_RTOL;
    //! @brief Maximal number of cuts in outer approximation of univariate terms
    unsigned SANDWICH_MAXCUT;
    //! @brief Rule for outer-approximation of nonlinear convex/concave terms
    SANDWICH SANDWICH_RULE;
    //! @brief Absolute tolerance in fractional terms to prevent division by zero
    double FRACTIONAL_ATOL;
    //! @brief Relative tolerance in fractional terms to prevent division by zero
    double FRACTIONAL_RTOL;
    //! @brief Rule for piecewise linear cuts of nonlinear convex/concave terms
    REFINE BREAKPOINT_TYPE;
    //! @brief Absolute tolerance in adding breakpoints in piecewise linear cuts
    double BREAKPOINT_ATOL;
    //! @brief Relative tolerance in adding breakpoints in piecewise linear cuts
    double BREAKPOINT_RTOL;
    //! @brief Whether or not to scale variables in DC decomposition of bilinear/fractional terms
    bool DCDECOMP_SCALE;
  } options;

  //! @brief Retreive reference to set of DAG variables in polytopic image
  t_Vars& Vars()
    { return _Vars; }

  //! @brief Retreive const reference to set of DAG variables in polytopic image
  const t_Vars& Vars() const
    { return _Vars; }
  
  //! @brief Retreive reference to set of auxiliary variables in polytopic image
  t_Aux& Aux()
    { return _Aux; }

  //! @brief Retreive reference to set of bilinear terms in polytopic image
  t_Bilin& Bilin()
    { return _Bilin; }

  //! @brief Retreive const reference to set of bilinear terms in polytopic image
  const t_Bilin& Bilin() const
    { return _Bilin; }

  //! @brief Retreive reference to set of cuts in polytopic image
  t_Cuts& Cuts()
    { return _Cuts; }

  //! @brief Reset polytopic image (all participating variables and cuts)
  void reset()
    { _erase_cuts(); _erase_vars(); _erase_aux(); _erase_bilin(); }

  //! @brief Append new relaxation cut w/ 1 variable
  typename t_Cuts::iterator add_cut
    ( const typename PolCut<T>::TYPE type, const double b,
      const PolVar<T>&X1, const double a1 )
    { return _append_cut( 0, type, b, X1, a1 ); }
  //! @brief Append new relaxation cut w/ 2 variables
  typename t_Cuts::iterator add_cut
    ( const typename PolCut<T>::TYPE type, const double b,
      const PolVar<T>&X1, const double a1,
      const PolVar<T>&X2, const double a2 )
    { return _append_cut( 0, type, b, X1, a1, X2, a2 ); }
  //! @brief Append new relaxation cut w/ 3 variables
  typename t_Cuts::iterator add_cut
    ( const typename PolCut<T>::TYPE type, const double b,
      const PolVar<T>&X1, const double a1,
      const PolVar<T>&X2, const double a2,
      const PolVar<T>&X3, const double a3 )
    { return _append_cut( 0, type, b, X1, a1, X2, a2, X3, a3 ); }
  //! @brief Append new relaxation cut w/ <a>n</a> variables
  typename t_Cuts::iterator add_cut
    ( const typename PolCut<T>::TYPE type, const double b,
      const unsigned n, const PolVar<T>*X, const double*a )
    { return _append_cut( 0, type, b, n, X, a ); }
  //! @brief Append new relaxation cut w/ <a>n</a> variables
  typename t_Cuts::iterator add_cut
    ( const typename PolCut<T>::TYPE type, const double b,
      const unsigned n, const PolVar<T>*X, const double*a,
      const PolVar<T>&X1, const double a1 )
    { return _append_cut( 0, type, b, n, X, a, X1, a1 ); }

  //! @brief Erase cut with iterator <a>itcut</a> from set of cuts
  void erase_cut
    ( typename t_Cuts::iterator itcut )
    { return _Cuts.erase( itcut ); }
  //! @brief Erase all cuts and auxiliary variables
  void reset_cuts
    ()
    { _reset_vars(); _erase_aux(); _erase_cuts(); }

/*
  //! @brief Add constraint to factorable program
  PolCut<T>* add_constraint
    ( const PolVar<T>&lhsVar, const typename PolImg<T>::CTRTYPE type,
      const PolVar<T>&rhsVar );
  //! @brief Remove constraint from factorable program
  virtual bool remove_constraint
    ( FPOp<T>* constr )
    { return _erase_operation( constr ); }
*/
};

template <typename T>
inline PolVar<T>*
PolImg<T>::_append_var
( FFVar*var, const T&range, const bool cont )
{
  auto itVar = _Vars.find( var );
  if( itVar != _Vars.end() ){
    itVar->second->_update( range, cont );
    return itVar->second;
  }
  PolVar<T>* pVar = new PolVar<T>;
  pVar->_set( this, *var, range, cont, _Vars.size() );
  _Vars.insert( std::make_pair( var, pVar ) );
  return pVar;
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
inline PolBilin<T>*
PolImg<T>::_append_bilin
( const PolVar<T>*varR )
{
  // Determine if bilinear term already exists
  auto itb = _Bilin.find( &varR->var() );
  if( itb != _Bilin.end() ){
    //itVar->second->_update( range, cont );
    return itb->second;
  }

  // Append bilinear or fractional term
  //const PolVar<T>*pvarR = _Vars.find( const_cast<FFVar*>(&varR->var()) )->second;
  //PolBilin<T>* pBilin = new PolBilin<T>( pvarR );
  PolBilin<T>* pBilin = new PolBilin<T>( varR );
  _Bilin.insert( std::make_pair( &varR->var(), pBilin ) );
  return pBilin;
}

template <typename T>
inline void
PolImg<T>::_decompose_bilin
( FFOp*op, PolBilin<T>* pBilin, const PolVar<T>&varR, const PolVar<T>&var1,
  const PolVar<T>&var2 )
{
  if( options.BREAKPOINT_TYPE == options.NONE ) return;

  // Add cuts for DC decomposition
  FFGraph* dag = varR._var.dag();

  if( options.DCDECOMP_SCALE ){
    double s1 = Op<T>::diam(var2.range());
    //double s1 = options.DCDECOMP_SCALE? 1./Op<T>::diam(var1.range()): 1.;
    if( pBilin->scal1 ){
      pBilin->scal1->_var.num() = s1;
      pBilin->scal1->_range = s1;
    }
    else{ 
      FFVar scal1( dag, s1 );
      FFVar* pscal1 = *dag->Vars().find(&scal1);
      PolVar<T> polscal1( this, *pscal1 );
      pBilin->scal1 = _Vars.find( &scal1 )->second;
    }

    double s2 = Op<T>::diam(var1.range());
    //double s2 = options.DCDECOMP_SCALE? 1./Op<T>::diam(var2.range()): 1.;
    if( pBilin->scal2 ){
      pBilin->scal2->_var.num() = s2;
      pBilin->scal2->_range = s2;
    }
    else{ 
      FFVar scal2( dag, s2 );
      FFVar* pscal2 = *dag->Vars().find(&scal2);
      PolVar<T> polscal2( this, *pscal2 );
      pBilin->scal2 = _Vars.find( &scal2 )->second;
    }

    FFVar scalvar1( pBilin->scal1->_var * var1._var );
    dag->curOp() = scalvar1.ops().first;
    PolVar<T> polscalvar1( *pBilin->scal1 * var1 );

    FFVar scalvar2( pBilin->scal2->_var * var2._var );
    dag->curOp() = scalvar2.ops().first;
    PolVar<T> polscalvar2( *pBilin->scal2 * var2 );

    FFVar sumvar( scalvar1 + scalvar2 );
    dag->curOp() = sumvar.ops().first;
    PolVar<T> polsumvar( polscalvar1 + polscalvar2 );

    FFVar subvar( scalvar1 - scalvar2 );
    dag->curOp() = subvar.ops().first;
    PolVar<T> polsubvar( polscalvar1 - polscalvar2 );

    FFVar sqrvar1( sqr( sumvar ) );
    dag->curOp() = sqrvar1.ops().first;
    PolVar<T> polsqrvar1( sqr( polsumvar ) );

    FFVar sqrvar2( sqr( subvar ) );
    dag->curOp() = sqrvar2.ops().first;
    PolVar<T> polsqrvar2( sqr( polsubvar ) );

    _append_cut( op, PolCut<T>::EQ, 0., varR, 4.*s1*s2, polsqrvar1, -1., polsqrvar2, 1. );
  }

  else{
    FFVar sumvar( var1._var + var2._var );
    dag->curOp() = sumvar.ops().first;
    PolVar<T> polsumvar( var1 + var2 );

    FFVar subvar( var1._var - var2._var );
    dag->curOp() = subvar.ops().first;
    PolVar<T> polsubvar( var1 - var2 );

    FFVar sqrvar1( sqr( sumvar ) );
    dag->curOp() = sqrvar1.ops().first;
    PolVar<T> polsqrvar1( sqr( polsumvar ) );

    FFVar sqrvar2( sqr( subvar ) );
    dag->curOp() = sqrvar2.ops().first;
    PolVar<T> polsqrvar2( sqr( polsubvar ) );

    _append_cut( op, PolCut<T>::EQ, 0., varR, 4., polsqrvar1, -1., polsqrvar2, 1. );
  }
  return;
}

template <typename T> inline void
PolImg<T>::_erase_bilin
()
{
  for( auto itv = _Bilin.begin(); itv != _Bilin.end(); ++itv )
    delete itv->second;
  _Bilin.clear();
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
PolImg<T>::_append_cut
( FFOp*op, const typename PolCut<T>::TYPE type,
  const double b, const PolVar<T>&X1, const double a1 )
{
  PolCut<T>* pCut = new PolCut<T>( op, type, b, X1, a1 );
  return _Cuts.insert( pCut );
}

template <typename T>
inline typename PolImg<T>::t_Cuts::iterator
PolImg<T>::_append_cut
( FFOp*op, const typename PolCut<T>::TYPE type,
  const double b, const PolVar<T>&X1, const double a1,
  const PolVar<T>&X2, const double a2 )
{
  PolCut<T>* pCut = new PolCut<T>( op, type, b, X1, a1, X2, a2 );
  return _Cuts.insert( pCut );
}

template <typename T>
inline typename PolImg<T>::t_Cuts::iterator
PolImg<T>::_append_cut
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
PolImg<T>::_append_cut
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
PolImg<T>::_append_cut
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

  out << ( P._Bilin.empty()? "\nNO BILINEAR OR FRACTIONAL TERM": "\nBILINEAR AND FRACTIONAL TERMS:\n" );
  for( auto itv=P._Bilin.begin(); itv!=P._Bilin.end(); ++itv ){
    out << "  " << itv->second->var->name();
  }
  out << std::endl;

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
  if( Var1._var.cst() )
    return Var2 ^ Var1.var().num().val();
  else if( Var2._var.cst() )
    return Var1 ^ Var2.var().num().val();

  if( Var1._img != Var2._img )
    throw typename PolImg<T>::Exceptions( PolImg<T>::Exceptions::ENVMIS );
  FFGraph* dag = Var1._var.dag();
  if( !dag || !dag->curOp() )
    throw typename PolImg<T>::Exceptions( PolImg<T>::Exceptions::NOTALLOWED );
  FFVar* pFFVarR = dag->curOp()->pres;
  T IVarR;
  if( !Op<T>::inter( IVarR, Var1._range, Var2._range ) )
    throw typename PolImg<T>::Exceptions( PolImg<T>::Exceptions::INTER );
  PolVar<T>* pVarR = Var1._img->_append_var( pFFVarR, IVarR, true );
  if( pVarR->cuts() ) return *pVarR;
  pVarR->set_cuts();

  Var1._img->_append_cut( pFFVarR->ops().first, PolCut<T>::EQ, 0., *pVarR, 1., Var1, -1. );
  Var1._img->_append_cut( pFFVarR->ops().first, PolCut<T>::EQ, 0., *pVarR, 1., Var2, -1. );
  return *pVarR;
}

template <typename T>
inline PolVar<T>
operator^
( const PolVar<T>&Var1, const double Cst2 )
{
  if( Var1._var.cst() )
    throw typename PolImg<T>::Exceptions( PolImg<T>::Exceptions::INTER );

  FFGraph* dag = Var1._var.dag();
  if( !dag || !dag->curOp() )
    throw typename PolImg<T>::Exceptions( PolImg<T>::Exceptions::NOTALLOWED );
  FFVar* pFFVarR = dag->curOp()->pres;
  T IVarR;
  if( !Op<T>::inter( IVarR, Var1._range, Cst2 ) )
    throw typename PolImg<T>::Exceptions( PolImg<T>::Exceptions::INTER );
  PolVar<T>* pVarR = Var1._img->_append_var( pFFVarR, IVarR, true );
  if( pVarR->cuts() ) return *pVarR;
  pVarR->set_cuts();

  Var1._img->_append_cut( pFFVarR->ops().first, PolCut<T>::EQ, 0., *pVarR, 1., Var1, -1. );
  Var1._img->_append_cut( pFFVarR->ops().first, PolCut<T>::EQ, Cst2, *pVarR, 1. );
  return *pVarR;
}

template <typename T>
inline PolVar<T>
operator^
( const double Cst1, const PolVar<T>&Var2 )
{
  return Var2 ^ Cst1;
}

template <typename T>
inline PolVar<T>
operator+
( const PolVar<T>&Var )
{
  return Var; 
}

template <typename T>
inline PolVar<T>
operator+
( const PolVar<T>&Var1, const PolVar<T>&Var2 )
{
  if( Var1._var.cst() )
    return Var2 + Var1.var().num().val();
  else if( Var2._var.cst() )
    return Var1 + Var2.var().num().val();

  if( Var1._img != Var2._img )
    throw typename PolImg<T>::Exceptions( PolImg<T>::Exceptions::ENVMIS );
  FFGraph* dag = Var1._var.dag();
  if( !dag || !dag->curOp() )
    throw typename PolImg<T>::Exceptions( PolImg<T>::Exceptions::NOTALLOWED );
  FFVar* pFFVarR = dag->curOp()->pres;
  PolVar<T>* pVarR = Var1._img->_append_var( pFFVarR, Var1._range + Var2._range, true );
  if( pVarR->cuts() ) return *pVarR;
  pVarR->set_cuts();

  Var1._img->_append_cut( pFFVarR->ops().first, PolCut<T>::EQ, 0., *pVarR, 1., Var1, -1., Var2, -1. );
  return *pVarR;
}

template <typename T>
inline PolVar<T>
operator+
( const PolVar<T>&Var1, const double Cst2 )
{
  if( Var1._var.cst() ) return Var1._var.num().val() + Cst2;

  FFGraph* dag = Var1._var.dag();
  if( !dag || !dag->curOp() )
    throw typename PolImg<T>::Exceptions( PolImg<T>::Exceptions::NOTALLOWED );
  FFVar* pFFVarR = dag->curOp()->pres;
  PolVar<T>* pVarR = Var1._img->_append_var( pFFVarR, Var1._range + Cst2, true );
  if( pVarR->cuts() ) return *pVarR;
  pVarR->set_cuts();

  Var1._img->_append_cut( pFFVarR->ops().first, PolCut<T>::EQ, Cst2, *pVarR, 1., Var1, -1. );
  return *pVarR;
}

template <typename T>
inline PolVar<T>
operator+
( const double Cst1, const PolVar<T>&Var2 )
{
  return Var2 + Cst1;
}

template <typename T>
inline PolVar<T>
operator-
( const PolVar<T>&Var1 )
{
  if( Var1._var.cst() ) return -Var1._var.num().val();

  FFGraph* dag = Var1._var.dag();
  if( !dag || !dag->curOp() )
    throw typename PolImg<T>::Exceptions( PolImg<T>::Exceptions::NOTALLOWED );
  FFVar* pFFVarR = dag->curOp()->pres;
  PolVar<T>* pVarR = Var1._img->_append_var( pFFVarR, -Var1._range, true );
  if( pVarR->cuts() ) return *pVarR;
  pVarR->set_cuts();

  Var1._img->_append_cut( pFFVarR->ops().first, PolCut<T>::EQ, 0, *pVarR, 1., Var1, 1. );
  return *pVarR;
}

template <typename T>
inline PolVar<T>
operator-
( const PolVar<T>&Var1, const PolVar<T>&Var2 )
{
  if( Var1._var.cst() )
    return Var1.var().num().val() - Var2;
  else if( Var2._var.cst() )
    return Var1 - Var2.var().num().val();

  if( Var1._img != Var2._img )
    throw typename PolImg<T>::Exceptions( PolImg<T>::Exceptions::ENVMIS );
  FFGraph* dag = Var1._var.dag();
  if( !dag || !dag->curOp() )
    throw typename PolImg<T>::Exceptions( PolImg<T>::Exceptions::NOTALLOWED );
  FFVar* pFFVarR = dag->curOp()->pres;
  PolVar<T>* pVarR = Var1._img->_append_var( pFFVarR, Var1._range - Var2._range, true );
  if( pVarR->cuts() ) return *pVarR;
  pVarR->set_cuts();

  Var1._img->_append_cut( pFFVarR->ops().first, PolCut<T>::EQ, 0., *pVarR, 1., Var1, -1., Var2, 1. );
  return *pVarR;
}

template <typename T>
inline PolVar<T>
operator-
( const PolVar<T>&Var1, const double Cst2 )
{
  if( Var1._var.cst() ) return Var1._var.num().val() - Cst2;

  FFGraph* dag = Var1._var.dag();
  if( !dag || !dag->curOp() )
    throw typename PolImg<T>::Exceptions( PolImg<T>::Exceptions::NOTALLOWED );
  FFVar* pFFVarR = dag->curOp()->pres;
  PolVar<T>* pVarR = Var1._img->_append_var( pFFVarR, Var1._range - Cst2, true );
  if( pVarR->cuts() ) return *pVarR;
  pVarR->set_cuts();

  Var1._img->_append_cut( pFFVarR->ops().first, PolCut<T>::EQ, -Cst2, *pVarR, 1., Var1, -1. );
  return *pVarR;
}

template <typename T>
inline PolVar<T>
operator-
( const double Cst1, const PolVar<T>&Var2 )
{
  if( Var2._var.cst() ) return Cst1 - Var2._var.num().val();

  FFGraph* dag = Var2._var.dag();
  if( !dag || !dag->curOp() )
    throw typename PolImg<T>::Exceptions( PolImg<T>::Exceptions::NOTALLOWED );
  FFVar* pFFVarR = dag->curOp()->pres;
  PolVar<T>* pVarR = Var2._img->_append_var( pFFVarR, Cst1 - Var2._range, true );
  if( pVarR->cuts() ) return *pVarR;
  pVarR->set_cuts();

  Var2._img->_append_cut( pFFVarR->ops().first, PolCut<T>::EQ, Cst1, *pVarR, 1., Var2, 1. );
  return *pVarR;
}

template <typename T>
inline PolVar<T>
operator*
( const PolVar<T>&Var1, const PolVar<T>&Var2 )
{
  if( Var1._var.cst() )
    return Var2 * Var1.var().num().val();
  else if( Var2._var.cst() )
    return Var1 * Var2.var().num().val();

  if( Var1._img != Var2._img )
    throw typename PolImg<T>::Exceptions( PolImg<T>::Exceptions::ENVMIS );
  FFGraph* dag = Var1._var.dag();
  if( !dag || !dag->curOp() )
    throw typename PolImg<T>::Exceptions( PolImg<T>::Exceptions::NOTALLOWED );
  FFVar* pFFVarR = dag->curOp()->pres;
  PolVar<T>* pVarR = Var1._img->_append_var( pFFVarR, Var1._range * Var2._range, true );
  if( pVarR->cuts() ) return *pVarR;
  pVarR->set_cuts();

  Var1._img->_append_cut( pFFVarR->ops().first, PolCut<T>::GE,
      -Op<T>::u(Var1._range)*Op<T>::u(Var2._range),
      *pVarR, 1., Var1, -Op<T>::u(Var2._range), Var2, -Op<T>::u(Var1._range) );
  Var1._img->_append_cut( pFFVarR->ops().first, PolCut<T>::GE,
      -Op<T>::l(Var1._range)*Op<T>::l(Var2._range),
      *pVarR, 1., Var1, -Op<T>::l(Var2._range), Var2, -Op<T>::l(Var1._range) );
  Var1._img->_append_cut( pFFVarR->ops().first, PolCut<T>::LE,
      -Op<T>::u(Var1._range)*Op<T>::l(Var2._range),
      *pVarR, 1., Var1, -Op<T>::l(Var2._range), Var2, -Op<T>::u(Var1._range) );
  Var1._img->_append_cut( pFFVarR->ops().first, PolCut<T>::LE,
      -Op<T>::l(Var1._range)*Op<T>::u(Var2._range),
      *pVarR, 1., Var1, -Op<T>::u(Var2._range), Var2, -Op<T>::l(Var1._range) );

  PolBilin<T>* pBilin = Var1._img->_append_bilin( pVarR );
  Var1._img->_decompose_bilin( pFFVarR->ops().first, pBilin, *pVarR, Var1, Var2 );
  //std::cout << *Var1._img; std::cout << "PAUSE "; int dum; std::cin >> dum;
  return *pVarR;
}

template <typename T>
inline PolVar<T>
operator*
( const PolVar<T>&Var1, const double Cst2 )
{
  if( Var1._var.cst() ) return Var1._var.num().val() * Cst2;

  FFGraph* dag = Var1._var.dag();
  if( !dag || !dag->curOp() )
    throw typename PolImg<T>::Exceptions( PolImg<T>::Exceptions::NOTALLOWED );
  FFVar* pFFVarR = dag->curOp()->pres;
  PolVar<T>* pVarR = Var1._img->_append_var( pFFVarR, Var1._range * Cst2, true );
  if( pVarR->cuts() ) return *pVarR;
  pVarR->set_cuts();

  Var1._img->_append_cut( pFFVarR->ops().first, PolCut<T>::EQ, 0., *pVarR, 1., Var1, -Cst2 );
  //std::cout << "operator*( const PolVar<T>&Var1, const double Cst2 ):\n";
  //std::cout << *Var1._img;
  return *pVarR;
}

template <typename T>
inline PolVar<T>
operator*
( const double Cst1, const PolVar<T>&Var2 )
{
  return Var2 * Cst1;
}

template <typename T>
inline PolVar<T>
operator/
( const PolVar<T>&Var1, const PolVar<T>&Var2 )
{
  if( Var1._var.cst() )
    return Var1.var().num().val() / Var2;
  else if( Var2._var.cst() )
    return Var1 / Var2.var().num().val();

  if( Var1._img != Var2._img )
    throw typename PolImg<T>::Exceptions( PolImg<T>::Exceptions::ENVMIS );
  FFGraph* dag = Var1._var.dag();
  if( !dag || !dag->curOp() )
    throw typename PolImg<T>::Exceptions( PolImg<T>::Exceptions::NOTALLOWED );

  FFVar* pFFVarR = dag->curOp()->pres;
  PolVar<T>* pVarR = Var1._img->_append_var( pFFVarR, Var1._range / Var2._range, true );
  if( pVarR->cuts() ) return *pVarR;
  pVarR->set_cuts();

  Var1._img->_append_cut( pFFVarR->ops().first, PolCut<T>::GE,
      -Op<T>::u(pVarR->_range)*Op<T>::u(Var2._range),
      Var1, 1., *pVarR, -Op<T>::u(Var2._range), Var2, -Op<T>::u(pVarR->_range) );
  Var1._img->_append_cut( pFFVarR->ops().first, PolCut<T>::GE,
      -Op<T>::l(pVarR->_range)*Op<T>::l(Var2._range),
      Var1, 1., *pVarR, -Op<T>::l(Var2._range), Var2, -Op<T>::l(pVarR->_range) );
  Var1._img->_append_cut( pFFVarR->ops().first, PolCut<T>::LE,
      -Op<T>::u(pVarR->_range)*Op<T>::l(Var2._range),
      Var1, 1., *pVarR, -Op<T>::l(Var2._range), Var2, -Op<T>::u(pVarR->_range) );
  Var1._img->_append_cut( pFFVarR->ops().first, PolCut<T>::LE,
      -Op<T>::l(pVarR->_range)*Op<T>::u(Var2._range),
      Var1, 1., *pVarR, -Op<T>::u(Var2._range), Var2, -Op<T>::l(pVarR->_range) );

  PolBilin<T>* pBilin = Var1._img->_append_bilin( pVarR );
  Var1._img->_decompose_bilin( pFFVarR->ops().first, pBilin, Var1, *pVarR, Var2 );
  return *pVarR;
}

template <typename T>
inline PolVar<T>
operator/
( const double Cst1, const PolVar<T>&Var2 )
{
  //return inv( Var2 ) * Cst1;
  if( Var2._var.cst() ) return Cst1/Var2.var().num().val();

  FFGraph* dag = Var2._var.dag();
  if( !dag || !dag->curOp() )
    throw typename PolImg<T>::Exceptions( PolImg<T>::Exceptions::NOTALLOWED );

  FFVar* pFFVarR = dag->curOp()->pres;
  FFOp* pOp = pFFVarR->ops().first;
  PolVar<T>* pVarR = Var2._img->_append_var( pFFVarR, Cst1/Var2._range, true );
  if( pVarR->cuts() ) return *pVarR;
  pVarR->set_cuts();

  // Constant term
  //if( Var2._var.cst() ){
  //  Var2._img->_append_cut( pOp, PolCut<T>::EQ, Cst1/Var2._var.num().val(), *pVarR, 1. );
  //  return *pVarR;
  //}

  struct loc{
    static std::pair<double,double> scalinv
      ( const double x, const double*rusr, const int*iusr )
      { return std::make_pair( *rusr/x, -*rusr/(x*x) ); }
  };
  // -- Convex Case
  if( (Cst1 >= 0. && Op<T>::l(Var2._range) > 0.) || (Cst1 <= 0. && Op<T>::u(Var2._range) < 0.) ){
    Var2._img->_semilinear_cuts( pOp, Var2, Op<T>::l(Var2._range), Op<T>::u(Var2._range),
      *pVarR, PolCut<T>::LE, loc::scalinv, &Cst1, 0 );
    Var2._img->_sandwich_cuts( pOp, Var2, Op<T>::l(Var2._range), Op<T>::u(Var2._range),
      *pVarR, Op<T>::l(pVarR->_range), Op<T>::u(pVarR->_range), PolCut<T>::GE, loc::scalinv, &Cst1, 0 );
  }
  // -- Concave Case
  else{
    Var2._img->_semilinear_cuts( pOp, Var2, Op<T>::l(Var2._range), Op<T>::u(Var2._range),
      *pVarR, PolCut<T>::GE, loc::scalinv, &Cst1, 0 );
    Var2._img->_sandwich_cuts( pOp, Var2, Op<T>::l(Var2._range), Op<T>::u(Var2._range),
      *pVarR, Op<T>::l(pVarR->_range), Op<T>::u(pVarR->_range), PolCut<T>::LE, loc::scalinv, &Cst1, 0 );
  }
  return *pVarR;
}

template <typename T>
inline PolVar<T>
operator/
( const PolVar<T>&Var1, const double Cst2 )
{
  if( Cst2 == 0. ) 
    throw typename PolImg<T>::Exceptions( PolImg<T>::Exceptions::DIV );
  return Var1 * ( 1. / Cst2 );
}

template <typename T>
inline PolVar<T>
inv
( const PolVar<T>&Var1 )
{
  return pow( Var1, PolVar<T>(-1) );
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
      const std::vector<double>& subint = X.create_subdiv( XL, XU );
      const unsigned nsubint = subint.size();
      if( nsubint > 2 ){
        for( ; NBCUTS<nsubint; NBCUTS++ ){
          _linearization_cut( pOp, subint[NBCUTS], X, XL, XU, Y, YL, YU, sense, f, rpar, ipar );
          if( !NBCUTS ) continue;
          std::pair<double,double> worst = _distmax( f, subint[NBCUTS-1], subint[NBCUTS], rpar, ipar );
          OA.push( OAsub( subint[NBCUTS-1], subint[NBCUTS], worst.first, worst.second ) );
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
  while( OA.top().gap() > dtol && NBCUTS < options.SANDWICH_MAXCUT ){
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
  _append_cut( pOp, sense, Fref.first-Xref*Fref.second, Y, 1., X, -Fref.second );
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

template <typename T>
inline void
PolImg<T>::_semilinear_cuts
( FFOp*pOp, const PolVar<T>&X, const double XL, const double XU, const PolVar<T>&Y,
  const typename PolCut<T>::TYPE sense, p_dUniv f, const double*rpar, const int*ipar )
{
  switch( options.BREAKPOINT_TYPE ){
   case options.BIN:{
    const std::vector<double>& subint = X.create_subdiv( XL, XU );
    const unsigned nsubint = subint.size()-1;
    if( nsubint > 1 ){
      // Represent variable range using linear binary transformation
      const std::vector< PolVar<T> >& subvar = X.BIN_subdiv( pOp );
      // Append semilinear cuts
      double coef[nsubint];
      double rhs = f( subint[0], rpar, ipar ).first;
      for( unsigned isub=0; isub<nsubint; isub++ )
        coef[isub] = f( subint[isub], rpar, ipar ).first - f( subint[isub+1], rpar, ipar ).first;
      _append_cut( pOp, sense, rhs, nsubint, subvar.data(), coef, Y, 1. );
    }
    break;
   }
   case options.SOS2:{
    const std::vector<double>& subint = X.create_subdiv( XL, XU );
    const unsigned nsubint = subint.size();
    if( nsubint > 2 ){
      // Represent variable range using SOS2 transformation
      const std::vector< PolVar<T> >& subvar = X.SOS2_subdiv( pOp );
      // Append semilinear cuts
      double coef[nsubint];
      for( unsigned isub=0; isub<nsubint; isub++ )
        coef[isub] = -f( subint[isub], rpar, ipar ).first;
      _append_cut( pOp, sense, 0., nsubint, subvar.data(), coef, Y, 1. );
    }
    break;
   }
   default:
    break;
  }
  double dX = XU-XL, YL = f(XL,rpar,ipar).first, dY = f(XU,rpar,ipar).first-YL;
  _append_cut( pOp, sense, dX*YL-dY*XL, Y, dX, X, -dY );
}

template <typename T>
inline PolVar<T>
exp
( const PolVar<T>&Var1 )
{
  if( Var1._var.cst() ) return std::exp( Var1.var().num().val() );

  FFGraph* dag = Var1._var.dag();
#ifdef MC__POLIMG_CHECK
  if( !dag || !dag->curOp() )
    throw typename PolImg<T>::Exceptions( PolImg<T>::Exceptions::NOTALLOWED );
#endif

  FFVar* pFFVarR = dag->curOp()->pres;
  FFOp* pOp = pFFVarR->ops().first;
  PolVar<T>* pVarR = Var1._img->_append_var( pFFVarR, Op<T>::exp( Var1._range ), true );
  if( pVarR->cuts() ) return *pVarR;
  pVarR->set_cuts();

  // Constant term
  //if( Var1._var.cst() ){
  //  Var1._img->_append_cut( pOp, PolCut<T>::EQ, std::exp(Var1._var.num().val()), *pVarR, 1. );
  //  return *pVarR;
  //}

  // Variable term
  struct loc{ static std::pair<double,double> exp
    ( const double x, const double*rusr, const int*iusr )
    { return std::make_pair( std::exp(x), std::exp(x) ); }
  };
  Var1._img->_semilinear_cuts( pOp, Var1, Op<T>::l(Var1._range),
    Op<T>::u(Var1._range), *pVarR, PolCut<T>::LE, loc::exp );
  Var1._img->_sandwich_cuts( pOp, Var1, Op<T>::l(Var1._range),
    Op<T>::u(Var1._range), *pVarR, Op<T>::l(pVarR->_range), Op<T>::u(pVarR->_range),
    PolCut<T>::GE, loc::exp );

  return *pVarR;
}

template <typename T>
inline PolVar<T>
log
( const PolVar<T>&Var1 )
{
  if( Var1._var.cst() ) return std::log( Var1.var().num().val() );

  FFGraph* dag = Var1._var.dag();
#ifdef MC__POLIMG_CHECK
  if( !dag || !dag->curOp() )
    throw typename PolImg<T>::Exceptions( PolImg<T>::Exceptions::NOTALLOWED );
    //return std::log( Var1.var().num().val() );
#endif

  FFVar* pFFVarR = dag->curOp()->pres;
  FFOp* pOp = pFFVarR->ops().first;
  PolVar<T>* pVarR = Var1._img->_append_var( pFFVarR, Op<T>::log( Var1._range ), true );
  if( pVarR->cuts() ) return *pVarR;
  pVarR->set_cuts();

  // Constant term
  //if( Var1._var.cst() ){
  //  Var1._img->_append_cut( pOp, PolCut<T>::EQ, std::log(Var1._var.num().val()), *pVarR, 1. );
  //  return *pVarR;
  //}

  // Variable term
  struct loc{
    static std::pair<double,double> exp
      ( const double x, const double*rusr, const int*iusr )
      { return std::make_pair( std::exp(x), std::exp(x) ); }
    static std::pair<double,double> log
      ( const double x, const double*rusr, const int*iusr )
      { return std::make_pair( std::log(x), 1/x ); }
  };
  Var1._img->_semilinear_cuts( pOp, Var1, Op<T>::l(Var1._range),
    Op<T>::u(Var1._range), *pVarR, PolCut<T>::GE, loc::log );
  Var1._img->_sandwich_cuts( pOp, *pVarR, Op<T>::l(pVarR->_range),
    Op<T>::u(pVarR->_range), Var1, Op<T>::l(Var1._range), Op<T>::u(Var1._range),
    PolCut<T>::GE, loc::exp );

  return *pVarR;
}

template <typename T>
inline PolVar<T>
sqr
( const PolVar<T>&Var1 )
{
  if( Var1._var.cst() ) return mc::sqr( Var1.var().num().val() );

  FFGraph* dag = Var1._var.dag();
#ifdef MC__POLIMG_CHECK
  if( !dag || !dag->curOp() )
    throw typename PolImg<T>::Exceptions( PolImg<T>::Exceptions::NOTALLOWED );
    //return mc::sqr( Var1.var().num().val() );
#endif

  FFVar* pFFVarR = dag->curOp()->pres;
  FFOp* pOp = pFFVarR->ops().first;
  PolVar<T>* pVarR = Var1._img->_append_var( pFFVarR, Op<T>::sqr( Var1._range ), true );
  if( pVarR->cuts() ) return *pVarR;
  pVarR->set_cuts();

  // Constant term
  //if( Var1._var.cst() ){
  //  Var1._img->_append_cut( pOp, PolCut<T>::EQ, mc::sqr(Var1._var.num().val()), *pVarR, 1. );
  //  return *pVarR;
  //}

  // Variable term
  struct loc{ static std::pair<double,double> sqr
    ( const double x, const double*rusr, const int*iusr )
    { return std::make_pair( x*x, 2.*x ); }
  };
  Var1._img->_semilinear_cuts( pOp, Var1, Op<T>::l(Var1._range),
    Op<T>::u(Var1._range), *pVarR, PolCut<T>::LE, loc::sqr );
  Var1._img->_sandwich_cuts( pOp, Var1, Op<T>::l(Var1._range),
    Op<T>::u(Var1._range), *pVarR, Op<T>::l(pVarR->_range), Op<T>::u(pVarR->_range),
    PolCut<T>::GE, loc::sqr );

  return *pVarR;
}

template <typename T>
inline PolVar<T>
sqrt
( const PolVar<T>&Var1 )
{
  if( Var1._var.cst() ) return std::sqrt( Var1.var().num().val() );

  FFGraph* dag = Var1._var.dag();
#ifdef MC__POLIMG_CHECK
  if( !dag || !dag->curOp() )
    throw typename PolImg<T>::Exceptions( PolImg<T>::Exceptions::NOTALLOWED );
    //return mc::sqr( Var1.var().num().val() );
#endif

  FFVar* pFFVarR = dag->curOp()->pres;
  FFOp* pOp = pFFVarR->ops().first;
  PolVar<T>* pVarR = Var1._img->_append_var( pFFVarR, Op<T>::sqrt( Var1._range ), true );
  if( pVarR->cuts() ) return *pVarR;
  pVarR->set_cuts();

  // Constant term
  //if( Var1._var.cst() ){
  //  Var1._img->_append_cut( pOp, PolCut<T>::EQ, std::sqrt(Var1._var.num().val()), *pVarR, 1. );
  //  return *pVarR;
  //}

  // Variable term
  struct loc{ 
    static std::pair<double,double> sqr
      ( const double x, const double*rusr, const int*iusr )
      { return std::make_pair( x*x, 2.*x ); }
    static std::pair<double,double> sqrt
      ( const double x, const double*rusr, const int*iusr )
      { return std::make_pair( std::sqrt(x), 1/(2*std::sqrt(x)) ); }
  };
  Var1._img->_semilinear_cuts( pOp, Var1, Op<T>::l(Var1._range),
    Op<T>::u(Var1._range), *pVarR, PolCut<T>::GE, loc::sqrt );
  Var1._img->_sandwich_cuts( pOp, *pVarR, Op<T>::l(pVarR->_range),
    Op<T>::u(pVarR->_range), Var1, Op<T>::l(Var1._range), Op<T>::u(Var1._range),
    PolCut<T>::GE, loc::sqr );

  return *pVarR;
}

template <typename T>
inline PolVar<T>
pow
( const PolVar<T>&Var1, const PolVar<T>&Var2 )
{
#ifdef MC__POLIMG_CHECK
  if( !Var2_.var.cst() && Var2_.var.num().t != FFNum::INT )
    throw typename PolImg<T>::Exceptions( PolImg<T>::Exceptions::NOTALLOWED );
#endif
  const int iExp = Var2._var.num().n;
#ifdef MC__POLIMG_CHECK
  if( iExp == 0 || iExp == 1 || iExp == 2 )
    throw typename PolImg<T>::Exceptions( PolImg<T>::Exceptions::NOTALLOWED );
#endif
  if( Var1._var.cst() ) return std::pow( Var1.var().num().val(), iExp );

  FFGraph* dag = Var1._var.dag();
#ifdef MC__POLIMG_CHECK
  if( !dag || !dag->curOp() )
    throw typename PolImg<T>::Exceptions( PolImg<T>::Exceptions::NOTALLOWED );
#endif

  FFVar* pFFVarR = dag->curOp()->pres;
  FFOp* pOp = pFFVarR->ops().first;
  PolVar<T>* pVarR = Var1._img->_append_var( pFFVarR, Op<T>::pow( Var1._range, iExp ), true );
  if( pVarR->cuts() ) return *pVarR;
  pVarR->set_cuts();

  // Constant term
  //if( Var1._var.cst() ){
  //  Var1._img->_append_cut( pOp, PolCut<T>::EQ, std::pow(Var1._var.num().val(),iExp), *pVarR, 1. );
  //  return *pVarR;
  //}

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
    Var1._img->_semilinear_cuts( pOp, Var1, Op<T>::l(Var1._range), Op<T>::u(Var1._range),
      *pVarR, PolCut<T>::LE, loc::pow, 0, &iExp );
    Var1._img->_sandwich_cuts( pOp, Var1, Op<T>::l(Var1._range), Op<T>::u(Var1._range),
      *pVarR, Op<T>::l(pVarR->_range), Op<T>::u(pVarR->_range), PolCut<T>::GE, loc::pow, 0, &iExp );
  }

  // Positive odd exponent term
  else if( iExp > 0 && Var1._img->options.ROOT_USE ){
    // -- Convex Portion
    if( Op<T>::l(Var1._range) >= 0. ){
      Var1._img->_semilinear_cuts( pOp, Var1, Op<T>::l(Var1._range), Op<T>::u(Var1._range),
        *pVarR, PolCut<T>::LE, loc::pow, 0, &iExp );
      Var1._img->_sandwich_cuts( pOp, Var1, Op<T>::l(Var1._range), Op<T>::u(Var1._range),
        *pVarR, Op<T>::l(pVarR->_range), Op<T>::u(pVarR->_range), PolCut<T>::GE, loc::pow, 0, &iExp );
    }
    // -- Concave Portion
    else if( Op<T>::u(Var1._range) <= 0. ){
      Var1._img->_semilinear_cuts( pOp, Var1, Op<T>::l(Var1._range), Op<T>::u(Var1._range),
        *pVarR, PolCut<T>::GE, loc::pow, 0, &iExp );
      Var1._img->_sandwich_cuts( pOp, Var1, Op<T>::l(Var1._range), Op<T>::u(Var1._range),
        *pVarR, Op<T>::l(pVarR->_range), Op<T>::u(pVarR->_range), PolCut<T>::LE, loc::pow, 0, &iExp );
    }
    // -- Nonconvex/Nonconcave Portion
    else{
      switch( Var1._img->options.BREAKPOINT_TYPE ){
        case PolImg<T>::Options::BIN:
        case PolImg<T>::Options::SOS2:{
          const unsigned nsubint = Var1.create_subdiv( Op<T>::l(Var1._range),
            Op<T>::u(Var1._range) ).size();
          if( nsubint > 2 ){
            struct dc{
              static std::pair<double,double> pow1
                ( const double x, const double*rusr, const int*iusr )
                { return std::make_pair( x<0?std::pow(x,*iusr):0., x<0?*iusr*std::pow(x,*iusr-1):0. ); }
              static std::pair<double,double> pow2
                ( const double x, const double*rusr, const int*iusr )
                { return std::make_pair( x>0?std::pow(x,*iusr):0., x>0?*iusr*std::pow(x,*iusr-1):0. ); }
            };
            PolVar<T>* pVar3 = Var1._img->_append_aux( Op<T>::pow( Var1._range, iExp ), true );
            PolVar<T>* pVar4 = Var1._img->_append_aux( Op<T>::pow( Var1._range, iExp ), true );
            Var1._img->_semilinear_cuts( pOp, Var1, Op<T>::l(Var1._range), Op<T>::u(Var1._range),
              *pVar3, PolCut<T>::GE, dc::pow1, 0, &iExp );
            Var1._img->_sandwich_cuts( pOp, Var1, Op<T>::l(Var1._range), Op<T>::u(Var1._range),
              *pVar3, Op<T>::l(pVarR->_range), 0., PolCut<T>::LE, dc::pow1, 0, &iExp );
            Var1._img->_semilinear_cuts( pOp, Var1, Op<T>::l(Var1._range), Op<T>::u(Var1._range),
              *pVar4, PolCut<T>::LE, dc::pow2, 0, &iExp );
              Var1._img->_sandwich_cuts( pOp, Var1, Op<T>::l(Var1._range), Op<T>::u(Var1._range),
              *pVar4, 0., Op<T>::u(pVarR->_range), PolCut<T>::GE, dc::pow2, 0, &iExp );
            Var1._img->_append_cut( pOp, PolCut<T>::EQ, 0., *pVarR, 1., *pVar3, -1., *pVar4, -1. );
            break;
          }
          // No break in order to append other "normal" cuts
        }
        case PolImg<T>::Options::NONE: default:{
          struct fct{ static std::pair<double,double> powoddfunc
            ( const double x, const double*rusr, const int*iusr )
            { return std::make_pair(
                ((*iusr-1)*x-(*iusr)*(*rusr))*std::pow(x,*iusr-1) + std::pow(*rusr,*iusr),
                (*iusr)*(*iusr-1)*(x-(*rusr))*std::pow(x,*iusr-2) ); }
          };
          double xJcc = Op<T>::u(Var1._range);
          xJcc = Var1._img->_newton( Op<T>::l(Var1._range), Op<T>::l(Var1._range), 0.,
            fct::powoddfunc, Var1._img->options.ROOT_TOL, Var1._img->options.ROOT_MAXIT, &xJcc, &iExp );
          if( mc::isequal(xJcc,Op<T>::l(Var1._range)) ){
            Var1._img->_semilinear_cuts( pOp, Var1, Op<T>::l(Var1._range), Op<T>::u(Var1._range),
              *pVarR, PolCut<T>::LE, loc::pow, 0, &iExp );
          }
          else{
            Var1._img->_sandwich_cuts( pOp, Var1, Op<T>::l(Var1._range), xJcc, *pVarR,
              Op<T>::l(pVarR->_range), Op<T>::u(pVarR->_range), PolCut<T>::LE, loc::pow, &xJcc, &iExp );
          }

          double xJcv = Op<T>::l(Var1._range);
          xJcv = Var1._img->_newton( Op<T>::u(Var1._range), 0., Op<T>::u(Var1._range),
            fct::powoddfunc, Var1._img->options.ROOT_TOL, Var1._img->options.ROOT_MAXIT, &xJcv, &iExp );
          if( mc::isequal(xJcv,Op<T>::u(Var1._range)) ){
            Var1._img->_semilinear_cuts( pOp, Var1, Op<T>::l(Var1._range), Op<T>::u(Var1._range),
              *pVarR, PolCut<T>::GE, loc::pow, 0, &iExp );
          }
          else{
            Var1._img->_sandwich_cuts( pOp, Var1, xJcv, Op<T>::u(Var1._range), *pVarR,
              Op<T>::l(pVarR->_range), Op<T>::u(pVarR->_range), PolCut<T>::GE, loc::pow, &xJcv, &iExp );
          }
          break;
        }
      }
    }
  }

  // Negative exponent term
  else if( iExp < 0 ){
    // -- Convex Case
    if( !(iExp%2) || Op<T>::l(Var1._range) > 0. ){
      Var1._img->_semilinear_cuts( pOp, Var1, Op<T>::l(Var1._range), Op<T>::u(Var1._range),
        *pVarR, PolCut<T>::LE, loc::pow, 0, &iExp );
      Var1._img->_sandwich_cuts( pOp, Var1, Op<T>::l(Var1._range), Op<T>::u(Var1._range),
        *pVarR, Op<T>::l(pVarR->_range), Op<T>::u(pVarR->_range), PolCut<T>::GE, loc::pow, 0, &iExp );
    }
    // -- Concave Case
    else{
      Var1._img->_semilinear_cuts( pOp, Var1, Op<T>::l(Var1._range), Op<T>::u(Var1._range),
        *pVarR, PolCut<T>::GE, loc::pow, 0, &iExp );
      Var1._img->_sandwich_cuts( pOp, Var1, Op<T>::l(Var1._range), Op<T>::u(Var1._range),
        *pVarR, Op<T>::l(pVarR->_range), Op<T>::u(pVarR->_range), PolCut<T>::LE, loc::pow, 0, &iExp );
    }
  }

  else
    throw typename PolImg<T>::Exceptions( PolImg<T>::Exceptions::NOTALLOWED );
    
  //std::cout << "pow( const PolVar<T>&Var1, const PolVar<T>&Var2 ):\n";
  //std::cout << *Var1._img;
  return *pVarR;
}

template <typename T>
inline PolVar<T>
cheb
( const PolVar<T>&Var1, const PolVar<T>&Var2 )
{
#ifdef MC__POLIMG_CHECK
  if( !Var2_.var.cst() && Var2_.var.num().t != FFNum::INT && Var2_.var.num().n < 2 )
    throw typename PolImg<T>::Exceptions( PolImg<T>::Exceptions::NOTALLOWED );
#endif
  const unsigned iOrd = Var2._var.num().n;
  if( Var1._var.cst() ) return mc::cheb( Var1.var().num().val(), iOrd );

  FFGraph* dag = Var1._var.dag();
#ifdef MC__POLIMG_CHECK
  if( !dag || !dag->curOp() )
    throw typename PolImg<T>::Exceptions( PolImg<T>::Exceptions::NOTALLOWED );
#endif

  FFVar* pFFVarR = dag->curOp()->pres;
  //FFOp* pOp = pFFVarR->ops().first;
#ifdef MC__POLIMG_DEBUG_CHEB
  std::cout << "X = " << Var1._range << std::endl;
  std::cout << "Cheb(" << Var1._range << "," << iOrd << ") = " << Op<T>::cheb( Var1._range, iOrd ) << std::endl;
#endif
  PolVar<T>* pVarR = Var1._img->_append_var( pFFVarR, Op<T>::cheb( Var1._range, iOrd ), true );
  if( pVarR->cuts() ) return *pVarR;
  pVarR->set_cuts();

  // Constant term
  //if( Var1._var.cst() ){
  //  Var1._img->_append_cut( pOp, PolCut<T>::EQ, mc::cheb(Var1._var.num().val(),iOrd), *pVarR, 1. );
  //  return *pVarR;
  //}

  //throw typename PolImg<T>::Exceptions( PolImg<T>::Exceptions::NOTALLOWED );
  return *pVarR;
}

} // namespace mc

#include "mcop.hpp"

namespace mc
{

//! @brief Specialization of the mc::Op templated structure to allow usage of type mc::PolVar as template type in other MC++ classes, such as mc::Taylor
template <> template<typename T> struct Op< mc::PolVar<T> >
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
  static PV xlog(const PV& x){ return x*mc::log(x); }
  static PV fabs(const PV& x)
    { throw std::runtime_error("operation not permitted"); }
//  { return mc::fabs(x); }
  static PV exp (const PV& x){ return mc::exp(x);  }
  static PV sin (const PV& x)
    { throw std::runtime_error("operation not permitted"); }
//  { return mc::sin(x);  }
  static PV cos (const PV& x)
    { throw std::runtime_error("operation not permitted"); }
//  { return mc::cos(x);  }
  static PV tan (const PV& x)
    { throw std::runtime_error("operation not permitted"); }
//  { return mc::tan(x);  }
  static PV asin(const PV& x)
    { throw std::runtime_error("operation not permitted"); }
//  { return mc::asin(x); }
  static PV acos(const PV& x)
    { throw std::runtime_error("operation not permitted"); }
//  { return mc::acos(x); }
  static PV atan(const PV& x)
    { throw std::runtime_error("operation not permitted"); }
//  { return mc::atan(x); }
  static PV erf (const PV& x)
    { throw std::runtime_error("operation not permitted"); }
//  { return mc::erf(x);  }
  static PV erfc(const PV& x)
    { throw std::runtime_error("operation not permitted"); }
//  { return mc::erfc(x); }
  static PV fstep(const PV& x)
    { throw std::runtime_error("operation not permitted"); }
//  { return mc::fstep(x); }
  static PV bstep(const PV& x)
    { throw std::runtime_error("operation not permitted"); }
//  { return mc::bstep(x); }
  static PV min (const PV& x, const PV& y)
    { throw std::runtime_error("operation not permitted"); }
//  { return mc::min(x,y);  }
  static PV max (const PV& x, const PV& y)
    { throw std::runtime_error("operation not permitted"); }
//  { return mc::max(x,y);  }
  static PV arh (const PV& x, const double k){ return mc::exp(-k/x); }
  static PV pow (const PV& x, const PV& y){ return mc::pow(x,y); }
  static PV cheb(const PV& x, const PV& y){ return mc::cheb(x,y); }
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
