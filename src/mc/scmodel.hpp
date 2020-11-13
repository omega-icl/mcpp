// Copyright (C) 2009-2016 Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the Eclipse Public License.

/*!
\page page_CHEBYSHEV Chebyshev Model Arithmetic for Factorable Functions
\author Jai Rajyaguru, Mario E. Villanueva, Beno&icirc;t Chachuat

A \f$q\f$th-order Chebyshev model of a Lipschitz-continuous function \f$f:\mathbb{R}^n\to\mathbb{R}\f$ on the domain \f$D\f$, consists of a \f$q^{\rm th}\f$-order multivariate polynomial \f$\mathcal P\f$ in Chebyshev basis , plus a remainder term \f$\mathcal R\f$, so that
\f{align*}
  f({x}) \in \mathcal P({x}-\hat{x}) \oplus \mathcal R, \quad \forall {x}\in D.
\f}
The polynomial part \f$\mathcal P\f$ is propagated symbolically and accounts for functional dependencies. The remainder term \f$\mathcal R\f$, on the other hand, is traditionally computed using interval analysis [Brisebarre & Joldes, 2010]; see figure below. More generally, convex/concave bounds or an ellipsoidal enclosure can be computed for the remainder term of vector-valued functions too. In particular, it can be established that the remainder term has convergence order (no less than) \f$q+1\f$ with respect to the diameter of the domain set \f$D\f$ under mild conditions [Bompadre <I>et al.</I>, 2012].

<CENTER><TABLE BORDER=0>
<TR>
<TD>\image html Chebyshev_model.png</TD>
</TR>
</TABLE></CENTER>

The classes mc::SCModel and mc::SCVar provide an implementation of Chebyshev model arithmetic. We note that mc::SCModel / mc::SCVar is <b>not a verified implementation</b> in the sense that rounding errors are not accounted for in propagating the coefficients in the multivariate polynomial part, which are treated as floating-point numbers.

The implementation of mc::SCModel and mc::SCVar relies on the operator/function overloading mechanism of C++. This makes the computation of Chebyshev models both simple and intuitive, similar to computing function values in real arithmetics or function bounds in interval arithmetic (see \ref page_INTERVAL). Moreover, mc::SCVar can be used as the template parameter of other available types in MC++; for instance, mc::SCVar can be used in order to propagate the underlying interval bounds in mc::McCormick. Likewise, mc::SCVar can be used as the template parameter of the types fadbad::F, fadbad::B and fadbad::T of <A href="http://www.fadbad.com/fadbad.html">FADBAD++</A> for computing Chebyshev models of either the partial derivatives or the Chebyshev coefficients of a factorable function (see \ref sec_CHEBYSHEV_fadbad).

mc::SCModel and mc::SCVar themselves are templated in the type used to propagate bounds on the remainder term. By default, mc::SCModel and mc::SCVar can be used with the non-verified interval type mc::Interval of MC++. For reliability, however, it is strongly recommended to use verified interval arithmetic such as <A href="http://www.ti3.tu-harburg.de/Software/PROFILEnglisch.html">PROFIL</A> (header file <tt>mcprofil.hpp</tt>) or <A href="http://www.math.uni-wuppertal.de/~xsc/software/filib.html">FILIB++</A> (header file <tt>mcfilib.hpp</tt>). As already noted, convex/concave bounds on the remainder term can also be propagated by using the type mc::McCormick of MC++, thereby enabling McCormick-Chebyshev models.

As well as propagating Chebyshev models for factorable functions, mc::SCModel and mc::SCVar provide support for computing bounds on the Chebyshev model range (multivariate polynomial part). We note that computing exact bounds for multivariate polynomials is a hard problem in general. Instead, a number of computationally tractable, yet typically conservative, bounding approaches are implemented in mc::SCModel and mc::SCVar, which include:
- Bounding every monomial term independently and adding these bounds;
- Bounding the first- and diagonal second-order terms exactly and adding bounds for the second-order off-diagonal and higher-order terms computed independently [Lin & Stadtherr, 2007];
- Bounding the terms up to order 2 based on an eigenvalue decomposition of the corresponding Hessian matrix and adding bounds for the higher-order terms computed independently;
- Expressing the multivariate polynomial in Bernstein basis, thereby providing bounds as the minimum/maximum among all Bernstein coefficients [Lin & Rokne, 1995; 1996].
.

Examples of Chebyshev models (blue lines) constructed with mc::SCModel and mc::SCVar are shown on the figure below for the factorable function \f$f(x)=x \exp(-x^2)\f$ (red line) for \f$x\in [-0.5,1]\f$. Also shown on these plots are the interval bounds computed from the Chebyshev models.

<CENTER><TABLE BORDER=0>
<TR>
<TD>\image html CM-1D.png</TD>
</TR>
</TABLE></CENTER>


\section sec_CHEBYSHEV_I How do I compute a Chebyshev model with interval remainder bound of a factorable function?

Suppose we want to compute a 4th-order Chebyshev model for the real-valued function \f$f(x,y)=x\exp(x+y^2)-y^2\f$ with \f$(x,y)\in [1,2]\times[0,1]\f$. For simplicity, bounds on the remainder terms are computed using the default interval type mc::Interval here:

\code
      #include "interval.hpp"
      #include "cmodel.hpp"
      typedef mc::Interval I;
      typedef mc::SCModel<I> CM;
      typedef mc::SCVar<I> CV;
\endcode

First, the number of independent variables in the factorable function (\f$x\f$ and \f$y\f$ here) as well as the order of the Chebyshev model (4th order here) are specified by defining an mc::SCModel object as:

\code
      CM mod( 2, 4 );
\endcode

Next, the variables \f$x\f$ and \f$y\f$ are defined as follows:

\code
      CV X( &mod, 0, I(1.,2.) );
      CV Y( &mod, 1, I(0.,1.) );
\endcode

Essentially, the first line means that <tt>X</tt> is a variable of class mc::SCVar, participating in the Chebyshev model <tt>mod</tt>, belonging to the interval \f$[1,2]\f$, and having index 0 (indexing in C/C++ start at 0 by convention!). The same holds for the Chebyshev variable <tt>Y</tt>, participating in the model <tt>mod</tt>, belonging to the interval \f$[0,1]\f$, and having index 1.

Having defined the variables, a Chebyshev model of \f$f(x,y)=x\exp(x+y^2)-y^2\f$ on \f$[1,2]\times[0,1]\f$ at the mid-point \f$(\frac{3}{2},\frac{1}{2})\f$ is simply computed as:

\code
      CV F = X*exp(X+pow(Y,2))-pow(Y,2);
\endcode

This model can be displayed to the standard output as:

\code
      std::cout << "f Chebyshev model: " << F << std::endl;
\endcode

which produces the following output:

\verbatim
f Chebyshev model: 
   a0    =  8.38199e+00     0  0
   a1    =  1.90755e+00     0  1
   a2    =  3.59621e+00     1  0
   a3    =  7.47482e-01     0  2
   a4    =  9.00782e-01     1  1
   a5    =  6.30186e-01     2  0
   a6    =  1.56945e-01     0  3
   a7    =  3.35238e-01     1  2
   a8    =  1.55141e-01     2  1
   a9    =  6.67468e-02     3  0
   a10   =  3.49519e-02     0  4
   a11   =  6.58449e-02     1  3
   a12   =  6.04330e-02     2  2
   a13   =  1.80397e-02     3  1
   a14   =  5.41191e-03     4  0
   R     =  [ -2.09182e+00 :  2.22652e+00 ]
   B     =  [ -1.02564e+01 :  3.93973e+01 ]
\endverbatim

<tt>a0</tt>,...,<tt>a14</tt> refer to the coefficients of the monomial terms in the Chebyshev model, with the corresponding variable orders given in the subsequent columns. The remainder term as well as the Chebyshev model range estimator are reported next.

Other operations involve retreiving the remainder bound, centering the remainder term in a Chebyshev model, or computing the value of its polynomial part at a given point:

\code
      I B = F.B();
      F.C();
      double x[2] = { 0.5, 1.5 };
      double Pval = F.P( x );
\endcode

See the documentations of mc::SCModel and mc::SCVar for a complete list of member functions. 


\section sec_CHEBYSHEV_fct Which functions are overloaded for Chebyshev model arithmetic?

mc::SCVar overloads the usual functions <tt>exp</tt>, <tt>log</tt>, <tt>sqr</tt>, <tt>sqrt</tt>, <tt>pow</tt>, <tt>inv</tt>, <tt>cos</tt>, <tt>sin</tt>, <tt>tan</tt>, <tt>acos</tt>, <tt>asin</tt>, <tt>atan</tt>. Unlike mc::Interval and mc::McCormick, the functions <tt>min</tt>, <tt>max</tt> and <tt>fabs</tt> are not overloaded in mc::SCVar as they are nonsmooth. Moreover, mc::SCVar defines the following functions:
- <tt>inter(x,y,z)</tt>, computing a Chebyshev model of the intersection \f$x = y\cap z\f$ of two Chebyshev models and returning true/false if the intersection is nonempty/empty. With Chebyshev models \f$\mathcal P_y\oplus\mathcal R_y\f$ and \f$\mathcal P_z\oplus\mathcal R_z\f$, this intersection is computed as follows:
\f{align*}
  \mathcal P_{x} =\ & (1-\eta) \mathcal P_y^{\rm C} + \eta \mathcal P_z^{\rm C}\\
  \mathcal R_{x} =\ & [\mathcal R_y^{\rm C}\oplus\eta\mathcal{B}(\mathcal P_y^{\rm C}-\mathcal P_z^{\rm C})] \cap [\mathcal R_z^{\rm C}\oplus (1-\eta)\mathcal{B}(\mathcal P_z^{\rm C}-\mathcal P_y^{\rm C})]\,.
\f}
with \f$\mathcal{B}(\cdot)\f$ the Chebyshev model range bounder, and \f$\eta\f$ a real scalar in \f$[0,1]\f$. Choosing \f$\eta=1\f$ amounts to setting the polynomial part \f$\mathcal P_{x}\f$ as \f$\mathcal P_y\f$, whereas \f$\eta=0\f$ sets \f$\mathcal P_{x}\f$ as \f$\mathcal P_z\f$. The parameter \f$\eta\f$ can be defined in mc::SCModel::Options::REF_POLY.
- <tt>hull(x,y)</tt>, computing a Chebyshev model of the union \f$x = y\cup z\f$ of two Chebyshev models. With Chebyshev models \f$\mathcal P_y\oplus\mathcal R_y\f$ and \f$\mathcal P_z\oplus\mathcal R_z\f$, this union is computed as follows:
\f{align*}
  \mathcal P_{x} =\ & (1-\eta) \mathcal P_y^{\rm C} + \eta \mathcal P_z^{\rm C}\\
  \mathcal R_{x} =\ & {\rm hull}\{\mathcal R_y^{\rm C}\oplus\eta\mathcal{B}(\mathcal P_y^{\rm C}-\mathcal P_z^{\rm C}), \mathcal R_z^{\rm C}\oplus (1-\eta)\mathcal{B}(\mathcal P_z^{\rm C}-\mathcal P_y^{\rm C})\}\,.
\f}
with \f$\mathcal{B}(\cdot)\f$ and \f$\eta\f$ as previously.


\section sec_CHEBYSHEV_opt How are the options set for the computation of a Chebyshev model?

The class mc::SCModel has a public member called mc::SCModel::options that can be used to set/modify the options; e.g.,

\code
      model.options.BOUNDER_TYPE = CM::Options::EIGEN;
      model.options.SCALE_VARIABLES = true;
\endcode

The available options are the following:

<TABLE border="1">
<CAPTION><EM>Options in mc::SCModel::Options: name, type and description</EM></CAPTION>
     <TR><TH><b>Name</b>  <TD><b>Type</b><TD><b>Default</b>
         <TD><b>Description</b>
     <TR><TH><tt>BOUNDER_TYPE</tt> <TD><tt>mc::SCModel::Options::BOUNDER</tt> <TD>mc::SCModel::Options::LSB
         <TD>Chebyshev model range bounder.
     <TR><TH><tt>BOUNDER_ORDER</tt> <TD><tt>unsigned int</tt> <TD>0
         <TD>Order of Bernstein polynomial for Chebyshev model range bounding, when mc::SCModel::options::BOUNDER_TYPE = mc::SCModel::options::BERNSTEIN is selected. Only values greater than the actual Chebyshev model order are accounted for; see [Lin & Rokne, 1996].
     <TR><TH><tt>REF_POLY</tt> <TD><tt>double</tt> <TD>0.
         <TD>Scalar in \f$[0,1]\f$ related to the choice of the polynomial part in the overloaded functions mc::inter and mc::hull (see \ref sec_CHEBYSHEV_fct). A value of 0. amounts to selecting the polynomial part of the left operand, whereas a value of 1. selects the right operand.
     <TR><TH><tt>DISPLAY_DIGITS</tt> <TD><tt>unsigned int</tt> <TD>5
         <TD>Number of digits in output stream for Chebyshev model coefficients.
</TABLE>


\section sec_CM_err Errors What errors can I encounter during computation of a Chebyshev model?

Errors are managed based on the exception handling mechanism of the C++ language. Each time an error is encountered, a class object of type mc::SCModel::Exceptions is thrown, which contains the type of error. It is the user's responsibility to test whether an exception was thrown during the computation of a Chebyshev model, and then make the appropriate changes. Should an exception be thrown and not caught by the calling program, the execution will abort.

Possible errors encountered during the computation of a Chebyshev model are:

<TABLE border="1">
<CAPTION><EM>Errors during the Computation of a Chebyshev Model</EM></CAPTION>
     <TR><TH><b>Number</b> <TD><b>Description</b>
     <TR><TH><tt>1</tt> <TD>Division by zero
     <TR><TH><tt>2</tt> <TD>Failed to compute eigenvalue decomposition in range bounder SCModel::Options::EIGEN
     <TR><TH><tt>3</tt> <TD>Failed to compute the maximum gap between a univariate term and its Bernstein model
     <TR><TH><tt>-1</tt> <TD>Number of variable in Chebyshev model must be nonzero
     <TR><TH><tt>-2</tt> <TD>Failed to construct Chebyshev variable
     <TR><TH><tt>-3</tt> <TD>Chebyshev model bound does not intersect with bound in template parameter arithmetic
     <TR><TH><tt>-4</tt> <TD>Operation between Chebyshev variables linked to different Chebyshev models
     <TR><TH><tt>-5</tt> <TD>Maximum size of Chebyshev model reached (monomials indexed as unsigned int)
     <TR><TH><tt>-33</tt> <TD>Feature not yet implemented in mc::SCModel
</TABLE>

Moreover, exceptions may be thrown by the template parameter class itself.


\section sec_CM_refs References

- Brisebarre, N., and M. Joldes, <A href="http://hal.archives-ouvertes.fr/docs/00/48/17/37/PDF/RRLIP2010-13.pdf">Chebyshev Interpolation Polynomial-based Tools for Rigorous Computing</A>, <i>Research Report No RR2010-13</i>, Ecole Normale Sup&eaccute;rieure de Lyon, Unit&eaccute; Mixte de Recherche CNRS-INRIA-ENS LYON-UCBL No 5668, 2010
- T Dzetkulic, <A HREF="http://dx.doi.org/10.1007%2Fs11075-014-9889-x">Rigorous integration of non-linear ordinary differential equations in Chebyshev basis</A>. <I>Numerical Algorithms</I> <B>69</B>(1):183–205, 2015
J Rajyaguru, ME Villanueva, B Houska, B Chachuat, <A href="">Chebyshev Model Arithmetic for Factorable Functions</A>, <I>Journal of Global Optimization</I> <B>68</B>:413-438, 2017
- J Rajyaguru, ME Villanueva, B Houska, B Chachuat, <A href="https://doi.org/10.1007/s10898-016-0474-9">Chebyshev Model Arithmetic for Factorable Functions</A>, <I>Journal of Global Optimization</I> <B>68</B>:413-438, 2017
- LN Trefethen,<A HREF="http://www.chebfun.org/ATAP/"><I>Approximation Theory and Approximation Practice</I></A>. SIAM, Philadelphia (PA), 2013
.
*/

#ifndef MC__SCMODEL_H
#define MC__SCMODEL_H

#include "spolymodel.hpp"
#include "mcop.hpp"

#undef  MC__SCMODEL_DEBUG
#undef  MC__SCMODEL_DEBUG_SCALE
#define MC__SCMODEL_CHECK
#undef  MC__SCMODEL_CHECK_PMODEL
#undef  MC__SCVAR_DEBUG_EXP
#undef  MC__SCVAR_DEBUG_BERNSTEIN
#define MC__SCMODEL_USE_PROD

//#undef MC__SCVAR_FABS_SQRT
//#undef MC__SCVAR_FORCE_REM_DERIV
//#undef MC__SCVAR_SPARSE_PRODUCT_FULL

namespace mc
{

template <typename T> class SCVar;

//! @brief C++ class for the computation of sparse Chebyshev models for factorable function: environment
////////////////////////////////////////////////////////////////////////
//! mc::SCModel is a C++ class for definition of Chebyshev model
//! environment in sparse format. Propagation of Chebyshev models for
//! factorable functions is via the C++ class mc::SCVar. The template
//! parameter corresponds to the type
//! used to propagate the remainder bound.
////////////////////////////////////////////////////////////////////////
template <typename T>
class SCModel
////////////////////////////////////////////////////////////////////////
{
  friend class SCVar<T>;
  template <typename U> friend class SCModel;
  template <typename U> friend SCVar<U> pow
    ( SCVar<U> const&, int const );
  template <typename U> friend SCVar<U> sqr
    ( SCVar<U> const& );

  typedef std::map< SPolyMon, double, lt_SPolyMon > t_coefmon;
  typedef std::set< unsigned > t_var;
  
protected:

  //!brief Unit ball in T arithmetic
  static T _TOne;

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

  //! @brief Coefficients in Chebyshev interpolant of univariate functions
  std::vector<double> _coefinterp;
      
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
  std::vector<double>& _resize_coefinterp
    ()
    { _coefinterp.resize( _maxord + options.INTERP_EXTRA + 1 );
      return _coefinterp; }

public:
  /** @addtogroup SCHEBYSHEV Chebyshev Model Arithmetic for Factorable Functions
   *  @{
   */
  //! @brief Constructor of Sparse Chebyshev model environment for maximal order <tt>maxord</tt>
  SCModel
    ( const unsigned maxord=3, const unsigned nvarres=0 )
    : _maxord( maxord ), _maxvar( 0 )
    { _coefinterp.reserve( _maxord + options.INTERP_EXTRA + 1 );
      if( !nvarres ) return;
      _bndvar.reserve( nvarres );
      _refvar.reserve( nvarres );
      _scalvar.reserve( nvarres ); }

  //! @brief Destructor of Sparse Chebyshev model environment
  ~SCModel()
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
/*
  //! @brief Polynomial range bounder using specified bounder <a>type</a> with basis functions <a>bndbasis</a> in U arithmetic and monomial coefficients <a>coefmon</a> in C arithmetic
  template <typename C, typename U> U get_bound
    ( const C*coefmon, const U*const*bndbasis, const U*bndrem, const int type,
      const std::set<unsigned>&ndxmon=std::set<unsigned>() )
    { return( bndrem? _polybound( coefmon, bndbasis, type, ndxmon ) + *bndrem:
                      _polybound( coefmon, bndbasis, type, ndxmon ) ); }
*/
  //! @brief Exceptions of mc::SCModel
  class Exceptions
  {
  public:
    //! @brief Enumeration type for SCModel exception handling
    enum TYPE{
      DIV=1,	//!< Division by zero scalar
      INV,	//!< Inverse operation with zero in range
      LOG,	//!< Log operation with non-positive numbers in range
      SQRT,	//!< Square-root operation with negative numbers in range
      TAN,	//!< Tangent operation with (k+1/2)·PI in range
      ACOS,	//!< Sine/Cosine inverse operation with range outside [-1,1]
      EIGEN,	//!< Failed to compute eigenvalue decomposition in range bounder SCModel::Options::EIGEN
      INIT=-1,	//!< Failed to construct Chebyshev variable
      INCON=-2, //!< Chebyshev model bound does not intersect with bound in template parameter arithmetic
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
      case TAN:
        return "mc::SCModel\t Tangent operation with (k+1/2)·PI in range";
      case ACOS:
        return "mc::SCModel\t Sine/Cosine inverse operation with range outside [-1,1]";
      case EIGEN:
        return "mc::SCModel\t Range bounder with eigenvalue decomposition failed";
      case INIT:
        return "mc::SCModel\t Chebyshev variable initialization failed";
      case INCON:
        return "mc::SCModel\t Inconsistent bounds with template parameter arithmetic";
      case SCMODEL:
        return "mc::SCModel\t Operation between Chebyshev variables in different Chebyshev model environment not allowed";
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

  //! @brief Options of mc::SCModel
  struct Options
  {
    //! @brief Constructor of mc::SCModel::Options
    Options():
      INTERP_EXTRA(0), INTERP_THRES(1e2*machprec()), BOUNDER_TYPE(LSB),
      BOUNDER_ORDER(0), MIXED_IA(true), MIN_FACTOR(0.), REF_POLY(0.),
      DISPLAY_DIGITS(5)
      {}
    //! @brief Copy constructor of mc::SCModel::Options
    template <typename U> Options
      ( U&options )
      : INTERP_EXTRA( options.INTERP_EXTRA ),
        INTERP_THRES( options.INTERP_THRES ),
        BOUNDER_TYPE( options.BOUNDER_TYPE ),
        BOUNDER_ORDER( options.BOUNDER_ORDER ),
        MIN_FACTOR( options.MIN_FACTOR ),
        MIXED_IA( options.MIXED_IA ),
        REF_POLY(options.REF_POLY),
        DISPLAY_DIGITS(options.DISPLAY_DIGITS)
      {}
    //! @brief Assignment of mc::SCModel::Options
    template <typename U> Options& operator =
      ( U&options ){
        INTERP_EXTRA     = options.INTERP_EXTRA;
        INTERP_THRES     = options.INTERP_THRES;
        BOUNDER_TYPE     = (BOUNDER)options.BOUNDER_TYPE;
        BOUNDER_ORDER    = options.BOUNDER_ORDER;
        MIN_FACTOR       = options.MIN_FACTOR;
        MIXED_IA         = options.MIXED_IA;
        REF_POLY         = options.REF_POLY;
        DISPLAY_DIGITS   = options.DISPLAY_DIGITS;
        return *this;
      }
    //! @brief Chebyshev model range bounder option
    enum BOUNDER{
      NAIVE=0,	//!< Naive polynomial range bounder
      LSB,	//!< Lin & Stadtherr range bounder
      EIGEN,	//!< Eigenvalue decomposition-based bounder
      BERNSTEIN //!< Bernstein range bounder
    };
    //! @brief Extra terms in chebyshev interpolation of univariates: 0-Chebyshev interpolation of order MAXORD; extra terms allow approximation of Chebyshev truncated series
    unsigned INTERP_EXTRA;
    //! @brief Threshold for coefficient values in Chebyshev expansion for bounding of transcendental univariates
    double INTERP_THRES;
    //! @brief Chebyshev model range bounder - See \ref sec_CHEBYSHEV_opt
    BOUNDER BOUNDER_TYPE;
    //! @brief Order of Bernstein polynomial for Chebyshev model range bounding (no less than Chebyshev model order!). Only if mc::SCModel::options::BOUNDER_TYPE is set to mc::SCModel::options::BERNSTEIN.
    unsigned BOUNDER_ORDER;
    //! @brief Array of Chebyshev model range bounder names (for display)
    static const std::string BOUNDER_NAME[5];
    //! @brief Whether to intersect internal bounds with underlying bounds in T arithmetics
    bool MIXED_IA;
    //! @brief Threshold for monomial coefficients below which the term is removed and appended to the remainder term.
    double MIN_FACTOR;
    //! @brief Scalar in \f$[0,1]\f$ related to the choice of the polynomial part in the overloaded functions mc::inter and mc::hull (see \ref sec_CHEBYSHEV_fct). A value of 0. amounts to selecting the polynomial part of the left operand, whereas a value of 1. selects the right operand.
    double REF_POLY;
    //! @brief Number of digits in output stream for Chebyshev model coefficients.
    unsigned DISPLAY_DIGITS;
  } options;
  /** @} */

private:

  //! @brief Get Chebyshev basis functions in U arithmetic for variable <a>bndvar</a>
  template <typename U> static U* _get_bndpow
    ( const unsigned maxord, U const& bndvar, double const ref, double const scal );

  //! @brief Get Chebyshev basis functions in U arithmetic for [-1,1] scaled variable <a>bndvar</a>
  template <typename U> static U* _get_bndpow
    ( const unsigned maxord, U const& bndvar );

  //! @brief Prototype real-valued function for interpolation
  typedef double (puniv)
    ( double const x );

  //! @brief Construct Chebyshev interpolating polynomial coefficient <a>coefmon</a> for univariate <a>f</a>
  static void _interpolation
    ( std::vector<double>&coefmon, const unsigned maxord, T const& X, puniv f );

  //! @brief Construct Chebyshev interpolating polynomial coefficient <a>coefmon</a> for univariate <a>f</a>
  static void _interpolation
    ( std::vector<double>&coefmon, double const TOL, unsigned& nord,
      T const& X, puniv f );

  //! @brief Apply Chebyshev composition to variable <a>CVI</a> using the coefficients <a>coefmon</a> of the outer function
  template <typename U> static SCVar<T> _composition
    ( std::vector<U> const& coefouter, unsigned const maxord,
      SCVar<T> const& CVinner );

  //! @brief Recursive calculation of nonnegative integer powers
  SCVar<T> _intpow
    ( SCVar<T> const& CV, int const n ) const;

  //! @brief Polynomial range bounder - Naive approach
  template <typename C, typename U> U _polybound_naive
    ( std::map<SPolyMon,C,lt_SPolyMon> const& coefmon, U const* const* bndbasis,
      unsigned const minord=0 ) const;

  //! @brief Polynomial range bounder - Lin & Stadtherr approach
  template <typename C, typename U> U _polybound_LSB
    ( std::map<SPolyMon,C,lt_SPolyMon> const& coefmon, U const* const* bndbasis ) const;

  //! @brief Polynomial range bounder - eigenvalue decomposition approach
  template <typename C, typename U> U _polybound_eigen
    ( std::map<SPolyMon,C,lt_SPolyMon> const& coefmon, U const* const* bndbasis ) const;

  //! @brief Polynomial range bounder - Bernstein approach
  template <typename C, typename U> U _polybound_bernstein
    ( std::map<SPolyMon,C,lt_SPolyMon> const& coefmon, U const* const* bndbasis ) const;

  //! @brief Compute Bernstein coefficient for variable with exponents <tt>jexp</tt>, given coefficients in Chebyshev form <tt>coefmon</tt>, maximum order <tt>maxord</tt> and transformation coefficients <a>trmat</a>
  template <typename C> C _coef_bernstein
    ( std::map<SPolyMon,C,lt_SPolyMon> const& coefmon, unsigned const* jexp, double const* trmat ) const;

  //! @brief Compute Bernstein tranformation coefficient into Bernstein basis function <a>j</a> with order <a>n</a> from Chebyshev basis function <a>k</a>
  double _transform_bernstein
    ( const unsigned j, const unsigned k, const unsigned n ) const;

  //! @brief Polynomial range bounder using specified bounder <a>type</a> with basis functions <a>bndbasis</a> in U arithmetic and monomial coefficients <a>coefmon</a> in C arithmetic
  template <typename C, typename U> U _polybound
    ( std::map<SPolyMon,C,lt_SPolyMon> const& coefmon, U const *const* bndbasis,
      int const type );

  //! @brief Recursive product of univariate Chebyshev polynomials
  void _sprod1D
    ( std::map<unsigned,t_coefmon> const& sp1map,
      std::map<unsigned,t_coefmon> const& sp2map,
      t_coefmon& coefmon, double& coefrem, std::set<unsigned> const& ndxvar,
      typename t_var::const_iterator itvar ) const;

  //! @brief Scaling of Chebyshev coefficient maps
  void _sscal1D
    ( const t_coefmon& coefmon0, double const&coefscal,
      t_coefmon& coefmon ) const;

  //! @brief Lifting of Chebyshev coefficient maps
  void _slift1D
    ( const t_coefmon& coefmon0, double const& dscal,
      t_coefmon& coefmon ) const;

  //! @brief Lifting of Chebyshev coefficient maps
  void _slift1D
    ( t_coefmon const& coefmon0, double const& dscal,
      t_coefmon& coefmon, double& coefrem,
      typename t_var::const_iterator itvar, unsigned const ndxord ) const;

  //! @brief Lifting of Chebyshev coefficient maps
  void _slift1D
    ( const t_coefmon& coefmon0, double const& dscal,
      double& coefrem ) const;

  //! @brief Build 1D vector of Chebyshev coefficients
  void _svec1D
    ( typename t_var::const_iterator itvar,
      std::pair<SPolyMon,double> const& mon,
      std::map<unsigned,t_coefmon>& mapspoly ) const;

  //! @brief Display of recursive univariate Chebyshev polynomials
  void _sdisp1D
    ( const std::vector<SCVar<T>>& coefmon, const unsigned ndxvar,
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
    ( unsigned const ndxvar, std::pair<SPolyMon,double> const& coefmon,
      std::vector<SCVar<T>>& vec ) const;

  //! @brief Product of multivariate Chebyshev polynomials in sparse format
  void _sprod
    ( SCVar<T> const& CV1, SCVar<T> const& CV2, t_coefmon& coefmon,
      double&coefrem ) const;

  //! @brief Squaring of multivariate Chebyshev polynomials in sparse format
  void _ssqr
    ( SCVar<T> const& CV, t_coefmon& coefmon, double&coefrem ) const;
};

template <typename T>
inline
const std::string SCModel<T>::Options::BOUNDER_NAME[5]
  = { "NAIVE", "LSB", "EIGEN", "BERNSTEIN", "HYBRID" };

template <typename T>
inline
T SCModel<T>::_TOne
  = 2.*Op<T>::zeroone()-1.;

//! @brief C++ class for Chebyshev model computation of factorable function - Chebyshev model propagation
////////////////////////////////////////////////////////////////////////
//! mc::SCVar is a C++ class for propagation of Chebyshev models through
//! factorable functions. The template parameter corresponds to the
//! type used in computing the remainder bound.
////////////////////////////////////////////////////////////////////////
template <typename T>
class SCVar: public SPolyVar<T>
////////////////////////////////////////////////////////////////////////
{
  template <typename U> friend class SCVar;
  template <typename U> friend class SCModel;

  template <typename U> friend SCVar<U> operator-
    ( const SCVar<U>& );
  template <typename U> friend SCVar<U> operator*
    ( const SCVar<U>&, const SCVar<U>& );
  template <typename U> friend std::ostream& operator<<
    ( std::ostream&, const SCVar<U>& );

  template <typename U> friend U funcptr
    ( const U, const unsigned  );
  template <typename U> friend void interpolation
    ( double*, const SCVar<U>&, const unsigned  );
  template <typename U> friend SCVar<U> composition
    ( double const*, const SCVar<U>& );

  template <typename U> friend SCVar<U> inv
    ( const SCVar<U>& );
  template <typename U> friend SCVar<U> sqr
    ( const SCVar<U>& );
  template <typename U> friend SCVar<U> sqrt
    ( const SCVar<U>& );
  template <typename U> friend SCVar<U> exp
    ( const SCVar<U>& );
  template <typename U> friend SCVar<U> log
    ( const SCVar<U>& );
  template <typename U> friend SCVar<U> xlog
    ( const SCVar<U>& );
  template <typename U> friend SCVar<U> pow
    ( const SCVar<U>&, const int );
  template <typename U> friend SCVar<U> pow
    ( const SCVar<U>&, double const );
  template <typename U> friend SCVar<U> pow
    ( double const, const SCVar<U>& );
  template <typename U> friend SCVar<U> pow
    ( const SCVar<U>&, const SCVar<U>& );
  template <typename U> friend SCVar<U> prod
    ( const unsigned int, const SCVar<U>* );
  template <typename U> friend SCVar<U> monom
    ( const unsigned int, const SCVar<U>*, const unsigned* );
  template <typename U> friend SCVar<U> cheb
    ( const SCVar<U>&, const unsigned );
  template <typename U> friend SCVar<U> cos
    ( const SCVar<U>& );
  template <typename U> friend SCVar<U> sin
    ( const SCVar<U>& );
  template <typename U> friend SCVar<U> tan
    ( const SCVar<U>& );
  template <typename U> friend SCVar<U> acos
    ( const SCVar<U>& );
  template <typename U> friend SCVar<U> asin
    ( const SCVar<U>& );
  template <typename U> friend SCVar<U> atan
    ( const SCVar<U>& );
  template <typename U> friend SCVar<U> cosh
    ( const SCVar<U>& );
  template <typename U> friend SCVar<U> sinh
    ( const SCVar<U>& );
  template <typename U> friend SCVar<U> tanh
    ( const SCVar<U>& );
  template <typename U> friend SCVar<U> fabs
    ( const SCVar<U>& );
  template <typename U> friend SCVar<U> hull
    ( const SCVar<U>&, const SCVar<U>& );
  template <typename U> friend bool inter
    ( SCVar<U>&, const SCVar<U>&, const SCVar<U>& );

public:

  typedef std::map< SPolyMon, double, lt_SPolyMon > t_coefmon;
  typedef std::set< unsigned > t_var;
  
  using SPolyVar<T>::_TOne;
  using SPolyVar<T>::_ndxvar;
  using SPolyVar<T>::_coefmon;
  using SPolyVar<T>::_bndrem;
  using SPolyVar<T>::_bndT;
  using SPolyVar<T>::_set_bndT;
  using SPolyVar<T>::_unset_bndT;
  using SPolyVar<T>::_bndpol;
  using SPolyVar<T>::_set_bndpol;
  using SPolyVar<T>::_unset_bndpol;
  using SPolyVar<T>::_reinit;
  using SPolyVar<T>::_set;
  using SPolyVar<T>::set;
  using SPolyVar<T>::nord;
  using SPolyVar<T>::nvar;
  using SPolyVar<T>::bound;
  using SPolyVar<T>::_center;
  using SPolyVar<T>::display;
  
private:
  //! @brief Pointer to Chebyshev model environment
  SCModel<T> *_CM;

  //! @brief Original bounds on variable <tt>ivar</tt>
  T const&  _bndvar
    ( const unsigned ivar ) const
    { return _CM->_bndvar[ivar]; };

  //! @brief Reference point for variable <tt>ivar</tt> in Chebyshev model
  double _refvar
    ( const unsigned ivar ) const
    { return _CM->_refvar[ivar]; };

  //! @brief Scaling for variable <tt>ivar</tt> in Cheyshev model
  double _scalvar
    ( const unsigned ivar ) const
    { return _CM->_scalvar[ivar]; };

  //! @brief Set Chebyshev variable with index <a>ix</a> (starting from 0) and bounded by <a>X</a>
  SCVar<T>& _set
    ( const unsigned ivar, T const& X, const bool updMod=true );

  //! @brief Product of multivariate Chebyshev polynomials in sparse format
  void _sprod
    ( const SCVar<T>& CV1, const SCVar<T>& CV2, t_coefmon& coefmon,
      double& coefrem ) const
    { return _CM->_sprod( CV1, CV2, coefmon, coefrem ); }

  //! @brief Squaring of multivariate Chebyshev polynomials in sparse format
  void _ssqr
    ( const SCVar<T>& CV, t_coefmon& coefmon, double& coefrem ) const
    { return _CM->_ssqr( CV, coefmon, coefrem ); }

  //! @brief Array of Chebyshev interpolant coefficients
  std::vector<double>& _coefinterp() const
    { return _CM->_resize_coefinterp(); };

public:
  /** @addtogroup CHEBYSHEV Chebyshev Model Arithmetic for Factorable Functions
   *  @{
   */
  //! @brief Get pointer to linked Chebyshev model environment
  SCModel<T>* env() const
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

  //! @brief Constructor of Chebyshev variable for a real scalar
  SCVar
    ( SCModel<T>*CM );

  //! @brief Constructor of Chebyshev variable for a real scalar
  SCVar
    ( double const d=0., SCModel<T>*CM=0 );

  //! @brief Constructor of Chebyshev variable for a remainder bound
  SCVar
    ( T const& B, SCModel<T>*CM=0 );

  //! @brief Constructor of Chebyshev variable with index <a>ix</a> (starting from 0) and bounded by <a>X</a>
  SCVar
    ( SCModel<T>*CM, const unsigned ix, T const& X );
/*
  //! @brief Copy constructor of Chebyshev variable in different Chebyshev model environment (with implicit type conversion)
  template <typename U> SCVar
    ( SCModel<T>*&CM, SCVar<U> const& CV );

  //! @brief Copy constructor of Chebyshev variable in different Chebyshev model environment (with explicit type conversion as given by class member function <a>method</a>)
  template <typename U> SCVar
    ( SCModel<T>*&CM, SCVar<U> const& CV, T const&  (U::*method)() const );

  //! @brief Copy constructor of Chebyshev variable in different Chebyshev model environment (with explicit type conversion as given by non-class member function <a>method</a>)
  template <typename U> SCVar
    ( SCModel<T>*&CM, SCVar<U> const& CV, T (*method)( U const&  ) );
*/
  //! @brief Copy constructor of Chebyshev variable
  SCVar
    ( SCVar<T> const& CV )
    : SPolyVar<T>( CV ), _CM( CV._CM )
    {
#ifdef  MC__SCMODEL_CHECK_PMODEL
      if( _CM != dynamic_cast< SCModel<T>* >( SPolyVar<T>::_CM ) ) assert( false );
#endif
    }

  //! @brief Destructor of Chebyshev variable
  ~SCVar()
    {}

  //! @brief Set Chebyshev variable with index <a>ix</a> (starting from 0) and bounded by <a>X</a>
  SCVar<T>& set
    ( SCModel<T>*CM, const unsigned ix, T const& X )
    { set( CM ); _set( ix, X ); return *this; }

  //! @brief Set environment in Chebyshev model
  SCVar<T>& set
    ( SCModel<T>*CM )
    { _CM = CM; return *this; }
/*
  //! @brief Set Chebyshev model environment in Chebyshev variable to <tt>CM</tt>
  SCVar<T>& set
    ( SCModel<T>*CM, const bool reset=false )
    { if( CM != _CM ){ _CM = CM; SPolyVar<T>::set( _CM, reset ); } return *this; }
*/
  //! @brief Retreive bound on variable using default bounder in U arithmetic
  template <typename U>
  U bound
    ( const U*const*bndbasis, U const& bndrem ) const
    { return _polybound( bndbasis ) + bndrem; }

  //! @brief Retreive bound on variable using bounder <a>type</a> in U arithmetic
  template <typename U>
  U bound
    ( const U*const*bndbasis, U const& bndrem, const int type ) const
    { return _polybound( bndbasis, type ) + bndrem; }

  //! @brief Retreive bound on all terms with (total) order <tt>minord</tt> or higher in polynomial part
  T bndord
    ( const unsigned minord )
    const
    { if( !_CM ) return !minord && !_coefmon.empty() && !_coefmon.begin()->first.tord? _coefmon.begin()->second: 0.;
      return _CM->_polybound_naive( _coefmon, (T const* const*)0, minord ); }

  //! @brief Evaluate polynomial part at <tt>x</tt>
  double polynomial
    ( double const*x ) const;

  //! @brief Shortcut to mc::SCVar::polynomial
  double P
    ( double const*x ) const
    { return polynomial( x ); }

  //! @brief Return new Chebyshev variable with same multivariate polynomial part but zero remainder
  SCVar<T> polynomial
    ()
    const
    { SCVar<T> var = *this; var._bndrem = 0.; return var; }

  //! @brief Shortcut to mc::SCVar::polynomial
  SCVar<T> P
    ()
    const
    { return polynomial(); }

  //! @brief Center remainder term of Chebyshev variable
  SCVar<T>& center
    ()
    { _center(); return *this; }

  //! @brief Shortcut to mc::SCVar::center
  SCVar<T>& C
    ()
    { return center(); }

  //! @brief Get coefficient of constant term in Chebyshev variable. The value of this coefficient is reset to 0 if <tt>reset=true</tt>, otherwise it is left unmodified (default).
  double constant
    ( bool const reset=false );

  //! @brief Get coefficients of linear term for variable <tt>ivar</tt> in Chebyshev variable. The value of this coefficient is reset to 0 if <tt>reset=true</tt>, otherwise it is left unmodified (default).
  double linear
    ( unsigned const ivar, bool const reset=false );

  //! @brief Simplify the model by appending coefficient less than TOL as well as terms with order greater than or equal to ORD to the remainder term
  SCVar<T>& simplify
    ( double const TOL=0e0, int const TORD=-1 );
  //! @brief Scale coefficients in Chebyshev variable for the modified range <a>X</a> of variable i
  SCVar<T>& scale
    ( unsigned const i, T const& X );
  //! @brief Scale coefficients in Chebyshev variable for the modified variable ranges <a>X</a>
  SCVar<T>& scale
    ( T const* X );
  //! @brief Rescale coefficients in Chebyshev variable for their original variable ranges
  SCVar<T>::t_coefmon unscale
    ()
    const;    
  //! @brief Return new coefficient map in monomial basis representation
  t_coefmon to_monomial
    ( bool const scaled=false )
    const;

 /** @} */

  SCVar<T>& operator=
    ( SCVar<T> const& );
  SCVar<T>& operator=
    ( double const );
  SCVar<T>& operator=
    ( T const&  );
  template <typename U> SCVar<T>& operator+=
    ( SCVar<U> const& );
  template <typename U> SCVar<T>& operator+=
    ( U const&  );
  SCVar<T>& operator+=
    ( double const );
  template <typename U> SCVar<T>& operator-=
    ( const SCVar<U>& );
  template <typename U> SCVar<T>& operator-=
    ( U const&  );
  SCVar<T>& operator-=
    ( double const );
  SCVar<T>& operator*=
    ( const SCVar<T>& );
  SCVar<T>& operator*=
    ( double const );
  SCVar<T>& operator*=
    ( T const&  );
  SCVar<T>& operator/=
    ( SCVar<T> const& );
  SCVar<T>& operator/=
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
    { return _polybound( (const T*const*)0, type ); }

  //! @brief Polynomial range bounder using default bounder
  T _polybound
    () const
    { return _polybound( _CM? _CM->options.BOUNDER_TYPE: 0 ); }

  //! @brief Prototype real-valued function for interpolation
  typedef double (puniv)
    ( double const x );

  //! @brief Construct Chebyshev interpolating polynomial coefficient <a>coefmon</a> for univariate <a>f</a>
  void _interpolation
    ( std::vector<double>& coefmon, puniv f ) const
    { SCModel<T>::_interpolation( coefmon, maxord()+_CM->options.INTERP_EXTRA, bound(), f ); }

  //! @brief Construct Chebyshev interpolating polynomial coefficient <a>coefmon</a> for univariate <a>f</a>
  void _interpolation
    ( std::vector<double>& coefmon, double const TOL, unsigned& nord, puniv f ) const
    { SCModel<T>::_interpolation( coefmon, TOL, nord, bound(), f ); }

  //! @brief Apply Chebyshev composition to variable <a>CVI</a> using the coefficients <a>coefmon</a> of the outer function
  template <typename U> SCVar<T> _composition
    ( std::vector<U> const& coefouter ) const
    { return SCModel<T>::_composition( coefouter, maxord(), *this ); }

  //! @brief Scale current variable in order for its range to be within [-1,1], with <a>c</a> and <a>w</a> respectively the center and width, respectively, of the orginal variable range
  SCVar<T> _rescale
    ( double const w, double const c ) const
    { return( !isequal(w,0.)? (*this-c)/w: c ); }

  //! @brief Scale variable <a>X</a>
  template <typename U> static
  U _rescale
    ( U& X, double const w, double const c )
    { return( !isequal(w,0.)? (X-c)/w: c ); }

  //! @brief Return an array of Chebyshev variables representing the coefficients in the univariate polynomial for variable <a>ivar</a> only
  SCVar<T>* _single
    ( unsigned const ivar );

  //! @brief Cancel zero entries in coefficient map <a>coefmon</a>
  void _simplify
    ( t_coefmon& coefmon )
    const;

  //! @brief Scale coefficient map <a>coefmon</a> for the modified range <a>X</a> of variable i
  void _scale
    ( unsigned const i, T const& X, typename SCVar<T>::t_coefmon& coefmon )
    const;

  //! @brief Convert CHebyshev basis into monomial basis for variable <a>ivar</a> in coefficient map <a>coefmon</a>
  void _to_monomial
    ( unsigned const ivar, typename SCVar<T>::t_coefmon& coefmon )
    const;
};

////////////////////////////////// SCModel //////////////////////////////////////

template <typename T>
template <typename U>
inline U*
SCModel<T>::_get_bndpow
( unsigned const maxord, U const& bndvar, double const ref, double const scal )
{
  U *bndcheb = new U[maxord+1];
  bndcheb[0] = 1;
  if( maxord ) bndcheb[1] = ( bndvar - ref ) / scal;
  for( unsigned i=2; i<=maxord; i++ )
    bndcheb[i] = Op<U>::cheb( bndcheb[1], i );
 return bndcheb;
}

template <typename T>
template <typename U>
inline U*
SCModel<T>::_get_bndpow
( unsigned const maxord, U const& bndvar )
{
  U* bndcheb = new U[maxord+1];
  bndcheb[0] = 1;
  if( maxord ) bndcheb[1] = bndvar;
  for( unsigned i=2; i<=maxord; i++ )
    bndcheb[i] = Op<U>::cheb( bndcheb[1], i );
 return bndcheb;
}

template <typename T>
template <typename U>
inline U**
SCModel<T>::get_basis
( unsigned const maxord, U const* bndvar, bool const scaled ) const
{
  assert ( bndvar );
  U** bndcheb = new U*[_maxvar];
  for( unsigned i=0; i<_maxvar; i++){
    if( _setvar.find(i) == _setvar.end() ){
      bndcheb[i] = nullptr;
      continue;
    }
    bndcheb[i] = scaled? _get_bndpow( maxord, bndvar[i] ):
                         _get_bndpow( maxord, bndvar[i], _refvar[i], _scalvar[i] );
  }
  return bndcheb;
}

template <typename T>
template <typename U>
inline void
SCModel<T>::get_bndmon
( std::map<SPolyMon,U,lt_SPolyMon>& bndmon, U const* bndvar,
  bool const scaled )
const
{
  if( bndmon.empty() ) return;

  unsigned const nord = bndmon.rbegin()->first.tord;
  U** basis = get_basis( nord, bndvar, scaled );

#ifdef MC__SCMODEL_USE_PROD
  std::vector<U> Umon; Umon.reserve(_maxvar);
  auto it = bndmon.begin();
  if( !it->first.tord ){ it->second = 1; ++it; }
  for( ; it!=bndmon.end(); ++it ){
    Umon.clear();
    for( auto& [ivar,iord] : it->first.second )
      Umon.push_back( basis[ivar][iord] );
    it->second = prod( Umon.size(), Umon.data() );
  }
#else
  auto it = bndmon.begin();
  if( !it->first.tord ){ it->second = 1; ++it; }
  for( ; it!=bndmon.end(); ++it ){
    it->second = 1;
    for( auto& [ivar,iord] : it->first.second )
      it->second *= basis[ivar][iord];
  }
#endif

  for( unsigned j : _setvar ) delete[] basis[j];
  delete[] basis;
}

template <typename T>
inline void
SCModel<T>::_interpolation
( std::vector<double>& coefmon, unsigned const maxord, T const& X, puniv f )
{
  double b( Op<T>::mid(X) ), a( Op<T>::u(X)-b ), x[maxord+1], fx[maxord+1];
  double mulconst( PI/(2.*double(maxord+1)) );
  for( unsigned i(0); i<=maxord; i++ ){
    x[i]  = std::cos(mulconst*(2.*double(i)+1.));
    fx[i] = f( a*x[i]+b );
  }

  switch( maxord ){
  case 0:
    coefmon[0] = fx[0];
    return;
  case 1:
    coefmon[0] = 0.5 * ( fx[0] + fx[1] );
    coefmon[1] = ( fx[1] - fx[0] ) / ( x[1] - x[0] );
    return;
  default:
    for( unsigned i(0); i<=maxord; i++ ){
      double mulconst2( std::cos(mulconst*double(i)) ),
             mulconst3( 4*std::pow(mulconst2,2)-2 ), 
             b0( 0 ), b1( 0 );
      b0 = fx[maxord];
      b1 = fx[maxord-1] + mulconst3*b0;
      for( unsigned j=maxord-2; j>1; j-=2 ){
        b0 = fx[j] + mulconst3*b1 - b0;
        b1 = fx[j-1] + mulconst3*b0 - b1;
      }
      if( !(maxord%2) )
        b0 = fx[0] + mulconst3*b1 - b0 - b1;
      else{
        b0 = fx[1] + mulconst3*b1 - b0;
        b0 = fx[0] + mulconst3*b0 - b1 - b0;
      }
      coefmon[i] = 2./double(maxord+1)*mulconst2*b0;
    }
    coefmon[0] *=0.5;
    return;
  }
}

template <typename T>
inline void
SCModel<T>::_interpolation
( std::vector<double>& coefmon, double const TOL, unsigned& nord,
  T const& X, puniv f )
{
  coefmon.resize( nord+1 );
  _interpolation( coefmon, nord, X, f );
  for( ; std::fabs(coefmon[nord])>TOL || (nord && std::fabs(coefmon[nord-1])>TOL); ){
    nord*=2;
    coefmon.resize( nord+1 );
    _interpolation( coefmon, nord, X, f );
#ifdef MC__CVAR_DEBUG_TANH
    for( unsigned i=0; i<=nord; i++ )
      std::cout << "a[" << i << "] = " << coefmon[i] << std::endl;
    { int dum; std::cout << "ENTER <1> TO CONTINUE"; std::cin >> dum; }
#endif
  }
}

template <typename T>
template <typename U>
inline SCVar<T>
SCModel<T>::_composition
( std::vector<U> const& coefouter, unsigned const maxord, SCVar<T> const& CVinner )
{
  //composition based on http://en.wikipedia.org/wiki/Clenshaw_algorithm#Special_case_for_Chebyshev_series
  if( !maxord )
    return coefouter[0];

  else if( maxord == 1 )
    return CVinner * coefouter[1] + coefouter[0];

  SCVar<T> CVinnerx2 = 2. * CVinner;
#ifdef MC__SCMODEL_DEBUG_COMPOSITION
  std::cout << "CVinner:" << CVinner;
  std::cout << "CVinnerx2:" << CVinnerx2;
#endif
  SCVar<T> CV1 = coefouter[maxord];
#ifdef MC__SCMODEL_DEBUG_COMPOSITION
  std::cout << "CV1:" << CV1;
#endif
  SCVar<T> CV2 = coefouter[maxord-1] + CVinnerx2 * CV1;
#ifdef MC__SCMODEL_DEBUG_COMPOSITION
  std::cout << "CV2:" << CV2;
#endif

  for( unsigned i=maxord-2; i>1; i-=2 ){
#ifdef MC__SCMODEL_DEBUG_COMPOSITION
    std::cout << "CVinnerx2 * CV2" << CVinnerx2 * CV2;
    std::cout << "coefouter[i]:" << coefouter[i];
    std::cout << "coefouter[i] + CVinnerx2 * CV2" << coefouter[i] + CVinnerx2 * CV2;
#endif
    CV1 = coefouter[i]   + CVinnerx2 * CV2 - CV1;
#ifdef MC__SCMODEL_DEBUG_COMPOSITION
    std::cout << "CV1:" << CV1;
#endif
    CV2 = coefouter[i-1] + CVinnerx2 * CV1 - CV2;
#ifdef MC__SCMODEL_DEBUG_COMPOSITION
    std::cout << "CV2:" << CV2;
#endif
  }

  if( !(maxord%2) ){
#ifdef MC__SCMODEL_DEBUG_COMPOSITION
    CV1 = coefouter[0] + CVinner * CV2 - CV1;
    std::cout << "CV1:" << CV1;
    return CV1;
#else
    return coefouter[0] + CVinner * CV2 - CV1;
#endif
  }

  CV1 = coefouter[1] + CVinnerx2 * CV2 - CV1;
#ifdef MC__SCMODEL_DEBUG_COMPOSITION
  CV2 = coefouter[0] + CVinner * CV1 - CV2;
  std::cout << "CV2:" << CV2;
  return CV2;
#else
    return coefouter[0] + CVinner * CV1 - CV2;
#endif

//if( !(maxord%2) ){
//  CV1 = coefouter[0] + CVinner * CV2 - CV1;
//SCVar<T> tmp = coefouter[1] + CVinnerx2 * CV2 - CV1;
//std::cout << "coefouter[1] + CVinnerx2 * CV2 - CV1:" << tmp;
//for( unsigned i=0; i<tmp.nmon(); i++ ) std::cout << tmp._coefmon[i] << std::endl;
//  CV1 = coefouter[1] + CVinnerx2 * CV2 - CV1;
//for( unsigned i=0; i<CV1.nmon(); i++ ) std::cout << CV1._coefmon[i] << std::endl;
//std::cout << "CV1:" << CV1;
//std::cout << "CVinner * CV1:" << CVinner * CV1;
//  return coefouter[0] + CVinner * CV1 - CV2;
}

template <typename T>
inline SCVar<T>
SCModel<T>::_intpow
( SCVar<T> const& CV, int const n ) const
{
  if( n == 0 ) return 1.;
  else if( n == 1 ) return CV;
  return n%2 ? sqr( _intpow( CV, n/2 ) ) * CV : sqr( _intpow( CV, n/2 ) );
}

template <typename T>
inline void
SCModel<T>::_sscal1D
( t_coefmon const& coefmon0, double const& dscal, t_coefmon& coefmon )
const
{
  if( isequal(dscal,0.) ) return;
  coefmon = coefmon0;
  if( isequal(dscal,1.) ) return;
  for( auto& [mon,coef] : coefmon )
    coef *= dscal;
  return;
}

template <typename T>
inline void
SCModel<T>::_slift1D
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
      itmon->second += (isequal(dscal,1.)? coef0: coef0*dscal);
    }
  }
}

template <typename T>
inline void
SCModel<T>::_slift1D
( t_coefmon const& coefmon0, double const& dscal, t_coefmon& coefmon,
  double& coefrem, typename t_var::const_iterator itvar, unsigned const ndxord )
const
{
  for( auto const& [mon0,coef0] : coefmon0 ){
    SPolyMon mon = mon0; // local copy for modification
    mon.tord += ndxord;
    if( mon.tord > _maxord ){ // append to remainder coefficient if total order too large
      coefrem += std::fabs( isequal(dscal,1.)? coef0: coef0*dscal );
      continue;
    }
#ifdef MC__SPOLYEXPR_DEBUG_SPROD
    std::cout << "mon: " << mon.display(1) << "  (" << *itvar << "," << ndxord << ")" << std::endl;
#endif
    assert( mon.expr.insert( std::make_pair( *itvar, ndxord ) ).second );
    auto [itmon,ins] = coefmon.insert( std::make_pair( mon, isequal(dscal,1.)? coef0: coef0*dscal ) );
    if( !ins ) itmon->second += isequal(dscal,1.)? coef0: coef0*dscal;
  }
}

template <typename T>
inline void
SCModel<T>::_slift1D
( t_coefmon const& coefmon0, double const& dscal, double& coefrem )
const
{
  for( auto const& [mon0,coef0] : coefmon0 )
    coefrem += isequal(dscal,1.)? std::fabs(coef0): std::fabs(coef0)*dscal;
}

template <typename T>
inline void
SCModel<T>::_sdisp1D
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
SCModel<T>::_sdisp1D
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
SCModel<T>::_sdisp1D
( std::vector<SCVar<T>> const& vec, unsigned const ivar, 
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
SCModel<T>::_sprod1D
( std::map<unsigned,t_coefmon> const& sp1map,
  std::map<unsigned,t_coefmon> const& sp2map,
  t_coefmon& coefmon, double& coefrem, std::set<unsigned> const& ndxvar,
  typename t_var::const_iterator itvar )
const
{
  // construct product matrix of polynomial coefficients
  auto itvarnext = itvar;
  if( !ndxvar.empty() ) std::advance( itvarnext, 1 );//++itvarnext;
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
      _sprod1D( sp11map, sp22map, it12->second, coefrem, ndxvar, itvarnext );
    }
  }
  
  // construct 1D product result and augment remainder as appropriate
  coefmon.clear();
  for( auto const& [ndx12,coefmon12] : sp12map ){
    auto const& [ndx1,ndx2] = ndx12;
    // Product involving two constant terms
    if( !ndx1 && !ndx2 )
      _slift1D( coefmon12, 1., coefmon );
    // Product involving exactly one constant term
    else if( !ndx1 || !ndx2 )
      _slift1D( coefmon12, 1., coefmon, coefrem, itvar, ndx1+ndx2 );
    // Product between non-constant Chebyshev basis functions
    else{
      if( ndx1+ndx2 <= _maxord )
        _slift1D( coefmon12, .5, coefmon, coefrem, itvar, ndx1+ndx2 );
      else
        _slift1D( coefmon12, .5, coefrem );
      if( ndx1 == ndx2 )
        _slift1D( coefmon12, .5, coefmon );
      else if( ndx1 >= ndx2 )
        _slift1D( coefmon12, .5, coefmon, coefrem, itvar, ndx1-ndx2 );
      else
        _slift1D( coefmon12, .5, coefmon, coefrem, itvar, ndx2-ndx1 );
    }
  }
#ifdef MC__POLYMODEL_DEBUG_SPROD
  _sdisp1D( coefmon, "Prod: " );
  std::endl;
#endif
}

template <typename T>
inline void
SCModel<T>::_svec1D
( typename t_var::const_iterator itvar, std::pair<SPolyMon,double> const& mon,
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
SCModel<T>::_svec1Dfull
( unsigned const ivar, std::pair<SPolyMon,double> const& coefmon,
  std::vector<SCVar<T>>& vec )
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

// ==> account for sparse format

template <typename T>
template <typename C, typename U>
inline
U
SCModel<T>::_polybound_eigen
( std::map<SPolyMon,C,lt_SPolyMon> const& coefmon, U const* const* bndbasis )
const
{
  return _polybound_naive( coefmon, bndbasis );
/*
  if( _maxord < 2 ) return _polybound_naive( coefmon, bndbasis );

  static double const TOL = 1e-8;

  U bndpol = (ndxmon.empty() || ndxmon.find(0)!=ndxmon.end())? coefmon[0]: 0.;
  //U bndpol = coefmon[0];
  if( _maxord == 1 ) bndpol += bndord[1];

  else if( _maxord > 1 ){
    double*Umat = new double[_nvar*_nvar];
    for( unsigned i=_posord[2]; i<_posord[3]; i++ ){
      unsigned i1=0, i2=_nvar;
      C Ci = (ndxmon.empty() || ndxmon.find(i)!=ndxmon.end())? coefmon[i]: 0.;
      const unsigned*iexp=_expmon+i*_nvar;
      for( ; i1<_nvar; i1++ ) if( iexp[i1] ) break;
      if( iexp[i1] == 2 ){
        Umat[_nvar*i1+i1] = 2.*Ci;
        bndpol -= Ci;
        continue;
      }
      for( i2=i1+1; i2<_nvar; i2++ ) if( iexp[i2] ) break;
      Umat[_nvar*i1+i2] = 0.;
      Umat[_nvar*i2+i1] = Ci/2.;
    }
#ifdef MC__POLYMODEL_DEBUG_POLYBOUND
    display( _nvar, _nvar, Umat, _nvar, "Matrix U", std::cout );
#endif
    double*Dmat = mc::dsyev_wrapper( _nvar, Umat, true );
    if( !Dmat ){
      delete[] Umat;
      return _polybound_LSB( coefmon, bndord, bndbasis );
    }

    for( unsigned i=0; i<_nvar; i++ ){
      double linaux = 0.;
      U bndaux(0.);
      for( unsigned k=0; k<_nvar; k++ ){
        C Ck = (ndxmon.empty() || ndxmon.find(_nvar-k)!=ndxmon.end())? coefmon[_nvar-k]: 0.;
        linaux += Umat[i*_nvar+k] * Ck;
        bndaux += Umat[i*_nvar+k] * bndbasis[k][1];
      }
#ifdef MC__POLYMODEL_DEBUG_POLYBOUND
      std::cout << i << ": LINAUX = " << linaux
                << "  BNDAUX = " << bndaux << std::endl;
#endif
      if( std::fabs(Dmat[i]) > TOL )
        bndpol += Dmat[i] * Op<U>::sqr( linaux/Dmat[i]/2. + bndaux )
                - linaux*linaux/Dmat[i]/4.;
      else     
        bndpol += linaux * bndaux + Dmat[i] * Op<U>::sqr( bndaux );
#ifdef MC__POLYMODEL_DEBUG_POLYBOUND
        std::cout << "BNDPOL: " << bndpol << std::endl;
#endif
    }
    delete[] Umat;
    delete[] Dmat;
  }
#ifdef MC__POLYMODEL_DEBUG_POLYBOUND
  int tmp; std::cin >> tmp;
#endif

  for( unsigned i=3; i<=_maxord; i++ ) bndpol += bndord[i];
  return bndpol;
*/
}

// ==> account for sparse format

template <typename T>
template <typename C, typename U>
inline
U
SCModel<T>::_polybound_LSB
( std::map<SPolyMon,C,lt_SPolyMon> const& coefmon, U const* const* bndbasis )
const
{
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
      bndpol += it2->second * ( bndbasis? bndbasis[ie2->first][1] * bndbasis[(++ie2)->first][1]: _TOne );
      continue;
    }
    SPolyMon explin( 1, std::map<unsigned,unsigned>() );
    explin.expr.insert( std::make_pair( ie2->first, 1 ) ); 
    it1 = coeflin.find( explin );
    if( it1 != coeflin.end() && std::fabs(it2->second) > TOL ){
       bndpol += (2.*it2->second) * Op<U>::sqr( (bndbasis? bndbasis[ie2->first][1]: _TOne)
                 + it1->second/(it2->second*4.) ) - it2->second - it1->second*it1->second/8./it2->second;
       coeflin.erase( it1 );
    }
    else if( it1 != coeflin.end() ){
      bndpol += it2->second * ( bndbasis? bndbasis[ie2->first][2]: _TOne )
              + it1->second * ( bndbasis? bndbasis[ie2->first][1]: _TOne );
       coeflin.erase( it1 );
    }
    else
      bndpol += it2->second * ( bndbasis? bndbasis[ie2->first][2]: _TOne );
  }

  // Remaining linear terms
  for( it1=coeflin.begin(); it1!=coeflin.end(); ++it1 ){
    auto ie1 = it1->first.expr.begin();
    bndpol += it1->second * ( bndbasis? bndbasis[ie1->first][1]: _TOne );
  }

  // Thrid and higher-order terms
  if( coefmon.rbegin()->first.tord > 2 )
    bndpol += _polybound_naive( coefmon, bndbasis, 3 );

  return bndpol;
}
/*
template <typename T> template <typename C, typename U> inline U
SCModel<T>::_polybound_bernstein
( const C*coefmon, const U*bndord, const U*const*bndbasis,
  const std::set<unsigned>&ndxmon )
{
  // Expand binomial coefficient and exponent arrays if needed
#ifdef MC__SCVAR_DEBUG_BERNSTEIN
  std::cout << "binom max: " << _binom_size.first << "  "
            << _binom_size.second << std::endl;
#endif
  const unsigned maxord = (options.BOUNDER_ORDER>_maxord? 
    options.BOUNDER_ORDER: _maxord );
  const poly_size maxmon = std::pow(maxord+1,_nvar);
  _ext_expmon( maxord, true );
  _ext_binom( 2*_maxord );
#ifdef MC__SCVAR_DEBUG_BERNSTEIN
  std::cout << "binom max: " << _binom_size.first << "  "
            << _binom_size.second << std::endl;
#endif

  // Compute transformation matrix
  double*trmat = new double[(maxord+1)*(_maxord+1)];
#ifdef MC__SCVAR_DEBUG_BERNSTEIN
  std::cout << "trmat:\n" << std::scientific << std::setprecision(5);
#endif
  for( unsigned j=0, jk=0; j<=maxord; j++ ){
    for( unsigned k=0; k<=_maxord; k++, jk++ ){
      trmat[jk] = _transform_bernstein( j, k, maxord );
#ifdef MC__SCVAR_DEBUG_BERNSTEIN
      std::cout << "  " << trmat[jk];
#endif
    }
#ifdef MC__SCVAR_DEBUG_BERNSTEIN
    std::cout << std::endl;
#endif
  }

  // Compute min/max amongst all Bernstein coefficients
  U bndpol = (ndxmon.empty() || ndxmon.find(0)!=ndxmon.end())? coefmon[0]: 0.;
  //U bndpol = coefmon[0];
#ifdef MC__SCVAR_DEBUG_BERNSTEIN
  std::cout << "\n0:  " << bndpol << std::endl;
#endif
  for( poly_size jmon=0; jmon<maxmon; jmon++ ){ // Loop over Bernstein terms
    const unsigned*jexp = _expmon + jmon*_nvar;
    C coefbern = _coef_bernstein( coefmon, jexp, trmat, ndxmon );
    bndpol = Op<U>::hull( bndpol, coefbern );
#ifdef MC__SCVAR_DEBUG_BERNSTEIN
    std::cout << jmon << " ["; 
    for( unsigned ivar=0; ivar<_nvar; ivar++ )
      std::cout << std::setw(3) << jexp[ivar];
    std::cout << "] : " << coefbern << " : " << bndpol << std::endl;
#endif
  }

  delete[] trmat;
  return bndpol;
}

template <typename T> inline double
SCModel<T>::_transform_bernstein
( const unsigned j, const unsigned k, const unsigned n ) const
{
  double Mjk = 0.;
  for( unsigned i=(j+k<=n?0:j+k-n); i<=(j<=k?j:k); i++ ){
    double Mjki = (double)_get_binom(2*k,2*i) * (double)_get_binom(n-k,j-i);
    (k+i)%2? Mjk -= Mjki: Mjk += Mjki;
    //std::cout << k << "  " << i << " (k-i)%2: " << (k-i)%2 << std::endl;
  }
  Mjk /= (double)_get_binom(n,j);
  return Mjk;
}

template <typename T> template <typename C> inline C
SCModel<T>::_coef_bernstein
( const std::map< t_expmon, C >& coefmon, const unsigned*jexp,
  double const*trmat ) const
{
  // Compute bernstein coefficient with variables indices <tt>jexp</tt>
  C coefbern = (ndxmon.empty() || ndxmon.find(0)!=ndxmon.end())? coefmon[0]: 0.;
#ifdef MC__SCVAR_DEBUG_BERNSTEIN
    std::cout << "  [  0]  " << coefmon[0] << "  " << coefbern << std::endl;
#endif
  // Sparse loop over Chebyshev terms
  for( auto it=ndxmon.begin(); it!=ndxmon.end(); ++it ){
    const unsigned*kexp = _expmon + (*it)*_nvar;
    C termbern = coefmon[(*it)]; // Chebyshev coefficient for (*it)
#ifdef MC__SCVAR_DEBUG_BERNSTEIN
    std::cout << "  ["; 
    for( unsigned ivar=0; ivar<_nvar; ivar++ )
      std::cout << std::setw(3) << kexp[ivar];
    std::cout << "]";
#endif
    for( unsigned ivar=0; ivar<_nvar; ivar++ )
      termbern *= trmat[jexp[ivar]*(_maxord+1)+kexp[ivar]];
    coefbern += termbern;      
#ifdef MC__SCVAR_DEBUG_BERNSTEIN
    std::cout << "  " << termbern << "  " << coefbern << std::endl;
#endif
  }  
  // Dense loop over Chebyshev terms
  for( unsigned kmon=1; ndxmon.empty() && kmon<_nmon; kmon++ ){
    const unsigned*kexp = _expmon + kmon*_nvar;
    C termbern = coefmon[kmon]; // Chebyshev coefficient for kmon
#ifdef MC__SCVAR_DEBUG_BERNSTEIN
    std::cout << "  ["; 
    for( unsigned ivar=0; ivar<_nvar; ivar++ )
      std::cout << std::setw(3) << kexp[ivar];
    std::cout << "]";
#endif
    for( unsigned ivar=0; ivar<_nvar; ivar++ )
      termbern *= trmat[jexp[ivar]*(_maxord+1)+kexp[ivar]];
    coefbern += termbern;      
#ifdef MC__SCVAR_DEBUG_BERNSTEIN
    std::cout << "  " << termbern << "  " << coefbern << std::endl;
#endif
  }
#ifdef MC__SCVAR_DEBUG_BERNSTEIN2
  std::cout << std::endl;
#endif
  return coefbern;
}
*/
template <typename T>
template <typename C, typename U>
inline
U
SCModel<T>::_polybound_naive
( std::map<SPolyMon,C,lt_SPolyMon> const& coefmon, U const* const* bndbasis,
  unsigned const minord )
const
{
  // Empty model
  if( coefmon.empty() || coefmon.rbegin()->first.tord < minord ) return 0.;
  auto it = coefmon.lower_bound( SPolyMon( minord, std::map<unsigned,unsigned>() ) );
  // Polynomial bounding in T arithmetic
  if( !bndbasis ){
    double bndcst = 0., bndcoef = 0.;
    if( !it->first.tord ){ bndcst = it->second; ++it; }
    for( ; it!=coefmon.end(); ++it )
      bndcoef += std::fabs( it->second );
    return bndcoef * _TOne + bndcst;
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
SCModel<T>::_polybound
( std::map<SPolyMon,C,lt_SPolyMon> const& coefmon, U const* const* bndbasis,
  int const type )
{
  switch( type ){
  //case Options::BERNSTEIN:
  // return _polybound_bernstein( coefmon, bndbasis );
  case Options::EIGEN:
    return _polybound_eigen( coefmon, bndbasis );
  case Options::LSB:
    return _polybound_LSB( coefmon, bndbasis );
  case Options::NAIVE: default:
    return _polybound_naive( coefmon, bndbasis );
  }
}

////////////////////////////////// SCVar ///////////////////////////////////////

template <typename T>
inline
SCVar<T>&
SCVar<T>::operator=
( SCVar<T> const& var )
{
  _CM = var._CM;
  _set( var ); // defined in spolymodel.hpp
  return *this;
}

template <typename T>
inline
SCVar<T>::SCVar
( SCModel<T>*CM )
: SPolyVar<T>(), _CM( CM )
{
  _bndrem = 0.;
  _set_bndpol( 0. );
  if( _CM->options.MIXED_IA ) _set_bndT( 0. );
}

template <typename T>
inline
SCVar<T>::SCVar
( double const d, SCModel<T>*CM )
: SPolyVar<T>(), _CM( CM )
{
  if( isequal( d, 0. ) ) return;
  _coefmon.insert( std::make_pair( SPolyMon(), d ) );
  _bndrem = 0.;
  _set_bndpol( d );
  if( !_CM || _CM->options.MIXED_IA ) _set_bndT( d );
}

template <typename T>
inline
SCVar<T>&
SCVar<T>::operator=
( double const d )
{
  _reinit();
  _CM = 0;
  _coefmon.insert( std::make_pair( SPolyMon(), d ) );
  _bndrem = 0.;
  _set_bndpol( d );
  _set_bndT( d );
  return *this;
}

template <typename T>
inline
SCVar<T>::SCVar
( T const& B, SCModel<T>*CM )
: SPolyVar<T>(), _CM( CM )
{
  double const midB = Op<T>::mid(B);
  _coefmon.insert( std::make_pair( SPolyMon(), midB ) );
  _bndrem = B - midB;
  _set_bndpol( midB );
  if( !_CM || _CM->options.MIXED_IA ) _set_bndT( B );
}

template <typename T>
inline
SCVar<T>&
SCVar<T>::operator=
( T const& B )
{
  _reinit();
  _CM = 0;
  double const midB = Op<T>::mid(B);
  _coefmon.insert( std::make_pair( SPolyMon(), midB ) );
  _bndrem = B - midB;
  _set_bndpol( midB );
  _set_bndT( B );
  return *this;
}

template <typename T>
inline
SCVar<T>::SCVar
( SCModel<T>* CM, unsigned const ivar, T const& X )
: _CM( CM )
{
  _set( ivar, X );
}

template <typename T>
inline
SCVar<T>&
SCVar<T>::_set
( unsigned const i, T const& X, bool const updMod )
{
  if( !_CM ) throw typename SCModel<T>::Exceptions( SCModel<T>::Exceptions::INIT );  

  // Keep data for variable #ivar in model environment
  if( updMod ) _CM->_set( i, X );

  // Populate model variable
  _ndxvar.clear();
  _ndxvar.insert(i);
  _coefmon.clear();
  _coefmon.insert( std::make_pair( SPolyMon(),_refvar(i) ) );
  if( _CM->_maxord && !isequal(_scalvar(i),0.) ){
    _coefmon.insert( std::make_pair( SPolyMon(i), _scalvar(i) ) );
    _set_bndpol( _bndvar(i) );
    _bndrem = 0.;
  }
  else{
    _set_bndpol( _refvar(i) );
    _bndrem = _bndvar(i) - _refvar(i);
  }
  //std::cout << "coefmon size: " << _coefmon.size() << std::endl;

  // Interval bounds
  if( _CM->options.MIXED_IA ) _set_bndT( _bndvar(i) );
  else                        _unset_bndT();
  if( _CM->options.MIN_FACTOR >= 0. ) simplify( _CM->options.MIN_FACTOR );
  return *this;
}

template <typename T>
inline
double
SCVar<T>::polynomial
( double const* x ) const
{
  // -> Is there a multivariate version of Clenshaw? Use recursively?

  double Pval = 0.;
  for( auto& [mon,coef] : _coefmon ){
    double val = coef;
    for( auto& [ivar,iord] : mon.expr )
      val *= isequal(_scalvar(ivar),0.)? _refvar(ivar): 
             mc::cheb((x[ivar]-_refvar(ivar))/_scalvar(ivar),iord);
    Pval += val;
  }
  return Pval;
}

template <typename T>
inline
double
SCVar<T>::constant( const bool reset )
{
  auto it_0 = ( _coefmon.empty() || _coefmon.begin()->first.tord? _coefmon.end(): _coefmon.begin() );
  double const coefcst = ( it_0 == _coefmon.end()? 0.: it_0->second );
  if( reset && it_0 != _coefmon.end() ){
    _coefmon.erase( it_0 );
    if( _bndpol ) *_bndpol -= coefcst;
    if( _bndT )   *_bndT -= coefcst;
  }
  return coefcst;
}

template <typename T>
inline
double
SCVar<T>::linear
( const unsigned i, const bool reset )
{
  if( i>=nvar() || !nord() ) return 0.;
  auto it_i = ( _coefmon.empty()? _coefmon.end(): _coefmon.find( SPolyMon( i ) ) );
  double const coeflin = ( it_i == _coefmon.end() || isequal(_scalvar(i),0.)? 0.: it_i->second/_scalvar(i) );
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
SCVar<T>::_scale
( unsigned const ivar, T const& Xivar, typename SCVar<T>::t_coefmon& coefmon )
const
{
  // Nothing to do if model _CM is NULL, i is outside of variable range,
  // or variable range X did not change
  if( !_CM || _ndxvar.find( ivar ) == _ndxvar.end()
   || ( isequal( Op<T>::l(Xivar), Op<T>::l(_bndvar(ivar)) )
     && isequal( Op<T>::u(Xivar), Op<T>::u(_bndvar(ivar)) ) ) ) return;

  // Get coefficients in univariate polynomial representation w.r.t variable #ivar
  std::vector<SCVar<T>> veccoef( nord()+1, SCVar<T>(_CM) );
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
  SCVar<T> CVivar( _CM ); CVivar._set( ivar, Xivar, false );
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
SCVar<T>&
SCVar<T>::scale
( unsigned const i, T const& X )
{
  _scale( i, X, _coefmon );
  _unset_bndpol();
  _unset_bndT();
  return *this;
}

template <typename T>
inline
SCVar<T>&
SCVar<T>::scale
( T const* X )
{
  // Return *this if null pointer to model _CM or variable ranges X
  if( !X || !_CM ) return *this;
  for( unsigned i : _ndxvar ) scale( i, X[i] );
  if( _CM && _CM->options.MIN_FACTOR >= 0. ) simplify( _CM->options.MIN_FACTOR );
  return *this;
}

template <typename T>
inline
typename SCVar<T>::t_coefmon
SCVar<T>::unscale
()
const
{
  // Return *this if null pointer to model _CM or variable ranges X
  if( !_CM ) return _coefmon;
  t_coefmon coefmon = _coefmon;
  for( unsigned i : _ndxvar ) _scale( i, T(-1,1), coefmon );
  if( _CM && _CM->options.MIN_FACTOR >= 0. ) simplify( _CM->options.MIN_FACTOR );
  return coefmon;
}

template <typename T>
inline
SCVar<T>&
SCVar<T>::simplify
( double const TOL, int const TORD )
{
  if( _coefmon.empty() ) return *this;
  for( auto it=_coefmon.begin(); it!=_coefmon.end(); ){
    // Eliminate any terms with zero coefficient
    if( it->second == 0. ){
      it = _coefmon.erase( it );
      continue;
    }
    // Eliminate any non-constant terms with small enough coefficient or large enough total order
    if( ( TOL > 0e0 && it->first.tord && std::fabs(it->second) <= TOL )
     || ( TORD >= 0 && (int)it->first.tord > TORD ) ){
      _bndrem += it->second * _TOne;
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
void
SCVar<T>::_simplify
( typename SCVar<T>::t_coefmon& coefmon )
const
{
  if( coefmon.empty() ) return;
  for( auto it=coefmon.begin(); it!=coefmon.end(); ){
    // Eliminate any terms with zero coefficient
    if( it->second == 0. ){
      it = coefmon.erase( it );
      continue;
    }
    ++it; // only increment if current element was not erased
  }
}

template <typename T>
inline
void
SCVar<T>::_to_monomial
( unsigned const ivar, typename SCVar<T>::t_coefmon& coefmon )
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
typename SCVar<T>::t_coefmon
SCVar<T>::to_monomial
( bool const scaled )
const
{
  if( !_CM || _coefmon.empty() || !nord() ) return _coefmon;
  //std::cout << *this;

  // Convert to monomial form
  t_coefmon coefmon = _coefmon;
  if( !scaled ){
    for( unsigned i : _ndxvar ){
      _scale( i, T(-1,1), coefmon );
      //_simplify( coefmon );
      std::cout << display( coefmon, 1 );
    }
  }
  
  for( unsigned i : _ndxvar ){
//    if( !scaled ){
//      _scale( i, T(-1,1), coefmon );
//      _simplify( coefmon );
//    }
    _to_monomial( i, coefmon );
    _simplify( coefmon );
  }
  
  return coefmon;
}

template <typename T>
inline
std::ostream&
operator<<
( std::ostream& out, SCVar<T> const& CV )
{
  //std::cout << "*NORD = " << CV.nord() << " " << CV._coefmon.rbegin()->first.display(1) << std::endl;
  unsigned IDISP = CV._CM? CV._CM->options.DISPLAY_DIGITS: 5;
  out << std::endl << std::scientific << std::setprecision(IDISP)
      << std::right;

  // Sparse multivariate polynomial
  for( auto const& [mon,coef] : CV._coefmon ){
    out << std::right << std::setw(IDISP+7) << coef << "  " << std::setw(2) << mon.tord << "  ";
    for( auto& [ivar,iord] : mon.expr )
      out << "T" << iord << "[" << ivar << "] ";
    out << std::endl;
  }

  // Remainder term
  out << std::right << "   R     =  " << CV._bndrem
      << std::endl;

  // Range bounder
  out << std::right << "   B     =  " << CV.B()
      << std::endl;

  // Index set
  out << std::right << "   I     =  {";
  for( auto const& ndx : CV._ndxvar )
    out << std::right << " " << ndx;
  out << std::right << " }" << std::endl;

  return out;
}

template <typename T>
inline
SCVar<T>
operator+
( SCVar<T> const& CV )
{
  return CV;
}

template <typename T>
template <typename U>
inline
SCVar<T>&
SCVar<T>::operator+=
( SCVar<U> const& CV )
{
//  if( CV._CM && !_CM ){
//    SCVar<T> CV2( *this );
//    *this = CV;
//    return *this += CV2;
//  }

  if( _CM && CV._CM && _CM != CV._CM )
    throw typename SCModel<T>::Exceptions( SCModel<T>::Exceptions::SCMODEL );
  if( CV._CM && !_CM ) _CM = CV._CM;

  _ndxvar.insert( CV._ndxvar.begin(), CV._ndxvar.end() );
  for( auto&& [mon,coef] : CV._coefmon ){
    //std::cout << mon << ": " << coef << std::endl;
    // No warm-start with map::insert unfortunately...
    auto [itmon, ins] = _coefmon.insert( std::make_pair( mon, coef ) );
    //if( !ins) std::cout << itmon->first << ": " << itmon->second << std::endl;
    if( !ins ) itmon->second += coef;
  }

  _bndrem += CV._bndrem;
  _unset_bndpol();
  if( _bndT && CV._bndT ) *_bndT += *CV._bndT;
  else _unset_bndT();
  if( _CM && _CM->options.MIN_FACTOR >= 0. ) simplify( _CM->options.MIN_FACTOR );

  return *this;
}

template <typename T, typename U>
inline SCVar<T>
operator+
( SCVar<T> const& CV1, SCVar<U> const& CV2 )
{
  if( CV1.nmon() >= CV2.nmon() ){
    SCVar<T> CV3( CV1 );
    CV3 += CV2;
    return CV3;
  }
  
  SCVar<T> CV3( CV2 );
  CV3 += CV1;
  return CV3;
}

template <typename T>
inline
SCVar<T>&
SCVar<T>::operator +=
( double const c )
{
  if( isequal( c, 0. ) ) return *this;
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
SCVar<T>
operator+
( SCVar<T> const& CV1, double const c )
{
  SCVar<T> CV3( CV1 );
  CV3 += c;
  return CV3;
}

template <typename T>
inline
SCVar<T>
operator+
( double const c, SCVar<T> const& CV2 )
{
  SCVar<T> CV3( CV2 );
  CV3 += c;
  return CV3;
}

template <typename T>
template <typename U>
inline
SCVar<T>&
SCVar<T>::operator+=
( U const& I )
{
  _bndrem += I;
  _center();
  if( _bndT ) *_bndT += I;
  return *this;
}

template <typename T>
inline
SCVar<T>
operator+
( SCVar<T> const& CV1, T const& I )
{
  SCVar<T> CV3( CV1 );
  CV3 += I;
  return CV3;
}

template <typename T>
inline
SCVar<T>
operator+
( T const& I, SCVar<T> const& CV2 )
{
  SCVar<T> CV3( CV2 );
  CV3 += I;
  return CV3;
}

template <typename T>
inline
SCVar<T>
operator-
( SCVar<T> const& CV )
{
  SCVar<T> CV2;
  CV2.set( CV._CM );
  CV2._ndxvar = CV._ndxvar;
  for( auto& [mon,coef] : CV._coefmon )
    CV2._coefmon.insert( std::make_pair( mon, -coef ) );
  CV2._bndrem = - CV._bndrem;
  if( CV._bndpol ) CV2._set_bndpol( - *CV._bndpol );
  if( CV._bndT )   CV2._set_bndT( - *CV._bndT );
  return CV2;
}

template <typename T>
template <typename U>
inline
SCVar<T>&
SCVar<T>::operator-=
( SCVar<U> const& CV )
{
//  if( CV._CM && !_CM ){
//    SCVar<T> CV2( *this );
//    *this = -CV;
//    return *this += CV2;
//  }

  if( _CM && CV._CM && _CM != CV._CM )
    throw typename SCModel<T>::Exceptions( SCModel<T>::Exceptions::SCMODEL );
  if( CV._CM && !_CM ) _CM = CV._CM;

  _ndxvar.insert( CV._ndxvar.begin(), CV._ndxvar.end() );
  for( auto& [mon,coef] : CV._coefmon ){
    // No warm-start with map::insert unfortunately...
    auto [itmon, ins] = _coefmon.insert( std::make_pair( mon, -coef ) );
    if( !ins ) itmon->second -= coef;
  }

  _bndrem -= CV._bndrem;
  _unset_bndpol();
  if( _bndT && CV._bndT ) *_bndT -= *CV._bndT;
  else _unset_bndT();
  if( _CM && _CM->options.MIN_FACTOR >= 0. ) simplify( _CM->options.MIN_FACTOR );

  return *this;
}

template <typename T, typename U>
inline
SCVar<T>
operator-
( SCVar<T> const& CV1, SCVar<U> const& CV2 )
{
  if( CV1.nmon() >= CV2.nmon() ){
    SCVar<T> CV3( CV1 );
    CV3 -= CV2;
    return CV3;
  }

  SCVar<T> CV3( -CV2 );
  CV3 += CV1;
  return CV3;
}

template <typename T>
inline
SCVar<T>&
SCVar<T>::operator-=
( double const c )
{
  if( isequal( c, 0. ) ) return *this;
  if( _coefmon.empty() || _coefmon.begin()->first.tord ) 
    _coefmon.insert( std::make_pair( SPolyMon(), -c ) );
  else
    _coefmon.begin()->second -= c;
  if( _bndpol ) *_bndpol -= c;
  if( _bndT )   *_bndT -= c;
  return *this;
}

template <typename T>
inline
SCVar<T>
operator-
( SCVar<T> const& CV1, double const c )
{
  SCVar<T> CV3( CV1 );
  CV3 -= c;
  return CV3;
}

template <typename T>
inline
SCVar<T>
operator-
( double const c, SCVar<T> const& CV2 )
{
  SCVar<T> CV3( -CV2 );
  CV3 += c;
  return CV3;
}

template <typename T>
template <typename U>
inline
SCVar<T>&
SCVar<T>::operator-=
( U const& I )
{
  *_bndrem -= I;
  _center();
  if( _bndT ) *_bndT -= I;
  return *this;
}

template <typename T>
inline SCVar<T>
operator-
( SCVar<T> const& CV1, T const& I )
{
  SCVar<T> CV3( CV1 );
  CV3 -= I;
  return CV3;
}

template <typename T>
inline
SCVar<T>
operator-
( T const& I, SCVar<T> const& CV2 )
{
  SCVar<T> CV3( -CV2 );
  CV3 += I;
  return CV3;
}

template <typename T>
inline
SCVar<T>&
SCVar<T>::operator*=
( SCVar<T> const& CV )
{
//  if( CV._CM && !_CM ){
//    SCVar<T> CV2( *this );
//    *this = -CV;
//    return *this += CV2;
//  }

  if( _CM && CV._CM && _CM != CV._CM )
    throw typename SCModel<T>::Exceptions( SCModel<T>::Exceptions::SCMODEL );
  if( CV._CM && !_CM ) _CM = CV._CM;

  // Remainder propagation
  T R1 = bound() * CV._bndrem + CV._polybound() * _bndrem;
  T R2 = _polybound() * CV._bndrem + CV.bound() * _bndrem;
  double coefrem( 0. );

  // Coefficient maps for first participating variable
  _ndxvar.insert( CV._ndxvar.begin(), CV._ndxvar.end() );
  auto itvar = _ndxvar.begin();
  std::map<unsigned,t_coefmon> sp1map, sp2map;
  for( auto&& mon : _coefmon ) _CM->_svec1D( itvar, mon, sp1map );
#ifdef MC__SCMODEL_DEBUG_SPROD
  _CM->_sdisp1D( sp1map, itvar, "Poly #1: " );
#endif
  for( auto&& mon : CV._coefmon ) _CM->_svec1D( itvar, mon, sp2map );
#ifdef MC__SCMODEL_DEBUG_SPROD
  _CM->_sdisp1D( sp2map, itvar, "Poly #2: " );
#endif

  // Recursive product of univariate Chebyshev polynomials
  _CM->_sprod1D( sp1map, sp2map, _coefmon, coefrem, _ndxvar, itvar );
#ifdef MC__SCMODEL_DEBUG_SPROD
  std::cout << "coefrem = " << coefrem << std::endl;
#endif

  // Remainder propagation
  if( !Op<T>::inter( _bndrem, R1, R2) )
    _bndrem = ( Op<T>::diam(R1) < Op<T>::diam(R2)? R1: R2 );
#ifdef MC__SCMODEL_DEBUG_SPROD
  std::cout << "bndrem = " << _bndrem << std::endl;
#endif
  _bndrem += coefrem * _TOne;
  
  _unset_bndpol();
  if( _bndT && CV._bndT ) *_bndT *= *CV._bndT;
  else _unset_bndT();
  if( _CM && _CM->options.MIN_FACTOR >= 0. ) simplify( _CM->options.MIN_FACTOR );

  return *this;
}

template <typename T>
inline
SCVar<T>
operator*
( SCVar<T> const& CV1, SCVar<T> const& CV2 )
{
  if( &CV1 == &CV2 ) return sqr( CV1 );
  if( CV1.nmon() >= CV2.nmon() ){
    SCVar<T> CV3( CV1 );
    CV3 *= CV2;
    return CV3;
  }

  SCVar<T> CV3( CV2 );
  CV3 *= CV1;
  return CV3;
}

template <typename T>
inline
SCVar<T>&
SCVar<T>::operator*=
( double const c )
{
  if( isequal( c, 0. ) ){ *this = 0.; return *this; }
  if( isequal( c, 1. ) ) return *this;
  for( auto& [mon,coef] : _coefmon ) coef *= c;
  _bndrem *= c;
  if( _bndpol ) *_bndpol *= c;
  if( _bndT ) *_bndT *= c;
  return *this;
}

template <typename T>
inline
SCVar<T>
operator*
( SCVar<T> const& CV1, double const c )
{
  SCVar<T> CV3( CV1 );
  CV3 *= c;
  return CV3;
}

template <typename T>
inline
SCVar<T>
operator*
( double const c, SCVar<T> const& CV2 )
{
  SCVar<T> CV3( CV2 );
  CV3 *= c;
  return CV3;
}

template <typename T>
inline
SCVar<T>&
SCVar<T>::operator*=
( T const& I )
{
  double const Imid = Op<T>::mid(I);
  T Icur = bound();
  for( auto& [mon,coef] : _coefmon ) coef *= Imid;
  _bndrem *= Imid;
  _bndrem += ( I - Imid ) * Icur;
  _unset_bndpol();
  if( _bndT ) *_bndT *= I;
  return *this;
}

template <typename T> inline SCVar<T>
operator*
( SCVar<T> const& CV1, T const& I )
{
  SCVar<T> CV3( CV1 );
  CV3 *= I;
  return CV3;
}

template <typename T>
inline
SCVar<T>
operator*
( T const& I, SCVar<T> const& CV2 )
{
  SCVar<T> CV3( CV2 );
  CV3 *= I;
  return CV3;
}

template <typename T>
inline
SCVar<T>
sqr
( SCVar<T> const& CV )
{
  if( !CV._CM ) return Op<T>::sqr( CV.B() );

  // Coefficient maps for first participating variable
  auto itvar = CV._ndxvar.begin();
  std::map<unsigned,typename SCVar<T>::t_coefmon> sp1map, sp2map;
  for( auto&& mon : CV._coefmon ) CV._CM->_svec1D( itvar, mon, sp1map );
#ifdef MC__SPOLYEXPR_DEBUG_SQR
  CV._CM->_sdisp1D( sp1map, itvar, "Var: " );
#endif

  // Recursive product of univariate Chebyshev polynomials
  SCVar<T> CVSQR;
  CVSQR.set( CV._CM );
  double coefrem( 0. );
  CVSQR._ndxvar = CV._ndxvar;
  CV._CM->_sprod1D( sp1map, sp1map, CVSQR._coefmon, coefrem, CVSQR._ndxvar, itvar );

  // Remainder propagation
  CVSQR._bndrem = coefrem * CV._TOne + ( CV.bound() + CV._polybound() ) * CV._bndrem;
  
  // Bound propagation
  if( CVSQR._CM->options.MIXED_IA ) CVSQR._set_bndT( Op<T>::sqr(CV.bound()) );
  if( CVSQR._CM->options.MIN_FACTOR >= 0. ) CVSQR.simplify( CVSQR._CM->options.MIN_FACTOR );
  return CVSQR;
}

template <typename T>
inline
SCVar<T>&
SCVar<T>::operator/=
( SCVar<T> const& CV )
{
  *this *= inv(CV);
  return *this;
}

template <typename T> inline SCVar<T>
operator/
( SCVar<T> const& CV1, SCVar<T> const& CV2 )
{
  return CV1 * inv(CV2);
}

template <typename T> inline SCVar<T>&
SCVar<T>::operator/=
( double const c )
{
  if( isequal( c, 0. ) )
    throw typename SCModel<T>::Exceptions( SCModel<T>::Exceptions::DIV );
  if( isequal( c, 1. ) ) return *this;
  *this *= (1./c);
  return *this;
}

template <typename T> inline SCVar<T>
operator/
( SCVar<T> const& CV, double const c )
{
  if ( isequal( c, 0. ))
    throw typename SCModel<T>::Exceptions( SCModel<T>::Exceptions::DIV );
  if( isequal( c, 1. ) ) return CV;
  return CV * (1./c);
}

template <typename T> inline SCVar<T>
operator/
( double const c, SCVar<T> const& CV )
{
  if( isequal( c, 0. ) ) return 0.;
  if( isequal( c, 1. ) ) return inv(CV);
  return inv(CV) * c;
}

template <typename T> inline SCVar<T>
inv
( SCVar<T> const& CV )
{
  if( !CV._CM )
    return SCVar<T>( Op<T>::inv( CV.B() ) );
  if ( Op<T>::l(CV.B()) <= 0. && Op<T>::u(CV.B()) >= 0. )
    throw typename SCModel<T>::Exceptions( SCModel<T>::Exceptions::INV );

  SCVar<T> CVI( CV._CM ), CV2( CV._CM );
  auto& coefmon = CV._coefinterp();
  CV._interpolation( coefmon, mc::inv );

  double m(Op<T>::mid(CV.B())), r(Op<T>::u(CV.B())-m), rem;
#ifndef MC__SCVAR_FORCE_REM_DERIV
  double ub(0), lb(0);
  for (unsigned i(0); i<=CV.maxord(); i++) {
    ub += coefmon[i];
    lb += std::pow(-1.,i)*coefmon[i];
  }
  rem = std::max(std::fabs(mc::inv(r+m)-ub), std::fabs(mc::inv(m-r)-lb));
#else
  rem = 4.*std::pow(r,double(CV.maxord()+2));
#endif

  CVI = CV._rescale(r,m);
  CV2 = CVI._composition( coefmon );
  CV2 += CV._TOne * rem;

  CV2._ndxvar = CV._ndxvar;
  if( CV._CM->options.MIXED_IA ) CV2._set_bndT( Op<T>::inv( CV.B() ) );
  if( CV._CM->options.MIN_FACTOR >= 0. ) CV2.simplify( CV._CM->options.MIN_FACTOR );
  return CV2;
}

template <typename T> inline SCVar<T>
sqrt
( SCVar<T> const& CV )
{
  if( !CV._CM )
    return SCVar<T>( Op<T>::sqrt( CV.B() ) );
  if ( Op<T>::l(CV.B()) < 0. )
    throw typename SCModel<T>::Exceptions( SCModel<T>::Exceptions::SQRT );

  SCVar<T> CVI( CV._CM ), CV2( CV._CM );
  auto& coefmon = CV._coefinterp();
  CV._interpolation( coefmon, std::sqrt );

  double b(Op<T>::mid(CV.B())), a(Op<T>::u(CV.B())-b), rem, ub(0), lb(0);
  for (unsigned i(0); i<=CV.maxord(); i++) {
    ub += coefmon[i];
    lb += i%2? -coefmon[i]: coefmon[i];
  }
  rem = std::max(std::fabs(std::sqrt(a+b)-ub), std::fabs(std::sqrt(b-a)-lb));

  CVI = CV._rescale(a,b);
  CV2 = CVI._composition( coefmon );
  CV2 += CV._TOne * rem;

  CV2._ndxvar = CV._ndxvar;
  if( CV._CM->options.MIXED_IA ) CV2._set_bndT( Op<T>::sqrt( CV.B() ) );
  if( CV._CM->options.MIN_FACTOR >= 0. ) CV2.simplify( CV._CM->options.MIN_FACTOR );
  return CV2;
}

template <typename T> inline SCVar<T>
exp
( SCVar<T> const& CV )
{ 
  if( !CV._CM )
    return SCVar<T>( Op<T>::exp( CV.B() ) );

  SCVar<T> CVI( CV._CM ), CV2( CV._CM );
  auto& coefmon = CV._coefinterp();
  CV._interpolation( coefmon, std::exp );

  double m(Op<T>::mid(CV.B())), r(Op<T>::u(CV.B())-m), rem;
#ifndef MC__SCVAR_FORCE_REM_DERIV
  double ub(0), lb(0);
  for (unsigned i(0); i<=CV.maxord(); i++) {
    ub += coefmon[i];
    lb += std::pow(-1.,i)*coefmon[i];
  }
  rem = std::max(std::fabs(std::exp(r+m)-ub), std::fabs(std::exp(m-r)-lb));
#else
  double fact(1);
  for (unsigned i(1); i<=CV.maxord()+1; i++) fact *= double(i);
  double M = Op<T>::abs(Op<T>::exp(CV.B()));
  rem = 2.*M*std::pow(r/2.,double(CV.maxord()+1))/fact;
#endif

  CVI = CV._rescale(r,m);
  CV2 = CVI._composition( coefmon );
  CV2 += CV._TOne * rem;

  CV2._ndxvar = CV._ndxvar;
  if( CV._CM->options.MIXED_IA ) CV2._set_bndT( Op<T>::exp( CV.B() ) );
  if( CV._CM->options.MIN_FACTOR >= 0. ) CV2.simplify( CV._CM->options.MIN_FACTOR );
  return CV2;
}

template <typename T> inline SCVar<T>
log
( SCVar<T> const& CV )
{
  if( !CV._CM )
    return SCVar<T>( Op<T>::log( CV.B() ) );
  if ( Op<T>::l(CV.B()) <= 0. )
    throw typename SCModel<T>::Exceptions( SCModel<T>::Exceptions::LOG );

  SCVar<T> CVI( CV._CM ), CV2( CV._CM );
  auto& coefmon = CV._coefinterp();
  CV._interpolation( coefmon, std::log );

  double b(Op<T>::mid(CV.B())), a(Op<T>::u(CV.B())-b), rem, ub(0), lb(0);
  for (unsigned i(0); i<=CV.maxord(); i++) {
    ub += coefmon[i];
    lb += std::pow(-1.,i)*coefmon[i];
  }
  rem = std::max(std::fabs(std::log(a+b)-ub), std::fabs(std::log(b-a)-lb));

  CVI = CV._rescale(a,b);
  CV2 = CVI._composition( coefmon );
  CV2 += CV._TOne * rem;

  CV2._ndxvar = CV._ndxvar;
  if( CV._CM->options.MIXED_IA ) CV2._set_bndT( Op<T>::log( CV.B() ) );
  if( CV._CM->options.MIN_FACTOR >= 0. ) CV2.simplify( CV._CM->options.MIN_FACTOR );
  return CV2;
}

template <typename T> inline SCVar<T>
xlog
( SCVar<T> const& CV )
{
#ifdef MC__SCVAR_NOINTERP_REM
  return CV * log( CV );
#endif

  if( !CV._CM )
    return SCVar<T>( Op<T>::xlog( CV.B() ) );
  if ( Op<T>::l(CV.B()) <= 0. )
    throw typename SCModel<T>::Exceptions( SCModel<T>::Exceptions::LOG );

  SCVar<T> CVI( CV._CM ), CV2( CV._CM );
  auto& coefmon = CV._coefinterp();
  unsigned nord = CV.nord()+2;
  double TOL = CV._CM->options.INTERP_THRES;
  CV._interpolation( coefmon, TOL, nord, mc::xlog );
  double rem = 2*TOL;
  for( unsigned iord=CV.nord()+1; iord<=nord; iord++ ){
    rem += std::fabs( coefmon[iord] );
#ifdef MC__SCVAR_DEBUG_XLOG
    std::cout << "a[" << iord << "] = " << coefmon[iord] << std::endl;
#endif
  }
#ifdef MC__SCVAR_DEBUG_XLOG
  { int dum; std::cout << "ENTER <1> TO CONTINUE"; std::cin >> dum; }
#endif

  double const m(Op<T>::mid(CV.B())), r(Op<T>::u(CV.B())-m);
  CVI = CV._rescale( r, m );
  CV2 = CVI._composition( coefmon );
  CV2 += CV._TOne * rem;

  CV2._ndxvar = CV._ndxvar;
  if( CV._CM->options.MIXED_IA ) CV2._set_bndT( Op<T>::xlog( CV.B() ) );
  if( CV._CM->options.MIN_FACTOR >= 0. ) CV2.simplify( CV._CM->options.MIN_FACTOR );
  return CV2;
}

template <typename T> inline SCVar<T>
pow
( SCVar<T> const& CV, const int n )
{
  if( !CV._CM )
    return SCVar<T>( Op<T>::pow( CV.B(), n ) );

  if( n < 0 ) return pow( inv( CV ), -n );
  SCVar<T> CV2( CV._CM->_intpow( CV, n ) );
  if( CV._CM->options.MIXED_IA ) CV2._set_bndT( Op<T>::pow( CV.B(), n ) );
  if( CV._CM->options.MIN_FACTOR >= 0. ) CV2.simplify( CV._CM->options.MIN_FACTOR );
  return CV2;
}

template <typename T> inline SCVar<T>
pow
( const SCVar<T> &CV, double const a )
{
  return exp( a * log( CV ) );
}

template <typename T> inline SCVar<T>
pow
( const SCVar<T> &CV1, const SCVar<T> &CV2 )
{
  return exp( CV2 * log( CV1 ) );
}

template <typename T> inline SCVar<T>
pow
( double const a, const SCVar<T> &CV )
{
  return exp( CV * std::log( a ) );
}

template <typename T> inline SCVar<T>
prod
( const unsigned n, const SCVar<T>*CV )
{
  switch( n ){
   case 0:  return 1.;
   case 1:  return CV[0];
   default: return CV[0] * prod( n-1, CV+1 );
  }
}

template <typename T> inline SCVar<T>
monom
( const unsigned n, const SCVar<T>*CV, const unsigned*k )
{
  switch( n ){
   case 0:  return 1.;
   case 1:  return pow( CV[0], (int)k[0] );
   default: return pow( CV[0], (int)k[0] ) * monom( n-1, CV+1, k+1 );
  }
}

template <typename T> inline SCVar<T>
cheb
( const SCVar<T> &CV, const unsigned n )
{
  switch( n ){
    case 0:  return 1.;
    case 1:  return CV;
    default: break;
  }
  SCVar<T> CV2( 2.*(CV*cheb(CV,n-1))-cheb(CV,n-2) );
  if( CV._CM->options.MIXED_IA ) CV2._set_bndT( Op<T>::cheb( CV.B(), n ) );
  if( CV._CM->options.MIN_FACTOR >= 0. ) CV2.simplify( CV._CM->options.MIN_FACTOR );
  return CV2;
}

template <typename T> inline SCVar<T>
cos
( const SCVar<T> &CV )
{
  if( !CV._CM )
    return SCVar<T>( Op<T>::cos( CV.B() ) );

  SCVar<T> CVI( CV._CM ), CV2( CV._CM );
  auto& coefmon = CV._coefinterp();

#ifdef MC__SCVAR_NOINTERP_REM
  CV._interpolation( coefmon, std::cos );
  double m(Op<T>::mid(CV.B())), r(Op<T>::u(CV.B())-m), rem, fact(1);
  for (unsigned i(1); i<=CV.nord()+1; i++) fact *= double(i);
  double M = CV.nord()%2? Op<T>::abs(Op<T>::cos(CV.B())):
                          Op<T>::abs(Op<T>::sin(CV.B()));
  rem = 2.*M*std::pow(r/2.,double(CV.nord()+1))/fact;
#else
  unsigned nord = CV.maxord()+2;
  double TOL = CV._CM->options.INTERP_THRES;
  CV._interpolation( coefmon, TOL, nord, std::cos );
  double rem = 2*TOL;
  for( unsigned iord=CV.maxord()+1; iord<=nord; iord++ ){
    rem += std::fabs( coefmon[iord] );
#ifdef MC__SCVAR_DEBUG_COS
    std::cout << "a[" << iord << "] = " << coefmon[iord] << std::endl;
#endif
  }
#ifdef MC__SCVAR_DEBUG_COS
  { int dum; std::cout << "ENTER <1> TO CONTINUE"; std::cin >> dum; }
#endif
  double const m(Op<T>::mid(CV.B())), r(Op<T>::u(CV.B())-m);  
#endif

  CVI = CV._rescale(r,m);
  CV2 = CVI._composition( coefmon );
  CV2 += CV._TOne * rem;

  CV2._ndxvar = CV._ndxvar;
  if( CV._CM->options.MIXED_IA ) CV2._set_bndT( Op<T>::cos( CV.B() ) );
  if( CV._CM->options.MIN_FACTOR >= 0. ) CV2.simplify( CV._CM->options.MIN_FACTOR );
  return CV2;
}

template <typename T> inline SCVar<T>
sin
( const SCVar<T> &CV )
{
  return cos( CV - PI/2. );
}

template <typename T> inline SCVar<T>
tan
( const SCVar<T> &CV )
{
#ifdef MC__SCVAR_NOINTERP_REM
  return sin( CV ) / cos( CV );
#endif

  if( !CV._CM )
    return SCVar<T>( Op<T>::tan( CV.B() ) );
  if ( Op<T>::l(Op<T>::cos(CV.B())) <= 0. && Op<T>::u(Op<T>::cos(CV.B())) >= 0. )
    throw typename SCModel<T>::Exceptions( SCModel<T>::Exceptions::TAN );

  SCVar<T> CVI( CV._CM ), CV2( CV._CM );
  auto& coefmon = CV._coefinterp();
  unsigned nord = CV.maxord()+2;
  double TOL = CV._CM->options.INTERP_THRES;
  CV._interpolation( coefmon, TOL, nord, std::tan );
  double rem = 2*TOL;
  for( unsigned iord=CV.maxord()+1; iord<=nord; iord++ ){
    rem += std::fabs( coefmon[iord] );
#ifdef MC__SCVAR_DEBUG_TAN
    std::cout << "a[" << iord << "] = " << coefmon[iord] << std::endl;
#endif
  }
#ifdef MC__SCVAR_DEBUG_TAN
  { int dum; std::cout << "ENTER <1> TO CONTINUE"; std::cin >> dum; }
#endif

  double const m(Op<T>::mid(CV.B())), r(Op<T>::u(CV.B())-m);
  CVI = CV._rescale( r, m );
  CV2 = CVI._composition( coefmon );
  CV2 += CV._TOne * rem;

  CV2._ndxvar = CV._ndxvar;
  if( CV._CM->options.MIXED_IA ) CV2._set_bndT( Op<T>::tan( CV.B() ) );
  if( CV._CM->options.MIN_FACTOR >= 0. ) CV2.simplify( CV._CM->options.MIN_FACTOR );
  return CV2;
}

template <typename T> inline SCVar<T>
acos
( const SCVar<T> &CV )
{
  if( !CV._CM )
    return SCVar<T>( Op<T>::acos( CV.B() ) );
  if ( Op<T>::l(CV.B()) < -1. || Op<T>::u(CV.B()) > 1. )
    throw typename SCModel<T>::Exceptions( SCModel<T>::Exceptions::ACOS );

  SCVar<T> CVI( CV._CM ), CV2( CV._CM );
  auto& coefmon = CV._coefinterp();
  unsigned nord = CV.maxord()+2;
  double TOL = CV._CM->options.INTERP_THRES;
  CV._interpolation( coefmon, TOL, nord, std::acos );
  double rem = 2*TOL;
  for( unsigned iord=CV.maxord()+1; iord<=nord; iord++ ){
    rem += std::fabs( coefmon[iord] );
#ifdef MC__SCVAR_DEBUG_ACOS
    std::cout << "a[" << iord << "] = " << coefmon[iord] << std::endl;
#endif
  }
#ifdef MC__SCVAR_DEBUG_ACOS
  { int dum; std::cout << "ENTER <1> TO CONTINUE"; std::cin >> dum; }
#endif

  double const m(Op<T>::mid(CV.B())), r(Op<T>::u(CV.B())-m);
  CVI = CV._rescale( r, m );
  CV2 = CVI._composition( coefmon );
  CV2 += CV._TOne * rem;

  CV2._ndxvar = CV._ndxvar;
  if( CV._CM->options.MIXED_IA ) CV2._set_bndT( Op<T>::acos( CV.B() ) );
  if( CV._CM->options.MIN_FACTOR >= 0. ) CV2.simplify( CV._CM->options.MIN_FACTOR );
  return CV2;
}

template <typename T> inline SCVar<T>
asin
( const SCVar<T> &CV )
{
  return PI/2. - acos( CV );
}

template <typename T> inline SCVar<T>
atan
( const SCVar<T> &CV )
{
#ifdef MC__SCVAR_NOINTERP_REM
  return asin( CV / sqrt( sqr( CV ) + 1. ) );
#endif

  if( !CV._CM )
    return SCVar<T>( Op<T>::atan( CV.B() ) );

  SCVar<T> CVI( CV._CM ), CV2( CV._CM );
  auto& coefmon = CV._coefinterp();
  unsigned nord = CV.maxord()+2;
  double TOL = CV._CM->options.INTERP_THRES;
  CV._interpolation( coefmon, TOL, nord, std::atan );
  double rem = 2*TOL;
  for( unsigned iord=CV.maxord()+1; iord<=nord; iord++ ){
    rem += std::fabs( coefmon[iord] );
#ifdef MC__SCVAR_DEBUG_ATAN
    std::cout << "a[" << iord << "] = " << coefmon[iord] << std::endl;
#endif
  }
#ifdef MC__SCVAR_DEBUG_ATAN
  { int dum; std::cout << "ENTER <1> TO CONTINUE"; std::cin >> dum; }
#endif

  double const m(Op<T>::mid(CV.B())), r(Op<T>::u(CV.B())-m);
  CVI = CV._rescale( r, m );
  CV2 = CVI._composition( coefmon );
  CV2 += CV._TOne * rem;

  CV2._ndxvar = CV._ndxvar;
  if( CV._CM->options.MIXED_IA ) CV2._set_bndT( Op<T>::atan( CV.B() ) );
  if( CV._CM->options.MIN_FACTOR >= 0. ) CV2.simplify( CV._CM->options.MIN_FACTOR );
  return CV2;
}

template <typename T> inline SCVar<T>
sinh
( const SCVar<T> &CV )
{
#ifdef MC__SCVAR_NOINTERP_REM
  return 0.5*(mc::exp(CV)-mc::exp(-CV));
#endif

  if( !CV._CM )
    return SCVar<T>( Op<T>::sinh( CV.B() ) );

  SCVar<T> CVI( CV._CM ), CV2( CV._CM );
  auto& coefmon = CV._coefinterp();
  unsigned nord = CV.maxord()+2;
  double TOL = CV._CM->options.INTERP_THRES;
  CV._interpolation( coefmon, TOL, nord, std::sinh );
  double rem = 2*TOL;
  for( unsigned iord=CV.maxord()+1; iord<=nord; iord++ ){
    rem += std::fabs( coefmon[iord] );
#ifdef MC__SCVAR_DEBUG_SINH
    std::cout << "a[" << iord << "] = " << coefmon[iord] << std::endl;
#endif
  }
#ifdef MC__SCVAR_DEBUG_SINH
  { int dum; std::cout << "ENTER <1> TO CONTINUE"; std::cin >> dum; }
#endif

  double const m(Op<T>::mid(CV.B())), r(Op<T>::u(CV.B())-m);
  CVI = CV._rescale( r, m );
  CV2 = CVI._composition( coefmon );
  CV2 += CV._TOne * rem;

  CV2._ndxvar = CV._ndxvar;
  if( CV._CM->options.MIXED_IA ) CV2._set_bndT( Op<T>::sinh( CV.B() ) );
  if( CV._CM->options.MIN_FACTOR >= 0. ) CV2.simplify( CV._CM->options.MIN_FACTOR );
  return CV2;
}

template <typename T> inline SCVar<T>
cosh
( const SCVar<T> &CV )
{
#ifdef MC__SCVAR_NOINTERP_REM
  return 0.5*(mc::exp(CV)+mc::exp(-CV));
#endif

  if( !CV._CM )
    return SCVar<T>( Op<T>::cosh( CV.B() ) );

  SCVar<T> CVI( CV._CM ), CV2( CV._CM );
  auto& coefmon = CV._coefinterp();
  unsigned nord = CV.maxord()+2;
  double TOL = CV._CM->options.INTERP_THRES;
  CV._interpolation( coefmon, TOL, nord, std::cosh );
  double rem = 2*TOL;
  for( unsigned iord=CV.maxord()+1; iord<=nord; iord++ ){
    rem += std::fabs( coefmon[iord] );
#ifdef MC__SCVAR_DEBUG_COSH
    std::cout << "a[" << iord << "] = " << coefmon[iord] << std::endl;
#endif
  }
#ifdef MC__SCVAR_DEBUG_COSH
  { int dum; std::cout << "ENTER <1> TO CONTINUE"; std::cin >> dum; }
#endif

  double const m(Op<T>::mid(CV.B())), r(Op<T>::u(CV.B())-m);
  CVI = CV._rescale( r, m );
  CV2 = CVI._composition( coefmon );
  CV2 += CV._TOne * rem;

  CV2._ndxvar = CV._ndxvar;
  if( CV._CM->options.MIXED_IA ) CV2._set_bndT( Op<T>::cosh( CV.B() ) );
  if( CV._CM->options.MIN_FACTOR >= 0. ) CV2.simplify( CV._CM->options.MIN_FACTOR );
  return CV2;
}

template <typename T> inline SCVar<T>
tanh
( const SCVar<T> &CV )
{
#ifdef MC__SCVAR_NOINTERP_REM
  return (mc::exp(2*CV)-1)/(mc::exp(2*CV)+1);
#endif

  if( !CV._CM )
    return SCVar<T>( Op<T>::tanh( CV.B() ) );

  SCVar<T> CVI( CV._CM ), CV2( CV._CM );
  auto& coefmon = CV._coefinterp();
  unsigned nord = CV.maxord()+2;
  double TOL = CV._CM->options.INTERP_THRES;
  CV._interpolation( coefmon, TOL, nord, std::tanh );
  double rem = 2*TOL;
  for( unsigned iord=CV.maxord()+1; iord<=nord; iord++ ){
    rem += std::fabs( coefmon[iord] );
#ifdef MC__SCVAR_DEBUG_TANH
    std::cout << "a[" << iord << "] = " << coefmon[iord] << std::endl;
#endif
  }
#ifdef MC__SCVAR_DEBUG_TANH
  { int dum; std::cout << "ENTER <1> TO CONTINUE"; std::cin >> dum; }
#endif

  double const m(Op<T>::mid(CV.B())), r(Op<T>::u(CV.B())-m);
  CVI = CV._rescale( r, m );
  CV2 = CVI._composition( coefmon );
  CV2 += CV._TOne* rem;

  CV2._ndxvar = CV._ndxvar;
  if( CV._CM->options.MIXED_IA ) CV2._set_bndT( Op<T>::tanh( CV.B() ) );
  if( CV._CM->options.MIN_FACTOR >= 0. ) CV2.simplify( CV._CM->options.MIN_FACTOR );
  return CV2;
}

template <typename T>
inline
SCVar<T>
fabs
( const SCVar<T> &CV )
{
  if( !CV._CM )
    return SCVar<T>( Op<T>::fabs( CV.B() ) );
  if ( Op<T>::l(CV.B()) >= 0. )
    return CV;
  if ( Op<T>::u(CV.B()) <= 0. )
    return -CV;

#ifdef MC__SCVAR_FABS_SQRT
  SCVar<T> CV2( sqrt( sqr(CV) ) );
  if( CV._CM->options.MIXED_IA ) CV2._set_bndT( Op<T>::fabs( CV.B() ) );
  return CV2;

#else
  if( CV.maxord() < 2 ){
    SCVar<T> CV2( sqrt( sqr(CV) ) );
    if( CV._CM->options.MIXED_IA ) CV2._set_bndT( Op<T>::fabs( CV.B() ) );
    return CV2;
  }

  SCVar<T> CVI( CV._CM ), CV2( CV._CM );
  auto& coefmon = CV._coefinterp();
  CV._interpolation( coefmon, std::fabs );

  double m(Op<T>::mid(CV.B())), r(Op<T>::u(CV.B())-m), rem, fact(1);
  for (unsigned i(1); i<=CV.maxord()+1; i++) fact *= double(i);
  rem = 2./mc::PI*r/double(CV.maxord()-1);
  //rem = 4.*std::pow(a/2.,double(CV.maxord()+1))/fact;

  CVI = CV._rescale(r,m);
  CV2 = CVI._composition( coefmon );
  CV2 += CV._TOne * rem;

  CV2._ndxvar = CV._ndxvar;
  if( CV._CM->options.MIXED_IA ) CV2._set_bndT( Op<T>::fabs( CV.B() ) );
  if( CV._CM->options.MIN_FACTOR >= 0. ) CV2.simplify( CV._CM->options.MIN_FACTOR );
  return CV2;
#endif
}

template <typename T>
inline
SCVar<T>
hull
( SCVar<T> const& CV1, SCVar<T> const& CV2 )
{
  // Neither operands associated to SCModel -- Make intersection in T type     
  if( !CV1._CM && !CV2._CM ){
    T R1 = CV1.B();
    T R2 = CV2.B();
    return Op<T>::hull(R1, R2);
  }

  // First operand not associated to SCModel
  else if( !CV1._CM )
    return hull( CV2, CV1 );

  // Second operand not associated to SCModel
  else if( !CV2._CM ){
    SCVar<T> CVR = CV1.P();
    return CVR + Op<T>::hull( CV1.R(), CV2._coefmon[0]+CV2._bndrem-CVR.B() );
  }

  // SCModel for first and second operands are inconsistent
  else if( CV1._CM != CV2._CM )
    throw typename SCModel<T>::Exceptions( SCModel<T>::Exceptions::SCMODEL );

  // Perform union
  SCVar<T> CV1C( CV1 ), CV2C( CV2 );
  double const eta = CV1._CM->options.REF_POLY;
  T R1C = CV1C.C().R(), R2C = CV2C.C().R(); 
  CV1C.set(T(0.));
  CV2C.set(T(0.));
  T BCVD = (CV1C-CV2C).B();
  return (1.-eta)*CV1C + eta*CV2C + Op<T>::hull( R1C+eta*BCVD, R2C+(eta-1.)*BCVD );
}

template <typename T>
inline
bool
inter
( SCVar<T>&CVR, SCVar<T> const& CV1, SCVar<T> const& CV2 )
{
  // Neither operands associated to SCModel -- Make intersection in T type     
  if( !CV1._CM && !CV2._CM ){
    T R1 = CV1.B();
    T R2 = CV2.B();
    T RR( 0. );
    bool flag = Op<T>::inter(RR, R1, R2);
    CVR = RR;
    return flag;
  }

  // First operand not associated to SCModel
  else if( !CV1._CM )
    return inter( CVR, CV2, CV1 );

  // Second operand not associated to SCModel
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
    CVR._center();

    if( CVR._CM->options.MIXED_IA ) CVR._set_bndT( BR );
    else CVR._unset_bndT();
    return true;
  }

  // SCModel for first and second operands are inconsistent
  else if( CV1._CM != CV2._CM )
    throw typename SCModel<T>::Exceptions( SCModel<T>::Exceptions::SCMODEL );

  // First intersect in T arithmetic
  T BR;
  if( CV1._CM->options.MIXED_IA && !Op<T>::inter( BR, CV1.B(), CV2.B() ) )
    return false;

  // Perform intersection in PM arithmetic
  SCVar<T> CV1C( CV1 ), CV2C( CV2 );
  double const eta = CV1._CM->options.REF_POLY;
  T R1C = CV1C.C().R(), R2C = CV2C.C().R(); 
  CV1C.set(T(0.));
  CV2C.set(T(0.));
  CVR = (1.-eta)*CV1C + eta*CV2C;
  CV1C -= CV2C;
  T BCVD = CV1C.B();
  if( !Op<T>::inter( CVR._bndrem, R1C+eta*BCVD, R2C+(eta-1.)*BCVD ) )
    return false;
  CVR._center();

  if( CVR._CM->options.MIXED_IA ) CVR._set_bndT( BR );
  else CVR._unset_bndT();
  if( CVR._CM->options.MIN_FACTOR >= 0. ) CVR.simplify( CVR._CM->options.MIN_FACTOR );
  return true;
}

} // namespace mc

#include "mcfadbad.hpp"

namespace fadbad
{

//! @brief Specialization of the structure fadbad::Op for use of the type mc::SCVar of MC++ as a template parameter of the classes fadbad::F, fadbad::B and fadbad::T of FADBAD++
template< typename T > struct Op< mc::SCVar<T> >
{ 
  typedef mc::SCVar<T> CV;
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

//! @brief C++ structure for specialization of the mc::Op templated structure for use of mc::SCVar in DAG evaluation and as template parameter in other MC++ types
template< typename T > struct Op< mc::SCVar<T> >
{
  typedef mc::SCVar<T> CV;
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
  static CV erf (const CV& x) { throw typename mc::SCModel<T>::Exceptions( SCModel<T>::Exceptions::UNDEF ); }
  static CV erfc(const CV& x) { throw typename mc::SCModel<T>::Exceptions( SCModel<T>::Exceptions::UNDEF ); }
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

