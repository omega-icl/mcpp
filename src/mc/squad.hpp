// Copyright (C) 2020 Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the Eclipse Public License.

/*!
\page page_SQUAD Decomposition of Sparse Polynomials into Quadratic Forms
\author Benoit Chachuat, Tanuj Karia & OMEGA Research Group (http://www3.imperial.ac.uk/environmentenergyoptimisation)
\date 2020
\bug No known bugs.

The class mc::SQuad defined in <tt>scquad.hpp</tt> enables the reformulation of sparse Chebyshev models into a set of quadratic forms via the introduction of auxiliary variables.

\section sec_SQUAD_process How Do I Reformulate a (Set of) Sparse Chebyshev Model(s)?

For illustration, suppose we want to process the factorable function \f${\bf f}:\mathbb{R}^3\to\mathbb{R}^2\f$ defined by
\f{align*}
  {\bf f}(x_0,x_1,x_2) = \left(\begin{array}{c} \left(x_0+x_1^2-2\cdot x_2\right)^3\\ 2\cdot x_1^2-1\end{array}\right)
\f}

The decomposition requires the header file <tt>quadexpr.hpp</tt> to be included:

\code
      #include "squad.hpp"
\endcode

A DAG of the factorable function \f${\bf f}\f$ is first created:

\code
      mc::FFGraph DAG;
      const unsigned NX = 3, NF = 2;
      mc::FFVar X[NX], F[NF];
      for( unsigned i(0); i<NX; i++ ) X[i].set( &DAG );
      F[0] = pow( X[0] + sqr( X[1] ) - 2 * X[2], 3 );
      F[1] = 2 * sqr( X[1] ) - 1;
\endcode

Then, a sparse multivariate polynomial representation of \f${\bf f}\f$ is generated:

\code
  mc::SPolyExpr SPX[NX], SPF[NF];
  for( unsigned i(0); i<NX; i++ ) SPX[i] = X[i];
  DAG.eval( NF, F, SPF, NX, X, SPX );
  
  std::cout << "Sparse polynomial expressions: " << std::endl << std::endl;
  for( unsigned i(0); i<NF; i++ ) std::cout << "P[" << i+1 << "] =" << SPF[i] << std::endl;
\endcode

The last two line display the resulting sparse polynomial expressions:

\verbatim
    Sparse polynomial expressions: 

    P[1] =
       1.00000e+00   X0^3
      -6.00000e+00   X0^2·X2
       1.20000e+01   X0·X2^2
      -8.00000e+00   X2^3
       3.00000e+00   X0^2·X1^2
      -1.20000e+01   X0·X1^2·X2
       1.20000e+01   X1^2·X2^2
       3.00000e+00   X0·X1^4
      -6.00000e+00   X1^4·X2
       1.00000e+00   X1^6

    P[2] =
      -1.00000e+00   1
       2.00000e+00   X1^2
\endverbatim

Next, an environment <a>mc::SQuad</a> is defined for reformulating sparse multivariate polynomials into quadratic form by calling the method <a>mc::SQuad::process</a>:

\code
      mc::SQuad QFF( &DAG );
      mc::SQuad::options.REDUC = mc::SQuad::Options::ALL;
      QFF.process( NF, SPF );
      std::cout << "Sparse quadratic forms: " << QFF << std::endl;
\endcode

The second line sets the options so that all of the reduction quadratic forms are generated. The final line displays the resulting quadratic forms:

\verbatim
Sparse Chebyshev quadratic form:1.62289e-01 CPU-sec

  Monomials: [ 1 T1[0] T1[1] T1[2] T2[0] T1[0]·T1[1] T1[0]·T1[2] T2[1] T1[1]·T1[2] T2[2] T3[1] T4[1] ]

  Quadratic form for P[0]:
    2.81250e-01    ( 1 ; 1 )
    1.25000e-01    ( 1 ; T1[0] )
   -6.25000e+00    ( 1 ; T1[2] )
    0.00000e+00    ( 1 ; T2[0] )
   -6.00000e+00    ( 1 ; T1[0]·T1[2] )
    4.68750e-01    ( 1 ; T2[1] )
    0.00000e+00    ( 1 ; T2[2] )
    1.87500e-01    ( 1 ; T4[1] )
    5.00000e-01    ( T1[0] ; T2[0] )
   -6.00000e+00    ( T1[0] ; T1[0]·T1[2] )
    3.75000e-01    ( T1[0] ; T4[1] )
    3.00000e+00    ( T1[1] ; T1[0]·T1[1] )
    1.20000e+01    ( T1[2] ; T1[0]·T1[2] )
   -3.00000e+00    ( T1[2] ; T2[1] )
   -4.00000e+00    ( T1[2] ; T2[2] )
   -7.50000e-01    ( T1[2] ; T4[1] )
    3.00000e+00    ( T1[0]·T1[1] ; T1[0]·T1[1] )
   -6.00000e+00    ( T1[0]·T1[2] ; T2[1] )
    1.20000e+01    ( T1[1]·T1[2] ; T1[1]·T1[2] )
    6.25000e-02    ( T3[1] ; T3[1] )

  Quadratic form for P[1]:
    0.00000e+00    ( 1 ; 1 )
    1.00000e+00    ( 1 ; T2[1] )

  Auxiliary quadratic form #1:
   -1.00000e+00    ( 1 ; 1 )
   -1.00000e+00    ( 1 ; T2[1] )
    2.00000e+00    ( T1[1] ; T1[1] )

  Auxiliary quadratic form #2:
   -1.00000e+00    ( 1 ; T1[1] )
   -1.00000e+00    ( 1 ; T3[1] )
    2.00000e+00    ( T1[1] ; T2[1] )

  Auxiliary quadratic form #3:
        0.00000e+00    ( 1 ; 1 )
        1.00000e+00    ( 1 ; T2[1] )
        1.84467e+19    ( T1[1] ; T3[1] )
        2.00000e+00    ( T2[1] ; T2[1] )

      Auxiliary quadratic form #4:
       -1.00000e+00    ( 1 ; 1 )
       -1.00000e+00    ( 1 ; T4[1] )
        2.00000e+00    ( T2[1] ; T2[1] )

      Auxiliary quadratic form #5:
        1.00000e+00    ( 1 ; 1 )
       -1.00000e+00    ( 1 ; T1[1] )
        1.00000e+00    ( 1 ; T3[1] )
        1.84467e+19    ( T1[1] ; T4[1] )
        2.00000e+00    ( T2[1] ; T3[1] )

      Auxiliary quadratic form #6:
        0.00000e+00    ( 1 ; 1 )
        1.00000e+00    ( 1 ; T2[1] )
        1.84467e+19    ( T2[1] ; T4[1] )
        2.00000e+00    ( T3[1] ; T3[1] )

      Auxiliary quadratic form #7:
       -1.00000e+00    ( 1 ; T1[1]·T1[2] )
        1.00000e+00    ( T1[1] ; T1[2] )

      Auxiliary quadratic form #8:
        1.00000e+00    ( 1 ; T1[2] )
        1.84467e+19    ( T1[1] ; T1[1]·T1[2] )
        1.00000e+00    ( T1[2] ; T2[1] )

      Auxiliary quadratic form #9:
       -1.00000e+00    ( 1 ; T1[0]·T1[2] )
        1.00000e+00    ( T1[0] ; T1[2] )

      Auxiliary quadratic form #10:
       -1.00000e+00    ( 1 ; T1[0]·T1[1] )
        1.00000e+00    ( T1[0] ; T1[1] )

      Auxiliary quadratic form #11:
       -1.00000e+00    ( 1 ; 1 )
       -1.00000e+00    ( 1 ; T2[2] )
        2.00000e+00    ( T1[2] ; T1[2] )

      Auxiliary quadratic form #12:
       -1.00000e+00    ( 1 ; T1[0] )
        1.84467e+19    ( T1[0] ; T2[2] )
        2.00000e+00    ( T1[2] ; T1[0]·T1[2] )

      Auxiliary quadratic form #13:
       -1.00000e+00    ( 1 ; T1[1] )
        1.84467e+19    ( T1[1] ; T2[2] )
        2.00000e+00    ( T1[2] ; T1[1]·T1[2] )

      Auxiliary quadratic form #14:
       -1.00000e+00    ( 1 ; T1[0]·T1[1] )
        1.84467e+19    ( T1[0]·T1[1] ; T2[2] )
        2.00000e+00    ( T1[0]·T1[2] ; T1[1]·T1[2] )

      Auxiliary quadratic form #15:
       -1.00000e+00    ( 1 ; 1 )
       -1.00000e+00    ( 1 ; T2[1] )
       -1.00000e+00    ( 1 ; T2[2] )
        1.84467e+19    ( T2[1] ; T2[2] )
        4.00000e+00    ( T1[1]·T1[2] ; T1[1]·T1[2] )

      Auxiliary quadratic form #16:
       -1.00000e+00    ( 1 ; 1 )
       -1.00000e+00    ( 1 ; T2[0] )
        2.00000e+00    ( T1[0] ; T1[0] )
\endverbatim

These results show that a total of 10 monomials participate in the quadratic forms: The constant monomial \f$1\f$; the participating variables \f$x_0,x_1,x_2\f$; and 6 lifted monomials \f$x_3:=x_0^2, x_4:=x_0\cdot x_1, x_5:=x_0\cdot x_2, x_6:=x_1^2, x_7:=x_2^2, x_8:=x_1^4\f$. A reformulation of the factorable function \f${\bf f}\f$ in terms of these monomials is given by:
\f{align*}
  {\bf f}(x_0,\ldots,x_8) = \left(\begin{array}{c} x_0\cdot x_3 + 12\cdot x_0\cdot x_7 + 3\cdot x_0\cdot x_8 - 6\cdot x_2\cdot x_3 - 8\cdot x_2\cdot x_7 - 6\cdot x_2\cdot x_8 + 3\cdot x_4^2 - 12\cdot x_5\cdot x_6 + 12\cdot x_6\cdot x_7 + x_6\cdot x_8\\ 2\cdot x_1^2-1\end{array}\right)
\f}
with the following 11 reduction constraints are also generated:
\f{align*}
  \left(\begin{array}{c} x_3-x_0^2\\ x_4-x_0\cdot x_1\\ x_5-x_0\cdot x_2\\ x_6-x_1^2\\ x_7-x_2^2\\ x_8-x_1^4\\ x_0\cdot x_4 - x_1\cdot x_3\\ x_0\cdot x_5 - x_2\cdot x_3\\ x_1\cdot x_5 - x_2\cdot x_4\\ x_0\cdot x_6 - x_1\cdot x_4\\ x_3\cdot x_6 - x_4^2\end{array}\right) = {\bf 0}
\f}

\section sec_SQUAD_refs References

- J Rajyaguru, ME Villanueva, B Houska, B Chachuat, <A href="https://doi.org/10.1007/s10898-016-0474-9">Chebyshev Model Arithmetic for Factorable Functions</A>, <I>Journal of Global Optimization</I> <B>68</B>:413-438, 2017
- P Rumschinski, S Borchers, S Bosio, R Weismantel, R Findeisen, <A href="http://www.biomedcentral.com/1752-0509/4/69">Set-base dynamical parameter estimation and model invalidation for biochemical reaction networks</A>, <i>BMC Systems Biology</i> <b>4</b>:69, 2010
- <A href="https://en.wikipedia.org/w/index.php?title=Sum-of-squares_optimization&oldid=929488354">Sum-of-squares optimization</A>, <i>Wikipedia, accessed: 27-Jan-2020
.
*/

#ifndef MC__SQUAD_H
#define MC__SQUAD_H

#include <list>
//#include "scmodel.hpp"
#include "spolyexpr.hpp"

#define MC__SQUAD_CHECK
#undef  MC__SQUAD_PROCESS_DEBUG

namespace mc
{
//! @brief C++ structure for ordering of monomial pairs
struct lt_SQuad
{
  // Comparison operator
  bool operator
    ()
    ( std::pair< SPolyMon const*, SPolyMon const* > const& pMon1,
      std::pair< SPolyMon const*, SPolyMon const* > const& pMon2 )
    const
    {
      // Order based on first monomial first
      if( lt_SPolyMon()( *pMon1.first, *pMon2.first ) ) return true;
      if( lt_SPolyMon()( *pMon2.first, *pMon1.first ) ) return false;
      // Order based on second monomial next
      if( lt_SPolyMon()( *pMon1.second, *pMon2.second ) ) return true;
      if( lt_SPolyMon()( *pMon2.second, *pMon1.second ) ) return false;
      // Pairs are identical on reaching this point
      return false;
    }
};

//! @brief C++ class for reformulation of sparse Chebyshev models as a set of quadratic forms
////////////////////////////////////////////////////////////////////////
//! mc::SQuad is a C++ class for reformulation of sparse Chebyshev
//! models as a set of quadratic forms via the introduction of
//! auxiliary variables (lifting). Sparse Chebyshev models are
//! computed using the classes mc::SCModel and mc::SCVar. 
////////////////////////////////////////////////////////////////////////
class SQuad
////////////////////////////////////////////////////////////////////////
{
  friend  std::ostream& operator<< ( std::ostream&, SQuad const& );

public:

  typedef std::set< SPolyMon, lt_SPolyMon > t_SPolyMon;
  typedef std::set< SPolyMon const*, lt_pSPolyMon > t_pSPolyMon;
  typedef std::map< SPolyMon, double, lt_SPolyMon > t_SPolyMonCoef;
  typedef std::pair< SPolyMon const*, SPolyMon const* > key_SQuad;
  typedef std::map< key_SQuad, double, lt_SQuad > t_SQuad;
  
  //! @brief Default Constructor
  SQuad
    ()
    {}

  //! @brief Destructor
  virtual ~SQuad
    ()
    { _reset(); }
  
  //! @brief Process the <a>nPol</a> sparse polynomials in array <a>pPol</a>
  double process
    ( unsigned const nSPol, SQuad::t_SPolyMonCoef const* pSPol, bool const check=false );

  //! @brief Process the sparse polynomial <a>Pol</a>
  double process
    ( SQuad::t_SPolyMonCoef const& SPol, bool const check=false );

  //! @brief Decompose the quadratic expression <a>mat</a> into separable expressions
  std::list< t_SQuad > separate
    ( t_SQuad const& mat )
    const;

  //! @brief Factorize the quadratic expression <a>mat</a> using eigenvalue decomposition
  std::multimap< double, t_SPolyMonCoef > factorize
    ( SQuad::t_SQuad const& mat )
    const;

  //! @brief Generate positive semi-definite cuts to tighten the quadratic reformulation
  void tighten
    ();

  //! @brief Retreive reference to vector of sparse coefficient matrices defining the main quadratic forms
  std::vector<t_SQuad> const& MatFct
    ()
    const
    { return _MatFct; }

  //! @brief Retreive reference to vector of sparse coefficient matrices defining the auxiliary quadratic forms
  std::vector<t_SQuad> const& MatRed
    ()
    const
    { return _MatRed; }

  //! @brief Retreive reference to vector of sparse coefficient matrices defining the positive semi-definite cuts
  std::vector<t_SQuad> const& MatPSD
    ()
    const
    { return _MatPSD; }

  //! @brief Retreive reference to set of monomials in quadratic forms
  t_SPolyMon const& SetMon
    ()
    const
    { return _SetMon; }

  //! @brief Reset quadratic form expressions
  void reset
    ()
    { _reset(); }
    
  //! @brief Exceptions of mc::SQuad
  class Exceptions
  {
   public:
    //! @brief Enumeration type for SQuad exception handling
    enum TYPE{
      INTERNAL = -33  //!< Internal error
    };
    //! @brief Constructor for error <a>ierr</a>
    Exceptions( TYPE ierr ) : _ierr( ierr ){}
    //! @brief Error flag
    int ierr(){ return _ierr; }
    //! @brief Error description
    std::string what(){
      switch( _ierr ){
       case INTERNAL:
       default:
        return "mc::SQuad\t Internal error";
      }
    }
   private:
    TYPE _ierr;
  };

  //! @brief Options of mc::SQuad
  static struct Options
  {
    //! @brief Constructor
    Options():
      BASIS(CHEB), REDUC(true), DISPLEN(5)
      {}
    //! @brief Assignment of mc::SQuad::Options
    Options& operator=
      ( Options& opt ){
        BASIS   = opt.BASIS;
        REDUC   = opt.REDUC;
        DISPLEN = opt.DISPLEN;
        return *this;
      }
    //! @brief Available basis representations
    enum BASIS_TYPE{
      MONOM=0,	//!< Monomial basis
      CHEB	    //!< Chebyshev basis
    };
    //! @brief Basis representation of the quadratic form
    int BASIS;
    //! @brief Whether to search for and append extra reduction constraints
    bool REDUC;
    //! @brief Number of digits in output stream for sparse polynomial coefficients
    unsigned DISPLEN;
  } options;

protected:

  //! @brief Set of monomials in quadratic forms
  t_SPolyMon _SetMon;

  //! @brief Vector of sparse coefficient matrices defining the main quadratic forms
  std::vector<t_SQuad> _MatFct;

  //! @brief Vector of sparse coefficient matrices defining the auxiliary quadratic forms
  std::vector<t_SQuad> _MatRed;

  //! @brief Vector of sparse coefficient matrices defining the postiive semi-defninite cuts
  std::vector<t_SQuad> _MatPSD;

  //! @brief Reorder entries in a monomial pair
  std::pair< SPolyMon const*, SPolyMon const* > _reorder
    ( std::pair< SPolyMon const*, SPolyMon const* > pMon )
    const;

  //! @brief Create set of monomials for Chebyshev product between <a>mon1</a> and <a>mon2</a>
  std::set< SPolyMon, lt_SPolyMon > _prodmon
    ( SPolyMon const& mon1, SPolyMon const& mon2 )
    const;

  //! @brief Insert new entry in matrix <a>mat</a> corresponding to the mononial <a>pMon</a> in a product with the constant monomial and with corresponding coefficient <a>coef</a>
  bool _insert
    ( t_SQuad& mat, SPolyMon const* pMon, double const coef );

  //! @brief Insert new entry in matrix <a>mat</a> corresponding to the product between mononials <a>pLMon</a> and <a>pRMon</a> with corresponding coefficient <a>coef</a>
  bool _insert
    ( t_SQuad& mat, SPolyMon const* pLMon, SPolyMon const* pRMon,
      double const coef );

  //! @brief Insert new entry in matrix <a>mat</a> corresponding to the product between mononials <a>pLMon</a> and <a>pRMon</a> with corresponding coefficient <a>coef</a>, and insert the associated low-order terms in monomial map <a>coefmon</a>
  bool _insert
    ( t_SQuad& mat, std::map< SPolyMon, double, lt_SPolyMon >& coefmon, 
      SPolyMon const* pLMon, SPolyMon const* pRMon, double const coef );
      
  //! @brief Decompose <a>mon</a> into product of two monomials in <a>_SetMon</a> by expanding <a>_SetMon</a> as necessary
  std::pair< SPolyMon const*, SPolyMon const* > _decompose
    ( SPolyMon const& mon );
    
  //! @brief Search for <a>mon</a> in <a>_SetMon</a> and append it to <a>_SetMon</a> if missing along reduction constraints in <a>_MatRed</a>
  typename t_SPolyMon::iterator _subexpression
    ( SPolyMon const& mon );

  //! @brief Search for extra reduction constraints for <a>mon</a> and append them to <a>_MatRed</a>
  void _reduction
    ( SPolyMon const& mon );
    
  //! @brief Check quadratic form expressions by comparing with Chebyshev model
  //template <typename T>
  double _check
    ( unsigned const nSPol, SQuad::t_SPolyMonCoef const* pSPol )
    //( unsigned const nscv, SCVar<T> const* pscv )
    const;

  //! @brief Check if a term exist in current quadratic forms 
  bool _find
    ( SPolyMon const* pMon1, SPolyMon const* pMon2 ) 
    const;

  //! @brief Check possible decompositions using existing monomial in _SetMon
  void _candidates
    ( std::set< std::pair< SPolyMon const*, SPolyMon const* >, lt_SQuad >& CandidateDec,
      SPolyMon const& mon )
    const;

  //! @brief Check connections between a quadratic entry and a quadratic form
  bool _isconnected
    ( std::pair< SPolyMon const*, SPolyMon const* > const& entry, 
      t_SQuad const& mat )
    const;
      
private:

  //! @brief Reset the quadratic form
  void _reset
    ();
};

inline SQuad::Options SQuad::options;

////////////////////////////////////////////////////////////////////////

inline std::ostream&
operator<<
( std::ostream& out, SQuad const& quad )
{
  const int BASIS = quad.options.BASIS;
  const unsigned DISPLEN = quad.options.DISPLEN;
  out << std::scientific << std::setprecision(DISPLEN)
      << std::right;

  // Output set of monomials
  out << std::endl << "  Monomials: [ ";
  for( auto const& mon : quad._SetMon )
    out << mon.display( BASIS ) << " ";
  out << "]" << std::endl;

  // Output quadratic forms
  unsigned count = 0;
  for( auto const& mat : quad._MatFct ){
    out << std::endl << "  Quadratic form for P[" << count++ << "]:" << std::endl;
    for( auto const& term : mat )
      out << "   " << std::right << std::setw(DISPLEN+7) << term.second << "   "
          << " ( " << term.first.first->display( BASIS )
          << " ; " << term.first.second->display( BASIS ) << " )"
          << std::endl;
  }
  count = 0;
  for( auto const& mat : quad._MatRed ){
    out << std::endl << "  Auxiliary quadratic form #" << ++count << ":" << std::endl;
    for( auto const& term : mat )
      out << "   " << std::right << std::setw(DISPLEN+7) << term.second << "   "
          << " ( " << term.first.first->display( BASIS )
          << " ; " << term.first.second->display( BASIS ) << " )"
          << std::endl;
  }
  count = 0;
  for( auto const& mat : quad._MatPSD ){
    out << std::endl << "  Positive semi-definite cut #" << ++count << ":" << std::endl;
    for( auto const& term : mat )
      out << "   " << std::right << std::setw(DISPLEN+7) << term.second << "   "
          << " ( " << term.first.first->display( BASIS )
          << " ; " << term.first.second->display( BASIS ) << " )"
          << std::endl;
  }
  return out;
}

inline bool
SQuad::_find
( SPolyMon const* pMon1, SPolyMon const* pMon2 ) 
const
{
  for( auto const& mat : _MatFct )
    if( mat.count( std::make_pair( pMon1, pMon2 ) ) ) return true;
  for( auto const& mat : _MatRed )
    if( mat.count( std::make_pair( pMon1, pMon2 ) ) ) return true;
  return false;
}


inline void
SQuad::tighten
()
{
  _MatPSD.clear();
  std::list< std::pair< SPolyMon const*, SPolyMon const* > > psdlist;
  auto itmon1 = _SetMon.cbegin();
  for( ; itmon1 != _SetMon.cend(); ++itmon1 ){

    // check square monomial *itmon1 participating
    if( itmon1->tord && !_find( &*itmon1, &*itmon1 ) ) continue;
    auto itmon2 = itmon1;
    for( ++itmon2; itmon2 != _SetMon.cend(); ++itmon2 ){

      // check square monomial *itmon2 and cross-term *itmon1.*itmon2 participating
      if( !_find( &*itmon1, &*itmon2 ) || !_find( &*itmon2, &*itmon2 ) ) continue;

      // add positive-definite cuts
      _MatPSD.push_back( t_SQuad() );
      auto& mat1 = _MatPSD.back();
      assert( _insert( mat1, &*itmon1, &*itmon1,  1. )
           && _insert( mat1, &*itmon2, &*itmon2,  1. )
           && _insert( mat1, &*itmon1, &*itmon2, -2. ) );
      if( options.BASIS == Options::MONOM && !itmon1->gcexp()%2 && !itmon2->gcexp()%2 )
        continue;
      _MatPSD.push_back( t_SQuad() );
      auto& mat2 = _MatPSD.back();
      assert( _insert( mat2, &*itmon1, &*itmon1,  1. )
           && _insert( mat2, &*itmon2, &*itmon2,  1. )
           && _insert( mat2, &*itmon1, &*itmon2,  2. ) );
    }
  }
}


inline std::list< SQuad::t_SQuad >
SQuad::separate
( SQuad::t_SQuad const& mat )
const
{
  std::list< t_SQuad > submat;
  for( auto itmat : mat ){
    std::pair< SPolyMon const*, SPolyMon const* > ijmon = itmat.first;
    double coef = itmat.second;
    auto itsubmat = submat.end();
//    // Search for alternative decompositions for terms multiplying monomial '1'
//    // NEED CORRECTION IN ORDER TO ACCOUNT THAT DECOMPOSED TERM MAY ALREADY BE PRESENT IN MAT!
//    if( options.BASIS == Options::MONOM && !ijmon.first->tord && ijmon.second->tord ){
//      std::set< std::pair< SPolyMon const*, SPolyMon const* >, lt_SQuad > CandidateDec;
//      _candidates( CandidateDec, *ijmon.second );
//      ijmon.first  = CandidateDec.crbegin()->first;
//      ijmon.second = CandidateDec.crbegin()->second;
//    }
    // Search connections between entry <a>ijmon</a> and existing terms
    for( auto it=submat.begin(); it!=submat.end(); ){
      // not connected to term
      if( !_isconnected( ijmon, *it ) ){
        ++it;
      }
      // merge terms if more than one connection 
      else if( itsubmat != submat.end() ){
        itsubmat->insert( it->begin(), it->end() );
        it = submat.erase( it );
      }
      // first connection to term
      else{
        (*it)[ijmon] = coef;
        itsubmat = it;
        ++it;
      }
    }
    if( itsubmat != submat.end() ) continue;
    // Create new independent term  if no connections
    submat.push_back( t_SQuad() );
    submat.back()[ijmon] = coef;
  }
  
#ifdef MC__SQUAD_DEBUG_SEPARATE
  unsigned i = 0;
  for( auto const& mat : submat ){
    std::cout << std::endl << "  Separable terms #" << ++i << ":" << std::endl;
    for( auto const& term : mat )
      std::cout << "   " << std::right << std::setw(options.DISPLEN+7) << term.second << "   "
                << " ( " << term.first.first->display( options.BASIS )
                << " ; " << term.first.second->display( options.BASIS ) << " )"
                << std::endl;
  }
#endif
  return submat;
}

inline std::multimap< double, std::map< SPolyMon, double, lt_SPolyMon > >
SQuad::factorize
( SQuad::t_SQuad const& mat )
const
{
  // Populate sparse symmetric coefficient matrix
  CPPL::dssmatrix coefmat;
  std::map< SPolyMon, unsigned, lt_SPolyMon > indexmap;
  unsigned index = 0;
  for( auto const& [ijmon,coef] : mat ){
    auto itimon = indexmap.find( *ijmon.first );
    if( itimon == indexmap.end() ){
      itimon = indexmap.insert( std::make_pair( *ijmon.first, index++ ) ).first;
      coefmat.stretch( 1 );
    }
    // Diagonal entry
    if( ijmon.first == ijmon.second ){
      coefmat.put( itimon->second, itimon->second, coef );
    }
    // Off-diagonal entry
    else{
      auto itjmon = indexmap.find( *ijmon.second );
      if( itjmon == indexmap.end() ){
        itjmon = indexmap.insert( std::make_pair( *ijmon.second, index++ ) ).first;
        coefmat.stretch( 1 );
      }
      coefmat.put( itimon->second, itjmon->second, coef/2e0 );
    }
  }
  
  // Perform eigenvalue decomposition
  std::multimap< double, std::map< SPolyMon, double, lt_SPolyMon > > eigdec;
  std::vector<double> eigval;
  std::vector<CPPL::dcovector> eigvec;
  CPPL::dsymatrix dcoefmat = coefmat.to_dsymatrix();
  if( dcoefmat.dsyev( eigval, eigvec ) ) return eigdec;

  // Populate decomposition map
  auto itval = eigval.begin();
  auto itvec = eigvec.begin();
  for( ; itval != eigval.end(); ++itval, ++itvec ){
    std::map< SPolyMon, double, lt_SPolyMon > eigterm;
    for( auto const& [mon,index] : indexmap ){
      //if( isequal( itvec->array[index], 0. ) ) continue;
      eigterm[mon] = itvec->array[index];
    }
    eigdec.insert( std::make_pair( *itval, eigterm ) );
  }

#ifdef MC__SQUAD_DEBUG_SEPARATE
  unsigned i = 0;
  for( auto const& [eigval,eigterm] : eigdec ){
    std::cout << std::endl << "  Eigen-direction #" << ++i << ": " << std::scientific 
              << std::right << std::setw(options.DISPLEN+7) << eigval << std::endl;
    for( auto const& [mon,coord] : eigterm )
      std::cout << "  " << mon.display(options.BASIS) << ": "
                << std::right << std::setw(options.DISPLEN+7) << coord << std::endl;
  }
#endif
  return eigdec;
}

inline bool
SQuad::_isconnected
( std::pair< SPolyMon const*, SPolyMon const* > const& entry, 
  SQuad::t_SQuad const& mat )
const
{
  for( auto const& [ijmon,coef] : mat ){
    if( entry.first->tord && ( entry.first == ijmon.first || entry.first == ijmon.second ) )
    //if( entry.first->tord && ( entry.first->inter( *ijmon.first ) || entry.first->inter( *ijmon.second ) ) )
      return true;
    if( entry.second->tord && ( entry.second == ijmon.first || entry.second == ijmon.second ) )
    //if( entry.second->tord && ( entry.second->inter( *ijmon.first ) || entry.second->inter( *ijmon.second ) ) )
      return true;
  }
  return false;  
}

//template <typename T>
//inline double
//SQuad::_check
//( unsigned const nscv, SCVar<T> const* pscv )
//const
//{
//  // Determine max monomial order
//  SCModel<T>* pscm = nullptr;
//  for( unsigned i=0; !pscm && i<nscv; i++ )
//    pscm = pscv[i].env(); 

//  // Process entries in _SetMon
//  std::map< SPolyMon, SCVar<T>, lt_SPolyMon > mapscv; 
//  for( auto const& mon : _SetMon ){
//    if( mon.tord == 0 ){
//      mapscv[mon] = SCVar<T>( 1. );
//    }
//    else if( mon.tord == 1 ){
//      unsigned const ivar = mon.expr.begin()->first;
//      mapscv[mon] = SCVar<T>( pscm, ivar, T(-1,1) );
//    }
//    else{
//      SCVar<T> scv( 1. );
//      for( auto const& [ivar,iord] : mon.expr )
//        switch( options.BASIS ){
//          case Options::CHEB:  scv *= cheb( mapscv[SPolyMon(ivar)], iord ); break;
//          case Options::MONOM: scv *= pow( mapscv[SPolyMon(ivar)], (int)iord );  break;
//        }
//      mapscv[mon] = scv;
//    }
//  }

//  // Check entries in _MatFct
//  double sumdiff = 0e0;
//  unsigned ifct = 0;
//  for( auto const& mat : _MatFct ){
//    SCVar<T> scv( 0. );
//    for( auto const& [ijmon,coef] : mat )
//      scv += coef * mapscv[*ijmon.first] * mapscv[*ijmon.second];
//#ifdef MC__SQUAD_DEBUG_CHECK
//    std::cout << "\n  Chebyshev model for quadratic form of P[" << ifct << "]:" << scv;
//#endif
//    sumdiff += Op<T>::abs( ( scv - pscv[ifct++] ).bndpol() );
//  }

//  // Check entries in _MatRed
//#ifdef MC__SQUAD_DEBUG_CHECK
//  unsigned ired = 1;
//#endif
//  for( auto const& mat : _MatRed ){
//    SCVar<T> scv( 0. );
//    for( auto const& [ijmon,coef] : mat )
//      scv += coef * mapscv[*ijmon.first] * mapscv[*ijmon.second];
//#ifdef MC__SQUAD_DEBUG_CHECK
//    std::cout << " \n Chebyshev model for auxiliary quadratic form #" << ired++ << scv;
//#endif
//    sumdiff += Op<T>::abs( scv.bndpol() );
//  }
//  
//  return sumdiff;
//}

inline double
SQuad::_check
( unsigned const nSPol, SQuad::t_SPolyMonCoef const* pSPol )
const
{ 
  assert( nSPol <= _MatFct.size() );
  SPolyExpr::options.BASIS = options.BASIS;
  double sumdiff = 0e0;
  
  // Process entries in _SetMon
  FFGraph dag;
  std::map< unsigned, FFVar > mapvar;
  std::map< SPolyMon, SPolyExpr, lt_SPolyMon > mapmon; 
  for( auto const& mon : _SetMon ){
    if( mon.tord == 0 ){
      mapmon[mon] = SPolyExpr( 1e0 );
    }
    else if( mon.tord == 1 ){
      unsigned const ivar = mon.expr.begin()->first;
      mapvar[ivar] = FFVar( &dag );
      mapmon[mon] = std::make_pair( FFMon( mapvar[ivar] ), 1e0 );
    }
    else{
      FFMon newmon;
      for( auto const& [ivar,iord] : mon.expr )
        assert( newmon.expr.insert( std::make_pair( &mapvar[ivar], iord ) ).second );
      newmon.tord = mon.tord;
      mapmon[mon] = std::make_pair( newmon, 1e0 );
    }
#ifdef MC__SQUAD_DEBUG_CHECK
    std::cout << mapmon[mon] << std::endl;
#endif
  }

  // Check entries in _MatFct
  unsigned ifct = 0;
  for( auto const& mat : _MatFct ){
    if( ifct < _MatFct.size() - nSPol ){
      ifct++;
      continue;
    }

    SPolyExpr spe( 0e0 );
    for( auto const& [mon,coef] : pSPol[ifct] ){
      if( mon.tord == 0 ){
        spe += coef;
        continue;
      }
      auto itmon = mapmon.find( mon );
      if( itmon != mapmon.end() ){
        spe += coef * itmon->second;
      }
      else{ 
        FFMon newmon;
        for( auto const& [ivar,iord] : mon.expr )
          assert( newmon.expr.insert( std::make_pair( &mapvar[ivar], iord ) ).second );
        newmon.tord = mon.tord;
        spe += SPolyExpr( std::make_pair( newmon, coef ) );
      }
    }

    for( auto const& [ijmon,coef] : mat )
      spe -= coef * mapmon[*ijmon.first] * mapmon[*ijmon.second];
#ifdef MC__SQUAD_DEBUG_CHECK
    std::cout << "\n  Quadratic form of P[" << ifct << "]:" << spe;
#endif
    ifct++;

    for( auto const& [mon,coef] : spe.mapmon() )
      sumdiff += std::fabs( coef ); 
  }

  // Check entries in _MatRed
#ifdef MC__SQUAD_DEBUG_CHECK
  unsigned ired = 1;
#endif
  for( auto const& mat : _MatRed ){
    SPolyExpr spe( 0e0 );
    for( auto const& [ijmon,coef] : mat )
      spe += coef * mapmon[*ijmon.first] * mapmon[*ijmon.second];
#ifdef MC__SQUAD_DEBUG_CHECK
    std::cout << " \n Auxiliary quadratic form #" << ired++ << spe;
#endif
    for( auto const& [mon,coef] : spe.mapmon() )
      sumdiff += std::fabs( coef );
  }
  
  return sumdiff;
}

inline double
SQuad::process
( unsigned const nSPol, SQuad::t_SPolyMonCoef const* pSPol, bool const check )
{
  for( unsigned i=0; i<nSPol; i++ )
    process( pSPol[i], false );
  return( check? _check( nSPol, pSPol ): 0. );
}

inline double
SQuad::process
( SQuad::t_SPolyMonCoef const& SPol, bool const check )
{
  // Append entry in <a>MatFct</a>
  _MatFct.push_back( t_SQuad() );
  auto& mat = _MatFct.back();

  // Always insert constant monomial 1
  _SetMon.insert( SPolyMon() );

  // Local copy and conversion to desired basis
  auto coefmon = SPol;

  // Iterate through monomial terms
  for( auto it=coefmon.crbegin(); it!=coefmon.crend(); ++it ){
    const auto& [mon,coef] = *it;

    // Insert variables participating in monomial
    for( auto const& [ivar,iord] : mon.expr )
      _SetMon.insert( SPolyMon( ivar ) );

    // Monomial already in _SetMon
    auto itmon = _SetMon.find( mon );
    if( itmon != _SetMon.end() ){
#ifdef MC__SQUAD_DEBUG_DECOMP
  std::cout << "Decomposed: " << itmon->display(options.BASIS)
            << " = " << _SetMon.cbegin()->display(options.BASIS)
            << " · " << itmon->display(options.BASIS)
            << std::endl;
#endif
      bool ins = _insert( mat, &(*itmon), coef );  
#ifdef MC__SQUAD_CHECK
      assert( ins );
#endif
      continue;
    }
      
    // Mononial needs further decomposition
    auto const& [plmon,prmon] = _decompose( mon );
#ifdef MC__SQUAD_CHECK
    assert( plmon && prmon );
#endif
#ifdef MC__SQUAD_DEBUG_DECOMP
  std::cout << "Decomposed: " << mon.display(options.BASIS)
            << " = " << plmon->display(options.BASIS)
            << " · " << prmon->display(options.BASIS)
            << std::endl;
#endif
    bool ins = _insert( mat, coefmon, plmon, prmon, coef );
#ifdef MC__SQUAD_CHECK
    assert( ins );
#endif
  }

  return( check? _check( 1, &SPol ): 0. );
}

inline std::pair< SPolyMon const*, SPolyMon const* >
SQuad::_reorder
( std::pair< SPolyMon const*, SPolyMon const* > pMon )
const
{
  if( lt_SPolyMon()( *pMon.second, *pMon.first ) )
    std::swap( pMon.first, pMon.second );
  return pMon;
}

inline std::set< SPolyMon, lt_SPolyMon >
SQuad::_prodmon
( SPolyMon const& mon1, SPolyMon const& mon2 )
const
{
#ifdef MC__SQUAD_DEBUG_PRODMON
  std::cout << "mon1:\n" << mon1.display(options.BASIS) << std::endl;
  std::cout << "mon2:\n" << mon2.display(options.BASIS) << std::endl;
#endif
  std::set<SPolyMon,lt_SPolyMon> prodmon;
  switch( options.BASIS ){
   // Monomial basis representation
   case Options::MONOM:
     prodmon.insert( mon1+mon2 );
     break;

   // Chebyshev basis representation
   case Options::CHEB:
    prodmon.insert( SPolyMon() );
    for( auto const& [ivar,iord] : (mon1+mon2).expr ){
      std::set<SPolyMon,lt_SPolyMon> prodmon2;
      auto&& it1 = mon1.expr.find( ivar );
      auto&& it2 = mon2.expr.find( ivar );
      if( it1 != mon1.expr.end() && it2 != mon2.expr.end() ){
        unsigned const& iord1 = it1->second;
        unsigned const& iord2 = it2->second;
        assert( iord1 && iord2 );
        for( auto const& mon3 : prodmon ){
          prodmon2.insert( mon3 + SPolyMon( iord1+iord2, {std::make_pair( ivar, iord1+iord2 )} ) );
          if( iord1 > iord2 )
            prodmon2.insert( mon3 + SPolyMon( iord1-iord2, {std::make_pair( ivar, iord1-iord2 )} ) );
          else if( iord1 < iord2 )
            prodmon2.insert( mon3 + SPolyMon( iord2-iord1, {std::make_pair( ivar, iord2-iord1 )} ) );
          else
            prodmon2.insert( mon3 );
        }
      }
      else if( it1 != mon1.expr.end() ){
        unsigned const& iord1 = it1->second;
        assert( iord1 );
        for( auto const& mon3 : prodmon )
          prodmon2.insert( mon3 + SPolyMon( iord1, {std::make_pair( ivar, iord1 )} ) );
      }    
      else{
        unsigned const& iord2 = it2->second;
        assert( iord2 );
        for( auto const& mon3 : prodmon )
          prodmon2.insert( mon3 + SPolyMon( iord2, {std::make_pair( ivar, iord2 )} ) );
      }    
      std::swap( prodmon, prodmon2 );
    }
    break;
  }
#ifdef MC__SQUAD_DEBUG_PRODMON
  std::cout << "monprod:\n";
  for( auto const& mon : prodmon )
    std::cout << mon.display(options.BASIS) << std::endl;
#endif
  return prodmon;
}

inline bool
SQuad::_insert
( t_SQuad& mat, SPolyMon const* pMon, double const coef )
{
  // New entry in quadratic form as product with constant monomial
  auto [itmat,ins] = mat.insert( std::make_pair( std::make_pair( &(*_SetMon.cbegin()), pMon ), coef ) );
  //if( !ins ) itmat->second += coef;
  //if( itmat->second == 0. ) mat.erase( itmat );
  return ins;
}

inline bool
SQuad::_insert
( t_SQuad& mat, SPolyMon const* pLMon, SPolyMon const* pRMon,
  double const coef )
{
  // New entry in quadratic form
  auto [itmat,ins] = mat.insert( std::make_pair( _reorder( std::make_pair( pLMon, pRMon ) ), coef ) );
  //if( !ins ) itmat->second += coef;
  //if( itmat->second == 0. ) mat.erase( itmat );
  return ins;
}

inline bool
SQuad::_insert
( t_SQuad& mat, std::map< SPolyMon, double, lt_SPolyMon >& coefmon,
  SPolyMon const* pLMon, SPolyMon const* pRMon, double const coef )
{
  // Extra lower-order terms generate by Chebyshev product
  auto&& prodmon = _prodmon( *pLMon, *pRMon );
  unsigned const nprod = prodmon.size();
  auto&& itmon = prodmon.crbegin(); 
  for( ++itmon; itmon != prodmon.crend(); ++itmon ){
    auto [itcmon,ins] = coefmon.insert( std::make_pair( *itmon, -coef ) );
    if( !ins ) itcmon->second -= coef; 
    if( itcmon->second == 0. ) coefmon.erase( itcmon );
  }
  
  // New entry in quadratic form
  auto [itmat,ins] = mat.insert( std::make_pair( _reorder( std::make_pair( pLMon, pRMon ) ), coef * nprod ) );
  //if( !ins ) itmat->second += coef;
  //if( itmat->second == 0. ) mat.erase( itmat );

  return ins;
}

inline void
SQuad::_candidates
( std::set< std::pair< SPolyMon const*, SPolyMon const* >, lt_SQuad >& CandidateDec,
  SPolyMon const& mon )
const
{
  for( auto&& mon2 : _SetMon ){
    if( !mon2.subseteq( mon ) ) continue;
    auto&& itmon3 = _SetMon.find( mon - mon2 );
    if( itmon3 == _SetMon.end() || lt_SPolyMon()( *itmon3, mon2 ) ) continue;
#ifdef MC__SQUAD_DEBUG_DECOMP
    std::cout << "Candidate: " << mon.display(options.BASIS)
              << " = " << mon2.display(options.BASIS)
              << " · " << itmon3->display(options.BASIS)
              << std::endl;
#endif
    CandidateDec.insert( std::make_pair( &mon2, &(*itmon3) ) );
  }
}

inline std::pair< SPolyMon const*, SPolyMon const* >
SQuad::_decompose
( SPolyMon const& mon )
{
  // Possible decompositions using existing monomial in _SetMon
  std::set< std::pair< SPolyMon const*, SPolyMon const* >, lt_SQuad > CandidateDec;
  _candidates( CandidateDec, mon );
  
  // Case 1: Monomial can be decomposed in terms of existing monomial in _SetMon
  if( !CandidateDec.empty() ){
#ifdef MC__SQUAD_DEBUG_DECOMP
    std::cout << "Decomposed: " << mon.display(options.BASIS)
              << " = " << CandidateDec.rbegin()->first->display(options.BASIS)
              << " · " << CandidateDec.rbegin()->second->display(options.BASIS)
              << std::endl;
#endif
    // Prefered candidate is based on order defined by lt_SQuad
    // **unless** monomial already present 'as is'
    if( !CandidateDec.begin()->first->tord ) return *CandidateDec.begin();
    return *CandidateDec.rbegin();
  }

  // Case 2: Monomial is linear in all of the variables (multilinear)
  if( mon.gexp() == 1 ){
    t_pSPolyMon CandidateMon;
    for( auto&& mon2 : _SetMon ){
      if( !mon2.subseteq( mon ) ) continue;
#ifdef MC__SQUAD_CHECK
      assert( _SetMon.find( mon - mon2 ) == _SetMon.end() ); // Covered by Case 1
#endif
      CandidateMon.insert( &mon2 );
#ifdef MC__SQUAD_DEBUG_DECOMP
      std::cout << "Candidate: " << (mon-mon2).display(options.BASIS)
                << " = " << mon.display(options.BASIS)
                << " / " << mon2.display(options.BASIS)
                << std::endl;
#endif
    }
#ifdef MC__SQUAD_CHECK
    assert( !CandidateMon.empty() ); // _SetMon comprises the participating variables
#endif

    // Case 2a: Use existing non-trivial monomial component in _SetMon
    if( (*CandidateMon.rbegin())->tord > 1 ){
      SPolyMon const* pmon2 = *CandidateMon.rbegin();
      SPolyMon mon3( mon - *pmon2 );
#ifdef MC__SQUAD_DEBUG_DECOMP
      std::cout << "Decomposed: " << mon.display(options.BASIS)
                << " = " << pmon2->display(options.BASIS)
                << " · " << mon3.display(options.BASIS)
                << std::endl;
#endif
      // Decompose mon3
      auto&& itmon3 = _subexpression( mon3 );
      //return _reorder( std::make_pair( pmon2, &(*itmon3) ) );
      return std::make_pair( pmon2, &(*itmon3) );
    }

    // Case 2b: Split monomial into two monomials of similar total order
    unsigned count = 0;
    SPolyMon mon2;
    for( auto const& [ivar,iord] : mon.expr ){
#ifdef MC__SQUAD_CHECK
      assert( iord == 1 );
#endif
      mon2 += SPolyMon( ivar );
      if( ++count >= mon.tord / 2 + mon.tord % 2 ) break;
    }
    SPolyMon mon3( mon - mon2 );
#ifdef MC__SQUAD_DEBUG_DECOMP
    std::cout << "Decomposed: " << mon.display(options.BASIS)
              << " = " << mon2.display(options.BASIS)
              << " · " << mon3.display(options.BASIS)
              << std::endl;
#endif
    // Decompose mon2 and mon3
    auto&& itmon2 = _subexpression( mon2 );
    auto&& itmon3 = _subexpression( mon3 );
    //return _reorder( std::make_pair( &(*itmon2), &(*itmon3) ) );
    return std::make_pair( &(*itmon2), &(*itmon3) );
  }

  // Case 3: Monomial is linear in some of the variables
  if( mon.lexp() == 1 && mon.gexp() > 1  ){
    SPolyMon mon2;
    for( auto&& [ivar,iord] : mon.expr )
      if( iord == 1 )
        mon2 += SPolyMon( ivar );
#ifdef MC__SQUAD_DEBUG_DECOMP
    std::cout << "Decomposed: " << mon.display(options.BASIS)
              << " = " << mon2.display(options.BASIS)
              << " · " << (mon-mon2).display(options.BASIS)
              << std::endl;
#endif
    SPolyMon mon3( mon - mon2 );
    // Decompose mon2 and mon3
    auto&& itmon2 = _subexpression( mon2 );
    auto&& itmon3 = _subexpression( mon3 );
    //return _reorder( std::make_pair( &(*itmon2), &(*itmon3) ) );
    return std::make_pair( &(*itmon2), &(*itmon3) );
  }

  // Case 4: Monomial has even partial order in all of the variables
  if( !(mon.gcexp() % 2) ){
    SPolyMon mon2 = mon / 2;
    // Decompose mon2
    auto&& itmon2 = _subexpression( mon2 );
    //return _reorder( std::make_pair( &(*itmon2), &(*itmon2) ) );   
    return std::make_pair( &(*itmon2), &(*itmon2) );
  }

  // Case 5: Monomial has partial order >1 in all of the variables w/ some odd partial order
  SPolyMon mon2 = mon / 2;
  SPolyMon mon3( mon - mon2 );
  // Decompose mon2 and mon3
  auto&& itmon2 = _subexpression( mon2 );
  auto&& itmon3 = _subexpression( mon3 );
  //return _reorder( std::make_pair( &(*itmon2), &(*itmon3) ) );
  return std::make_pair( &(*itmon2), &(*itmon3) );
}

inline typename SQuad::t_SPolyMon::iterator
SQuad::_subexpression
( SPolyMon const& mon )
{
  // Monomial mon already in _SetMon
  auto itmon0 = _SetMon.find( mon );
  if( itmon0 != _SetMon.end() ) return itmon0;

  // Perform further decomposition of monomial <a>mon</a>
  auto const& [plmon,prmon] = _decompose( mon );
#ifdef MC__SQUAD_CHECK
  assert( plmon && prmon );
#endif

  // Append new reduction constraint for <a>mon</a>
  _MatRed.push_back( t_SQuad() );
  auto& mat = _MatRed.back();
  auto itmon = _SetMon.insert( mon ).first;
  std::map< SPolyMon, double, lt_SPolyMon > coefmon;
  bool ins = _insert( mat, &(*itmon), -1. )
          && _insert( mat, coefmon, plmon, prmon, 1. );
#ifdef MC__SQUAD_CHECK
    assert( ins );
#endif
  for( auto it=coefmon.crbegin(); it!=coefmon.crend(); ++it ){
    auto const& [monlow,coeflow] = *it;
    auto itmonlow = _subexpression( monlow );
    ins = _insert( mat, &(*itmonlow), coeflow );
#ifdef MC__SQUAD_CHECK
    assert( ins );
#endif
  }
  
  // Search for extra reduction constraints for <a>mon</a>
  if( options.REDUC ) _reduction( *itmon );
   
  return itmon;
}

inline void
SQuad::_reduction
( SPolyMon const& mon )
{
  // Search for extra reduction constraints for <a>mon</a>
  // that don't add additional low-order monomials
  auto itmon2 = _SetMon.begin();
  for( ++itmon2; itmon2 != _SetMon.end(); ++itmon2 ){
  
    // Low-order terms generated by Chebyshev product
    auto&& prodmon12 = _prodmon( mon, *itmon2 );
    bool monmis = false;
    auto itprodmon12 = prodmon12.crbegin(); 
    // Check that all low-order terms are present -> OR INTRODUCE THEM?
    for( ++itprodmon12; !monmis && itprodmon12 != prodmon12.crend(); ++itprodmon12 )
       if( _SetMon.find( *itprodmon12 ) == _SetMon.end() ) 
         //_subexpression( *itprodmon12 );
         monmis = true;
    if( monmis ) continue;

    // Find alternative decompositions
    SPolyMon montot( mon + *itmon2 );
    auto itmon3 = itmon2;
    for( ++itmon3; itmon3 != _SetMon.end(); ++itmon3 ){
      if( !itmon3->subseteq( montot ) ) continue;
      auto itmon4 = _SetMon.find( montot - *itmon3 );
      if( itmon4 == _SetMon.end() || lt_SPolyMon()( *itmon4, *itmon3 ) ) continue;

      // Low-order terms generated by Chebyshev product
      auto&& prodmon34 = _prodmon( *itmon3, *itmon4 );
      bool monmis = false;
      auto itprodmon34 = prodmon34.crbegin(); 
      // Check that all low-order terms are present -> OR INTRODUCE THEM?
      for( ++itprodmon34; !monmis && itprodmon34 != prodmon34.crend(); ++itprodmon34 )
         if( _SetMon.find( *itprodmon34 ) == _SetMon.end() ) 
           //_subexpression( *itprodmon34 );
           monmis = true;
      if( monmis ) continue;
#ifdef MC__SQUAD_DEBUG_REDUC
      std::cout << "Reduction: " << mon.display(options.BASIS)
                << " · " << itmon2->display(options.BASIS)
                << " == " << itmon3->display(options.BASIS)
                << " · " << itmon4->display(options.BASIS)
                << std::endl;
#endif

      // Populate reduction constraint
      _MatRed.push_back( t_SQuad() );
      auto& mat = _MatRed.back();
      // Insert product terms
      bool ins = _insert( mat, &mon,       &(*itmon2), -(double)prodmon12.size() )
              && _insert( mat, &(*itmon3), &(*itmon4),  (double)prodmon34.size() );
#ifdef MC__SQUAD_CHECK
      assert( ins );
#endif
      // Insert low-order terms 
      itprodmon12 = prodmon12.crbegin(); 
      for( ++itprodmon12; itprodmon12 != prodmon12.crend(); ++itprodmon12 ){
        ins = _insert( mat, &(*_SetMon.find( *itprodmon12 )),  1. );
#ifdef MC__SQUAD_CHECK
        assert( ins );
#endif
      }
      itprodmon34 = prodmon34.crbegin(); 
      for( ++itprodmon34; itprodmon34 != prodmon34.crend(); ++itprodmon34 ){
        ins = _insert( mat, &(*_SetMon.find( *itprodmon34 )), -1. );
#ifdef MC__SQUAD_CHECK
        assert( ins );
#endif
      }
    }
  }
}

inline void
SQuad::_reset
()
{
  _SetMon.clear();
  _MatFct.clear();
  _MatRed.clear();
}

} // namespace mc

#endif
