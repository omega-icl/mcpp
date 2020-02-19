// Copyright (C) 2020 Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the Eclipse Public License.

/*!
\page page_QUADEXPR Decomposition of Sparse Mutlivariate Polynomials into Quadratic Forms
\author Benoit Chachuat, Tanuj Karia & OMEGA Research Group (http://www3.imperial.ac.uk/environmentenergyoptimisation)
\date 2020
\bug No known bugs.

The class mc::QuadEnv defined in <tt>quadform.hpp</tt> enables the reformulation of sparse polynomial expressions into a set of quadratic forms via the introduction of auxiliary variables.

\section sec_QEXPR_process How Do I Reformulate a (Set of) Sparse Mutlivariate Polynomial(s)?

For illustration, suppose we want to process the factorable function \f${\bf f}:\mathbb{R}^3\to\mathbb{R}^2\f$ defined by
\f{align*}
  {\bf f}(x_0,x_1,x_2) = \left(\begin{array}{c} \left(x_0+x_1^2-2\cdot x_2\right)^3\\ 2\cdot x_1^2-1\end{array}\right)
\f}

The decomposition requires the header file <tt>quadexpr.hpp</tt> to be included:

\code
      #include "quadexpr.hpp"
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

Next, an environment <a>mc::QuadEnv</a> is defined for reformulating sparse multivariate polynomials into quadratic form by calling the method <a>mc::QuadEnv::process</a>:

\code
      mc::QuadEnv QFF( &DAG );
      mc::QuadEnv::options.REDUC = mc::QuadEnv::Options::ALL;
      QFF.process( NF, SPF );
      std::cout << "Sparse quadratic forms: " << QFF << std::endl;
\endcode

The second line sets the options so that all of the reduction quadratic forms are generated. The final line displays the resulting quadratic forms:

\verbatim
    Sparse quadratic form in monomial basis: 

      Monomials: [ 1 X0 X1 X2 X0^2 X0·X1 X0·X2 X1^2 X2^2 X1^4 ]

      Quadratic form for P[1]:
        1.00000e+00    ( X0 , X0^2 )
        1.20000e+01    ( X0 , X2^2 )
        3.00000e+00    ( X0 , X1^4 )
       -6.00000e+00    ( X2 , X0^2 )
       -8.00000e+00    ( X2 , X2^2 )
       -6.00000e+00    ( X2 , X1^4 )
        3.00000e+00    ( X0·X1 , X0·X1 )
       -1.20000e+01    ( X0·X2 , X1^2 )
        1.20000e+01    ( X1^2 , X2^2 )
        1.00000e+00    ( X1^2 , X1^4 )

      Quadratic form for P[2]:
       -1.00000e+00    ( 1 , 1 )
        2.00000e+00    ( X1 , X1 )

      Auxiliary quadratic form #1:
        1.00000e+00    ( 1 , X0^2 )
       -1.00000e+00    ( X0 , X0 )

      Auxiliary quadratic form #2:
        1.00000e+00    ( 1 , X2^2 )
       -1.00000e+00    ( X2 , X2 )

      Auxiliary quadratic form #3:
        1.00000e+00    ( 1 , X0·X1 )
       -1.00000e+00    ( X0 , X1 )

      Auxiliary quadratic form #4:
        1.00000e+00    ( X0 , X0·X1 )
       -1.00000e+00    ( X1 , X0^2 )

      Auxiliary quadratic form #5:
        1.00000e+00    ( 1 , X0·X2 )
       -1.00000e+00    ( X0 , X2 )

      Auxiliary quadratic form #6:
        1.00000e+00    ( X0 , X0·X2 )
       -1.00000e+00    ( X2 , X0^2 )

      Auxiliary quadratic form #7:
        1.00000e+00    ( X1 , X0·X2 )
       -1.00000e+00    ( X2 , X0·X1 )

      Auxiliary quadratic form #8:
        1.00000e+00    ( 1 , X1^2 )
       -1.00000e+00    ( X1 , X1 )

      Auxiliary quadratic form #9:
        1.00000e+00    ( X0 , X1^2 )
       -1.00000e+00    ( X1 , X0·X1 )

      Auxiliary quadratic form #10:
        1.00000e+00    ( X0^2 , X1^2 )
       -1.00000e+00    ( X0·X1 , X0·X1 )

      Auxiliary quadratic form #11:
        1.00000e+00    ( 1 , X1^4 )
       -1.00000e+00    ( X1^2 , X1^2 )
\endverbatim

These results show that a total of 10 monomials participate in the quadratic forms: The constant monomial \f$1\f$; the participating variables \f$x_0,x_1,x_2\f$; and 6 lifted monomials \f$x_3:=x_0^2, x_4:=x_0\cdot x_1, x_5:=x_0\cdot x_2, x_6:=x_1^2, x_7:=x_2^2, x_8:=x_1^4\f$. A reformulation of the factorable function \f${\bf f}\f$ in terms of these monomials is given by:
\f{align*}
  {\bf f}(x_0,\ldots,x_8) = \left(\begin{array}{c} x_0\cdot x_3 + 12\cdot x_0\cdot x_7 + 3\cdot x_0\cdot x_8 - 6\cdot x_2\cdot x_3 - 8\cdot x_2\cdot x_7 - 6\cdot x_2\cdot x_8 + 3\cdot x_4^2 - 12\cdot x_5\cdot x_6 + 12\cdot x_6\cdot x_7 + x_6\cdot x_8\\ 2\cdot x_1^2-1\end{array}\right)
\f}
with the following 11 reduction constraints are also generated:
\f{align*}
  \left(\begin{array}{c} x_3-x_0^2\\ x_4-x_0\cdot x_1\\ x_5-x_0\cdot x_2\\ x_6-x_1^2\\ x_7-x_2^2\\ x_8-x_1^4\\ x_0\cdot x_4 - x_1\cdot x_3\\ x_0\cdot x_5 - x_2\cdot x_3\\ x_1\cdot x_5 - x_2\cdot x_4\\ x_0\cdot x_6 - x_1\cdot x_4\\ x_3\cdot x_6 - x_4^2\end{array}\right) = {\bf 0}
\f}

\section sec_QEXPR_refs References

- Rumschinski, P., Borchers, S., Bosio, S., Weismantel, R., Findeisen, R., <A href="http://www.biomedcentral.com/1752-0509/4/69">Set-base dynamical parameter estimation and
model invalidation for biochemical reaction networks</A>, <i>MBMC Systems Biology</i>, <b>4</b>:69, 2010
- <A href="https://en.wikipedia.org/w/index.php?title=Sum-of-squares_optimization&oldid=929488354">Sum-of-squares optimization</A>, <i>Wikipedia, accessed: 27-Jan-2020
.
*/

// TO DO:
// - Quadratic forms in Chebyshev basis

#ifndef MC__QUADEXPR_H
#define MC__QUADEXPR_H

//#include "sparseexpr.hpp"
#include "spolyexpr.hpp"

#define MC__QUADENV_CHECK
#undef  MC__QUADENV_PROCESS_DEBUG

namespace mc
{
//! @brief C++ structure for ordering of monomial pairs
struct lt_FFQuad
{
  // Comparison operator
  bool operator
    ()
    ( std::pair< FFMon const*, FFMon const* > const& pMon1,
      std::pair< FFMon const*, FFMon const* > const& pMon2 )
    const
    {
      // Order based on first monomial first
      if( lt_FFMon()( *pMon1.first, *pMon2.first ) ) return true;
      if( lt_FFMon()( *pMon2.first, *pMon1.first ) ) return false;
      // Order based on second monomial next
      if( lt_FFMon()( *pMon1.second, *pMon2.second ) ) return true;
      if( lt_FFMon()( *pMon2.second, *pMon1.second ) ) return false;
      // Pairs are identical on reaching this point
      return false;
    }
};

//! @brief C++ class for reformulation of sparse polynomial expressions as a set of quadratic forms
////////////////////////////////////////////////////////////////////////
//! mc::QuadEnv is a C++ class for reformulation of sparse polynomial
//! expressions as a set of quadratic forms via the introduction of
//! auxiliary variables (lifting)
////////////////////////////////////////////////////////////////////////
class QuadEnv
//: public SparseEnv
////////////////////////////////////////////////////////////////////////
{
  friend  std::ostream& operator<< ( std::ostream&, QuadEnv const& );

public:

  typedef std::set< FFMon, lt_FFMon > t_FFMon;
  typedef std::set< FFMon const*, lt_pFFMon > t_pFFMon;
  typedef std::map< std::pair< FFMon const*, FFMon const* >, double, lt_FFQuad > t_FFQuad;
  typedef std::map< FFMon const*, FFVar, lt_pFFMon > t_pFFVar;
  typedef std::map< const FFVar*, const FFVar*, lt_FFVar > t_FFAux;
  
  //! @brief Default Constructor
  QuadEnv
    ( FFGraph* dag )
    : _dag( dag )
    {}

  //! @brief Destructor
  virtual ~QuadEnv
    ()
    { _reset(); }
  
  // Retreive pointer to DAG
  FFGraph* dag
    ()
    const
    { return _dag; };

  //! @brief Process the <a>nPol</a> sparse polynomials in array <a>pPol</a>
  void process
    ( unsigned const nPol, SPolyExpr const* pPol, bool const add2dag=true );

  //! @brief Process the sparse polynomial <a>Pol</a>
  void process
    ( SPolyExpr const& Pol, bool const add2dag=true );

  //! @brief Retreive reference to vector of sparse coefficient matrices defining the main quadratic forms
  std::vector<t_FFQuad> const& MatFct
    ()
    const
    { return _MatFct; }

  //! @brief Retreive reference to vector of sparse coefficient matrices defining the auxiliary quadratic forms
  std::vector<t_FFQuad> const& MatRed
    ()
    const
    { return _MatRed; }

  //! @brief Retreive reference to set of monomials in quadratic forms
  t_FFMon const& SetMon
    ()
    const
    { return _SetMon; }

  //! @brief Retreive reference to map of monomials and their transcription in DAG
  t_pFFVar const& MapMon
    ()
    const
    { return _MapMon; }

  //! @brief Retreive reference to mapping between monomial transcriptions and auxiliaries in DAG
  t_FFAux& MapAux
    ()
    { return _MapAux; }

  //! @brief Retreive reference to vector of quadratic expressions transcribed in DAG
  std::vector<FFVar>& VecQExpr
    ()
    { return _VecQExpr; }

  //! @brief Reset quadratic form expressions
  void reset
    ()
    { _reset(); }

  //! @brief Exceptions of mc::QuadEnv
  class Exceptions
  {
   public:
    //! @brief Enumeration type for QuadEnv exception handling
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
        return "mc::QuadEnv\t Internal error";
      }
    }
   private:
    TYPE _ierr;
  };

  //! @brief Options of mc::QuadEnv
  static struct Options
  {
    //! @brief Constructor
    Options():
      REDUC(ALL)
      {}
    //! @brief Assignment of mc::QuadEnv::Options
    Options& operator=
      ( Options& opt ){
        REDUC = opt.REDUC;
        return *this;
      }
    //! @brief Alternatives in handling of reduction constraints
    enum REDUC_TYPE{
      ONE=0,	 //!< Append a single reduction constraint
      ALL        //!< Append all possible reduction constraints
    };
    //! @brief Handling of reduction constraints
    REDUC_TYPE REDUC;  typedef std::map< const FFVar*, const FFVar*, lt_FFVar > t_Aux;
  } options;  typedef std::map< const FFVar*, const FFVar*, lt_FFVar > t_Aux;

protected:

  //! @brief pointer to underlying dag
  FFGraph* _dag;

  //! @brief Set of monomials in quadratic forms
  t_FFMon _SetMon;

  //! @brief Vector of sparse coefficient matrices defining the main quadratic forms
  std::vector<t_FFQuad> _MatFct;

  //! @brief Vector of sparse coefficient matrices defining the auxiliary quadratic forms
  std::vector<t_FFQuad> _MatRed;

  //! @brief Map of monomials to DAG variables
  t_pFFVar _MapMon;

  //! @brief Map of existing DAG auxiliaries to new DAG variables
  t_FFAux _MapAux;

  //! @brief Vector of quadratic expressions in DAG
  std::vector<FFVar> _VecQExpr;

  //! @brief Reorder entries in a monomial pair
  std::pair< FFMon const*, FFMon const* > _reorder
    ( std::pair< FFMon const*, FFMon const* > pMon )
    const
    {
      if( lt_FFMon()( *pMon.second, *pMon.first ) )
        std::swap( pMon.first, pMon.second );
      return pMon;
    }

  //! @brief Decompose a monomial into product of two monomials in _SetMon by expanding _SetMon if necessary
  std::pair< FFMon const*, FFMon const* > _decompose
    ( FFMon const& mon );
    
  //! @brief Decompose a monomial into product of two monomials in _SetMon by expanding _SetMon if necessary
  typename t_FFMon::iterator _subexpression
    ( FFMon const& mon );

  //! @brief Trancribe quadratic expression into DAG
  FFVar const* _add_to_dag
    ( t_FFQuad const& QExpr )
    const;

  //! @brief Trancribe monomial term into DAG
  FFVar const* _add_to_dag
    ( FFMon const& Mon )
    const;

private:

  //! @brief Reset the quadratic form in _Mon and _Mat
  void _reset
    ();
};

inline QuadEnv::Options QuadEnv::options;

////////////////////////////////////////////////////////////////////////

inline std::ostream&
operator<<
( std::ostream & out, QuadEnv const& env )
{
  const unsigned DISPLEN = SPolyExpr::options.DISPLEN;
  out << std::endl << std::scientific << std::setprecision(DISPLEN)
      << std::right;

  // Output set of monomials
  out << std::endl << "  Monomials: [ ";
  for( auto&& mon : env._SetMon )
    out << mon << " ";
  out << "]" << std::endl;

  // Output quadratic forms
  unsigned count = 0;
  for( auto&& mat : env._MatFct ){
    out << std::endl << "  Quadratic form for P[" << ++count << "]:" << std::endl;
    for( auto&& term : mat )
      out << "   " << std::right << std::setw(DISPLEN+7) << term.second << "   "
          << " ( " << *term.first.first << " , " << *term.first.second << " )"
          << std::endl;
  }
  count = 0;
  for( auto&& mat : env._MatRed ){
    out << std::endl << "  Auxiliary quadratic form #" << ++count << ":" << std::endl;
    for( auto&& term : mat )
      out << "   " << std::right << std::setw(DISPLEN+7) << term.second << "   "
          << " ( " << *term.first.first << " , " << *term.first.second << " )"
          << std::endl;
  }
  return out;
}

inline void
QuadEnv::process
( unsigned const nPol, SPolyExpr const* pPol, bool const add2dag )
{
  for( unsigned i=0; i<nPol; i++ )
    process( pPol[i], add2dag );

  // No transcription in DAG if <a>addtodag</a> is false
  if( !add2dag ) return;
}

inline void
QuadEnv::process
( SPolyExpr const& spoly, bool const add2dag )
{
  // Append entry in <a>MatFct</a>
  _MatFct.push_back( t_FFQuad() );
  auto&& mat = _MatFct.back();

  // Always insert constant monomial 1
  _SetMon.insert( FFMon() );

  // Iterate through monomial terms
  for( auto&& term : spoly.mapmon() ){
    // Insert variables participating in monomial
    for( auto&& varpow : term.first.expr )
      _SetMon.insert( FFMon( *varpow.first ) );
    // Check for product combination in existing monomials
#ifdef MC__QUADENV_CHECK
    auto&& dec = _decompose( term.first );
    assert( dec.first && dec.second );
    auto res = mat.insert( std::make_pair( dec, term.second ) );
    assert( res.second );
#else
    _decompose( term.first );
    mat.insert( std::make_pair( dec, term.second ) );
#endif
  }

  // No transcription in DAG if <a>addtodag</a> is false
  if( !add2dag ) return;
  
  _MapMon.clear();
  _MapAux.clear();
  for( auto&& mon : _SetMon ){
    FFVar const* paux = _add_to_dag( mon );
    if( !paux ) continue;
    _MapMon.insert( std::make_pair( &mon, *paux ) );
      
    // Nothing to do if operand is a DAG constant or leaf variable
    FFVar const* pvar = paux;
    if( paux->ops().first->type != FFOp::VAR
     && paux->ops().first->type != FFOp::CNST )
      pvar = new FFVar( _dag );
    _MapAux.insert( std::make_pair( paux, pvar ) );
  }

  _VecQExpr.clear(); _VecQExpr.reserve( _MatFct.size() + _MatRed.size() );
  for( auto&& qexpr : _MatFct )
    _VecQExpr.push_back( *_add_to_dag( qexpr ) );
  for( auto&& qexpr : _MatRed )
    _VecQExpr.push_back( *_add_to_dag( qexpr ) );
}

inline FFVar const*
QuadEnv::_add_to_dag
( FFMon const& Mon )
const
{
  FFVar prod = 1., one( _dag, 1. );
  for( auto&& term : Mon.expr ){
    switch( SPolyExpr::options.BASIS ){
     case SPolyExpr::Options::MONOM:
      prod *= pow( *term.first, (int)term.second );
      break;
     case SPolyExpr::Options::CHEB:
      prod *= cheb( *term.first, term.second );
      break;
    }
  }
  auto itprod = _dag->Vars().find( &prod );
  if( itprod != _dag->Vars().end() )
    //itprod = _dag->Vars().find( &one );
    return *itprod;
  return 0;
}

inline FFVar const*
QuadEnv::_add_to_dag
( t_FFQuad const& QExpr )
const
{
  FFVar QVar = 0.;
  for( auto&& qterm : QExpr ){
    auto ita1 = _MapMon.find( qterm.first.first );
    auto ita2 = _MapMon.find( qterm.first.second );
    if( ita1 == _MapMon.end() ){
      if( ita2 == _MapMon.end() )
        QVar += qterm.second;
      else{
        FFVar const* pVar2 = _MapAux.find( &ita2->second )->second;
        QVar += qterm.second * *pVar2;
      }
    }
    else{
      FFVar const* pVar1 = _MapAux.find( &ita1->second )->second;
      if( ita2 == _MapMon.end() )
        QVar += qterm.second * *pVar1;
      else{
        FFVar const* pVar2 = _MapAux.find( &ita2->second )->second;
        QVar += qterm.second * ( *pVar1 * *pVar2 );
      }
    }

//    FFVar Aux1 = _MapMon.find( qterm.first.first )->second;
//    FFVar const* pVar1 = _MapAux.find( &Aux1 )->second;
//    FFVar Aux2 = _MapMon.find( qterm.first.second )->second;
//    FFVar const* pVar2 = _MapAux.find( &Aux2 )->second;
//    QVar += qterm.second * ( *pVar1 * *pVar2 );
  }
  auto itvar = _dag->Vars().find( &QVar );
#ifdef MC__QUADENV_CHECK
    assert( itvar != _dag->Vars().end() );
#endif
  return *itvar;
}

inline std::pair< FFMon const*, FFMon const* >
QuadEnv::_decompose
( FFMon const& mon )
{
  std::set< std::pair< FFMon const*, FFMon const* >, lt_FFQuad > CandidateDec;
  for( auto&& mon2 : _SetMon ){
    if( !(mon2 <= mon) ) continue;
    auto&& itmon3 = _SetMon.find( mon - mon2 );
    if( itmon3 == _SetMon.end() || lt_FFMon()( *itmon3, mon2 ) ) continue;
#ifdef MC__QUADENV_PROCESS_DEBUG
    std::cout << "Candidate: " << mon << " = " << mon2 << " · " << *itmon3 << std::endl;
#endif
    CandidateDec.insert( std::make_pair( &mon2, &(*itmon3) ) );
  }
  
  // Case 1: Monomial can be decomposed in terms of existing monomial in _SetMon
  if( !CandidateDec.empty() ){
#ifdef MC__QUADENV_PROCESS_DEBUG
    std::cout << "Decomposed: " << mon << " = " << *CandidateDec.rbegin()->first << " · " << *CandidateDec.rbegin()->second << std::endl;
#endif
    // Prefered candidate isbased on order defined by lt_FFQuad
    return _reorder( *CandidateDec.rbegin() );
  }

  // Case 2: Monomial is linear in all of the variables (multilinear)
  if( mon.gexp() == 1 ){
    t_pFFMon CandidateMon;
    for( auto&& mon2 : _SetMon ){
      if( !(mon2 <= mon) ) continue;
#ifdef MC__QUADENV_CHECK
      assert( _SetMon.find( mon - mon2 ) == _SetMon.end() ); // Covered by Case 1
#endif
      CandidateMon.insert( &mon2 );
#ifdef MC__QUADENV_PROCESS_DEBUG
      std::cout << "Candidate: " << mon - mon2 << " = " << mon << " / " << mon2 << std::endl;
#endif
    }
#ifdef MC__QUADENV_CHECK
    assert( !CandidateMon.empty() ); // _SetMon comprises the participating variables
#endif

    // Case 2a: Use existing non-trivial monomial component in _SetMon
    if( (*CandidateMon.rbegin())->tord > 1 ){
      FFMon const* pmon2 = *CandidateMon.rbegin();
      FFMon mon3( mon - *pmon2 );
#ifdef MC__QUADENV_PROCESS_DEBUG
      std::cout << "Decomposed: " << mon << " = " << *pmon2 << " · " << mon3 << std::endl;
#endif
      // Decompose mon3
      auto&& itmon3 = _subexpression( mon3 );
      return _reorder( std::make_pair( pmon2, &(*itmon3) ) );
    }

    // Case 2b: Split monomial into two monomials of similar total order
    unsigned count = 0;
    FFMon mon2;
    for( auto&& varpow : mon.expr ){
      mon2 += FFMon( *varpow.first );
      if( ++count >= mon.tord / 2 + mon.tord % 2 ) break;
    }
    FFMon mon3( mon - mon2 );
#ifdef MC__QUADENV_PROCESS_DEBUG
    std::cout << "Decomposed: " << mon << " = " << mon2 << " · " << mon3 << std::endl;
#endif
    // Decompose mon2 and mon3
    auto&& itmon2 = _subexpression( mon2 );
    auto&& itmon3 = _subexpression( mon3 );
    return _reorder( std::make_pair( &(*itmon2), &(*itmon3) ) );
  }

  // Case 3: Monomial is linear in some of the variables
  if( mon.lexp() == 1 && mon.gexp() > 1  ){
    FFMon mon2;
    for( auto&& varpow : mon.expr )
      if( varpow.second == 1 )
        mon2 += FFMon( *varpow.first );
#ifdef MC__QUADENV_PROCESS_DEBUG
    std::cout << "Decomposed: " << mon << " = " << mon2 << " · " << mon-mon2 << std::endl;
#endif
    FFMon mon3( mon - mon2 );
    // Decompose mon2 and mon3
    auto&& itmon2 = _subexpression( mon2 );
    auto&& itmon3 = _subexpression( mon3 );
    return _reorder( std::make_pair( &(*itmon2), &(*itmon3) ) );
  }

  // Case 4: Monomial has even partial order in all of the variables
  if( !(mon.gcexp() % 2) ){
    FFMon mon2 = mon / 2;
    // Decompose mon2
    auto&& itmon2 = _subexpression( mon2 );
    return _reorder( std::make_pair( &(*itmon2), &(*itmon2) ) );   
  }

  // Case 5: Monomial has partial order >1 in all of the variables w/ some odd partial order
  FFMon mon2 = mon / 2;
  FFMon mon3( mon - mon2 );
  // Decompose mon2 and mon3
  auto&& itmon2 = _subexpression( mon2 );
  auto&& itmon3 = _subexpression( mon3 );
  return _reorder( std::make_pair( &(*itmon2), &(*itmon3) ) );
}

inline typename QuadEnv::t_FFMon::iterator
QuadEnv::_subexpression
( FFMon const& mon )
{
  // Monomial mon already in _SetMon
  auto&& itmon = _SetMon.find( mon );
  if( itmon != _SetMon.end() )
    return itmon;

  // Monomial mon needs further decomposition
#ifdef MC__QUADENV_CHECK
  auto&& dec = _decompose( mon );
  assert( dec.first && dec.second );
#else
  _decompose( mon );
#endif
  itmon = _SetMon.insert( mon ).first;

  // Append new reduction constraints in <a>MatRed</a>
  auto&& itmon2 = _SetMon.begin();
  for( auto&& mon2 : _SetMon ){
    FFMon montot( *itmon + mon2 );
    auto itmon3 = itmon2;
    for( ++itmon3; itmon3 != _SetMon.end(); ++itmon3 ){
      if( !(*itmon3 <= montot) ) continue;
      auto&& itmon4 = _SetMon.find( montot - *itmon3 );
      if( itmon4 == _SetMon.end() || lt_FFMon()( *itmon4, *itmon3 ) || itmon == itmon4 ) continue;
#ifdef MC__QUADENV_PROCESS_DEBUG
      std::cout << "Reduction: " << mon << " · " << mon2 << " = " << *itmon3 << " · " << *itmon4 << std::endl;
#endif
      _MatRed.push_back( t_FFQuad() );
      auto&& mat = _MatRed.back();
#ifdef MC__QUADENV_CHECK
      assert( mat.insert( std::make_pair( _reorder( std::make_pair( &(*itmon),  &(*itmon2) ) ),  1 ) ).second );
      assert( mat.insert( std::make_pair( _reorder( std::make_pair( &(*itmon3), &(*itmon4) ) ), -1 ) ).second );
#else
      mat.insert( std::make_pair( _reorder( std::make_pair( &(*itmon),  &(*itmon2) ) ),  1 ) );
      mat.insert( std::make_pair( _reorder( std::make_pair( &(*itmon3), &(*itmon4) ) ), -1 ) );
#endif
      // Return after append a single reduction constraint?
      if( options.REDUC == Options::ONE ) return itmon;
    }
    ++itmon2;
  }
  
  return itmon;
}

inline void
QuadEnv::_reset
()
{
  _SetMon.clear();
  _MatFct.clear();
  _MatRed.clear();
  _MapMon.clear();
  _MapAux.clear();
  _VecQExpr.clear();
}

} // namespace mc

#endif
