// Copyright (C) 2018 Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the Eclipse Public License.

/*!
\page page_SPEXPR Manipulation of Sparse Factorable Expressions
\author Benoit Chachuat & OMEGA Research Group (http://www3.imperial.ac.uk/environmentenergyoptimisation)
\date 2018
\bug No known bugs.

The class mc::SparseEnv defined in <tt>sparseexpr.hpp</tt> enables the processing of factorable expressions into sparse polynomial and transcendental subexpressions via the introduction of auxiliary variables.

\section sec_SPEXPR_process How Do I Process a Factorable Expression?

For illustration, suppose we want to process the factorable function \f${\bf f}:\mathbb{R}^2\to\mathbb{R}^2\f$ defined by
\f{align*}
  {\bf f}(x_0,x_1) = \left(\begin{array}{c} \left(x_0+\frac{1}{x_1^2}\right)^3\\ \exp\left(2\cdot x_1^2-1\right)\end{array}\right)
\f}

The decomposition requires the header file <tt>sparseexpr.hpp</tt> to be included:

\code
      #include "sparseexpr.hpp"
\endcode

A DAG of the factorable function is first created:

\code
      mc::FFGraph DAG;
      const unsigned NX = 2, NF = 2;
      mc::FFVar X[NX];
      for( unsigned i(0); i<NX; i++ ) X[i].set( &DAG );
      mc::FFVar F[NF];
      F[0] = pow( X[0] + 1 / sqr( X[1] ), 3 );
      F[1] = exp( 2 * sqr( X[1] ) - 1 );
      std::cout << DAG;
\endcode

The last line displays the following information about the factorable function DAG:

\verbatim
    DAG VARIABLES:
      X0     => { Z3 }
      X1     => { Z0 }

    DAG INTERMEDIATES:
      Z0    <=  SQR( X1 )       => { Z2 Z7 }
      Z2    <=  Z1 / Z0         => { Z3 }
      Z3    <=  X0 + Z2         => { Z5 }
      Z5    <=  POW( Z3, Z4 )   => { }
      Z7    <=  Z0 * Z6         => { Z9 }
      Z9    <=  Z7 + Z8         => { Z10 }
      Z10   <=  EXP( Z9 )       => { }
      Z4    <=  3(I)            => { Z5 }
      Z8    <=  -1(D)           => { Z9 }
      Z1    <=  1(D)            => { Z2 }
      Z6    <=  2(D)            => { Z7 }
\endverbatim

Next, an environment <a>mc::SparseEnv</a> is defined for manipulating the factorable expressions in <a>DAG</a>. The factorable expressions are processed into sparse polynomial and transcendental subexpressions by calling the method <a>mc::SparseEnv::process</a>:

\code
      mc::SparseEnv SPE( &DAG );
      SPE.process( NF, F );
\endcode

The resulting participating variables in the processed expressions, the lifted auxiliary variables, and the resulting expressions can be retreived and displayed as follows:

\code
      std::cout << std::endl << SPE.Var().size() << " participating variables: ";
      for( auto&& var : SPE.Var() ) std::cout << var << " ";
      std::cout << std::endl;
      std::cout << std::endl << SPE.Aux().size() << " auxiliary variables: ";
      for( auto&& aux : SPE.Aux() ) std::cout << *aux.first << "->" << *aux.second << " ";
      std::cout << std::endl;
      std::cout << std::endl << SPE.Poly().size() << " polynomial constraints: " << std::endl;
      for( auto&& expr : SPE.Poly() ) DAG.output( DAG.subgraph( 1, &expr ) );
      std::cout << std::endl;
      std::cout << std::endl << SPE.Trans().size() << " transcendental constraints: " << std::endl;
      for( auto&& expr : SPE.Trans() ) DAG.output( DAG.subgraph( 1, &expr ) );
\endcode

The following information is displayed in this instance:

\verbatim
    7 participating variables: X0 X1 X2 X3 X4 X5 X6 

    5 lifted variables: Z0->X2 Z2->X3 Z5->X6 Z9->X4 Z10->X5 

    4 lifted expressions: 

    FACTORS IN SUBGRAPH:
      X2    <=  VARIABLE
      X1    <=  VARIABLE
      Z0    <=  SQR( X1 )	
      Z11   <=  X2 - Z0	

    FACTORS IN SUBGRAPH:
      X2    <=  VARIABLE
      X3    <=  VARIABLE
      Z12   <=  X2 * X3	
      Z8    <=  -1(D)	
      Z13   <=  Z12 + Z8	

    FACTORS IN SUBGRAPH:
      X4    <=  VARIABLE
      X1    <=  VARIABLE
      Z0    <=  SQR( X1 )	
      Z6    <=  2(D)	
      Z7    <=  Z0 * Z6	
      Z8    <=  -1(D)	
      Z9    <=  Z7 + Z8	
      Z14   <=  X4 - Z9	

    FACTORS IN SUBGRAPH:
      X6    <=  VARIABLE
      X0    <=  VARIABLE
      Z17   <=  3(D)	
      Z18   <=  X0 * Z17	
      X3    <=  VARIABLE
      Z19   <=  SQR( X3 )	
      Z20   <=  Z18 * Z19	
      Z21   <=  SQR( X0 )	
      Z22   <=  Z21 * Z17	
      Z23   <=  X3 * Z22	
      Z24   <=  Z20 + Z23	
      Z4    <=  3(I)	
      Z25   <=  POW( X0, Z4 )
      Z26   <=  Z24 + Z25	
      Z27   <=  POW( X3, Z4 )
      Z28   <=  Z26 + Z27	
      Z29   <=  X6 - Z28	

    1 transcendental constraints: 

    FACTORS IN SUBGRAPH:
      X5	<=  VARIABLE
      X4	<=  VARIABLE
      Z15	<=  EXP( X4 )	
      Z16	<=  X5 - Z15	

\endverbatim

These results show that 5 auxiliary variables have been added to the DAG, \f$x_2,\ldots,x_6\f$. These variables can be determined from the following implicit equations in terms of the original variables \f$x_0,x_1\f$:
\f{align*}
  \left\{\begin{array}{rcl} 0 & = & x_2 - x_1^2\\ 0 & = & x_2\cdot x_3 - 1\\ 0 & = & 2\cdot x_1 - x_4 - 1\\ 0 & = & \exp(x_4) - x_5 \\ 0 & = & x_0^3 + 3\cdot x_0^2\cdot x_3 + 3\cdot x_0\cdot x_3^2 + x_3^3 - x_6  \end{array}\right.
\f}
Finally, the original vector-valued function \f${\bf f}(x_0,x_1)\f$ is equal to \f$(x_6,x_5)^{\sf T}\f$.
*/

// TO DO:
// - Documentation
// - Split subsets of polynomial and transcendental expressions
// - Enforce reusing of intermediate variables in polynomial expressions - as a follow-up reduction step?

#ifndef MC__SPARSEEXPR_H
#define MC__SPARSEEXPR_H

#include "spolyexpr.hpp"

#define MC__SPARSEENV_CHECK
#undef  MC__SPARSEENV_PROCESS_DEBUG

namespace mc
{

class SparseExpr;

//! @brief C++ class for reformulation of a factorable function in sparse polynomial/rational/transcendental subexpressions
////////////////////////////////////////////////////////////////////////
//! mc::SparseEnv is a C++ class for reformulation of a factorable
//! function in sparse polynomial/rational/transcendental subexpressions
////////////////////////////////////////////////////////////////////////
class SparseEnv
////////////////////////////////////////////////////////////////////////
{
  friend class SparseExpr;
  friend  std::ostream& operator<< ( std::ostream&, const SparseEnv& );
  friend  SparseExpr inv ( const SparseExpr& );
  friend  SparseExpr exp ( const SparseExpr& );
  friend  SparseExpr log ( const SparseExpr& );
  friend  SparseExpr xlog( const SparseExpr& );
  friend  SparseExpr sqrt( const SparseExpr& );
  friend  SparseExpr sqr ( const SparseExpr& );
  friend  SparseExpr pow ( const SparseExpr&, const int );  
  friend  SparseExpr pow ( const SparseExpr&, const double );  
  friend  SparseExpr cheb( const SparseExpr&, const unsigned );  
  friend  SparseExpr prod( const unsigned, const SparseExpr* );  
  friend  SparseExpr cos ( const SparseExpr& );
  friend  SparseExpr sin ( const SparseExpr& );
  friend  SparseExpr tan ( const SparseExpr& );
  friend  SparseExpr acos( const SparseExpr& );
  friend  SparseExpr asin( const SparseExpr& );
  friend  SparseExpr atan( const SparseExpr& );
  friend  SparseExpr cosh( const SparseExpr& );
  friend  SparseExpr sinh( const SparseExpr& );
  friend  SparseExpr tanh( const SparseExpr& );
  friend  SparseExpr fabs( const SparseExpr& );
  friend  SparseExpr erf( const SparseExpr& );
  friend  SparseExpr fstep( const SparseExpr& );
  friend  SparseExpr max( const SparseExpr&, const SparseExpr& );  
  friend  SparseExpr min( const SparseExpr&, const SparseExpr& );  
  friend  SparseExpr lmtd( const SparseExpr&, const SparseExpr& );  
  friend  SparseExpr rlmtd( const SparseExpr&, const SparseExpr& );  

public:

  typedef std::list< std::pair< const FFOp*, std::vector<const SparseExpr*> > > t_Interm;
  typedef std::map< const FFVar*, const FFVar*, lt_FFVar > t_Aux;
  typedef std::vector< const FFVar* > t_Expr;

  //! @brief Default Constructor
  SparseEnv
    ( FFGraph* dag )
    : _dag( dag )
    {}

  //! @brief Destructor
  virtual ~SparseEnv
    ()
    { _reset(); }
  
  // Retreive pointer to DAG
  FFGraph* dag
    ()
    const
    { return _dag; };

  //! @brief Retreive reference to vector of new DAG polynmial constraints
  std::vector<FFVar> Poly
    ()
    { return _Poly; }

  //! @brief Retreive reference to vector of new DAG transcendental constraints
  std::vector<FFVar> Trans
    ()
    { return _Trans; }

  //! @brief Retreive reference to vector of independent DAG variables participating in expressions
  std::vector<FFVar> Var
    ()
    { return _Var; }

  //! @brief Retreive reference to mapping between existing DAG auxiliaries and new DAG variables
  t_Aux& Aux
    ()
    { return _Aux; }

  //! @brief Retreive reference to intermediate sparse expressions
  t_Interm& Interm
    ()
    { return _Interm; }

  //! @brief Reset sparse intermediate expressions
  void reset
    ()
    { _reset(); }

  //! @brief Process the <a>nDep</a> dependents in array <a>pDep</a>
  void process
    ( const unsigned nDep, const FFVar*pDep, const bool add2dag=true );

  // List of auxiliary variables corresponding to intermediate SparseExpr
  // Auxiliairies should correspond to already existing DAG variables! -> use dag->curop() and see _append_var in ImgPol

  //! @brief Exceptions of mc::SparseExpr
  class Exceptions
  {
   public:
    //! @brief Enumeration type for SparseExpr exception handling
    enum TYPE{
      DAGERR=0,       //!< Operation involving a sparse expression linked to a different DAG
      ENVERR,         //!< Operation between sparse expressions linked to different environments or without an environment
      INTERNAL=-33    //!< Internal error
    };
    //! @brief Constructor for error <a>ierr</a>
    Exceptions( TYPE ierr ) : _ierr( ierr ){}
    //! @brief Error flag
    int ierr(){ return _ierr; }
    //! @brief Error description
    std::string what(){
      switch( _ierr ){
      case DAGERR:
        return "mc::SparseEnv\t Operation involving a sparse expression linked to a different DAG is not allowed";
      case ENVERR:
        return "mc::SparseEnv\t Operation between sparse rational expressions linked to different environments or without an environment is not allowed";
      case INTERNAL:
      default:
        return "mc::SparseEnv\t Internal error";
      }
    }
   private:
    TYPE _ierr;
  };

  //! @brief Options of mc::SparseEnv
  static struct Options
  {
    //! @brief Constructor
    Options():
      LIFTDIV(true)
      {}
    //! @brief Assignment of mc::SparseEnv::Options
    Options& operator=
      ( Options& opt ){
        LIFTDIV = opt.LIFTDIV;
        return *this;
      }
    //! @brief Whether to lift division terms using auxiliary variables (default: false)
    bool LIFTDIV;
  } options;

private:
  //! @brief pointer to underlying dag
  FFGraph* _dag;

  //! @brief Map of DAG variables and sparse expressions
  t_Interm _Interm;

  //! @brief Map of existing DAG auxiliaries to new DAG variables
  t_Aux _Aux;

  //! @brief Vector of new DAG polynomial constraints
  std::vector<FFVar> _Poly;

  //! @brief Vector of new DAG transcendental constraints
  std::vector<FFVar> _Trans;

  //! @brief Vector of independent DAG variables participating in expressions
  std::vector<FFVar> _Var;

  //! @brief Vector of independent sparse variables
  std::vector<SparseExpr> _SPVar;

  //! @brief Vector of dependent sparse expressions
  std::vector<SparseExpr> _SPDep;

  //! @brief Add new intermediate uni- or bi-variate expression in _Interm
  void _append_interm
    ( const FFOp*op, const SparseExpr*var1, const SparseExpr*var2=0 );

//  //! @brief Add new intermediate n-variate expression in _Interm
//  void _append_interm
//    ( const FFOp*op, const unsigned nvars, const SparseExpr**pvars );

  //! @brief Insert an auxiliary variable corresponding to a rational/polynomial expression into DAG
  const FFVar* _insert_expr
    ( const FFVar*oper, const SparseExpr*expr );

  //! @brief Insert a non-rational operation into DAG via the introduction of auxiliary variables
  void _insert_expr
    ( const FFOp*pOp, std::vector<const FFVar*>&vAux );

  //! @brief Transcribe sparse rational/polynomial expression (SparseExpr) into DAG
  FFVar _SPolyExpr_to_FFVar
    ( const SPolyExpr&expr );

  //! @brief Find auxiliary variable in DAG and return corresponding new DAG variable (or NULL if undefined)
  const FFVar* _find_aux
    ( const FFVar*aux );

  //! @brief Lift current univariate operation
  SparseExpr _lift_univariate_term
    ( const SparseExpr&var );

  //! @brief Lift current bivariate operation
  SparseExpr _lift_bivariate_term
    ( const SparseExpr&var1, const SparseExpr&var2 );

  //! @brief Erase all entries in _Interm
  void _reset
    ();
};

SparseEnv::Options SparseEnv::options;

//! @brief C++ class for sparse rational function representation and arithmetic
////////////////////////////////////////////////////////////////////////
//! mc::SparseExpr is a C++ class for arithmetic manipulation in sparse
//! functions contaning polynomial, rational or transcendental
//! subexpressions
////////////////////////////////////////////////////////////////////////
class SparseExpr
////////////////////////////////////////////////////////////////////////
{
private:
  //! @brief Pointer to sparse rational function environment
  SparseEnv *_env;

protected:
  // numerator sparse polynomial
  SPolyExpr _numer;

  // Denominator sparse polynomial
  SPolyExpr _denom;

  //! @brief Set rational expression equal to <a>var</a>
  SparseExpr& _set
    ( const SparseExpr& var );

  //! @brief Set rational expression equal to constant <a>d</a>
  SparseExpr& _set
    ( const double d );

  //! @brief Set rational expression equal to variable <a>x</a>
  SparseExpr& _set
    ( const FFVar& x );

public:
  /** @ingroup FFunc
   *  @{
   */
  //! @brief Default constructor of sparse rational expression
  SparseExpr
    ()
    {}

  //! @brief Constructor of sparse rational expression as constant
  SparseExpr
    ( const double d )
    : _env( 0 )
    { _set( d ); }

  //! @brief Constructor of sparse rational expression as variable
  SparseExpr
    ( SparseEnv* env, const FFVar& x)
    : _env( env )
    { if( env->dag() != x.dag() )
        throw typename SparseEnv::Exceptions( SparseEnv::Exceptions::DAGERR );
      _set( x ); }

  //! @brief Copy constructor of sparse rational expression
  SparseExpr
    ( const SparseExpr& var )
    { _set( var ); }

  //! @brief Constructor of sparse rational expression as variable
  SparseExpr
    ( SparseEnv* env, const SPolyExpr& n, const SPolyExpr& d )
    : _env( env ), _numer( n ), _denom( d )
    {}

  //! @brief Destructor of sparse rational expression
  virtual ~SparseExpr()
    {}

  //! @brief Initialize variable in sparse rational envrionment <a>env</a> corresponding to DAG variable <a>x</a>
  SparseExpr& set
    ( SparseEnv* env, const FFVar& x )
    { _env = env;
      if( env->dag() != x.dag() )
        throw typename SparseEnv::Exceptions( SparseEnv::Exceptions::DAGERR );
      _set( x ); return *this; }

  //! @brief Overloaded operator '=' for sparse rational expression
  SparseExpr& operator=
    ( const SparseExpr& var )
    { _set( var ); return *this; }

  //! @brief Overloaded operator '=' for constant 
  SparseExpr& operator=
    ( const double d )
    { _env = 0; _set( d ); return *this; }

  //! @brief Overloaded operator '+=' for sparse rational function
  SparseExpr& operator+=
    ( const SparseExpr& var );

  //! @brief Overloaded operator '-=' for sparse rational function
  SparseExpr& operator-=
    ( const SparseExpr& var );

  //! @brief Overloaded operator '*=' for sparse rational function
  SparseExpr& operator*=
    ( const SparseExpr& var );

  //! @brief Overloaded operator '/=' for sparse rational function
  SparseExpr& operator/=
    ( const SparseExpr& var );

  // Sparse rational polynomial environment
  SparseEnv* env
    ()
    const
    { return _env; };

  // numerator sparse polynomial
  const SPolyExpr& numer
    ()
    const
    { return _numer; };

  // Denominator sparse polynomial
  const SPolyExpr& denom
    ()
    const
    { return _denom; };
  /** @} */
};

////////////////////////////////////////////////////////////////////////

inline std::ostream&
operator<<
( std::ostream&out, const SparseExpr&var )
{
  out << std::endl
      << "NUMERATOR:"   << var.numer()
      //<< std::endl
      << "DENOMINATOR:" << var.denom();
  return out;
}

inline std::ostream&
operator<<
( std::ostream&out, const SparseEnv&env)
{
  unsigned count = 0;
  for( auto&& expr : env._Interm ){
    out << std::endl << "Intermediate #" << ++count << ": "
        << *(expr.first->pres) << " = " << *(expr.first) << std::endl;
    unsigned pos = 0;
    for( auto&& oper : expr.second )
      out << "Operand " << *expr.first->pops[pos++] << ": " << *oper;
  }
  return out;
}

inline void
SparseEnv::process
( const unsigned nDep, const FFVar*pDep, const bool add2dag )
{
  // Reset intermediate / auxiliary arrays
  _reset();

  // Update participating variables in _Var
  auto lOps = _dag->subgraph( nDep, pDep );
  for( auto&& Op : lOps ){
    if( Op->type != FFOp::VAR ) continue;
    _Var.push_back( *Op->pres );
    _SPVar.push_back( SparseExpr( this, *Op->pres ) );
  }

#ifdef MC__SPARSEENV_PROCESS_DEBUG
  std::cout << std::endl << _Var.size() << " Original Variables: ";
  for( auto&& var : _Var ) std::cout << var << " ";
  std::cout << std::endl;
#endif

  // Process DAG dependents
  _SPDep.resize( nDep );
  _dag->eval( lOps, nDep, pDep, _SPDep.data(), _Var.size(), _Var.data(), _SPVar.data() );

#ifdef MC__SPARSEENV_PROCESS_DEBUG
  std::cout << *this;
  for( unsigned i=0; i<nDep; i++ )
    std::cout << std::endl << "SPDep[" << i << "]:" << _SPDep[i];
#endif

  // No transcription in DAG if <a>addtodag</a> is false
  if( !add2dag ) return;

  // Insert auxiliary expressions into DAG
  for( auto&& expr : _Interm ){
    // Insert all operands and their defining expressions
    std::vector<const FFVar*> vAux;
    auto itSV = expr.second.begin();
    for( auto&& oper : expr.first->pops ){
#ifdef MC__SPARSEENV_CHECK
      if( itSV == expr.second.end() )
        throw typename SparseEnv::Exceptions( SparseEnv::Exceptions::INTERNAL );
#endif
      vAux.push_back( _insert_expr( oper, *itSV ) );
      ++itSV;
    }
    // Insert operation result and their defining expressions
    _insert_expr( expr.first, vAux );
  }

  // Insert terminal expressions into DAG
  for( unsigned i=0; i<nDep; i++ )
    // Insert operation result and defining expression
    _insert_expr( pDep+i, &_SPDep.at(i) );

#ifdef MC__SPARSEENV_PROCESS_DEBUG
  std::cout << std::endl << _Aux.size() << " Auxiliary Variables: ";
  for( auto&& aux : _Aux ) std::cout << *aux.first << "->" << *aux.second << " ";
  std::cout << std::endl;
  std::cout << std::endl << _Poly.size() << " Polynomial Constraints: " << std::endl;
  for( auto&& expr : _Poly ) _dag->output( _dag->subgraph( 1, &expr ) );
  std::cout << std::endl;
  std::cout << std::endl << _Trans.size() << " Transcendental Constraints: " << std::endl;
  for( auto&& expr : _Trans ) _dag->output( _dag->subgraph( 1, &expr ) );
#endif
}

inline const FFVar*
SparseEnv::_insert_expr
( const FFVar*oper, const SparseExpr*expr )
{ 
  // Nothing to do if operand is a DAG constant or leaf variable
  if( oper->ops().first->type == FFOp::VAR || oper->ops().first->type == FFOp::CNST )
    return oper;

  // Nothing to do if DAG auxiliary was already made a DAG variable
  auto itv = _Aux.find( oper );
  if( itv != _Aux.end() )
    return itv->second;

  // Append new DAG variable in _Aux and defining polynomial constraint in _Poly 
  auto newvar = new FFVar( _dag );
  _Aux.insert( std::make_pair( oper, newvar ) );
  _Var.push_back( *newvar );
  _Poly.push_back( *newvar * _SPolyExpr_to_FFVar( expr->denom() ) - _SPolyExpr_to_FFVar( expr->numer() ) );
  return newvar;
}

inline void
SparseEnv::_insert_expr
( const FFOp*pOp, std::vector<const FFVar*>&vAux )
{ 
#ifdef MC__SPARSEENV_CHECK
  // Throw exception if DAG auxiliary was already made a DAG variable
  if( _Aux.find( pOp->pres ) != _Aux.end() )
    throw typename SparseEnv::Exceptions( SparseEnv::Exceptions::INTERNAL );
#endif

  // Append new DAG variable in _Aux
  auto newvar = new FFVar( _dag );
  _Aux.insert( std::make_pair( pOp->pres, newvar ) );
  _Var.push_back( *newvar );

  // Append defining expression in _Poly or _Transc 
  switch( pOp->type ){
   case FFOp::INV:   return _Poly.push_back( *newvar * *vAux.at(1) - vAux.at(0)->num().val() );
   case FFOp::DIV:   return _Poly.push_back( *newvar * *vAux.at(1) - *vAux.at(0) );
   case FFOp::SQRT:  return _Poly.push_back( sqr( *newvar ) - *vAux.at(0) );
   case FFOp::EXP:   return _Trans.push_back( *newvar - exp( *vAux.at(0) ) );
   case FFOp::LOG:   return _Trans.push_back( exp( *newvar ) - *vAux.at(0) );
   case FFOp::XLOG:  return _Trans.push_back( *newvar - xlog( *vAux.at(0) ) );
   case FFOp::DPOW:  return _Trans.push_back( *newvar - pow( *vAux.at(0), vAux.at(1)->num().val() ) );
   case FFOp::LMTD:  return _Trans.push_back( *newvar - lmtd( *vAux.at(0), *vAux.at(1) ) );
   case FFOp::RLMTD: return _Trans.push_back( *newvar - rlmtd( *vAux.at(0), *vAux.at(1) ) );
   case FFOp::COS:   return _Trans.push_back( *newvar - cos( *vAux.at(0) ) );
   case FFOp::SIN:   return _Trans.push_back( *newvar - sin( *vAux.at(0) ) );
   case FFOp::TAN:   return _Trans.push_back( *newvar - tan( *vAux.at(0) ) );
   case FFOp::ACOS:  return _Trans.push_back( *newvar - acos( *vAux.at(0) ) );
   case FFOp::ASIN:  return _Trans.push_back( *newvar - asin( *vAux.at(0) ) );
   case FFOp::ATAN:  return _Trans.push_back( *newvar - atan( *vAux.at(0) ) );
   case FFOp::COSH:  return _Trans.push_back( *newvar - cosh( *vAux.at(0) ) );
   case FFOp::SINH:  return _Trans.push_back( *newvar - sinh( *vAux.at(0) ) );
   case FFOp::TANH:  return _Trans.push_back( *newvar - tanh( *vAux.at(0) ) );
   case FFOp::ERF:   return _Trans.push_back( *newvar - erf( *vAux.at(0) ) );
   case FFOp::FABS:  return _Trans.push_back( *newvar - fabs( *vAux.at(0) ) );
   case FFOp::FSTEP: return _Trans.push_back( *newvar - fstep( *vAux.at(0) ) );
   case FFOp::MINF:  return _Trans.push_back( *newvar - min( *vAux.at(0), *vAux.at(1) ) );
   case FFOp::MAXF:  return _Trans.push_back( *newvar - max( *vAux.at(0), *vAux.at(1) ) );
   case FFOp::INTER: return _Trans.push_back( *newvar - inter( *vAux.at(0), *vAux.at(1) ) );
   case FFOp::VAR:
   case FFOp::CNST:
   case FFOp::SHIFT:
   case FFOp::PLUS:
   case FFOp::NEG:
   case FFOp::MINUS:
   case FFOp::SCALE:
   case FFOp::TIMES:
   case FFOp::SQR:  
   case FFOp::IPOW:
   case FFOp::CHEB:
   case FFOp::PROD:
   default:          throw typename SparseEnv::Exceptions( SparseEnv::Exceptions::INTERNAL );
  }
}

inline const FFVar*
SparseEnv::_find_aux
( const FFVar*aux )
{
  if( aux->ops().first->type == FFOp::VAR ) return aux;
  auto it = _Aux.find( aux );
  if( it != _Aux.end() ) return it->second;
  return (FFVar*)0;
}

inline FFVar
SparseEnv::_SPolyExpr_to_FFVar
( const SPolyExpr&expr )
{
  FFVar var = 0.;
  for( auto it=expr.mapmon().begin(); it!=expr.mapmon().end(); ++it ){
    FFVar prodmon = it->second;
    for( auto ie=it->first.second.begin(); ie!=it->first.second.end(); ++ie ){
      const FFVar*oper = _find_aux( ie->first );
      if( !oper ) 
        throw typename SparseEnv::Exceptions( SparseEnv::Exceptions::INTERNAL );
      switch( SPolyExpr::options.BASIS ){
       case SPolyExpr::Options::MONOM:
        prodmon *= pow( *oper, (int)ie->second );
        break;
       case SPolyExpr::Options::CHEB:
        prodmon *= cheb( *oper, ie->second );
        break;
      }
    }
    if( it->first.second.empty() ) var = prodmon;
    else                           var+= prodmon;
  }
  // ADD OPTION TO RETURN A PRODMON OR MONOM TERM?
  return var;
}

inline void
SparseEnv::_reset
()
{
  for( auto&& expr : _Interm )
    for( auto&& oper : expr.second )
      delete oper;
  _Interm.clear();

  for( auto&& aux : _Aux )
    delete aux.second;
  _Aux.clear();

  _Poly.clear();
  _Trans.clear();
  _Var.clear();
  _SPVar.clear();
  _SPDep.clear();
}

inline void
SparseEnv::_append_interm
( const FFOp*op, const SparseExpr*var1, const SparseExpr*var2 )
{
#ifdef MC__SPARSEENV_CHECK
  if( !var1 )
    throw typename SparseEnv::Exceptions( SparseEnv::Exceptions::INTERNAL );
#endif
  std::vector<const SparseExpr*> vops;
  vops.push_back( new SparseExpr( *var1 ) );
  if( var2 ) vops.push_back( new SparseExpr( *var2 ) );
  _Interm.push_back( std::make_pair( op, vops ) );
}

//inline bool
//SparseEnv::_append_interm
//( const FFOp*op, const unsigned nvars, const SparseExpr**pvars )
//{
//#ifdef MC__SPARSEENV_CHECK
//  if( !nvars || !pvars )
//    throw typename SparseEnv::Exceptions( SparseEnv::Exceptions::INTERNAL );
//#endif
//  std::vector<const SparseExpr*> vops( pvars, pvars+nvars );
//  _Interm.push_back( std::make_pair( op, vops ) );
//}

inline SparseExpr
SparseEnv::_lift_univariate_term
( const SparseExpr&var )
{
#ifdef MC__SPARSENV_CHECK
  if( !_dag->curOp() )
    throw typename SparseEnv::Exceptions( SparseEnv::Exceptions::INTERNAL );
#endif

  // Append new intermediate expression and assert that same operation was not previously appended
  _append_interm( _dag->curOp(), &var );
  //if( !_append_interm( _dag->curOp(), &var ) )
  //  throw typename SparseEnv::Exceptions( SparseEnv::Exceptions::INTERNAL );

  return SparseExpr( this, *(_dag->curOp()->pres) );
}

inline SparseExpr
SparseEnv::_lift_bivariate_term
( const SparseExpr&var1, const SparseExpr&var2 )
{
#ifdef MC__SPARSENV_CHECK
  if( !_dag->curOp() )
    throw typename SparseEnv::Exceptions( SparseEnv::Exceptions::INTERNAL );
#endif

  // Append new intermediate expression and assert that same operation was not previously appended
  _append_interm( _dag->curOp(), &var1, &var2 );
  //if( !_append_interm( _dag->curOp(), &var1, &var2 ) )
  //  throw typename SparseEnv::Exceptions( SparseEnv::Exceptions::INTERNAL );

  return SparseExpr( this, *(_dag->curOp()->pres) );
}

inline SparseExpr&
SparseExpr::_set
( const double d )
{
  _numer = d;
  _denom = 1;
  return *this;
}

inline SparseExpr&
SparseExpr::_set
( const FFVar&x )
{
  _numer = x;
  _denom = 1;
  return *this;
}

inline SparseExpr&
SparseExpr::_set
( const SparseExpr&var )
{
  if( this == &var ) return *this;
  _env = var._env;
  _numer = var._numer;
  _denom = var._denom;
  return *this;
}

inline SparseExpr
operator+
( const SparseExpr&var )
{
  return var;
}

inline SparseExpr&
SparseExpr::operator+=
( const SparseExpr&var )
{
  if( !_env )
    _env = var._env;
  else if( var._env && _env != var._env )
    throw typename SparseEnv::Exceptions( SparseEnv::Exceptions::ENVERR );

  _numer *= var._denom;
  _numer += var._numer * _denom;
  _denom *= var._denom;
  return *this;
}

inline SparseExpr
operator+
( const SparseExpr&var1, const SparseExpr&var2 )
{
  SparseExpr var3( var1 );
  var3 += var2;
  return var3;
}

inline SparseExpr
operator-
( const SparseExpr&var )
{
  return SparseExpr( var.env(), -var.numer(), var.denom() );
}

inline SparseExpr&
SparseExpr::operator-=
( const SparseExpr&var )
{
  if( !_env )
    _env = var._env;
  else if( var._env && _env != var._env )
    throw typename SparseEnv::Exceptions( SparseEnv::Exceptions::ENVERR );

  _numer *= var._denom;
  _numer -= var._numer * _denom;
  _denom *= var._denom;
  return *this;
}

inline SparseExpr
operator-
( const SparseExpr&var1, const SparseExpr&var2 )
{
  SparseExpr var3( var1 );
  var3 -= var2;
  return var3;
}

inline SparseExpr&
SparseExpr::operator*=
( const SparseExpr&var )
{
  if( !_env )
    _env = var._env;
  else if( var._env && _env != var._env )
    throw typename SparseEnv::Exceptions( SparseEnv::Exceptions::ENVERR );

  _numer *= var._numer;
  _denom *= var._denom;
  return *this;
}

inline SparseExpr
operator*
( const SparseExpr&var1, const SparseExpr&var2 )
{
  SparseExpr var3( var1 );
  var3 *= var2;
  return var3;
}

inline SparseExpr&
SparseExpr::operator/=
( const SparseExpr&var )
{
  if( !_env )
    _env = var._env;
  else if( var._env && _env != var._env )
    throw typename SparseEnv::Exceptions( SparseEnv::Exceptions::ENVERR );

  if( SparseEnv::options.LIFTDIV )
    *this = _env->_lift_bivariate_term( *this, var );
  else{
    _numer *= var._denom;
    _denom *= var._numer;
  }
  return *this;
}

inline SparseExpr
operator/
( const SparseExpr&var1, const SparseExpr&var2 )
{
  SparseExpr var3( var1 );
  var3 /= var2;
  return var3;
}

inline SparseExpr
sqr
( const SparseExpr&var )
{
  return SparseExpr( var.env(), sqr( var.numer() ), sqr( var.denom() ) );
}

inline SparseExpr
inv
( const SparseExpr&var )
{
#ifdef MC__SPARSENV_CHECK
  if( !var.env() )
    throw typename SparseEnv::Exceptions( SparseEnv::Exceptions::ENVERR );
#endif
  if( SparseEnv::options.LIFTDIV )
    return var.env()->_lift_univariate_term( var );
  return SparseExpr( var.env(), var.denom(), var.numer() );
}

inline SparseExpr
pow
( const SparseExpr&var, const int n )
{
  if( n < 0 ) return pow( inv( var ), -n );
  switch( n ){
   case 0:  return 1.;
   case 1:  return var;
   case 2:  return sqr( var );
   default: return SparseExpr( var.env(), pow( var.numer(), (unsigned)n ), pow( var.denom(), (unsigned)n ) );
  }
}

inline SparseExpr
cheb
( const SparseExpr&var, const unsigned n )
{
  switch( n ){
   case 0:  return 1.;
   case 1:  return var;
   case 2:  return 2 * sqr( var ) - 1;
   default: return 2 * var * cheb( var, n-1 ) - cheb( var, n-2 );
  }
}

inline SparseExpr
prod
(const unsigned int nvars, const SparseExpr*pvars)
{
  switch( nvars ){
   case 0:  return 1.;
   case 1:  return pvars[0];
   default: return pvars[0] * prod( nvars-1, pvars+1 );
  }
}

inline SparseExpr
exp
( const SparseExpr&var )
{
#ifdef MC__SPARSENV_CHECK
  if( !var.env() )
    throw typename SparseEnv::Exceptions( SparseEnv::Exceptions::ENVERR );
#endif
  return var.env()->_lift_univariate_term( var );
}

inline SparseExpr
log
( const SparseExpr&var )
{
#ifdef MC__SPARSENV_CHECK
  if( !var.env() )
    throw typename SparseEnv::Exceptions( SparseEnv::Exceptions::ENVERR );
#endif
  return var.env()->_lift_univariate_term( var );
}

inline SparseExpr
xlog
( const SparseExpr&var )
{
#ifdef MC__SPARSENV_CHECK
  if( !var.env() )
    throw typename SparseEnv::Exceptions( SparseEnv::Exceptions::ENVERR );
#endif
  return var.env()->_lift_univariate_term( var );
}

inline SparseExpr
sqrt
( const SparseExpr&var )
{
#ifdef MC__SPARSENV_CHECK
  if( !var.env() )
    throw typename SparseEnv::Exceptions( SparseEnv::Exceptions::ENVERR );
#endif
  return var.env()->_lift_univariate_term( var );
}

inline SparseExpr
cos
( const SparseExpr&var )
{
#ifdef MC__SPARSENV_CHECK
  if( !var.env() )
    throw typename SparseEnv::Exceptions( SparseEnv::Exceptions::ENVERR );
#endif
  return var.env()->_lift_univariate_term( var );
}

inline SparseExpr
sin
( const SparseExpr&var )
{
#ifdef MC__SPARSENV_CHECK
  if( !var.env() )
    throw typename SparseEnv::Exceptions( SparseEnv::Exceptions::ENVERR );
#endif
  return var.env()->_lift_univariate_term( var );
}

inline SparseExpr
tan
( const SparseExpr&var )
{
#ifdef MC__SPARSENV_CHECK
  if( !var.env() )
    throw typename SparseEnv::Exceptions( SparseEnv::Exceptions::ENVERR );
#endif
  return var.env()->_lift_univariate_term( var );
}

inline SparseExpr
acos
( const SparseExpr&var )
{
#ifdef MC__SPARSENV_CHECK
  if( !var.env() )
    throw typename SparseEnv::Exceptions( SparseEnv::Exceptions::ENVERR );
#endif
  return var.env()->_lift_univariate_term( var );
}

inline SparseExpr
asin
( const SparseExpr&var )
{
#ifdef MC__SPARSENV_CHECK
  if( !var.env() )
    throw typename SparseEnv::Exceptions( SparseEnv::Exceptions::ENVERR );
#endif
  return var.env()->_lift_univariate_term( var );
}

inline SparseExpr
atan
( const SparseExpr&var )
{
#ifdef MC__SPARSENV_CHECK
  if( !var.env() )
    throw typename SparseEnv::Exceptions( SparseEnv::Exceptions::ENVERR );
#endif
  return var.env()->_lift_univariate_term( var );
}

inline SparseExpr
cosh
( const SparseExpr&var )
{
#ifdef MC__SPARSENV_CHECK
  if( !var.env() )
    throw typename SparseEnv::Exceptions( SparseEnv::Exceptions::ENVERR );
#endif
  return var.env()->_lift_univariate_term( var );
}

inline SparseExpr
sinh
( const SparseExpr&var )
{
#ifdef MC__SPARSENV_CHECK
  if( !var.env() )
    throw typename SparseEnv::Exceptions( SparseEnv::Exceptions::ENVERR );
#endif
  return var.env()->_lift_univariate_term( var );
}

inline SparseExpr
tanh
( const SparseExpr&var )
{
#ifdef MC__SPARSENV_CHECK
  if( !var.env() )
    throw typename SparseEnv::Exceptions( SparseEnv::Exceptions::ENVERR );
#endif
  return var.env()->_lift_univariate_term( var );
}

inline SparseExpr
fabs
( const SparseExpr&var )
{
#ifdef MC__SPARSENV_CHECK
  if( !var.env() )
    throw typename SparseEnv::Exceptions( SparseEnv::Exceptions::ENVERR );
#endif
  return var.env()->_lift_univariate_term( var );
}

inline SparseExpr
erf
( const SparseExpr&var )
{
#ifdef MC__SPARSENV_CHECK
  if( !var.env() )
    throw typename SparseEnv::Exceptions( SparseEnv::Exceptions::ENVERR );
#endif
  return var.env()->_lift_univariate_term( var );
}

inline SparseExpr
fstep
( const SparseExpr&var )
{
#ifdef MC__SPARSENV_CHECK
  if( !var.env() )
    throw typename SparseEnv::Exceptions( SparseEnv::Exceptions::ENVERR );
#endif
  return var.env()->_lift_univariate_term( var );
}

inline SparseExpr
max
( const SparseExpr&var1, const SparseExpr&var2 )
{
  if( var1.env() && var2.env() && var1.env() != var2.env() )
    throw typename SparseEnv::Exceptions( SparseEnv::Exceptions::ENVERR );

  return var1.env()->_lift_bivariate_term( var1, var2 );
}

inline SparseExpr
min
( const SparseExpr&var1, const SparseExpr&var2 )
{
  if( var1.env() && var2.env() && var1.env() != var2.env() )
    throw typename SparseEnv::Exceptions( SparseEnv::Exceptions::ENVERR );

  return var1.env()->_lift_bivariate_term( var1, var2 );
}

inline SparseExpr
lmtd
( const SparseExpr&var1, const SparseExpr&var2 )
{
  if( var1.env() && var2.env() && var1.env() != var2.env() )
    throw typename SparseEnv::Exceptions( SparseEnv::Exceptions::ENVERR );

  return var1.env()->_lift_bivariate_term( var1, var2 );
}

inline SparseExpr
rlmtd
( const SparseExpr&var1, const SparseExpr&var2 )
{
  if( var1.env() && var2.env() && var1.env() != var2.env() )
    throw typename SparseEnv::Exceptions( SparseEnv::Exceptions::ENVERR );

  return var1.env()->_lift_bivariate_term( var1, var2 );
}

} // namespace mc

//#include "mcop.hpp"

namespace mc
{

//! @brief Specialization of the structure mc::Op to allow usage of the type mc::Interval for DAG evaluation or as a template parameter in other MC++ classes
template <> struct Op<mc::SparseExpr>
{
  typedef mc::SparseExpr T;
  static T point( const double c ) { return T(c); }
  static T zeroone() { throw std::runtime_error("operation not permitted"); }
  static void I(T& x, const T&y) { x = y; }
  static double l(const T& x) { throw std::runtime_error("operation not permitted"); }
  static double u(const T& x) { throw std::runtime_error("operation not permitted"); }
  static double abs (const T& x) { throw std::runtime_error("operation not permitted");  }
  static double mid (const T& x) { throw std::runtime_error("operation not permitted");  }
  static double diam(const T& x) { throw std::runtime_error("operation not permitted"); }
  static T inv (const T& x) { return mc::inv(x);  }
  static T sqr (const T& x) { return mc::sqr(x);  }
  static T sqrt(const T& x) { return mc::sqrt(x); }
  static T exp (const T& x) { return mc::exp(x);  }
  static T log (const T& x) { return mc::log(x);  }
  static T xlog(const T& x) { return mc::xlog(x); }
  static T lmtd(const T& x, const T& y) { return lmtd(x,y); }
  static T rlmtd(const T& x, const T& y) { return rlmtd(x,y); }
  static T fabs(const T& x) { return mc::fabs(x); }
  static T sin (const T& x) { return mc::sin(x);  }
  static T cos (const T& x) { return mc::cos(x);  }
  static T tan (const T& x) { return mc::tan(x);  }
  static T asin(const T& x) { return mc::asin(x); }
  static T acos(const T& x) { return mc::acos(x); }
  static T atan(const T& x) { return mc::atan(x); }
  static T sinh(const T& x) { return mc::sinh(x); }
  static T cosh(const T& x) { return mc::cosh(x); }
  static T tanh(const T& x) { return mc::tanh(x); }
  static T erf (const T& x) { return mc::erf(x);  }
  static T erfc(const T& x) { throw std::runtime_error("operation not permitted"); }//return mc::erfc(x); }
  static T fstep(const T& x) { return mc::fstep(x); }
  static T bstep(const T& x) { throw std::runtime_error("operation not permitted"); }//return mc::bstep(x); }
  static T min (const T& x, const T& y) { return mc::min(x,y);  }
  static T max (const T& x, const T& y) { return mc::max(x,y);  }
  static T arh (const T& x, const double k) { throw std::runtime_error("operation not permitted"); }//return mc::exp( - a / x ); }
  template <typename EXP> static T pow(const T& x, const EXP& y) { return mc::pow(x,y); }
  static T cheb (const T& x, const unsigned n) { return mc::cheb(x,n); }
  static T prod (const unsigned int n, const T* x) { return mc::prod(n,x); }
  static T monom (const unsigned int n, const T* x, const unsigned* k) { throw std::runtime_error("operation not permitted"); }//return mc::monom(n,x,k); }
  static T hull(const T& x, const T& y) { throw std::runtime_error("operation not permitted"); }
  static bool inter(T& xIy, const T& x, const T& y) { throw std::runtime_error("operation not permitted"); }
  static bool eq(const T& x, const T& y) { throw std::runtime_error("operation not permitted"); }
  static bool ne(const T& x, const T& y) { throw std::runtime_error("operation not permitted"); }
  static bool lt(const T& x, const T& y) { throw std::runtime_error("operation not permitted");  }
  static bool le(const T& x, const T& y) { throw std::runtime_error("operation not permitted"); }
  static bool gt(const T& x, const T& y) { throw std::runtime_error("operation not permitted");  }
  static bool ge(const T& x, const T& y) { throw std::runtime_error("operation not permitted"); }
};

} // namespace mc

#endif
