// Copyright (C) Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the Eclipse Public License.

/*!
\page page_SElim Variable Elimination in Factorable Equality Constraints 
\author Benoit Chachuat & Dominik Bongartz
\date 2023
\bug No known bugs.

The classes mc::SElimEnv and mc::SELIM defined in <tt>selim.hpp</tt> enable the elimination of (a subset of) variables from factorable expressions through expoiting solvable (invertible) equality constraints.

\section sec_SElim_process How do I eliminate variables using factorable equality constraints?

For illustration, consider the factorable function \f${\bf f}:\mathbb{R}^2\to\mathbb{R}^2\f$ defined by
\f{align*}
  {\bf f}(x_0,x_1) = \left(\begin{array}{c} \left(x_0+\frac{1}{x_1^2}\right)^3\\ \exp\left(2\cdot x_1^2-1\right)\end{array}\right)
\f}

The lifting requires the header file <tt>selim.hpp</tt> to be included:

\code
      #include "selim.hpp"
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

The last line displays the following information about the DAG:

\verbatim
    DAG VARIABLES:
      X0     => { Z3 }
      X1     => { Z0 }

    DAG INTERMEDIATES:
      Z0    <=  SQR( X1 )       => { Z2 Z7 }
      Z2    <=  Z1 / Z0         => { Z3 }
      Z3    <=  X0 + Z2         => { Z5 }
      Z5    <=  POW( Z3, Z4 )   => { }
      Z7    <=  Z0 x Z6         => { Z9 }
      Z9    <=  Z7 + Z8         => { Z10 }
      Z10   <=  EXP( Z9 )       => { }
      Z4    <=  3(I)            => { Z5 }
      Z8    <=  -1(D)           => { Z9 }
      Z1    <=  1(D)            => { Z2 }
      Z6    <=  2(D)            => { Z7 }
\endverbatim

Next, an environment <a>mc::SElimEnv</a> is defined for lifting the factorable expressions in <a>DAG</a>. The method <a>mc::SElimEnv::process</a> decomposes the factorable expressions recursively into sparse polynomial and transcendental subexpressions:

\code
      mc::SElimEnv<mc::FFGraph<>> SPE( &DAG );
      SPE.process( NF, F );
\endcode

The resulting participating variables in the processed expressions, the lifted variables, and the resulting subexpressions can be retreived and displayed as follows:

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

    5 auxiliary variables: Z0->X2 Z2->X3 Z5->X6 Z9->X4 Z10->X5

    4 polynomial constraints: 

    FACTORS IN SUBGRAPH:
      X2    <=  VARIABLE
      X1    <=  VARIABLE
      Z0    <=  SQR( X1 )	
      Z11   <=  X2 - Z0	

    FACTORS IN SUBGRAPH:
      X2    <=  VARIABLE
      X3    <=  VARIABLE
      Z12   <=  X2 x X3	
      Z8    <=  -1(D)	
      Z13   <=  Z12 + Z8	

    FACTORS IN SUBGRAPH:
      X4    <=  VARIABLE
      X1    <=  VARIABLE
      Z0    <=  SQR( X1 )	
      Z6    <=  2(D)	
      Z7    <=  Z0 x Z6	
      Z8    <=  -1(D)	
      Z9    <=  Z7 + Z8	
      Z14   <=  X4 - Z9	

    FACTORS IN SUBGRAPH:
      X6    <=  VARIABLE
      X0    <=  VARIABLE
      Z17   <=  3(D)	
      Z18   <=  X0 x Z17	
      X3    <=  VARIABLE
      Z19   <=  SQR( X3 )	
      Z20   <=  Z18 x Z19	
      Z21   <=  SQR( X0 )	
      Z22   <=  Z21 x Z17	
      Z23   <=  X3 x Z22	
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
  \left\{\begin{array}{rcl} 0 & = & x_2 - x_1^2\\ 0 & = & x_2\cdot x_3 - 1\\ 0 & = & 2\cdot x_1^2 - x_4 - 1\\ 0 & = & \exp(x_4) - x_5 \\ 0 & = & x_0^3 + 3\cdot x_0^2\cdot x_3 + 3\cdot x_0\cdot x_3^2 + x_3^3 - x_6  \end{array}\right.
\f}
Finally, the original vector-valued function \f${\bf f}(x_0,x_1)\f$ is equal to \f$(x_6,x_5)^{\sf T}\f$.
*/

// TO DO:
// - Complete documentation

#ifndef MC__SELIM_H
#define MC__SELIM_H

#include "ffunc.hpp"
#include "ffinv.hpp"

#define MC__SELIM_CHECK

namespace mc
{

template <typename DAG> class SElimExpr;

//! @brief Environment for variable elimination in factorable equality constraints
////////////////////////////////////////////////////////////////////////
//! mc::SElimEnv is a C++ class defining the environment for variable
//! elimination from factorable expressions through expoiting solvable
//! equality constraints
////////////////////////////////////////////////////////////////////////
template < typename DAG >
class SElimEnv
////////////////////////////////////////////////////////////////////////
{
  friend class SElimExpr<DAG>;
  template <typename D> friend  std::ostream& operator<< ( std::ostream&, SElimEnv<D> const& );
  template <typename D> friend  SElimExpr<D> inv ( SElimExpr<D> const& );
  template <typename D> friend  SElimExpr<D> exp ( SElimExpr<D> const& );
  template <typename D> friend  SElimExpr<D> log ( SElimExpr<D> const& );
  template <typename D> friend  SElimExpr<D> xlog( SElimExpr<D> const& );
  template <typename D> friend  SElimExpr<D> sqrt( SElimExpr<D> const& );
  template <typename D> friend  SElimExpr<D> sqr ( SElimExpr<D> const& );
  template <typename D> friend  SElimExpr<D> pow ( SElimExpr<D> const&, int const );  
  template <typename D> friend  SElimExpr<D> pow ( SElimExpr<D> const&, double const& );  
  template <typename D> friend  SElimExpr<D> cheb( SElimExpr<D> const&, const unsigned );  
  template <typename D> friend  SElimExpr<D> prod( const unsigned, SElimExpr<D> const* );  
  template <typename D> friend  SElimExpr<D> cos ( SElimExpr<D> const& );
  template <typename D> friend  SElimExpr<D> sin ( SElimExpr<D> const& );
  template <typename D> friend  SElimExpr<D> tan ( SElimExpr<D> const& );
  template <typename D> friend  SElimExpr<D> acos( SElimExpr<D> const& );
  template <typename D> friend  SElimExpr<D> asin( SElimExpr<D> const& );
  template <typename D> friend  SElimExpr<D> atan( SElimExpr<D> const& );
  template <typename D> friend  SElimExpr<D> cosh( SElimExpr<D> const& );
  template <typename D> friend  SElimExpr<D> sinh( SElimExpr<D> const& );
  template <typename D> friend  SElimExpr<D> tanh( SElimExpr<D> const& );
  template <typename D> friend  SElimExpr<D> fabs( SElimExpr<D> const& );
  template <typename D> friend  SElimExpr<D> erf( SElimExpr<D> const& );
  template <typename D> friend  SElimExpr<D> fstep( SElimExpr<D> const& );
  template <typename D> friend  SElimExpr<D> max( SElimExpr<D> const&, SElimExpr<D> const& );  
  template <typename D> friend  SElimExpr<D> min( SElimExpr<D> const&, SElimExpr<D> const& );  
  template <typename D> friend  SElimExpr<D> lmtd( SElimExpr<D> const&, SElimExpr<D> const& );  
  template <typename D> friend  SElimExpr<D> rlmtd( SElimExpr<D> const&, SElimExpr<D> const& );  

public:

  typedef std::list< std::pair< FFOp const*, std::vector<SElimExpr<DAG> const*> > > t_Interm;
  typedef std::map< FFVar const*, FFVar const*, lt_FFVar > t_Aux;
  typedef std::vector< FFVar const* > t_Expr;
  typedef SPoly<FFVar const*, lt_FFVar> SPolyExpr;

  //! @brief Default Constructor
  SElimEnv
    ( DAG* dag=nullptr )
    : _dag( dag )
    {}

  //! @brief Destructor
  virtual ~SElimEnv
    ()
    { _reset(); }
  
  // Retreive pointer to DAG
  DAG* dag
    ()
    const
    { return _dag; };

  //! @brief Retreive reference to vector of new DAG rational constraints
  std::vector<FFVar>& Poly
    ()
    { return _Poly; }

  //! @brief Retreive reference to vector of new DAG transcendental constraints
  std::vector<FFVar>& Trans
    ()
    { return _Trans; }

  //! @brief Retreive reference to vector of independent DAG variables participating in expressions
  std::vector<FFVar>& Var
    ()
    { return _Var; }

  //! @brief Retreive reference to mapping between existing DAG auxiliaries and new DAG variables
  t_Aux& Aux
    ()
    { return _Aux; }

  //! @brief Retreive reference to intermediate expressions
  t_Interm& Interm
    ()
    { return _Interm; }

  //! @brief Set DAG environment
  void set
    (  DAG* dag )
    { _dag = dag; _reset(); }

  //! @brief Reset intermediate expressions
  void reset
    ()
    { _reset(); }

  //! @brief Process the dependents in set <a>sDep</a>
  void process
    ( std::set<unsigned> const& ndxDep, FFVar const* pDep, bool const add2dag=true );

  //! @brief Process the <a>nDep</a> dependents in array <a>pDep</a>
  void process
    ( unsigned const nDep, FFVar const* pDep, const bool add2dag=true );

  //! @brief Exceptions of mc::SElimExpr
  class Exceptions
  {
   public:
    //! @brief Enumeration type for SElimExpr exception handling
    enum TYPE{
      DAGERR=0,       //!< Operation involving a factorable expression linked to a different DAG
      ENVERR,         //!< Operation between factorable expressions linked to different environments or without an environment
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
        return "mc::SElimEnv\t Operation involving a factorable expression linked to a different DAG is not allowed";
      case ENVERR:
        return "mc::SElimEnv\t Operation between factorable expressions linked to different environments or without an environment is not allowed";
      case INTERNAL:
      default:
        return "mc::SElimEnv\t Internal error";
      }
    }
   private:
    TYPE _ierr;
  };

  //! @brief Options of mc::SElimEnv
  static struct Options
  {
    //! @brief Constructor
    Options():
      LIFTDIV( true ), LIFTIPOW( false )//, LIFTUPOL( false )
      {}
    //! @brief Assignment of mc::SElimEnv<DAG>::Options
    Options& operator=
      ( Options& opt ){
        LIFTDIV  = opt.LIFTDIV;
        LIFTIPOW = opt.LIFTIPOW;
        //LIFTUPOL = opt.LIFTUPOL;
        return *this;
      }
    //! @brief Whether to lift division terms using auxiliary variables (default: true)
    bool LIFTDIV;
    //! @brief Whether to lift integral power terms using auxiliary variables (default: true)
    bool LIFTIPOW;
    //! @brief Whether to lift univariate polynomials using auxiliary variables (default: true)
    //bool LIFTUPOL;
  } options;

protected:
  //! @brief pointer to underlying dag
  DAG* _dag;

  //! @brief Map of DAG variables and sparse expressions
  t_Interm _Interm;

  //! @brief Map of existing DAG auxiliaries to new DAG variables
  t_Aux _Aux;

  //! @brief Vector of new DAG polynomial/rational constraints
  std::vector<FFVar> _Poly;

  //! @brief Vector of new DAG transcendental constraints
  std::vector<FFVar> _Trans;

  //! @brief Vector of independent DAG variables participating in expressions
  std::vector<FFVar> _Var;

  //! @brief Vector of independent sparse variables
  std::vector<SElimExpr<DAG>> _SPVar;

  //! @brief Vector of dependent sparse expressions
  std::vector<SElimExpr<DAG>> _SPDep;

  //! @brief Add new intermediate uni- or bi-variate expression in _Interm
  void _append_interm
    ( FFOp const* op, SElimExpr<DAG> const* var1, SElimExpr<DAG> const* var2=nullptr );

//  //! @brief Add new intermediate n-variate expression in _Interm
//  void _append_interm
//    ( FFOp const* op, unsigned const nvars, SElimExpr<DAG> const*const* pvars );

  //! @brief Insert an auxiliary variable corresponding to a rational/polynomial expression into DAG
  FFVar const* _insert_expr
    ( FFVar const* oper, SElimExpr<DAG> const* expr );

  //! @brief Insert a non-rational operation into DAG via the introduction of auxiliary variables
  void _insert_expr
    ( FFOp const* pOp, std::vector<const FFVar*>& vAux );

  //! @brief Transcribe sparse rational/polynomial expression (SElimExpr) into DAG
  FFVar _SPolyExpr_to_FFVar
    ( SPolyExpr const& expr );

  //! @brief Find auxiliary variable in DAG and return corresponding new DAG variable (or NULL if undefined)
  FFVar const* _find_aux
    ( FFVar const* aux );

  //! @brief Lift current univariate operation
  SElimExpr<DAG> _lift_univariate_term
    ( SElimExpr<DAG> const& var );

  //! @brief Lift current bivariate operation
  SElimExpr<DAG> _lift_bivariate_term
    ( SElimExpr<DAG> const& var1, SElimExpr<DAG> const& var2 );

  //! @brief Erase all entries in _Interm
  void _reset
    ();
};


template <typename DAG>
inline typename SElimEnv<DAG>::Options SElimEnv<DAG>::options;

//! @brief Arithmetic for variable elimination in factorable equality constraints
////////////////////////////////////////////////////////////////////////
//! mc::SElimExpr is a C++ class implementing an arithmetic for variable
//! elimination from factorable expressions through expoiting solvable
//! equality constraints
////////////////////////////////////////////////////////////////////////
template < typename DAG >
class SElimExpr
////////////////////////////////////////////////////////////////////////
{
public:

  typedef SPoly<FFVar const*, lt_FFVar> SPolyExpr;

private:

  //! @brief Pointer to recursive decomposition environment
  SElimEnv<DAG> *_env;

protected:

  // numerator sparse polynomial
  SPolyExpr _numer;

  // Denominator sparse polynomial
  SPolyExpr _denom;

  //! @brief Initialize sparse rational expression with existing expression
  SElimExpr<DAG>& _set
    ( SElimExpr<DAG> const& var );

  //! @brief Initialize sparse expression as constant
  SElimExpr<DAG>& _set
    ( double const& d );

  //! @brief Initialize sparse rational expression as DAG variable
  SElimExpr<DAG>& _set
    ( FFVar const& x );

public:

  //! @brief Default constructor of sparse rational expression
  SElimExpr
    ()
    {}

  //! @brief Constructor of sparse rational expression as constant
  SElimExpr
    ( double const& d )
    : _env( 0 )
    { _set( d ); }

  //! @brief Constructor of sparse rational expression as DAG variable
  SElimExpr
    ( SElimEnv<DAG>* env, FFVar const& x)
    : _env( env )
    { if( env->dag() != x.dag() )
        throw typename SElimEnv<DAG>::Exceptions( SElimEnv<DAG>::Exceptions::DAGERR );
      _set( x ); }

  //! @brief Copy constructor of sparse rational expression
  SElimExpr
    ( SElimExpr<DAG> const& var )
    { _set( var ); }

  //! @brief Constructor of sparse rational expression
  SElimExpr
    ( SElimEnv<DAG>* env, SPolyExpr const& n, SPolyExpr const& d )
    : _env( env ), _numer( n ), _denom( d )
    {}

  //! @brief Destructor of sparse rational expression
  virtual ~SElimExpr()
    {}

  //! @brief Initialize variable in sparse rational envrionment <a>env</a> corresponding to DAG variable <a>x</a>
  SElimExpr<DAG>& set
    ( SElimEnv<DAG>* env, FFVar const& x )
    { _env = env;
      if( env->dag() != x.dag() )
        throw typename SElimEnv<DAG>::Exceptions( SElimEnv<DAG>::Exceptions::DAGERR );
      _set( x ); return *this; }

  //! @brief Overloaded operator '=' for sparse rational expression
  SElimExpr<DAG>& operator=
    ( SElimExpr<DAG> const& var )
    { _set( var ); return *this; }

  //! @brief Overloaded operator '=' for constant 
  SElimExpr<DAG>& operator=
    ( double const& d )
    { _env = 0; _set( d ); return *this; }

  //! @brief Overloaded operator '+=' for sparse rational function
  SElimExpr<DAG>& operator+=
    ( SElimExpr<DAG> const& var );

  //! @brief Overloaded operator '-=' for sparse rational function
  SElimExpr<DAG>& operator-=
    ( SElimExpr<DAG> const& var );

  //! @brief Overloaded operator '*=' for sparse rational function
  SElimExpr<DAG>& operator*=
    ( SElimExpr<DAG> const& var );

  //! @brief Overloaded operator '/=' for sparse rational function
  SElimExpr<DAG>& operator/=
    ( SElimExpr<DAG> const& var );

  // Sparse rational polynomial environment
  SElimEnv<DAG>* env
    ()
    const
    { return _env; };

  // numerator sparse polynomial
  SPolyExpr const& numer
    ()
    const
    { return _numer; };

  // Denominator sparse polynomial
  SPolyExpr const& denom
    ()
    const
    { return _denom; };
};

////////////////////////////////////////////////////////////////////////

template < typename DAG >
inline std::ostream&
operator<<
( std::ostream& out, SElimExpr<DAG> const& var )
{
  out << std::endl
      << "NUMERATOR:"   << var.numer()
      << "DENOMINATOR:" << var.denom();
  return out;
}

template < typename DAG >
inline std::ostream&
operator<<
( std::ostream& out, SElimEnv<DAG> const& env)
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

template < typename DAG >
inline void
SElimEnv<DAG>::process
( std::set<unsigned> const& ndxDep, FFVar const* pDep, bool const add2dag )
{
  if( ndxDep.empty() ) return; // Nothing to do!
  std::vector<FFVar> vpDep;//( sDep.begin(), sDep.end() );
  vpDep.reserve( ndxDep.size() );
  for( unsigned const& i : ndxDep ) vpDep.push_back( pDep[i] );
  process( ndxDep.size(), vpDep.data(), add2dag );
}

template < typename DAG >
inline void
SElimEnv<DAG>::process
( unsigned const nDep, FFVar const* pDep, bool const add2dag )
{
  // Reset intermediate / auxiliary arrays
  _reset();

  // Update participating variables in _Var
  auto sgDep = _dag->subgraph( nDep, pDep );
  for( auto&& Op : sgDep.l_op ){
    if( Op->type != FFOp::VAR ) continue;
    _Var.push_back( *Op->pres );
    _SPVar.push_back( SElimExpr( this, *Op->pres ) );
  }

#ifdef MC__SELIM_DEBUG_PROCESS
  std::cout << std::endl << _Var.size() << " Original Variables: ";
  for( auto&& var : _Var ) std::cout << var << " ";
  std::cout << std::endl;
#endif

  // Process DAG dependents
  _SPDep.resize( nDep );
  _dag->eval( sgDep, nDep, pDep, _SPDep.data(), _Var.size(), _Var.data(), _SPVar.data() );

#ifdef MC__SELIM_DEBUG_PROCESS
  std::cout << *this;
  for( unsigned i=0; i<nDep; i++ )
    std::cout << std::endl << "SPDep[" << i << "]:" << _SPDep[i];
#endif

  // No transcription in DAG if <a>add2dag</a> is false
  if( !add2dag ) return;

  // Insert auxiliary expressions into DAG
  for( auto&& expr : _Interm ){
    // Insert all operands and their defining expressions
    std::vector<FFVar const*> vAux;
    auto itSV = expr.second.begin();
    for( auto&& operand : expr.first->pops ){
#ifdef MC__SELIM_CHECK
      if( itSV == expr.second.end() )
        throw typename SElimEnv<DAG>::Exceptions( SElimEnv<DAG>::Exceptions::INTERNAL );
#endif
      vAux.push_back( _insert_expr( operand, *itSV ) );
      ++itSV;
    }
    // Insert operation result and their defining expressions
    _insert_expr( expr.first, vAux );
  }

  // Insert terminal expressions into DAG
  for( unsigned i=0; i<nDep; i++ )
    // Insert operation result and defining expression
//    _Dep.push_back( *_insert_expr( pDep+i, &_SPDep.at(i) ) );
    _insert_expr( pDep+i, &_SPDep.at(i) );

#ifdef MC__SELIM_DEBUG_PROCESS
  std::cout << std::endl << _Aux.size() << " Auxiliary Variables: ";
  for( auto&& aux : _Aux ) std::cout << *aux.first << "->" << *aux.second << " ";
  std::cout << std::endl;
//  std::cout << std::endl << _Dep.size() << " Dependent Variables: ";
//  for( auto&& dep : _Dep ) std::cout << dep << " ";
//  std::cout << std::endl;
  std::cout << std::endl << _Poly.size() << " Polynomial Constraints: " << std::endl;
  for( auto&& expr : _Poly ) _dag->output( _dag->subgraph( 1, &expr ) );
  std::cout << std::endl;
  std::cout << std::endl << _Trans.size() << " Transcendental Constraints: " << std::endl;
  for( auto&& expr : _Trans ) _dag->output( _dag->subgraph( 1, &expr ) );
#endif
}

template < typename DAG >
inline FFVar const*
SElimEnv<DAG>::_insert_expr
( FFVar const* var, SElimExpr<DAG> const* expr )
{ 
  auto itdagvar = _dag->Vars().find( const_cast<FFVar*>(var) );
#ifdef MC__SELIM_CHECK
  assert( itdagvar != _dag->Vars().end() );
#endif

  // Nothing to do if operand is a DAG constant or leaf variable
  if( var->ops().first->type == FFOp::VAR || var->ops().first->type == FFOp::CNST )
    return *itdagvar;

  // Nothing to do if DAG auxiliary was already made a DAG variable
  auto itv = _Aux.find( *itdagvar );
  if( itv != _Aux.end() )
    return itv->second;

  // Append new DAG variable in _Aux and defining polynomial constraint in _Poly 
#ifdef MC__SELIM_DEBUG_PROCESS
  std::cout <<std::endl << "operand: " << **itdagvar << std::endl;
#endif
  FFVar newvar( _dag );
  auto itnewvar = _dag->Vars().find( &newvar );
#ifdef MC__SELIM_CHECK
  assert( itnewvar != _dag->Vars().end() );
#endif
#ifdef MC__SELIM_DEBUG_PROCESS
  std::cout << "paired with new DAG variable: " << **itnewvar << std::endl;
#endif
  _Aux.insert( std::make_pair( *itdagvar, *itnewvar ) );
  _Var.push_back( **itnewvar );

  FFVar polyctr = **itnewvar * _SPolyExpr_to_FFVar( expr->denom() ) - _SPolyExpr_to_FFVar( expr->numer() );
  auto itpolyctr = _dag->Vars().find( &polyctr );
#ifdef MC__SELIM_CHECK
  assert( itpolyctr != _dag->Vars().end() );
#endif
#ifdef MC__SELIM_DEBUG_PROCESS
  std::cout << "defined by DAG subexpression: ";
  _dag->output( _dag->subgraph( 1, *itpolyctr ) );
#endif
  _Poly.push_back( **itpolyctr );
  return *itnewvar;
}

template < typename DAG >
inline void
SElimEnv<DAG>::_insert_expr
( FFOp const* pOp, std::vector<FFVar const*>& vAux )
{ 
#ifdef MC__SELIM_CHECK
  // Throw exception if DAG auxiliary was already made a DAG variable
  if( _Aux.find( pOp->pres ) != _Aux.end() )
    throw typename SElimEnv<DAG>::Exceptions( SElimEnv<DAG>::Exceptions::INTERNAL );
#endif
#ifdef MC__SELIM_DEBUG_PROCESS
  std::cout << std::endl << "operand: " << *pOp->pres << std::endl;
#endif

  // Append new DAG variable in _Aux
  FFVar newvar( _dag );
  auto itnewvar = _dag->Vars().find( &newvar );
#ifdef MC__SELIM_CHECK
  assert( itnewvar != _dag->Vars().end() );
#endif
#ifdef MC__SELIM_DEBUG_PROCESS
  std::cout << "paired with new DAG variable: " << **itnewvar << std::endl;
#endif
  _Aux.insert( std::make_pair( pOp->pres, *itnewvar ) );
  _Var.push_back( **itnewvar );


  // Append defining expression in _Poly or _Transc 
  switch( pOp->type ){
   case FFOp::SQR:   return _Poly.push_back( **itnewvar - sqr( *vAux.at(0) ) );
   case FFOp::IPOW:  return _Poly.push_back( **itnewvar - pow( *vAux.at(0), vAux.at(1)->num().n ) );
   case FFOp::CHEB:  return _Poly.push_back( **itnewvar - cheb( *vAux.at(0), vAux.at(1)->num().n ) );
   case FFOp::SQRT:  return _Poly.push_back( sqr( **itnewvar ) - *vAux.at(0) );
   case FFOp::INV:   return _Poly.push_back( **itnewvar * *vAux.at(1) - vAux.at(0)->num().val() );
   case FFOp::DIV:   return _Poly.push_back( **itnewvar * *vAux.at(1) - *vAux.at(0) );
   case FFOp::EXP:   return _Trans.push_back( **itnewvar - exp( *vAux.at(0) ) );
   case FFOp::LOG:   return _Trans.push_back( exp( **itnewvar ) - *vAux.at(0) );
   case FFOp::XLOG:  return _Trans.push_back( **itnewvar - xlog( *vAux.at(0) ) );
   case FFOp::DPOW:  return _Trans.push_back( **itnewvar - pow( *vAux.at(0), vAux.at(1)->num().val() ) );
   case FFOp::COS:   return _Trans.push_back( **itnewvar - cos( *vAux.at(0) ) );
   case FFOp::SIN:   return _Trans.push_back( **itnewvar - sin( *vAux.at(0) ) );
   case FFOp::TAN:   return _Trans.push_back( **itnewvar - tan( *vAux.at(0) ) );
   case FFOp::ACOS:  return _Trans.push_back( cos( **itnewvar ) - *vAux.at(0) );
   case FFOp::ASIN:  return _Trans.push_back( sin( **itnewvar ) - *vAux.at(0) );
   case FFOp::ATAN:  return _Trans.push_back( tan( **itnewvar ) - *vAux.at(0) );
   case FFOp::COSH:  return _Trans.push_back( **itnewvar - cosh( *vAux.at(0) ) );
   case FFOp::SINH:  return _Trans.push_back( **itnewvar - sinh( *vAux.at(0) ) );
   case FFOp::TANH:  return _Trans.push_back( **itnewvar - tanh( *vAux.at(0) ) );
   case FFOp::ERF:   return _Trans.push_back( **itnewvar - erf( *vAux.at(0) ) );
   case FFOp::FABS:  return _Trans.push_back( **itnewvar - fabs( *vAux.at(0) ) );
   case FFOp::FSTEP: return _Trans.push_back( **itnewvar - fstep( *vAux.at(0) ) );
   case FFOp::MINF:  return _Trans.push_back( **itnewvar - min( *vAux.at(0), *vAux.at(1) ) );
   case FFOp::MAXF:  return _Trans.push_back( **itnewvar - max( *vAux.at(0), *vAux.at(1) ) );
   case FFOp::INTER: return _Trans.push_back( **itnewvar - inter( *vAux.at(0), *vAux.at(1) ) );
   case FFOp::VAR:
   case FFOp::CNST:
   case FFOp::SHIFT:
   case FFOp::PLUS:
   case FFOp::NEG:
   case FFOp::MINUS:
   case FFOp::SCALE:
   case FFOp::TIMES:
   case FFOp::PROD:
   default:          throw typename SElimEnv<DAG>::Exceptions( SElimEnv<DAG>::Exceptions::INTERNAL );
  }
}

template < typename DAG >
inline FFVar const*
SElimEnv<DAG>::_find_aux
( FFVar const* aux )
{
  if( aux->ops().first->type == FFOp::VAR ) return aux;
  auto it = _Aux.find( aux );
  if( it != _Aux.end() ) return it->second;
  return (FFVar*)0;
}

template < typename DAG >
inline FFVar
SElimEnv<DAG>::_SPolyExpr_to_FFVar
( SPolyExpr const& expr )
{
  FFVar var = 0.;
  for( auto it=expr.mapmon().begin(); it!=expr.mapmon().end(); ++it ){
    FFVar prodmon = it->second;
    for( auto ie=it->first.expr.begin(); ie!=it->first.expr.end(); ++ie ){
      const FFVar*oper = _find_aux( ie->first );
      if( !oper ) 
        throw typename SElimEnv<DAG>::Exceptions( SElimEnv<DAG>::Exceptions::INTERNAL );
      switch( SPolyExpr::options.BASIS ){
       case SPolyExpr::Options::MONOM:
        prodmon *= pow( *oper, (int)ie->second );
        break;
       case SPolyExpr::Options::CHEB:
        prodmon *= cheb( *oper, ie->second );
        break;
      }
    }
    if( it->first.expr.empty() ) var = prodmon;
    else                         var+= prodmon;
  }
  // ADD OPTION TO RETURN A PRODMON OR MONOM TERM?
  return var;
}

template < typename DAG >
inline void
SElimEnv<DAG>::_reset
()
{
  for( auto&& expr : _Interm )
    for( auto&& oper : expr.second )
      delete oper;
  _Interm.clear();

  _Aux.clear();
  _Poly.clear();
  _Trans.clear();
  _Var.clear();
//  _Dep.clear();
  _SPVar.clear();
  _SPDep.clear();
}

template < typename DAG >
inline void
SElimEnv<DAG>::_append_interm
( FFOp const* op, SElimExpr<DAG> const* var1, SElimExpr<DAG> const* var2 )
{
#ifdef MC__SELIM_CHECK
  if( !var1 )
    throw typename SElimEnv<DAG>::Exceptions( SElimEnv<DAG>::Exceptions::INTERNAL );
#endif
  std::vector<SElimExpr<DAG> const*> vops;
  vops.push_back( new SElimExpr( *var1 ) );
  if( var2 ) vops.push_back( new SElimExpr( *var2 ) );
  _Interm.push_back( std::make_pair( op, vops ) );
}

template < typename DAG >
inline SElimExpr<DAG>
SElimEnv<DAG>::_lift_univariate_term
( SElimExpr<DAG> const& var )
{
#ifdef MC__SELIM_CHECK
  if( !_dag->curOp() )
    throw typename SElimEnv<DAG>::Exceptions( SElimEnv<DAG>::Exceptions::INTERNAL );
#endif

  // Append new intermediate expression and assert that same operation was not previously appended
  _append_interm( _dag->curOp(), &var );
  return SElimExpr( this, *(_dag->curOp()->pres) );
}

template < typename DAG >
inline SElimExpr<DAG>
SElimEnv<DAG>::_lift_bivariate_term
( SElimExpr<DAG> const& var1, SElimExpr<DAG> const& var2 )
{
#ifdef MC__SELIM_CHECK
  if( !_dag->curOp() )
    throw typename SElimEnv<DAG>::Exceptions( SElimEnv<DAG>::Exceptions::INTERNAL );
#endif

  // Append new intermediate expression and assert that same operation was not previously appended
  _append_interm( _dag->curOp(), &var1, &var2 );
  return SElimExpr( this, *(_dag->curOp()->pres) );
}

template < typename DAG >
inline SElimExpr<DAG>&
SElimExpr<DAG>::_set
( double const& d )
{
  _numer = d;
  _denom = 1;
  return *this;
}

template < typename DAG >
inline SElimExpr<DAG>&
SElimExpr<DAG>::_set
( FFVar const& x )
{
  _numer.var( &x );
  _denom = 1;
  return *this;
}

template < typename DAG >
inline SElimExpr<DAG>&
SElimExpr<DAG>::_set
( SElimExpr<DAG> const& var )
{
  if( this == &var ) return *this;
  _env = var._env;
  _numer = var._numer;
  _denom = var._denom;
  return *this;
}

template < typename DAG >
inline SElimExpr<DAG>
operator+
( SElimExpr<DAG> const& var )
{
  return var;
}

template < typename DAG >
inline SElimExpr<DAG>&
SElimExpr<DAG>::operator+=
( SElimExpr<DAG> const& var )
{
  if( !_env )
    _env = var._env;
  else if( var._env && _env != var._env )
    throw typename SElimEnv<DAG>::Exceptions( SElimEnv<DAG>::Exceptions::ENVERR );

  _numer *= var._denom;
  _numer += var._numer * _denom;
  _denom *= var._denom;
  return *this;
}

template < typename DAG >
inline SElimExpr<DAG>
operator+
( SElimExpr<DAG> const& var1, SElimExpr<DAG> const& var2 )
{
  SElimExpr<DAG> var3( var1 );
  var3 += var2;
  return var3;
}

template < typename DAG >
inline SElimExpr<DAG>
operator+
( SElimExpr<DAG> const& var1, double const& cst2 )
{
  SElimExpr<DAG> var3( var1 );
  var3 += cst2;
  return var3;
}

template < typename DAG >
inline SElimExpr<DAG>
operator-
( SElimExpr<DAG> const& var )
{
  return SElimExpr( var.env(), -var.numer(), var.denom() );
}

template < typename DAG >
inline SElimExpr<DAG>&
SElimExpr<DAG>::operator-=
( SElimExpr<DAG> const& var )
{
  if( !_env )
    _env = var._env;
  else if( var._env && _env != var._env )
    throw typename SElimEnv<DAG>::Exceptions( SElimEnv<DAG>::Exceptions::ENVERR );

  _numer *= var._denom;
  _numer -= var._numer * _denom;
  _denom *= var._denom;
  return *this;
}

template < typename DAG >
inline SElimExpr<DAG>
operator-
( SElimExpr<DAG> const& var1, SElimExpr<DAG> const& var2 )
{
  SElimExpr<DAG> var3( var1 );
  var3 -= var2;
  return var3;
}

template < typename DAG >
inline SElimExpr<DAG>
operator-
( SElimExpr<DAG> const& var1, double const& cst2 )
{
  SElimExpr<DAG> var3( var1 );
  var3 -= cst2;
  return var3;
}

template < typename DAG >
inline SElimExpr<DAG>&
SElimExpr<DAG>::operator*=
( SElimExpr<DAG> const& var )
{
  if( !_env )
    _env = var._env;
  else if( var._env && _env != var._env )
    throw typename SElimEnv<DAG>::Exceptions( SElimEnv<DAG>::Exceptions::ENVERR );

  _numer *= var._numer;
  _denom *= var._denom;
  return *this;
}

template < typename DAG >
inline SElimExpr<DAG>
operator*
( SElimExpr<DAG> const& var1, SElimExpr<DAG> const& var2 )
{
  SElimExpr<DAG> var3( var1 );
  var3 *= var2;
  return var3;
}

template < typename DAG >
inline SElimExpr<DAG>
operator*
( SElimExpr<DAG> const& var1, double const& cst2 )
{
  SElimExpr<DAG> var3( var1 );
  var3 *= cst2;
  return var3;
}

template < typename DAG >
inline SElimExpr<DAG>&
SElimExpr<DAG>::operator/=
( SElimExpr<DAG> const& var )
{
  if( !_env )
    _env = var._env;
  else if( var._env && _env != var._env )
    throw typename SElimEnv<DAG>::Exceptions( SElimEnv<DAG>::Exceptions::ENVERR );

  if( SElimEnv<DAG>::options.LIFTDIV )
    *this = _env->_lift_bivariate_term( *this, var );
  else{
    _numer *= var._denom;
    _denom *= var._numer;
  }
  return *this;
}

template < typename DAG >
inline SElimExpr<DAG>
operator/
( SElimExpr<DAG> const& var1, SElimExpr<DAG> const& var2 )
{
  SElimExpr<DAG> var3( var1 );
  var3 /= var2;
  return var3;
}

template < typename DAG >
inline SElimExpr<DAG>
operator/
( double const& cst1, SElimExpr<DAG> const& var2 )
{
  SElimExpr<DAG> var3( cst1 );
  var3 /= var2;
  return var3;
}

template < typename DAG >
inline SElimExpr<DAG>
sqr
( SElimExpr<DAG> const& var )
{
#ifdef MC__SELIM_CHECK
  if( !var.env() )
    throw typename SElimEnv<DAG>::Exceptions( SElimEnv<DAG>::Exceptions::ENVERR );
#endif
  if( SElimEnv<DAG>::options.LIFTIPOW )
    return var.env()->_lift_univariate_term( var );
  return SElimExpr<DAG>( var.env(), sqr( var.numer() ), sqr( var.denom() ) );
}

template < typename DAG >
inline SElimExpr<DAG>
inv
( SElimExpr<DAG> const& var )
{
#ifdef MC__SELIM_CHECK
  if( !var.env() )
    throw typename SElimEnv<DAG>::Exceptions( SElimEnv<DAG>::Exceptions::ENVERR );
#endif
  if( SElimEnv<DAG>::options.LIFTDIV )
    return var.env()->_lift_bivariate_term( 1, var );
  return SElimExpr<DAG>( var.env(), var.denom(), var.numer() );
}

template < typename DAG >
inline SElimExpr<DAG>
pow
( SElimExpr<DAG> const& var, double const& a )
{
#ifdef MC__SELIM_CHECK
  if( !var.env() )
    throw typename SElimEnv<DAG>::Exceptions( SElimEnv<DAG>::Exceptions::ENVERR );
#endif
  return var.env()->_lift_bivariate_term( var, a );
}

template < typename DAG >
inline SElimExpr<DAG>
pow
( SElimExpr<DAG> const& var, int const n )
{
  if( n < 0 ) return pow( inv( var ), -n );
  switch( n ){
   case 0:  return 1.;
   case 1:  return var;
   case 2:  return sqr( var );
   default: if( SElimEnv<DAG>::options.LIFTIPOW ) return var.env()->_lift_bivariate_term( var, n );
            return SElimExpr( var.env(), pow( var.numer(), (unsigned)n ), pow( var.denom(), (unsigned)n ) );
  }
}

template < typename DAG >
inline SElimExpr<DAG>
cheb
( SElimExpr<DAG> const& var, unsigned const n )
{
  switch( n ){
   case 0:  return 1.;
   case 1:  return var;
   case 2:  return sqr( var ) * 2. - 1.;
   default: if( SElimEnv<DAG>::options.LIFTIPOW ) return var.env()->_lift_bivariate_term( var, n );
            return var * cheb( var, n-1 ) * 2. - cheb( var, n-2 );
  }
}

template < typename DAG >
inline SElimExpr<DAG>
prod
( unsigned int const nvars, SElimExpr<DAG> const* pvars )
{
  switch( nvars ){
   case 0:  return 1.;
   case 1:  return pvars[0];
   default: return pvars[0] * prod( nvars-1, pvars+1 );
  }
}

template < typename DAG >
inline SElimExpr<DAG>
monom
( unsigned int const nvars, SElimExpr<DAG> const* pvars, unsigned const* k, bool const chebbasis=false )
{
  switch( nvars ){
   case 0:  return 1.;
   case 1:  return chebbasis? cheb( pvars[0], k[0] ): pow( pvars[0], (int)k[0] );
   default: return ( chebbasis? cheb( pvars[0], k[0] ): pow( pvars[0], (int)k[0] ) ) * monom( nvars-1, pvars+1, k+1 );
  }
}

template < typename DAG >
inline SElimExpr<DAG>
exp
( SElimExpr<DAG> const& var )
{
#ifdef MC__SELIM_CHECK
  if( !var.env() )
    throw typename SElimEnv<DAG>::Exceptions( SElimEnv<DAG>::Exceptions::ENVERR );
#endif
  return var.env()->_lift_univariate_term( var );
}

template < typename DAG >
inline SElimExpr<DAG>
log
( SElimExpr<DAG> const& var )
{
#ifdef MC__SELIM_CHECK
  if( !var.env() )
    throw typename SElimEnv<DAG>::Exceptions( SElimEnv<DAG>::Exceptions::ENVERR );
#endif
  return var.env()->_lift_univariate_term( var );
}

template < typename DAG >
inline SElimExpr<DAG>
xlog
( SElimExpr<DAG> const& var )
{
#ifdef MC__SELIM_CHECK
  if( !var.env() )
    throw typename SElimEnv<DAG>::Exceptions( SElimEnv<DAG>::Exceptions::ENVERR );
#endif
  return var.env()->_lift_univariate_term( var );
}

template < typename DAG >
inline SElimExpr<DAG>
sqrt
( SElimExpr<DAG> const& var )
{
#ifdef MC__SELIM_CHECK
  if( !var.env() )
    throw typename SElimEnv<DAG>::Exceptions( SElimEnv<DAG>::Exceptions::ENVERR );
#endif
  return var.env()->_lift_univariate_term( var );
}

template < typename DAG >
inline SElimExpr<DAG>
cos
( SElimExpr<DAG> const& var )
{
#ifdef MC__SELIM_CHECK
  if( !var.env() )
    throw typename SElimEnv<DAG>::Exceptions( SElimEnv<DAG>::Exceptions::ENVERR );
#endif
  return var.env()->_lift_univariate_term( var );
}

template < typename DAG >
inline SElimExpr<DAG>
sin
( SElimExpr<DAG> const& var )
{
#ifdef MC__SELIM_CHECK
  if( !var.env() )
    throw typename SElimEnv<DAG>::Exceptions( SElimEnv<DAG>::Exceptions::ENVERR );
#endif
  return var.env()->_lift_univariate_term( var );
}

template < typename DAG >
inline SElimExpr<DAG>
tan
( SElimExpr<DAG> const& var )
{
#ifdef MC__SELIM_CHECK
  if( !var.env() )
    throw typename SElimEnv<DAG>::Exceptions( SElimEnv<DAG>::Exceptions::ENVERR );
#endif
  return var.env()->_lift_univariate_term( var );
}

template < typename DAG >
inline SElimExpr<DAG>
acos
( SElimExpr<DAG> const& var )
{
#ifdef MC__SELIM_CHECK
  if( !var.env() )
    throw typename SElimEnv<DAG>::Exceptions( SElimEnv<DAG>::Exceptions::ENVERR );
#endif
  return var.env()->_lift_univariate_term( var );
}

template < typename DAG >
inline SElimExpr<DAG>
asin
( SElimExpr<DAG> const& var )
{
#ifdef MC__SELIM_CHECK
  if( !var.env() )
    throw typename SElimEnv<DAG>::Exceptions( SElimEnv<DAG>::Exceptions::ENVERR );
#endif
  return var.env()->_lift_univariate_term( var );
}

template < typename DAG >
inline SElimExpr<DAG>
atan
( SElimExpr<DAG> const& var )
{
#ifdef MC__SELIM_CHECK
  if( !var.env() )
    throw typename SElimEnv<DAG>::Exceptions( SElimEnv<DAG>::Exceptions::ENVERR );
#endif
  return var.env()->_lift_univariate_term( var );
}

template < typename DAG >
inline SElimExpr<DAG>
cosh
( SElimExpr<DAG> const& var )
{
#ifdef MC__SELIM_CHECK
  if( !var.env() )
    throw typename SElimEnv<DAG>::Exceptions( SElimEnv<DAG>::Exceptions::ENVERR );
#endif
  return var.env()->_lift_univariate_term( var );
}

template < typename DAG >
inline SElimExpr<DAG>
sinh
( SElimExpr<DAG> const& var )
{
#ifdef MC__SELIM_CHECK
  if( !var.env() )
    throw typename SElimEnv<DAG>::Exceptions( SElimEnv<DAG>::Exceptions::ENVERR );
#endif
  return var.env()->_lift_univariate_term( var );
}

template < typename DAG >
inline SElimExpr<DAG>
tanh
( SElimExpr<DAG> const& var )
{
#ifdef MC__SELIM_CHECK
  if( !var.env() )
    throw typename SElimEnv<DAG>::Exceptions( SElimEnv<DAG>::Exceptions::ENVERR );
#endif
  return var.env()->_lift_univariate_term( var );
}

template < typename DAG >
inline SElimExpr<DAG>
fabs
( SElimExpr<DAG> const& var )
{
#ifdef MC__SELIM_CHECK
  if( !var.env() )
    throw typename SElimEnv<DAG>::Exceptions( SElimEnv<DAG>::Exceptions::ENVERR );
#endif
  return var.env()->_lift_univariate_term( var );
}

template < typename DAG >
inline SElimExpr<DAG>
erf
( SElimExpr<DAG> const& var )
{
#ifdef MC__SELIM_CHECK
  if( !var.env() )
    throw typename SElimEnv<DAG>::Exceptions( SElimEnv<DAG>::Exceptions::ENVERR );
#endif
  return var.env()->_lift_univariate_term( var );
}

template < typename DAG >
inline SElimExpr<DAG>
fstep
( SElimExpr<DAG> const& var )
{
#ifdef MC__SELIM_CHECK
  if( !var.env() )
    throw typename SElimEnv<DAG>::Exceptions( SElimEnv<DAG>::Exceptions::ENVERR );
#endif
  return var.env()->_lift_univariate_term( var );
}

template < typename DAG >
inline SElimExpr<DAG>
max
( SElimExpr<DAG> const& var1, SElimExpr<DAG> const& var2 )
{
  if( var1.env() && var2.env() && var1.env() != var2.env() )
    throw typename SElimEnv<DAG>::Exceptions( SElimEnv<DAG>::Exceptions::ENVERR );

  return var1.env()->_lift_bivariate_term( var1, var2 );
}

template < typename DAG >
inline SElimExpr<DAG>
min
( SElimExpr<DAG> const& var1, SElimExpr<DAG> const& var2 )
{
  if( var1.env() && var2.env() && var1.env() != var2.env() )
    throw typename SElimEnv<DAG>::Exceptions( SElimEnv<DAG>::Exceptions::ENVERR );

  return var1.env()->_lift_bivariate_term( var1, var2 );
}

template < typename DAG >
inline SElimExpr<DAG>
lmtd
( SElimExpr<DAG> const& var1, SElimExpr<DAG> const& var2 )
{
  if( var1.env() && var2.env() && var1.env() != var2.env() )
    throw typename SElimEnv<DAG>::Exceptions( SElimEnv<DAG>::Exceptions::ENVERR );

  return var1.env()->_lift_bivariate_term( var1, var2 );
}

template < typename DAG >
inline SElimExpr<DAG>
rlmtd
( SElimExpr<DAG> const& var1, SElimExpr<DAG> const& var2 )
{
  if( var1.env() && var2.env() && var1.env() != var2.env() )
    throw typename SElimEnv<DAG>::Exceptions( SElimEnv<DAG>::Exceptions::ENVERR );

  return var1.env()->_lift_bivariate_term( var1, var2 );
}

} // namespace mc

//#include "mcop.hpp"

namespace mc
{

//! @brief Specialization of the structure mc::Op to allow usage of the type mc::SElimExpr for DAG evaluation or as a template parameter in other MC++ classes
template <typename DAG> struct Op<mc::SElimExpr<DAG>>
{
  typedef mc::SElimExpr<DAG> T;
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
