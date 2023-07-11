// Copyright (C) Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the Eclipse Public License.

/*!
\page page_SLIFT Recursive Decomposition in Factorable Expressions 
\author Benoit Chachuat
\date 2023
\bug No known bugs.

The classes mc::SLiftEnv and mc::SLiftVar defined in <tt>slift.hpp</tt> enable the recursive decomposition of factorable expressions as a collection of polynomial/rational and transcendental subexpressions via the introduction of auxiliary variables.

\section sec_SLift_process How do I decompose a factorable expression?

For illustration, consider the factorable function \f${\bf f}:\mathbb{R}^2\to\mathbb{R}^2\f$ defined by
\f{align*}
  {\bf f}(x_0,x_1) = \left(\begin{array}{c} \left(x_0+\frac{1}{x_1^2}\right)^3\\ \exp\left(2\cdot x_1^2-1\right)\end{array}\right)
\f}

The lifting requires the header file <tt>slift.hpp</tt> to be included:

\code
      #include "slift.hpp"
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
      V0	 => { Z3 }
      V1	 => { Z0 }

    DAG INTERMEDIATES:
      Z0	<=  SQR( V1 )		 => { Z2 Z7 }
      Z2	<=  Z1 / Z0		 => { Z3 }
      Z3	<=  V0 + Z2		 => { Z5 }
      Z5	<=  IPOW( Z3, Z4 )	 => { }
      Z7	<=  Z0 x Z6		 => { Z9 }
      Z9	<=  Z7 + Z8		 => { Z10 }
      Z10	<=  EXP( Z9 )		 => { }
      Z8	<=  -1(I)		 => { Z9 }
      Z1	<=  1(I)		 => { Z2 }
      Z6	<=  2(I)		 => { Z7 }
      Z4	<=  3(I)		 => { Z5 }
\endverbatim

Next, an environment <a>mc::SLiftEnv</a> is defined for lifting the factorable expressions in <a>DAG</a>. The method <a>mc::SLiftEnv::process</a> decomposes the factorable expressions recursively into sparse polynomial and transcendental subexpressions:

\code
      mc::SLiftEnv<mc::FFGraph<>> SPE( &DAG );
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
    6 participating variables: V0 V1 V2 V3 V4 V5 

    2 dependent expressions: 

    OPERATIONS IN SUBGRAPH:
      V0	<-  VARIABLE
      Z4	<-  3(I)	
      Z17	<<  IPOW( V0, Z4 )
      V3	<-  VARIABLE
      Z18	<<  SQR( V0 )	
      Z19	<<  V3 x Z18	
      Z20	<<  Z19 x Z4	
      Z21	<<  Z17 + Z20	
      Z22	<<  SQR( V3 )	
      Z23	<<  V0 x Z22	
      Z24	<<  Z23 x Z4	
      Z25	<<  Z21 + Z24	
      Z26	<<  IPOW( V3, Z4 )
      Z27	<<  Z25 + Z26	
    DEPENDENTS IN SUBGRAPH:
      0:  Z27
    WORK ARRAY SIZE: 14

    OPERATIONS IN SUBGRAPH:
      V5	<<  VARIABLE
    DEPENDENTS IN SUBGRAPH:
      0:  V5
    WORK ARRAY SIZE: 1

    3 polynomial constraints: 

    OPERATIONS IN SUBGRAPH:
      V2	<<  VARIABLE
      V1	<<  VARIABLE
      Z0	<<  SQR( V1 )	
      Z11	<<  V2 - Z0	
    DEPENDENTS IN SUBGRAPH:
      0:  Z11
    WORK ARRAY SIZE: 4

    OPERATIONS IN SUBGRAPH:
      V2	<<  VARIABLE
      V3	<<  VARIABLE
      Z12	<<  V2 x V3	
      Z8	<<  -1(I)	
      Z13	<<  Z12 + Z8	
    DEPENDENTS IN SUBGRAPH:
      0:  Z13
    WORK ARRAY SIZE: 5

    OPERATIONS IN SUBGRAPH:
      V4	<<  VARIABLE
      V1	<<  VARIABLE
      Z0	<<  SQR( V1 )	
      Z6	<<  2(I)	
      Z7	<<  Z0 x Z6	
      Z8	<<  -1(I)	
      Z9	<<  Z7 + Z8	
      Z14	<<  V4 - Z9	
    DEPENDENTS IN SUBGRAPH:
      0:  Z14
    WORK ARRAY SIZE: 8

    1 transcendental constraints: 

    OPERATIONS IN SUBGRAPH:
      V5	<<  VARIABLE
      V4	<<  VARIABLE
      Z15	<<  EXP( V4 )	
      Z16	<<  V5 - Z15	
    DEPENDENTS IN SUBGRAPH:
      0:  Z16
    WORK ARRAY SIZE: 4

    4 auxiliary variables: Z0->V2 Z2->V3 Z9->V4 Z10->V5 
\endverbatim

These results show that 4 auxiliary variables have been added to the DAG, \f$x_2,\ldots,x_5\f$. These variables can be determined from the following (possibly implicit) equations in terms of both original and lifted variables \f$x_0,\ldots,x_5\f$:
\f{align*}
  \left\{\begin{array}{rcl} 0 & = & x_2 - x_1^2\\ 0 & = & x_2\cdot x_3 - 1\\ 0 & = & x_4 - 2\cdot x_1^2 + 1\\ 0 & = & \exp(x_4) - x_5 \end{array}\right.
\f}
Finally, the original vector-valued function \f${\bf f}\f$ is now given by:
\f{align*}
  {\bf f}(x_0,\ldots,x_5) = \left(\begin{array}{c} x_0^3 + 3\cdot x_0^2\cdot x_3 + 3\cdot x_0\cdot x_3^2 + x_3^3\\ x_5
  
  \left(x_0+\frac{1}{x_1^2}\right)^3\\ \exp\left(2\cdot x_1^2-1\right)\end{array}\right)
\f}
*/

// TO DO:
// - Complete documentation
// - Enforce reusing of intermediate variables in polynomial expressions - as a follow-up reduction step?

#ifndef MC__SLIFT_H
#define MC__SLIFT_H

#include <list>
#include "ffunc.hpp"
#include "spoly.hpp"

#define MC__SLIFT_CHECK

namespace mc
{

class SLiftBase;

//! @brief Arithmetic for recursive decomposition of factorable expressions into polynomial/rational and transcendental subexpressions
////////////////////////////////////////////////////////////////////////
//! mc::SLiftVar is a C++ class implementing an arithmetic for
//! recursive decomposition of factorable expressions into polynomial/
//! rational and transcendental subexpressions
////////////////////////////////////////////////////////////////////////
class SLiftVar
////////////////////////////////////////////////////////////////////////
{
public:

  typedef SPoly<FFVar const*, lt_FFVar> t_poly;

private:

  //! @brief Pointer to recursive decomposition environment
  SLiftBase *_env;

protected:

  // numerator sparse polynomial
  t_poly _numer;

  // Denominator sparse polynomial
  t_poly _denom;

  //! @brief Initialize sparse rational expression with existing expression
  SLiftVar& _set
    ( SLiftVar const& var );

  //! @brief Initialize sparse expression as constant
  SLiftVar& _set
    ( double const& d );

  //! @brief Initialize sparse rational expression as DAG variable
  SLiftVar& _set
    ( FFVar const& x );

public:

  //! @brief Default constructor of sparse rational expression
  SLiftVar
    ()
    {}

  //! @brief Constructor of sparse rational expression as constant
  SLiftVar
    ( double const& d )
    : _env( nullptr )
    { _set( d ); }

  //! @brief Constructor of sparse rational expression as DAG variable
  SLiftVar
    ( SLiftBase* env, FFVar const& x);
    
  //! @brief Copy constructor of sparse rational expression
  SLiftVar
    ( SLiftVar const& var );

  //! @brief Constructor of sparse rational expression
  SLiftVar
    ( SLiftBase* env, t_poly const& n, t_poly const& d );

  //! @brief Destructor of sparse rational expression
  virtual ~SLiftVar()
    {}

  //! @brief Initialize variable in sparse rational envrionment <a>env</a> corresponding to DAG variable <a>x</a>
  SLiftVar& set
    ( SLiftBase* env, FFVar const& x );

  //! @brief Overloaded operator '=' for sparse rational expression
  SLiftVar& operator=
    ( SLiftVar const& var )
    { _set( var ); return *this; }

  //! @brief Overloaded operator '=' for constant 
  SLiftVar& operator=
    ( double const& d )
    { _env = nullptr; _set( d ); return *this; }

  //! @brief Overloaded operator '+=' for sparse rational function
  SLiftVar& operator+=
    ( SLiftVar const& var );

  //! @brief Overloaded operator '-=' for sparse rational function
  SLiftVar& operator-=
    ( SLiftVar const& var );

  //! @brief Overloaded operator '*=' for sparse rational function
  SLiftVar& operator*=
    ( SLiftVar const& var );

  //! @brief Overloaded operator '/=' for sparse rational function
  SLiftVar& operator/=
    ( SLiftVar const& var );

  // Sparse rational polynomial environment
  SLiftBase* env
    ()
    const
    { return _env; };

  // numerator sparse polynomial
  t_poly const& numer
    ()
    const
    { return _numer; };

  // Denominator sparse polynomial
  t_poly const& denom
    ()
    const
    { return _denom; };
};

//! @brief C++ base class for recursive decomposition of factorable expressions into polynomial/rational and transcendental subexpressions
////////////////////////////////////////////////////////////////////////
//! mc::SLiftBase is a C++ base class defining the environment for
//! recursive decomposition of factorable expressions into polynomial/
//! rational and transcendental subexpressions
////////////////////////////////////////////////////////////////////////
class SLiftBase
////////////////////////////////////////////////////////////////////////
{
  friend class SLiftVar;
  friend  std::ostream& operator<< ( std::ostream&, SLiftBase const& );
  friend  SLiftVar inv ( SLiftVar const& );
  friend  SLiftVar exp ( SLiftVar const& );
  friend  SLiftVar log ( SLiftVar const& );
  friend  SLiftVar xlog( SLiftVar const& );
  friend  SLiftVar sqrt( SLiftVar const& );
  friend  SLiftVar sqr ( SLiftVar const& );
  friend  SLiftVar pow ( SLiftVar const&, int const );  
  friend  SLiftVar pow ( SLiftVar const&, double const& );  
  friend  SLiftVar cheb( SLiftVar const&, const unsigned );  
  friend  SLiftVar prod( const unsigned, SLiftVar const* );  
  friend  SLiftVar cos ( SLiftVar const& );
  friend  SLiftVar sin ( SLiftVar const& );
  friend  SLiftVar tan ( SLiftVar const& );
  friend  SLiftVar acos( SLiftVar const& );
  friend  SLiftVar asin( SLiftVar const& );
  friend  SLiftVar atan( SLiftVar const& );
  friend  SLiftVar cosh( SLiftVar const& );
  friend  SLiftVar sinh( SLiftVar const& );
  friend  SLiftVar tanh( SLiftVar const& );
  friend  SLiftVar fabs( SLiftVar const& );
  friend  SLiftVar erf( SLiftVar const& );
  friend  SLiftVar fstep( SLiftVar const& );
  friend  SLiftVar max( SLiftVar const&, SLiftVar const& );  
  friend  SLiftVar min( SLiftVar const&, SLiftVar const& );  
  friend  SLiftVar lmtd( SLiftVar const&, SLiftVar const& );  
  friend  SLiftVar rlmtd( SLiftVar const&, SLiftVar const& );  

public:

  typedef std::list< std::pair< FFOp const*, std::vector<SLiftVar const*> > > t_Interm;
  typedef std::map< FFVar const*, FFVar const*, lt_FFVar > t_Aux;
  typedef std::vector< FFVar const* > t_Expr;
  
  typedef SMon<FFVar const*, lt_FFVar> t_mon;
  typedef SPoly<FFVar const*, lt_FFVar> t_poly;

  //! @brief Default Constructor
  SLiftBase
    ( FFBase* dag=nullptr )
    : _dag( dag )
    {}

  //! @brief Destructor
  virtual ~SLiftBase
    ()
    { _reset(); }
  
  // Retreive pointer to DAG
  FFBase* dag
    ()
    const
    { return _dag; };

  //! @brief Retreive reference to vector of depedent DAG subexpressions
  std::vector<FFVar>& Dep
    ()
    { return _Dep; }

  //! @brief Retreive reference to vector of auxiliary DAG polynomial constraints
  std::vector<FFVar>& Poly
    ()
    { return _Poly; }

  //! @brief Retreive reference to vector of auxiliary DAG non-polynomial constraints
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
    (  FFBase* dag )
    { _dag = dag; _reset(); }

  //! @brief Reset intermediate expressions
  void reset
    ()
    { _reset(); }

  //! @brief Transcribe sparse polynomial into DAG
  FFVar insert_dag
    ( t_poly const& poly, bool const useprod=false, bool const dagaux=false );

  //! @brief Transcribe sparse monomial into DAG
  FFVar insert_dag
    ( t_mon const& mon, bool const useprod=false, bool const dagaux=false );

  //! @brief Lift an external operation
  void lift
    ( unsigned const nres, SLiftVar* pres,
      unsigned const nvar, SLiftVar const* pvar );

  //! @brief Exceptions of mc::SLiftVar
  class Exceptions
  {
   public:
    //! @brief Enumeration type for SLiftVar exception handling
    enum TYPE{
      DAGERR=0,       //!< Operation involving a factorable expression linked to a different DAG
      ENVERR,         //!< Operation between factorable expressions linked to different environments or without an environment
      EXTERNAL,       //!< Invalid external operation
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
        return "mc::SLiftBase\t Operation involving a factorable expression linked to a different DAG is not allowed";
      case ENVERR:
        return "mc::SLiftBase\t Operation between factorable expressions linked to different environments or without an environment is not allowed";
      case EXTERNAL:
        return "mc::SLiftBase\t Invalid external operation";
      case INTERNAL:
      default:
        return "mc::SLiftBase\t Internal error";
      }
    }
   private:
    TYPE _ierr;
  };

  //! @brief Options of mc::SLiftBase
  struct Options
  {
    //! @brief Constructor
    Options():
      LIFTDIV( true ), LIFTIPOW( false )//, LIFTUPOL( false )
      {}
    //! @brief Assignment of mc::SLiftBase::Options
    Options& operator=
      ( Options const& opt ){
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

private:

  //! @brief pointer to underlying dag
  FFBase* _dag;

protected:

  //! @brief Map of DAG variables and sparse expressions
  t_Interm _Interm;

  //! @brief Map of existing DAG auxiliaries to new DAG variables
  t_Aux _Aux;

  //! @brief Vector of dependent DAG subexpressions
  std::vector<FFVar> _Dep;

  //! @brief Vector of auxiliary DAG polynomial/rational constraints
  std::vector<FFVar> _Poly;

  //! @brief Vector of auxiliary DAG transcendental constraints
  std::vector<FFVar> _Trans;

  //! @brief Vector of independent DAG variables participating in expressions
  std::vector<FFVar> _Var;

  //! @brief Vector of independent sparse variables
  std::vector<SLiftVar> _SPVar;

  //! @brief Vector of dependent sparse expressions
  std::vector<SLiftVar> _SPDep;

  //! @brief Add new intermediate uni- or bi-variate expression in _Interm
  void _append_interm
    ( SLiftVar const* var1, SLiftVar const* var2=nullptr );

  //! @brief Insert an auxiliary variable corresponding to a rational/polynomial expression into DAG
  FFVar const* _insert_expr
    ( FFVar const* oper, SLiftVar const* expr, bool const isdep );

  //! @brief Insert a non-rational operation into DAG via the introduction of auxiliary variables
  void _insert_expr
    ( FFOp const* pOp, std::vector<const FFVar*>& vAux );

  //! @brief Insert an external operation into DAG via the introduction of auxiliary variables
  template <typename ExtOp, typename... NextOps>
  void _insert_expr_external
    ( FFOp const* pOp, std::vector<FFVar const*>& vAux,
      ExtOp op, std::tuple<NextOps...> ops );

  //! @brief Find auxiliary variable in DAG and return corresponding new DAG variable (or NULL if undefined)
  FFVar const* _find_aux
    ( FFVar const* aux );

  //! @brief Lift current univariate operation
  SLiftVar _lift_univariate_term
    ( SLiftVar const& var );

  //! @brief Lift current bivariate operation
  SLiftVar _lift_bivariate_term
    ( SLiftVar const& var1, SLiftVar const& var2 );

  //! @brief Erase all entries in _Interm
  void _reset
    ();
};

////////////////////////////////////////////////////////////////////////
//! mc::SLiftEnv is a C++ class defining the environment for recursive
//! decomposition of factorable expressions into polynomial/rational
//! and transcendental subexpressions
////////////////////////////////////////////////////////////////////////
template <typename... ExtOps>
class SLiftEnv
: public SLiftBase
////////////////////////////////////////////////////////////////////////
{
  typedef std::tuple<ExtOps...> t_ExtOps;
  typedef typename remove_first_type< t_ExtOps >::type t_NextExtOps;
  typedef typename first_type_of< ExtOps... >::type FirstExtOps;

public:

  //! @brief Default Constructor
  SLiftEnv
    ( FFGraph<ExtOps...>* dag=nullptr )
    : SLiftBase( dag ),
      _dag( dag )
    {}
  
  // Retreive pointer to DAG
  FFGraph<ExtOps...>* dag
    ()
    const
    { return _dag; };

  //! @brief Process the dependents in set <a>sDep</a>
  void process
    ( std::set<unsigned> const& ndxDep, FFVar const* pDep, bool const add2dag=true );

  //! @brief Process the <a>nDep</a> dependents in array <a>pDep</a>
  void process
    ( unsigned const nDep, FFVar const* pDep, const bool add2dag=true );

  //! @brief Set DAG environment
  void set
    ( FFGraph<ExtOps...>* dag )
    { SLiftBase::set( dag ); _dag = dag; }

protected:

  //! @brief pointer to underlying dag
  FFGraph<ExtOps...>* _dag;
};

////////////////////////////////////////////////////////////////////////

inline std::ostream&
operator<<
( std::ostream& out, SLiftVar const& var )
{
  out << std::endl
      << "NUMERATOR:"   << var.numer()
      << "DENOMINATOR:" << var.denom();
  return out;
}

inline std::ostream&
operator<<
( std::ostream& out, SLiftBase const& env)
{
  unsigned count = 0;

  if( env._Dep.empty() ){
    for( auto const& expr : env._Interm ){
      out << std::endl << "INTERMEDIATE #" << ++count << ": ";
      for( auto const& var : expr.first->varout )
        out << *var << " ";
      out << " = " << *(expr.first) << std::endl;
      unsigned pos = 0;
      for( auto&& oper : expr.second )
        out << "Operand " << *expr.first->varin[pos++] << ": " << *oper;
    }
    
    count = 0;
    for( auto const& spdep : env._SPDep )
      out << std::endl << "DEPENDENT #" << ++count << ": " << spdep;
    return out;
  }

  std::cout << std::endl
            << env._Aux.size() << " AUXILIARY VARIABLE";
  if( env._Aux.size() > 1 ) std::cout << "S";
  std::cout << ":";
  for( auto&& aux : env._Aux )
    std::cout << " " << *aux.first << "->" << *aux.second;
  std::cout << std::endl;

  count = 0;
  for( auto const& expr : env._Dep ){
    std::ostringstream ext; 
    ext << " OF DEPENDENT EXPRESSION #" << ++count;
    env._dag->output( env._dag->subgraph( 1, &expr ), ext.str(), out );
  }

  count = 0;
  for( auto const& expr : env._Poly ){
    std::ostringstream ext; 
    ext << " OF AUXILIARY POLYNOMIAL CONSTRAINT #" << ++count;
    env._dag->output( env._dag->subgraph( 1, &expr ), ext.str(), out );
  }

  count = 0;
  for( auto const& expr : env._Trans ){
    std::ostringstream ext; 
    ext << " OF AUXILIARY NON-POLYNOMIAL CONSTRAINT #" << ++count;
    env._dag->output( env._dag->subgraph( 1, &expr ), ext.str(), out );
  
  }
  return out;
}

template <typename... ExtOps>
inline void
SLiftEnv<ExtOps...>::process
( std::set<unsigned> const& ndxDep, FFVar const* pDep, bool const add2dag )
{
  if( ndxDep.empty() ) return; // Nothing to do!
  std::vector<FFVar> vpDep;//( sDep.begin(), sDep.end() );
  vpDep.reserve( ndxDep.size() );
  for( unsigned const& i : ndxDep ) vpDep.push_back( pDep[i] );
  process( ndxDep.size(), vpDep.data(), add2dag );
}

template <typename... ExtOps>
inline void
SLiftEnv<ExtOps...>::process
( unsigned const nDep, FFVar const* pDep, bool const add2dag )
{
  // Reset intermediate / auxiliary arrays
  _reset();

  // Update participating variables in _Var
  auto sgDep = _dag->subgraph( nDep, pDep );
  for( auto const& Op : sgDep.l_op ){
    if( Op->type != FFOp::VAR ) continue;
    _Var.push_back( *Op->varout[0] );
    _SPVar.push_back( SLiftVar( this, *Op->varout[0] ) );
  }

#ifdef MC__SLIFT_DEBUG_PROCESS
  std::cout << std::endl << _Var.size() << " Original Variables: ";
  for( auto const& var : _Var ) std::cout << var << " ";
  std::cout << std::endl;
#endif

  // Process DAG dependents
  _SPDep.resize( nDep );
  _dag->eval( sgDep, nDep, pDep, _SPDep.data(), _Var.size(), _Var.data(), _SPVar.data() );

#ifdef MC__SLIFT_DEBUG_PROCESS
  std::cout << *this;
  for( unsigned i=0; i<nDep; i++ )
    std::cout << std::endl << "SPDep[" << i << "]:" << _SPDep[i];
#endif

  // No transcription in DAG if <a>add2dag</a> is false
  if( !add2dag ) return;

  // Insert auxiliary expressions into DAG
  for( auto const& expr : _Interm ){
    // Insert all operands and their defining expressions
    std::vector<FFVar const*> vAux;
    auto itSV = expr.second.begin();
    for( auto const& operand : expr.first->varin ){
#ifdef MC__SLIFT_CHECK
      if( itSV == expr.second.end() )
        throw Exceptions( Exceptions::INTERNAL );
#endif
      vAux.push_back( _insert_expr( operand, *itSV, false ) );
      ++itSV;
    }
    // Insert operation result and their defining expressions
    FFOp const* pOp = expr.first;
    if( pOp->type < FFOp::EXTERN )
      _insert_expr( pOp, vAux );
    else if( !sizeof...(ExtOps) )
      throw Exceptions( Exceptions::EXTERNAL );
    else
      _insert_expr_external( pOp, vAux, FirstExtOps(), t_NextExtOps() );
  }

  // Insert terminal expressions into DAG
  // Separate two cases:
  // - _SPDep(i) is a polynomial function -> no need for an auxiliary variable
  // - _SPDep(i) is a rational function   -> need to introduce an auxiliary
  for( unsigned i=0; i<nDep; i++ )
    // Insert operation result and defining expression
    _Dep.push_back( *_insert_expr( pDep+i, &_SPDep.at(i), true ) );

#ifdef MC__SLIFT_DEBUG_PROCESS
  std::cout << std::endl << _Aux.size() << " Auxiliary Variables: ";
  for( auto const& aux : _Aux ) std::cout << *aux.first << "->" << *aux.second << " ";
  std::cout << std::endl;
  std::cout << std::endl << _Dep.size() << " Dependent Expressions: ";
  for( auto const& expr : _Dep ) _dag->output( _dag->subgraph( 1, &expr ) );
  std::cout << std::endl;
  std::cout << std::endl << _Poly.size() << " Auxiliary Polynomial Constraints: " << std::endl;
  for( auto const& expr : _Poly ) _dag->output( _dag->subgraph( 1, &expr ) );
  std::cout << std::endl;
  std::cout << std::endl << _Trans.size() << " Auxiliary Transcendental Constraints: " << std::endl;
  for( auto const& expr : _Trans ) _dag->output( _dag->subgraph( 1, &expr ) );
#endif
}

inline FFVar const*
SLiftBase::_insert_expr
( FFVar const* var, SLiftVar const* expr, bool const isdep )
{ 
  auto itdagvar = _dag->Vars().find( const_cast<FFVar*>(var) );
#ifdef MC__SLIFT_CHECK
  assert( itdagvar != _dag->Vars().end() );
#endif

  // Nothing to do if operand is a DAG constant or leaf variable
  if( var->opdef().first->type == FFOp::VAR || var->opdef().first->type == FFOp::CNST )
    return *itdagvar;

  // Nothing to do if DAG auxiliary was already made a DAG variable
  auto itv = _Aux.find( *itdagvar );
  if( itv != _Aux.end() )
    return itv->second;

  // If polynomial subexpression, add subexpression to DAG without defining new polynomial constraint in _Poly 
  if( isdep && !expr->denom().maxord() ){
    FFVar polyexpr = insert_dag( expr->numer() * expr->denom().coefmon() );
#ifdef MC__SLIFT_CHECK
    auto itdagaux = _dag->Vars().find( &polyexpr );
    assert( itdagvar != _dag->Vars().end() );
    return *itdagaux;
#else
    return _dag->Vars().find( &polyexpr );
#endif  
  }

  // Otherwise rational subexpression, append new DAG variable in _Aux and define new polynomial constraint in _Poly 
#ifdef MC__SLIFT_DEBUG_PROCESS
  std::cout <<std::endl << "operand: " << **itdagvar << std::endl;
#endif
  FFVar newvar( _dag );
  auto itnewvar = _dag->Vars().find( &newvar );
#ifdef MC__SLIFT_CHECK
  assert( itnewvar != _dag->Vars().end() );
#endif
#ifdef MC__SLIFT_DEBUG_PROCESS
  std::cout << "paired with new DAG variable: " << **itnewvar << std::endl;
#endif
  _Aux.insert( std::make_pair( *itdagvar, *itnewvar ) );
  _Var.push_back( **itnewvar );

  FFVar polyctr = **itnewvar * insert_dag( expr->denom() ) - insert_dag( expr->numer() );
  auto itpolyctr = _dag->Vars().find( &polyctr );
#ifdef MC__SLIFT_CHECK
  assert( itpolyctr != _dag->Vars().end() );
#endif
#ifdef MC__SLIFT_DEBUG_PROCESS
  std::cout << "defined by DAG subexpression: ";
  _dag->output( _dag->subgraph( 1, *itpolyctr ) );
#endif
  _Poly.push_back( **itpolyctr );
  return *itnewvar;
}

inline void
SLiftBase::_insert_expr
( FFOp const* pOp, std::vector<FFVar const*>& vAux )
{
  assert( pOp->varout.size() == 1 );
  FFVar const* pres = pOp->varout[0];
#ifdef MC__SLIFT_CHECK
  // Throw exception if DAG auxiliary was already made a DAG variable
  if( _Aux.find( pres ) != _Aux.end() )
    throw Exceptions( Exceptions::INTERNAL );
#endif
#ifdef MC__SLIFT_DEBUG_PROCESS
  std::cout << std::endl << "operand: " << *pres << std::endl;
#endif

  // Append new DAG variable in _Aux
  FFVar newvar( _dag );
  auto itnewvar = _dag->Vars().find( &newvar );
#ifdef MC__SLIFT_CHECK
  assert( itnewvar != _dag->Vars().end() );
#endif
#ifdef MC__SLIFT_DEBUG_PROCESS
  std::cout << "paired with new DAG variable: " << **itnewvar << std::endl;
#endif
  _Aux.insert( std::make_pair( pres, *itnewvar ) );
  _Var.push_back( **itnewvar );

  // Append defining expression in _Poly or _Transc 
  switch( pOp->type ){
   case FFOp::SQR:   return _Poly.push_back( _Var.back() - sqr( *vAux.at(0) ) );
   case FFOp::IPOW:  return _Poly.push_back( _Var.back() - pow( *vAux.at(0), vAux.at(1)->num().n ) );
   case FFOp::CHEB:  return _Poly.push_back( _Var.back() - cheb( *vAux.at(0), vAux.at(1)->num().n ) );
   case FFOp::SQRT:  return _Poly.push_back( sqr( _Var.back() ) - *vAux.at(0) );
   case FFOp::INV:   return _Poly.push_back( _Var.back() * *vAux.at(1) - vAux.at(0)->num().val() );
   case FFOp::DIV:   return _Poly.push_back( _Var.back() * *vAux.at(1) - *vAux.at(0) );
   case FFOp::EXP:   return _Trans.push_back( _Var.back() - exp( *vAux.at(0) ) );
   //case FFOp::LOG:   return _Trans.push_back( _Var.back() - log( *vAux.at(0) ) );
   case FFOp::LOG:   return _Trans.push_back( exp( _Var.back() ) - *vAux.at(0) );
   case FFOp::XLOG:  return _Trans.push_back( _Var.back() - xlog( *vAux.at(0) ) );
   case FFOp::DPOW:  return _Trans.push_back( _Var.back() - pow( *vAux.at(0), vAux.at(1)->num().val() ) );
   case FFOp::COS:   return _Trans.push_back( _Var.back() - cos( *vAux.at(0) ) );
   case FFOp::SIN:   return _Trans.push_back( _Var.back() - sin( *vAux.at(0) ) );
   //case FFOp::TAN:   return _Trans.push_back( _Var.back() - tan( *vAux.at(0) ) );
   case FFOp::TAN:   return _Trans.push_back( atan( _Var.back() ) - *vAux.at(0) );
   case FFOp::ACOS:  return _Trans.push_back( cos( _Var.back() ) - *vAux.at(0) );
   case FFOp::ASIN:  return _Trans.push_back( sin( _Var.back() ) - *vAux.at(0) );
   //case FFOp::ATAN:  return _Trans.push_back( tan( _Var.back() ) - *vAux.at(0) );
   case FFOp::ATAN:  return _Trans.push_back( _Var.back() - atan( *vAux.at(0) ) );
   case FFOp::COSH:  return _Trans.push_back( _Var.back() - cosh( *vAux.at(0) ) );
   case FFOp::SINH:  return _Trans.push_back( _Var.back() - sinh( *vAux.at(0) ) );
   case FFOp::TANH:  return _Trans.push_back( _Var.back() - tanh( *vAux.at(0) ) );
   case FFOp::ERF:   return _Trans.push_back( _Var.back() - erf( *vAux.at(0) ) );
   case FFOp::FABS:  return _Trans.push_back( _Var.back() - fabs( *vAux.at(0) ) );
   case FFOp::FSTEP: return _Trans.push_back( _Var.back() - fstep( *vAux.at(0) ) );
   case FFOp::MINF:  return _Trans.push_back( _Var.back() - min( *vAux.at(0), *vAux.at(1) ) );
   case FFOp::MAXF:  return _Trans.push_back( _Var.back() - max( *vAux.at(0), *vAux.at(1) ) );
   case FFOp::INTER: return _Trans.push_back( _Var.back() - inter( *vAux.at(0), *vAux.at(1) ) );
   case FFOp::VAR:
   case FFOp::CNST:
   case FFOp::SHIFT:
   case FFOp::PLUS:
   case FFOp::NEG:
   case FFOp::MINUS:
   case FFOp::SCALE:
   case FFOp::TIMES: // Throw exception since polynomial elements
   case FFOp::PROD:  throw Exceptions( Exceptions::INTERNAL );
   default:          // Should not reach this point
                     throw Exceptions( Exceptions::INTERNAL );
  }
}

template <typename ExtOp, typename... NextOps>
inline void
SLiftBase::_insert_expr_external
( FFOp const* pOp, std::vector<FFVar const*>& vAux,
  ExtOp op, std::tuple<NextOps...> ops )
{
  if( !pOp || pOp->type < FFOp::EXTERN )
    throw Exceptions( Exceptions::INTERNAL );

  // Current operation matches ExtOp type
  if( typeid(*pOp) == typeid(op) ){

    for( auto const& pres : pOp->varout ){
#ifdef MC__SLIFT_CHECK
      // Throw exception if DAG auxiliary was already made a DAG variable
      if( _Aux.find( pres ) != _Aux.end() )
        throw Exceptions( Exceptions::INTERNAL );
#endif
#ifdef MC__SLIFT_DEBUG_PROCESS
      std::cout << std::endl << "operand: " << *pres << std::endl;
#endif

      // Append new DAG variable in _Aux
      FFVar newvar( _dag );
      auto itnewvar = _dag->Vars().find( &newvar );
#ifdef MC__SLIFT_CHECK
      assert( itnewvar != _dag->Vars().end() );
#endif
#ifdef MC__SLIFT_DEBUG_PROCESS
      std::cout << "paired with new DAG variable: " << **itnewvar << std::endl;
#endif
      _Aux.insert( std::make_pair( pres, *itnewvar ) );
      _Var.push_back( **itnewvar );
    }

    std::vector<FFVar> vVar( vAux.size() ), vRes(pOp->varout.size());
    std::vector<unsigned> mVar( vAux.size(), 0 );
    for( unsigned i=0; i<vAux.size(); ++i )
      vVar[i] = *vAux.at(i);
    op.eval( vRes.size(), vRes.data(), vVar.size(), vVar.data(), mVar.data() );
    FFVar* pRes = _Var.data() + _Var.size() - pOp->varout.size();
    for( unsigned j=0; j<vRes.size(); ++j )
      _Trans.push_back( pRes[j] - vRes[j] );
    return;
  }
  
  // No more external operations to peel off parameter pack
  if( !sizeof...( NextOps ) )
    throw Exceptions( Exceptions::EXTERNAL );

  // Separate first element off parameter pack and recursive call
  typedef std::tuple<NextOps...> t_NextOps;
  typedef typename remove_first_type< t_NextOps >::type t_NextNextOps;
  typedef typename first_type_of< NextOps... >::type FirstNextOps;
  _insert_expr_external( pOp, vAux, FirstNextOps(), t_NextNextOps() );
}

inline FFVar const*
SLiftBase::_find_aux
( FFVar const* aux )
{
  if( aux->opdef().first->type == FFOp::VAR ) return aux;
  auto it = _Aux.find( aux );
  if( it != _Aux.end() ) return it->second;
  return (FFVar*)0;
}

inline FFVar
SLiftBase::insert_dag
( t_poly const& pol, bool const useprod, bool const dagaux )
{
  assert( !pol.mapmon().empty() );

  FFVar var;
  bool first = true;
  for( auto const& [mon,coef] : pol.mapmon() ){
    if( first ){
      var  = coef * insert_dag( mon, useprod, dagaux );
      first = false;
    }
    else if( coef >= 0 )
      var += coef * insert_dag( mon, useprod, dagaux );
    else
      var -= (-coef) * insert_dag( mon, useprod, dagaux );
  }
  return var;
}

inline FFVar
SLiftBase::insert_dag
( t_mon const& mon, bool const useprod, bool const dagaux )
{
  if( useprod ){
    std::vector<FFVar> pvar;
    pvar.reserve( mon.expr.size() );
  for( auto const& [var,ord] : mon.expr ){
      FFVar const* oper = nullptr;
      if( dagaux ) oper = *_dag->Vars().find( const_cast<FFVar*>(var) );
      else         oper = _find_aux( var );
      //FFVar const* oper = _find_aux( var );
      if( oper == nullptr )
        throw Exceptions( Exceptions::INTERNAL );
      switch( t_poly::options.BASIS ){
       case t_poly::Options::MONOM:
        pvar.push_back( pow( *oper, (int)ord ) );
        break;
       case t_poly::Options::CHEB:
        pvar.push_back( cheb( *oper, ord ) );
        break;
      }
    }
    return prod( pvar.size(), pvar.data() );
  }
    
  FFVar prodmon = 1;
  for( auto const& [var,ord] : mon.expr ){
    //std::cout << "[var,ord] = " << *var << "^" << ord << std::endl;
    FFVar const* oper = nullptr;
    if( dagaux ) oper = *_dag->Vars().find( const_cast<FFVar*>(var) );
    else         oper = _find_aux( var );
    //FFVar const* oper = _find_aux( var );
    if( oper == nullptr )
      throw Exceptions( Exceptions::INTERNAL );
    switch( t_poly::options.BASIS ){
     case t_poly::Options::MONOM:
      prodmon *= pow( *oper, (int)ord );
      break;
     case t_poly::Options::CHEB:
      prodmon *= cheb( *oper, ord );
      break;
    }
  }
  return prodmon;
}

inline void
SLiftBase::_reset
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
  _Dep.clear();
  _SPVar.clear();
  _SPDep.clear();
}

inline void
SLiftBase::lift
( unsigned const nres, SLiftVar* pres,
  unsigned const nvar, SLiftVar const* pvar )
{
#ifdef MC__SLIFT_CHECK
  if( !nvar || !pvar || !nres || !pres )
    throw Exceptions( Exceptions::EXTERNAL );
#endif
  FFOp const* pOp = _dag->curOp();
  std::vector<SLiftVar const*> vops;
  vops.reserve( nvar );
  for( unsigned i=0; i<nvar; ++i )
    vops.push_back( new SLiftVar( pvar[i] ) );
  _Interm.push_back( std::make_pair( pOp, vops ) );
  for( unsigned j=0; j<nres; ++j )
    pres[j] = SLiftVar( this, *(pOp->varout[j]) );
}

inline void
SLiftBase::_append_interm
( SLiftVar const* var1, SLiftVar const* var2 )
{
#ifdef MC__SLIFT_CHECK
  if( !var1 ) throw Exceptions( Exceptions::INTERNAL );
#endif
  FFOp const* pOp = _dag->curOp();
  std::vector<SLiftVar const*> vops;
  vops.push_back( new SLiftVar( *var1 ) );
  if( var2 ) vops.push_back( new SLiftVar( *var2 ) );
  _Interm.push_back( std::make_pair( pOp, vops ) );
}

inline SLiftVar
SLiftBase::_lift_univariate_term
( SLiftVar const& var )
{
  FFOp const* pOp = _dag->curOp();
#ifdef MC__SLIFT_CHECK
  if( !pOp ) throw Exceptions( Exceptions::INTERNAL );
#endif

  // Append new intermediate expression and assert that same operation was not previously appended
  assert( pOp->varout.size() == 1 );
  _append_interm( &var );
  return SLiftVar( this, *(pOp->varout[0]) );
}

inline SLiftVar
SLiftBase::_lift_bivariate_term
( SLiftVar const& var1, SLiftVar const& var2 )
{
  FFOp const* pOp = _dag->curOp();
#ifdef MC__SLIFT_CHECK
  if( !pOp ) throw Exceptions( Exceptions::INTERNAL );
#endif

  // Append new intermediate expression and assert that same operation was not previously appended
  assert( pOp->varout.size() == 1 );
  _append_interm( &var1, &var2 );
  return SLiftVar( this, *(pOp->varout[0]) );
}

inline
SLiftVar::SLiftVar
( SLiftBase* env, FFVar const& x)
: _env( env )
{
  if( env->dag() != x.dag() )
  throw typename SLiftBase::Exceptions( SLiftBase::Exceptions::DAGERR );
  _set( x );
}

inline
SLiftVar::SLiftVar
( SLiftVar const& var )
{
  _set( var );
}

inline
SLiftVar::SLiftVar
( SLiftBase* env, t_poly const& n, t_poly const& d )
: _env( env ), _numer( n ), _denom( d )
{}

inline SLiftVar&
SLiftVar::set
( SLiftBase* env, FFVar const& x )
{
  _env = env;
  if( env->dag() != x.dag() )
  throw typename SLiftBase::Exceptions( SLiftBase::Exceptions::DAGERR );
  _set( x );
  return *this;
}

inline SLiftVar&
SLiftVar::_set
( double const& d )
{
  _numer = d;
  _denom = 1;
  return *this;
}

inline SLiftVar&
SLiftVar::_set
( FFVar const& x )
{
  _numer.var( &x );
  _denom = 1;
  return *this;
}

inline SLiftVar&
SLiftVar::_set
( SLiftVar const& var )
{
  if( this == &var ) return *this;
  _env = var._env;
  _numer = var._numer;
  _denom = var._denom;
  return *this;
}

inline SLiftVar
operator+
( SLiftVar const& var )
{
  return var;
}

inline SLiftVar&
SLiftVar::operator+=
( SLiftVar const& var )
{
  if( !_env )
    _env = var._env;
  else if( var._env && _env != var._env )
    throw typename SLiftBase::Exceptions( SLiftBase::Exceptions::ENVERR );

  _numer *= var._denom;
  _numer += var._numer * _denom;
  _denom *= var._denom;
  return *this;
}

inline SLiftVar
operator+
( SLiftVar const& var1, SLiftVar const& var2 )
{
  SLiftVar var3( var1 );
  var3 += var2;
  return var3;
}

inline SLiftVar
operator+
( SLiftVar const& var1, double const& cst2 )
{
  SLiftVar var3( var1 );
  var3 += cst2;
  return var3;
}

inline SLiftVar
operator-
( SLiftVar const& var )
{
  return SLiftVar( var.env(), -var.numer(), var.denom() );
}

inline SLiftVar&
SLiftVar::operator-=
( SLiftVar const& var )
{
  if( !_env )
    _env = var._env;
  else if( var._env && _env != var._env )
    throw typename SLiftBase::Exceptions( SLiftBase::Exceptions::ENVERR );

  _numer *= var._denom;
  _numer -= var._numer * _denom;
  _denom *= var._denom;
  return *this;
}

inline SLiftVar
operator-
( SLiftVar const& var1, SLiftVar const& var2 )
{
  SLiftVar var3( var1 );
  var3 -= var2;
  return var3;
}

inline SLiftVar
operator-
( SLiftVar const& var1, double const& cst2 )
{
  SLiftVar var3( var1 );
  var3 -= cst2;
  return var3;
}

inline SLiftVar&
SLiftVar::operator*=
( SLiftVar const& var )
{
  if( !_env )
    _env = var._env;
  else if( var._env && _env != var._env )
    throw typename SLiftBase::Exceptions( SLiftBase::Exceptions::ENVERR );

  _numer *= var._numer;
  _denom *= var._denom;
  return *this;
}

inline SLiftVar
operator*
( SLiftVar const& var1, SLiftVar const& var2 )
{
  SLiftVar var3( var1 );
  var3 *= var2;
  return var3;
}

inline SLiftVar
operator*
( SLiftVar const& var1, double const& cst2 )
{
  SLiftVar var3( var1 );
  var3 *= cst2;
  return var3;
}

inline SLiftVar&
SLiftVar::operator/=
( SLiftVar const& var )
{
  if( !_env )
    _env = var._env;
  else if( var._env && _env != var._env )
    throw typename SLiftBase::Exceptions( SLiftBase::Exceptions::ENVERR );

  if( _env->options.LIFTDIV )
    *this = _env->_lift_bivariate_term( *this, var );
  else{
    _numer *= var._denom;
    _denom *= var._numer;
  }
  return *this;
}

inline SLiftVar
operator/
( SLiftVar const& var1, SLiftVar const& var2 )
{
  SLiftVar var3( var1 );
  var3 /= var2;
  return var3;
}

inline SLiftVar
operator/
( double const& cst1, SLiftVar const& var2 )
{
  SLiftVar var3( cst1 );
  var3 /= var2;
  return var3;
}

inline SLiftVar
sqr
( SLiftVar const& var )
{
#ifdef MC__SLIFT_CHECK
  if( !var.env() )
    throw typename SLiftBase::Exceptions( SLiftBase::Exceptions::ENVERR );
#endif
  if( var.env()->options.LIFTIPOW )
    return var.env()->_lift_univariate_term( var );
  return SLiftVar( var.env(), sqr( var.numer() ), sqr( var.denom() ) );
}

inline SLiftVar
inv
( SLiftVar const& var )
{
#ifdef MC__SLIFT_CHECK
  if( !var.env() )
    throw typename SLiftBase::Exceptions( SLiftBase::Exceptions::ENVERR );
#endif
  if( var.env()->options.LIFTDIV )
    return var.env()->_lift_bivariate_term( 1, var );
  return SLiftVar( var.env(), var.denom(), var.numer() );
}

inline SLiftVar
pow
( SLiftVar const& var, double const& a )
{
#ifdef MC__SLIFT_CHECK
  if( !var.env() )
    throw typename SLiftBase::Exceptions( SLiftBase::Exceptions::ENVERR );
#endif
  return var.env()->_lift_bivariate_term( var, a );
}

inline SLiftVar
pow
( SLiftVar const& var, int const n )
{
  if( n < 0 ) return pow( inv( var ), -n );
  switch( n ){
   case 0:  return 1.;
   case 1:  return var;
   case 2:  return sqr( var );
   default: if( var.env()->options.LIFTIPOW ) return var.env()->_lift_bivariate_term( var, n );
            return SLiftVar( var.env(), pow( var.numer(), (unsigned)n ), pow( var.denom(), (unsigned)n ) );
  }
}

inline SLiftVar
cheb
( SLiftVar const& var, unsigned const n )
{
  switch( n ){
   case 0:  return 1.;
   case 1:  return var;
   case 2:  return sqr( var ) * 2. - 1.;
   default: if( var.env()->options.LIFTIPOW ) return var.env()->_lift_bivariate_term( var, n );
            return var * cheb( var, n-1 ) * 2. - cheb( var, n-2 );
  }
}

inline SLiftVar
prod
( unsigned int const nvars, SLiftVar const* pvars )
{
  switch( nvars ){
   case 0:  return 1.;
   case 1:  return pvars[0];
   default: return pvars[0] * prod( nvars-1, pvars+1 );
  }
}

inline SLiftVar
monom
( unsigned int const nvars, SLiftVar const* pvars, unsigned const* k, bool const chebbasis=false )
{
  switch( nvars ){
   case 0:  return 1.;
   case 1:  return chebbasis? cheb( pvars[0], k[0] ): pow( pvars[0], (int)k[0] );
   default: return ( chebbasis? cheb( pvars[0], k[0] ): pow( pvars[0], (int)k[0] ) ) * monom( nvars-1, pvars+1, k+1 );
  }
}

inline SLiftVar
exp
( SLiftVar const& var )
{
#ifdef MC__SLIFT_CHECK
  if( !var.env() )
    throw typename SLiftBase::Exceptions( SLiftBase::Exceptions::ENVERR );
#endif
  return var.env()->_lift_univariate_term( var );
}

inline SLiftVar
log
( SLiftVar const& var )
{
#ifdef MC__SLIFT_CHECK
  if( !var.env() )
    throw typename SLiftBase::Exceptions( SLiftBase::Exceptions::ENVERR );
#endif
  return var.env()->_lift_univariate_term( var );
}

inline SLiftVar
xlog
( SLiftVar const& var )
{
#ifdef MC__SLIFT_CHECK
  if( !var.env() )
    throw typename SLiftBase::Exceptions( SLiftBase::Exceptions::ENVERR );
#endif
  return var.env()->_lift_univariate_term( var );
}

inline SLiftVar
sqrt
( SLiftVar const& var )
{
#ifdef MC__SLIFT_CHECK
  if( !var.env() )
    throw typename SLiftBase::Exceptions( SLiftBase::Exceptions::ENVERR );
#endif
  return var.env()->_lift_univariate_term( var );
}

inline SLiftVar
cos
( SLiftVar const& var )
{
#ifdef MC__SLIFT_CHECK
  if( !var.env() )
    throw typename SLiftBase::Exceptions( SLiftBase::Exceptions::ENVERR );
#endif
  return var.env()->_lift_univariate_term( var );
}

inline SLiftVar
sin
( SLiftVar const& var )
{
#ifdef MC__SLIFT_CHECK
  if( !var.env() )
    throw typename SLiftBase::Exceptions( SLiftBase::Exceptions::ENVERR );
#endif
  return var.env()->_lift_univariate_term( var );
}

inline SLiftVar
tan
( SLiftVar const& var )
{
#ifdef MC__SLIFT_CHECK
  if( !var.env() )
    throw typename SLiftBase::Exceptions( SLiftBase::Exceptions::ENVERR );
#endif
  return var.env()->_lift_univariate_term( var );
}

inline SLiftVar
acos
( SLiftVar const& var )
{
#ifdef MC__SLIFT_CHECK
  if( !var.env() )
    throw typename SLiftBase::Exceptions( SLiftBase::Exceptions::ENVERR );
#endif
  return var.env()->_lift_univariate_term( var );
}

inline SLiftVar
asin
( SLiftVar const& var )
{
#ifdef MC__SLIFT_CHECK
  if( !var.env() )
    throw typename SLiftBase::Exceptions( SLiftBase::Exceptions::ENVERR );
#endif
  return var.env()->_lift_univariate_term( var );
}

inline SLiftVar
atan
( SLiftVar const& var )
{
#ifdef MC__SLIFT_CHECK
  if( !var.env() )
    throw typename SLiftBase::Exceptions( SLiftBase::Exceptions::ENVERR );
#endif
  return var.env()->_lift_univariate_term( var );
}

inline SLiftVar
cosh
( SLiftVar const& var )
{
#ifdef MC__SLIFT_CHECK
  if( !var.env() )
    throw typename SLiftBase::Exceptions( SLiftBase::Exceptions::ENVERR );
#endif
  return var.env()->_lift_univariate_term( var );
}

inline SLiftVar
sinh
( SLiftVar const& var )
{
#ifdef MC__SLIFT_CHECK
  if( !var.env() )
    throw typename SLiftBase::Exceptions( SLiftBase::Exceptions::ENVERR );
#endif
  return var.env()->_lift_univariate_term( var );
}

inline SLiftVar
tanh
( SLiftVar const& var )
{
#ifdef MC__SLIFT_CHECK
  if( !var.env() )
    throw typename SLiftBase::Exceptions( SLiftBase::Exceptions::ENVERR );
#endif
  return var.env()->_lift_univariate_term( var );
}

inline SLiftVar
fabs
( SLiftVar const& var )
{
#ifdef MC__SLIFT_CHECK
  if( !var.env() )
    throw typename SLiftBase::Exceptions( SLiftBase::Exceptions::ENVERR );
#endif
  return var.env()->_lift_univariate_term( var );
}

inline SLiftVar
erf
( SLiftVar const& var )
{
#ifdef MC__SLIFT_CHECK
  if( !var.env() )
    throw typename SLiftBase::Exceptions( SLiftBase::Exceptions::ENVERR );
#endif
  return var.env()->_lift_univariate_term( var );
}

inline SLiftVar
fstep
( SLiftVar const& var )
{
#ifdef MC__SLIFT_CHECK
  if( !var.env() )
    throw typename SLiftBase::Exceptions( SLiftBase::Exceptions::ENVERR );
#endif
  return var.env()->_lift_univariate_term( var );
}

inline SLiftVar
max
( SLiftVar const& var1, SLiftVar const& var2 )
{
  if( var1.env() && var2.env() && var1.env() != var2.env() )
    throw typename SLiftBase::Exceptions( SLiftBase::Exceptions::ENVERR );

  return var1.env()->_lift_bivariate_term( var1, var2 );
}

inline SLiftVar
min
( SLiftVar const& var1, SLiftVar const& var2 )
{
  if( var1.env() && var2.env() && var1.env() != var2.env() )
    throw typename SLiftBase::Exceptions( SLiftBase::Exceptions::ENVERR );

  return var1.env()->_lift_bivariate_term( var1, var2 );
}

inline SLiftVar
lmtd
( SLiftVar const& var1, SLiftVar const& var2 )
{
  if( var1.env() && var2.env() && var1.env() != var2.env() )
    throw typename SLiftBase::Exceptions( SLiftBase::Exceptions::ENVERR );

  return var1.env()->_lift_bivariate_term( var1, var2 );
}

inline SLiftVar
rlmtd
( SLiftVar const& var1, SLiftVar const& var2 )
{
  if( var1.env() && var2.env() && var1.env() != var2.env() )
    throw typename SLiftBase::Exceptions( SLiftBase::Exceptions::ENVERR );

  return var1.env()->_lift_bivariate_term( var1, var2 );
}

} // namespace mc

//#include "mcop.hpp"

namespace mc
{

//! @brief Specialization of the structure mc::Op to allow usage of the type mc::SLiftVar for DAG evaluation or as a template parameter in other MC++ classes
template <> struct Op<mc::SLiftVar>
{
  typedef mc::SLiftVar T;
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
