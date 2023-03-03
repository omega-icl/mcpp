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
// - Split subsets of polynomial and transcendental expressions
// - Enforce reusing of intermediate variables in polynomial expressions - as a follow-up reduction step?

#ifndef MC__SLIFT_H
#define MC__SLIFT_H

#include <list>
#include "ffunc.hpp"
#include "spoly.hpp"

#define MC__SLIFT_CHECK

namespace mc
{

template <typename DAG> class SLiftVar;

//! @brief Environment for recursive decomposition of factorable expressions into polynomial/rational and transcendental subexpressions
////////////////////////////////////////////////////////////////////////
//! mc::SLiftEnv is a C++ class defining the environment for recursive
//! decomposition of factorable expressions into polynomial/rational
//! and transcendental subexpressions
////////////////////////////////////////////////////////////////////////
template < typename DAG >
class SLiftEnv
////////////////////////////////////////////////////////////////////////
{
  friend class SLiftVar<DAG>;
  template <typename D> friend  std::ostream& operator<< ( std::ostream&, SLiftEnv<D> const& );
  template <typename D> friend  SLiftVar<D> inv ( SLiftVar<D> const& );
  template <typename D> friend  SLiftVar<D> exp ( SLiftVar<D> const& );
  template <typename D> friend  SLiftVar<D> log ( SLiftVar<D> const& );
  template <typename D> friend  SLiftVar<D> xlog( SLiftVar<D> const& );
  template <typename D> friend  SLiftVar<D> sqrt( SLiftVar<D> const& );
  template <typename D> friend  SLiftVar<D> sqr ( SLiftVar<D> const& );
  template <typename D> friend  SLiftVar<D> pow ( SLiftVar<D> const&, int const );  
  template <typename D> friend  SLiftVar<D> pow ( SLiftVar<D> const&, double const& );  
  template <typename D> friend  SLiftVar<D> cheb( SLiftVar<D> const&, const unsigned );  
  template <typename D> friend  SLiftVar<D> prod( const unsigned, SLiftVar<D> const* );  
  template <typename D> friend  SLiftVar<D> cos ( SLiftVar<D> const& );
  template <typename D> friend  SLiftVar<D> sin ( SLiftVar<D> const& );
  template <typename D> friend  SLiftVar<D> tan ( SLiftVar<D> const& );
  template <typename D> friend  SLiftVar<D> acos( SLiftVar<D> const& );
  template <typename D> friend  SLiftVar<D> asin( SLiftVar<D> const& );
  template <typename D> friend  SLiftVar<D> atan( SLiftVar<D> const& );
  template <typename D> friend  SLiftVar<D> cosh( SLiftVar<D> const& );
  template <typename D> friend  SLiftVar<D> sinh( SLiftVar<D> const& );
  template <typename D> friend  SLiftVar<D> tanh( SLiftVar<D> const& );
  template <typename D> friend  SLiftVar<D> fabs( SLiftVar<D> const& );
  template <typename D> friend  SLiftVar<D> erf( SLiftVar<D> const& );
  template <typename D> friend  SLiftVar<D> fstep( SLiftVar<D> const& );
  template <typename D> friend  SLiftVar<D> max( SLiftVar<D> const&, SLiftVar<D> const& );  
  template <typename D> friend  SLiftVar<D> min( SLiftVar<D> const&, SLiftVar<D> const& );  
  template <typename D> friend  SLiftVar<D> lmtd( SLiftVar<D> const&, SLiftVar<D> const& );  
  template <typename D> friend  SLiftVar<D> rlmtd( SLiftVar<D> const&, SLiftVar<D> const& );  

public:

  typedef std::list< std::pair< FFOp const*, std::vector<SLiftVar<DAG> const*> > > t_Interm;
  typedef std::map< FFVar const*, FFVar const*, lt_FFVar > t_Aux;
  typedef std::vector< FFVar const* > t_Expr;
  
  typedef SMon<FFVar const*, lt_FFVar> t_mon;
  typedef SPoly<FFVar const*, lt_FFVar> t_poly;

  //! @brief Default Constructor
  SLiftEnv
    ( DAG* dag=nullptr )
    : _dag( dag )
    {}

  //! @brief Destructor
  virtual ~SLiftEnv
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

  //! @brief Transcribe sparse polynomial into DAG
  FFVar insert_dag
    ( t_poly const& poly, bool const useprod=false, bool const dagaux=false );

  //! @brief Transcribe sparse monomial into DAG
  FFVar insert_dag
    ( t_mon const& mon, bool const useprod=false, bool const dagaux=false );

  //! @brief Exceptions of mc::SLiftVar
  class Exceptions
  {
   public:
    //! @brief Enumeration type for SLiftVar exception handling
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
        return "mc::SLiftEnv\t Operation involving a factorable expression linked to a different DAG is not allowed";
      case ENVERR:
        return "mc::SLiftEnv\t Operation between factorable expressions linked to different environments or without an environment is not allowed";
      case INTERNAL:
      default:
        return "mc::SLiftEnv\t Internal error";
      }
    }
   private:
    TYPE _ierr;
  };

  //! @brief Options of mc::SLiftEnv
  static struct Options
  {
    //! @brief Constructor
    Options():
      LIFTDIV( true ), LIFTIPOW( false )//, LIFTUPOL( false )
      {}
    //! @brief Assignment of mc::SLiftEnv<DAG>::Options
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
  std::vector<SLiftVar<DAG>> _SPVar;

  //! @brief Vector of dependent sparse expressions
  std::vector<SLiftVar<DAG>> _SPDep;

  //! @brief Add new intermediate uni- or bi-variate expression in _Interm
  void _append_interm
    ( FFOp const* op, SLiftVar<DAG> const* var1, SLiftVar<DAG> const* var2=nullptr );

//  //! @brief Add new intermediate n-variate expression in _Interm
//  void _append_interm
//    ( FFOp const* op, unsigned const nvars, SLiftVar<DAG> const*const* pvars );

  //! @brief Insert an auxiliary variable corresponding to a rational/polynomial expression into DAG
  FFVar const* _insert_expr
    ( FFVar const* oper, SLiftVar<DAG> const* expr );

  //! @brief Insert a non-rational operation into DAG via the introduction of auxiliary variables
  void _insert_expr
    ( FFOp const* pOp, std::vector<const FFVar*>& vAux );

  //! @brief Find auxiliary variable in DAG and return corresponding new DAG variable (or NULL if undefined)
  FFVar const* _find_aux
    ( FFVar const* aux );

  //! @brief Lift current univariate operation
  SLiftVar<DAG> _lift_univariate_term
    ( SLiftVar<DAG> const& var );

  //! @brief Lift current bivariate operation
  SLiftVar<DAG> _lift_bivariate_term
    ( SLiftVar<DAG> const& var1, SLiftVar<DAG> const& var2 );

  //! @brief Erase all entries in _Interm
  void _reset
    ();
};


template <typename DAG>
inline typename SLiftEnv<DAG>::Options SLiftEnv<DAG>::options;

//! @brief Arithmetic for recursive decomposition of factorable expressions into polynomial/rational and transcendental subexpressions
////////////////////////////////////////////////////////////////////////
//! mc::SLiftVar is a C++ class implementing an arithmetic for
//! recursive decomposition of factorable expressions into polynomial/
//! rational and transcendental subexpressions
////////////////////////////////////////////////////////////////////////
template < typename DAG >
class SLiftVar
////////////////////////////////////////////////////////////////////////
{
public:

  typedef SPoly<FFVar const*, lt_FFVar> t_poly;

private:

  //! @brief Pointer to recursive decomposition environment
  SLiftEnv<DAG> *_env;

protected:

  // numerator sparse polynomial
  t_poly _numer;

  // Denominator sparse polynomial
  t_poly _denom;

  //! @brief Initialize sparse rational expression with existing expression
  SLiftVar<DAG>& _set
    ( SLiftVar<DAG> const& var );

  //! @brief Initialize sparse expression as constant
  SLiftVar<DAG>& _set
    ( double const& d );

  //! @brief Initialize sparse rational expression as DAG variable
  SLiftVar<DAG>& _set
    ( FFVar const& x );

public:

  //! @brief Default constructor of sparse rational expression
  SLiftVar
    ()
    {}

  //! @brief Constructor of sparse rational expression as constant
  SLiftVar
    ( double const& d )
    : _env( 0 )
    { _set( d ); }

  //! @brief Constructor of sparse rational expression as DAG variable
  SLiftVar
    ( SLiftEnv<DAG>* env, FFVar const& x)
    : _env( env )
    { if( env->dag() != x.dag() )
        throw typename SLiftEnv<DAG>::Exceptions( SLiftEnv<DAG>::Exceptions::DAGERR );
      _set( x ); }

  //! @brief Copy constructor of sparse rational expression
  SLiftVar
    ( SLiftVar<DAG> const& var )
    { _set( var ); }

  //! @brief Constructor of sparse rational expression
  SLiftVar
    ( SLiftEnv<DAG>* env, t_poly const& n, t_poly const& d )
    : _env( env ), _numer( n ), _denom( d )
    {}

  //! @brief Destructor of sparse rational expression
  virtual ~SLiftVar()
    {}

  //! @brief Initialize variable in sparse rational envrionment <a>env</a> corresponding to DAG variable <a>x</a>
  SLiftVar<DAG>& set
    ( SLiftEnv<DAG>* env, FFVar const& x )
    { _env = env;
      if( env->dag() != x.dag() )
        throw typename SLiftEnv<DAG>::Exceptions( SLiftEnv<DAG>::Exceptions::DAGERR );
      _set( x ); return *this; }

  //! @brief Overloaded operator '=' for sparse rational expression
  SLiftVar<DAG>& operator=
    ( SLiftVar<DAG> const& var )
    { _set( var ); return *this; }

  //! @brief Overloaded operator '=' for constant 
  SLiftVar<DAG>& operator=
    ( double const& d )
    { _env = 0; _set( d ); return *this; }

  //! @brief Overloaded operator '+=' for sparse rational function
  SLiftVar<DAG>& operator+=
    ( SLiftVar<DAG> const& var );

  //! @brief Overloaded operator '-=' for sparse rational function
  SLiftVar<DAG>& operator-=
    ( SLiftVar<DAG> const& var );

  //! @brief Overloaded operator '*=' for sparse rational function
  SLiftVar<DAG>& operator*=
    ( SLiftVar<DAG> const& var );

  //! @brief Overloaded operator '/=' for sparse rational function
  SLiftVar<DAG>& operator/=
    ( SLiftVar<DAG> const& var );

  // Sparse rational polynomial environment
  SLiftEnv<DAG>* env
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

////////////////////////////////////////////////////////////////////////

template < typename DAG >
inline std::ostream&
operator<<
( std::ostream& out, SLiftVar<DAG> const& var )
{
  out << std::endl
      << "NUMERATOR:"   << var.numer()
      << "DENOMINATOR:" << var.denom();
  return out;
}

template < typename DAG >
inline std::ostream&
operator<<
( std::ostream& out, SLiftEnv<DAG> const& env)
{
  unsigned count = 0;

  if( env._Poly.empty() ){
    for( auto const& expr : env._Interm ){
      out << std::endl << "INTERMEDIATE #" << ++count << ": "
          << *(expr.first->pres) << " = " << *(expr.first) << std::endl;
      unsigned pos = 0;
      for( auto&& oper : expr.second )
        out << "Operand " << *expr.first->pops[pos++] << ": " << *oper;
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
  for( auto&& expr : env._Poly ){
    std::ostringstream ext; 
    ext << " OF POLYNOMIAL CONSTRAINT #" << ++count;
    env._dag->output( env._dag->subgraph( 1, &expr ), ext.str(), out );
  }

  count = 0;
  for( auto&& expr : env._Trans ){
    std::ostringstream ext; 
    ext << " OF NON-POLYNOMIAL CONSTRAINT #" << ++count;
    env._dag->output( env._dag->subgraph( 1, &expr ), ext.str(), out );
  
  }
  return out;
}

template < typename DAG >
inline void
SLiftEnv<DAG>::process
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
SLiftEnv<DAG>::process
( unsigned const nDep, FFVar const* pDep, bool const add2dag )
{
  // Reset intermediate / auxiliary arrays
  _reset();

  // Update participating variables in _Var
  auto sgDep = _dag->subgraph( nDep, pDep );
  for( auto const& [Op,mov] : sgDep.l_op ){
    if( Op->type != FFOp::VAR ) continue;
    _Var.push_back( *Op->pres );
    _SPVar.push_back( SLiftVar( this, *Op->pres ) );
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
    for( auto const& operand : expr.first->pops ){
#ifdef MC__SLIFT_CHECK
      if( itSV == expr.second.end() )
        throw typename SLiftEnv<DAG>::Exceptions( SLiftEnv<DAG>::Exceptions::INTERNAL );
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

#ifdef MC__SLIFT_DEBUG_PROCESS
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
SLiftEnv<DAG>::_insert_expr
( FFVar const* var, SLiftVar<DAG> const* expr )
{ 
  auto itdagvar = _dag->Vars().find( const_cast<FFVar*>(var) );
#ifdef MC__SLIFT_CHECK
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

template < typename DAG >
inline void
SLiftEnv<DAG>::_insert_expr
( FFOp const* pOp, std::vector<FFVar const*>& vAux )
{ 
#ifdef MC__SLIFT_CHECK
  // Throw exception if DAG auxiliary was already made a DAG variable
  if( _Aux.find( pOp->pres ) != _Aux.end() )
    throw typename SLiftEnv<DAG>::Exceptions( SLiftEnv<DAG>::Exceptions::INTERNAL );
#endif
#ifdef MC__SLIFT_DEBUG_PROCESS
  std::cout << std::endl << "operand: " << *pOp->pres << std::endl;
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
   default:          throw Exceptions( Exceptions::INTERNAL );
  }
}

template < typename DAG >
inline FFVar const*
SLiftEnv<DAG>::_find_aux
( FFVar const* aux )
{
  if( aux->ops().first->type == FFOp::VAR ) return aux;
  auto it = _Aux.find( aux );
  if( it != _Aux.end() ) return it->second;
  return (FFVar*)0;
}

template < typename DAG >
inline FFVar
SLiftEnv<DAG>::insert_dag
( t_poly const& pol, bool const useprod, bool const dagaux )
{
  FFVar var;
  bool first = true;
  for( auto const& [mon,coef] : pol.mapmon() ){
    if( first ){
      var  = coef * insert_dag( mon, useprod, dagaux );
      first = false;
    }
    else if( coef >= 0 )
      var += coef * insert_dag( mon, useprod,dagaux );
    else
      var -= (-coef) * insert_dag( mon, useprod,dagaux );
  }
  return var;
}

template < typename DAG >
inline FFVar
SLiftEnv<DAG>::insert_dag
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

template < typename DAG >
inline void
SLiftEnv<DAG>::_reset
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
SLiftEnv<DAG>::_append_interm
( FFOp const* op, SLiftVar<DAG> const* var1, SLiftVar<DAG> const* var2 )
{
#ifdef MC__SLIFT_CHECK
  if( !var1 )
    throw typename SLiftEnv<DAG>::Exceptions( SLiftEnv<DAG>::Exceptions::INTERNAL );
#endif
  std::vector<SLiftVar<DAG> const*> vops;
  vops.push_back( new SLiftVar( *var1 ) );
  if( var2 ) vops.push_back( new SLiftVar( *var2 ) );
  _Interm.push_back( std::make_pair( op, vops ) );
}

template < typename DAG >
inline SLiftVar<DAG>
SLiftEnv<DAG>::_lift_univariate_term
( SLiftVar<DAG> const& var )
{
#ifdef MC__SLIFT_CHECK
  if( !_dag->curOp() )
    throw typename SLiftEnv<DAG>::Exceptions( SLiftEnv<DAG>::Exceptions::INTERNAL );
#endif

  // Append new intermediate expression and assert that same operation was not previously appended
  _append_interm( _dag->curOp(), &var );
  return SLiftVar( this, *(_dag->curOp()->pres) );
}

template < typename DAG >
inline SLiftVar<DAG>
SLiftEnv<DAG>::_lift_bivariate_term
( SLiftVar<DAG> const& var1, SLiftVar<DAG> const& var2 )
{
#ifdef MC__SLIFT_CHECK
  if( !_dag->curOp() )
    throw typename SLiftEnv<DAG>::Exceptions( SLiftEnv<DAG>::Exceptions::INTERNAL );
#endif

  // Append new intermediate expression and assert that same operation was not previously appended
  _append_interm( _dag->curOp(), &var1, &var2 );
  return SLiftVar( this, *(_dag->curOp()->pres) );
}

template < typename DAG >
inline SLiftVar<DAG>&
SLiftVar<DAG>::_set
( double const& d )
{
  _numer = d;
  _denom = 1;
  return *this;
}

template < typename DAG >
inline SLiftVar<DAG>&
SLiftVar<DAG>::_set
( FFVar const& x )
{
  _numer.var( &x );
  _denom = 1;
  return *this;
}

template < typename DAG >
inline SLiftVar<DAG>&
SLiftVar<DAG>::_set
( SLiftVar<DAG> const& var )
{
  if( this == &var ) return *this;
  _env = var._env;
  _numer = var._numer;
  _denom = var._denom;
  return *this;
}

template < typename DAG >
inline SLiftVar<DAG>
operator+
( SLiftVar<DAG> const& var )
{
  return var;
}

template < typename DAG >
inline SLiftVar<DAG>&
SLiftVar<DAG>::operator+=
( SLiftVar<DAG> const& var )
{
  if( !_env )
    _env = var._env;
  else if( var._env && _env != var._env )
    throw typename SLiftEnv<DAG>::Exceptions( SLiftEnv<DAG>::Exceptions::ENVERR );

  _numer *= var._denom;
  _numer += var._numer * _denom;
  _denom *= var._denom;
  return *this;
}

template < typename DAG >
inline SLiftVar<DAG>
operator+
( SLiftVar<DAG> const& var1, SLiftVar<DAG> const& var2 )
{
  SLiftVar<DAG> var3( var1 );
  var3 += var2;
  return var3;
}

template < typename DAG >
inline SLiftVar<DAG>
operator+
( SLiftVar<DAG> const& var1, double const& cst2 )
{
  SLiftVar<DAG> var3( var1 );
  var3 += cst2;
  return var3;
}

template < typename DAG >
inline SLiftVar<DAG>
operator-
( SLiftVar<DAG> const& var )
{
  return SLiftVar( var.env(), -var.numer(), var.denom() );
}

template < typename DAG >
inline SLiftVar<DAG>&
SLiftVar<DAG>::operator-=
( SLiftVar<DAG> const& var )
{
  if( !_env )
    _env = var._env;
  else if( var._env && _env != var._env )
    throw typename SLiftEnv<DAG>::Exceptions( SLiftEnv<DAG>::Exceptions::ENVERR );

  _numer *= var._denom;
  _numer -= var._numer * _denom;
  _denom *= var._denom;
  return *this;
}

template < typename DAG >
inline SLiftVar<DAG>
operator-
( SLiftVar<DAG> const& var1, SLiftVar<DAG> const& var2 )
{
  SLiftVar<DAG> var3( var1 );
  var3 -= var2;
  return var3;
}

template < typename DAG >
inline SLiftVar<DAG>
operator-
( SLiftVar<DAG> const& var1, double const& cst2 )
{
  SLiftVar<DAG> var3( var1 );
  var3 -= cst2;
  return var3;
}

template < typename DAG >
inline SLiftVar<DAG>&
SLiftVar<DAG>::operator*=
( SLiftVar<DAG> const& var )
{
  if( !_env )
    _env = var._env;
  else if( var._env && _env != var._env )
    throw typename SLiftEnv<DAG>::Exceptions( SLiftEnv<DAG>::Exceptions::ENVERR );

  _numer *= var._numer;
  _denom *= var._denom;
  return *this;
}

template < typename DAG >
inline SLiftVar<DAG>
operator*
( SLiftVar<DAG> const& var1, SLiftVar<DAG> const& var2 )
{
  SLiftVar<DAG> var3( var1 );
  var3 *= var2;
  return var3;
}

template < typename DAG >
inline SLiftVar<DAG>
operator*
( SLiftVar<DAG> const& var1, double const& cst2 )
{
  SLiftVar<DAG> var3( var1 );
  var3 *= cst2;
  return var3;
}

template < typename DAG >
inline SLiftVar<DAG>&
SLiftVar<DAG>::operator/=
( SLiftVar<DAG> const& var )
{
  if( !_env )
    _env = var._env;
  else if( var._env && _env != var._env )
    throw typename SLiftEnv<DAG>::Exceptions( SLiftEnv<DAG>::Exceptions::ENVERR );

  if( SLiftEnv<DAG>::options.LIFTDIV )
    *this = _env->_lift_bivariate_term( *this, var );
  else{
    _numer *= var._denom;
    _denom *= var._numer;
  }
  return *this;
}

template < typename DAG >
inline SLiftVar<DAG>
operator/
( SLiftVar<DAG> const& var1, SLiftVar<DAG> const& var2 )
{
  SLiftVar<DAG> var3( var1 );
  var3 /= var2;
  return var3;
}

template < typename DAG >
inline SLiftVar<DAG>
operator/
( double const& cst1, SLiftVar<DAG> const& var2 )
{
  SLiftVar<DAG> var3( cst1 );
  var3 /= var2;
  return var3;
}

template < typename DAG >
inline SLiftVar<DAG>
sqr
( SLiftVar<DAG> const& var )
{
#ifdef MC__SLIFT_CHECK
  if( !var.env() )
    throw typename SLiftEnv<DAG>::Exceptions( SLiftEnv<DAG>::Exceptions::ENVERR );
#endif
  if( SLiftEnv<DAG>::options.LIFTIPOW )
    return var.env()->_lift_univariate_term( var );
  return SLiftVar<DAG>( var.env(), sqr( var.numer() ), sqr( var.denom() ) );
}

template < typename DAG >
inline SLiftVar<DAG>
inv
( SLiftVar<DAG> const& var )
{
#ifdef MC__SLIFT_CHECK
  if( !var.env() )
    throw typename SLiftEnv<DAG>::Exceptions( SLiftEnv<DAG>::Exceptions::ENVERR );
#endif
  if( SLiftEnv<DAG>::options.LIFTDIV )
    return var.env()->_lift_bivariate_term( 1, var );
  return SLiftVar<DAG>( var.env(), var.denom(), var.numer() );
}

template < typename DAG >
inline SLiftVar<DAG>
pow
( SLiftVar<DAG> const& var, double const& a )
{
#ifdef MC__SLIFT_CHECK
  if( !var.env() )
    throw typename SLiftEnv<DAG>::Exceptions( SLiftEnv<DAG>::Exceptions::ENVERR );
#endif
  return var.env()->_lift_bivariate_term( var, a );
}

template < typename DAG >
inline SLiftVar<DAG>
pow
( SLiftVar<DAG> const& var, int const n )
{
  if( n < 0 ) return pow( inv( var ), -n );
  switch( n ){
   case 0:  return 1.;
   case 1:  return var;
   case 2:  return sqr( var );
   default: if( SLiftEnv<DAG>::options.LIFTIPOW ) return var.env()->_lift_bivariate_term( var, n );
            return SLiftVar( var.env(), pow( var.numer(), (unsigned)n ), pow( var.denom(), (unsigned)n ) );
  }
}

template < typename DAG >
inline SLiftVar<DAG>
cheb
( SLiftVar<DAG> const& var, unsigned const n )
{
  switch( n ){
   case 0:  return 1.;
   case 1:  return var;
   case 2:  return sqr( var ) * 2. - 1.;
   default: if( SLiftEnv<DAG>::options.LIFTIPOW ) return var.env()->_lift_bivariate_term( var, n );
            return var * cheb( var, n-1 ) * 2. - cheb( var, n-2 );
  }
}

template < typename DAG >
inline SLiftVar<DAG>
prod
( unsigned int const nvars, SLiftVar<DAG> const* pvars )
{
  switch( nvars ){
   case 0:  return 1.;
   case 1:  return pvars[0];
   default: return pvars[0] * prod( nvars-1, pvars+1 );
  }
}

template < typename DAG >
inline SLiftVar<DAG>
monom
( unsigned int const nvars, SLiftVar<DAG> const* pvars, unsigned const* k, bool const chebbasis=false )
{
  switch( nvars ){
   case 0:  return 1.;
   case 1:  return chebbasis? cheb( pvars[0], k[0] ): pow( pvars[0], (int)k[0] );
   default: return ( chebbasis? cheb( pvars[0], k[0] ): pow( pvars[0], (int)k[0] ) ) * monom( nvars-1, pvars+1, k+1 );
  }
}

template < typename DAG >
inline SLiftVar<DAG>
exp
( SLiftVar<DAG> const& var )
{
#ifdef MC__SLIFT_CHECK
  if( !var.env() )
    throw typename SLiftEnv<DAG>::Exceptions( SLiftEnv<DAG>::Exceptions::ENVERR );
#endif
  return var.env()->_lift_univariate_term( var );
}

template < typename DAG >
inline SLiftVar<DAG>
log
( SLiftVar<DAG> const& var )
{
#ifdef MC__SLIFT_CHECK
  if( !var.env() )
    throw typename SLiftEnv<DAG>::Exceptions( SLiftEnv<DAG>::Exceptions::ENVERR );
#endif
  return var.env()->_lift_univariate_term( var );
}

template < typename DAG >
inline SLiftVar<DAG>
xlog
( SLiftVar<DAG> const& var )
{
#ifdef MC__SLIFT_CHECK
  if( !var.env() )
    throw typename SLiftEnv<DAG>::Exceptions( SLiftEnv<DAG>::Exceptions::ENVERR );
#endif
  return var.env()->_lift_univariate_term( var );
}

template < typename DAG >
inline SLiftVar<DAG>
sqrt
( SLiftVar<DAG> const& var )
{
#ifdef MC__SLIFT_CHECK
  if( !var.env() )
    throw typename SLiftEnv<DAG>::Exceptions( SLiftEnv<DAG>::Exceptions::ENVERR );
#endif
  return var.env()->_lift_univariate_term( var );
}

template < typename DAG >
inline SLiftVar<DAG>
cos
( SLiftVar<DAG> const& var )
{
#ifdef MC__SLIFT_CHECK
  if( !var.env() )
    throw typename SLiftEnv<DAG>::Exceptions( SLiftEnv<DAG>::Exceptions::ENVERR );
#endif
  return var.env()->_lift_univariate_term( var );
}

template < typename DAG >
inline SLiftVar<DAG>
sin
( SLiftVar<DAG> const& var )
{
#ifdef MC__SLIFT_CHECK
  if( !var.env() )
    throw typename SLiftEnv<DAG>::Exceptions( SLiftEnv<DAG>::Exceptions::ENVERR );
#endif
  return var.env()->_lift_univariate_term( var );
}

template < typename DAG >
inline SLiftVar<DAG>
tan
( SLiftVar<DAG> const& var )
{
#ifdef MC__SLIFT_CHECK
  if( !var.env() )
    throw typename SLiftEnv<DAG>::Exceptions( SLiftEnv<DAG>::Exceptions::ENVERR );
#endif
  return var.env()->_lift_univariate_term( var );
}

template < typename DAG >
inline SLiftVar<DAG>
acos
( SLiftVar<DAG> const& var )
{
#ifdef MC__SLIFT_CHECK
  if( !var.env() )
    throw typename SLiftEnv<DAG>::Exceptions( SLiftEnv<DAG>::Exceptions::ENVERR );
#endif
  return var.env()->_lift_univariate_term( var );
}

template < typename DAG >
inline SLiftVar<DAG>
asin
( SLiftVar<DAG> const& var )
{
#ifdef MC__SLIFT_CHECK
  if( !var.env() )
    throw typename SLiftEnv<DAG>::Exceptions( SLiftEnv<DAG>::Exceptions::ENVERR );
#endif
  return var.env()->_lift_univariate_term( var );
}

template < typename DAG >
inline SLiftVar<DAG>
atan
( SLiftVar<DAG> const& var )
{
#ifdef MC__SLIFT_CHECK
  if( !var.env() )
    throw typename SLiftEnv<DAG>::Exceptions( SLiftEnv<DAG>::Exceptions::ENVERR );
#endif
  return var.env()->_lift_univariate_term( var );
}

template < typename DAG >
inline SLiftVar<DAG>
cosh
( SLiftVar<DAG> const& var )
{
#ifdef MC__SLIFT_CHECK
  if( !var.env() )
    throw typename SLiftEnv<DAG>::Exceptions( SLiftEnv<DAG>::Exceptions::ENVERR );
#endif
  return var.env()->_lift_univariate_term( var );
}

template < typename DAG >
inline SLiftVar<DAG>
sinh
( SLiftVar<DAG> const& var )
{
#ifdef MC__SLIFT_CHECK
  if( !var.env() )
    throw typename SLiftEnv<DAG>::Exceptions( SLiftEnv<DAG>::Exceptions::ENVERR );
#endif
  return var.env()->_lift_univariate_term( var );
}

template < typename DAG >
inline SLiftVar<DAG>
tanh
( SLiftVar<DAG> const& var )
{
#ifdef MC__SLIFT_CHECK
  if( !var.env() )
    throw typename SLiftEnv<DAG>::Exceptions( SLiftEnv<DAG>::Exceptions::ENVERR );
#endif
  return var.env()->_lift_univariate_term( var );
}

template < typename DAG >
inline SLiftVar<DAG>
fabs
( SLiftVar<DAG> const& var )
{
#ifdef MC__SLIFT_CHECK
  if( !var.env() )
    throw typename SLiftEnv<DAG>::Exceptions( SLiftEnv<DAG>::Exceptions::ENVERR );
#endif
  return var.env()->_lift_univariate_term( var );
}

template < typename DAG >
inline SLiftVar<DAG>
erf
( SLiftVar<DAG> const& var )
{
#ifdef MC__SLIFT_CHECK
  if( !var.env() )
    throw typename SLiftEnv<DAG>::Exceptions( SLiftEnv<DAG>::Exceptions::ENVERR );
#endif
  return var.env()->_lift_univariate_term( var );
}

template < typename DAG >
inline SLiftVar<DAG>
fstep
( SLiftVar<DAG> const& var )
{
#ifdef MC__SLIFT_CHECK
  if( !var.env() )
    throw typename SLiftEnv<DAG>::Exceptions( SLiftEnv<DAG>::Exceptions::ENVERR );
#endif
  return var.env()->_lift_univariate_term( var );
}

template < typename DAG >
inline SLiftVar<DAG>
max
( SLiftVar<DAG> const& var1, SLiftVar<DAG> const& var2 )
{
  if( var1.env() && var2.env() && var1.env() != var2.env() )
    throw typename SLiftEnv<DAG>::Exceptions( SLiftEnv<DAG>::Exceptions::ENVERR );

  return var1.env()->_lift_bivariate_term( var1, var2 );
}

template < typename DAG >
inline SLiftVar<DAG>
min
( SLiftVar<DAG> const& var1, SLiftVar<DAG> const& var2 )
{
  if( var1.env() && var2.env() && var1.env() != var2.env() )
    throw typename SLiftEnv<DAG>::Exceptions( SLiftEnv<DAG>::Exceptions::ENVERR );

  return var1.env()->_lift_bivariate_term( var1, var2 );
}

template < typename DAG >
inline SLiftVar<DAG>
lmtd
( SLiftVar<DAG> const& var1, SLiftVar<DAG> const& var2 )
{
  if( var1.env() && var2.env() && var1.env() != var2.env() )
    throw typename SLiftEnv<DAG>::Exceptions( SLiftEnv<DAG>::Exceptions::ENVERR );

  return var1.env()->_lift_bivariate_term( var1, var2 );
}

template < typename DAG >
inline SLiftVar<DAG>
rlmtd
( SLiftVar<DAG> const& var1, SLiftVar<DAG> const& var2 )
{
  if( var1.env() && var2.env() && var1.env() != var2.env() )
    throw typename SLiftEnv<DAG>::Exceptions( SLiftEnv<DAG>::Exceptions::ENVERR );

  return var1.env()->_lift_bivariate_term( var1, var2 );
}

} // namespace mc

//#include "mcop.hpp"

namespace mc
{

//! @brief Specialization of the structure mc::Op to allow usage of the type mc::SLiftVar for DAG evaluation or as a template parameter in other MC++ classes
template <typename DAG> struct Op<mc::SLiftVar<DAG>>
{
  typedef mc::SLiftVar<DAG> T;
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
