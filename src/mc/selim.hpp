// Copyright (C) Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the Eclipse Public License.

/*!
\page page_SELIM Variable Elimination in Factorable Expressions
\author Benoit Chachuat & Dominik Bongartz
\date 2023
\bug No known bugs.

The classes mc::SElimEnv and mc::SElimVar defined in <tt>selim.hpp</tt> enable the elimination of a subset of variables from factorable expressions through expoiting solvable (invertible) equality constraints. 

Given a set of equality constraints \f${\bf f}({\bf x}) = {\bf 0}\f$, we week a partition of the variable set \f${\bf x} =: [{\bf x},  = {\bf 0}\f$
The elimination proceeds in 3 steps:
-# Identify which equality constraints can be inverted analytically for which variables.

-# Determine which variables to eliminate using which equality constraints in order to maximize the number of eliminated variables, based on: (i) a mixed-integer programming (MIP) formulation; or (ii) a greedy heuristic inspired by the method of Hernandez & Sargent (1979) to reorder the equations then select the variables through reordering the incidence matrix to bordered lower-triangular form.

-# Construct expressions for the eliminated variables, so they can be eliminated through composition.
.

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
#include "slift.hpp"

#if defined(MC__USE_GUROBI)
 #include "gurobi_c++.h"
 extern "C"{
  #include <fenv.h>
  int fedisableexcept( int );
 }
#endif

#define MC__SELIM_CHECK
//#define MC__SELIM_DEBUG_PROCESS
//#define MC__SELIM_DEBUG_INSERT

namespace mc
{

//! @brief Environment for variable elimination in factorable equality constraints
////////////////////////////////////////////////////////////////////////
//! mc::SElimEnv is a C++ class defining the environment for variable
//! elimination from factorable expressions through expoiting solvable
//! equality constraints
////////////////////////////////////////////////////////////////////////
template <typename... ExtOps>
class SElimEnv:
  protected virtual SLiftEnv<ExtOps...>
////////////////////////////////////////////////////////////////////////
{
  using SLiftEnv<ExtOps...>::_dag;
  using SLiftEnv<ExtOps...>::_Interm;
  using SLiftEnv<ExtOps...>::_SPDep;

  using SLiftEnv<ExtOps...>::dag;
  using SLiftEnv<ExtOps...>::insert_dag;

  template <typename... Ops> friend std::ostream& operator<<
    ( std::ostream&, const SElimEnv<Ops...>& );

public:

  typedef std::tuple< std::vector<FFVar>, std::vector<FFVar>, std::vector<FFVar> > t_VarElim;
  typedef SLiftEnv<ExtOps...> t_lift;
  typedef typename SLiftVar::t_poly t_poly;
  typedef typename t_poly::t_mon t_mon;

  //! @brief Default Constructor
  SElimEnv
    ( FFGraph<ExtOps...>* dag=nullptr )
    : SLiftEnv<ExtOps...>(dag)
    {
#if defined(MC__USE_GUROBI)
      _GRBenv   = new GRBEnv();
      _GRBmodel = nullptr;
#endif
    }
    
  //! @brief Destructor
  virtual ~SElimEnv
    ()
    {
      _reset();
#if defined(MC__USE_GUROBI)
      delete _GRBmodel;
      delete _GRBenv;
#endif
    }

  //! @brief Retreive tuple of vectors <ELIMINATED VARIABLE, INVERTED CONSTRAINT, INVERTED EXPRESSION> with entries in a feasible order of elimination
  t_VarElim& VarElim
    ()
    { return _VarElim; }

  //! @brief Set DAG environment
  void set
    (  FFGraph<ExtOps...>* dag )
    { SLiftEnv<ExtOps...>::set( dag );
      _reset(); }

  //! @brief Reset intermediate expressions
  void reset
    ()
    { SLiftEnv<ExtOps...>::_reset();
      _reset(); }

  //! @brief Process the <a>ndxCtr</a> equality constraints in array <a>pDep</a> 
  void process
    ( std::set<unsigned> const& ndxCtr, FFVar const* pCtr, 
      std::map<FFVar const*,double,lt_FFVar> const& wVar=std::map<FFVar const*,double,lt_FFVar>(),
      bool const add2dag=true );

  //! @brief Process the <a>nCtr</a> equality constraints in array <a>pDep</a>
  void process
    ( unsigned const nCtr, FFVar const* pCtr, 
      std::map<FFVar const*,double,lt_FFVar> const& wVar=std::map<FFVar const*,double,lt_FFVar>(),
      const bool add2dag=true );

  //! @brief Exceptions of mc::SElimVar
  class Exceptions
  {
   public:
    //! @brief Enumeration type for SElimVar exception handling
    enum TYPE{
      MIPERR=0,       //!< Call to MIP solver disabled
      INVERT,         //!< Elimnation process failed
      INTERNAL=-33    //!< Internal error
    };
    //! @brief Constructor for error <a>ierr</a>
    Exceptions( TYPE ierr ) : _ierr( ierr ){}
    //! @brief Error flag
    int ierr(){ return _ierr; }
    //! @brief Error description
    std::string what(){
      switch( _ierr ){
      case MIPERR:
        return "mc::SElimEnv\t Mixed-integer programming solver is disabled";
      case INVERT:
        return "mc::SElimEnv\t Internal error during elimination process";
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
#if defined(MC__USE_GUROBI)
      LPALGO( LPALGO_DEFAULT ), LPPRESOLVE(-1),
      LPFEASTOL(1e-7), LPOPTIMTOL(1e-7),
      MIPRELGAP(1e-3), MIPABSGAP(1e-3), MIPTHREADS(0),
      MIPCONCURRENT(1), MIPFOCUS(0), MIPHEURISTICS(0.2),
      MIPNUMFOCUS(0), MIPDISPLEVEL(1), MIPOUTPUTFILE(""),
      MIPTIMELIMIT(600),
#endif
      MULTMAX(0), ELIMLIN(true), ELIMMLIN(true),
      ELIMNLIN( {FFInv::Options::INV,FFInv::Options::SQRT,FFInv::Options::EXP,
                  FFInv::Options::LOG,FFInv::Options::RPOW} )
      {}
    //! @brief Assignment of mc::SElimEnv<ExtOps...>::Options
    Options& operator=
      ( Options& opt ){
#if defined(MC__USE_GUROBI)
        LPALGO          = opt.LPALGO;
        LPPRESOLVE      = opt.LPPRESOLVE;
        LPFEASTOL       = opt.LPFEASTOL;
        LPOPTIMTOL      = opt.LPOPTIMTOL;
        MIPRELGAP       = opt.MIPRELGAP;
        MIPABSGAP       = opt.MIPABSGAP;
        MIPTHREADS      = opt.MIPTHREADS;
        MIPCONCURRENT   = opt.MIPCONCURRENT;
        MIPFOCUS        = opt.MIPFOCUS;
        MIPHEURISTICS   = opt.MIPHEURISTICS;
        MIPNUMFOCUS     = opt.MIPNUMFOCUS;
        MIPDISPLEVEL    = opt.MIPDISPLEVEL;
        MIPOUTPUTFILE   = opt.MIPOUTPUTFILE;
        MIPTIMELIMIT    = opt.MIPTIMELIMIT;
#endif
        MULTMAX  = opt.MULTMAX;
        ELIMLIN  = opt.ELIMLIN;
        ELIMMLIN = opt.ELIMMLIN;
        ELIMNLIN = opt.ELIMNLIN;
        return *this;
      }
#if defined(MC__USE_GUROBI)
    //! @brief LP algorithm used by MIP solver
    int LPALGO;
    //! @brief LP presolve strategy in MIP solver
    int LPPRESOLVE;
    //! @brief Tolerance on LP feasibility in MIP solver
    double LPFEASTOL;
     //! @brief Tolerance on LP optimality in MIP solver
    double LPOPTIMTOL;
    //! @brief Tolerance on relative gap in MIP solver
    double MIPRELGAP;
    //! @brief Tolerance on absolute gap in MIP solver
    double MIPABSGAP;
    //! @brief Number of threads used by MIP solver - default value of 0 allows to use all available threads
    int MIPTHREADS;
    //! @brief Number of independent MIP solves in parallel
    int MIPCONCURRENT;
    //! @brief MIP high-level solution strategy - default value of 0 seeks to a balance between finding new feasible solutions and proving that the current solution is optimal
    int MIPFOCUS;
    //! @brief Fraction of time spent in MIP heuristics - default value of 0.2 aims to spend 20% of runtime on heuristics
    double MIPHEURISTICS;
    //! @brief Degree to which the code attempts to detect and manage numerical issues - default setting of 0 makes an automatic choice, with a slight preference for speed 
    int MIPNUMFOCUS;
    //! @brief Display level for MIP solver
    int MIPDISPLEVEL;
    //! @brief Name of output file for MIP model
    std::string MIPOUTPUTFILE;
    //! @brief Maximum MIP run time (seconds)
    double MIPTIMELIMIT;
  //! @brief Default option for LP solver
    static const int LPALGO_DEFAULT = -1;
#endif
    //! @brief Maximal multiplicity of eliminated variables
    unsigned MULTMAX;
    //! @brief Whether to invert linear operations
    bool ELIMLIN;
    //! @brief Whether to invert multilinear operations
    bool ELIMMLIN;
    //! @brief Set of invertible nonlinear operations
    std::set<FFInv::Options::NLINV> ELIMNLIN;
  } options;


protected:

  //! @brief Map between eliminated DAG variables and coresponding equality constraint
  t_VarElim _VarElim;

  //! @brief Independent DAG variables participating in equality constraints
  std::vector<FFVar> _Var;

  //! @brief Independent DAG variables with associated weights/priorities
  std::map<unsigned,double> _VarWeight;

  //! @brief Candidate variables
  std::set<unsigned> _ndxVar;
  
  //! @brief Variable multiplicity
  std::map<unsigned,unsigned> _mapVarMult;

  //! @brief Candidate constraints from which to eliminate a given variable
  std::map<unsigned,std::set<unsigned>> _mapVarCtrCand;

  //! @brief Candidate variables to eliminate a given variable
  std::map<unsigned,std::set<unsigned>> _mapVarElVar;

  //! @brief Equality constraints
  std::vector<FFVar> _Ctr;

  //! @brief Candidate constraints
  std::set<unsigned> _ndxCtr;

  //! @brief Participating variables in a given constraint
  std::map<unsigned,std::set<unsigned>> _mapCtrVar;

  //! @brief Candidate variables to eliminate from a given constraint
  std::map<unsigned,std::set<unsigned>> _mapCtrVarCand;

  //! @brief Invertibility information for DAG variables
  std::vector<FFInv> _IVar;

  //! @brief Invertibility information for equality constraints
  std::vector<FFInv> _ICtr;
  
#if defined(MC__USE_GUROBI)
  //! @brief whether the MIP solver has sent an exception
  bool _MIPexcpt;
  
  //! @brief MIP environment
  GRBEnv* _GRBenv;
  
  //! @brief MIP model
  GRBModel* _GRBmodel;
  
  //! @brief Binary variables indicating selected equality constraints
  std::map<unsigned,GRBVar> _MIP_ctr;
  
  //! @brief Binary variables indicating eliminated variables
  std::map<unsigned,GRBVar> _MIP_var;
  
  //! @brief Continuous variables for MTZ subtour elimination
  std::map<unsigned,GRBVar> _MIP_varMTZ;

  //! @brief Binary variables matching equality constraints and eliminated variables
  std::map<unsigned,std::map<unsigned,GRBVar>> _MIP_ctrvar;

  //! @brief Binary variables matching eliminated variable with other participating variables
  std::map<unsigned,std::map<unsigned,GRBVar>> _MIP_varvar;

  //! @brief Optimize the quadratic expressions for minimum number of auxiliaries 
  void _MIP_optimize
    ();

  //! @brief Solve MIP optimization model
  void _MIP_solve
    ();

  //! @brief Encode MIP optimization model for maximal elimination
  void _MIP_encode
    ();

  //! @brief Decode MIP optimal solution for maximal elimination
  void _MIP_decode
    ();

  //! @brief Reset variable vectors in MIP optimization model
  void _MIP_reset
    ();

  //! @brief Set options in MIP optimization model
  void _MIP_options
    ();
#endif

  //! @brief Insert an auxiliary variable corresponding to operand of inverse operation into DAG
  FFVar const* _insert_expr
    ( FFVar& var, t_poly const& polyindep, t_poly const& polydep, bool const useprod,
      int const order );

  //! @brief Insert an auxiliary variable corresponding to inverse operation into DAG
  std::pair< FFVar const*, SLiftVar const* > _insert_expr
    ( long const ndxVarEl, FFVar& var, std::vector<SLiftVar const*> const& SPVar,
      FFOp const* pOp, t_poly const& polyindep, t_poly const& polydep, bool const useprod,
      int const order );

  //! @brief Check dependency in variable ndxVarEl in a lifted variable
  std::set< FFVar const*, lt_FFVar > _dep_expr
    ( long const ndxVarEl, t_poly const& numer, t_poly const& denom );

  //! @brief Return pointer to intrenal DAG variable of expression
  FFVar const* _ptr_expr
    ( FFVar& var );

  //! @brief Erase all entries in _Interm
  void _reset
    ();
};

template <typename... ExtOps>
inline typename SElimEnv<ExtOps...>::Options SElimEnv<ExtOps...>::options;

#if defined(MC__USE_GUROBI)
template <typename... ExtOps>
inline int const SElimEnv<ExtOps...>::Options::LPALGO_DEFAULT;
#endif

template <typename... ExtOps>
inline std::ostream&
operator<<
( std::ostream& out, SElimEnv<ExtOps...> const& env)
{
  auto const& [vvar,vctr,vaux] = env._VarElim;
  std::cout << std::endl
            << vvar.size() << " VARIABLES MAY BE ELIMINATED" << std::endl;

  if( vaux.empty() ){
    for( auto itvar = vvar.begin(), itctr = vctr.begin();
         itvar != vvar.end(); ++itvar, ++itctr )
      out << "VARIABLE " << *itvar << " USING CONSTRAINT" << *itctr << "=0";
    return out;
  }

  for( auto itvar = vvar.begin(), itctr = vctr.begin(), itaux = vaux.begin();
       itvar != vvar.end(); ++itvar, ++itctr, ++itaux ){
    std::ostringstream ext; 
    ext << " OF " << *itvar << " USING " << *itctr << "=0";
    env._dag->output( env._dag->subgraph( 1, &*itaux ), ext.str(), out );
  }
  return out;
}

template <typename... ExtOps>
inline void
SElimEnv<ExtOps...>::_reset
()
{
  _Var.clear();
  _IVar.clear();
  _VarWeight.clear();
  _ndxVar.clear();
  _mapVarCtrCand.clear();
  _mapVarElVar.clear();
  _mapVarMult.clear();

  _Ctr.clear();
  _ICtr.clear();
  _ndxCtr.clear();
  _mapCtrVar.clear();
  _mapCtrVarCand.clear();

  std::get<0>( _VarElim ).clear();
  std::get<1>( _VarElim ).clear();
  std::get<2>( _VarElim ).clear();
}

template <typename... ExtOps>
inline void
SElimEnv<ExtOps...>::process
( std::set<unsigned> const& ndxCtr, FFVar const* pCtr,
  std::map<FFVar const*,double,lt_FFVar> const& wVar,
  bool const add2dag )
{
  if( ndxCtr.empty() ) return; // Nothing to do!
  std::vector<FFVar> vCtr;
  vCtr.reserve( ndxCtr.size() );
  for( unsigned const& i : ndxCtr ) vCtr.push_back( pCtr[i] );
  process( ndxCtr.size(), vCtr.data(), wVar, add2dag );
}

template <typename... ExtOps>
inline void
SElimEnv<ExtOps...>::process
( unsigned const nCtr, FFVar const* pCtr,
  std::map<FFVar const*,double,lt_FFVar> const& wVar,
  bool const add2dag )
{
  // Reset internal variables
  _reset();

  // Update participating variables in _Var
  auto sgCtr = _dag->subgraph( nCtr, pCtr );
  unsigned v=0;
  for( auto const& Op : sgCtr.l_op ){
    if( Op->type != FFOp::VAR ) continue;
    _Var.push_back( *Op->varout[0] );
    _IVar.push_back( FFInv().indep( v ) );
    auto itv = wVar.find( Op->varout[0] ); 
    _VarWeight[v] = ( itv != wVar.end() ? itv->second : 1. );
    ++v;
  }
#ifdef MC__SELIM_DEBUG_PROCESS
  std::cout << std::endl << _Var.size() << " Original Variables: ";
  for( auto const& var : _Var ) std::cout << var << " ";
  std::cout << std::endl;
#endif

  // Find out variable invertibility in equality constraints
  _Ctr.assign( pCtr, pCtr+nCtr );
  _ICtr.resize( nCtr );
  FFInv::options.INVOP = options.ELIMNLIN;
  //std::cout << "FFInv::options.INVOP: " << FFInv::options.INVOP.size() << std::endl;
  _dag->eval( sgCtr, nCtr, pCtr, _ICtr.data(), _Var.size(), _Var.data(), _IVar.data() );

  // Create variable multiplicity map !!!THIS IS NOT TRULY MULTIPLICITY - MULTIPLE OCCURENCES COULD BE WITHIN SINGLE FUNCTION
  for( unsigned j=0; j<_ICtr.size(); ++j )
    for( auto const& [i,type] : _ICtr[j].inv() )
      if( _mapVarMult.count(i) ) _mapVarMult[i]++;
      else                       _mapVarMult[i] = 1;
#ifdef MC__SELIM_DEBUG_PROCESS
  std::cout << std::endl << "Participating variable multiplicity:" << std::endl;
  for( auto const& [i,m] : _mapVarMult ){
    std::cout << i << " (" << _Var[i] << "):" << m;
    std::cout << std::endl;
  }
#endif

  // Create candidate variable and constraint sets and maps
  for( unsigned j=0; j<_ICtr.size(); ++j ){
#ifdef MC__SELIM_DEBUG_PROCESS
    std::cout << j << ": " << _ICtr[j] << std::endl;
#endif
    bool isinvert = false;
    for( auto const& [i,type] : _ICtr[j].inv() ){
      // Exclude variables based on weight and/or multiplicity
      if( _VarWeight[i] <= 0.
       || (options.MULTMAX > 0 && _mapVarMult[i] > options.MULTMAX) ) continue;
      switch( type ){
        case FFInv::L:
          if( options.ELIMLIN ){
	    _ndxVar.insert( i );
	    isinvert = true;
	  }
          break;
	case FFInv::S:
          if( options.ELIMMLIN ){
	    _ndxVar.insert( i );
	    isinvert = true;
	  }
          break;
	case FFInv::N:
          _ndxVar.insert( i );
	  isinvert = true;
          break;
	case FFInv::U:
        default:
	  break;
      }
    }
    if( isinvert ) _ndxCtr.insert( j );
  }
#ifdef MC__SELIM_DEBUG_PROCESS
  std::cout << std::endl << "Candidate variables for elimination:";
  for( unsigned const& v : _ndxVar ) std::cout << " " << v;
  std::cout << std::endl;
  std::cout << std::endl << "Candidate constraints for elimination:";
  for( unsigned const& e : _ndxCtr ) std::cout << " " << e;
  std::cout << std::endl;
#endif

  for( unsigned const& i : _ndxVar )
    _mapVarElVar[i] = std::set<unsigned>();
  for( unsigned const& j : _ndxCtr ){
    for( auto const& [i,type] : _ICtr[j].inv() ){
      if( _ndxVar.find( i ) == _ndxVar.end() ) continue;
      _mapCtrVar[j].insert( i );
      bool insert = false;
      switch( type ){
        case FFInv::L: if( options.ELIMLIN )  insert = true; break;
	case FFInv::S: if( options.ELIMMLIN ) insert = true; break;
	case FFInv::N: insert = true; break;
	case FFInv::U: default: break;
      }
      if( insert ){
        _mapCtrVarCand[j].insert( i );
        _mapVarCtrCand[i].insert( j );
	for( auto const& [k,type] : _ICtr[j].inv() )
	  if( k != i && _ndxVar.find( k ) != _ndxVar.end() )
	    _mapVarElVar[i].insert( k );
      }
    }
  }
#ifdef MC__SELIM_DEBUG_PROCESS
  std::cout << std::endl << "Participating variables in each constraint:" << std::endl;
  for( auto const& [e,ndx] : _mapCtrVar ){
    std::cout << e << " (" << _Ctr[e] << "):";
    for( auto const& v : ndx ) std::cout << " " << v;
    std::cout << std::endl;
  }
  std::cout << std::endl << "Candidate variables in each constraint:" << std::endl;
  for( auto const& [j,ndx] : _mapCtrVarCand ){
    std::cout << j << " (" << _Ctr[j] << "):";
    for( auto const& i : ndx ) std::cout << " " << i;
    std::cout << std::endl;
  }
  std::cout << std::endl << "Candidate variables to substitute each variable:" << std::endl;
  for( auto const& [i,ndx] : _mapVarElVar ){
    std::cout << i << " (" << _Var[i] << "):";
    for( auto const& k : ndx ) std::cout << " " << k << " (" << _Var[k] << ")";
    std::cout << std::endl;
  }
  std::cout << std::endl << "Candidate constraints for variable elimination:" << std::endl;
  for( auto const& [i,ndx] : _mapVarCtrCand ){
    std::cout << i << " (" << _Var[i] << "):";
    for( auto const& j : ndx ) std::cout << " " << j;
    std::cout << std::endl;
  }
#endif

  // Determine an optimal elimination set
#if defined(MC__USE_GUROBI)
  _MIP_optimize();
#else
  throw Exceptions( Exceptions::MIPERR );
#endif

  auto& [vvar,vctr,vaux] = _VarElim;
#ifdef MC__SELIM_DEBUG_PROCESS
  std::cout << std::endl << "Eliminated variables and corresponding constraints:" << std::endl;
  for( auto itvar = vvar.cbegin(), itctr = vctr.cbegin();
       itvar != vvar.cend(); ++itvar, ++itctr )
    std::cout << *itvar << " <- " << *itctr << std::endl;
#endif
  
  // No transcription in DAG if <a>add2dag</a> is false
  if( !add2dag ) return;

  t_poly::options.BASIS = t_poly::Options::MONOM;
  t_lift::options.LIFTDIV  = true;
  t_lift::options.LIFTIPOW = false;


#ifdef MC__SELIM_CHECK
  assert( vvar.size() == vctr.size() );
#endif
  vaux.resize( vvar.size() );
  for( auto itvarel = vvar.begin(), itctr = vctr.begin(), itaux = vaux.begin();
       itvarel != vvar.end(); ++itvarel, ++itctr, ++itaux ){
    t_lift::process( 1, &*itctr, false );
#ifdef MC__SELIM_DEBUG_PROCESS
    std::cout << *(t_lift*)this;
#endif
    // Initialize lifted variable to current dependent
    auto const* spvar = &_SPDep.front();
#ifdef MC__SELIM_DEBUG_PROCESS
      std::cout << "First dependent variable: " << *spvar << std::endl;
#endif

    // Initialize lifted variable and loop over intermediate operations
    for( FFVar VarEl = 0.; ; ){

      // Check active dependence of spvar
      auto const& denom = spvar->denom();
      auto const& numer = spvar->numer();
      long ndxVarEl = itvarel->id().second;
      auto&& sdep = _dep_expr( ndxVarEl, numer, denom );
      //std::cout << "HERE2\n";
      
      if( sdep.empty() )
        throw Exceptions( Exceptions::INVERT );

      // Case multiple dependents
      else if( sdep.size() > 1 ){
        // Check all dependences are separably linear
        t_mon monvarel( &*itvarel );
	std::vector<t_mon> vmondep; vmondep.reserve( sdep.size() );
	for( auto pdep : sdep ) vmondep.push_back( t_mon( pdep ) );
        t_poly polyindep;
        FFVar funcdep, funcindep;
        for( auto const& [mon,coef] : numer.mapmon() ){
	  auto itdep = sdep.cbegin();
	  for( auto const& mondep : vmondep ){
            if( !mondep.subseteq( mon ) ){
	      ++itdep;
	      continue;
	    }

	    // Eliminated variable directly participates in monomial
	    if( monvarel.subseteq( mon ) )
              funcdep += insert_dag( mon - monvarel, false, true ) * coef;

	    // Eliminated variable indirectly participates in monomial
            else{
              bool found_Op = false;
              for( auto const& [pOp,SPVar] : _Interm ){
                if( pOp->varout[0]->id().first  != FFVar::AUX ) throw Exceptions( Exceptions::INVERT );
                if( pOp->varout[0]->id().second != (*itdep)->id().second ) continue;
		//std::cout << "OPERATION: " << *pOp << std::endl;
		if( pOp->type != FFOp::DIV ) throw Exceptions( Exceptions::INVERT );
                assert( SPVar.size() == 2 ); // exactly two operands
                found_Op = true;
		// Check denominator is independent
                int ndep2 = _dep_expr( ndxVarEl, SPVar.at(1)->numer(), SPVar.at(1)->denom() ).size();
		assert( ndep2 == 0 );
		// Check numerator is dependent
                int ndep1 = _dep_expr( ndxVarEl, SPVar.at(0)->numer(), SPVar.at(0)->denom() ).size();
		assert( ndep1 == 1 );
                // Split numerator into dependent and independent parts
                t_poly polydep, polyindep;
                for( auto const& [monaux,coefaux] : SPVar.at(0)->numer().mapmon() ){
                  if( monvarel.subseteq( monaux ) )
                    polydep += std::make_pair( monaux - monvarel, coefaux );
                  else
                    polyindep -= std::make_pair( mon, coef );
                }
                // Append dependent and independent parts
		FFVar vardenom = insert_dag( SPVar.at(1)->numer(), false, true );
                if( !polydep.mapmon().empty() )
                  funcdep += insert_dag( polydep, false, true ) / vardenom;
                if( !polyindep.mapmon().empty() )
                  funcindep -= insert_dag( polyindep, false, true ) / vardenom;
		break;
	      }
              if( !found_Op ) throw Exceptions( Exceptions::INVERT );

	    }
	    break;
	  }
	  if( itdep == sdep.cend() )
            funcindep -= insert_dag( mon, false, true ) * coef;
        }

        // Update eliminated variable expression and return pointer to internal DAG variable
        VarEl += funcindep;
        VarEl /= funcdep;
        *itaux = *_ptr_expr( VarEl );
        break; // EXIT THE INNER FOR LOOP ON INTERMEDIATE VARIABLES	
      }
      
      // Split numerator into dependent and independent parts
      FFVar const* pdep = *sdep.cbegin();
      t_mon mondep( pdep );
      unsigned orddep = 0;
      t_poly polydep, polyindep;
      for( auto const& [mon,coef] : numer.mapmon() ){
        // Detect if variable has order >1 in monomial
        t_mon mondepord;
        for( ; (mondepord+mondep).subseteq( mon ); orddep++ ) mondepord += mondep;
        if( mondepord.tord )
          polydep   += std::make_pair( mon - mondepord, coef );
        else
          polyindep -= std::make_pair( mon, coef );
      }
#ifdef MC__SELIM_DEBUG_PROCESS
      std::cout << "Dependent variable: " << *pdep << "^" << orddep << std::endl;
      std::cout << "Dependent part:" << polydep;
      std::cout << "Independent part:" << polyindep;
#endif

      // Case unique dependent is a DAG variable
      if( pdep->id().first == FFVar::VAR ){
        // Create new DAG expression and copy pointer to expression in _VarElim
        *itaux = *_insert_expr( VarEl, polyindep, polydep, false, orddep );
#ifdef MC__SELIM_DEBUG_PROCESS
        std::cout << "INSERTED DAG EXPRESSION:";
        std::ostringstream ext; 
        ext << " OF " << *pdep;
        _dag->output( _dag->subgraph( 1, &*itaux ), ext.str() );
#endif
        break; // EXIT THE INNER FOR LOOP ON INTERMEDIATE VARIABLES
      }

      // Case dependent is not a DAG variable, identify intermediate operation
      bool found_Op = false;
      for( auto const& [pOp,SPVar] : _Interm ){
        if( pOp->varout[0]->id().first  != FFVar::AUX ) throw Exceptions( Exceptions::INVERT );
        if( pOp->varout[0]->id().second != pdep->id().second ) continue;
        found_Op = true;
        // Create new DAG expression and update inverted expression
        auto pexpr = _insert_expr( ndxVarEl, VarEl, SPVar, pOp, polyindep, polydep, false, orddep );
        *itaux = *pexpr.first;
        spvar  = pexpr.second;
#ifdef MC__SELIM_DEBUG_PROCESS
      std::cout << "Next dependent variable: " << *spvar << std::endl;
#endif
        break;
      }
      if( !found_Op ) throw Exceptions( Exceptions::INVERT );
    }
  }
}

template <typename... ExtOps>
inline std::set< FFVar const*, lt_FFVar >
SElimEnv<ExtOps...>::_dep_expr
( long const ndxVarEl, t_poly const& numer, t_poly const& denom )
{
  // Check denominator is 1
  if( denom.maxord() || denom.coefmon( t_mon() ) != 1e0 )
    throw Exceptions( Exceptions::INVERT );
	
  // Check numerator active dependence in ndxVarEl
  std::set< FFVar const*, lt_FFVar > sdep;
#ifdef MC__SELIM_DEBUG_PROCESS
  std::cout << "DEPENDENTS:";
#endif
  for( auto const& pvar : numer.setvar() ){
    if( pvar->dep().dep( ndxVarEl ).first ){
      sdep.insert( pvar );
#ifdef MC__SELIM_DEBUG_PROCESS
      std::cout << " " << *pvar;
#endif
    }
  }
#ifdef MC__SELIM_DEBUG_PROCESS
  std::cout << std::endl;
#endif
  return sdep;
}

template <typename... ExtOps>
inline std::pair< FFVar const*, SLiftVar const* >
SElimEnv<ExtOps...>::_insert_expr
( long const ndxVarEl, FFVar& var, std::vector<SLiftVar const*> const& SPVar,
  FFOp const* pOp, t_poly const& polyindep, t_poly const& polydep, bool const useprod,
  int const order )
{
  //std::cout << *pOp << std::endl;
  unsigned ndxsp = 0; // variable index for SLiftVar dependence in eliminated variable
  
  // Insert inverse operation in DAG
  switch( pOp->type ){
   case FFOp::IPOW:  _insert_expr( var, polyindep, polydep, useprod, order );
                     var = pow( var, 1./(double)pOp->varin.at(1)->num().n ); break;
   case FFOp::DPOW:  _insert_expr( var, polyindep, polydep, useprod, order );
                     var = pow( var, 1./pOp->varin.at(1)->num().val() );     break;

   case FFOp::SQR:   _insert_expr( var, polyindep, polydep, useprod, order );
                     var = sqrt( var ); break;
   case FFOp::SQRT:  _insert_expr( var, polyindep, polydep, useprod, order );
                     var = sqr( var );  break;
   case FFOp::EXP:   _insert_expr( var, polyindep, polydep, useprod, order );
                     var = log( var );  break;
   case FFOp::LOG:   _insert_expr( var, polyindep, polydep, useprod, order );
                     var = exp( var );  break;
   case FFOp::COS:   _insert_expr( var, polyindep, polydep, useprod, order );
                     var = acos( var ); break;
   case FFOp::SIN:   _insert_expr( var, polyindep, polydep, useprod, order );
                     var = asin( var ); break;
   case FFOp::TAN:   _insert_expr( var, polyindep, polydep, useprod, order );
                     var = atan( var ); break;
   case FFOp::ACOS:  _insert_expr( var, polyindep, polydep, useprod, order );
                     var = cos( var );  break;
   case FFOp::ASIN:  _insert_expr( var, polyindep, polydep, useprod, order );
                     var = sin( var );  break;
   case FFOp::ATAN:  _insert_expr( var, polyindep, polydep, useprod, order );
                     var = tan( var );  break;

   case FFOp::INV:   assert( pOp->varin.at(0)->cst() );
                     _insert_expr( var, polyindep, polydep*pOp->varin.at(0)->num().val(), useprod, -order );
                     ndxsp = 1;
                     break;
		     
   case FFOp::DIV: { assert( SPVar.size() == 2 ); // exactly two operands
                     int ndep1 = _dep_expr( ndxVarEl, SPVar.at(0)->numer(), SPVar.at(0)->denom() ).size();
                     int ndep2 = _dep_expr( ndxVarEl, SPVar.at(1)->numer(), SPVar.at(1)->denom() ).size();
		     assert( ndep1 + ndep2 == 1 ); // exactly one dependence in numerator or denominator
                     // numerator depends on variable ndxVarEl
                     if( ndep1 ){
                       _insert_expr( var, polyindep, polydep, useprod, order );
		       var *= *pOp->varin.at(1);
                       ndxsp = 0;
		     }
                     // denominator depends on variable ndxVarEl
                     else{
                       _insert_expr( var, polyindep, polydep, useprod, -order );
		       var *= *pOp->varin.at(0);
                       ndxsp = 1;
		     }
                     break;
		   }

   case FFOp::CHEB:
   case FFOp::XLOG:
   case FFOp::COSH:
   case FFOp::SINH:
   case FFOp::TANH:
   case FFOp::ERF:
   case FFOp::FABS:
   case FFOp::FSTEP:
   case FFOp::MINF:
   case FFOp::MAXF:
   case FFOp::INTER:
   case FFOp::VAR:
   case FFOp::CNST:
   case FFOp::SHIFT:
   case FFOp::PLUS:
   case FFOp::NEG:
   case FFOp::MINUS:
   case FFOp::SCALE:
   case FFOp::TIMES:
   case FFOp::PROD:
   default:          throw Exceptions( Exceptions::INVERT );
  }

  // Return pointer to internal DAG variable
  return std::make_pair( _ptr_expr( var ), SPVar.at( ndxsp ) );
}

template <typename... ExtOps>
inline FFVar const*
SElimEnv<ExtOps...>::_insert_expr
( FFVar& var, t_poly const& polyindep, t_poly const& polydep, bool const useprod,
  int const order )
{
  //std::cout << "ENTERING _insert_expr\n";
  if( !order ) throw Exceptions( Exceptions::INVERT );
  bool const dagaux = true;

  // Case inverse operation is needed
  if( order < 0 ){
    if( !polyindep.mapmon().empty() ){
      var += insert_dag( polyindep, useprod, dagaux );
    }
    var =  insert_dag( polydep, useprod, dagaux ) / var;
    // Variable has order greater than 1
    if( order < -1 ) var = pow( var, -1./(double)order );
  }
  
  // Case dependent term is constant
  else if( !polydep.maxord() ){
    double const cst = polydep.coefmon( t_mon() );
    // Constant dependent term is one
    if( cst == 1e0 ){
      if( !polyindep.mapmon().empty() )
        var += insert_dag( polyindep, useprod, dagaux );
    }
    // Constant dependent term is negative one
    else if( cst == -1e0 ){
      if( !polyindep.mapmon().empty() )
        var = insert_dag( - polyindep, useprod, dagaux ) - var;
      else
        var = -var;
    }
    else{
      if( !polyindep.mapmon().empty() )
        var += insert_dag( polyindep, useprod, dagaux );
      var /=  cst;
    }
    // Variable has order greater than 1
    if( order > 1 ) var = pow( var, 1./(double)order );
  }

  // Case dependent term is monomial
  else if( polydep.nmon() == 1 && order == 1 ){

    // Split numerator between reduced and non-reduced parts first
    auto const& [mondep,coefdep] = *polydep.mapmon().cbegin();
    t_poly polyindepred, polyindepnred;
    for( auto const& [mon,coef] : polyindep.mapmon() ){
      if( mondep.subseteq( mon ) )
        polyindepred += std::make_pair( mon-mondep, coef/coefdep );
      else
        polyindepnred += std::make_pair( mon, coef );
    }
    if( !polyindepnred.mapmon().empty() )
      var += insert_dag( polyindepnred, useprod, dagaux );
    var /= insert_dag( polydep, useprod, dagaux );
    if( !polyindepred.mapmon().empty() )
      var += insert_dag( polyindepred, useprod, dagaux );
  }
  
  // Case dependent term is linear
  // General case
  else{
    var += insert_dag( polyindep, useprod, dagaux );
    var /= insert_dag( polydep, useprod, dagaux );
    // Variable has order greater than 1
    if( order > 1 ) var = pow( var, 1./(double)order );
  }

  // Return pointer to internal DAG variable
  return _ptr_expr( var );
}

template <typename... ExtOps>
inline FFVar const*
SElimEnv<ExtOps...>::_ptr_expr
( FFVar& var )
{
  // Return pointer to internal DAG variable
  auto itvar = _dag->Vars().find( &var );
  if( itvar != _dag->Vars().end() ){
#ifdef MC__SELIM_DEBUG_INSERT
    std::cout << "INSERTED DAG EXPRESSION:";
    _dag->output( _dag->subgraph( 1, *itvar ) );
#endif
    return *itvar;
  }
  
#ifdef MC__SELIM_CHECK
  assert( var.cst() );
#endif
  return _dag->add_constant( var.num().val() );
}

#if defined(MC__USE_GUROBI)
template <typename... ExtOps>
inline void
SElimEnv<ExtOps...>::_MIP_optimize
()
{
  _MIPexcpt = false;

  try{
    // Run MIP optimization for a minimal representation
    _MIP_reset();
    _MIP_encode();
    //_MIP_initialize();
    _MIP_options();
    _MIP_solve();
    _MIP_decode();
  }
  
  catch(GRBException& e){
    if( options.MIPDISPLEVEL )
      std::cout << "Error code = " << e.getErrorCode() << std::endl
                << e.getMessage() << std::endl;
    _MIPexcpt = true;
  }
}

template <typename... ExtOps>
inline void
SElimEnv<ExtOps...>::_MIP_solve
()
{
  _GRBmodel->update();
  if( options.MIPOUTPUTFILE != "" )
    _GRBmodel->write( options.MIPOUTPUTFILE );
  fedisableexcept(FE_ALL_EXCEPT);
  _GRBmodel->set( GRB_IntAttr_ModelSense, -1 );
  _GRBmodel->optimize();
}

template <typename... ExtOps>
inline void
SElimEnv<ExtOps...>::_MIP_encode
()
{
  // Define Gurobi continuous and binary variables  
  for( auto const& [v,vset] : _mapVarElVar ){
    // This also defines the objective function
    _MIP_var[v]    = _GRBmodel->addVar( 0., 1., _VarWeight[v], GRB_BINARY );
    _MIP_varMTZ[v] = _GRBmodel->addVar( 0., _ndxVar.size()-1., 0., GRB_CONTINUOUS );
    for( auto const& vv : vset )
      _MIP_varvar[v][vv] = _GRBmodel->addVar( 0., 1., 0., GRB_BINARY );
  }

  for( auto const& [e,vset] : _mapCtrVarCand ){
    _MIP_ctr[e] = _GRBmodel->addVar( 0., 1., 0., GRB_BINARY );
#ifdef MC__SELIM_CHECK
    assert( !_mapCtrVarCand[e].empty() );
#endif
    for( auto const& v : vset ){
#ifdef MC__SELIM_DEBUG_MIP
      std::cout << "Define: _MIP_ctrvar[" << e << "][" << v << "]\n";
#endif
      _MIP_ctrvar[e][v] = _GRBmodel->addVar( 0., 1., 0., GRB_BINARY );
    }
  }

  for( auto const& [v,eset] : _mapVarCtrCand ){
    GRBLinExpr sumctrvar;
    for( auto const& e : eset ){
      sumctrvar += _MIP_ctrvar[e][v];
#ifdef MC__SELIM_DEBUG_MIP
      std::cout << "Use: _MIP_ctrvar[" << e << "][" << v << "]\n";
#endif
    }
    _GRBmodel->addConstr( sumctrvar, GRB_EQUAL, _MIP_var[v] );
  }

  for( auto const& [v,vset] : _mapVarElVar ){
    for( auto const& vv : vset )
      _GRBmodel->addConstr( _MIP_varvar[v][vv], GRB_LESS_EQUAL, _MIP_var[v] );
  }

  for( auto const& [e,vset] : _mapCtrVarCand ){
    GRBLinExpr sumctrvar;
    for( auto const& v : vset ){
#ifdef MC__SELIM_DEBUG_MIP
      std::cout << "Use: _MIP_ctrvar[" << e << "][" << v << "]\n";
#endif
      sumctrvar += _MIP_ctrvar[e][v];
    }
    _GRBmodel->addConstr( sumctrvar, GRB_EQUAL, _MIP_ctr[e] );
  }

  for( auto const& [e,vset] : _mapCtrVarCand ){
    for( auto const& v : vset ){
#ifdef MC__SELIM_DEBUG_MIP
      std::cout << "Use: _MIP_ctrvar[" << e << "][" << v << "]\n";
#endif
      for( auto const& vv : _mapVarElVar[v] ){
        if( _mapCtrVar[e].find( vv ) == _mapCtrVar[e].end() )
	  _GRBmodel->addConstr( 1. - _MIP_varvar[v][vv], GRB_GREATER_EQUAL, _MIP_ctrvar[e][v] );
        else
	  _GRBmodel->addConstr( _MIP_varvar[v][vv], GRB_GREATER_EQUAL, _MIP_ctrvar[e][v] );
      }
    }
  }
  
  // MTZ subtour elimination
  for( auto const& [v,vset] : _mapVarElVar ){
    for( auto const& vv : vset )
      _GRBmodel->addConstr( _MIP_varMTZ[v] - _MIP_varMTZ[vv], GRB_LESS_EQUAL,
                            _ndxVar.size() * (1.-_MIP_varvar[v][vv]) - 1. );
  }
}

template <typename... ExtOps>
inline void
SElimEnv<ExtOps...>::_MIP_decode
()
{
  // Get order of variable elimination
  std::multimap<double,unsigned> VarMTZ;
  for( auto const& [v,grbvar] : _MIP_var ){
#ifdef MC__SELIM_DEBUG_MIP
    std::cout << "Variable " << v << " (" << _Var[v] << "): "
              << _MIP_var[v].get(GRB_DoubleAttr_X) << std::endl;
#endif
    if( grbvar.get(GRB_DoubleAttr_X) < 0.9 ) continue;
    VarMTZ.insert( std::make_pair( _MIP_varMTZ[v].get(GRB_DoubleAttr_X), v ) );
#ifdef MC__SELIM_DEBUG_MIP
    std::cout << "Variable " << v << " (" << _Var[v] << ") eliminated - MTZ index: "
              << _MIP_varMTZ[v].get(GRB_DoubleAttr_X) << std::endl;
#endif
  }
  std::map<unsigned,unsigned> OrdElim;
  unsigned count=VarMTZ.size();
  for( auto const& [u,v] : VarMTZ )
    OrdElim[v] = --count;
  
  // Add entries in _VarElim for eliminated variables and corresponding constraints
  auto& [vvar,vctr,vaux] = _VarElim;
  vvar.resize( OrdElim.size() );
  vctr.resize( OrdElim.size() );
  //_VarElim.resize( OrdElim.size() );
  for( auto const& [e,vmap] : _MIP_ctrvar ){
    if( _MIP_ctr[e].get(GRB_DoubleAttr_X) < 0.9 ) continue;
#ifdef MC__SELIM_DEBUG_MIP
    std::cout << "Constraint " << e << " (" << _Ctr[e] << ") used for elimination\n";
#endif
    for( auto const& [v,grbvar] : vmap ){
      if( grbvar.get(GRB_DoubleAttr_X) < 0.9 ) continue;
      auto itdagvar = _dag->Vars().find( const_cast<FFVar*>(&_Var[v]) );
#ifdef MC__SELIM_CHECK
      assert( itdagvar != _dag->Vars().end() );
#endif
      auto itdagctr = _dag->Vars().find( const_cast<FFVar*>(&_Ctr[e]) );
#ifdef MC__SELIM_CHECK
      assert( itdagctr != _dag->Vars().end() );
#endif
      vvar[OrdElim[v]] = **itdagvar;
      vctr[OrdElim[v]] = **itdagctr;
      //_VarElim[OrdElim[v]] = std::make_tuple( *itdagvar, *itdagctr, nullptr );
    }
  }
}

template <typename... ExtOps>
inline void
SElimEnv<ExtOps...>::_MIP_options
()
{
  _GRBmodel->getEnv().set( GRB_IntParam_OutputFlag,        options.MIPDISPLEVEL );
  _GRBmodel->getEnv().set( GRB_IntParam_Method,            options.LPALGO );
  _GRBmodel->getEnv().set( GRB_IntParam_Presolve,          options.LPPRESOLVE );
  _GRBmodel->getEnv().set( GRB_IntParam_Threads,           options.MIPTHREADS );
  _GRBmodel->getEnv().set( GRB_IntParam_ConcurrentMIP,     options.MIPCONCURRENT );
  _GRBmodel->getEnv().set( GRB_IntParam_MIPFocus,          options.MIPFOCUS );
  _GRBmodel->getEnv().set( GRB_IntParam_NumericFocus,      options.MIPNUMFOCUS );
  _GRBmodel->getEnv().set( GRB_DoubleParam_Heuristics,     options.MIPHEURISTICS );
  _GRBmodel->getEnv().set( GRB_DoubleParam_FeasibilityTol, options.LPFEASTOL );
  _GRBmodel->getEnv().set( GRB_DoubleParam_OptimalityTol,  options.LPOPTIMTOL );
  _GRBmodel->getEnv().set( GRB_DoubleParam_MIPGap,         options.MIPRELGAP );
  _GRBmodel->getEnv().set( GRB_DoubleParam_MIPGapAbs,      options.MIPABSGAP );
  _GRBmodel->getEnv().set( GRB_DoubleParam_TimeLimit,      options.MIPTIMELIMIT );
}

template <typename... ExtOps>
inline void
SElimEnv<ExtOps...>::_MIP_reset
()
{
  _MIP_var.clear();
  _MIP_varMTZ.clear();
  _MIP_ctr.clear();
  _MIP_ctrvar.clear();
  _MIP_varvar.clear();
  delete _GRBmodel;
  _GRBmodel = new GRBModel( *_GRBenv );
}
#endif // #if defined(MC__USE_GUROBI)

} // namespace mc

#endif
