// Copyright (C) Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the Eclipse Public License.

/*!
\page page_SELIM Variable Elimination in Factorable Expressions
\author Benoit Chachuat & Dominik Bongartz
\date 2023
\bug No known bugs.

The classes mc::SElimEnv and mc::SElimVar defined in <tt>selim.hpp</tt> enable the elimination of a subset of variables from factorable expressions through expoiting solvable (invertible) equality constraints. 


\section sec_SELIM_algo What is the algorithm used for the elimination?

Given a set of equality constraints \f${\bf f}({\bf x}) = {\bf 0}\f$ with \f${\bf f}:\mathbb{R}^n\to\mathbb{R}^m\f$ and \f$n\geq m\f$, we seek to partition the variable set into \f${\bf x} =: [{\bf x}_{\rm indep},{\bf x}_{\rm dep}]\f$, where \f${\bf x}_{\rm indep}\in\mathbb{R}^{n_{\rm indep}}\f$ and \f${\bf x}_{\rm dep}\in\mathbb{R}^{n_{\rm dep}}\f$ are the independent and dependent decisions, respectively; and to construct an explicit mapping \f${\bf g}:\mathbb{R}^n\to\mathbb{R}^m\f$ such that \f${\bf x}_{\rm dep} = {\bf g}({\bf x}_{\rm indep})\f$.

This elimination proceeds in 3 steps:
-# Identify which equality constraints can be inverted analytically in terms of which variables.

-# Determine which variables to eliminate using which equality constraints, e.g. to maximize the number of eliminated variables, based on: (i) a mixed-integer programming (MIP) formulation; or (ii) a greedy heuristic inspired by the method of <A href="https://doi.org/10.1016/0098-1354(79)80057-5">Hernandez & Sargent (1979)</A> to reorder the equations then select the variables through reordering the incidence matrix to bordered lower-triangular form.

-# Construct expressions for the eliminated variables, so they can be eliminated through composition.
.

The MIP formulation solves the following problem using Gurobi:
\f{align*}
\displaystyle\min_{\boldsymbol{v},\boldsymbol{c},\boldsymbol{d},\boldsymbol{e},\boldsymbol{z}}\ & \sum_{i=1}^{\bar{V}} \omega_i v_i\\
\displaystyle\text{s.t.}\ \ \ 
& \sum_{k\in E_i} e_{i,k} = v_i,\ \ i=1\ldots\bar{V}\\
& \sum_{i\in C_k} e_{i,k} = c_k,\ \ k=1\ldots\bar{C}\\
& d_{i,j} \leq v_i,\ \ i=1\ldots\bar{V},\ j\in D_i\\
& d_{i,j} \geq e_{i,k},\ \ k=1\ldots\bar{C},\ i\in C_k,\ j\in D_i \cap V_k\\
& d_{i,j} + e_{i,k} \leq 1,\ \ k=1\ldots\bar{C},\ i\in C_k,\ j\notin D_i \cap V_k\\
& z_i - z_j \leq (1-d_{i,j})\bar{V}-1,\ \ i=1\ldots\bar{V},\ j\in D_i\\
& v_i, d_{i,j}, c_k, e_{i,k}\in\{0,1\},\ \ z_i\in [0,\bar{V}-1]
\f}
with the following sets, parameters and variables:
- \f$E_i\subset\{1,\ldots,\bar{C}\}\f$ subset of candidate constraints for variable \f$i\f$ elimination
- \f$D_i\subset\{1,\ldots,\bar{V}\}\f$ subset of candidate variables for variable \f$i\f$ elimination
- \f$V_k\subset\{1,\ldots,\bar{V}\}\f$ subset of variables participating in constraint \f$k\f$
- \f$C_k\subset\{1,\ldots,\bar{V}\}\f$ subset of candidate variables for elimination participating in constraint \f$k\f$
- \f$\omega_i\in\{0,1\}\f$ weight coefficient for variable \f$i=1\ldots \bar{V}\f$
- \f$v_i\in\{0,1\}\f$ encoding whether variable \f$i=1\ldots \bar{V}\f$ is selected for elimination
- \f$c_k\in\{0,1\}\f$ encoding whether constraint \f$k=1\ldots \bar{C}\f$ is selected for elimination
- \f$d_{i,j}\in\{0,1\}\f$ encoding whether variable \f$j\in D_i\f$ is used to eliminate variable \f$i=1\ldots \bar{V}\f$
- \f$e_{i,k}\in\{0,1\}\f$ encoding whether constraint \f$k\in E_i\f$ is selected to eliminate variable \f$i=1\ldots \bar{V}\f$
- \f$z_i\in[0,\bar{V}-1]\f$ encoding precedence in MTZ subtour elimination
.


\section sec_SElim_process How do I eliminate variables using factorable equality constraints?

For illustration, consider the following pair of nonlinear equality constraints in three variables \f$x_0,x_1,x_2\f$:
\f{align*}
  0 &= \frac{3x_0(x_2)^2}{x_1} - 2x_0x_1 - x_0 - 1\\
  0 &= \frac{2}{x_1} + \frac{3}{x_2} - 1
\f}

The elimination requires the header file <tt>selim.hpp</tt> to be included:

\code
      #include "selim.hpp"
\endcode

A DAG of the factorable function defining the equality constraints is first created:

\code
      mc::FFGraph DAG;
      const unsigned NX = 3, NF = 2;
      mc::FFVar X[NX];
      for( unsigned i(0); i<NX; i++ ) X[i].set( &DAG );
      mc::FFVar F[NF];
      F[0] = ( 3. * X[0] * sqr( X[2] ) ) / X[1] - 2. * X[0] * X[1] - X[0] - 1;
      F[1] = 2./X[1] + 3./X[2] - 1.;
      std::cout << DAG;
\endcode

The last line displays the following information about the DAG:

\verbatim
    DAG VARIABLES:
      V0	 => { Z1 Z5 Z9 }
      V1	 => { Z2 Z7 Z13 }
      V2	 => { Z3 Z12 }

    DAG INTERMEDIATES:
      Z1	<=  V0 x Z0		 => { Z2 }
      Z2	<=  V1 x Z1		 => { Z8 }
      Z3	<=  SQR( V2 )		 => { Z6 }
      Z5	<=  V0 x Z4		 => { Z6 }
      Z6	<=  Z3 x Z5		 => { Z7 }
      Z7	<=  Z6 / V1		 => { Z8 }
      Z8	<=  Z7 - Z2		 => { Z9 }
      Z9	<=  Z8 - V0		 => { Z11 }
      Z11	<=  Z9 + Z10		 => { }
      Z12	<=  Z4 / V2		 => { Z14 }
      Z13	<=  Z0 / V1		 => { Z14 }
      Z14	<=  Z12 + Z13		 => { Z15 }
      Z15	<=  Z14 + Z10		 => { }
      Z10	<=  -1(I)		 => { Z11 Z15 }
      Z0	<=  2(I)		 => { Z1 Z13 }
      Z4	<=  3(I)		 => { Z5 Z12 }
\endverbatim

Next, an environment <a>mc::SElimEnv</a> is defined and the method <a>mc::SElimEnv::process</a> is invoked to perform the elimination:

\code
      mc::SElimEnv SPE( &DAG );
      SPE.process( NF, F );
      std::cout << SPE;
\endcode

The following information is displayed in this instance:

\verbatim
    2 VARIABLES MAY BE ELIMINATED

    OPERATIONS IN SUBGRAPH OF V2 USING Z15=0:
      Z4	<<  3(I)	
      Z0	<<  2(I)	
      V1	<<  VARIABLE
      Z13	<<  Z0 / V1	
      Z16	<<  - Z13	
      Z17	<<  1(I)	
      Z18	<<  Z16 + Z17	
      Z19	<<  Z4 / Z18	
    DEPENDENTS IN SUBGRAPH OF V2 USING Z15=0:
      0:  Z19
    WORK ARRAY SIZE: 8

    OPERATIONS IN SUBGRAPH OF V0 USING Z11=0:
      Z17	<<  1(I)	
      V2	<<  VARIABLE
      Z3	<<  SQR( V2 )	
      Z4	<<  3(I)	
      Z20	<<  Z3 x Z4	
      V1	<-  VARIABLE
      Z21	<<  Z20 / V1	
      Z10	<<  -1(I)	
      Z22	<<  Z21 + Z10	
      Z23	<<  -2(I)	
      Z24	<<  V1 x Z23	
      Z25	<<  Z22 + Z24	
      Z26	<<  Z17 / Z25	
    DEPENDENTS IN SUBGRAPH OF V0 USING Z11=0:
      0:  Z26
    WORK ARRAY SIZE: 13
\endverbatim

These results show that both variables \f$x_0\f$ and \f$x_2\f$ can be eliminated through inverting both equality constraints.
- The expression of variable \f$x_2\f$ in terms of \f$x_1\f$ is given by:
\f{align*}
  x_2 &= \frac{3}{\displaystyle 1-\frac{2}{x_1}}
\f}
- The expression of variable \f$x_0\f$ in terms of \f$x_1\f$ and \f$x_2\f$---where \f$x_2\f$ could be eliminated recursively in terms of \f$x_1\f$ using the previous expression---is given by:
\f{align*}
  x_0 &= \frac{1}{\displaystyle 3\frac{(x_2)^2}{x_1} - 2x_1 - 1}
\f}
.

The method mc::SElimEnv::VarElim allows retreiving a 3-tuple of vectors of (i) the eliminated variables, (ii) the original constraint used for the elimination, and (iii) the inverted constraint expression, with entries in a feasible order of elimination (lower-triangular form). 


\section sec_SELIM_refs References

- Hernandez, R., Sargent, R. W. H., <A href="https://doi.org/10.1016/0098-1354(79)80057-5">A new algorithm for process flowsheeting</A>, <I>Computers & Chemical Engineering</I>, <b>3</b>(1-4):363-371, 2015
*/

// TO DO:
// - Complete documentation

#ifndef MC__SELIM_H
#define MC__SELIM_H

#include "ffunc.hpp"
#include "ffdep.hpp"
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

//! @brief Base class for variable elimination in factorable equality constraints
////////////////////////////////////////////////////////////////////////
//! mc::SElimBase is a C++ base class defining the environment for
//! variable elimination from factorable expressions through expoiting
//! solvable equality constraints
////////////////////////////////////////////////////////////////////////
class SElimBase
////////////////////////////////////////////////////////////////////////
{
public:

  //! @brief Default Constructor
  SElimBase
    ()
    {
#if defined(MC__USE_GUROBI)
      _GRBenv   = new GRBEnv();
      _GRBmodel = nullptr;
#endif
    }
    
  //! @brief Destructor
  virtual ~SElimBase
    ()
    {
#if defined(MC__USE_GUROBI)
      delete _GRBmodel;
      delete _GRBenv;
#endif
    }

  //! @brief Options of mc::SElimBase
  struct Options
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
      DISPFULL(false)
      {}
    //! @brief Assignment of mc::SElimBase::Options
    Options& operator=
      ( Options const& opt ){
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
        DISPFULL        = opt.DISPFULL;
#endif
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
    //! @brief Whether to display a full subgraph of the lifted expressions (default: false)
    bool DISPFULL;
  } options;


protected:

  //! @brief Candidate variables
  std::set<unsigned> _ndxVar;

  //! @brief Independent DAG variables with associated weights/priorities
  std::map<unsigned,double> _VarWeight;

  //! @brief Candidate constraints from which to eliminate a given variable
  std::map<unsigned,std::set<unsigned>> _mapVarCtrCand;

  //! @brief Candidate variables to eliminate a given variable
  std::map<unsigned,std::set<unsigned>> _mapVarElVar;

  //! @brief Candidate constraints
  std::set<unsigned> _ndxCtr;

  //! @brief Participating variables in a given constraint
  std::map<unsigned,std::set<unsigned>> _mapCtrVar;

  //! @brief Candidate variables to eliminate from a given constraint
  std::map<unsigned,std::set<unsigned>> _mapCtrVarCand;

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
    ( Options const& opt );

  //! @brief Solve MIP optimization model
  void _MIP_solve
    ();

  //! @brief Encode MIP optimization model for maximal elimination
  void _MIP_encode
    ();

  //! @brief Reset variable vectors in MIP optimization model
  void _MIP_reset
    ();

  //! @brief Set options in MIP optimization model
  void _MIP_options
    ();
#endif
};

#if defined(MC__USE_GUROBI)
inline int const SElimBase::Options::LPALGO_DEFAULT;
//inline SElimBase::Options SElimBase::options;

inline
void
SElimBase::_MIP_optimize
( Options const& opt )
{
  _MIPexcpt = false;
  options   = opt;

  try{
    // Run MIP optimization for a minimal representation
    _MIP_reset();
    _MIP_encode();
    _MIP_options();
    _MIP_solve();
  }
  
  catch(GRBException& e){
    if( options.MIPDISPLEVEL )
      std::cout << "Error code = " << e.getErrorCode() << std::endl
                << e.getMessage() << std::endl;
    _MIPexcpt = true;
  }
}

inline
void
SElimBase::_MIP_solve
()
{
  _GRBmodel->update();
  if( options.MIPOUTPUTFILE != "" )
    _GRBmodel->write( options.MIPOUTPUTFILE );
  fedisableexcept(FE_ALL_EXCEPT);
  _GRBmodel->set( GRB_IntAttr_ModelSense, -1 );
  _GRBmodel->optimize();
}

inline
void
SElimBase::_MIP_encode
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

inline
void
SElimBase::_MIP_options
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

inline
void
SElimBase::_MIP_reset
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

//! @brief Environment for variable elimination in factorable equality constraints
////////////////////////////////////////////////////////////////////////
//! mc::SElimEnv is a C++ class defining the environment for variable
//! elimination from factorable expressions through expoiting solvable
//! equality constraints
////////////////////////////////////////////////////////////////////////
class SElimEnv
: public SElimBase,
  protected virtual SLiftEnv
////////////////////////////////////////////////////////////////////////
{
  using SLiftEnv::_dag;
  using SLiftEnv::_OpLift;
  using SLiftEnv::_SPDep;

  using SLiftEnv::dag;
  using SLiftEnv::insert_dag;

  friend std::ostream& operator<<
    ( std::ostream&, const SElimEnv& );

public:

  typedef std::tuple< std::vector<FFVar>, std::vector<FFVar>, std::vector<FFVar> > t_VarElim;
  typedef SLiftEnv t_lift;
  typedef typename SLiftVar::t_poly t_poly;
  typedef typename t_poly::t_mon t_mon;

  //! @brief Default Constructor
  SElimEnv
    ( FFGraph* dag=nullptr )
    : SLiftEnv( dag ),
      SElimBase()
    {}
    
  //! @brief Destructor
  virtual ~SElimEnv
    ()
    { _reset(); }

  //! @brief Retreive vector of partipating variables in equality constraints
  std::vector<FFVar> const& Var
    ()
    const
    { return _Var; }


  //! @brief Retreive tuple of vectors <ELIMINATED VARIABLE, INVERTED CONSTRAINT, INVERTED EXPRESSION> with entries in a feasible order of elimination
  t_VarElim const& VarElim
    ()
    const
    { return _VarElim; }

  //! @brief Set DAG environment
  void set
    (  FFGraph* dag )
    {
      SLiftEnv::set( dag );
      _reset();
    }

  //! @brief Reset intermediate expressions
  void reset
    ()
    {
      SLiftEnv::_reset();
      _reset();
    }

  //! @brief Process the <a>ndxCtr</a> equality constraints in array <a>pCtr</a> 
  void process
    ( std::set<unsigned> const& ndxCtr, FFVar const* pCtr, 
      std::map<FFVar const*,double,lt_FFVar> const& wVar=std::map<FFVar const*,double,lt_FFVar>(),
      bool const add2dag=true );

  //! @brief Process the <a>nCtr</a> equality constraints in array <a>pCtr</a>
  void process
    ( unsigned const nCtr, FFVar const* pCtr, 
      std::map<FFVar const*,double,lt_FFVar> const& wVar=std::map<FFVar const*,double,lt_FFVar>(),
      const bool add2dag=true );

  //! @brief Process the <a>ndxCtr</a> equality constraints in vector <a>vCtr</a> 
  void process
    ( std::set<unsigned> const& ndxCtr, std::vector<FFVar> const& vCtr, 
      std::map<FFVar const*,double,lt_FFVar> const& wVar=std::map<FFVar const*,double,lt_FFVar>(),
      bool const add2dag=true )
    { process( ndxCtr, vCtr.data(), wVar, add2dag ); }

  //! @brief Process the <a>nCtr</a> equality constraints in vector <a>vCtr</a>
  void process
    ( std::vector<FFVar> const& vCtr, 
      std::map<FFVar const*,double,lt_FFVar> const& wVar=std::map<FFVar const*,double,lt_FFVar>(),
      const bool add2dag=true )
    { process( vCtr.size(), vCtr.data(), wVar, add2dag ); }

  //! @brief Exceptions of mc::SElimEnv
  class Exceptions
  {
   public:
    //! @brief Enumeration type for SElimVar exception handling
    enum TYPE{
      MIPERR=0,       //!< Call to MIP solver disabled
      INVERT,         //!< Elimination process failed
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
  struct Options
  : public SElimBase::Options
  {
    //! @brief Constructor
    Options
    ()
    : SElimBase::Options(),
      MULTMAX(0), ELIMLIN(true), ELIMMLIN(true),
      ELIMNLIN( { FFInv::Options::INV, FFInv::Options::SQRT, FFInv::Options::EXP,
                  FFInv::Options::LOG, FFInv::Options::RPOW } )
      {
        SLIFT.KEEPFACT = true;
        SLIFT.LIFTDIV  = true;
        SLIFT.LIFTIPOW = false;
      }
    //! @brief Assignment of mc::SElimEnv::Options
    Options& operator=
      ( Options const& opt ){
        SElimBase::Options::operator=( opt );
        MULTMAX  = opt.MULTMAX;
        ELIMLIN  = opt.ELIMLIN;
        ELIMMLIN = opt.ELIMMLIN;
        ELIMNLIN = opt.ELIMNLIN;
        SLIFT    = opt.SLIFT;
        return *this;
      }
    //! @brief Maximal multiplicity of eliminated variables
    unsigned MULTMAX;
    //! @brief Whether to invert linear operations
    bool ELIMLIN;
    //! @brief Whether to invert multilinear operations
    bool ELIMMLIN;
    //! @brief Set of invertible nonlinear operations
    std::set<FFInv::Options::NLINV> ELIMNLIN;
    //! @brief Options for factorable function decomposition using mc::SLiftEnv
    typename SLiftEnv::Options SLIFT;
  } options;


protected:

  //! @brief Independent DAG variables participating in equality constraints
  std::vector<FFVar> _Var;

  //! @brief Equality constraints
  std::vector<FFVar> _Ctr;

  //! @brief Dependency information for DAG variables
  std::vector<FFDep> _DVar;

  //! @brief Invertibility information for DAG variables
  std::vector<FFInv> _IVar;

  //! @brief Invertibility information for equality constraints
  std::vector<FFInv> _ICtr;
  
  //! @brief Variable multiplicity
  std::map<unsigned,unsigned> _mapVarMult;

  //! @brief Tuple of eliminated DAG variables and corresponding equality constraint and auxiliary variable
  t_VarElim _VarElim;

#if defined(MC__USE_GUROBI)
  //! @brief Decode MIP optimal solution for maximal elimination
  void _MIP_decode
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
    ( long const ndxVarEl, t_poly const& numer, t_poly const& denom, bool toplevel=false );

  //! @brief Return pointer to intrenal DAG variable of expression
  FFVar const* _ptr_expr
    ( FFVar& var );

  //! @brief Erase private storage data
  void _reset
    ();
};

inline
std::ostream&
operator<<
( std::ostream& out, SElimEnv const& env)
{
  auto const& [vvar,vctr,vaux] = env._VarElim;
  std::cout << std::endl
            << vvar.size() << " VARIABLES CAN BE ELIMINATED" << std::endl;

  if( vaux.empty() ){
    for( auto itvar = vvar.begin(), itctr = vctr.begin();
         itvar != vvar.end(); ++itvar, ++itctr )
      out << "VARIABLE " << *itvar << " USING CONSTRAINT " << *itctr << "=0" << std::endl;
    return out;
  }

  else if( !env.options.DISPFULL ){
    FFExpr::options.LANG = FFExpr::Options::DAG;
    auto sgExpr = env._dag->subgraph( vaux.size(), vaux.data() );
    auto vExpr  = FFExpr::subgraph( env._dag, sgExpr );
    unsigned iaux = 0;
    for( auto itvar = vvar.begin(), itctr = vctr.begin();
         itvar != vvar.end(); ++itvar, ++itctr, ++iaux )
      out << "  0 = " << std::left << std::setw(6) <<  *itctr << "-> " << *itvar << " = " << vExpr[iaux] << std::endl;
  }

  else{
    for( auto itvar = vvar.begin(), itctr = vctr.begin(), itaux = vaux.begin();
         itvar != vvar.end(); ++itvar, ++itctr, ++itaux ){
      std::ostringstream ext; 
      ext << " OF " << *itvar << " USING " << *itctr << "=0";
      env._dag->output( env._dag->subgraph( 1, &*itaux ), ext.str(), out );
    }
  }
  
  return out;
}

inline
void
SElimEnv::_reset
()
{
  _Var.clear();
  _DVar.clear();
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

inline
void
SElimEnv::process
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

inline
void
SElimEnv::process
( unsigned const nCtr, FFVar const* pCtr,
  std::map<FFVar const*,double,lt_FFVar> const& wVar,
  bool const add2dag )
{
  // Reset internal variables
  _reset();

  // Pass options to base class
  SElimBase::options = options;

  // Update participating variables in _Var
  auto sgCtr = _dag->subgraph( nCtr, pCtr );
  unsigned v=0;
  for( auto const& Op : sgCtr.l_op ){
    if( Op->type != FFOp::VAR ) continue;
    _Var.push_back( *Op->varout[0] );
    _DVar.push_back( FFDep().indep( Op->varout[0]->id().second ) );
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
  _MIP_optimize( options );
  if( _MIPexcpt ) throw Exceptions( Exceptions::MIPERR );
  _MIP_decode();
#else
  throw Exceptions( Exceptions::MIPERR );
#endif

  // Decode optimization results
  auto& [vvar,vctr,vaux] = _VarElim;
#ifdef MC__SELIM_DEBUG_PROCESS
  std::cout << std::endl << vvar.size() << " eliminated variables and corresponding constraints:" << std::endl;
  for( auto itvar = vvar.cbegin(), itctr = vctr.cbegin();
       itvar != vvar.cend(); ++itvar, ++itctr ){
    std::cout << *itvar << " <- " << *itctr;
    auto sgAux = _dag->subgraph( 1, &*itctr );
    auto vAux  = FFExpr::subgraph( _dag, sgAux );
    std::cout << " = " << vAux[0] << std::endl;
  }
  { std::cout << "PAUSED, ENTER <1> TO CONTINUE "; int dum; std::cin >> dum; }
#endif
  
  // No transcription in DAG if <a>add2dag</a> is false
  if( !add2dag ) return;

  t_poly::options.BASIS = t_poly::Options::MONOM;
  t_lift::options = options.SLIFT;

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
    bool toplevel = true;
    for( FFVar VarEl = 0.; ; ){

      // Check active dependence of spvar
      auto const& denom = spvar->denom();
      auto const& numer = spvar->numer();
      long ndxVarEl = itvarel->id().second;
      auto&& sdep = _dep_expr( ndxVarEl, numer, denom, toplevel );
      toplevel = false;
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
              for( auto const& [pOp,SPVar] : _OpLift ){
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
        for( unsigned q=0; (mondepord+mondep).subseteq( mon ); ++q )
          mondepord += mondep;
        // Make sure order is consistent throughout expression
        if( orddep && mondepord.tord && orddep != mondepord.tord ) throw Exceptions( Exceptions::INVERT );
        if( orddep < mondepord.tord ) orddep = mondepord.tord; 
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
      for( auto const& [pOp,SPVar] : _OpLift ){
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

inline
std::set< FFVar const*, lt_FFVar >
SElimEnv::_dep_expr
( long const ndxVarEl, t_poly const& numer, t_poly const& denom, bool toplevel )
{
  // Check denominator is 1
  if( !toplevel && ( denom.maxord() || denom.coefmon( t_mon() ) != 1e0 ) )
    //std::cerr << "SElimEnv::_dep_expr ** warning: denominator is not 1" << std::endl;
    throw Exceptions( Exceptions::INVERT );

  // Check numerator active dependence in ndxVarEl
  std::set< FFVar const*, lt_FFVar > sdep;
#ifdef MC__SELIM_DEBUG_PROCESS
  std::cout << "DEPENDENTS (" << ndxVarEl << "):";
#endif
  FFDep dvar;
  for( auto const& pvar : numer.setvar() ){
    _dag->eval( 1, pvar, &dvar, _Var.size(), _Var.data(), _DVar.data() );
    if( dvar.dep( ndxVarEl ).first ){
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

inline
std::pair< FFVar const*, SLiftVar const* >
SElimEnv::_insert_expr
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

inline
FFVar const*
SElimEnv::_insert_expr
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

inline
FFVar const*
SElimEnv::_ptr_expr
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
inline
void
SElimEnv::_MIP_decode
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
#endif // #if defined(MC__USE_GUROBI)

} // namespace mc

#endif
