// Copyright (C) Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the Eclipse Public License.

/*!
\page page_SELIM Variable Elimination in Factorable Expressions
\author Benoit Chachuat & Dominik Bongartz
\date 2023
\bug No known bugs.

The classes mc::SElimEnv and mc::SElimVar defined in <tt>selim.hpp</tt> enable the elimination of (a subset of) variables from factorable expressions through expoiting solvable (invertible) equality constraints. 

Given a set of equality constraints \f${\bf f}({\bf x}) = {\bf 0}\f$, we week a partition of the variable set \f${\bf x} =: [{\bf x},  = {\bf 0}\f$
The elimination proceeds in 3 steps:
- 1. Identify which equality constraints can be solved analytically for which variables.

2. Determine the order to process the equality constraints to minimize the number of remaining variables based on the method of Hernandez & Sargent (1979).

3. Select which variables to eliminate through reordering the incidence matrix to bordered lower-triangular form.

4. Rearrange the optimization problem accordingly, and introduce new (possibly nonlinear) inequality constraints to account for the variable bounds of the eliminated variables 𝑥𝑒𝑙. 
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
#define MC__SELIM_DEBUG_PROCESS
#define MC__SELIM_DEBUG_INSERT

namespace mc
{

//template <typename DAG> class SElimVar;

//! @brief Environment for variable elimination in factorable equality constraints
////////////////////////////////////////////////////////////////////////
//! mc::SElimEnv is a C++ class defining the environment for variable
//! elimination from factorable expressions through expoiting solvable
//! equality constraints
////////////////////////////////////////////////////////////////////////
template < typename DAG >
class SElimEnv:
  protected virtual SLiftEnv<DAG>
////////////////////////////////////////////////////////////////////////
{
/*
  friend class SElimVar<DAG>;
  template <typename D> friend  std::ostream& operator<< ( std::ostream&, SElimEnv<D> const& );
  template <typename D> friend  SElimVar<D> inv ( SElimVar<D> const& );
  template <typename D> friend  SElimVar<D> exp ( SElimVar<D> const& );
  template <typename D> friend  SElimVar<D> log ( SElimVar<D> const& );
  template <typename D> friend  SElimVar<D> xlog( SElimVar<D> const& );
  template <typename D> friend  SElimVar<D> sqrt( SElimVar<D> const& );
  template <typename D> friend  SElimVar<D> sqr ( SElimVar<D> const& );
  template <typename D> friend  SElimVar<D> pow ( SElimVar<D> const&, int const );  
  template <typename D> friend  SElimVar<D> pow ( SElimVar<D> const&, double const& );  
  template <typename D> friend  SElimVar<D> cheb( SElimVar<D> const&, const unsigned );  
  template <typename D> friend  SElimVar<D> prod( const unsigned, SElimVar<D> const* );  
  template <typename D> friend  SElimVar<D> cos ( SElimVar<D> const& );
  template <typename D> friend  SElimVar<D> sin ( SElimVar<D> const& );
  template <typename D> friend  SElimVar<D> tan ( SElimVar<D> const& );
  template <typename D> friend  SElimVar<D> acos( SElimVar<D> const& );
  template <typename D> friend  SElimVar<D> asin( SElimVar<D> const& );
  template <typename D> friend  SElimVar<D> atan( SElimVar<D> const& );
  template <typename D> friend  SElimVar<D> cosh( SElimVar<D> const& );
  template <typename D> friend  SElimVar<D> sinh( SElimVar<D> const& );
  template <typename D> friend  SElimVar<D> tanh( SElimVar<D> const& );
  template <typename D> friend  SElimVar<D> fabs( SElimVar<D> const& );
  template <typename D> friend  SElimVar<D> erf( SElimVar<D> const& );
  template <typename D> friend  tSElimVar<D> fstep( SElimVar<D> const& );
  template <typename D> friend  SElimVar<D> max( SElimVar<D> const&, SElimVar<D> const& );  
  template <typename D> friend  SElimVar<D> min( SElimVar<D> const&, SElimVar<D> const& );  
  template <typename D> friend  SElimVar<D> lmtd( SElimVar<D> const&, SElimVar<D> const& );  
  template <typename D> friend  SElimVar<D> rlmtd( SElimVar<D> const&, SElimVar<D> const& );  
*/

  using SLiftEnv<DAG>::_dag;
  using SLiftEnv<DAG>::_Interm;
  using SLiftEnv<DAG>::_SPDep;
  using SLiftEnv<DAG>::insert_dag;

public:

  typedef std::map< FFVar const*, FFVar const*, lt_FFVar > t_VarElim;
  typedef SLiftEnv<DAG> t_lift;
  typedef typename SLiftVar<DAG>::t_poly t_poly;
  typedef typename t_poly::t_mon t_mon;

  //! @brief Default Constructor
  SElimEnv
    ( DAG* dag=nullptr )
    : SLiftEnv<DAG>(dag)
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
    
  // Retreive pointer to DAG
  DAG* dag
    ()
    const
    { return _dag; };

  //! @brief Retreive map between eliminated DAG variables and coresponding equality constraint
  t_VarElim& VarElim
    ()
    { return _VarElim; }

  //! @brief Set DAG environment
  void set
    (  DAG* dag )
    { SLiftEnv<DAG>::set( dag );
      _reset(); }

  //! @brief Reset intermediate expressions
  void reset
    ()
    { SLiftEnv<DAG>::_reset();
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
      DAGERR=0,       //!< Operation involving a factorable expression linked to a different DAG
      ENVERR,         //!< Operation between factorable expressions linked to different environments or without an environment
      MIPERR,         //!< Call to MIP solver disabled
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
      case DAGERR:
        return "mc::SElimEnv\t Operation involving a factorable expression linked to a different DAG is not allowed";
      case ENVERR:
        return "mc::SElimEnv\t Operation between factorable expressions linked to different environments or without an environment is not allowed";
      case MIPERR:
        return "mc::SElimEnv\t Mixed-integer programming solver is disabled";
      case INVERT:
        return "mc::SElimEnv\t Internal error during elimination process";
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
#if defined(MC__USE_GUROBI)
      LPALGO( LPALGO_DEFAULT ), LPPRESOLVE(-1),
      LPFEASTOL(1e-7), LPOPTIMTOL(1e-7),
      MIPRELGAP(1e-3), MIPABSGAP(1e-3), MIPTHREADS(0),
      MIPCONCURRENT(1), MIPFOCUS(0), MIPHEURISTICS(0.2),
      MIPNUMFOCUS(0), MIPDISPLEVEL(1), MIPOUTPUTFILE(""),
      MIPTIMELIMIT(600),
#endif
      ELIM_LIN(true), ELIM_MLIN(true),
      ELIM_NLIN( {FFInv::Options::INV,FFInv::Options::SQRT,FFInv::Options::EXP,
                  FFInv::Options::LOG,FFInv::Options::RPOW} )
      {}
    //! @brief Assignment of mc::SElimEnv<DAG>::Options
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
        ELIM_LIN  = opt.ELIM_LIN;
        ELIM_MLIN = opt.ELIM_MLIN;
        ELIM_NLIN = opt.ELIM_NLIN;
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
    //! @brief Whether to invert linear operations
    bool ELIM_LIN;
    //! @brief Whether to invert multilinear operations
    bool ELIM_MLIN;
    //! @brief Set of invertible nonlinear operations
    std::set<FFInv::Options::NLINV> ELIM_NLIN;
  } options;


protected:
  //! @brief pointer to underlying dag
  //DAG* _dag;

  //! @brief Map between eliminated DAG variables and coresponding equality constraint
  t_VarElim _VarElim;

  //! @brief Independent DAG variables participating in equality constraints
  std::vector<FFVar> _Var;

  //! @brief Independent DAG variables with associated weights/priorities
  std::map<unsigned,double> _VarWeight;

  //! @brief Candidate variables
  std::set<unsigned> _ndxVar;

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
#endif

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

  //! @brief Insert an auxiliary variable corresponding to operand of inverse operation into DAG
  FFVar const* _insert_expr
    ( FFVar& var, t_poly const& polyindep, t_poly const& polydep, bool const useprod,
      bool const inverse = false );

  //! @brief Insert an auxiliary variable corresponding to inverse operation into DAG
  FFVar const* _insert_expr
    ( FFVar& var, FFOp const* op, t_poly const& polyindep, t_poly const& polydep,
      bool const useprod );

  //! @brief Erase all entries in _Interm
  void _reset
    ();
};

template <typename DAG>
inline typename SElimEnv<DAG>::Options SElimEnv<DAG>::options;

#if defined(MC__USE_GUROBI)
template <typename DAG>
inline int const SElimEnv<DAG>::Options::LPALGO_DEFAULT;
#endif
/*
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
*/

template < typename DAG >
inline void
SElimEnv<DAG>::_reset
()
{
  _Var.clear();
  _IVar.clear();
  _VarWeight.clear();
  _ndxVar.clear();
  _mapVarCtrCand.clear();
  _mapVarElVar.clear();
  
  _Ctr.clear();
  _ICtr.clear();
  _ndxCtr.clear();
  _mapCtrVar.clear();
  _mapCtrVarCand.clear();

  _VarElim.clear();
}

template < typename DAG >
inline void
SElimEnv<DAG>::process
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

template < typename DAG >
inline void
SElimEnv<DAG>::process
( unsigned const nCtr, FFVar const* pCtr,
  std::map<FFVar const*,double,lt_FFVar> const& wVar,
  bool const add2dag )
{
  // Reset internal variables
  _reset();

  // Update participating variables in _Var
  auto sgCtr = _dag->subgraph( nCtr, pCtr );
  unsigned v=0;
  for( auto&& Op : sgCtr.l_op ){
    if( Op->type != FFOp::VAR ) continue;
    _Var.push_back( *Op->pres );
    _IVar.push_back( FFInv().indep( v ) );
    auto itv = wVar.find( Op->pres ); 
    _VarWeight[v] = ( itv != wVar.end() ? itv->second : 1. );
    ++v;
  }
#ifdef MC__SELIM_DEBUG_PROCESS
  std::cout << std::endl << _Var.size() << " Original Variables: ";
  for( auto&& var : _Var ) std::cout << var << " ";
  std::cout << std::endl;
#endif

  // Find out variable invertibility in equality constraints
  _Ctr.assign( pCtr, pCtr+nCtr );
  _ICtr.resize( nCtr );
  FFInv::options.INVOP = options.ELIM_NLIN;
  _dag->eval( sgCtr, nCtr, pCtr, _ICtr.data(), _Var.size(), _Var.data(), _IVar.data() );

  // Create candidate variable and constraint sets and maps
  for( unsigned j=0; j<_ICtr.size(); ++j ){
#ifdef MC__SELIM_DEBUG_PROCESS
    std::cout << j << ": " << _ICtr[j] << std::endl;
#endif
    bool isinvert = false;
    for( auto const& [i,type] : _ICtr[j].inv() ){
      switch( type ){
        case FFInv::L:
          if( options.ELIM_LIN ){
	    _ndxVar.insert( i );
	    isinvert = true;
	  }
          break;
	case FFInv::S:
          if( options.ELIM_MLIN ){
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
        case FFInv::L: if( options.ELIM_LIN )  insert = true; break;
	case FFInv::S: if( options.ELIM_MLIN ) insert = true; break;
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
  _MIP_optimize();

#ifdef MC__SELIM_DEBUG_PROCESS
  std::cout << std::endl << "Eliminated variables and corresponding constraints:" << std::endl;
  for( auto const& [pvar,pctr] : _VarElim )
    std::cout << *pvar << " <- " << *pctr << std::endl;
#endif
  
  // No transcription in DAG if <a>add2dag</a> is false
  //if( !add2dag ) return;

  t_poly::options.BASIS = t_poly::Options::MONOM;
  t_lift::options.LIFTDIV  = true;
  t_lift::options.LIFTIPOW = false;

  for( auto&& [pvarel,pctr] : _VarElim ){
    t_lift::process( 1, pctr, false );
#ifdef MC__SELIM_DEBUG_PROCESS
    std::cout << *this;
#endif
    // Initialize lifted variable to current dependent
    auto const* spvar = &_SPDep.front();

    // Initialize lifted variable and loop over intermediate operations
    for( FFVar VarEl = 0.; ; ){

      // Check denominator is 1
      auto const& denom = spvar->denom();
      if( denom.maxord() || denom.coefmon( t_mon() ) != 1e0 )
        throw Exceptions( Exceptions::INVERT );

      // Check numerator dependence in pvarel
      auto const& numer = spvar->numer();
      long ndxvarel = pvarel->id().second;
      FFVar const* pdep = nullptr;
      for( auto const& pvar : numer.setvar() ){
        if( pvar->dep().dep( ndxvarel ).first ){
          if( pdep ) throw Exceptions( Exceptions::INVERT );
          pdep = pvar;
        }  
      }
      if( pdep == nullptr ) throw Exceptions( Exceptions::INVERT );

      // Split numerator into dependent and independent parts
      t_mon mondep( pdep );
      t_poly polydep, polyindep;
      for( auto const& [mon,coef] : numer.mapmon() ){
        if( mondep.subseteq( mon ) )
          polydep += std::make_pair( mon - mondep, coef );
        else
          polyindep -= std::make_pair( mon, coef );
      }
#ifdef MC__SELIM_DEBUG_PROCESS
      std::cout << "Dependent variable: " << *pdep << std::endl;
      std::cout << "Dependent part:" << polydep;
      std::cout << "Independent part:" << polyindep;
#endif

      // Case dependent is a DAG variable
      if( pdep->id().first == FFVar::VAR ){
        // Create new DAG expression and copy pointer to expression in _VarElim
        pctr = _insert_expr( VarEl, polyindep, polydep, false );
        break; // EXIT THE INNER FOR LOOP ON INTERMEDIATE VARIABLES
      }

      // Case dependent is not a DAG variable, identify intermediate operation
      bool found_Op = false;
      for( auto const& [pOp,SPVar] : _Interm ){
        if( pOp->pres->id().first  != FFVar::AUX ) throw Exceptions( Exceptions::INVERT );
        if( pOp->pres->id().second != pdep->id().second ) continue;
        found_Op = true;
        // Create new DAG expression and update pctr with pointer to expression
        pctr  = _insert_expr( VarEl, pOp, polyindep, polydep, false );
	spvar = SPVar.at(0);
        break;
      }
      if( !found_Op ) throw Exceptions( Exceptions::INVERT );
    }
  }
}

template < typename DAG >
inline FFVar const*
SElimEnv<DAG>::_insert_expr
( FFVar& var, FFOp const* pOp, t_poly const& polyindep, t_poly const& polydep, bool const useprod )
{
  // Insert inverse operation in DAG
  switch( pOp->type ){
   case FFOp::IPOW:  _insert_expr( var, polyindep, polydep, useprod );
                     var = pow( var, 1./(double)pOp->pops.at(1)->num().n ); break;
   case FFOp::DPOW:  _insert_expr( var, polyindep, polydep, useprod );
                     var = pow( var, 1./pOp->pops.at(1)->num().val() );     break;

   case FFOp::INV:   _insert_expr( var, polyindep, polydep, useprod, true );
                     break;
		     
   case FFOp::SQR:   _insert_expr( var, polyindep, polydep, useprod );
                     var = sqrt( var ); break;
   case FFOp::SQRT:  _insert_expr( var, polyindep, polydep, useprod );
                     var = sqr( var );  break;
   case FFOp::EXP:   _insert_expr( var, polyindep, polydep, useprod );
                     var = log( var );  break;
   case FFOp::LOG:   _insert_expr( var, polyindep, polydep, useprod );
                     var = exp( var );  break;
   case FFOp::COS:   _insert_expr( var, polyindep, polydep, useprod );
                     var = acos( var ); break;
   case FFOp::SIN:   _insert_expr( var, polyindep, polydep, useprod );
                     var = asin( var ); break;
   case FFOp::TAN:   _insert_expr( var, polyindep, polydep, useprod );
                     var = atan( var ); break;
   case FFOp::ACOS:  _insert_expr( var, polyindep, polydep, useprod );
                     var = cos( var );  break;
   case FFOp::ASIN:  _insert_expr( var, polyindep, polydep, useprod );
                     var = sin( var );  break;
   case FFOp::ATAN:  _insert_expr( var, polyindep, polydep, useprod );
                     var = tan( var );  break;

   case FFOp::CHEB:
   case FFOp::DIV:
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

template < typename DAG >
inline FFVar const*
SElimEnv<DAG>::_insert_expr
( FFVar& var, t_poly const& polyindep, t_poly const& polydep, bool const useprod,
  bool const inverse )
{
  // Case inverse operation is needed
  if( inverse ){
    var += insert_dag( polyindep, useprod );
    var =  insert_dag( polydep, useprod ) / var;
  }
  
  // Case dependent term is constant
  else if( !polydep.maxord() ){
    double const cst = polydep.coefmon( t_mon() );
    // Constant dependent term is one
    if( cst == 1e0 ){
      var += insert_dag( polyindep, useprod );
    }
    // Constant dependent term is negative one
    else if( cst == -1e0 ){
      var = insert_dag( - polyindep, useprod ) - var;
    }
    else{
      var += insert_dag( polyindep, useprod );
      var *=  cst;
    }
  }

  // Case dependent term is monomial
  else if( polydep.nmon() == 1 ){
    // Check if common divider
    bool common_div = true;
    t_mon const& mondep = polydep.mapmon().cbegin()->first;
    for( auto const& [mon,coef] : polyindep.mapmon() ){
      if( mondep.subseteq( mon ) ) continue;
      common_div = false;
      break;
    }
    // If common divider, simplify polynomial division first
    if( common_div ){
      t_poly polyindepred;
      auto it = polyindep.mapmon().cbegin();
      for( ; it!=polyindep.mapmon().cend(); ++it ){
        auto const& [mon,coef] = *it;
        polyindepred += std::make_pair( mon-mondep, coef );
      }
      var /= insert_dag( mondep, useprod );
      var += insert_dag( polyindepred, useprod );
    }
    // No common divider
    else{
      var += insert_dag( polyindep, useprod );
      var /= insert_dag( mondep, useprod );
    }
  }
  
  // Case dependent term is linear
  // General case
  else{
    var += insert_dag( polyindep, useprod );
    var /= insert_dag( polydep, useprod );
  }

  // Return pointer to internal DAG variable
  auto itvar = _dag->Vars().find( &var );
#ifdef MC__SELIM_CHECK
  assert( itvar != _dag->Vars().end() );
#endif
#ifdef MC__SELIM_DEBUG_INSERT
  std::cout << "INSERTED DAG EXPRESSION:";
  _dag->output( _dag->subgraph( 1, *itvar ) );
#endif
  return *itvar;
}

//#if defined(MC__USE_GUROBI)
template <typename DAG>
inline void
SElimEnv<DAG>::_MIP_optimize
()
{
  _MIPexcpt = false;

  try{
    // Run MIP optimization for a minimal representation
    _MIP_reset();
    _MIP_encode();
    //_MIP_initialize( minOrd );
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

template <typename DAG>
inline void
SElimEnv<DAG>::_MIP_solve
()
{
#if defined(MC__USE_GUROBI)
  _GRBmodel->update();
//#ifdef MC__SELIM_DEBUG_MIP
//  _MIP_display( false );
//#endif
  if( options.MIPOUTPUTFILE != "" )
    _GRBmodel->write( options.MIPOUTPUTFILE );
  fedisableexcept(FE_ALL_EXCEPT);
  _GRBmodel->set( GRB_IntAttr_ModelSense, -1 );
  _GRBmodel->optimize();
  //if( options.MIPDISPLEVEL )
  //  std::cout << "  #eliminated variables: " << _GRBmodel->get( GRB_DoubleAttr_ObjVal ) << std::endl;
  //_MIP_display( true );
#else
  throw Exceptions( Exceptions::MIPERR );
#endif
}

template <typename DAG>
inline void
SElimEnv<DAG>::_MIP_encode
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

template <typename DAG>
inline void
SElimEnv<DAG>::_MIP_decode
()
{
  // Add entries in _VarElim for eliminated variables and corresponding constraints
#ifdef MC__SELIM_DEBUG_MIP
  for( auto const& [v,grbvar] : _MIP_var ){
    if( grbvar.get(GRB_DoubleAttr_X) < 0.9 ) continue;
    std::cout << "Variable " << v << " eliminated\n";
  }
#endif
  for( auto const& [e,vmap] : _MIP_ctrvar ){
    if( _MIP_ctr[e].get(GRB_DoubleAttr_X) < 0.9 ) continue;
#ifdef MC__SELIM_DEBUG_MIP
    std::cout << "Constraint " << e << " used for elimination\n";
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
      _VarElim[*itdagvar] = *itdagctr;
    }
  }
}

template < typename DAG >
inline void
SElimEnv<DAG>::_MIP_options
()
{
#if defined(MC__USE_GUROBI)
  // Gurobi options
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
#endif
}

template < typename DAG >
inline void
SElimEnv<DAG>::_MIP_reset
()
{
  _MIP_var.clear();
  _MIP_varMTZ.clear();
  _MIP_ctr.clear();
  _MIP_ctrvar.clear();
  _MIP_varvar.clear();
#if defined(MC__USE_GUROBI)
  delete _GRBmodel;
  _GRBmodel = new GRBModel( *_GRBenv );
#endif
}
//#endif

} // namespace mc

#endif
