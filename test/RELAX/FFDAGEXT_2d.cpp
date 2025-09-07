#define TEST_EXP2	    // <-- select test function here
const int NGRID = 40;	    // <-- select discretization here
#define SAVE_RESULTS    // <-- specify whether to save results to file

////////////////////////////////////////////////////////////////////////

#include <fstream>
#include <iomanip>

#ifdef MC__USE_PROFIL
 #include "mcprofil.hpp"
 typedef INTERVAL I;
#else
 #ifdef MC__USE_FILIB
  #include "mcfilib.hpp"
  typedef filib::interval<double,filib::native_switched,filib::i_mode_extended> I;
 #else
  #ifdef MC__USE_BOOST
   #include "mcboost.hpp"
   typedef boost::numeric::interval_lib::save_state<boost::numeric::interval_lib::rounded_transc_opp<double>> T_boost_round;
   typedef boost::numeric::interval_lib::checking_base<double> T_boost_check;
   typedef boost::numeric::interval_lib::policies<T_boost_round,T_boost_check> T_boost_policy;
   typedef boost::numeric::interval<double,T_boost_policy> I;
  #else
   #include "interval.hpp"
   typedef mc::Interval I;
  #endif
 #endif
#endif

#include "ffdagext.hpp"

#if defined( MC__USE_GUROBI )
 #include "gurobi_c++.h"
 extern "C"{
  #include <fenv.h>
  int fedisableexcept( int );
 }
 typedef std::map< const mc::PolVar<I>*, GRBVar, mc::lt_PolVar<I> > t_GRBVar;
 t_GRBVar DAGVars;
 typedef std::map< const mc::PolCut<I>*, GRBConstr, mc::lt_PolCut<I> > t_GRBCut;
 t_GRBCut DAGCuts;
#endif

//using namespace mc;

////////////////////////////////////////////////////////////////////////

#if defined( TEST_PROD )
const double X0L   = -.5; // <-- range lower bound
const double X0U   =  1.5; // <-- range upper bound
const double X1L   = -1.; // <-- range lower bound
const double X1U   =  0.5; // <-- range upper bound
template <class T>
T myfunc
( const T*x )
{
  return x[0]*x[1];
}

#elif defined( TEST_PROD2 )
const double X0L   = -.5; // <-- range lower bound
const double X0U   =  1.5; // <-- range upper bound
const double X1L   = -1.; // <-- range lower bound
const double X1U   =  0.5; // <-- range upper bound
template <class T>
T myfunc
( const T*x )
{
  return 0.25 * ( sqr( x[0]+x[1] ) - sqr( x[0]-x[1] ) );
}

#elif defined( TEST_PROD3 )
const double X0L   =  0.1; // <-- range lower bound
const double X0U   =  1.5; // <-- range upper bound
const double X1L   =  0.2; // <-- range lower bound
const double X1U   =  1.0; // <-- range upper bound
template <class T>
T myfunc
( const T*x )
{
  return exp( log(x[0]) + log(x[1]) );
}

#elif defined( TEST_FRAC )
const double X0L   = -.5; // <-- range lower bound
const double X0U   =  .5; // <-- range upper bound
const double X1L   =  0.5; // <-- range lower bound
const double X1U   =  1.5; // <-- range upper bound
template< class T >
T myfunc
( const T*x )
{
  return x[0]/x[1];
}

#elif defined( TEST_POLY )
const double X0L   = -1; // <-- range lower bound
const double X0U   =  1; // <-- range upper bound
const double X1L   = -1; // <-- range lower bound
const double X1U   =  1; // <-- range upper bound
template <class T>
T myfunc
( const T*x )
{
  return x[0] * ( x[1] + sqr(x[0]) + 1 );
}

#elif defined( TEST_FABS )
const double X0L   = -1.; // <-- range lower bound
const double X0U   =  2.; // <-- range upper bound
const double X1L   = -2.; // <-- range lower bound
const double X1U   =  1.; // <-- range upper bound
template< class T >
T myfunc
( const T*x )
{
  return sqrt(fabs(x[0]*x[1]));
}

#elif defined( TEST_EXP )
const double X0L   = -.5; // <-- range lower bound
const double X0U   =  .5; // <-- range upper bound
const double X1L   =  0.5; // <-- range lower bound
const double X1U   =  1.2; // <-- range upper bound
template <class T>
T myfunc
( const T*x )
{
  return exp( 1 - x[0]/x[1] - x[0]*log(x[1]) );
}

#elif defined( TEST_EXP2 )
const double X0L   = -.5; // <-- range lower bound
const double X0U   =  .5; // <-- range upper bound
const double X1L   =  0.5; // <-- range lower bound
const double X1U   =  1.2; // <-- range upper bound
template <class T>
T myfunc
( const T*x )
{
  return x[0]*x[1]*(x[0]*(exp(x[0])-exp(-x[0]))-x[1]*(exp(x[1])-exp(-x[1])));
}

#elif defined( TEST_EXP3 )
const double X0L   = -.5; // <-- range lower bound
const double X0U   =  .5; // <-- range upper bound
const double X1L   =  0.5; // <-- range lower bound
const double X1U   =  1.2; // <-- range upper bound
template <class T>
T myfunc
( const T*x )
{
  return x[0]*exp(x[0]+pow(x[1],2))-pow(x[1],2);
}

#elif defined( TEST_EXP4 )
const double X0L   = -.5; // <-- range lower bound
const double X0U   =  .5; // <-- range upper bound
const double X1L   =  0.5; // <-- range lower bound
const double X1U   =  1.2; // <-- range upper bound
template <class T>
T myfunc
( const T*x )
{
  return pow(x[0]*exp(fabs(x[0])/x[1]),3);
}

#elif defined( TEST_FSTEP )
const double X0L   = -1.; // <-- range lower bound
const double X0U   =  .5; // <-- range upper bound
const double X1L   = -.5; // <-- range lower bound
const double X1U   =  1.; // <-- range upper bound
template< class T >
T myfunc
( const T*x )
{
  return fstep(x[0]*x[1])*mc::sqr(x[0]) + bstep(x[0]*x[1])*pow(x[1],3);
}

#elif defined( TEST_TRIG )
const double X0L   = -1.; // <-- range lower bound
const double X0U   =  2.; // <-- range upper bound
const double X1L   = -2.; // <-- range lower bound
const double X1U   =  1.; // <-- range upper bound
template< class T >
T myfunc
( const T*x )
{
  return 1.+x[0]-sin(2.*x[0]+3.*x[1])-cos(3.*x[0]-5.*x[1]);
}

#elif defined( TEST_INV )
const double X0L   = 0.; // <-- range lower bound
const double X0U   = 7.; // <-- range upper bound
const double X1L   = 0.; // <-- range lower bound
const double X1U   = 7.; // <-- range upper bound
template< class T >
T myfunc
( const T*x )
{
  return -1./(pow(x[0]-4.,2) + pow(x[1]-4.,2) + 0.1)
         -1./(pow(x[0]-1.,2) + pow(x[1]-1.,2) + 0.2)
         -1./(pow(x[0]-8.,2) + pow(x[1]-8.,2) + 0.2);
}

#elif defined( TEST_MAX )
const double X0L   = -2.; // <-- range lower bound
const double X0U   =  1.; // <-- range upper bound
const double X1L   = -1.; // <-- range lower bound
const double X1U   =  2.; // <-- range upper bound
template< class T >
T myfunc
( const T*x )
{
  T f[3] = { -pow(x[0]*x[1]-2.,2), -pow(x[0]*x[1],2), -pow(x[0]*x[1]+2.,2) };
  return max( (unsigned)3, f );
}

#elif defined( TEST_DISJ )
const double X0L   = -2.; // <-- range lower bound
const double X0U   =  1.; // <-- range upper bound
const double X1L   = -1.; // <-- range lower bound
const double X1U   =  2.; // <-- range upper bound
template< class T >
T myfunc
( const T*x )
{
  return mc::Op<T>::max( x[0], x[1] ) + mc::Op<T>::min( 2*x[0], x[1] ) + fstep( x[0] - x[1] );
}

#endif

////////////////////////////////////////////////////////////////////////

#if defined( MC__USE_GUROBI )

void append_cut
( GRBModel &model, const mc::PolCut<I>*pCut )
{
  GRBVar* VarSOS = 0;
  double* WeiSOS = 0;
  int TypSOS = 0;
  if( pCut->type() == mc::PolCut<I>::SOS1 || pCut->type() == mc::PolCut<I>::SOS2 ){
    VarSOS = new GRBVar[pCut->nvar()];
    WeiSOS = new double[pCut->nvar()];
    TypSOS = ( pCut->type() == mc::PolCut<I>::SOS1? GRB_SOS_TYPE1: GRB_SOS_TYPE2 );
  }

  GRBLinExpr lhs;
  for( unsigned k=0; k<pCut->nvar(); k++ ){
    auto ivar = DAGVars.find( pCut->var()+k );
    if( ivar==DAGVars.end() ) throw std::runtime_error("variable not found");
    GRBVar&Var = ivar->second;
    lhs += GRBLinExpr( Var, pCut->coef()[k] );
    if( pCut->type() == mc::PolCut<I>::SOS1 || pCut->type() == mc::PolCut<I>::SOS2 ){
      VarSOS[k] = Var; WeiSOS[k] = (double)k/(double)pCut->nvar();
    }
  }

  switch( pCut->type() ){
    case mc::PolCut<I>::SOS1:
    case mc::PolCut<I>::SOS2:
      model.addSOS( VarSOS, WeiSOS, pCut->nvar(), TypSOS );
      delete [] VarSOS;
      delete [] WeiSOS;
    case mc::PolCut<I>::EQ:
      DAGCuts.insert( std::make_pair( pCut, model.addConstr( lhs,
        GRB_EQUAL, pCut->rhs() ) ) );
      break;
    case mc::PolCut<I>::LE:
      DAGCuts.insert( std::make_pair( pCut, model.addConstr( lhs,
        GRB_LESS_EQUAL, pCut->rhs() ) ) );
      break;
    case mc::PolCut<I>::GE:
      DAGCuts.insert( std::make_pair( pCut, model.addConstr( lhs,
        GRB_GREATER_EQUAL, pCut->rhs() ) ) );
      break;
    default:
      throw std::runtime_error("cut type not found");
  }
}

#endif

////////////////////////////////////////////////////////////////////////

int main()
{
#ifdef SAVE_RESULTS
  std::ofstream res( "FFDAGEXT_2d.out", std::ios_base::out );
  res << std::scientific << std::setprecision(5) << std::right;
#endif

  // Create DAG
  mc::FFGraph DAG;
  std::vector<mc::FFVar> X{ mc::FFVar( &DAG ), mc::FFVar( &DAG ) };

  mc::FFVar F = myfunc( X.data() );
  mc::DAGEXT<I> DAGF( &DAG, X, {F} );
  mc::FFDAGEXT<I> OpF;
  OpF.options.RELAX  = { OpF.options.PWCS, OpF.options.AUX };//MC };//INT };
  OpF.options.PWCDIV = 32;
  OpF.options.PWCREL = 1;
  OpF.options.PWCSUP.USE_SHADOW = 1;
  OpF.options.PWCSHADOW = 1;
  OpF.options.PWLINI = 4;
  OpF.options.PWLREL = 1;
  OpF.options.PWLSUP.USE_SHADOW = 1;
  OpF.options.PWLSHADOW = 1;
  std::vector<mc::FFVar> FF{ OpF( 0, X, &DAGF, 1 ) }; // with DAG copy

  auto SgF = DAG.subgraph( FF );
  DAG.output( SgF, " OF F" );
  auto StrF = mc::FFExpr::subgraph( &DAG, SgF );
  std::cout << "F: " << StrF[0] << std::endl;

  mc::PolImg<I> PolEnv;
  PolEnv.options.BREAKPOINT_TYPE = mc::PolImg<I>::Options::CONT;//BIN;//SOS2;
  PolEnv.options.AGGREG_LQ       = true;
  PolEnv.options.BREAKPOINT_RTOL =
  PolEnv.options.BREAKPOINT_ATOL = 0e0;
  //PolEnv.options.ALLOW_DISJ = { mc::FFOp::FSTEP, mc::FFOp::MAXF, mc::FFOp::MINF, mc::FFOp::FABS };
  //PolEnv.options.ALLOW_NLIN      = { mc::FFOp::TANH, mc::FFOp::EXP  };
  PolEnv.options.SANDWICH_MAXCUT = 6;

  std::vector<I> IX{ I(X0L,X0U), I(X1L,X1U) };
  std::vector<mc::PolVar<I>> PolX{ mc::PolVar<I>( &PolEnv, X[0], IX[0] ), mc::PolVar<I>( &PolEnv, X[1], IX[1] ) }, PolF;
  DAG.eval( SgF, FF, PolF, X, PolX );

  PolEnv.generate_cuts( PolF );
  std::cout << "\n Polyhedral relaxation:" << PolEnv << std::endl;

  return 0;
  
#if defined( MC__USE_GUROBI )
  try{ 
    GRBEnv GRBenv;
    GRBModel GRBmodel( GRBenv );
    GRBmodel.getEnv().set( GRB_DoubleParam_FeasibilityTol, 1e-9 );
    GRBmodel.getEnv().set( GRB_DoubleParam_OptimalityTol,  1e-9 );
    GRBmodel.getEnv().set( GRB_IntParam_OutputFlag,        0    );
    
    auto itx0 = PolEnv.Vars().find( &PolX[0].var() );
    auto itx1 = PolEnv.Vars().find( &PolX[1].var() );
    auto itobj = PolEnv.Vars().find( &PolF[0].var() );
    auto jtx0 = DAGVars.end(), jtx1 = DAGVars.end(), jtobj = DAGVars.end();

    for( auto itv=PolEnv.Vars().begin(); itv!=PolEnv.Vars().end(); ++itv ){
      GRBVar DAGVar = GRBmodel.addVar( mc::Op<I>::l(itv->second->range()),
          mc::Op<I>::u(itv->second->range()), 0.0, GRB_CONTINUOUS,
          itv->second->name() );
      auto jtv = DAGVars.insert( std::make_pair( itv->second, DAGVar ) );
      if( itv == itx0 ) jtx0 = jtv.first;
      if( itv == itx1 ) jtx1 = jtv.first;
      if( itv == itobj ) jtobj = jtv.first;
    }
    
    for( auto itv=PolEnv.Aux().begin(); itv!=PolEnv.Aux().end(); ++itv ){
      switch( (*itv)->id().first ){
       case mc::PolVar<I>::AUXCONT:{
        GRBVar DAGVar = GRBmodel.addVar( mc::Op<I>::l((*itv)->range()),
          mc::Op<I>::u((*itv)->range()), 0.0, GRB_CONTINUOUS, (*itv)->name() );
        DAGVars.insert( std::make_pair( *itv, DAGVar ) );
        break;
       }
       case mc::PolVar<I>::AUXINT:{
        GRBVar DAGVar = GRBmodel.addVar( mc::Op<I>::l((*itv)->range()),
          mc::Op<I>::u((*itv)->range()), 0.0, GRB_INTEGER, (*itv)->name() );
        DAGVars.insert( std::make_pair( *itv, DAGVar ) );
        break;
       }
       default:
        break;
      }
    }
    
    GRBmodel.update();
    for( auto itc=PolEnv.Cuts().begin(); itc!=PolEnv.Cuts().end(); ++itc )
      append_cut( GRBmodel, *itc );

    for( unsigned iX0=0; iX0<NGRID; iX0++ ){ 
     for( unsigned iX1=0; iX1<NGRID; iX1++ ){ 

      fedisableexcept(FE_ALL_EXCEPT);

      double X_val[2] = { X0L+iX0*(X0U-X0L)/(NGRID-1.), X1L+iX1*(X1U-X1L)/(NGRID-1.) };
      double Z_val = myfunc( X_val);

      GRBmodel.update();
      GRBConstr cut0val = GRBmodel.addConstr( GRBLinExpr( jtx0->second,  1. ), GRB_EQUAL, X_val[0] );
      GRBConstr cut1val = GRBmodel.addConstr( GRBLinExpr( jtx1->second,  1. ), GRB_EQUAL, X_val[1] );
      GRBmodel.setObjective( GRBLinExpr( jtobj->second,  1. ) );

      GRBmodel.set( GRB_IntAttr_ModelSense, 1 ); // MIN:1, MAX:-1
      GRBmodel.update();
      GRBmodel.write( "FFDAGEXT_2d.lp" );
      GRBmodel.optimize();
      double Zcv = GRBmodel.get( GRB_DoubleAttr_ObjVal );

      GRBmodel.set( GRB_IntAttr_ModelSense, -1 ); // MIN:1, MAX:-1
      GRBmodel.update();
      GRBmodel.optimize();
      double Zcc = GRBmodel.get( GRB_DoubleAttr_ObjVal );

      GRBmodel.remove( cut0val );
      GRBmodel.remove( cut1val );

#ifdef SAVE_RESULTS
      res << std::setw(14) << X_val[0] << std::setw(14) << X_val[1] << std::setw(14) << Z_val
          << std::setw(14) << mc::Op<I>::l(PolF[0].range())
          << std::setw(14) << mc::Op<I>::u(PolF[0].range())
          << std::setw(14) << Zcv << std::setw(14) << Zcc
          << std::endl;
#endif
     }
#ifdef SAVE_RESULTS
     res << std::endl;
#endif
    }
  }
  catch(GRBException& ex) {
    std::cerr << "Error code = " << ex.getErrorCode() << std::endl;
    std::cerr << ex.getMessage() << std::endl;
  }
#endif

  return 0;
}





