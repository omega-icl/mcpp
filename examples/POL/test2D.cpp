#define TEST_EXP	// <-- select test function here
const int NX = 40;	// <-- select discretization here
#define SAVE_RESULTS    // <-- specify whether to save results to file

////////////////////////////////////////////////////////////////////////

#include <fstream>
#include <iomanip>

#include "interval.hpp"
typedef mc::Interval  I;

#ifdef USE_CMODEL
  #include "cmodel.hpp"
  typedef mc::CModel<I> PM;
  typedef mc::CVar<I> PV;
#else
  #include "tmodel.hpp"
  typedef mc::TModel<I> PM;
  typedef mc::TVar<I> PV;
#endif

#include "polimage.hpp"
#include "gurobi_c++.h"

typedef std::map< const mc::PolVar<I>*, GRBVar, mc::lt_PolVar<I> > t_GRBVar;
t_GRBVar DAGVars;
typedef std::map< const mc::PolCut<I>*, GRBConstr, mc::lt_PolCut<I> > t_GRBCut;
t_GRBCut DAGCuts;

extern "C"{
  #include <fenv.h>
  int fedisableexcept( int );
}

//using namespace mc;

////////////////////////////////////////////////////////////////////////

#if defined( TEST_EXP )
const double X0L   = -.5; // <-- range lower bound
const double X0U   =  .5; // <-- range upper bound
const double X1L   =  0.5; // <-- range lower bound
const double X1U   =  1.2; // <-- range upper bound
template <class T>
T myfunc
( const T*x )
{
  //return x[0]/x[1];
  return exp( 1 - x[0]/x[1] - x[0]*log(x[1]) );
}

#endif

////////////////////////////////////////////////////////////////////////

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
  }
}

////////////////////////////////////////////////////////////////////////

int main()
{

#ifdef SAVE_RESULTS
  std::ofstream res( "POLIMG-2D.out", std::ios_base::out );
  res << std::scientific << std::setprecision(5) << std::right;
#endif

  try{ 
    I IX[2] = { I(X0L,X0U), I(X1L,X1U) };

    mc::FFGraph DAG;
    mc::FFVar X[2]; X[0].set( &DAG ); X[1].set( &DAG );

    mc::FFVar F = myfunc( X );
    std::cout << DAG;
    //return 0;

    mc::PolImg<I> PolEnv;
    PolEnv.options.SANDWICH_MAXCUT = 5;
    mc::PolVar<I> X_Pol[2], F_Pol;
    X_Pol[0].set( &PolEnv, X[0], IX[0] );
    X_Pol[1].set( &PolEnv, X[1], IX[1] );
    DAG.eval( 1, &F, &F_Pol, 2, X, X_Pol );
    std::cout << PolEnv;

    //DAG.output( DAG.subgraph( 1, &F ) );
    //std::ofstream o_F( "F.dot", std::ios_base::out );
    //DAG.dot_script( 1, &F, o_F );
    //o_F.close();


    GRBEnv GRBenv;
    GRBModel GRBmodel( GRBenv );
    //GRBmodel.getEnv().set( GRB_IntParam_Method,            FPMethod );
    //GRBmodel.getEnv().set( GRB_IntParam_OutputFlag,        FPRelax<T>::options.SOLVER_DISPLAY );
    GRBmodel.getEnv().set( GRB_DoubleParam_FeasibilityTol, 1e-9 );
    GRBmodel.getEnv().set( GRB_DoubleParam_OptimalityTol,  1e-9 );
    //GRBmodel.getEnv().set( GRB_DoubleParam_MIPGap,         FPRelax<T>::options.MILP_RELGAP );
    //GRBmodel.getEnv().set( GRB_DoubleParam_MIPGapAbs,      FPRelax<T>::options.MILP_ABSGAP );
    //GRBmodel.getEnv().set( GRB_IntParam_Presolve, -1 );

    auto itx0 = PolEnv.Vars().find( &X_Pol[0].var() );
    auto itx1 = PolEnv.Vars().find( &X_Pol[1].var() );
    auto itobj = PolEnv.Vars().find( &F_Pol.var() );
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

    for( int iX0=0; iX0<NX; iX0++ ){ 
     for( int iX1=0; iX1<NX; iX1++ ){ 

      fedisableexcept(FE_ALL_EXCEPT);

      double X_val[2] = { X0L+iX0*(X0U-X0L)/(NX-1.), X1L+iX1*(X1U-X1L)/(NX-1.) };
      double Z_val = myfunc( X_val);

      GRBmodel.update();
      GRBConstr cut0val = GRBmodel.addConstr( GRBLinExpr( jtx0->second,  1. ), GRB_EQUAL, X_val[0] );
      GRBConstr cut1val = GRBmodel.addConstr( GRBLinExpr( jtx1->second,  1. ), GRB_EQUAL, X_val[1] );
      GRBmodel.setObjective( GRBLinExpr( jtobj->second,  1. ) );

      GRBmodel.set( GRB_IntAttr_ModelSense, 1 ); // MIN:1, MAX:-1
      GRBmodel.update();
      //GRBmodel.write( "POLIMG.lp" );
      GRBmodel.optimize();
      double Zcv = GRBmodel.get( GRB_DoubleAttr_ObjVal );
      //return 0;

      GRBmodel.set( GRB_IntAttr_ModelSense, -1 ); // MIN:1, MAX:-1
      GRBmodel.update();
      GRBmodel.optimize();
      double Zcc = GRBmodel.get( GRB_DoubleAttr_ObjVal );

      GRBmodel.remove( cut0val );
      GRBmodel.remove( cut1val );

#ifdef SAVE_RESULTS
      res << std::setw(14) << X_val[0] << std::setw(14) << X_val[1] << std::setw(14) << Z_val
          << std::setw(14) << mc::Op<I>::l(F_Pol.range())
          << std::setw(14) << mc::Op<I>::u(F_Pol.range())
          << std::setw(14) << Zcv << std::setw(14) << Zcc
          << std::endl;
#endif
     }
#ifdef SAVE_RESULTS
     res << std::endl;
#endif
    }

#ifdef USE_POLYMOD
    delete[] Xrcheb[0];
    delete[] Xrcheb;
#endif
  }
  catch(GRBException& ex) {
    std::cerr << "Error code = " << ex.getErrorCode() << std::endl;
    std::cerr << ex.getMessage() << std::endl;
  }

  return 0;
}





