#define TEST_EXP       // <-- select test function here
const int NX = 200;     // <-- select discretization here
const int NE = 5;       // <-- select polynomial model expansion here
#define SAVE_RESULTS    // <-- specify whether to save results to file
#undef  USE_POLYMOD     // <-- specify whether to use a Chebyshev expansion before relaxation
#define  ADD_BREAKPOINT  // <-- specify whether to add breakpoints to the variables
const int NDIV = 2;     // <-- select number of breakpoints
#define USE_CMODEL	    // <-- Use Chebyshev models?
#define USE_MIP         // <-- specify whether to use piecewise-linear cuts

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

using namespace mc;

////////////////////////////////////////////////////////////////////////

#if defined( TEST_FABS )
const double XL   = -1.;	// <-- range lower bound
const double XU   =  2.;	// <-- range upper bound
template <class T>
T myfunc
( const T&x )
{
  return fabs(sqr(x)-1);
}

#elif defined( TEST_XLOG )
const double XL   =  .2;	// <-- range lower bound
const double XU   =  3.;	// <-- range upper bound
template <class T>
T myfunc
( const T&x )
{
  //return xlog(x);
  return sqr(x);
}

#elif defined( TEST_DPOW )
const double XL   =  .2;	// <-- range lower bound
const double XU   =  3.;	// <-- range upper bound
template <class T>
T myfunc
( const T&x )
{
  return pow(x,1.6);
}

#elif defined( TEST_2NORM )
const double XL   = -2.;	// <-- range lower bound
const double XU   =  1.;	// <-- range upper bound
template <class T>
T myfunc
( const T&x )
{
  return sqrt(fabs(x));
}

#elif defined( TEST_HILL )
const double XL   =  0.05;	// <-- range lower bound
const double XU   =  5.;	// <-- range upper bound
template <class T>
T myfunc
( const T&x )
{
  return 1./(1.+pow(x,-6));
}

#elif defined( TEST_HILL2 )
const double XL   =  0.6;	// <-- range lower bound
const double XU   =  1.2;	// <-- range upper bound
template <class T>
T myfunc
( const T&x )
{
  return pow(x,3)/(1+pow(x,3));
}

#elif defined( TEST_EXP )
const double XL   = -2.;	// <-- range lower bound
const double XU   =  1.;	// <-- range upper bound
template <class T>
T myfunc
( const T&x )
{
  return exp(x);
  //return x*exp(-pow(x,2));
}

#elif defined( TEST_EXP2 )
const double XL   = -2.;	// <-- range lower bound
const double XU   =  1.;	// <-- range upper bound
template <class T>
T myfunc
( const T&x )
{
  return x*exp(-fabs(pow(x,3)));
}

#elif defined( TEST_ERF )
const double XL   = -3.;	// <-- range lower bound
const double XU   =  3.;	// <-- range upper bound
template <class T>
T myfunc
( const T&x )
{
  return erf(x);
}

#elif defined( TEST_ERFC )
const double XL   = -2.;	// <-- range lower bound
const double XU   =  3.;	// <-- range upper bound
template <class T>
T myfunc
( const T&x )
{
  return x*erfc(-pow(x,2));
}

#elif defined( TEST_TRIG )
const double XL   =  mc::PI/2.;	// <-- range lower bound
const double XU   =  5.*mc::PI/4.;	// <-- range upper bound
template <class T>
T myfunc
( const T&x )
{
  return sin(x);
}

#elif defined( TEST_TRIG2 )
const double XL   =  mc::PI/6.;	// <-- range lower bound
const double XU   =  mc::PI/3.;	// <-- range upper bound
template <class T>
T myfunc
( const T&x )
{
  return cos(pow(x,2))*sin(pow(x,-3));
}

#elif defined( TEST_TRIG3 )
const double XL   =  -mc::PI;	// <-- range lower bound
const double XU   =  mc::PI;	// <-- range upper bound
template <class T>
T myfunc
( const T&x )
{
  return tan(cos(x*atan(x)));
}

#elif defined( TEST_TRIG4 )
const double XL   = -0.95;	// <-- range lower bound
const double XU   =  0.95;	// <-- range upper bound
template <class T>
T myfunc
( const T&x )
{
  return asin(x)+acos(x);
}

#elif defined( TEST_XLOG )
const double XL   =  0.1;	// <-- range lower bound
const double XU   =  0.9;	// <-- range upper bound
template <class T>
T myfunc
( const T&x )
{
  return pow(-fabs(x-0.5)-xlog(x),3);
}

#elif defined( TEST_MIN1 )
const double XL   = -2;	// <-- range lower bound
const double XU   =  1;	// <-- range upper bound
using std::min;
template <class T>
T myfunc
( const T&x )
{
  return min(pow(x-1,2), 1.);
}

#elif defined( TEST_MIN2 )
const double XL   = -.8;	// <-- range lower bound
const double XU   =  .8;	// <-- range upper bound
template <class T>
T myfunc
( const T&x )
{
  T m[3] = { pow(x+.5,2), pow(x,2), pow(x-.5,2) };
  return min(3,m);
}

#elif defined( TEST_DISC )
const double XL   =  2.5;	// <-- range lower bound
const double XU   =  4.5;	// <-- range upper bound
template <class T>
T myfunc
( const T&x )
{
  return fstep(4-x)*(fstep(x-3)*(exp(4-x)+3-(fstep(3-x)*(-sqr(x-2.5)+4)))
                   +(fstep(3-x)*(-sqr(x-2.5)+4))-(2*x-7))+(2*x-7);
}

#elif defined( TEST_LTCOND )
const double XL   =  0.;	// <-- range lower bound
const double XU   =  1.;	// <-- range upper bound
template <class T>
T myfunc
( const T&x )
{
  return ltcond( sqr(x)-0.25, pow(x-0.5,3), pow(-x-0.25,3) );
}

#elif defined( TEST_GTCOND )
const double XL   =  0.;	// <-- range lower bound
const double XU   =  1.;	// <-- range upper bound
template <class T>
T myfunc
( const T&x )
{
  return gtcond( sqr(x)-.25, sqr(x)-0.8, 0.1-sqr(x) );
}

#elif defined( TEST_CHEB )
const double XL   = -1.;	// <-- range lower bound
const double XU   =  1.;	// <-- range upper bound
template <class T>
T myfunc
( const T&x )
{
  return cheb(x,2)+0.5*cheb(x,3)-0.3*cheb(x,4);
}

#elif defined( TEST_INTER )
const double XL   =  0.;	// <-- range lower bound
const double XU   =  1.;	// <-- range upper bound
bool inter( double&x, const double&y, const double&z )
{ x = std::max(y,z); return true; }
template <class T>
T myfunc
( const T&x )
{
  T xinter;
  assert( inter( xinter, x*pow(x,2), pow(x,2) ) );
  return xinter;
}

#elif defined( TEST_HULL )
const double XL   =  -1.;	// <-- range lower bound
const double XU   =   .5;	// <-- range upper bound
bool hull( const double&y, const double&z )
{ return 0.5*(y+z); }
template <class T>
T myfunc
( const T&x )
{
  return hull( sqrt(fabs(x-.5)), sqrt(fabs(x+.5)) );
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
      break;
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
  std::ofstream res( "LINREL_1d.out", std::ios_base::out );
  res << std::scientific << std::setprecision(5) << std::right;
#endif

  I IX = { XL, XU };

  mc::FFGraph DAG;
  mc::FFVar X( &DAG );

#ifndef USE_POLYMOD
  mc::FFVar F = myfunc( X );
  std::cout << DAG;
#ifdef SAVE_RESULTS
  DAG.output( DAG.subgraph( 1, &F ) );
  std::ofstream ofdag( "LINREL_1d.dot", std::ios_base::out );
  DAG.dot_script( 1, &F, ofdag );
  ofdag.close();
#endif

  mc::PolImg<I> PolEnv;
  PolEnv.options.AGGREG_LQ = true;
#ifndef USE_MIPRELAX
  PolEnv.options.RELAX_DISC = 1;
#else
  PolEnv.options.RELAX_DISC = 2;
#endif
  PolEnv.options.SANDWICH_MAXCUT = 6;
  mc::PolVar<I> X_Pol( &PolEnv, X, IX ), F_Pol;
#ifndef ADD_BREAKPOINT
  DAG.eval( 1, &F, &F_Pol, 1, &X, &X_Pol );
#else
  PolEnv.options.BREAKPOINT_TYPE = mc::PolImg<I>::Options::BIN;//SOS2;//NONE;
  DAG.eval( 1, &F, &F_Pol, 1, &X, &X_Pol );
  // Add breakpoints to all variables in DAG
  for( auto&& var : PolEnv.Vars() ){
    for( unsigned i=0; i<NDIV; i++ ){
      double pt = mc::Op<I>::l(var.second->range())
                + mc::Op<I>::diam(var.second->range())*(i+1.)/(NDIV+1.);
      var.second->add_breakpt( pt );
    }
  }
#endif

#else
  // Compute polynomial model
  PM PMenv( 1, NE );
  PV PMX( &PMenv, 0, IX );
  PV PMF = myfunc( PMX );
  std::cout << PMF;

  // Add polynomial model to DAG
  mc::FFVar Y[2];
  Y[0] = X;
  Y[1].set( &DAG );
#ifdef USE_CMODEL
  mc::FFVar**Xrbasis = PMenv.get_basis( NE, Y );
#else
  mc::FFVar**Xrbasis = PMenv.get_basis( Y );
#endif
  mc::FFVar F = PMenv.get_bound( PMF.coefmon().second, Xrbasis, &Y[1], PM::Options::NAIVE );
  std::cout << DAG;

  // Relax polynomial model range
  mc::PolImg<I> PolEnv;
  PolEnv.options.SANDWICH_MAXCUT = 6;
  PolEnv.options.AGGREG_LQ = true;
  mc::PolVar<I> X_Pol, Y_Pol[2], F_Pol;
  Y_Pol[0].set( &PolEnv, Y[0], IX );
  Y_Pol[1].set( &PolEnv, Y[1], PMF.R() );
  X_Pol = Y_Pol[0];
  DAG.eval( 1, &F, &F_Pol, 2, Y, Y_Pol );
#endif

  PolEnv.generate_cuts( 1, &F_Pol, true );
  std::cout << PolEnv;

#if defined( MC__USE_GUROBI )
  try{
    GRBEnv GRBenv;
    GRBModel GRBmodel( GRBenv );
    GRBmodel.getEnv().set( GRB_DoubleParam_FeasibilityTol, 1e-9 );
    GRBmodel.getEnv().set( GRB_DoubleParam_OptimalityTol,  1e-9 );
    GRBmodel.getEnv().set( GRB_IntParam_OutputFlag,        0    );

    auto itvar = PolEnv.Vars().find( &X_Pol.var() );
    auto itobj = PolEnv.Vars().find( &F_Pol.var() );
    auto jtvar = DAGVars.end(), jtobj = DAGVars.end();
    for( auto itv=PolEnv.Vars().begin(); itv!=PolEnv.Vars().end(); ++itv ){
      //std::cout << itv->second->name() << ": " << itv->second->has_cuts() << std::endl;
      GRBVar DAGVar = GRBmodel.addVar( mc::Op<I>::l(itv->second->range()),
            mc::Op<I>::u(itv->second->range()), 0.0, GRB_CONTINUOUS,
            itv->second->name() );
      auto jtv = DAGVars.insert( std::make_pair( itv->second, DAGVar ) );
      if( itv == itvar ) jtvar = jtv.first;
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

    for( int iX=0; iX<NX; iX++ ){
      fedisableexcept(FE_ALL_EXCEPT);

      double X_val = XL+iX*(XU-XL)/(NX-1.);
      double Z_val = myfunc( X_val);

      GRBmodel.update();
      GRBConstr cutval = GRBmodel.addConstr( GRBLinExpr( jtvar->second,  1. ), GRB_EQUAL, X_val );
      GRBmodel.setObjective( GRBLinExpr( jtobj->second,  1. ) );
      GRBmodel.set( GRB_IntAttr_ModelSense, 1 ); // MIN:1, MAX:-1

      GRBmodel.set( GRB_IntAttr_ModelSense, 1 ); // MIN:1, MAX:-1
      GRBmodel.update();
      GRBmodel.write( "LINREL_1d.lp" ); //return;
      GRBmodel.optimize();
      double Zcv = GRBmodel.get( GRB_DoubleAttr_ObjVal );

      GRBmodel.set( GRB_IntAttr_ModelSense, -1 ); // MIN:1, MAX:-1
      GRBmodel.update();
      GRBmodel.optimize();
      double Zcc = GRBmodel.get( GRB_DoubleAttr_ObjVal );

      GRBmodel.remove( cutval );

#ifdef SAVE_RESULTS
      res << std::setw(14) << X_val << std::setw(14) << Z_val
          << std::setw(14) << mc::Op<I>::l(F_Pol.range())
          << std::setw(14) << mc::Op<I>::u(F_Pol.range())
          << std::setw(14) << Zcv << std::setw(14) << Zcc
          << std::endl;
#endif
    }
  }
  
  catch(GRBException& ex) {
    std::cerr << "Error code = " << ex.getErrorCode() << std::endl;
    std::cerr << ex.getMessage() << std::endl;
  }
#endif

#ifdef USE_POLYMOD
  delete[] Xrbasis[0];
  delete[] Xrbasis;
#endif

  return 0;
}

