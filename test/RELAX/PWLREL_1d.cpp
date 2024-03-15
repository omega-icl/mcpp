#define TEST_EXP2       // <-- select test function here
const int NX = 100;     // <-- select discretization here
#define SAVE_RESULTS    // <-- specify whether to save results to file
const int NBKPT = 6;    // <-- select number of breakpoints

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
  //return exp(x-1);
  return x*exp(-pow(x,2));
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
      DAGCuts.insert( std::make_pair( pCut, model.addConstr( lhs, GRB_EQUAL, pCut->rhs() ) ) );
      break;
    case mc::PolCut<I>::LE:
      DAGCuts.insert( std::make_pair( pCut, model.addConstr( lhs, GRB_LESS_EQUAL, pCut->rhs() ) ) );
      break;
    case mc::PolCut<I>::GE:
      DAGCuts.insert( std::make_pair( pCut, model.addConstr( lhs, GRB_GREATER_EQUAL, pCut->rhs() ) ) );
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
  std::ofstream res( "PWLREL_1d.out", std::ios_base::out );
  res << std::scientific << std::setprecision(5) << std::right;
#endif

  mc::FFGraph DAG;
  mc::FFVar X( &DAG );
  mc::FFVar F = myfunc( X );
  std::cout << DAG;
#ifdef SAVE_RESULTS
  DAG.output( DAG.subgraph( 1, &F ) );
  std::ofstream ofdag( "PWLREL_1d.dot", std::ios_base::out );
  DAG.dot_script( 1, &F, ofdag );
  ofdag.close();
#endif

  I IX = { XL, XU };
  I IF;
  DAG.eval( 1, &F, &IF, 1, &X, &IX );

  mc::PolImg<I> PolEnv;
  mc::PolVar<I> PX( &PolEnv, X, IX );
  mc::PolVar<I> PF( &PolEnv, F, IF );

  PolEnv.options.BREAKPOINT_TYPE = mc::PolImg<I>::Options::BIN;//SOS2;//NONE;
  std::vector<double> Xk(NBKPT+2), Fk(NBKPT+2);
  Xk[0] = XL;
  for( unsigned i=1; i<NBKPT+1; i++ )
    Xk[i] = XL + (XU-XL)*i/(NBKPT+1.);
  Xk[NBKPT+1] = XU;
  for( unsigned i=0; i<NBKPT+2; i++ )
    DAG.eval( 1, &F, &Fk[i], 1, &X, &Xk[i] );
  PolEnv.add_semilinear_cuts( nullptr, NBKPT+2, PX, Xk.data(), PF, Fk.data(), mc::PolCut<I>::EQ ); 
  std::cout << PolEnv;

#if defined( MC__USE_GUROBI )
  try{
    GRBEnv GRBenv;
    GRBModel GRBmodel( GRBenv );
    GRBmodel.getEnv().set( GRB_DoubleParam_FeasibilityTol, 1e-9 );
    GRBmodel.getEnv().set( GRB_DoubleParam_OptimalityTol,  1e-9 );
    GRBmodel.getEnv().set( GRB_IntParam_OutputFlag,        0    );

    auto itvar = PolEnv.Vars().find( &PX.var() );
    auto itobj = PolEnv.Vars().find( &PF.var() );
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
      GRBmodel.write( "PWLREL_1d.lp" ); //return;
      GRBmodel.optimize();
      double Zcv = GRBmodel.get( GRB_DoubleAttr_ObjVal );

      GRBmodel.set( GRB_IntAttr_ModelSense, -1 ); // MIN:1, MAX:-1
      GRBmodel.update();
      GRBmodel.optimize();
      double Zcc = GRBmodel.get( GRB_DoubleAttr_ObjVal );

      GRBmodel.remove( cutval );

#ifdef SAVE_RESULTS
      res << std::setw(14) << X_val << std::setw(14) << Z_val
          << std::setw(14) << mc::Op<I>::l(PF.range())
          << std::setw(14) << mc::Op<I>::u(PF.range())
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

  return 0;
}

