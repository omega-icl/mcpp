////////////////////////////////////////////////////////////////////////
#define AD		// <-- select model here
unsigned int NTE = 30;	// <-- select order of Taylor expansion here
unsigned int NTM = 5;	// <-- select Taylor model order here
#undef CREATE_DOT	// <-- specify whether to save the DAG constructs to dot files
#undef SHOW_RESULTS     // <-- specify whether to show the DAG constructs and evaluation results
#undef USE_PROFIL	// <-- specify to use PROFIL for interval arithmetic
#undef USE_FILIB	// <-- specify to use FILIB++ for interval arithmetic
////////////////////////////////////////////////////////////////////////

#define _DEBUG
#define _TRACE
#include "ffunc.hpp"
#include <fstream>
#include <iomanip>
#include <sys/time.h>

#ifdef USE_PROFIL
  #include "mcprofil.hpp"
  typedef INTERVAL I;
#else
  #ifdef USE_FILIB
    #include "mcfilib.hpp"
    typedef filib::interval<double> I;
  #else
    #include "interval.hpp"
    typedef mc::Interval I;
  #endif
#endif

#include "ellimage.hpp"
typedef mc::EllImg<I>   EI    ;
typedef mc::EllVar<I>   EV    ;
//typedef mc::Ellipsoid   E     ;
typedef CPPL::dcovector dcv   ;
typedef CPPL::dsymatrix dsm   ;

#include "tmodel.hpp"
typedef mc::TModel<I> TM;
typedef mc::TVar<I> TV;
typedef mc::TModel<EV> TME;
typedef mc::TVar<EV> TVE;

#include "cmodel.hpp"
typedef mc::CModel<I> CM;
typedef mc::CVar<I> CV;
typedef mc::CModel<EV> CME;
typedef mc::CVar<EV> CVE;

#include "odetaylor.hpp"

#if defined( LV )
////////////////////////////////////////////////////////////////////////
//                                         _
//             dx1dt = p*x1*(1-x2)          |  t in (0,tf]
//             dx2dt = p*x2*(x1-1)         _|
//             x1(0) = 1.2, x2(0) = 1.1
//             pmin <= p <= pmax
////////////////////////////////////////////////////////////////////////

const unsigned int NP = 1;
const unsigned int NX = 2;
const double PVAL[NP] = { 3. }; 
const      I PDOM[NP] = { I(2.95,3.05) }; 

class IVP : public mc::ODESTRUCT
{
  public:
    IVP(): mc::ODESTRUCT( NP, NX )
    {}

  template <typename TX, typename TP, typename TT>
  TX RHS
    ( const unsigned int ix, const TP*p, const TX*x, const TT&t, const unsigned int is )
    {
      assert( ix < NX );
      switch( ix ){
        case 0: return p[0]*x[0]*(1.-x[1]);
        case 1: return p[0]*x[1]*(x[0]-1.);
        default: throw std::runtime_error("invalid index");
      }
    }

  template <typename T>
  T IC
    ( const unsigned int ix, const T*p )
    {
      using mc::sqr;
      assert( ix < NX );
      switch( ix ){
        case 0: return 1.2+sqr(p[0]-3.);
        case 1: return 1.1+(p[0]-3.);
        default: throw std::runtime_error("invalid index");
      }
    }
};

#elif defined( AD )
////////////////////////////////////////////////////////////////////////
//  ANAEROBIC DIGESTION MODEL FROM THE PAPER BY BERNARD ET AL (2001)
////////////////////////////////////////////////////////////////////////
const unsigned int NP = 2;    // Number of parameters
const unsigned int NX = 6;    // Number of states
const double PVAL[NP] = { 1., 1. }; 
const      I PDOM[NP] = { I(0.95,1.05), I(0.95,1.05) }; 

namespace AD_param{
const double mumax1  = 1.2e0;	// /day
const double KS1     = 7.1e0;	// g(COD)/L
const double mumax2  = 0.74e0;	// /day
const double KS2     = 9.28e0;	// mmol/L
const double KI2     = 256e0;	// mmol/L
const double alpha   = 0.5;	// -
const double kLa     = 19.8e0;	// /day
const double k1      = 42.14e0;	// g(COD)/g(VSS)
const double k2      = 116.5e0;	// mmol/g(VSS)
const double k3      = 268.0e0;	// mmol/g(VSS)
const double k4      =  50.6e0;	// mmol/g(VSS)
const double k5      = 343.6e0;	// mmol/g(VSS)
const double k6      = 453.0e0;	// mmol/g(VSS)
const double kVFA    = 64e-3;	// g(COD)/mmol
const double Ka      = 1.5e-2;	// mmol/L
const double Kb      = 6.5e-4;	// mmol/L
const double KH      = 16e0;	// mmol/L/atm
const double S1in    = 5e0;	// g(COD)/L
const double S2in    = 80e0;	// mmol/L
const double Cin     = 0e0;	// mmol/L
const double Zin     = 50e0;	// mmol/L
const double PT      = 1e0;	// atm
const double X10     = .5e0;	// g(VSS)/L
const double X20     = 1e0;	// g(VSS)/L
const double S10     = 1e0;	// g(COD)/L
const double S20     = 5e0;	// mmol/L
const double C0      = 40e0;	// mmol/L
const double Z0      = 50e0;	// mmol/L
const double D       = 0.4e0;	// /day
}
using namespace AD_param;

class IVP : public mc::ODESTRUCT
{
  public:
    IVP(): mc::ODESTRUCT( NP, NX )
    {}

  template <typename TX, typename TP, typename TT>
  TX RHS
  (const unsigned int ix, const TP*p, const TX*x, const TT&t, const unsigned int is)
  {
    assert(ix < NX );
    TX mu1  = mumax1 * x[2] / ( x[2] + KS1 );
    TX mu2  = mumax2 * x[3] / ( x[3] + KS2 + x[3]*x[3] / KI2 );
    TX phi  = x[5] + x[3] - x[4] + KH * PT + k6 / kLa * mu2 * x[1];
    TX PCO2 = ( phi - sqrt( phi*phi - 4e0 * KH * PT * ( x[5] + x[3] - x[4] ) ) ) / ( 2e0 * KH );
    TX qCO2 = kLa * ( x[5] + x[3] - x[4] - KH * PCO2 );
    switch( ix ){
      case 0: return ( mu1 - alpha * D ) * x[0];
      case 1: return ( mu2 - alpha * D ) * x[1];
      case 2: return D * ( S1in - x[2] ) - k1 * mu1 * x[0];
      case 3: return D * ( S2in - x[3] ) + k2 * mu1 * x[0] - k3 * mu2 * x[1];
      case 4: return D * ( Zin - x[4] );
      case 5: return D * ( Cin - x[5] ) + k4 * mu1 * x[0] + k5 * mu2 * x[1] - qCO2;
      default: throw std::runtime_error("invalid size");
    }
  }

  template <typename T>
  T IC
  ( const unsigned int ix, const T*p )
  {
    assert( ix < NX );
    switch(ix){
      case 0: return X10 * p[0];
      case 1: return X20 * p[1];
      case 2: return S10;
      case 3: return S20;
      case 4: return Z0;
      case 5: return C0;
      default: throw std::runtime_error("invalid size");
    }
  }
};
#endif

double cpuclock()
{
  timeval time;
  gettimeofday(&time, 0) ;
  return time.tv_sec + time.tv_usec*1e-6;
}

int test_F_ODE
( const bool eval=true )
{
  double cputime = -cpuclock();

  // Construction and Display of ODE RHS
  mc::FFGraph FF;
  mc::FFVar t, P[NP], X[NX], F[NX];
  for( unsigned int i=0; i<NP; i++ ) P[i].set( &FF );
  for( unsigned int i=0; i<NX; i++ ) X[i].set( &FF );
  t.set( &FF );
  for( unsigned int i=0; i<NX; i++ ) F[i] = IVP().RHS( i, P, X, t, 0 );
#ifdef SHOW_RESULTS
  std::cout << FF;
#endif

  std::vector<const mc::FFVar*> vX, vF;
  for( unsigned int i=0; i<NX; i++ ) vX.push_back(&X[i]);
  for( unsigned int i=0; i<NX; i++ ) vF.push_back(&F[i]);
  std::list<const mc::FFOp*> F_op = FF.subgraph( vF );
#ifdef SHOW_RESULTS
  FF.output( F_op );
#endif
#ifdef CREATE_DOT
  std::ofstream o_F( "F.dot", std::ios_base::out );
  FF.dot_script( vF, o_F );
  o_F.close();
#endif

  // Construction and Display of Jacobian matrix ODE RHS
  std::vector<const mc::FFVar*> vdFdX = FF.FAD( vF, vX );
#ifdef SHOW_RESULTS
  std::cout << FF;
#endif
  std::list<const mc::FFOp*> dFdX_op = FF.subgraph( vdFdX );
#ifdef SHOW_RESULTS
  FF.output( dFdX_op );
#endif
#ifdef CREATE_DOT
  std::ofstream odFdX( "dFdX.dot", std::ios_base::out );
  FF.dot_script( vdFdX, odFdX );
  odFdX.close();
#endif

  // Display of ODE RHS and Jacobian
  std::vector<const mc::FFVar*> vF_dFdX = vF;
  vF_dFdX.insert( vF_dFdX.end(), vdFdX.begin(), vdFdX.end() );
#ifdef SHOW_RESULTS
  std::cout << FF;
#endif
  std::list<const mc::FFOp*> F_dFdX_op = FF.subgraph( vF_dFdX );
#ifdef SHOW_RESULTS
  FF.output( F_dFdX_op );
#endif
#ifdef CREATE_DOT
  std::ofstream oF_dFdX( "F_dFdX.dot", std::ios_base::out );
  FF.dot_script( vF_dFdX, oF_dFdX );
  oF_dFdX.close();
#endif

  // statistics
  cputime += cpuclock();
  std::cout << std::setw(10) << F_op.size()
            << std::setw(10) << dFdX_op.size()
            << std::setw(10) << F_dFdX_op.size()
            << std::fixed << std::setprecision(4)
            << std::setw(10) << cputime;
  if( !eval ){ std::cout << std::endl; return 0; }

  // Evaluate RHS Jacobian in real arithmetic
  double cpu_deval = -cpuclock();
  double XVAL[NX], tVAL(0.);
  for( unsigned int i=0; i<NX; i++ ) XVAL[i] = IVP().IC( i, PVAL );
  std::vector< std::pair<const mc::FFVar*,double> > dVar;
  for( unsigned int i=0; i<NX; i++ ) dVar.push_back( std::make_pair( X+i, XVAL[i] ) );
  for( unsigned int i=0; i<NP; i++ ) dVar.push_back( std::make_pair( P+i, PVAL[i] ) );
  dVar.push_back( std::make_pair( &t, tVAL ) );
  std::vector<double> vdFdX_d = FF.eval( dFdX_op, vdFdX, dVar );
#ifdef SHOW_RESULTS
  for( unsigned i=0; i<vdFdX_d.size(); i++ )
    std::cout << std::scientific << std::setprecision(5)
              << "  dFdX(" << i << ") = " << vdFdX_d[i] << std::endl;
#endif
  cpu_deval += cpuclock();

  // Evaluate TE coefficients in interval arithmetic
  double cpu_Ieval = -cpuclock();
  I XI[NX], tI(0.);
  for( unsigned int i=0; i<NX; i++ ) XI[i] = IVP().IC( i, PDOM );
  std::vector< std::pair<const mc::FFVar*,I> > IVar;
  for( unsigned int i=0; i<NX; i++ ) IVar.push_back( std::make_pair( X+i, XI[i] ) );
  for( unsigned int i=0; i<NP; i++ ) IVar.push_back( std::make_pair( P+i, PDOM[i] ) );
  IVar.push_back( std::make_pair( &t, tI ) );
  std::vector<I> vdFdX_I = FF.eval( dFdX_op, vdFdX, IVar );
#ifdef SHOW_RESULTS
  for( unsigned i=0; i<vdFdX_I.size(); i++ )
    std::cout << std::scientific << std::setprecision(5)
              << "  dFdX(" << i << ") = " << vdFdX_I[i] << std::endl;
#endif
  cpu_Ieval += cpuclock();

  // Evaluate TE coefficients in Taylor model arithmetic
  double cpu_TVeval = -cpuclock();
  TV PTV[NP], XTV[NX], tTV(0.);
  TM TM_env( NP, NTM );
  for( unsigned int i=0; i<NP; i++ ) PTV[i].set( &TM_env, i, PDOM[i] );
  for( unsigned int i=0; i<NX; i++ ) XTV[i] = IVP().IC( i, PTV );
  std::vector< std::pair<const mc::FFVar*,TV> > TVVar;
  for( unsigned int i=0; i<NX; i++ ) TVVar.push_back( std::make_pair( X+i, XTV[i] ) );
  for( unsigned int i=0; i<NP; i++ ) TVVar.push_back( std::make_pair( P+i, PTV[i] ) );
  TVVar.push_back( std::make_pair( &t, tTV ) );
  std::vector<TV> vdFdX_TV = FF.eval( dFdX_op, vdFdX, TVVar );
#ifdef SHOW_RESULTS
  for( unsigned i=0; i<TVF_TE.size(); i++ )
    std::cout << std::scientific << std::setprecision(5)
              << "  dFdX(" << i << ") = " << vdFdX_TV[i] << std::endl;
#endif
  cpu_TVeval += cpuclock();

  // Evaluate TE coefficients in Chebyshev model arithmetic
  double cpu_CVeval = -cpuclock();
  CV PCV[NP], XCV[NX], tCV(0.);
  CM CM_env( NP, NTM );
  for( unsigned int i=0; i<NP; i++ ) PCV[i].set( &CM_env, i, PDOM[i] );
  for( unsigned int i=0; i<NX; i++ ) XCV[i] = IVP().IC( i, PCV );
  std::vector< std::pair<const mc::FFVar*,CV> > CVVar;
  for( unsigned int i=0; i<NX; i++ ) CVVar.push_back( std::make_pair( X+i, XCV[i] ) );
  for( unsigned int i=0; i<NP; i++ ) CVVar.push_back( std::make_pair( P+i, PCV[i] ) );
  CVVar.push_back( std::make_pair( &t, tCV ) );
  std::vector<CV> vdFdX_CV = FF.eval( dFdX_op, vdFdX, CVVar );
#ifdef SHOW_RESULTS
  for( unsigned i=0; i<CVF_TE.size(); i++ )
    std::cout << std::scientific << std::setprecision(5)
              << "  dFdX(" << i << ") = " << vdFdX_CV[i] << std::endl;
#endif
  cpu_CVeval += cpuclock();

/*
  // Evaluate function in fadbad::F type with Taylor model arithmetic
  double cpu_FTVeval = -cpuclock();
  //TV PTV[NP], XTV[NX], tTV(0.);
  //TM TM_env( NP, NTM );
  //for( unsigned int i=0; i<NP; i++ ) PTV[i].set( &TM_env, i, PDOM[i] );
  //for( unsigned int i=0; i<NX; i++ ) XTV[i] = IVP().IC( i, PTV );

  fadbad::F<TV> PFTV[NP], XFTV[NX], tFTV(tTV);
  for( unsigned int i=0; i<NP; i++ ) PFTV[i] = PTV[i];
  for( unsigned int i=0; i<NX; i++ ){ XFTV[i] = XTV[i]; XFTV[i].diff( i, NX ); }

  std::vector< std::pair< const mc::FFVar*,fadbad::F<TV> > > FTVVar;
  for( unsigned int i=0; i<NX; i++ ) FTVVar.push_back( std::make_pair( X+i, XFTV[i] ) );
  for( unsigned int i=0; i<NP; i++ ) FTVVar.push_back( std::make_pair( P+i, PFTV[i] ) );
  FTVVar.push_back( std::make_pair( &t, tFTV ) );
  std::vector< fadbad::F<TV> > FTVF = FF.eval( F_op, vF, FTVVar );

#ifdef SHOW_RESULTS
  for( unsigned i=0; i<FTVF.size(); i++ )
    for( unsigned j=0; j<NX.size(); j++ )
      std::cout << std::scientific << std::setprecision(5)
                << "  dFdX(" << i << "," << j << ") = " << FTVF[i]->d(j) << std::endl;
#endif
  cpu_FTVeval += cpuclock();
*/
  // statistics
  cputime += cpuclock();
  std::cout << std::fixed << std::setprecision(6)
            << std::setw(10) << cpu_deval
            << std::setw(10) << cpu_Ieval
            << std::setw(10) << cpu_TVeval
            << std::setw(10) << cpu_CVeval
//            << std::setw(10) << cpu_FTVeval
            << std::endl;

  return 0;
}

int AD_F_ODE()
{
  // Evaluation of ODE RHS Jacobian in real arithmetic
  double cpu_deval = -cpuclock();
  double X_d[NX], t_d(0.);
  for( unsigned int i=0; i<NX; i++ ) X_d[i] = IVP().IC( i, PVAL );
  fadbad::F<double> P_Fd[NP], X_Fd[NX], t_Fd(t_d);
  for( unsigned int i=0; i<NP; i++ ) P_Fd[i] = PVAL[i];
  for( unsigned int i=0; i<NX; i++ ){ X_Fd[i] = X_d[i]; X_Fd[i].diff( i, NX ); }
  fadbad::F<double> F_Fd[NX];
  for( unsigned int i=0; i<NX; i++ ) F_Fd[i] = IVP().RHS( i, P_Fd, X_Fd, t_Fd, 0 );
#ifdef SHOW_RESULTS
  for( unsigned i=0; i<NX; i++ )
    for( unsigned j=0; j<NX; j++ )
      std::cout << std::scientific << std::setprecision(5)
                << "  dFdX(" << i << "," << j << ") = " << F_Fd[i].d(j) << std::endl;
#endif
  cpu_deval += cpuclock();

  // Evaluation of ODE RHS Jacobian in interval arithmetic
  double cpu_Ieval = -cpuclock();
  I X_I[NX], t_I(0.);
  for( unsigned int i=0; i<NX; i++ ) X_I[i] = IVP().IC( i, PDOM );
  fadbad::F<I> P_FI[NP], X_FI[NX], t_FI(t_I);
  for( unsigned int i=0; i<NP; i++ ) P_FI[i] = PDOM[i];
  for( unsigned int i=0; i<NX; i++ ){ X_FI[i] = X_I[i]; X_FI[i].diff( i, NX ); }
  fadbad::F<I> F_FI[NX];
  for( unsigned int i=0; i<NX; i++ ) F_FI[i] = IVP().RHS( i, P_FI, X_FI, t_FI, 0 );
#ifdef SHOW_RESULTS
  for( unsigned i=0; i<NX; i++ )
    for( unsigned j=0; j<NX; j++ )
      std::cout << std::scientific << std::setprecision(5)
                << "  dFdX(" << i << "," << j << ") = " << F_FI[i].d(j) << std::endl;
#endif
  cpu_Ieval += cpuclock();

  // Evaluation of ODE RHS Jacobian in Taylor model arithmetic
  double cpu_TVeval = -cpuclock();
  TV P_TV[NP], X_TV[NX], t_TV(0.);
  TM TM_env( NP, NTM );
  for( unsigned int i=0; i<NP; i++ ) P_TV[i].set( &TM_env, i, PDOM[i] );
  for( unsigned int i=0; i<NX; i++ ) X_TV[i] = IVP().IC( i, P_TV );
  fadbad::F<TV> P_FTV[NP], X_FTV[NX], t_FTV(t_TV);
  for( unsigned int i=0; i<NP; i++ ) P_FTV[i] = P_TV[i];
  for( unsigned int i=0; i<NX; i++ ){ X_FTV[i] = X_TV[i]; X_FTV[i].diff( i, NX ); }
  fadbad::F<TV> F_FTV[NX];
  for( unsigned int i=0; i<NX; i++ ) F_FTV[i] = IVP().RHS( i, P_FTV, X_FTV, t_FTV, 0 );
#ifdef SHOW_RESULTS
  for( unsigned i=0; i<NX; i++ )
    for( unsigned j=0; j<NX; j++ )
      std::cout << std::scientific << std::setprecision(5)
                << "  dFdX(" << i << "," << j << ") = " << F_FTV[i].d(j) << std::endl;
#endif
  cpu_TVeval += cpuclock();

  // Evaluation of ODE RHS Jacobian in Chebyshev model arithmetic
  double cpu_CVeval = -cpuclock();
  CV P_CV[NP], X_CV[NX], t_CV(0.);
  CM CM_env( NP, NTM );
  for( unsigned int i=0; i<NP; i++ ) P_CV[i].set( &CM_env, i, PDOM[i] );
  for( unsigned int i=0; i<NX; i++ ) X_CV[i] = IVP().IC( i, P_CV );
  fadbad::F<CV> P_FCV[NP], X_FCV[NX], t_FCV(t_CV);
  for( unsigned int i=0; i<NP; i++ ) P_FCV[i] = P_CV[i];
  for( unsigned int i=0; i<NX; i++ ){ X_FCV[i] = X_CV[i]; X_FCV[i].diff( i, NX ); }
  fadbad::F<CV> F_FCV[NX];
  for( unsigned int i=0; i<NX; i++ ) F_FCV[i] = IVP().RHS( i, P_FCV, X_FCV, t_FCV, 0 );
#ifdef SHOW_RESULTS
  for( unsigned i=0; i<NX; i++ )
    for( unsigned j=0; j<NX; j++ )
      std::cout << std::scientific << std::setprecision(5)
                << "  dFdX(" << i << "," << j << ") = " << F_FCV[i].d(j) << std::endl;
#endif
  cpu_CVeval += cpuclock();

  // statistics
  std::cout << std::fixed << std::setprecision(6)
            << std::setw(10) << cpu_deval
            << std::setw(10) << cpu_Ieval
            << std::setw(10) << cpu_TVeval
            << std::setw(10) << cpu_CVeval
            << std::endl;

  return 0;
}

int test_TE_ODE
( const unsigned int TEORDER, const bool eval=true )
{
  double cputime = -cpuclock();

  // Construction and Display of ODE RHS
  mc::FFGraph FF;
  mc::FFVar t, P[NP], X[NX], F[NX];
  for( unsigned int i=0; i<NP; i++ ) P[i].set( &FF );
  for( unsigned int i=0; i<NX; i++ ) X[i].set( &FF );
  t.set( &FF );
  for( unsigned int i=0; i<NX; i++ ) F[i] = IVP().RHS( i, P, X, t, 0 );
#ifdef SHOW_RESULTS
  std::cout << FF;
#endif

  std::vector<const mc::FFVar*> vX, vF;
  for( unsigned int i=0; i<NX; i++ ) vX.push_back(&X[i]);
  for( unsigned int i=0; i<NX; i++ ) vF.push_back(&F[i]);
  std::list<const mc::FFOp*> F_op = FF.subgraph( vF );
#ifdef SHOW_RESULTS
  FF.output( F_op );
#endif
#ifdef CREATE_DOT
  std::ofstream o_F( "F.dot", std::ios_base::out );
  FF.dot_script( vF, o_F );
  o_F.close();
#endif

  // Construction and Display of Taylor coefficients of ODE solutions
  std::vector<const mc::FFVar*> vF_TE = FF.TAD( TEORDER, vF, vX );
#ifdef SHOW_RESULTS
  std::cout << FF;
#endif
  std::list<const mc::FFOp*> F_TE_op = FF.subgraph( vF_TE );
#ifdef SHOW_RESULTS
  FF.output( F_TE_op );
#endif
#ifdef CREATE_DOT
  std::ofstream oF_TE( "F_TE.dot", std::ios_base::out );
  FF.dot_script( vF_TE, oF_TE );
  oF_TE.close();
#endif

  // Construction and Display of Jacobian matrix of Taylor coefficients of ODE solutions
  std::vector<const mc::FFVar*> vF_FTE = FF.FAD( vF_TE, vX );
#ifdef SHOW_RESULTS
  std::cout << FF;
#endif
  std::list<const mc::FFOp*> F_FTE_op = FF.subgraph( vF_FTE );
#ifdef SHOW_RESULTS
  FF.output( F_FTE_op );
#endif
#ifdef CREATE_DOT
  std::ofstream oF_FTE( "F_FTE.dot", std::ios_base::out );
  FF.dot_script( vF_FTE, oF_FTE );
  oF_FTE.close();
#endif

  // Display of Taylor coefficients of ODE solutions and their Jacobian matrix
  std::vector<const mc::FFVar*> vF_TE_FTE = vF_TE;
  vF_TE_FTE.insert( vF_TE_FTE.end(), vF_FTE.begin(), vF_FTE.end() );
#ifdef SHOW_RESULTS
  std::cout << FF;
#endif
  std::list<const mc::FFOp*> F_TE_FTE_op = FF.subgraph( vF_TE_FTE );
#ifdef SHOW_RESULTS
  FF.output( F_TE_FTE_op );
#endif
#ifdef CREATE_DOT
  std::ofstream oF_TE_FTE( "F_TE_FTE.dot", std::ios_base::out );
  FF.dot_script( vF_TE_FTE, oF_TE_FTE );
  oF_TE_FTE.close();
#endif

  // statistics
  cputime += cpuclock();
  std::cout << std::setw(4) << TEORDER
            << std::setw(10) << F_op.size()
            << std::setw(10) << F_TE_op.size()
            << std::setw(10) << F_FTE_op.size()
            << std::setw(10) << F_TE_FTE_op.size()
            << std::fixed << std::setprecision(4)
            << std::setw(10) << cputime;
  if( !eval ){ std::cout << std::endl; return 0; }

  // Evaluate TE coefficients in real arithmetic
  double cpu_deval = -cpuclock();
  double XVAL[NX], tVAL(0.);
  for( unsigned int i=0; i<NX; i++ ) XVAL[i] = IVP().IC( i, PVAL );
  std::vector< std::pair<const mc::FFVar*,double> > dVar;
  for( unsigned int i=0; i<NX; i++ ) dVar.push_back( std::make_pair( X+i, XVAL[i] ) );
  for( unsigned int i=0; i<NP; i++ ) dVar.push_back( std::make_pair( P+i, PVAL[i] ) );
  dVar.push_back( std::make_pair( &t, tVAL ) );
  std::vector<double> dF_TE = FF.eval( F_TE_op, vF_TE, dVar );
#ifdef SHOW_RESULTS
  for( unsigned i=0; i<dF_TE.size(); i++ )
    std::cout << std::scientific << std::setprecision(5)
              << "  F_TE(" << i << ") = " << dF_TE[i] << std::endl;
#endif
  cpu_deval += cpuclock();

  // Evaluate TE coefficients in interval arithmetic
  double cpu_Ieval = -cpuclock();
  I XI[NX], tI(0.);
  for( unsigned int i=0; i<NX; i++ ) XI[i] = IVP().IC( i, PDOM );
  std::vector< std::pair<const mc::FFVar*,I> > IVar;
  for( unsigned int i=0; i<NX; i++ ) IVar.push_back( std::make_pair( X+i, XI[i] ) );
  for( unsigned int i=0; i<NP; i++ ) IVar.push_back( std::make_pair( P+i, PDOM[i] ) );
  IVar.push_back( std::make_pair( &t, tI ) );
  std::vector<I> IF_TE = FF.eval( F_TE_op, vF_TE, IVar );
#ifdef SHOW_RESULTS
  for( unsigned i=0; i<IF_TE.size(); i++ )
    std::cout << std::scientific << std::setprecision(5)
              << "  F_TE(" << i << ") = " << IF_TE[i] << std::endl;
#endif
  cpu_Ieval += cpuclock();

  // Evaluate TE coefficients in Taylor model arithmetic
  double cpu_TVeval = -cpuclock();
  TV PTV[NP], XTV[NX], tTV(0.);
  TM TM_env( NP, NTM );
  for( unsigned int i=0; i<NP; i++ ) PTV[i].set( &TM_env, i, PDOM[i] );
  for( unsigned int i=0; i<NX; i++ ) XTV[i] = IVP().IC( i, PTV );
  std::vector< std::pair<const mc::FFVar*,TV> > TVVar;
  for( unsigned int i=0; i<NX; i++ ) TVVar.push_back( std::make_pair( X+i, XTV[i] ) );
  for( unsigned int i=0; i<NP; i++ ) TVVar.push_back( std::make_pair( P+i, PTV[i] ) );
  TVVar.push_back( std::make_pair( &t, tTV ) );
  std::vector<TV> TVF_TE = FF.eval( F_TE_op, vF_TE, TVVar );
  //std::vector<TV> TVF_TE = FF.eval( vF_TE, TVVar );
#ifdef SHOW_RESULTS
  for( unsigned i=0; i<TVF_TE.size(); i++ )
    std::cout << std::scientific << std::setprecision(5)
              << "  F_TE(" << i << ") = " << TVF_TE[i] << std::endl;
#endif
  cpu_TVeval += cpuclock();

  // Evaluate function in fadbad::T type with Taylor model arithmetic
  double cpu_TTVeval = -cpuclock();
  //TV PTV[NP], XTV[NX], tTV(0.);
  //TM TM_env( NP, NTM );
  //for( unsigned int i=0; i<NP; i++ ) PTV[i].set( &TM_env, i, PDOM[i] );
  //for( unsigned int i=0; i<NX; i++ ) XTV[i] = IVP().IC( i, PTV );

  fadbad::T<TV> PTTV[NP], XTTV[NX], tTTV(tTV);
  for( unsigned int i=0; i<NP; i++ ) PTTV[i] = PTV[i];
  for( unsigned int i=0; i<NX; i++ ) XTTV[i] = XTV[i];

  std::vector< std::pair< const mc::FFVar*,fadbad::T<TV> > > TTVVar;
  for( unsigned int i=0; i<NX; i++ ) TTVVar.push_back( std::make_pair( X+i, XTTV[i] ) );
  for( unsigned int i=0; i<NP; i++ ) TTVVar.push_back( std::make_pair( P+i, PTTV[i] ) );
  TTVVar.push_back( std::make_pair( &t, tTTV ) );
  std::vector< fadbad::T<TV> > TTVF = FF.eval( F_op, vF, TTVVar );

  for( unsigned int q=0; q<TEORDER; q++ ){
    for( unsigned int i=0; i<TTVF.size(); i++ ){
      // Evaluate q-th Taylor coefficient for f[i]
      TTVF[i].eval(q);
      // Set result as (q+1)-th Taylor coefficient for x[i]
      XTTV[i][q+1] = TTVF[i][q] / double(q+1);
    }
  }
#ifdef SHOW_RESULTS
  for( unsigned int q=0; q<TEORDER; q++ )
    for( unsigned i=0; i<TTVF.size(); i++ )
      std::cout << std::scientific << std::setprecision(5)
                << "  F_TE(" << q << "," << i << ") = " << TTVF[i][q] << std::endl;
#endif
  cpu_TTVeval += cpuclock();

  // statistics
  cputime += cpuclock();
  std::cout << std::fixed << std::setprecision(6)
            << std::setw(10) << cpu_deval
            << std::setw(10) << cpu_Ieval
            << std::setw(10) << cpu_TVeval
            << std::setw(10) << cpu_TTVeval
            << std::endl;

  return 0;
}

int AD_TE_ODE
( const unsigned int TEORDER )
{
  // Evaluate function in real arithmetic using AD
  mc::ODETAYLOR<double,IVP> TE_double;
  double cpu_deval = -cpuclock();
  double XVAL[NX], tVAL(0.);
  for( unsigned int i=0; i<NX; i++ ) XVAL[i] = IVP().IC( i, PVAL );
  double* dF_TE[TEORDER+1];
  for( unsigned int q=0; q<=TEORDER; q++ ) dF_TE[q] = new double[NX];
  TE_double.T_expand( 0, tVAL, XVAL, PVAL, TEORDER, dF_TE );
#ifdef SHOW_RESULTS
  for( unsigned q=0; q<=TEORDER; q++ )
    for( unsigned i=0; i<NX; i++ )
      std::cout << std::scientific << std::setprecision(5)
                << "  F_TE(" << q << "," << i << ") = " << dF_TE[q][i] << std::endl;
#endif
  for( unsigned int q=0; q<=TEORDER; q++ ) delete[] dF_TE[q];
  cpu_deval += cpuclock();

  // Evaluate function in interval arithmetic using AD
  mc::ODETAYLOR<I,IVP> TE_I;
  double cpu_Ieval = -cpuclock();
  I XI[NX]; double tI(0.);
  for( unsigned int i=0; i<NX; i++ ) XI[i] = IVP().IC( i, PDOM );
  I* IF_TE[TEORDER+1];
  for( unsigned int q=0; q<=TEORDER; q++ ) IF_TE[q] = new I[NX];
  TE_I.T_expand( 0, tI, XI, PDOM, TEORDER, IF_TE );
#ifdef SHOW_RESULTS
  for( unsigned q=0; q<=TEORDER; q++ )
    for( unsigned i=0; i<NX; i++ )
      std::cout << std::scientific << std::setprecision(5)
                << "  F_TE(" << q << "," << i << ") = " << IF_TE[q][i] << std::endl;
#endif
  for( unsigned int q=0; q<=TEORDER; q++ ) delete[] IF_TE[q];
  cpu_Ieval += cpuclock();

  // Evaluate function in Taylor model arithmetic using AD
  mc::ODETAYLOR<TV,IVP> TE_TV;
  double cpu_TVeval = -cpuclock();
  TV PTV[NP], XTV[NX]; double tTV(0.);
  TM TM_env( NP, NTM );
  for( unsigned int i=0; i<NP; i++ ) PTV[i].set( &TM_env, i, PDOM[i] );
  for( unsigned int i=0; i<NX; i++ ) XTV[i] = IVP().IC( i, PTV );
  TV *TVF_TE[TEORDER+1], *TVF_FTE[TEORDER+1];
  for( unsigned int q=0; q<=TEORDER; q++ )
    { TVF_TE[q] = new TV[NX]; TVF_FTE[q] = new TV[NX*NX]; }
  TE_TV.T_expand( 0, tTV, XTV, PTV, TEORDER, TVF_TE );
  //TE_TV.TF_expand( 0, tTV, XTV, PTV, TEORDER, TVF_TE, TVF_FTE );
#ifdef SHOW_RESULTS
  for( unsigned q=0; q<=TEORDER; q++ )
    for( unsigned i=0; i<NX; i++ )
      std::cout << std::scientific << std::setprecision(5)
                << "  F_TE(" << q << "," << i << ") = " << TVF_TE[q][i] << std::endl;
#endif
  for( unsigned int q=0; q<=TEORDER; q++ )
    { delete[] TVF_TE[q]; delete[] TVF_FTE[q]; }
  cpu_TVeval += cpuclock();

  // statistics
  std::cout << std::setw(4) << TEORDER
            << std::fixed << std::setprecision(6)
            << std::setw(10) << cpu_deval
            << std::setw(10) << cpu_Ieval
            << std::setw(10) << cpu_TVeval
            << std::endl;

  return 0;
}

int test_TAD
( const unsigned int TEORDER )
{
  mc::FFGraph FF;
  mc::FFVar T( &FF );
  mc::FFVar X( &FF ); std::vector<const mc::FFVar*> vX; vX.push_back(&X);
  mc::FFVar F = T*T;  std::vector<const mc::FFVar*> vF; vF.push_back(&F);
  std::cout << FF;
  std::vector<const mc::FFVar*> vF_TE = FF.TAD( TEORDER, vF, vX, &T );
  std::cout << FF;

  std::list<const mc::FFOp*> FTE_op = FF.subgraph( vF_TE );
  FF.output( FTE_op );

  std::ofstream o_FTE( "FTE.dot", std::ios_base::out );
  FF.dot_script( vF_TE, o_FTE );
  o_FTE.close();

  return 0;
}

int test_eval()
{
      // DAG environment
      mc::FFGraph FF;

      // Independent variables
      const unsigned int NX = 4;
      mc::FFVar X[NX];
      std::vector<const mc::FFVar*> v_X;
      for( unsigned int i=0; i<NX; i++ ){
        X[i].set( &FF );
        v_X.push_back(&X[i]);
      }

      // Dependent variables
      const unsigned int NF = 2;
      mc::FFVar F[NF]
        = { X[2]*X[3]-X[0],
            X[0]*pow(exp(X[2]*X[3])+3.,4)+X[1] };
      std::vector<const mc::FFVar*> v_F;
      for( unsigned int j=0; j<NF; j++ )
        v_F.push_back( &F[j] );

      // DAG of second-order derivatives
      std::vector<const mc::FFVar*> v_d2FdX2 = FF.FAD( FF.FAD( v_F, v_X ), v_X );

      // Evaluation in interval arithmetic
      I IX[NX] = { I(0,0.5), I(1,2), I(-1,-0.8), I(0.5,1) };
      std::vector< std::pair<const mc::FFVar*,I> > v_IX;
      for( unsigned int i=0; i<NX; i++ )
        v_IX.push_back( std::make_pair( &X[i], IX[i] ) );
      std::vector<I> v_Id2FdX2 = FF.eval( v_d2FdX2, v_IX );

      // Display results
      for( unsigned i=0, k=0; i<v_Id2FdX2.size(); i++ ){
        std::cout << "  d2FdX2(" << i << ") = " << v_Id2FdX2[i] << std::endl;
        if( ++k == NX ){ std::cout << std::endl; k = 0; }
      }

      // Evaluation in 5th-order Chebyshev model arithmetic
      TM TM_env( NX, 5 );
      TV TVX[NX];
      std::vector< std::pair<const mc::FFVar*,TV> > v_TVX;
      for( unsigned int i=0; i<NX; i++ ){
        TVX[i].set( &TM_env, i, IX[i] );
        v_TVX.push_back( std::make_pair( &X[i], TVX[i] ) );
      }
      std::vector<TV> v_TVd2FdX2 = FF.eval( v_d2FdX2, v_TVX );

      // Display results
      for( unsigned i=0, k=0; i<v_TVd2FdX2.size(); i++ ){
        std::cout << "  d2FdX2(" << i << ") = " << v_TVd2FdX2[i].B() << std::endl;
        if( ++k == NX ){ std::cout << std::endl; k = 0; }
      }

  return 0;
}

int test_DAG()
{
  const unsigned int NX = 4, NF = 2;
  mc::FFGraph FF;
  mc::FFVar X[NX];

  // DAG construct and use
  for( unsigned int i=0; i<NX; i++ ) X[i].set( &FF );
  mc::FFVar F[NF]
    //= { X[1]*sqr(X[0]) };//,
    //    X[1]*X[2]*X[3] };
    //= { 2*X[2]*X[3],
    //    2*X[0]*X[1] };
    = { X[2]*X[3]-X[0],
        X[0]*pow(exp(X[2]*X[3])+3.,4)+X[1] };
  std::cout << FF;

  std::list<const mc::FFOp*> F_op  = FF.subgraph( NF, F );    FF.output( F_op );
  std::ofstream o_F( "F.dot", std::ios_base::out );
  FF.dot_script( NF, F, o_F );
  o_F.close();

  std::list<const mc::FFOp*> F0_op = FF.subgraph( 1, F );     FF.output( F0_op );
  std::ofstream o_F0( "F0.dot", std::ios_base::out );
  FF.dot_script( 1, F, o_F0 );
  o_F0.close();

  for( unsigned i=0; i<NF; i++ )
    std::cout << "Variable dependence of F[" << i << "]: " << F[i].dep() << std::endl;

  // 1st-order derivative DAG construct
  std::vector<const mc::FFVar*> v_X, v_F;
  for( unsigned int i=0; i<NX; i++ ) v_X.push_back( &X[i] );
  for( unsigned int j=0; j<NF; j++ ) v_F.push_back( &F[j] );
//  std::vector<const mc::FFVar*> v_dFdX = FF.FAD( v_F, v_X );
  std::vector<const mc::FFVar*> v_dFdX = FF.BAD( v_F, v_X );
  std::cout << FF;

  std::vector<const mc::FFVar*> v_F_dFdX = v_F; v_F_dFdX.insert( v_F_dFdX.end(), v_dFdX.begin(), v_dFdX.end() );
  std::list<const mc::FFOp*> op_F_dFdX = FF.subgraph( v_F_dFdX );
  std::ofstream o_F_dFdX( "F_dFdX.dot", std::ios_base::out );
  FF.dot_script( v_F_dFdX, o_F_dFdX );
  o_F_dFdX.close();
      
  std::vector<const mc::FFVar*> v_dF1dX3; v_dF1dX3.insert( v_dF1dX3.end(), v_dFdX[NX+3] );
  std::list<const mc::FFOp*> op_dF1dX3 = FF.subgraph( v_dF1dX3 );
  std::ofstream o_dF1dX3( "dF1dX3.dot", std::ios_base::out );
  FF.dot_script( v_dF1dX3, o_dF1dX3 );
  o_dF1dX3.close();

  return 0;
/*
  std::vector<const mc::FFVar*> l_F_dFdX;
  for( unsigned int i=0, ij=0; i<NF; i++ ){
    l_F_dFdX.push_back( vF[i] );
    for( unsigned int j=0; j<NX; j++, ij++ )
      l_F_dFdX.push_back( vdFdX[ij] );
  }
  std::list<const mc::FFOp*> F_dFdX_op = FF.subgraph( l_F_dFdX );
  FF.output( F_dFdX_op );

  std::ofstream o_F_dFdX( "F_dFdX.dot", std::ios_base::out );
  FF.dot_script( l_F_dFdX, o_F_dFdX );
  o_F_dFdX.close();

  for( unsigned i=0; i<NF*NX; i++ )
    std::cout << "Variable dependence of dFdX[" << i << "]: " << vdFdX[i]->dep() << std::endl;

  // 2nd-order derivative DAG construct
  std::vector<const mc::FFVar*> vd2FdX2 = FF.FAD( vdFdX, vX );
  std::cout << FF;

  std::vector<const mc::FFVar*> l_F_dFdX_d2Fdx2;
  for( unsigned int i=0, ij=0, ijk=0; i<NF; i++ ){
    l_F_dFdX_d2Fdx2.push_back( vF[i] );
    for( unsigned int j=0; j<NX; j++, ij++ ){
      l_F_dFdX_d2Fdx2.push_back( vdFdX[ij] );
      for( unsigned int k=0; k<NX; k++, ijk++ )
        l_F_dFdX_d2Fdx2.push_back( vd2FdX2[ijk] );
    }
  }
  std::list<const mc::FFOp*> F_dFdX_d2Fdx2_op = FF.subgraph( l_F_dFdX_d2Fdx2 );
  FF.output( F_dFdX_d2Fdx2_op );

  std::ofstream o_F_dFdX_d2Fdx2( "F_dFdX_d2Fdx2.dot", std::ios_base::out );
  FF.dot_script( l_F_dFdX_d2Fdx2, o_F_dFdX_d2Fdx2 );
  o_F_dFdX_d2Fdx2.close();

  for( unsigned i=0; i<NF*NX*NX; i++ )
    std::cout << "Variable dependence of d2FdX2[" << i << "]: " << vd2FdX2[i]->dep() << std::endl;

/*
  // Combination with fadbad::F for DAG of derivatives
  // FF.clear();
  fadbad::F<mc::FFVar> F_X[NX];
  for( int i=0; i<NX; i++ ){
    F_X[i].x() = X[i].set( &FF );
    F_X[i].diff( i, NX );
  }
  fadbad::F<mc::FFVar> F_F[NF]
    = { F_X[2]*F_X[3]-F_X[0],
        F_X[0]*pow(exp(F_X[2]*F_X[3]*2)+3.,3)+F_X[1] };
  std::cout << FF;
  std::cout << "Variable dependence of F[0]: " << F_F[0].x().dep() << std::endl;
  std::cout << "Variable dependence of F[1]: " << F_F[1].x().dep() << std::endl;
  for( int i=0; i<NX; i++ ){
    std::cout << "Variable dependence of dFdX[0," << i << "]: " << F_F[0].d(i).dep() << std::endl;
    std::cout << "Variable dependence of dFdX[1," << i << "]: " << F_F[1].d(i).dep() << std::endl;
  }

  std::vector<const mc::FFVar*> l_F;
  for( int j=0; j<NF; j++ ) l_F.push_back( &F_F[j].x() );
  std::list<const mc::FFOp*> F_op = FF.subgraph( l_F );
  FF.output( F_op );

  std::ofstream o_F( "F.dot", std::ios_base::out );
  FF.dot_script( l_F, o_F );
  o_F.close();

  std::vector<const mc::FFVar*> l_dFdX;
  for( int i=0; i<NX; i++ ){
    for( int j=0; j<NF; j++ ) l_dFdX.push_back( &F_F[j].d(i) );
  }
  std::list<const mc::FFOp*> dFdX_op = FF.subgraph( l_dFdX );
  FF.output( dFdX_op );

  std::ofstream o_dFdX( "dFdX.dot", std::ios_base::out );
  FF.dot_script( l_dFdX, o_dFdX );
  o_dFdX.close();

  std::vector<const mc::FFVar*> l_F_dFdX;
  for( int j=0; j<NF; j++ ) l_F_dFdX.push_back( &F_F[j].x() );
  for( int i=0; i<NX; i++ ){
    for( int j=0; j<NF; j++ ) l_F_dFdX.push_back( &F_F[j].d(i) );
  }
  std::list<const mc::FFOp*> F_dFdX_op = FF.subgraph( l_F_dFdX );
  FF.output( F_dFdX_op );

  std::ofstream o_F_dFdX( "F_dFdX.dot", std::ios_base::out );
  FF.dot_script( l_F_dFdX, o_F_dFdX );
  o_F_dFdX.close();

  // Combination with fadbad::B for DAG of derivatives
  // FF.clear();
  fadbad::B<mc::FFVar> B_X[NX];
  for( unsigned i=0; i<NX; i++ ){
    B_X[i] = X[i].set( &FF );
  }
  fadbad::B<mc::FFVar> Z0 = sqr(B_X[0]);
  fadbad::B<mc::FFVar> B_F[NF]
    = { B_X[1]*Z0 };
    //= { B_X[1]*B_X[2],
    //    2.*B_X[0] };
    //= { B_X[2]*B_X[3]-B_X[0],
    //    B_X[0]*pow(exp(B_X[2]*B_X[3]*2)+3.,4)+B_X[1] };

  //return 0;
  //std::cout << "Z0.m_rc: " << Z0.getBTypeNameHV()->m_rc;
  //std::cout << "F0.m_rc: " << B_F[0].getBTypeNameHV()->m_rc;

  for( unsigned j=0; j<NF; j++ )
    B_F[j].diff( j, NF );

  std::cout << FF;
  for( unsigned j=0; j<NF; j++ ){
    std::cout << "Variable dependence of F[" << j << "]: " << B_F[j].x().dep() << std::endl;
    for( unsigned i=0; i<NX; i++ ){
      std::cout << "Variable dependence of dFdX[" << j << "," << i << "]: " << B_X[i].d(j).dep() << std::endl;
    }
  }

  std::vector<const mc::FFVar*> l_F;
  for( unsigned j=0; j<NF; j++ ) l_F.push_back( &B_F[j].x() );
  std::list<const mc::FFOp*> F_op = FF.subgraph( l_F );
  FF.output( F_op );

  std::ofstream o_F( "F_2.dot", std::ios_base::out );
  FF.dot_script( l_F, o_F );
  o_F.close();

  std::vector<const mc::FFVar*> l_dFdX;
  for( unsigned i=0; i<NX; i++ ){
    for( unsigned j=0; j<NF; j++ ) l_dFdX.push_back( &B_X[i].d(j) );
  }
  std::list<const mc::FFOp*> dFdX_op = FF.subgraph( l_dFdX );
  FF.output( dFdX_op );

  std::ofstream o_dFdX( "dFdX_2.dot", std::ios_base::out );
  FF.dot_script( l_dFdX, o_dFdX );
  o_dFdX.close();

  std::vector<const mc::FFVar*> l_F_dFdX;
  for( unsigned j=0; j<NF; j++ ) l_F_dFdX.push_back( &B_F[j].x() );
  for( unsigned i=0; i<NX; i++ ){
    for( unsigned j=0; j<NF; j++ ) l_F_dFdX.push_back( &B_X[i].d(j) );
  }
  std::list<const mc::FFOp*> F_dFdX_op = FF.subgraph( l_F_dFdX );
  FF.output( F_dFdX_op );

  std::ofstream o_F_dFdX( "F_dFdX_2.dot", std::ios_base::out );
  FF.dot_script( l_F_dFdX, o_F_dFdX );
  o_F_dFdX.close();
*/
  return 0;
}

int test_DAG2()
{
  // Create DAG for f1,f2
  const unsigned int NX = 2, NF = 2;
  mc::FFGraph DAG;
  mc::FFVar X[NX];
  for( unsigned int i=0; i<NX; i++ ) X[i].set( &DAG );
  mc::FFVar F[NF] = { X[0]*exp(X[1]), pow(X[1],3)*sqrt(X[0]) };
  std::cout << DAG;

  // Display DAG for f1,f2
  std::list<const mc::FFOp*> F_op  = DAG.subgraph( NF, F );    DAG.output( F_op );
  std::ofstream o_F( "F.dot", std::ios_base::out );
  DAG.dot_script( NF, F, o_F );
  o_F.close();

  // Display DAG for f1 alone
  std::list<const mc::FFOp*> F0_op = DAG.subgraph( 1, F );     DAG.output( F0_op );
  std::ofstream o_F0( "F0.dot", std::ios_base::out );
  DAG.dot_script( 1, F, o_F0 );
  o_F0.close();

  // Evaluate DAG for f1,f2 in double arithmetic
  double F_d[NF], X_d[NX] = { 2., -1. };
  DAG.eval( F_op, NF, F, F_d, NX, X, X_d );
  std::cout <<"\nDouble Value:\n";
  for( unsigned int i=0; i<NX; i++ ) std::cout << "F[" << i << "] = " << F_d[i] << std::endl;

  // Evaluate DAG for f1,f2 in ellipsoidal arithmetic
  dcv X_c(NX);    dsm X_Q(NX);
  X_c(0) =  2.;   X_Q(0,0) = 0.1;
  X_c(1) = -1.;   X_Q(1,0) = 0.04; X_Q(1,1) = 0.1; 
  EI X_EI( X_Q, X_c, DAG.depmap( F_op, NF, F, NX, X ) );
  EV F_EV[NF], X_EV[NX];
  for( unsigned i=0; i<NX; ++i ) X_EV[i].set( X_EI, i );
  DAG.eval( F_op, NF, F, F_EV, NX, X, X_EV );
  EI F_EI = X_EI.get( NF, F_EV ) ;
  std::cout <<"\nEllipsoidal Image:" << F_EI << std::endl;

/*
  // Evaluate DAG for f1,f2 in Taylor model with ellipsoidal arithmetic
  X_EI.set( X_Q, X_c );
  for( unsigned i=0; i<NX; ++i ) X_EV[i].set( X_EI, i );
  std::cout << X_EV[0] << std::endl;
  TME TM_env( NX, 1 );
  TVE X_TVE[NX], F_TVE[NF];
  for( unsigned int i=0; i<NX; i++ )
     X_TVE[i].set( &TM_env, i, X_EV[i] );
  //X_EI.output();
  //std::cout << X_TVE[0];
  DAG.eval( F_op, NF, F, F_TVE, NX, X, X_TVE );
  EI RF_EI = X_EI.get( NF, F_EV ) ;
  std::cout <<"\nEllipsoidal Image:" << F_EI << std::endl;
*/

  return 0;
}

int test_comp()
{
  mc::FFGraph DAG;
  mc::FFVar X, Y, F, G;
  X.set( &DAG );
  Y.set( &DAG );
  F = exp(X);
  G = sqr(Y)+F;
  std::cout << DAG;

  std::ofstream o_comp0( "compose0.dot", std::ios_base::out );
  DAG.dot_script( 1, &G, o_comp0 );
  o_comp0.close();

  const mc::FFVar* GoF = DAG.compose( 1, &G, 1, &Y, &F );
  std::cout << DAG;

  const mc::FFVar* dGoFdX = DAG.FAD( 1, GoF, 1, &X );
  std::cout << DAG;

  std::ofstream o_comp1( "compose1.dot", std::ios_base::out );
  DAG.dot_script( 1, GoF, o_comp1 );
  o_comp1.close();

  std::ofstream o_comp2( "compose2.dot", std::ios_base::out );
  DAG.dot_script( 1, dGoFdX, o_comp2 );
  o_comp2.close();

  delete[] GoF;
  delete[] dGoFdX;
  return 0;
}

int test_sparseder()
{
  mc::FFGraph DAG;
  const unsigned NX = 5, NF = 2;
  mc::FFVar X[NX];
  for( unsigned i(0); i<NX; i++ )  X[i].set( &DAG );
  mc::FFVar F[NF] = { sqrt(X[0])*X[1], X[3]+2. };
  std::cout << DAG;

  auto dFdX = DAG.SBAD( NF, F, NX, X );
  //std::cout << DAG;
  std::cout << "\nNumber of non-zero elements in Jacobian matrix: "
            << std::get<0>(dFdX) << std::endl;
  for( unsigned ie=0; ie<std::get<0>(dFdX); ie++ )
    std::cout << "(" << std::get<1>(dFdX)[ie] << "," << std::get<2>(dFdX)[ie] << ") "
              << std::get<3>(dFdX)[ie] << std::endl;

  DAG.output( DAG.subgraph( std::get<0>(dFdX), std::get<3>(dFdX) ) );

  std::ofstream o_dFdX( "sparseFAD.dot", std::ios_base::out );
  DAG.dot_script( std::get<0>(dFdX), std::get<3>(dFdX), o_dFdX );
  o_dFdX.close();

  auto d2FdX2 = DAG.SFAD( std::get<0>(dFdX), std::get<3>(dFdX), NX, X );
  //std::cout << DAG;
  std::cout << "\nNumber of non-zero elements in Hessian matrix: "
            << std::get<0>(d2FdX2) << std::endl;
  for( unsigned ie=0; ie<std::get<0>(d2FdX2); ie++ )
    std::cout << "(" << std::get<1>(d2FdX2)[ie] << "," << std::get<2>(d2FdX2)[ie] << ") "
              << std::get<3>(d2FdX2)[ie] << std::endl;

  DAG.output( DAG.subgraph( std::get<0>(d2FdX2), std::get<3>(d2FdX2) ) );

  std::ofstream o_d2FdX2( "sparseFAD2.dot", std::ios_base::out );
  DAG.dot_script( std::get<0>(d2FdX2), std::get<3>(d2FdX2), o_d2FdX2 );
  o_d2FdX2.close();

  delete[] std::get<1>(dFdX);
  delete[] std::get<2>(dFdX);
  delete[] std::get<3>(dFdX);

  delete[] std::get<1>(d2FdX2);
  delete[] std::get<2>(d2FdX2);
  delete[] std::get<3>(d2FdX2);
  return 0;
}

int test_dirder()
{
  mc::FFGraph DAG;
  mc::FFVar X[2], F, D[2];
  X[0].set( &DAG );
  X[1].set( &DAG );
  F = sqrt(X[0])*X[1];
  D[0].set( &DAG, 1. );
  D[1].set( &DAG, 1. );
  //D[0] = 1.;
  //D[1] = 1.;
  std::cout << DAG;

  const mc::FFVar* dFdXxD = DAG.FAD( 1, &F, 2, X, D );
  //std::cout << DAG;
  DAG.output( DAG.subgraph( 1, dFdXxD ) );

  std::ofstream o_dFdXxD( "directional.dot", std::ios_base::out );
  DAG.dot_script( 1, dFdXxD, o_dFdXxD );
  o_dFdXxD.close();

  delete[] dFdXxD;
  return 0;
}

int test_const()
{
  mc::FFGraph DAG;
  std::vector<mc::FFVar> X(2);
  mc::FFVar F, PI;
  X[0].set( &DAG );
  X[1].set( &DAG );
  PI.set( &DAG, mc::PI );
  F = exp(X[0])*X[1]/PI;
  std::cout << DAG;

  std::ofstream o_cst1( "constant.dot", std::ios_base::out );
  DAG.dot_script( 1, &F, o_cst1 );
  o_cst1.close();

  // Evaluate DAG for F in interval arithmetic
  I X_I[2] = { I(-1.,1.), I(0.,2.) }, F_I;
  DAG.eval( 1, &F, &F_I, 2, X.data(), X_I );
  std::cout <<"\nInterval bound:\n";
  std::cout << "F = " << F_I << std::endl;

  X[1].set(1);
  DAG.eval( 1, &F, &F_I, 1, X.data(), X_I );
  std::cout <<"\nInterval bound:\n";
  std::cout << "F = " << F_I << std::endl;

  //X[1].unset();
  DAG.eval( 1, &F, &F_I, 2, X.data(), X_I );
  std::cout <<"\nInterval bound:\n";
  std::cout << "F = " << F_I << std::endl;

  const mc::FFVar* dFdX = DAG.FAD( 1, &F, 1, X.data() );
  std::cout << DAG;

  std::ofstream o_cst2( "compose2.dot", std::ios_base::out );
  DAG.dot_script( 1, dFdX, o_cst2 );
  o_cst2.close();

  delete[] dFdX;
  return 0;
}

int test_badiff0()
{
  const unsigned NX = 1, NF = 1;
  double X[NX];
  for( unsigned i(0); i<NX; i++ )  X[i] = 2.;
  //double F[NF] = { X[0] };

  fadbad::B<double> X_B = X[0];
  fadbad::B<double> F_B = X_B;
  F_B.diff(0,1);
  //F_B = 0.;
  double dFdX = X_B.d(0);
  std::cout << dFdX << std::endl;
  return 0;
}

int test_badiff()
{
  mc::FFGraph DAG;
  const unsigned NX = 2, NF = 2;
  mc::FFVar X[NX];
  for( unsigned i(0); i<NX; i++ )  X[i].set( &DAG );
  mc::FFVar F[NF] = { X[0], 1. };//exp(X[1])*X[0] };
  std::cout << DAG;

  fadbad::B<mc::FFVar> X0_B = X[0];
  fadbad::B<mc::FFVar> X1_B = X[1];
  fadbad::B<mc::FFVar> F0_B = X0_B;
  fadbad::B<mc::FFVar> F1_B = 1.;//X1_B;//exp(X1_B)*X0_B;
  F0_B.diff(0,2);
  F1_B.diff(1,2);
  //F_B = 0.;
  X0_B.d(0);  
  X1_B.d(0);
  std::cout << X0_B.d(0) << std::endl;
  std::cout << X0_B.d(1) << std::endl;
  std::cout << X1_B.d(0) << std::endl;
  std::cout << X1_B.d(1) << std::endl;
  //return 0;

  fadbad::B<mc::FFVar> X_B[NX], F_B[NF];
  X_B[0] = X[0]; X_B[1] = X[1];
  DAG.eval( NF, F, F_B, NX, X, X_B );
  F_B[0].diff(0,NF);
  F_B[1].diff(1,NF);
  std::cout << X_B[0].d(0) << std::endl;
  std::cout << X_B[0].d(1) << std::endl;
  std::cout << X_B[1].d(0) << std::endl;
  std::cout << X_B[1].d(1) << std::endl;
  //return 0;
/*
  std::vector< fadbad::B<mc::FFVar> > X_B(NX);
  X_B[0] = X[0]; X_B[1] = X[1];
  std::vector< std::pair< const mc::FFVar*,fadbad::B<mc::FFVar> > > vX_B;
  vX_B.push_back( std::make_pair( &X[0], X_B[0] ) );
  vX_B.push_back( std::make_pair( &X[1], X_B[1] ) );
  std::vector< const mc::FFVar* > vF;
  vF.push_back( &F[0] );
  vF.push_back( &F[1] );
  std::vector< fadbad::B<mc::FFVar> > vF_B = DAG.eval( vF, vX_B );
  vF_B[0].diff(0,NF);
  vF_B[1].diff(1,NF);
  std::cout << vX_B[0].second.d(0) << std::endl;
  std::cout << vX_B[0].second.d(1) << std::endl;
  std::cout << vX_B[1].second.d(0) << std::endl;
  std::cout << vX_B[1].second.d(1) << std::endl;
  //return 0;
*/
  auto dFdX = DAG.SBAD( NF, F, NX, X );
  std::cout << DAG;
  std::cout << "\nNumber of non-zero elements in Jacobian matrix: "
            << std::get<0>(dFdX) << std::endl;
  for( unsigned ie=0; ie<std::get<0>(dFdX); ie++ )
    std::cout << "(" << std::get<1>(dFdX)[ie] << "," << std::get<2>(dFdX)[ie] << ") "
              << std::get<3>(dFdX)[ie] << std::endl;

  DAG.output( DAG.subgraph( std::get<0>(dFdX), std::get<3>(dFdX) ) );

  std::ofstream o_dFdX( "sparseBAD.dot", std::ios_base::out );
  DAG.dot_script( std::get<0>(dFdX), std::get<3>(dFdX), o_dFdX );
  o_dFdX.close();

  delete[] std::get<1>(dFdX);
  delete[] std::get<2>(dFdX);
  delete[] std::get<3>(dFdX);

  return 0;
}

int main()
{
  try{
    //test_DAG();
    //test_DAG2();
    //test_TAD(10);
    //test_comp();
    //test_dirder();
    test_sparseder();
    test_badiff();
    //test_eval();

    //for( unsigned i=0; i<10; i++ ) AD_F_ODE();
    //for( unsigned i=0; i<10; i++ ) test_F_ODE( true );
    //AD_TE_ODE( NTE );
    //test_TE_ODE( NTE, true );
    //for( unsigned q=0; q<=NTE; q++ )
    //  test_TE_ODE( q, true );
    //for( unsigned q=0; q<=NTE; q++ )
    //  AD_TE_ODE( q );
  }

#ifndef USE_PROFIL
#ifndef USE_FILIB
  catch( I::Exceptions &eObj ){
    std::cerr << "Error " << eObj.ierr()
              << " in natural interval extension:" << std::endl
              << eObj.what() << std::endl
              << "Aborts." << std::endl;
    return eObj.ierr();
  }
#endif
#endif
  catch( TM::Exceptions &eObj ){
    std::cerr << "Error " << eObj.ierr()
              << " in Taylor model computation:" << std::endl
              << eObj.what() << std::endl
              << "Aborts." << std::endl;
    return eObj.ierr();
  }
  catch( mc::FFGraph::Exceptions &eObj ){
    std::cerr << "Error " << eObj.ierr()
              << " in factorable function manipulation:" << std::endl
              << eObj.what() << std::endl
              << "Aborts." << std::endl;
    return eObj.ierr();
  }
}

