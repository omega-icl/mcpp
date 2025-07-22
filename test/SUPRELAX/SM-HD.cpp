#define TEST_SHEKEL	// <-- select test function here
#define SAVE_RESULTS    // <-- specify whether to save results to file
#undef  ANALYSE_RATE    // <-- specify whether to analyse rate of convergence
#undef  ANALYSE_TIME    // <-- specify whether to analyse computational time
////////////////////////////////////////////////////////////////////////

#include <fstream>
#include <iomanip>
#include <chrono>

#include <boost/random/sobol.hpp>
#include <boost/random/uniform_01.hpp>
#include <boost/random/variate_generator.hpp>

#include "mcboost.hpp"
typedef boost::numeric::interval_lib::save_state<boost::numeric::interval_lib::rounded_transc_opp<double>> T_boost_round;
typedef boost::numeric::interval_lib::checking_base<double> T_boost_check;
typedef boost::numeric::interval_lib::policies<T_boost_round,T_boost_check> T_boost_policy;
typedef boost::numeric::interval<double,T_boost_policy> I;

#include "mccormick.hpp"
typedef mc::McCormick<I> MC;

#include "supmodel.hpp"
#include "pwcu.hpp"
typedef mc::SupModel<mc::PWCU> PWCSM;
typedef mc::SupVar<mc::PWCU> PWCSV;
#include "pwlu.hpp"
typedef mc::SupModel<mc::PWLU> PWLSM;
typedef mc::SupVar<mc::PWLU> PWLSV;

using namespace std;
using namespace mc;

////////////////////////////////////////////////////////////////////////

#if defined( TEST_EXP )
unsigned const NX = 4;
const double XL[NX] = { -2, -1, -2, -1 };	// <-- X range lower bound
const double XU[NX] = {  1,  2,  1,  2 };	// <-- X range upper bound
const double XREF[NX] = { -0.5, 0.5, -0.5, 0.5 };		// <-- X reference point

template <class T>
T myfunc
( T const* x )
{
  T w = x[0],
    z = x[0]*(exp(x[0])-exp(-x[0]));
  for( unsigned i=1; i<NX; ++i ){
    w *= x[i];
    if( i%2 ) z -= x[i]*(exp(x[i])-exp(-x[i]));
    else      z += x[i]*(exp(x[i])-exp(-x[i]));
  }
  return w*z;
}

#elif defined( TEST_SHEKEL )
unsigned const NX = 6;
unsigned const M  = 10;

const double B[10]     = {0.1, 0.2, 0.2, 0.4, 0.4, 0.6, 0.3, 0.7, 0.5, 0.5};
const double C[NX][10] = {{4.0, 1.0, 8.0, 6.0, 3.0, 2.0, 5.0, 8.0, 6.0, 7.0},
                          {4.0, 1.0, 8.0, 6.0, 7.0, 9.0, 3.0, 1.0, 2.0, 3.6},
                          {4.0, 1.0, 8.0, 6.0, 3.0, 2.0, 5.0, 8.0, 6.0, 7.0},
                          {4.0, 1.0, 8.0, 6.0, 7.0, 9.0, 3.0, 1.0, 2.0, 3.6},
                          {4.0, 1.0, 8.0, 6.0, 3.0, 2.0, 5.0, 8.0, 6.0, 7.0},
                          {4.0, 1.0, 8.0, 6.0, 7.0, 9.0, 3.0, 1.0, 2.0, 3.6}};

const double XL[6] = {  0,  0,  0,  0,  0,  0 };	// <-- X range lower bound
const double XU[6] = { 10, 10, 10, 10, 10, 10 };	// <-- X range upper bound
const double XREF[6] = { 4, 4, 4, 4, 4, 4 };		// <-- X reference point

template <class T>
T myfunc
( T const* x )
{
  T z = 0.;
  for( unsigned i=0; i<M; ++i ){
    T w = 0.;
    for( unsigned j=0; j<NX; ++j )
      w += pow(x[j]-C[j][i],2);
    z += 1./(w+B[i]);
  }
  return z;
}

#elif defined( TEST_CCPP )
#include "ffexpr.hpp"
#include "gamsio.hpp"

unsigned const NX = 5;
const double XL[5] = { 0.01, 0.3,   5, 2480, 0.01 };	// <-- X range lower bound
const double XU[5] = { 0.5,  10., 100, 3750, 0.2  };	// <-- X range upper bound
const double XREF[5] = { 0.0200, 4.5263, 25.3575, 3643.3781, 0.0328 };		// <-- X reference point

template <class T>
T myfunc
( T const* x )
{
  static FFGraph* DAG = nullptr;
  static FFSubgraph GF;
  static std::vector<FFVar> X, F;

  if( !DAG ){
    static GAMSIO model;
    model.read( "CS2_LCOE.gms", 0, 0 );
    //model.read( "CS2_WNET.gms", 0, 0 );
    assert( model.var().size() == NX );
    DAG = model.dag();
    X   = model.var();
    F   = std::get<1>( model.obj() );
    GF  = DAG->subgraph( F );
    
    vector<mc::FFExpr> EX(NX), EF;
    for( unsigned int i=0; i<NX; i++ ) EX[i].set( X[i] );
    DAG->eval( GF, F, EF, X, EX );
    std::cout << "Expression of F: " << EF[0] << std::endl;
  }

  static std::vector<T> wkF;
  static T f;
  DAG->eval( GF, wkF, 1, F.data(), &f, NX, X.data(), x ); 
  return f;
}

#endif

////////////////////////////////////////////////////////////////////////
int main()
////////////////////////////////////////////////////////////////////////
{
/*
  GAMSIO model;
  model.read( "CS2_WNET.gms", 0, 1 );

  vector<mc::FFExpr> EX(model.var().size()), EF;
  for( unsigned int i=0; i<model.var().size(); i++ ) EX[i].set( model.var()[i] );
  model.dag()->eval( std::get<1>(model.obj()), EF, model.var(), EX );
  std::cout << "Expression of F: " << EF[0] << std::endl;

  vector<double> DX(model.var().size()), DF;
  DX = { 0.0200, 4.5263, 25.3575, 3643.3781, 0.0328 };
  model.dag()->eval( std::get<1>(model.obj()), DF, model.var(), DX );
  std::cout << "Value of F: " << DF[0] << std::endl;
  
  return 0;
*/
  try{
    MC::options.ENVEL_USE   = 1;
    MC::options.ENVEL_MAXIT = 100;
    MC::options.ENVEL_TOL   = 1e-12;
    MC::options.MVCOMP_USE  = 1;

    PWLSM pwlmod( NX );
    pwlmod.options.PROD_METH      = PWLSM::Options::PARTIAL;//FULL;//LOG;//NONE;
    pwlmod.options.PROD_CUT       = 0;
    pwlmod.options.REF_WEIGHT     = 0.5;
    pwlmod.options.MAX_SUBDIV     = 0;//8;//16;
    pwlmod.options.USE_SHADOW     = 0;
    pwlmod.options.DISPLAY_SHADOW = 1;
    pwlmod.options.DISPLAY_DIGITS = 10;

    PWCSM pwcmod( NX );
    pwcmod.options = pwlmod.options;

    // Calculate superposition relaxations
    vector<PWLSV> PWLSVX( NX ), PWLSVF( 1 );
    int const NPWL = 16;	// <-- select initial variable partition >=1
    for( unsigned i=0; i<NX; ++i ) PWLSVX[i].set( pwlmod, i, I(XL[i],XU[i]), NPWL );
    vector<double> DF( 1 );
    PWLSVF[0] = myfunc( PWLSVX.data() );
    cout << PWLSVF[0];

    vector<I> XBND0(NX);
    for( unsigned i=0; i<NX; ++i ) XBND0[i] = I(XL[i],XU[i]);
    
#ifdef ANALYSE_RATE
    ofstream ratefile( "SM-2D_rate.out", ios_base::out );
    ratefile << scientific << setprecision(5) << right;

    auto const& red = [=]( const I& bnd, const double& ref, const double& r ){ return r*bnd + (1-r)*ref; };
    auto const& min = [=]( const double& x, const double& y ){ return x<y?x:y; };
    auto const& max = [=]( const double& x, const double& y ){ return x>y?x:y; };

    size_t const NGRID  = pow(2,22);
    double const rate   = 0.95;
    double const rhomin = 1e-5;
    DF[0] = myfunc( XREF );//.data() );
    for( double rho = 1; rho > rhomin; rho *= rate ){
      double min_F = DF[0],           max_F = DF[0];
      double min_MCFcv = DF[0],       max_MCFcc = DF[0];
      double distmax_MCFcv = 0.,      distmax_MCFcc = 0.;
      double distmax_PWC16SVFu = 0.,  distmax_PWC16SVFo = 0.;
      double distmax_PWC32SVFu = 0.,  distmax_PWC32SVFo = 0.;
      double distmax_PWC64SVFu = 0.,  distmax_PWC64SVFo = 0.;
      double distmax_PWC128SVFu = 0., distmax_PWC128SVFo = 0.;
      double distmax_PWL1SVFu = 0.,   distmax_PWL1SVFo = 0.;
      double distmax_PWL2SVFu = 0.,   distmax_PWL2SVFo = 0.;
      double distmax_PWL4SVFu = 0.,   distmax_PWL4SVFo = 0.;
      double distmax_PWL8SVFu = 0.,   distmax_PWL8SVFo = 0.;
      double distmax_PWL16SVFu = 0.,  distmax_PWL16SVFo = 0.;
      double min_PWC16SVFu = DF[0],   max_PWC16SVFo = DF[0];
      double min_PWC32SVFu = DF[0],   max_PWC32SVFo = DF[0];
      double min_PWC64SVFu = DF[0],   max_PWC64SVFo = DF[0];
      double min_PWC128SVFu = DF[0],  max_PWC128SVFo = DF[0];
      double min_PWL1SVFu = DF[0],    max_PWL1SVFo = DF[0];
      double min_PWL2SVFu = DF[0],    max_PWL2SVFo = DF[0];
      double min_PWL4SVFu = DF[0],    max_PWL4SVFo = DF[0];
      double min_PWL8SVFu = DF[0],    max_PWL8SVFo = DF[0];
      double min_PWL16SVFu = DF[0],   max_PWL16SVFo = DF[0];

      vector<I> XBND(NX);
      cout << rho;
      for( unsigned i=0; i<NX; ++i ){
        XBND[i] = red( XBND0[i], XREF[i], rho );
        cout << " " << XBND[i];
      }
      cout << std::endl;

      vector<PWCSV> PWC16SVX(NX), PWC16SVF( 1 );
      for( unsigned i=0; i<NX; ++i )
        PWC16SVX[i].set( pwcmod, i, XBND[i], 16 );
      try{
        PWC16SVF[0] = myfunc( PWC16SVX.data() );
      }
      catch(...){
        PWC16SVF[0] = 0./0.;      
      }
      //cout << " PWC16SVX ";

      vector<PWCSV> PWC32SVX(NX), PWC32SVF( 1 );
      for( unsigned i=0; i<NX; ++i )
        PWC32SVX[i].set( pwcmod, i, XBND[i], 32 );
      try{
        PWC32SVF[0] = myfunc( PWC32SVX.data() );
      }
      catch(...){
        PWC32SVF[0] = 0./0.;      
      }

      vector<PWCSV> PWC64SVX(NX), PWC64SVF( 1 );
      for( unsigned i=0; i<NX; ++i )
        PWC64SVX[i].set( pwcmod, i, XBND[i], 64 );
      try{
        PWC64SVF[0] = myfunc( PWC64SVX.data() );
      }
      catch(...){
        PWC64SVF[0] = 0./0.;      
      }

      vector<PWCSV> PWC128SVX(NX), PWC128SVF( 1 );
      for( unsigned i=0; i<NX; ++i )
        PWC128SVX[i].set( pwcmod, i, XBND[i], 128 );
      try{
        PWC128SVF[0] = myfunc( PWC128SVX.data() );
      }
      catch(...){
        PWC128SVF[0] = 0./0.;      
      }

      vector<PWLSV> PWL1SVX(NX), PWL1SVF( 1 );
      for( unsigned i=0; i<NX; ++i )
        PWL1SVX[i].set( pwlmod, i, XBND[i], 1 );
      try{
        PWL1SVF[0] = myfunc( PWL1SVX.data() );
      }
      catch(...){
        PWL1SVF[0] = 0./0.;      
      }
      //cout << " PWL1SVX ";

      vector<PWLSV> PWL2SVX(NX), PWL2SVF( 1 );
      for( unsigned i=0; i<NX; ++i )
        PWL2SVX[i].set( pwlmod, i, XBND[i], 2 );
      PWL2SVX[1].set( pwlmod, 1, XBND[1], 2 );
      try{
        PWL2SVF[0] = myfunc( PWL2SVX.data() );
      }
      catch(...){
        PWL2SVF[0] = 0./0.;      
      }
      //cout << " PWL2SVX ";

      vector<PWLSV> PWL4SVX(NX), PWL4SVF( 1 );
      for( unsigned i=0; i<NX; ++i )
        PWL4SVX[i].set( pwlmod, i, XBND[i], 4 );
      try{
        PWL4SVF[0] = myfunc( PWL4SVX.data() );
      }
      catch(...){
        PWL4SVF[0] = 0./0.;      
      }
      //cout << " PWL4SVX ";

      vector<PWLSV> PWL8SVX(NX), PWL8SVF( 1 );
      for( unsigned i=0; i<NX; ++i )
        PWL8SVX[i].set( pwlmod, i, XBND[i], 8 );
      try{
        PWL8SVF[0] = myfunc( PWL8SVX.data() );
      }
      catch(...){
        PWL8SVF[0] = 0./0.;      
      }
      //cout << " PWL8SVX " << endl;

      vector<PWLSV> PWL16SVX(NX), PWL16SVF( 1 );
      for( unsigned i=0; i<NX; ++i )
        PWL16SVX[i].set( pwlmod, i, XBND[i], 16 );
      try{
        PWL16SVF[0] = myfunc( PWL16SVX.data() );
      }
      catch(...){
        PWL16SVF[0] = 0./0.;      
      }
      //cout << " PWL8SVX " << endl;

      typedef boost::random::sobol_engine< boost::uint_least64_t, 64u > sobol64;
      typedef boost::variate_generator< sobol64, boost::uniform_01< double > > qrgen;
      sobol64 eng( NX );
      qrgen gen( eng, boost::uniform_01<double>() );
      gen.engine().seed( 0 );

      vector<double> DX(NX);
      map<unsigned,double> mDX;
      for( unsigned s=0; s<NGRID; s++ ){
        for( unsigned i=0; i<NX; ++i )
          mDX[i] = DX[i] = Op<I>::l(XBND[i]) + Op<I>::diam(XBND[i]) * gen();
        vector<double> DF( 1 );
        DF[0] = myfunc( DX.data() );
        min_F = min( min_F, DF[0] );
        max_F = max( max_F, DF[0] );

        vector<MC> MCX(NX), MCF( 1 );
        for( unsigned i=0; i<NX; ++i ) MCX[i].I( XBND[i] ).c( DX[i] );
        MCF[0] = myfunc( MCX.data() );
        min_MCFcv = min( min_MCFcv, MCF[0].cv() );
        max_MCFcc = max( max_MCFcc, MCF[0].cc() );
        distmax_MCFcv = max( distmax_MCFcv, DF[0] - MCF[0].cv() );
        distmax_MCFcc = max( distmax_MCFcc, MCF[0].cc() - DF[0] );

        double const PWC16SVFu = PWC16SVF[0].uval( mDX );
        double const PWC16SVFo = PWC16SVF[0].oval( mDX );
        distmax_PWC16SVFu = max( distmax_PWC16SVFu, DF[0] - PWC16SVFu );
        distmax_PWC16SVFo = max( distmax_PWC16SVFo, PWC16SVFo - DF[0] );
        min_PWC16SVFu = min( min_PWC16SVFu, PWC16SVFu );
        max_PWC16SVFo = max( max_PWC16SVFo, PWC16SVFo );

        double const PWC32SVFu = PWC32SVF[0].uval( mDX );
        double const PWC32SVFo = PWC32SVF[0].oval( mDX );
        distmax_PWC32SVFu = max( distmax_PWC32SVFu, DF[0] - PWC32SVFu );
        distmax_PWC32SVFo = max( distmax_PWC32SVFo, PWC32SVFo - DF[0] );
        min_PWC32SVFu = min( min_PWC32SVFu, PWC32SVFu );
        max_PWC32SVFo = max( max_PWC32SVFo, PWC32SVFo );

        double const PWC64SVFu = PWC64SVF[0].uval( mDX );
        double const PWC64SVFo = PWC64SVF[0].oval( mDX );
        distmax_PWC64SVFu = max( distmax_PWC64SVFu, DF[0] - PWC64SVFu );
        distmax_PWC64SVFo = max( distmax_PWC64SVFo, PWC64SVFo - DF[0] );
        min_PWC64SVFu = min( min_PWC64SVFu, PWC64SVFu );
        max_PWC64SVFo = max( max_PWC64SVFo, PWC64SVFo );

        double const PWC128SVFu = PWC128SVF[0].uval( mDX );
        double const PWC128SVFo = PWC128SVF[0].oval( mDX );
        distmax_PWC128SVFu = max( distmax_PWC128SVFu, DF[0] - PWC128SVFu );
        distmax_PWC128SVFo = max( distmax_PWC128SVFo, PWC128SVFo - DF[0] );
        min_PWC128SVFu = min( min_PWC128SVFu, PWC128SVFu );
        max_PWC128SVFo = max( max_PWC128SVFo, PWC128SVFo );

        double const PWL1SVFu = PWL1SVF[0].uval( mDX );
        double const PWL1SVFo = PWL1SVF[0].oval( mDX );
        distmax_PWL1SVFu = max( distmax_PWL1SVFu, DF[0] - PWL1SVFu );
        distmax_PWL1SVFo = max( distmax_PWL1SVFo, PWL1SVFo - DF[0] );
        min_PWL1SVFu = min( min_PWL1SVFu, PWL1SVFu );
        max_PWL1SVFo = max( max_PWL1SVFo, PWL1SVFo );

        double const PWL2SVFu = PWL2SVF[0].uval( mDX );
        double const PWL2SVFo = PWL2SVF[0].oval( mDX );
        distmax_PWL2SVFu = max( distmax_PWL2SVFu, DF[0] - PWL2SVFu );
        distmax_PWL2SVFo = max( distmax_PWL2SVFo, PWL2SVFo - DF[0] );
        min_PWL2SVFu = min( min_PWL2SVFu, PWL2SVFu );
        max_PWL2SVFo = max( max_PWL2SVFo, PWL2SVFo );

        double const PWL4SVFu = PWL4SVF[0].uval( mDX );
        double const PWL4SVFo = PWL4SVF[0].oval( mDX );
        distmax_PWL4SVFu = max( distmax_PWL4SVFu, DF[0] - PWL4SVFu );
        distmax_PWL4SVFo = max( distmax_PWL4SVFo, PWL4SVFo - DF[0] );
        min_PWL4SVFu = min( min_PWL4SVFu, PWL4SVFu );
        max_PWL4SVFo = max( max_PWL4SVFo, PWL4SVFo );

        double const PWL8SVFu = PWL8SVF[0].uval( mDX );
        double const PWL8SVFo = PWL8SVF[0].oval( mDX );
        distmax_PWL8SVFu = max( distmax_PWL8SVFu, DF[0] - PWL8SVFu );
        distmax_PWL8SVFo = max( distmax_PWL8SVFo, PWL8SVFo - DF[0] );
        min_PWL8SVFu = min( min_PWL8SVFu, PWL8SVFu );
        max_PWL8SVFo = max( max_PWL8SVFo, PWL8SVFo );

        double const PWL16SVFu = PWL16SVF[0].uval( mDX );
        double const PWL16SVFo = PWL16SVF[0].oval( mDX );
        distmax_PWL16SVFu = max( distmax_PWL16SVFu, DF[0] - PWL16SVFu );
        distmax_PWL16SVFo = max( distmax_PWL16SVFo, PWL16SVFo - DF[0] );
        min_PWL16SVFu = min( min_PWL16SVFu, PWL16SVFu );
        max_PWL16SVFo = max( max_PWL16SVFo, PWL16SVFo );
      }
      
      std::cout << "  minF: " << min_F << "  maxF: " << max_F
                << "  minPWL8SVFu: " << min_PWL8SVFu << " =?= " << PWL8SVF[0].l()
                << "  maxPWL8SVFo: " << max_PWL8SVFo << " =?= " << PWL8SVF[0].u()
                << std::endl;
      
      ratefile << scientific << setprecision(5) << right
               << setw(14) << rho;
      for( unsigned i=0; i<NX; ++i )
        ratefile << setw(14) << mc::Op<I>::l(XBND[i]) << setw(14) << mc::Op<I>::u(XBND[i]);
      ratefile << setw(14) << max( distmax_MCFcv,            distmax_MCFcc )
               << setw(14) << max( min_F - min_MCFcv,        max_MCFcc - max_F )
               << setw(14) << max( distmax_PWC16SVFu,        distmax_PWC16SVFo )
               << setw(14) << max( min_F - min_PWC16SVFu,    max_PWC16SVFo - max_F )
               << setw(14) << max( distmax_PWC32SVFu,        distmax_PWC32SVFo )
               << setw(14) << max( min_F - min_PWC32SVFu,    max_PWC32SVFo - max_F )
               << setw(14) << max( distmax_PWC64SVFu,        distmax_PWC64SVFo )
               << setw(14) << max( min_F - min_PWC64SVFu,    max_PWC64SVFo - max_F )
               << setw(14) << max( distmax_PWC128SVFu,       distmax_PWC128SVFo )
               << setw(14) << max( min_F - min_PWC128SVFu,   max_PWC128SVFo - max_F )
               << setw(14) << max( distmax_PWL1SVFu,         distmax_PWL1SVFo )
               << setw(14) << max( min_F - min_PWL1SVFu,     max_PWL1SVFo - max_F )
               << setw(14) << max( distmax_PWL2SVFu,         distmax_PWL2SVFo )
               << setw(14) << max( min_F - min_PWL2SVFu,     max_PWL2SVFo - max_F )
               << setw(14) << max( distmax_PWL4SVFu,         distmax_PWL4SVFo )
               << setw(14) << max( min_F - min_PWL4SVFu,     max_PWL4SVFo - max_F )
               << setw(14) << max( distmax_PWL8SVFu,         distmax_PWL8SVFo )
               << setw(14) << max( min_F - min_PWL8SVFu,     max_PWL8SVFo - max_F )
               << setw(14) << max( distmax_PWL16SVFu,         distmax_PWL16SVFo )
               << setw(14) << max( min_F - min_PWL16SVFu,     max_PWL16SVFo - max_F )
               << endl;
    }
#endif

#ifdef ANALYSE_TIME
    std::chrono::time_point<std::chrono::system_clock> start;
    std::chrono::microseconds walltime;
    size_t NREPEAT = 100000;

    vector<MC> MCX(NX), MCF( 1 );
    for( unsigned i=0; i<NX; ++i ) MCX[i].I( XBND0[i] ).c( XREF[i] );

    start = std::chrono::system_clock::now();
    for( unsigned i=0; i<NREPEAT; ++i )
       MCF[0] = myfunc( MCX.data() );
    walltime = std::chrono::duration_cast<std::chrono::microseconds>( std::chrono::system_clock::now() - start );
    std::cout << "McCormick walltime: " << (walltime.count() * 1e-6) / (double)NREPEAT << std::endl;

    vector<MC> MCsubX(NX), MCsubF( 1 );
    for( unsigned i=0; i<NX; ++i ) MCsubX[i].I( XBND0[i] ).c( XREF[i] ).sub( NX, i );

    start = std::chrono::system_clock::now();
    for( unsigned i=0; i<NREPEAT; ++i )
      MCsubF[0] = myfunc( MCsubX.data() );
    walltime = std::chrono::duration_cast<std::chrono::microseconds>( std::chrono::system_clock::now() - start );
    std::cout << "McCormick subgradient walltime: " << (walltime.count() * 1e-6) / (double)NREPEAT << std::endl;

    NREPEAT = 50000;

    std::vector<PWCSV> PWC16SVX(NX), PWC16SVF( 1 );
    for( unsigned i=0; i<NX; ++i )
      PWC16SVX[i].set( pwcmod, i, XBND0[i], 16 );

    start = std::chrono::system_clock::now();
    for( unsigned i=0; i<NREPEAT; ++i )
      PWC16SVF[0] = myfunc( PWC16SVX.data() );
    walltime = std::chrono::duration_cast<std::chrono::microseconds>( std::chrono::system_clock::now() - start );
    std::cout << "Superposition PWC16 walltime: " << (walltime.count() * 1e-6) / (double)NREPEAT << std::endl;

    std::vector<PWCSV> PWC32SVX(NX), PWC32SVF( 1 );
    for( unsigned i=0; i<NX; ++i )
      PWC32SVX[i].set( pwcmod, i, XBND0[i], 32 );

    start = std::chrono::system_clock::now();
    for( unsigned i=0; i<NREPEAT; ++i )
      PWC32SVF[0] = myfunc( PWC32SVX.data() );
    walltime = std::chrono::duration_cast<std::chrono::microseconds>( std::chrono::system_clock::now() - start );
    std::cout << "Superposition PWC32 walltime: " << (walltime.count() * 1e-6) / (double)NREPEAT << std::endl;

    std::vector<PWCSV> PWC64SVX(NX), PWC64SVF( 1 );
    for( unsigned i=0; i<NREPEAT; ++i )
      PWC64SVX[i].set( pwcmod, i, XBND0, 64 );

    start = std::chrono::system_clock::now();
    for( unsigned i=0; i<NREPEAT; ++i )
      PWC64SVF[0] = myfunc( PWC64SVX.data() );
    walltime = std::chrono::duration_cast<std::chrono::microseconds>( std::chrono::system_clock::now() - start );
    std::cout << "Superposition PWC64 walltime: " << (walltime.count() * 1e-6) / (double)NREPEAT << std::endl;

    std::vector<PWCSV> PWC128SVX(NX), PWC128SVF( 1 );
    for( unsigned i=0; i<NREPEAT; ++i )
      PWC128SVX[i].set( pwcmod, i, XBND0, 128 );

    start = std::chrono::system_clock::now();
    for( unsigned i=0; i<NREPEAT; ++i )
      PWC128SVF[0] = myfunc( PWC128SVX.data() );
    walltime = std::chrono::duration_cast<std::chrono::microseconds>( std::chrono::system_clock::now() - start );
    std::cout << "Superposition PWC128 walltime: " << (walltime.count() * 1e-6) / (double)NREPEAT << std::endl;

    NREPEAT = 50000;

    std::vector<PWLSV> PWL1SVX(NX), PWL1SVF( 1 );
    for( unsigned i=0; i<NREPEAT; ++i )
      PWL1SVX[i].set( pwlmod, i, XBND0, 1 );

    start = std::chrono::system_clock::now();
    for( unsigned i=0; i<NREPEAT; ++i )
      PWL1SVF[0] = myfunc( PWL1SVX.data() );
    walltime = std::chrono::duration_cast<std::chrono::microseconds>( std::chrono::system_clock::now() - start );
    std::cout << "Superposition PWL1 walltime: " << (walltime.count() * 1e-6) / (double)NREPEAT << std::endl;

    std::vector<PWLSV> PWL2SVX(NX), PWL2SVF( 1 );
    for( unsigned i=0; i<NREPEAT; ++i )
      PWL2SVX[i].set( pwlmod, i, XBND0, 2 );

    start = std::chrono::system_clock::now();
    for( unsigned i=0; i<NREPEAT; ++i )
      PWL2SVF[0] = myfunc( PWL2SVX.data() );
    walltime = std::chrono::duration_cast<std::chrono::microseconds>( std::chrono::system_clock::now() - start );
    std::cout << "Superposition PWL2 walltime: " << (walltime.count() * 1e-6) / (double)NREPEAT << std::endl;

    std::vector<PWLSV> PWL4SVX(NX), PWL4SVF( 1 );
    for( unsigned i=0; i<NREPEAT; ++i )
      PWL4SVX[i].set( pwlmod, i, XBND0, 4 );

    start = std::chrono::system_clock::now();
    for( unsigned i=0; i<NREPEAT; ++i )
      PWL4SVF[0] = myfunc( PWL4SVX.data() );
    walltime = std::chrono::duration_cast<std::chrono::microseconds>( std::chrono::system_clock::now() - start );
    std::cout << "Superposition PWL4 walltime: " << (walltime.count() * 1e-6) / (double)NREPEAT << std::endl;

    std::vector<PWLSV> PWL8SVX(NX), PWL8SVF( 1 );
    for( unsigned i=0; i<NREPEAT; ++i )
      PWL8SVX[i].set( pwlmod, i, XBND0, 8 );

    start = std::chrono::system_clock::now();
    for( unsigned i=0; i<NREPEAT; ++i )
      PWL8SVF[0] = myfunc( PWL8SVX.data() );
    walltime = std::chrono::duration_cast<std::chrono::microseconds>( std::chrono::system_clock::now() - start );
    std::cout << "Superposition PWL8 walltime: " << (walltime.count() * 1e-6) / (double)NREPEAT << std::endl;
#endif
  }

  catch( MC::Exceptions &eObj ){
    std::cerr << "Error " << eObj.ierr()
              << " in McCormick relaxation:" << std::endl
	      << eObj.what() << std::endl
              << "Aborts." << std::endl;
    return eObj.ierr();
  }

  catch( PWCSM::Exceptions &eObj ){
    std::cerr << "Error " << eObj.ierr()
              << " in superposition relaxation:" << std::endl
	      << eObj.what() << std::endl
              << "Aborts." << std::endl;
    return eObj.ierr();
  }

  catch( PWLSM::Exceptions &eObj ){
    std::cerr << "Error " << eObj.ierr()
              << " in superposition relaxation:" << std::endl
	      << eObj.what() << std::endl
              << "Aborts." << std::endl;
    return eObj.ierr();
  }

  return 0;
}
