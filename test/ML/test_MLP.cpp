#define PEAK_RELU_40L4       // <-- select test function here
#undef SAVE_RESULTS    // <-- specify whether to save results to file
#undef ANALYSE_RATE    // <-- specify whether to analyse rate of convergence
#undef ANALYSE_TIME    // <-- specify whether to analyse computational time

#undef MC__FFMLP_CHECK
#undef MC__FFMLP_DEBUG

////////////////////////////////////////////////////////////////////////

#include <fstream>
#include <iomanip>
#include <chrono>

#if defined( MC__USE_PROFIL )
 #include "mcprofil.hpp"
 typedef INTERVAL I;
#elif defined( MC__USE_FILIB )
 #include "mcfilib.hpp"
 typedef filib::interval<double,filib::native_switched,filib::i_mode_extended> I;
#elif defined( MC__USE_BOOST )
 #include "mcboost.hpp"
 typedef boost::numeric::interval_lib::save_state<boost::numeric::interval_lib::rounded_transc_opp<double>> T_boost_round;
 typedef boost::numeric::interval_lib::checking_base<double> T_boost_check;
 typedef boost::numeric::interval_lib::policies<T_boost_round,T_boost_check> T_boost_policy;
 typedef boost::numeric::interval<double,T_boost_policy> I;
#else
 #include "interval.hpp"
 typedef mc::Interval I;
#endif

#include "ffmlp.hpp"

typedef mc::McCormick<I> MC;
typedef mc::SupModel<mc::PWCU> PWCSM;
typedef mc::SupVar<mc::PWCU> PWCSV;
typedef mc::SupModel<mc::PWLU> PWLSM;
typedef mc::SupVar<mc::PWLU> PWLSV;

#if defined( TEST_RELU )
size_t const NX = 2;
double const X1L = -2., X1U = 2.; // [-10:10]
double const X2L = -2., X2U = 2.; // [-10:10]
std::vector<std::vector<std::vector<double>>> const MLPCOEF =
{ { { 0.0, 0.2, 0.3 },
    { -3.0, 0.5, -0.2 },
    { 0.0, 0.2, -0.4 },
    { 0.0, -0.5, 0.0 }           },
  { { 0.0, 1.0, -1.0, 1.0, 1.0 } } };

#elif defined( PEAK_RELU_30L1 )
size_t const NX = 2;
double const X1L = -3., X1U = 3.;
double const X2L = -3., X2U = 3.;
#include "peak_ReLU_30L1.hpp"

#elif defined( PEAK_RELU_40L4 )
size_t const NX = 2;
//double const X1L = -1., X1U = 1.;
//double const X2L = -1., X2U = 1.;
double const X1L = -3., X1U = 3.;
double const X2L = -3., X2U = 3.;
#include "peak_ReLU_40L4.hpp"
#endif

///////////////////////////////////////////////////////////////////////////////

int main()
{
  try{
    // Create ANN
    mc::MLP<I> ANN;
    ANN.options.RELU2ABS  = false;
    unsigned l=0;
    for( auto const& layer : MLPCOEF )
      ANN.append_data( layer, (++l)<MLPCOEF.size()? ANN.RELU: ANN.LINEAR );

    // Create DAG
    mc::FFGraph DAG;
    DAG.options.MAXTHREAD = 0;
    std::vector<mc::FFVar> X = DAG.add_vars( NX, "X" );
    mc::FFMLP<I> OpMLP;
    OpMLP.options.RELAX  = { OpMLP.options.AUX }; //PWLS PWCS MC AUX INT
    OpMLP.options.PWCDIV = 32;
    OpMLP.options.PWCREL = 1;
    OpMLP.options.PWCSUP.USE_SHADOW = 1;
    OpMLP.options.PWCSHADOW = 1;
    OpMLP.options.PWLINI = 1;
    OpMLP.options.PWLREL = 1;
    OpMLP.options.PWLSUP.USE_SHADOW = 1;
    OpMLP.options.PWLSHADOW = 1;
    std::vector<mc::FFVar> Y{ OpMLP( 0, X, &ANN, OpMLP.COPY ) };

    // Create and display subgraph
    auto SgY = DAG.subgraph( Y );
    //DAG.output( SgY, " OF Y" );
    auto StrY = mc::FFExpr::subgraph( &DAG, SgY );
    std::cout << "Y: " << StrY[0] << std::endl;
    
    // Evaluate DAG
    std::vector<double> DX{ 1, 1 }, DY(1);
    ANN.eval( DX.data(), DY.data() );
    std::cout << "Y(" << DX[0] << "," << DX[1] << ") = " << DY[0] << std::endl;
    DAG.eval( SgY, Y, DY, X, DX );
    std::cout << "Y(" << DX[0] << "," << DX[1] << ") = " << DY[0] << std::endl;

    // Vectorized evaluation
    std::vector<std::vector<double>> vDX{ { 1, 1 }, { -1, -1 }, { 2, 2 }, { -2, -2 } }, vDY;
    DAG.veval( SgY, Y, vDY, X, vDX );
    for( size_t i=0; i<vDX.size(); ++i )
      std::cout << "Y(" << vDX[i][0] << "," << vDX[i][1] << ") = " << vDY[i][0] << std::endl;

    // Differentiate and evaluate DAG
    std::vector<mc::FFVar> dYdX = DAG.FAD( Y, X );
    auto SgdYdX = DAG.subgraph( dYdX );
    auto StrdYdX = mc::FFExpr::subgraph( &DAG, SgdYdX );
    for( unsigned i=0; i<StrdYdX.size(); ++i )
      std::cout << "dYdX[" << i << "]: " << StrdYdX[i] << std::endl;

    std::vector<double> DdYdX;
    DAG.eval( SgdYdX, dYdX, DdYdX, X, DX );
    for( unsigned i=0; i<StrdYdX.size(); ++i )
      std::cout << "dYdX[" << i << "](" << DX[0] << "," << DX[1] << ") = " << DdYdX[i] << std::endl;

    // McCormick relaxation
    std::vector<I> IX{ I(X1L,X1U), I(X2L,X2U) };
    std::vector<MC> MCX{ MC( IX[0], DX[0] ), MC( IX[1], DX[1] ) }, MCY;
    DAG.eval( SgY, Y, MCY, X, MCX );
    std::cout << "Y = " << MCY[0] << std::endl;

    // Constraint propagation
    std::vector<I> IY{ I(-10.,10.) }, Iwk;
    //ANN.reval( IX.data(), IY.data() );
    DAG.reval( SgY, Iwk, Y, IY, X, IX, 1E20*I(-1,1), 10, 0. );
    std::cout << "Y(" << DX[0] << "," << DX[1] << ") in " << IY[0] << std::endl;

    // Polyhedral relaxation
    mc::PolImg<I> PolEnv;
    PolEnv.options.BREAKPOINT_TYPE = mc::PolImg<I>::Options::BIN;//CONT;//BIN;//SOS2;
    PolEnv.options.AGGREG_LQ       = true;
    //PolEnv.options.BREAKPOINT_RTOL =
    //PolEnv.options.BREAKPOINT_ATOL = 0e0;
    //PolEnv.options.ALLOW_DISJ = { mc::FFOp::FSTEP, mc::FFOp::MAXF, mc::FFOp::MINF, mc::FFOp::FABS };
    //PolEnv.options.ALLOW_NLIN      = { mc::FFOp::TANH, mc::FFOp::EXP  };
    //PolEnv.options.SANDWICH_MAXCUT = 6;
    std::vector<mc::PolVar<I>> PolX{ mc::PolVar<I>( &PolEnv, X[0], IX[0] ), mc::PolVar<I>( &PolEnv, X[1], IX[1] ) }, PolY;
    DAG.eval( SgY, Y, PolY, X, PolX );
    PolEnv.generate_cuts( PolY );
    std::cout << "\n Polyhedral relaxation:" << PolEnv << std::endl;

    // Piecewise-constant superposition relaxation
    PWCSM pwcmod( NX );
    pwcmod.options.PROD_METH  = PWCSM::Options::PARTIAL;//FULL;//LOG;//NONE;
    pwcmod.options.PROD_CUT   = 0;
    pwcmod.options.USE_SHADOW = 0;//OpMLP.options.PWCSHADOW;
    std::vector<PWCSV> PWCSVX{ PWCSV( pwcmod, 0, I(X1L,X1U), OpMLP.options.PWCDIV ),
                               PWCSV( pwcmod, 1, I(X2L,X2U), OpMLP.options.PWCDIV ) }, PWCSVY;
    DAG.eval( SgY, Y, PWCSVY, X, PWCSVX );
    std::cout << "Y = " << PWCSVY[0] << std::endl;

    // Piecewise-linear superposition relaxation
    PWLSM pwlmod( NX );
    pwlmod.options.PROD_METH  = PWLSM::Options::PARTIAL;//PARTIAL;//FULL;//LOG;//NONE;
    pwlmod.options.PROD_CUT   = 0;
    pwlmod.options.USE_SHADOW = 0;//OpMLP.options.PWLSHADOW;
    pwlmod.options.MAX_SUBDIV = 0;//16;
    std::vector<PWLSV> PWLSVX{ PWLSV( pwlmod, 0, I(X1L,X1U), OpMLP.options.PWLINI ),
                               PWLSV( pwlmod, 1, I(X2L,X2U), OpMLP.options.PWLINI ) }, PWLSVY;
    DAG.eval( SgY, Y, PWLSVY, X, PWLSVX );
    std::cout << "Y = " << PWLSVY[0] << std::endl;

#ifdef SAVE_RESULTS
    // Repeated calculations at grid points
    std::ofstream resfile( "test_MLP.out", std::ios_base::out );
    resfile << std::scientific << std::setprecision(5) << std::right;

    int const NPTS = 100;
    for( int iX1=0; iX1<NPTS; iX1++ ){
     for( int iX2=0; iX2<NPTS; iX2++ ){
       std::vector<double> DX{ X1L+iX1*(X1U-X1L)/(NPTS-1.), X2L+iX2*(X2U-X2L)/(NPTS-1.) };
       std::vector<double> DY;
       DAG.eval( SgY, Y, DY, X, DX );
 
       std::vector<MC> MCX{ MC( I(X1L,X1U), DX[0] ), MC( I(X2L,X2U), DX[1] ) }, MCY;
       DAG.eval( SgY, Y, MCY, X, MCX );

       // Calculate relaxations + propagate all subgradient components
       resfile << std::setw(14) << DX[0] << std::setw(14) << DX[1] << std::setw(14) << DY[0]
               << std::setw(14) << MCY[0].l()  << std::setw(14) << MCY[0].u()
               << std::setw(14) << MCY[0].cv() << std::setw(14) << MCY[0].cc()
               << std::setw(14) << PWLSVY[0].l() << std::setw(14) << PWLSVY[0].u()
               << std::setw(14) << PWLSVY[0].uval({{0,DX[0]},{1,DX[1]}}) << std::setw(14) << PWLSVY[0].oval({{0,DX[0]},{1,DX[1]}})
               << std::endl;

     }
     resfile << std::endl;
    }
#endif

    std::vector<I> const XBND0{ I(X1L,X1U), I(X2L,X2U) };
    std::map<unsigned,double> XREF{ {0,-1.497271e-02}, {1,1.617562e+00} }; // Global max
    //std::map<unsigned,double> XREF{ {0,2.034074e-01}, {1,-1.637012e+00} }; // Global min

#ifdef ANALYSE_RATE
    auto const& red = [=]( const I& bnd, const double& ref, const double& r ){ return r*bnd + (1-r)*ref; };
    auto const& min = [=]( const double& x, const double& y ){ return x<y?x:y; };
    auto const& max = [=]( const double& x, const double& y ){ return x>y?x:y; };

    size_t const NGRID  = 100;
    double const rate   = 0.95;
    double const rhomin = 1e-5;
    for( double rho = 1.; rho > rhomin; rho *= rate ){
      double min_Y = DY[0],           max_Y = DY[0];
      double min_MCYcv = DY[0],       max_MCYcc = DY[0];
      double distmax_MCYcv = 0.,      distmax_MCYcc = 0.;
      double distmax_PWC16SVYu = 0.,  distmax_PWC16SVYo = 0.;
      double distmax_PWC32SVYu = 0.,  distmax_PWC32SVYo = 0.;
      double distmax_PWC64SVYu = 0.,  distmax_PWC64SVYo = 0.;
      double distmax_PWC128SVYu = 0., distmax_PWC128SVYo = 0.;
      double distmax_PWL1SVYu = 0.,   distmax_PWL1SVYo = 0.;

      std::vector<I> XBND{ red( XBND0[0], XREF[0], rho ), red( XBND0[1], XREF[1], rho ) };

      std::vector<PWCSV> PWC16SVX(NX), PWC16SVY;
      PWC16SVX[0].set( pwcmod, 0, XBND[0], 16 );
      PWC16SVX[1].set( pwcmod, 1, XBND[1], 16 );
      DAG.eval( SgY, Y, PWC16SVY, X, PWC16SVX );

      std::vector<PWCSV> PWC32SVX(NX), PWC32SVY;
      PWC32SVX[0].set( pwcmod, 0, XBND[0], 32 );
      PWC32SVX[1].set( pwcmod, 1, XBND[1], 32 );
      DAG.eval( SgY, Y, PWC32SVY, X, PWC32SVX );

      std::vector<PWCSV> PWC64SVX(NX), PWC64SVY;
      PWC64SVX[0].set( pwcmod, 0, XBND[0], 64 );
      PWC64SVX[1].set( pwcmod, 1, XBND[1], 64 );
      DAG.eval( SgY, Y, PWC64SVY, X, PWC64SVX );

      std::vector<PWCSV> PWC128SVX(NX), PWC128SVY;
      PWC128SVX[0].set( pwcmod, 0, XBND[0], 128 );
      PWC128SVX[1].set( pwcmod, 1, XBND[1], 128 );
      DAG.eval( SgY, Y, PWC128SVY, X, PWC128SVX );

      std::vector<PWLSV> PWL1SVX(NX), PWL1SVY;
      PWL1SVX[0].set( pwlmod, 0, XBND[0], 1 );
      PWL1SVX[1].set( pwlmod, 1, XBND[1], 1 );
      DAG.eval( SgY, Y, PWL1SVY, X, PWL1SVX );

      for( unsigned iX1=0; iX1<NGRID; iX1++ ){
        for( unsigned iX2=0; iX2<NGRID; iX2++ ){
          std::vector<double> DX{ mc::Op<I>::l(XBND[0])+iX1*mc::Op<I>::diam(XBND[0])/(NGRID-1.),
                                  mc::Op<I>::l(XBND[1])+iX2*mc::Op<I>::diam(XBND[1])/(NGRID-1.) };
          std::vector<double> DY;
          DAG.eval( SgY, Y, DY, X, DX );
 
          min_Y = min( min_Y, DY[0] );
          max_Y = max( max_Y, DY[0] );

          std::vector<MC> MCX{ MC( XBND[0], DX[0] ), MC( XBND[1], DX[1] ) }, MCY;
          DAG.eval( SgY, Y, MCY, X, MCX );

          min_MCYcv = min( min_MCYcv, MCY[0].cv() );
          max_MCYcc = max( max_MCYcc, MCY[0].cc() );

          distmax_MCYcv = max( distmax_MCYcv, DY[0] - MCY[0].cv() );
          distmax_MCYcc = max( distmax_MCYcc, MCY[0].cc() - DY[0] );

          double const PWC16SVYu = PWC16SVY[0].uval({{0,DX[0]},{1,DX[1]}});
          double const PWC16SVYo = PWC16SVY[0].oval({{0,DX[0]},{1,DX[1]}});

          distmax_PWC16SVYu = max( distmax_PWC16SVYu, DY[0] - PWC16SVYu );
          distmax_PWC16SVYo = max( distmax_PWC16SVYo, PWC16SVYo - DY[0] );

          double const PWC32SVYu = PWC32SVY[0].uval({{0,DX[0]},{1,DX[1]}});
          double const PWC32SVYo = PWC32SVY[0].oval({{0,DX[0]},{1,DX[1]}});

          distmax_PWC32SVYu = max( distmax_PWC32SVYu, DY[0] - PWC32SVYu );
          distmax_PWC32SVYo = max( distmax_PWC32SVYo, PWC32SVYo - DY[0] );

          double const PWC64SVYu = PWC64SVY[0].uval({{0,DX[0]},{1,DX[1]}});
          double const PWC64SVYo = PWC64SVY[0].oval({{0,DX[0]},{1,DX[1]}});

          distmax_PWC64SVYu = max( distmax_PWC64SVYu, DY[0] - PWC64SVYu );
          distmax_PWC64SVYo = max( distmax_PWC64SVYo, PWC64SVYo - DY[0] );

          double const PWC128SVYu = PWC128SVY[0].uval({{0,DX[0]},{1,DX[1]}});
          double const PWC128SVYo = PWC128SVY[0].oval({{0,DX[0]},{1,DX[1]}});

          distmax_PWC128SVYu = max( distmax_PWC128SVYu, DY[0] - PWC128SVYu );
          distmax_PWC128SVYo = max( distmax_PWC128SVYo, PWC128SVYo - DY[0] );

          double const PWL1SVYu = PWL1SVY[0].uval({{0,DX[0]},{1,DX[1]}});
          double const PWL1SVYo = PWL1SVY[0].oval({{0,DX[0]},{1,DX[1]}});

          distmax_PWL1SVYu = max( distmax_PWL1SVYu, DY[0] - PWL1SVYu );
          distmax_PWL1SVYo = max( distmax_PWL1SVYo, PWL1SVYo - DY[0] );
        }
      }
      
      std::cout << std::scientific << std::setprecision(5) << std::right
                << std::setw(14) << rho
                << std::setw(14) << mc::Op<I>::l(XBND[0]) << std::setw(14) << mc::Op<I>::u(XBND[0])
                << std::setw(14) << mc::Op<I>::l(XBND[1]) << std::setw(14) << mc::Op<I>::u(XBND[1])
                //<< std::setw(14) << min_MCYcv << std::setw(14) << max_MCYcc
                << std::setw(14) << max( distmax_MCYcv,            distmax_MCYcc )
                << std::setw(14) << max( min_Y - min_MCYcv,        max_MCYcc - max_Y )
                << std::setw(14) << max( distmax_PWC16SVYu,        distmax_PWC16SVYo )
                << std::setw(14) << max( min_Y - PWC16SVY[0].l(),  PWC16SVY[0].u() - max_Y )
                << std::setw(14) << max( distmax_PWC32SVYu,        distmax_PWC32SVYo )
                << std::setw(14) << max( min_Y - PWC32SVY[0].l(),  PWC32SVY[0].u() - max_Y )
                << std::setw(14) << max( distmax_PWC64SVYu,        distmax_PWC64SVYo )
                << std::setw(14) << max( min_Y - PWC64SVY[0].l(),  PWC64SVY[0].u() - max_Y )
                << std::setw(14) << max( distmax_PWC128SVYu,       distmax_PWC128SVYo )
                << std::setw(14) << max( min_Y - PWC128SVY[0].l(), PWC128SVY[0].u() - max_Y )
                << std::setw(14) << max( distmax_PWL1SVYu,         distmax_PWL1SVYo )
                << std::setw(14) << max( min_Y - PWL1SVY[0].l(),   PWL1SVY[0].u() - max_Y )
                << std::endl;
    }
#endif

#ifdef ANALYSE_TIME
    std::chrono::time_point<std::chrono::system_clock> start;
    std::chrono::microseconds walltime;
    size_t NREPEAT = 1000;

    MCX[0].c( XREF[0] );
    MCX[1].c( XREF[1] );

    start = std::chrono::system_clock::now();
    std::vector<MC> MCwk;
    for( unsigned i=0; i<NREPEAT; ++i )
      DAG.eval( SgY, MCwk, Y, MCY, X, MCX );
    walltime = std::chrono::duration_cast<std::chrono::microseconds>( std::chrono::system_clock::now() - start );
    std::cout << "McCormick walltime: " << (walltime.count() * 1e-6) / (double)NREPEAT << std::endl;

    std::vector<MC> MCsubX{ MC( XBND0[0], XREF[0] ), MC( XBND0[1], XREF[1] ) }, MCsubY;
    MCsubX[0].sub( 2, 0 );
    MCsubX[1].sub( 2, 1 );

    start = std::chrono::system_clock::now();
    for( unsigned i=0; i<NREPEAT; ++i )
      DAG.eval( SgY, Y, MCsubY, X, MCsubX );
    walltime = std::chrono::duration_cast<std::chrono::microseconds>( std::chrono::system_clock::now() - start );
    std::cout << "McCormick subgradient walltime: " << (walltime.count() * 1e-6) / (double)NREPEAT << std::endl;

    std::vector<PWCSV> PWC16SVX(NX), PWC16SVY, PWCwk;
    PWC16SVX[0].set( pwcmod, 0, I(X1L,X1U), 16 );
    PWC16SVX[1].set( pwcmod, 1, I(X2L,X2U), 16 );
    DAG.eval( SgY, PWCwk, Y, PWC16SVY, X, PWC16SVX );

    start = std::chrono::system_clock::now();
    for( unsigned i=0; i<NREPEAT; ++i )
      DAG.eval( SgY, Y, PWC16SVY, X, PWC16SVX );
    walltime = std::chrono::duration_cast<std::chrono::microseconds>( std::chrono::system_clock::now() - start );
    std::cout << "Superposition PWC16 walltime: " << (walltime.count() * 1e-6) / (double)NREPEAT << std::endl;

    std::vector<PWCSV> PWC32SVX(NX), PWC32SVY;
    PWC32SVX[0].set( pwcmod, 0, I(X1L,X1U), 32 );
    PWC32SVX[1].set( pwcmod, 1, I(X2L,X2U), 32 );
    DAG.eval( SgY, PWCwk, Y, PWC32SVY, X, PWC32SVX );

    start = std::chrono::system_clock::now();
    for( unsigned i=0; i<NREPEAT; ++i )
      DAG.eval( SgY, Y, PWC32SVY, X, PWC32SVX );
    walltime = std::chrono::duration_cast<std::chrono::microseconds>( std::chrono::system_clock::now() - start );
    std::cout << "Superposition PWC32 walltime: " << (walltime.count() * 1e-6) / (double)NREPEAT << std::endl;

    std::vector<PWCSV> PWC64SVX(NX), PWC64SVY;
    PWC64SVX[0].set( pwcmod, 0, I(X1L,X1U), 64 );
    PWC64SVX[1].set( pwcmod, 1, I(X2L,X2U), 64 );
    DAG.eval( SgY, PWCwk, Y, PWC64SVY, X, PWC64SVX );

    start = std::chrono::system_clock::now();
    for( unsigned i=0; i<NREPEAT; ++i )
      DAG.eval( SgY, Y, PWC64SVY, X, PWC64SVX );
    walltime = std::chrono::duration_cast<std::chrono::microseconds>( std::chrono::system_clock::now() - start );
    std::cout << "Superposition PWC64 walltime: " << (walltime.count() * 1e-6) / (double)NREPEAT << std::endl;

    std::vector<PWCSV> PWC128SVX(NX), PWC128SVY;
    PWC128SVX[0].set( pwcmod, 0, I(X1L,X1U), 128 );
    PWC128SVX[1].set( pwcmod, 1, I(X2L,X2U), 128 );
    DAG.eval( SgY, PWCwk, Y, PWC128SVY, X, PWC128SVX );

    start = std::chrono::system_clock::now();
    for( unsigned i=0; i<NREPEAT; ++i )
      DAG.eval( SgY, PWCwk, Y, PWC128SVY, X, PWC128SVX );
    walltime = std::chrono::duration_cast<std::chrono::microseconds>( std::chrono::system_clock::now() - start );
    std::cout << "Superposition PWC128 walltime: " << (walltime.count() * 1e-6) / (double)NREPEAT << std::endl;

    //NREPEAT = 2000;
    std::vector<PWLSV> PWLwk;
    start = std::chrono::system_clock::now();
    for( unsigned i=0; i<NREPEAT; ++i )
      DAG.eval( SgY, PWLwk, Y, PWLSVY, X, PWLSVX );
    walltime = std::chrono::duration_cast<std::chrono::microseconds>( std::chrono::system_clock::now() - start );
    std::cout << "Superposition PWL walltime: " << (walltime.count() * 1e-6) / (double)NREPEAT << std::endl;
#endif

  }

  catch( mc::FFBase::Exceptions &eObj ){
    std::cerr << "Error " << eObj.ierr()
              << " in factorable function manipulation:" << std::endl
              << eObj.what() << std::endl
              << "Aborts." << std::endl;
    return eObj.ierr();
  }

#if !defined(MC__USE_PROFIL) && !defined(MC__USE_FILIB) && !defined(MC__USE_BOOST)
  catch( I::Exceptions &eObj ){
    std::cerr << "Error " << eObj.ierr()
              << " in natural interval extension:" << std::endl
	      << eObj.what() << std::endl
              << "Aborts." << std::endl;
    return eObj.ierr();
  }
#endif

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

