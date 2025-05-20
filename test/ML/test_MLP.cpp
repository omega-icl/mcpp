#define PEAK_RELU_40L4       // <-- select test function here
#define  SAVE_RESULTS    // <-- specify whether to save results to file
#define ANALYSE_RATE    // <-- specify whether to analyse rate of convergence
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
#include "ffunc.hpp"

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
double const X1L = -3., X1U = 3.;
double const X2L = -3., X2U = 3.;
#include "peak_ReLU_40L4.hpp"
#endif

///////////////////////////////////////////////////////////////////////////////

int main()
{
  try{
    // Create ANN
    mc::MLP<I> f;
    f.options.RELAX     = mc::MLP<I>::Options::MC;//AUX;//MCISM;
    f.options.SUPDIV    = 64;
    f.options.INIPWL    = 1;
    f.options.SUPCONT   = true;
    f.options.SUPSHADOW = false;
    //f.options.CUTSHADOW = false;
    f.options.RELU2ABS  = false;
    unsigned l=0;
    for( auto const& layer : MLPCOEF )
      f.append_data( layer, (++l)<MLPCOEF.size()? mc::MLP<I>::Options::RELU: mc::MLP<I>::Options::LINEAR );

    // Create DAG
    mc::FFGraph dag;
    std::vector<mc::FFVar> X( NX );
    for( unsigned int i=0; i<NX; i++ ) X[i].set( &dag );
    mc::FFMLP<I> MLP;
    mc::FFVar F = MLP( 0, NX, X.data(), &f );
    auto GF = dag.subgraph( 1, &F );
    std::cout << dag;

    // Evaluate DAG
    std::vector<double> DX{ 1, 1 }, DF;
    dag.eval( GF, {F}, DF, X, DX );
    std::cout << "F = " << DF[0] << std::endl;

    // Relax DAG
    std::vector<MC> MCX{ MC( I(X1L,X1U), 1 ), MC( I(X2L,X2U), 1 ) }, MCF;
    dag.eval( GF, {F}, MCF, X, MCX );
    std::cout << "F = " << MCF[0] << std::endl;

    PWCSM pwcmod( NX );
    pwcmod.options.PROD_METH  = PWCSM::Options::PARTIAL;//FULL;//LOG;//NONE;
    pwcmod.options.PROD_CUT   = 0;
    pwcmod.options.USE_SHADOW = f.options.SUPSHADOW;
    std::vector<PWCSV> PWCSVX{ PWCSV( pwcmod, 0, I(X1L,X1U), f.options.SUPDIV ),
                               PWCSV( pwcmod, 1, I(X2L,X2U), f.options.SUPDIV ) }, PWCSVF;
    dag.eval( GF, {F}, PWCSVF, X, PWCSVX );
    std::cout << "F = " << PWCSVF[0] << std::endl;

    PWLSM pwlmod( NX );
    pwlmod.options.PROD_METH  = PWLSM::Options::PARTIAL;//PARTIAL;//FULL;//LOG;//NONE;
    pwlmod.options.PROD_CUT   = 0;
    pwlmod.options.USE_SHADOW = f.options.SUPSHADOW;
    std::vector<PWLSV> PWLSVX{ PWLSV( pwlmod, 0, I(X1L,X1U), f.options.INIPWL ),
                               PWLSV( pwlmod, 1, I(X2L,X2U), f.options.INIPWL ) }, PWLSVF;
    dag.eval( GF, {F}, PWLSVF, X, PWLSVX );
    std::cout << "F = " << PWLSVF[0] << std::endl;

#ifdef SAVE_RESULTS
    // Repeated calculations at grid points
    std::ofstream resfile( "test_MLP.out", std::ios_base::out );
    resfile << std::scientific << std::setprecision(5) << std::right;

    int const NPTS = 20;
    for( int iX1=0; iX1<NPTS; iX1++ ){
     for( int iX2=0; iX2<NPTS; iX2++ ){
       std::vector<double> DX{ X1L+iX1*(X1U-X1L)/(NPTS-1.), X2L+iX2*(X2U-X2L)/(NPTS-1.) };
       std::vector<double> DF;
       dag.eval( GF, {F}, DF, X, DX );
 
       std::vector<MC> MCX{ MC( I(X1L,X1U), DX[0] ), MC( I(X2L,X2U), DX[1] ) }, MCF;
       dag.eval( GF, {F}, MCF, X, MCX );

       // Calculate relaxations + propagate all subgradient components
       resfile << std::setw(14) << DX[0] << std::setw(14) << DX[1] << std::setw(14) << DF[0]
               << std::setw(14) << MCF[0].l()  << std::setw(14) << MCF[0].u()
               << std::setw(14) << MCF[0].cv() << std::setw(14) << MCF[0].cc()
               << std::setw(14) << PWLSVF[0].l() << std::setw(14) << PWLSVF[0].u()
               << std::setw(14) << PWLSVF[0].uval({{0,DX[0]},{1,DX[1]}}) << std::setw(14) << PWLSVF[0].oval({{0,DX[0]},{1,DX[1]}})
               << std::endl;

     }
     resfile << std::endl;
    }
#endif

#ifdef ANALYSE_TIME
    std::chrono::time_point<std::chrono::system_clock> start;
    std::chrono::microseconds walltime;
    size_t NREPEAT = 10000;

    start = std::chrono::system_clock::now();
    for( unsigned i=0; i<NREPEAT; ++i )
      dag.eval( GF, {F}, MCF, X, MCX );
    walltime = std::chrono::duration_cast<std::chrono::microseconds>( std::chrono::system_clock::now() - start );
    std::cout << "McCormick walltime: " << (walltime.count() * 1e-6) / (double)NREPEAT << std::endl;

    std::vector<MC> MCsubX{ MC( I(X1L,X1U), 1 ), MC( I(X2L,X2U), 1 ) }, MCsubF;
    MCsubX[0].sub( 2, 0 );
    MCsubX[1].sub( 2, 1 );

    start = std::chrono::system_clock::now();
    for( unsigned i=0; i<NREPEAT; ++i )
      dag.eval( GF, {F}, MCsubF, X, MCsubX );
    walltime = std::chrono::duration_cast<std::chrono::microseconds>( std::chrono::system_clock::now() - start );
    std::cout << "McCormick subgradient walltime: " << (walltime.count() * 1e-6) / (double)NREPEAT << std::endl;

    std::vector<PWCSV> PWC16SVX(NX), PWC16SVF;
    PWC16SVX[0].set( pwcmod, 0, I(X1L,X1U), 16 );
    PWC16SVX[1].set( pwcmod, 1, I(X2L,X2U), 16 );
    dag.eval( GF, {F}, PWC16SVF, X, PWC16SVX );

    start = std::chrono::system_clock::now();
    for( unsigned i=0; i<NREPEAT; ++i )
      dag.eval( GF, {F}, PWC16SVF, X, PWC16SVX );
    walltime = std::chrono::duration_cast<std::chrono::microseconds>( std::chrono::system_clock::now() - start );
    std::cout << "Superposition PWC16 walltime: " << (walltime.count() * 1e-6) / (double)NREPEAT << std::endl;

    std::vector<PWCSV> PWC32SVX(NX), PWC32SVF;
    PWC32SVX[0].set( pwcmod, 0, I(X1L,X1U), 32 );
    PWC32SVX[1].set( pwcmod, 1, I(X2L,X2U), 32 );
    dag.eval( GF, {F}, PWC32SVF, X, PWC32SVX );

    start = std::chrono::system_clock::now();
    for( unsigned i=0; i<NREPEAT; ++i )
      dag.eval( GF, {F}, PWC32SVF, X, PWC32SVX );
    walltime = std::chrono::duration_cast<std::chrono::microseconds>( std::chrono::system_clock::now() - start );
    std::cout << "Superposition PWC32 walltime: " << (walltime.count() * 1e-6) / (double)NREPEAT << std::endl;

    std::vector<PWCSV> PWC64SVX(NX), PWC64SVF;
    PWC64SVX[0].set( pwcmod, 0, I(X1L,X1U), 64 );
    PWC64SVX[1].set( pwcmod, 1, I(X2L,X2U), 64 );
    dag.eval( GF, {F}, PWC64SVF, X, PWC64SVX );

    start = std::chrono::system_clock::now();
    for( unsigned i=0; i<NREPEAT; ++i )
      dag.eval( GF, {F}, PWC64SVF, X, PWC64SVX );
    walltime = std::chrono::duration_cast<std::chrono::microseconds>( std::chrono::system_clock::now() - start );
    std::cout << "Superposition PWC64 walltime: " << (walltime.count() * 1e-6) / (double)NREPEAT << std::endl;

    std::vector<PWCSV> PWC128SVX(NX), PWC128SVF;
    PWC128SVX[0].set( pwcmod, 0, I(X1L,X1U), 128 );
    PWC128SVX[1].set( pwcmod, 1, I(X2L,X2U), 128 );
    dag.eval( GF, {F}, PWC128SVF, X, PWC128SVX );

    start = std::chrono::system_clock::now();
    for( unsigned i=0; i<NREPEAT; ++i )
      dag.eval( GF, {F}, PWC128SVF, X, PWC128SVX );
    walltime = std::chrono::duration_cast<std::chrono::microseconds>( std::chrono::system_clock::now() - start );
    std::cout << "Superposition PWC128 walltime: " << (walltime.count() * 1e-6) / (double)NREPEAT << std::endl;

    NREPEAT = 2000;
    start = std::chrono::system_clock::now();
    for( unsigned i=0; i<NREPEAT; ++i )
      dag.eval( GF, {F}, PWLSVF, X, PWLSVX );
    walltime = std::chrono::duration_cast<std::chrono::microseconds>( std::chrono::system_clock::now() - start );
    std::cout << "Superposition PWL walltime: " << (walltime.count() * 1e-6) / (double)NREPEAT << std::endl;
#endif

#ifdef ANALYSE_RATE
    std::vector<I> const XBND0{ I(X1L,X1U), I(X2L,X2U) };
    std::map<unsigned,double> XREF{ {0,-1.497271e-02}, {1,1.617562e+00} }; // Global max
    //std::map<unsigned,double> XREF{ {0,2.034074e-01}, {1,-1.637012e+00} }; // Global min
    auto const& red = [=]( const I& bnd, const double& ref, const double& r ){ return r*bnd + (1-r)*ref; };
    auto const& min = [=]( const double& x, const double& y ){ return x<y?x:y; };
    auto const& max = [=]( const double& x, const double& y ){ return x>y?x:y; };

    size_t const NGRID  = 100;
    double const rate   = 0.95;
    double const rhomin = 1e-5;
    for( double rho = 1.; rho > rhomin; rho *= rate ){
      double min_F = DF[0],           max_F = DF[0];
      double min_MCFcv = DF[0],       max_MCFcc = DF[0];
      double distmax_MCFcv = 0.,      distmax_MCFcc = 0.;
      double distmax_PWC16SVFu = 0.,  distmax_PWC16SVFo = 0.;
      double distmax_PWC32SVFu = 0.,  distmax_PWC32SVFo = 0.;
      double distmax_PWC64SVFu = 0.,  distmax_PWC64SVFo = 0.;
      double distmax_PWC128SVFu = 0., distmax_PWC128SVFo = 0.;
      double distmax_PWL1SVFu = 0.,   distmax_PWL1SVFo = 0.;

      std::vector<I> XBND{ red( XBND0[0], XREF[0], rho ), red( XBND0[1], XREF[1], rho ) };

      std::vector<PWCSV> PWC16SVX(NX), PWC16SVF;
      PWC16SVX[0].set( pwcmod, 0, XBND[0], 16 );
      PWC16SVX[1].set( pwcmod, 1, XBND[1], 16 );
      dag.eval( GF, {F}, PWC16SVF, X, PWC16SVX );

      std::vector<PWCSV> PWC32SVX(NX), PWC32SVF;
      PWC32SVX[0].set( pwcmod, 0, XBND[0], 32 );
      PWC32SVX[1].set( pwcmod, 1, XBND[1], 32 );
      dag.eval( GF, {F}, PWC32SVF, X, PWC32SVX );

      std::vector<PWCSV> PWC64SVX(NX), PWC64SVF;
      PWC64SVX[0].set( pwcmod, 0, XBND[0], 64 );
      PWC64SVX[1].set( pwcmod, 1, XBND[1], 64 );
      dag.eval( GF, {F}, PWC64SVF, X, PWC64SVX );

      std::vector<PWCSV> PWC128SVX(NX), PWC128SVF;
      PWC128SVX[0].set( pwcmod, 0, XBND[0], 128 );
      PWC128SVX[1].set( pwcmod, 1, XBND[1], 128 );
      dag.eval( GF, {F}, PWC128SVF, X, PWC128SVX );

      std::vector<PWLSV> PWL1SVX(NX), PWL1SVF;
      PWL1SVX[0].set( pwlmod, 0, XBND[0], 1 );
      PWL1SVX[1].set( pwlmod, 1, XBND[1], 1 );
      dag.eval( GF, {F}, PWL1SVF, X, PWL1SVX );

      for( unsigned iX1=0; iX1<NGRID; iX1++ ){
        for( unsigned iX2=0; iX2<NGRID; iX2++ ){
          std::vector<double> DX{ mc::Op<I>::l(XBND[0])+iX1*mc::Op<I>::diam(XBND[0])/(NGRID-1.),
                                  mc::Op<I>::l(XBND[1])+iX2*mc::Op<I>::diam(XBND[1])/(NGRID-1.) };
          std::vector<double> DF;
          dag.eval( GF, {F}, DF, X, DX );
 
          min_F = min( min_F, DF[0] );
          max_F = max( max_F, DF[0] );

          std::vector<MC> MCX{ MC( XBND[0], DX[0] ), MC( XBND[1], DX[1] ) }, MCF;
          dag.eval( GF, {F}, MCF, X, MCX );

          min_MCFcv = min( min_MCFcv, MCF[0].cv() );
          max_MCFcc = max( max_MCFcc, MCF[0].cc() );

          distmax_MCFcv = max( distmax_MCFcv, DF[0] - MCF[0].cv() );
          distmax_MCFcc = max( distmax_MCFcc, MCF[0].cc() - DF[0] );

          double const PWC16SVFu = PWC16SVF[0].uval({{0,DX[0]},{1,DX[1]}});
          double const PWC16SVFo = PWC16SVF[0].oval({{0,DX[0]},{1,DX[1]}});

          distmax_PWC16SVFu = max( distmax_PWC16SVFu, DF[0] - PWC16SVFu );
          distmax_PWC16SVFo = max( distmax_PWC16SVFo, PWC16SVFo - DF[0] );

          double const PWC32SVFu = PWC32SVF[0].uval({{0,DX[0]},{1,DX[1]}});
          double const PWC32SVFo = PWC32SVF[0].oval({{0,DX[0]},{1,DX[1]}});

          distmax_PWC32SVFu = max( distmax_PWC32SVFu, DF[0] - PWC32SVFu );
          distmax_PWC32SVFo = max( distmax_PWC32SVFo, PWC32SVFo - DF[0] );

          double const PWC64SVFu = PWC64SVF[0].uval({{0,DX[0]},{1,DX[1]}});
          double const PWC64SVFo = PWC64SVF[0].oval({{0,DX[0]},{1,DX[1]}});

          distmax_PWC64SVFu = max( distmax_PWC64SVFu, DF[0] - PWC64SVFu );
          distmax_PWC64SVFo = max( distmax_PWC64SVFo, PWC64SVFo - DF[0] );

          double const PWC128SVFu = PWC128SVF[0].uval({{0,DX[0]},{1,DX[1]}});
          double const PWC128SVFo = PWC128SVF[0].oval({{0,DX[0]},{1,DX[1]}});

          distmax_PWC128SVFu = max( distmax_PWC128SVFu, DF[0] - PWC128SVFu );
          distmax_PWC128SVFo = max( distmax_PWC128SVFo, PWC128SVFo - DF[0] );

          double const PWL1SVFu = PWL1SVF[0].uval({{0,DX[0]},{1,DX[1]}});
          double const PWL1SVFo = PWL1SVF[0].oval({{0,DX[0]},{1,DX[1]}});

          distmax_PWL1SVFu = max( distmax_PWL1SVFu, DF[0] - PWL1SVFu );
          distmax_PWL1SVFo = max( distmax_PWL1SVFo, PWL1SVFo - DF[0] );
        }
      }
      
      std::cout << std::scientific << std::setprecision(5) << std::right
                << std::setw(14) << rho
                << std::setw(14) << mc::Op<I>::l(XBND[0]) << std::setw(14) << mc::Op<I>::u(XBND[0])
                << std::setw(14) << mc::Op<I>::l(XBND[1]) << std::setw(14) << mc::Op<I>::u(XBND[1])
                //<< std::setw(14) << min_MCFcv << std::setw(14) << max_MCFcc
                << std::setw(14) << max( distmax_MCFcv,            distmax_MCFcc )
                << std::setw(14) << max( min_F - min_MCFcv,        max_MCFcc - max_F )
                << std::setw(14) << max( distmax_PWC16SVFu,        distmax_PWC16SVFo )
                << std::setw(14) << max( min_F - PWC16SVF[0].l(),  PWC16SVF[0].u() - max_F )
                << std::setw(14) << max( distmax_PWC32SVFu,        distmax_PWC32SVFo )
                << std::setw(14) << max( min_F - PWC32SVF[0].l(),  PWC32SVF[0].u() - max_F )
                << std::setw(14) << max( distmax_PWC64SVFu,        distmax_PWC64SVFo )
                << std::setw(14) << max( min_F - PWC64SVF[0].l(),  PWC64SVF[0].u() - max_F )
                << std::setw(14) << max( distmax_PWC128SVFu,       distmax_PWC128SVFo )
                << std::setw(14) << max( min_F - PWC128SVF[0].l(), PWC128SVF[0].u() - max_F )
                << std::setw(14) << max( distmax_PWL1SVFu,         distmax_PWL1SVFo )
                << std::setw(14) << max( min_F - PWL1SVF[0].l(),   PWL1SVF[0].u() - max_F )
                << std::endl;
    }
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

