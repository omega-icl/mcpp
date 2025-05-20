#define TEST_EXP	// <-- select test function here
#undef  USE_DAG        // <-- specify to evaluate via a DAG of the function
#define SAVE_RESULTS   // <-- specify whether to save results to file
#define  ANALYSE_RATE    // <-- specify whether to analyse rate of convergence
#undef  ANALYSE_TIME    // <-- specify whether to analyse computational time
////////////////////////////////////////////////////////////////////////

#include <fstream>
#include <iomanip>
#include <chrono>

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

#ifdef USE_DAG
 #include "ffunc.hpp"
#endif

using namespace std;
using namespace mc;

////////////////////////////////////////////////////////////////////////

#if defined( TEST_FABS )
const double XL   = -2.;	// <-- X range lower bound
const double XU   =  1.;	// <-- X range upper bound
const double YL   = -1.;	// <-- Y range lower bound
const double YU   =  2.;	// <-- Y range upper bound
template <class T>
T myfunc
( const T&x, const T&y )
{
  //return (-x)+2*y;
  //return log(x+y+4) - exp(x+y);
  //return pow(x+y,2);
  //return x*y;
  return sqrt(fabs(x-y));
}

#elif defined( TEST_EXP )
const double XL   = -2;	// <-- X range lower bound
const double XU   =  1; // <-- X range upper bound
const double YL   = -1; // <-- Y range lower bound
const double YU   =  2;	// <-- Y range upper bound
template <class T>
T myfunc
( const T&x, const T&y )
{
  //return pow(x+exp(x)+2,2);
  //return x*exp(x);
  //return x*(exp(x)-exp(-x))-y*(exp(y)-exp(-y));
  return x*y*(x*(exp(x)-exp(-x))-y*(exp(y)-exp(-y)));
}

#elif defined( TEST_PEAK )
const double XL   = -3.; // <-- X range lower bound
const double XU   =  3.; // <-- X range upper bound
const double YL   = -3.; // <-- Y range lower bound
const double YU   =  3.; // <-- Y range upper bound
template <class T>
T myfunc
( const T&x, const T&y )
{
  //return exp(-pow(x,2)-pow(y,2));
  //return pow((x/5-pow(x,3)-pow(y,5))-exp(-pow(x,2)-pow(y,2)),2);
  //return (x/5-pow(x,3)-pow(y,5))*exp(-pow(x,2)-pow(y,2));
  //return 3*pow(1-x,2)*exp(-pow(x,2)-pow(y+1,2));// - 10*(x/5-pow(x,3)-pow(y,5))*exp(-pow(x,2)-pow(y,2));
  return 3*pow(1-x,2)*exp(-pow(x,2)-pow(y+1,2))-10*(x/5-pow(x,3)-pow(y,5))*exp(-pow(x,2)-pow(y,2))-(1/3)*exp(-pow(x+1,2)-pow(y,2));
}

#elif defined( TEST_EXP1 )
const double XL   = 1.;	// <-- X range lower bound
const double XU   = 2.;	// <-- X range upper bound
const double YL   = 0.;	// <-- Y range lower bound
const double YU   = 1.;	// <-- Y range upper bound
template <class T>
T myfunc
( const T&x, const T&y )
{
  return x*exp(x+sqr(y))-sqr(y);
}

#elif defined( TEST_EXP2 )
const double XL   = -1.;	// <-- X range lower bound
const double XU   =  1.;	// <-- X range upper bound
const double YL   =  0.5;	// <-- Y range lower bound
const double YU   =  2.;	// <-- Y range upper bound
template <class T>
T myfunc
( const T&x, const T&y )
{
  return fabs(x)*exp(-fabs(x)/y);
}

#elif defined( TEST_DIV )
const double XL   = -2.;	// <-- X range lower bound
const double XU   =  0.;	// <-- X range upper bound
const double YL   = -2.;	// <-- Y range lower bound
const double YU   =  0.;	// <-- Y range upper bound
template <class T>
T myfunc
( const T&x, const T&y )
{
  return +1./(pow(x-1.,3)+pow(y-1.,3)+0.1)
         -1./(pow(x-2.,2)+pow(y-3.,4)+0.2)
         +1./(pow(x-3.,3)+pow(y-2.,1)+0.2);
}

#elif defined( TEST_TRIG )
const double XL   = -1.;	// <-- X range lower bound
const double XU   =  1.;	// <-- X range upper bound
const double YL   = -2.;	// <-- Y range lower bound
const double YU   =  2.;	// <-- Y range upper bound
template <class T>
T myfunc
( const T&x, const T&y )
{
  return 1.+x-sin(2.*x+3.*y)-cos(3.*x-5.*y);
}

#elif defined( TEST_RELU )
const double XL   = -3;	// <-- X range lower bound
const double XU   =  3;	// <-- X range upper bound
const double YL   = -3;	// <-- Y range lower bound
const double YU   =  3;	// <-- Y range upper bound
template <class T>
T myfunc
( const T&x, const T&y )
{
  return relu(0.2*x+0.3*y) - relu(0.5*x-0.2*y-3) + relu(0.2*x-0.4*y) + relu(-0.5*x);
}

#elif defined( TEST_TANH )
const double XL   =  -1.;	// <-- X range lower bound
const double XU   =   2.;	// <-- X range upper bound
const double YL   =  -2.;	// <-- Y range lower bound
const double YU   =   1.;	// <-- Y range upper bound
template <class T>
T myfunc
( const T&x, const T&y )
{
  return tanh(x+y);
}

#elif defined( TEST_NORM )
const double XL   = -1.;	// <-- X range lower bound
const double XU   =  1.;	// <-- X range upper bound
const double YL   = -1.;	// <-- Y range lower bound
const double YU   =  1.;	// <-- Y range upper bound
template <class T>
T myfunc
( const T&x, const T&y )
{
  return sqrt(pow(x,2)+pow(y,2));
}
#endif

////////////////////////////////////////////////////////////////////////
int main()
////////////////////////////////////////////////////////////////////////
{
  try{
    MC::options.ENVEL_USE   = 1;
    MC::options.ENVEL_MAXIT = 100;
    MC::options.ENVEL_TOL   = 1e-12;
    MC::options.MVCOMP_USE  = 1;

    PWLSM pwlmod( 2 );
    pwlmod.options.PROD_METH      = PWLSM::Options::PARTIAL;//FULL;//LOG;//NONE;
    pwlmod.options.PROD_CUT       = 0;
    pwlmod.options.REF_WEIGHT     = 0.5;
    pwlmod.options.MAX_SUBDIV     = 0;//8;//16;
    pwlmod.options.USE_SHADOW     = 0;
    pwlmod.options.DISPLAY_SHADOW = 1;
    pwlmod.options.DISPLAY_DIGITS = 10;

    PWCSM pwcmod( 2 );
    pwcmod.options = pwlmod.options;

    // Calculate superposition relaxations
    int const NPWL = 8;	// <-- select initial variable partition >=1
    vector<PWLSV> PWLSVX{ PWLSV( pwlmod, 0, I(XL,XU), NPWL ), PWLSV( pwlmod, 1, I(YL,YU), NPWL ) },
                  PWLSVF( 1 );
    vector<double> DF( 1 );

#ifdef USE_DAG
    // Construct DAG representation of the factorable function
    FFGraph DAG;
    vector<mc::FFVar> X( 2 );
    for( unsigned int i=0; i<2; i++ ) X[i].set( &dag );
    FFVar F = myfunc( X[0], X[1] );
    auto GF = DAG.subgraph( 1, &F );
#ifdef SAVE_RESULTS
    DAG.output( GF );
    ofstream ofdag( "SM-2D.dot", ios_base::out );
    DAG.dot_script( 1, &F, ofdag );
    ofdag.close();
#endif
    DAG.eval( GF, {F}, PWLSVF, X, PWLSVX );

#else
    PWLSVF[0] = myfunc( PWLSVX[0], PWLSVX[1] );
#endif
    cout << PWLSVF[0];

#ifdef SAVE_RESULTS
    ofstream resfile( "SM-2D.out", ios_base::out );
    resfile << scientific << setprecision(5) << right;

    // Repeated calculations at grid points
    int const NPTS = 32;
    for( int iX=0; iX<NPTS; iX++ ){
     for( int iY=0; iY<NPTS; iY++ ){
       vector<double> DX{ XL+iX*(XU-XL)/(NPTS-1.), YL+iY*(YU-YL)/(NPTS-1.) };
#ifdef USE_DAG
       DAG.eval( GF, {F}, DF, X, DX );
#else
       DF[0] = myfunc( DX[0], DX[1] );
#endif

       // Calculate relaxations + propagate all subgradient components
       resfile << setw(14) << DX[0] << setw(14) << DX[1] << setw(14) << DF[0]
               << setw(14) << PWLSVF[0].l() << setw(14) << PWLSVF[0].u()
               << setw(14) << PWLSVF[0].uval({{0,DX[0]},{1,DX[1]}})
               << setw(14) << PWLSVF[0].oval({{0,DX[0]},{1,DX[1]}})
               << endl;
     }
     resfile << endl;
    }

    resfile.close();
#endif

#ifdef ANALYSE_RATE
    ofstream ratefile( "SM-2D_rate.out", ios_base::out );
    ratefile << scientific << setprecision(5) << right;

    vector<I> const XBND0{ I(XL,XU), I(YL,YU) };
#ifdef TEST_PEAK
    map<unsigned,double> XREF{ {0,-1.059997e-02}, {1,1.580344e+00} }; // peak maximum point
    //map<unsigned,double> XREF{ {0,2.288469e-01}, {1,-1.626050e+00} }; // peak minimum point
#else
    map<unsigned,double> XREF{ {0,Op<I>::mid(XBND0[0])}, {1,Op<I>::mid(XBND0[1])} }; // mid-point
#endif
    auto const& red = [=]( const I& bnd, const double& ref, const double& r ){ return r*bnd + (1-r)*ref; };
    auto const& min = [=]( const double& x, const double& y ){ return x<y?x:y; };
    auto const& max = [=]( const double& x, const double& y ){ return x>y?x:y; };

    size_t const NGRID  = 128;
    double const rate   = 0.95;
    double const rhomin = 1e-5;
    DF[0] = myfunc( XREF[0], XREF[1] );
    for( double rho = 1.; rho > rhomin; rho *= rate ){
    //for( double rho = pow(rate,45); rho > rhomin; rho *= rate ){
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

      vector<I> XBND = { red( XBND0[0], XREF[0], rho ), red( XBND0[1], XREF[1], rho ) };
      cout << rho << " " << XBND[0] << " " << XBND[1];

      vector<PWCSV> PWC16SVX(2), PWC16SVF( 1 );
      PWC16SVX[0].set( pwcmod, 0, XBND[0], 16 );
      PWC16SVX[1].set( pwcmod, 1, XBND[1], 16 );
      PWC16SVF[0] = myfunc( PWC16SVX[0], PWC16SVX[1] );

      vector<PWCSV> PWC32SVX(2), PWC32SVF( 1 );
      PWC32SVX[0].set( pwcmod, 0, XBND[0], 32 );
      PWC32SVX[1].set( pwcmod, 1, XBND[1], 32 );
      PWC32SVF[0] = myfunc( PWC32SVX[0], PWC32SVX[1] );
      PWC32SVF[0] = myfunc( PWC32SVX[0], PWC32SVX[1] );

      vector<PWCSV> PWC64SVX(2), PWC64SVF( 1 );
      PWC64SVX[0].set( pwcmod, 0, XBND[0], 64 );
      PWC64SVX[1].set( pwcmod, 1, XBND[1], 64 );
      PWC64SVF[0] = myfunc( PWC64SVX[0], PWC64SVX[1] );

      vector<PWCSV> PWC128SVX(2), PWC128SVF( 1 );
      PWC128SVX[0].set( pwcmod, 0, XBND[0], 128 );
      PWC128SVX[1].set( pwcmod, 1, XBND[1], 128 );
      PWC128SVF[0] = myfunc( PWC128SVX[0], PWC128SVX[1] );

      vector<PWLSV> PWL1SVX(2), PWL1SVF( 1 );
      PWL1SVX[0].set( pwlmod, 0, XBND[0], 1 );
      PWL1SVX[1].set( pwlmod, 1, XBND[1], 1 );
      PWL1SVF[0] = myfunc( PWL1SVX[0], PWL1SVX[1] );
      cout << " PWL1SVX ";

      vector<PWLSV> PWL2SVX(2), PWL2SVF( 1 );
      PWL2SVX[0].set( pwlmod, 0, XBND[0], 2 );
      PWL2SVX[1].set( pwlmod, 1, XBND[1], 2 );
      PWL2SVF[0] = myfunc( PWL2SVX[0], PWL2SVX[1] );
      cout << " PWL2SVX ";

      vector<PWLSV> PWL4SVX(2), PWL4SVF( 1 );
      PWL4SVX[0].set( pwlmod, 0, XBND[0], 4 );
      PWL4SVX[1].set( pwlmod, 1, XBND[1], 4 );
      PWL4SVF[0] = myfunc( PWL4SVX[0], PWL4SVX[1] );
      cout << " PWL4SVX ";

      vector<PWLSV> PWL8SVX(2), PWL8SVF( 1 );
      PWL8SVX[0].set( pwlmod, 0, XBND[0], 8 );
      PWL8SVX[1].set( pwlmod, 1, XBND[1], 8 );
      PWL8SVF[0] = myfunc( PWL8SVX[0], PWL8SVX[1] );
      cout << " PWL8SVX " << endl;

      for( unsigned iX1=0; iX1<NGRID; iX1++ ){
        for( unsigned iX2=0; iX2<NGRID; iX2++ ){
          vector<double> DX{ mc::Op<I>::l(XBND[0])+iX1*mc::Op<I>::diam(XBND[0])/(NGRID-1.),
                                  mc::Op<I>::l(XBND[1])+iX2*mc::Op<I>::diam(XBND[1])/(NGRID-1.) };
          vector<double> DF( 1 );
          DF[0] = myfunc( DX[0], DX[1] );

          min_F = min( min_F, DF[0] );
          max_F = max( max_F, DF[0] );

          vector<MC> MCX{ MC( XBND[0], DX[0] ), MC( XBND[1], DX[1] ) }, MCF( 1 );
          MCF[0] = myfunc( MCX[0], MCX[1] );

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

          double const PWL2SVFu = PWL2SVF[0].uval({{0,DX[0]},{1,DX[1]}});
          double const PWL2SVFo = PWL2SVF[0].oval({{0,DX[0]},{1,DX[1]}});

          distmax_PWL2SVFu = max( distmax_PWL2SVFu, DF[0] - PWL2SVFu );
          distmax_PWL2SVFo = max( distmax_PWL2SVFo, PWL2SVFo - DF[0] );

          double const PWL4SVFu = PWL4SVF[0].uval({{0,DX[0]},{1,DX[1]}});
          double const PWL4SVFo = PWL4SVF[0].oval({{0,DX[0]},{1,DX[1]}});

          distmax_PWL4SVFu = max( distmax_PWL4SVFu, DF[0] - PWL4SVFu );
          distmax_PWL4SVFo = max( distmax_PWL4SVFo, PWL4SVFo - DF[0] );

          double const PWL8SVFu = PWL8SVF[0].uval({{0,DX[0]},{1,DX[1]}});
          double const PWL8SVFo = PWL8SVF[0].oval({{0,DX[0]},{1,DX[1]}});

          distmax_PWL8SVFu = max( distmax_PWL8SVFu, DF[0] - PWL8SVFu );
          distmax_PWL8SVFo = max( distmax_PWL8SVFo, PWL8SVFo - DF[0] );
        }
      }
      
      ratefile << scientific << setprecision(5) << right
               << setw(14) << rho
               << setw(14) << mc::Op<I>::l(XBND[0]) << setw(14) << mc::Op<I>::u(XBND[0])
               << setw(14) << mc::Op<I>::l(XBND[1]) << setw(14) << mc::Op<I>::u(XBND[1])
               //<< setw(14) << min_MCFcv << setw(14) << max_MCFcc
               << setw(14) << max( distmax_MCFcv,            distmax_MCFcc )
               << setw(14) << max( min_F - min_MCFcv,        max_MCFcc - max_F )
               << setw(14) << max( distmax_PWC16SVFu,        distmax_PWC16SVFo )
               << setw(14) << max( min_F - PWC16SVF[0].l(),  PWC16SVF[0].u() - max_F )
               << setw(14) << max( distmax_PWC32SVFu,        distmax_PWC32SVFo )
               << setw(14) << max( min_F - PWC32SVF[0].l(),  PWC32SVF[0].u() - max_F )
               << setw(14) << max( distmax_PWC64SVFu,        distmax_PWC64SVFo )
               << setw(14) << max( min_F - PWC64SVF[0].l(),  PWC64SVF[0].u() - max_F )
               << setw(14) << max( distmax_PWC128SVFu,       distmax_PWC128SVFo )
               << setw(14) << max( min_F - PWC128SVF[0].l(), PWC128SVF[0].u() - max_F )
               << setw(14) << max( distmax_PWL1SVFu,         distmax_PWL1SVFo )
               << setw(14) << max( min_F - PWL1SVF[0].l(),   PWL1SVF[0].u() - max_F )
               << setw(14) << max( distmax_PWL2SVFu,         distmax_PWL2SVFo )
               << setw(14) << max( min_F - PWL2SVF[0].l(),   PWL2SVF[0].u() - max_F )
               << setw(14) << max( distmax_PWL4SVFu,         distmax_PWL4SVFo )
               << setw(14) << max( min_F - PWL4SVF[0].l(),   PWL4SVF[0].u() - max_F )
               << setw(14) << max( distmax_PWL8SVFu,         distmax_PWL8SVFo )
               << setw(14) << max( min_F - PWL8SVF[0].l(),   PWL8SVF[0].u() - max_F )
               << endl;
    }
#endif
  }

#ifdef USE_DAG
  catch( mc::FFBase::Exceptions &eObj ){
    std::cerr << "Error " << eObj.ierr()
              << " in factorable function manipulation:" << std::endl
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
