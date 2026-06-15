#define TEST_TANH	// <-- select test function here
#undef  USE_DAG        // <-- specify to evaluate via a DAG of the function
#define SAVE_RESULTS   // <-- specify whether to save results to file
////////////////////////////////////////////////////////////////////////

#include <fstream>
#include <iomanip>
#include <chrono>

#include "mcboost.hpp"
typedef boost::numeric::interval_lib::save_state<boost::numeric::interval_lib::rounded_transc_opp<double>> T_boost_round;
typedef boost::numeric::interval_lib::checking_base<double> T_boost_check;
typedef boost::numeric::interval_lib::policies<T_boost_round,T_boost_check> T_boost_policy;
typedef boost::numeric::interval<double,T_boost_policy> I;

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

#if defined( TEST_TANH )
const double XL   = -2.;	// <-- X range lower bound
const double XU   =  1.;	// <-- X range upper bound
template <class T>
T myfunc
( const T&x )
{
  return tanh(pow(exp(x),2)-pow(x,3));
}

#elif defined( TEST_ATAN )
const double XL   = -2.;	// <-- X range lower bound
const double XU   =  1.;	// <-- X range upper bound
template <class T>
T myfunc
( const T&x )
{
  //return tanh(x);
  //return atan(x);
  //return erf(x);
  return erfc(x);
}

#elif defined( TEST_ACOS )
const double XL   = PI/2.;	// <-- X range lower bound
const double XU   =  15.*PI/2.;	// <-- X range upper bound
template <class T>
T myfunc
( const T&x )
{
  return cos(x);
  //return pow(x,2);
  //return pow(x,3);
  //return erf(x);
  //return erfc(x);
  //return acos(x);
  //return asin(x);
  //return acos( erf(x) );
}

#elif defined( TEST_CUBIC )
const double XL   = -2.;	// <-- X range lower bound
const double XU   =  1.;	// <-- X range upper bound
template <class T>
T myfunc
( const T&x )
{
  return pow(x,3);
}

#endif

////////////////////////////////////////////////////////////////////////
int main()
////////////////////////////////////////////////////////////////////////
{
  size_t NX = 1;

  try{
    // Calculate superposition relaxations w/ continuous piecewise linear univariates
    PWLSM pwlmod( 1 );
    pwlmod.options.PROD_METH      = PWLSM::Options::PARTIAL;//FULL;//LOG;//NONE;
    pwlmod.options.PROD_CUT       = 0;
    pwlmod.options.REF_WEIGHT     = 0.5;
    pwlmod.options.MAX_SUBDIV     = 0;//8;//16;
    pwlmod.options.USE_CVXCCV     = 1;
    pwlmod.options.USE_SHADOW     = 0;
    pwlmod.options.DISPLAY_SHADOW = 1;
    pwlmod.options.DISPLAY_DIGITS = 10;

    int const NPWL = 16;	// <-- select initial variable partition >=1
    vector<PWLSV> PWLSVX{ PWLSV( pwlmod, 0, I(XL,XU), NPWL ) },
                  PWLSVF( NX );

    // Calculate superposition relaxations w/ piecewise constant univariates
    PWCSM pwcmod( NX );
    pwcmod.options = pwlmod.options;
    mc::PWCU::options.SLOPEUSE = 2;

    int const NPWC = 16;	// <-- select variable partition
    vector<PWCSV> PWCSVX{ PWCSV( pwcmod, 0, I(XL,XU), NPWC ) },
                  PWCSVF( NX );

#ifdef USE_DAG
    // Construct DAG representation of the factorable function
    FFGraph DAG;
    vector<mc::FFVar> X( NX );
    for( unsigned int i=0; i<NX; i++ ) X[i].set( &dag );
    FFVar F = myfunc( X[0] );
    auto GF = DAG.subgraph( 1, &F );
#ifdef SAVE_RESULTS
    DAG.output( GF );
    ofstream ofdag( "SM-2D.dot", ios_base::out );
    DAG.dot_script( 1, &F, ofdag );
    ofdag.close();
#endif
    DAG.eval( GF, {F}, PWLSVF, X, PWLSVX );
    DAG.eval( GF, {F}, PWCSVF, X, PWCSVX );
#else
    PWLSVF[0] = myfunc( PWLSVX[0] );
    PWCSVF[0] = myfunc( PWCSVX[0] );
#endif
    cout << PWLSVF[0];
    cout << PWCSVF[0];

#ifdef SAVE_RESULTS
    ofstream resfile( "SM-1D.out", ios_base::out );
    resfile << scientific << setprecision(5) << right;

    // Repeated calculations at grid points
    int const NPTS = 200;
    std::vector<double> DF( NX );
    for( int iX=0; iX<NPTS; iX++ ){
      vector<double> DX{ XL+iX*(XU-XL)/(NPTS-1.) };
#ifdef USE_DAG
      DAG.eval( GF, {F}, DF, X, DX );
#else
      DF[0] = myfunc( DX[0] );
#endif

      // Output relaxations
      resfile << setw(14) << DX[0] << setw(14) << DF[0]
              << setw(14) << PWLSVF[0].l() << setw(14) << PWLSVF[0].u()
              << setw(14) << PWLSVF[0].uval({{0,DX[0]}})
              << setw(14) << PWLSVF[0].oval({{0,DX[0]}})
              << setw(14) << PWCSVF[0].l() << setw(14) << PWCSVF[0].u()
              << setw(14) << PWCSVF[0].uval({{0,DX[0]}})
              << setw(14) << PWCSVF[0].oval({{0,DX[0]}})
              << endl;
    }
    resfile.close();
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
