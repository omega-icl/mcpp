#define TEST_TRIG	// <-- select test function here
#define USE_DAG        // <-- specify whether to use a DAG or operator overloading
#define SAVE_RESULTS   // <-- specify whether to save results to file
const int NX = 50;	// <-- select X discretization here
const int NY = 50;	// <-- select Y discretization here

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
  return sqrt(fabs(x-y));
}

#elif defined( TEST_EXP )
const double XL   = -2.;	// <-- X range lower bound
const double XU   =  1.;	// <-- X range upper bound
const double YL   = -1.;	// <-- Y range lower bound
const double YU   =  2.;	// <-- Y range upper bound
template <class T>
T myfunc
( const T&x, const T&y )
{
  return x*y*(x*(exp(x)-exp(-x))-y*(exp(y)-exp(-y)));
}

#elif defined( TEST_EXP2 )
const double XL   = -1.;	// <-- X range lower bound
const double XU   =  1.;	// <-- X range upper bound
const double YL   =  2.;	// <-- Y range lower bound
const double YU   =  0.5;	// <-- Y range upper bound
template <class T>
T myfunc
( const T&x, const T&y )
{
  return sqrt(fabs(x)*exp(-fabs(x)/y));
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

#elif defined( TEST_DISC )
const double XL   = 0.5;	// <-- X range lower bound
const double XU   = 1.5;	// <-- X range upper bound
const double YL   = 0.5;	// <-- Y range lower bound
const double YU   = 1.5;	// <-- Y range upper bound
template <class T>
T myfunc
( const T&x, const T&y )
{
  return fstep(1-x*y)*(2*(x+y)-exp(x*y+1)-0.5*sin(6*y-1)*sqr(x))+0.5*sin(6*y-1)*sqr(x);
}

#elif defined( TEST_GTCOND )
const double XL   = 0.5;	// <-- X range lower bound
const double XU   = 1.5;	// <-- X range upper bound
const double YL   = 0.5;	// <-- Y range lower bound
const double YU   = 1.5;	// <-- Y range upper bound
template <class T>
T myfunc
( const T&x, const T&y )
{
  return gtcond( 1-x*y, x/sqr(y), pow(x*y,-3) );
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

#elif defined( TEST_LMTD )
const double XL   =  0.1;	// <-- X range lower bound
const double XU   =  2.;	// <-- X range upper bound
const double YL   =  0.1;	// <-- Y range lower bound
const double YU   =  2.;	// <-- Y range upper bound
template <class T>
T myfunc
( const T&x, const T&y )
{
  return lmtd(x,y);
}

#elif defined( TEST_RLMTD )
const double XL   =  0.1;	// <-- X range lower bound
const double XU   =  2.;	// <-- X range upper bound
const double YL   =  0.1;	// <-- Y range lower bound
const double YU   =  2.;	// <-- Y range upper bound
template <class T>
T myfunc
( const T&x, const T&y )
{
  return rlmtd(x,y);
}

#endif

////////////////////////////////////////////////////////////////////////
int main()
////////////////////////////////////////////////////////////////////////
{  
  cout << "INTERVAL LIBRARY: "; 
#ifdef MC__USE_PROFIL
  cout << "PROFIL/BIAS" << endl;
#else
 #ifdef MC__USE_FILIB
  cout << "FILIB++" << endl;
 #else
  #ifdef MC__USE_BOOST
  cout << "BOOST" << endl;
  #else
  cout << "MC++ NON-VERIFIED" << endl;
  #endif
 #endif
#endif

#ifdef SAVE_RESULTS
  ofstream res( "INT-2D.out", ios_base::out );
  res << scientific << setprecision(5) << right;
#endif

  try{
#ifdef USE_DAG
    // Construct DAG representation of the factorable function
    FFGraph DAG;
    FFVar X( &DAG );
    FFVar Y( &DAG );
    FFVar F = myfunc( X, Y );
    auto GF = DAG.subgraph( 1, &F );
#ifdef SAVE_RESULTS
    DAG.output( GF );
    ofstream ofdag( "INT-2D.dot", ios_base::out );
    DAG.dot_script( 1, &F, ofdag );
    ofdag.close();
#endif
#endif

    // Calculate interval bounds
    I IX( XL, XU );
    I IY( YL, YU );
#ifdef USE_DAG
    I IF;
    DAG.eval( GF, 1, &F, &IF, 1, &X, &IX, 1, &Y, &IY );
#else
    IF = myfunc( IX, IY );
#endif
    cout << endl
         << "DOMAIN: " << IX << " x " << IY << endl
         << "BOUNDS: " << IF << endl;

    // Repeated calculations at grid points
    for( int iX=0; iX<NX; iX++ ){ 
     for( int iY=0; iY<NY; iY++ ){ 
      double DX = XL+iX*(XU-XL)/(NX-1.);
      double DY = YL+iY*(YU-YL)/(NY-1.);
#ifdef USE_DAG
      double DF;
      DAG.eval( GF, 1, &F, &DF, 1, &X, &DX, 1, &Y, &DY );
#else
      double DF = myfunc( DX, DY );
#endif
#ifdef SAVE_RESULTS
      res << setw(14) << DX
          << setw(14) << DY
          << setw(14) << DF
          << setw(14) << Op<I>::l(IF)
          << setw(14) << Op<I>::u(IF)
          << endl;
#endif
    }
#ifdef SAVE_RESULTS
    res << endl;
#endif
   }
  }
  
#if !defined(MC__USE_PROFIL) && !defined(MC__USE_FILIB) && !defined(MC__USE_BOOST)
  catch( I::Exceptions &eObj ){
    cerr << "Error " << eObj.ierr()
         << " in natural interval extension:" << endl
	 << eObj.what() << endl
         << "Aborts." << endl;
    return eObj.ierr();
  }
#endif
  catch( FFGraph::Exceptions &eObj ){
    cerr << "Error " << eObj.ierr()
         << " in DAG evaluation:" << endl
	 << eObj.what() << endl
         << "Aborts." << endl;
    return eObj.ierr();
  }
#ifdef SAVE_RESULTS
  res.close();
#endif
  return 0;
}

