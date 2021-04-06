#define TEST_TRIG	// <-- select test function here
const int NX = 50;	// <-- select X discretization here
const int NY = 50;	// <-- select Y discretization here
#define SAVE_RESULTS    // <-- specify whether to save results to file
#undef  USE_PROFIL      // <-- specify to use PROFIL for interval arithmetic
#undef  USE_FILIB        // <-- specify to use FILIB++ for interval arithmetic
#undef  USE_BOOST        // <-- specify to use BOOST for interval arithmetic

////////////////////////////////////////////////////////////////////////

#include <fstream>
#include <iomanip>

#ifdef USE_PROFIL
 #include "mcprofil.hpp"
 typedef INTERVAL I;
#else
 #ifdef USE_FILIB
  #include "mcfilib.hpp"
  typedef filib::interval<double> I;
 #else
  #ifdef USE_BOOST
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

#include "ffunc.hpp"

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
#ifdef SAVE_RESULTS
  ofstream res( "INT-2D.out", ios_base::out );
  res << scientific << setprecision(5) << right;
#endif

  try{
    // Construct DAG representation of the factorable function
    FFGraph DAG;
    FFVar X( &DAG );
    FFVar Y( &DAG );
    FFVar F = myfunc( X, Y );
#ifdef SAVE_RESULTS
    DAG.output( DAG.subgraph( 1, &F ) );
    ofstream ofdag( "INT-2D.dot", ios_base::out );
    DAG.dot_script( 1, &F, ofdag );
    ofdag.close();
#endif

    // Calculate interval bounds
    I Xint( XL, XU );
    I Yint( YL, YU );
    I Fint;
    DAG.eval( 1, &F, &Fint, 1, &X, &Xint, 1, &Y, &Yint );
    cout << "Domain: " << Xint << " x " << Yint << endl;
    cout << "Bounds: " << Fint << endl;

    // Repeated calculations at grid points
    for( int iX=0; iX<NX; iX++ ){ 
     for( int iY=0; iY<NY; iY++ ){ 

      double Xval = XL+iX*(XU-XL)/(NX-1.);
      double Yval = YL+iY*(YU-YL)/(NY-1.);
      double Fval;
      DAG.eval( 1, &F, &Fval, 1, &X, &Xval, 1, &Y, &Yval );

#ifdef SAVE_RESULTS
      res << setw(14) << Xval
          << setw(14) << Yval
          << setw(14) << Fval
          << setw(14) << Op<I>::l(Fint)
          << setw(14) << Op<I>::u(Fint)
          << endl;
#endif
    }
#ifdef SAVE_RESULTS
    res << endl;
#endif
   }
  }
#ifndef USE_PROFIL
#ifndef USE_FILIB
#ifndef USE_BOOST
  catch( I::Exceptions &eObj ){
    cerr << "Error " << eObj.ierr()
         << " in natural interval extension:" << endl
	 << eObj.what() << endl
         << "Aborts." << endl;
    return eObj.ierr();
  }
#endif
#endif
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

