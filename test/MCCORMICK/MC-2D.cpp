#define TEST_EXP	// <-- select test function here
#define USE_DAG        // <-- specify to evaluate via a DAG of the function
#define SAVE_RESULTS   // <-- specify whether to save results to file
const int NX = 32;	// <-- select X discretization here
const int NY = 32;	// <-- select Y discretization here
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

#include "mccormick.hpp"
typedef mc::McCormick<I> MC;

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

#elif defined( TEST_PEAK )
const double XL   = -3.;	// <-- X range lower bound
const double XU   =  3.;	// <-- X range upper bound
const double YL   = -3.;	// <-- Y range lower bound
const double YU   =  3.;	// <-- Y range upper bound
template <class T>
T myfunc
( const T&x, const T&y )
{
  return 3*pow(1-x,2)*exp(-pow(x,2)-pow(y+1,2))-10*(x/5-pow(x,3)-pow(y,5))*exp(-pow(x,2)-pow(y,2))-(1/3)*exp(-pow(x+1,2)-pow(y,2));
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

#elif defined( TEST_RELU )
const double XL   = -3;	// <-- X range lower bound
const double XU   =  3;	// <-- X range upper bound
const double YL   = -3; // <-- Y range lower bound
const double YU   =  3;	// <-- Y range upper bound
template <typename T>
T relu
( T const& x )
{
  return Op<T>::max( x, 0. );
}
template <class T>
T myfunc
( const T&x, const T&y )
{
  //return relu(x+y);
  //return - relu(0.5*x-0.2*y);
  return relu(0.2*x+0.3*y) - relu(0.5*x-0.2*y-3) + relu(0.2*x-0.4*y) + relu(-0.5*x);
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
  ofstream allsub( "MC-2D.out", ios_base::out );
  ofstream dirsub( "MC2-2D.out", ios_base::out );
  allsub << scientific << setprecision(5) << right;
  dirsub << scientific << setprecision(5) << right;
#endif

  // <-- set options here -->
  MC::options.ENVEL_USE   = true;
  MC::options.ENVEL_MAXIT = 100;
  MC::options.ENVEL_TOL   = 1e-12;
  MC::options.MVCOMP_USE  = true;

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
    ofstream ofdag( "MC-2D.dot", ios_base::out );
    DAG.dot_script( 1, &F, ofdag );
    ofdag.close();
#endif
#endif

#ifdef SAVE_RESULTS
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

       // Calculate relaxations + propagate all subgradient components
       MC MCX( I(XL,XU), DX );
       MC MCY( I(YL,YU), DY );
       MCX.sub( 2, 0 );
       MCY.sub( 2, 1 );
#ifdef USE_DAG
       MC MCF;
       DAG.eval( GF, 1, &F, &MCF, 1, &X, &MCX, 1, &Y, &MCY );
#else
       MC MCF = myfunc( MCX, MCY );
#endif
       if( !iX && !iY ) std::cout << MCF << std::endl;

       allsub << setw(14) << DX << setw(14) << DY << setw(14) << DF
              << setw(14) << MCF.l() << setw(14) <<  MCF.u()
              << setw(14) << MCF.cv() << setw(14) << MCF.cc()
              << setw(14) << MCF.cvsub(0) << setw(14) << MCF.ccsub(0)
              << setw(14) << MCF.cvsub(1) << setw(14) << MCF.ccsub(1)
              << endl;

       // Calculate relaxations + propagate directional subgradient
       const double dir[2] = { 1, -1 };
       MCX.sub( 1, &dir[0], &dir[0] );
       MCY.sub( 1, &dir[1], &dir[1] );
#ifdef USE_DAG
       DAG.eval( GF, 1, &F, &MCF, 1, &X, &MCX, 1, &Y, &MCY );
#else
       MCF = myfunc( MCX, MCY );
#endif
       dirsub << setw(14) << DX << setw(14) << DY << setw(14) << DF
              << setw(14) << MCF.l() << setw(14) <<  MCF.u()
              << setw(14) << MCF.cv() << setw(14) << MCF.cc()
              << setw(14) << MCF.cvsub(0) << setw(14) << MCF.ccsub(0)
              << endl;
      }
      allsub << endl;
      dirsub << endl;
    }
#endif
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
  catch( MC::Exceptions &eObj ){
    cerr << "Error " << eObj.ierr()
         << " in McCormick relaxation:" << endl
	 << eObj.what() << endl
         << "Aborts." << endl;
    return eObj.ierr();
  }

#ifdef SAVE_RESULTS
  allsub.close();
  dirsub.close();
#endif
  return 0;
}
