#define TEST_DISC       // <-- select test function here
#define USE_DAG         // <-- specify to evaluate via a DAG of the function
#define SAVE_RESULTS    // <-- specify whether to save results to file
const int NX = 500;	 // <-- select discretization here
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
const double XL   = -2.;	// <-- range lower bound
const double XU   =  1.;	// <-- range upper bound
const double XREF =  0.;	// <-- linearization point
template <class T>
T myfunc
( const T&x )
{
  return -x*fabs(x);
}

#elif defined( TEST_SQRT )
const double XL   = -2.;	// <-- range lower bound
const double XU   =  1.;	// <-- range upper bound
const double XREF =  0.;	// <-- linearization point
template <class T>
T myfunc
( const T&x )
{
  return sqrt(fabs(x));
}

#elif defined( TEST_EXP )
const double XL   = -2.;	// <-- range lower bound
const double XU   =  1.;	// <-- range upper bound
const double XREF =  0.;	// <-- linearization point
template <class T>
T myfunc
( const T&x )
{
  return x*exp(-pow(x,2));
}

#elif defined( TEST_EXP2 )
const double XL   = -2.;	// <-- range lower bound
const double XU   =  1.;	// <-- range upper bound
const double XREF = -1.;	// <-- linearization point
template <class T>
T myfunc
( const T&x )
{
  return x*exp(-fabs(pow(x,3)));
}

#elif defined( TEST_ERF )
const double XL   = -10.;	// <-- range lower bound
const double XU   =  10.;	// <-- range upper bound
const double XREF = -2.;	// <-- linearization point
template <class T>
T myfunc
( const T&x )
{
  return erf(x);
}

#elif defined( TEST_ERFC )
const double XL   = -2.;	// <-- range lower bound
const double XU   =  3.;	// <-- range upper bound
const double XREF =  0.;	// <-- linearization point
template <class T>
T myfunc
( const T&x )
{
  return x*erfc(-pow(x,2));
}

#elif defined( TEST_HYP )
const double XL   = -3;	    // <-- range lower bound
const double XU   =  2.;	// <-- range upper bound
const double XREF =  0.;	// <-- linearization point
template <class T>
T myfunc
( const T&x )
{
  return tanh(x);
}

#elif defined( TEST_HYP2 )
const double XL   = -3;	    // <-- range lower bound
const double XU   =  2.5;	// <-- range upper bound
const double XREF =  0.;	// <-- linearization point
template <class T>
T myfunc
( const T&x )
{
  return sqr(sinh(x))-sqr(cosh(x));
}

#elif defined( TEST_TRIG )
const double XL   = -4*PI;	// <-- range lower bound
const double XU   =  PI/3.;	// <-- range upper bound
const double XREF =  0.;	// <-- linearization point
template <class T>
T myfunc
( const T&x )
{
  return tan(cos(x*atan(x)));
}

#elif defined( TEST_TRIG2 )
const double XL   =  PI/6.;	// <-- range lower bound
const double XU   =  PI/3.;	// <-- range upper bound
const double XREF =  PI/4.;	// <-- linearization point
template <class T>
T myfunc
( const T&x )
{
  return cos(pow(x,2))*sin(pow(x,-3));
}

#elif defined( TEST_XLOG )
const double XL   =  0.1;	// <-- range lower bound
const double XU   =  0.9;	// <-- range upper bound
const double XREF =  0.5;	// <-- linearization point
template <class T>
T myfunc
( const T&x )
{
  return pow(-fabs(x-0.5)-xlog(x),3);
}

#elif defined( TEST_MIN1 )
const double XL   = -2;	// <-- range lower bound
const double XU   =  1;	// <-- range upper bound
const double XREF =  0.;	// <-- linearization point
template <class T>
T myfunc
( const T&x )
{
  //T m[2] = { -x, x };
  //T m[2] = { pow(x-1,2), pow(x+1,2) };
  //return min(2,m);
  return max( pow(x-1,2), pow(x+1,2) );
}

#elif defined( TEST_MIN2 )
const double XL   = -.8;	// <-- range lower bound
const double XU   =  .8;	// <-- range upper bound
const double XREF =  0.;	// <-- linearization point
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
const double XREF =  3.;	// <-- linearization point
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
const double XREF =  .5;	// <-- linearization point
template <class T>
T myfunc
( const T&x )
{
  return ltcond( sqr(x)-0.25, pow(x-0.5,3), pow(-x-0.25,3) );
}

#elif defined( TEST_GTCOND )
const double XL   =  0.;	// <-- range lower bound
const double XU   =  1.;	// <-- range upper bound
const double XREF =  .5;	// <-- linearization point
template <class T>
T myfunc
( const T&x )
{
  return gtcond( sqr(x)-.25, sqr(x)-0.8, 0.1-sqr(x) );
}

#elif defined( TEST_CHEB )
const double XL   = -1.;	// <-- range lower bound
const double XU   =  1.;	// <-- range upper bound
const double XREF =  .5;	// <-- linearization point
template <class T>
T myfunc
( const T&x )
{
  return cheb(x,2)+0.5*cheb(x,3)-0.3*cheb(x,4);
  //return cheb(x,2)+0.5*cheb(x,3)-0.3*cheb(x,4);
}

#elif defined( TEST_PROD )
const double XL   = -1.;	// <-- range lower bound
const double XU   =  1.;	// <-- range upper bound
const double XREF =  .5;	// <-- linearization point
template <class T>
T myfunc
( const T&x )
{
  T m[3] = { x, x+.5, x-.5 };
  return prod( 3, m );
}

#elif defined( TEST_INTER )
const double XL   =  0.;	// <-- range lower bound
const double XU   =  1.;	// <-- range upper bound
const double XREF =  .5;	// <-- linearization point
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
const double XREF =  -.5;	// <-- linearization point
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
  ofstream res( "MC-1D.out", ios_base::out );
  res << scientific << setprecision(5) << right;
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
    FFVar F = myfunc( X );
    auto GF = DAG.subgraph( 1, &F );
#ifdef SAVE_RESULTS
    DAG.output( GF );
    ofstream ofdag( "MC-1D.dot", ios_base::out );
    DAG.dot_script( 1, &F, ofdag );
    ofdag.close();
#endif
#endif

    // Calculate relaxations & subgradient at point XREF
    MC MCX( I(XL,XU), XREF );
    MCX.sub(1,0);
#ifdef USE_DAG
    MC MCFREF;
    DAG.eval( GF, 1, &F, &MCFREF, 1, &X, &MCX );
#else
    MC MCFREF = myfunc( MCX );
#endif
    cout << "RELAXATION AT REFERENCE POINT:\n" << MCFREF << endl;

    // Repeated calculations at grid points
#ifdef SAVE_RESULTS
    for( int iX=0; iX<NX; iX++ ){ 
      // Calculate original function
      double DX = XL+iX*(XU-XL)/(NX-1.);
#ifdef USE_DAG
      double DF;
      DAG.eval( GF, 1, &F, &DF, 1, &X, &DX );
#else
      double DF = myfunc( DX );
#endif
      // Calculate relaxations + propagate subgradient component
      MCX.c( DX );
      MCX.sub( 1, 0 );
#ifdef USE_DAG
      MC MCF;
      DAG.eval( GF, 1, &F, &MCF, 1, &X, &MCX );
#else
      MC MCF = myfunc( MCX );
#endif
      res << setw(14) << DX           << setw(14) << DF
          << setw(14) << MCF.l()      << setw(14) << MCF.u()
          << setw(14) << MCF.cv()     << setw(14) << MCF.cc()
          << setw(14) << MCF.cvsub(0) << setw(14) << MCF.ccsub(0)
          << setw(14) << MCFREF.cv()+MCFREF.cvsub(0)*(DX-XREF)
          << setw(14) << MCFREF.cc()+MCFREF.ccsub(0)*(DX-XREF)
          << endl;
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
  res.close();
#endif
  return 0;
}

