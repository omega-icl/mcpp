#define TEST_FABS       // <-- select test function here
const int NX = 500;	    // <-- select discretization here
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
const double XL   = -2.;	// <-- range lower bound
const double XU   =  1.;	// <-- range upper bound
template <class T>
T myfunc
( const T&x )
{
  return -x*fabs(x);
}

#elif defined( TEST_SQRT )
const double XL   = -2.;	// <-- range lower bound
const double XU   =  1.;	// <-- range upper bound
template <class T>
T myfunc
( const T&x )
{
  return sqrt(fabs(x));
}

#elif defined( TEST_EXP )
const double XL   = -2.;	// <-- range lower bound
const double XU   =  1.;	// <-- range upper bound
template <class T>
T myfunc
( const T&x )
{
  return x*exp(-pow(x,2));
}

#elif defined( TEST_EXP2 )
const double XL   = -2.;	// <-- range lower bound
const double XU   =  1.;	// <-- range upper bound
template <class T>
T myfunc
( const T&x )
{
  return x*exp(-fabs(pow(x,3)));
}

#elif defined( TEST_ERF )
const double XL   = -10.;	// <-- range lower bound
const double XU   =  10.;	// <-- range upper bound
template <class T>
T myfunc
( const T&x )
{
  return erf(x);
}

#elif defined( TEST_ERFC )
const double XL   = -2.;	// <-- range lower bound
const double XU   =  3.;	// <-- range upper bound
template <class T>
T myfunc
( const T&x )
{
  return x*erfc(-pow(x,2));
}

#elif defined( TEST_HYP )
const double XL   = -3;	    // <-- range lower bound
const double XU   =  2.;	// <-- range upper bound
template <class T>
T myfunc
( const T&x )
{
  return tanh(x);
}

#elif defined( TEST_HYP2 )
const double XL   = -3;	    // <-- range lower bound
const double XU   =  2.5;	// <-- range upper bound
template <class T>
T myfunc
( const T&x )
{
  return sqr(sinh(x))-sqr(cosh(x));
}

#elif defined( TEST_TRIG )
const double XL   = -4*PI;	// <-- range lower bound
const double XU   =  PI/3.;	// <-- range upper bound
template <class T>
T myfunc
( const T&x )
{
  return tan(cos(x*atan(x)));
}

#elif defined( TEST_TRIG2 )
const double XL   =  PI/6.;	// <-- range lower bound
const double XU   =  PI/3.;	// <-- range upper bound
template <class T>
T myfunc
( const T&x )
{
  return cos(pow(x,2))*sin(pow(x,-3));
}

#elif defined( TEST_XLOG )
const double XL   =  0.1;	// <-- range lower bound
const double XU   =  0.6;	// <-- range upper bound
template <class T>
T myfunc
( const T&x )
{
  return pow(-fabs(x-0.5)-xlog(x),3);
  //return xlog(x);
}

#elif defined( TEST_MIN1 )
const double XL   = -2;	// <-- range lower bound
const double XU   =  1;	// <-- range upper bound
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
template <class T>
T myfunc
( const T&x )
{
  return fstep(4-x)*(fstep(x-3)*(exp(4-x)+3-(fstep(3-x)*(-sqr(x-2.5)+4)))
                   +(fstep(3-x)*(-sqr(x-2.5)+4))-(2*x-7))+(2*x-7);
}

#elif defined( TEST_CHEB )
const double XL   = -1.;	// <-- range lower bound
const double XU   =  1.;	// <-- range upper bound
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
template <class T>
T myfunc
( const T&x )
{
  T m[3] = { x, x+.5, x-.5 };
  return prod( 3, m );
}

#elif defined( TEST_HULL )
const double XL   =  -1.;	// <-- range lower bound
const double XU   =   .5;	// <-- range upper bound
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
#ifdef SAVE_RESULTS
  ofstream res( "INT-1D.out", ios_base::out );
  res << scientific << setprecision(5) << right;
#endif

  try{ 
    // Construct DAG representation of the factorable function
    FFGraph DAG;
    FFVar X( &DAG );
    FFVar F = myfunc( X );
#ifdef SAVE_RESULTS
    DAG.output( DAG.subgraph( 1, &F ) );
    ofstream ofdag( "INT-1D.dot", ios_base::out );
    DAG.dot_script( 1, &F, ofdag );
    ofdag.close();
#endif

    // Calculate interval bounds
    I Xint( XL, XU );
    I Fint;
    DAG.eval( 1, &F, &Fint, 1, &X, &Xint );
    cout << "Domain: " << Xint << endl;
    cout << "Bounds: " << Fint << endl;

    // Repeated calculations at grid points
    for( int iX=0; iX<NX; iX++ ){ 

      double Xval = XL+iX*(XU-XL)/(NX-1.);
      double Fval;
      DAG.eval( 1, &F, &Fval, 1, &X, &Xval );

#ifdef SAVE_RESULTS
      res << setw(14) << Xval
          << setw(14) << Fval
          << setw(14) << Op<I>::l(Fint)
          << setw(14) << Op<I>::u(Fint)
          << endl;
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

