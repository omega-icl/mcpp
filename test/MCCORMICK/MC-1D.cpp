#define TEST_EXP2       // <-- select test function here
const int NX = 500;	    // <-- select discretization here
#define SAVE_RESULTS    // <-- specify whether to save results to file
#undef USE_PROFIL      // <-- specify to use PROFIL for interval arithmetic
#define USE_FILIB        // <-- specify to use FILIB++ for interval arithmetic
#define USE_DAG          // <-- specify to evaluate via a DAG of the function

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
    #include "interval.hpp"
    typedef mc::Interval I;
  #endif
#endif

#include "mccormick.hpp"
typedef mc::McCormick<I> MC;

using namespace std;
using namespace mc;

////////////////////////////////////////////////////////////////////////

#if defined( TEST_FABS )
const double XL   = -2.;	// <-- range lower bound
const double XU   =  1.;	// <-- range upper bound
const double Xref =  0.;	// <-- linearization point
template <class T>
T myfunc
( const T&x )
{
  return -x*fabs(x);
}

#elif defined( TEST_SQRT )
const double XL   = -2.;	// <-- range lower bound
const double XU   =  1.;	// <-- range upper bound
const double Xref =  0.;	// <-- linearization point
template <class T>
T myfunc
( const T&x )
{
  return sqrt(fabs(x));
}

#elif defined( TEST_EXP )
const double XL   = -2.;	// <-- range lower bound
const double XU   =  1.;	// <-- range upper bound
const double Xref =  0.;	// <-- linearization point
template <class T>
T myfunc
( const T&x )
{
  return x*exp(-pow(x,2));
}

#elif defined( TEST_EXP2 )
const double XL   = -2.;	// <-- range lower bound
const double XU   =  1.;	// <-- range upper bound
const double Xref = -1.;	// <-- linearization point
template <class T>
T myfunc
( const T&x )
{
  return x*exp(-fabs(pow(x,3)));
}

#elif defined( TEST_ERF )
const double XL   = -10.;	// <-- range lower bound
const double XU   =  10.;	// <-- range upper bound
const double Xref = -2.;	// <-- linearization point
template <class T>
T myfunc
( const T&x )
{
  return erf(x);
}

#elif defined( TEST_ERFC )
const double XL   = -2.;	// <-- range lower bound
const double XU   =  3.;	// <-- range upper bound
const double Xref =  0.;	// <-- linearization point
template <class T>
T myfunc
( const T&x )
{
  return x*erfc(-pow(x,2));
}

#elif defined( TEST_HYP )
const double XL   = -3;	    // <-- range lower bound
const double XU   =  2.;	// <-- range upper bound
const double Xref =  0.;	// <-- linearization point
template <class T>
T myfunc
( const T&x )
{
  return tanh(x);
}

#elif defined( TEST_HYP2 )
const double XL   = -3;	    // <-- range lower bound
const double XU   =  2.5;	// <-- range upper bound
const double Xref =  0.;	// <-- linearization point
template <class T>
T myfunc
( const T&x )
{
  return sqr(sinh(x))-sqr(cosh(x));
}

#elif defined( TEST_TRIG )
const double XL   = -4*PI;	// <-- range lower bound
const double XU   =  PI/3.;	// <-- range upper bound
const double Xref =  0.;	// <-- linearization point
template <class T>
T myfunc
( const T&x )
{
  return tan(cos(x*atan(x)));
}

#elif defined( TEST_TRIG2 )
const double XL   =  PI/6.;	// <-- range lower bound
const double XU   =  PI/3.;	// <-- range upper bound
const double Xref =  PI/4.;	// <-- linearization point
template <class T>
T myfunc
( const T&x )
{
  return cos(pow(x,2))*sin(pow(x,-3));
}

#elif defined( TEST_XLOG )
const double XL   =  0.1;	// <-- range lower bound
const double XU   =  0.9;	// <-- range upper bound
const double Xref =  0.5;	// <-- linearization point
template <class T>
T myfunc
( const T&x )
{
  return pow(-fabs(x-0.5)-xlog(x),3);
}

#elif defined( TEST_MIN )
const double XL   = -.8;	// <-- range lower bound
const double XU   =  .8;	// <-- range upper bound
const double Xref =  0.;	// <-- linearization point
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
const double Xref =  3.;	// <-- linearization point
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
const double Xref =  .5;	// <-- linearization point
template <class T>
T myfunc
( const T&x )
{
  return ltcond( sqr(x)-0.25, pow(x-0.5,3), pow(-x-0.25,3) );
}

#elif defined( TEST_GTCOND )
const double XL   =  0.;	// <-- range lower bound
const double XU   =  1.;	// <-- range upper bound
const double Xref =  .5;	// <-- linearization point
template <class T>
T myfunc
( const T&x )
{
  return gtcond( sqr(x)-.25, sqr(x)-0.8, 0.1-sqr(x) );
}

#elif defined( TEST_CHEB )
const double XL   = -1.;	// <-- range lower bound
const double XU   =  1.;	// <-- range upper bound
const double Xref =  .5;	// <-- linearization point
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
const double Xref =  .5;	// <-- linearization point
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
const double Xref =  .5;	// <-- linearization point
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
const double Xref =  -.5;	// <-- linearization point
bool hull( const double&y, const double&z )
{ return 0.5*(y+z); }
template <class T>
T myfunc
( const T&x )
{
  return hull( sqrt(fabs(x-.5)), sqrt(fabs(x+.5)) );
}

#endif

#ifndef USE_DAG
////////////////////////////////////////////////////////////////////////
int main()
////////////////////////////////////////////////////////////////////////
{  

  // <-- set options here -->
  MC::options.ENVEL_USE   = true;
  MC::options.ENVEL_MAXIT = 100;
  MC::options.ENVEL_TOL   = 1e-12;
  MC::options.MVCOMP_USE  = true;

#ifdef SAVE_RESULTS
  ofstream res( "MC-1D.out", ios_base::out );
  res << scientific << setprecision(5) << right;
#endif

  try{ 

    // Calculate relaxations & subgradient at point Xref
    MC Xrel( I(XL,XU), Xref );
    Xrel.sub(1,0);
    MC Zref = myfunc( Xrel );
    cout << "Relaxation at reference point:\n" << Zref << endl;

    // Repeated calculations at grid points
    for( int iX=0; iX<NX; iX++ ){ 

      double Xval = XL+iX*(XU-XL)/(NX-1.);
      double Zval = myfunc( Xval );

      MC Xrel( I(XL,XU), Xval );

      // Calculate relaxations + propagate subgradient component
      Xrel.sub(1,0);
      MC Zrel = myfunc( Xrel );

#ifdef SAVE_RESULTS
      res << setw(14) << Xval          << setw(14) << Zval
          << setw(14) << Zrel.l()      << setw(14) << Zrel.u()
          << setw(14) << Zrel.cv()     << setw(14) << Zrel.cc()
          << setw(14) << Zrel.cvsub(0) << setw(14) << Zrel.ccsub(0)
          << setw(14) << Zref.cv()+Zref.cvsub(0)*(Xval-Xref)
          << setw(14) << Zref.cc()+Zref.ccsub(0)*(Xval-Xref)
          << endl;
#endif
    }

  }
#ifndef USE_PROFIL
#ifndef USE_FILIB
  catch( I::Exceptions &eObj ){
    cerr << "Error " << eObj.ierr()
         << " in natural interval extension:" << endl
	 << eObj.what() << endl
         << "Aborts." << endl;
    return eObj.ierr();
  }
#endif
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

#else

#include "ffunc.hpp"

////////////////////////////////////////////////////////////////////////
int main()
////////////////////////////////////////////////////////////////////////
{  

  // <-- set options here -->
  MC::options.ENVEL_USE   = true;
  MC::options.ENVEL_MAXIT = 100;
  MC::options.ENVEL_TOL   = 1e-12;
  MC::options.MVCOMP_USE  = true;

#ifdef SAVE_RESULTS
  ofstream res( "MC-1D.out", ios_base::out );
  res << scientific << setprecision(5) << right;
#endif

  try{ 
    // Construct DAG representation of the factorable function
    FFGraph DAG;
    FFVar X( &DAG );
    FFVar Z = myfunc( X );
#ifdef SAVE_RESULTS
    DAG.output( DAG.subgraph( 1, &Z ) );
    ofstream ofdag( "MC-1D.dot", ios_base::out );
    DAG.dot_script( 1, &Z, ofdag );
    ofdag.close();
#endif

    // Calculate relaxations & subgradient at point Xref
    MC Xrel( I(XL,XU), Xref );
    Xrel.sub(1,0); 
    MC Zref;
    DAG.eval( 1, &Z, &Zref, 1, &X, &Xrel );
    cout << "Relaxation at reference point:\n" << Zref << endl;

    // Repeated calculations at grid points
    for( int iX=0; iX<NX; iX++ ){ 

      double Xval = XL+iX*(XU-XL)/(NX-1.);
      double Zval = myfunc( Xval );

      MC Xrel( I(XL,XU), Xval );

      // Calculate relaxations + propagate subgradient component
      Xrel.sub(1,0);
      MC Zrel;
      DAG.eval( 1, &Z, &Zrel, 1, &X, &Xrel );

#ifdef SAVE_RESULTS
      res << setw(14) << Xval          << setw(14) << Zval
          << setw(14) << Zrel.l()      << setw(14) << Zrel.u()
          << setw(14) << Zrel.cv()     << setw(14) << Zrel.cc()
          << setw(14) << Zrel.cvsub(0) << setw(14) << Zrel.ccsub(0)
          << setw(14) << Zref.cv()+Zref.cvsub(0)*(Xval-Xref)
          << setw(14) << Zref.cc()+Zref.ccsub(0)*(Xval-Xref)
          << endl;
#endif
    }

  }
#ifndef USE_PROFIL
#ifndef USE_FILIB
  catch( I::Exceptions &eObj ){
    cerr << "Error " << eObj.ierr()
         << " in natural interval extension:" << endl
	 << eObj.what() << endl
         << "Aborts." << endl;
    return eObj.ierr();
  }
#endif
#endif
  catch( MC::Exceptions &eObj ){
    cerr << "Error " << eObj.ierr()
         << " in McCormick relaxation:" << endl
	 << eObj.what() << endl
         << "Aborts." << endl;
    return eObj.ierr();
  }
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
#endif
