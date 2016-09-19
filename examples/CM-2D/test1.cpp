#define TEST_EXP	// <-- select test function here
const int NTE = 5;	// <-- select Taylor expansion order here
const int NX = 50;	// <-- select X discretization here
const int NY = 50;	// <-- select Y discretization here
#define SAVE_RESULTS    // <-- specify whether to save results to file
#undef USE_PROFIL	// <-- specify to use PROFIL for interval arithmetic
#undef USE_FILIB	// <-- specify to use FILIB++ for interval arithmetic
#define USE_SPARSE	// <-- specify to use sparse models
#undef  MC__CVAR_SPARSE_PRODUCT_NAIVE
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

#include "cmodel.hpp"
typedef mc::CModel<I> CM;
typedef mc::CVar<I> CV;

using namespace std;
using namespace mc;

////////////////////////////////////////////////////////////////////////

#if defined( TEST_MULT )
const double XL   = -1.;	// <-- X range lower bound
const double XU   =  2.;	// <-- X range upper bound
const double YL   = -2.;	// <-- Y range lower bound
const double YU   =  1.;	// <-- Y range upper bound
template <class T>
T myfunc
( const T&x, const T&y )
{
  return -pow(x+y,3)*x*y;
  //return -(sqr(x+y)*(x+y))*x*y;
}

#elif defined( TEST_HILL )
const double XL   =  0.1;	// <-- X range lower bound
const double XU   =  1.;	// <-- X range upper bound
const double Xref =  5.;	// <-- X ref point for McCormick
const double YL   =  1.;	// <-- Y range lower bound
const double YU   =  4.;	// <-- Y range upper bound
const double Yref =  2.5;	// <-- Y ref point for McCormick
template <class T>
T myfunc
( const T&x, const T&y )
{
  return pow(x,y)/(0.5+pow(x,y));
  //return 1./(1.+0.5*pow(x,-y));
}

#elif defined( TEST_EXP )
const double XL   =  -2.;	// <-- X range lower bound
const double XU   =   0.;	// <-- X range upper bound
const double Xref =  -1.;	// <-- X ref point for McCormick
const double YL   =  0.;	// <-- Y range lower bound
const double YU   =  2.;	// <-- Y range upper bound
const double Yref =  1.;	// <-- Y ref point for McCormick
template <class T>
T myfunc
( const T&x, const T&y )
{
  return exp(x*y);
}

#elif defined( TEST_EXP1 )
const double XL   =  1.;	// <-- X range lower bound
const double XU   =  2.;	// <-- X range upper bound
const double Xref =  1.5;	// <-- X ref point for McCormick
const double YL   =  0.;	// <-- Y range lower bound
const double YU   =  1.;	// <-- Y range upper bound
const double Yref =  0.5;	// <-- Y ref point for McCormick
template <class T>
T myfunc
( const T&x, const T&y )
{
  return x*exp(pow(y,2))-pow(y,2);
}

#elif defined( TEST_EXP2 )
const double XL   = -2.;	// <-- X range lower bound
const double XU   =  2.;	// <-- X range upper bound
const double Xref =  0.;	// <-- X ref point for McCormick
const double YL   = -1.;	// <-- Y range lower bound
const double YU   =  1.;	// <-- Y range upper bound
const double Yref =  0.;	// <-- Y ref point for McCormick
template <class T>
T myfunc
( const T&x, const T&y )
{
  return x*y*(x*(exp(x)-exp(-x))-y*(exp(y)-exp(-y)));
}

#elif defined( TEST_INV )
const double XL   = -2.;	// <-- X range lower bound
const double XU   =  0.;	// <-- X range upper bound
const double Xref = -1.;	// <-- X ref point for McCormick
const double YL   =  1.;	// <-- Y range lower bound
const double YU   =  3.;	// <-- Y range upper bound
const double Yref =  2.;	// <-- Y ref point for McCormick
template <class T>
T myfunc
( const T&x, const T&y )
{
  return -1./(pow(x-4.,2)+pow(y-4.,2)+0.1)
         -1./(pow(x-1.,2)+pow(y-1.,2)+0.2)
         -1./(pow(x-8.,2)+pow(y-8.,2)+0.2);
}

#elif defined( TEST_INV2 )
double XL   = -2.;	// <-- X range lower bound
double XU   =  0.;	// <-- X range upper bound
double Xref = -2.;	// <-- X ref point for McCormick
double YL   = -1.;	// <-- Y range lower bound
double YU   =  0.;	// <-- Y range upper bound
double Yref = -1.;	// <-- Y ref point for McCormick
template <class T>
T myfunc
( const T&x, const T&y )
{
  return 1./(((x-1.)*(x-1.)*(x-1.)+(y-1.)*(y-1.)*(y-1.)+0.1));
//  return +1./(pow(x-1.,3)+pow(y-1.,3)+0.1);
//  return +1./(pow(x-1.,3)+pow(y-1.,3)+0.1)
//         -1./(pow(x-2.,2)+pow(y-3.,4)+0.2)
//         +1./(pow(x-3.,3)+pow(y-2.,1)+0.2);
}

#elif defined( TEST_TRIG )
const double XL   = -0.7;	// <-- X range lower bound
const double XU   =  0.7;	// <-- X range upper bound
const double Xref =  0.;	// <-- X ref point for McCormick
const double YL   = -0.7;	// <-- Y range lower bound
const double YU   =  0.7;	// <-- Y range upper bound
const double Yref =  0.;	// <-- Y ref point for McCormick
template <class T>
T myfunc
( const T&x, const T&y )
{
  return 1.+x-sin(2.*x+3.*y)-cos(3.*x-5.*y);
}

#elif defined( TEST_TRIG2 )
const double XL   = -0.6;	// <-- X range lower bound
const double XU   =  0.6;	// <-- X range upper bound
const double Xref =  0.;	// <-- X ref point for McCormick
const double YL   = -0.6;	// <-- Y range lower bound
const double YU   =  0.6;	// <-- Y range upper bound
const double Yref =  0.;	// <-- Y ref point for McCormick
template <class T>
T myfunc
( const T&x, const T&y )
{
  return exp(-sqr(x)-sqr(y)) * cos(2*PI*(sqr(x)+sqr(y)));
}

#elif defined( TEST_NORM )
const double XL   =  0.5;	// <-- X range lower bound
const double XU   =  2.;	// <-- X range upper bound
const double Xref =  1.;	// <-- X ref point for McCormick
const double YL   =  0.5;	// <-- Y range lower bound
const double YU   =  2.;	// <-- Y range upper bound
const double Yref =  1.;	// <-- Y ref point for McCormick
template <class T>
T myfunc
( const T&x, const T&y )
{
  return sqrt(pow(x,2)+pow(y,2));
}
#endif

////////////////////////////////////////////////////////////////////////

int main()
{

#ifdef SAVE_RESULTS
  ofstream res( "test1.out", ios_base::out );
  res << std::scientific << std::setprecision(5) << std::right;
#endif

  try{ 

    // Define Chebyshev model environment
#ifdef USE_SPARSE
    CM mod( 2, NTE, true );
#else
    CM mod( 2, NTE );
#endif

    // <-- set options here -->
    mod.options.BOUNDER_TYPE = CM::Options::NAIVE;//LSB;//BERNSTEIN;
    //mod.options.BOUNDER_ORDER = 20;
    mod.options.MIXED_IA = true;
    
    // Define variables X and Y, and evaluate Chebyshev model
    CV CVX( &mod, 0, I(XL,XU) );
    CV CVY( &mod, 1, I(YL,YU) );
    CV CVF = myfunc( CVX, CVY );
    std::cout << "\nChebyshev model of f(x,y):" << CVF;
    //std::cout << "\nBernstein bounder:" << CVF.bound(CM::Options::BERNSTEIN) << std::endl;
/*
    XL = XL+0.25*(XU-XL); YU = YU-0.5*(YU-YL);
    I IXY[2] = { I(XL,XU), I(YL,YU) };
    CVF = CVF.scale( IXY );
    std::cout << CVF;

    CVX.set( &mod, 0, I(XL,XU) );
    CVY.set( &mod, 1, I(YL,YU) );
    CVF = myfunc( CVX, CVY );
    std::cout << "\nChebyshev model of f(x,y):" << CVF;
*/
    // Repeated calculations at grid points (for display)
    for( int iX=0; iX<NX; iX++ ){ 
      for( int iY=0; iY<NY; iY++ ){ 

        double DXY[2] = { XL+iX*(XU-XL)/(NX-1.), YL+iY*(YU-YL)/(NY-1.) };
        double DF = myfunc( DXY[0], DXY[1] );
        I IF = CVF.P( DXY ) + CVF.R();

#ifdef SAVE_RESULTS
        res << std::setw(14) << DXY[0] << std::setw(14) << DXY[1] << std::setw(14) << DF
            << std::setw(14) << CVF.P( DXY ) << std::setw(14) << IF.l()  << setw(14) << IF.u()
            << std::setw(14) << CVF.B().l()  << std::setw(14) << CVF.B().u()
            << std::endl;
#endif
      }
      res << endl;
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
  catch( CM::Exceptions &eObj ){
    cerr << "Error " << eObj.ierr()
         << " in Chebyshev model computation:" << endl
	 << eObj.what() << endl
         << "Aborts." << endl;
    return eObj.ierr();
  }

#ifdef SAVE_RESULTS
  res.close();
#endif
  return 0;
}
