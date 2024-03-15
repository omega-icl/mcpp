#define TEST_EXP3	// <-- select test function here
const int NTE = 5;	// <-- select expansion order here
const int NX = 20;	// <-- select X discretization here
const int NY = 20;	// <-- select Y discretization here
#define SAVE_RESULTS    // <-- specify whether to save results to file
#define USE_PROFIL	// <-- specify to use PROFIL for interval arithmetic
#define USE_SPARSE	// <-- specify to use sparse models
#define USE_INTCOEF	// <-- specify to use interval coefficients
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

#ifndef USE_SPARSE
  #include "cmodel.hpp"
  typedef mc::CModel<I> CM;
  typedef mc::CVar<I> CV;
#else
  #ifndef USE_INTCOEF
    #include "scmodel.hpp"
    typedef mc::SCModel<I> CM;
    typedef mc::SCVar<I> CV;
  #else
    #include "sicmodel.hpp"
    typedef mc::SICModel<I> CM;
    typedef mc::SICVar<I> CV; 
  #endif
#endif

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

#elif defined( TEST_EXP0 )
const double XL   =  -2.;	// <-- X range lower bound
const double XU   =   0.;	// <-- X range upper bound
const double Xref =  -1.;	// <-- X ref point for McCormick
const double YL   =   0.;	// <-- Y range lower bound
const double YU   =   2.;	// <-- Y range upper bound
const double Yref =   1.;	// <-- Y ref point for McCormick
template <class T>
T myfunc
( const T&x, const T&y )
{
  return exp(x*y);
}

#elif defined( TEST_EXP1 )
const double XL   =  -1.;	// <-- X range lower bound
const double XU   =   1.;	// <-- X range upper bound
const double Xref =   0.;	// <-- X ref point for McCormick
const double YL   =  -1.;	// <-- Y range lower bound
const double YU   =   1.;	// <-- Y range upper bound
const double Yref =   0.;	// <-- Y ref point for McCormick
template <class T>
T myfunc
( const T&x, const T&y )
{
  return (exp(x)-sqr(y))*x*y;
}

#elif defined( TEST_EXP2 )
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

#elif defined( TEST_EXP3 )
const double XL   = -2.;	// <-- X range lower bound
const double XU   =  2.;	// <-- X range upper bound
const double Xref =  0.;	// <-- X ref point for McCormick
const double YL   = 0.1;	// <-- Y range lower bound
const double YU   =  1.;	// <-- Y range upper bound
const double Yref =  0.;	// <-- Y ref point for McCormick
template <class T>
T myfunc
( const T&x, const T&y )
{
  using mc::inv;
  return inv(y);
  //return x*y*(x*(exp(x)-exp(-x))-y*(exp(y)-exp(-y)));
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

#elif defined( TEST_TRIG3 )
const double XL   =  0.;	// <-- X range lower bound
const double XU   =  10.;	// <-- X range upper bound
const double Xref =  5.;	// <-- X ref point for McCormick
const double YL   =  0.;	// <-- Y range lower bound
const double YU   =  20.;	// <-- Y range upper bound
const double Yref =  10.;	// <-- X ref point for McCormick
template <class T>
T myfunc
( const T&x, const T&y )
{
  return exp( sin(x) + sin(y)*cos(y) );
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

#elif defined( TEST_GP )
const double XL   = -3.;	// <-- X range lower bound
const double XU   =  3.;	// <-- X range upper bound
const double Xref =  1.;	// <-- X ref point for McCormick
const double YL   = -3.;	// <-- Y range lower bound
const double YU   =  3.;	// <-- Y range upper bound
const double Yref =  1.;	// <-- Y ref point for McCormick
template <class T>
T myfunc
( const T&x, const T&y )
{
  return 1.86360571672641*(0.14220168012508*(1 + 2.23606797749979*sqrt(
     64.2558879895505*sqr(0.194030125231034 - 0.166666666666667*x) + 
     0.451453304154821*sqr((-0.404841797459896) - 0.166666666666667*y)) + 
     1.66666666666667*(64.2558879895505*sqr(0.194030125231034 - 
     0.166666666666667*x) + 0.451453304154821*sqr((-0.404841797459896) - 
     0.166666666666667*y)))*exp(-2.23606797749979*sqrt(64.2558879895505*sqr(
     0.194030125231034 - 0.166666666666667*x) + 0.451453304154821*sqr((-
     0.404841797459896) - 0.166666666666667*y))) + 0.0468045532937668*(1 + 
     2.23606797749979*sqrt(64.2558879895505*sqr(0.268324435572052 - 
     0.166666666666667*x) + 0.451453304154821*sqr((-0.220830234095203) - 
     0.166666666666667*y)) + 1.66666666666667*(64.2558879895505*sqr(
     0.268324435572052 - 0.166666666666667*x) + 0.451453304154821*sqr((-
     0.220830234095203) - 0.166666666666667*y)))*exp(-2.23606797749979*sqrt(
     64.2558879895505*sqr(0.268324435572052 - 0.166666666666667*x) + 
     0.451453304154821*sqr((-0.220830234095203) - 0.166666666666667*y))) - 
     0.436827958615331*(1 + 2.23606797749979*sqrt(64.2558879895505*sqr((-
     0.0505822037784427) - 0.166666666666667*x) + 0.451453304154821*sqr(
     0.190399166851636 - 0.166666666666667*y)) + 1.66666666666667*(
     64.2558879895505*sqr((-0.0505822037784427) - 0.166666666666667*x) + 
     0.451453304154821*sqr(0.190399166851636 - 0.166666666666667*y)))*exp(-
     2.23606797749979*sqrt(64.2558879895505*sqr((-0.0505822037784427) - 
     0.166666666666667*x) + 0.451453304154821*sqr(0.190399166851636 - 
     0.166666666666667*y))) + 0.0405780260973915*(1 + 2.23606797749979*sqrt(
     64.2558879895505*sqr(0.335653234585858 - 0.166666666666667*x) + 
     0.451453304154821*sqr(0.319274707641782 - 0.166666666666667*y)) + 
     1.66666666666667*(64.2558879895505*sqr(0.335653234585858 - 
     0.166666666666667*x) + 0.451453304154821*sqr(0.319274707641782 - 
     0.166666666666667*y)))*exp(-2.23606797749979*sqrt(64.2558879895505*sqr(
     0.335653234585858 - 0.166666666666667*x) + 0.451453304154821*sqr(
     0.319274707641782 - 0.166666666666667*y))) + 0.0383717635737363*(1 + 
     2.23606797749979*sqrt(64.2558879895505*sqr(0.416115262788179 - 
     0.166666666666667*x) + 0.451453304154821*sqr((-0.0829399801155143) - 
     0.166666666666667*y)) + 1.66666666666667*(64.2558879895505*sqr(
     0.416115262788179 - 0.166666666666667*x) + 0.451453304154821*sqr((-
     0.0829399801155143) - 0.166666666666667*y)))*exp(-2.23606797749979*sqrt(
     64.2558879895505*sqr(0.416115262788179 - 0.166666666666667*x) + 
     0.451453304154821*sqr((-0.0829399801155143) - 0.166666666666667*y))) + 
     0.226121620987036*(1 + 2.23606797749979*sqrt(64.2558879895505*sqr((-
     0.324344059222606) - 0.166666666666667*x) + 0.451453304154821*sqr(
     0.0775173811415889 - 0.166666666666667*y)) + 1.66666666666667*(
     64.2558879895505*sqr((-0.324344059222606) - 0.166666666666667*x) + 
     0.451453304154821*sqr(0.0775173811415889 - 0.166666666666667*y)))*exp(-
     2.23606797749979*sqrt(64.2558879895505*sqr((-0.324344059222606) - 
     0.166666666666667*x) + 0.451453304154821*sqr(0.0775173811415889 - 
     0.166666666666667*y))) + 0.055301550115541*(1 + 2.23606797749979*sqrt(
     64.2558879895505*sqr((-0.289414019566429) - 0.166666666666667*x) + 
     0.451453304154821*sqr((-0.337636568589209) - 0.166666666666667*y)) + 
     1.66666666666667*(64.2558879895505*sqr((-0.289414019566429) - 
     0.166666666666667*x) + 0.451453304154821*sqr((-0.337636568589209) - 
     0.166666666666667*y)))*exp(-2.23606797749979*sqrt(64.2558879895505*sqr((-
     0.289414019566429) - 0.166666666666667*x) + 0.451453304154821*sqr((-
     0.337636568589209) - 0.166666666666667*y))) - 0.220359443423157*(1 + 
     2.23606797749979*sqrt(64.2558879895505*sqr((-0.120888900038076) - 
     0.166666666666667*x) + 0.451453304154821*sqr((-0.125103261913789) - 
     0.166666666666667*y)) + 1.66666666666667*(64.2558879895505*sqr((-
     0.120888900038076) - 0.166666666666667*x) + 0.451453304154821*sqr((-
     0.125103261913789) - 0.166666666666667*y)))*exp(-2.23606797749979*sqrt(
     64.2558879895505*sqr((-0.120888900038076) - 0.166666666666667*x) + 
     0.451453304154821*sqr((-0.125103261913789) - 0.166666666666667*y))) + 
     0.0439612252066741*(1 + 2.23606797749979*sqrt(64.2558879895505*sqr(
     0.0826904704977199 - 0.166666666666667*x) + 0.451453304154821*sqr(
     0.458984274331446 - 0.166666666666667*y)) + 1.66666666666667*(
     64.2558879895505*sqr(0.0826904704977199 - 0.166666666666667*x) + 
     0.451453304154821*sqr(0.458984274331446 - 0.166666666666667*y)))*exp(-
     2.23606797749979*sqrt(64.2558879895505*sqr(0.0826904704977199 - 
     0.166666666666667*x) + 0.451453304154821*sqr(0.458984274331446 - 
     0.166666666666667*y))) + 0.0560980731767207*(1 + 2.23606797749979*sqrt(
     64.2558879895505*sqr((-0.435192808594306) - 0.166666666666667*x) + 
     0.451453304154821*sqr(0.220400711043696 - 0.166666666666667*y)) + 
     1.66666666666667*(64.2558879895505*sqr((-0.435192808594306) - 
     0.166666666666667*x) + 0.451453304154821*sqr(0.220400711043696 - 
     0.166666666666667*y)))*exp(-2.23606797749979*sqrt(64.2558879895505*sqr((-
     0.435192808594306) - 0.166666666666667*x) + 0.451453304154821*sqr(
     0.220400711043696 - 0.166666666666667*y))));
}

#endif

////////////////////////////////////////////////////////////////////////

int main()
{

#ifdef SAVE_RESULTS
  ofstream res( "CM-2D.out", ios_base::out );
  res << std::scientific << std::setprecision(5) << std::right;
#endif

  try{ 

    // Define Chebyshev model environment
#ifndef USE_SPARSE
    CM mod( 2, NTE );
#else
  #ifndef USE_INTCOEF
    CM mod( NTE );
    mod.options.BASIS        = 1; //0;//1;
    mod.options.REMEZ_USE    = true;//false;
  #else
    CM mod( NTE );
    mod.options.BASIS        = 1; //0;//1;
    mod.options.REMEZ_USE    = true;//false;
    mod.options.HOT_SPLIT    = CM::Options::FULL;//NONE;
    mod.options.REMEZ_USE    = true;//false; 
  #endif
#endif

    // <-- set options here -->
    mod.options.BOUNDER_TYPE = CM::Options::LSB;//NAIVE;
    mod.options.MIXED_IA     = true;//false;
    
    // Define variables X and Y, and evaluate Chebyshev model
    CV CVX( &mod, 0, I(XL,XU) );
    CV CVY( &mod, 1, I(YL,YU) );
    CV CVF = myfunc( CVX, CVY );
    std::cout << "\nChebyshev model of f(x,y):" << CVF;

    // Repeated calculations at grid points (for display)
    for( int iX=0; iX<NX; iX++ ){ 
      for( int iY=0; iY<NY; iY++ ){ 

        double DXY[2] = { XL+iX*(XU-XL)/(NX-1.), YL+iY*(YU-YL)/(NY-1.) };
        double DF = myfunc( DXY[0], DXY[1] );
#if defined (USE_SPARSE) && defined (USE_INTCOEF)
        I IF = CVF.P( std::map<unsigned,double>({std::make_pair(0,DXY[0]),std::make_pair(1,DXY[1])}) );
        double PF = Op<I>::mid( IF );
#elif defined (USE_SPARSE)
        double PF = CVF.P( std::map<unsigned,double>({std::make_pair(0,DXY[0]),std::make_pair(1,DXY[1])}) );
        I IF = PF + CVF.R();
#else
        double PF = CVF.P( DXY );
        I IF = PF + CVF.R();
#endif

#ifdef SAVE_RESULTS
        res << std::setw(14) << DXY[0] << std::setw(14) << DXY[1] << std::setw(14) << DF
            << std::setw(14) << PF << std::setw(14) << Op<I>::l(IF)  << setw(14) << Op<I>::u(IF)
            << std::setw(14) << Op<I>::l(CVF.B())  << std::setw(14) << Op<I>::u(CVF.B())
            << std::endl;
#endif
      }
      res << endl;
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
