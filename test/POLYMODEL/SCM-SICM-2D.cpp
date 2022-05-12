#define TEST_EXP2            // <-- select test function here
const int NTE = 5;           // <-- select Taylor expansion order here
const unsigned NREP = 2000;   // <-- select repetition for acurate timing 
const int NX = 25;	      // <-- select X discretization here
const int NY = 25;	      // <-- select Y discretization here
#define SAVE_RESULTS         // <-- specify whether to save results to file
#define USE_PROFIL           // <-- specify to use PROFIL for interval arithmetic
#undef USE_FILIB             // <-- specify to use FILIB++ for interval arithmetic
////////////////////////////////////////////////////////////////////////
//#define MC__SCMODEL_DEBUG_SPROD
//#define MC__SICMODEL_DEBUG_SPROD
//#define MC__SICMODEL_DEBUG_COMPOSITION
//#define MC__SCMODEL_DEBUG_COMPOSITION
//#define MC__SICMODEL_DEBUG_SLIFT
//#define MC__SICMODEL_DEBUG_SPROD

#include <fstream>
#include <iomanip>

#include "mctime.hpp"
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

#include "scmodel.hpp"
typedef mc::SCModel<I> SCM;
typedef mc::SCVar<I> SCV;

#include "sicmodel.hpp"
typedef mc::SICModel<I> SICM;
typedef mc::SICVar<I> SICV;

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
const double YL   =  1.;	// <-- Y range lower bound
const double YU   =  4.;	// <-- Y range upper bound
template <class T>
T myfunc
( const T&x, const T&y )
{
  return pow(x,y)/(0.5+pow(x,y));
  //return 1./(1.+0.5*pow(x,-y));
}

#elif defined( TEST_EXP0 )
const double XL   =  -1.;	// <-- X range lower bound
const double XU   =  1.;	// <-- X range upper bound
const double YL   =  -1.;	// <-- Y range lower bound
const double YU   =  1.;	// <-- Y range upper bound
template <class T>
T myfunc
( const T&x, const T&y )
{
  //return exp(sqrt(x+y));
  return cos(exp(x+y));
}

#elif defined( TEST_EXP1 )
const double XL   =  1.;	// <-- X range lower bound
const double XU   =  2.;	// <-- X range upper bound
const double YL   =  0.;	// <-- Y range lower bound
const double YU   =  1.;	// <-- Y range upper bound
template <class T>
T myfunc
( const T&x, const T&y )
{
  //return exp(x)*exp(y);
  return x*exp(pow(y,2))-pow(y,2);
}

#elif defined( TEST_EXP2 )
const double XL   = -2.;	// <-- X range lower bound
const double XU   =  2.;	// <-- X range upper bound
const double YL   = -1.;	// <-- Y range lower bound
const double YU   =  1.;	// <-- Y range upper bound
template <class T>
T myfunc
( const T&x, const T&y )
{
  return x*y*(x*(exp(x)-exp(-x))-y*(exp(y)-exp(-y)));
}

#elif defined( TEST_INV )
const double XL   = -2.;	// <-- X range lower bound
const double XU   =  0.;	// <-- X range upper bound
const double YL   =  1.;	// <-- Y range lower bound
const double YU   =  3.;	// <-- Y range upper bound
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
double YL   = -1.;	// <-- Y range lower bound
double YU   =  0.;	// <-- Y range upper bound
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
const double YL   = -0.7;	// <-- Y range lower bound
const double YU   =  0.7;	// <-- Y range upper bound
template <class T>
T myfunc
( const T&x, const T&y )
{
  return 1.+x-sin(2.*x+3.*y)-cos(3.*x-5.*y);
}

#elif defined( TEST_TRIG2 )
const double XL   = -0.6;	// <-- X range lower bound
const double XU   =  0.6;	// <-- X range upper bound
const double YL   = -0.6;	// <-- Y range lower bound
const double YU   =  0.6;	// <-- Y range upper bound
template <class T>
T myfunc
( const T&x, const T&y )
{
  return exp(-sqr(x)-sqr(y)) * cos(2*PI*(sqr(x)+sqr(y)));
}

#elif defined( TEST_NORM )
const double XL   =  0.5;	// <-- X range lower bound
const double XU   =  2.;	// <-- X range upper bound
const double YL   =  0.5;	// <-- Y range lower bound
const double YU   =  2.;	// <-- Y range upper bound
template <class T>
T myfunc
( const T&x, const T&y )
{
  return sqrt(pow(x,2)+pow(y,2));
}

#elif defined( TEST_GP )
const double XL   = -1.;	// <-- X range lower bound
const double XU   =  1.;	// <-- X range upper bound
const double YL   = -1.;	// <-- Y range lower bound
const double YU   =  1.;	// <-- Y range upper bound
template <class T>
T myfunc
( const T&x, const T&y )
{
  return (-2.23606797749979*exp(64.2558879895505*sqr(
     0.194030125231034 - 0.166666666666667*x) + 0.451453304154821*sqr((-
     0.404841797459896) - 0.166666666666667*y)));
//  return 0.14220168012508*(1 + 2.23606797749979*sqrt(
//     64.2558879895505*sqr(0.194030125231034 - 0.166666666666667*x) + 
//     0.451453304154821*sqr((-0.404841797459896) - 0.166666666666667*y)) + 
//     1.66666666666667*(64.2558879895505*sqr(0.194030125231034 - 
//     0.166666666666667*x) + 0.451453304154821*sqr((-0.404841797459896) - 
//     0.166666666666667*y)));
//  return 0.14220168012508*(1 + 2.23606797749979*sqrt(
//     64.2558879895505*sqr(0.194030125231034 - 0.166666666666667*x) + 
//     0.451453304154821*sqr((-0.404841797459896) - 0.166666666666667*y)) + 
//     1.66666666666667*(64.2558879895505*sqr(0.194030125231034 - 
//     0.166666666666667*x) + 0.451453304154821*sqr((-0.404841797459896) - 
//     0.166666666666667*y)))*exp(-2.23606797749979*sqrt(64.2558879895505*sqr(
//     0.194030125231034 - 0.166666666666667*x) + 0.451453304154821*sqr((-
//     0.404841797459896) - 0.166666666666667*y)));
//  return 1.86360571672641*(0.14220168012508*(1 + 2.23606797749979*sqrt(
//     64.2558879895505*sqr(0.194030125231034 - 0.166666666666667*x) + 
//     0.451453304154821*sqr((-0.404841797459896) - 0.166666666666667*y)) + 
//     1.66666666666667*(64.2558879895505*sqr(0.194030125231034 - 
//     0.166666666666667*x) + 0.451453304154821*sqr((-0.404841797459896) - 
//     0.166666666666667*y)))*exp(-2.23606797749979*sqrt(64.2558879895505*sqr(
//     0.194030125231034 - 0.166666666666667*x) + 0.451453304154821*sqr((-
//     0.404841797459896) - 0.166666666666667*y))) + 0.0468045532937668*(1 + 
//     2.23606797749979*sqrt(64.2558879895505*sqr(0.268324435572052 - 
//     0.166666666666667*x) + 0.451453304154821*sqr((-0.220830234095203) - 
//     0.166666666666667*y)) + 1.66666666666667*(64.2558879895505*sqr(
//     0.268324435572052 - 0.166666666666667*x) + 0.451453304154821*sqr((-
//     0.220830234095203) - 0.166666666666667*y)))*exp(-2.23606797749979*sqrt(
//     64.2558879895505*sqr(0.268324435572052 - 0.166666666666667*x) + 
//     0.451453304154821*sqr((-0.220830234095203) - 0.166666666666667*y))) - 
//     0.436827958615331*(1 + 2.23606797749979*sqrt(64.2558879895505*sqr((-
//     0.0505822037784427) - 0.166666666666667*x) + 0.451453304154821*sqr(
//     0.190399166851636 - 0.166666666666667*y)) + 1.66666666666667*(
//     64.2558879895505*sqr((-0.0505822037784427) - 0.166666666666667*x) + 
//     0.451453304154821*sqr(0.190399166851636 - 0.166666666666667*y)))*exp(-
//     2.23606797749979*sqrt(64.2558879895505*sqr((-0.0505822037784427) - 
//     0.166666666666667*x) + 0.451453304154821*sqr(0.190399166851636 - 
//     0.166666666666667*y))) + 0.0405780260973915*(1 + 2.23606797749979*sqrt(
//     64.2558879895505*sqr(0.335653234585858 - 0.166666666666667*x) + 
//     0.451453304154821*sqr(0.319274707641782 - 0.166666666666667*y)) + 
//     1.66666666666667*(64.2558879895505*sqr(0.335653234585858 - 
//     0.166666666666667*x) + 0.451453304154821*sqr(0.319274707641782 - 
//     0.166666666666667*y)))*exp(-2.23606797749979*sqrt(64.2558879895505*sqr(
//     0.335653234585858 - 0.166666666666667*x) + 0.451453304154821*sqr(
//     0.319274707641782 - 0.166666666666667*y))) + 0.0383717635737363*(1 + 
//     2.23606797749979*sqrt(64.2558879895505*sqr(0.416115262788179 - 
//     0.166666666666667*x) + 0.451453304154821*sqr((-0.0829399801155143) - 
//     0.166666666666667*y)) + 1.66666666666667*(64.2558879895505*sqr(
//     0.416115262788179 - 0.166666666666667*x) + 0.451453304154821*sqr((-
//     0.0829399801155143) - 0.166666666666667*y)))*exp(-2.23606797749979*sqrt(
//     64.2558879895505*sqr(0.416115262788179 - 0.166666666666667*x) + 
//     0.451453304154821*sqr((-0.0829399801155143) - 0.166666666666667*y))) + 
//     0.226121620987036*(1 + 2.23606797749979*sqrt(64.2558879895505*sqr((-
//     0.324344059222606) - 0.166666666666667*x) + 0.451453304154821*sqr(
//     0.0775173811415889 - 0.166666666666667*y)) + 1.66666666666667*(
//     64.2558879895505*sqr((-0.324344059222606) - 0.166666666666667*x) + 
//     0.451453304154821*sqr(0.0775173811415889 - 0.166666666666667*y)))*exp(-
//     2.23606797749979*sqrt(64.2558879895505*sqr((-0.324344059222606) - 
//     0.166666666666667*x) + 0.451453304154821*sqr(0.0775173811415889 - 
//     0.166666666666667*y))) + 0.055301550115541*(1 + 2.23606797749979*sqrt(
//     64.2558879895505*sqr((-0.289414019566429) - 0.166666666666667*x) + 
//     0.451453304154821*sqr((-0.337636568589209) - 0.166666666666667*y)) + 
//     1.66666666666667*(64.2558879895505*sqr((-0.289414019566429) - 
//     0.166666666666667*x) + 0.451453304154821*sqr((-0.337636568589209) - 
//     0.166666666666667*y)))*exp(-2.23606797749979*sqrt(64.2558879895505*sqr((-
//     0.289414019566429) - 0.166666666666667*x) + 0.451453304154821*sqr((-
//     0.337636568589209) - 0.166666666666667*y))) - 0.220359443423157*(1 + 
//     2.23606797749979*sqrt(64.2558879895505*sqr((-0.120888900038076) - 
//     0.166666666666667*x) + 0.451453304154821*sqr((-0.125103261913789) - 
//     0.166666666666667*y)) + 1.66666666666667*(64.2558879895505*sqr((-
//     0.120888900038076) - 0.166666666666667*x) + 0.451453304154821*sqr((-
//     0.125103261913789) - 0.166666666666667*y)))*exp(-2.23606797749979*sqrt(
//     64.2558879895505*sqr((-0.120888900038076) - 0.166666666666667*x) + 
//     0.451453304154821*sqr((-0.125103261913789) - 0.166666666666667*y))) + 
//     0.0439612252066741*(1 + 2.23606797749979*sqrt(64.2558879895505*sqr(
//     0.0826904704977199 - 0.166666666666667*x) + 0.451453304154821*sqr(
//     0.458984274331446 - 0.166666666666667*y)) + 1.66666666666667*(
//     64.2558879895505*sqr(0.0826904704977199 - 0.166666666666667*x) + 
//     0.451453304154821*sqr(0.458984274331446 - 0.166666666666667*y)))*exp(-
//     2.23606797749979*sqrt(64.2558879895505*sqr(0.0826904704977199 - 
//     0.166666666666667*x) + 0.451453304154821*sqr(0.458984274331446 - 
//     0.166666666666667*y))) + 0.0560980731767207*(1 + 2.23606797749979*sqrt(
//     64.2558879895505*sqr((-0.435192808594306) - 0.166666666666667*x) + 
//     0.451453304154821*sqr(0.220400711043696 - 0.166666666666667*y)) + 
//     1.66666666666667*(64.2558879895505*sqr((-0.435192808594306) - 
//     0.166666666666667*x) + 0.451453304154821*sqr(0.220400711043696 - 
//     0.166666666666667*y)))*exp(-2.23606797749979*sqrt(64.2558879895505*sqr((-
//     0.435192808594306) - 0.166666666666667*x) + 0.451453304154821*sqr(
//     0.220400711043696 - 0.166666666666667*y))));
}
#endif

////////////////////////////////////////////////////////////////////////

int main()
{ 
 {
  CM modCM( 2, NTE );
  modCM.options.BOUNDER_TYPE = CM::Options::LSB;//NAIVE;//LSB;
  modCM.options.MIXED_IA = true;//false;

  CV X( &modCM, 0, I(XL,XU) );
  CV Y( &modCM, 1, I(YL,YU) );
  CV F = myfunc( X, Y );
  std::cout << F;
  double tStart = mc::userclock();
  for( unsigned i=0; i<NREP; i++ )
    F = myfunc( X, Y );
  std::cout << "\nChebyshev model (dense implementation):" << (mc::userclock()-tStart)/(double)NREP << " CPU-sec\n";
 }
 {
  SCM modSCM( NTE );
  modSCM.options.REMEZ_USE = true;//false;
  modSCM.options.BOUNDER_TYPE = SCM::Options::LSB; //NAIVE;
  modSCM.options.MIXED_IA = false;//true;//false;

  SCV X( &modSCM, 0, I(XL,XU) );
  SCV Y( &modSCM, 1, I(YL,YU) );
  SCV F = myfunc( X, Y );
  std::cout << F;
  double tStart = mc::userclock();
  for( unsigned i=0; i<NREP; i++ )
    F = myfunc( X, Y );
  std::cout << "\nChebyshev model (sparse implementation):" << (mc::userclock()-tStart)/(double)NREP << " CPU-sec\n";

#ifdef SAVE_RESULTS
  ofstream res( "SCM-2D.out", ios_base::out );
  res << std::scientific << std::setprecision(5) << std::right;
  // Repeated calculations at grid points (for display)
  for( int iX=0; iX<NX; iX++ ){ 
    for( int iY=0; iY<NY; iY++ ){ 
      double DXY[2] = { XL+iX*(XU-XL)/(NX-1.), YL+iY*(YU-YL)/(NY-1.) };
      double DF = myfunc( DXY[0], DXY[1] );
      double PF = F.P( DXY );
      I IF = PF + F.R();
      res << std::setw(14) << DXY[0] << std::setw(14) << DXY[1]
          << std::setw(14) << DF << std::setw(14) << PF
          << std::setw(14) << Op<I>::l(IF) << std::setw(14) << Op<I>::u(IF)
          << std::endl;
    }
    res << endl;
  }
#endif
 }
 {
  SICM modSICM( NTE );
  
  modSICM.options.REMEZ_USE    = true;
  modSICM.options.INTERP_EXTRA = 0;//2*NTE;
  modSICM.options.BOUNDER_TYPE = SICM::Options::LSB; //NAIVE;
  modSICM.options.HOT_SPLIT    = SICM::Options::FULL;
  modSICM.options.MIXED_IA     = false;//true;//false;

  SICV X( &modSICM, 0, I(XL,XU) );
  SICV Y( &modSICM, 1, I(YL,YU) );
  SICV F = myfunc( X, Y );
  std::cout << F;
  double tStart = mc::userclock();
  for( unsigned i=0; i<NREP; i++ )
    F = myfunc( X, Y );
  std::cout << "\nChebyshev model (sparse interval implementation):" << (mc::userclock()-tStart)/(double)NREP << " CPU-sec\n";

#ifdef SAVE_RESULTS
  ofstream res( "SICM-2D.out", ios_base::out );
  res << std::scientific << std::setprecision(5) << std::right;
  // Repeated calculations at grid points (for display)
  for( int iX=0; iX<NX; iX++ ){ 
    for( int iY=0; iY<NY; iY++ ){ 
      double DXY[2] = { XL+iX*(XU-XL)/(NX-1.), YL+iY*(YU-YL)/(NY-1.) };
      double DF = myfunc( DXY[0], DXY[1] );
      double PF = F.P( DXY );
      I IF = F.IP( DXY );
      res << std::setw(14) << DXY[0] << std::setw(14) << DXY[1]
          << std::setw(14) << DF << std::setw(14) << PF
          << std::setw(14) << Op<I>::l(IF) << std::setw(14) << Op<I>::u(IF)
          << std::endl;
    }
    res << endl;
  }
#endif
 }

  return 0;
} 

