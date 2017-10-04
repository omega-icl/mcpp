#define TEST_EXP1        // <-- select test function here
const int NTE = 3;      // <-- select Taylor expansion order here
#undef USE_PROFIL       // <-- specify to use PROFIL for interval arithmetic
#undef USE_FILIB        // <-- specify to use FILIB++ for interval arithmetic
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

#include "scmodel.hpp"
typedef mc::SCModel<I> SCM;
typedef mc::SCVar<I> SCV;

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
const double XL   =  -2.;	// <-- X range lower bound
const double XU   =   0.;	// <-- X range upper bound
const double YL   =  0.;	// <-- Y range lower bound
const double YU   =  2.;	// <-- Y range upper bound
template <class T>
T myfunc
( const T&x, const T&y )
{
  return exp(x*y);
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
#endif

////////////////////////////////////////////////////////////////////////

int main()
{ 
  const unsigned NVAR = 2; 
 {
  CM modCM( NVAR, NTE );
  modCM.options.BOUNDER_TYPE = CM::Options::LSB;
  modCM.options.MIXED_IA = true;

  CV X( &modCM, 0, I(XL,XU) );
  CV Y( &modCM, 1, I(YL,YU) );
  CV F = myfunc( X, Y );
  std::cout << "\nChebyshev model (dense implementation):\n";
  std::cout << F;
 }
 {
  SCM modSCM( NTE );
  modSCM.options.BOUNDER_TYPE = SCM::Options::LSB;
  modSCM.options.MIXED_IA = true;

  SCV X( &modSCM, 0, I(XL,XU) );
  SCV Y( &modSCM, 1, I(YL,YU) );
  SCV F = myfunc( X, Y );
  std::cout << "\nChebyshev model (sparse implementation):\n";
  std::cout << F;
 }
  return 0;
} 

