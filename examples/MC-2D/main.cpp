#define TEST_TRIG	// <-- select test function here
const int NX = 50;	// <-- select X discretization here
const int NY = 50;	// <-- select Y discretization here
#define SAVE_RESULTS    // <-- specify whether to save results to file
#undef USE_PROFIL	// <-- specify to use PROFIL for interval arithmetic
#undef USE_FILIB	// <-- specify to use FILIB++ for interval arithmetic

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

#endif

////////////////////////////////////////////////////////////////////////
int main()
////////////////////////////////////////////////////////////////////////
{  

  MC::options.ENVEL_USE=true;
  MC::options.ENVEL_MAXIT=100;
  MC::options.ENVEL_TOL=1e-12;
  MC::options.MVCOMP_USE=true;

#ifdef SAVE_RESULTS
  ofstream allsub( "MC-2D.out", ios_base::out );
  ofstream dirsub( "MC2-2D.out", ios_base::out );
  allsub << scientific << setprecision(5) << right;
  dirsub << scientific << setprecision(5) << right;
#endif

  try{ 

    for( int iX=0; iX<NX; iX++ ){ 
      for( int iY=0; iY<NY; iY++ ){

        double Xval = XL+iX*(XU-XL)/(NX-1.);
        double Yval = YL+iY*(YU-YL)/(NY-1.);
        double Zval = myfunc( Xval, Yval );

        MC Xrel( I(XL,XU), Xval );
        MC Yrel( I(YL,YU), Yval );

        // Calculate relaxations + propagate all subgradient components
        Xrel.sub(2,0);
        Yrel.sub(2,1);
        MC Zrel = myfunc( Xrel, Yrel );

#ifdef SAVE_RESULTS
        allsub << setw(14) << Xval << setw(14) << Yval << setw(14) << Zval
               << setw(14) << Zrel.l() << setw(14) <<  Zrel.u()
               << setw(14) << Zrel.cv() << setw(14) << Zrel.cc()
               << setw(14) << Zrel.cvsub(0) << setw(14) << Zrel.ccsub(0)
               << setw(14) << Zrel.cvsub(1) << setw(14) << Zrel.ccsub(1)
               << endl;
#endif
        
        // Calculate relaxations + propagate directional subgradient
        const double dir[2] = {1,-1};
        Xrel.sub(1,&dir[0],&dir[0]);
        Yrel.sub(1,&dir[1],&dir[1]);
        Zrel = myfunc( Xrel, Yrel );

#ifdef SAVE_RESULTS
        dirsub << setw(14) << Xval << setw(14) << Yval << setw(14) << Zval
               << setw(14) << Zrel.l() << setw(14) <<  Zrel.u()
               << setw(14) << Zrel.cv() << setw(14) << Zrel.cc()
               << setw(14) << Zrel.cvsub(0) << setw(14) << Zrel.ccsub(0)
               << endl;
#endif
      }
#ifdef SAVE_RESULTS
      allsub << endl;
      dirsub << endl;
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
  allsub.close();
  dirsub.close();
#endif
  return 0;
}
