#define TEST_MULT	// <-- select test function here
const int NTE = 3;	// <-- select Taylor expansion order here
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

#include "tmodel.hpp"
typedef mc::TModel<MC> TMMC;
typedef mc::TVar<MC> TVMC;

#include "mcfadbad.hpp"
typedef fadbad::F<TVMC> FTVMC;

using namespace std;
using namespace mc;

////////////////////////////////////////////////////////////////////////

#if defined( TEST_MULT )
const double XL   = -1.;	// <-- X range lower bound
const double XU   =  2.;	// <-- X range upper bound
const double Xref =  0.;	// <-- X ref point for McCormick
const double YL   = -2.;	// <-- Y range lower bound
const double YU   =  1.;	// <-- Y range upper bound
const double Yref =  0.;	// <-- Y ref point for McCormick
template <class T>
T myfunc
( const T&x, const T&y )
{
  return -pow(x+y,3)*x*y;
}

#elif defined( TEST_EXP )
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
  return x*exp(x+pow(y,2))-pow(y,2);
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
const double XL   = -2.;	// <-- X range lower bound
const double XU   =  0.;	// <-- X range upper bound
const double Xref = -1.;	// <-- X ref point for McCormick
const double YL   = -2.;	// <-- Y range lower bound
const double YU   =  0.;	// <-- Y range upper bound
const double Yref = -1.;	// <-- Y ref point for McCormick
template <class T>
T myfunc
( const T&x, const T&y )
{
  return +1./(pow(x-1.,3)+pow(y-1.,3)+0.1)
         -1./(pow(x-2.,2)+pow(y-3.,4)+0.2)
         +1./(pow(x-3.,3)+pow(y-2.,1)+0.2);
}

#elif defined( TEST_TRIG )
const double XL   = -0.5;	// <-- X range lower bound
const double XU   =  0.5;	// <-- X range upper bound
const double Xref =  0.;	// <-- X ref point for McCormick
const double YL   = -0.5;	// <-- Y range lower bound
const double YU   =  0.5;	// <-- Y range upper bound
const double Yref =  0.;	// <-- Y ref point for McCormick
template <class T>
T myfunc
( const T&x, const T&y )
{
  return 1.+x-sin(2.*x+3.*y)-cos(3.*x-5.*y);
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
////////////////////////////////////////////////////////////////////////
{  

#ifdef SAVE_RESULTS
  ofstream res( "TM-2D.out", ios_base::out );
  res << std::scientific << std::setprecision(5) << std::right;
#endif

  try{ 

    // Define Taylor model environment
    TMMC mod( 2, NTE );

    // <-- set options here -->
    MC::options.MVCOMP_USE = true;
    mod.options.BOUNDER_TYPE = TMMC::Options::BERNSTEIN;
    mod.options.BOUNDER_ORDER = 0;
    mod.options.BERNSTEIN_USE = true;

    // Define variables X and Y, and evaluate Taylor model
    TVMC TVMCX( &mod, 0, MC(I(XL,XU),Xref).sub(2,0) );
    TVMC TVMCY( &mod, 1, MC(I(YL,YU),Yref).sub(2,1) );
    TVMC TVMCF = myfunc( TVMCX, TVMCY );
    std::cout << "\nMcCormick-Taylor model of f(x,y):" << TVMCF;

    // Evaluate Taylor models of first derivatives using FADBAD++
    FTVMC FTVMCX = TVMCX;
    FTVMCX.diff(0,2);
    FTVMC FTVMCY = TVMCY;
    FTVMCY.diff(1,2);
    FTVMC FTVMCF = myfunc( FTVMCX, FTVMCY );
    std::cout << "\nMcCormick-Taylor model of df/dx: " << FTVMCF.d(0) << std::endl;
    std::cout << "\nMcCormick-Taylor model of df/dy: " << FTVMCF.d(1) << std::endl;

    // Repeated calculations at grid points (for display)
    for( int iX=0; iX<NX; iX++ ){ 
      for( int iY=0; iY<NY; iY++ ){ 

        double DXY[2] = { XL+iX*(XU-XL)/(NX-1.), YL+iY*(YU-YL)/(NY-1.) };
        double DF = myfunc( DXY[0], DXY[1] );
        TVMCX.set( &mod, 0, MC(I(XL,XU),DXY[0]).sub(2,0) );
        TVMCY.set( &mod, 1, MC(I(YL,YU),DXY[1]).sub(2,1) );
        TVMCF = myfunc( TVMCX, TVMCY );
        MC MCF = TVMCF.P( DXY ) + TVMCF.R();

#ifdef SAVE_RESULTS
      res << std::setw(14) << DXY[0] << std::setw(14) << DXY[1] << std::setw(14) << DF
          << std::setw(14) << MCF.l()  << std::setw(14) << MCF.u()
          << std::setw(14) << MCF.cv() << std::setw(14) << MCF.cc()
          << std::setw(14) << TVMCF.B().l()  << std::setw(14) << TVMCF.B().u()
          << std::setw(14) << TVMCF.B().cv() << std::setw(14) << TVMCF.B().cc()
          << std::setw(14) << TVMCF.B().cvsub(0) << std::setw(14) << TVMCF.B().ccsub(0)
          << std::setw(14) << TVMCF.B().cvsub(1) << std::setw(14) << TVMCF.B().ccsub(1)
          << std::endl;
#endif
     }
     res << std::endl;
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
  catch( TMMC::Exceptions &eObj ){
    cerr << "Error " << eObj.ierr()
         << " in Taylor model computation:" << endl
	 << eObj.what() << endl
         << "Aborts." << endl;
    return eObj.ierr();
  }

#ifdef SAVE_RESULTS
  res.close();
#endif
  return 0;
}
