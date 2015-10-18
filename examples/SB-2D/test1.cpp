#define TEST_CHEB	// <-- select test function here
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

#include "specbnd.hpp"
typedef mc::Specbnd<I> SB;


#include "mccormick.hpp"
typedef mc::McCormick<I> MC;
typedef mc::Specbnd<MC> SBMC;

#include "mcfadbad.hpp"
typedef fadbad::F<double> FD;
typedef fadbad::F<FD> FFD;
typedef fadbad::B<FD> BFD;

typedef fadbad::F<I> FI;
typedef fadbad::F<FI> FFI;
typedef fadbad::B<FI> BFI;

typedef fadbad::F<SB> FSB;

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

#elif defined( TEST_CHEB )
const double XL   = -1.;	// <-- range lower bound
const double XU   =  1.;	// <-- range upper bound
const double Xref =  0.;	// <-- X ref point for McCormick
const double YL   = -1.;	// <-- Y range lower bound
const double YU   =  1.;	// <-- Y range upper bound
const double Yref =  0.;	// <-- Y ref point for McCormick
template <class T>
T myfunc
( const T&x, const T&y )
{
  return cheb(x+y,2)+cheb(x+y,3);
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

#elif defined( TEST_TRIG2 )
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
  return tan(x*y);
}

#elif defined( TEST_TRIG3 )
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
  return asin(x)*acos(y);
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
  ofstream res;
#endif

   try{ 

#ifdef SAVE_RESULTS
    ofstream res;
    res.open( "test1.out", ios_base::out );
    res << std::scientific << std::setprecision(5) << std::right;
#endif

    // Compute spectral bounds using eigenvalue arithmetic
    SB SBX( I(XL,XU), 0, 2 );
    SB SBY( I(YL,YU), 1, 2 );
    SB SBF = myfunc( SBX, SBY );
    std::cout << "\nSpectral bound (eigenvalue arithmetic):\n" << SBF;

    // Compute spectrum at reference point (forward-forward)
    FD FXref = Xref; FXref.diff(0,2);
    FD FYref = Yref; FYref.diff(1,2);
    FFD FFXref = FXref; FFXref.diff(0,2);
    FFD FFYref = FYref; FFYref.diff(1,2);
    FFD FFFref = myfunc( FFXref, FFYref );
    std::pair<double,double> specF = SB::spectrum( FFFref );
    std::cout << "\nSpectrum at reference point: " << I(specF.first,specF.second)
              << std::endl;
/*
    // Compute spectral bounds from interval Hessian matrix (forward-forward)
    FI FX = I(XL,XU); FX.diff(0,2);
    FI FY = I(YL,YU); FY.diff(1,2);
    FFI FFX = FX; FFX.diff(0,2);
    FFI FFY = FY; FFY.diff(1,2);
    FFI FFF = myfunc( FFX, FFY );

    SB::options.HESSBND = SB::Options::GERSHGORIN;
    std::pair<double,double> spbndG = SB::spectral_bound( FFF );
    std::cout << "\nSpectral bound (Gershgorin, forward-forward): "
              << I(spbndG.first,spbndG.second) << std::endl;

    SB::options.HESSBND = SB::Options::HERTZROHN;
    std::pair<double,double> spbndHR = SB::spectral_bound( FFF );
    std::cout << "\nSpectral bound (Hertz&Rohn, forward-forward): "
              << I(spbndHR.first,spbndHR.second) << std::endl;

    // Compute spectral bounds of first derivatives using FADBAD++
    FSB FSBX = SBX;
    FSBX.diff(0,2);
    FSB FSBY = SBY;
    FSBY.diff(1,2);
    FSB FSBF = myfunc( FSBX, FSBY );
    std::cout << "\nSpectral bounds of df/dx:\n" << FSBF.d(0) << std::endl;
    std::cout << "\nSpectral bounds of df/dy:\n" << FSBF.d(1) << std::endl;

#ifdef SAVE_RESULTS
    // Repeated calculations at grid points (for display)
    for( int iX=0; iX<NX; iX++ ){ 
      for( int iY=0; iY<NY; iY++ ){ 

        double XY[2] = { XL+iX*(XU-XL)/(NX-1.), YL+iY*(YU-YL)/(NY-1.) };
        FD FXY[2] = { XY[0], XY[1] }; FXY[0].diff(0,2); FXY[1].diff(1,2);
        BFD BFXY[2] = { FXY[0], FXY[1] };
        BFD BFF = myfunc( BFXY[0], BFXY[1] );
        BFF.diff(0,1);
        std::pair<double,double> spec = SB::spectrum( BFXY );
        SBMC SBMCX( MC(I(XL,XU),XY[0]).sub(2,0), 0, 2 );
        SBMC SBMCY( MC(I(YL,YU),XY[1]).sub(2,1), 1, 2 );
        SBMC SBMCF = myfunc( SBMCX, SBMCY );

        res << std::setw(14) << XY[0] << std::setw(14) << XY[1]
            << std::setw(14) << spec.first << std::setw(14) << spec.second
            << std::setw(14) << SBMCF.SI().l()  << std::setw(14) << SBMCF.SI().u()
            << std::setw(14) << SBMCF.SI().cv() << std::setw(14) << SBMCF.SI().cc()
            << std::setw(14) << BFF.x().x()
            << std::setw(14) << SBMCF.I().l()  << std::setw(14) << SBMCF.I().u()
            << std::setw(14) << SBMCF.I().cv() << std::setw(14) << SBMCF.I().cc()
	    << std::setw(14) << BFF.x().d(0)
            << std::setw(14) << SBMCF.FI().deriv(0).l()  << std::setw(14) << SBMCF.FI().deriv(0).u()
            << std::setw(14) << SBMCF.FI().deriv(0).cv() << std::setw(14) << SBMCF.FI().deriv(0).cc()
	    << std::setw(14) << BFF.x().d(1)
            << std::setw(14) << SBMCF.FI().deriv(1).l()  << std::setw(14) << SBMCF.FI().deriv(1).u()
            << std::setw(14) << SBMCF.FI().deriv(1).cv() << std::setw(14) << SBMCF.FI().deriv(1).cc()
            << std::setw(14) << spbndG.first << std::setw(14) << spbndG.second
            << std::setw(14) << spbndHR.first << std::setw(14) << spbndHR.second
            << std::endl;
      }
      res << endl;
    }
#endif
*/
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
/*
  catch( MC::Exceptions &eObj ){
    cerr << "Error " << eObj.ierr()
         << " in McCormick relaxation:" << endl
	 << eObj.what() << endl
         << "Aborts." << endl;
    return eObj.ierr();
  }
  catch( SB::Exceptions &eObj ){
    cerr << "Error " << eObj.ierr()
         << " in spectral bound computation:" << endl
	 << eObj.what() << endl
         << "Aborts." << endl;
    return eObj.ierr();
  }
  catch( SBMC::Exceptions &eObj ){
    cerr << "Error " << eObj.ierr()
         << " in spectral bound computation:" << endl
	 << eObj.what() << endl
         << "Aborts." << endl;
    return eObj.ierr();
  }
*/
#ifdef SAVE_RESULTS
  res.close();
#endif
  return 0;
}
