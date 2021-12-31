#define TEST_TRIG      // <-- select test function here
#define USE_DAG        // <-- specify to evaluate via a DAG of the function
#undef  USE_BAD        // <-- specify to differentiate via backward AD
#define SAVE_RESULTS   // <-- specify whether to save results to file
const int NX = 40;	// <-- select X discretization here
const int NY = 40;	// <-- select Y discretization here
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

#include "specbnd.hpp"
typedef mc::Specbnd<I> SBI;

#ifdef USE_DAG
 #include "ffunc.hpp"
#else
#include "mcfadbad.hpp"
 typedef fadbad::F<double> FD;
 typedef fadbad::F<I> FI;
 #ifndef USE_BAD
  typedef fadbad::F<FD> FFD;
  typedef fadbad::F<FI> FFI;
 #else
  typedef fadbad::B<FD> BFD;
  typedef fadbad::B<FI> BFI;
 #endif
 typedef fadbad::F<SBI> FSBI;
#endif

using namespace std;
using namespace mc;

////////////////////////////////////////////////////////////////////////

#if defined( TEST_POLY )
const double XL   = -1.;	// <-- X range lower bound
const double XU   =  2.;	// <-- X range upper bound
const double XREF =  0.;	// <-- X ref point for McCormick
const double YL   = -2.;	// <-- Y range lower bound
const double YU   =  1.;	// <-- Y range upper bound
const double YREF =  0.;	// <-- Y ref point for McCormick
template <class T>
T myfunc
( const T&x, const T&y )
{
  return -pow(x+y,3)*x*y;
}

#elif defined( TEST_EXP )
const double XL   =  1.;	// <-- X range lower bound
const double XU   =  2.;	// <-- X range upper bound
const double XREF =  1.5;	// <-- X ref point for McCormick
const double YL   =  0.;	// <-- Y range lower bound
const double YU   =  1.;	// <-- Y range upper bound
const double YREF =  0.5;	// <-- Y ref point for McCormick
template <class T>
T myfunc
( const T&x, const T&y )
{
  return x*exp(x+pow(y,2))-pow(y,2);
}

#elif defined( TEST_EXP2 )
const double XL   = -2.;	// <-- X range lower bound
const double XU   =  2.;	// <-- X range upper bound
const double XREF =  0.;	// <-- X ref point for McCormick
const double YL   = -1.;	// <-- Y range lower bound
const double YU   =  1.;	// <-- Y range upper bound
const double YREF =  0.;	// <-- Y ref point for McCormick
template <class T>
T myfunc
( const T&x, const T&y )
{
  return x*y*(x*(exp(x)-exp(-x))-y*(exp(y)-exp(-y)));
}

#elif defined( TEST_INV )
const double XL   = -2.;	// <-- X range lower bound
const double XU   =  0.;	// <-- X range upper bound
const double XREF = -1.;	// <-- X ref point for McCormick
const double YL   =  1.;	// <-- Y range lower bound
const double YU   =  3.;	// <-- Y range upper bound
const double YREF =  2.;	// <-- Y ref point for McCormick
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
const double XREF = -1.;	// <-- X ref point for McCormick
const double YL   = -2.;	// <-- Y range lower bound
const double YU   =  0.;	// <-- Y range upper bound
const double YREF = -1.;	// <-- Y ref point for McCormick
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
const double XREF =  0.;	// <-- X ref point for McCormick
const double YL   = -1.;	// <-- Y range lower bound
const double YU   =  1.;	// <-- Y range upper bound
const double YREF =  0.;	// <-- Y ref point for McCormick
template <class T>
T myfunc
( const T&x, const T&y )
{
  return cheb(x+y,2)+cheb(x+y,3);
}

#elif defined( TEST_TRIG )
const double XL   = -0.5;	// <-- X range lower bound
const double XU   =  0.5;	// <-- X range upper bound
const double XREF =  0.;	// <-- X ref point for McCormick
const double YL   = -0.5;	// <-- Y range lower bound
const double YU   =  0.5;	// <-- Y range upper bound
const double YREF =  0.;	// <-- Y ref point for McCormick
template <class T>
T myfunc
( const T&x, const T&y )
{
  return 1.+x-sin(2.*x+3.*y)-cos(3.*x-5.*y);
}

#elif defined( TEST_TRIG2 )
const double XL   = -0.5;	// <-- X range lower bound
const double XU   =  0.5;	// <-- X range upper bound
const double XREF =  0.;	// <-- X ref point for McCormick
const double YL   = -0.5;	// <-- Y range lower bound
const double YU   =  0.5;	// <-- Y range upper bound
const double YREF =  0.;	// <-- Y ref point for McCormick
template <class T>
T myfunc
( const T&x, const T&y )
{
  return tan(x*y);
}

#elif defined( TEST_TRIG3 )
const double XL   = -0.5;	// <-- X range lower bound
const double XU   =  0.5;	// <-- X range upper bound
const double XREF =  0.;	// <-- X ref point for McCormick
const double YL   = -0.5;	// <-- Y range lower bound
const double YU   =  0.5;	// <-- Y range upper bound
const double YREF =  0.;	// <-- Y ref point for McCormick
template <class T>
T myfunc
( const T&x, const T&y )
{
  return asin(x)*acos(y);
}

#elif defined( TEST_NORM )
const double XL   =  0.5;	// <-- X range lower bound
const double XU   =  2.;	// <-- X range upper bound
const double XREF =  1.;	// <-- X ref point for McCormick
const double YL   =  0.5;	// <-- Y range lower bound
const double YU   =  2.;	// <-- Y range upper bound
const double YREF =  1.;	// <-- Y ref point for McCormick
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
  ofstream res;
  res.open( "SB-2D.out", ios_base::out );
  res << std::scientific << std::setprecision(5) << std::right;
#endif

  try{ 
#ifdef USE_DAG
    // Construct DAG representation of factorable function
    FFGraph DAG;
    FFVar X( &DAG );
    FFVar Y( &DAG );
    FFVar F = myfunc( X, Y );
    auto GF = DAG.subgraph( 1, &F );
#ifdef SAVE_RESULTS
    DAG.output( GF );
    ofstream ofdag( "SBF-2D.dot", ios_base::out );
    DAG.dot_script( 1, &F, ofdag );
    ofdag.close();
#endif
    // Construct DAG representation of Hessian matrix of factorable function
    FFVar *DF = DAG.FAD( 1, &F, 1, &X, 1, &Y );
#ifndef USE_BAD
    FFVar *D2F= DAG.FAD( 2, DF, 1, &X, 1, &Y );
#else
    FFVar *D2F= DAG.BAD( 2, DF, 1, &X, 1, &Y );
#endif
    auto GD2F = DAG.subgraph( 2*2, D2F );
#ifdef SAVE_RESULTS
    DAG.output( GD2F );
    ofstream od2fdag( "SBD2F-2D.dot", ios_base::out );
    DAG.dot_script( 2*2, D2F, od2fdag );
    ofdag.close();
#endif
#endif

    // Compute spectrum at reference point
#ifdef USE_DAG
    double DD2F[2*2];
    DAG.eval( 2*2, D2F, DD2F, 1, &X, &XREF, 1, &Y, &YREF );
    pair<double,double> specF = SBI::spectrum( 2, DD2F );
#else
    FD FXREF = XREF; FXREF.diff(0,2);
    FD FYREF = YREF; FYREF.diff(1,2);
#ifndef USE_BAD
    FFD FFXREF = FXREF; FFXREF.diff(0,2);
    FFD FFYREF = FYREF; FFYREF.diff(1,2);
    FFD FFFREF = myfunc( FFXREF, FFYREF );
    pair<double,double> specF = SBI::spectrum( FFFREF );
#else
    BFD BFXYREF[2] = { FXREF, FYREF };
    BFD BFFREF = myfunc( BFXYREF[0], BFXYREF[1] );
    BFFREF.diff(0,1);
    pair<double,double> specF = SBI::spectrum( BFXYREF );  
#endif
#endif
    cout << "\nSPECTRUM @REFERENCE POINT: " << I( specF.first, specF.second ) << endl;

    // Compute spectral interval inclusion using eigenvalue arithmetic
    I IX( XL, XU );
    I IY( YL, YU );
    SBI SBX( IX, 0, 2 );
    SBI SBY( IY, 1, 2 );
#ifdef USE_DAG
    SBI SBF;
    DAG.eval( 1, &F, &SBF, 1, &X, &SBX, 1, &Y, &SBY );
#else
    SBI SBF = myfunc( SBX, SBY );
#endif
    cout << "\nSPECTRAL BOUND (EIGENVALUE ARITHMETIC): " << SBF << endl;

    // Compute spectral bounds from interval Hessian matrix (forward-forward)
#ifdef USE_DAG
    I ID2F[2*2];
    DAG.eval( 2*2, D2F, ID2F, 1, &X, &IX, 1, &Y, &IY );

    SBI::options.HESSBND = SBI::Options::GERSHGORIN;
    pair<double,double> spbndG = SBI::spectral_bound( 2, ID2F );

    SBI::options.HESSBND = SBI::Options::HERTZROHN;
    pair<double,double> spbndHR = SBI::spectral_bound( 2, ID2F );
#else
    FI FX = IX; FX.diff(0,2);
    FI FY = IY; FY.diff(1,2);
#ifndef USE_BAD
    FFI FFX = FX; FFX.diff(0,2);
    FFI FFY = FY; FFY.diff(1,2);
    FFI FFF = myfunc( FFX, FFY );

    SBI::options.HESSBND = SBI::Options::GERSHGORIN;
    std::pair<double,double> spbndG = SBI::spectral_bound( FFF );

    SBI::options.HESSBND = SBI::Options::HERTZROHN;
    std::pair<double,double> spbndHR = SBI::spectral_bound( FFF );
#else
    BFI BFXY[2] = { FX, FY };
    BFI BFF = myfunc( BFXY[0], BFXY[1] );
    BFF.diff(0,1);

    SBI::options.HESSBND = SBI::Options::GERSHGORIN;
    std::pair<double,double> spbndG = SBI::spectral_bound( BFXY );

    SBI::options.HESSBND = SBI::Options::HERTZROHN;
    std::pair<double,double> spbndHR = SBI::spectral_bound( BFXY );
#endif
#endif
    cout << "\nSPECTRAL BOUND (GERSHGORIN): " << I( spbndG.first, spbndG.second) << endl
         << "\nSPECTRAL BOUND (HERTZ&ROHN): " << I( spbndHR.first, spbndHR.second) << endl;

    // Compute spectral bounds of first derivatives using FADBAD++
#ifdef USE_DAG
    SBI SBDF[2];
    DAG.eval( 2, DF, SBDF, 1, &X, &SBX, 1, &Y, &SBY );
    cout << "\nSPECTRAL BOUND OF DF/DX (EIGENVALUE ARITHMETIC): " << SBDF[0] << endl
         << "\nSPECTRAL BOUND OF DF/DY (EIGENVALUE ARITHMETIC): " << SBDF[1] << endl;
#else
    FSBI FSBX = SBX;
    FSBX.diff(0,2);
    FSBI FSBY = SBY;
    FSBY.diff(1,2);
    FSBI FSBF = myfunc( FSBX, FSBY );
    cout << "\nSPECTRAL BOUND OF DF/DX (EIGENVALUE ARITHMETIC): " << FSBF.d(0) << endl
         << "\nSPECTRAL BOUND OF DF/DY (EIGENVALUE ARITHMETIC): " << FSBF.d(1) << endl;
#endif

    // Repeated calculations at grid points
#ifdef SAVE_RESULTS
    for( int iX=0; iX<NX; iX++ ){
     for( int iY=0; iY<NY; iY++ ){
       double DX = XL+iX*(XU-XL)/(NX-1.);
       double DY = YL+iY*(YU-YL)/(NY-1.);
#ifdef USE_DAG
       DAG.eval( 2*2, D2F, DD2F, 1, &X, &DX, 1, &Y, &DY );
       specF = SBI::spectrum( 2, DD2F );
#else
       FD FX = DX; FX.diff(0,2);
       FD FY = DY; FY.diff(1,2);
#ifndef USE_BAD
       FFD FFX = FX; FFX.diff(0,2);
       FFD FFY = FY; FFY.diff(1,2);
       FFD FFF = myfunc( FFX, FFY );
       specF = SBI::spectrum( FFF );
#else
       BFD BFXY[2] = { FX, FY };
       BFD BFF = myfunc( BFXY[0], BFXY[1] );
       BFF.diff(0,1);
       specF = SBI::spectrum( BFXYREF );  
#endif
#endif
       res << std::setw(14) << DX << std::setw(14) << DY
           << std::setw(14) << specF.first << std::setw(14) << specF.second
           << std::setw(14) << Op<I>::l(SBF.SI()) << std::setw(14) << Op<I>::u(SBF.SI())
           << std::setw(14) << spbndG.first << std::setw(14) << spbndG.second
           << std::setw(14) << spbndHR.first << std::setw(14) << spbndHR.second
           << std::endl;
      }
      res << endl;
    }
#endif

#ifdef USE_DAG
    delete[] DF;
    delete[] D2F;
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
  catch( SBI::Exceptions &eObj ){
    cerr << "Error " << eObj.ierr()
         << " in spectral bound computation:" << endl
	 << eObj.what() << endl
         << "Aborts." << endl;
    return eObj.ierr();
  }

#ifdef SAVE_RESULTS
  res.close();
#endif
  return 0;
}
