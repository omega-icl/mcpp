#define TEST_EXP	    // <-- select test function here
const int NTE = 4;	    // <-- select expansion order here
const int NX = 100;	    // <-- select X discretization here
#define SAVE_RESULTS    // <-- specify whether to save results to file
#define TEST_CVG        // <-- specify whether to save results to file
#define USE_PROFIL   	// <-- specify to use PROFIL for interval arithmetic
#undef  USE_FILIB	    // <-- specify to use FILIB++ for interval arithmetic
#define  USE_SPARSE      // <-- specify whether to use sparse Chebyshev models
#undef  MC__CVAR_SPARSE_PRODUCT_NAIVE
#undef  MC__POLYMODEL_DEBUG_SPROD
#undef  MC__CVAR_DEBUG_INTERPOLATION
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

#ifdef USE_SPARSE
  #include "scmodel.hpp"
  typedef mc::SCModel<I> CM;
  typedef mc::SCVar<I> CV;
  typedef mc::SCModel<MC> CMMC;
  typedef mc::SCVar<MC> CVMC;
#else
  #include "cmodel.hpp"
  typedef mc::CModel<I> CM;
  typedef mc::CVar<I> CV;
  typedef mc::CModel<MC> CMMC;
  typedef mc::CVar<MC> CVMC;
#endif

using namespace std;
using namespace mc;

////////////////////////////////////////////////////////////////////////

#if defined( TEST_POW )
const double XL   = -1.;	// <-- X range lower bound
const double XU   =  1.;	// <-- X range upper bound
const double Xref =  0.;	// <-- X ref point for McCormick
template <class T>
T myfunc
( const T&x )
{
  return 1. - 5.*x - pow(x,2)/2. + pow(x,3)/3.;
  //return 1. - 5.*x - x*x/2. + x*x*x/3. + x*x*x*x*x/9.;
}

#elif defined( TEST_DPOW )
const double XL   =  500;	// <-- X range lower bound
const double XU   =  3000;	// <-- X range upper bound
const double Xref =  1500;	// <-- X ref point for McCormick
template <class T>
T myfunc
( const T&x )
{
  return pow(x,0.6);
}

#elif defined( TEST_XLOG )
const double XL   =  -1;	// <-- X range lower bound
const double XU   =   1;	// <-- X range upper bound
const double Xref =   0;	// <-- X ref point for McCormick
template <class T>
T myfunc
( const T&x )
{
  return exp(sqr(x));//xlog(x);
}

#elif defined( TEST_POLY )
const double XL   = -0.125; //-.5;	// <-- X range lower bound
const double XU   = 0.5; // 1.;	// <-- X range upper bound
const double Xref =  0.;	// <-- X ref point for McCormick
template <class T>
T myfunc
( const T&x )
{
  return x*(x-1.1)*(x+2)*(x+2.2)*(x+2.5)*(x+3.)*sin(1.7*x+0.5);
  //return pow(x,6)+8.6*pow(x,5)+24.33*pow(x,4)+17.2*pow(x,3)-28.27*pow(x,2)-36.3*x;
}

#elif defined( TEST_EXP )
const double XL   = -2.;	// <-- X range lower bound
const double XU   =  2.;	// <-- X range upper bound
const double Xref =  0.;	// <-- X ref point for McCormick
template <class T>
T myfunc
( const T&x )
{
  return exp(x);
}

#elif defined( TEST_EXP1 )
const double XL   = 0.;	// <-- X range lower bound
const double XU   = 1.;	// <-- X range upper bound
const double Xref =  0.;	// <-- X ref point for McCormick
template <class T>
T myfunc
( const T&x )
{
  //return sqrt(x*x-x+.5)/sqrt(x*x+.5);
  return x*exp(-pow(x,2));
}

#elif defined( TEST_EXP2 )
const double XL   = -1.;	// <-- X range lower bound
const double XU   =  1.;	// <-- X range upper bound
const double Xref =  0.;	// <-- X ref point for McCormick
template <class T>
T myfunc
( const T&x )
{
  return exp(x)*pow(x+2,-3);
}

#elif defined( TEST_MIX )
const double XL   =  1.;	// <-- X range lower bound
const double XU   =  2.9;	// <-- X range upper bound
const double Xref = 1.5;	// <-- X ref point for McCormick
template <class T>
T myfunc
( const T&x )
{
  return log(cos(x));
  //return exp(x)*sqrt(x)*log(x)*cos(x);
  //return inv(x);
}

#elif defined( TEST_CHEB )
const double XL   = -1.;	// <-- range lower bound
const double XU   =  1.;	// <-- range upper bound
const double Xref = 0.5;	// <-- X ref point for McCormick
template <class T>
T myfunc
( const T&x )
{
  return cheb(x,2)+cheb(x,3);
}

#elif defined( TEST_INV )
const double XL   =  .2;	// <-- X range lower bound
const double XU   = 1.5;	// <-- X range upper bound
const double Xref =  1.;	// <-- X ref point for McCormick
template <class T>
T myfunc
( const T&x )
{
  return sqrt(1./x);
}

#elif defined( TEST_HYP )
const double XL   = -6.0;	// <-- X range lower bound
const double XU   = 8.0;	// <-- X range upper bound
const double Xref =  1.;	// <-- X ref point for McCormick
template <class T>
T myfunc
( const T&x )
{
  return tanh(x);//-sinh(x)/cosh(x);
}

#elif defined( TEST_TAN )
const double XL   = -1.2;	// <-- X range lower bound
const double XU   =  1.5;	// <-- X range upper bound
const double Xref =  0.;	// <-- X ref point for McCormick
template <class T>
T myfunc
( const T&x )
{
  return tan(x);//-sin(x)/cos(x);
}

#elif defined( TEST_ARC )
const double XL   = -.8;	// <-- X range lower bound
const double XU   =  .5;	// <-- X range upper bound
const double Xref =  1.;	// <-- X ref point for McCormick
template <class T>
T myfunc
( const T&x )
{
  return atan(x)-acos(x);
}

#elif defined( TEST_FABS )
const double XL   = -PI/2.;	// <-- range lower bound
const double XU   =  PI/4.;	// <-- range upper bound
const double Xref =  0.;	// <-- X ref point for McCormick
template <class T>
T myfunc
( const T&x )
{
  return fabs(x);
  //return sin(fabs(x));
}

#elif defined( TEST_TRIG )
const double XL   = -2*PI;	// <-- range lower bound
const double XU   = PI; 	// <-- range upper bound
const double Xref = PI; 	// <-- X ref point for McCormick
template <class T>
T myfunc
( const T&x )
{
  //return cos(x*atan(x));
  //return acos(x);
  return sin(x);
}

#elif defined( TEST_TRIG2 )
const double XL   = PI/6.;	// <-- X range lower bound
const double XU   = PI/3.;	// <-- X range upper bound
const double Xref = PI/4.;	// <-- X ref point for McCormick
template <class T>
T myfunc
( const T&x )
{
  return sin(pow(x,-3))*cos(sqrt(x));
}
#endif

////////////////////////////////////////////////////////////////////////
int main()
////////////////////////////////////////////////////////////////////////
{

#ifdef SAVE_RESULTS
  ofstream res( "CM-1D.out", ios_base::out );
  res << std::scientific << std::setprecision(5) << std::right;
#endif

  try{ 

    // Define Chebyshev model environment
#ifdef USE_SPARSE
    CM modCM( NTE );
    CMMC modCMMC( NTE );
#else
    CM modCM( 1, NTE );
    CMMC modCMMC( 1, NTE );
#endif

    // <-- set options here -->
    MC::options.MVCOMP_USE = true;
    //modCM.options.BOUNDER_TYPE = CM::Options::BERNSTEIN;
    //modCM.options.BOUNDER_ORDER = 40;
    modCM.options.DISPLAY_DIGITS = 10;
    //modCM.options.INTERP_EXTRA   = 1000;
    modCM.options.REMEZ_USE      = false;//true;
    modCM.options.MIXED_IA       = true;//false;
    modCM.options.BASIS          = 1; //0;//1;
    modCMMC.options = modCM.options;

    // Define variable X, and evaluate Chebyshev model
    //CVMC CVMCX( &mod, 0, MC(I(XL,XU),Xref).sub(1,0) );
    CV CVX( &modCM, 0, I(XL,XU) );
    std::cout << "\nChebyshev model of x:" << CVX;
    CV CVF = myfunc( CVX );
    std::cout << "\nChebyshev model of f(x):" << CVF;
    //std::cout << "\nBernstein bounder:" << CVMCF.bound(CMMC::Options::BERNSTEIN) << std::endl;

    // Define variable X, and evaluate Chebyshev-McCormick model
    //CVMC CVMCX( &mod, 0, MC(I(XL,XU),Xref).sub(1,0) );
    CVMC CVMCX( &modCMMC, 0, MC(I(XL,XU),Xref) );
    std::cout << "\nMcCormick-Chebyshev model of x:" << CVMCX;
    CVMC CVMCF = myfunc( CVMCX );
    std::cout << "\nMcCormick-Chebyshev model of f(x):" << CVMCF;

    // Repeated calculations at grid points (for display)
    for( int iX=0; iX<NX; iX++ ){ 

      double DX = XL+iX*(XU-XL)/(NX-1.);
      double DF = myfunc( DX );
      //CVMCX.set( &mod, 0, MC(I(XL,XU),DX).sub(1,0) );
      CVMCX.set( &modCMMC, 0, MC(I(XL,XU),DX) );
      CVMCF = myfunc( CVMCX );
      MC MCBF = CVMCF.B();
      double PF = CVMCF.P( &DX );
      MC MCF = PF + CVMCF.R();
      //std::cout << "CVMCF: " << CVMCF << std::endl;

      double a = CVMCF.linear(0,true);
      //std::cout << "CVMCF - lin: " << CVMCF << std::endl;
      MC MCRF = CVMCF.B();
      //std::cout << "MCRF: " << MCRF << std::endl;

#ifdef SAVE_RESULTS
      const unsigned IPREC = 9;
      res << std::scientific << std::setprecision(IPREC)
          << std::setw(IPREC+9) << DX << std::setw(IPREC+9) << std::setw(IPREC+9) << DF
          << std::setw(IPREC+9) << PF << std::setw(IPREC+9) << MCF.l()  << setw(IPREC+9) << MCF.u()
          << std::setw(IPREC+9) << MCF.cv() << setw(IPREC+9) << MCF.cc()
          << std::setw(IPREC+9) << MCBF.l()  << std::setw(IPREC+9) << MCBF.u()
          << std::setw(IPREC+9) << MCBF.cv() << std::setw(IPREC+9) << MCBF.cc()
          //<< std::setw(IPREC+9) << MCBF.cvsub(0) << std::setw(IPREC+9) << MCBF.ccsub(0)
          << std::setw(IPREC+9) << a*(DX-(XL+XU)/2.)+MCRF.l()  << std::setw(IPREC+9) << a*(DX-(XL+XU)/2.)+MCRF.u()
          << std::endl;
#endif
    }

#ifdef TEST_CVG
    // Repeated calculations for convergence test at range mid-point
    const unsigned N = 20;
    const double r = 2.;
    I Xred = I(XL,XU);
    CV CVXred( &modCM, 0, Xred );
    CV CVFred = myfunc( CVXred );
    std::cout << "\n  Range        Error         Order\n";
    for( unsigned i=0; i<N; i++ ){
      CV CVXred0 = CVXred;
      CV CVFred0 = CVFred;
      Xred = (Xred-Op<I>::mid(Xred))/r + Op<I>::mid(Xred);
      CVXred.set( &modCM, 0, Xred );
      CVFred = myfunc( CVXred );
      std::cout << Op<I>::diam(Xred) << "  " << Op<I>::diam(CVFred.R()) << "  "
                << (log(Op<I>::diam(CVFred0.R()))-log(Op<I>::diam(CVFred.R())))
                  /(log(Op<I>::diam(CVXred0.B()))-log(Op<I>::diam(CVXred.B()))) << std::endl;
    }
#endif

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
  catch( CMMC::Exceptions &eObj ){
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
