#define TEST_TRIG	// <-- select test function here
const int NTE = 2;	// <-- select Chebyshev expansion order here
const int NX = 200;	// <-- select X discretization here
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

#include "cmodel.hpp"
typedef mc::CModel<I> CM;
typedef mc::CVar<I> CV;
typedef mc::CModel<MC> CMMC;
typedef mc::CVar<MC> CVMC;

//#include "mcfadbad.hpp"
//typedef fadbad::F<CVMC> FCVMC;

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

#elif defined( TEST_ABS )
const double XL   = -1.;	// <-- X range lower bound
const double XU   =  1.;	// <-- X range upper bound
const double Xref =  0.;	// <-- X ref point for McCormick
template <class T>
T myfunc
( const T&x )
{
  return sqrt(sqr(x))*x;
}

#elif defined( TEST_POLY )
const double XL   = -5.55556e-02; //-0.125; //-.5;	// <-- X range lower bound
const double XU   = -4.68750e-02; // 0.5; // 1.;	// <-- X range upper bound
const double Xref =  0.;	// <-- X ref point for McCormick
template <class T>
T myfunc
( const T&x )
{
  return x*(x-1.1)*(x+2)*(x+2.2)*(x+2.5)*(x+3.)*sin(1.7*x+0.5);
  //return pow(x,6)+8.6*pow(x,5)+24.33*pow(x,4)+17.2*pow(x,3)-28.27*pow(x,2)-36.3*x;
}

#elif defined( TEST_EXP )
const double XL   = -1.5;	// <-- X range lower bound
const double XU   =  1.;	// <-- X range upper bound
const double Xref =  0.;	// <-- X ref point for McCormick
template <class T>
T myfunc
( const T&x )
{
  //return exp(-pow(x,-2));//fabs(x);
  //return sqrt(x*x-x+.5)/sqrt(x*x+.5);
  return x*exp(-pow(x,2));
  //return exp(-pow(x,2));
  //return exp(pow(x,3));
  //return -pow(x,2);
  //return x*exp(-x);
  //return x*pow(x,3);
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
  //return 1./pow(x+2,3);
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
  //return x;
  //return 1./x;
  return sqrt(1./x);
}

#elif defined( TEST_TRIG )
const double XL   =  0.;	// <-- range lower bound
const double XU   =  PI;	// <-- range upper bound
const double Xref = PI/3.;	// <-- X ref point for McCormick
template <class T>
T myfunc
( const T&x )
{
  //return cos(x*atan(x));
  //return acos(x);
  return x*sin(x);
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

  //CM modCM2( 4, 3 );//, true );
  //CV p[4];
  //for( unsigned i=0; i<4; i++ ) p[i].set( &modCM2, i, I(1,2) );
  //std::cout << p[0];
/*  CV x = -p[0]+p[1]-2*p[2];
  std::cout << x;
  x += p[0];
  std::cout << x;
  x -= p[1];
  std::cout << x;
  CV y = p[0]*p[1]*p[2]*p[3];
  std::cout << y;
  CV z = p[2]*p[3];
  std::cout << z;
  z *= p[0]*p[1];
  std::cout << z;
  z *= z;
  std::cout << z;
  z *= y;
  y += z;
  std::cout << y;*/
  //CV w = exp(p[1]/inv(p[0]));
  //std::cout << w;
  //return 0;
/*
  for( unsigned iord=0; iord<5; iord++ ){
    for( unsigned ivar=1; ivar<20; ivar++ )
      CM modCM2( ivar, iord );
    std::cout << std::endl;
  }

  CM modCM2( 4, 3 );
  CV p[4];
  for( unsigned i=0; i<4; i++ ) p[i].set( &modCM2, i, I(1,5) );
  CV f = (p[0]*p[3])*(p[0]+p[1]+p[2])+p[2];
  std::cout << f;
  return 0;
*/
#ifdef SAVE_RESULTS
  ofstream res( "CM-1D.out", ios_base::out );
  res << std::scientific << std::setprecision(5) << std::right;
#endif

  try{ 

    // Define Chebyshev model environment
    CM modCM( 1, NTE );
    CMMC modCMMC( 1, NTE );

    // <-- set options here -->
    MC::options.MVCOMP_USE = false;
    //modCM.options.BOUNDER_TYPE = CM::Options::BERNSTEIN;
    //modCM.options.BOUNDER_ORDER = 40;
    modCM.options.DISPLAY_DIGITS = 10;
    modCM.options.INTERP_EXTRA   = 0;//100;
    modCM.options.MIXED_IA       = true;
    modCMMC.options = modCM.options;

    // Define variable X, and evaluate Taylor model
    //CVMC CVMCX( &mod, 0, MC(I(XL,XU),Xref).sub(1,0) );
    CV CVX( &modCM, 0, I(XL,XU) );
    std::cout << "\nChebyshev model of x:" << CVX;
    CV CVF = myfunc( CVX );
    std::cout << "\nChebyshev model of f(x):" << CVF;
    //std::cout << "\nBernstein bounder:" << CVMCF.bound(CMMC::Options::BERNSTEIN) << std::endl;

    // Define variable X, and evaluate Taylor model
    //CVMC CVMCX( &mod, 0, MC(I(XL,XU),Xref).sub(1,0) );
    CVMC CVMCX( &modCMMC, 0, MC(I(XL,XU),Xref) );
    std::cout << "\nMcCormick-Chebyshev model of x:" << CVMCX;
    CVMC CVMCF = myfunc( CVMCX );
    std::cout << "\nMcCormick-Chebyshev model of f(x):" << CVMCF;
    //std::cout << "\nBernstein bounder:" << CVMCF.bound(CMMC::Options::BERNSTEIN) << std::endl;

//return 0;

    // Evaluate Taylor models of first derivatives using FADBAD++
    //FCVMC FCVMCX = CVMCX;
    //FCVMCX.diff(0,1);
    //FCVMC FCVMCF = myfunc( FCVMCX );
    //std::cout << "\nMcCormick-Chebyshev model of df/dx: " << FCVMCF.d(0) << std::endl;

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
      std::cout << "CVMCF: " << CVMCF << std::endl;

      double a = CVMCF.linear(0,true);
      std::cout << "CVMCF - lin: " << CVMCF << std::endl;
      MC MCRF = CVMCF.B();
      std::cout << "MCRF: " << MCRF << std::endl;

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
/*
    // Repeated calculations for convergence test at range mid-point
    const unsigned N = 20;
    const double r = 2.;
    I Xred = I(XL,XU);
    CV CVXred( &mod, 0, Xred );
    CV CVFred = myfunc( CVXred );
    for( unsigned i=0; i<N; i++ ){
      CV CVXred0 = CVXred;
      CV CVFred0 = CVFred;
      Xred = (Xred-Op<I>::mid(Xred))/r + Op<I>::mid(Xred);
      CVXred.set( &mod, 0, Xred );
      CVFred = myfunc( CVXred );
      std::cout << Op<I>::diam(Xred) << "  " << Op<I>::diam(CVFred.R()) << "  "
                << (log(Op<I>::diam(CVFred0.R()))-log(Op<I>::diam(CVFred.R())))
                  /(log(Op<I>::diam(CVXred0.B()))-log(Op<I>::diam(CVXred.B()))) << std::endl;
    }
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
