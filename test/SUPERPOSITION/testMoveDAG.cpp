#define TEST_TRIG	    // <-- select test function here
#define TEST_MOVE     // <-- enable the display of function calls (either copy or move)
#define MC__FISA_VERIFY // enable the verification of resultant ISM/SeISM enclosures
#define MC__FISA_DISPLAY_FISM  // enable displays (it needs to enable MC__FISA_VERIFY)

const unsigned long long int NSUB = 8;	// <-- select partition size here
const int NX = 256;	    // <-- select X discretization here
const int NY = 256;	    // <-- select Y discretization here
#undef USE_PROFIL	    // <-- specify to use PROFIL for interval arithmetic
#ifndef MC__USE_FILIB
  #define MC__USE_FILIB 
#endif  
//#undef MC__USE_FILIB
#undef MC__ISMODEL_TRACE 
////////////////////////////////////////////////////////////////////////

#include <fstream>    
#include <iomanip>  
 
#ifdef USE_PROFIL
  #include "mcprofil.hpp"
  typedef INTERVAL I;
#else
 #ifdef MC__USE_FILIB 
  #include "mcfilib.hpp"
  //typedef filib::interval<double,filib::native_switched,filib::i_mode_extended> I;
  typedef filib::interval<double, filib::rounding_strategy::native_switched, filib::interval_mode::i_mode_extended> I;
 #else
  #include "interval.hpp"
  typedef mc::Interval I; 
 #endif
#endif
  
#include "ismodel.hpp"
#include "ffunc.hpp"
typedef mc::ISModel<I> ISM;
typedef mc::ISVar<I> ISV;



#include <chrono>
 




using namespace mc;
////////////////////////////////////////////////////////////////////////

#if defined( TEST_MULT )
const double XL   =  1.;	// <-- X range lower bound
const double XU   =  2.;	// <-- X range upper bound
const double YL   =  0.;	// <-- Y range lower bound
const double YU   =  1.;	// <-- Y range upper bound
template <class T>
T myfunc
( const T&x, const T&y )
{
  return x*y;
}


#elif defined( TEST_EXP )
const double XL   =  1.;	// <-- X range lower bound
const double XU   =  2.;	// <-- X range upper bound
const double YL   =  0.;	// <-- Y range lower bound
const double YU   =  1.;	// <-- Y range upper bound
template <class T>
T myfunc
( const T&x, const T&y )
{
  using mc::sqr;
  return x*exp(x+sqr(y))-sqr(y);
}

#elif defined( TEST_TRIG )
template <class T>
T myTanh( const T&x)
{
  return 1-2./(exp(2.*x)+1);
} 
 
// const double XL   =  -2.1;	// <-- X range lower bound
// const double XU   =  -1.5;	  // <-- X range upper bound  
// const double YL   =  0.1;	// <-- Y range lower bound
// const double YU   =  2.7;	// <-- Y range upper bound
  
const double XL   =  -7.0;//7.8096536739679;	// <-- X range lower bound 
const double XU   =  1.0;//7.8234458656;	  // <-- X range upper bound
const double YL   =  -1.0;//7.0453216264833;	// <-- Y range lower bound
const double YU   =  3.0;//7.0918463404166;	// <-- Y range upper bound

// const double XL   =  -0.25;	// <-- X range lower bound
// const double XU   =  1.0;	  // <-- X range upper bound
// const double YL   =  0.5;	// <-- Y range lower bound
// const double YU   =  1.5;	// <-- Y range upper bound

// const double XL   =  0.;	// <-- X range lower bound
// const double XU   =  10.0;	  // <-- X range upper bound
// const double YL   =  0.;	// <-- Y range lower bound
// const double YU   =  20.;	// <-- Y range upper bound

// const double XL   =  3.79811943170321387;	// <-- X range lower bound
// const double XU   =  5.76074027115023846;	  // <-- X range upper bound
// const double YL   =  1.68107199008463071;	// <-- Y range lower bound
// const double YU   =  9.19320424169311146;	// <-- Y range upper bound


// const double XL   =  1.53758454701880032;	// <-- X range lower bound
// const double XU   =  1.56955970468480599;	  // <-- X range upper bound

// // const double YL   =  1.9622593320e1;
// // const double YU   =  1.9636673921e1;
// const double YL   =  10.186691405608883;	// <-- Y range lower bound
// const double YU   =  10.2336714768485955;	// <-- Y range upper bound

// const double XL   =  0.523573010650866;	// <-- X range lower bound
// const double XU   =  2.61801964293893;	  // <-- X range upper bound
// const double YL   =  5.;	// <-- Y range lower bound
// const double YU   =  10.;	// <-- Y range upper bound 
   
// const double XL   =  -PI;	// <-- X range lower bound
// const double XU   =  PI;	  // <-- X range upper bound
// const double YL   =  -PI;	// <-- Y range lower bound
// const double YU   =  PI;	// <-- Y range upper bound      

// template <class T>
// std::vector<T> myfunc
// ( const std::vector<T>& x)
// { 
//   std::vector<T> output;
//   output.resize(2); 
//   output[0] = -exp(sin(x[0])+sin(x[1])*cos(x[1]));
//   output[1] = exp(sin(x[0])+sin(x[1])*cos(x[1]));
//   return output; 
// } 

const double x1L   =  78.;	// <-- X range lower bound
const double x1U   =  102.;	  // <-- X range upper bound  
const double x3L   =  27.;	// <-- Y range lower bound
const double x3U   =  45.;	// <-- Y range upper bound 
const double x4L   =  27.;	// <-- X range lower bound 
const double x4U   =  45.;	  // <-- X range upper bound  
const double x5L   =  27.;	// <-- Y range lower bound
const double x5U   =  45.;	// <-- Y range upper bound 


template <class T>  
T myfunc 
( const T&x1, const T&x3, const T&x4, const T&x5) 
{
//  return sqr(0.5*inv(x5)+0.5*(inv(x3)*2275.1327 - 0.2668*x1 - 0.40584*x4));
  return inv(x5)*(inv(x3)*2275.1327 - 0.2668*x1 - 0.40584*x4);
  //return  (-sqr(0.5*(inv(x3)*2275.1327 - 0.2668*x1 - 0.40584*x4)-0.5*inv(x5))); 
 
}  
 
    
template <class T>  
T myfunc 
( const T&x, const T&y) 
{    
  //return ReLU(x+y+100);  
  //return 3*(1-x)*(1-x)*exp(-x*x-(y+1)*(y+1)) -  10*(x/5. - x*x*x - y*y*y*y*y)*exp(-x*x - y*y) - exp(-(x + 1)*(x + 1) - y*y)/3.;
  //return myTanh(x+y);
  //return max(x+y,0.);//(x+y)*(x+y)*(x+y);
  //return sqr(x)+sqr(y)-1.;//sqr(x)+sqr(y);//pow(x,4) + pow(y,4);//sqr(x)+sqr(y);//pow(x,2) + sqr(y);//x*y;//x*y+1;//pow(x+y,3);//sqr(x+y)*(x+y);//pow(x+y,3);//xlog(x+y);
  //return -exp(sin(x)+sqr(0.5*sin(y)+0.5*cos(y))-sqr(0.5*sin(y)-0.5*cos(y)));
  double cst = 0.1;
  return pow(exp(tanh(tanh(tanh(tanh(-x+(y*0.1)+x-x-cst))))),2.);//-exp(sin(x)+sin(y)*cos(y)); 
  //return sin(y)*cos(y);
  //return x*y;    
  ///return cos(x)+sin(y) + x + y;   
  //return 3.0*x+y-2.6*x;  
  //return -(sin(y)*cos(y));   
//  return sqr(x)+sqr(y)-100.;
  //return -(sin(x)+sin(y)*cos(y));  
  //return -(sin(x)+sin(2.0*y)/2.0);  
  //return exp(sin(x)+sin(2.0*y)/2.0);   
  //return exp(x+y);
  //return sqr(x)*(log(y)+exp(-y))-x*pow(log(y)+exp(-y),3);
  //return x*(x*(log(y)+exp(-y))-pow(log(y)+exp(-y),3));
  //return ((log(y)+exp(-y))*(x-sqr(log(y)+exp(-y))))*x;
  //return tanh(x+y);//sqr( sqr(x)*(log(y)+exp(-y))-x*(log(y)+exp(-y))) * ( sqr(x)*(log(y)+exp(-y))-x*(log(y)+exp(-y)) );  
  //return 1./(x+4);
  //return exp(cos(exp(3*x*x*x)+1)+y)-cos(exp(3*y*y*y)+1);
} 

 
#endif 

////////////////////////////////////////////////////////////////////////

int main()
{
 

 

  try{ 

    ISM mod( 2, NSUB );
    mod.options.ASYREM_USE = true; // non-constant remainder
    mod.options.DCDEC_USE = true;
    mod.options.SCALING_TYPE = ISM::Options::NONE;    
    mod.options.INTERSECTION_USE = true;   
    mod.options.SLOPE_USE = true;   
    I IX(XL,XU);
    I IY(YL,YU);
  
    ISV ISX( &mod, 0, IX );
    ISV ISY( &mod, 1, IY ); 

    auto start = std::chrono::steady_clock::now();
    auto ISF = myfunc( ISX, ISY);
    auto end = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed_seconds = end-start;
    std::cout << "elapsed time: " << elapsed_seconds.count() << "s\n";


    std::cout << "Interval superposition model of f:\n" << ISF << std::endl;


/*
    ISM mod4( 4, NSUB );
    mod4.options.ASYREM_USE = true; // non-constant remainder
    mod4.options.DCDEC_USE = true;
    mod4.options.SCALING_TYPE = ISM::Options::FULL;    
    mod4.options.INTERSECTION_USE = true;   
    mod4.options.SLOPE_USE = true;   
    I Ix1(x1L,x1U);  
    I Ix3(x3L,x3U);    
    I Ix4(x4L,x4U);
    I Ix5(x5L,x5U);    
    ISV ISx1( &mod4, 0, Ix1 );
    ISV ISx3( &mod4, 1, Ix3 );
    ISV ISx4( &mod4, 2, Ix4 );
    ISV ISx5( &mod4, 3, Ix5 );    
    auto ISF4 = myfunc(ISx1,ISx3,ISx4,ISx5);
   std::cout << "Interval superposition model of f:\n" << ISF4 << std::endl;    
    std::cout << " lb " << Op<ISV>::l(ISF4) << std::endl; 
*/

#ifdef MC__FISA_DISPLAY_FISM  
    std::ofstream ofile( "test2_ism.out", std::ofstream::out );
    ISF.display( 0, ofile );
    ofile.close();
#endif
    ///////////////////////////


    // auto IF  = myfunc( IX, IY);
    //auto ISF = myfunc( ISX, ISY);

    mc::FFGraph DAG; 
    mc::FFVar X( &DAG), Y(&DAG);
    mc::FFVar F = myfunc( X, Y );

    // Evaluate in interval arithmetic  
    // std::vector<I> IWK; 
    // DAG.eval( IWK, 1, &F, &IF, 1, &X, &IX, 1, &Y, &IY );
    // std::cout << "Interval inclusion of f:\n" << IF << std::endl;

    // Evaluate in interval superposition arithmetic
    std::vector<ISV> ISWK;

    start = std::chrono::steady_clock::now();
    DAG.eval( ISWK, 1, &F, &ISF, 1, &X, &ISX, 1, &Y, &ISY );
    end = std::chrono::steady_clock::now();
    elapsed_seconds = end-start;
    std::cout << "elapsed time: " << elapsed_seconds.count() << "s\n";
    std::cout << "Interval superposition model of f:\n" << ISF << std::endl;



    // I ttt(-1,1);
    // std::cout << "Interval" << ttt.l() << std::endl; 
   
    /////////////////////////// 
 
      
    // Repeated calculations at grid points (for display)   
        
#ifdef MC__FISA_VERIFY    
    std::ofstream ofile2( "test2_fct.out", std::ofstream::out );
    for( int iX=0; iX<NX; iX++ ){ 
      for( int iY=0; iY<NY; iY++ ){
        double DXY[2] = { ((double)(XL+iX*(XU-XL)/(NX-1.))), ((double)(YL+iY*(YU-YL)/(NY-1.))) };
        double DF = myfunc(DXY[0],DXY[1]);//( DXY[0]+ DXY[1] )*( DXY[0]+ DXY[1] ); 
        //double DF = std::max( DXY[0] + DXY[1],4.); 
        I BF = ISF.eval( DXY );
        if( ((double)mc::Op<I>::l(BF)) - DF > 1e-12 || DF - ((double)(mc::Op<I>::u(BF)))  >  1e-12)
          std::cout << std::setprecision(18) << mc::Op<I>::l(BF) << " < " << DF << " < "<< mc::Op<I>::u(BF) << "\n"
               << "at " << DXY[0] << "," << DXY[1] << std::endl;
#ifdef MC__FISA_DISPLAY_FISM                 
        ofile2 << std::setw(14) << DXY[0] << std::setw(14) << DXY[1] << std::setw(14) << DF
               << std::setw(14) << mc::Op<I>::l(BF) << std::setw(14) << mc::Op<I>::u(BF) 
               << std::endl;
#endif
      } 
      ofile2 << std::endl; 
    }
    ofile2.close();
#endif    



  } // end: try
  
#if !defined(MC__USE_PROFIL) && !defined(MC__USE_FILIB) 
  catch( I::Exceptions &eObj ){
    std::cerr << "Error " << eObj.ierr()
         << " in natural interval extension:" << endl
	 << eObj.what() << endl
         << "Aborts." << endl;
    return eObj.ierr();
  }
#endif
  catch( ISM::Exceptions &eObj ){
    std::cerr << "Error " << eObj.ierr()
              << " in interval superposition model:" << std::endl
              << eObj.what() << std::endl
              << "Aborts." << std::endl;
    return eObj.ierr();
  }

  return 0;
}

