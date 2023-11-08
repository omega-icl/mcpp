#define TEST_RELU	    // <-- select test function here
#define SAVE_RESULTS    // <-- specify whether to save results to file
#undef  USE_OPTIM       // <-- specify whether to use optimization
const int NSUB = 8;	    // <-- select partition size here
const int NX = 200;	    // <-- select X discretization here
const int NY = 200;	    // <-- select Y discretization here
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


#include "univarmodels.hpp" 
#include "asmodel.hpp"
typedef mc::ASModel<I> ASM;
typedef mc::ASVar<I>   ASV;

#ifdef USE_OPTIM
 #include "gurobi_c++.h"
#endif

using namespace std;
using namespace mc;

////////////////////////////////////////////////////////////////////////

#if defined( TEST_MULT )
const double XL   =  -20.;	// <-- X range lower bound
const double XU   =  20.;	// <-- X range upper bound
const double YL   =  -1.;	// <-- Y range lower bound
const double YU   =  1.;	// <-- Y range upper bound
template <class T>
T myfunc
( const T&x, const T&y )
{
  return x*y;
}

#elif defined( TEST_SQRT )
const double XL   =  1.;	// <-- X range lower bound
const double XU   =  2.;	// <-- X range upper bound
const double YL   =  0.5;	// <-- Y range lower bound
const double YU   =  1.;	// <-- Y range upper bound
template <class T>
T myfunc
( const T&x, const T&y )
{
  return sqrt(x+y);
}

#elif defined( TEST_SQR )
const double XL   =  -2.;	// <-- X range lower bound
const double XU   =  10.;	// <-- X range upper bound
const double YL   =  -2.;	// <-- Y range lower bound
const double YU   =  10.;	// <-- Y range upper bound
template <class T>
T myfunc
( const T&x, const T&y )
{
  using mc::sqr;
  return sqr(x+y);
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
  return x*exp(x+sqr(y));
}

#elif defined( TEST_XLOG )
const double XL   =  0.2;	// <-- X range lower bound
const double XU   =  2.;	// <-- X range upper bound
const double YL   =  0.;	// <-- Y range lower bound
const double YU   =  1.;	// <-- Y range upper bound
template <class T>
T myfunc
( const T&x, const T&y )
{
  return xlog(x+y);
}

#elif defined( TEST_EXP2 )
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
const double XL   =  0.;	// <-- X range lower bound
const double XU   =  10.;	// <-- X range upper bound
const double YL   =  0.;	// <-- Y range lower bound
const double YU   =  20.;	// <-- Y range upper bound
template <class T>
T myfunc
( const T&x, const T&y )
{
  return exp( sin(x) + sin(y)*cos(y) );
}

#elif defined( TEST_POW )
const double XL   =  -2.;	// <-- X range lower bound
const double XU   =   2.;	// <-- X range upper bound
const double YL   =  -2.;	// <-- Y range lower bound
const double YU   =   2.;	// <-- Y range upper bound
template <class T>
T myfunc
( const T&x, const T&y )
{
  return pow(x+y,3);
}

#elif defined( TEST_POW2 )
const double XL   =  -2.;	// <-- X range lower bound
const double XU   =   2.;	// <-- X range upper bound
const double YL   =  -2.;	// <-- Y range lower bound
const double YU   =   2.;	// <-- Y range upper bound
template <class T>
T myfunc
( const T&x, const T&y )
{
  return pow(x,3)+pow(y,3);
}

#elif defined( TEST_INV )
const double XL   = -0.6;	// <-- X range lower bound
const double XU   = -0.5;	// <-- X range upper bound
const double YL   = -0.5;	// <-- Y range lower bound
const double YU   = -0.2;	// <-- Y range upper bound
template <class T>
T myfunc
( const T&x, const T&y )
{
  return 1./(x+y);
}

#elif defined( TEST_MAX )
const double XL   =  -2.;	// <-- X range lower bound
const double XU   =   2.;	// <-- X range upper bound
const double YL   =  -2.;	// <-- Y range lower bound
const double YU   =   2.;	// <-- Y range upper bound
template <class T>
T myfunc
( const T&x, const T&y )
{
  return max(x+y,0.);
}

#elif defined( TEST_TANH )
const double XL   =  -1.;	// <-- X range lower bound
const double XU   =   2.;	// <-- X range upper bound
const double YL   =  -2.;	// <-- Y range lower bound
const double YU   =   1.;	// <-- Y range upper bound
template <class T>
T myfunc
( const T&x, const T&y )
{
  return tanh(x+y);
}

#elif defined( TEST_TANH2 )  
const double XL   =  -1.;	// <-- X range lower bound
const double XU   =   2.;	// <-- X range upper bound
const double YL   =  -2.;	// <-- Y range lower bound
const double YU   =   1.;	// <-- Y range upper bound
template <class T> 
T myfunc
( const T&x, const T&y )
{
  return tanh(x)+tanh(y);
}

#elif defined( TEST_RELU )
const double XL   =  -1.;	// <-- X range lower bound  
const double XU   =  2.;	// <-- X range upper bound
const double YL   =  -2.;	// <-- Y range lower bound
const double YU   =  1.;	// <-- Y range upper bound
template <class T>
T myfunc
( const T&x, const T&y )
{
  return  max(max(x+y+1,0.) + -max(x+y,0.)  -0.3 ,0.); // (-max(x+y+1,0.)) + (-max(x+y,0.));// max((-max(x+y+1,0.)) + (-max(x+y,0.))+7,0.0);//;max(x+y+1,0.); //,0.);
}  
#endif

////////////////////////////////////////////////////////////////////////   
   
int main() 
{
  try{          

    ASM mod( 2, NSUB );        
    mod.options.ASYREM_USE = true;    
    mod.options.DCDEC_USE  = true;   
    mod.options.SHADOW_USE  = true;     
    mc::UnivarPWLE<double>::nbpsMax = 1000;       
    I IX( XL, XU );        
    I IY( YL, YU );       
    ASV ISX( &mod, 0, IX );       
    ASV ISY( &mod, 1, IY );
    ISX.debug_check_overNunder_flags();          
    ISY.debug_check_overNunder_flags();       
    auto ISZ1 = -relu(ISX + ISY);     
    auto ISZ2 = relu(ISX + ISY + 1);        
    auto ISZ3 =  relu(ISZ2 + ISZ1 - 0.3); // ISZ2;// - ISZ1;   

/*
    auto outC = ISZ3.C();
    for(unsigned jj = 0; jj < 8; jj++){
      std::cout << outC[0][jj];      
    }
    std::cout << std::endl;
    for(unsigned jj = 0; jj < 8; jj++){
      std::cout << outC[1][jj];      
    }
    std::cout << std::endl;
    std::cout << ISZ3 << std::endl;
*/ 
    //std::cout << ISZ2 + ISZ1 << std::endl;       
/*       
    auto IF  = myfunc( IX, IY );     
    std::cout << "Interval inclusion of f:\n" << IF << std::endl;  
     
    auto ISF = myfunc( ISX, ISY );    
    std::cout << "Interval superposition model of f:\n" << ASF << std::endl;
*/  
#ifdef USE_OPTIM
    I IDX[NSUB], IDY[NSUB], IDF[NSUB][NSUB];  
    for( int iX=0; iX<NSUB; iX++ )  
      IDX[iX] = XL+I(iX,iX+1)*(XU-XL)/NSUB; 
    for( int iY=0; iY<NSUB; iY++ ) 
      IDY[iY] = YL+I(iY,iY+1)*(YU-YL)/NSUB; 
    for( int iX=0; iX<NSUB; iX++ )  
      for( int iY=0; iY<NSUB; iY++ )
        IDF[iX][iY] = myfunc( IDX[iX], IDY[iY] );

    GRBEnv GRBenv;
    GRBModel GRBmodel( GRBenv );
    GRBmodel.getEnv().set( GRB_DoubleParam_FeasibilityTol, 1e-9 );
    GRBmodel.getEnv().set( GRB_DoubleParam_OptimalityTol,  1e-9 ); 

    GRBVar FXL[NSUB], FXU[NSUB], FYL[NSUB], FYU[NSUB];  
    for( int iX=0; iX<NSUB; iX++ ){
      FXL[iX] = GRBmodel.addVar( -GRB_INFINITY, GRB_INFINITY, 0.0, GRB_CONTINUOUS );
      FXU[iX] = GRBmodel.addVar( -GRB_INFINITY, GRB_INFINITY, 0.0, GRB_CONTINUOUS );
    }
    for( int iY=0; iY<NSUB; iY++ ){
      FYL[iY] = GRBmodel.addVar( -GRB_INFINITY, GRB_INFINITY, 0.0, GRB_CONTINUOUS );
      FYU[iY] = GRBmodel.addVar( -GRB_INFINITY, GRB_INFINITY, 0.0, GRB_CONTINUOUS );
    } 
    //GRBVar ERR = GRBmodel.addVar( 0, GRB_INFINITY, 1.0, GRB_CONTINUOUS );    
    //GRBmodel.update(); 
    GRBVar ERR[NSUB][NSUB];    
    for( int iX=0; iX<NSUB; iX++ )
      for( int iY=0; iY<NSUB; iY++ )
        ERR[iX][iY] = GRBmodel.addVar( 0, GRB_INFINITY, 1./(NSUB*NSUB), GRB_CONTINUOUS );
    GRBmodel.update();

    for( int iX=0; iX<NSUB; iX++ ) 
      GRBmodel.addConstr( FXL[iX], GRB_LESS_EQUAL, FXU[iX] );
    for( int iY=0; iY<NSUB; iY++ )
      GRBmodel.addConstr( FYL[iY], GRB_LESS_EQUAL, FYU[iY] );
    for( int iX=0; iX<NSUB; iX++ ){ 
      for( int iY=0; iY<NSUB; iY++ ){
        GRBmodel.addConstr( FXL[iX] + FYL[iY], GRB_LESS_EQUAL,    Op<I>::l( IDF[iX][iY] ) );
        GRBmodel.addConstr( FXL[iX] + FYL[iY], GRB_GREATER_EQUAL, Op<I>::l( IDF[iX][iY] ) - ERR[iX][iY] ); //ERR );
        GRBmodel.addConstr( FXU[iX] + FYU[iY], GRB_GREATER_EQUAL, Op<I>::u( IDF[iX][iY] ) );
        GRBmodel.addConstr( FXU[iX] + FYU[iY], GRB_LESS_EQUAL,    Op<I>::u( IDF[iX][iY] ) + ERR[iX][iY] ); //ERR );
        GRBmodel.addConstr( FXL[iX] + FYL[iY], GRB_GREATER_EQUAL, Op<I>::l( IF ) );
        GRBmodel.addConstr( FXU[iX] + FYU[iY], GRB_LESS_EQUAL,    Op<I>::u( IF ) );
      }
    }

    GRBmodel.set( GRB_IntAttr_ModelSense, 1 ); // MIN:1, MAX:-1
    GRBmodel.update();
    GRBmodel.write( "ISMOpt.lp" ); 
    GRBmodel.optimize();
    //double ERRMIN = GRBmodel.get( GRB_DoubleAttr_ObjVal ); 

    ISF.set( &mod );
    ISF.C()[0].resize( NSUB );
    for( int iX=0; iX<NSUB; iX++ )
      ISF.C()[0][iX] = I( FXL[iX].get(GRB_DoubleAttr_X), FXU[iX].get(GRB_DoubleAttr_X) );
    ISF.C()[1].resize( NSUB );
    for( int iY=0; iY<NSUB; iY++ )
      ISF.C()[1][iY] = I( FYL[iY].get(GRB_DoubleAttr_X), FYU[iY].get(GRB_DoubleAttr_X) );
    ISF.ndep() += 2;
    std::cout << "Optimal interval superposition model of f:\n" << ISF << std::endl;
#endif 

   
#ifdef SAVE_RESULTS
    std::ofstream ofile( "test2_ism.out", std::ofstream::out );  
    ISZ3.display( true,0, ofile );    
    ofile.close(); 
  
    // Repeated calculations at grid points (for display) 
    std::ofstream ofile2( "test2_fct.out", std::ofstream::out ); 
    for( int iX=0; iX<NX; iX++ ){
      for( int iY=0; iY<NY; iY++ ){
        double DXY[2] = { XL+iX*(XU-XL)/(NX-1.), YL+iY*(YU-YL)/(NY-1.) };
        double DF = myfunc( DXY[0], DXY[1] );
        I BF = ISZ3.eval( DXY );
        if ((DF < mc::Op<I>::l(BF)-MC__ASM_COMPUTATION_TOL ) || DF > mc::Op<I>::u(BF) + MC__ASM_COMPUTATION_TOL) std::cout<< mc::Op<I>::l(BF) << " "<< DF << " " << mc::Op<I>::u(BF)<< std::endl;
        ofile2 << std::setw(14) << DXY[0] << std::setw(14) << DXY[1] << std::setw(14) << DF
               << std::setw(14) << mc::Op<I>::l(BF) << std::setw(14) << mc::Op<I>::u(BF) 
               << std::endl;
      }
      ofile2 << std::endl;
    }    
    ofile2.close();  
#endif
  } // end: try   
    
#ifndef USE_PROFIL 
#ifndef MC__USE_FILIB 
  catch( I::Exceptions &eObj ){
    std::cerr << "Error " << eObj.ierr()
              << " in natural interval extension:" << std::endl
              << eObj.what() << std::endl
              << "Aborts." << std::endl;
    return eObj.ierr();
  }
#endif
#endif
  catch( ASM::Exceptions &eObj ){
    std::cerr << "Error " << eObj.ierr()
              << " in interval superposition model:" << std::endl
              << eObj.what() << std::endl
              << "Aborts." << std::endl;
    return eObj.ierr();
  }
#ifdef USE_OPTIM
  catch(GRBException& ex) {
    std::cerr << "Error code = " << ex.getErrorCode() << std::endl;
    std::cerr << ex.getMessage() << std::endl;
  }
#endif

  return 0;
}

