#define TEST_ANN	    // <-- select test function here
const int NSUB = 48;	// <-- select partition size here
const int NX = 200;	    // <-- select X discretization here
const int NY = 200;	    // <-- select Y discretization here
#undef USE_PROFIL	    // <-- specify to use PROFIL for interval arithmetic
#undef MC__ISMODEL_TRACE
////////////////////////////////////////////////////////////////////////

#include <fstream>

#ifdef USE_PROFIL
  #include "mcprofil.hpp"
  typedef INTERVAL I;
#else
  #include "interval.hpp"
  typedef mc::Interval I;
#endif

#include "ismodel.hpp"
#include "ffunc.hpp"
typedef mc::ISModel<I> ISM;
typedef mc::ISVar<I> ISV;

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

#elif defined( TEST_SQRT )
const double XL   =  1.;	// <-- X range lower bound
const double XU   =  2.;	// <-- X range upper bound
const double YL   =  0.;	// <-- Y range lower bound
const double YU   =  1.;	// <-- Y range upper bound
template <class T>
T myfunc
( const T&x, const T&y )
{
  return sqrt(x+y);
}

#elif defined( TEST_SQR )
const double XL   =  -10.;	// <-- X range lower bound
const double XU   =  10.;	// <-- X range upper bound
const double YL   =  -10.;	// <-- Y range lower bound
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
  return x*exp(x+sqr(y))-sqr(y);
}

#elif defined( TEST_TRIG )
const double XL   =  -140.;	// <-- X range lower bound
const double XU   =  140.;	// <-- X range upper bound
const double YL   =  -80.;	// <-- Y range lower bound
const double YU   =  80.;	// <-- Y range upper bound
template <class T>
T myfunc
( const T&x, const T&y , const T&z)
{
  return ReLU(x+y+100);//exp(cos(exp(3*x*x*x)+1)+y)-cos(exp(3*y*y*y)+1);
  
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
#endif

////////////////////////////////////////////////////////////////////////

std::vector<std::vector<double>> getWeights(std::string filename){

  std::ifstream inputParas (filename);
  //inputParas.open("w1.csv",std::ios::out)
  if (!inputParas.is_open()) throw std::runtime_error("Could not open file");
	else {
    std::vector<std::vector<double>> w1;
    char comma;
    std::string _linestr;
      
    while(std::getline(inputParas,_linestr)){
      std::vector<double> w1i;
      double _weight;
      std::istringstream _line(_linestr);
    
      while (_line >> _weight){      // read a number  
        w1i.push_back(_weight);  // store the number 
        _line >> comma;         // read and discard the comma
      }
      w1.push_back(w1i);
    }
		inputParas.close(); 
	  // std::cout << w1[1][1] <<std::endl;
    return w1;
  }
}

double ReLU(double x){
  return std::max(x,0.);
}

#if defined( TEST_ANN )
template <class T>
T myTanh( const T&x)
{
  return 1-2./(exp(2.*x)+1);
}

const double XL   =  0.2;	// <-- X range lower bound
const double XU   =  0.4;	// <-- X range upper bound
const double YL   =  -1.8;	// <-- Y range lower bound
const double YU   =  -1.6;	// <-- Y range upper bound
template <class T>
T myfunc
( const T&x, const T&y)
{
  std::vector<std::vector<double>> w1 (getWeights("w1Tanh.csv"));
  std::vector<std::vector<double>> w2 (getWeights("w2Tanh.csv"));
  std::vector<std::vector<double>> b1 (getWeights("b1Tanh.csv"));
  std::vector<std::vector<double>> b2 (getWeights("b2Tanh.csv"));  
  std::vector<T> V1;
  V1.push_back(myTanh(w1[0][0]*x+w1[0][1]*y+b1[0][0])*w2[0][0]+b2[0][0]);
  for(int i=1;i<47;i++){
    V1[0] += myTanh(w1[i][0]*x+w1[i][1]*y+b1[i][0])*w2[0][i];}

  //for(int i=0;i<47;i++){
  //  V1.push_back(ReLU(w1[i][0]*x+w1[i][1]*y));
  //}
  //V1.push_back(ReLU(w1[46][0]*x+w1[46][1]*y)*w2[0][46]);
  //for(int i=0;i<46;i++)
  //  V1[47] += V1[i]*w2[0][i];

  //return V1[47];//exp(cos(exp(3*x*x*x)+1)+y)-cos(exp(3*y*y*y)+1);
  return V1[0];
}
#endif




int main()
{
  try{ 

    ISM mod( 3, NSUB );
    I IX(XL,XU);
    I IY(YL,YU);
    ISV ISX( &mod, 0, IX );
    ISV ISY( &mod, 1, IY );


    auto IF  = myfunc( IX, IY); 
    std::cout << "Interval inclusion of f:\n" << IF << std::endl;
  
    auto ISF = myfunc( ISX, ISY);
    std::cout << "Interval superposition model of f:\n" << ISF << std::endl;
    //auto&& mat = ISF.C();


    // For debuging reading weights from a csv file. 
    //std::vector<std::vector<double>> w1 (getWeights("w2.csv"));
    //std::cout << w1[0][1] <<std::endl;


    std::ofstream ofile( "test2_ism.out", std::ofstream::out );
    ISF.display( 0, ofile );
    ofile.close();



    ///////////////////////////
    /*
    mc::FFGraph DAG;
    mc::FFVar X( &DAG), Y(&DAG);
    mc::FFVar F = myfunc( X, Y );

    // Evaluate in interval arithmetic
    std::vector<I> IWK;
    DAG.eval( IWK, 1, &F, &IF, 1, &X, &IX, 1, &Y, &IY );
    std::cout << "Interval inclusion of f:\n" << IF << std::endl;

    // Evaluate in interval superposition arithmetic
    std::vector<ISV> ISWK;
    DAG.eval( ISWK, 1, &F, &ISF, 1, &X, &ISX, 1, &Y, &ISY );
    std::cout << "Interval superposition model of f:\n" << ISF << std::endl;
     */
    // I ttt(-1,1);
    // std::cout << "Interval" << ttt.l() << std::endl;

    ///////////////////////////

    
    // Repeated calculations at grid points (for display)
    std::ofstream ofile2( "test2_fct.out", std::ofstream::out );
    for( int iX=0; iX<NX; iX++ ){
      for( int iY=0; iY<NY; iY++ ){
        double DXY[2] = { XL+iX*(XU-XL)/(NX-1.), YL+iY*(YU-YL)/(NY-1.) };
        double DF = myfunc( DXY[0], DXY[1]);
        I BF = ISF.eval( DXY );
        ofile2 << std::setw(14) << DXY[0] << std::setw(14) << DXY[1] << std::setw(14) << DF
               << std::setw(14) << mc::Op<I>::l(BF) << std::setw(14) << mc::Op<I>::u(BF) 
               << std::endl;
      }
      ofile2 << std::endl;
    }
    ofile2.close();

    ///////////////////////////
    /*
    // Repeated calculations for different orders
    std::ofstream ofile3( "test2_div.out", std::ofstream::out );
    for( unsigned int q=0; q<16; q++ ){
      ISM mod( 2, pow(2,q) );
      I IX(XL,XU);
      I IY(YL,YU);
      ISV ISX( &mod, 0, IX );
      ISV ISY( &mod, 1, IY );

      auto ISF = myfunc( ISX, ISY );
      ofile3 << std::setw(5) << pow(2,q)
             << std::setw(14) << mc::Op<I>::l(ISF.B())
             << std::setw(14) << mc::Op<I>::u(ISF.B())
             << std::endl;
    }
    ofile3.close();
  */  
  } // end: try
  
#ifndef USE_PROFIL
  catch( I::Exceptions &eObj ){
    std::cerr << "Error " << eObj.ierr()
              << " in natural interval extension:" << std::endl
              << eObj.what() << std::endl
              << "Aborts." << std::endl;
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

