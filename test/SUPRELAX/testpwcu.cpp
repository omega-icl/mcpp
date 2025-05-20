#define MC__PWCU_CHECK
//#define MC__PWCU_DEBUG
#define MC__SUPMODEL_TRACE

#include <iostream>

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

#include "pwcu.hpp"
#include "supmodel.hpp"
typedef mc::SupModel<mc::PWCU> SM;
typedef mc::SupVar<mc::PWCU> SV;

using namespace mc;

void
doxygen_supmodel
()
{
  size_t const N = 8;
  SM mod( 2 );
  SV X( mod, 0, I(1.,2.), N );
  SV Y( mod, 1, I(0.,1.), N );

  SV F = X*exp(X+pow(Y,2))-pow(Y,2);

  std::cout << "Superposition relaxation of f:\n" << F;

  std::vector<mc::PWCU> const& Fuest = F.uest();
  std::vector<mc::PWCU> const& Foest = F.oest();

  double const& Fl = F.l();
  double const& Fu = F.u();

  std::cout << "uest(1.5,0.5) = " << F.uval({{0,1.5},{1,0.5}}) << std::endl;
  std::cout << "oest(1.5,0.5) = " << F.oval({{0,1.5},{1,0.5}}) << std::endl;

  std::cout << std::scientific << std::setprecision(5) << std::endl;
  for( double x=1.; x<=2.; x+=1./20.-DBL_EPSILON*10. ){
   for( double y=0.; y<=1.; y+=1./20.-DBL_EPSILON*10. )
    std::cout << std::setw(13) << x
              << std::setw(13) << y
              << std::setw(13) << x*std::exp(x+std::pow(y,2))-std::pow(y,2)
              << std::setw(13) << F.uval({{0,x},{1,y}})
              << std::setw(13) << F.oval({{0,x},{1,y}})
              << std::endl; 
   std::cout << std::endl;
  }
  std::cout << std::endl;

  double const tol = 1e-5;
  for( size_t i=0; i<=N; ++i ){

   double x = i? 1.+i/(double)N*(2.-1.)-tol: 1.;

   for( size_t j=0; j<=N; ++j ){
    double y = j? 0.+j/(double)N*(1.-0.)-tol: 0.;
    std::cout << std::setw(13) << x
              << std::setw(13) << y
              << std::setw(13) << F.uval({{0,x},{1,y}})
              << std::setw(13) << F.oval({{0,x},{1,y}})
              << std::endl;
    if( j == N || j == 0 ) continue;
    
    y = 0.+j/(double)N*(1.-0.)+tol;
    std::cout << std::setw(13) << x
              << std::setw(13) << y
              << std::setw(13) << F.uval({{0,x},{1,y}})
              << std::setw(13) << F.oval({{0,x},{1,y}})
              << std::endl;
   }

   std::cout << std::endl;
   if( i == N || i == 0 ) continue;

   x = 1.+i/(double)N*(2.-1.)+tol;
   
   for( size_t j=0; j<=N; ++j ){
    double y = j? 0.+j/(double)N*(1.-0.)-tol: 0.;
    std::cout << std::setw(13) << x
              << std::setw(13) << y
              << std::setw(13) << F.uval({{0,x},{1,y}})
              << std::setw(13) << F.oval({{0,x},{1,y}})
              << std::endl;
    if( j == N || j == 0 ) continue;
    
    y = 0.+j/(double)N*(1.-0.)+tol;
    std::cout << std::setw(13) << x
              << std::setw(13) << y
              << std::setw(13) << F.uval({{0,x},{1,y}})
              << std::setw(13) << F.oval({{0,x},{1,y}})
              << std::endl;
   }

   std::cout << std::endl;
  }
}
/*
void
test_supmodel
()
{
  SM mod( 2 );
  
  SV C( 2. );
  SV X( mod, 0, I(0.,1.) );
  SV Y( mod, 1, mc::PWCU( 0., {0.2,0.6,0.2} ) );
  std::cout << mod;
  std::cout << X << Y;

  X -= 1.;
  std::cout << X;

  X += X;
  std::cout << X;

  X -= Y;
  std::cout << X;

  auto Z = exp(X);
  std::cout << Z;
}

const double XL   = -1.;	// <-- range lower bound
const double XU   =  2.;	// <-- range upper bound
template <class T>
T f_1d
( const T&x )
{
  //return -exp(x);
  //return x*exp(-pow(x,2)+1);
  return x*exp(-fabs(pow(x,3)));
}

void
test_supmodel_1d
()
{
  SM mod( 1 );

  //SV X( mod, 0, I(XL,XU) );
  SV X( mod, 0, mc::PWCU( XL, std::list<double>(8,(XU-XL)/8) ) );

  auto Z = f_1d(X);
  std::cout << Z;
 
  std::cout << std::scientific << std::setprecision(5) << std::endl;
  for( double x=XL; x<=XU; x+=(XU-XL)/100-DBL_EPSILON*10. )
    std::cout << std::setw(13) << x
              << std::setw(13) << f_1d(x)
              << std::setw(13) << Z.l({{0,x}})
              << std::setw(13) << Z.u({{0,x}}) 
              << std::endl; 
}

//const double XL   = -2.;	// <-- X range lower bound
//const double XU   =  1.5;	// <-- X range upper bound
//const double YL   = -0.3;	// <-- Y range lower bound
//const double YU   =  0.2;	// <-- Y range upper bound
//template <class T>
//T f_2d
//( T const& x, T const& y )
//{
//  //return exp(exp(x+y));
//  //return 1./exp(x+y);
//  //return tanh(x+y);
//  //return sqr(x+y);
//  //return sqr(x+y)-sqr(x-y);
//  return x*y;
//}

//const double XL   =  0.2;	// <-- X range lower bound
//const double XU   =  1.5;	// <-- X range upper bound
//const double YL   =  0.;	// <-- Y range lower bound
//const double YU   =  1.;	// <-- Y range upper bound
//template <class T>
//T f_2d
//( T const& x, T const& y )
//{
//  //return 1./(x+y);
//  //return log(x+y);
//  //return sqrt(x+y);
//  //return pow(x+y,0.3);
//  //return sqrt(log(x+y));
//  //return pow(x+y,3);
//  //return x/y;
//  //using std::min;
//  //return min(x-y,0.);
//  //return fabs(x-y);
//  //return xlog(x+y);
//  using std::max;
//  return max(x,y);
//}

//const double XL   =  -2.;	// <-- X range lower bound
//const double XU   =  -1.;	// <-- X range upper bound
//const double YL   =  -1.5;	// <-- Y range lower bound
//const double YU   =  -0.8;	// <-- Y range upper bound
//template <class T>
//T f_2d
//( T const& x, T const& y )
//{
//  //return 1./(x+y);
//  //return exp(1./(x+y));
//  //return pow(x+y,-3);
//  return tanh(x+y);
//}

//const double XL   =  -2.;	// <-- X range lower bound
//const double XU   =  1.;	// <-- X range upper bound
//const double YL   =  -1.5;	// <-- Y range lower bound
//const double YU   =  0.8;	// <-- Y range upper bound
//template <class T>
//T f_2d
//( T const& x, T const& y )
//{
//  return pow(x+y,3);
//}

//const double XL   =  1.;	// <-- X range lower bound
//const double XU   =  2.;	// <-- X range upper bound
//const double YL   =  0.;	// <-- Y range lower bound
//const double YU   =  1.;	// <-- Y range upper bound
//template <class T>
//T f_2d
//( T const& x, T const& y )
//{
//  return x*exp(x+sqr(y))-sqr(y);
//}

const double XL   = -1.;	// <-- X range lower bound
const double XU   =  1.;	// <-- X range upper bound
const double YL   =  0.5;	// <-- Y range lower bound
const double YU   =  2.;	// <-- Y range upper bound
template <class T>
T f_2d
( T const& x, T const& y )
{
  //return fabs(x);
  //return exp(-x/y);
  //return x*exp(-x/y);
  //return exp(-fabs(x)/y);
  //return fabs(x)*exp(-fabs(x)/y);
  using std::max;
  return max(fabs(x)*exp(-fabs(x)/y),0.);
}



void
test_supmodel_2d
()
{
  SM mod( 2 );
  mod.options.PROD_METH  = SM::Options::PARTIAL;//FULL;//LOG;//NONE;
  mod.options.PROD_CUT   = 0;
  mod.options.REF_WEIGHT = 0.5;
  //SV X( mod, 0, I(XL,XU) );
  //SV X( mod, 0, mc::PWCU( XL, {0.5*(XU-XL),0.5*(XU-XL)} ) );
  SV X( mod, 0, mc::PWCU( XL, std::list<double>(5, 0.2*(XU-XL)) ) );
  //SV Y( mod, 1, mc::PWCU( YL, {0.5*(YU-YL),0.5*(YU-YL)} ) );
  SV Y( mod, 1, mc::PWCU( YL, std::list<double>(5, 0.2*(YU-YL)) ) );

  auto Z = f_2d(X,Y);
  std::cout << Z;

  std::cout << std::scientific << std::setprecision(5) << std::endl;
  for( double x=XL; x<=XU; x+=(XU-XL)/20-DBL_EPSILON*10. ){
    for( double y=YL; y<=YU; y+=(YU-YL)/20-DBL_EPSILON*10. )
    std::cout << std::setw(13) << x
              << std::setw(13) << y
              << std::setw(13) << f_2d(x,y)
              << std::setw(13) << Z.l({{0,x},{1,y}})
              << std::setw(13) << Z.u({{0,x},{1,y}})
              << std::endl; 
   std::cout << std::endl;
 }
}
*/
void
test_linear
()
{
  size_t const N = 8;
  PWCU X( -1., 1., N, 1 );
  std::cout << "X:" << X << std::endl;

  PWCU Y = X;
  Y += X;
  std::cout << "Y+X:" << Y+X << std::endl;

  PWCU Z = X;
  std::cout << "Z+Y:" << Z+Y << std::endl;
  std::cout << "Z/2.-1:" << Z/2.-1. << std::endl;

  return;
}

void
test_exp
()
{
  size_t const N = 8;
  double XL = -1., XU = 1.;
  PWCU expX( XL, XU, N );
  std::cout << "X:" << expX << std::endl;

  auto const& f = [=]( const double& x ){ return std::exp(x); };
  auto const& df = [=]( const double& x ){ return std::exp(x); };

  std::cout << "exp(X):" << expX.compose( f, df, 1, 1, 1 ) << std::endl;

  std::cout << std::scientific << std::setprecision(5) << std::endl;
  for( double x=XL; x<=XU; x+=(XU-XL)/100-DBL_EPSILON*10. )
    std::cout << std::setw(13) << x
              << std::setw(13) << f(x)
              << std::setw(13) << expX.l(x)
              << std::setw(13) << expX.u(x)
              << std::endl;

  return;
}

void
test_log
()
{
  size_t const N = 8;
  double XL = 0.5, XU = 5.;
  PWCU logX( XL, XU, N );
  std::cout << "X:" << logX << std::endl;

  auto const& f = [=]( const double& x ){ return std::log(x); };
  auto const& df = [=]( const double& x ){ return 1./x; };

  std::cout << "log(X):" << logX.compose( f, df, 1, 0, 1 ) << std::endl;

  std::cout << std::scientific << std::setprecision(5) << std::endl;
  for( double x=XL; x<=XU; x+=(XU-XL)/100-DBL_EPSILON*10. )
    std::cout << std::setw(13) << x
              << std::setw(13) << f(x)
              << std::setw(13) << logX.l(x)
              << std::setw(13) << logX.u(x)
              << std::endl;

  return;
}

void
test_sqr
()
{
  size_t const N = 8;
  double XL = -1., XU = 1.;
  PWCU sqr1X( XL, XU, N );
  std::cout << "X:" << sqr1X << std::endl;

  PWCU sqr2X = sqr1X;

  auto const& f1  = [=]( const double& x ){ return std::pow(x>0?x:0,2); };
  auto const& df1 = [=]( const double& x ){ return 2*(x>0?x:0); };

  auto const& f2  = [=]( const double& x ){ return std::pow(x<0?x:0,2); };
  auto const& df2 = [=]( const double& x ){ return 2*(x<0?x:0); };
  
  std::cout << "sqr1(X):" << sqr1X.compose( f1, df1, 1, 1, 1 ) << std::endl;
  std::cout << "sqr2(X):" << sqr2X.compose( f2, df2, 1, 1, 0 ) << std::endl;
  std::cout << "sqr(X):" << sqr1X + sqr2X << std::endl;

  PWCU sqrX = sqr1X + sqr2X;

  std::cout << std::scientific << std::setprecision(5) << std::endl;
  for( double x=XL; x<=XU; x+=(XU-XL)/100-DBL_EPSILON*100. )
    std::cout << std::setw(13) << x
              << std::setw(13) << f1(x)+f2(x)
              << std::setw(13) << sqr1X.l(x)
              << std::setw(13) << sqr1X.u(x)
              << std::setw(13) << sqr2X.l(x)
              << std::setw(13) << sqr2X.u(x)
              << std::setw(13) << sqrX.l(x)
              << std::setw(13) << sqrX.u(x)
              << std::endl;

  return;
}

void
test_cub
()
{
  size_t const N = 8;
  double XL = -1., XU = 1.;
  PWCU cub1X( XL, XU, N );
  std::cout << "X:" << cub1X << std::endl;

  PWCU cub2X = cub1X;

  auto const& f1  = [=]( const double& x ){ return std::pow(x>0?x:0,3); };
  auto const& df1 = [=]( const double& x ){ return 3*std::pow(x>0?x:0,2); };

  auto const& f2  = [=]( const double& x ){ return std::pow(x<0?x:0,3); };
  auto const& df2 = [=]( const double& x ){ return 3*std::pow(x<0?x:0,2); };
  
  std::cout << "cub1(X):" << cub1X.compose( f1, df1, 1, 1, 1 ) << std::endl;
  std::cout << "cub2(X):" << cub2X.compose( f2, df2, 1, 0, 1 ) << std::endl;
  std::cout << "cub(X):" << cub1X + cub2X << std::endl;

  PWCU cubX = cub1X + cub2X;

  std::cout << std::scientific << std::setprecision(5) << std::endl;
  for( double x=XL; x<=XU; x+=(XU-XL)/100-DBL_EPSILON*100. )
    std::cout << std::setw(13) << x
              << std::setw(13) << f1(x)+f2(x)
              << std::setw(13) << cub1X.l(x)
              << std::setw(13) << cub1X.u(x)
              << std::setw(13) << cub2X.l(x)
              << std::setw(13) << cub2X.u(x)
              << std::setw(13) << cubX.l(x)
              << std::setw(13) << cubX.u(x)
              << std::endl;

  return;
}

int
main
( int argc, char* argv[] )
{
  doxygen_supmodel();
  //test_supmodel();
  //test_supmodel_1d();
  //test_supmodel_2d();
  //test_linear();
  //test_exp();
  //test_log();
  //test_sqr();
  //test_cub();

  return 0;
}
