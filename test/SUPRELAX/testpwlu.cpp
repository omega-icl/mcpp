#define MC__PWLU_CHECK
//#define MC__PWLU_DEBUG
#define MC__SUPMODEL_TRACE

#include <iostream>
#include "pwlu.hpp"

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

#include "supmodel.hpp"
typedef mc::SupModel<mc::PWLU> SM;
typedef mc::SupVar<mc::PWLU> SV;

using namespace mc;

void
doxygen_supmodel
()
{
  SM mod( 2 );
  SV X( mod, 0, I(1.,2.), 8 );
  SV Y( mod, 1, I(0.,1.), 8 );

  SV F = X*exp(X+pow(Y,2))-pow(Y,2);

  std::cout << "Suposition relaxation of f:\n" << F;

  std::vector<mc::PWLU> const& Fuest = F.uest();
  std::vector<mc::PWLU> const& Foest = F.oest();

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
}

void
test_supmodel
()
{
  SM mod( 2 );
  
  SV C( 2. );
  SV X( mod, 0, I(0.,1.) );
  SV Y( mod, 1, mc::PWLU( 0., {0.2,0.6,0.2} ) );
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
/*
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
  SV X( mod, 0, mc::PWLU( XL, std::list<double>(8,(XU-XL)/8) ) );

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
*/
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
  //SV X( mod, 0, mc::PWLU( XL, {0.5*(XU-XL),0.5*(XU-XL)} ) );
  SV X( mod, 0, mc::PWLU( XL, std::list<double>(5, 0.2*(XU-XL)) ) );
  //SV Y( mod, 1, mc::PWLU( YL, {0.5*(YU-YL),0.5*(YU-YL)} ) );
  SV Y( mod, 1, mc::PWLU( YL, std::list<double>(5, 0.2*(YU-YL)) ) );

  auto Z = f_2d(X,Y);
  std::cout << Z;

  std::cout << std::scientific << std::setprecision(5) << std::endl;
  for( double x=XL; x<=XU; x+=(XU-XL)/20-DBL_EPSILON*10. ){
    for( double y=YL; y<=YU; y+=(YU-YL)/20-DBL_EPSILON*10. )
    std::cout << std::setw(13) << x
              << std::setw(13) << y
              << std::setw(13) << f_2d(x,y)
              << std::setw(13) << Z.uval({{0,x},{1,y}})
              << std::setw(13) << Z.oval({{0,x},{1,y}})
              << std::endl; 
   std::cout << std::endl;
 }
}

void
test_reduce
()
{
  PWLU Q;
  Q.options.reset();

  Q.xL() = -3.0000000000000000e+00;
  Q.dx() = { 3.5714285714285715e-01, 3.9285721258070888e-01, 3.4999993027643400e-01, 1.8125207858698206e-01, 2.1874792141301802e-01, 3.3333333333333331e-01, 4.1666666666666685e-02, 1.0077995278366680e-01, 2.7422004721633320e-01, 2.5000000000000000e-01, 1.2500000000000000e-01, 3.4972241417178418e-02, 3.4002775858282158e-01, 3.7500000000000000e-01, 1.2499622524063553e-01, 2.5000377475936447e-01, 2.7422004721633320e-01, 1.0077995278366680e-01, 4.1666666666666685e-02, 3.3333333333333331e-01, 2.1874792141301791e-01, 1.5625207858698209e-01, 2.5000000000000022e-02, 3.4999999999999998e-01, 3.9285708324239554e-01, 3.5714291675760457e-01 };
  Q.yL() = 7.1121599991958528e+03;
  Q.dy() = { -7.8542713315829869e+03, -7.8542713315829869e+03, -1.7454716427862743e+03, -1.7454716427862743e+03, -9.0369359989084478e+02, -1.6673694819625516e+02, -1.6673694819625516e+02, -1.6673694819625516e+02, -5.9826235145381204e+01, 8.2371153214518031e+00, 6.1329652304636841e+00, -3.6552526741069751e+00, -2.6492546004474962e+00, 4.4057752852412122e+00, -5.6709355792903091e+00, -6.2822504762894589e+00, 7.5949536536794454e+01, 7.5949536536794454e+01, 1.7715816574258304e+02, 1.7715816574258304e+02, 1.0783571940853521e+03, 1.0783571940853521e+03, 1.9206334744114920e+03, 1.9206334744114920e+03, 7.8542656554719788e+03, 7.8542656554719788e+03 };

  //Q.xL() = Q.yL() = -1;

  //Q.dx() = { 9.62963e-02, 1.03704e-01, 9.52381e-02, 1.98095e-01 };//, 1.95556e-01, 1.77778e-01, 1.33333e-01, 2.00000e-01, 2.00000e-01, 2.00000e-01, 2.00000e-01, 2.00000e-01 };
  //Q.dy() = { 3.00000e+00, 1.92000e+00, 1.92000e+00, 1.08000e+00 };//, 4.80000e-01, 1.20000e-01, 0.00000e+00, 4.00000e-02, 2.80000e-01, 7.60000e-01, 1.48000e+00, 2.44000e+00 };

  //Q.dx() = { 9.62963e-02, 1.03704e-01, 9.52381e-02, 1.04762e-01, 9.33333e-02, 1.95556e-01, 1.77778e-01, 1.33333e-01, 2.00000e-01, 2.00000e-01, 2.00000e-01, 2.00000e-01, 2.00000e-01 };
  //Q.dy() = { 3.00000e+00, 1.92000e+00, 1.92000e+00, 1.08000e+00, 1.08000e+00, 4.80000e-01, 1.20000e-01, 0.00000e+00, 4.00000e-02, 2.80000e-01, 7.60000e-01, 1.48000e+00, 2.44000e+00 };

  //Q.dx() = { 9.62963e-02, 1.03704e-01, 9.52381e-02, 1.04762e-01, 9.33333e-02, 1.06667e-01, 8.88889e-02, 1.11111e-01, 6.66667e-02, 1.33333e-01, 2.00000e-01, 2.00000e-01, 2.00000e-01, 2.00000e-01, 2.00000e-01 };
  //Q.dy() = { 3.00000e+00, 1.92000e+00, 1.92000e+00, 1.08000e+00, 1.08000e+00, 4.80000e-01, 4.80000e-01, 1.20000e-01, 1.20000e-01, 0.00000e+00, 4.00000e-02, 2.80000e-01, 7.60000e-01, 1.48000e+00, 2.44000e+00 };

  //Q.dx() = { 9.62963e-02, 1.98942e-01, 1.04762e-01, 9.33333e-02, 1.06667e-01, 8.88889e-02, 1.77778e-01, 1.33333e-01, 2.00000e-01, 2.00000e-01, 2.00000e-01, 2.00000e-01, 2.00000e-01 };
  //Q.dy() = { 3.00000e+00, 1.92000e+00, 1.08000e+00, 1.08000e+00, 4.80000e-01, 4.80000e-01, 1.20000e-01, 0.00000e+00, 4.00000e-02, 2.80000e-01, 7.60000e-01, 1.48000e+00, 2.44000e+00 };

  //Q.dx() = { 1.75, 0.583333, 1.66667, 0.666667, 1.33333 };
  //Q.dy() = { 0,    -2,       1,       -2,       1 };

  //Q.dx() = { 1.75, 0.25, 1, 1, 0.666667, 1.33333 };
  //Q.dy() = { 0,    -2,   0, 1, -2,       1 };

  //Q.dx() = { 1, 0.5, 0.5, 1, 1, 0.666667, 1.33333 };
  //Q.dy() = { 0, 1,   -2,  0, 1, -2,       1 };

  //Q.dx() = { 1, 0.5, 0.5, 1, 0.5, 0.5, 1, 1 };
  //Q.dy() = { 0, 1,  -2,   0, 1,   -2,  0, 1 };

  //Q.dy() = { 0, 1, 2, 1 };
  //Q.dx() = { 1, 1, 1 };
  //Q.dy() = { 1, -2, 0 };
  //Q.dy() = { -1, 2, 0 };
  std::cout << Q << std::endl;
  
  PWLU Qu = Q, Qo = Q;
  //Qu.reduce( 1, 4 );
  //std::cout << Qu << std::endl;
  Qo.reduce( 0, 16 );
  std::cout << Qo << std::endl;

  std::cout << std::scientific << std::setprecision(5) << std::endl;
  for( double x=Q.xL(); x<=Q.xU(); x+=(Q.xU()-Q.xL())/300-DBL_EPSILON*10. )
    std::cout << std::setw(13) << x
              << std::setw(13) << Q.val(x)
              //<< std::setw(13) << Qu.val(x)
              << std::setw(13) << Qo.val(x)
              << std::endl;

  return;
}

void
test_merge
()
{
  PWLU Q;
  Q.options.DISPNUM=16;
  Q.options.BKPTATOL=1e-1;
  Q.dx() = { 1e-2, 1e-2, 1, 1e-2, 1e-2, 1e-2, 2e-2, 1, 1e-2, 1e-2 };
  Q.dy() = { 0, 4, 1, 2, 0, 2, 0, 1, 0, 2 };
  std::cout << Q << std::endl;
  
  PWLU Qu = Q, Qo = Q;
  Qu.clean( 1 ).merge( 1 );
  std::cout << Qu << std::endl;
  Qo.clean( 0 ).merge( 0 );
  std::cout << Qo << std::endl;

  std::cout << std::scientific << std::setprecision(5) << std::endl;
  for( double x=Q.xL(); x<=Q.xU(); x+=(Q.xU()-Q.xL())/500-DBL_EPSILON*10. )
    std::cout << std::setw(13) << x
              << std::setw(13) << Q.val(x)
              << std::setw(13) << Qu.val(x)
              << std::setw(13) << Qo.val(x)
              << std::endl;

  return;
}

void
test_linear
()
{
  PWLU X( -1., 1. );
  std::cout << "X:" << X << std::endl;

  PWLU Y = X;
  Y.insert( 0.5 );
  Y.insert( 0.55 );
  Y.insert( 0.85 );
  std::cout << "Y=X:" << Y << std::endl;

  Y += X;
  std::cout << "Y+X:" << Y+X << std::endl;

  PWLU Z = X;
  Z.insert( {0.4,0.45,0.6} );
  std::cout << "Z=X:" << Z << std::endl;

  std::cout << "Z+Y:" << Z+Y << std::endl;

  std::cout << "Z/2.-1:" << Z/2.-1. << std::endl;

  return;
}

void
test_sqr
()
{
  PWLU X( -1., 1. );
  std::cout << "X:" << X << std::endl;

  auto const& f = [=]( const double& x ){ return x*x; };
  auto const& df = [=]( const double& x ){ return 2*x; };

  PWLU sqrX_u = X; //sqrX_u.insert( 0.5 );
  PWLU sqrX_o = X; //sqrX_o.insert( 0.5 );

  std::cout << "sqr(X) oest:" << sqrX_o.max( 0. ).compose( f, df, 0, 1 ) << std::endl;
  std::cout << "sqr(X) uest:" << sqrX_u.max( 0. ).compose( f, df, 1, 1 ) << std::endl;

  std::cout << std::scientific << std::setprecision(5) << std::endl;
  for( double x=X.l(); x<=X.u(); x+=(X.u()-X.l())/200-DBL_EPSILON*10. )
    std::cout << std::setw(13) << x
              << std::setw(13) << f(x>0?x:0)
              << std::setw(13) << sqrX_u.val(x)
              << std::setw(13) << sqrX_o.val(x)
              << std::endl;

  return;
}

void
test_sqr2
()
{
  PWLU X( -1., 1. );
  std::cout << "X:" << X << std::endl;

  auto const& f = [=]( const double& x ){ return x*x; };
  auto const& df = [=]( const double& x ){ return 2*x; };

  X.insert({-0.5,0.,0.5});
  //X.insert({-0.8,-0.6,-0.4,-0.2,0.,0.2,0.4,0.6,0.8});
  PWLU sqrX_u = X; //sqrX_u.insert( 0.5 );
  PWLU sqrX_o = X; //sqrX_o.insert( 0.5 );
  
  std::cout << "sqr(X) oest:" << sqrX_o.compose( f, df, 0, 1 ) << std::endl;
  std::cout << "sqr(X) uest:" << sqrX_u.compose( f, df, 1, 1 ) << std::endl;

//  std::cout << "max(sqr(X),0.2) oest:" << sqrX_o.max( 0.25 ) << std::endl;
//  std::cout << "max(sqr(X),0.2) uest:" << sqrX_u.max( 0.25 ) << std::endl;

//  std::cout << "sqr(X) oest:" << sqrX_o.reduce(0) << std::endl;
//  std::cout << "sqr(X) uest:" << sqrX_u.reduce(1) << std::endl;

  std::cout << "min(sqr(X),0.2) oest:" << sqrX_o.min( 0.4 ).reduce(0,20) << std::endl;
  std::cout << "min(sqr(X),0.2) uest:" << sqrX_u.min( 0.4 ).reduce(1,20) << std::endl;

  std::cout << std::scientific << std::setprecision(5) << std::endl;
  for( double x=X.l(); x<=X.u(); x+=(X.u()-X.l())/500-DBL_EPSILON*10. )
    std::cout << std::setw(13) << x
              << std::setw(13) << sqrX_u.val(x)
              << std::setw(13) << sqrX_o.val(x)
              << std::endl;

  return;
}

void
test_cub
()
{
  PWLU X( -1., 1. );
  std::cout << "X:" << X << std::endl;

  auto const& f1  = [=]( const double& x ){ return std::pow(x>0?x:0,3); };
  auto const& df1 = [=]( const double& x ){ return 3*std::pow(x>0?x:0,2); };

  auto const& f2  = [=]( const double& x ){ return std::pow(x<0?x:0,3); };
  auto const& df2 = [=]( const double& x ){ return 3*std::pow(x<0?x:0,2); };

  //X.insert({-0.5,0.,0.5});
  X.insert({-0.8,-0.6,-0.4,-0.2,0.,0.2,0.4,0.6,0.8});
  PWLU cub1X_u = X;
  PWLU cub1X_o = X;
  PWLU cub2X_u = X;
  PWLU cub2X_o = X;
  
  std::cout << "cub1(X) oest:" << cub1X_o.compose( f1, df1, 0, 1 ) << std::endl;
  std::cout << "cub1(X) uest:" << cub1X_u.compose( f1, df1, 1, 1 ) << std::endl;

  std::cout << "cub2(X) oest:" << cub2X_o.compose( f2, df2, 0, 0 ) << std::endl;
  std::cout << "cub2(X) uest:" << cub2X_u.compose( f2, df2, 1, 0 ) << std::endl;

  std::cout << "cub(X) oest:" << cub1X_o + cub2X_o << std::endl;
  std::cout << "cub(X) uest:" << cub1X_u + cub2X_u << std::endl;

  PWLU cubX_u = cub1X_u + cub2X_u;
  PWLU cubX_o = cub1X_o + cub2X_o;
  cubX_u.reduce( 1, 5 );
  cubX_o.reduce( 0, 5 );

  std::cout << std::scientific << std::setprecision(5) << std::endl;
  for( double x=X.l(); x<=X.u(); x+=(X.u()-X.l())/500-DBL_EPSILON*10. )
    std::cout << std::setw(13) << x
              << std::setw(13) << x*x*x
              << std::setw(13) << cub1X_u.val(x)
              << std::setw(13) << cub1X_o.val(x)
              << std::setw(13) << cub2X_u.val(x)
              << std::setw(13) << cub2X_o.val(x)
              << std::setw(13) << cub1X_u.val(x) + cub2X_u.val(x)
              << std::setw(13) << cub1X_o.val(x) + cub2X_o.val(x)
              << std::setw(13) << cubX_u.val(x)
              << std::setw(13) << cubX_o.val(x)
              << std::endl;

  return;
}

int
main
( int argc, char* argv[] )
{
  //doxygen_supmodel();
  //test_supmodel();
  //test_supmodel_1d();
  //test_supmodel_2d();
  test_reduce();
  //test_merge();
  //test_linear();
  //test_sqr();
  //test_sqr2();
  //test_cub();

  return 0;
}
