#define TEST_1         // <-- specify to evaluate via a DAG of the function
#define  USE_DAG        // <-- specify to evaluate via a DAG of the function
#define SAVE_RESULTS   // <-- specify whether to save results to file
const int NRAD = 30;	// <-- select radius discretization here
const int NPHI = 400;	// <-- select angle discretization here
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

#include "ellimage.hpp"
typedef mc::Ellipsoid   E     ;
typedef mc::EllImg<I>   EI    ;
typedef mc::EllVar<I>   EV    ;
typedef CPPL::dcovector dcv   ;
typedef CPPL::dsymatrix dsm   ;
typedef CPPL::dgematrix dgm   ;

#ifdef USE_DAG
 #include "ffunc.hpp"
#endif

////////////////////////////////////////////////////////////////////////

#if defined( TEST_0 ) 
const int NX = 2;	
const int NF = 2;	
dcv q0(NX);
dsm Q0(NX);
template <class T>
T f
( const T*x, const long i )
{
 using mc::sqr;
  switch( i ){
    case 0:  return x[0];
    case 1:  return erf(x[0]);
    //case 1:  return (sqr(x[0]+x[1])-sqr(x[0]-x[1]))/4.;
    default: throw std::runtime_error("invalid index");
  }
}

#elif defined( TEST_1 ) 
const int NX = 2;	
const int NF = 2;	
dcv q0(NX);
dsm Q0(NX);
template <class T>
T f
( const T*x, const long i )
{
  switch( i ){
    case 0:  return x[0]*x[1];
    case 1:  return x[1]*x[0];
    default: throw std::runtime_error("invalid index");
  }
}

#elif defined( TEST_2)
const int NX = 2;	
const int NF = 2;	
dcv q0(NX);
dsm Q0(NX);
template <class T>
T f
( const T*x, const long i )
{
  using mc::sqr;
  switch( i ){
    //case 0:  return pow(x[0]*x[1],2) + log(x[1]);
    case 0:  return sqrt(x[0]+x[1]) + x[0]*x[1];
    case 1:  return pow(x[0]-x[1],2) + 3.*x[1];
    default: throw std::runtime_error("invalid index");
  }
}

#elif defined( TEST_3)
const int NX = 2;	
const int NF = 2;	
dcv q0(NX);
dsm Q0(NX);
template <class T>
T f
( const T*x, const long i )
{
  switch( i ){
    case 0:  return log(x[0]) + pow(x[1],2);
    case 1:  return sin(x[0]) - cos(x[1]);
    default: throw std::runtime_error("invalid index");
  }
}

#elif defined( TEST_4 )
const int NX = 2;	
const int NF = 2;	
dcv q0(NX);
dsm Q0(NX);
template <class T>
T f
( const T*x, const long i )
{
  switch( i ){
    case 0:  return  pow(x[0]*x[1],2) + log(x[1]);
    case 1:  return  sin(x[0]) * cos(x[1]);
    default: throw std::runtime_error("invalid index");
  }
}

#elif defined( TEST_5 )
const int NX = 2;	
const int NF = 2;	
dcv q0(NX);
dsm Q0(NX);
template <class T>
T f
( const T*x, const long i )
{
  using mc::sqr;
  switch( i ){
    case 0:  return  x[0] + x[1];
    case 1:  return  x[0] - x[1];
    default: throw std::runtime_error("invalid index");
  }
}

#endif

////////////////////////////////////////////////////////////////////////

void test_ellintersection_ia()
{
  const unsigned n=2, m=4;
  dgm A(m,n); dcv b(m);
  A(0,0) = -1.; A(0,1) =  0.; b(0) = 0.; // x >= 0
  A(1,0) =  1.; A(1,1) =  0.; b(1) = 1.; // x <= 1
  A(2,0) =  2.; A(2,1) = -1.; b(2) = 1.; // y >= 2x-1 
  A(3,0) = -2.; A(3,1) =  1.; b(3) = 1.; // y <= 2x+1 

  std::vector< std::pair<dcv,double> > HP;
  for( unsigned i=0; i<m; i++ )
    HP.push_back( std::make_pair( t(A.row(i)), b(i) ) );

  mc::Ellipsoid Eia( mc::ellintersection_ia( HP, 1e-5 ) );

  std::pair<CPPL::dcovector,CPPL::dgematrix> res = Eia.eigQ();
  std::cout << "Eigenvalues:\n" << res.first << std::endl;
  std::cout << "Eigenvectors:\n" << res.second << std::endl; 
  std::cout << "Square Root:\n" << Eia.sqrtQ() << std::endl;
}

void test_intersection_ea()
{
  const unsigned n=2;
  //center         shape
  dcv q1(n);      dsm Q1(n);
  q1(0) = 0.;     Q1(0,0) = 1.;
  q1(1) = 0.;     Q1(1,0) = 0.; Q1(1,1) = 1.; 
  mc::Ellipsoid E1(Q1,q1);
  std::cout << "E1:" << E1 << std::endl;

  //center         shape
  dcv q2(n);      dsm Q2(n);
  q2(0) = 1.;     Q2(0,0) = 1.;
  q2(1) = 1.;     Q2(1,0) = 0.; Q2(1,1) = 1.; 
  mc::Ellipsoid E2(Q2,q2);
  std::cout << "E2:" << E2 << std::endl;

  mc::Ellipsoid Eea1( mc::intersection_ea( E1, E2 ) );
  std::cout << "E1_inter_E2:" << Eea1 << std::endl;

  //hyperplane x1+x2 <= 0
  dcv a1(n); double b1(0.);
  a1(0) = 1.;   a1(1) = 1.;
  std::pair<dcv,double> HP1( a1, b1 );
  std::cout << "HP1:\n" << a1 << " <= " << b1 << std::endl << std::endl;

  mc::Ellipsoid Eea2( mc::intersection_ea( E1, HP1 ) );
  std::cout << "E1_inter_HP1:" << Eea2 << std::endl;
}

////////////////////////////////////////////////////////////////////////

int main()
{
  //test_ellintersection_ia();
  //test_intersection_ea();
  //return 0;

  #if defined( TEST_0 ) || defined( TEST_1 ) || defined( TEST_2 ) || defined( TEST_3 ) || defined( TEST_4 ) || defined( TEST_5 )
  //center         shape
//  q0(0) = 1.;     Q0(0,0) = 3.;
//  q0(1) = 1.;     Q0(1,0) = 2.; Q0(1,1) = 3.; 
  q0(0) = 3.;     Q0(0,0) = 2.;
  q0(1) = 4.;     Q0(1,0) = 1.; Q0(1,1) = 2.; 
  #endif
  mc::Ellipsoid EX(Q0,q0);
  
#ifdef USE_DAG
  mc::FFGraph FF;
  mc::FFVar X[NX];
  for( unsigned int i=0; i<NX; i++ ) X[i].set( &FF );
  mc::FFVar F[NF] = { f( X, 0 ), f( X, 1 ) };
#endif
  
  // Construct ellipsoidal image
  EI img;
  E::options.PSDCHK = false;
  img.options.PREALLOC = 0;
  img.options.DCPROD_USE = false;

  // Set independent variables in the ellipsoidal model
  img.set( EX.Q(), EX.c() );
  EV EVX[NX];
  for( long i=0; i<NX; ++i ) EVX[i].set( img, i ); 

  // Evaluate dependent variables
  EV EVF[NF];
#ifdef USE_DAG
  FF.eval( NF, F, EVF, NX, X, EVX );
#else
  for( long i=0; i<NF; ++i ) EVF[i] = f( EVX, i );
#endif

  // Print lifted ellipsoid
  std::cout << "\nLifted ellipsoid:\n" << img << std::endl;
 
  // Project lifted ellipsoidal image
  E EF = img.get( NF, EVF ) ;
  std::cout <<"\nEllipsoid image:\n" << EF << std::endl;

#ifdef SAVE_RESULTS
  // Parameters for sampling
  I      radius    = I(0.,1.);
  I      phi       = I(0.,2.*mc::PI);
  double discr     = (mc::Op<I>::u(radius) - mc::Op<I>::l(radius) )/(double ) NRAD;
  double discphi   = (mc::Op<I>::u(phi)    - mc::Op<I>::l(phi)    )/(double ) NPHI;
  double phiEr_sam = mc::Op<I>::l(phi);

  std::ostringstream name1;
  name1 << "Fsam.out" ; 
  std::ofstream fout(name1.str().c_str() , std::ios_base::out );
  fout << std::scientific << std::setprecision(8) << std::right;

  std::ostringstream name2;
  name2 << "Esam.out"; 
  std::ofstream eout(name2.str().c_str() , std::ios_base::out );
  eout << std::scientific << std::setprecision(8) << std::right;    

  CPPL::dcovector xy_vec(NF), z_sam(NF), E_sam(NF), Fimg(NF) ;
  while(phiEr_sam <= mc::Op<I>::u(phi)){
    xy_vec(0) = mc::Op<I>::u(radius) * cos(phiEr_sam);
    xy_vec(1) = mc::Op<I>::u(radius) * sin(phiEr_sam);
    E_sam     = EF.c() + EF.sqrtQ()*xy_vec;
    //print
    eout  <<   E_sam(0) << '\t' <<  E_sam(1) << std::endl;

    double  rad_sam = mc::Op<I>::u(radius);
    while( rad_sam >= mc::Op<I>::l(radius) +discr ){	
      double  phi_sam = mc::Op<I>::l(phi);
      while(phi_sam <= mc::Op<I>::u(phi) ){
        xy_vec(0) = rad_sam * cos(phi_sam);
        xy_vec(1) = rad_sam * sin(phi_sam);
        z_sam     = img.c() + img.sqrtQ()*xy_vec;
        // Evaluate dependent variables
        double z_d[2] ={z_sam(0),z_sam(1)};
        for( long i=0; i<NX; ++i ) Fimg(i) = f( z_d, i );
        //print
        if( phiEr_sam == mc::Op<I>::l(phi) ) fout <<  Fimg(0) << '\t' << Fimg(1) << std::endl;
        phi_sam += discphi;
      }
      rad_sam -= discr;
    }
    // Evaluate at center
    double z_c[2] ={img.c()(0), img.c()(1)};
    for( long i=0; i<NX; ++i ) Fimg(i) = f( z_c , i );
    //print evaluation at center
    if( phiEr_sam == mc::Op<I>::l(phi) ) fout <<  Fimg(0) << '\t' << Fimg(1) << std::endl;

    phiEr_sam += discphi;
  }
#endif

  // Construct interval enclosure
  I IX[NX], IF[NF];
  for( long i=0; i<NX; ++i ) IX[i] = I(img.l(i),img.u(i)); 
#ifdef USE_DAG
  FF.eval( NF, F, IF, NX, X, IX );
#else
  for( long i=0; i<NF; ++i ) IF[i] = f( IX, i );
#endif
  std::cout << "\nInterval enclosure:\n";
  for( long i=0; i<NF; ++i ) std::cout << "f[" << i << "] in " << IF[i] << std::endl;

#ifdef SAVE_RESULTS
  std::ostringstream name3;
  name3 << "Isam.out" ; 
  std::ofstream Iout(name3.str().c_str() , std::ios_base::out );
  Iout << std::scientific << std::setprecision(8) << std::right;
  Iout << mc::Op<I>::l(IF[0]) << "  " << mc::Op<I>::l(IF[1]) <<std::endl;
  Iout << mc::Op<I>::u(IF[0]) << "  " << mc::Op<I>::l(IF[1]) <<std::endl;
  Iout << mc::Op<I>::u(IF[0]) << "  " << mc::Op<I>::u(IF[1]) <<std::endl;
  Iout << mc::Op<I>::l(IF[0]) << "  " << mc::Op<I>::u(IF[1]) <<std::endl;
  Iout << mc::Op<I>::l(IF[0]) << "  " << mc::Op<I>::l(IF[1]) <<std::endl;
#endif

  return 0;
}





