////////////////////////////////////////////////////////////////////////
#define TEST_2		// <-- select model here
#define  SAVE_RESULTS      // <-- specify whether to save results to file for plotting 
unsigned long samples_r   = 30;    // <-- select number of samples keep one order of magnitude less than phi
unsigned long samples_phi = 400;   // <-- select number of samples
#undef USE_PROFIL	// <-- specify to use PROFIL for interval arithmetic
#undef USE_FILIB	// <-- specify to use FILIB++ for interval arithmetic
////////////////////////////////////////////////////////////////////////

#define _DEBUG
#define _TRACE
#include "ffunc.hpp"
#include <fstream>
#include <iomanip>
#include <sys/time.h>

#include "ellimage.hpp"

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

#include "ffunc.hpp"
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <sys/time.h>

#include "ellimage.hpp"
#include "interval.hpp"
typedef mc::Interval    I     ;
typedef mc::Ellipsoid   E     ;
typedef mc::EllImg<I>   EI    ;
typedef mc::EllVar<I>   EV    ;
typedef CPPL::dcovector dcv   ;
typedef CPPL::dsymatrix dsm   ;
typedef CPPL::dgematrix dgm   ;


#if defined( TEST_1 ) 
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
    case 1:  return sqr(x[0]-x[1]) + 3*x[1];
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
  using mc::sqr;
  switch( i ){
    case 0:  return log(x[0]) + sqr(x[1]);
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
  using mc::sqr;
  switch( i ){
    case 0:  return  sqr(x[0]*x[1]) + log(x[1]);
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

int main()
{
  //test_ellintersection_ia();
  //test_intersection_ea();
  //return 0;

  #if defined( TEST_1 ) || defined( TEST_2 ) || defined( TEST_3 ) || defined( TEST_4 ) || defined( TEST_5 )
  //center         shape
  q0(0) = 3.;     Q0(0,0) = 2.;
  q0(1) = 4.;     Q0(1,0) = 1.; Q0(1,1) = 2.; 
  #endif
  
  // Construct ellipsoidal domain
  EI Ex;
  E::options.PSDCHK = false;
  Ex.set( Q0, q0 );//, depmap );
  Ex.options.PREALLOC = 0;
  Ex.options.CHEBUSE  = false;

  I Ix[NX], If[NF];
  for( long i=0; i<NX; ++i ) Ix[i] = I(Ex.l(i),Ex.u(i)); 
  for( long i=0; i<NF; ++i ) If[i] = f( Ix, i );
  std::cout << "\nInterval enclosure:\n";
  for( long i=0; i<NF; ++i ) std::cout << "f[" << i << "] in " << If[i] << std::endl;

  // Set independent variables in the ellipsoidal model
  EV X[NX];
  for( long i=0; i<NX; ++i ) X[i].set( Ex, i ); 

  // Evaluate dependent variables
  EV F[NF];
  for( long i=0; i<NF; ++i ) F[i] = f( X, i );

  // Print lifted ellipsoid
  std::cout << "\nLifted ellipsoidal image:\n";
  Ex.output();
 
  // Project lifted ellipsoidal image
  //EV XF[NX+NF] = { EX[0], EX[1], EF[0] };
  //EI Ef = Ex.get( NX+NF, XF ) ;
  EI Ef = Ex.get( NF, F ) ;
  std::cout <<"\nEllipsoidal Image:\n" << Ef << std::endl;

  #ifdef SAVE_RESULTS
  // Parameters for sampling
  I      radius    = I(0.,1.);
  I      phi       = I(0.,2.*mc::PI);
  double discr     = (mc::Op<I>::u(radius) - mc::Op<I>::l(radius) )/(double ) samples_r;
  double discphi   = (mc::Op<I>::u(phi)    - mc::Op<I>::l(phi)    )/(double ) samples_phi;
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
    E_sam     = Ef.c() + Ef.sqrtQ()*xy_vec;
    //print
    eout  <<   E_sam(0) << '\t' <<  E_sam(1) << std::endl;

    double  rad_sam = mc::Op<I>::u(radius);
    while( rad_sam >= mc::Op<I>::l(radius) +discr ){	
      double  phi_sam = mc::Op<I>::l(phi);
      while(phi_sam <= mc::Op<I>::u(phi) ){
        xy_vec(0) = rad_sam * cos(phi_sam);
        xy_vec(1) = rad_sam * sin(phi_sam);
        z_sam     = Ex.c() + Ex.sqrtQ()*xy_vec;
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
    double z_c[2] ={Ex.c()(0), Ex.c()(1)};
    for( long i=0; i<NX; ++i ) Fimg(i) = f( z_c , i );
    //print evaluation at center
    if( phiEr_sam == mc::Op<I>::l(phi) ) fout <<  Fimg(0) << '\t' << Fimg(1) << std::endl;

    phiEr_sam += discphi;
  }

  std::ostringstream name3;
  name3 << "Isam.out" ; 
  std::ofstream Iout(name3.str().c_str() , std::ios_base::out );
  Iout << std::scientific << std::setprecision(8) << std::right;
  Iout << mc::Op<I>::l(If[0]) << "  " << mc::Op<I>::l(If[1]) <<std::endl;
  Iout << mc::Op<I>::u(If[0]) << "  " << mc::Op<I>::l(If[1]) <<std::endl;
  Iout << mc::Op<I>::u(If[0]) << "  " << mc::Op<I>::u(If[1]) <<std::endl;
  Iout << mc::Op<I>::l(If[0]) << "  " << mc::Op<I>::u(If[1]) <<std::endl;
  Iout << mc::Op<I>::l(If[0]) << "  " << mc::Op<I>::l(If[1]) <<std::endl;
  #endif
  
  
  return 0;
}





