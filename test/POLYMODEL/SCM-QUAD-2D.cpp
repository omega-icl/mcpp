const int NTE = 6;      // <-- select Taylor expansion order here
#undef USE_PROFIL       // <-- specify to use PROFIL for interval arithmetic
#undef USE_FILIB        // <-- specify to use FILIB++ for interval arithmetic
#undef  MC__SCMODEL_DEBUG_COMPOSITION
#undef  MC__SQUAD_DEBUG_PRODMON
#undef  MC__SQUAD_DEBUG_DECOMP
#undef  MC__SQUAD_DEBUG_REDUC
#define MC__SQUAD_DEBUG_SEPARATE
////////////////////////////////////////////////////////////////////////

#include <fstream>
#undef MC__SCMODEL_DEBUG_SPROD
#include <iomanip>

#include "mctime.hpp"
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
#include "scmodel.hpp"
typedef mc::SCModel<I> SCM;
typedef mc::SCVar<I> SCV;
#include "squad.hpp"

using namespace std;
using namespace mc;

////////////////////////////////////////////////////////////////////////

int main()
{ 
  SCM modSCM( NTE );
  modSCM.options.BOUNDER_TYPE = SCM::Options::LSB;
  modSCM.options.MIXED_IA = true;//false;

  double tStart = mc::userclock();

  unsigned const NX = 20, NF = 1;
  SCV X[NX], F[NF];
  for( unsigned i=0; i<NX; i++ ) X[i].set( &modSCM, i, I(-1,1) );

  // sum {i in 1..N-2} (x[i]+x[i+1]+x[N])^4 + (x[1]-x[2])^2 + (x[N-1]+x[N])^2;
  //F[0] = pow( X[0]-X[1], 2 ) + pow( X[NX-2]+X[NX-1], 2 );
  //for( unsigned i=0; i<NX-2; i++ ) F[0] += pow( X[i]+X[i+1]+X[NX-1], 4 ); 

  // set J{i in 1..N} := {j in 1..N : j != i && max(1,i-5) ≤ j ≤ min(N,i+1)  max(1,i-ml) <= j <= min(N,i+mu) };
  // sum {i in 1..N} ( x[i]*(2+5*x[i]^2) + 1 - sum {j in J[i]} x[j]*(1+x[j]) );
  F[0] = 0.;
  for( int i=0; i<(int)NX; i++ ){
    SCV S = X[i]*(2+5*pow(X[i],2)) + 1;
    for( int j=0; j<(int)NX; j++ ){
      if( j == i || j < i-5 || j > i+1 ) continue;
      S -= X[j]*(1+X[j]);
    }
    F[0] += pow(S,2);
  }

//  unsigned const NX = 3, NF = 2;
//  SCV X[NX], F[NF];
//  for( unsigned i=0; i<NX; i++ ) X[i].set( &modSCM, i, I(-1,1) );
//  F[0] = pow( X[0], 6 );
//  F[1] = pow( X[0], 3 ) * X[1];
//  //F[1] = pow( X[0] + sqr( X[1] ) - 2 * X[2], 3 );
//  //F[0] = pow( X[1] - 2, 3 ) - 1;
//  //F[1] = exp( X[0] * X[1] ) + 2 * sqr( X[2] ) - 1;
//  //F[0] = sqr( X[0] ) + pow( X[1], 3 );
  std::cout << "\nSparse Chebyshev model:" << mc::userclock()-tStart << " CPU-sec\n";
  for( unsigned i=0; i<NF; i++ ) std::cout << F[i];

  SQuad quad;
  quad.options.REDUC = true;
  quad.options.ORDER = SQuad::Options::DEC;//INC;// 
  quad.options.BASIS = SQuad::Options::MONOM;//CHEB;//
  tStart = mc::userclock();
  for( unsigned i=0; i<NF; i++ ){
    //auto coefmon = F[i].coefmon();
    //std::cout << F[i].display( coefmon, SQuad::Options::CHEB );
    //double viol = quad.process( coefmon, SQuad::Options::CHEB, true );
    auto coefmon = F[i].to_monomial( true );//false );//true );
    std::cout << F[i].display( coefmon, SQuad::Options::MONOM );
    double viol = quad.process( coefmon, SQuad::Options::MONOM, true );
    std::cout << "\nSparse quadratic form: " << mc::userclock()-tStart << " CPU-sec\n"
              << "(violation: " << viol << ")\n"
              << quad;
  }
/*
  auto&& quadlist = quad.separate( quad.MatFct().back() );
  for( auto const& quadterm : quadlist )
    auto&& eigendec = quad.factorize( quadterm );
*/            
  return 0;
} 

