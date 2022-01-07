#define TEST_DOXYGEN   // <-- select test example
#define USE_CHEB       // <-- whether to perform the decomposition in Chebyshev basis
////////////////////////////////////////////////////////////////////////

#include "interval.hpp"
typedef mc::Interval I;

#include "scmodel.hpp"
typedef mc::SCModel<I> SCM;
typedef mc::SCVar<I> SCV;

#include "mctime.hpp"
#include "squad.hpp"

using namespace std;
using namespace mc;

////////////////////////////////////////////////////////////////////////
int main()
////////////////////////////////////////////////////////////////////////
{  
  try{
    // Construction of multivariate polynomials
#if defined( TEST_DOXYGEN )
    SCM modSCM( 6 );
#if defined( USE_CHEB)
    modSCM.options.BASIS = SCM::Options::CHEB;
#else
    modSCM.options.BASIS = SCM::Options::MONOM;
#endif
    unsigned const NX = 3, NP = 2;
    SCV X[NX], P[NP];
    for( unsigned i=0; i<NX; i++ )
      X[i].set( &modSCM, i, I(-1,1) );
    P[0] = pow( X[0] + sqr( X[1] ) - 2 * X[2], 3 );
    P[1] = 2 * sqr( X[1] ) - 1;

#elif defined( TEST_NONDQUAR )
    SCM modSCM( 4 );
#if defined( USE_CHEB)
    modSCM.options.BASIS = SCM::Options::CHEB;
#else
    modSCM.options.BASIS = SCM::Options::MONOM;
#endif
    unsigned const NX = 32, NP = 1;
    SCV X[NX], P[NP];
    for( unsigned i=0; i<NX; i++ )
      X[i].set( &modSCM, i, I(-1,1) );
    // sum {i in 1..N-2} (x[i]+x[i+1]+x[N])^4 + (x[1]-x[2])^2 + (x[N-1]+x[N])^2;
    P[0] = pow( X[0]-X[1], 2 ) + pow( X[NX-2]+X[NX-1], 2 );
    for( unsigned i=0; i<NX-2; i++ )
      P[0] += pow( X[i]+X[i+1]+X[NX-1], 4 ); 

#elif defined( TEST_BROYDENBAND )
    SCM modSCM( 6 );
#if defined( USE_CHEB)
    modSCM.options.BASIS = SCM::Options::CHEB;
#else
    modSCM.options.BASIS = SCM::Options::MONOM;
#endif
    unsigned const NX = 60, NP = 1;
    SCV X[NX], P[NP];
    for( unsigned i=0; i<NX; i++ )
      X[i].set( &modSCM, i, I(-1,1) );
    // set J{i in 1..N} := {j in 1..N : j != i && max(1,i-5) ≤ j ≤ min(N,i+1)  max(1,i-ml) <= j <= min(N,i+mu) };
    // sum {i in 1..N} ( x[i]*(2+5*x[i]^2) + 1 - sum {j in J[i]} x[j]*(1+x[j]) );
    P[0] = 0.;
    for( int i=0; i<(int)NX; i++ ){
      SCV S = X[i]*(2+5*pow(X[i],2)) + 1;
      for( int j=0; j<(int)NX; j++ ){
        if( j == i || j < i-5 || j > i+1 ) continue;
        S -= X[j]*(1+X[j]);
      }
      P[0] += pow(S,2);
    }
#endif
    cout << "\nSparse multivariate polynomials:\n";
    for( unsigned i=0; i<NP; i++ )
      cout << P[i];

    // Quadratization of multivariate polynomials
    double tStart = mc::userclock();
    double viol = 0.;
    SQuad SQF;
    SQuad::options.ORDER = mc::SQuad::Options::DEC;
    SQuad::options.REDUC = true;
#if defined( USE_CHEB)
    SQuad::options.BASIS = mc::SQuad::Options::CHEB;
    for( unsigned i=0; i<NP; i++ )
      viol = SQF.process( P[i].coefmon(), mc::SQuad::Options::CHEB, i+1==NP?true:false );
#else
    SQuad::options.BASIS = mc::SQuad::Options::MONOM;
    for( unsigned i=0; i<NP; i++ )
      viol = SQF.process( P[i].coefmon(), mc::SQuad::Options::MONOM, i+1==NP?true:false );
#endif
    std::cout << "\nSparse quadratic forms: " << mc::userclock()-tStart << " CPU-sec\n"
              << "(discrepancy: " << viol << ")\n"
              << SQF;
  }

  catch( I::Exceptions &eObj ){
    cerr << "Error " << eObj.ierr()
         << " in natural interval extension:" << endl
	 << eObj.what() << endl
         << "Aborts." << endl;
    return eObj.ierr();
  }

  catch( SCM::Exceptions &eObj ){
    cerr << "Error " << eObj.ierr()
         << " in Chebyshev model arithmetic:" << endl
	 << eObj.what() << endl
         << "Aborts." << endl;
    return eObj.ierr();
  }

  catch( SQuad::Exceptions &eObj ){
    cerr << "Error " << eObj.ierr()
         << " in quadratization:" << endl
	 << eObj.what() << endl
         << "Aborts." << endl;
    return eObj.ierr();
  }

  return 0;
}
