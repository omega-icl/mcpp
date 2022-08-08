#define TEST_TUNCPHD30  // <-- select test example
#undef  USE_CHEB       // <-- whether to perform the decomposition in Chebyshev basis
#define USE_DAG        // <-- whether to define a DAG of the expressions
#define USE_OPTIM      // <-- whether to optimize the quadratic form
////////////////////////////////////////////////////////////////////////

#include "mctime.hpp"
#include "spoly.hpp"
#include "squad.hpp"
#if defined( USE_DAG )
  #include "ffunc.hpp"
  typedef mc::SPoly<mc::FFVar const*,mc::lt_FFVar> t_SPoly;
  typedef mc::SQuad<mc::FFVar const*,mc::lt_FFVar> t_SQuad;
#else
  typedef mc::SPoly<> t_SPoly;
  typedef mc::SQuad<> t_SQuad;
#endif

using namespace std;
using namespace mc;

////////////////////////////////////////////////////////////////////////
int main()
////////////////////////////////////////////////////////////////////////
{  
  try{
    // Construction of multivariate polynomials
#if defined( USE_CHEB )
  t_SPoly::options.BASIS = t_SPoly::Options::CHEB;
#else
  t_SPoly::options.BASIS = t_SPoly::Options::MONOM;
#endif
#if defined( TEST_DOXYGEN )
    //unsigned const NX = 1, NP = 1;
    unsigned const NX = 3, NP = 2;
    t_SPoly X[NX], P[NP];
#if defined( USE_DAG )
    FFGraph DAG;
    FFVar DAGX[NX];
    for( unsigned i=0; i<NX; i++ ) X[i].var( &DAGX[i].set( &DAG ) );
#else
    for( unsigned i=0; i<NX; i++ ) X[i].var( i );
#endif
    //P[0] = pow( X[0], 4 ) + pow( X[0], 6 );
    P[0] = pow( X[0] + sqr( X[1] ) - 2 * X[2], 3 );
    P[1] = 2 * sqr( X[1] ) - 1;

#elif defined( TEST_POP )
    unsigned const NX = 4, NP = 3;
    t_SPoly X[NX], P[NP];
#if defined( USE_DAG )
    FFGraph DAG;
    FFVar DAGX[NX];
    for( unsigned i=0; i<NX; i++ ) X[i].var( &DAGX[i].set( &DAG ) );
#else
    for( unsigned i=0; i<NX; i++ ) X[i].var( i );
#endif
    //P[0] = (X[0]*X[3])*X[1];
    //P[1] = (X[0]*X[3])*X[1]*X[2];
    P[0] = (X[0]*X[3])*(X[0]+X[1]+X[2])+X[2];
    P[1] = (X[0]*X[3])*X[1]*X[2]-25;
    P[2] = sqr(X[0])+sqr(X[1])+sqr(X[2])+sqr(X[3])-40;

#elif defined( TEST_TUNCPHD30 )
    unsigned const NX = 6, NP = 3;
    t_SPoly X[NX], P[NP];
#if defined( USE_DAG )
    FFGraph DAG;
    FFVar DAGX[NX];
    for( unsigned i=0; i<NX; i++ ) X[i].var( &DAGX[i].set( &DAG ) );
#else
    for( unsigned i=0; i<NX; i++ ) X[i].var( i );
#endif
    P[0] = 0.0204*X[0]*X[3]*(X[0]+X[1]+X[2]) + 0.0187*X[1]*X[2]*(X[0]+1.57*X[1]+X[3]) + 0.0607*X[0]*X[3]*X[4]*X[4]*(X[0]+X[1]+X[2]) + 0.0437*X[1]*X[2]*X[5]*X[5]*(X[0]+1.57*X[1]+X[3]);
    P[1] = 0.001*X[0]*X[1]*X[2]*X[3]*X[4]*X[5] - 2.07;
    P[2] = 0.00062*X[0]*X[3]*X[4]*X[4]*(X[0]+X[1]+X[2])+0.00058*X[1]*X[2]*X[5]*X[5]*(X[0]+1.57*X[1]+X[3])-1;

#elif defined( TEST_NONDQUAR )
    unsigned const NX = 8, NP = 1;
    t_SPoly X[NX], P[NP];
#if defined( USE_DAG )
    FFGraph DAG;
    FFVar DAGX[NX];
    for( unsigned i=0; i<NX; i++ ) X[i].var( &DAGX[i].set( &DAG ) );
#else
    for( unsigned i=0; i<NX; i++ ) X[i].var( i );
#endif
    // sum {i in 1..N-2} (x[i]+x[i+1]+x[N])^4 + (x[1]-x[2])^2 + (x[N-1]+x[N])^2;
    P[0] = pow( X[0]-X[1], 2 ) + pow( X[NX-2]+X[NX-1], 2 );
    for( unsigned i=0; i<NX-2; i++ )
      P[0] += pow( X[i]+X[i+1]+X[NX-1], 4 ); 

#elif defined( TEST_BROYDENBAND )
    unsigned const NX = 60, NP = 1;
    t_SPoly X[NX], P[NP];
#if defined( USE_DAG )
    FFGraph DAG;
    FFVar DAGX[NX];
    for( unsigned i=0; i<NX; i++ ) X[i].var( &DAGX[i].set( &DAG ) );
#else
    for( unsigned i=0; i<NX; i++ ) X[i].var( i );
#endif
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

    t_SQuad SQF;
    t_SQuad::options.ORDER = t_SQuad::Options::INC;
    t_SQuad::options.REDUC = false;
    t_SQuad::options.MIPOUTPUTFILE = "quad.lp";
#if defined( USE_CHEB )
    t_SQuad::options.BASIS = t_SQuad::Options::CHEB;
    viol = SQF.process( NP, P, &t_SPoly::mapmon, t_SQuad::Options::CHEB, true );
#else
    t_SQuad::options.BASIS = t_SQuad::Options::MONOM;
    viol = SQF.process( NP, P, &t_SPoly::mapmon, t_SQuad::Options::MONOM, true );
#endif
    std::cout << "\nSparse quadratic forms: " << mc::userclock()-tStart << " CPU-sec\n"
              << "(discrepancy: " << viol << ")\n"
              << SQF;
#if defined( USE_OPTIM )
    SQF.optimize( 2, true );
    viol = SQF.check( NP, P, &t_SPoly::mapmon, t_SQuad::Options::MONOM );
    std::cout << "\nSparse quadratic forms: " << mc::userclock()-tStart << " CPU-sec\n"
              << "(discrepancy: " << viol << ")\n"
              << SQF;
#endif
  }

  catch( t_SPoly::Exceptions &eObj ){
    cerr << "Error " << eObj.ierr()
         << " in polynomial expression:" << endl
	 << eObj.what() << endl
         << "Aborts." << endl;
    return eObj.ierr();
  }

  catch( t_SQuad::Exceptions &eObj ){
    cerr << "Error " << eObj.ierr()
         << " in quadratization:" << endl
	 << eObj.what() << endl
         << "Aborts." << endl;
    return eObj.ierr();
  }

  return 0;
}
