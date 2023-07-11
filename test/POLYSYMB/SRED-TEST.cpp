#define TEST_QUAD1 // <-- select test example
#undef USE_DAG        // <-- select test example
////////////////////////////////////////////////////////////////////////

#include "mctime.hpp"
#include "spoly.hpp"
#include "sred.hpp"

#if defined( USE_DAG )
  #include "ffunc.hpp"
  typedef mc::SPoly<mc::FFVar const*,mc::lt_FFVar> t_SPoly;
  typedef mc::SRed<mc::FFVar const*,mc::lt_FFVar> t_SRed;
#else
  typedef mc::SPoly<> t_SPoly;
  typedef mc::SRed<> t_SRed;
#endif

using namespace std;
using namespace mc;

////////////////////////////////////////////////////////////////////////
int main()
////////////////////////////////////////////////////////////////////////
{  
  try{
    // Construction of multivariate polynomials
    t_SPoly::options.BASIS = t_SPoly::Options::MONOM;
#if defined( TEST_DOXYGEN )
    const unsigned NX = 5, NF = 3, NC = 3;
    t_SPoly X[NX], &FT = X[0], &F1 = X[1], &F2 = X[2], &x1 = X[3], &x2 = X[4];
#if defined( USE_DAG )
    FFGraph DAG;
    FFVar DAGX[NX];
    for( unsigned i=0; i<NX; i++ ) X[i].var( &DAGX[i].set( &DAG ) );
#else
    for( unsigned i=0; i<NX; i++ ) X[i].var( i );
#endif
    t_SPoly F[NF];
    F[0] = x1 * FT - F1;
    F[1] = x2 * FT - F2;
    //F[2] = x1 + x2 - 1.;
    F[2] = x1 * FT + x2 * FT - FT;

#elif defined( TEST_QUAD1 )
    const unsigned NX = 3, NF = 2, NC = 2;
    t_SPoly X[NX];
#if defined( USE_DAG )
    FFGraph DAG;
    FFVar DAGX[NX];
    for( unsigned i=0; i<NX; i++ ) X[i].var( &DAGX[i].set( &DAG ) );
#else
    for( unsigned i=0; i<NX; i++ ) X[i].var( i );
#endif
    t_SPoly F[NF];
    // X0^2 + X0·X1 - X0 = 0 && X0·X1 + X1^2 - X1 = 0 => X0^2 - X1^2 - X0 + X1 = 0
    F[0] = X[0] + X[1] - 1.;
    F[1] = sqr(X[0]) - sqr(X[1]) - X[2];

#elif defined( TEST_QUAD2 )
    const unsigned NX = 7, NF = 5, NC = 5;
    t_SPoly X[NX];
#if defined( USE_DAG )
    FFGraph DAG;
    FFVar DAGX[NX];
    for( unsigned i=0; i<NX; i++ ) X[i].var( &DAGX[i].set( &DAG ) );
#else
    for( unsigned i=0; i<NX; i++ ) X[i].var( i );
#endif
    t_SPoly F[NF];
    F[0] =   X[0]          + 2*X[2]          +   X[4] +   X[5] - 1;
    F[1] = 2*X[0] -   X[1]          +   X[3]          + 3*X[5] - 2;
    F[2] =            X[1]          + 6*X[3] + 2*X[4] - 3*X[5] + 1;
    F[3] = 2*X[0]                   +   X[3] + 3*X[4]          - 1;
    F[4] = X[6] - X[0]*X[0] - X[1]*X[1] - 3*(X[3]*X[3]) - X[4]*X[4] - 2*(X[5]*X[5])
         - 2*(X[0]*X[3]) - X[0]*X[4] - 2*(X[0]*X[5]) - X[1]*X[3] + X[1]*X[4]
         - 2*(X[1]*X[5]) + X[2]*X[3] - 4*(X[2]*X[4]) - 3*(X[2]*X[5]) - 6*(X[3]*X[4])
         - 9*(X[3]*X[5]) - X[4]*X[5];

#elif defined( TEST_TUNCPHD29 )
    const unsigned NX = 12, NF = 14, NC = 5;
    t_SPoly X[NX];
#if defined( USE_DAG )
    FFGraph DAG;
    FFVar DAGX[NX];
    for( unsigned i=0; i<NX; i++ ) X[i].var( &DAGX[i].set( &DAG ) );
#else
    for( unsigned i=0; i<NX; i++ ) X[i].var( i );
#endif
    t_SPoly F[NF];
    F[0] = 1.22*X[3] - X[0] - X[4];
    F[1] = 98e3*X[2]*X[10] - X[5];
    F[2] = (X[1]+X[4])*X[11] - X[7];
    F[3] = (X[3]*X[8] + 1e3*X[2])*X[10] - 1;
    F[4] = X[0]*X[11] - 1;
    F[5] = 5.04*X[0] + 0.035*X[1] +10*X[2] + 3.36*X[4] - 0.063*X[3]*X[6];
    F[6] = 35.82 - 0.222*X[9] - 0.9*X[8];
    F[7] = -133 + 3*X[6] - 0.99*X[9];
    F[8] = -(35.82-0.222*X[9]-0.9*X[8]) + 0.211111*X[8];
    F[9] = -(-133+3*X[6]-0.99*X[9]) + 0.020101*X[9];
    F[10] = 1.12*X[0] + 0.13167*X[0]*X[7] - 0.00667*X[0]*X[7]*X[7] - 0.99*X[3];
    F[11] = 57.425 + 1.098*X[7] - 0.038*X[7]*X[7] + 0.325*X[5] - 0.99*X[6];
    F[12] = -(1.12*X[0] + 0.13167*X[0]*X[7] - 0.00667*X[0]*X[7]*X[7] - 0.99*X[3]) + 0.020101*X[3];
    F[13] = -(57.425 + 1.098*X[7] - 0.038*X[7]*X[7] + .325*X[5] - 0.99*X[6]) + 0.020101*X[6];
#endif

    t_SRed SRF;
    //t_SRed::options.ORDER = 3;
    //t_SRed::options.NODIV = 0;
    //t_SRed::options.MIPOUTPUTFILE = "red.lp";
    //t_SRed::options.MIPNUMFOCUS = 3;
    //t_SRed::options.LPFEASTOL = 1e-7;
    //t_SRed::options.MIPTIMELIMIT  = 1800;
    //t_SRed::options.MIPCONCURRENT = 4;
    //t_SRed::options.MIPFOCUS = 0;
    //t_SRed::options.MIPHEURISTICS = 0.2;

    SRF.set_monomials( NF, F, &t_SPoly::mapmon );
    SRF.search_reductions( NC, F, &t_SPoly::mapmon );//, false );
    cout << SRF;
  }

  catch( t_SPoly::Exceptions &eObj ){
    cerr << "Error " << eObj.ierr()
         << " in polynomial expression:" << endl
	 << eObj.what() << endl
         << "Aborts." << endl;
    return eObj.ierr();
  }

  catch( t_SRed::Exceptions &eObj ){
    cerr << "Error " << eObj.ierr()
         << " in reduction constraint:" << endl
	 << eObj.what() << endl
         << "Aborts." << endl;
    return eObj.ierr();
  }

  return 0;
}
