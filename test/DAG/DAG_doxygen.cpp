////////////////////////////////////////////////////////////////////////

#include <fstream>
#include <iomanip>

#include "ffunc.hpp"
#include "interval.hpp"
#include "scmodel.hpp"

typedef mc::Interval I;
typedef mc::SCModel<I> SCM;
typedef mc::SCVar<I> SCV;

I IINF = 1e20 * I(-1,1);

///////////////////////////////////////////////////////////////////////////////

int test_build()
{
  std::cout << "\n==============================================\ntest_build:\n";

  // How Do I Construct the DAG of a Factorable Function?

  mc::FFGraph DAG;
  const unsigned int NX = 4;
  mc::FFVar X[NX];
  for( unsigned int i=0; i<NX; i++ ) X[i].set( &DAG );

  const unsigned int NF = 2;
  mc::FFVar F[NF]
    = { X[2]*X[3]-X[0],
        X[0]*pow(exp(X[2]*X[3])+3.,4)+X[1] };
  std::cout << DAG;

  DAG.output( DAG.subgraph( NF, F ), " F" );
  DAG.output( DAG.subgraph( 1, &F[0] ), " F0" );

  std::ofstream o_F( "F.dot", std::ios_base::out );
  DAG.dot_script( NF, F, o_F );
  o_F.close();

  std::ofstream o_F0( "F0.dot", std::ios_base::out );
  DAG.dot_script( 1, F, o_F0 );
  o_F0.close();

  // How Do I Obtain the DAG of a Factorable Function's Derivatives?

  const mc::FFVar* dFdX_FAD = DAG.FAD( NF, F, NX, X );
  std::cout << DAG;

  std::ofstream o_dFdX_FAD( "dFdX_FAD.dot", std::ios_base::out );
  DAG.dot_script( NX*NF, dFdX_FAD, o_dFdX_FAD );
  o_dFdX_FAD.close();
  std::ofstream o_dF1dX3_FAD( "dF1dX3_FAD.dot", std::ios_base::out );
  DAG.dot_script( 1, &dFdX_FAD[NX+3], o_dF1dX3_FAD );
  o_dF1dX3_FAD.close();
  delete[] dFdX_FAD;

  const mc::FFVar* dFdX_BAD = DAG.BAD( NF, F, NX, X );
  std::ofstream o_dFdX_BAD( "dFdX_BAD.dot", std::ios_base::out );
  DAG.dot_script( NX*NF, dFdX_BAD, o_dFdX_BAD );
  o_dFdX_BAD.close();
  std::ofstream o_dF1dX3_BAD( "dF1dX3_BAD.dot", std::ios_base::out );
  DAG.dot_script( 1, &dFdX_BAD[NX+3], o_dF1dX3_BAD );
  o_dF1dX3_BAD.close();
  delete[] dFdX_BAD;

  return 0;
}

///////////////////////////////////////////////////////////////////////////////

int test_eval()
{
  std::cout << "\n==============================================\ntest_eval:\n";

  // DAG environment
  mc::FFGraph DAG;
  // Independent variables and derivative direction
  const unsigned int NX = 4;
  mc::FFVar X[NX], D[NX] = { 0., 1., 1., 0. };
  for( unsigned int i=0; i<NX; i++ ) X[i].set( &DAG );
  // Dependent variables
  const unsigned int NF = 2;
  mc::FFVar F[NF]
    = { X[2]*X[3]-X[0],
        X[0]*pow(exp(X[2]*X[3])+3.,4)+X[1] };
  // DAG of second-order derivatives
  const mc::FFVar* dFdXdir = DAG.DFAD( NF, F, NX, X, D );

  DAG.output( DAG.subgraph( 1, &dFdXdir[0] ), " dF(0)dX·D" );
  DAG.output( DAG.subgraph( 1, &dFdXdir[1] ), " dF(1)dX·D" );

  // Evaluation in interval arithmetic
  try{
    I IX[NX] = { I(0,0.5), I(1,2), I(-1,-0.8), I(0.5,1) }, IdFdXdir[NF];
    std::vector<I> IWK;
    DAG.eval( IWK, NF, dFdXdir, IdFdXdir, NX, X, IX );
    // Display results
    for( unsigned i=0; i<NF; i++ )
      std::cout << "  dF("<< i << ")dX·D = " << IdFdXdir[i] << std::endl;
  }
  catch(...){
    std::cout << "\nInterval forward propagation failed\n";
  }

  // Evaluation in 3rd-order Chebyshev model arithmetic
  try{
    const unsigned ORD = 3;
    SCM CMenv( ORD );
    SCV CMX[NX], CMdFdXdir[NF];
    I IX[NX] = { I(0,0.5), I(1,2), I(-1,-0.8), I(0.5,1) };
    for( unsigned i=0; i<NX; i++ ) CMX[i].set( &CMenv, i, IX[i] );
    std::vector<SCV> SCVWK;
    DAG.eval( SCVWK, NF, dFdXdir, CMdFdXdir, NX, X, CMX );
    // Display results
    for( unsigned i=0; i<NF; i++ )
      std::cout << "  dF("<< i << ")dX·D = " << CMdFdXdir[i] << std::endl;
  }
  catch(...){
    std::cout << "\nSparse Chebyshev model forward propagation failed\n";
  }

  // Forward/backward evaluation in interval arithmetic
  try{
    I IX[NX] = { I(0,0.5), I(1,2), I(-1,-0.8), I(0.5,1) },
      IdFdXdir[NF] = { I(0.,1.), I(0.,5.) };
    std::vector<I> IWK;
    int flag = DAG.reval( IWK, NF, dFdXdir, IdFdXdir, NX, X, IX, IINF );
    std::cout << "\nDAG interval evaluation w/ forward/backward passes:\n";
    // Display results
    std::cout << "FLAG = " << flag << std::endl;
    for( unsigned i=0; i<NX; i++ )
      std::cout << "  X(" << i << ") = " << IX[i] << std::endl;
    for( unsigned i=0; i<NF; i++ )
      std::cout << "  dF("<< i << ")dX·D = " << IdFdXdir[i] << std::endl;
  }
  catch(...){
    std::cout << "\nInterval forward/backward propagation failed\n";
  }

  delete[] dFdXdir;
  return 0;
}

///////////////////////////////////////////////////////////////////////////////

int main()
{
  try{
    test_build();
    test_eval();
  }
  catch( mc::FFBase::Exceptions &eObj ){
    std::cerr << "Error " << eObj.ierr()
              << " in factorable function manipulation:" << std::endl
              << eObj.what() << std::endl
              << "Aborts." << std::endl;
    return eObj.ierr();
  }
}

