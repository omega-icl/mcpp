////////////////////////////////////////////////////////////////////////

#include <fstream>
#include <iomanip>

#include "ffunc.hpp"
#include "interval.hpp"
#include "scmodel.hpp"

typedef mc::Interval I;
I IINF = 1e20 * I(-1,1);

typedef mc::SCModel<I,mc::FFVar*,mc::lt_FFVar> SCM;
typedef mc::SCVar<I,mc::FFVar*,mc::lt_FFVar> SCV;

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
        X[0]*pow(exp(X[2]*X[3])+3.1,4)+X[1] };
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
    DAG.eval( NF, dFdXdir, IdFdXdir, NX, X, IX );
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
    for( unsigned i=0; i<NX; i++ ) CMX[i].set( &CMenv, &X[i], IX[i] );
    DAG.eval( NF, dFdXdir, CMdFdXdir, NX, X, CMX );
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
      IdFdXdir[NF] = { I(0.6,0.9), I(2.,5.) };
    int flag = DAG.reval( NF, dFdXdir, IdFdXdir, NX, X, IX, IINF );
    std::cout << "\nDAG interval evaluation w/ " << flag << " forward/backward passes:\n";
    // Display results
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

namespace mc
{
class FFnorm2
: public FFOp
{
public:
  // Constructors
  FFnorm2
    ()
    : FFOp( (int)EXTERN )
    {}

  // Functor
  FFVar& operator()
    ( unsigned const nVar, FFVar const* pVar )
    const
    {
      auto dep = FFDep();
      for( unsigned i=0; i<nVar; ++i ) dep += pVar[i].dep();
      dep.update( FFDep::TYPE::N );
      return insert_external_operation( *this, dep, nVar, pVar );
    }

  // Evaluation overloads
  template< typename T > void eval
    ( T& vRes, unsigned const nVar, T const* vVar )
    {
      switch( nVar ){
        case 0: vRes = T( 0. ); break;
        case 1: vRes = vVar[0]; break;
        default: vRes = Op<T>::sqr( vVar[0] );
                 for( unsigned i=1; i<nVar; ++i ) vRes += Op<T>::sqr( vVar[i] );
                 vRes = Op<T>::sqrt( vRes ); break;
      }
    }
  void eval
    ( FFVar& vRes, unsigned const nVar, FFVar const* pVar )
    const
    {
      vRes = operator()( nVar, pVar );
    }

  // Properties
  std::string name
    ()
    const
    { return "NORM2"; }
  //! @brief Return whether or not operation is commutative
  bool commutative
    ()
    const
    { return true; }
};
}

int test_extern()
{
  std::cout << "\n==============================================\ntest_extern:\n";

  // DAG environment
  mc::FFGraph<mc::FFnorm2> DAG;
  const unsigned int NX = 3;
  mc::FFVar X[NX];
  for( unsigned int i=0; i<NX; i++ ) X[i].set( &DAG );
  mc::FFnorm2 norm2;
  mc::FFVar F = norm2( NX, X );
  std::cout << DAG;

  // Evaluation in interval arithmetic
  try{
    I IX[NX] = { I(0,0.5), I(1,2), I(-1,-0.8) }, IF;
    DAG.eval( 1, &F, &IF, NX, X, IX );
    // Display results
    std::cout << "  " << F << " = " << IF << std::endl;
  }
  catch(...){
    std::cout << "\nInterval evaluation failed\n";
  }

  // Differentiation
  const mc::FFVar* dFFAD = DAG.FAD( 1, &F, NX, X );
  std::cout << DAG;
  DAG.output( DAG.subgraph( NX, dFFAD ), " dnorm2_FAD" );

  std::ofstream o_FFAD( "dnorm2_FAD.dot", std::ios_base::out );
  DAG.dot_script( NX, dFFAD, o_FFAD );
  o_FFAD.close();
  delete[] dFFAD;

  const mc::FFVar* dFBAD = DAG.BAD( 1, &F, NX, X );
  std::cout << DAG;
  DAG.output( DAG.subgraph( NX, dFBAD ), " dnorm2_BAD" );

  std::ofstream o_FBAD( "dnorm2_BAD.dot", std::ios_base::out );
  DAG.dot_script( NX, dFBAD, o_FBAD );
  o_FBAD.close();
  delete[] dFBAD;

  return 0;
}

///////////////////////////////////////////////////////////////////////////////

int main()
{
  try{
    test_build();
    test_eval();
    test_extern();
  }
  catch( mc::FFBase::Exceptions &eObj ){
    std::cerr << "Error " << eObj.ierr()
              << " in factorable function manipulation:" << std::endl
              << eObj.what() << std::endl
              << "Aborts." << std::endl;
    return eObj.ierr();
  }
}

