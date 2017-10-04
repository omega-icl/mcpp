////////////////////////////////////////////////////////////////////////
#undef USE_PROFIL	// <-- specify to use PROFIL for interval arithmetic
#undef USE_FILIB	// <-- specify to use FILIB++ for interval arithmetic
////////////////////////////////////////////////////////////////////////

#include <fstream>
#include <iomanip>

#include "mctime.hpp"
#include "ffunc.hpp"

///////////////////////////////////////////////////////////////////////////////

int test_fadiff1()
{
  std::cout << "\n==============================================\ntest_fadiff1:\n";

  // Create DAG
  const unsigned NX = 4, NF = 2;
  mc::FFGraph DAG;
  mc::FFVar X[NX];
  for( unsigned int i=0; i<NX; i++ ) X[i].set( &DAG );
  mc::FFVar F[NF] = { X[2]*X[3]-2./(X[1]+X[2]),
                      X[0]*pow(exp(X[2]*X[1])+3.,4)+tanh(X[3]) };
  std::cout << DAG;

  std::ofstream o_F( "fadiff1_F.dot", std::ios_base::out );
  DAG.dot_script( NF, F, o_F );
  o_F.close();

  // Forward AD
  const mc::FFVar* dFdX_FAD = DAG.FAD( NF, F, NX, X );
  std::cout << DAG;

  DAG.output( DAG.subgraph( NX*NF, dFdX_FAD ), " dFdX_FAD" );
  std::ofstream o_dFdX_FAD( "fadiff1_dFdX_FAD.dot", std::ios_base::out );
  DAG.dot_script( NX*NF, dFdX_FAD, o_dFdX_FAD );
  o_dFdX_FAD.close();
      
  DAG.output( DAG.subgraph( 1, &dFdX_FAD[NX+3] ), " dF1dX3_FAD" );
  std::ofstream o_dF1dX3_FAD( "fadiff1_dF1dX3_FAD.dot", std::ios_base::out );
  DAG.dot_script( 1, &dFdX_FAD[NX+3], o_dF1dX3_FAD );
  o_dF1dX3_FAD.close();

  delete[] dFdX_FAD;

  // Backward AD
  const mc::FFVar* dFdX_BAD = DAG.BAD( NF, F, NX, X );
  std::cout << DAG;

  DAG.output( DAG.subgraph( NX*NF, dFdX_BAD ), " dFdX_BAD" );
  std::ofstream o_dFdX_BAD( "fadiff1_dFdX_BAD.dot", std::ios_base::out );
  DAG.dot_script( NX*NF, dFdX_BAD, o_dFdX_BAD );
  o_dFdX_BAD.close();
      
  DAG.output( DAG.subgraph( 1, &dFdX_BAD[NX+3] ), " dF1dX3_BAD" );
  std::ofstream o_dF1dX3_BAD( "fadiff1_dF1dX3_BAD.dot", std::ios_base::out );
  DAG.dot_script( 1, &dFdX_BAD[NX+3], o_dF1dX3_BAD );
  o_dF1dX3_BAD.close();

  delete[] dFdX_BAD;

  return 0;
}

///////////////////////////////////////////////////////////////////////////////

int test_fadiff2()
{
  std::cout << "\n==============================================\ntest_fadiff2:\n";

  // Create DAG
  const unsigned NX = 2, NF = 2;
  mc::FFGraph DAG;
  mc::FFVar X[NX];
  for( unsigned int i=0; i<NX; i++ ) X[i].set( &DAG );
  mc::FFVar F[NF] = { sqrt(X[0])*exp(X[1])*X[0],
                      pow(X[1],3)*sqrt(X[0])     };
  std::cout << DAG;
  std::ofstream o_F( "fadiff2_F.dot", std::ios_base::out );
  DAG.dot_script( NF, F, o_F );
  o_F.close();

  // Forward AD
  const mc::FFVar* dFdX = DAG.FAD( NF, F, NX, X, true );
  DAG.output( DAG.subgraph( NX*NF, dFdX ), " dFdX (non-recursive)" );
  std::ofstream o_dFdX( "fadiff2_dFdX.dot", std::ios_base::out );
  DAG.dot_script( NF*NX, dFdX, o_dFdX );
  o_dFdX.close();
  delete[] dFdX;

  // Forward AD recursive
  const mc::FFVar* dFdX_recur = DAG.FAD( NF, F, 1, &X[0], 1, &X[1], true );
  DAG.output( DAG.subgraph( NX*NF, dFdX_recur ), " dFdX (recursive)" );
  std::ofstream o_dFdX_recur( "fadiff2_dFdX_recur.dot", std::ios_base::out );
  DAG.dot_script( NF*NX, dFdX_recur, o_dFdX_recur );
  o_dFdX_recur.close();
  delete[] dFdX_recur;

  return 0;
}

///////////////////////////////////////////////////////////////////////////////

int test_fadiff_directional()
{
  std::cout << "\n==============================================\ntest_fadiff_directional:\n";

  // Create DAG
  mc::FFGraph DAG;
  mc::FFVar X[2], F[2], D[2];
  X[0].set( &DAG );
  X[1].set( &DAG );
  F[0] = sqrt(X[0])*X[1];
  F[1] = X[0]*X[1];
  D[0] = 1.;
  D[1] = X[0];
  std::cout << DAG;

  // Directional forward AD
  const mc::FFVar* dFdXxD = DAG.DFAD( 2, F, 2, X, D );
  DAG.output( DAG.subgraph( 1, dFdXxD ), " dF0dX·D (non-recursive)" );
  DAG.output( DAG.subgraph( 1, dFdXxD+1), " dF1dX·D (non-recursive)" );
  std::ofstream o_dFdXxD( "fadiff_directional.dot", std::ios_base::out );
  DAG.dot_script( 1, dFdXxD, o_dFdXxD );
  o_dFdXxD.close();
  delete[] dFdXxD;

  // Directional forward AD recursive
  const mc::FFVar* dFdXxD_recur = DAG.DFAD( 2, F, 1, &X[0], &D[0], 1, &X[1], &D[1] );
  DAG.output( DAG.subgraph( 1, dFdXxD_recur ), " dF0dX·D (recursive)" );
  DAG.output( DAG.subgraph( 1, dFdXxD_recur+1), " dF1dX·D (recursive)" );
  std::ofstream o_dFdXxD_recur( "fadiff_directional_recur.dot", std::ios_base::out );
  DAG.dot_script( 1, dFdXxD_recur, o_dFdXxD_recur );
  o_dFdXxD_recur.close();
  delete[] dFdXxD_recur;

  return 0;
}

///////////////////////////////////////////////////////////////////////////////

int test_gradient_sparse()
{
  std::cout << "\n==============================================\ntest_gradient_sparse:\n";

  // Create DAG
  mc::FFGraph DAG;
  const unsigned NX = 2, NF = 3;
  mc::FFVar X[NX], D[NX] = { 1., 2. };
  for( unsigned i(0); i<NX; i++ )  X[i].set( &DAG );
  mc::FFVar F[NF] = { 0.5*X[0], X[0]*X[1], sqrt(X[0])*exp(X[1])*X[0]+1. };
  std::cout << DAG;

  // Sparse forward AD
  //auto dFdX_FAD = DAG.SDFAD( NF, F, NX-1, X, D, 1, X+NX-1, D+NX-1 );
  //auto dFdX_FAD = DAG.SFAD( NF, F, NX, X );
  auto dFdX_FAD = DAG.SFAD( NF, F, NX-1, X, 1, X+NX-1 );
  std::cout << "\nNon-zero Jacobian elements (FAD, non-recursive): "
            << std::get<0>(dFdX_FAD) << std::endl;
  for( unsigned ie=0; ie<std::get<0>(dFdX_FAD); ie++ )
    std::cout << "(" << std::get<1>(dFdX_FAD)[ie] << "," << std::get<2>(dFdX_FAD)[ie] << ") "
                     << std::get<3>(dFdX_FAD)[ie] << std::endl;
  DAG.output( DAG.subgraph( std::get<0>(dFdX_FAD), std::get<3>(dFdX_FAD) ) );
  std::ofstream o_dFdX_FAD( "gradient_sparse_FAD.dot", std::ios_base::out );
  DAG.dot_script( std::get<0>(dFdX_FAD), std::get<3>(dFdX_FAD), o_dFdX_FAD );
  o_dFdX_FAD.close();
  delete[] std::get<1>(dFdX_FAD);
  delete[] std::get<2>(dFdX_FAD);
  delete[] std::get<3>(dFdX_FAD);

  // Sparse backward AD
  //auto dFdX_BAD = DAG.SBAD( NF, F, NX, X );
  auto dFdX_BAD = DAG.SBAD( NF, F, NX-1, X, 1, X+NX-1 );
  std::cout << "\nNon-zero Jacobian elements (BAD, non-recursive): "
            << std::get<0>(dFdX_BAD) << std::endl;
  for( unsigned ie=0; ie<std::get<0>(dFdX_BAD); ie++ )
    std::cout << "(" << std::get<1>(dFdX_BAD)[ie] << "," << std::get<2>(dFdX_BAD)[ie] << ") "
                     << std::get<3>(dFdX_BAD)[ie] << std::endl;
  DAG.output( DAG.subgraph( std::get<0>(dFdX_BAD), std::get<3>(dFdX_BAD) ) );
  std::ofstream o_dFdX_BAD( "gradient_sparse_BAD.dot", std::ios_base::out );
  DAG.dot_script( std::get<0>(dFdX_BAD), std::get<3>(dFdX_BAD), o_dFdX_BAD );
  o_dFdX_BAD.close();

  delete[] std::get<1>(dFdX_BAD);
  delete[] std::get<2>(dFdX_BAD);
  delete[] std::get<3>(dFdX_BAD);

  return 0;
}

///////////////////////////////////////////////////////////////////////////////

int test_hessian_sparse()
{
  std::cout << "\n==============================================\ntest_hessian_sparse:\n";

  mc::FFGraph DAG;
  const unsigned NX = 5, NF = 1;
  mc::FFVar X[NX];
  for( unsigned i(0); i<NX; i++ )  X[i].set( &DAG );
  mc::FFVar F[NF] = { sqrt(X[0])*X[1]*(X[3]+2.) };
  std::cout << DAG;

  auto dFdX = DAG.SBAD( NF, F, NX, X );
  //std::cout << DAG;
  std::cout << "\nNumber of non-zero elements in Jacobian matrix: "
            << std::get<0>(dFdX) << std::endl;
  for( unsigned ie=0; ie<std::get<0>(dFdX); ie++ )
    std::cout << "(" << std::get<1>(dFdX)[ie] << "," << std::get<2>(dFdX)[ie] << ") "
              << std::get<3>(dFdX)[ie] << std::endl;

  DAG.output( DAG.subgraph( std::get<0>(dFdX), std::get<3>(dFdX) ), " dFdX" );
  std::ofstream o_dFdX( "hessian_sparse_BAD.dot", std::ios_base::out );
  DAG.dot_script( std::get<0>(dFdX), std::get<3>(dFdX), o_dFdX );
  o_dFdX.close();

  auto d2FdX2 = DAG.SFAD( std::get<0>(dFdX), std::get<3>(dFdX), NX, X, true );
  //std::cout << DAG;
  std::cout << "\nNumber of non-zero elements in Hessian matrix: "
            << std::get<0>(d2FdX2) << std::endl;
  for( unsigned ie=0; ie<std::get<0>(d2FdX2); ie++ )
    std::cout << "(" << std::get<1>(d2FdX2)[ie] << "," << std::get<2>(d2FdX2)[ie] << ") "
              << std::get<3>(d2FdX2)[ie] << std::endl;

  DAG.output( DAG.subgraph( std::get<0>(d2FdX2), std::get<3>(d2FdX2) ), " d2FdX2" );
  std::ofstream o_d2FdX2( "hessian_sparse_FAD+BAD.dot", std::ios_base::out );
  DAG.dot_script( std::get<0>(d2FdX2), std::get<3>(d2FdX2), o_d2FdX2 );
  o_d2FdX2.close();

  delete[] std::get<1>(dFdX);
  delete[] std::get<2>(dFdX);
  delete[] std::get<3>(dFdX);

  delete[] std::get<1>(d2FdX2);
  delete[] std::get<2>(d2FdX2);
  delete[] std::get<3>(d2FdX2);
  return 0;
}

///////////////////////////////////////////////////////////////////////////////

int test_tadiff()
{
  mc::FFGraph DAG;
  mc::FFVar T( &DAG );
  mc::FFVar X( &DAG ); 
  mc::FFVar F = sqr(X)+exp(T);
  std::cout << DAG;

  const unsigned NTE = 5;
  const mc::FFVar* F_TAD = DAG.TAD( NTE, 1, &F, 1, &X, &T );
  std::cout << DAG;

  DAG.output( DAG.subgraph( NTE+1, F_TAD ), " F_TAD" );
  std::ofstream o_TAD( "tadiff.dot", std::ios_base::out );
  DAG.dot_script( NTE+1, F_TAD, o_TAD );
  o_TAD.close();

  delete[] F_TAD;
  return 0;
}

///////////////////////////////////////////////////////////////////////////////

int main()
{
  try{
    test_fadiff1();
    test_fadiff2();
    test_fadiff_directional();
    test_gradient_sparse();
    test_hessian_sparse();
    test_tadiff();
  }
  catch( mc::FFGraph::Exceptions &eObj ){
    std::cerr << "Error " << eObj.ierr()
              << " in factorable function manipulation:" << std::endl
              << eObj.what() << std::endl
              << "Aborts." << std::endl;
    return eObj.ierr();
  }
}

