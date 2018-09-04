////////////////////////////////////////////////////////////////////////
#undef USE_PROFIL	// <-- specify to use PROFIL for interval arithmetic
#undef USE_FILIB	// <-- specify to use FILIB++ for interval arithmetic
////////////////////////////////////////////////////////////////////////

#include <fstream>
#include <iomanip>

#include "mctime.hpp"
#include "ffunc.hpp"

///////////////////////////////////////////////////////////////////////////////

int test_compose1()
{
  std::cout << "\n==============================================\ntest_compose1:\n";

  mc::FFGraph DAG;
  mc::FFVar X, Y, F, G;
  X.set( &DAG );
  Y.set( &DAG );
  F = exp(X);
  G = sqr(Y)+F;
  std::cout << DAG;

  std::ofstream o_comp0( "compose1_0.dot", std::ios_base::out );
  DAG.dot_script( 1, &G, o_comp0 );
  o_comp0.close();

  const mc::FFVar* GoF = DAG.compose( 1, &G, 1, &Y, &F );
  std::cout << DAG;

  std::ofstream o_comp1( "compose1_1.dot", std::ios_base::out );
  DAG.dot_script( 1, GoF, o_comp1 );
  o_comp1.close();

  delete[] GoF;
  return 0;
}

///////////////////////////////////////////////////////////////////////////////

int test_compose2()
{
  std::cout << "\n==============================================\ntest_compose2:\n";

  mc::FFGraph DAG;
  mc::FFVar X, Y, Z, F, G;
  X.set( &DAG );
  Y.set( &DAG );
  Z.set( &DAG );
  F = exp(X);
  G = sqr(Y)+Z;
  std::cout << DAG;

  std::ofstream o_comp0( "compose2_0.dot", std::ios_base::out );
  DAG.dot_script( 1, &G, o_comp0 );
  o_comp0.close();

  const mc::FFVar* GoF = DAG.compose( 1, &G, 1, &Y, &F, 1, &Z, &F );
  std::cout << DAG;

  std::ofstream o_comp1( "compose2_1.dot", std::ios_base::out );
  DAG.dot_script( 1, GoF, o_comp1 );
  o_comp1.close();

  delete[] GoF;
  return 0;
}

///////////////////////////////////////////////////////////////////////////////

int test_compose3()
{
  std::cout << "\n==============================================\ntest_compose3:\n";

  mc::FFGraph DAG;
  mc::FFVar X, Y, Z, F, G;
  X.set( &DAG );
  Y.set( &DAG );
  Z.set( &DAG );
  F = exp(X);
  G = sqr(Y)+F;
  std::cout << DAG;

  std::ofstream o_comp0( "compose3_0.dot", std::ios_base::out );
  DAG.dot_script( 1, &G, o_comp0 );
  o_comp0.close();

  const mc::FFVar* GoF = DAG.compose( 1, &G, 1, &F, &Z );
  std::cout << DAG;

  std::ofstream o_comp1( "compose3_1.dot", std::ios_base::out );
  DAG.dot_script( 1, GoF, o_comp1 );
  o_comp1.close();

  delete[] GoF;
  return 0;
}

///////////////////////////////////////////////////////////////////////////////

int main()
{
  try{
//    test_compose1();
//    test_compose2();
    test_compose3();
  }
  catch( mc::FFGraph::Exceptions &eObj ){
    std::cerr << "Error " << eObj.ierr()
              << " in factorable function manipulation:" << std::endl
              << eObj.what() << std::endl
              << "Aborts." << std::endl;
    return eObj.ierr();
  }
}

