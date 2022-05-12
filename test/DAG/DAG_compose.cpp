#include <fstream>
#include <iomanip>

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

int test_insert()
{
  std::cout << "\n==============================================\ntest_insert:\n";

  mc::FFGraph DAG1;
  mc::FFVar X1[3], F1[3];
  for( unsigned i=0; i<3; i++ )
    X1[i].set( &DAG1 );
  F1[0] = exp(X1[1]);
  F1[1] = sqr(X1[0])-exp(X1[1])*X1[0];
  F1[2] = sqr(X1[2])-exp(X1[1])*X1[2];
  std::cout << DAG1;

  mc::FFGraph DAG2;
  mc::FFVar F2[2];
  DAG2.insert( &DAG1, 1, F1, F2 );
  std::cout << DAG2;
  DAG2.insert( &DAG1, 1, F1+1, F2+1 );
  std::cout << DAG2;

  mc::FFGraph DAG3;
  mc::FFVar F3[3];
  DAG3.insert( &DAG1, 3, F1, F3 );
  std::cout << DAG3;

  return 0;
}

///////////////////////////////////////////////////////////////////////////////

int main()
{
  try{
    test_compose1();
    test_compose2();
    test_compose3();
    test_insert();
  }
  catch( mc::FFGraph::Exceptions &eObj ){
    std::cerr << "Error " << eObj.ierr()
              << " in factorable function manipulation:" << std::endl
              << eObj.what() << std::endl
              << "Aborts." << std::endl;
    return eObj.ierr();
  }
}

