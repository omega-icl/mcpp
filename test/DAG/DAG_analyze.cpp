#undef MC__HSL_USE

#include <fstream>
#include <sstream>
#include <iomanip>

#include "mctime.hpp"
#include "ffunc.hpp"
#include "rltred.hpp"

///////////////////////////////////////////////////////////////////////////////

int test_dependency()
{
  std::cout << "\n==============================================\ntest_dependency:\n";

  const int NX = 4;
  mc::FFDep X[NX];
  for( int i=0; i<NX; i++ ) X[i].indep(i);

  const int NF = 2;
  mc::FFDep F[NF] = { X[2]*X[3]+X[0]/X[2],
                      X[0]*pow(exp(X[2]-X[3]),2)+X[1] };

  std::map<int,int> F0_dep = F[0].dep();
  std::map<int,int> F1_dep = F[1].dep();

  std::cout << "Variable dependence of F[0]: " << F[0] << std::endl;
  std::cout << "Variable dependence of F[1]: " << F[1] << std::endl;

  return 0;
}

///////////////////////////////////////////////////////////////////////////////

int test_rltred1()
{
  std::cout << "\n==============================================\ntest_rltred1:\n";

  mc::FFGraph DAG;
  const unsigned NX = 5, NF = 3;
  mc::FFVar X[NX], &FT = X[0], &F1 = X[1], &F2 = X[2], &x1 = X[3], &x2 = X[4];
  for( unsigned i(0); i<NX; i++ )  X[i].set( &DAG );
  mc::FFVar F[NF];
  F[0] = x1 * FT - F1;
  F[1] = x2 * FT - F2;
  F[2] = x1 + x2 -1.;
  std::cout << DAG;

  mc::RLTRED RRLT( &DAG );
  RRLT.options.DISPLAY = 1;
  RRLT.options.LEVEL   = mc::RLTRED::Options::PRIMSEQ;
  RRLT.options.NODIV   = false;

  RRLT.search( NF, F );

  auto FRED = RRLT.constraints();
  for( auto it=FRED.begin(); it!=FRED.end(); ++it ){
    std::ostringstream ostr; ostr << " OF REDUCTION CONSTRAINT " << **it;
    DAG.output( DAG.subgraph( 1, *it ), ostr.str() );
  }

  F[2] = F1 + F2 - FT;
  RRLT.search( NF, F );

  FRED = RRLT.constraints();
  for( auto it=FRED.begin(); it!=FRED.end(); ++it ){
    std::ostringstream ostr; ostr << " OF REDUCTION CONSTRAINT " << **it;
    DAG.output( DAG.subgraph( 1, *it ), ostr.str() );
  }

  return 0;
}

///////////////////////////////////////////////////////////////////////////////

int test_rltred2()
{
  std::cout << "\n==============================================\ntest_rltred2:\n";

  mc::FFGraph DAG;
  const unsigned NX = 2, NF = 2;
  mc::FFVar X[NX];
  for( unsigned i(0); i<NX; i++ ) X[i].set( &DAG );
  mc::FFVar F[NF];
  F[0] = X[0] + X[1] - 1.;
  F[1] = sqr(X[0]) - sqr(X[1]);
  std::cout << DAG;

  mc::RLTRED RRLT( &DAG );
  RRLT.options.DISPLAY = 1;
  RRLT.options.LEVEL   = mc::RLTRED::Options::FULLSEQ;
  RRLT.options.NODIV   = false;
  RRLT.search( NF, F );

  auto FRED = RRLT.constraints();
  for( auto it=FRED.begin(); it!=FRED.end(); ++it ){
    std::ostringstream ostr; ostr << " OF REDUCTION CONSTRAINT " << **it;
    DAG.output( DAG.subgraph( 1, *it ), ostr.str() );
  }

  return 0;
}

///////////////////////////////////////////////////////////////////////////////

int test_rltred3()
{
  std::cout << "\n==============================================\ntest_rltred3:\n";

  mc::FFGraph DAG;
  const unsigned NX = 7, NF = 5;
  std::vector<mc::FFVar> X(NX);
  for( unsigned i(0); i<NX; i++ ) X[i].set( &DAG );
  std::vector<mc::FFVar> F(NF);
  F[0] =   X[0]          + 2*X[2]          +   X[4] +   X[5] - 1;
  F[1] = 2*X[0] -   X[1]          +   X[3]          + 3*X[5] - 2;
  F[2] =            X[1]          + 6*X[3] + 2*X[4] - 3*X[5] + 1;
  F[3] = 2*X[0]                   +   X[3] + 3*X[4]          - 1;
  F[4] = X[6] - X[0]*X[0] - X[1]*X[1] - 3*(X[3]*X[3]) - X[4]*X[4] - 2*(X[5]*X[5])
       - 2*(X[0]*X[3]) - X[0]*X[4] - 2*(X[0]*X[5]) - X[1]*X[3] + X[1]*X[4]
       - 2*(X[1]*X[5]) + X[2]*X[3] - 4*(X[2]*X[4]) - 3*(X[2]*X[5]) -6*(X[3]*X[4])
       - 9*(X[3]*X[5]) - X[4]*X[5];
  std::cout << DAG;

  mc::RLTRED RRLT( &DAG );
  RRLT.options.DISPLAY = 1;
  RRLT.options.LEVEL   = mc::RLTRED::Options::PRIMSEQ;
  RRLT.options.NODIV   = false;
  RRLT.search( NF, F.data() );
  auto FRED = RRLT.constraints();
  for( auto it=FRED.begin(); it!=FRED.end(); ++it ){
    std::ostringstream ostr; ostr << " OF REDUCTION CONSTRAINT " << **it;
    DAG.output( DAG.subgraph( 1, *it ), ostr.str() );
    break;
  }

#ifdef MC__HSL_USE
  for( auto it=FRED.begin(); it!=FRED.end(); ++it ) F.push_back(**it);
  const unsigned NFRED = FRED.size();
  std::vector<int> IP(NF+NFRED), IQ(NX), IPROF(NX), IFLAG(3);
  DAG.MC33( NF, F.data(), NX, X.data(), IP.data(), IQ.data(),
            IPROF.data(), IFLAG.data(), true );
  DAG.MC33( NF+NFRED, F.data(), NX, X.data(), IP.data(), IQ.data(),
            IPROF.data(), IFLAG.data(), true );
#endif
  return 0;
}

///////////////////////////////////////////////////////////////////////////////

int main()
{
  try{
    test_dependency();
    test_rltred1();
    test_rltred2();
    test_rltred3();
  }
  catch( mc::FFGraph::Exceptions &eObj ){
    std::cerr << "Error " << eObj.ierr()
              << " in factorable function manipulation:" << std::endl
              << eObj.what() << std::endl
              << "Aborts." << std::endl;
    return eObj.ierr();
  }
}

