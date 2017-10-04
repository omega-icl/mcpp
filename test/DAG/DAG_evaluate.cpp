////////////////////////////////////////////////////////////////////////
#define USE_PROFIL	// <-- specify to use PROFIL for interval arithmetic
#undef  USE_FILIB	// <-- specify to use FILIB++ for interval arithmetic
#undef  MC__FFUNC_CPU_EVAL
////////////////////////////////////////////////////////////////////////

#include <fstream>
#include <iomanip>

#include "mctime.hpp"
#include "ffunc.hpp"

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

#include "cmodel.hpp"
typedef mc::CModel<I> CM;
typedef mc::CVar<I> CV;

#include "scmodel.hpp"
typedef mc::SCModel<I> SCM;
typedef mc::SCVar<I> SCV;

///////////////////////////////////////////////////////////////////////////////

int test_eval1()
{
  std::cout << "\n==============================================\ntest_eval1:\n";

  // Create DAG
  const unsigned NX = 4, NF = 2;
  mc::FFGraph DAG;
  mc::FFVar X[NX];
  for( unsigned int i=0; i<NX; i++ ) X[i].set( &DAG );
  mc::FFVar F[NF] = { X[2]*X[3]-2./(X[1]+X[2]),
                      X[0]/pow(exp(X[2]*X[1])+3.,3)+tanh(X[3]) };
  std::cout << DAG;

  std::ofstream o_F( "eval1_F.dot", std::ios_base::out );
  DAG.dot_script( NF, F, o_F );
  o_F.close();

  double cputime;
  const unsigned NREP=1000000;

  // Evaluate with doubles, no parameter pack
  std::list<const mc::FFOp*> F_op  = DAG.subgraph( NF, F );
  double dX[NX] = { -1., -1., 2., 3. };
  for( unsigned i=0; i<NX; i++ ) std::cout << "X[" << i << "] = " << dX[i] << std::endl;

  double dF[NF];
  cputime = -mc::cpuclock();
  for( unsigned i=0; i<NREP; i++ ){
    dF[0] = dX[2]*dX[3]-2./(dX[1]+dX[2]);
    dF[1] = dX[0]/pow(exp(dX[2]*dX[1])+3.,3)+tanh(dX[3]);
  }
  cputime += mc::cpuclock();
  std::cout << "\nCompiled evaluation: " << (cputime/=NREP) << " CPU-sec\n";
  for( unsigned i=0; i<NF; i++ ) std::cout << "F[" << i << "] = " << dF[i] << std::endl;

  cputime = -mc::cpuclock();
  for( unsigned i=0; i<NREP; i++ )
    DAG.eval( F_op, NF, F, dF, NX, X, dX );
  cputime += mc::cpuclock();
  std::cout << "\nDAG evaluation - no preallocation, no variadic template: " << (cputime/=NREP) << " CPU-sec\n";
  for( unsigned i=0; i<NF; i++ ) std::cout << "F[" << i << "] = " << dF[i] << std::endl;

  cputime = -mc::cpuclock();
  for( unsigned i=0; i<NREP; i++ )
    DAG.eval( F_op, NF, F, dF, 1, &X[0], &dX[0], 1, &X[1], &dX[1], 1, &X[2], &dX[2], 1, &X[3], &dX[3] );
  cputime += mc::cpuclock();
  std::cout << "\nDAG evaluation - no preallocation, with variadic template: " << (cputime/=NREP) << " CPU-sec\n";
  for( unsigned i=0; i<NF; i++ ) std::cout << "F[" << i << "] = " << dF[i] << std::endl;

  std::vector<double> WK;
  cputime = -mc::cpuclock();
  for( unsigned i=0; i<NREP; i++ )
    DAG.eval( F_op, WK, NF, F, dF, NX, X, dX );
  cputime += mc::cpuclock();
  std::cout << "\nDAG evaluation - with preallocation, no variadic template: " << (cputime/=NREP) << " CPU-sec\n";
  for( unsigned i=0; i<NF; i++ ) std::cout << "F[" << i << "] = " << dF[i] << std::endl;

  cputime = -mc::cpuclock();
  for( unsigned i=0; i<NREP; i++ )
    DAG.eval( F_op, WK, NF, F, dF, 1, &X[0], &dX[0], 1, &X[1], &dX[1], 1, &X[2], &dX[2], 1, &X[3], &dX[3] );
  cputime += mc::cpuclock();
  std::cout << "\nDAG evaluation - with preallocation, with variadic template: " << (cputime/=NREP) << " CPU-sec\n";
  for( unsigned i=0; i<NF; i++ ) std::cout << "F[" << i << "] = " << dF[i] << std::endl;

  return 0;
}

///////////////////////////////////////////////////////////////////////////////

int test_eval2()
{
  std::cout << "\n==============================================\ntest_fadiff2:\n";

  // Create DAG
  const unsigned int NX = 2, NF = 2;
  mc::FFGraph DAG;
  mc::FFVar X[NX];
  for( unsigned int i=0; i<NX; i++ ) X[i].set( &DAG );
  mc::FFVar F[NF] = { sqrt(X[0])*exp(X[1])*X[0],
                      pow(X[1],3)*sqrt(X[0]) };
  std::cout << DAG;

  std::list<const mc::FFOp*> F_op  = DAG.subgraph( NF, F );
  DAG.output( F_op );
  std::ofstream o_F( "eval2_F.dot", std::ios_base::out );
  DAG.dot_script( NF, F, o_F );
  o_F.close();

  double cputime;
  const unsigned NREP=10000;

  // Evaluate in interval arithmetic
  I IX[NX] = { I(1.,2.), I(2.,3.) }, IF[2];
  for( unsigned i=0; i<NX; i++ ) std::cout << "X[" << i << "] = " << IX[i] << std::endl;
  cputime = -mc::cpuclock();
  for( unsigned i=0; i<NREP; i++ )
    DAG.eval( F_op, NF, F, IF, NX, X, IX );
  cputime += mc::cpuclock();
  std::cout << "\nDAG interval evaluation - no preallocation, no variadic template: " << (cputime/=NREP) << " CPU-sec\n";
  for( unsigned i=0; i<NF; i++ ) std::cout << "F[" << i << "] = " << IF[i] << std::endl;

  for( unsigned NTE=1; NTE<=5; NTE++ ){
    SCM modSCM( NTE );
    //modSCM.options.BOUNDER_TYPE = SCM::Options::LSB;
    modSCM.options.MIXED_IA = false;
    SCV SCX[NX] = {  SCV( &modSCM, 0, IX[0] ), SCV( &modSCM, 1, IX[1] ) }, SCF[2];
    cputime = -mc::cpuclock();
    for( unsigned i=0; i<NREP; i++ )
      DAG.eval( F_op, NF, F, SCF, NX, X, SCX );
    cputime += mc::cpuclock();
    std::cout << "\nDAG " << NTE << "th-order sparse Chebyshev model evaluation - no preallocation, no variadic template: " << (cputime/=NREP) << " CPU-sec\n";
    for( unsigned i=0; i<NF; i++ ) std::cout << "F[" << i << "] = " << SCF[i].R() << std::endl;
  }

  for( unsigned NTE=1; NTE<=5; NTE++ ){
    CM modCM( NX, NTE );
    //modSCM.options.BOUNDER_TYPE = SCM::Options::LSB;
    modCM.options.MIXED_IA = false;
    CV CX[NX] = {  CV( &modCM, 0, IX[0] ), CV( &modCM, 1, IX[1] ) }, CF[2];
    cputime = -mc::cpuclock();
    for( unsigned i=0; i<NREP; i++ )
      DAG.eval( F_op, NF, F, CF, NX, X, CX );
    cputime += mc::cpuclock();
    std::cout << "\nDAG " << NTE << "th-order dense Chebyshev model evaluation - no preallocation, no variadic template: " << (cputime/=NREP) << " CPU-sec\n";
    for( unsigned i=0; i<NF; i++ ) std::cout << "F[" << i << "] = " << CF[i].R() << std::endl;
  }

  return 0;
}

///////////////////////////////////////////////////////////////////////////////

int main()
{
  try{
    test_eval1();
    test_eval2();
  }
  catch( mc::FFGraph::Exceptions &eObj ){
    std::cerr << "Error " << eObj.ierr()
              << " in factorable function manipulation:" << std::endl
              << eObj.what() << std::endl
              << "Aborts." << std::endl;
    return eObj.ierr();
  }
#ifndef USE_PROFIL
#ifndef USE_FILIB
  catch( I::Exceptions &eObj ){
    std::cerr << "Error " << eObj.ierr()
              << " in interval arithmetic:" << std::endl
              << eObj.what() << std::endl
              << "Aborts." << std::endl;
    return eObj.ierr();
  }
#endif
#endif
  catch( SCM::Exceptions &eObj ){
    std::cerr << "Error " << eObj.ierr()
              << " in sparse Chebyshev model arithmetic:" << std::endl
              << eObj.what() << std::endl
              << "Aborts." << std::endl;
    return eObj.ierr();
  }
  catch( CM::Exceptions &eObj ){
    std::cerr << "Error " << eObj.ierr()
              << " in dense Chebyshev model arithmetic:" << std::endl
              << eObj.what() << std::endl
              << "Aborts." << std::endl;
    return eObj.ierr();
  }
}

