#undef MC__SELIM_DEBUG_PROCESS

#include <fstream>
#include <sstream>
#include <iomanip>

#include "ffunc.hpp"
#include "ffinv.hpp"
#include "ffexpr.hpp"
#include "slift.hpp"
#include "selim.hpp"

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

///////////////////////////////////////////////////////////////////////////////

int test_dep1()
{
  std::cout << "\n==============================================\ntest_dep1:\n";

  const int NX = 4;
  mc::FFDep X[NX];
  for( int i=0; i<NX; i++ ) X[i].indep(i);

  const int NF = 2;
  mc::FFDep F[NF] = { X[2]*X[3]+X[0]/X[2],
                      sqr(X[0])*(exp(X[2])+X[3])+X[1] };
  
  std::cout << "Dependence structure of F[0]: " << F[0] << std::endl;
  std::cout << "Dependence structure of F[1]: " << F[1] << std::endl;

  return 0;
}

///////////////////////////////////////////////////////////////////////////////

int test_dep2()
{
  std::cout << "\n==============================================\ntest_dep2:\n";

  // Create DAG
  const unsigned NX = 4, NF = 1;
  mc::FFGraph DAG;
  mc::FFVar X[NX];
  for( unsigned int i=0; i<NX; i++ ) X[i].set( &DAG );
  //mc::FFVar F[NF] = { X[0]*sqr(X[3])};
  mc::FFVar F[NF] = { X[2]*X[3] + sqr(X[1]) + X[0]*sqr(X[3])};
  //mc::FFVar F[NF] = { X[2]*X[3]+X[0]/X[2],
  //                    X[0]*pow(exp(X[2]-X[3]),2)+X[1] };
  std::cout << DAG;

  std::ofstream o_F( "dep2_F.dot", std::ios_base::out );
  DAG.dot_script( NF, F, o_F );
  o_F.close();

  // Evaluate with dependents
  auto F_op  = DAG.subgraph( NF, F );
  mc::FFDep depX[NX], depF[NF];
  for( unsigned int i=0; i<NX; i++ ) depX[i].indep(i);
  DAG.eval( F_op, NF, F, depF, NX, X, depX );

  for( unsigned int i=0; i<NF; i++ )
    std::cout << "Dependence structure of F[" << i << "]: " << depF[i] << std::endl;

  return 0;
}

///////////////////////////////////////////////////////////////////////////////

int test_expr1()
{
  std::cout << "\n==============================================\ntest_expr1:\n";

  // Create DAG
  const unsigned NX = 4, NF = 2;
  mc::FFGraph DAG;
  mc::FFVar X[NX];
  for( unsigned int i=0; i<NX; i++ ) X[i].set( &DAG );
  mc::FFVar F[NF] = { -1 + log( X[2]/(-8*X[3]) - 2 * sqr(X[1]) + X[0] ) * sqr(X[3]) - 2.,
                      X[2]*X[3]+X[0]/X[2] };
  std::cout << DAG;

  std::ofstream o_F( "expr1_F.dot", std::ios_base::out );
  DAG.dot_script( NF, F, o_F );
  o_F.close();

  // Evaluate with dependents
  auto Fop  = DAG.subgraph( NF, F );
  mc::FFExpr EX[NX], EF[NF];
  for( unsigned int i=0; i<NX; i++ ) EX[i].set( X[i] );
  DAG.eval( Fop, NF, F, EF, NX, X, EX );

  for( unsigned int i=0; i<NF; i++ )
    std::cout << "Expression of F[" << i << "]: " << EF[i] << std::endl;

  return 0;
}

///////////////////////////////////////////////////////////////////////////////

int test_inv1()
{
  std::cout << "\n==============================================\ntest_inv1:\n";

  const int NX = 4;
  mc::FFInv X[NX];
  for( int i=0; i<NX; i++ ) X[i].indep(i);

  const int NF = 2;
  mc::FFInv F[NF] = { X[2]*X[3]+X[0]/X[2],
                      sqr(X[0])*(exp(X[2])+X[3])+X[1] };
  
  std::cout << "Invertible structure of F[0]: " << F[0] << std::endl;
  std::cout << "Invertible structure of F[1]: " << F[1] << std::endl;

  return 0;
}

///////////////////////////////////////////////////////////////////////////////

int test_inv2()
{
  std::cout << "\n==============================================\ntest_inv2:\n";

  // Create DAG
  const unsigned NX = 4, NF = 1;
  mc::FFGraph DAG;
  mc::FFVar X[NX];
  for( unsigned int i=0; i<NX; i++ ) X[i].set( &DAG );
  mc::FFVar F[NF] = { X[2]*X[3] + sqr(X[1]) + X[0]*sqr(X[3])};
  std::cout << DAG;

  std::ofstream o_F( "inv2_F.dot", std::ios_base::out );
  DAG.dot_script( NF, F, o_F );
  o_F.close();

  // Evaluate with dependents
  auto F_op  = DAG.subgraph( NF, F );
  mc::FFInv invX[NX], invF[NF];
  for( unsigned int i=0; i<NX; i++ ) invX[i].indep(i);
  DAG.eval( F_op, NF, F, invF, NX, X, invX );

  for( unsigned int i=0; i<NF; i++ )
    std::cout << "Invertible structure of F[" << i << "]: " << invF[i] << std::endl;

  return 0;
}

///////////////////////////////////////////////////////////////////////////////

int test_spoly1()
{
  std::cout << "\n==============================================\ntest_spoly1:\n";

  mc::FFGraph DAG;
  const unsigned NX = 4;
  mc::FFVar X[NX];
  typedef mc::SPoly< mc::FFVar const*, mc::lt_FFVar > t_SPoly;
  t_SPoly SPX[NX];
  for( unsigned i(0); i<NX; i++ ){
    X[i].set( &DAG );
    SPX[i].var( &X[i] );
  }
  //mc::FFVar F = pow( X[0] - sqr( X[1] ) - 1, 3 );
  //mc::FFVar F = pow( 2 + X[0] + X[1] - X[2] + (-X[3]) + X[0] - 3, 4 );
  mc::FFVar F = sqr( sqr( 2 + X[0] + X[1] - X[2] + (-X[3]) + X[0] - 3 ) );
  t_SPoly SPF;

  SPF.options.BASIS = t_SPoly::Options::MONOM;
  DAG.eval( 1, &F, &SPF, NX, X, SPX );
  std::cout << std::endl << "Sparse polynomial expression in monomial basis: " << SPF << std::endl;

  SPF.options.BASIS = t_SPoly::Options::CHEB;
  DAG.eval( 1, &F, &SPF, NX, X, SPX );
  std::cout << std::endl << "Sparse polynomial expression in chebyshev basis: " << SPF << std::endl;

  return 0;
}

///////////////////////////////////////////////////////////////////////////////

int test_slift0()
{
  std::cout << "\n==============================================\ntest_slift0:\n";

  mc::FFGraph DAG;
  const unsigned NX = 2, NF = 2;
  mc::FFVar X[NX];
  for( unsigned i(0); i<NX; i++ ) X[i].set( &DAG );
  mc::FFVar F[NF];
  F[0] = pow( X[0] + 1 / sqr( X[1] ), 3 );
  F[1] = exp( 2 * sqr( X[1] ) - 1 );
  std::cout << DAG;

  //mc::SLiftVar::t_poly::options.BASIS = mc::SLiftVar::t_poly::Options::MONOM;
  mc::SLiftEnv SPE( &DAG );
  //SPE.options.LIFTDIV = true;//false;//
  //SPE.options.LIFTIPOW = false;//
  SPE.process( NF, F );
  //std::cout << SPE;

  std::cout << std::endl << SPE.Var().size() << " participating variables: ";
  for( auto&& var : SPE.Var() ) std::cout << var << " ";
  std::cout << std::endl;
  std::cout << std::endl << SPE.Dep().size() << " dependent expressions: " << std::endl;;
  for( auto&& expr : SPE.Dep() ) DAG.output( DAG.subgraph( 1, &expr ) );
  std::cout << std::endl << SPE.Poly().size() << " polynomial constraints: " << std::endl;
  for( auto&& expr : SPE.Poly() ) DAG.output( DAG.subgraph( 1, &expr ) );
  std::cout << std::endl << SPE.Trans().size() << " transcendental constraints: " << std::endl;
  for( auto&& expr : SPE.Trans() ) DAG.output( DAG.subgraph( 1, &expr ) );
  std::cout << std::endl << SPE.Aux().size() << " auxiliary variables: ";
  for( auto&& aux : SPE.Aux() ) std::cout << *aux.first << "->" << *aux.second << " ";
  std::cout << std::endl;
  
  return 0;
}

///////////////////////////////////////////////////////////////////////////////

int test_slift1()
{
  std::cout << "\n==============================================\ntest_slift1:\n";

  mc::FFGraph DAG;
  const unsigned NX = 4, NF = 1;
  mc::FFVar X[NX];
  for( unsigned i(0); i<NX; i++ ) X[i].set( &DAG );
  mc::FFVar F[NF];
  //F[0] = pow( X[0] + X[1], 3 );
  //F[0] = 1/X[0] + 1/X[1] - 1;
  //F[0] = pow( X[0] + 1 / sqr( X[1] ), 3 );
  //F[0] = 4*sqr(X[0]) - 2.1*pow(X[0],4) + pow(X[0],6)/3 + X[0]*X[1] - 4*sqr(X[1]) + 4*pow(X[1],4);
  //F[0] = 250*exp(X[2])*X[0] + 250*pow(X[3],0.6)*X[1];
  F[0] = exp( 2 * sqr( X[1] ) - 1 );
  std::cout << DAG;

  mc::SLiftVar::t_poly::options.BASIS = mc::SLiftVar::t_poly::Options::MONOM;
  mc::SLiftEnv SPL( &DAG );
  SPL.options.LIFTDIV = true;//false;//
  SPL.options.LIFTIPOW = false;//

  SPL.process( NF, F, true );//false );
  std::cout << SPL;
  
  return 0;
}

///////////////////////////////////////////////////////////////////////////////

int test_selim0()
{
  std::cout << "\n==============================================\ntest_selim0:\n";

  mc::FFGraph DAG;
  const unsigned NX = 3, NF = 2;
  mc::FFVar X[NX];
  for( unsigned i(0); i<NX; i++ ) X[i].set( &DAG );
  mc::FFVar F[NF];
  F[0] = ( 3. * X[0] * sqr( X[2] ) ) / X[1] - 2. * X[0] * X[1] - X[0] - 1;
  F[1] = 2./X[1] + 3./X[2] - 1.;
  std::cout << DAG;

  mc::SElimEnv SPE( &DAG );
  SPE.options.ELIMMLIN      = true;
  SPE.options.MIPDISPLEVEL  = 1;
  SPE.options.MIPOUTPUTFILE = "test_selim0.lp";

  SPE.process( NF, F );
  std::cout << SPE;
  
  return 0;
}

///////////////////////////////////////////////////////////////////////////////

int test_selim1()
{
  std::cout << "\n==============================================\ntest_selim1:\n";

  mc::FFGraph DAG;
  const unsigned NX = 5, NF = 4;
  mc::FFVar X[NX];
  for( unsigned i(0); i<NX; i++ ) X[i].set( &DAG );
  mc::FFVar F[NF];
  F[0] = X[0] + X[1] * sqr( X[3] ) + X[1] * exp( X[2] );
  F[1] = exp( 2 * sqrt( X[3] - X[1] ) ) - 2;
  F[2] = sqr( X[0] ) + sqrt( X[1] ) - 1;
  F[3] = sum( NX, X ) - 1;
  std::cout << DAG;

  mc::SElimEnv SPE( &DAG );
  SPE.options.ELIMMLIN      = true;
  SPE.options.MIPDISPLEVEL  = 0;
  SPE.options.MIPOUTPUTFILE = "test_selim1.lp";

  SPE.process( NF, F );
  std::cout << SPE;
  
  return 0;
}

///////////////////////////////////////////////////////////////////////////////

int test_selim2()
{
  std::cout << "\n==============================================\ntest_selim2:\n";

  mc::FFGraph DAG;
  const unsigned NX = 4, NF = 1;
  mc::FFVar X[NX];
  for( unsigned i(0); i<NX; i++ ) X[i].set( &DAG );
  mc::FFVar F[NF];
  F[0] = pow( X[0] + 1 / sqr( X[1] ), 3 );
  std::cout << DAG;

  mc::SLiftVar::t_poly::options.BASIS = mc::SLiftVar::t_poly::Options::MONOM;
  mc::SLiftEnv SPL( &DAG );
  SPL.options.LIFTDIV = true;//false;//
  SPL.options.LIFTIPOW = true; //false;//

  SPL.process( NF, F );
  std::cout << SPL;

  std::vector<mc::FFVar> Flift( SPL.Poly() );
  Flift.insert( Flift.end(), SPL.Trans().begin(), SPL.Trans().end() );

  mc::SElimEnv SPE( &DAG );
  SPE.options.ELIMMLIN      = true;
  SPE.options.MIPDISPLEVEL  = 1;
  SPE.options.MIPOUTPUTFILE = "test_selim2.lp";

  SPE.process( Flift.size(), Flift.data() );
  std::cout << SPE;

  return 0;
}

///////////////////////////////////////////////////////////////////////////////

int test_selim3()
{
  std::cout << "\n==============================================\ntest_selim3:\n";

  mc::FFGraph DAG;
  const unsigned NX = 4, NF = 1;
  mc::FFVar X[NX];
  for( unsigned i(0); i<NX; i++ ) X[i].set( &DAG );
  mc::FFVar F[NF];
  F[0] = (X[1]+1)*X[3] - 1;
  //F[0] = (X[1]*X[2] + 1e3*X[0])*X[3] - 1;
  std::cout << DAG;

  mc::SElimEnv SPE( &DAG );
  SPE.options.ELIMMLIN      = true;
  SPE.options.MIPDISPLEVEL  = 0;
  SPE.options.MIPOUTPUTFILE = "test_selim3.lp";

  SPE.process( NF, F );
  std::cout << SPE;
  
  return 0;
}

///////////////////////////////////////////////////////////////////////////////

int main()
{
  try{
    test_expr1();
    test_dep1();
    test_inv1();
    test_dep2();
    test_inv2();
    test_spoly1();
    test_slift0();
    test_slift1();
    test_selim0();
    test_selim1();
    test_selim2();
    test_selim3();
  }
  catch( mc::FFBase::Exceptions &eObj ){
    std::cerr << "Error " << eObj.ierr()
              << " in factorable function manipulation:" << std::endl
              << eObj.what() << std::endl
              << "Aborts." << std::endl;
    return eObj.ierr();
  }
  catch( mc::SLiftEnv<mc::FFGraph<>>::Exceptions &eObj ){
    std::cerr << "Error " << eObj.ierr()
              << " in variable lifting manipulation:" << std::endl
              << eObj.what() << std::endl
              << "Aborts." << std::endl;
    return eObj.ierr();
  }
  catch( mc::SElimEnv<mc::FFGraph<>>::Exceptions &eObj ){
    std::cerr << "Error " << eObj.ierr()
              << " in variable elimination manipulation:" << std::endl
              << eObj.what() << std::endl
              << "Aborts." << std::endl;
    return eObj.ierr();
  }
}

