////////////////////////////////////////////////////////////////////////
#undef USE_PROFIL	// <-- specify to use PROFIL for interval arithmetic
#undef USE_FILIB	// <-- specify to use FILIB++ for interval arithmetic
////////////////////////////////////////////////////////////////////////

#include <fstream>
#include <iomanip>

#include "mctime.hpp"
#include "ffunc.hpp"
#include "rltred.hpp"

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

#include "scmodel.hpp"
typedef mc::SCModel<I> SCM;
typedef mc::SCVar<I> SCV;

///////////////////////////////////////////////////////////////////////////////

int test_DAG()
{
  const unsigned int NX = 4, NF = 2;
  mc::FFGraph DAG;
  mc::FFVar X[NX];

  // DAG construct and use
  for( unsigned int i=0; i<NX; i++ ) X[i].set( &DAG );
  mc::FFVar F[NF]
    = { X[2]*X[3]-X[0],
        X[0]*pow(exp(X[2]*X[3])+3.,4)+X[1] };
  std::cout << DAG;

  DAG.output( DAG.subgraph( NF, F ) );
  std::ofstream o_F( "F.dot", std::ios_base::out );
  DAG.dot_script( NF, F, o_F );
  o_F.close();

  DAG.output( DAG.subgraph( 1, &F[0] ) );
  std::ofstream o_F0( "F0.dot", std::ios_base::out );
  DAG.dot_script( 1, F, o_F0 );
  o_F0.close();

  for( unsigned i=0; i<NF; i++ )
    std::cout << "\nVariable dependence of F[" << i << "]: "
              << F[i].dep() << std::endl;

  // 1st-order derivative DAG construct using FAD
  const mc::FFVar* dFdX_FAD = DAG.FAD( NF, F, NX, X );
  std::cout << DAG;

  DAG.output( DAG.subgraph( NX*NF, dFdX_FAD ) );
  std::ofstream o_dFdX_FAD( "dFdX_FAD.dot", std::ios_base::out );
  DAG.dot_script( NX*NF, dFdX_FAD, o_dFdX_FAD );
  o_dFdX_FAD.close();
      
  DAG.output( DAG.subgraph( 1, &dFdX_FAD[NX+3] ) );
  std::ofstream o_dF1dX3_FAD( "dF1dX3_FAD.dot", std::ios_base::out );
  DAG.dot_script( 1, &dFdX_FAD[NX+3], o_dF1dX3_FAD );
  o_dF1dX3_FAD.close();


  // 1st-order derivative DAG construct using BAD
  const mc::FFVar* dFdX_BAD = DAG.BAD( NF, F, NX, X );
  std::cout << DAG;

  DAG.output( DAG.subgraph( NX*NF, dFdX_BAD ) );
  std::ofstream o_dFdX_BAD( "dFdX_BAD.dot", std::ios_base::out );
  DAG.dot_script( NX*NF, dFdX_BAD, o_dFdX_BAD );
  o_dFdX_BAD.close();
      
  DAG.output( DAG.subgraph( 1, &dFdX_BAD[NX+3] ) );
  std::ofstream o_dF1dX3_BAD( "dF1dX3_BAD.dot", std::ios_base::out );
  DAG.dot_script( 1, &dFdX_BAD[NX+3], o_dF1dX3_BAD );
  o_dF1dX3_BAD.close();

  delete[] dFdX_FAD;
  delete[] dFdX_BAD;

  return 0;
}

///////////////////////////////////////////////////////////////////////////////

int test_TAD()
{
  mc::FFGraph DAG;
  mc::FFVar T( &DAG );
  mc::FFVar X( &DAG ); 
  mc::FFVar F = sqr(X);
  std::cout << DAG;

  const unsigned NTE = 5;
  const mc::FFVar* F_TAD = DAG.TAD( NTE, 1, &F, 1, &X, &T );
  std::cout << DAG;

  DAG.output( DAG.subgraph( NTE+1, F_TAD ) );
  std::ofstream o_TAD( "FTE.dot", std::ios_base::out );
  DAG.dot_script( NTE+1, F_TAD, o_TAD );
  o_TAD.close();

  delete[] F_TAD;
  return 0;
}

///////////////////////////////////////////////////////////////////////////////

int test_eval()
{
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
  const mc::FFVar* dFdXdir = DAG.FAD( NF, F, NX, X, D );

  // Evaluation in interval arithmetic
  I IX[NX] = { I(0,0.5), I(1,2), I(-1,-0.8), I(0.5,1) }, IdFdXdir[NF];
  DAG.eval( NF, dFdXdir, IdFdXdir, NX, X, IX );

  // Display results
  for( unsigned i=0; i<NF; i++ )
    std::cout << "  dF("<< i << ")dX·D = " << IdFdXdir[i] << std::endl;

  // Evaluation in 3rd-order Chebyshev model arithmetic
  const unsigned ORD = 3;
  SCM CMenv( ORD );
  SCV CMX[NX], CMdFdXdir[NF];
  for( unsigned i=0; i<NX; i++ ) CMX[i].set( &CMenv, i, IX[i] );
  DAG.eval( NF, dFdXdir, CMdFdXdir, NX, X, CMX );

  // Display results
  for( unsigned i=0; i<NF; i++ )
    std::cout << "  dF("<< i << ")dX·D = " << CMdFdXdir[i] << std::endl;

  return 0;
}

///////////////////////////////////////////////////////////////////////////////

int test_DAG2()
{
  // Create DAG for f1,f2
  const unsigned int NX = 2, NF = 2;
  mc::FFGraph DAG;
  mc::FFVar X[NX];
  for( unsigned int i=0; i<NX; i++ ) X[i].set( &DAG );
  mc::FFVar F[NF] = { sqrt(X[0])*exp(X[1])*X[0], pow(X[1],3)*sqrt(X[0]) };
  std::cout << DAG;

  // Display DAG for f1,f2
  std::list<const mc::FFOp*> F_op  = DAG.subgraph( NF, F );    DAG.output( F_op );
  std::ofstream o_F( "F.dot", std::ios_base::out );
  DAG.dot_script( NF, F, o_F );
  o_F.close();

  // Derivatives
  const mc::FFVar* dFdX = DAG.FAD( NF, F, NX, X );
  std::ofstream o_dFdX( "dFdX.dot", std::ios_base::out );
  DAG.dot_script( NF*NX, dFdX, o_dFdX );

  // Evaluate
  double dX[NX] = { 1., 2. }, dF[2];
  DAG.eval( NF, F, dF, NX, X, dX );
  std::cout << "F1 = " << dF[0] << std::endl;
  std::cout << "F2 = " << dF[1] << std::endl;

  double ddFdX[NX*NF];
  DAG.eval( NF*NX, dFdX, ddFdX, NX, X, dX );
  std::cout << "dF1dX1 = " << ddFdX[0] << std::endl;
  std::cout << "dF2dX1 = " << ddFdX[2] << std::endl;
  std::cout << "dF1dX2 = " << ddFdX[1] << std::endl;
  std::cout << "dF2dX2 = " << ddFdX[3] << std::endl;

  // Evaluate
  I IX[NX] = { I(1.,2.), I(2.,3.) }, IF[2];
  DAG.eval( NF, F, IF, NX, X, IX );
  std::cout << "F1 = " << IF[0] << std::endl;
  std::cout << "F2 = " << IF[1] << std::endl;

/*
  // Display DAG for f1 alone
  std::list<const mc::FFOp*> F0_op = DAG.subgraph( 1, F );     DAG.output( F0_op );
  std::ofstream o_F0( "F0.dot", std::ios_base::out );
  DAG.dot_script( 1, F, o_F0 );
  o_F0.close();

  // Evaluate DAG for f1,f2 in double arithmetic
  double F_d[NF], X_d[NX] = { 2., -1. };
  DAG.eval( F_op, NF, F, F_d, NX, X, X_d );
  std::cout <<"\nDouble Value:\n";
  for( unsigned int i=0; i<NX; i++ ) std::cout << "F[" << i << "] = " << F_d[i] << std::endl;

  // Evaluate DAG for f1,f2 in ellipsoidal arithmetic
  dcv X_c(NX);    dsm X_Q(NX);
  X_c(0) =  2.;   X_Q(0,0) = 0.1;
  X_c(1) = -1.;   X_Q(1,0) = 0.04; X_Q(1,1) = 0.1; 
  EI X_EI( X_Q, X_c, DAG.depmap( F_op, NF, F, NX, X ) );
  EV F_EV[NF], X_EV[NX];
  for( unsigned i=0; i<NX; ++i ) X_EV[i].set( X_EI, i );
  DAG.eval( F_op, NF, F, F_EV, NX, X, X_EV );
  EI F_EI = X_EI.get( NF, F_EV ) ;
  std::cout <<"\nEllipsoidal Image:" << F_EI << std::endl;


  // Evaluate DAG for f1,f2 in Taylor model with ellipsoidal arithmetic
  X_EI.set( X_Q, X_c );
  for( unsigned i=0; i<NX; ++i ) X_EV[i].set( X_EI, i );
  std::cout << X_EV[0] << std::endl;
  TME TM_env( NX, 1 );
  TVE X_TVE[NX], F_TVE[NF];
  for( unsigned int i=0; i<NX; i++ )
     X_TVE[i].set( &TM_env, i, X_EV[i] );
  //X_EI.output();
  //std::cout << X_TVE[0];
  DAG.eval( F_op, NF, F, F_TVE, NX, X, X_TVE );
  EI RF_EI = X_EI.get( NF, F_EV ) ;
  std::cout <<"\nEllipsoidal Image:" << F_EI << std::endl;
*/

  return 0;
}

///////////////////////////////////////////////////////////////////////////////

int test_composition()
{
  mc::FFGraph DAG;
  mc::FFVar X, Y, F, G;
  X.set( &DAG );
  Y.set( &DAG );
  F = exp(X);
  G = sqr(Y)+F;
  std::cout << DAG;

  std::ofstream o_comp0( "compose0.dot", std::ios_base::out );
  DAG.dot_script( 1, &G, o_comp0 );
  o_comp0.close();

  const mc::FFVar* GoF = DAG.compose( 1, &G, 1, &Y, &F );
  std::cout << DAG;

  std::ofstream o_comp1( "compose1.dot", std::ios_base::out );
  DAG.dot_script( 1, GoF, o_comp1 );
  o_comp1.close();

  delete[] GoF;
  return 0;
}

///////////////////////////////////////////////////////////////////////////////

int test_FAD_sparse()
{
  mc::FFGraph DAG;
  const unsigned NX = 5, NF = 2;
  mc::FFVar X[NX];
  for( unsigned i(0); i<NX; i++ )  X[i].set( &DAG );
  mc::FFVar F[NF] = { sqrt(X[0])*X[1], X[3]+2. };
  std::cout << DAG;

  auto dFdX = DAG.SBAD( NF, F, NX, X );
  //std::cout << DAG;
  std::cout << "\nNumber of non-zero elements in Jacobian matrix: "
            << std::get<0>(dFdX) << std::endl;
  for( unsigned ie=0; ie<std::get<0>(dFdX); ie++ )
    std::cout << "(" << std::get<1>(dFdX)[ie] << "," << std::get<2>(dFdX)[ie] << ") "
              << std::get<3>(dFdX)[ie] << std::endl;

  DAG.output( DAG.subgraph( std::get<0>(dFdX), std::get<3>(dFdX) ) );

  std::ofstream o_dFdX( "sparseFAD.dot", std::ios_base::out );
  DAG.dot_script( std::get<0>(dFdX), std::get<3>(dFdX), o_dFdX );
  o_dFdX.close();

  auto d2FdX2 = DAG.SFAD( std::get<0>(dFdX), std::get<3>(dFdX), NX, X );
  //std::cout << DAG;
  std::cout << "\nNumber of non-zero elements in Hessian matrix: "
            << std::get<0>(d2FdX2) << std::endl;
  for( unsigned ie=0; ie<std::get<0>(d2FdX2); ie++ )
    std::cout << "(" << std::get<1>(d2FdX2)[ie] << "," << std::get<2>(d2FdX2)[ie] << ") "
              << std::get<3>(d2FdX2)[ie] << std::endl;

  DAG.output( DAG.subgraph( std::get<0>(d2FdX2), std::get<3>(d2FdX2) ) );

  std::ofstream o_d2FdX2( "sparseFAD2.dot", std::ios_base::out );
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

int test_FAD_directional()
{
  mc::FFGraph DAG;
  mc::FFVar X[2], F[2], D[2];
  X[0].set( &DAG );
  X[1].set( &DAG );
  F[0] = sqrt(X[0])*X[1];
  F[1] = X[0]*X[1];
  //D[0].set( &DAG, 1. );
  //D[1].set( &DAG, 1. );
  //F = X[0]*X[1];
  D[0] = 1.;
  D[1] = 1.;
  std::cout << DAG;

  const mc::FFVar* dFdXxD = DAG.FAD( 2, F, 2, X, D );
  std::cout << DAG;
  DAG.output( DAG.subgraph( 1, dFdXxD ) );
  DAG.output( DAG.subgraph( 1, dFdXxD+1 ) );

  std::ofstream o_dFdXxD( "directional.dot", std::ios_base::out );
  DAG.dot_script( 1, dFdXxD, o_dFdXxD );
  o_dFdXxD.close();

  delete[] dFdXxD;
  return 0;
}

///////////////////////////////////////////////////////////////////////////////

int test_const()
{
  mc::FFGraph DAG;
  std::vector<mc::FFVar> X(2);
  mc::FFVar F, PI;
  X[0].set( &DAG );
  X[1].set( &DAG );
  PI.set( &DAG, mc::PI );
  F = exp(X[0])*X[1]/PI;
  std::cout << DAG;

  std::ofstream o_cst1( "constant.dot", std::ios_base::out );
  DAG.dot_script( 1, &F, o_cst1 );
  o_cst1.close();

  // Evaluate DAG for F in interval arithmetic
  I X_I[2] = { I(-1.,1.), I(0.,2.) }, F_I;
  DAG.eval( 1, &F, &F_I, 2, X.data(), X_I );
  std::cout <<"\nInterval bound:\n";
  std::cout << "F = " << F_I << std::endl;

  X[1].set(1);
  DAG.eval( 1, &F, &F_I, 1, X.data(), X_I );
  std::cout <<"\nInterval bound:\n";
  std::cout << "F = " << F_I << std::endl;

  //X[1].unset();
  DAG.eval( 1, &F, &F_I, 2, X.data(), X_I );
  std::cout <<"\nInterval bound:\n";
  std::cout << "F = " << F_I << std::endl;

  const mc::FFVar* dFdX = DAG.FAD( 1, &F, 1, X.data() );
  std::cout << DAG;

  std::ofstream o_cst2( "compose2.dot", std::ios_base::out );
  DAG.dot_script( 1, dFdX, o_cst2 );
  o_cst2.close();

  delete[] dFdX;
  return 0;
}

///////////////////////////////////////////////////////////////////////////////

int test_BAD()
{
  mc::FFGraph DAG;
  const unsigned NX = 2, NF = 2;
  mc::FFVar X[NX];
  for( unsigned i(0); i<NX; i++ )  X[i].set( &DAG );
  mc::FFVar F[NF] = { X[0], 1. };//exp(X[1])*X[0] };
  std::cout << DAG;

  fadbad::B<mc::FFVar> X0_B = X[0];
  fadbad::B<mc::FFVar> X1_B = X[1];
  fadbad::B<mc::FFVar> F0_B = X0_B;
  fadbad::B<mc::FFVar> F1_B = 1.;//X1_B;//exp(X1_B)*X0_B;
  F0_B.diff(0,2);
  F1_B.diff(1,2);
  //F_B = 0.;
  X0_B.d(0);  
  X1_B.d(0);
  std::cout << X0_B.d(0) << std::endl;
  std::cout << X0_B.d(1) << std::endl;
  std::cout << X1_B.d(0) << std::endl;
  std::cout << X1_B.d(1) << std::endl;
  //return 0;

  fadbad::B<mc::FFVar> X_B[NX], F_B[NF];
  X_B[0] = X[0]; X_B[1] = X[1];
  DAG.eval( NF, F, F_B, NX, X, X_B );
  F_B[0].diff(0,NF);
  F_B[1].diff(1,NF);
  std::cout << X_B[0].d(0) << std::endl;
  std::cout << X_B[0].d(1) << std::endl;
  std::cout << X_B[1].d(0) << std::endl;
  std::cout << X_B[1].d(1) << std::endl;
  //return 0;

  auto dFdX = DAG.SBAD( NF, F, NX, X );
  std::cout << DAG;
  std::cout << "\nNumber of non-zero elements in Jacobian matrix: "
            << std::get<0>(dFdX) << std::endl;
  for( unsigned ie=0; ie<std::get<0>(dFdX); ie++ )
    std::cout << "(" << std::get<1>(dFdX)[ie] << "," << std::get<2>(dFdX)[ie] << ") "
              << std::get<3>(dFdX)[ie] << std::endl;

  DAG.output( DAG.subgraph( std::get<0>(dFdX), std::get<3>(dFdX) ) );

  std::ofstream o_dFdX( "sparseBAD.dot", std::ios_base::out );
  DAG.dot_script( std::get<0>(dFdX), std::get<3>(dFdX), o_dFdX );
  o_dFdX.close();

  delete[] std::get<1>(dFdX);
  delete[] std::get<2>(dFdX);
  delete[] std::get<3>(dFdX);

  return 0;
}

///////////////////////////////////////////////////////////////////////////////

int test_rltred1()
{
  mc::FFGraph DAG;
  const unsigned NX = 5, NF = 3;
  mc::FFVar X[NX], &FT = X[0], &F1 = X[1], &F2 = X[2], &x1 = X[3], &x2 = X[4];
  for( unsigned i(0); i<NX; i++ )  X[i].set( &DAG );
  mc::FFVar F[NF];
  F[0] = x1 * FT - F1;
  F[1] = x2 * FT - F2;
  //F[2] = x1 + x2 -1.;
  F[2] = F1 + F2 - FT;
  std::cout << DAG;

  mc::RLTRED RRLT( &DAG );
  RRLT.options.DISPLAY = 1;
  RRLT.options.LEVEL   = mc::RLTRED::Options::FULLSEQ;
  RRLT.options.NODIV   = false;
  RRLT.search( NF, F );

  return 0;
}

///////////////////////////////////////////////////////////////////////////////

int test_rltred2()
{
  mc::FFGraph DAG;
  const unsigned NX = 2, NF = 2;
  mc::FFVar X[NX];
  for( unsigned i(0); i<NX; i++ ) X[i].set( &DAG );
  mc::FFVar F[NF];
  F[0] = X[0] + X[1] - 1.;
  F[1] = sqr(X[0]) - sqr(X[1]);
  std::cout << DAG;

  mc::RLTRED RRLT( &DAG );
  RRLT.options.DISPLAY = 2;
  RRLT.options.LEVEL   = mc::RLTRED::Options::FULLSEQ;
  RRLT.options.NODIV   = false;
  RRLT.search( NF, F );

  return 0;
}

///////////////////////////////////////////////////////////////////////////////

int main()
{
  try{
    //test_DAG();
    //test_DAG2();
    //test_TAD();
    //test_comp();
    //test_dirder();
    //test_sparseder();
    //test_badiff();
    test_eval();
    //test_rltred1();
    //test_rltred2();
  }

#ifndef USE_PROFIL
#ifndef USE_FILIB
  catch( I::Exceptions &eObj ){
    std::cerr << "Error " << eObj.ierr()
              << " in natural interval extension:" << std::endl
              << eObj.what() << std::endl
              << "Aborts." << std::endl;
    return eObj.ierr();
  }
#endif
#endif
  catch( SCM::Exceptions &eObj ){
    std::cerr << "Error " << eObj.ierr()
              << " in Taylor model computation:" << std::endl
              << eObj.what() << std::endl
              << "Aborts." << std::endl;
    return eObj.ierr();
  }
  catch( mc::FFGraph::Exceptions &eObj ){
    std::cerr << "Error " << eObj.ierr()
              << " in factorable function manipulation:" << std::endl
              << eObj.what() << std::endl
              << "Aborts." << std::endl;
    return eObj.ierr();
  }
}

