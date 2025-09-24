////////////////////////////////////////////////////////////////////////
#undef  MC__FFUNC_CPU_EVAL
#undef  MC__SCMODEL_TRACE
#undef  MC__CMODEL_TRACE
#undef  MC__USE_THREADLOCAL
#undef MC__VEVAL_DEBUG
////////////////////////////////////////////////////////////////////////

#include <fstream>
#include <iomanip>
#include <chrono>

#include "mctime.hpp"
#include "ffunc.hpp"

#ifdef MC__USE_PROFIL
 #include "mcprofil.hpp"
 typedef INTERVAL I;
#else
 #ifdef MC__USE_FILIB
  #include "mcfilib.hpp"
  typedef filib::interval<double,filib::native_switched,filib::i_mode_extended> I;
 #else
  #ifdef MC__USE_BOOST
   #include "mcboost.hpp"
   typedef boost::numeric::interval_lib::save_state<boost::numeric::interval_lib::rounded_transc_opp<double>> T_boost_round;
   typedef boost::numeric::interval_lib::checking_base<double> T_boost_check;
   typedef boost::numeric::interval_lib::policies<T_boost_round,T_boost_check> T_boost_policy;
   typedef boost::numeric::interval<double,T_boost_policy> I;
  #else
   #include "interval.hpp"
   typedef mc::Interval I;
  #endif
 #endif
#endif

#include "cmodel.hpp"
typedef mc::CModel<I> CM;
typedef mc::CVar<I> CV;

#include "scmodel.hpp"
typedef mc::SCModel<I> SCM;
typedef mc::SCVar<I> SCV;

I IINF = 1e20 * I(-1,1);

///////////////////////////////////////////////////////////////////////////////

int test_eval1()
{
  std::cout << "\n==============================================\ntest_eval1:\n";

  // Create DAG
  const unsigned NX = 4, NF = 2;
  mc::FFGraph DAG;
  std::vector<mc::FFVar> X(NX);
  for( auto& Xi : X ) Xi.set( &DAG );
  std::vector<mc::FFVar> F{
    X[2]*X[3]-2./(X[1]+X[2]),
    X[0]/pow(exp(X[2]*X[1])+3.,3)+tanh(X[3])
  };
  std::cout << DAG;

  std::ofstream o_F( "eval1_F.dot", std::ios_base::out );
  DAG.dot_script( F, o_F );
  o_F.close();

  double cputime;
  const unsigned NREP=1000;

  // Evaluate with doubles, no parameter pack
  auto&& F_op  = DAG.subgraph( F );
  std::vector<double> dX{ -1., -1., 2., 3. };
  for( unsigned i=0; i<NX; i++ ) std::cout << "X[" << i << "] = " << dX[i] << std::endl;

  std::vector<double> dF(NF);
  cputime = -mc::cpuclock();
  for( unsigned i=0; i<NREP; i++ ){
    dF[0] = dX[2]*dX[3]-2./(dX[1]+dX[2]);
    dF[1] = dX[0]/pow(exp(dX[2]*dX[1])+3.,3)+tanh(dX[3]);
  }
  cputime += mc::cpuclock();
  std::cout << "\nCompiled evaluation: " << (cputime/=NREP) << " CPU-sec\n";
  for( unsigned i=0; i<NF; i++ ) std::cout << "F[" << i << "] = " << dF[i] << std::endl;

  std::vector<double> WK;
  cputime = -mc::cpuclock();
  for( unsigned i=0; i<NREP; i++ )
    DAG.eval( F_op, WK, F, dF, X, dX );
  cputime += mc::cpuclock();
  std::cout << "\nDAG evaluation, w/ preallocation, w/o variadic template: " << (cputime/=NREP) << " CPU-sec\n";
  for( unsigned i=0; i<NF; i++ ) std::cout << "F[" << i << "] = " << dF[i] << std::endl;

  cputime = -mc::cpuclock();
  for( unsigned i=0; i<NREP; i++ )
    DAG.eval( F_op, WK, NF, &F[0], &dF[0], 1, &X[0], &dX[0], 1, &X[1], &dX[1], 1, &X[2], &dX[2], 1, &X[3], &dX[3] );
  cputime += mc::cpuclock();
  std::cout << "\nDAG evaluation, w/ preallocation, w/ variadic template: " << (cputime/=NREP) << " CPU-sec\n";
  for( unsigned i=0; i<NF; i++ ) std::cout << "F[" << i << "] = " << dF[i] << std::endl;

  return 0;
}

///////////////////////////////////////////////////////////////////////////////

int test_eval2()
{
  std::cout << "\n==============================================\ntest_eval2:\n";

  // Create DAG
  const unsigned NX = 4, NF = 2;
  mc::FFGraph DAG;
  DAG.options.USEMOVE = true;
  std::vector<mc::FFVar> X(NX);
  for( auto& Xi : X ) Xi.set( &DAG );
  std::vector<mc::FFVar> F{ X[2]*X[3]-2./(X[1]+X[2]),
                            X[0]/pow(exp(X[2]*X[1])+3.,3)+tanh(X[3]) };
//  std::vector<mc::FFVar> F{ sqrt(X[0])*exp(X[1])*X[0],
//                            pow(X[1],3)*sqrt(X[0]) };
  std::cout << DAG;

  auto F_op  = DAG.subgraph( F );
  DAG.output( F_op );
  std::ofstream o_F( "eval2_F.dot", std::ios_base::out );
  DAG.dot_script( F, o_F );
  o_F.close();

  double cputime;
  const unsigned NREP=1;//000;
  const unsigned NTEMIN=1, NTEMAX=7;

  // Evaluate in interval arithmetic
  std::vector<I> IWK;
  std::vector<I> IX{ I(-1.1,-0.9), I(-1.1, -0.9), I(1.6,2.4), I(2.5,3.5) }, IF;
  for( unsigned i=0; i<NX; i++ ) std::cout << "X[" << i << "] = " << IX[i] << std::endl;
  cputime = -mc::cpuclock();
  for( unsigned i=0; i<NREP; i++ )
    DAG.eval( F_op, IWK, F, IF, X, IX );
  cputime += mc::cpuclock();
  std::cout << "\nDAG interval evaluation, w/ preallocation, w/o variadic template: " << (cputime/=NREP) << " CPU-sec\n";
  for( unsigned i=0; i<NF; i++ ) std::cout << "F[" << i << "] = " << IF[i] << std::endl;

  for( unsigned NTE=NTEMIN; NTE<=NTEMAX; NTE++ ){
    SCM modSCM( NTE );
    //modSCM.options.BOUNDER_TYPE = SCM::Options::LSB;
    modSCM.options.REMEZ_USE = false;
    modSCM.options.MIXED_IA = false;
    std::vector<SCV> SCX{  SCV( &modSCM, 0, IX[0] ), SCV( &modSCM, 1, IX[1] ), SCV( &modSCM, 2, IX[2] ), SCV( &modSCM, 3, IX[3] ) }, SCF(NF);
    cputime = -mc::cpuclock();
    for( unsigned i=0; i<NREP; i++ ){
      SCF[0] = SCX[2]*SCX[3]-2./(SCX[1]+SCX[2]);
      SCF[1] = SCX[0]/pow(exp(SCX[2]*SCX[1])+3.,3)+tanh(SCX[3]);
    }
    cputime += mc::cpuclock();
    std::cout << "\nDAG " << NTE << "th-order sparse Chebyshev model evaluation, w/ compilation: " << (cputime/=NREP) << " CPU-sec\n";
    for( unsigned i=0; i<NF; i++ ) std::cout << "F[" << i << "] = " << SCF[i].P().B() << " +/- " << SCF[i].R() << std::endl;
  }

  for( unsigned NTE=NTEMIN; NTE<=NTEMAX; NTE++ ){
    std::vector<SCV> SCWK;
    SCM modSCM( NTE );
    //modSCM.options.BOUNDER_TYPE = SCM::Options::LSB;
    modSCM.options.REMEZ_USE = false;
    modSCM.options.MIXED_IA = false;
    std::vector<SCV> SCX{  SCV( &modSCM, 0, IX[0] ), SCV( &modSCM, 1, IX[1] ), SCV( &modSCM, 2, IX[2] ), SCV( &modSCM, 3, IX[3] ) }, SCF;
    cputime = -mc::cpuclock();
    for( unsigned i=0; i<NREP; i++ )
      DAG.eval( F_op, SCWK, F, SCF, X, SCX );
    cputime += mc::cpuclock();
    std::cout << "\nDAG " << NTE << "th-order sparse Chebyshev model evaluation, w/ preallocation, w/o variadic template: " << (cputime/=NREP) << " CPU-sec\n";
    for( unsigned i=0; i<NF; i++ ) std::cout << "F[" << i << "] = " << SCF[i].P().B() << " +/- " << SCF[i].R() << std::endl;
  }

  for( unsigned NTE=NTEMIN; NTE<=NTEMAX; NTE++ ){
    CM modCM( NX, NTE );
    //modSCM.options.BOUNDER_TYPE = SCM::Options::LSB;
    modCM.options.MIXED_IA = false;
    std::vector<CV> CX{  CV( &modCM, 0, IX[0] ), CV( &modCM, 1, IX[1] ), CV( &modCM, 2, IX[2] ), CV( &modCM, 3, IX[3] ) }, CF(NF);
    cputime = -mc::cpuclock();
    for( unsigned i=0; i<NREP; i++ ){
      CF[0] = CX[2]*CX[3]-2./(CX[1]+CX[2]);
      CF[1] = CX[0]/pow(exp(CX[2]*CX[1])+3.,3)+tanh(CX[3]);
    }
    cputime += mc::cpuclock();
    std::cout << "\nDAG " << NTE << "th-order dense Chebyshev model evaluation, w/ compilation: " << (cputime/=NREP) << " CPU-sec\n";
    for( unsigned i=0; i<NF; i++ ) std::cout << "F[" << i << "] = " << CF[i].P().B() << " +/- " << CF[i].R() << std::endl;
  }

  for( unsigned NTE=NTEMIN; NTE<=NTEMAX; NTE++ ){
    std::vector<CV> CWK;
    CM modCM( NX, NTE );
    //modSCM.options.BOUNDER_TYPE = SCM::Options::LSB;
    modCM.options.MIXED_IA = false;
    std::vector<CV> CX{  CV( &modCM, 0, IX[0] ), CV( &modCM, 1, IX[1] ), CV( &modCM, 2, IX[2] ), CV( &modCM, 3, IX[3] ) }, CF;
    cputime = -mc::cpuclock();
    for( unsigned i=0; i<NREP; i++ )
      DAG.eval( F_op, CWK, F, CF, X, CX );
    cputime += mc::cpuclock();
    std::cout << "\nDAG " << NTE << "th-order dense Chebyshev model evaluation, w/ preallocation, w/o variadic template: " << (cputime/=NREP) << " CPU-sec\n";
    for( unsigned i=0; i<NF; i++ ) std::cout << "F[" << i << "] = " << CF[i].P().B() << " +/- " << CF[i].R() << std::endl;
  }

  return 0;
}

///////////////////////////////////////////////////////////////////////////////

int test_eval3()
{
  std::cout << "\n==============================================\ntest_eval3:\n";

  // Create DAG for the residual r = ( y - 0.5 * x / exp(p) )^2
  mc::FFGraph DAG;
  mc::FFVar X( &DAG ), P( &DAG ), Y( &DAG );
  mc::FFVar F = sqr( Y - 0.5 * X * exp( P ) );
  std::cout << DAG;

  auto F_op  = DAG.subgraph( 1, &F );
  DAG.output( F_op );
  std::ofstream o_F( "eval3_F.dot", std::ios_base::out );
  DAG.dot_script( 1, &F, o_F );
  o_F.close();

  // Compute bounds on residual sum-of-squares by specializing the DAG at 3 different data points
  std::vector< double > xdat = { 3, 2, 1 };
  std::vector< double > ydat = { 9, 5, 2 };
  I IP(-1.,1.), IF( 0. );  
  std::vector<I> IWK;
  double const scaladd = 1.;
  for( size_t k=0; k<xdat.size(); ++k ){
    // Assign constant values to the variables X and Y
    X.set( xdat[k] );
    Y.set( ydat[k] );
    // Evaluate the current residual - the final 'true' argument is to append the result to IF instead of overwriting IF
    DAG.eval( F_op, IWK, 1, &F, &IF, 1, &P, &IP, &scaladd );
    std::cout << "\nSSE bound: " << IF << std::endl;
  }

  return 0;
}

///////////////////////////////////////////////////////////////////////////////

int test_reval1()
{
  std::cout << "\n==============================================\ntest_reval1:\n";

  // Create DAG for: h(z,p)=z2+zp+4
  const unsigned NX = 2, NF = 1;
  mc::FFGraph DAG;
  std::vector<mc::FFVar> X(NX);
  for( auto& Xi : X ) Xi.set( &DAG );
  std::vector<mc::FFVar> F{ sqr(X[0])+X[0]*X[1]+4 };
  auto SGF  = DAG.subgraph( F );
  DAG.output( SGF );

  std::ofstream o_F( "reval1_F.dot", std::ios_base::out );
  DAG.dot_script( F, o_F );
  o_F.close();

  double cputime;
  const unsigned NREP=1000;

  // Evaluate in interval arithmetic, both forward and backward
  std::vector<I> IWKF;
  std::vector<I> IX{ I(-0.8,-0.3), I(6.,9.) }, IF;
  std::cout << "\nInterval hull:\n";
  std::cout << "X[0] = "
            << I(-mc::Op<I>::l(IX[1])/2.+sqrt(mc::sqr(mc::Op<I>::l(IX[1])/2.)-4),
                 -mc::Op<I>::u(IX[1])/2.+sqrt(mc::sqr(mc::Op<I>::u(IX[1])/2.)-4))
            << std::endl;
  std::cout << "X[1] = " << IX[1] << std::endl;  

  cputime = -mc::cpuclock();
  for( unsigned i=0; i<NREP; i++ ){
    DAG.eval( SGF, IWKF, F, IF, X, IX );
  }
  cputime += mc::cpuclock();
  std::cout << "\nDAG interval evaluation w/ forward pass only: "
            << (cputime/=NREP) << " CPU-sec\n";
  for( unsigned i=0; i<NX; i++ ) std::cout << "X[" << i << "] = " << IX[i] << std::endl;
  for( unsigned i=0; i<NF; i++ ) std::cout << "F[" << i << "] = " << IF[i] << std::endl;

  IF[0] = 0.;
  for( unsigned maxpass=1; maxpass<=10; ++maxpass ){
    cputime = -mc::cpuclock();
    int flag = 0;
    for( unsigned i=0; i<NREP; i++ ){
      std::vector<I> IX{ I(-0.8,-0.3), I(6.,9.) };
      flag = DAG.reval( SGF, IWKF, F, IF, X, IX, IINF, maxpass );
    }
    cputime += mc::cpuclock();
    std::cout << "\nDAG interval evaluation w/ " << flag << " forward/backward passes: "
              << (cputime/=NREP) << " CPU-sec\n";
    std::vector<I> IX{ I(-0.8,-0.3), I(6.,9.) };
    DAG.reval( SGF, IWKF, F, IF, X, IX, IINF, maxpass );
    for( unsigned i=0; i<NX; i++ ) std::cout << "X[" << i << "] = " << IX[i] << std::endl;
    for( unsigned i=0; i<NF; i++ ) std::cout << "F[" << i << "] = " << IF[i] << std::endl;
  }

  return 0;
}

///////////////////////////////////////////////////////////////////////////////

int test_reval2()
{
  std::cout << "\n==============================================\ntest_reval2:\n";

  // Example DAG from Schichl & Neumaier (2005)
  const unsigned NX = 3, NF = 3;
  mc::FFGraph DAG;
  mc::FFVar X[NX];
  for( unsigned int i=0; i<NX; i++ ) X[i].set( &DAG );
  mc::FFVar F[NF] = { (4*X[0]-X[1]*X[2])*(X[0]*X[1]+X[2]),
                      sqr(X[0])+sqr(X[1])+X[0]*X[1]+X[1]*X[2]+X[1],
                      exp(X[0]*X[1]+X[1]*X[2]+X[1]+sqrt(X[2])) };
  auto SGF  = DAG.subgraph( NF, F );
  DAG.output( SGF );

  std::ofstream o_F( "reval2_F.dot", std::ios_base::out );
  DAG.dot_script( NF, F, o_F );
  o_F.close();

  // Evaluate in interval arithmetic, both forward and backward
  double INF = 1e20;
  std::vector<I> IWKF;
  try{
    I IX[NX] = { I(-1,1), I(-1,1), I(0,4) },
      IF[NF];
    DAG.eval( SGF, IWKF, NF, F, IF, NX, X, IX );
    std::cout << "\nDAG interval evaluation w/ forward pass only:\n";
    for( unsigned i=0; i<NX; i++ ) std::cout << "X[" << i << "] = " << IX[i] << std::endl;
    for( unsigned i=0; i<NF; i++ ) std::cout << "F[" << i << "] = " << IF[i] << std::endl;
  }
  catch(...){
    std::cout << "\nDAG interval evaluation w/ forward pass only: FAILED\n";
  }

  try{
    I IX[NX] = { I(-1,1), I(-1,1), I(0,4) },
      IF[NF] = { I(-INF,-10), I(-INF,INF), I(0,2) };
    int flag = DAG.reval( SGF, IWKF, NF, F, IF, NX, X, IX, IINF );
    std::cout << "\nDAG interval evaluation w/ forward/backward passes:\n";
    for( unsigned i=0; i<NX; i++ ) std::cout << "X[" << i << "] = " << IX[i] << std::endl;
    for( unsigned i=0; i<NF; i++ ) std::cout << "F[" << i << "] = " << IF[i] << std::endl;
    std::cout << "FLAG = " << flag << std::endl;
  }
  catch(...){
    std::cout << "\nDAG interval evaluation w/ forward/backward passes: FAILED\n";
  }

  return 0;
}

///////////////////////////////////////////////////////////////////////////////

int test_reval3()
{
  std::cout << "\n==============================================\ntest_reval3:\n";

  // Example constraint propagation from Wikipedia (https://en.wikipedia.org/wiki/Interval_propagation)
  const unsigned NX = 7, NF = 4;
  mc::FFGraph DAG;
  mc::FFVar X[NX], &E=X[0], &J=X[1], &U1=X[2], &U2=X[3], &P=X[4], &R1=X[5], &R2=X[6];
  for( unsigned int i=0; i<NX; i++ ) X[i].set( &DAG );
  mc::FFVar F[NF] = { P - E * J,
                      U1 - R1 * J,
                      U2 - R2 * J,
                      E - U1 - U2 };
  auto SGF  = DAG.subgraph( NF, F );
  DAG.output( SGF );

  std::ofstream o_F( "reval3_F.dot", std::ios_base::out );
  DAG.dot_script( NF, F, o_F );
  o_F.close();

  // Evaluate in interval arithmetic, both forward and backward
  double INF = 1e20;
  std::vector<I> IWKF;
  try{
    I IX[NX] = { I(23,26), I(4,8), I(10,11), I(14,17), I(124,130), I(0,INF), I(0,INF) },
      IF[NF];
    DAG.eval( SGF, IWKF, NF, F, IF, NX, X, IX );
    std::cout << "\nDAG interval evaluation w/ forward pass only:\n";
    for( unsigned i=0; i<NX; i++ ) std::cout << "X[" << i << "] = " << IX[i] << std::endl;
    for( unsigned i=0; i<NF; i++ ) std::cout << "F[" << i << "] = " << IF[i] << std::endl;
  }
  catch(...){
    std::cout << "\nDAG interval evaluation w/ forward pass only: FAILED\n";
  }

  try{
    I IX[NX] = { I(23,26), I(4,8), I(10,11), I(14,17), I(124,130), I(0,INF), I(0,INF) },
      IF[NF] = { I(0), I(0), I(0), I(0) };
    int flag = DAG.reval( SGF, IWKF, NF, F, IF, NX, X, IX, IINF );
    std::cout << "\nDAG interval evaluation w/ forward/backward passes:\n";
    for( unsigned i=0; i<NX; i++ ) std::cout << "X[" << i << "] = " << IX[i] << std::endl;
    for( unsigned i=0; i<NF; i++ ) std::cout << "F[" << i << "] = " << IF[i] << std::endl;
    std::cout << "FLAG = " << flag << std::endl;
  }
  catch(...){
    std::cout << "\nDAG interval evaluation w/ forward/backward passes: FAILED\n";
  }

  return 0;
}

///////////////////////////////////////////////////////////////////////////////

int test_reval4()
{
  std::cout << "\n==============================================\ntest_reval4:\n";

  //No. stages, feed location, No. components 
  const unsigned nc = 3, ns = 3, pf = 1;
  //const unsigned nc = 3, ns = 5, pf = 3;
  //const unsigned nc = 3, ns = 7, pf = 4;
  //const unsigned nc = 3, ns = 10, pf = 5;
  //const unsigned nc = 3, ns = 20, pf = 10;

  mc::FFGraph DAG;
  const unsigned NF = ns*(nc*2+1), NP = NF+1;
  mc::FFVar P[NP], F[NF];
  for( unsigned int i=0; i<NP; i++ ) P[i].set( &DAG );

  /*
    ################### DATA SECTION ###############################################

  */
  // Model Parameters
  double Patm = 1e5;				//Atmospheric pressure
  double V    = 1.38;				//Vapour flow rate
  /*
  	let a[1] := 23.4832; 	let a[2] := 20.5110; 	let a[3] := 20.9064;
	let b[1] := -3634.01;	let b[2] := -2664.30;	let b[3] := -3096.52;
	let c[1] := -33.768;	let c[2] := -79.483;	let c[3] := -53.668;
	
	let r[1,1] := 0.0;	let r[1,2] :=  0.7411;	let r[1,3] :=  0.9645;	
	let r[2,1] := -1.0250;	let r[2,2] := 0.0;	let r[2,3] := -1.4350;
	let r[3,1] := -0.9645;	let r[3,2] :=  2.7470;	let r[3,3] := 0.0;

	let s[1,1] := 0.0;	let s[1,2] := -477.00;	let s[1,3] := -903.1024;
	let s[2,1] :=  72.78;	let s[2,2] := 0.0;	let s[2,3] :=  768.20;
	let s[3,1] := -140.9995;let s[3,2] := -1419.0;	let s[3,3] := 0.0;
  */
  double a[nc] = { 23.4832,  20.5110, 20.9064 }; 
  double b[nc] = { -3634.01, -2664.30, -3096.52}; 
  double c[nc] = { -33.768, -79.483, -53.668}; 
  double r[nc*nc] = { 0.0,  0.7411, 0.9645, -1.0250, 0.0, -1.4350, -0.9645, 2.7470, 0.0};
  double s[nc*nc] = { 0.0, -477.00, -903.1024, 72.78, 0.0, 768.20, -140.9995, -1419.0, 0.0};
  /*
    for {i in 1..C, j in 0..N}
	let f[i,j] := 0.0;

	let f[1, N_F] := 0.4098370;
	let f[2, N_F] := 0.01229769;
	let f[3, N_F] := 0.06090665;
  */
  double f[nc*(ns+1)];   
  for( unsigned i = 0; i < nc; ++i ){
  	for( unsigned j = 0; j < ns+1; ++j ) f[i*ns+j] = 0.0;
  }
  f[ 0*ns+pf ] = 0.4098370;
  f[ 1*ns+pf ] = 0.01229769;
  f[ 2*ns+pf ] = 0.06090665;
  /*
    for {j in 0..N} let F[j] := sum{i in 1..C} f[i,j];
 */
  double FF[ns+1], sumd;
  for( unsigned j = 0; j < ns+1; ++j ){
  	sumd = 0.0;
  	for( unsigned i = 0; i < nc; ++i ) sumd += f[ i*ns + j ];
  	FF[j] = sumd;
  } 

  /*
    ############### VARIABLES ######################################################

  */
  mc::FFVar  D  = *P;
  /*
 	var x{i in 1..C, j in 1..N} >= x_L[i,j], <= x_U[i,j], := x_0[i,j];
    var K{i in 1..C, j in 1..N} >= K_L[i,j], <= K_U[i,j], := K_0[i,j];
    var T{j in 1..N} >= 336.3, <= 383.4, := T_0[j];
  */
  mc::FFVar *X  = P+1;
  mc::FFVar *K  = X+nc*ns;
  mc::FFVar *Te = K+nc*ns;

  /*
    ####### DEFINED VARIABLES ######################################################

  */
  mc::FFVar B, L[ns+1], p[nc*ns], Lambda[nc*nc*ns], sum_xLambda[nc*ns], G[nc*ns], sum;
  /*
  	param B := F[N_F] - D;
  */
  B = FF[pf] - D;
  /*
  	param L{j in 0..N} = V - D + sum{k in 0..j} F[k];
  */
  for( unsigned j = 0; j < ns+1 ; ++j ){
    sumd = 0.0;
    for( unsigned k = 0; k <= j; ++k )
      sumd += FF[k];
    L[j] = V - D + sumd ;	
  }
  /*
  	var p{i in 1..C, j in 1..N} = exp(a[i]+b[i]/(T[j]+c[i]));
  */
  for( unsigned i = 0; i < nc; ++i  )
    for( unsigned j = 0; j < ns; ++j )
      p[i*ns+j] = exp( a[i] + b[i] / (Te[j]+c[i]) );  	
  /*
    var rcp_T{j in 1..N} = 1.0/T[j];
  	var Lambda{i1 in 1..C, i2 in 1..C, j in 1..N} = exp(r[i1,i2]+s[i1,i2]*rcp_T[j]);
  */
  for( unsigned i1 = 0; i1 < nc; ++i1 )
    for( unsigned i2 = 0; i2 < nc; ++i2 )
      for( unsigned j = 0; j < ns; ++j )
        Lambda[nc*ns*i1+ns*i2+j] = exp(r[nc*i1+i2] + s[nc*i1+i2]/Te[j]);
  /*
  	var sum_xLambda{i in 1..C, j in 1..N} = sum{i1 in 1..C} (x[i1,j]*Lambda[i,i1,j]);
  */
  for( unsigned i = 0; i < nc; ++i  )
    for( unsigned j = 0; j < ns; ++j ){
      sum = 0.;
      for( unsigned i1 = 0; i1 < nc; ++i1 )
        sum += X[i1*ns+j] * Lambda[i*nc*ns+i1*ns+j];	
      sum_xLambda[i*ns+j] = sum;
    }
  /*
    var rcp_sum_xLambda{i in 1..C, j in 1..N} = 1.0/sum_xLambda[i,j];
    var gamma{i in 1..C, j in 1..N} = exp( -log(sum_xLambda[i,j]) + 1.0
                - (sum{i2 in 1..C} (x[i2,j]*Lambda[i2,i,j]*rcp_sum_xLambda[i2,j])) );
  */
  for( unsigned i = 0; i < nc; ++i  )
    for( unsigned j = 0; j < ns; ++j ){
      sum = 0.;
      for( unsigned i2 = 0; i2 < nc; ++i2 )
        sum += X[i2*ns+j] * Lambda[nc*ns*i2+ns*i+j] / sum_xLambda[i2*ns+j];
      G[i*ns+j] = exp( 1.0 - sum - log( sum_xLambda[i*ns+j] ) );
    }

  /*
    ############## EQUATIONS #######################################################

  */
  unsigned ieq=0;
  /*
  	AUXILIARY EQUATIONS
  	E_aux_K{j in 1..N, i in 1..C}: 	K[i,j] - gamma[i,j]*(p[i,j]/P) = 0.0;  
  */
  for( unsigned j = 0; j < ns; ++j  ){
    for( unsigned i = 0; i < nc; ++i ){
      F[ ieq++ ] = K[i*ns+j] - G[i*ns+j] * ( p[i*ns+j] / Patm );
    }
  }
  /*
  	MATERIAL BALANCES
  	M_tot{i in 1..C}: D*(K[i,1]*x[i,1]) + B*x[i,N] - f[i,N_F] = 0.0;
  */
  for( unsigned i = 0; i < nc; ++i  ){
    F[ ieq++ ] = D*(K[i*ns+0]*X[i*ns+0]) + B*X[i*ns+(ns-1)] - f[i*ns+pf];
  } 
  /* 
  	NOTE THE UNUSUAL FORMULATION
	M_eq{j in 1..N-1, i in 1..C}:
	L[j]*x[i,j] + sum{i1 in j+1..N} f[i,i1] - B*x[i,N] - V*(K[i,j+1]*x[i,j+1]) = 0.0;
  */
  for( unsigned j = 0; j < ns-1; ++j  ){
    for( unsigned i = 0; i < nc; ++i ){
      sum = 0.; 
      for( unsigned i2 = j+1; i2 < ns; ++i2 )
        sum += f[ i*ns+i2 ];
      F[ ieq++ ] = L[j]*X[i*ns+j] + sum - B*X[i*ns+(ns-1)] - V*( K[i*ns+j+1]*X[i*ns+j+1] );
    } 
  } 
  /*
  	SUMMATION EQUATIONS
    S_x_eq{j in 1..N}:
    sum{i in 1..C} x[i,j] - 1.0 = 0.0;
  */
  for( unsigned j = 0; j < ns; ++j  ){
    sum = -1.;
    for( unsigned i = 0; i < nc; ++i )
      sum += X[i*ns+j];
    F[ ieq++ ] = sum;
  }

  std::vector<I> IWKF;
  auto SGF  = DAG.subgraph( NF, F );

  /*
    ############## VARIABLE BOUNDS #################################################

  */
  // Variable bounds
//  I Ip[NP], &ID = *Ip, *IX = Ip+1, *IK = IX+nc*ns, *ITe = IK+nc*ns;
//  //ID = I( 0.455, 0.455 );
//  //ID = I( 0.45, 0.46 );
//  ID = I( 0.44, 0.47 );
//  for( unsigned j=0; j<ns; ++j ){
//    IX[j]  = I( 0.0001, 0.9999 ); IX[ns+j] = I( 0.0001, 0.9999 ); IX[2*ns+j] = I( 0.0001, 0.9999 );
//    IK[j]  = I( 0.97 , 40.52  );  IK[ns+j] = I( 0.2445, 1.317  ); IK[2*ns+j] = I( 0.2745, 1.975  );
//    ITe[j] = I( 350., 351.  );//I( 336.3, 383.4  );
//  }

  const I Ip_3stg[NP] = { 
    I(4.51882712332260e-01,4.53400900000000e-01),
    I(8.51605100000000e-01,8.52364800000000e-01),
    I(7.98175092047527e-01,8.00297694717604e-01),
    I(4.41986225216501e-01,4.59569700000000e-01),
    I(2.95642700000000e-02,2.99233900000000e-02),
    I(5.11901000000000e-02,5.19930600000000e-02),
    I(1.57389600000000e-01,1.62274478001944e-01),
    I(1.18065000000000e-01,1.18472600000000e-01),
    I(1.48509700000000e-01,1.49843600000000e-01),
    I(3.82601600000000e-01,3.89211950000000e-01),
    I(1.02686900000000e+00,1.02749100000000e+00),
    I(1.07442300000000e+00,1.07664100000000e+00),
    I(1.75775700000000e+00,1.82370300000000e+00),
    I(5.51961100000000e-01,5.52927300000000e-01),
    I(4.90793000000000e-01,4.93025500000000e-01),
    I(3.05665600000000e-01,3.09595538816592e-01),
    I(9.15561400000000e-01,9.17969100000000e-01),
    I(7.68530500000000e-01,7.73717500000000e-01),
    I(3.65186900000000e-01,3.74117500000000e-01),
    I(3.36703700000000e+02,3.36707700000000e+02),
    I(3.36960000000000e+02,3.36970800000000e+02),
    I(3.38874100000000e+02,3.39011800000000e+02)
  };

  /*
    ############## CONSTRAINT PROPAGATION ##########################################

  */
  // Compute forward bounds on equality constraints
  std::cout << "\nDAG interval evaluation using forward pass only:\n";
  std::cout << "Ip =" << std::endl;
  std::vector<I> vIp( Ip_3stg, Ip_3stg+NP ), vIf( NF, I(0.) );
  for( unsigned i=0; i<vIp.size(); i++ )
    std::cout << "  " << P[i] << " = " << vIp[i] << std::endl;
  DAG.eval( SGF, IWKF, NF, F, vIf.data(), NP, P, vIp.data() );
  std::cout << "If =" << std::endl;
  for( unsigned i=0; i<vIf.size(); i++ )
    std::cout << "  " << F[i] << " = " << vIf[i] << std::endl;

  // Compute forward/backward bounds on equality constraints
  vIp.assign( Ip_3stg, Ip_3stg+NP ); vIf.assign( NF, I(0.) );
  std::cout << "\nDAG interval evaluation using forward/backward passes:";
  int flag = 0;
  try{ flag = DAG.reval( SGF, IWKF, NF, F, vIf.data(), NP, P, vIp.data(), IINF ); }
  catch(...){ std::cout << "Failure\n"; }
  if( flag < 0 ) std::cout << " infeasible\n";
  else std::cout << " " << flag << " passes\n";
  for( unsigned i=0; i<vIp.size(); i++ )
    std::cout << "  " << P[i] << " = " << vIp[i] << std::endl;
  std::cout << "If =" << std::endl;
  for( unsigned i=0; i<vIf.size(); i++ )
    std::cout << "  " << F[i] << " = " << vIf[i] << std::endl;

  return 0;
}

///////////////////////////////////////////////////////////////////////////////

int test_veval1()
{
  std::cout << "\n==============================================\ntest_veval1:\n";

  // Create DAG
  const unsigned NX = 4, NF = 2;
  mc::FFGraph DAG;
  std::vector<mc::FFVar> X(NX);
  for( auto& Xi : X ) Xi.set( &DAG );
  std::vector<mc::FFVar> F{
    X[2]*X[3]-2./(X[1]+X[2]),
    X[0]/pow(exp(X[2]*X[1])+3.,3)+tanh(X[3])
  };
  auto&& F_op  = DAG.subgraph( F );
  std::cout << DAG;

  const unsigned NREP=1000000;//000;

  // Evaluate with doubles, no parameter pack
  std::vector<double> dX{ -1., -1., 2., 3. };
  for( unsigned i=0; i<NX; i++ ) std::cout << "X[" << i << "] = " << dX[i] << std::endl;
  std::vector<std::vector<double>> v_dX( NREP, dX );
  static double const one = 2.;
    
  for( DAG.options.MAXTHREAD = 1; DAG.options.MAXTHREAD <= 12; ++DAG.options.MAXTHREAD ){
    std::vector<std::vector<double>> v_dF;
    std::vector<double> dwk;
    std::vector<mc::FFGraph::Worker<double>> thwk;
    auto starttime = std::chrono::system_clock::now();
    DAG.veval( F_op, dwk, thwk, F, v_dF, X, v_dX );
    auto walltime = std::chrono::duration_cast<std::chrono::microseconds>( std::chrono::system_clock::now() - starttime );
    std::cout << "\nvectorized DAG evaluation on " << DAG.options.MAXTHREAD << " threads: " << walltime.count()*1e-6 << " sec\n";
    starttime = std::chrono::system_clock::now();
    DAG.veval( F_op, dwk, thwk, F, v_dF, X, v_dX );
    walltime = std::chrono::duration_cast<std::chrono::microseconds>( std::chrono::system_clock::now() - starttime );
    std::cout << "\nvectorized DAG evaluation on " << DAG.options.MAXTHREAD << " threads: " << walltime.count()*1e-6 << " sec\n";
    starttime = std::chrono::system_clock::now();
    std::cout << "v_dF: " << v_dF.front()[0] << "  " << v_dF.front()[1] << std::endl;
    std::cout << "v_dF: " << v_dF.back()[0] << "  " << v_dF.back()[1] << std::endl;
    DAG.veval( F_op, dwk, thwk, F, v_dF, X, v_dX, &one );
    walltime = std::chrono::duration_cast<std::chrono::microseconds>( std::chrono::system_clock::now() - starttime );
    std::cout << "\nvectorized DAG evaluation on " << DAG.options.MAXTHREAD << " threads: " << walltime.count()*1e-6 << " sec\n";
    std::cout << "v_dF: " << v_dF.front()[0] << "  " << v_dF.front()[1] << std::endl;
    std::cout << "v_dF: " << v_dF.back()[0] << "  " << v_dF.back()[1] << std::endl;
  }
  
  //for( auto dF : v_dF ){
  //  for( auto dFi : dF ) std::cout << " " << dFi;
  //  std::cout << std::endl;
  //} 
  //for( unsigned i=0; i<NF; i++ ) std::cout << "F[" << i << "] = " << dF[i] << std::endl;

  return 0;
}

///////////////////////////////////////////////////////////////////////////////

int main()
{
  try{
//    test_eval1();
//    test_eval2();
//    test_eval3();
//    test_reval1();
//    test_reval2();
//    test_reval3();
//    test_reval4();
    test_veval1();
  }
  catch( mc::FFBase::Exceptions &eObj ){
    std::cerr << "Error " << eObj.ierr()
              << " in factorable function manipulation:" << std::endl
              << eObj.what() << std::endl
              << "Aborts." << std::endl;
    return eObj.ierr();
  }
#if !defined(MC__USE_PROFIL) && !defined(MC__USE_FILIB) && !defined(MC__USE_BOOST)
  catch( I::Exceptions &eObj ){
    std::cerr << "Error " << eObj.ierr()
         << " in natural interval extension:" <<    std::endl
	 << eObj.what() << std::endl
         << "Aborts." << std::endl;
    return eObj.ierr();
  }
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

