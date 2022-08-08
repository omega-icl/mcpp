#include <fstream>
#include <sstream>
#include <iomanip>

#include "ffunc.hpp"
#include "rltred.hpp"
#include "sparseexpr.hpp"

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
  
  std::cout << "Variable dependence of F[0]: " << F[0] << std::endl;
  std::cout << "Variable dependence of F[1]: " << F[1] << std::endl;

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

  std::cout << "Variable dependence of F[0]: " << depF[0] << std::endl;
  //std::cout << "Variable dependence of F[1]: " << depF[1] << std::endl;

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

  mc::RLTRed RRLT( &DAG );
  RRLT.options.METHOD        = mc::RLTRed::Options::ILP;
  RRLT.options.LEVEL         = mc::RLTRed::Options::FULLSIM;
  RRLT.options.NODIV         = false;
  RRLT.options.DISPLEVEL     = 1;
  RRLT.options.MIPDISPLEVEL  = 0;
  RRLT.options.MIPOUTPUTFILE = "rltred1.lp";

  RRLT.search( NF, F );

  std::cout << DAG;
  auto FRED = RRLT.constraints();
  for( auto it=FRED.begin(); it!=FRED.end(); ++it ){
    std::ostringstream ostr; ostr << " OF REDUCTION CONSTRAINT " << **it;
    DAG.output( DAG.subgraph( 1, *it ), ostr.str() );
  }

  RRLT.stats.display();

  F[2] = F1 + F2 - FT;
  std::cout << DAG;

  RRLT.search( NF, F );

  std::cout << DAG;
  auto FRED2 = RRLT.constraints();
  for( auto it=FRED2.begin(); it!=FRED2.end(); ++it ){
    std::ostringstream ostr; ostr << " OF REDUCTION CONSTRAINT " << **it;
    DAG.output( DAG.subgraph( 1, *it ), ostr.str() );
  }

  RRLT.stats.display();

  return 0;
}

///////////////////////////////////////////////////////////////////////////////

int test_rltred2()
{
  std::cout << "\n==============================================\ntest_rltred2:\n";

  mc::FFGraph DAG;
  const unsigned NX = 3, NF = 2;
  mc::FFVar X[NX];
  for( unsigned i(0); i<NX; i++ ) X[i].set( &DAG );
  mc::FFVar F[NF];
  // X0^2 + X0·X1 - X0 = 0 && X0·X1 + X1^2 - X1 = 0 => X0^2 - X1^2 - X0 + X1 = 0
  F[0] = X[0] + X[1] - 1.;
  F[1] = sqr(X[0]) - sqr(X[1]) - X[2];
  std::cout << DAG;

  mc::RLTRed RRLT( &DAG );
  RRLT.options.METHOD        = mc::RLTRed::Options::ILP;
  RRLT.options.LEVEL         = mc::RLTRed::Options::FULLSIM;
  RRLT.options.NODIV         = false;
  RRLT.options.DISPLEVEL     = 2;
  RRLT.options.MIPDISPLEVEL  = 1;
  RRLT.options.MIPOUTPUTFILE = "rltred2.lp";

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
       - 2*(X[1]*X[5]) + X[2]*X[3] - 4*(X[2]*X[4]) - 3*(X[2]*X[5]) - 6*(X[3]*X[4])
       - 9*(X[3]*X[5]) - X[4]*X[5];
  std::cout << DAG;

  mc::RLTRed RRLT( &DAG );
  RRLT.options.METHOD        = mc::RLTRed::Options::ILP;
  RRLT.options.LEVEL         = mc::RLTRed::Options::FULLSIM;
  RRLT.options.NODIV         = false;
  RRLT.options.DISPLEVEL     = 2;
  RRLT.options.MIPDISPLEVEL  = 1;
  RRLT.options.MIPOUTPUTFILE = "rltred3.lp";

  RRLT.search( NF, F.data() );

  auto FRED = RRLT.constraints();
  for( auto it=FRED.begin(); it!=FRED.end(); ++it ){
    std::ostringstream ostr; ostr << " OF REDUCTION CONSTRAINT " << **it;
    DAG.output( DAG.subgraph( 1, *it ), ostr.str() );
    //break;
  }

  RRLT.stats.display();

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

int test_rltred4()
{
  std::cout << "\n==============================================\ntest_rltred4:\n";

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

  /*
    ############## REDUCTION CONSTRAINTS ###########################################
  */
  mc::RLTRed RRLT( &DAG );
  RRLT.options.METHOD        = mc::RLTRed::Options::ILP;
  RRLT.options.LEVEL         = mc::RLTRed::Options::FULLSIM;
  RRLT.options.NODIV         = false;
  RRLT.options.DISPLEVEL     = 1;
  RRLT.options.MIPDISPLEVEL  = 1;
  RRLT.options.MIPOUTPUTFILE = "rltred4.lp";

  RRLT.search( NF, F );

  auto FRED = RRLT.constraints();
  for( auto it=FRED.begin(); it!=FRED.end(); ++it ){
    std::ostringstream ostr; ostr << " OF REDUCTION CONSTRAINT " << **it;
    DAG.output( DAG.subgraph( 1, *it ), ostr.str() );
  }

  RRLT.stats.display();

//  RRLT.search( SPE.Expr().size(), SPE.Expr().data() );

//  auto FRED2 = RRLT.constraints();
//  for( auto it=FRED2.begin(); it!=FRED2.end(); ++it ){
//    std::ostringstream ostr; ostr << " OF REDUCTION CONSTRAINT " << **it;
//    DAG.output( DAG.subgraph( 1, *it ), ostr.str() );
//  }

//  RRLT.stats.display();

  return 0;
}

///////////////////////////////////////////////////////////////////////////////

template <class U>
U ACTMOD
( const U&x, const unsigned n, const double*a )
{
  if( !n ) return 0.;
  U act = a[n-1];
  for( unsigned i=n-1; i>0; i-- ){ act *= x; act += a[i-1]; }
  return act;
}

int test_rltred5()
{
  std::cout << "\n==============================================\ntest_rltred5:\n";

  mc::FFGraph DAG;
  const unsigned NCOMP = 2, NSTG = 1, NOUT = 2;
  const unsigned NP = NSTG*(5+5*NCOMP) + (NSTG+NOUT)*(1+2*NSTG) + NOUT*(1+NCOMP) + 1;
  const unsigned NAUX = (1+2*NSTG)*(NSTG+NOUT)*NCOMP;
  mc::FFVar p[NP+NAUX];
  for( unsigned i=0; i<NP; i++ ) p[i].set( &DAG );
  std::vector<mc::FFVar> EQ;
  /*
    ################### DATA SECTION ##############################################
  */
  enum COMP { C7=0, C16 };
  const double MVOL[NCOMP]  = { 1.482e-4,  2.997e-4   }, // [m3/mol]
               PERM[NCOMP]  = { 3.5177e0,  0.54340e0 },  // [kg/m2·bar·h]
               CFEED[NCOMP] = { 0.778,     0.222     };  // [wt%]
  const double ACT[NCOMP*4] = { 0.7906, 0.5099, -0.3958, 0.0953,
                                0.6757, 0.3681, 0.2066, -0.2504 };
  const double R      = 8.314e-5,  // [m3·bar/K·mol]
               T      = 298e0,     // [K]
               FFEED  = 1.0124,    // [kg/day]
               ACELL  = 14e-4;     // [m2]
  /*
    ############### VARIABLES ######################################################
  */
  mc::FFVar *FFMIX = p, *FRMIX = FFMIX+NSTG+NOUT, *FPMIX = FRMIX+NSTG*(NSTG+NOUT),
            *FFMEM = FPMIX+NSTG*(NSTG+NOUT), *FRMEM = FFMEM+NSTG, *FPMEM = FRMEM+NSTG, *SCMEM = FPMEM+NSTG,
            *CFMEM = SCMEM+NSTG, *CRMEM = CFMEM+NSTG*NCOMP, *CPMEM = CRMEM+NSTG*NCOMP,
            *JCMEM = CPMEM+NSTG*NCOMP, *ARATIO = JCMEM+NSTG*NCOMP, *DPMEM = ARATIO+NSTG*NCOMP, 
            *FOUT = DPMEM+NSTG, *COUT = FOUT+NOUT, *SCOUT = COUT+NOUT*NCOMP,
            *CFMIX = SCOUT+1, *CRMIX = CFMIX+(NSTG+NOUT)*NCOMP, *CPMIX = CRMIX+NSTG*(NSTG+NOUT)*NCOMP;

  /*
    ############## EQUATIONS #######################################################
  */
  // Feed splitter
  mc::FFVar BALFFEED = FFEED;
  for( unsigned stg=0; stg<NSTG+NOUT; ++stg )
    BALFFEED -= FFMIX[stg];
  EQ.push_back( BALFFEED );
  for( unsigned comp=0; comp<NCOMP; ++comp ){
    for( unsigned stg=0; stg<NSTG+NOUT; ++stg )
      //EQ.push_back( CFMIX[stg*NCOMP+comp] - CFEED[comp] );
      CFMIX[stg*NCOMP+comp] = CFEED[comp];
  }

  // Membrane retentate and permeate splitters
  for( unsigned stg=0; stg<NSTG; ++stg ){
    mc::FFVar BALRSPLIT = FRMEM[stg], BALPSPLIT = FPMEM[stg];
    for( unsigned mix=0; mix<NSTG+NOUT; ++mix ){
      BALRSPLIT -= FRMIX[stg*(NSTG+NOUT)+mix];
      BALPSPLIT -= FPMIX[stg*(NSTG+NOUT)+mix];
    }
    EQ.push_back( BALRSPLIT );
    EQ.push_back( BALPSPLIT );
    for( unsigned mix=0; mix<NSTG+NOUT; ++mix )
      for( unsigned comp=0; comp<NCOMP; ++comp ){
        //EQ.push_back( CRMIX[(stg*(NSTG+NOUT)+mix)*NCOMP+comp] - CRMEM[stg*NCOMP+comp] );
        //EQ.push_back( CPMIX[(stg*(NSTG+NOUT)+mix)*NCOMP+comp] - CPMEM[stg*NCOMP+comp] );
        CRMIX[(stg*(NSTG+NOUT)+mix)*NCOMP+comp] = CRMEM[stg*NCOMP+comp];
        CPMIX[(stg*(NSTG+NOUT)+mix)*NCOMP+comp] = CPMEM[stg*NCOMP+comp];
      }
  }

  // Membrane modules
  for( unsigned stg=0; stg<NSTG; ++stg ){
    mc::FFVar BALCFMEM = -1., BALCPMEM = -1., BALCRMEM = -1., JTMEM = 0.;
    for( unsigned comp=0; comp<NCOMP; ++comp ){
      BALCFMEM += CFMEM[stg*NCOMP+comp];
      BALCPMEM += CPMEM[stg*NCOMP+comp];
      BALCRMEM += CRMEM[stg*NCOMP+comp];
      JTMEM    += JCMEM[stg*NCOMP+comp];
    }
    EQ.push_back( FFMEM[stg] - FPMEM[stg] - FRMEM[stg] );
    EQ.push_back( SCMEM[stg] * FFMEM[stg] - FPMEM[stg] );
    EQ.push_back( FPMEM[stg] - ACELL * JTMEM * 24. * DPMEM[stg] );
    for( unsigned comp=0; comp<NCOMP; ++comp ){
      if( comp )
        EQ.push_back( CFMEM[stg*NCOMP+comp] * FFMEM[stg] - CPMEM[stg*NCOMP+comp] * FPMEM[stg] - CRMEM[stg*NCOMP+comp] * FRMEM[stg] );
      else
        EQ.push_back( BALCRMEM );
      EQ.push_back( ARATIO[stg*NCOMP+comp] * ACTMOD( CRMEM[stg*NCOMP+comp], 4, ACT+4*comp ) - ACTMOD( CPMEM[stg*NCOMP+comp], 4, ACT+4*comp ) );
      EQ.push_back( JCMEM[stg*NCOMP+comp] - PERM[comp] * ( CRMEM[stg*NCOMP+comp] - ARATIO[stg*NCOMP+comp] * CPMEM[stg*NCOMP+comp] * exp( - MVOL[comp] * DPMEM[stg] / R / T ) ) );
      EQ.push_back( JCMEM[stg*NCOMP+comp] - CPMEM[stg*NCOMP+comp] * JTMEM );
    }
  }

  // Membrane retentate and output mixers
  for( unsigned mix=0; mix<NSTG+NOUT; ++mix ){
    mc::FFVar BALCOMIX = -1.;
    for( unsigned comp=0; comp<NCOMP; ++comp )
      BALCOMIX += (mix<NSTG? CFMEM[mix*NCOMP+comp]: COUT[(mix-NSTG)*NCOMP+comp]);
    mc::FFVar BALFMIX = (mix<NSTG? FFMEM[mix]: FOUT[mix-NSTG]) - FFMIX[mix];
    for( unsigned stg=0; stg<NSTG; ++stg )
      BALFMIX -= FPMIX[stg*(NSTG+NOUT)+mix] + FRMIX[stg*(NSTG+NOUT)+mix];
    EQ.push_back( BALFMIX );
    for( unsigned comp=0; comp<NCOMP; ++comp ){
      if( comp ){
        mc::FFVar BALCMIX = (mix<NSTG? FFMEM[mix] * CFMEM[mix*NCOMP+comp]: FOUT[mix-NSTG] * COUT[(mix-NSTG)*NCOMP+comp]) - FFMIX[mix] * CFMIX[mix*NCOMP+comp];
        for( unsigned stg=0; stg<NSTG; ++stg ){
          BALCMIX -= FPMIX[stg*(NSTG+NOUT)+mix] * CPMIX[(stg*(NSTG+NOUT)+mix)*NCOMP+comp] + FRMIX[stg*(NSTG+NOUT)+mix] * CRMIX[(stg*(NSTG+NOUT)+mix)*NCOMP+comp];
        }
        EQ.push_back( BALCMIX );
      }
      else
        EQ.push_back( BALCOMIX );
    }
  }

  // Overall balances
  EQ.push_back( *SCOUT * FFEED - FOUT[0] );

  /*
    ############## REDUCTION CONSTRAINTS ###########################################
  */
  mc::RLTRed RRLT( &DAG );
  RRLT.options.METHOD        = mc::RLTRed::Options::ILP;
  RRLT.options.LEVEL         = mc::RLTRed::Options::PRIMSIM;
  RRLT.options.NODIV         = false;
  RRLT.options.DISPLEVEL     = 1;
  RRLT.options.MIPDISPLEVEL  = 1;
  RRLT.options.MIPOUTPUTFILE = "rltred5.lp";

  RRLT.search( EQ.size(), EQ.data() );

  auto FRED = RRLT.constraints();
  for( auto it=FRED.begin(); it!=FRED.end(); ++it ){
    std::ostringstream ostr; ostr << " OF REDUCTION CONSTRAINT " << **it;
    DAG.output( DAG.subgraph( 1, *it ), ostr.str() );
  }

  RRLT.stats.display();

  return 0;
}

///////////////////////////////////////////////////////////////////////////////

int test_rltred6()
{
  std::cout << "\n==============================================\ntest_rltred6:\n";

  mc::FFGraph DAG;
  const unsigned NX = 10, NF = 2;
  std::vector<mc::FFVar> X(NX);
  for( unsigned i(0); i<NX; i++ ) X[i].set( &DAG );
  std::vector<mc::FFVar> F(NF);
  F[0] = 1;
  for( unsigned i(0); i<NX; i++ ) F[0] -= X[i];
  F[1] = X[0]*X[8] + X[0]*X[9] + X[1]*X[9] + X[0]*X[4] + X[3]*X[6];
  for( unsigned i(0); i<NX-1; i++ ) F[1] += X[i]*X[i+1];
  for( unsigned i(0); i<NX-2; i++ ) F[1] += X[i]*X[i+2];
  std::cout << DAG;

  mc::RLTRed RRLT( &DAG );
  RRLT.options.METHOD        = mc::RLTRed::Options::ILP;
  RRLT.options.LEVEL         = mc::RLTRed::Options::FULLSIM;
  RRLT.options.NODIV         = false;
  RRLT.options.DISPLEVEL     = 1;
  RRLT.options.MIPDISPLEVEL  = 1;
  RRLT.options.MIPOUTPUTFILE = "rltred6.lp";

  RRLT.search( NF, F.data() );

  auto FRED = RRLT.constraints();
  for( auto it=FRED.begin(); it!=FRED.end(); ++it ){
    std::ostringstream ostr; ostr << " OF REDUCTION CONSTRAINT " << **it;
    DAG.output( DAG.subgraph( 1, *it ), ostr.str() );
    //break;
  }

  RRLT.stats.display();

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

int test_sexpr1()
{
  std::cout << "\n==============================================\ntest_sexpr1:\n";

  mc::FFGraph DAG;
  const unsigned NX = 4, NF = 1;
  mc::FFVar X[NX];
  for( unsigned i(0); i<NX; i++ ) X[i].set( &DAG );
  mc::FFVar F[NF];
  //F[0] = pow( X[0] + X[1], 3 );
  F[0] = pow( X[0] + 1 / sqr( X[1] ), 3 );
  //F[0] = 4*sqr(X[0]) - 2.1*pow(X[0],4) + pow(X[0],6)/3 + X[0]*X[1] - 4*sqr(X[1]) + 4*pow(X[1],4);
  //F[0] = 250*exp(X[2])*X[0] + 250*pow(X[3],0.6)*X[1];
  //F[1] = exp( 2 * sqr( X[1] ) - 1 );
  std::cout << DAG;

  mc::SparseExpr<mc::FFGraph<>>::SPolyExpr::options.BASIS = mc::SparseExpr<mc::FFGraph<>>::SPolyExpr::Options::MONOM;
  mc::SparseEnv<mc::FFGraph<>> SPE( &DAG );
  SPE.options.LIFTDIV = true;//false;//
  SPE.options.LIFTIPOW = true; //false;//

  SPE.process( NF, F );

  std::cout << std::endl << SPE.Var().size() << " participating variables: ";
  for( auto&& var : SPE.Var() ) std::cout << var << " ";
  std::cout << std::endl;
  std::cout << std::endl << SPE.Aux().size() << " auxiliary variables: ";
  for( auto&& aux : SPE.Aux() ) std::cout << *aux.first << "->" << *aux.second << " ";
  std::cout << std::endl;
  std::cout << std::endl << SPE.Poly().size() << " polynomial constraints: " << std::endl;
  for( auto&& expr : SPE.Poly() ) DAG.output( DAG.subgraph( 1, &expr ) );
  //std::cout << std::endl;
  std::cout << std::endl << SPE.Trans().size() << " transcendental constraints: " << std::endl;
  for( auto&& expr : SPE.Trans() ) DAG.output( DAG.subgraph( 1, &expr ) );

  return 0;
}

///////////////////////////////////////////////////////////////////////////////

int main()
{
  try{
//    test_dep1();
//    test_dep2();
//    test_rltred1();
//    test_rltred2();
//    test_rltred3();
//    test_rltred4();
//    test_rltred5();
//    test_rltred6();
    test_spoly1();
    test_sexpr1();
  }
  catch( mc::FFBase::Exceptions &eObj ){
    std::cerr << "Error " << eObj.ierr()
              << " in factorable function manipulation:" << std::endl
              << eObj.what() << std::endl
              << "Aborts." << std::endl;
    return eObj.ierr();
  }
  catch( mc::SparseEnv<mc::FFGraph<>>::Exceptions &eObj ){
    std::cerr << "Error " << eObj.ierr()
              << " in sparse expression manipulation:" << std::endl
              << eObj.what() << std::endl
              << "Aborts." << std::endl;
    return eObj.ierr();
  }
}

