#undef USE_SCMODEL
#include <iomanip>
#include "polimage.hpp"
#include "interval.hpp"
typedef mc::Interval I;
#ifdef USE_SCMODEL
  #include "scmodel.hpp"
  typedef mc::SCModel<I> SCM;
  typedef mc::SCVar<I> SCV;
#endif

int main()
{
  mc::FFGraph DAG;
  const unsigned NX = 5; mc::FFVar X[NX]; for( unsigned i=0; i<NX; i++ ) X[i].set( &DAG );
  //const unsigned NF = 2; mc::FFVar F[NF]; F[0] = log(X[0])+pow(X[1],2); F[1] = sin(X[0])-cos(X[1]);
  //const unsigned NF = 1; mc::FFVar F[NF]; F[0] = ( X[0] + X[0] ) * cheb( X[1], 3 ) * 3;
  //const unsigned NF = 1; mc::FFVar F[NF]; F[0] = prod( NX, X );
  //const unsigned NF = 1; mc::FFVar F[NF]; F[0] = X[0] + X[0]*sqr(X[1]) + pow(X[2],3)*X[3]*X[4];
  const unsigned NF = 1; mc::FFVar F[NF]; F[0] = fstep( X[0] ) + max( X[0], X[1] );
  std::cout << DAG;

  I IX[NX] = { I(-2,2), I(-3,3), I(-4,4), I(-1,1), I(-1,1) };
#ifndef USE_SCMODEL
  mc::PolImg<I> Env;
  mc::PolVar<I> POLX[NX]; for( unsigned i=0; i<NX; i++ ) POLX[i].set( &Env, X[i], IX[i] );
  mc::PolVar<I> POLF[NF]; 
#else
  mc::PolImg<SCV> Env;
  SCM SCMenv( 5, NX );
  SCMenv.options.MIXED_IA = true;
  SCMenv.options.MIN_FACTOR = mc::machprec();
  SCV SCX[2]; for( unsigned i=0; i<NX; i++ ) SCX[i].set( &SCMenv, i, IX[i] );
  mc::PolVar<SCV> POLX[NX]; for( unsigned i=0; i<NX; i++ ) POLX[i].set( &Env, X[i], SCX[i] );
  mc::PolVar<SCV> POLF[NF];
#endif
  Env.options.AGGREG_LQ       = 1;
  Env.options.ALLOW_QUAD      = 1;
  //Env.options.ALLOW_DISJ      = { mc::FFOp::FSTEP, mc::FFOp::MAXF };
  //Env.options.ALLOW_NLIN      = { mc::FFOp::SQR, mc::FFOp::IPOW, mc::FFOp::CHEB, mc::FFOp::FSTEP, mc::FFOp::MAXF };
  Env.options.SANDWICH_MAXCUT = 5;

  DAG.eval( NF, F, POLF, NX, X, POLX );
  Env.generate_cuts( NF, POLF );
  std::cout << Env;

  return 0;
}
