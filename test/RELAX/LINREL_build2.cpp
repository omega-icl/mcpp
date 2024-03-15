#include <iomanip>

#include "polimage.hpp"
#include "interval.hpp"

typedef mc::Interval I;

int main()
{
  mc::FFGraph DAG;
  const unsigned NX = 2; mc::FFVar X[NX]; for( unsigned i=0; i<NX; i++ ) X[i].set( &DAG );
  const unsigned NF = 1; mc::FFVar F[NF]; F[0] = sqr(X[0])+X[0]*X[1]+4;
  auto SGF  = DAG.subgraph( NF, F );
  DAG.output( SGF );

  I IX[NX] = { I(-0.8,-0.3), I(6.,9.) };

  mc::PolImg<I> Env;
  mc::PolVar<I> POLX[NX]; for( unsigned i=0; i<NX; i++ ) POLX[i].set( &Env, X[i], IX[i] );
  mc::PolVar<I> POLF[NF]; 
  std::vector<mc::PolVar<I>> POLWKF;

  Env.options.AGGREG_LQ       = 1;
  Env.options.ALLOW_QUAD      = 0;
  Env.options.SANDWICH_MAXCUT = 5;

  DAG.eval( SGF, POLWKF, NF, F, POLF, NX, X, POLX );
  Env.generate_cuts( NF, POLF );
  std::cout << "\nCuts w/ forward pass only:\n";
  std::cout << Env;


  I IF[NF] = { 0. }, IINF(-1e20,1e20);
  std::vector<I> IWKF;
  int flag = DAG.reval( SGF, IWKF, NF, F, IF, NX, X, IX, IINF, 10 );
  std::cout << "\nDAG interval evaluation w/ " << flag << " forward/backward passes:\n";
  for( unsigned i=0; i<NX; i++ ) std::cout << "X[" << i << "] = " << IX[i] << std::endl;
  for( unsigned i=0; i<NF; i++ ) std::cout << "F[" << i << "] = " << IF[i] << std::endl;
  for( auto && wk : IWKF ) std::cout << wk << std::endl;

  Env.reset_cuts();
  DAG.eval( SGF, POLWKF, NF, F, POLF, NX, X, POLX );
  for( unsigned i=0; i<POLWKF.size(); i++ ) POLWKF[i].update( IWKF[i] );
  Env.generate_cuts( NF, POLF );
  std::cout << "\nCuts w/ " << flag << " forward/backward passes:\n";
  std::cout << Env;

  return 0;
}
