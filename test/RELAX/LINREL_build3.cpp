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
  const unsigned NX = 3; mc::FFVar X[NX]; for( unsigned i=0; i<NX; i++ ) X[i].set( &DAG );
  const unsigned NF = 2; mc::FFVar F[NF]; F[0] = X[1] - 2*X[2]; F[1] = sqr(X[0])+X[0]*X[1]+4; 
  auto SGF  = DAG.subgraph( NF, F );
  DAG.output( SGF );
  auto SGF1 = DAG.subgraph( 1, F+1 );
  DAG.output( SGF1 );

  I IX[NX] = { I(-0.8,-0.3), I(-1e20,1e20), I(3.,4.) };

  mc::PolImg<I> Env;
  mc::PolVar<I> POLX[NX]; for( unsigned i=0; i<NX; i++ ) POLX[i].set( &Env, X[i], IX[i] );
  mc::PolVar<I> POLF[NF]; 
  std::vector<mc::PolVar<I>> POLWKF1;

  Env.options.AGGREG_LQ = true;
  Env.options.RELAX_QUAD = true;
  Env.options.SANDWICH_MAXCUT = 5;

  DAG.eval( SGF1, POLWKF1, 1, F+1, POLF+1, NX, X, POLX );
  Env.generate_cuts( 1, POLF+1 );
  std::cout << "\nCuts w/ forward pass only:\n";
  std::cout << Env;

  I IF[NF] = { 0., 0. }, IINF(-1e20,1e20);
  std::vector<I> IWKF, IWKF1;
  int flag = DAG.reval( SGF, IWKF, NF, F, IF, NX, X, IX, IINF, 10 );
  std::cout << "\nDAG interval evaluation w/ " << flag << " forward/backward passes:\n";
  for( unsigned i=0; i<NX; i++ ) std::cout << "X[" << i << "] = " << IX[i] << std::endl;
  for( unsigned i=0; i<NF; i++ ) std::cout << "F[" << i << "] = " << IF[i] << std::endl << std::endl;
  for( auto && wk : IWKF ) std::cout << wk << std::endl;
  std::cout << std::endl;
  DAG.wkextract( SGF1, IWKF1, SGF, IWKF );
  for( auto && wk : IWKF1 ) std::cout << wk << std::endl;

  Env.reset_cuts();
  DAG.eval( SGF1, POLWKF1, 1, F+1, POLF+1, NX, X, POLX );
  for( unsigned i=0; i<POLWKF1.size(); i++ ) POLWKF1[i].update( IWKF1[i] );
  Env.generate_cuts( 1, POLF+1 );
  std::cout << "\nCuts w/ " << flag << " forward/backward passes:\n";
  std::cout << Env;

  return 0;
}
