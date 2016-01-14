#include <iomanip>

#include "interval.hpp"
typedef mc::Interval  I;

#include "cmodel.hpp"
typedef mc::CModel<I> PM;
typedef mc::CVar<I> PV;

#include "polimage.hpp"

int main()
{
  mc::FFGraph DAG;
  mc::FFVar X[2]; X[0].set( &DAG ); X[1].set( &DAG );
  mc::FFVar F[2]; F[0] = log(X[0])+pow(X[1],2); F[1] = sin(X[0])-cos(X[1]); 
  std::cout << DAG;

  mc::PolImg<I> Env;
  Env.options.AGGREG_LIN = true;
  Env.options.SANDWICH_MAXCUT = 5;
  I IX[2] = { I(1,5), I(2,6) };
  mc::PolVar<I> PX[2]; PX[0].set( &Env, X[0], IX[0] ); PX[1].set( &Env, X[1], IX[1] );
  mc::PolVar<I> PF[2]; DAG.eval( 2, F, PF, 2, X, PX );
  std::cout << Env;

  Env.generate_cuts( 2, PF );
  std::cout << Env;

  return 0;
}
/*
int main()
{
  mc::FFGraph DAG;
  mc::FFVar X[2]; X[0].set( &DAG ); X[1].set( &DAG );
  mc::FFVar F[2]; F[0] = log(X[0])+pow(X[1],2); F[1] = sin(X[0])-cos(X[1]); 
  std::cout << DAG;

  PM PMenv( 2, 5, true );
  PMenv.options.MIXED_IA = true;
  //I  IX[2] = { I(1,5), I(2,6) };
  PV IX[2]; IX[0].set( &PMenv, 0, I(1,5) ); IX[1].set( &PMenv, 1, I(2,6) );

  mc::PolImg<PV> Env;
  Env.options.AGGREG_LIN = true;
  Env.options.SANDWICH_MAXCUT = 5;
  mc::PolVar<PV> PX[2]; PX[0].set( &Env, X[0], IX[0] ); PX[1].set( &Env, X[1], IX[1] );
  mc::PolVar<PV> PF[2]; DAG.eval( 2, F, PF, 2, X, PX );
  std::cout << Env;

  Env.generate_cuts( 2, PF );
  std::cout << Env;

  return 0;
}
*/
