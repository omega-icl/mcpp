#undef USE_SCMODEL
#include <iomanip>

#include "interval.hpp"
#include "polimage.hpp"
typedef mc::Interval I;
#ifdef USE_SCMODEL
#include "scmodel.hpp"
typedef mc::SCModel<I> SCM;
typedef mc::SCVar<I> SCV;
#endif

int
main()
{
  mc::FFGraph DAG;
  const unsigned NX        = 2;
  std::vector<mc::FFVar> X = DAG.add_vars(NX, "X");
  std::vector<mc::FFVar> F = {X[0] * sqr(exp(X[0]) - X[1])};
  std::cout << DAG;

  std::vector<I> IX{I(-2, 1), I(-1, 2)};
  mc::PolImg<I> Env;
  std::vector<mc::PolVar<I>> POLX(NX);
  for (unsigned i = 0; i < NX; i++) POLX[i].set(&Env, X[i], IX[i]);
  std::vector<mc::PolVar<I>> POLF;

  // Env.options.AGGREG_LQ       = 1;
  // Env.options.ALLOW_QUAD      = 1;
  // Env.options.ALLOW_DISJ      = { mc::FFOp::FSTEP, mc::FFOp::MAXF };
  Env.options.ALLOW_NLIN = {mc::FFOp::EXP};
  // Env.options.SANDWICH_MAXCUT = 5;

  DAG.eval(F, POLF, X, POLX);
  Env.generate_cuts(POLF);
  std::cout << Env;

  return 0;
}
