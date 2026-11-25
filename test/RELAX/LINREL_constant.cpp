// =============================================================================
// LINREL_constant.cpp
// Functions that do not depend on the variables, e.g. the zero objective of a
// feasibility problem: cuts are generated for the other dependents of a
// polyhedral image (PolImg::generate_cuts), and the zero polynomial is
// transcribed into the DAG as the constant 0 (SLiftEnv::insert_dag).
//
// The process exits with code 0 if every test passes, 1 otherwise.
// =============================================================================

#include <iostream>
#include <string>

#include "interval.hpp"
typedef mc::Interval I;
#include "polimage.hpp"
#include "slift.hpp"

static int npass = 0, nfail = 0;
static void
check(bool const ok, std::string const& what)
{
  std::cout << (ok ? "  PASS  " : "  FAIL  ") << what << std::endl;
  (ok ? npass : nfail)++;
}

int
main()
{
  mc::FFGraph DAG;
  mc::FFVar X[2];
  for (auto& x : X) x.set(&DAG);
  I IX[2] = {I(-1, 2), I(0, 3)};

  // Cuts with a constant dependent, as many as without it
  size_t ncuts[2] = {0, 0};
  for (unsigned withcst = 0; withcst < 2; ++withcst)
  {
    mc::FFVar F[2] = {X[0] * X[1] + exp(X[0]), 0.};
    mc::PolImg<I> Env;
    mc::PolVar<I> POLX[2], POLF[2];
    for (unsigned i = 0; i < 2; i++) POLX[i].set(&Env, X[i], IX[i]);
    DAG.eval(1 + withcst, F, POLF, 2, X, POLX);
    if (withcst) POLF[1] = mc::PolVar<I>(0.);
    try
    {
      Env.generate_cuts(1 + withcst, POLF);
      ncuts[withcst] = Env.Cuts().size();
    }
    catch (...)
    {
      ncuts[withcst] = 0;
    }
  }
  check(ncuts[0] > 0 && ncuts[1] == ncuts[0],
        "PolImg::generate_cuts with a constant dependent: " + std::to_string(ncuts[1]) +
            " cuts, as without it");

  // Zero polynomial
  mc::SLiftEnv SPE(&DAG);
  mc::SPoly<mc::FFVar const*, mc::lt_FFVar> P0;
  mc::FFVar Z;
  bool thrown = false;
  try
  {
    Z = SPE.insert_dag(P0);
  }
  catch (...)
  {
    thrown = true;
  }
  check(!thrown && Z.cst() && Z.num().val() == 0., "SLiftEnv::insert_dag of the zero polynomial: the constant 0");

  std::cout << npass << " passed, " << nfail << " failed" << std::endl;
  return nfail ? 1 : 0;
}
