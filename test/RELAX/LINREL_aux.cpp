// Relaxation cuts for squares and products whose operands are auxiliary
// variables or constants of a polyhedral image, rather than DAG variables:
//   W  = X0 - 2*X1     (auxiliary variable defined by a linear cut)
//   S  = W^2           (append_cuts_SQR)
//   P1 = W*X0, P2 = W*W, P3 = W*3, P4 = 3*X1, P5 = 2*3   (append_cuts_TIMES)
// The cuts are checked on a grid of points: the program returns the number of
// violated cuts, for each type of breakpoint and with or without quadratic
// cuts.
#include <iomanip>
#include <map>

#include "interval.hpp"
#include "polimage.hpp"
typedef mc::Interval I;
typedef mc::PolVar<I> PV;
typedef mc::PolCut<I> PC;

int
main()
{
  mc::FFGraph DAG;
  const unsigned NX = 2;
  mc::FFVar X[NX];
  for (unsigned i = 0; i < NX; i++) X[i].set(&DAG);
  I IX[NX] = {I(-1, 2), I(0, 3)};

  const mc::PolImg<I>::Options::REFINE BKPT[4] = {
      mc::PolImg<I>::Options::NONE, mc::PolImg<I>::Options::CONT,
      mc::PolImg<I>::Options::BIN, mc::PolImg<I>::Options::SOS2};
  unsigned nfail = 0;

  for (unsigned icase = 0; icase < 8; icase++)
  {
    mc::PolImg<I> Env;
    Env.options.BREAKPOINT_TYPE = BKPT[icase / 2];
    Env.options.ALLOW_QUAD      = icase % 2;
    Env.options.SANDWICH_MAXCUT = 6;

    PV POLX[NX];
    for (unsigned i = 0; i < NX; i++) POLX[i].set(&Env, X[i], IX[i]);

    // Auxiliary variable W = X0 - 2*X1
    const I IW = IX[0] - 2. * IX[1];
    PV W(&Env, IW);
    Env.add_cut(nullptr, PC::EQ, 0., W, 1., POLX[0], -1., POLX[1], 2.);

    // Square and products involving the auxiliary variable or constants
    PV S(&Env, mc::Op<I>::sqr(IW)), P1(&Env, IW * IX[0]), P2(&Env, IW * IW),
        P3(&Env, IW * 3.), P4(&Env, 3. * IX[1]), P5(&Env, I(6.));
    Env.append_cuts_SQR(S, W);
    Env.append_cuts_TIMES(P1, W, POLX[0]);
    Env.append_cuts_TIMES(P2, W, W);
    Env.append_cuts_TIMES(P3, W, PV(3.));
    Env.append_cuts_TIMES(P4, PV(3.), POLX[1]);
    Env.append_cuts_TIMES(P5, PV(2), PV(3.));
    if (!icase) std::cout << Env;

    // Check the cuts on a grid
    const unsigned NG = 12;
    unsigned nviol    = 0;
    for (unsigned i0 = 0; i0 <= NG; i0++)
      for (unsigned i1 = 0; i1 <= NG; i1++)
      {
        const double x0 = IX[0].l() + i0 * (IX[0].u() - IX[0].l()) / NG,
                     x1 = IX[1].l() + i1 * (IX[1].u() - IX[1].l()) / NG,
                     w  = x0 - 2. * x1;
        std::map<PV::t_idVar, double> val;
        val[POLX[0].id()] = x0;
        val[POLX[1].id()] = x1;
        val[W.id()]       = w;
        val[S.id()]       = w * w;
        val[P1.id()]      = w * x0;
        val[P2.id()]      = w * w;
        val[P3.id()]      = w * 3.;
        val[P4.id()]      = 3. * x1;
        val[P5.id()]      = 6.;

        for (auto const& pCut : Env.Cuts())
        {
          double lhs = 0.;
          for (unsigned k = 0; k < pCut->nvar(); k++)
            lhs += pCut->coef()[k] * val.at(pCut->var()[k].id());
          for (unsigned k = 0; k < pCut->nqvar(); k++)
            lhs += pCut->qcoef()[k] * val.at(pCut->qvar1()[k].id()) *
                   val.at(pCut->qvar2()[k].id());
          const double tol = 1e-10 * (1. + std::fabs(pCut->rhs()));
          switch (pCut->type())
          {
            case PC::EQ:
              if (std::fabs(lhs - pCut->rhs()) > tol) ++nviol;
              break;
            case PC::LE:
              if (lhs > pCut->rhs() + tol) ++nviol;
              break;
            case PC::GE:
              if (lhs < pCut->rhs() - tol) ++nviol;
              break;
            default:
              ++nviol;  // no special set expected without breakpoints
              break;
          }
        }
      }

    std::cout << "BREAKPOINT_TYPE=" << Env.options.BREAKPOINT_TYPE
              << "  ALLOW_QUAD=" << Env.options.ALLOW_QUAD << ": "
              << Env.Cuts().size() << " cuts, " << Env.Aux().size()
              << " auxiliaries, " << nviol << " violations" << std::endl;
    nfail += nviol;
  }

  // An auxiliary variable of another image is rejected
  {
    mc::PolImg<I> Env1, Env2;
    PV W1(&Env1, I(-1., 1.)), S2(&Env2, I(0., 1.));
    bool caught = false;
    try
    {
      Env2.append_cuts_SQR(S2, W1);
    }
    catch (mc::PolImg<I>::Exceptions& e)
    {
      caught = (e.ierr() == mc::PolImg<I>::Exceptions::ENVMIS);
    }
    std::cout << "Operand from another image: "
              << (caught ? "rejected" : "NOT rejected") << std::endl;
    if (!caught) ++nfail;
  }

  return nfail;
}
