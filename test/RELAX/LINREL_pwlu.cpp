// =============================================================================
// LINREL_pwlu.cpp
// Polyhedral relaxations of piecewise-linear univariates (mc::FFPWLU), for
// every encoding (HULL, CONT, BIN, SOS2), with the functions of an operation
// on one shared subdivision or on one each (Options::SEPARATE):
//  1. Validity, without solver: every cut holds at each point (x, f(x)) of the
//     graph, for x on a grid of the operand range. The variables of the
//     subdivision of an operand are found from the cut linking them to it
//     (CONT, BIN: x = x_0 + sum_k (x_{k+1}-x_k) d_k; SOS2: x = sum_k x_k w_k)
//     and take the values that encode x; the binaries of BIN get ranges from
//     the cuts in which they are the only unknown.
//  2. Tightness of HULL, without solver: at each x, the bounds on f(x) implied
//     by the cuts are the convex and concave envelopes (mc::PWLU::hull).
//  3. Tightness of every encoding, with Gurobi (MC__USE_GUROBI): for x fixed,
//     the minimum and maximum of each result over the relaxation are f(x) for
//     BIN and SOS2, and the envelopes for HULL and CONT.
// Cases: four functions of one operand on ranges inside, across and outside
// their breakpoints, a single point, bounds at breakpoints; and a two-layer
// Kolmogorov-Arnold network whose second layer acts on intermediate variables.
// The program returns the number of failed cases.
// =============================================================================

#include <cmath>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include "interval.hpp"
typedef mc::Interval I;
#include "ffpwlu.hpp"
typedef mc::PolVar<I> PV;
typedef mc::PolCut<I> PC;
typedef mc::FFPWLU<I> OP;

#if defined(MC__USE_GUROBI)
#include "gurobi_c++.h"
#endif

// Values of the auxiliaries of the subdivisions of the operands, from the cuts
// linking them to the operands (taken in order, so that the output cut of a
// first layer, which has the same form, is not mistaken for the linking cut of
// the second); the binaries of BIN are left unknown
static void
subdiv_values(mc::PolImg<I> const& Env, std::vector<PV::t_idVar> const& operands,
              std::map<PV::t_idVar, double>& val)
{
  std::set<PV::t_idVar> aux;
  for (auto const* pa : Env.Aux()) aux.insert(pa->id());
  for (auto const& idx : operands)
  {
    double const x = val.at(idx);
    for (auto const& pCut : Env.Cuts())
    {
      if (pCut->type() != PC::EQ || pCut->nvar() < 3) continue;
      int ix = -1;
      bool ok = true;
      for (unsigned k = 0; ok && k < pCut->nvar(); ++k)
      {
        auto const id = pCut->var()[k].id();
        if (id == idx)
          ix = k;
        else
          ok = aux.count(id) && !val.count(id);
      }
      if (!ok || ix != (int)pCut->nvar() - 1) continue;
      unsigned const n = pCut->nvar() - 1;
      double const a   = pCut->coef()[ix];
      std::vector<double> knot;
      if (a == 1.)  // CONT, BIN: sum_k (x_k-x_{k+1}) d_k + x = x_0
      {
        knot.push_back(pCut->rhs());
        for (unsigned k = 0; k < n; ++k) knot.push_back(knot.back() - pCut->coef()[k]);
        for (unsigned k = 0; k < n; ++k)
          val[pCut->var()[k].id()] =
              std::min(1., std::max(0., (x - knot[k]) / (knot[k + 1] - knot[k])));
      }
      else if (a == -1. && pCut->rhs() == 0.)  // SOS2: sum_k x_k w_k - x = 0
      {
        for (unsigned k = 0; k < n; ++k) knot.push_back(pCut->coef()[k]);
        size_t j = 0;
        while (j + 2 < knot.size() && x > knot[j + 1]) ++j;
        for (unsigned k = 0; k < n; ++k) val[pCut->var()[k].id()] = 0.;
        val[pCut->var()[j].id()]     = (knot[j + 1] - x) / (knot[j + 1] - knot[j]);
        val[pCut->var()[j + 1].id()] = (x - knot[j]) / (knot[j + 1] - knot[j]);
      }
    }
  }
}

// Number of continuous auxiliaries without a value: the subdivision variables
// not found from a linking cut
static unsigned
unfound(mc::PolImg<I> const& Env, std::map<PV::t_idVar, double> const& val)
{
  unsigned n = 0;
  for (auto const* pa : Env.Aux())
    if (pa->id().first == PV::AUXCONT && !val.count(pa->id())) ++n;
  return n;
}

// Largest violation of the cuts at the values in val, the other auxiliaries
// ranging over the bounds derived from the cuts with one unknown; -1 if a cut
// has unknowns that are not auxiliaries
static double
violation(mc::PolImg<I> const& Env, std::map<PV::t_idVar, double> const& val)
{
  std::map<PV::t_idVar, std::pair<double, double>> rng;
  for (auto const* pa : Env.Aux())
    if (!val.count(pa->id())) rng[pa->id()] = {mc::Op<I>::l(pa->range()), mc::Op<I>::u(pa->range())};

  struct Split
  {
    double known, scale;
    std::vector<std::pair<PV::t_idVar, double>> unk;
  };
  auto split = [&](PC const& cut) {
    Split sp{0., 1. + std::fabs(cut.rhs()), {}};
    for (unsigned k = 0; k < cut.nvar(); ++k)
    {
      auto itv = val.find(cut.var()[k].id());
      if (itv != val.end())
      {
        sp.known += cut.coef()[k] * itv->second;
        sp.scale += std::fabs(cut.coef()[k] * itv->second);
      }
      else if (cut.coef()[k] != 0.)
        sp.unk.push_back({cut.var()[k].id(), cut.coef()[k]});
    }
    return sp;
  };

  for (unsigned pass = 0; pass < 5; ++pass)
    for (auto const& pCut : Env.Cuts())
    {
      if (pCut->type() == PC::SOS2) continue;
      Split const sp = split(*pCut);
      if (sp.unk.size() != 1 || !rng.count(sp.unk[0].first)) continue;
      auto& [lo, hi]   = rng[sp.unk[0].first];
      double const a   = sp.unk[0].second, r = (pCut->rhs() - sp.known) / a;
      bool const upper = (pCut->type() == PC::LE) == (a > 0.);
      if (pCut->type() == PC::EQ || upper) hi = std::min(hi, r);
      if (pCut->type() == PC::EQ || !upper) lo = std::max(lo, r);
    }

  double worst = 0.;
  for (auto const& pCut : Env.Cuts())
  {
    if (pCut->type() == PC::SOS2) continue;  // holds by construction of the weights
    Split const sp = split(*pCut);
    double smin = 0., smax = 0.;
    for (auto const& [id, a] : sp.unk)
    {
      auto itr = rng.find(id);
      if (itr == rng.end()) return -1.;
      double const lo = itr->second.first, hi = itr->second.second;
      if (lo > hi) worst = std::max(worst, (lo - hi) / sp.scale);
      smin += std::min(a * lo, a * hi);
      smax += std::max(a * lo, a * hi);
    }
    double viol = 0.;
    switch (pCut->type())
    {
      case PC::EQ: viol = std::max(sp.known + smin - pCut->rhs(), pCut->rhs() - sp.known - smax); break;
      case PC::LE: viol = sp.known + smin - pCut->rhs(); break;
      case PC::GE: viol = pCut->rhs() - sp.known - smax; break;
      default: break;
    }
    worst = std::max(worst, viol / sp.scale);
  }
  return worst;
}

int
main()
{
  using mc::FFVar;
  using mc::PWLU;
  int const ENC[4]       = {OP::Options::HULL, OP::Options::CONT, OP::Options::BIN, OP::Options::SOS2};
  char const* ENCNAME[4] = {"HULL", "CONT", "BIN", "SOS2"};
  unsigned nfail = 0, nrun = 0;

  // Functions of one operand: a hat, a ReLU form (merged breakpoints, a
  // constant term), a steep piece between breakpoints 1e-4 apart, an affine
  // function
  std::vector<PWLU> const F = {
      PWLU({-1., 0., 2.}, {1., 0., 1.}, 0., 0.5),
      PWLU::relu(0.1, -0.5, {1.2, -0.4, 0.7, 0.3}, {1., -2., 0., 2.}, {0.3, 0.2, 0.5, 0.6}),
      PWLU({-0.2, -0.1999, 0.4}, {0., 3., -2.}, -1., 0.),
      PWLU({0.}, {0.25}, 1.5, 1.5)};
  std::vector<std::pair<double, double>> const ranges = {
      {-1.5, 2.5}, {-0.25, 0.45}, {-0.2, 0.4}, {0.5, 1.5}, {-3., -2.}, {0.3, 0.3}};

  std::cout << "Four functions of one operand" << std::endl;
  for (auto const& [xl, xu] : ranges)
  {
    std::ostringstream msg;
    for (unsigned ie = 0; ie < 4; ++ie)
      for (bool sep : {false, true})
      {
        mc::FFGraph DAG;
        FFVar X(&DAG);
        OP Op;
        Op.options.ENCODING = ENC[ie];
        Op.options.SEPARATE = sep;
        std::vector<FFVar> Y = Op(X, F);
        mc::PolImg<I> Env;
        PV PX(&Env, X, I(xl, xu)), PY[4];
        DAG.eval(4, Y.data(), PY, 1, &X, &PX);
        Env.generate_cuts(4, PY);
        ++nrun;

        double worst = 0., gap = 0.;
        unsigned const NG = (xu > xl ? 200 : 0);
        for (unsigned ig = 0; ig <= NG; ++ig)
        {
          double const x = (NG ? xl + (xu - xl) * ig / NG : xl);
          std::map<PV::t_idVar, double> val;
          val[PX.id()] = x;
          for (int j = 0; j < 4; ++j) val[PY[j].id()] = F[j].val(x, true);
          subdiv_values(Env, {PX.id()}, val);
          double const v = (unfound(Env, val) ? -1. : violation(Env, val));
          worst          = (v < 0. || worst < 0. ? -1. : std::max(worst, v));

          // Bounds on each result implied by the HULL cuts at x, against the envelopes
          if (ENC[ie] != OP::Options::HULL) continue;
          for (int j = 0; j < 4; ++j)
          {
            double lo = -INFINITY, up = INFINITY;
            for (auto const& pCut : Env.Cuts())
            {
              // cuts in Y_j and X only
              double ay = 0., ax = 0.;
              bool other = false;
              for (unsigned k = 0; k < pCut->nvar(); ++k)
              {
                auto const id = pCut->var()[k].id();
                if (id == PY[j].id())
                  ay += pCut->coef()[k];
                else if (id == PX.id())
                  ax += pCut->coef()[k];
                else
                  other = true;
              }
              if (other || ay == 0. || pCut->type() == PC::SOS2) continue;
              double const r   = (pCut->rhs() - ax * x) / ay;
              bool const upper = (pCut->type() == PC::LE) == (ay > 0.);
              if (pCut->type() == PC::EQ || upper) up = std::min(up, r);
              if (pCut->type() == PC::EQ || !upper) lo = std::max(lo, r);
            }
            PWLU const R = F[j].restrict(xl, xu);
            double const hl = (xu > xl ? R.hull(true).val(x, true) : F[j].val(x, true)),
                         hu = (xu > xl ? R.hull(false).val(x, true) : F[j].val(x, true));
            gap = std::max({gap, std::fabs(lo - hl) / (1. + std::fabs(hl)), std::fabs(up - hu) / (1. + std::fabs(hu))});
          }
        }
        if (worst < 0. || worst > 1e-9 || gap > 1e-9)
        {
          ++nfail;
          msg << "  " << ENCNAME[ie] << (sep ? "/sep" : "") << ": "
              << (worst < 0. ? std::string("subdivision variables not found")
                             : "violation " + std::to_string(worst) + ", envelope gap " + std::to_string(gap))
              << ";";
        }

#if defined(MC__USE_GUROBI)
        // Minimum and maximum of each result for x fixed at a few points
        try
        {
          GRBEnv env(true);
          env.set(GRB_IntParam_OutputFlag, 0);
          env.start();
          double gapg = 0.;
          for (unsigned ig = 0; ig <= 4; ++ig)
          {
            double const x = xl + (xu - xl) * ig / 4.;
            GRBModel model(env);
            model.getEnv().set(GRB_DoubleParam_FeasibilityTol, 1e-9);
            model.getEnv().set(GRB_DoubleParam_OptimalityTol, 1e-9);
            model.getEnv().set(GRB_DoubleParam_IntFeasTol, 1e-9);
            model.getEnv().set(GRB_DoubleParam_MIPGap, 0.);
            std::map<PV::t_idVar, GRBVar> var;
            for (auto const& [id, pv] : Env.Vars())
              var[pv->id()] = model.addVar(mc::Op<I>::l(pv->range()), mc::Op<I>::u(pv->range()), 0., GRB_CONTINUOUS);
            for (auto const* pa : Env.Aux())
              var[pa->id()] = model.addVar(mc::Op<I>::l(pa->range()), mc::Op<I>::u(pa->range()), 0.,
                                           pa->id().first == PV::AUXINT ? GRB_BINARY : GRB_CONTINUOUS);
            var[PX.id()].set(GRB_DoubleAttr_LB, x);
            var[PX.id()].set(GRB_DoubleAttr_UB, x);
            for (auto const& pCut : Env.Cuts())
            {
              GRBLinExpr lhs;
              std::vector<GRBVar> sv;
              std::vector<double> sw;
              for (unsigned k = 0; k < pCut->nvar(); ++k)
              {
                lhs += pCut->coef()[k] * var.at(pCut->var()[k].id());
                sv.push_back(var.at(pCut->var()[k].id()));
                sw.push_back(k + 1.);
              }
              switch (pCut->type())
              {
                case PC::EQ: model.addConstr(lhs, GRB_EQUAL, pCut->rhs()); break;
                case PC::LE: model.addConstr(lhs, GRB_LESS_EQUAL, pCut->rhs()); break;
                case PC::GE: model.addConstr(lhs, GRB_GREATER_EQUAL, pCut->rhs()); break;
                case PC::SOS2: model.addSOS(sv.data(), sw.data(), sv.size(), GRB_SOS_TYPE2); break;
                default: break;
              }
            }
            for (int j = 0; j < 4; ++j)
            {
              PWLU const R  = F[j].restrict(xl, xu);
              bool const lp = (ENC[ie] == OP::Options::HULL || ENC[ie] == OP::Options::CONT) && xu > xl;
              double const ref[2] = {lp ? R.hull(true).val(x, true) : F[j].val(x, true),
                                     lp ? R.hull(false).val(x, true) : F[j].val(x, true)};
              for (int s = 0; s < 2; ++s)
              {
                // minimum of Y, then of -Y
                model.setObjective(GRBLinExpr(var.at(PY[j].id()), s ? -1. : 1.));
                model.optimize();
                double const obj = (model.get(GRB_IntAttr_Status) == GRB_OPTIMAL
                                        ? (s ? -1. : 1.) * model.get(GRB_DoubleAttr_ObjVal) : NAN);
                gapg = std::max(gapg, std::isnan(obj) ? INFINITY : std::fabs(obj - ref[s]) / (1. + std::fabs(ref[s])));
              }
            }
          }
          if (gapg > 1e-6)
          {
            ++nfail;
            msg << "  " << ENCNAME[ie] << (sep ? "/sep" : "") << ": Gurobi bounds off by " << gapg << ";";
          }
        }
        catch (GRBException& e)
        {
          ++nfail;
          msg << "  " << ENCNAME[ie] << (sep ? "/sep" : "") << ": Gurobi error " << e.getMessage() << ";";
        }
#endif
      }
    std::ostringstream name;
    name << "X in [" << xl << "," << xu << "]";
    std::cout << (msg.str().empty() ? "  ok    " : "  FAIL  ") << std::left << std::setw(22) << name.str() << msg.str()
              << std::endl;
  }

  // Two-layer Kolmogorov-Arnold network: U_j = A_j(X1) + B_j(X2), Z = C(U_0) + D(U_1)
  std::cout << "Two-layer network, second layer on intermediate variables" << std::endl;
  {
    std::vector<PWLU> const A = {PWLU::relu(0.2, 0.5, {1., -0.8}, {2., -1.5}, {0.4, 0.3}),
                                 PWLU({-0.5, 0.2, 0.6}, {0.3, -0.2, 0.5}, -0.4, 1.)};
    std::vector<PWLU> const B = {PWLU({0.5, 1.5}, {0., 1.}, 0.2, -0.3),
                                 PWLU::relu(-0.1, 0.3, {0.6, 0.9}, {-1., 1.}, {1.2, -0.8})};
    PWLU const C              = PWLU({-0.3, 0.4, 1.1}, {0.5, -0.4, 0.8}, -1., 0.3);
    PWLU const D              = PWLU::relu(0., 0.4, {1., -2.}, {1., 1.}, {0.1, -0.6});
    std::ostringstream msg;
    for (unsigned ie = 0; ie < 4; ++ie)
      for (bool sep : {false, true})
      {
        mc::FFGraph DAG;
        FFVar X1(&DAG), X2(&DAG);
        OP Op;
        Op.options.ENCODING = ENC[ie];
        Op.options.SEPARATE = sep;
        std::vector<FFVar> YA = Op(X1, A), YB = Op(X2, B);
        FFVar U[2] = {YA[0] + YB[0], YA[1] + YB[1]};
        FFVar Z    = Op(U[0], C) + Op(U[1], D);
        mc::PolImg<I> Env;
        Env.options.AGGREG_LQ = false;
        PV PX1(&Env, X1, I(-1., 1.)), PX2(&Env, X2, I(0., 2.)), PZ;
        PV PX[2] = {PX1, PX2};
        FFVar X[2] = {X1, X2};
        DAG.eval(1, &Z, &PZ, 2, X, PX);
        Env.generate_cuts(1, &PZ);
        ++nrun;

        // Every DAG variable of the image at the point, then the subdivisions
        std::vector<FFVar> dvar;
        std::vector<PV const*> pvar;
        for (auto const& [id, pv] : Env.Vars())
        {
          dvar.push_back(pv->var());
          pvar.push_back(pv);
        }
        std::vector<PV::t_idVar> operands = {PX1.id(), PX2.id()};
        for (auto const& [id, pv] : Env.Vars())
          if (pv->var().id() == U[0].id() || pv->var().id() == U[1].id()) operands.push_back(pv->id());
        std::vector<double> dval(dvar.size());
        double worst = 0.;
        for (unsigned i1 = 0; i1 <= 20; ++i1)
          for (unsigned i2 = 0; i2 <= 20; ++i2)
          {
            double const x[2] = {-1. + 2. * i1 / 20., 2. * i2 / 20.};
            DAG.eval(dvar.size(), dvar.data(), dval.data(), 2, X, x);
            std::map<PV::t_idVar, double> val;
            for (size_t k = 0; k < dvar.size(); ++k) val[pvar[k]->id()] = dval[k];
            subdiv_values(Env, operands, val);
            double const v = (unfound(Env, val) ? -1. : violation(Env, val));
            worst          = (v < 0. || worst < 0. ? -1. : std::max(worst, v));
          }
        if (worst < 0. || worst > 1e-9)
        {
          ++nfail;
          msg << "  " << ENCNAME[ie] << (sep ? "/sep" : "") << ": "
              << (worst < 0. ? std::string("subdivision variables not found") : "violation " + std::to_string(worst)) << ";";
        }
      }
    std::cout << (msg.str().empty() ? "  ok    " : "  FAIL  ") << std::left << std::setw(22) << "X1 x X2 grid 21x21"
              << msg.str() << std::endl;
  }

  std::cout << nrun << " relaxations, " << nfail << " failed" << std::endl;
  return nfail;
}
