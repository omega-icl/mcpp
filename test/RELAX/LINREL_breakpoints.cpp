// Polyhedral relaxations of univariate operations, and of sums of them on the
// same operand, when the operand has breakpoints: every cut must hold at each
// point (x, f(x)) of the graph, for x on a grid of the operand range. The
// variables of the subdivision of an operand take the values that encode x:
//   CONT and BIN: x = x_0 + sum_k (x_{k+1}-x_k) l_k, with l_k = 1 below the
//                 interval of x, its fraction in it, 0 above; the binaries
//                 of BIN are 1 below the interval of x, 0 from it on;
//   SOS2:         x = sum_k w_k x_k, with the weights of the two knots around
//   x.
// The intermediate variables of the DAG take their values at x. Auxiliary
// variables that x does not determine (e.g. those of the decomposition of a
// convexo-concave function into convex and concave parts) are checked through
// the cuts: ranges are derived for them from the cuts in which they are the
// only unknown, and every other cut must hold for some values in these ranges
// (binary auxiliaries are relaxed to [0,1] there). Covers the
// tangent points of the outer approximation (taken from the breakpoints in
// the part of the range being approximated) and the secants of sin and cos,
// and checks the other operations for regressions.
// The program returns the number of failed cases.
#include <cmath>
#include <functional>
#include <iomanip>
#include <map>
#include <string>
#include <vector>

#include "interval.hpp"
#include "polimage.hpp"
typedef mc::Interval I;
typedef mc::PolVar<I> PV;
typedef mc::PolCut<I> PC;
typedef mc::PolImg<I>::Options::REFINE REFINE;

struct Case
{
  std::string name;
  std::function<mc::FFVar(mc::FFVar const&)> op;  // expression in X
  double xl, xu;                                  // operand range
  std::vector<double> bkpt;                       // breakpoints of the operand
};

// Values of the variables of the subdivision of a variable, at the value v
static void
subdiv_values(PV const& var, double const v, REFINE const type,
              std::map<PV::t_idVar, double>& val)
{
  auto const& sub  = var.subdiv();
  auto const& knot = sub.first;
  auto const& svar = sub.second;
  if (svar.empty() || knot.size() < 2) return;
  size_t const NINT = knot.size() - 1;
  size_t j          = 0;  // interval holding v
  while (j + 1 < NINT && v > knot[j + 1]) ++j;
  switch (type)
  {
    case mc::PolImg<I>::Options::CONT:
    case mc::PolImg<I>::Options::BIN:
      for (size_t k = 0; k < NINT && k < svar.size(); ++k)
      {
        double const lk =
            (v - knot[k]) / (knot[k + 1] - knot[k]);  // fraction of interval k
        val[svar[k].id()] = std::min(1., std::max(0., lk));
      }
      for (size_t k = 0; NINT + k < svar.size(); ++k)  // binaries of BIN
        val[svar[NINT + k].id()] = (k < j ? 1. : 0.);
      break;
    case mc::PolImg<I>::Options::SOS2:
      for (size_t k = 0; k < svar.size(); ++k) val[svar[k].id()] = 0.;
      val[svar[j].id()]     = (knot[j + 1] - v) / (knot[j + 1] - knot[j]);
      val[svar[j + 1].id()] = (v - knot[j]) / (knot[j + 1] - knot[j]);
      break;
    default:
      break;
  }
}

int
main()
{
  using mc::FFVar;

  // clang-format off
  const std::vector<Case> cases = {
    { "sin, X=[9,25]",        [](FFVar const& x){ return sin(x); },  9., 25., { 15.264, 20. } },
    { "sin, X=[2.7,7.5]",     [](FFVar const& x){ return sin(x); },  2.7, 7.5, { 4.5793 } },
    { "sin, X=[0.5,4]",       [](FFVar const& x){ return sin(x); },  0.5, 4., { 1.2, 3.6 } },
    { "sin, X=[-1,1]",        [](FFVar const& x){ return sin(x); },  -1., 1., { -0.4, 0.7 } },
    { "sin, X=[-7,-1]",       [](FFVar const& x){ return sin(x); },  -7., -1., { -5., -2.2 } },
    { "sin, X=[0.5,2.5]",     [](FFVar const& x){ return sin(x); },  0.5, 2.5, { 1., 2. } },
    { "cos, X=[1,12]",        [](FFVar const& x){ return cos(x); },  1., 12., { 2., 7.5 } },
    { "cos, X=[0.5,4]",       [](FFVar const& x){ return cos(x); },  0.5, 4., { 1., 3.5 } },
    { "cos, X=[-5,2]",        [](FFVar const& x){ return cos(x); },  -5., 2., { -3., 0.5 } },
    { "cos, X=[2,4.5]",       [](FFVar const& x){ return cos(x); },  2., 4.5, { 2.5, 4. } },
    { "tan, X=[-1,1.2]",      [](FFVar const& x){ return tan(x); },  -1., 1.2, { -0.5, 0.8 } },
    { "tanh, X=[-2,3]",       [](FFVar const& x){ return tanh(x); }, -2., 3., { -1.5, 0.3, 2. } },
    { "atan, X=[-2,3]",       [](FFVar const& x){ return atan(x); }, -2., 3., { -1.5, 0.3, 2. } },
    { "erf, X=[-2,1]",        [](FFVar const& x){ return erf(x); },  -2., 1., { -1., 0.5 } },
    { "pow(x,3), X=[-2,3]",   [](FFVar const& x){ return pow(x, 3); }, -2., 3., { -1.5, 0.3, 2. } },
    { "pow(x,3), X=[-1,3]",   [](FFVar const& x){ return pow(x, 3); }, -1., 3., { -0.5, 2. } },
    { "pow(x,5), X=[-1.5,1]", [](FFVar const& x){ return pow(x, 5); }, -1.5, 1., { -1., 0.5 } },
    { "pow(x,4), X=[-2,1.5]", [](FFVar const& x){ return pow(x, 4); }, -2., 1.5, { -1., 0.5 } },
    { "pow(x,-2), X=[0.5,4]", [](FFVar const& x){ return pow(x, -2); }, 0.5, 4., { 1., 3. } },
    { "pow(x,1.5), X=[0.5,4]",[](FFVar const& x){ return pow(x, 1.5); }, 0.5, 4., { 1., 3. } },
    { "sqr, X=[-1,2]",        [](FFVar const& x){ return sqr(x); },  -1., 2., { 0.3, 1.2 } },
    { "exp, X=[-1,2]",        [](FFVar const& x){ return exp(x); },  -1., 2., { 0.3, 1.2 } },
    { "log, X=[0.5,4]",       [](FFVar const& x){ return log(x); },  0.5, 4., { 1., 3. } },
    { "sqrt, X=[0.5,4]",      [](FFVar const& x){ return sqrt(x); }, 0.5, 4., { 1., 3. } },
    { "inv, X=[0.5,4]",       [](FFVar const& x){ return inv(x); },  0.5, 4., { 1., 3. } },
    { "xlog, X=[0.2,4]",      [](FFVar const& x){ return xlog(x); }, 0.2, 4., { 1., 3. } },
    { "fabs, X=[-2,3]",       [](FFVar const& x){ return fabs(x); }, -2., 3., { -1., 1.5 } },
    { "cosh, X=[-2,1]",       [](FFVar const& x){ return cosh(x); }, -2., 1., { -1., 0.5 } },
    { "sinh, X=[-2,1]",       [](FFVar const& x){ return sinh(x); }, -2., 1., { -1., 0.5 } },
    // several operations on the same operand, which share its subdivision
    { "exp(x)+sin(x), X=[0.5,4]",        [](FFVar const& x){ return exp(x) + sin(x); }, 0.5, 4., { 1.2, 3.6 } },
    { "sin(x)+cos(x)+x^3, X=[-2,6]",     [](FFVar const& x){ return sin(x) + cos(x) + pow(x, 3); }, -2., 6., { -1., 1.5, 4. } },
    { "sin(x)+sin(10x/3), X=[2.7,7.5]",  [](FFVar const& x){ return sin(x) + sin(10. * x / 3.); }, 2.7, 7.5, { 4.5793, 5.1457 } },
  };
  // clang-format on

  const REFINE TYPES[4] = {
      mc::PolImg<I>::Options::NONE, mc::PolImg<I>::Options::CONT,
      mc::PolImg<I>::Options::BIN, mc::PolImg<I>::Options::SOS2};
  const char* TNAME[4] = {"NONE", "CONT", "BIN", "SOS2"};
  const unsigned NG    = 400;
  unsigned nfail = 0, nrun = 0;

  for (auto const& c : cases)
  {
    std::ostringstream msg;
    for (unsigned it = 0; it < 4; ++it)
    {
      for (unsigned withbk = 0; withbk < 2; ++withbk)
      {
        mc::FFGraph DAG;
        FFVar X(&DAG);
        FFVar F = c.op(X);

        mc::PolImg<I> Env;
        Env.options.BREAKPOINT_TYPE = TYPES[it];
        Env.options.AGGREG_LQ       = false;
        Env.options.SANDWICH_MAXCUT = 6;
        PV PX(&Env, X, I(c.xl, c.xu)), PF;
        if (withbk)
          for (double b : c.bkpt) PX.add_breakpt(b);
        DAG.eval(1, &F, &PF, 1, &X, &PX);
        Env.generate_cuts(1, &PF);
        ++nrun;

        // DAG variables of the image, evaluated together at each x
        std::vector<FFVar> dvar;
        std::vector<PV const*> pvar;
        for (auto const& [id, pv] : Env.Vars())
        {
          dvar.push_back(pv->var());
          pvar.push_back(pv);
        }
        std::vector<double> dval(dvar.size());

        unsigned nviol = 0, nunk = 0;
        double worst = 0., xworst = 0.;
        auto record = [&](double const viol, double const x)
        {
          ++nviol;
          if (viol > worst)
          {
            worst  = viol;
            xworst = x;
          }
        };
        for (unsigned ig = 0; ig <= NG; ++ig)
        {
          double const x = c.xl + (c.xu - c.xl) * ig / NG;
          DAG.eval(dvar.size(), dvar.data(), dval.data(), 1, &X, &x);
          std::map<PV::t_idVar, double> val;
          for (size_t k = 0; k < dvar.size(); ++k) val[pvar[k]->id()] = dval[k];
          for (size_t k = 0; k < dvar.size(); ++k)
            subdiv_values(*pvar[k], dval[k], TYPES[it], val);

          // Ranges of the auxiliary variables that x does not determine
          std::map<PV::t_idVar, std::pair<double, double>> rng;
          for (auto const* pa : Env.Aux())
            if (!val.count(pa->id()))
              rng[pa->id()] = {mc::Op<I>::l(pa->range()),
                               mc::Op<I>::u(pa->range())};

          // Known part, unknown terms and scale of a cut
          struct Split
          {
            double known, scale;
            std::vector<std::pair<PV::t_idVar, double>> unk;
          };
          auto split = [&](PC const& cut)
          {
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

          // Ranges from the cuts with one unknown, a few passes
          for (unsigned pass = 0; pass < 5; ++pass)
            for (auto const& pCut : Env.Cuts())
            {
              if (pCut->type() == PC::SOS2 || pCut->nqvar()) continue;
              Split const sp = split(*pCut);
              if (sp.unk.size() != 1 || !rng.count(sp.unk[0].first)) continue;
              auto& [lo, hi]   = rng[sp.unk[0].first];
              double const a   = sp.unk[0].second,
                           r   = (pCut->rhs() - sp.known) / a;
              bool const upper = (pCut->type() == PC::LE) == (a > 0.);
              if (pCut->type() == PC::EQ || upper) hi = std::min(hi, r);
              if (pCut->type() == PC::EQ || !upper) lo = std::max(lo, r);
            }

          for (auto const& pCut : Env.Cuts())
          {
            if (pCut->type() == PC::SOS2) continue;  // holds by construction
            if (pCut->nqvar())
            {
              if (!ig) ++nunk;
              continue;
            }
            Split const sp = split(*pCut);
            double smin = 0., smax = 0.;
            bool unchecked = false, empty = false;
            for (auto const& [id, a] : sp.unk)
            {
              auto itr = rng.find(id);
              if (itr == rng.end())
              {
                unchecked = true;  // neither determined by x nor auxiliary
                break;
              }
              double const lo = itr->second.first, hi = itr->second.second;
              if (lo > hi + 1e-8 * sp.scale)
              {
                record(lo - hi, x);  // no value of the auxiliary fits
                empty = true;
                break;
              }
              smin += std::min(a * lo, a * hi);
              smax += std::max(a * lo, a * hi);
            }
            if (unchecked && !ig) ++nunk;
            if (unchecked || empty) continue;
            double viol = 0.;
            switch (pCut->type())
            {
              case PC::EQ:
                viol = std::max(sp.known + smin - pCut->rhs(),
                                pCut->rhs() - sp.known - smax);
                break;
              case PC::LE:
                viol = sp.known + smin - pCut->rhs();
                break;
              case PC::GE:
                viol = pCut->rhs() - sp.known - smax;
                break;
              default:
                break;
            }
            if (viol > 1e-8 * sp.scale) record(viol, x);
          }
        }
        if (nviol || nunk)
        {
          ++nfail;
          msg << "  " << TNAME[it] << (withbk ? "+bkpt" : "") << ": ";
          if (nviol)
            msg << nviol << " violations, worst " << std::setprecision(3)
                << worst << " at x = " << xworst;
          if (nunk) msg << nunk << " cut(s) not checked";
          msg << ";";
        }
      }
    }
    std::cout << (msg.str().empty() ? "  ok    " : "  FAIL  ") << std::left
              << std::setw(34) << c.name << msg.str() << std::endl;
  }

  std::cout << cases.size() << " cases, " << nrun
            << " relaxations (4 breakpoint types, with and without "
               "breakpoints), "
            << nfail << " failed" << std::endl;
  return nfail;
}
