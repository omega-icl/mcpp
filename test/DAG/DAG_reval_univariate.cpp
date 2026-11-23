// Backward propagation (FFGraph::reval) through univariate operations, and
// through a sum of periodic terms: after propagation, the range of the operand
// must contain every x of its initial range whose image lies in the range set
// for the result. The admissible x are sampled on a fine grid. Where the
// operand stays on one branch of the operation (one monotonic piece), the
// propagated range must also be the hull of the admissible x, up to the grid
// spacing; an empty set of admissible x must be reported infeasible.
// Covers the inversion of sin, cos and tan outside their principal branch and
// that of fstep, and checks the other univariate operations for regressions.
// The program returns the number of failed cases.
#undef MC__REVAL_DEBUG

#include <cmath>
#include <functional>
#include <iomanip>
#include <string>
#include <vector>

#include "ffunc.hpp"

#ifdef MC__USE_PROFIL
#include "mcprofil.hpp"
typedef INTERVAL I;
#else
#ifdef MC__USE_FILIB
#include "mcfilib.hpp"
typedef filib::interval<double, filib::native_switched, filib::i_mode_extended>
    I;
#else
#ifdef MC__USE_BOOST
#include "mcboost.hpp"
typedef boost::numeric::interval_lib::save_state<
    boost::numeric::interval_lib::rounded_transc_opp<double> >
    T_boost_round;
typedef boost::numeric::interval_lib::checking_base<double> T_boost_check;
typedef boost::numeric::interval_lib::policies<T_boost_round, T_boost_check>
    T_boost_policy;
typedef boost::numeric::interval<double, T_boost_policy> I;
#else
#include "interval.hpp"
typedef mc::Interval I;
#endif
#endif
#endif

struct Case
{
  std::string name;
  std::function<mc::FFVar(mc::FFVar const&)> op;  // expression in X
  std::function<double(double)> f;                // the same in double
  double xl, xu;                                  // initial operand range
  double fl, fu;                                  // range set for the result
  bool tight;  // operand range expected to be the hull of the admissible x
};

int
main()
{
  using mc::FFVar;
  const double INF = 1e30;

  // clang-format off
  const std::vector<Case> cases = {
    // sin: operand range across several branches (no reduction expected)
    { "sin, X=[1,5]",                 [](FFVar const& x){ return sin(x); }, [](double x){ return std::sin(x); }, 1., 5., -INF, INF, true },
    { "sin, X=[9,12]",                [](FFVar const& x){ return sin(x); }, [](double x){ return std::sin(x); }, 9., 12., -INF, INF, true },
    { "sin, X=[0,25]",                [](FFVar const& x){ return sin(x); }, [](double x){ return std::sin(x); }, 0., 25., -INF, INF, true },
    { "sin, X=[1,5], F=[0,0.5]",      [](FFVar const& x){ return sin(x); }, [](double x){ return std::sin(x); }, 1., 5., 0., 0.5, false },
    // sin: operand range on one branch
    { "sin, X=[-1,1], F=[0,0.5]",     [](FFVar const& x){ return sin(x); }, [](double x){ return std::sin(x); }, -1., 1., 0., 0.5, true },
    { "sin, X=[2,4], F=[0,0.5]",      [](FFVar const& x){ return sin(x); }, [](double x){ return std::sin(x); }, 2., 4., 0., 0.5, true },
    { "sin, X=[8,10.5], F=[-0.5,0]",  [](FFVar const& x){ return sin(x); }, [](double x){ return std::sin(x); }, 8., 10.5, -0.5, 0., true },
    { "sin, X=[-4,-2], F=[0,0.5]",    [](FFVar const& x){ return sin(x); }, [](double x){ return std::sin(x); }, -4., -2., 0., 0.5, true },
    { "sin, X=[1.6,4.7], F=[-.5,.5]", [](FFVar const& x){ return sin(x); }, [](double x){ return std::sin(x); }, 1.6, 4.7, -0.5, 0.5, true },
    { "sin, X=[-7.8,-4.8], F=[.2,.9]",[](FFVar const& x){ return sin(x); }, [](double x){ return std::sin(x); }, -7.8, -4.8, 0.2, 0.9, true },
    { "sin, X=[2,3], F=[0.95,1]",     [](FFVar const& x){ return sin(x); }, [](double x){ return std::sin(x); }, 2., 3., 0.95, 1., true },
    // cos
    { "cos, X=[1,5]",                 [](FFVar const& x){ return cos(x); }, [](double x){ return std::cos(x); }, 1., 5., -INF, INF, true },
    { "cos, X=[9,25]",                [](FFVar const& x){ return cos(x); }, [](double x){ return std::cos(x); }, 9., 25., -INF, INF, true },
    { "cos, X=[0.5,3], F=[0,0.5]",    [](FFVar const& x){ return cos(x); }, [](double x){ return std::cos(x); }, 0.5, 3., 0., 0.5, true },
    { "cos, X=[3.5,6], F=[0,0.5]",    [](FFVar const& x){ return cos(x); }, [](double x){ return std::cos(x); }, 3.5, 6., 0., 0.5, true },
    { "cos, X=[-3,-0.5], F=[0,0.5]",  [](FFVar const& x){ return cos(x); }, [](double x){ return std::cos(x); }, -3., -0.5, 0., 0.5, true },
    { "cos, X=[-6,-3.5], F=[0,0.5]",  [](FFVar const& x){ return cos(x); }, [](double x){ return std::cos(x); }, -6., -3.5, 0., 0.5, true },
    { "cos, X=[6.3,9.4], F=[-.3,.6]", [](FFVar const& x){ return cos(x); }, [](double x){ return std::cos(x); }, 6.3, 9.4, -0.3, 0.6, true },
    { "cos, X=[-9.4,-6.3], F=[-.3,.6]",[](FFVar const& x){ return cos(x); }, [](double x){ return std::cos(x); }, -9.4, -6.3, -0.3, 0.6, true },
    { "cos, X=[4,5], F=[0.5,1]",      [](FFVar const& x){ return cos(x); }, [](double x){ return std::cos(x); }, 4., 5., 0.5, 1., true },
    // tan
    { "tan, X=[-1,1], F=[0,1]",       [](FFVar const& x){ return tan(x); }, [](double x){ return std::tan(x); }, -1., 1., 0., 1., true },
    { "tan, X=[2.5,4], F=[0,1]",      [](FFVar const& x){ return tan(x); }, [](double x){ return std::tan(x); }, 2.5, 4., 0., 1., true },
    { "tan, X=[-4,-2.5], F=[0,1]",    [](FFVar const& x){ return tan(x); }, [](double x){ return std::tan(x); }, -4., -2.5, 0., 1., true },
    { "tan, X=[5,7.5], F=[-2,3]",     [](FFVar const& x){ return tan(x); }, [](double x){ return std::tan(x); }, 5., 7.5, -2., 3., true },
    // fstep
    { "fstep, X=[2,5]",               [](FFVar const& x){ return fstep(x); }, [](double x){ return x >= 0. ? 1. : 0.; }, 2., 5., -INF, INF, true },
    { "fstep, X=[-5,-2]",             [](FFVar const& x){ return fstep(x); }, [](double x){ return x >= 0. ? 1. : 0.; }, -5., -2., -INF, INF, true },
    { "fstep, X=[-2,5]",              [](FFVar const& x){ return fstep(x); }, [](double x){ return x >= 0. ? 1. : 0.; }, -2., 5., -INF, INF, true },
    { "fstep, X=[-2,5], F=[1,1]",     [](FFVar const& x){ return fstep(x); }, [](double x){ return x >= 0. ? 1. : 0.; }, -2., 5., 1., 1., true },
    { "fstep, X=[-2,5], F=[0,0]",     [](FFVar const& x){ return fstep(x); }, [](double x){ return x >= 0. ? 1. : 0.; }, -2., 5., 0., 0., true },
    { "fstep, X=[2,5], F=[0,0]",      [](FFVar const& x){ return fstep(x); }, [](double x){ return x >= 0. ? 1. : 0.; }, 2., 5., 0., 0., true },
    // sums of periodic terms
    { "sin(x)+sin(10x/3), X=[2.7,7.5], F<=-1.8", [](FFVar const& x){ return sin(x) + sin(10. * x / 3.); },
      [](double x){ return std::sin(x) + std::sin(10. * x / 3.); }, 2.7, 7.5, -INF, -1.8, false },
    { "cos(x)+cos(3x), X=[-5,9], F>=1.5",        [](FFVar const& x){ return cos(x) + cos(3. * x); },
      [](double x){ return std::cos(x) + std::cos(3. * x); }, -5., 9., 1.5, INF, false },
    // other operations: no regression
    { "pow(x,-3), X=[1,4]",           [](FFVar const& x){ return pow(x, -3); }, [](double x){ return std::pow(x, -3); }, 1., 4., -INF, INF, true },
    { "pow(x,-2), X=[1,4], F=[.1,.5]",[](FFVar const& x){ return pow(x, -2); }, [](double x){ return std::pow(x, -2); }, 1., 4., 0.1, 0.5, true },
    { "pow(x,1.5), X=[1,4], F=[2,5]", [](FFVar const& x){ return pow(x, 1.5); }, [](double x){ return std::pow(x, 1.5); }, 1., 4., 2., 5., true },
    { "pow(x,3), X=[-2,3], F=[-1,2]", [](FFVar const& x){ return pow(x, 3); }, [](double x){ return x * x * x; }, -2., 3., -1., 2., false },
    { "pow(x,4), X=[-2,3], F=[1,2]",  [](FFVar const& x){ return pow(x, 4); }, [](double x){ return x * x * x * x; }, -2., 3., 1., 2., true },
    { "sqr, X=[-2,3], F=[1,2]",       [](FFVar const& x){ return sqr(x); }, [](double x){ return x * x; }, -2., 3., 1., 2., true },
    { "fabs, X=[-2,3], F=[1,2]",      [](FFVar const& x){ return fabs(x); }, [](double x){ return std::fabs(x); }, -2., 3., 1., 2., true },
    { "cosh, X=[-2,3], F=[1,2]",      [](FFVar const& x){ return cosh(x); }, [](double x){ return std::cosh(x); }, -2., 3., 1., 2., false },
    { "exp, X=[-2,3], F=[1,2]",       [](FFVar const& x){ return exp(x); }, [](double x){ return std::exp(x); }, -2., 3., 1., 2., true },
    { "log, X=[0.5,4], F=[0,1]",      [](FFVar const& x){ return log(x); }, [](double x){ return std::log(x); }, 0.5, 4., 0., 1., true },
    { "sqrt, X=[0.5,4], F=[1,1.5]",   [](FFVar const& x){ return sqrt(x); }, [](double x){ return std::sqrt(x); }, 0.5, 4., 1., 1.5, true },
    { "asin, X=[-1,1], F=[0,1]",      [](FFVar const& x){ return asin(x); }, [](double x){ return std::asin(x); }, -1., 1., 0., 1., true },
    { "acos, X=[-1,1], F=[0,1]",      [](FFVar const& x){ return acos(x); }, [](double x){ return std::acos(x); }, -1., 1., 0., 1., true },
    { "atan, X=[-3,3], F=[0,1]",      [](FFVar const& x){ return atan(x); }, [](double x){ return std::atan(x); }, -3., 3., 0., 1., true },
  };
  // clang-format on

  const unsigned NSAM = 400000;
  unsigned nfail      = 0;
  std::cout << std::setprecision(6);
  for (auto const& c : cases)
  {
    mc::FFGraph DAG;
    FFVar X(&DAG);
    FFVar F  = c.op(X);
    auto sgF = DAG.subgraph(1, &F);

    std::vector<I> wk;
    I IX(c.xl, c.xu), IF(c.fl, c.fu);
    int flag = 0;
    bool exc = false;
    try
    {
      flag =
          DAG.reval(sgF, wk, 1, &F, &IF, 1, &X, &IX, I(-INF, INF), 20, 1e-10);
    }
    catch (...)
    {
      exc = true;
    }

    // Hull of the admissible x, sampled
    double tl = INF, tu = -INF;
    for (unsigned k = 0; k <= NSAM; k++)
    {
      const double x = c.xl + (c.xu - c.xl) * k / NSAM, f = c.f(x);
      if (f < c.fl || f > c.fu) continue;
      tl = std::min(tl, x);
      tu = std::max(tu, x);
    }
    const bool empty  = (tl > tu);
    const double gtol = 2. * (c.xu - c.xl) / NSAM + 1e-9;

    bool ok = !exc;
    if (ok && flag < 0)
      ok = empty;
    else if (ok && !empty)
      ok = (mc::Op<I>::l(IX) <= tl + 1e-9 && mc::Op<I>::u(IX) >= tu - 1e-9);
    bool tight = true;
    if (ok && c.tight && flag >= 0 && !empty)
      tight = (mc::Op<I>::l(IX) >= tl - gtol && mc::Op<I>::u(IX) <= tu + gtol);
    if (ok && c.tight && flag >= 0 && empty) tight = false;
    if (!ok || !tight) ++nfail;

    std::cout << (!ok      ? "  FAIL   "
                  : !tight ? "  LOOSE  "
                           : "  ok     ")
              << std::left << std::setw(42) << c.name << std::right;
    if (exc)
      std::cout << "exception";
    else if (flag < 0)
      std::cout << "infeasible";
    else
      std::cout << "X -> [" << mc::Op<I>::l(IX) << "," << mc::Op<I>::u(IX)
                << "]";
    if (empty)
      std::cout << "   no admissible x";
    else
      std::cout << "   admissible x in [" << tl << "," << tu << "]";
    std::cout << std::endl;
  }

  std::cout << cases.size() << " cases, " << nfail << " failed" << std::endl;
  return nfail;
}
