// Constraint propagation (FFGraph::reval) when an operation returns a bound
// that is not a number in the arithmetic used: such a result carries no
// information, as when the operation throws. The propagation must not report
// infeasibility for it, and must leave the ranges it cannot justify unchanged.
// The program returns the number of failed checks.
#undef MC__REVAL_DEBUG

#include <cmath>
#include <iomanip>
#include <limits>
#include <string>
#include <vector>

#include "fflin.hpp"
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

static unsigned nfail = 0;

static void
check(bool const ok, std::string const& what)
{
  if (!ok) ++nfail;
  std::cout << (ok ? "  ok     " : "  FAIL   ") << what << std::endl;
}

static double
lo(I const& x)
{
  return mc::Op<I>::l(x);
}

static double
up(I const& x)
{
  return mc::Op<I>::u(x);
}

int
main()
{
  using mc::FFVar;
  const double INF = 1e30, NaN = std::numeric_limits<double>::quiet_NaN();
  const I IINF(-INF, INF);

  // A range whose bounds are not numbers, as left by an operation that failed
  const I INaN(NaN, NaN);
  check(mc::FFOp::_has_nan(INaN) && !mc::FFOp::_has_nan(I(1., 2.)) &&
            !mc::FFOp::_has_nan(IINF),
        "_has_nan: bounds that are not numbers, against proper ranges");

  // Intersection: a range with such bounds is ignored
  {
    I r(0.);
    bool ok = mc::FFOp::_inter(r, I(1., 3.), INaN);
    check(ok && lo(r) == 1. && up(r) == 3., "_inter: second argument ignored");
    r  = I(0.);
    ok = mc::FFOp::_inter(r, INaN, I(1., 3.));
    check(ok && lo(r) == 1. && up(r) == 3., "_inter: first argument ignored");
    r  = I(1., 3.);
    ok = mc::FFOp::_inter(r, INaN, r);
    check(ok && lo(r) == 1. && up(r) == 3., "_inter: in place");
    ok = mc::FFOp::_inter(r, I(1., 3.), I(2., 5.));
    check(ok && lo(r) == 2. && up(r) == 3., "_inter: two proper ranges");
    check(!mc::FFOp::_inter(r, I(1., 2.), I(3., 5.)),
          "_inter: two disjoint ranges are infeasible");
  }

  // F = X + c with c not a number: nothing is known of F, X is left as it is
  {
    mc::FFGraph DAG;
    FFVar X(&DAG), C(&DAG, NaN);
    FFVar F  = X + C;
    auto sgF = DAG.subgraph(1, &F);
    std::vector<I> wk;
    I IX(2., 3.), IF(0., 1.);
    int const flag =
        DAG.reval(sgF, wk, 1, &F, &IF, 1, &X, &IX, IINF, 20, 1e-10);
    check(flag >= 0, "X + NaN in [0,1]: not reported infeasible");
    check(lo(IX) == 2. && up(IX) == 3. && lo(IF) == 0. && up(IF) == 1.,
          "X + NaN in [0,1]: ranges of X and F unchanged");
  }

  // F = exp(X) / exp(Z) with X and Z in [700,720]: both terms reach the
  // largest number, and infinity over infinity, or times zero, has no bounds
  // in some arithmetics; F = 1 at X = Z
  {
    mc::FFGraph DAG;
    FFVar X(&DAG), Z(&DAG);
    FFVar F  = exp(X) / exp(Z);
    auto sgF = DAG.subgraph(1, &F);
    std::vector<I> wk;
    FFVar V[2] = {X, Z};
    I IV[2]    = {I(700., 720.), I(700., 720.)}, IF(0.5, 2.);
    int flag   = -1;
    try
    {
      flag = DAG.reval(sgF, wk, 1, &F, &IF, 2, V, IV, IINF, 20, 1e-10);
    }
    catch (...)
    {
    }
    check(flag >= 0, "exp(X)/exp(Z) in [0.5,2], overflow: not infeasible");
    check(lo(IV[0]) <= 700. && up(IV[0]) >= 720. && lo(IV[1]) <= 700. &&
              up(IV[1]) >= 720.,
          "exp(X)/exp(Z) in [0.5,2], overflow: X = Z in [700,720] all kept");
  }

  // Linear sum X0 + 0*X1 + 1 with X1 unbounded: zero times infinity has no
  // bounds in some arithmetics; X0 = F - 1 is in [0.2,0.4]
  {
    mc::FFGraph DAG;
    std::vector<FFVar> X = DAG.add_vars(2);
    mc::FFLin<I> Sum;
    FFVar F  = Sum(X, std::vector<double>({1., 0.}), 1.);
    auto sgF = DAG.subgraph(1, &F);
    std::vector<I> wk;
    const double inf = std::numeric_limits<double>::infinity();
    std::vector<I> IX({I(0., 1.), I(-inf, inf)});
    I IF(1.2, 1.4);
    int flag = -1;
    try
    {
      flag = DAG.reval(sgF, wk, 1, &F, &IF, 2, X.data(), IX.data(), IINF, 20,
                       1e-10);
    }
    catch (...)
    {
    }
    check(flag >= 0,
          "X0 + 0*X1 + 1 in [1.2,1.4], X1 unbounded: not infeasible");
    check(lo(IX[0]) <= 0.2 + 1e-12 && up(IX[0]) >= 0.4 - 1e-12 &&
              lo(IX[0]) >= 0. && up(IX[0]) <= 1.,
          "X0 + 0*X1 + 1 in [1.2,1.4], X1 unbounded: X0 keeps [0.2,0.4]");
  }

  std::cout << nfail << " failed" << std::endl;
  return nfail;
}
