// =============================================================================
// DAG_pwlu.cpp
// Unit tests for mc::PWLU and the external operation mc::FFPWLU (continuous
// piecewise-linear univariates): values, ReLU form, derivatives (numeric
// forward and reverse modes, symbolic), interval ranges (exact), backward
// contraction, McCormick envelopes, polyhedral cuts, sharing in a DAG.
//
// The process exits with code 0 if every test passes, 1 otherwise.
// =============================================================================

#include <cmath>
#include <iostream>
#include <random>
#include <string>
#include <vector>

#include "interval.hpp"
typedef mc::Interval I;
#include "ffpwlu.hpp"
#include <map>

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
  std::mt19937 gen(1);
  std::uniform_real_distribution<double> U(-1., 1.);

  // Functions: a hat, a ReLU sum with a constant term and two terms with the
  // same breakpoint, a single affine function, a step-like steep function
  mc::PWLU const F0({-1., 0., 2.}, {1., 0., 1.}, 0., 0.5);
  mc::PWLU const F1 = mc::PWLU::relu(0.1, -0.5, {1.2, -0.4, 0.7, 0.3, 0.}, {1., -2., 0., 2., 3.},
                                     {0.3, 0.2, 0.5, 0.6, 1.});
  mc::PWLU const F2({0.}, {0.25}, 1.5, 1.5);
  mc::PWLU const F3({-0.2, -0.1999, 0.4}, {0., 3., -2.}, -1., 0.);
  std::vector<mc::PWLU> const F = {F0, F1, F2, F3};
  auto relu1 = [](double const x) {
    return 0.1 - 0.5 * x + 1.2 * std::max(x + 0.3, 0.) - 0.4 * std::max(-2. * x + 0.2, 0.) +
           0.7 * std::max(0.5, 0.) + 0.3 * std::max(2. * x + 0.6, 0.);
  };

  // PWLU
  {
    bool ok = true;
    for (double x : {-3., -1., -0.5, 0., 0.3, 2., 2.5}) ok = ok && F0.val(x, true) == (x <= -1. ? 1. : x <= 0. ? -x : x <= 2. ? 0.5 * x : 1. + 0.5 * (x - 2.));
    check(ok, "PWLU: values of a hat function, before, at and after its breakpoints, extended beyond its domain");
    double err = 0.;
    for (int k = 0; k < 1000; ++k)
    {
      double const x = 3. * U(gen);
      err = std::max(err, std::fabs(F1.val(x, true) - relu1(x)));
    }
    check(err < 1e-14 && F1.bkpts().size() == 2, "PWLU::relu: same function, breakpoints merged (" +
                                             std::to_string(F1.bkpts().size()) + "), max error " + std::to_string(err));
    check(std::fabs(F1.slope(-5.) - (-0.5 + (-0.4) * (-2.))) < 1e-15 && std::fabs(F1.slope(5.) - (-0.5 + 1.2 + 0.3 * 2.)) < 1e-15,
          "PWLU::relu: end slopes");
    bool thrown = false;
    try
    {
      mc::PWLU({0., 0.}, {1., 2.}, 0., 0.);
    }
    catch (...)
    {
      thrown = true;
    }
    check(thrown, "PWLU: breakpoints not increasing are refused");
  }

  // PWLU: restriction (exact range), envelopes (valid, tight at the vertices), breakpoints
  {
    bool range = true, env = true;
    for (int k = 0; k < 200; ++k)
    {
      double a = 3. * U(gen), b = 3. * U(gen);
      if (a > b) std::swap(a, b);
      for (auto const& f : F)
      {
        mc::PWLU const R = f.restrict(a, b), Hl = R.hull(true), Hu = R.hull(false);
        double lo = INFINITY, up = -INFINITY;
        for (int s = 0; s <= 400; ++s)
        {
          double const t = a + (b - a) * s / 400., ft = f.val(t, true);
          lo = std::min(lo, ft);
          up = std::max(up, ft);
          range = range && std::fabs(R.val(t, true) - ft) <= 1e-10 * (1. + std::fabs(ft));
          env = env && Hl.val(t, true) <= ft + 1e-10 && Hu.val(t, true) >= ft - 1e-10;
        }
        for (double t : f.bkpts())
          if (t > a && t < b)
          {
            lo = std::min(lo, f.val(t, true));
            up = std::max(up, f.val(t, true));
          }
        // up to the rounding of the values, accumulated over the segments (steep for F3)
        range = range && std::fabs(R.l() - lo) <= 1e-10 * (1. + std::fabs(lo)) &&
                std::fabs(R.u() - up) <= 1e-10 * (1. + std::fabs(up));
        // envelopes equal to the function at the ends of the interval and at their own breakpoints
        env = env && std::fabs(Hl.val(a, true) - f.val(a, true)) < 1e-10 && std::fabs(Hu.val(b, true) - f.val(b, true)) < 1e-10;
        for (double t : Hl.bkpts()) env = env && std::fabs(Hl.val(t, true) - f.val(t, true)) < 1e-10;
      }
    }
    check(range, "PWLU::restrict: same function, range equal to the extreme values over the interval");
    check(env, "PWLU::hull: convex and concave envelopes enclose the function, and touch it at their breakpoints");
    check(F0.bkpts() == std::vector<double>({-1., 0.}) && F2.bkpts().empty(),
          "PWLU::bkpts: breakpoints where the slope changes (none at 2 for the hat function, whose slope there stays 1/2)");
    bool thrown = false;
    try
    {
      F0.val(10.);
    }
    catch (...)
    {
      thrown = true;
    }
    check(thrown, "PWLU::val: no extension beyond the domain unless asked");
  }

  // DAG with the four functions of X, and two of them of another operation
  mc::FFGraph DAG;
  mc::FFVar X(&DAG);
  mc::FFPWLU<I> OpF;
  std::vector<mc::FFVar> Y = OpF(X, F);
  mc::FFVar Z = OpF(X, F0);
  std::vector<mc::FFVar> Y2 = mc::FFPWLU<I>()(X, F);
  check(Y2[0].id() == Y[0].id() && Y2[3].id() == Y[3].id() && Z.id() != Y[0].id(),
        "DAG: the same functions of the same operand share an operation");

  mc::FFVar G = sqr(Y[0]) + Y[1] * Y[3] + 2. * Y[2];
  auto gval = [&](double const x) { return F0.val(x, true) * F0.val(x, true) + F1.val(x, true) * F3.val(x, true) + 2. * F2.val(x, true); };

  // Values
  {
    double err = 0.;
    auto sg = DAG.subgraph(Y);
    std::vector<double> wk;
    for (int k = 0; k < 200; ++k)
    {
      double x = 3. * U(gen), y[4];
      DAG.eval(sg, wk, 4, Y.data(), y, 1, &X, &x);
      for (int j = 0; j < 4; ++j) err = std::max(err, std::fabs(y[j] - F[j].val(x, true)));
    }
    check(err == 0., "double: values of the four functions");
  }

  // Derivatives: numeric forward and reverse, symbolic forward and reverse
  {
    auto sgG = DAG.subgraph(1, &G);
    auto dGf = DAG.FAD({G}, {X});
    auto dGb = DAG.BAD({G}, {X});
    auto sgf = DAG.subgraph(dGf), sgb = DAG.subgraph(dGb);
    double errf = 0., errb = 0., errsf = 0., errsb = 0.;
    std::vector<double> wk;
    for (int k = 0; k < 200; ++k)
    {
      double const x = 3. * U(gen), h = 1e-7;
      double const ref = (gval(x + h) - gval(x - h)) / (2. * h);
      std::vector<mc::FADType<double>> Fwk;
      mc::FADType<double> FX = x, FG;
      FX.diff(0, 1);
      DAG.eval(sgG, Fwk, 1, &G, &FG, 1, &X, &FX);
      errf = std::max(errf, std::fabs(FG.d(0) - ref) / (1. + std::fabs(ref)));
      std::vector<mc::BADType<double>> Bwk;
      mc::BADType<double> BX = x, BG;
      DAG.eval(sgG, Bwk, 1, &G, &BG, 1, &X, &BX);
      Bwk.clear();
      BG.diff(0, 1);
      errb = std::max(errb, std::fabs(BX.d(0) - FG.d(0)) / (1. + std::fabs(FG.d(0))));
      double d;
      DAG.eval(sgf, wk, 1, dGf.data(), &d, 1, &X, &x);
      errsf = std::max(errsf, std::fabs(d - FG.d(0)) / (1. + std::fabs(FG.d(0))));
      DAG.eval(sgb, wk, 1, dGb.data(), &d, 1, &X, &x);
      errsb = std::max(errsb, std::fabs(d - FG.d(0)) / (1. + std::fabs(FG.d(0))));
    }
    check(errf < 1e-6, "forward mode: derivative against finite differences, error " + std::to_string(errf));
    check(errb < 1e-13, "reverse mode: same derivative as the forward mode");
    check(errsf < 1e-10 && errsb < 1e-10, "symbolic forward and reverse: same derivative");
    // right derivative at a breakpoint
    std::vector<mc::FADType<double>> Fwk;
    mc::FADType<double> FX = 0., FY;
    FX.diff(0, 1);
    auto sg0 = DAG.subgraph(1, &Y[0]);
    DAG.eval(sg0, Fwk, 1, &Y[0], &FY, 1, &X, &FX);
    check(FY.d(0) == 0.5, "forward mode: right derivative at a breakpoint");
  }

  // Interval: exact ranges on random boxes, and backward contraction
  {
    auto sg = DAG.subgraph(Y);
    std::vector<I> Iwk;
    bool exact = true, valid = true;
    for (int k = 0; k < 200; ++k)
    {
      double a = 3. * U(gen), b = 3. * U(gen);
      if (a > b) std::swap(a, b);
      I IX(a, b), IY[4];
      DAG.eval(sg, Iwk, 4, Y.data(), IY, 1, &X, &IX);
      for (int j = 0; j < 4; ++j)
      {
        double lo = INFINITY, up = -INFINITY;
        std::vector<double> pts = {a, b};
        for (double t : F[j].bkpts())
          if (t > a && t < b) pts.push_back(t);
        for (int s = 0; s <= 1000; ++s) pts.push_back(a + (b - a) * s / 1000.);
        for (double t : pts)
        {
          lo = std::min(lo, F[j].val(t, true));
          up = std::max(up, F[j].val(t, true));
        }
        // extreme values at the vertices, up to the rounding of the values in between
        exact = exact && std::fabs(mc::Op<I>::l(IY[j]) - lo) <= 1e-14 * (1. + std::fabs(lo)) &&
                std::fabs(mc::Op<I>::u(IY[j]) - up) <= 1e-14 * (1. + std::fabs(up));
      }
      (void)valid;
    }
    check(exact, "interval: ranges equal to the extreme values over the box");

    // Contraction of X in [-3,3] from F0(X) in [0,0.25]: X in [-1,0.5]
    mc::FFGraph DAG2;
    mc::FFVar X2(&DAG2);
    mc::FFVar Y0 = mc::FFPWLU<I>()(X2, F0);
    std::vector<I> IX2{I(-3., 3.)}, IY0{I(0., 0.25)};
    int const flag = DAG2.reval(std::vector<mc::FFVar>{Y0}, IY0, std::vector<mc::FFVar>{X2}, IX2, I(-1e20, 1e20));
    check(flag >= 0 && std::fabs(mc::Op<I>::l(IX2[0]) + 0.25) < 1e-12 && std::fabs(mc::Op<I>::u(IX2[0]) - 0.5) < 1e-12,
          "interval: backward contraction to the preimage hull [-0.25,0.5], got [" +
              std::to_string(mc::Op<I>::l(IX2[0])) + "," + std::to_string(mc::Op<I>::u(IX2[0])) + "]");
    std::vector<I> IX3{I(-3., 3.)}, IY3{I(5., 6.)};
    int const flag3 = DAG2.reval(std::vector<mc::FFVar>{Y0}, IY3, std::vector<mc::FFVar>{X2}, IX3, I(-1e20, 1e20));
    check(flag3 < 0, "interval: backward contraction detects an empty preimage");
  }

  // McCormick: envelopes valid, exact at the bounds, subgradients valid
  {
    typedef mc::McCormick<I> MC;
    auto sg = DAG.subgraph(Y);
    std::vector<MC> MCwk;
    bool valid = true, tight = true;
    for (int k = 0; k < 100; ++k)
    {
      double a = 3. * U(gen), b = 3. * U(gen);
      if (a > b) std::swap(a, b);
      double const xr = a + (b - a) * 0.5 * (1. + U(gen));
      MC MX(I(a, b), xr), MY[4];
      MX.sub(1, 0);
      DAG.eval(sg, MCwk, 4, Y.data(), MY, 1, &X, &MX);
      for (int j = 0; j < 4; ++j)
      {
        double const f = F[j].val(xr, true);
        valid = valid && MY[j].cv() <= f + 1e-12 && MY[j].cc() >= f - 1e-12;
        for (int s = 0; s <= 100; ++s)
        {
          double const t = a + (b - a) * s / 100.;
          valid = valid && MY[j].cv() + MY[j].cvsub(0) * (t - xr) <= F[j].val(t, true) + 1e-12 &&
                  MY[j].cc() + MY[j].ccsub(0) * (t - xr) >= F[j].val(t, true) - 1e-12;
        }
      }
      // at a bound of the box, both envelopes equal the function
      MC MA(I(a, b), a), MA2[4];
      MA.sub(1, 0);
      DAG.eval(sg, MCwk, 4, Y.data(), MA2, 1, &X, &MA);
      for (int j = 0; j < 4; ++j)
        tight = tight && std::fabs(MA2[j].cv() - F[j].val(a, true)) < 1e-12 && std::fabs(MA2[j].cc() - F[j].val(a, true)) < 1e-12;
    }
    check(valid, "McCormick: envelopes and their subgradient cuts enclose the functions on the box");
    check(tight, "McCormick: envelopes equal the functions at the bounds of the box");
  }

  // Chebyshev models: valid enclosures; expressions; dependencies
  {
    mc::SCModel<I> CM(4);
    mc::SCVar<I> CX(&CM, 0, I(-1.5, 2.5)), CY[4];
    auto sg = DAG.subgraph(Y);
    std::vector<mc::SCVar<I>> CMwk;
    DAG.eval(sg, CMwk, 4, Y.data(), CY, 1, &X, &CX);
    bool valid = true;
    for (int j = 0; j < 4; ++j)
      for (int s = 0; s <= 400; ++s)
      {
        double const t = -1.5 + 4. * s / 400.;
        valid = valid && mc::Op<I>::l(CY[j].B()) <= F[j].val(t, true) + 1e-12 && mc::Op<I>::u(CY[j].B()) >= F[j].val(t, true) - 1e-12;
      }
    check(valid, "Chebyshev models: bounds enclose the functions");

    auto expr = mc::FFExpr::subgraph(&DAG, sg);
    std::ostringstream os;
    os << expr[1];
    check(expr.size() == 4 && os.str().find("PWLU[") != std::string::npos, "expressions: " + os.str());

    std::vector<mc::FFDep> Dwk;
    mc::FFDep DX = 0, DY[4];
    DX.indep(0);
    DAG.eval(sg, Dwk, 4, Y.data(), DY, 1, &X, &DX);
    check(DY[0].dep().at(0) == mc::FFDep::TYPE::N && DY[2].dep().at(0) == mc::FFDep::TYPE::L,
          "dependencies: nonlinear, and linear for an affine function");
  }

  // PWLU::compose: exact composite with a piecewise-linear estimator
  {
    mc::PWLU const u(-1., 0.5, {0.6, 0.9, 1.}, {2., -1.5, 0.4});  // on [-1,1.5]
    bool ok = true;
    for (auto const& f : F)
    {
      mc::PWLU const c = mc::PWLU(u).compose(f);
      for (int s = 0; s <= 500; ++s)
      {
        double const t = -1. + 2.5 * s / 500.;
        ok = ok && std::fabs(c.val(t, true) - f.val(u.val(t, true), true)) <= 1e-10 * (1. + std::fabs(f.val(u.val(t, true), true)));
      }
    }
    check(ok, "PWLU::compose: equal to the composite, with the crossings of the breakpoints added");
  }

  // Superposition models: exact for a variable or an affine function of it (PWLU), exact range on
  // each cell (PWCU); valid for an operand depending on two variables
  {
    double gapL = 0., cellC = 0., viol = 0.;
    for (int aff = 0; aff < 2; ++aff)
    {
      mc::FFGraph DAG4;
      mc::FFVar X4(&DAG4), Z4(&DAG4);
      mc::FFVar A4 = (aff ? 0.8 * X4 + 0.3 : X4);
      std::vector<mc::FFVar> Y4 = mc::FFPWLU<I>()(A4, F);
      mc::SupModel<mc::PWLU> modL(2);
      mc::SupVar<mc::PWLU> SXL(modL, 0, I(-1.5, 2.5), 8), SYL[4];
      DAG4.eval(4, Y4.data(), SYL, 1, &X4, &SXL);
      mc::SupModel<mc::PWCU> modC(2);
      mc::SupVar<mc::PWCU> SXC(modC, 0, I(-1.5, 2.5), 16), SYC[4];
      DAG4.eval(4, Y4.data(), SYC, 1, &X4, &SXC);
      for (int j = 0; j < 4; ++j)
      {
        for (int s = 0; s <= 400; ++s)
        {
          double const x = -1.5 + 4. * s / 400., f = F[j].val(aff ? 0.8 * x + 0.3 : x, true);
          std::map<unsigned, double> const pt{{0, x}};
          gapL = std::max({gapL, std::fabs(SYL[j].uval(pt) - f), std::fabs(SYL[j].oval(pt) - f)});
          viol = std::max({viol, SYC[j].uval(pt) - f, f - SYC[j].oval(pt)});
        }
        for (int k = 0; k < 16; ++k)
        {
          double const a = -1.5 + 0.25 * k, b = a + 0.25;
          mc::PWLU const R = F[j].restrict(aff ? 0.8 * a + 0.3 : a, aff ? 0.8 * b + 0.3 : b);
          cellC = std::max({cellC, std::fabs(SYC[j].uest()[0].yL()[k] - R.l()),
                            std::fabs(SYC[j].oest()[0].yU()[k] - R.u())});
        }
      }
    }
    check(gapL < 1e-12, "superposition (PWLU): estimators equal to the function of a variable or of an affine function of it");
    check(cellC < 1e-12 && viol < 1e-12, "superposition (PWCU): exact range on each cell, estimators valid");

    // with MAX_SUBDIV below the number of segments of the composite: reduced, still valid
    {
      mc::FFGraph DAG6;
      mc::FFVar X6(&DAG6);
      mc::FFVar Y6 = mc::FFPWLU<I>()(X6, F[3]);
      mc::SupModel<mc::PWLU> mod6(1);
      mod6.options.MAX_SUBDIV = 2;
      mc::SupVar<mc::PWLU> SX6(mod6, 0, I(-1.5, 2.5), 8), SY6;
      DAG6.eval(1, &Y6, &SY6, 1, &X6, &SX6);
      double viol6 = 0., gap6 = 0.;
      for (int s = 0; s <= 400; ++s)
      {
        double const x = -1.5 + 4. * s / 400., f = F[3].val(x, true);
        std::map<unsigned, double> const pt{{0, x}};
        viol6 = std::max({viol6, SY6.uval(pt) - f, f - SY6.oval(pt)});
        gap6  = std::max({gap6, f - SY6.uval(pt), SY6.oval(pt) - f});
      }
      check(viol6 < 1e-10 && gap6 > 0. && SY6.uest()[0].size() <= 2 && SY6.oest()[0].size() <= 2,
            "superposition (PWLU): reduced to MAX_SUBDIV = 2 segments, estimators valid");
    }

    mc::FFGraph DAG5;
    mc::FFVar X5(&DAG5), Z5(&DAG5);
    std::vector<mc::FFVar> Y5 = mc::FFPWLU<I>()(0.6 * X5 + 0.4 * Z5, F);
    mc::SupModel<mc::PWLU> mod5(2);
    mc::SupVar<mc::PWLU> SV5[2] = {mc::SupVar<mc::PWLU>(mod5, 0, I(-1.5, 2.5), 8),
                                   mc::SupVar<mc::PWLU>(mod5, 1, I(-1., 1.), 8)}, SY5[4];
    mc::FFVar V5[2] = {X5, Z5};
    DAG5.eval(4, Y5.data(), SY5, 2, V5, SV5);
    double viol5 = 0.;
    for (int j = 0; j < 4; ++j)
      for (int a = 0; a <= 40; ++a)
        for (int b = 0; b <= 20; ++b)
        {
          double const x = -1.5 + 4. * a / 40., z = -1. + 2. * b / 20., f = F[j].val(0.6 * x + 0.4 * z, true);
          std::map<unsigned, double> const pt{{0, x}, {1, z}};
          viol5 = std::max({viol5, SY5[j].uval(pt) - f, f - SY5[j].oval(pt)});
        }
    check(viol5 < 1e-12, "superposition (PWLU): estimators valid for an operand depending on two variables");
  }

  // Polyhedral image: cuts generated for each encoding; HULL cuts valid on the graph
  {
    for (int enc : {mc::FFPWLU<I>::Options::HULL, mc::FFPWLU<I>::Options::CONT,
                    mc::FFPWLU<I>::Options::BIN, mc::FFPWLU<I>::Options::SOS2})
      for (bool sep : {false, true})
      {
        mc::FFGraph DAG3;
        mc::FFVar X3(&DAG3);
        mc::FFPWLU<I> Op3;
        Op3.options.ENCODING = enc;
        Op3.options.SEPARATE = sep;
        std::vector<mc::FFVar> Y3 = Op3(X3, F);
        mc::PolImg<I> Env;
        mc::PolVar<I> PX(&Env, X3, I(-1.5, 2.5)), PY[4];
        DAG3.eval(4, Y3.data(), PY, 1, &X3, &PX);
        Env.generate_cuts(4, PY);
        size_t ncut = Env.Cuts().size(), naux = Env.Aux().size();
        bool ok = ncut > 0;
        if (enc == mc::FFPWLU<I>::Options::HULL)
        {
          ok = ok && naux == 0;
          // every cut holds on the graph of the functions
          for (int s = 0; s <= 400; ++s)
          {
            double const t = -1.5 + 4. * s / 400.;
            for (auto const& cut : Env.Cuts())
            {
              double lhs = 0.;
              for (unsigned i = 0; i < cut->nvar(); ++i)
              {
                auto const& v = cut->var()[i];
                double val = (v.id() == PX.id() ? t : NAN);
                for (int j = 0; j < 4; ++j)
                  if (v.id() == PY[j].id()) val = F[j].val(t, true);
                lhs += cut->coef()[i] * val;
              }
              double const rhs = cut->rhs();
              switch (cut->type())
              {
                case mc::PolCut<I>::EQ: ok = ok && std::fabs(lhs - rhs) < 1e-12; break;
                case mc::PolCut<I>::LE: ok = ok && lhs <= rhs + 1e-12; break;
                case mc::PolCut<I>::GE: ok = ok && lhs >= rhs - 1e-12; break;
                default: break;
              }
            }
          }
        }
        else
          ok = ok && naux > 0;
        check(ok, "polyhedral image, encoding " + std::to_string(enc) + (sep ? ", separate" : ", shared") +
                      ": " + std::to_string(ncut) + " cuts, " + std::to_string(naux) + " auxiliaries");
      }
  }

  std::cout << npass << " passed, " << nfail << " failed" << std::endl;
  return nfail ? 1 : 0;
}
