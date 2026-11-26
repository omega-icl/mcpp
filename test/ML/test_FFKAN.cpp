// Copyright (C) Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the EPL-2.0 with GPL-2.0-or-later as a Secondary License; see the LICENSE file.
// SPDX-License-Identifier: EPL-2.0 OR GPL-2.0-or-later

// A Kolmogorov-Arnold network as an external operation of a DAG: mc::FFKAN, which is mc::FFMAP on
// an mc::KAN, as mc::FFMLP is mc::FFMAP on an mc::MLP. The operation is checked against the network
// itself (mc::KAN, see test_KAN.cpp): values, interval and McCormick bounds, the Jacobian, tangent
// and adjoint operations that the differentiation of the DAG inserts, the polyhedral relaxation
// with each encoding of the edges, the backward propagation of bounds, and the network that the
// DAG holds under each policy. The last checks are on an mc::MLP through the same class.
// The program returns the number of failed checks.

#include <cmath>
#include <iostream>
#include <random>
#include <sstream>
#include <string>
#include <type_traits>
#include <vector>

#include "ffkan.hpp"
#include "ffmlp.hpp"
#include "interval.hpp"
#include "mccormick.hpp"
#include "polimage.hpp"

typedef mc::Interval I;
typedef mc::McCormick<I> MC;
typedef mc::PolVar<I> PV;
using mc::FFGraph;
using mc::FFVar;
using mc::PWLU;

static_assert(std::is_same<mc::FFKAN<I>, mc::FFMAP<I, mc::KAN<I>>>::value, "FFKAN is FFMAP on a KAN");
static_assert(std::is_same<mc::FFMLP<I>, mc::FFMAP<I, mc::MLP<I>>>::value, "FFMLP is FFMAP on an MLP");
static_assert(std::is_same<mc::FFGradKAN<I>, mc::FFGradMAP<I, mc::KAN<I>>>::value, "FFGradKAN");
static_assert(std::is_same<mc::FFTanMLP<>, mc::FFTanMAP<I, mc::MLP<I>>>::value, "FFTanMLP");

static int npass = 0, nfail = 0;
static void
check(bool const ok, std::string const& what)
{
  std::cout << (ok ? "  ok    " : "  FAIL  ") << what << std::endl;
  (ok ? npass : nfail)++;
}

// Random piecewise-linear edge with n breakpoints in [-1,1]
static PWLU
edge(std::mt19937& gen, size_t const n)
{
  std::uniform_real_distribution<double> U(-1., 1.);
  std::vector<double> a(n), w(n), b(n);
  for (size_t r = 0; r < n; ++r)
  {
    a[r] = U(gen);
    w[r] = (U(gen) > 0. ? 1. : -1.) * (0.5 + std::fabs(U(gen)));
    b[r] = -w[r] * U(gen);
  }
  return PWLU::relu(0.3 * U(gen), U(gen), a, w, b);
}

static std::vector<std::vector<PWLU>>
layer(std::mt19937& gen, size_t const nout, size_t const nin, size_t const n)
{
  std::vector<std::vector<PWLU>> L(nout);
  for (auto& row : L)
    for (size_t i = 0; i < nin; ++i) row.push_back(edge(gen, n));
  return L;
}

static std::vector<double>
evaluate(FFGraph& DAG, std::vector<FFVar> const& F, std::vector<FFVar> const& X, std::vector<double> const& x)
{
  std::vector<double> f(F.size());
  if (!F.empty()) DAG.eval(F, f, X, x);
  return f;
}

// names of the external operations that F depends on
static std::string
operations(FFGraph& DAG, std::vector<FFVar> const& F)
{
  std::string names;
  auto sg = DAG.subgraph(F);
  for (auto const& op : sg.l_op)
  {
    if (op->type < mc::FFOp::EXTERN) continue;
    std::string const n = op->name();
    names += n.substr(0, n.find('[')) + " ";
  }
  return names;
}

static bool
has(std::string const& names, std::string const& name)
{
  return (" " + names).find(" " + name + " ") != std::string::npos;
}

static double
dist(std::vector<double> const& a, std::vector<double> const& b)
{
  if (a.size() != b.size()) return INFINITY;
  double d = 0.;
  for (size_t i = 0; i < a.size(); ++i) d = std::max(d, std::fabs(a[i] - b[i]) / (1. + std::fabs(b[i])));
  return d;
}

int
main()
{
  std::mt19937 gen(7);
  size_t const nin = 2, nout = 2;
  mc::KAN<I> NN;
  NN.append_data(layer(gen, 3, nin, 3));
  NN.append_data(layer(gen, nout, 3, 2));
  std::vector<std::vector<double>> const pts = {{0.3, -0.2}, {-0.7, 0.9}, {0.05, 0.45}, {0.8, -0.85}};
  std::vector<I> const IX = {I(-1., 1.), I(-0.5, 1.)};

  std::cout << "Kolmogorov-Arnold network as an external operation" << std::endl;
  FFGraph DAG;
  std::vector<FFVar> X(nin);
  for (auto& Xi : X) Xi.set(&DAG);
  mc::FFKAN<I> Op;
  std::vector<FFVar> const Y = Op(X, &NN);
  {
    std::ostringstream os;
    os << "inserted: " << Y.size() << " results, operations " << operations(DAG, Y);
    check(Y.size() == nout && operations(DAG, Y) == "KAN ", os.str());
  }

  // Values, alone and inside an expression
  {
    std::vector<FFVar> const F = {Y[0] * X[1] + exp(Y[1])};
    double err = 0., errf = 0.;
    for (auto const& x : pts)
    {
      std::vector<double> y(nout);
      NN.eval(x.data(), y.data());
      err  = std::max(err, dist(evaluate(DAG, Y, X, x), y));
      errf = std::max(errf, dist(evaluate(DAG, F, X, x), {y[0] * x[1] + std::exp(y[1])}));
    }
    std::ostringstream os;
    os << "values: error " << err << ", inside an expression " << errf;
    check(err < 1e-14 && errf < 1e-14, os.str());
  }

  // Interval and McCormick bounds: those of the network
  {
    std::vector<I> IY(nout), IYref(nout);
    DAG.eval(Y, IY, X, IX);
    NN.eval(IX.data(), IYref.data());
    bool same = true;
    for (size_t k = 0; k < nout; ++k)
      same = same && mc::Op<I>::l(IY[k]) == mc::Op<I>::l(IYref[k]) && mc::Op<I>::u(IY[k]) == mc::Op<I>::u(IYref[k]);
    check(same, "interval bounds: those of the network");

    std::vector<double> const xref = {0.2, 0.1};
    std::vector<MC> MX(nin), MY(nout), MYref(nout);
    for (size_t i = 0; i < nin; ++i) MX[i] = MC(IX[i], xref[i]).sub(nin, i);
    DAG.eval(Y, MY, X, MX);
    NN.eval(MX.data(), MYref.data());
    std::vector<double> y(nout);
    NN.eval(xref.data(), y.data());
    bool ok = true;
    for (size_t k = 0; k < nout; ++k)
      ok = ok && MY[k].cv() == MYref[k].cv() && MY[k].cc() == MYref[k].cc() && MY[k].cv() <= y[k] + 1e-12 &&
           y[k] <= MY[k].cc() + 1e-12 && MY[k].cvsub(0) == MYref[k].cvsub(0) && MY[k].ccsub(1) == MYref[k].ccsub(1);
    check(ok, "McCormick relaxations: those of the network, enclosing its value");
  }

  // Derivatives: the Jacobian operation, then the tangent and adjoint operations with the options
  // of the network that the DAG holds
  {
    auto jacerr = [&](std::vector<FFVar> const& dY) {
      double err = 0.;
      for (auto const& x : pts)
      {
        std::vector<double> J(nout * nin);
        NN.grad(x.data(), J.data());
        err = std::max(err, dist(evaluate(DAG, dY, X, x), J));
      }
      return err;
    };
    std::vector<FFVar> const dY = DAG.FAD(Y, X), bY = DAG.BAD(Y, X);
    std::ostringstream os;
    os << "Jacobian: forward " << jacerr(dY) << " (" << operations(DAG, dY) << "), backward " << jacerr(bY) << " ("
       << operations(DAG, bY) << ")";
    check(jacerr(dY) < 1e-13 && jacerr(bY) < 1e-13 && has(operations(DAG, dY), "GradKAN") &&
              has(operations(DAG, bY), "GradKAN"),
          os.str());

    Op.map().options.TANGENT = true;
    std::vector<FFVar> const tY = DAG.FAD(Y, X);
    std::ostringstream ost;
    ost << "tangent operation: " << jacerr(tY) << " (" << operations(DAG, tY) << ")";
    check(jacerr(tY) < 1e-13 && has(operations(DAG, tY), "TanKAN"), ost.str());
    Op.map().options.TANGENT = false;

    Op.map().options.ADJOINT = true;
    std::vector<FFVar> const aY = DAG.BAD(Y, X);
    std::ostringstream osa;
    osa << "adjoint operation: " << jacerr(aY) << " (" << operations(DAG, aY) << ")";
    check(jacerr(aY) < 1e-13 && has(operations(DAG, aY), "AdjKAN"), osa.str());
    Op.map().options.ADJOINT = false;
  }

  // Polyhedral relaxation through the expression tree of the network (RELAX = AUX), with each
  // encoding of the edges: ranges of the results, cuts and binary variables
  {
    Op.options.RELAX = {Op.options.AUX};
    typedef mc::FFPWLU<I>::Options ENC;
    int const enc[4]       = {ENC::HULL, ENC::CONT, ENC::BIN, ENC::SOS2};
    char const* encname[4] = {"HULL", "CONT", "BIN", "SOS2"};
    std::vector<I> IY(nout);
    NN.eval(IX.data(), IY.data());
    for (int ie = 0; ie < 4; ++ie)
    {
      Op.map().options.ENCODING = enc[ie];
      mc::PolImg<I> Env;
      std::vector<PV> PX, PY(nout);
      for (size_t i = 0; i < nin; ++i) PX.push_back(PV(&Env, X[i], IX[i]));
      DAG.eval(Y, PY, X, PX);
      Env.generate_cuts(PY);
      size_t nbin = 0, nsos = 0;
      for (auto const* pa : Env.Aux()) nbin += (pa->id().first == PV::AUXINT);
      for (auto const& pc : Env.Cuts()) nsos += (pc->type() == mc::PolCut<I>::SOS2);
      double erng = 0.;
      for (size_t k = 0; k < nout; ++k)
        erng = std::max({erng, std::fabs(mc::Op<I>::l(PY[k].range()) - mc::Op<I>::l(IY[k])),
                         std::fabs(mc::Op<I>::u(PY[k].range()) - mc::Op<I>::u(IY[k]))});
      std::ostringstream os;
      os << "polyhedral relaxation, " << encname[ie] << ": " << Env.Cuts().size() << " cuts, " << Env.Aux().size()
         << " auxiliary variables, " << nbin << " binary, " << nsos << " SOS2 sets, ranges off by " << erng;
      bool const ok = !Env.Cuts().empty() && erng < 1e-12 && (enc[ie] == ENC::BIN) == (nbin > 0) &&
                      (enc[ie] == ENC::SOS2) == (nsos > 0);
      check(ok, os.str());
    }
    Op.map().options.ENCODING = ENC::DEFAULT;
    Op.options.RELAX          = {Op.options.MC};
  }

  // Backward propagation of bounds on the results to the inputs: that of the network
  {
    std::vector<I> IY(nout), IXd = IX, IXn = IX, IYn(nout);
    NN.eval(IX.data(), IY.data());
    // upper tenth of the range of the first result
    IY[0] = I(mc::Op<I>::u(IY[0]) - 0.1 * mc::Op<I>::diam(IY[0]), mc::Op<I>::u(IY[0]));
    IYn   = IY;
    bool const feas = NN.reval(IXn.data(), IYn.data());
    auto sg         = DAG.subgraph(Y);
    std::vector<I> wk;
    int const flag = DAG.reval(sg, wk, Y.size(), Y.data(), IY.data(), X.size(), X.data(), IXd.data(), I(-1e30, 1e30),
                               NN.options.CPMAX, NN.options.CPTHRES);
    double d = 0.;
    for (size_t i = 0; i < nin; ++i)
      d = std::max({d, std::fabs(mc::Op<I>::l(IXd[i]) - mc::Op<I>::l(IXn[i])),
                    std::fabs(mc::Op<I>::u(IXd[i]) - mc::Op<I>::u(IXn[i]))});
    std::ostringstream os;
    os << "backward propagation: X1 in " << IXd[0] << ", X2 in " << IXd[1] << ", off those of the network by " << d;
    bool const contracted = mc::Op<I>::diam(IXd[0]) < mc::Op<I>::diam(IX[0]) || mc::Op<I>::diam(IXd[1]) < mc::Op<I>::diam(IX[1]);
    check(feas && flag >= 0 && d < 1e-12 && contracted, os.str());
  }

  // Network that the DAG holds: a copy by default, the network given with the policy SHALLOW
  {
    bool const copy = (&Op.map() != &NN) && (&static_cast<mc::FFKAN<I>*>(Y[0].opdef().first)->object() == &Op.map());
    FFGraph DAG2;
    std::vector<FFVar> X2(nin);
    for (auto& Xi : X2) Xi.set(&DAG2);
    mc::FFKAN<I> Op2;
    std::vector<FFVar> const Y2 = Op2(X2, &NN, Op2.SHALLOW);
    bool const shallow          = (&Op2.map() == &NN);
    double err                  = 0.;
    for (auto const& x : pts) err = std::max(err, dist(evaluate(DAG2, Y2, X2, x), evaluate(DAG, Y, X, x)));
    check(copy && shallow && err == 0., "network of the DAG: a copy by default, the network given with SHALLOW");
  }

  // The same class on a multilayer perceptron
  std::cout << "Multilayer perceptron through the same class" << std::endl;
  {
    mc::MLP<I> ANN;
    typedef std::vector<std::vector<double>> t_layer;
    ANN.append_data(t_layer{{0.1, 0.5, -0.3}, {-0.2, 0.8, 0.4}, {0.3, -0.6, 0.7}}, mc::MLP<I>::TANH);
    ANN.append_data(t_layer{{0.05, 1.2, -0.7, 0.4}});
    FFGraph DAGM;
    std::vector<FFVar> XM(2);
    for (auto& Xi : XM) Xi.set(&DAGM);
    mc::FFMLP<I> OpM;
    std::vector<FFVar> const YM = OpM(XM, &ANN), dYM = DAGM.FAD(YM, XM);
    double err = 0., errj = 0.;
    for (auto const& x : pts)
    {
      std::vector<double> y(1), J(2);
      ANN.grad(x.data(), J.data(), y.data());
      err  = std::max(err, dist(evaluate(DAGM, YM, XM, x), y));
      errj = std::max(errj, dist(evaluate(DAGM, dYM, XM, x), J));
    }
    std::ostringstream os;
    os << "FFMLP: operations " << operations(DAGM, YM) << "and " << operations(DAGM, dYM) << ", values " << err
       << ", Jacobian " << errj;
    check(operations(DAGM, YM) == "MLP " && has(operations(DAGM, dYM), "GradMLP") && err < 1e-14 && errj < 1e-13 &&
              &OpM.MLP() == &OpM.map() && &OpM.map() != &ANN,
          os.str());
  }

  std::cout << "=== Results: " << npass << "/" << npass + nfail << " passed ===" << std::endl;
  return nfail;
}
