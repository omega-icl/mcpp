// Copyright (C) Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the EPL-2.0 with GPL-2.0-or-later as a Secondary License; see the LICENSE file.
// SPDX-License-Identifier: EPL-2.0 OR GPL-2.0-or-later

// Kolmogorov-Arnold networks (mc::KAN): data, evaluation in double, interval, McCormick and
// superposition arithmetics, Jacobian, directional derivative and adjoint, expression tree,
// polyhedral image and backward propagation; and the networks read without Torch from
// dictionaries of torch.save (mc::NETIO::read_dict) by mc::KAN::read_data and mc::MLP::read_data,
// against their values and bounds computed in Python (dict_values.txt, from test_KAN_dict.py),
// with the files that must be refused.

#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <random>
#include <sstream>
#include <string>
#include <vector>

#include "ffkan.hpp"
#include "ffmlp.hpp"
#include "interval.hpp"

typedef mc::Interval I;
typedef mc::KAN<I> KAN;
using mc::PWLU;

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

int
main()
{
  std::mt19937 gen(5);
  std::uniform_real_distribution<double> U(-1., 1.);

  std::cout << "Data\n";
  KAN NN;
  auto const L0 = layer(gen, 4, 3, 5), L1 = layer(gen, 3, 4, 4), L2 = layer(gen, 2, 3, 6);
  bool ok = NN.append_data(L0) && NN.append_data(L1) && NN.append_data(L2);
  check(ok && NN.nin() == 3 && NN.nout() == 2 && NN.nlay() == 3,
        "three layers appended: 3 inputs, 2 outputs, " + std::to_string(NN.nwei()) + " segments");
  size_t const uid = NN.uid();
  check(!NN.append_data(layer(gen, 2, 3, 2)) && NN.uid() == uid && NN.nlay() == 3,
        "a layer with another number of inputs is refused, the network unchanged");
  {
    auto bad = layer(gen, 2, 2, 2);
    bad[1].pop_back();
    KAN N2;
    check(!N2.append_data(bad) && !N2.append_data({}) && !N2.set_data({}),
          "a ragged layer, an empty layer and an empty network are refused");
    KAN N3(NN), N4;
    N4 = NN;
    check(N3.uid() == uid && N4.uid() == uid && N3.nwei() == NN.nwei(),
          "a copy and an assignment keep the stamp of the network");
    N4.append_data(layer(gen, 1, 2, 2));
    check(N4.uid() != uid && N4.nout() == 1, "appending a layer draws a new stamp");
  }

  // Reference: the layers evaluated one edge at a time
  auto ref = [&](double const* x, double* y) {
    std::vector<double> u(x, x + 3), v;
    for (auto const* L : {&L0, &L1, &L2})
    {
      v.assign(L->size(), 0.);
      for (size_t j = 0; j < L->size(); ++j)
        for (size_t i = 0; i < u.size(); ++i) v[j] += (*L)[j][i].val(u[i], true);
      u = v;
    }
    y[0] = u[0];
    y[1] = u[1];
  };

  std::cout << "\nValues and derivatives\n";
  {
    double errv = 0., errj = 0., errfd = 0., errt = 0., erra = 0., errad = 0.;
    for (int k = 0; k < 200; ++k)
    {
      double x[3] = {U(gen), U(gen), U(gen)}, y[2], y0[2], Jf[6], Jb[6], yj[2];
      NN.eval(x, y);
      ref(x, y0);
      errv = std::max({errv, std::fabs(y[0] - y0[0]), std::fabs(y[1] - y0[1])});
      NN.options.AUTODIFF = KAN::Options::F;
      NN.grad(x, Jf, yj);  // 2 outputs for 3 inputs: backward by default
      NN.options.AUTODIFF = KAN::Options::B;
      NN.grad(x, Jb);
      // forward accumulation: the same network with its outputs summed into 4, more than its inputs
      for (int e = 0; e < 6; ++e) errj = std::max(errj, std::fabs(Jf[e] - Jb[e]));
      errv = std::max({errv, std::fabs(yj[0] - y0[0]), std::fabs(yj[1] - y0[1])});
      // finite differences, away from the breakpoints with high probability
      double const h = 1e-7;
      for (int i = 0; i < 3; ++i)
      {
        double xp[3] = {x[0], x[1], x[2]}, xm[3] = {x[0], x[1], x[2]}, yp[2], ym[2];
        xp[i] += h;
        xm[i] -= h;
        NN.eval(xp, yp);
        NN.eval(xm, ym);
        for (int o = 0; o < 2; ++o)
        {
          double const fd = (yp[o] - ym[o]) / (2. * h);
          if (std::fabs(fd - Jf[o * 3 + i]) < 1e-5 * (1. + std::fabs(fd))) continue;
          errfd += 1.;  // a breakpoint within h of the point: counted
        }
      }
      double d[3] = {U(gen), U(gen), U(gen)}, w[2] = {U(gen), U(gen)}, t[2], a[3];
      NN.tang(x, d, t);
      NN.adj(x, w, a);
      for (int o = 0; o < 2; ++o)
        errt = std::max(errt, std::fabs(t[o] - (Jf[o * 3] * d[0] + Jf[o * 3 + 1] * d[1] + Jf[o * 3 + 2] * d[2])));
      for (int i = 0; i < 3; ++i) erra = std::max(erra, std::fabs(a[i] - (Jf[i] * w[0] + Jf[3 + i] * w[1])));
      // forward and reverse modes of AD through the layers
      mc::FADType<double> FX[3], FY[2];
      mc::BADType<double> BX[3], BY[2];
      for (int i = 0; i < 3; ++i)
      {
        FX[i] = x[i];
        FX[i].diff(i, 3);
        BX[i] = x[i];
      }
      NN.eval(FX, FY);
      NN.eval(BX, BY);
      BY[0].diff(0, 2);
      BY[1].diff(1, 2);
      for (int o = 0; o < 2; ++o)
        for (int i = 0; i < 3; ++i)
          errad = std::max({errad, std::fabs(FY[o].d(i) - Jf[o * 3 + i]), std::fabs(BX[i].d(o) - Jf[o * 3 + i])});
    }
    check(errv < 1e-13, "values equal to the edges summed layer by layer");
    check(errj < 1e-12, "Jacobian: backward and forward accumulations agree");
    check(errfd <= 4., "Jacobian against finite differences (" + std::to_string((int)errfd) +
                           " of 1200 entries within 1e-7 of a breakpoint)");
    check(errt < 1e-12 && erra < 1e-12, "directional derivative and adjoint equal to the products with the Jacobian");
    check(errad < 1e-12, "forward and reverse modes of AD give the Jacobian");

    // a network with more outputs than inputs: forward accumulation by default
    KAN N5;
    N5.append_data(L0);
    double x[3] = {0.2, -0.4, 0.7}, J1[12], J2[12];
    N5.grad(x, J1);
    N5.options.AUTODIFF = KAN::Options::B;
    N5.grad(x, J2);
    double e5 = 0.;
    for (int e = 0; e < 12; ++e) e5 = std::max(e5, std::fabs(J1[e] - J2[e]));
    for (int j = 0; j < 4; ++j)
      for (int i = 0; i < 3; ++i) e5 = std::max(e5, std::fabs(J1[j * 3 + i] - L0[j][i].slope(x[i])));
    check(e5 < 1e-14, "Jacobian of a single layer: the slopes of its edges, forward and backward");
  }

  std::cout << "\nBounds and relaxations\n";
  {
    // Interval: valid, and exact for a single layer
    I IX[3] = {I(-0.8, 0.5), I(-0.2, 0.9), I(-1., 1.)}, IY[2];
    NN.eval(IX, IY);
    bool valid = true;
    typedef mc::McCormick<I> MC;
    MC MX[3], MY[2];
    double const xr[3] = {-0.1, 0.4, 0.3};
    for (int i = 0; i < 3; ++i) MX[i] = MC(IX[i], xr[i]).sub(3, i);
    NN.eval(MX, MY);
    bool mcvalid = true;
    for (int k = 0; k < 2000; ++k)
    {
      double x[3], y[2];
      for (int i = 0; i < 3; ++i) x[i] = mc::Op<I>::l(IX[i]) + mc::Op<I>::diam(IX[i]) * 0.5 * (1. + U(gen));
      NN.eval(x, y);
      for (int o = 0; o < 2; ++o)
      {
        valid = valid && y[o] >= mc::Op<I>::l(IY[o]) - 1e-12 && y[o] <= mc::Op<I>::u(IY[o]) + 1e-12;
        double cv = MY[o].cv(), cc = MY[o].cc();
        for (int i = 0; i < 3; ++i)
        {
          cv += MY[o].cvsub(i) * (x[i] - xr[i]);
          cc += MY[o].ccsub(i) * (x[i] - xr[i]);
        }
        mcvalid = mcvalid && cv <= y[o] + 1e-12 && cc >= y[o] - 1e-12;
      }
    }
    check(valid, "interval bounds enclose the network");
    check(mcvalid, "McCormick relaxations and their subgradient cuts enclose the network");

    KAN N1;
    N1.append_data(L0);
    I JY[4];
    N1.eval(IX, JY);
    double gap = 0.;
    for (int j = 0; j < 4; ++j)
    {
      double lo = 0., up = 0.;
      for (int i = 0; i < 3; ++i)
      {
        PWLU const R = L0[j][i].restrict(mc::Op<I>::l(IX[i]), mc::Op<I>::u(IX[i]));
        lo += R.l();
        up += R.u();
      }
      gap = std::max({gap, std::fabs(mc::Op<I>::l(JY[j]) - lo), std::fabs(mc::Op<I>::u(JY[j]) - up)});
    }
    check(gap < 1e-13, "interval bounds of a single layer: the exact range of the layer");

    // Superposition models: a single layer is a sum of univariates, hence exact
    mc::SupModel<PWLU> mod(3);
    mc::SupVar<PWLU> SX[3], SY[4];
    for (int i = 0; i < 3; ++i) SX[i].set(mod, i, IX[i], 4);
    N1.eval(SX, SY);
    double sgap = 0.;
    for (int k = 0; k < 500; ++k)
    {
      double x[3], y[4];
      std::map<unsigned, double> pt;
      for (int i = 0; i < 3; ++i) pt[i] = x[i] = mc::Op<I>::l(IX[i]) + mc::Op<I>::diam(IX[i]) * 0.5 * (1. + U(gen));
      N1.eval(x, y);
      for (int j = 0; j < 4; ++j) sgap = std::max({sgap, std::fabs(SY[j].uval(pt) - y[j]), std::fabs(SY[j].oval(pt) - y[j])});
    }
    check(sgap < 1e-12, "superposition estimators of a single layer equal to the layer");
    mc::SupVar<PWLU> SZ[2];
    NN.eval(SX, SZ);
    bool svalid = true;
    for (int k = 0; k < 500; ++k)
    {
      double x[3], y[2];
      std::map<unsigned, double> pt;
      for (int i = 0; i < 3; ++i) pt[i] = x[i] = mc::Op<I>::l(IX[i]) + mc::Op<I>::diam(IX[i]) * 0.5 * (1. + U(gen));
      NN.eval(x, y);
      for (int o = 0; o < 2; ++o) svalid = svalid && SZ[o].uval(pt) <= y[o] + 1e-10 && SZ[o].oval(pt) >= y[o] - 1e-10;
    }
    check(svalid, "superposition estimators enclose the network");
  }

  std::cout << "\nExpression tree\n";
  {
    mc::FFGraph* dag = NN.dag();
    auto const& X = NN.varin();
    auto const& Y = NN.varout();
    double errd = 0.;
    for (int k = 0; k < 50; ++k)
    {
      double x[3] = {U(gen), U(gen), U(gen)}, y[2], yd[2];
      NN.eval(x, y);
      dag->eval(2, Y.data(), yd, 3, X.data(), x);
      errd = std::max({errd, std::fabs(y[0] - yd[0]), std::fabs(y[1] - yd[1])});
    }
    auto const sg = dag->subgraph(Y);
    size_t next = 0;
    for (auto const* op : sg.l_op) next += (op->type == mc::FFOp::EXTERN);
    // one FFPWLU per input of a layer (3+4+3), one FFLin per output (4+3+2)
    check(errd < 1e-13 && next == 19, "tree equal to the network: " + std::to_string(next) + " operations, one per input and one per output of each layer");

    // Inserted in another DAG
    mc::FFGraph DAG;
    std::vector<mc::FFVar> Z = DAG.add_vars(3, "Z"), F(2);
    NN.eval(Z.data(), F.data());
    double x[3] = {0.3, -0.6, 0.1}, y[2], yd[2];
    NN.eval(x, y);
    DAG.eval(2, F.data(), yd, 3, Z.data(), x);
    check(std::fabs(y[0] - yd[0]) < 1e-13 && std::fabs(y[1] - yd[1]) < 1e-13, "network inserted in another DAG");

    // Polyhedral image: ranges of the outputs equal to the interval bounds, cuts for each encoding
    I IX[3] = {I(-0.8, 0.5), I(-0.2, 0.9), I(-1., 1.)}, IY[2];
    NN.eval(IX, IY);
    bool okp = true;
    std::string counts;
    for (int enc : {0, 1, 2, 3})
      for (bool sep : {false, true})
      {
        NN.options.ENCODING = enc;
        NN.options.SEPARATE = sep;
        mc::PolImg<I> Env;
        mc::PolVar<I> PX[3], PY[2];
        for (int i = 0; i < 3; ++i) PX[i].set(&Env, NN.varin()[i], IX[i]);
        NN.eval(PX, PY);
        Env.generate_cuts(2, PY);
        for (int o = 0; o < 2; ++o)
          okp = okp && std::fabs(mc::Op<I>::l(PY[o].range()) - mc::Op<I>::l(IY[o])) < 1e-12 &&
                std::fabs(mc::Op<I>::u(PY[o].range()) - mc::Op<I>::u(IY[o])) < 1e-12;
        okp = okp && !Env.Cuts().empty() && (enc == 0) == Env.Aux().empty();
        counts += " " + std::to_string(Env.Cuts().size()) + "/" + std::to_string(Env.Aux().size());
      }
    NN.options.reset();
    check(okp, "polyhedral image for each encoding, shared and separate (cuts/auxiliaries):" + counts);

    // Backward propagation: points with outputs in the range asked stay inside the contracted box
    I CX[3] = {IX[0], IX[1], IX[2]}, CY[2] = {I(mc::Op<I>::l(IY[0]), mc::Op<I>::mid(IY[0])), IY[1]};
    NN.options.CPMAX = 5;
    bool const feas = NN.reval(CX, CY);
    bool kept = true, tighter = false;
    for (int i = 0; i < 3; ++i) tighter = tighter || mc::Op<I>::diam(CX[i]) < mc::Op<I>::diam(IX[i]) - 1e-9;
    for (int k = 0; k < 5000; ++k)
    {
      double x[3], y[2];
      for (int i = 0; i < 3; ++i) x[i] = mc::Op<I>::l(IX[i]) + mc::Op<I>::diam(IX[i]) * 0.5 * (1. + U(gen));
      NN.eval(x, y);
      if (y[0] > mc::Op<I>::mid(IY[0])) continue;
      for (int i = 0; i < 3; ++i) kept = kept && x[i] >= mc::Op<I>::l(CX[i]) - 1e-9 && x[i] <= mc::Op<I>::u(CX[i]) + 1e-9;
    }
    check(feas && kept, std::string("backward propagation keeps the points with outputs in the range asked") +
                            (tighter ? " (box contracted)" : ""));
    // a layer of increasing edges: an upper bound on the output contracts the inputs from above
    {
      KAN N6;
      N6.append_data({{PWLU({0., 1.}, {0., 2.}, 0.5, 1.), PWLU({-1., 0.5}, {-1., 0.}, 1., 3.)}});
      N6.options.CPMAX = 5;
      I DX[2] = {I(-1., 2.), I(-2., 1.)}, DY[1] = {I(-10., -1.)};
      // output in [-10,-1] with f1 in [-0.5,3], f2 in [-2,1.5]: f1 <= 1 and f2 <= -0.5, i.e. x1 <= 0.5
      // and x2 <= -0.25
      bool const f6 = N6.reval(DX, DY);
      check(f6 && std::fabs(mc::Op<I>::u(DX[0]) - 0.5) < 1e-12 && std::fabs(mc::Op<I>::u(DX[1]) + 0.25) < 1e-12 &&
                mc::Op<I>::l(DX[0]) == -1. && mc::Op<I>::l(DX[1]) == -2.,
            "backward propagation contracts the inputs of increasing edges to [-1,0.5] and [-2,-0.25]");
    }
    I EX[3] = {IX[0], IX[1], IX[2]}, EY[2] = {I(mc::Op<I>::u(IY[0]) + 1., mc::Op<I>::u(IY[0]) + 2.), IY[1]};
    check(!NN.reval(EX, EY), "backward propagation detects outputs out of reach");
    NN.options.reset();
  }

  std::cout << "\nDictionaries of torch.save, read without Torch\n";
  {
    std::ifstream is("dict_values.txt");
    if (!is)
    {
      std::cout << "dict_values.txt not found: run test_KAN_dict.py\n";
      return 1;
    }
    std::string line, last;
    KAN kan, kanm;
    mc::MLP<I> mlp, mlpm;
    bool iskan = false, read = false;
    double err = 0., errm = 0.;
    auto report = [&]() {
      if (last.empty() || !read) return;
      std::ostringstream os;
      os << last << ": values, error " << err << " (raw), " << errm << " (scaled coordinates)";
      check(err < 1e-12 && errm < 1e-12, os.str());
      err = errm = 0.;
    };
    while (std::getline(is, line))
    {
      std::istringstream ls(line);
      std::string tag, name;
      ls >> tag >> name;
      if (name != last)
      {
        report();
        last  = name;
        iskan = (name.find("_kan_") != std::string::npos);
        // a network without activation in the file or bounds to recognise it from is read with
        // the activation given
        if (name == "dict_ann_plain.pt")
        {
          mc::NETIO::Dict D;
          std::string msg;
          bool const refused = !mlp.read_data(name) && !mc::NETIO::read_dict(name, D, msg);
          read               = mc::NETIO::read_dict(name, D, msg, mc::NETIO::RELU) && mlpm.set_data(D.ann) &&
                 mc::NETIO::fold_scaling(D.ann, D, msg) && mlp.set_data(D.ann);
          check(refused && read, name + ": refused where the activation is neither named nor given, read with it given");
        }
        else
        {
          // the network with the scalings of the file folded in (read_data), and the one in the
          // scaled coordinates of the file, from its layers as read_dict returns them
          mc::NETIO::Dict D;
          std::string msg;
          read = mc::NETIO::read_dict(name, D, msg);
          if (read && iskan)
          {
            std::vector<mc::KAN<I>::Layer> L(D.kan.size());
            for (size_t l = 0; l < L.size(); ++l)
            {
              auto const& K = D.kan[l];
              L[l].pwl      = mc::KAN<I>::relu_layer(K.nout, K.nin, K.R, K.c.data(), K.v.data(), K.a.data(),
                                                     K.w.data(), K.b.data());
            }
            read = kan.read_data(name) && kanm.set_data(L);
          }
          else if (read)
            read = mlp.read_data(name) && mlpm.set_data(D.ann);
          std::ostringstream os;
          os << name << ": read, ";
          if (iskan)
            os << kan.nin() << " inputs, " << kan.nlay() << " layers, " << kan.nout() << " outputs, "
               << kan.nwei() << " segments";
          else
            os << mlp.nin() << " inputs, " << mlp.nhid() << " hidden layers, " << mlp.nout() << " outputs";
          check(read, os.str());
        }
      }
      if (!read) continue;
      if (tag == "value")
      {
        size_t nin, nout;
        ls >> nin >> nout;
        std::vector<double> x(nin), z(nout), y(nout), v(nout), vm(nout), u(nin);
        for (auto& xi : x) ls >> xi;
        for (auto& zi : z) ls >> zi;
        for (auto& yi : y) ls >> yi;
        // inputs in the scaled coordinates, from the scalings of the file
        mc::NETIO::Dict D;
        std::string msg;
        mc::NETIO::read_dict(name, D, msg, mc::NETIO::RELU);
        for (size_t i = 0; i < nin; ++i)
          u[i] = (D.xlower.empty() ? x[i] : (x[i] - D.xlower[i]) / (D.xupper[i] - D.xlower[i]));
        if (iskan)
        {
          kan.eval(x.data(), v.data());
          kanm.eval(u.data(), vm.data());
        }
        else
        {
          mlp.eval(x.data(), v.data());
          mlpm.eval(u.data(), vm.data());
        }
        for (size_t k = 0; k < nout; ++k)
        {
          err  = std::max(err, std::fabs(v[k] - z[k]) / (1. + std::fabs(z[k])));
          errm = std::max(errm, std::fabs(vm[k] - y[k]) / (1. + std::fabs(y[k])));
        }
      }
      else if (tag == "bound")
      {
        // interval bounds on the scaled input box [0,1]^n against those of the file: exact
        // ranges of the edges for a KAN, natural interval bounds for a dense network
        size_t nout;
        ls >> nout;
        std::vector<double> lo(nout), up(nout);
        for (auto& l : lo) ls >> l;
        for (auto& u : up) ls >> u;
        size_t const nin = (iskan ? kanm.nin() : mlpm.nin());
        std::vector<I> IX(nin, I(0., 1.)), IY(nout);
        if (iskan)
          kanm.eval(IX.data(), IY.data());
        else
          mlpm.eval(IX.data(), IY.data());
        double eb = 0.;
        for (size_t k = 0; k < nout; ++k)
          eb = std::max({eb, std::fabs(mc::Op<I>::l(IY[k]) - lo[k]) / (1. + std::fabs(lo[k])),
                         std::fabs(mc::Op<I>::u(IY[k]) - up[k]) / (1. + std::fabs(up[k]))});
        std::ostringstream os;
        os << name << ": interval bounds equal to those of the file, error " << eb;
        check(eb < 1e-12, os.str());
      }
    }
    report();
  }

  std::cout << "\nFiles refused\n";
  {
    for (std::string const name : {"dict_bad_type.pt", "dict_bad_shape.pt", "dict_bad_missing.pt",
                                   "dict_bad_scaling.pt", "dict_bad_tensor.pt", "dict_bad_activ.pt",
                                   "dict_bad_truncated.pt", "dict_bad_deflated.pt", "dict_missing.pt",
                                   "torch_tanh.pt", "onnx_torch_tanh.onnx"})
    {
      mc::NETIO::Dict D;
      std::string msg;
      bool const refused = !mc::NETIO::read_dict(name, D, msg);
      check(refused && !msg.empty(), name + ": " + msg);
    }
    KAN kan;
    kan.read_data("dict_kan_two.pt");
    size_t const uid = kan.uid();
    check(!kan.read_data("dict_bad_shape.pt") && !kan.read_data("dict_ann_tanh.pt") && kan.uid() == uid &&
              kan.nin() == 3,
          "a refused file, or a dense network, leaves the KAN as it was");
    mc::MLP<I> mlp;
    mlp.read_data("dict_ann_tanh.pt");
    size_t const uid2 = mlp.uid();
    check(!mlp.read_data("dict_kan_two.pt") && mlp.uid() == uid2 && mlp.nin() == 3,
          "a Kolmogorov-Arnold network leaves the MLP as it was");
    check(mc::NETIO::is_torch_dict("dict_kan_two.pt") && !mc::NETIO::is_torch_dict("torch_tanh.pt") &&
              !mc::NETIO::is_torch_dict("onnx_torch_tanh.onnx") && !mc::NETIO::is_torch_dict("dict_missing.pt"),
          "a dictionary is told from a TorchScript Module, an ONNX file and a missing file");
  }

  // A dictionary read once, then given to the network: set_data(dict) is read_data without the file
  {
    mc::NETIO::Dict DK, DA;
    std::string msg;
    bool const read = mc::NETIO::read_dict("dict_kan_scaled.pt", DK, msg) && mc::NETIO::read_dict("dict_ann_relu.pt", DA, msg);
    KAN kanf, kand;
    mc::MLP<I> mlpf, mlpd;
    bool const set = read && kanf.read_data("dict_kan_scaled.pt") && kand.set_data(DK) && mlpf.read_data("dict_ann_relu.pt") &&
                     mlpd.set_data(DA);
    double d = 0.;
    if (set)
    {
      std::vector<double> x(kanf.nin(), 0.3), yf(kanf.nout()), yd(kanf.nout());
      for (size_t i = 0; i < x.size(); ++i) x[i] = DK.xlower.empty() ? 0.3 + 0.1 * i : DK.xlower[i] + (0.3 + 0.1 * i) * (DK.xupper[i] - DK.xlower[i]);
      kanf.eval(x.data(), yf.data());
      kand.eval(x.data(), yd.data());
      for (size_t k = 0; k < yf.size(); ++k) d = std::max(d, std::fabs(yf[k] - yd[k]));
      std::vector<double> u(mlpf.nin(), 0.4), zf(mlpf.nout()), zd(mlpf.nout());
      mlpf.eval(u.data(), zf.data());
      mlpd.eval(u.data(), zd.data());
      for (size_t k = 0; k < zf.size(); ++k) d = std::max(d, std::fabs(zf[k] - zd[k]));
    }
    check(set && d == 0., "set_data(dict): the network of read_data, for a KAN and for a dense network");
    size_t const nk = kand.nwei(), nm = mlpd.nin();
    std::string why;
    bool const refused = !kand.set_data(DA, &why) && !why.empty() && kand.nwei() == nk && !mlpd.set_data(DK) && mlpd.nin() == nm;
    check(refused, "set_data(dict) of the other type of network is refused: " + why);
    check(!kand.set_data({}) && kand.nwei() == nk, "set_data({}) remains a call with layers, refused");
  }

  std::cout << "\n=== Results: " << npass << "/" << npass + nfail << " passed ===\n";
  return nfail ? 1 : 0;
}
