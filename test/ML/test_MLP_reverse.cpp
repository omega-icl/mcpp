// =============================================================================
// test_MLP_reverse.cpp
// Unit tests for the evaluation of FFMLP and FFMLPREG operations in the
// reverse mode of AD (mc::BADType<double>): each output is one node whose
// partial derivatives are the Jacobian of the network. The gradient of a
// training loss in the weights, and of a network in its inputs, agree with the
// forward mode (mc::FADType<double>) and with symbolic differentiation.
//
// Build (assuming MC++ headers are on the include path):
//   g++ -std=c++17 -O2 -I<path-to-mc++> -o test_MLP_reverse test_MLP_reverse.cpp
//
// The process exits with code 0 if every test passes, 1 otherwise.
// =============================================================================

#include <algorithm>
#include <cmath>
#include <iostream>
#include <random>
#include <string>
#include <vector>

#include "interval.hpp"
typedef mc::Interval I;

#include "ffmlpreg.hpp"
#include "ffmlp.hpp"

static int npass = 0, nfail = 0;
static void
check(bool const ok, std::string const& what)
{
  std::cout << (ok ? "  PASS  " : "  FAIL  ") << what << std::endl;
  (ok ? npass : nfail)++;
}

// Largest difference between two vectors, relative to 1 + the largest entry
static double
reldiff(std::vector<double> const& a, std::vector<double> const& b)
{
  double d = 0., m = 0.;
  for (size_t i = 0; i < a.size(); ++i)
  {
    d = std::max(d, std::fabs(a[i] - b[i]));
    m = std::max(m, std::fabs(a[i]));
  }
  return d / (1. + m);
}

// Gradient of F in the variables X at x: reverse mode, forward mode, symbolic
// (FFGraph::BAD)
static void
gradients(mc::FFGraph& DAG, mc::FFVar const& F, std::vector<mc::FFVar> const& X,
          std::vector<double> const& x, double& f, std::vector<double>& gb,
          std::vector<double>& gf, std::vector<double>& gs)
{
  size_t const n = X.size();
  auto sg = DAG.subgraph(1, &F);

  std::vector<mc::BADType<double>> BX(n), Bwk;
  for (size_t i = 0; i < n; ++i) BX[i] = x[i];
  mc::BADType<double> BF;
  DAG.eval(sg, Bwk, 1, &F, &BF, n, X.data(), BX.data());
  Bwk.clear();
  BF.diff(0, 1);
  f = BF.val();
  gb.resize(n);
  for (size_t i = 0; i < n; ++i) gb[i] = BX[i].d(0);

  std::vector<mc::FADType<double>> FX(n), Fwk;
  for (size_t i = 0; i < n; ++i)
  {
    FX[i] = x[i];
    FX[i].diff(i, n);
  }
  mc::FADType<double> FF;
  DAG.eval(sg, Fwk, 1, &F, &FF, n, X.data(), FX.data());
  gf.resize(n);
  for (size_t i = 0; i < n; ++i) gf[i] = FF.d(i);

  auto G   = DAG.BAD({F}, X);
  auto sgG = DAG.subgraph(G);
  std::vector<double> dwk;
  gs.resize(n);
  DAG.eval(sgG, dwk, n, G.data(), gs.data(), n, X.data(), x.data());
}

int
main()
{
  std::mt19937 gen(1);
  std::uniform_real_distribution<double> U(-0.5, 0.5);

  // MLPREG: training loss over three samples, the inputs being DAG constants
  {
    mc::MLPREG<I> REG(2, {{5, mc::MLPREG<I>::TANH},
                          {4, mc::MLPREG<I>::SIGMOID},
                          {1, mc::MLPREG<I>::LINEAR}});
    mc::FFGraph DAG;
    auto W = DAG.add_vars(REG.nwei(), "W");
    mc::FFMLPREG<I> Op;
    double const xs[3][3] = {{0.1, -0.4, 0.3}, {-0.7, 0.2, -0.1}, {0.5, 0.6, 0.8}};
    mc::FFVar Loss = 0.;
    for (auto const& s : xs)
    {
      mc::FFVar Y = Op(0, {mc::FFVar(&DAG, s[0]), mc::FFVar(&DAG, s[1])}, W, &REG);
      Loss += sqr(Y - s[2]);
    }
    std::vector<double> w(REG.nwei());
    for (auto& wi : w) wi = U(gen);
    double f;
    std::vector<double> gb, gf, gs;
    gradients(DAG, Loss, W, w, f, gb, gf, gs);
    double f0;
    auto sg = DAG.subgraph(1, &Loss);
    std::vector<double> dwk;
    DAG.eval(sg, dwk, 1, &Loss, &f0, W.size(), W.data(), w.data());
    check(f == f0, "MLPREG loss: reverse-mode value equal to the value in double");
    check(reldiff(gb, gf) < 1e-13, "MLPREG loss: gradient in the weights, reverse mode = forward mode");
    check(reldiff(gb, gs) < 1e-13, "MLPREG loss: gradient in the weights, reverse mode = symbolic");
    // a second evaluation at another point
    for (auto& wi : w) wi = U(gen);
    gradients(DAG, Loss, W, w, f, gb, gf, gs);
    check(reldiff(gb, gf) < 1e-13, "MLPREG loss: again at another point");
  }

  // MLPREG: inputs and weights as variables, two outputs combined
  {
    mc::MLPREG<I> REG(3, {{6, mc::MLPREG<I>::RELU}, {2, mc::MLPREG<I>::LINEAR}});
    mc::FFGraph DAG;
    auto X = DAG.add_vars(3, "X");
    auto W = DAG.add_vars(REG.nwei(), "W");
    mc::FFMLPREG<I> Op;
    auto Y = Op(X, W, &REG);
    mc::FFVar F = Y[0] * exp(Y[1]) + X[2];
    std::vector<mc::FFVar> V(X);
    V.insert(V.end(), W.begin(), W.end());
    std::vector<double> v(V.size());
    for (auto& vi : v) vi = U(gen);
    double f;
    std::vector<double> gb, gf, gs;
    gradients(DAG, F, V, v, f, gb, gf, gs);
    check(reldiff(gb, gf) < 1e-13 && reldiff(gb, gs) < 1e-13,
          "MLPREG with two outputs: gradient in the inputs and the weights");
  }

  // MLP: gradient in the inputs
  {
    mc::MLP<I> NN;
    std::vector<std::vector<double>> L1(4, std::vector<double>(3)), L2(1, std::vector<double>(5));
    for (auto& r : L1) for (auto& c : r) c = U(gen);
    for (auto& r : L2) for (auto& c : r) c = U(gen);
    NN.append_data(L1, mc::MLP<I>::TANH);
    NN.append_data(L2, mc::MLP<I>::LINEAR);
    mc::FFGraph DAG;
    auto X = DAG.add_vars(2, "X");
    mc::FFMLP<I> Op;
    mc::FFVar F = sqr(Op(0, X, &NN)) + X[0];
    std::vector<double> x = {0.3, -0.2};
    double f;
    std::vector<double> gb, gf, gs;
    gradients(DAG, F, X, x, f, gb, gf, gs);
    check(reldiff(gb, gf) < 1e-13 && reldiff(gb, gs) < 1e-13, "MLP: gradient in the inputs");
  }

  std::cout << npass << " passed, " << nfail << " failed" << std::endl;
  return nfail ? 1 : 0;
}
