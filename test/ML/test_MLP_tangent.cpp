// =============================================================================
// test_MLP_tangent.cpp
// Unit tests for the directional derivative of multilayer perceptrons:
// mc::MLP::tang, the external operation mc::FFTanMLP, and the forward
// symbolic differentiation of mc::FFMLP through the tangent hook
// (FFOp::has_tangent / FFOp::tangent) when mc::MLP::Options::TANGENT is set,
// against the Jacobian operation mc::FFGradMLP that it replaces.
//
// Build (assuming MC++ headers are on the include path):
//   g++ -std=c++17 -O2 -I<path-to-mc++> -o test_MLP_tangent test_MLP_tangent.cpp
//
// Each test prints PASS or FAIL and a summary is printed at the end.
// The process exits with code 0 if every test passes, 1 otherwise.
// =============================================================================

#include <cmath>
#include <cstdio>
#include <functional>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "interval.hpp"
typedef mc::Interval I;

#include "ffmlp.hpp"

typedef mc::McCormick<I> MC;
using mc::FFGraph;
using mc::FFVar;

// -----------------------------------------------------------------------------
// Lightweight test harness
// -----------------------------------------------------------------------------
namespace test
{

static int total = 0, passed = 0, failed = 0;

struct Result
{
  bool ok;
  std::string msg;
};

inline bool
near(double a, double b, double tol = 1e-12)
{
  return std::fabs(a - b) <= tol * (1. + std::fabs(b));
}

inline bool
near(std::vector<double> const& a, std::vector<double> const& b, double tol = 1e-12)
{
  if (a.size() != b.size()) return false;
  for (size_t i = 0; i < a.size(); ++i)
    if (!near(a[i], b[i], tol)) return false;
  return true;
}

inline void
run(std::string const& name, std::function<Result()> fn)
{
  ++total;
  Result r;
  try
  {
    r = fn();
  }
  catch (std::exception const& e)
  {
    r = {false, std::string("exception: ") + e.what()};
  }
  catch (mc::FFBase::Exceptions& e)
  {
    r = {false, std::string("FFBase exception: ") + e.what()};
  }
  catch (...)
  {
    r = {false, "unknown exception"};
  }
  if (r.ok)
  {
    ++passed;
    std::printf("  PASS  %s\n", name.c_str());
  }
  else
  {
    ++failed;
    std::printf("  FAIL  %s -- %s\n", name.c_str(), r.msg.c_str());
  }
}

}  // namespace test

// -----------------------------------------------------------------------------
// Network with 2 inputs, hidden layers of 3 (tanh) and 2 (sigmoid) neurons, and 2 linear outputs;
// each neuron is a bias followed by the weights of the previous layer
// -----------------------------------------------------------------------------
static void
network(mc::MLP<I>& ANN, bool const tangent)
{
  ANN.append_data({{0.1, 0.8, -0.5}, {-0.3, 0.4, 0.9}, {0.2, -0.7, 0.6}}, ANN.TANH, true);
  ANN.append_data({{0.05, 1.2, -0.4, 0.3}, {-0.1, -0.6, 0.7, 0.5}}, ANN.SIGMOID);
  ANN.append_data({{0.3, 1.5, -2.0}, {-0.2, 0.5, 0.25}}, ANN.LINEAR);
  ANN.options.TANGENT = tangent;
}

// Network with the other activation functions: ReLU, sigmoid in its exponential form, and tanh outputs
static void
network2(mc::MLP<I>& ANN)
{
  ANN.append_data({{0.1, 0.8, -0.5}, {-0.3, 0.4, 0.9}, {0.2, -0.7, 0.6}}, ANN.RELU, true);
  ANN.append_data({{0.05, 1.2, -0.4, 0.3}, {-0.1, -0.6, 0.7, 0.5}}, ANN.SIGMOID);
  ANN.append_data({{0.3, 1.5, -2.0}, {-0.2, 0.5, 0.25}}, ANN.TANH);
  ANN.options.SIG2EXP = true;
}

static std::vector<double>
evaluate(FFGraph& DAG, std::vector<FFVar> const& F, std::vector<FFVar> const& Q,
         std::vector<double> const& q)
{
  std::vector<double> f(F.size());
  if (!F.empty()) DAG.eval(F, f, Q, q);
  return f;
}

// Names of the external operations in the subgraph of F
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

static std::vector<std::vector<double>> const points = {{0.3, -0.7}, {-1.2, 0.4}, {2.0, 1.5}};

int
main()
{
#if defined(MC__DAG_FADIFF) || defined(MC__DAG_BADIFF)
  // The tangent and adjoint hooks belong to the default differentiation engines: the legacy engines, in
  // FADBAD++ arithmetic, do not call them
  std::printf("test_MLP_tangent: skipped with MC__DAG_FADIFF or MC__DAG_BADIFF, the hooks belong to the default engines\n");
  return 0;
#endif

  using test::near;
  using test::Result;
  using test::run;

  std::printf("Directional derivative of multilayer perceptrons\n");

  run("MLP::tang against the Jacobian of MLP::grad", []() -> Result {
    mc::MLP<I> ANN;
    network(ANN, false);
    for (auto const& x : points)
    {
      std::vector<double> J(4), d = {0.5, -2.}, t(2), y(2), yref(2);
      ANN.grad(x.data(), J.data(), yref.data());
      ANN.tang(x.data(), d.data(), t.data(), y.data());
      if (!near(y, yref) || !near(t, {J[0] * d[0] + J[1] * d[1], J[2] * d[0] + J[3] * d[1]}))
        return {false, "tangent or value"};
    }
    return {true, ""};
  });

  run("MLP::tang with ReLU, exponential-form sigmoid and tanh output activations", []() -> Result {
    mc::MLP<I> ANN;
    network2(ANN);
    for (auto const& x : points)
    {
      std::vector<double> J(4), d = {0.5, -2.}, t(2), y(2), yref(2);
      ANN.grad(x.data(), J.data(), yref.data());
      ANN.tang(x.data(), d.data(), t.data(), y.data());
      if (!near(y, yref) || !near(t, {J[0] * d[0] + J[1] * d[1], J[2] * d[0] + J[3] * d[1]}))
        return {false, "tangent or value"};
      // enclosure in interval arithmetic
      std::vector<I> const IX = {x[0] + I(-0.1, 0.1), x[1] + I(-0.1, 0.1)}, ID = {I(d[0]), I(d[1])};
      std::vector<I> IT(2), IY(2);
      ANN.tang(IX.data(), ID.data(), IT.data(), IY.data());
      for (int k = 0; k < 2; ++k)
        if (t[k] < IT[k].l() - 1e-12 || t[k] > IT[k].u() + 1e-12 || y[k] < IY[k].l() - 1e-12 || y[k] > IY[k].u() + 1e-12)
          return {false, "interval bounds"};
    }
    return {true, ""};
  });

  run("FAD: same derivatives with Options::TANGENT, through FFTanMLP operations", []() -> Result {
    mc::MLP<I> ANN, ANNT;
    network(ANN, false);
    network(ANNT, true);
    FFGraph DAG;
    std::vector<FFVar> Q = DAG.add_vars(3);
    // inputs as expressions of the independents, outputs in further expressions
    std::vector<FFVar> const X = {Q[0] * Q[1], Q[1] - 2. * Q[2]};
    mc::FFMLP<I> Op, OpT;
    std::vector<FFVar> const Y = Op(X, &ANN), YT = OpT(X, &ANNT);
    std::vector<FFVar> const F = {Y[0] * Y[1] + Q[2], exp(Y[1])}, FT = {YT[0] * YT[1] + Q[2], exp(YT[1])};
    std::vector<FFVar> const dF = DAG.FAD(F, Q), dFT = DAG.FAD(FT, Q);
    std::string const ops = operations(DAG, dF), opsT = operations(DAG, dFT);
    if (ops.find("GradMLP") == std::string::npos || ops.find("TanMLP") != std::string::npos)
      return {false, "operations without the option: " + ops};
    if (opsT.find("TanMLP") == std::string::npos || opsT.find("GradMLP") != std::string::npos)
      return {false, "operations with the option: " + opsT};
    for (auto const& p : points)
    {
      std::vector<double> const q = {p[0], p[1], 0.6};
      if (!near(evaluate(DAG, dFT, Q, q), evaluate(DAG, dF, Q, q))) return {false, "derivatives differ"};
    }
    return {true, ""};
  });

  run("SFAD and DFAD: same derivatives with Options::TANGENT", []() -> Result {
    mc::MLP<I> ANN, ANNT;
    network(ANN, false);
    network(ANNT, true);
    FFGraph DAG;
    std::vector<FFVar> Q = DAG.add_vars(3), D = DAG.add_vars(3);
    std::vector<FFVar> const X = {Q[0] * Q[1], sqr(Q[1])};  // Q[2] does not reach the network
    mc::FFMLP<I> Op, OpT;
    std::vector<FFVar> const F = Op(X, &ANN), FT = OpT(X, &ANNT);
    auto [r0, c0, v0] = DAG.SFAD(F, Q);
    auto [r1, c1, v1] = DAG.SFAD(FT, Q);
    if (r0 != r1 || c0 != c1 || r1.size() != 4) return {false, "sparsity patterns differ"};
    std::vector<FFVar> const dir = DAG.DFAD(F, Q, D), dirT = DAG.DFAD(FT, Q, D);
    std::vector<FFVar> QD = Q;
    QD.insert(QD.end(), D.begin(), D.end());
    for (auto const& p : points)
    {
      std::vector<double> const q = {p[0], p[1], 0.6}, qd = {p[0], p[1], 0.6, 0.5, -1., 2.};
      if (!near(evaluate(DAG, v1, Q, q), evaluate(DAG, v0, Q, q))) return {false, "sparse derivatives differ"};
      if (!near(evaluate(DAG, dirT, QD, qd), evaluate(DAG, dir, QD, qd))) return {false, "directional derivatives differ"};
    }
    return {true, ""};
  });

  run("the Jacobian operation is kept where the sweeps through the network are many", []() -> Result {
    mc::MLP<I> ANN, ANNT;
    network(ANN, false);
    network(ANNT, true);
    FFGraph DAG;
    std::vector<FFVar> Q = DAG.add_vars(7), D = DAG.add_vars(7);
    mc::FFMLP<I> Op, OpT;
    // with 2 inputs the directional derivative is used for up to 5 sweeps: here 5, then 7 independents
    std::vector<FFVar> const X5 = {Q[0] + Q[1] * Q[2], Q[3] * Q[4]}, X7 = {Q[0] + Q[1] * Q[2], Q[3] * Q[4] + Q[5] * Q[6]};
    std::vector<FFVar> const F5 = OpT(X5, &ANNT), F7 = OpT(X7, &ANNT), G7 = Op(X7, &ANN);
    std::vector<FFVar> const dF5 = DAG.FAD(F5, Q), dF7 = DAG.FAD(F7, Q), dG7 = DAG.FAD(G7, Q), dirF7 = DAG.DFAD(F7, Q, D);
    std::string ops = operations(DAG, dF5);
    if (ops.find("TanMLP") == std::string::npos || ops.find("GradMLP") != std::string::npos) return {false, "5 sweeps: " + ops};
    ops = operations(DAG, dF7);
    if (ops.find("GradMLP") == std::string::npos || ops.find("TanMLP") != std::string::npos) return {false, "7 sweeps: " + ops};
    ops = operations(DAG, dirF7);
    if (ops.find("TanMLP") == std::string::npos || ops.find("GradMLP") != std::string::npos) return {false, "one direction: " + ops};
    std::vector<double> const q = {0.3, -0.7, 1.1, 0.4, -2., 0.8, 0.1};
    return {near(evaluate(DAG, dF7, Q, q), evaluate(DAG, dG7, Q, q)), "derivatives differ"};
  });

  run("BAD: the Jacobian operation either way", []() -> Result {
    mc::MLP<I> ANNT;
    network(ANNT, true);
    FFGraph DAG;
    std::vector<FFVar> Q = DAG.add_vars(2);
    mc::FFMLP<I> OpT;
    std::vector<FFVar> const FT = OpT({Q[0] * Q[1], Q[1]}, &ANNT);
    std::vector<FFVar> const dFb = DAG.BAD(FT, Q), dFf = DAG.FAD(FT, Q);
    std::string const ops = operations(DAG, dFb);
    if (ops.find("GradMLP") == std::string::npos || ops.find("TanMLP") != std::string::npos)
      return {false, "operations: " + ops};
    for (auto const& q : points)
      if (!near(evaluate(DAG, dFb, Q, q), evaluate(DAG, dFf, Q, q))) return {false, "derivatives differ"};
    return {true, ""};
  });

  run("FFTanMLP defined directly, and in a copy of the DAG", []() -> Result {
    mc::MLP<I> ANN;
    network(ANN, false);
    std::vector<FFVar> X2(2), D2(2), T2(2);
    FFGraph DAG2;
    {
      FFGraph DAG;
      std::vector<FFVar> X = DAG.add_vars(2), D = DAG.add_vars(2);
      mc::FFTanMLP<I> OpTan;
      std::vector<FFVar> const T = OpTan(X, D, &ANN);
      std::vector<FFVar> XD = X;
      XD.insert(XD.end(), D.begin(), D.end());
      std::vector<double> J(4);
      ANN.grad(points[0].data(), J.data());
      if (!near(evaluate(DAG, T, XD, {points[0][0], points[0][1], 0.5, -2.}), {J[0] * 0.5 - J[1] * 2., J[2] * 0.5 - J[3] * 2.}))
        return {false, "values"};
      DAG2.insert(&DAG, 2, X.data(), X2.data());
      DAG2.insert(&DAG, 2, D.data(), D2.data());
      DAG2.insert(&DAG, 2, T.data(), T2.data());
    }
    // the first DAG is gone: the copy holds its own network
    std::vector<FFVar> XD2 = X2;
    XD2.insert(XD2.end(), D2.begin(), D2.end());
    std::vector<double> J(4);
    ANN.grad(points[1].data(), J.data());
    return {near(evaluate(DAG2, T2, XD2, {points[1][0], points[1][1], 1., 3.}), {J[0] + J[1] * 3., J[2] + J[3] * 3.}), "values in the copy"};
  });

  run("interval and McCormick evaluation of the derivatives enclose their values", []() -> Result {
    mc::MLP<I> ANNT;
    network(ANNT, true);
    FFGraph DAG;
    std::vector<FFVar> Q = DAG.add_vars(2);
    mc::FFMLP<I> OpT;
    OpT.options.RELAX = {OpT.options.MC};
    std::vector<FFVar> const FT = OpT({Q[0] + Q[1], Q[0] * Q[1]}, &ANNT);
    std::vector<FFVar> const dFT = DAG.FAD(FT, Q);
    std::vector<I> const IQ = {I(-0.5, 1.), I(0.2, 0.9)};
    std::vector<I> IdF(dFT.size());
    DAG.eval(dFT, IdF, Q, IQ);
    std::vector<MC> MCdF(dFT.size());
    for (double const a : {0.1, 0.5, 0.9})
      for (double const b : {0.2, 0.7})
      {
        std::vector<double> const q = {IQ[0].l() + a * (IQ[0].u() - IQ[0].l()), IQ[1].l() + b * (IQ[1].u() - IQ[1].l())};
        std::vector<double> const d = evaluate(DAG, dFT, Q, q);
        std::vector<MC> const MCQ = {MC(IQ[0], q[0]), MC(IQ[1], q[1])};
        DAG.eval(dFT, MCdF, Q, MCQ);
        for (size_t k = 0; k < d.size(); ++k)
        {
          if (d[k] < IdF[k].l() - 1e-12 || d[k] > IdF[k].u() + 1e-12) return {false, "interval bounds"};
          if (d[k] < MCdF[k].cv() - 1e-12 || d[k] > MCdF[k].cc() + 1e-12) return {false, "McCormick relaxations"};
        }
      }
    return {true, ""};
  });

  run("a second differentiation is refused, as for FFGradMLP", []() -> Result {
    mc::MLP<I> ANNT;
    network(ANNT, true);
    FFGraph DAG;
    std::vector<FFVar> Q = DAG.add_vars(2);
    mc::FFMLP<I> OpT;
    std::vector<FFVar> const FT  = OpT(Q, &ANNT);
    std::vector<FFVar> const dFT = DAG.FAD(FT, Q);
    bool refused = false;
    try
    {
      std::vector<FFVar> const d2FT = DAG.FAD(dFT, Q);
    }
    catch (mc::FFBase::Exceptions& e)
    {
      refused = (e.ierr() == mc::FFBase::Exceptions::EXTERN);
    }
    std::vector<double> J(4);
    ANNT.grad(points[2].data(), J.data());
    return {refused && near(evaluate(DAG, dFT, Q, points[2]), J), "not refused, or DAG not usable afterwards"};
  });

  std::printf("\n=== Results: %d/%d passed ===\n", test::passed, test::total);
  return test::failed ? 1 : 0;
}
