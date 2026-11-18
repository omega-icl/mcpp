// =============================================================================
// test_MLP_adjoint.cpp
// Unit tests for the adjoint of multilayer perceptrons: mc::MLP::adj, the
// external operation mc::FFAdjMLP, and the backward symbolic differentiation
// of mc::FFMLP through the adjoint hook (FFOp::has_adjoint / FFOp::adjoint)
// when mc::MLP::Options::ADJOINT is set, against the Jacobian operation
// mc::FFGradMLP that it replaces.
//
// Build (assuming MC++ headers are on the include path):
//   g++ -std=c++17 -O2 -I<path-to-mc++> -o test_MLP_adjoint test_MLP_adjoint.cpp
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
network(mc::MLP<I>& ANN, bool const adjoint, bool const tangent = false)
{
  ANN.append_data({{0.1, 0.8, -0.5}, {-0.3, 0.4, 0.9}, {0.2, -0.7, 0.6}}, ANN.TANH, true);
  ANN.append_data({{0.05, 1.2, -0.4, 0.3}, {-0.1, -0.6, 0.7, 0.5}}, ANN.SIGMOID);
  ANN.append_data({{0.3, 1.5, -2.0}, {-0.2, 0.5, 0.25}}, ANN.LINEAR);
  ANN.options.ADJOINT = adjoint;
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
  std::printf("test_MLP_adjoint: skipped with MC__DAG_FADIFF or MC__DAG_BADIFF, the hooks belong to the default engines\n");
  return 0;
#endif

  using test::near;
  using test::Result;
  using test::run;

  std::printf("Adjoint of multilayer perceptrons\n");

  run("MLP::adj against the Jacobian of MLP::grad", []() -> Result {
    mc::MLP<I> ANN;
    network(ANN, false);
    for (auto const& x : points)
    {
      std::vector<double> J(4), w = {0.5, -2.}, a(2), y(2), yref(2);
      ANN.grad(x.data(), J.data(), yref.data());
      ANN.adj(x.data(), w.data(), a.data(), y.data());
      if (!near(y, yref) || !near(a, {J[0] * w[0] + J[2] * w[1], J[1] * w[0] + J[3] * w[1]}))
        return {false, "adjoint or value"};
      // a second call, and without the values
      ANN.adj(x.data(), w.data(), a.data());
      if (!near(a, {J[0] * w[0] + J[2] * w[1], J[1] * w[0] + J[3] * w[1]})) return {false, "second call"};
    }
    return {true, ""};
  });

  run("MLP::adj with ReLU, exponential-form sigmoid and tanh output activations", []() -> Result {
    mc::MLP<I> ANN;
    network2(ANN);
    for (auto const& x : points)
    {
      std::vector<double> J(4), w = {0.5, -2.}, a(2), y(2), yref(2);
      ANN.grad(x.data(), J.data(), yref.data());
      ANN.adj(x.data(), w.data(), a.data(), y.data());
      if (!near(y, yref) || !near(a, {J[0] * w[0] + J[2] * w[1], J[1] * w[0] + J[3] * w[1]}))
        return {false, "adjoint or value"};
      // enclosure in interval arithmetic
      std::vector<I> const IX = {x[0] + I(-0.1, 0.1), x[1] + I(-0.1, 0.1)}, IW = {I(w[0]), I(w[1])};
      std::vector<I> IA(2), IY(2);
      ANN.adj(IX.data(), IW.data(), IA.data(), IY.data());
      for (int k = 0; k < 2; ++k)
        if (a[k] < IA[k].l() - 1e-12 || a[k] > IA[k].u() + 1e-12 || y[k] < IY[k].l() - 1e-12 || y[k] > IY[k].u() + 1e-12)
          return {false, "interval bounds"};
    }
    return {true, ""};
  });

  run("BAD: same derivatives with Options::ADJOINT, through FFAdjMLP operations", []() -> Result {
    mc::MLP<I> ANN, ANNA;
    network(ANN, false);
    network(ANNA, true);
    FFGraph DAG;
    std::vector<FFVar> Q = DAG.add_vars(3);
    // inputs as expressions of the independents, outputs in further expressions
    std::vector<FFVar> const X = {Q[0] * Q[1], Q[1] - 2. * Q[2]};
    mc::FFMLP<I> Op, OpA;
    std::vector<FFVar> const Y = Op(X, &ANN), YA = OpA(X, &ANNA);
    std::vector<FFVar> const F = {Y[0] * Y[1] + Q[2], exp(Y[1])}, FA = {YA[0] * YA[1] + Q[2], exp(YA[1])};
    std::vector<FFVar> const dF = DAG.BAD(F, Q), dFA = DAG.BAD(FA, Q);
    std::string const ops = operations(DAG, dF), opsA = operations(DAG, dFA);
    if (ops.find("GradMLP") == std::string::npos || ops.find("AdjMLP") != std::string::npos)
      return {false, "operations without the option: " + ops};
    if (opsA.find("AdjMLP") == std::string::npos || opsA.find("GradMLP") != std::string::npos)
      return {false, "operations with the option: " + opsA};
    for (auto const& p : points)
    {
      std::vector<double> const q = {p[0], p[1], 0.6};
      if (!near(evaluate(DAG, dFA, Q, q), evaluate(DAG, dF, Q, q))) return {false, "derivatives differ"};
    }
    return {true, ""};
  });

  run("SBAD and DBAD: same derivatives with Options::ADJOINT", []() -> Result {
    mc::MLP<I> ANN, ANNA;
    network(ANN, false);
    network(ANNA, true);
    FFGraph DAG;
    std::vector<FFVar> Q = DAG.add_vars(3), W = DAG.add_vars(2);
    std::vector<FFVar> const X = {Q[0] * Q[1], sqr(Q[1])};  // Q[2] does not reach the network
    mc::FFMLP<I> Op, OpA;
    std::vector<FFVar> const F = Op(X, &ANN), FA = OpA(X, &ANNA);
    auto [r0, c0, v0] = DAG.SBAD(F, Q);
    auto [r1, c1, v1] = DAG.SBAD(FA, Q);
    if (r0 != r1 || c0 != c1 || r1.size() != 4) return {false, "sparsity patterns differ"};
    std::vector<FFVar> const dir = DAG.DBAD(F, W, Q), dirA = DAG.DBAD(FA, W, Q);
    std::vector<FFVar> QW = Q;
    QW.insert(QW.end(), W.begin(), W.end());
    for (auto const& p : points)
    {
      std::vector<double> const q = {p[0], p[1], 0.6}, qw = {p[0], p[1], 0.6, 0.5, -1.};
      if (!near(evaluate(DAG, v1, Q, q), evaluate(DAG, v0, Q, q))) return {false, "sparse derivatives differ"};
      if (!near(evaluate(DAG, dirA, QW, qw), evaluate(DAG, dir, QW, qw))) return {false, "weighted derivatives differ"};
    }
    return {true, ""};
  });

  run("the Jacobian operation is kept where the sweeps through the network are many", []() -> Result {
    mc::MLP<I> ANN, ANNA;
    network(ANN, false);
    network(ANNA, true);
    FFGraph DAG;
    std::vector<FFVar> Q = DAG.add_vars(2), W = DAG.add_vars(5);
    mc::FFMLP<I> Op, OpA;
    std::vector<FFVar> const X = {Q[0] * Q[1], Q[1]};
    std::vector<FFVar> const Y = OpA(X, &ANNA), Z = Op(X, &ANN);
    // with 2 inputs the adjoint is used for up to 4 sweeps: here 4, then 5 dependents
    std::vector<FFVar> const F5 = {Y[0], Y[1], Y[0] * Y[1], exp(Y[0]), sqr(Y[1]) + Q[0]}, F4(F5.begin(), F5.begin() + 4);
    std::vector<FFVar> const G5 = {Z[0], Z[1], Z[0] * Z[1], exp(Z[0]), sqr(Z[1]) + Q[0]};
    std::vector<FFVar> const dF4 = DAG.BAD(F4, Q), dF5 = DAG.BAD(F5, Q), dG5 = DAG.BAD(G5, Q), dirF5 = DAG.DBAD(F5, W, Q);
    std::string ops = operations(DAG, dF4);
    if (ops.find("AdjMLP") == std::string::npos || ops.find("GradMLP") != std::string::npos) return {false, "4 sweeps: " + ops};
    ops = operations(DAG, dF5);
    if (ops.find("GradMLP") == std::string::npos || ops.find("AdjMLP") != std::string::npos) return {false, "5 sweeps: " + ops};
    ops = operations(DAG, dirF5);
    if (ops.find("AdjMLP") == std::string::npos || ops.find("GradMLP") != std::string::npos) return {false, "one weighted sum: " + ops};
    for (auto const& q : points)
      if (!near(evaluate(DAG, dF5, Q, q), evaluate(DAG, dG5, Q, q))) return {false, "derivatives differ"};
    return {true, ""};
  });

  run("FAD: the Jacobian operation with Options::ADJOINT alone, FFTanMLP with both options", []() -> Result {
    mc::MLP<I> ANNA, ANNB;
    network(ANNA, true);
    network(ANNB, true, true);
    FFGraph DAG;
    std::vector<FFVar> Q = DAG.add_vars(2);
    mc::FFMLP<I> OpA, OpB;
    std::vector<FFVar> const FA = OpA({Q[0] * Q[1], Q[1]}, &ANNA), FB = OpB({Q[0] * Q[1], Q[1]}, &ANNB);
    std::vector<FFVar> const dFAf = DAG.FAD(FA, Q), dFAb = DAG.BAD(FA, Q), dFBf = DAG.FAD(FB, Q), dFBb = DAG.BAD(FB, Q);
    std::string ops = operations(DAG, dFAf);
    if (ops.find("GradMLP") == std::string::npos || ops.find("AdjMLP") != std::string::npos)
      return {false, "forward operations with ADJOINT alone: " + ops};
    ops = operations(DAG, dFBf);
    if (ops.find("TanMLP") == std::string::npos || ops.find("GradMLP") != std::string::npos || ops.find("AdjMLP") != std::string::npos)
      return {false, "forward operations with both options: " + ops};
    ops = operations(DAG, dFBb);
    if (ops.find("AdjMLP") == std::string::npos || ops.find("GradMLP") != std::string::npos || ops.find("TanMLP") != std::string::npos)
      return {false, "backward operations with both options: " + ops};
    for (auto const& q : points)
    {
      std::vector<double> const ref = evaluate(DAG, dFAf, Q, q);
      if (!near(evaluate(DAG, dFAb, Q, q), ref) || !near(evaluate(DAG, dFBf, Q, q), ref) || !near(evaluate(DAG, dFBb, Q, q), ref))
        return {false, "derivatives differ"};
    }
    return {true, ""};
  });

  run("FFAdjMLP defined directly, and in a copy of the DAG", []() -> Result {
    mc::MLP<I> ANN;
    network(ANN, false);
    std::vector<FFVar> X2(2), W2(2), A2(2);
    FFGraph DAG2;
    {
      FFGraph DAG;
      std::vector<FFVar> X = DAG.add_vars(2), W = DAG.add_vars(2);
      mc::FFAdjMLP<I> OpAdj;
      std::vector<FFVar> const A = OpAdj(X, W, &ANN);
      std::vector<FFVar> XW = X;
      XW.insert(XW.end(), W.begin(), W.end());
      std::vector<double> J(4);
      ANN.grad(points[0].data(), J.data());
      if (!near(evaluate(DAG, A, XW, {points[0][0], points[0][1], 0.5, -2.}), {J[0] * 0.5 - J[2] * 2., J[1] * 0.5 - J[3] * 2.}))
        return {false, "values"};
      DAG2.insert(&DAG, 2, X.data(), X2.data());
      DAG2.insert(&DAG, 2, W.data(), W2.data());
      DAG2.insert(&DAG, 2, A.data(), A2.data());
    }
    // the first DAG is gone: the copy holds its own network
    std::vector<FFVar> XW2 = X2;
    XW2.insert(XW2.end(), W2.begin(), W2.end());
    std::vector<double> J(4);
    ANN.grad(points[1].data(), J.data());
    return {near(evaluate(DAG2, A2, XW2, {points[1][0], points[1][1], 1., 3.}), {J[0] + J[2] * 3., J[1] + J[3] * 3.}), "values in the copy"};
  });

  run("interval and McCormick evaluation of the derivatives enclose their values", []() -> Result {
    mc::MLP<I> ANNA;
    network(ANNA, true);
    FFGraph DAG;
    std::vector<FFVar> Q = DAG.add_vars(2);
    mc::FFMLP<I> OpA;
    OpA.options.RELAX = {OpA.options.MC};
    std::vector<FFVar> const FA = OpA({Q[0] + Q[1], Q[0] * Q[1]}, &ANNA);
    std::vector<FFVar> const dFA = DAG.BAD(FA, Q);
    std::vector<I> const IQ = {I(-0.5, 1.), I(0.2, 0.9)};
    std::vector<I> IdF(dFA.size());
    DAG.eval(dFA, IdF, Q, IQ);
    std::vector<MC> MCdF(dFA.size());
    for (double const a : {0.1, 0.5, 0.9})
      for (double const b : {0.2, 0.7})
      {
        std::vector<double> const q = {IQ[0].l() + a * (IQ[0].u() - IQ[0].l()), IQ[1].l() + b * (IQ[1].u() - IQ[1].l())};
        std::vector<double> const d = evaluate(DAG, dFA, Q, q);
        std::vector<MC> const MCQ = {MC(IQ[0], q[0]), MC(IQ[1], q[1])};
        DAG.eval(dFA, MCdF, Q, MCQ);
        for (size_t k = 0; k < d.size(); ++k)
        {
          if (d[k] < IdF[k].l() - 1e-12 || d[k] > IdF[k].u() + 1e-12) return {false, "interval bounds"};
          if (d[k] < MCdF[k].cv() - 1e-12 || d[k] > MCdF[k].cc() + 1e-12) return {false, "McCormick relaxations"};
        }
      }
    // not degenerate: the bounds have some width
    return {IdF[0].u() > IdF[0].l(), "degenerate bounds"};
  });

  run("a second differentiation is refused, as for FFGradMLP", []() -> Result {
    mc::MLP<I> ANNA;
    network(ANNA, true);
    FFGraph DAG;
    std::vector<FFVar> Q = DAG.add_vars(2);
    mc::FFMLP<I> OpA;
    std::vector<FFVar> const FA  = OpA(Q, &ANNA);
    std::vector<FFVar> const dFA = DAG.BAD(FA, Q);
    int refused = 0;
    for (int mode = 0; mode < 2; ++mode) try
      {
        std::vector<FFVar> const d2FA = (mode ? DAG.BAD(dFA, Q) : DAG.FAD(dFA, Q));
      }
      catch (mc::FFBase::Exceptions& e)
      {
        if (e.ierr() == mc::FFBase::Exceptions::EXTERN) ++refused;
      }
    std::vector<double> J(4);
    ANNA.grad(points[2].data(), J.data());
    return {refused == 2 && near(evaluate(DAG, dFA, Q, points[2]), J), "not refused, or DAG not usable afterwards"};
  });

  std::printf("\n=== Results: %d/%d passed ===\n", test::passed, test::total);
  return test::failed ? 1 : 0;
}
