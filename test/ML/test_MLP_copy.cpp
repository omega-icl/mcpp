// =============================================================================
// test_MLP_copy.cpp
// Unit tests for the copy construction and the assignment of mc::MLP: the
// copy holds the network and the options of the original, builds its own
// expression tree, and does not depend on the original afterwards.
//
// Build (assuming MC++ headers are on the include path):
//   g++ -std=c++17 -O2 -I<path-to-mc++> -o test_MLP_copy test_MLP_copy.cpp
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
// Two networks of different sizes
// -----------------------------------------------------------------------------
typedef std::vector<std::vector<double>> layer;

// 2 inputs, 3 tanh neurons, 2 linear outputs
static void
network1(mc::MLP<I>& ANN)
{
  ANN.append_data(layer{{0.1, 0.8, -0.5}, {-0.3, 0.4, 0.9}, {0.2, -0.7, 0.6}}, ANN.TANH, true);
  ANN.append_data(layer{{0.3, 1.5, -2.0, 0.4}, {-0.2, 0.5, 0.25, -1.}}, ANN.LINEAR);
}

// 3 inputs, 2 sigmoid neurons, 4 ReLU neurons, 1 linear output
static void
network2(mc::MLP<I>& ANN)
{
  ANN.append_data(layer{{0.05, 1.2, -0.4, 0.3}, {-0.1, -0.6, 0.7, 0.5}}, ANN.SIGMOID, true);
  ANN.append_data(layer{{0.1, 0.8, -0.5}, {-0.3, 0.4, 0.9}, {0.2, -0.7, 0.6}, {0., 1., 1.}}, ANN.RELU);
  ANN.append_data(layer{{0.3, 1.5, -2.0, 0.4, 0.7}}, ANN.LINEAR);
}

// Everything an MLP evaluates at a point: values and Jacobian in double, interval bounds and
// McCormick relaxations on a box around the point, and the values through its expression tree
static std::vector<double>
evaluate(mc::MLP<I>& ANN)
{
  size_t const nin = ANN.nin(), nout = ANN.nout();
  std::vector<double> x(nin), y(nout), J(nout * nin), res;
  for (size_t i = 0; i < nin; ++i) x[i] = 0.3 - 0.4 * i;
  ANN.eval(x.data(), y.data());
  res.insert(res.end(), y.begin(), y.end());
  ANN.grad(x.data(), J.data());
  res.insert(res.end(), J.begin(), J.end());
  std::vector<I> IX(nin), IY(nout);
  std::vector<MC> MCX(nin), MCY(nout);
  for (size_t i = 0; i < nin; ++i)
  {
    IX[i]  = I(x[i] - 0.5, x[i] + 0.5);
    MCX[i] = MC(IX[i], x[i]);
  }
  ANN.eval(IX.data(), IY.data());
  ANN.eval(MCX.data(), MCY.data());
  for (size_t k = 0; k < nout; ++k) res.insert(res.end(), {IY[k].l(), IY[k].u(), MCY[k].cv(), MCY[k].cc()});
  // through the expression tree of the MLP
  FFGraph* dag = ANN.dag();
  res.push_back(dag ? dag->Vars().size() : 0.);
  return res;
}

int
main()
{
  using test::near;
  using test::Result;
  using test::run;

  std::printf("Copy and assignment of multilayer perceptrons\n");

  run("copy construction: same network and options, independent of the original", []() -> Result {
    std::vector<double> ref;
    mc::MLP<I>* pANN = new mc::MLP<I>;
    network1(*pANN);
    pANN->options.TANGENT = true;
    pANN->dag();
    ref = evaluate(*pANN);
    mc::MLP<I> ANN2(*pANN);
    if (ANN2.nin() != 2 || ANN2.nout() != 2 || !ANN2.options.TANGENT) return {false, "sizes or options"};
    if (evaluate(ANN2) != ref) return {false, "evaluation of the copy"};
    // the original changes, then goes: the copy is unaffected
    network2(*pANN);
    if (evaluate(ANN2) != ref) return {false, "copy changed with the original"};
    delete pANN;
    return {evaluate(ANN2) == ref, "copy after the original is gone"};
  });

  run("assignment: replaces a network of another size, and its expression tree", []() -> Result {
    mc::MLP<I> ANN2;
    network2(ANN2);
    ANN2.dag();
    std::vector<double> const old = evaluate(ANN2);
    std::vector<double> ref;
    {
      mc::MLP<I> ANN1;
      network1(ANN1);
      ANN1.options.ADJOINT = true;
      ANN1.dag();
      ref  = evaluate(ANN1);
      ANN2 = ANN1;
      if (ANN2.nin() != 2 || ANN2.nout() != 2 || !ANN2.options.ADJOINT) return {false, "sizes or options"};
      if (evaluate(ANN2) != ref || ref == old) return {false, "evaluation after the assignment"};
      if (evaluate(ANN1) != ref) return {false, "the original changed"};
      // the original changes: the assigned MLP is unaffected
      network2(ANN1);
      if (evaluate(ANN2) != ref) return {false, "assigned MLP changed with the original"};
    }
    return {evaluate(ANN2) == ref, "assigned MLP after the original is gone"};
  });

  run("assignment to itself, from an empty MLP, and in a chain", []() -> Result {
    mc::MLP<I> ANN1, ANN2, ANN3, empty;
    network1(ANN1);
    network2(ANN2);
    std::vector<double> const ref = evaluate(ANN1);
    mc::MLP<I>& self = ANN1;
    ANN1 = self;
    if (evaluate(ANN1) != ref) return {false, "assignment to itself"};
    ANN3 = ANN2 = ANN1;
    if (evaluate(ANN2) != ref || evaluate(ANN3) != ref) return {false, "chain"};
    ANN2 = empty;
    if (ANN2.nin() != 0 || ANN2.nout() != 0 || ANN2.dag() != nullptr) return {false, "from an empty MLP"};
    return {evaluate(ANN3) == ref, "after the others changed"};
  });

  run("an assigned MLP in a DAG: operation and derivatives", []() -> Result {
    mc::MLP<I> ANN1, ANN2;
    network1(ANN1);
    network2(ANN2);
    ANN2 = ANN1;
    FFGraph DAG;
    std::vector<FFVar> Q = DAG.add_vars(2);
    mc::FFMLP<I> Op1, Op2;
    std::vector<FFVar> const F1 = Op1({Q[0] * Q[1], Q[1]}, &ANN1), F2 = Op2({Q[0] * Q[1], Q[1]}, &ANN2);
    std::vector<FFVar> const dF1 = DAG.FAD(F1, Q), dF2 = DAG.BAD(F2, Q);
    std::vector<double> const q = {0.3, -0.7};
    std::vector<double> f1(2), f2(2), d1(4), d2(4);
    DAG.eval(F1, f1, Q, q);
    DAG.eval(F2, f2, Q, q);
    DAG.eval(dF1, d1, Q, q);
    DAG.eval(dF2, d2, Q, q);
    return {f1 == f2 && near(d1, d2), "values or derivatives"};
  });

  std::printf("\n=== Results: %d/%d passed ===\n", test::passed, test::total);
  return test::failed ? 1 : 0;
}
