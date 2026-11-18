// =============================================================================
// test_MLP_ident.cpp
// Unit tests for the identification of MLP operations in a DAG. An operation
// that holds a copy of a network (policy COPY, the default) is identified by
// the stamp of the network it copied, mc::MLP::uid(): a network changed after
// its insertion, or another network created at the address of the first, is
// another operation; an unchanged network, or a copy of it, is the same one.
// An operation that refers to the network itself (SHALLOW, TRANSFER) is
// identified by its address. The options are not part of the identification.
//
// Build (assuming MC++ headers are on the include path):
//   g++ -std=c++17 -O2 -I<path-to-mc++> -o test_MLP_ident test_MLP_ident.cpp
//
// Each test prints PASS or FAIL and a summary is printed at the end.
// The process exits with code 0 if every test passes, 1 otherwise.
// =============================================================================

#include <cmath>
#include <cstdio>
#include <functional>
#include <iostream>
#include <optional>
#include <set>
#include <stdexcept>
#include <string>
#include <thread>
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
// 2 inputs, 3 tanh neurons, nout linear outputs; the weights follow from seed
static void
network(mc::MLP<I>& ANN, size_t const nout, double const seed)
{
  std::vector<std::vector<double>> layer1(3, std::vector<double>(3)), layer2(nout, std::vector<double>(4));
  for (size_t i = 0; i < 3; ++i)
    for (size_t j = 0; j < 3; ++j) layer1[i][j] = std::cos(seed + 2. * i + 3. * j);
  for (size_t i = 0; i < nout; ++i)
    for (size_t j = 0; j < 4; ++j) layer2[i][j] = std::sin(seed + 5. * i + 7. * j);
  ANN.append_data(layer1, ANN.TANH, true);
  ANN.append_data(layer2, ANN.LINEAR);
}

// the network evaluated outside any DAG
static std::vector<double>
direct(mc::MLP<I> const& ANN, std::vector<double> const& x)
{
  mc::MLP<I> copy(ANN);
  std::vector<double> y(copy.nout());
  copy.eval(x.data(), y.data());
  return y;
}

static std::vector<double>
jacobian(mc::MLP<I> const& ANN, std::vector<double> const& x)
{
  mc::MLP<I> copy(ANN);
  std::vector<double> J(copy.nout() * copy.nin());
  copy.grad(x.data(), J.data());
  return J;
}

static std::vector<double>
evaluate(FFGraph& DAG, std::vector<FFVar> const& F, std::vector<FFVar> const& X, std::vector<double> const& x)
{
  std::vector<double> f(F.size());
  if (!F.empty()) DAG.eval(F, f, X, x);
  return f;
}

// number of external operations that F depends on
static unsigned
nops(FFGraph& DAG, std::vector<FFVar> const& F)
{
  unsigned n = 0;
  auto sg = DAG.subgraph(F);
  for (auto const& op : sg.l_op)
    if (op->type >= mc::FFOp::EXTERN) ++n;
  return n;
}

static bool
same(std::vector<FFVar> const& A, std::vector<FFVar> const& B)
{
  if (A.size() != B.size()) return false;
  for (size_t i = 0; i < A.size(); ++i)
    if (A[i].id() != B[i].id()) return false;
  return true;
}

static std::vector<double> const x0 = {0.3, -0.7};

int
main()
{
  using namespace test;
  std::printf("Identification of multilayer perceptrons in a DAG\n");

  run("MLP::uid: kept by copies and assignment, drawn again when the network changes", []() -> Result {
    mc::MLP<I> A, B;
    if (!A.uid() || A.uid() == B.uid()) return {false, "two new networks"};
    network(A, 2, 1.);
    size_t const uidA = A.uid();
    mc::MLP<I> C(A);
    B = A;
    if (C.uid() != uidA || B.uid() != uidA) return {false, "copy or assignment"};
    C.options.SIG2EXP = !C.options.SIG2EXP;
    if (C.uid() != uidA) return {false, "options"};
    std::set<size_t> seen = {uidA};
    // every change draws a stamp never seen, and leaves the copies alone
    A.append_data(std::vector<double>{1., 2., 3.}, A.LINEAR);
    if (!seen.insert(A.uid()).second) return {false, "append_data"};
    std::vector<std::pair<std::vector<std::vector<double>>, int>> const data = {{{{1., 2., 3.}}, mc::MLP<I>::TANH}};
    if (!A.set_data(data) || !seen.insert(A.uid()).second) return {false, "set_data"};
    A.reset_data();
    if (!seen.insert(A.uid()).second) return {false, "reset_data"};
    network(A, 2, 1.);  // the same weights as before: still another stamp
    if (!seen.insert(A.uid()).second) return {false, "same weights set again"};
    if (C.uid() != uidA || B.uid() != uidA) return {false, "copies changed"};
    // a refused change leaves the stamp
    size_t const uidB = B.uid();
    if (B.append_data(std::vector<double>{1., 2.}, B.LINEAR) || B.uid() != uidB) return {false, "refused append_data"};
    // stamps drawn from several threads are distinct
    std::vector<std::vector<size_t>> drawn(4);
    std::vector<std::thread> threads;
    for (size_t t = 0; t < 4; ++t)
      threads.emplace_back([&drawn, t]() {
        for (size_t k = 0; k < 500; ++k) drawn[t].push_back(mc::MLP<I>().uid());
      });
    for (auto& th : threads) th.join();
    for (auto const& d : drawn)
      for (size_t u : d)
        if (!seen.insert(u).second) return {false, "stamps from several threads"};
    return {true, ""};
  });

  run("COPY: a network changed after its insertion is another operation", []() -> Result {
    FFGraph DAG;
    std::vector<FFVar> X = DAG.add_vars(2);
    mc::MLP<I> ANN;
    network(ANN, 2, 1.);
    std::vector<double> const y1 = direct(ANN, x0);
    mc::FFMLP<I> Op;
    std::vector<FFVar> const Y1 = Op(X, &ANN);
    // other weights, then a further layer and another number of outputs
    network(ANN, 2, 2.);
    std::vector<double> const y2 = direct(ANN, x0);
    std::vector<FFVar> const Y2 = Op(X, &ANN);
    ANN.append_data(std::vector<double>{0.5, -1., 2.}, ANN.TANH);
    std::vector<double> const y3 = direct(ANN, x0);
    std::vector<FFVar> const Y3 = Op(X, &ANN);
    if (y1 == y2 || Y3.size() != 1) return {false, "test networks"};
    if (same(Y1, Y2)) return {false, "the changed network was taken for the first"};
    if (evaluate(DAG, Y1, X, x0) != y1) return {false, "values of the first"};
    if (evaluate(DAG, Y2, X, x0) != y2) return {false, "values after other weights"};
    if (evaluate(DAG, Y3, X, x0) != y3) return {false, "values after a further layer"};
    std::vector<FFVar> all = Y1;
    all.insert(all.end(), Y2.begin(), Y2.end());
    all.insert(all.end(), Y3.begin(), Y3.end());
    return {nops(DAG, all) == 3, "3 operations expected"};
  });

  run("COPY: networks succeeding each other at the same address are distinct operations", []() -> Result {
    FFGraph DAG;
    std::vector<FFVar> X = DAG.add_vars(2);
    std::optional<mc::MLP<I>> slot;
    void const* address = nullptr;
    std::vector<std::vector<FFVar>> Y;
    std::vector<std::vector<double>> y;
    for (size_t nout : {1, 4, 2})
    {
      slot.emplace();  // the previous network is destroyed, the next built in its place
      if (address && address != &*slot) return {false, "test: the address changed"};
      address = &*slot;
      network(*slot, nout, 1.);
      mc::FFMLP<I> Op;
      Y.push_back(Op(X, &*slot));
      y.push_back(direct(*slot, x0));
    }
    slot.reset();  // none of the networks exists any more
    for (size_t k = 0; k < Y.size(); ++k)
      if (Y[k].size() != y[k].size() || evaluate(DAG, Y[k], X, x0) != y[k])
        return {false, "network " + std::to_string(k)};
    return {true, ""};
  });

  run("COPY: an unchanged network, a copy of it and its outputs one by one are the same operation", []() -> Result {
    FFGraph DAG;
    std::vector<FFVar> X = DAG.add_vars(2), Z = DAG.add_vars(2);
    mc::MLP<I> ANN;
    network(ANN, 3, 1.);
    mc::FFMLP<I> Op, Op2;
    std::vector<FFVar> const Y = Op(X, &ANN);
    if (!same(Y, Op(X, &ANN)) || !same(Y, Op2(X, &ANN))) return {false, "inserted again"};
    std::vector<FFVar> const Yi = {Op2(0, X, &ANN), Op2(1, X, &ANN), Op2(2, X, &ANN)};
    if (!same(Y, Yi) || nops(DAG, Yi) != 1) return {false, "outputs one by one"};
    mc::MLP<I> const COPY(ANN);
    mc::MLP<I> ASSIGNED;
    ASSIGNED = ANN;
    mc::MLP<I> copy(COPY);
    if (!same(Y, Op(X, &copy)) || !same(Y, Op(X, &ASSIGNED))) return {false, "copy of the network"};
    // the options are not part of the identification: the operation found is the one inserted first
    ANN.options.SIG2EXP = !ANN.options.SIG2EXP;
    if (!same(Y, Op(X, &ANN))) return {false, "options changed"};
    // other operands: another operation of the same network
    std::vector<FFVar> const W = Op(Z, &ANN);
    std::vector<FFVar> both = Y;
    both.insert(both.end(), W.begin(), W.end());
    if (same(Y, W) || nops(DAG, both) != 2) return {false, "other operands"};
    std::vector<FFVar> XZ = X;
    XZ.insert(XZ.end(), Z.begin(), Z.end());
    std::vector<double> const z0 = {-1.2, 0.4};
    return {evaluate(DAG, W, XZ, {x0[0], x0[1], z0[0], z0[1]}) == direct(ANN, z0), "values"};
  });

  run("SHALLOW and TRANSFER: an operation is identified by its network, not by the content", []() -> Result {
    FFGraph DAG;
    std::vector<FFVar> X = DAG.add_vars(2);
    mc::MLP<I> ANN;
    network(ANN, 2, 1.);
    mc::MLP<I> TWIN(ANN);  // same content, same stamp, another object
    mc::FFMLP<I> Op;
    std::vector<FFVar> const Y = Op(X, &ANN, Op.SHALLOW), YT = Op(X, &TWIN, Op.SHALLOW);
    if (same(Y, YT)) return {false, "two objects taken for one"};
    if (!same(Y, Op(X, &ANN, Op.SHALLOW))) return {false, "inserted again"};
    // a copy in the same DAG is a third operation
    std::vector<FFVar> const YC = Op(X, &ANN);
    if (same(Y, YC) || same(YT, YC)) return {false, "copy taken for the object"};
    // the operation follows its object: a change is seen without inserting again, and inserting
    // again finds the same operation; the twin and the copy keep the first network
    std::vector<double> const y1 = direct(ANN, x0);
    network(ANN, 2, 2.);
    std::vector<double> const y2 = direct(ANN, x0);
    if (evaluate(DAG, Y, X, x0) != y2) return {false, "change of the object not seen"};
    if (!same(Y, Op(X, &ANN, Op.SHALLOW))) return {false, "inserted again after the change"};
    if (evaluate(DAG, YT, X, x0) != y1 || evaluate(DAG, YC, X, x0) != y1) return {false, "twin or copy changed"};
    // the copy of the changed network is a fourth operation
    std::vector<FFVar> const YC2 = Op(X, &ANN);
    if (same(YC, YC2) || evaluate(DAG, YC2, X, x0) != y2) return {false, "copy after the change"};
    // ownership handed to the DAG: identified by the object too
    mc::MLP<I>* pANN = new mc::MLP<I>;
    network(*pANN, 2, 3.);
    std::vector<double> const y3 = direct(*pANN, x0);
    std::vector<FFVar> const YO = Op(X, pANN, Op.TRANSFER);
    if (!same(YO, Op(X, pANN, Op.SHALLOW))) return {false, "transferred network inserted again"};
    return {evaluate(DAG, YO, X, x0) == y3, "values of the transferred network"};
  });

  run("derivatives: Jacobian, tangent and adjoint operations of a network and of its changed self", []() -> Result {
    for (bool const hooks : {false, true})
    {
      FFGraph DAG;
      std::vector<FFVar> X = DAG.add_vars(2);
      mc::MLP<I> ANN;
      network(ANN, 2, 1.);
      ANN.options.TANGENT = ANN.options.ADJOINT = hooks;
      mc::MLP<I> const OLD(ANN);
      mc::FFMLP<I> Op;
      std::vector<FFVar> const Y1 = Op(X, &ANN);
      network(ANN, 2, 2.);
      std::vector<FFVar> const Y2 = Op(X, &ANN);
      std::vector<double> const J1 = jacobian(OLD, x0), J2 = jacobian(ANN, x0);
      if (near(J1, J2, 1e-3)) return {false, "test networks"};
      // differentiated together, so that the derivative operations of both are in the DAG at once
      std::vector<FFVar> Y12 = Y1;
      Y12.insert(Y12.end(), Y2.begin(), Y2.end());
      std::vector<double> J12 = J1;
      J12.insert(J12.end(), J2.begin(), J2.end());
      if (!near(evaluate(DAG, DAG.FAD(Y12, X), X, x0), J12)) return {false, hooks ? "FAD with the hooks" : "FAD"};
      if (!near(evaluate(DAG, DAG.BAD(Y12, X), X, x0), J12)) return {false, hooks ? "BAD with the hooks" : "BAD"};
      // the derivative operations of the first network are those of a copy of it
      mc::MLP<I> old(OLD);
      if (!same(DAG.FAD(Y1, X), DAG.FAD(Op(X, &old), X))) return {false, "derivatives of a copy"};
    }
    // a network that its operation refers to, changed after a first differentiation: the derivative
    // operations hold copies of the network, and a second differentiation makes those of the changed one
    for (bool const hooks : {false, true})
    {
      FFGraph DAG;
      std::vector<FFVar> X = DAG.add_vars(2);
      mc::MLP<I> ANN;
      network(ANN, 2, 1.);
      ANN.options.TANGENT = ANN.options.ADJOINT = hooks;
      mc::FFMLP<I> Op;
      std::vector<FFVar> const Y = Op(X, &ANN, Op.SHALLOW);
      if (!near(evaluate(DAG, DAG.FAD(Y, X), X, x0), jacobian(ANN, x0))) return {false, "SHALLOW: FAD"};
      if (!near(evaluate(DAG, DAG.BAD(Y, X), X, x0), jacobian(ANN, x0))) return {false, "SHALLOW: BAD"};
      network(ANN, 2, 2.);
      if (!near(evaluate(DAG, DAG.FAD(Y, X), X, x0), jacobian(ANN, x0))) return {false, "SHALLOW: FAD after the change"};
      if (!near(evaluate(DAG, DAG.BAD(Y, X), X, x0), jacobian(ANN, x0))) return {false, "SHALLOW: BAD after the change"};
    }
    return {true, ""};
  });

  run("a copy of the DAG keeps the operations distinct, and finds the networks again", []() -> Result {
    mc::MLP<I> ANN;
    network(ANN, 2, 1.);
    mc::MLP<I> OLD(ANN);
    std::vector<double> const y1 = direct(ANN, x0);
    std::vector<FFVar> X2(2), Y12(4);
    FFGraph DAG2;
    {
      FFGraph DAG;
      std::vector<FFVar> X = DAG.add_vars(2);
      mc::FFMLP<I> Op;
      std::vector<FFVar> Y = Op(X, &ANN);
      network(ANN, 2, 2.);
      std::vector<FFVar> const Y2 = Op(X, &ANN);
      Y.insert(Y.end(), Y2.begin(), Y2.end());
      DAG2.insert(&DAG, 2, X.data(), X2.data());
      DAG2.insert(&DAG, 4, Y.data(), Y12.data());
    }
    std::vector<double> y12 = y1;
    std::vector<double> const y2 = direct(ANN, x0);
    y12.insert(y12.end(), y2.begin(), y2.end());
    if (evaluate(DAG2, Y12, X2, x0) != y12) return {false, "values in the copy"};
    mc::FFMLP<I> Op;
    if (!same(Op(X2, &OLD), {Y12[0], Y12[1]})) return {false, "first network not found in the copy"};
    return {same(Op(X2, &ANN), {Y12[2], Y12[3]}), "second network not found in the copy"};
  });

  std::printf("\n=== Results: %d/%d passed ===\n", test::passed, test::total);
  return test::failed ? 1 : 0;
}
