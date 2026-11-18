// =============================================================================
// DAG_dagext_home.cpp
// Unit tests for expression trees as external operations, mc::FFDAGEXT: their
// identification in a DAG and the home of an expression tree there. A tree
// that the DAG holds a copy of (policy COPY, the default) is identified by
// its stamp, mc::DAGEXT::uid(): another tree of the same shape, or the same
// object set to another tree, is another operation. The home of a tree in a
// DAG -- the tree, one copy per DAG, and the options of the operations on it
// -- is shared by every operation on that tree; its Jacobian is
// differentiated once, in a home that has the same options. The options may
// be changed after the insertion, and the operations follow.
//
// Build (assuming MC++ headers are on the include path):
//   g++ -std=c++17 -O2 -I<path-to-mc++> -o DAG_dagext_home DAG_dagext_home.cpp
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

#include "ffdagext.hpp"

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

typedef mc::DAGEXT<I> TREE;
typedef mc::FFDAGEXT<I> OPTREE;

// the operation that defines a variable, as an operation on an expression tree
static OPTREE*
opof(FFVar const& V)
{
  return static_cast<OPTREE*>(V.opdef().first);
}

// expression trees held by the external operations that F depends on
static std::set<TREE*>
trees(FFGraph& DAG, std::vector<FFVar> const& F)
{
  std::set<TREE*> all;
  auto sg = DAG.subgraph(F);
  for (auto const& op : sg.l_op)
    if (op->type >= mc::FFOp::EXTERN) all.insert(static_cast<OPTREE const*>(op)->object());
  return all;
}

static std::vector<double>
evaluate(FFGraph& DAG, std::vector<FFVar> const& F, std::vector<FFVar> const& X, std::vector<double> const& x)
{
  std::vector<double> f(F.size());
  if (!F.empty()) DAG.eval(F, f, X, x);
  return f;
}

static bool
same(std::vector<FFVar> const& A, std::vector<FFVar> const& B)
{
  if (A.size() != B.size()) return false;
  for (size_t i = 0; i < A.size(); ++i)
    if (A[i].id() != B[i].id()) return false;
  return true;
}

static std::vector<FFVar>
join(std::vector<FFVar> A, std::vector<FFVar> const& B)
{
  A.insert(A.end(), B.begin(), B.end());
  return A;
}

// polyhedral relaxation of Y over a box: numbers of cuts and of auxiliary variables
static std::pair<size_t, size_t>
relax(FFGraph& DAG, std::vector<FFVar> const& Y, std::vector<FFVar> const& X)
{
  std::vector<I> const IX = {I(-1., 1.), I(-0.5, 2.), I(0., 1.), I(-1., 0.5)};
  mc::PolImg<I> PolEnv;
  std::vector<mc::PolVar<I>> PolX, PolY(Y.size());
  for (size_t i = 0; i < X.size(); ++i) PolX.push_back(mc::PolVar<I>(&PolEnv, X[i], IX[i]));
  DAG.eval(Y, PolY, X, PolX);
  PolEnv.generate_cuts(PolY);
  return {PolEnv.Cuts().size(), PolEnv.Aux().size()};
}

// the tree u -> (u0 u1, exp(u0)), its values and its Jacobian by rows
static TREE
tree(FFGraph& G, std::vector<FFVar> const& U)
{
  return TREE(&G, U, {U[0] * U[1], exp(U[0])});
}
static std::vector<double>
values(std::vector<double> const& u)
{
  return {u[0] * u[1], std::exp(u[0])};
}
static std::vector<double>
jacobian(std::vector<double> const& u)
{
  return {u[1], u[0], std::exp(u[0]), 0.};
}

static std::vector<double> const x0 = {3., 5.};

int
main()
{
  using namespace test;
  std::printf("Expression trees as external operations: identification and home in a DAG\n");

  run("DAGEXT::uid: kept by copies, drawn again when the tree is set", []() -> Result {
    FFGraph G;
    std::vector<FFVar> U = G.add_vars(2);
    TREE A, B, D = tree(G, U);
    if (!A.uid() || A.uid() == B.uid() || A.uid() == D.uid()) return {false, "new objects"};
    TREE C(D);
    if (C.uid() != D.uid()) return {false, "copy"};
    size_t const uidD = D.uid();
    D.options.CPMAX = 5;
    if (D.uid() != uidD) return {false, "options"};
    D.set(&G, U, {U[0] + U[1]});
    if (D.uid() == uidD || C.uid() != uidD) return {false, "set"};
    TREE E(&G, U, {U[0] * U[1], exp(U[0])});  // the same expressions: still another stamp
    return {E.uid() != uidD && E.uid() != D.uid(), "same expressions"};
  });

  run("trees of the same shape, and a tree set again, are distinct operations", []() -> Result {
    FFGraph DAG, G1, G2;
    std::vector<FFVar> X = DAG.add_vars(2), U = G1.add_vars(2), V = G2.add_vars(2);
    TREE D1(&G1, U, {U[0] * U[1]}), D2(&G2, V, {V[0] + V[1]});
    OPTREE Op;
    std::vector<FFVar> const Y1 = Op(X, &D1), Y2 = Op(X, &D2);
    if (same(Y1, Y2)) return {false, "taken for the same operation"};
    if (evaluate(DAG, join(Y1, Y2), X, x0) != std::vector<double>{15., 8.}) return {false, "values"};
    // their Jacobians too: 2 entries each
    if (!near(evaluate(DAG, DAG.FAD(join(Y1, Y2), X), X, x0), {5., 3., 1., 1.})) return {false, "FAD"};
    if (!near(evaluate(DAG, DAG.BAD(join(Y1, Y2), X), X, x0), {5., 3., 1., 1.})) return {false, "BAD"};
    // the same object set to another tree
    D1.set(&G1, U, {U[0] - U[1]});
    std::vector<FFVar> const Y3 = Op(X, &D1);
    if (same(Y1, Y3) || evaluate(DAG, join(Y1, Y3), X, x0) != std::vector<double>{15., -2.}) return {false, "tree set again"};
    return {near(evaluate(DAG, DAG.FAD(join(Y1, Y3), X), X, x0), {5., 3., 1., -1.}), "FAD of the tree set again"};
  });

  run("an unchanged tree, a copy of it and its outputs one by one are the same operation", []() -> Result {
    FFGraph DAG, G;
    std::vector<FFVar> X = DAG.add_vars(2), U = G.add_vars(2);
    TREE D = tree(G, U);
    OPTREE Op, Op2;
    std::vector<FFVar> const Y = Op(X, &D);
    if (!same(Y, Op(X, &D)) || !same(Y, Op2(X, &D))) return {false, "inserted again"};
    if (!same(Y, {Op2(0, X, &D), Op2(1, X, &D)})) return {false, "outputs one by one"};
    TREE C(D);
    if (!same(Y, Op(X, &C))) return {false, "copy of the tree"};
    D.options.CPMAX = 7;  // not part of the identification
    return {same(Y, Op(X, &D)) && evaluate(DAG, Y, X, x0) == values(x0), "options changed, or values"};
  });

  run("one copy of the tree per DAG, and one Jacobian, for all the operations on it", []() -> Result {
    FFGraph DAG, G;
    std::vector<FFVar> X = DAG.add_vars(2), Z = DAG.add_vars(2), U = G.add_vars(2);
    TREE D = tree(G, U);
    OPTREE Op, Op2;
    std::vector<FFVar> const Y = Op(X, &D), W = Op2(Z, &D);
    TREE* const obj = Op.pDAGEXT();
    if (!obj || obj == &D || Op2.pDAGEXT() != obj || opof(Y[0])->object() != obj || opof(W[0])->object() != obj)
      return {false, "tree of the operations"};
    // Jacobian operations on either operands, by either mode: one further tree
    std::vector<FFVar> const dY = DAG.FAD(Y, X), dW = DAG.FAD(W, Z), bY = DAG.BAD(Y, X);
    std::set<TREE*> const all = trees(DAG, join(join(join(Y, W), join(dY, dW)), bY));
    if (all.size() != 2 || !all.count(obj)) return {false, std::to_string(all.size()) + " trees for the value and the Jacobian"};
    if (opof(dY[0])->object() == obj || opof(dY[0])->object() != opof(dW[0])->object()) return {false, "tree of the Jacobian operations"};
    // differentiated again: the same operations
    if (!same(dY, DAG.FAD(Y, X)) || !same(bY, DAG.BAD(Y, X))) return {false, "differentiated again"};
    std::vector<double> const z0 = {-1.2, 0.4};
    std::vector<FFVar> const XZ = join(X, Z);
    std::vector<double> const xz0 = {x0[0], x0[1], z0[0], z0[1]};
    if (!near(evaluate(DAG, dY, XZ, xz0), jacobian(x0)) || !near(evaluate(DAG, dW, XZ, xz0), jacobian(z0)) ||
        !near(evaluate(DAG, bY, XZ, xz0), jacobian(x0)))
      return {false, "derivative values"};
    // second derivatives: the Jacobian of the Jacobian tree, in a home of its own
    std::vector<FFVar> const ddY = DAG.FAD(dY, X);
    if (!near(evaluate(DAG, ddY, X, x0), {0., 1., 1., 0., std::exp(x0[0]), 0., 0., 0.})) return {false, "second derivatives"};
    if (trees(DAG, ddY).size() != 1 || trees(DAG, join(join(Y, dY), ddY)).size() != 3) return {false, "trees of the second derivatives"};
    // dense Jacobian: the same values
    OPTREE OpD(false);
    FFGraph DAGD;
    std::vector<FFVar> XD = DAGD.add_vars(2);
    std::vector<FFVar> const YD = OpD(XD, &D);
    return {near(evaluate(DAGD, DAGD.FAD(YD, XD), XD, x0), jacobian(x0)) && near(evaluate(DAGD, DAGD.BAD(YD, XD), XD, x0), jacobian(x0)),
            "dense Jacobian"};
  });

  run("options of the handle and of the tree of the home: changed after the insertion", []() -> Result {
    FFGraph DAG, G;
    DAG.options.MAXTHREAD = 0;
    std::vector<FFVar> X = DAG.add_vars(2), Z = DAG.add_vars(2), U = G.add_vars(2);
    TREE D(&G, U, {tanh(U[0] * U[1]) + sqr(U[0])});
    D.options.CPMAX = 3;
    OPTREE Op;
    Op.options.RELAX    = {Op.options.INT};
    Op.options.SCMODORD = 7;
    std::vector<FFVar> const Y = Op(X, &D), W = Op(Z, &D), dY = DAG.FAD(Y, X);
    // the options of the handle are those of the operations, the Jacobian ones included
    if (&opof(Y[0])->options != &Op.options || &opof(W[0])->options != &Op.options || &opof(dY[0])->options != &Op.options)
      return {false, "the operations do not share the options of the handle"};
    auto const rINT = relax(DAG, Y, X);
    Op.options.RELAX = {Op.options.AUX};
    auto const rAUX = relax(DAG, Y, X), rAUXW = relax(DAG, W, Z);
    if (rINT.second != 0 || rAUX.second == 0 || rAUXW != rAUX) return {false, "relaxation after the change"};
    // relaxed a second time, and back: nothing is kept from the previous options
    if (relax(DAG, Y, X) != rAUX) return {false, "relaxed again"};
    Op.options.RELAX = {Op.options.INT};
    if (relax(DAG, Y, X) != rINT) return {false, "options changed back"};
    // the tree of the Jacobian has the options of its parent, as of each differentiation
    if (opof(dY[0])->object()->options.CPMAX != 3) return {false, "options of the Jacobian tree"};
    D.options.CPMAX = 4;  // the original: no effect
    DAG.FAD(Y, X);
    if (opof(dY[0])->object()->options.CPMAX != 3) return {false, "options of the original taken"};
    Op.pDAGEXT()->options.CPMAX = 5;
    if (!same(dY, DAG.FAD(Y, X)) || opof(dY[0])->object()->options.CPMAX != 5) return {false, "options of the tree of the home"};
    // another AD strategy: the Jacobian is differentiated again, in a new home with the same options
    TREE* const jac = opof(dY[0])->object();
    Op.pDAGEXT()->options.AUTODIFF = TREE::Options::B;
    std::vector<FFVar> const bY = DAG.FAD(Y, X);
    if (opof(bY[0])->object() == jac || &opof(bY[0])->options != &Op.options) return {false, "Jacobian after a change of AD strategy"};
    return {near(evaluate(DAG, bY, X, x0), evaluate(DAG, dY, X, x0)), "derivative values"};
  });

  run("SHALLOW and TRANSFER: the home holds the tree given", []() -> Result {
    FFGraph DAG, G;
    std::vector<FFVar> X = DAG.add_vars(2), U = G.add_vars(2);
    TREE D = tree(G, U);
    OPTREE Op;
    std::vector<FFVar> const Y = Op(X, &D, Op.SHALLOW), dY = DAG.FAD(Y, X);
    if (Op.pDAGEXT() != &D || opof(Y[0])->object() != &D) return {false, "tree of the operation"};
    if (!same(Y, Op(X, &D, Op.SHALLOW)) || evaluate(DAG, Y, X, x0) != values(x0)) return {false, "inserted again, or values"};
    // the tree set again: the operation follows, and its next differentiation is that of the new tree
    D.set(&G, U, {U[0] - U[1], sqr(U[1])});
    if (!same(Y, Op(X, &D, Op.SHALLOW)) || evaluate(DAG, Y, X, x0) != std::vector<double>{-2., 25.}) return {false, "tree set again"};
    if (!near(evaluate(DAG, DAG.FAD(Y, X), X, x0), {1., -1., 0., 10.})) return {false, "FAD of the tree set again"};
    if (!near(evaluate(DAG, DAG.BAD(Y, X), X, x0), {1., -1., 0., 10.})) return {false, "BAD of the tree set again"};
    // a copy in the same DAG has another home
    std::vector<FFVar> const YC = Op(X, &D);
    if (same(Y, YC) || Op.pDAGEXT() == &D) return {false, "copy in the same DAG"};
    // a tree handed over to the DAG: held as given, deleted with its home
    TREE* pD = new TREE(&G, U, {U[0] * U[1], exp(U[0])});
    OPTREE OpT;
    std::vector<FFVar> const YT = OpT(X, pD, OpT.TRANSFER);
    if (OpT.pDAGEXT() != pD || !same(YT, OpT(X, pD, OpT.SHALLOW))) return {false, "tree handed over"};
    return {evaluate(DAG, YT, X, x0) == values(x0) && near(evaluate(DAG, DAG.FAD(YT, X), X, x0), jacobian(x0)), "values of the tree handed over"};
  });

  run("a copy of the DAG has homes of its own", []() -> Result {
    FFGraph G;
    std::vector<FFVar> U = G.add_vars(2);
    TREE D = tree(G, U);
    std::vector<FFVar> X2(2), Y2(2), dY2(4);
    FFGraph DAG2;
    OPTREE Op;
    Op.options.SCMODORD = 7;
    {
      FFGraph DAG;
      std::vector<FFVar> X = DAG.add_vars(2);
      std::vector<FFVar> const Y = Op(X, &D), dY = DAG.FAD(Y, X);
      DAG2.insert(&DAG, 2, X.data(), X2.data());
      DAG2.insert(&DAG, 2, Y.data(), Y2.data());
      DAG2.insert(&DAG, 4, dY.data(), dY2.data());
      // the tree and its Jacobian tree, copied; one set of options for both, copied too
      std::set<TREE*> const all = trees(DAG2, join(Y2, dY2));
      if (all.size() != 2 || all.count(Op.pDAGEXT()) || all.count(opof(dY[0])->object())) return {false, "trees of the copy"};
      if (&opof(Y2[0])->options == &Op.options || &opof(dY2[0])->options != &opof(Y2[0])->options) return {false, "options of the copy"};
      Op.options.SCMODORD = 9;
      if (opof(Y2[0])->options.SCMODORD != 7 || opof(dY2[0])->options.SCMODORD != 7) return {false, "change in the first DAG"};
    }
    // the first DAG is gone
    if (evaluate(DAG2, Y2, X2, x0) != values(x0) || !near(evaluate(DAG2, dY2, X2, x0), jacobian(x0))) return {false, "values in the copy"};
    // the handle finds the home of the tree in the copy; differentiation there gives the same values
    if (!same(Op(X2, &D), Y2) || Op.pDAGEXT() != opof(Y2[0])->object()) return {false, "tree not found in the copy"};
    return {near(evaluate(DAG2, DAG2.BAD(Y2, X2), X2, x0), jacobian(x0)), "differentiation in the copy"};
  });

  run("vector evaluation on several threads: values and derivatives, through the DAG copies of the threads", []() -> Result {
    FFGraph DAG, G;
    DAG.options.MAXTHREAD = 4;
    std::vector<FFVar> X = DAG.add_vars(2), U = G.add_vars(2);
    TREE D = tree(G, U);
    OPTREE Op;
    std::vector<FFVar> const Y = Op(X, &D);
    std::vector<FFVar> const F = join(join(Y, DAG.FAD(Y, X)), DAG.BAD(Y, X));
    std::vector<std::vector<double>> pts, vals;
    for (size_t k = 0; k < 64; ++k) pts.push_back({std::cos(1. + k), std::sin(2. * k)});
    for (int round = 0; round < 2; ++round)  // the second round reuses the copies of the DAG
    {
      DAG.veval(F, vals, X, pts);
      if (vals.size() != pts.size()) return {false, "number of results"};
      for (size_t k = 0; k < pts.size(); ++k)
      {
        std::vector<double> const y = values(pts[k]), J = jacobian(pts[k]);
        if (!near(std::vector<double>(vals[k].begin(), vals[k].begin() + 2), y)) return {false, "values at point " + std::to_string(k)};
        if (!near(std::vector<double>(vals[k].begin() + 2, vals[k].begin() + 6), J) ||
            !near(std::vector<double>(vals[k].begin() + 6, vals[k].end()), J))
          return {false, "derivatives at point " + std::to_string(k)};
      }
    }
    return {true, ""};
  });

  std::printf("\n=== Results: %d/%d passed ===\n", test::passed, test::total);
  return test::failed ? 1 : 0;
}
