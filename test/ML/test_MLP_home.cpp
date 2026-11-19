// =============================================================================
// test_MLP_home.cpp
// Unit tests for the home of a network in a DAG: the network -- one copy per
// DAG under policy COPY, the network given under SHALLOW and TRANSFER -- and
// the options of the operations on it, shared by every operation on that
// network in that DAG: on whichever operands, and the Jacobian, tangent and
// adjoint operations that differentiation makes. The options may be changed
// after the insertion, through the handle or through the network of the home,
// and the operations follow.
//
// Build (assuming MC++ headers are on the include path):
//   g++ -std=c++17 -O2 -I<path-to-mc++> -o test_MLP_home test_MLP_home.cpp
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

typedef mc::FFEXTERN<I, mc::MLP<I>> EXT;

// the operation that defines a variable, as an operation on a network
static EXT*
opof(FFVar const& V)
{
  return static_cast<EXT*>(V.opdef().first);
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

// networks held by the external operations that F depends on
static std::set<mc::MLP<I>*>
networks(FFGraph& DAG, std::vector<FFVar> const& F)
{
  std::set<mc::MLP<I>*> nets;
  auto sg = DAG.subgraph(F);
  for (auto const& op : sg.l_op)
    if (op->type >= mc::FFOp::EXTERN) nets.insert(&static_cast<EXT const*>(op)->object());
  return nets;
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

static std::vector<double> const x0 = {0.3, -0.7};

int
main()
{
  using namespace test;
  std::printf("Home of a multilayer perceptron in a DAG\n");

  run("COPY: one copy of the network per DAG, for all its operations and their derivatives", []() -> Result {
    FFGraph DAG;
    std::vector<FFVar> X = DAG.add_vars(2), Z = DAG.add_vars(2);
    mc::MLP<I> ANN;
    network(ANN, 2, 1.);
    ANN.options.TANGENT = true;  // tangent operations too
    mc::FFMLP<I> Op, Op2;
    std::vector<FFVar> const Y = Op(X, &ANN), W = Op(Z, &ANN);
    mc::MLP<I>* const net = &Op.MLP();
    if (!net || net == &ANN) return {false, "the handle does not point to a copy"};
    if (&opof(Y[0])->object() != net || &opof(W[0])->object() != net) return {false, "operations on other operands"};
    // a second handle, a copy of the network, the outputs one by one
    mc::MLP<I> COPY(ANN);
    if (&opof(Op2(1, X, &ANN))->object() != net || &Op2.MLP() != net) return {false, "second handle"};
    if (&opof(Op2(X, &COPY)[0])->object() != net) return {false, "copy of the network"};
    // Jacobian, tangent and adjoint operations
    std::vector<FFVar> const XZ = join(X, Z), F = {Y[0] * W[1], Y[1] + W[0]};
    std::vector<FFVar> const dF = DAG.BAD(F, XZ), tF = DAG.FAD(F, XZ);
    std::string const ops = operations(DAG, join(dF, tF));
    if (ops.find("GradMLP") == std::string::npos || ops.find("TanMLP") == std::string::npos)
      return {false, "test: derivative operations " + ops};
    std::set<mc::MLP<I>*> const nets = networks(DAG, join(join(F, dF), tF));
    if (nets.size() != 1 || *nets.begin() != net) return {false, "derivative operations: " + std::to_string(nets.size()) + " networks"};
    // derivatives are those of the network
    std::vector<double> const z0 = {-1.2, 0.4}, xz0 = {x0[0], x0[1], z0[0], z0[1]};
    std::vector<double> const y = direct(ANN, x0), w = direct(ANN, z0), Jy = jacobian(ANN, x0), Jw = jacobian(ANN, z0);
    std::vector<double> const ref = {Jy[0] * w[1], Jy[1] * w[1], y[0] * Jw[2], y[0] * Jw[3], Jy[2], Jy[3], Jw[0], Jw[1]};
    if (!near(evaluate(DAG, dF, XZ, xz0), ref) || !near(evaluate(DAG, tF, XZ, xz0), ref)) return {false, "derivative values"};
    // a changed network has another home, and the first one is left alone
    network(ANN, 2, 2.);
    std::vector<FFVar> const Y2 = Op(X, &ANN);
    if (&opof(Y2[0])->object() == net || &Op.MLP() != &opof(Y2[0])->object()) return {false, "changed network"};
    if (&opof(Y[0])->object() != net) return {false, "first network after the change"};
    // another DAG has its own copy
    FFGraph DAGB;
    std::vector<FFVar> XB = DAGB.add_vars(2);
    std::vector<FFVar> const YB = Op(XB, &COPY);
    if (&opof(YB[0])->object() == net || &opof(YB[0])->object() == &COPY) return {false, "another DAG"};
    return {evaluate(DAGB, YB, XB, x0) == y, "values in another DAG"};
  });

  run("options of the handle: those of its operations, before and after the insertion", []() -> Result {
    FFGraph DAG;
    DAG.options.MAXTHREAD = 0;
    std::vector<FFVar> X = DAG.add_vars(2), Z = DAG.add_vars(2);
    mc::MLP<I> ANN;
    network(ANN, 1, 1.);
    mc::FFMLP<I> Op;
    Op.options.RELAX    = {Op.options.INT};
    Op.options.SCMODORD = 7;
    std::vector<FFVar> const Y = Op(X, &ANN), W = Op(Z, &ANN), dY = DAG.FAD(Y, X);
    if (operations(DAG, dY).find("GradMLP") == std::string::npos) return {false, "test: no Jacobian operation"};
    if (&opof(Y[0])->options != &Op.options || &opof(W[0])->options != &Op.options || &opof(dY[0])->options != &Op.options)
      return {false, "the operations do not share the options of the handle"};
    if (opof(dY[0])->options.SCMODORD != 7) return {false, "Jacobian operation: options set before the insertion"};
    // the relaxation follows a change made after the insertion: bounds only, then auxiliary variables
    auto const rINT = relax(DAG, Y, X), rINTW = relax(DAG, W, Z);
    Op.options.RELAX = {Op.options.AUX};
    auto const rAUX = relax(DAG, Y, X), rAUXW = relax(DAG, W, Z);
    if (rINT.second != 0 || rAUX.second == 0 || rINTW != rINT || rAUXW != rAUX)
      return {false, "relaxation after the change: " + std::to_string(rINT.first) + "/" + std::to_string(rINT.second) + " then " +
                         std::to_string(rAUX.first) + "/" + std::to_string(rAUX.second)};
    // the same relaxation as with the option set before the insertion
    FFGraph DAGB;
    DAGB.options.MAXTHREAD = 0;
    std::vector<FFVar> XB = DAGB.add_vars(2);
    mc::FFMLP<I> OpB;
    OpB.options.RELAX = {OpB.options.AUX};
    if (relax(DAGB, OpB(XB, &ANN), XB) != rAUX) return {false, "differs from the option set before the insertion"};
    // and that handle was left alone
    if (OpB.options.SCMODORD == 7) return {false, "options shared between handles"};
    // operations made after the change: further operands, a new differentiation
    std::vector<FFVar> const dW = DAG.FAD(W, Z);
    Op.options.SCMODORD = 9;
    return {opof(dW[0])->options.SCMODORD == 9 && opof(dY[0])->options.SCMODORD == 9, "operations made later"};
  });

  run("options of the network of the home: changed after the insertion, for all its operations", []() -> Result {
    FFGraph DAG;
    std::vector<FFVar> X = DAG.add_vars(2), Z = DAG.add_vars(2);
    mc::MLP<I> ANN;
    network(ANN, 2, 1.);
    mc::FFMLP<I> Op;
    std::vector<FFVar> const Y = Op(X, &ANN), W = Op(Z, &ANN);
    // the Jacobian operation without the option
    std::vector<FFVar> const dY = DAG.FAD(Y, X);
    if (operations(DAG, dY).find("TanMLP") != std::string::npos) return {false, "test: tangent operation by default"};
    // the original is a different object: no effect
    ANN.options.TANGENT = true;
    if (operations(DAG, DAG.FAD(W, Z)).find("TanMLP") != std::string::npos) return {false, "option of the original taken"};
    // the network of the home: the next differentiation of any of its operations follows
    ANN.options.TANGENT        = false;
    Op.MLP().options.TANGENT = true;
    std::vector<FFVar> const tY = DAG.FAD(Y, X), tW = DAG.FAD(W, Z);
    if (operations(DAG, tY).find("TanMLP") == std::string::npos || operations(DAG, tW).find("TanMLP") == std::string::npos)
      return {false, "option of the home not taken: " + operations(DAG, tY) + "| " + operations(DAG, tW)};
    // the tangent operations see the options too, being on the same network
    if (!opof(tY[0])->object().options.TANGENT) return {false, "tangent operation on another network"};
    // derivatives already made are left as they are, and all agree
    if (operations(DAG, dY).find("TanMLP") != std::string::npos) return {false, "existing derivatives changed"};
    return {near(evaluate(DAG, tY, X, x0), evaluate(DAG, dY, X, x0)) && near(evaluate(DAG, tY, X, x0), jacobian(ANN, x0)), "values"};
  });

  run("a network inserted first by another handle keeps the options of that handle", []() -> Result {
    FFGraph DAG;
    DAG.options.MAXTHREAD = 0;
    std::vector<FFVar> X = DAG.add_vars(2), Z = DAG.add_vars(2);
    mc::MLP<I> ANN;
    network(ANN, 1, 1.);
    mc::FFMLP<I> Op, Op2;
    Op.options.RELAX  = {Op.options.INT};
    Op2.options.RELAX = {Op2.options.AUX};
    if (Op.attached() || &Op.home_options() != &Op.options) return {false, "handle before any insertion"};
    std::vector<FFVar> const Y = Op(X, &ANN), W = Op2(Z, &ANN);
    if (!Op.attached() || !Op2.attached()) return {false, "handles not attached"};
    if (&Op.home_options() != &Op.options) return {false, "first handle"};
    if (&Op2.home_options() != &Op.options || &opof(W[0])->options != &Op.options) return {false, "second handle"};
    if (relax(DAG, W, Z).second != 0) return {false, "relaxation with the options of the second handle"};
    // through the second handle all the same
    Op2.home_options().RELAX = {Op2.options.AUX};
    if (relax(DAG, Y, X).second == 0 || relax(DAG, W, Z).second == 0) return {false, "change through home_options"};
    // another network through the second handle: a home of its own, with the options of that handle
    mc::MLP<I> ANN2;
    network(ANN2, 1, 2.);
    std::vector<FFVar> const V = Op2(X, &ANN2);
    return {&opof(V[0])->options == &Op2.options && &Op2.home_options() == &Op2.options, "another network"};
  });

  run("SHALLOW and TRANSFER: the home holds the network given, also for the derivatives", []() -> Result {
    FFGraph DAG;
    std::vector<FFVar> X = DAG.add_vars(2);
    mc::MLP<I> ANN;
    network(ANN, 2, 1.);
    mc::FFMLP<I> Op;
    std::vector<FFVar> const Y = Op(X, &ANN, Op.SHALLOW), dY = DAG.FAD(Y, X), bY = DAG.BAD(Y, X);
    std::set<mc::MLP<I>*> const nets = networks(DAG, join(join(Y, dY), bY));
    if (&Op.MLP() != &ANN || nets.size() != 1 || *nets.begin() != &ANN) return {false, "network of the operations"};
    // a change of the network is seen by the value and by the derivatives already made
    network(ANN, 2, 2.);
    if (evaluate(DAG, Y, X, x0) != direct(ANN, x0)) return {false, "values after the change"};
    if (!near(evaluate(DAG, dY, X, x0), jacobian(ANN, x0)) || !near(evaluate(DAG, bY, X, x0), jacobian(ANN, x0)))
      return {false, "derivatives after the change"};
    // and so are its options: tangent operations at the next differentiation
    ANN.options.TANGENT = true;
    if (operations(DAG, DAG.FAD(Y, X)).find("TanMLP") == std::string::npos) return {false, "option of the network"};
    // a copy in the same DAG has another home, with its own network
    std::vector<FFVar> const YC = Op(X, &ANN);
    if (&Op.MLP() == &ANN || &opof(YC[0])->object() == &ANN) return {false, "copy in the same DAG"};
    // a network handed over to the DAG: held as given, and deleted with its home
    mc::MLP<I>* pANN = new mc::MLP<I>;
    network(*pANN, 2, 3.);
    std::vector<double> const y3 = direct(*pANN, x0), J3 = jacobian(*pANN, x0);
    mc::FFMLP<I> OpT;
    std::vector<FFVar> const YT = OpT(X, pANN, OpT.TRANSFER), dYT = DAG.FAD(YT, X);
    if (&OpT.MLP() != pANN || &opof(dYT[0])->object() != pANN) return {false, "network handed over"};
    return {evaluate(DAG, YT, X, x0) == y3 && near(evaluate(DAG, dYT, X, x0), J3), "values of the network handed over"};
  });

  run("a copy of the DAG has homes of its own: network and options", []() -> Result {
    mc::MLP<I> ANN;
    network(ANN, 2, 1.);
    std::vector<double> const y = direct(ANN, x0), J = jacobian(ANN, x0);
    std::vector<FFVar> X2(2), Y2(2), dY2(4);
    FFGraph DAG2;
    mc::FFMLP<I> Op;
    Op.options.SCMODORD = 7;
    mc::MLP<I>* net = nullptr;
    {
      FFGraph DAG;
      std::vector<FFVar> X = DAG.add_vars(2);
      std::vector<FFVar> const Y = Op(X, &ANN), dY = DAG.FAD(Y, X);
      net = &Op.MLP();
      DAG2.insert(&DAG, 2, X.data(), X2.data());
      DAG2.insert(&DAG, 2, Y.data(), Y2.data());
      DAG2.insert(&DAG, 4, dY.data(), dY2.data());
      // one home in the copy, for the value and for the Jacobian: another network, other options with the same values
      std::set<mc::MLP<I>*> const nets = networks(DAG2, join(Y2, dY2));
      if (nets.size() != 1 || *nets.begin() == net) return {false, "network of the copy"};
      if (&opof(Y2[0])->options == &Op.options || &opof(dY2[0])->options != &opof(Y2[0])->options) return {false, "options of the copy"};
      if (opof(Y2[0])->options.SCMODORD != 7) return {false, "values of the options of the copy"};
      // each DAG follows its own
      Op.options.SCMODORD = 9;
      if (opof(Y[0])->options.SCMODORD != 9 || opof(Y2[0])->options.SCMODORD != 7) return {false, "change in the first DAG"};
      opof(Y2[0])->options.SCMODORD      = 5;
      opof(Y2[0])->object().options.TANGENT = true;
      if (opof(Y[0])->options.SCMODORD != 9 || net->options.TANGENT) return {false, "change in the copy"};
    }
    // the first DAG is gone
    if (evaluate(DAG2, Y2, X2, x0) != y || !near(evaluate(DAG2, dY2, X2, x0), J)) return {false, "values in the copy"};
    // the handle finds the home of the network in the copy, and is then attached to it
    std::vector<FFVar> const Y2again = Op(X2, &ANN);
    if (Y2again[0].id() != Y2[0].id() || &Op.MLP() != &opof(Y2[0])->object()) return {false, "network not found in the copy"};
    return {Op.home_options().SCMODORD == 5 && operations(DAG2, DAG2.FAD(Y2again, X2)).find("TanMLP") != std::string::npos,
            "options of the home in the copy"};
  });

  run("homes from several threads, each with its own DAG", []() -> Result {
    mc::MLP<I> ANN;
    network(ANN, 2, 1.);
    std::vector<double> const y = direct(ANN, x0), J = jacobian(ANN, x0);
    std::vector<int> ok(4, 0);
    std::vector<std::thread> threads;
    for (size_t t = 0; t < 4; ++t)
      threads.emplace_back([&ANN, &y, &J, &ok, t]() {
        bool good = true;
        for (size_t k = 0; k < 20 && good; ++k)
        {
          FFGraph DAG;
          std::vector<FFVar> X = DAG.add_vars(2);
          mc::FFMLP<I> Op;
          mc::MLP<I> copy(ANN);
          std::vector<FFVar> const Y = Op(X, &copy), dY = DAG.FAD(Y, X);
          good = (evaluate(DAG, Y, X, x0) == y) && near(evaluate(DAG, dY, X, x0), J) && &Op.MLP() != &copy;
        }
        ok[t] = good;
      });
    for (auto& th : threads) th.join();
    return {ok[0] && ok[1] && ok[2] && ok[3], "values"};
  });

  run("MLP(): the network of the DAG by reference, from every class of MLP operation", []() -> Result {
    mc::FFMLP<I> Op;
    mc::FFGradMLP<I> OpG;
    mc::FFTanMLP<I> OpT;
    mc::FFAdjMLP<I> OpA;
    // before any insertion: an exception, and a null pointer from the deprecated accessor
    int thrown = 0;
    for (std::function<void()> f : std::vector<std::function<void()>>{[&]() { Op.MLP(); }, [&]() { OpG.MLP(); },
                                                                      [&]() { OpT.MLP(); }, [&]() { OpA.MLP(); }})
    {
      try { f(); }
      catch (std::runtime_error const&) { ++thrown; }
    }
    if (thrown != 4) return {false, "no exception before an insertion: " + std::to_string(4 - thrown) + " classes"};
#if defined(__GNUC__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#endif
    if (Op.pMLP()) return {false, "pMLP() before an insertion"};
    FFGraph DAG;
    std::vector<FFVar> X = DAG.add_vars(2);
    mc::MLP<I> ANN;
    network(ANN, 2, 1.);
    ANN.options.TANGENT = ANN.options.ADJOINT = true;
    std::vector<FFVar> const Y = Op(X, &ANN);
    if (Op.pMLP() != &Op.MLP()) return {false, "pMLP() and MLP() differ"};
#if defined(__GNUC__)
#pragma GCC diagnostic pop
#endif
    // the derivative operations of every class hold the same network
    std::vector<FFVar> const F = {Y[0] * Y[1]};
    std::vector<FFVar> const dJ = DAG.FAD(Y, X), dT = DAG.FAD(F, X), dA = DAG.BAD(F, X);
    std::string const ops = operations(DAG, join(join(dJ, dT), dA));
    if (ops.find("TanMLP") == std::string::npos || ops.find("AdjMLP") == std::string::npos) return {false, "test: operations " + ops};
    for (auto const& V : join(join(dJ, dT), dA))
      if (V.opdef().first && V.opdef().first->type >= mc::FFOp::EXTERN && &opof(V)->object() != &Op.MLP())
        return {false, "a derivative operation holds another network"};
    // const access
    mc::FFMLP<I> const& cOp = Op;
    return {&cOp.MLP() == &Op.MLP() && cOp.MLP().nout() == 2, "const access"};
  });

  run("a network of the DAG given other weights: its operations follow, the original gets another one", []() -> Result {
    FFGraph DAG;
    std::vector<FFVar> X = DAG.add_vars(2);
    mc::MLP<I> ANN;
    network(ANN, 2, 1.);
    std::vector<double> const y1 = direct(ANN, x0);
    mc::FFMLP<I> Op;
    std::vector<FFVar> const Y = Op(X, &ANN);
    mc::MLP<I>& net = Op.MLP();
    network(net, 2, 2.);  // the copy of the DAG, given other weights
    std::vector<double> const y2 = direct(net, x0);
    if (y1 == y2) return {false, "test networks"};
    if (evaluate(DAG, Y, X, x0) != y2) return {false, "the operation does not follow its network"};
    if (!near(evaluate(DAG, DAG.FAD(Y, X), X, x0), jacobian(net, x0))) return {false, "derivatives made after the change"};
    // the original, unchanged, inserted again: another operation, with a copy of the original
    mc::FFMLP<I> Op2;
    std::vector<FFVar> const Y1 = Op2(X, &ANN);
    if (same(Y, Y1) || &Op2.MLP() == &net || evaluate(DAG, Y1, X, x0) != y1) return {false, "the original inserted again"};
    // through the same handle, and again: found each time
    if (!same(Op(X, &ANN), Y1) || !same(Op2(X, &ANN), Y1)) return {false, "the original inserted a third time"};
    // a copy of the changed network is another operation too, with the same values
    mc::MLP<I> copy(net);
    std::vector<FFVar> const Y2 = Op(X, &copy);
    return {!same(Y2, Y) && !same(Y2, Y1) && evaluate(DAG, Y2, X, x0) == y2, "a copy of the changed network"};
  });

  run("vector evaluation on several threads: values and derivatives, through the DAG copies of the threads", []() -> Result {
    FFGraph DAG;
    DAG.options.MAXTHREAD = 4;
    std::vector<FFVar> X = DAG.add_vars(2);
    mc::MLP<I> ANN;
    network(ANN, 2, 1.);
    mc::MLP<I> ANNT(ANN);
    ANNT.options.TANGENT = true;
    mc::FFMLP<I> Op;
    std::vector<FFVar> const Y = Op(X, &ANN), YT = Op(X, &ANNT);
    std::vector<FFVar> const F = join(join(Y, DAG.FAD(Y, X)), DAG.FAD(YT, X));
    std::vector<std::vector<double>> pts, vals;
    for (size_t k = 0; k < 64; ++k) pts.push_back({std::cos(1. + k), std::sin(2. * k)});
    for (int round = 0; round < 2; ++round)  // the second round reuses the copies of the DAG
    {
      DAG.veval(F, vals, X, pts);
      if (vals.size() != pts.size()) return {false, "number of results"};
      for (size_t k = 0; k < pts.size(); ++k)
      {
        std::vector<double> const y = direct(ANN, pts[k]), J = jacobian(ANN, pts[k]);
        if (std::vector<double>(vals[k].begin(), vals[k].begin() + 2) != y) return {false, "values at point " + std::to_string(k)};
        if (!near(std::vector<double>(vals[k].begin() + 2, vals[k].begin() + 6), J) ||
            !near(std::vector<double>(vals[k].begin() + 6, vals[k].end()), J))
          return {false, "derivatives at point " + std::to_string(k)};
      }
    }
    return {true, ""};
  });

  run("RELU2ABS and SIG2EXP changed after a relaxation: the expression tree of the network follows", []() -> Result {
    // operations of each type in the expression tree of a network
    auto count = [](mc::MLP<I>& NN, int const type) {
      size_t n = 0;
      auto sg = NN.dag()->subgraph(NN.varout());
      for (auto const& op : sg.l_op) n += (op->type == type);
      return n;
    };
    mc::MLP<I> ANN;
    ANN.append_data({{0.1, 1., -0.5}, {-0.2, 0.3, 0.8}}, ANN.RELU, true);
    ANN.append_data({{0.0, 1., 1.}, {0.5, -1., 2.}}, ANN.SIGMOID);
    FFGraph DAG;
    std::vector<FFVar> X = DAG.add_vars(2);
    mc::FFMLP<I> Op;
    Op.options.RELAX = {Op.options.AUX};
    std::vector<FFVar> Y = Op(X, &ANN);
    auto const r0 = relax(DAG, Y, X);
    mc::MLP<I>& NN = Op.MLP();
    bool ok = count(NN, mc::FFOp::MAXF) == 2 && count(NN, mc::FFOp::FABS) == 0 &&
              count(NN, mc::FFOp::TANH) == 2 && count(NN, mc::FFOp::EXP) == 0;
    NN.options.RELU2ABS = true;
    NN.options.SIG2EXP  = true;
    auto const r1 = relax(DAG, Y, X);
    ok = ok && count(NN, mc::FFOp::MAXF) == 0 && count(NN, mc::FFOp::FABS) == 2 &&
         count(NN, mc::FFOp::TANH) == 0 && count(NN, mc::FFOp::EXP) == 2 && r1 != r0;
    // and back
    NN.options.RELU2ABS = false;
    NN.options.SIG2EXP  = false;
    ok = ok && relax(DAG, Y, X) == r0 && count(NN, mc::FFOp::MAXF) == 2;
    // the values are the same with either form
    std::vector<double> const x = {0.3, -0.7};
    std::vector<double> const y0 = evaluate(DAG, Y, X, x);
    NN.options.RELU2ABS = NN.options.SIG2EXP = true;
    ok = ok && near(evaluate(DAG, Y, X, x), y0);
    return {ok, "cuts " + std::to_string(r0.first) + " then " + std::to_string(r1.first)};
  });

  std::printf("\n=== Results: %d/%d passed ===\n", test::passed, test::total);
  return test::failed ? 1 : 0;
}
