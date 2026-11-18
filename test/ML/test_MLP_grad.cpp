// =============================================================================
// test_MLP_grad.cpp
// Unit tests for the Jacobian of mc::MLP in double arithmetic, mc::MLP::grad:
// a reverse pass through the layers where the outputs are no more than the
// inputs (or with Options::AUTODIFF = B), a forward pass otherwise. Against
// central differences, the tangent and adjoint passes, and the Jacobian
// operation mc::FFGradMLP in a DAG.
//
// Build (assuming MC++ headers are on the include path):
//   g++ -std=c++17 -O2 -I<path-to-mc++> -o test_MLP_grad test_MLP_grad.cpp
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
// Networks with given sizes {inputs, hidden layers..., outputs} and activation
// function of the hidden layers; the outputs are linear, or of that function too
// -----------------------------------------------------------------------------
static void
network(mc::MLP<I>& ANN, std::vector<size_t> const& sizes, int const activ, bool const outactiv = false)
{
  for (size_t l = 1; l < sizes.size(); ++l)
  {
    std::vector<std::vector<double>> layer(sizes[l], std::vector<double>(sizes[l - 1] + 1));
    for (size_t i = 0; i < sizes[l]; ++i)
      for (size_t j = 0; j <= sizes[l - 1]; ++j) layer[i][j] = 0.8 * std::cos(1. + 3. * l + 5. * i + 7. * j);
    ANN.append_data(layer, (l + 1 < sizes.size() || outactiv) ? activ : (int)ANN.LINEAR, l == 1);
  }
}

static std::vector<double>
point(size_t const n, int const k)
{
  std::vector<double> x(n);
  for (size_t i = 0; i < n; ++i) x[i] = 0.9 * std::cos(0.7 + 2.3 * k + 1.7 * i);
  return x;
}

// Jacobian by central differences
static std::vector<double>
differences(mc::MLP<I>& ANN, std::vector<double> const& x)
{
  size_t const nin = ANN.nin(), nout = ANN.nout();
  std::vector<double> J(nout * nin), yp(nout), ym(nout), xp = x, xm = x;
  double const h = 1e-6;
  for (size_t i = 0; i < nin; ++i)
  {
    xp[i] = x[i] + h;
    xm[i] = x[i] - h;
    ANN.eval(xp.data(), yp.data());
    ANN.eval(xm.data(), ym.data());
    for (size_t k = 0; k < nout; ++k) J[k * nin + i] = (yp[k] - ym[k]) / (2 * h);
    xp[i] = xm[i] = x[i];
  }
  return J;
}

// Jacobian by tangent passes (columns) and by adjoint passes (rows)
static std::vector<double>
passes(mc::MLP<I>& ANN, std::vector<double> const& x, bool const tangent)
{
  size_t const nin = ANN.nin(), nout = ANN.nout();
  std::vector<double> J(nout * nin), d(nin, 0.), t(nout), w(nout, 0.);
  if (tangent)
    for (size_t i = 0; i < nin; ++i)
    {
      d[i] = 1.;
      ANN.tang(x.data(), d.data(), t.data());
      for (size_t k = 0; k < nout; ++k) J[k * nin + i] = t[k];
      d[i] = 0.;
    }
  else
    for (size_t k = 0; k < nout; ++k)
    {
      w[k] = 1.;
      ANN.adj(x.data(), w.data(), J.data() + k * nin);
      w[k] = 0.;
    }
  return J;
}

// {inputs, hidden layers..., outputs}: more inputs than outputs, as many, fewer, a single one of each,
// and no hidden layer
static std::vector<std::vector<size_t>> const shapes = {{4, 6, 5, 2}, {3, 5, 3}, {2, 5, 4, 7}, {1, 4, 1}, {5, 6, 1},
                                                        {1, 6, 5}, {3, 2}, {2, 3}};

int
main()
{
  using test::near;
  using test::Result;
  using test::run;

  std::printf("Jacobian of multilayer perceptrons in double arithmetic\n");

  run("against central differences and the tangent and adjoint passes: all sizes and activation functions", []() -> Result {
    for (auto const& sizes : shapes)
      for (int const variant : {0, 1, 2, 3, 4, 5})
      {
        mc::MLP<I> ANN;
        int const activ = (variant == 0 ? (int)ANN.TANH : variant <= 2 ? (int)ANN.SIGMOID : variant <= 4 ? (int)ANN.RELU : (int)ANN.LINEAR);
        ANN.options.SIG2EXP  = (variant == 2);
        ANN.options.RELU2ABS = (variant == 4);
        network(ANN, sizes, activ, variant % 2);
        size_t const nin = ANN.nin(), nout = ANN.nout();
        std::string const id = std::to_string(nin) + " inputs, " + std::to_string(nout) + " outputs, variant " + std::to_string(variant);
        for (int k = 0; k < 3; ++k)
        {
          std::vector<double> const x = point(nin, k);
          std::vector<double> J(nout * nin), y(nout), yref(nout);
          ANN.eval(x.data(), yref.data());
          ANN.grad(x.data(), J.data(), y.data());
          // the values come out as eval computes them
          if (y != yref) return {false, id + ": values returned with the Jacobian"};
          if (!near(J, differences(ANN, x), 1e-6)) return {false, id + ": central differences"};
          if (!near(J, passes(ANN, x, true), 1e-12)) return {false, id + ": tangent passes"};
          if (!near(J, passes(ANN, x, false), 1e-12)) return {false, id + ": adjoint passes"};
          // without the values
          std::vector<double> J2(nout * nin);
          ANN.grad(x.data(), J2.data());
          if (J2 != J) return {false, id + ": without the values"};
        }
      }
    return {true, ""};
  });

  run("the reverse mode asked for with AUTODIFF = B agrees with the forward mode", []() -> Result {
    for (auto const& sizes : shapes)
    {
      mc::MLP<I> ANN;
      network(ANN, sizes, ANN.TANH);
      size_t const nin = ANN.nin(), nout = ANN.nout();
      for (int k = 0; k < 3; ++k)
      {
        std::vector<double> const x = point(nin, k);
        std::vector<double> JF(nout * nin), JB(nout * nin), yF(nout), yB(nout);
        ANN.options.AUTODIFF = ANN.options.F;
        ANN.grad(x.data(), JF.data(), yF.data());
        ANN.options.AUTODIFF = ANN.options.B;
        ANN.grad(x.data(), JB.data(), yB.data());
        if (yF != yB || !near(JB, JF, 1e-12)) return {false, std::to_string(nin) + " inputs, " + std::to_string(nout) + " outputs"};
        // with fewer inputs than outputs the two are different passes: not the same rounding errors
        if (nin < nout && sizes.size() > 2 && JB == JF) return {false, "the same pass for both modes"};
      }
    }
    return {true, ""};
  });

  run("neurons with a flat activation and weights that are zero", []() -> Result {
    for (auto const& sizes : std::vector<std::vector<size_t>>{{3, 4, 2}, {2, 4, 3}})
    {
      // ReLU neurons: the second never active, the third always, and the first input without weights
      mc::MLP<I> ANN;
      std::vector<std::vector<double>> layer1(4, std::vector<double>(sizes[0] + 1)), layer2(sizes[2], std::vector<double>(5));
      for (size_t i = 0; i < 4; ++i)
        for (size_t j = 0; j <= sizes[0]; ++j) layer1[i][j] = (j == 1 ? 0. : 0.8 * std::cos(1. + 5. * i + 7. * j));
      layer1[1][0] = -100.;
      layer1[2][0] = 100.;
      for (size_t i = 0; i < sizes[2]; ++i)
        for (size_t j = 0; j < 5; ++j) layer2[i][j] = 0.8 * std::cos(2. + 3. * i + 11. * j);
      ANN.append_data(layer1, ANN.RELU, true);
      ANN.append_data(layer2, ANN.LINEAR);
      size_t const nin = ANN.nin(), nout = ANN.nout();
      for (int k = 0; k < 3; ++k)
      {
        std::vector<double> const x = point(nin, k);
        std::vector<double> J(nout * nin);
        ANN.grad(x.data(), J.data());
        for (size_t i = 0; i < nout; ++i)
          if (J[i * nin] != 0.) return {false, "derivative in an input without weights"};
        if (!near(J, differences(ANN, x), 1e-6) || !near(J, passes(ANN, x, true), 1e-12)) return {false, "Jacobian"};
      }
    }
    return {true, ""};
  });

  run("the Jacobian operation in a DAG, and the Jacobian in interval arithmetic", []() -> Result {
    for (auto const& sizes : shapes)
    {
      mc::MLP<I> ANN;
      network(ANN, sizes, ANN.TANH);
      size_t const nin = ANN.nin(), nout = ANN.nout();
      FFGraph DAG;
      std::vector<FFVar> X = DAG.add_vars(nin);
      mc::FFMLP<I> Op;
      std::vector<FFVar> const Y = Op(X, &ANN), dYf = DAG.FAD(Y, X), dYb = DAG.BAD(Y, X);
      std::vector<double> const x = point(nin, 0);
      std::vector<double> J(nout * nin), df(nout * nin), db(nout * nin);
      ANN.grad(x.data(), J.data());
      DAG.eval(dYf, df, X, x);
      DAG.eval(dYb, db, X, x);
      if (df != J || db != J) return {false, "Jacobian operation"};
      // interval arithmetic still goes through forward AD: bounds that enclose the Jacobian
      std::vector<I> IX(nin), IJ(nout * nin);
      for (size_t i = 0; i < nin; ++i) IX[i] = I(x[i] - 0.05, x[i] + 0.05);
      ANN.grad(IX.data(), IJ.data());
      for (size_t k = 0; k < nout * nin; ++k)
        if (J[k] < IJ[k].l() - 1e-12 || J[k] > IJ[k].u() + 1e-12 || !(IJ[k].u() - IJ[k].l() < 10.)) return {false, "interval bounds"};
    }
    return {true, ""};
  });

  run("several threads on one MLP", []() -> Result {
    for (auto const& sizes : std::vector<std::vector<size_t>>{{4, 6, 5, 2}, {2, 5, 4, 7}})
    {
      mc::MLP<I> ANN;
      network(ANN, sizes, ANN.TANH);
      size_t const nin = ANN.nin(), nout = ANN.nout(), npts = 5;
      std::vector<std::vector<double>> pts(npts), ref(npts, std::vector<double>(nout * nin));
      for (size_t k = 0; k < npts; ++k)
      {
        pts[k] = point(nin, k);
        ANN.grad(pts[k].data(), ref[k].data());
      }
      std::vector<int> nbad(4, 0);
      std::vector<std::thread> threads;
      for (size_t ith = 0; ith < 4; ++ith)
        threads.emplace_back([&, ith]() {
          std::vector<double> J(nout * nin);
          for (size_t it = 0; it < 2000; ++it)
          {
            size_t const k = (it + ith) % npts;
            ANN.grad(pts[k].data(), J.data());
            if (J != ref[k]) ++nbad[ith];
          }
        });
      for (auto& th : threads) th.join();
      for (int const n : nbad)
        if (n) return {false, std::to_string(n) + " Jacobians differ in a thread"};
    }
    return {true, ""};
  });

  std::printf("\n=== Results: %d/%d passed ===\n", test::passed, test::total);
  return test::failed ? 1 : 0;
}
