// =============================================================================
// DAG_tangent.cpp
// Unit tests for the tangent hook of external operations, FFOp::has_tangent()
// and FFOp::tangent(): an operation that propagates tangents itself in the
// forward-mode symbolic differentiation of mc::FFGraph (FAD, SFAD, DFAD), in
// place of the partial derivatives returned by FFOp::deriv().
//
// The operation under test is y = ( x0 x1, sin(x0)+x2^2, exp(x1 x2) ), with
//  - FFToy     its value                                      [3 operands, 3 results]
//  - FFToyJac  its Jacobian, the result of deriv()            [3 operands, 9 results]
//  - FFToyTan  its directional derivative J(x) t, the result  [6 operands, 3 results]
//              of tangent(); it has a tangent hook of its own
//  - FFToyTan2 its second directional derivative D2y[t,s]     [9 operands, 3 results]
// The derivatives obtained through the hook are compared with those obtained
// through deriv() and with analytic values.
//
// Build (assuming MC++ headers are on the include path):
//   g++ -std=c++17 -O2 -I<path-to-mc++> -o DAG_tangent DAG_tangent.cpp
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

#include "ffunc.hpp"
#include "ffdep.hpp"

using mc::FFGraph;
using mc::FFOp;
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
// The function and its derivatives, in double arithmetic
// -----------------------------------------------------------------------------
namespace toy
{

inline void
val(double const* x, double* y)
{
  y[0] = x[0] * x[1];
  y[1] = std::sin(x[0]) + x[2] * x[2];
  y[2] = std::exp(x[1] * x[2]);
}

// J[k*3+i] = dy_k/dx_i
inline void
jac(double const* x, double* J)
{
  double const e = std::exp(x[1] * x[2]);
  J[0] = x[1];           J[1] = x[0];      J[2] = 0.;
  J[3] = std::cos(x[0]); J[4] = 0.;        J[5] = 2. * x[2];
  J[6] = 0.;             J[7] = x[2] * e;  J[8] = x[1] * e;
}

// J(x) t
inline void
tan(double const* x, double const* t, double* r)
{
  double const e = std::exp(x[1] * x[2]);
  r[0] = x[1] * t[0] + x[0] * t[1];
  r[1] = std::cos(x[0]) * t[0] + 2. * x[2] * t[2];
  r[2] = e * (x[2] * t[1] + x[1] * t[2]);
}

// D2y[t,s]
inline void
tan2(double const* x, double const* t, double const* s, double* r)
{
  double const e = std::exp(x[1] * x[2]);
  r[0] = t[0] * s[1] + t[1] * s[0];
  r[1] = -std::sin(x[0]) * t[0] * s[0] + 2. * t[2] * s[2];
  r[2] = e * ((x[2] * t[1] + x[1] * t[2]) * (x[2] * s[1] + x[1] * s[2]) + t[1] * s[2] + t[2] * s[1]);
}

}  // namespace toy

// -----------------------------------------------------------------------------
// The external operations
// -----------------------------------------------------------------------------

// What the operations record, and how FFToy differentiates
struct Trace
{
  static bool use_tangent;     // FFToy::has_tangent()
  static bool sparse_tangent;  // FFToy::tangent() leaves the tangent of y1 a structural zero where it is one
  static bool throw_tangent;   // FFToy::tangent() throws
  static int n_deriv;          // calls of FFToy::deriv()
  static int n_tangent;        // calls of FFToy::tangent()
  static int n_tangent2;       // calls of FFToyTan::tangent()
  static std::string pattern;  // operand tangents of the calls of FFToy::tangent(): 'T' given, '0' null
  static bool detached;        // FFToy::tangent() was given a tangent that is not attached to the DAG
  static unsigned max_sweeps;  // FFToy::has_tangent() declines beyond this number of sweeps
  static std::string asked;    // numbers of sweeps FFToy::has_tangent() was asked with
  static void
  reset()
  {
    n_deriv = n_tangent = n_tangent2 = 0;
    pattern.clear();
    asked.clear();
    detached = false;
  }
};
bool Trace::use_tangent    = true;
bool Trace::sparse_tangent = false;
bool Trace::throw_tangent  = false;
int Trace::n_deriv         = 0;
int Trace::n_tangent       = 0;
int Trace::n_tangent2      = 0;
std::string Trace::pattern;
bool Trace::detached = false;
unsigned Trace::max_sweeps = 1000;
std::string Trace::asked;

// Common part: evaluation in double, in FFVar (insertion in a DAG) and in FFDep
template <typename OP>
class FFToyBase : public FFOp
{
 public:
  FFToyBase() : FFOp(EXTERN) {}

  // Insert the operation with operands vVar and nRes results
  FFVar**
  insert(unsigned const nRes, std::vector<FFVar> const& vVar) const
  {
    data    = nullptr;
    owndata = false;
    return insert_external_operation(static_cast<OP const&>(*this), nRes, vVar.size(), vVar.data());
  }

  virtual void
  feval(std::type_info const& idU, unsigned const nRes, void* vRes, unsigned const nVar,
        void const* vVar, unsigned const* mVar) const
  {
    if (mc::same_type(idU, typeid(double)))
      return static_cast<OP const*>(this)->evald(static_cast<double*>(vRes),
                                                 static_cast<double const*>(vVar));
    if (mc::same_type(idU, typeid(FFVar)))
    {
      FFVar const* pVar = static_cast<FFVar const*>(vVar);
      FFVar** ppRes     = insert(nRes, std::vector<FFVar>(pVar, pVar + nVar));
      for (unsigned k = 0; k < nRes; ++k) static_cast<FFVar*>(vRes)[k] = *ppRes[k];
      return;
    }
    if (mc::same_type(idU, typeid(mc::FFDep)))
    {
      mc::FFDep* pRes       = static_cast<mc::FFDep*>(vRes);
      mc::FFDep const* pVar = static_cast<mc::FFDep const*>(vVar);
      pRes[0]               = 0;
      for (unsigned i = 0; i < nVar; ++i) pRes[0] += pVar[i];
      pRes[0].update(mc::FFDep::TYPE::N);
      for (unsigned k = 1; k < nRes; ++k) pRes[k] = pRes[0];
      return;
    }
    throw std::runtime_error("FFToy: no evaluation in this type");
  }

  bool
  commutative() const
  {
    return false;
  }
};

// Second directional derivative D2y[t,s]: operands (x, t, s)
class FFToyTan2 : public FFToyBase<FFToyTan2>
{
 public:
  void
  evald(double* r, double const* v) const
  {
    toy::tan2(v, v + 3, v + 6, r);
  }
  std::string
  name() const
  {
    return "ToyTan2";
  }
};

// Directional derivative J(x) t: operands (x, t).  Its own tangent along (dx, dt) is
// D2y[t,dx] + J(x) dt; it has no deriv(), so that backward mode through it is refused
class FFToyTan : public FFToyBase<FFToyTan>
{
 public:
  void
  evald(double* r, double const* v) const
  {
    toy::tan(v, v + 3, r);
  }
  std::string
  name() const
  {
    return "ToyTan";
  }
  bool
  has_tangent(unsigned const nSweep) const
  {
    return true;
  }
  void
  tangent(unsigned const nRes, FFVar const* vRes, unsigned const nVar, FFVar const* vVar,
          FFVar const* const* vTan, FFVar* vResTan) const
  {
    ++Trace::n_tangent2;
    bool dx = false, dt = false;
    std::vector<FFVar> x(vVar, vVar + 3), t(vVar + 3, vVar + 6), vdx(3, FFVar(0.)), vdt(3, FFVar(0.));
    for (unsigned i = 0; i < 3; ++i)
    {
      if (vTan[i])     { vdx[i] = *vTan[i];     dx = true; }
      if (vTan[3 + i]) { vdt[i] = *vTan[3 + i]; dt = true; }
    }
    if (dx)
    {
      std::vector<FFVar> v = x;
      v.insert(v.end(), t.begin(), t.end());
      v.insert(v.end(), vdx.begin(), vdx.end());
      FFVar** pp = FFToyTan2().insert(3, v);
      for (unsigned k = 0; k < 3; ++k) vResTan[k] = *pp[k];
    }
    if (dt)
    {
      std::vector<FFVar> v = x;
      v.insert(v.end(), vdt.begin(), vdt.end());
      FFVar** pp = FFToyTan().insert(3, v);
      for (unsigned k = 0; k < 3; ++k) vResTan[k] = (dx ? vResTan[k] + *pp[k] : *pp[k]);
    }
  }
};

// Jacobian: operands x, results J[k*3+i]
class FFToyJac : public FFToyBase<FFToyJac>
{
 public:
  void
  evald(double* r, double const* v) const
  {
    toy::jac(v, r);
  }
  std::string
  name() const
  {
    return "ToyJac";
  }
};

// Value: operands x
class FFToy : public FFToyBase<FFToy>
{
 public:
  void
  evald(double* r, double const* v) const
  {
    toy::val(v, r);
  }
  std::string
  name() const
  {
    return "Toy";
  }

  // Define the operation
  std::vector<FFVar>
  operator()(std::vector<FFVar> const& vVar) const
  {
    FFVar** pp = insert(3, vVar);
    return {*pp[0], *pp[1], *pp[2]};
  }

  // Partial derivatives: the Jacobian operation
  void
  deriv(unsigned const nRes, FFVar const* vRes, unsigned const nVar, FFVar const* vVar,
        FFVar** vDer) const
  {
    ++Trace::n_deriv;
    FFVar** pp = FFToyJac().insert(9, std::vector<FFVar>(vVar, vVar + nVar));
    for (unsigned k = 0; k < nRes; ++k)
      for (unsigned i = 0; i < nVar; ++i) vDer[k][i] = *pp[k * nVar + i];
  }

  // Tangents: the directional derivative operation
  bool
  has_tangent(unsigned const nSweep) const
  {
    if (Trace::use_tangent) Trace::asked += std::to_string(nSweep) + " ";
    return Trace::use_tangent && nSweep <= Trace::max_sweeps;
  }
  void
  tangent(unsigned const nRes, FFVar const* vRes, unsigned const nVar, FFVar const* vVar,
          FFVar const* const* vTan, FFVar* vResTan) const
  {
    ++Trace::n_tangent;
    if (!Trace::pattern.empty()) Trace::pattern += ' ';
    for (unsigned i = 0; i < nVar; ++i) Trace::pattern += (vTan[i] ? 'T' : '0');
    for (unsigned i = 0; i < nVar; ++i) Trace::detached = Trace::detached || (vTan[i] && !vTan[i]->dag());
    if (Trace::throw_tangent) throw std::runtime_error("FFToy::tangent");

    std::vector<FFVar> v(vVar, vVar + nVar);
    for (unsigned i = 0; i < nVar; ++i) v.push_back(vTan[i] ? *vTan[i] : FFVar(0.));
    FFVar** pp = FFToyTan().insert(3, v);
    for (unsigned k = 0; k < nRes; ++k) vResTan[k] = *pp[k];
    // y1 = sin(x0)+x2^2 does not depend on x1
    if (Trace::sparse_tangent && !vTan[0] && !vTan[2]) vResTan[1] = 0.;
  }
};

// -----------------------------------------------------------------------------
// Helpers
// -----------------------------------------------------------------------------

static std::vector<double>
evaluate(FFGraph& DAG, std::vector<FFVar> const& F, std::vector<FFVar> const& Q,
         std::vector<double> const& q)
{
  std::vector<double> f(F.size());
  if (!F.empty()) DAG.eval(F, f, Q, q);
  return f;
}

// Dense Jacobian of F in Q at q, by FAD with the tangent hook on or off
static std::vector<double>
fad(FFGraph& DAG, std::vector<FFVar> const& F, std::vector<FFVar> const& Q,
    std::vector<double> const& q, bool const hook)
{
  Trace::use_tangent = hook;
  Trace::reset();
  std::vector<FFVar> const dF = DAG.FAD(F, Q);
  Trace::use_tangent = true;
  return evaluate(DAG, dF, Q, q);
}

static std::vector<std::vector<double>> const points = {
    {0.3, -0.7, 1.1}, {-1.2, 0.4, 0.5}, {2.0, 1.5, -0.6}};

// -----------------------------------------------------------------------------
// Tests
// -----------------------------------------------------------------------------
int
main()
{
#if defined(MC__DAG_FADIFF) || defined(MC__DAG_BADIFF)
  // The tangent and adjoint hooks belong to the default differentiation engines: the legacy engines, in
  // FADBAD++ arithmetic, do not call them
  std::printf("DAG_tangent: skipped with MC__DAG_FADIFF or MC__DAG_BADIFF, the hooks belong to the default engines\n");
  return 0;
#endif

  using test::near;
  using test::Result;
  using test::run;

  std::printf("Tangent hook of external operations\n");

  run("value of the operation and of its derivative operations", []() -> Result {
    FFGraph DAG;
    std::vector<FFVar> X = DAG.add_vars(3);
    std::vector<FFVar> Y = FFToy()(X);
    std::vector<double> const x = points[0], y = evaluate(DAG, Y, X, x);
    double yref[3], h = 1e-6;
    toy::val(x.data(), yref);
    // directional derivatives against differences
    double t[3] = {0.5, -1., 2.}, s[3] = {1., 0.25, -0.5}, r[3], r2[3], xp[3], xm[3], yp[3], ym[3], rp[3], rm[3];
    toy::tan(x.data(), t, r);
    toy::tan2(x.data(), t, s, r2);
    for (int i = 0; i < 3; ++i) { xp[i] = x[i] + h * t[i]; xm[i] = x[i] - h * t[i]; }
    toy::val(xp, yp);
    toy::val(xm, ym);
    bool ok = near(y, {yref[0], yref[1], yref[2]});
    for (int k = 0; k < 3; ++k) ok = ok && near((yp[k] - ym[k]) / (2 * h), r[k], 1e-8);
    for (int i = 0; i < 3; ++i) { xp[i] = x[i] + h * s[i]; xm[i] = x[i] - h * s[i]; }
    toy::tan(xp, t, rp);
    toy::tan(xm, t, rm);
    for (int k = 0; k < 3; ++k) ok = ok && near((rp[k] - rm[k]) / (2 * h), r2[k], 1e-8);
    return {ok, "values or analytic derivatives of the test function"};
  });

  run("FAD: same derivatives through tangent() as through deriv()", []() -> Result {
    FFGraph DAG;
    std::vector<FFVar> Q = DAG.add_vars(3);
    std::vector<FFVar> Y = FFToy()({Q[0] * Q[1], Q[0] + 2. * Q[2], sqr(Q[1])});
    std::vector<FFVar> F = {Y[0] + Q[2] * Y[2], Y[1] * Y[1], exp(Y[0]) - Y[2]};
    for (auto const& q : points)
    {
      std::vector<double> const ref = fad(DAG, F, Q, q, false);
      int const nd = Trace::n_deriv, nt0 = Trace::n_tangent;
      std::vector<double> const der = fad(DAG, F, Q, q, true);
      if (nd != 1 || nt0 != 0) return {false, "deriv() path: " + std::to_string(nd) + " deriv, " + std::to_string(nt0) + " tangent calls"};
      if (Trace::n_deriv != 0 || Trace::n_tangent != 3)
        return {false, "tangent() path: " + std::to_string(Trace::n_deriv) + " deriv, " + std::to_string(Trace::n_tangent) + " tangent calls"};
      if (ref.size() != 9 || !near(der, ref)) return {false, "derivatives differ"};
    }
    return {true, ""};
  });

  run("FAD: analytic Jacobian, operands the independents themselves", []() -> Result {
    FFGraph DAG;
    std::vector<FFVar> X = DAG.add_vars(3);
    std::vector<FFVar> Y = FFToy()(X);
    for (auto const& x : points)
    {
      std::vector<double> J(9);
      toy::jac(x.data(), J.data());
      if (!near(fad(DAG, Y, X, x, true), J)) return {false, "tangent() path"};
      // the seeds are constants: they are handed over as constants of the DAG
      if (Trace::pattern != "T00 0T0 00T" || Trace::detached) return {false, "operand tangents: " + Trace::pattern};
      if (!near(fad(DAG, Y, X, x, false), J)) return {false, "deriv() path"};
    }
    return {true, ""};
  });

  run("tangent() receives the tangents of all operands, null for structural zeros", []() -> Result {
    FFGraph DAG;
    std::vector<FFVar> Q = DAG.add_vars(4);
    // x0 = q0 q1, x1 = 2 (constant), x2 = q1+q3; q2 does not reach the operation
    std::vector<FFVar> Y = FFToy()({Q[0] * Q[1], FFVar(2.), Q[1] + Q[3]});
    std::vector<FFVar> F = {Y[0] + Y[1] + Y[2] + Q[2]};
    std::vector<double> const q = {0.3, -0.7, 1.1, 0.4};
    std::vector<double> const ref = fad(DAG, F, Q, q, false);
    std::vector<double> const der = fad(DAG, F, Q, q, true);
    // one sweep per independent, restricted to its cone: q0 -> x0; q1 -> x0, x2; q2 -> none; q3 -> x2
    if (Trace::pattern != "T00 T0T 00T") return {false, "operand tangents: " + Trace::pattern};
    return {near(der, ref) && near(der[2], 1.), "derivatives differ"};
  });

  run("SFAD: same sparse Jacobian through tangent() as through deriv()", []() -> Result {
    FFGraph DAG;
    std::vector<FFVar> Q = DAG.add_vars(5);
    std::vector<FFVar> Y = FFToy()({Q[0] * Q[1], Q[0] + 2. * Q[2], sqr(Q[1])});
    std::vector<FFVar> F = {Y[0] + Q[4], Q[3] * Q[3], Y[1] * Y[2], Q[4] - Q[3]};
    std::vector<double> const q = {0.3, -0.7, 1.1, 0.4, -2.};
    Trace::use_tangent = false;
    auto [r0, c0, v0] = DAG.SFAD(F, Q);
    Trace::use_tangent = true;
    Trace::reset();
    auto [r1, c1, v1] = DAG.SFAD(F, Q);
    if (r0 != r1 || c0 != c1) return {false, "sparsity patterns differ"};
    if (Trace::n_tangent != 3) return {false, std::to_string(Trace::n_tangent) + " tangent calls for 3 independents in the cone"};
    // F0 in q0, q1, q2, q4; F1 in q3; F2 in q0, q1, q2; F3 in q3, q4
    if (r1.size() != 10) return {false, std::to_string(r1.size()) + " entries"};
    return {near(evaluate(DAG, v1, Q, q), evaluate(DAG, v0, Q, q)), "values differ"};
  });

  run("SFAD: a result tangent left at zero is a structural zero", []() -> Result {
    FFGraph DAG;
    std::vector<FFVar> Q = DAG.add_vars(2);
    // only x1 depends on q1, and y1 does not depend on x1
    std::vector<FFVar> Y = FFToy()({Q[0], Q[1], sqr(Q[0])});
    std::vector<double> const q = {0.3, -0.7};
    Trace::sparse_tangent = true;
    auto [r1, c1, v1] = DAG.SFAD(Y, Q);
    Trace::sparse_tangent = false;
    auto [r0, c0, v0] = DAG.SFAD(Y, Q);
    bool found1 = false, found0 = false;
    for (size_t k = 0; k < r1.size(); ++k) found1 = found1 || (r1[k] == 1 && c1[k] == 1);
    for (size_t k = 0; k < r0.size(); ++k) found0 = found0 || (r0[k] == 1 && c0[k] == 1);
    if (found1 || !found0 || r1.size() + 1 != r0.size()) return {false, "entry (y1, q1) of the pattern"};
    std::vector<double> const d1 = evaluate(DAG, v1, Q, q), d0 = evaluate(DAG, v0, Q, q);
    for (size_t k = 0, l = 0; k < r0.size(); ++k)
    {
      if (r0[k] == 1 && c0[k] == 1)
      {
        if (d0[k] != 0.) return {false, "value of the entry that is a structural zero"};
        continue;
      }
      if (!near(d1[l++], d0[k])) return {false, "values differ"};
    }
    return {true, ""};
  });

  run("DFAD: one call of tangent() for a direction", []() -> Result {
    FFGraph DAG;
    std::vector<FFVar> Q = DAG.add_vars(3), D = DAG.add_vars(3);
    std::vector<FFVar> Y = FFToy()({Q[0] * Q[1], Q[0] + 2. * Q[2], sqr(Q[1])});
    std::vector<FFVar> F = {Y[0] + Q[2] * Y[2], Y[1] * Y[1]};
    std::vector<FFVar> QD = Q;
    QD.insert(QD.end(), D.begin(), D.end());
    Trace::use_tangent = false;
    std::vector<FFVar> const ref = DAG.DFAD(F, Q, D);
    Trace::use_tangent = true;
    Trace::reset();
    std::vector<FFVar> const dir = DAG.DFAD(F, Q, D);
    if (Trace::n_tangent != 1 || Trace::n_deriv != 0 || Trace::pattern != "TTT")
      return {false, std::to_string(Trace::n_tangent) + " tangent calls, operand tangents " + Trace::pattern};
    for (auto const& q : points)
    {
      std::vector<double> qd = q;
      qd.insert(qd.end(), {0.5, -1., 2.});
      if (!near(evaluate(DAG, dir, QD, qd), evaluate(DAG, ref, QD, qd))) return {false, "directional derivatives differ"};
    }
    // constant direction with a zero component: its operand tangent is a structural zero
    std::vector<FFVar> const D0 = {FFVar(1.), FFVar(0.), FFVar(0.)};
    std::vector<FFVar> const Y0 = FFToy()(Q);
    Trace::reset();
    std::vector<FFVar> const dir0 = DAG.DFAD(Y0, Q, D0);
    std::vector<double> J(9);
    toy::jac(points[0].data(), J.data());
    return {Trace::pattern == "T00" && near(evaluate(DAG, dir0, Q, points[0]), {J[0], J[3], J[6]}),
            "constant direction: operand tangents " + Trace::pattern};
  });

  run("has_tangent() is asked with the number of sweeps, and may decline", []() -> Result {
    FFGraph DAG;
    std::vector<FFVar> Q = DAG.add_vars(5), D = DAG.add_vars(5);
    // the operation depends on q0, q1 and q3; F1 does not depend on the operation
    std::vector<FFVar> Y = FFToy()({Q[0] * Q[1], FFVar(2.), Q[1] + Q[3]});
    std::vector<FFVar> F = {Y[0] + Y[2] * Q[2], Q[4] * Q[2]};
    std::vector<double> const q = {0.3, -0.7, 1.1, 0.4, -2.};
    std::vector<double> const ref = fad(DAG, F, Q, q, false);
    std::vector<double> const der = fad(DAG, F, Q, q, true);
    // first with one sweep, then with the three independents that the operation depends on
    if (Trace::asked != "1 3 " || Trace::n_tangent != 3 || Trace::n_deriv != 0) return {false, "asked with " + Trace::asked};
    if (!near(der, ref)) return {false, "derivatives differ"};
    // declining for three sweeps: the partial derivatives instead
    Trace::max_sweeps = 2;
    std::vector<double> const der2 = fad(DAG, F, Q, q, true);
    Trace::max_sweeps = 1000;
    if (Trace::asked != "1 3 " || Trace::n_tangent != 0 || Trace::n_deriv != 1) return {false, "declining: asked with " + Trace::asked};
    if (!near(der2, ref)) return {false, "derivatives differ when declining"};
    // a directional derivative is one sweep, whatever the number of independents
    Trace::max_sweeps = 1;
    Trace::reset();
    std::vector<FFVar> const dir = DAG.DFAD(F, Q, D);
    Trace::max_sweeps = 1000;
    if (Trace::asked != "1 " || Trace::n_tangent != 1 || Trace::n_deriv != 0) return {false, "directional: asked with " + Trace::asked};
    // an operation that a single independent reaches is not asked again
    std::vector<FFVar> const Y1 = FFToy()({Q[4], sqr(Q[4]), FFVar(1.)});
    Trace::reset();
    std::vector<FFVar> const dY1 = DAG.FAD(Y1, Q);
    return {Trace::asked == "1 " && Trace::n_tangent == 1, "single independent: asked with " + Trace::asked};
  });

  run("repeated operands", []() -> Result {
    FFGraph DAG;
    std::vector<FFVar> Q = DAG.add_vars(2);
    std::vector<FFVar> Y = FFToy()({Q[0], Q[0], Q[1]});
    for (auto const& p : points)
    {
      std::vector<double> const q = {p[0], p[2]};
      std::vector<double> const ref = fad(DAG, Y, Q, q, false);
      std::vector<double> const der = fad(DAG, Y, Q, q, true);
      if (Trace::pattern != "TT0 00T") return {false, "operand tangents: " + Trace::pattern};
      double J[9], x[3] = {q[0], q[0], q[1]};
      toy::jac(x, J);
      if (!near(der, ref) || !near(der, {J[0] + J[1], J[2], J[3] + J[4], J[5], J[6] + J[7], J[8]}))
        return {false, "derivatives differ"};
    }
    return {true, ""};
  });

  run("operation outside the cone of every independent is not called", []() -> Result {
    FFGraph DAG;
    std::vector<FFVar> Q = DAG.add_vars(2), P = DAG.add_vars(3);
    std::vector<FFVar> Y = FFToy()(P);
    std::vector<FFVar> F = {Y[0] * Q[0] + Q[1]};
    std::vector<FFVar> QP = Q;
    QP.insert(QP.end(), P.begin(), P.end());
    Trace::reset();
    std::vector<FFVar> const dF = DAG.FAD(F, Q);
    if (Trace::n_tangent != 0 || Trace::n_deriv != 0) return {false, "operation called"};
    std::vector<double> const qp = {0.3, -0.7, 1.1, 0.4, -2.};
    return {near(evaluate(DAG, dF, QP, qp), {1.1 * 0.4, 1.}), "derivatives"};
  });

  run("second order: tangent() of the directional derivative operation", []() -> Result {
    FFGraph DAG;
    std::vector<FFVar> X = DAG.add_vars(3);
    std::vector<FFVar> Y = FFToy()(X);
    std::vector<FFVar> F = {Y[0] + Y[1] + Y[2]};
    Trace::reset();
    std::vector<FFVar> const dF  = DAG.FAD(F, X);
    std::vector<FFVar> const d2F = DAG.FAD(dF, X);
    // first pass: one FFToyTan operation per independent; second pass: each of them in each sweep
    if (Trace::n_tangent != 3 || Trace::n_tangent2 != 9 || Trace::n_deriv != 0)
      return {false, std::to_string(Trace::n_tangent) + " and " + std::to_string(Trace::n_tangent2) + " tangent calls"};
    for (auto const& x : points)
    {
      // Hessian of y0+y1+y2: entry (i,j) is the sum over k of D2y_k[e_i,e_j]
      std::vector<double> H(9);
      for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
        {
          double ei[3] = {0., 0., 0.}, ej[3] = {0., 0., 0.}, r[3];
          ei[i] = ej[j] = 1.;
          toy::tan2(x.data(), ei, ej, r);
          H[i * 3 + j] = r[0] + r[1] + r[2];
        }
      if (!near(evaluate(DAG, d2F, X, x), H)) return {false, "Hessian"};
    }
    return {true, ""};
  });

  run("second order through expressions of the independents", []() -> Result {
    FFGraph DAG;
    std::vector<FFVar> Q = DAG.add_vars(2);
    std::vector<FFVar> Y = FFToy()({Q[0] * Q[1], Q[0] + 2. * Q[1], sqr(Q[1])});
    std::vector<FFVar> F = {Y[0] * Y[2] + Y[1]};
    std::vector<FFVar> const dF  = DAG.FAD(F, Q);
    std::vector<FFVar> const d2F = DAG.FAD(dF, Q);
    // against central differences of the first derivatives
    double const h = 1e-5;
    for (auto const& p : points)
    {
      std::vector<double> const q = {p[0], p[1]}, H = evaluate(DAG, d2F, Q, q);
      for (int j = 0; j < 2; ++j)
      {
        std::vector<double> qp = q, qm = q;
        qp[j] += h;
        qm[j] -= h;
        std::vector<double> const gp = evaluate(DAG, dF, Q, qp), gm = evaluate(DAG, dF, Q, qm);
        for (int i = 0; i < 2; ++i)
          if (!near((gp[i] - gm[i]) / (2 * h), H[i * 2 + j], 1e-6)) return {false, "Hessian against differences"};
      }
      if (!near(H[1], H[2], 1e-10)) return {false, "Hessian not symmetric"};
    }
    return {true, ""};
  });

  run("BAD: backward mode still calls deriv()", []() -> Result {
    FFGraph DAG;
    std::vector<FFVar> Q = DAG.add_vars(3);
    std::vector<FFVar> Y = FFToy()({Q[0] * Q[1], Q[0] + 2. * Q[2], sqr(Q[1])});
    std::vector<FFVar> F = {Y[0] + Q[2] * Y[2], Y[1] * Y[1]};
    Trace::reset();
    std::vector<FFVar> const dFb = DAG.BAD(F, Q);
    if (Trace::n_deriv != 1 || Trace::n_tangent != 0)
      return {false, std::to_string(Trace::n_deriv) + " deriv, " + std::to_string(Trace::n_tangent) + " tangent calls"};
    for (auto const& q : points)
      if (!near(evaluate(DAG, dFb, Q, q), fad(DAG, F, Q, q, true))) return {false, "derivatives differ"};
    return {true, ""};
  });

  run("BAD through an operation with tangent() only is refused, and the DAG stays usable", []() -> Result {
    FFGraph DAG;
    std::vector<FFVar> X = DAG.add_vars(3);
    std::vector<FFVar> Y = FFToy()(X);
    std::vector<FFVar> const dY = DAG.FAD(Y, X);  // FFToyTan operations, which have no deriv()
    bool refused = false;
    try
    {
      std::vector<FFVar> const d2Y = DAG.BAD(dY, X);
    }
    catch (mc::FFBase::Exceptions& e)
    {
      refused = (e.ierr() == mc::FFBase::Exceptions::EXTERN);
    }
    if (!refused) return {false, "not refused"};
    std::vector<double> J(9);
    toy::jac(points[1].data(), J.data());
    return {near(evaluate(DAG, dY, X, points[1]), J) && near(fad(DAG, Y, X, points[1], true), J), "after the refusal"};
  });

  run("exception in tangent(): propagated, and the DAG stays usable", []() -> Result {
    FFGraph DAG;
    std::vector<FFVar> X = DAG.add_vars(3);
    std::vector<FFVar> Y = FFToy()(X);
    bool thrown = false;
    Trace::throw_tangent = true;
    try
    {
      std::vector<FFVar> const dY = DAG.FAD(Y, X);
    }
    catch (std::runtime_error const& e)
    {
      thrown = (std::string(e.what()) == "FFToy::tangent");
    }
    Trace::throw_tangent = false;
    if (!thrown) return {false, "not propagated"};
    std::vector<double> y(3), J(9);
    toy::val(points[2].data(), y.data());
    toy::jac(points[2].data(), J.data());
    return {near(evaluate(DAG, Y, X, points[2]), y) && near(fad(DAG, Y, X, points[2], true), J), "after the exception"};
  });

  run("default: an operation without tangent() is differentiated through deriv()", []() -> Result {
    // FFOp::has_tangent() is false unless overridden, and FFOp::tangent() throws
    FFToyJac Op;
    FFVar res, var, tan;
    FFVar const* ptan = &var;
    bool thrown       = false;
    try
    {
      Op.FFOp::tangent(1, &res, 1, &var, &ptan, &tan);
    }
    catch (mc::FFBase::Exceptions& e)
    {
      thrown = (e.ierr() == mc::FFBase::Exceptions::EXTERN);
    }
    return {!Op.has_tangent(1) && thrown, "defaults of the hook"};
  });

  std::printf("\n=== Results: %d/%d passed ===\n", test::passed, test::total);
  return test::failed ? 1 : 0;
}
