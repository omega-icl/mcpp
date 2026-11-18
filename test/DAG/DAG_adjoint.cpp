// =============================================================================
// DAG_adjoint.cpp
// Unit tests for the adjoint hook of external operations, FFOp::has_adjoint()
// and FFOp::adjoint(): an operation that propagates adjoints itself in the
// backward-mode symbolic differentiation of mc::FFGraph (BAD, SBAD, DBAD), in
// place of the partial derivatives returned by FFOp::deriv().
//
// The operation under test is y = ( x0 x1, sin(x0)+x2^2, exp(x1 x2) ), with
//  - FFToy     its value                                      [3 operands, 3 results]
//  - FFToyJac  its Jacobian, the result of deriv()            [3 operands, 9 results]
//  - FFToyAdj  its adjoint J(x)^T w, the result of adjoint(); [6 operands, 3 results]
//              it has tangent and adjoint hooks of its own
//  - FFToyTan  its directional derivative J(x) t, the result  [6 operands, 3 results]
//              of tangent(); it has an adjoint hook of its own
//  - FFToyAdj2 its second-order adjoint H_w(x) s, where H_w   [9 operands, 3 results]
//              is the Hessian of w^T y
// The tangent and adjoint operations are closed under both modes:
//   tangent of FFToyAdj(x,w) along (dx,dw) = FFToyAdj2(x,w,dx) + FFToyAdj(x,dw)
//   adjoint of FFToyAdj(x,w) for r         = FFToyAdj2(x,w,r) in x, FFToyTan(x,r) in w
//   adjoint of FFToyTan(x,t) for r         = FFToyAdj2(x,r,t) in x, FFToyAdj(x,r) in t
// The derivatives obtained through the hook are compared with those obtained
// through deriv() and with analytic values.
//
// Build (assuming MC++ headers are on the include path):
//   g++ -std=c++17 -O2 -I<path-to-mc++> -o DAG_adjoint DAG_adjoint.cpp
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

// J(x)^T w
inline void
adj(double const* x, double const* w, double* r)
{
  double const e = std::exp(x[1] * x[2]);
  r[0] = x[1] * w[0] + std::cos(x[0]) * w[1];
  r[1] = x[0] * w[0] + x[2] * e * w[2];
  r[2] = 2. * x[2] * w[1] + x[1] * e * w[2];
}

// H_w(x) s, with H_w the Hessian of w^T y
inline void
adj2(double const* x, double const* w, double const* s, double* r)
{
  double const e = std::exp(x[1] * x[2]), c = 1. + x[1] * x[2];
  r[0] = w[0] * s[1] - w[1] * std::sin(x[0]) * s[0];
  r[1] = w[0] * s[0] + w[2] * e * (x[2] * x[2] * s[1] + c * s[2]);
  r[2] = 2. * w[1] * s[2] + w[2] * e * (c * s[1] + x[1] * x[1] * s[2]);
}

}  // namespace toy

// -----------------------------------------------------------------------------
// The external operations
// -----------------------------------------------------------------------------

// What the operations record, and how FFToy differentiates
struct Trace
{
  static bool use_adjoint;     // FFToy::has_adjoint()
  static bool use_tangent;     // FFToy::has_tangent()
  static bool sparse_adjoint;  // FFToy::adjoint() leaves the contribution to x1 a structural zero where it is one
  static bool throw_adjoint;   // FFToy::adjoint() throws
  static int n_deriv;          // calls of FFToy::deriv()
  static int n_adjoint;        // calls of FFToy::adjoint()
  static int n_tangent;        // calls of FFToy::tangent()
  static int n_adjadj;         // calls of FFToyAdj::adjoint()
  static int n_adjtan;         // calls of FFToyAdj::tangent()
  static int n_tanadj;         // calls of FFToyTan::adjoint()
  static std::string pattern;  // result adjoints of the calls of FFToy::adjoint(): 'A' given, '0' null
  static bool detached;        // FFToy::adjoint() was given an adjoint that is not attached to the DAG
  static unsigned max_sweeps;  // FFToy::has_adjoint() declines beyond this number of sweeps
  static std::string asked;    // numbers of sweeps FFToy::has_adjoint() was asked with
  static void
  reset()
  {
    n_deriv = n_adjoint = n_tangent = n_adjadj = n_adjtan = n_tanadj = 0;
    pattern.clear();
    asked.clear();
    detached = false;
  }
};
bool Trace::use_adjoint    = true;
bool Trace::use_tangent    = false;
bool Trace::sparse_adjoint = false;
bool Trace::throw_adjoint  = false;
int Trace::n_deriv         = 0;
int Trace::n_adjoint       = 0;
int Trace::n_tangent       = 0;
int Trace::n_adjadj        = 0;
int Trace::n_adjtan        = 0;
int Trace::n_tanadj        = 0;
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

  // Insert the operation with the operands in v1, v2 and v3, a null pointer standing for zero
  FFVar**
  insert(unsigned const nRes, FFVar const* v1, FFVar const* const* v2 = nullptr,
         FFVar const* const* v3 = nullptr) const
  {
    std::vector<FFVar> v(v1, v1 + 3);
    for (auto const* p : {v2, v3})
      for (unsigned i = 0; p && i < 3; ++i) v.push_back(p[i] ? *p[i] : FFVar(0.));
    return insert(nRes, v);
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

// Pointers to the elements of an array, for the operand arguments of insert()
static std::vector<FFVar const*>
pointers(FFVar const* v, unsigned const n = 3)
{
  std::vector<FFVar const*> p(n);
  for (unsigned i = 0; i < n; ++i) p[i] = v + i;
  return p;
}

// Second-order adjoint H_w(x) s: operands (x, w, s).  No derivative of any kind
class FFToyAdj2 : public FFToyBase<FFToyAdj2>
{
 public:
  void
  evald(double* r, double const* v) const
  {
    toy::adj2(v, v + 3, v + 6, r);
  }
  std::string
  name() const
  {
    return "ToyAdj2";
  }
};

class FFToyAdj;

// Directional derivative J(x) t: operands (x, t).  Its adjoint for r is H_r(x) t in x and
// J(x)^T r in t; it has neither tangent() nor deriv()
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
  has_adjoint(unsigned const nSweep) const
  {
    return true;
  }
  void adjoint(unsigned const nRes, FFVar const* vRes, unsigned const nVar, FFVar const* vVar,
               FFVar const* const* vResAdj, FFVar* vVarAdj) const;
};

// Adjoint J(x)^T w: operands (x, w).  Its tangent along (dx, dw) is H_w(x) dx + J(x)^T dw; its
// adjoint for r is H_w(x) r in x and J(x) r in w; it has no deriv()
class FFToyAdj : public FFToyBase<FFToyAdj>
{
 public:
  void
  evald(double* r, double const* v) const
  {
    toy::adj(v, v + 3, r);
  }
  std::string
  name() const
  {
    return "ToyAdj";
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
    ++Trace::n_adjtan;
    bool const dx = (vTan[0] || vTan[1] || vTan[2]), dw = (vTan[3] || vTan[4] || vTan[5]);
    auto const w = pointers(vVar + 3);
    if (dx)
    {
      FFVar** pp = FFToyAdj2().insert(3, vVar, w.data(), vTan);
      for (unsigned i = 0; i < 3; ++i) vResTan[i] = *pp[i];
    }
    if (dw)
    {
      FFVar** pp = FFToyAdj().insert(3, vVar, vTan + 3);
      for (unsigned i = 0; i < 3; ++i) vResTan[i] = (dx ? vResTan[i] + *pp[i] : *pp[i]);
    }
  }
  bool
  has_adjoint(unsigned const nSweep) const
  {
    return true;
  }
  void
  adjoint(unsigned const nRes, FFVar const* vRes, unsigned const nVar, FFVar const* vVar,
          FFVar const* const* vResAdj, FFVar* vVarAdj) const
  {
    ++Trace::n_adjadj;
    auto const w = pointers(vVar + 3);
    FFVar** ppx  = FFToyAdj2().insert(3, vVar, w.data(), vResAdj);
    FFVar** ppw  = FFToyTan().insert(3, vVar, vResAdj);
    for (unsigned i = 0; i < 3; ++i)
    {
      vVarAdj[i]     = *ppx[i];
      vVarAdj[3 + i] = *ppw[i];
    }
  }
};

inline void
FFToyTan::adjoint(unsigned const nRes, FFVar const* vRes, unsigned const nVar, FFVar const* vVar,
                  FFVar const* const* vResAdj, FFVar* vVarAdj) const
{
  ++Trace::n_tanadj;
  auto const t = pointers(vVar + 3);
  FFVar** ppx  = FFToyAdj2().insert(3, vVar, vResAdj, t.data());
  FFVar** ppt  = FFToyAdj().insert(3, vVar, vResAdj);
  for (unsigned i = 0; i < 3; ++i)
  {
    vVarAdj[i]     = *ppx[i];
    vVarAdj[3 + i] = *ppt[i];
  }
}

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

  // Adjoints: the adjoint operation
  bool
  has_adjoint(unsigned const nSweep) const
  {
    if (Trace::use_adjoint) Trace::asked += std::to_string(nSweep) + " ";
    return Trace::use_adjoint && nSweep <= Trace::max_sweeps;
  }
  void
  adjoint(unsigned const nRes, FFVar const* vRes, unsigned const nVar, FFVar const* vVar,
          FFVar const* const* vResAdj, FFVar* vVarAdj) const
  {
    ++Trace::n_adjoint;
    if (!Trace::pattern.empty()) Trace::pattern += ' ';
    for (unsigned k = 0; k < nRes; ++k) Trace::pattern += (vResAdj[k] ? 'A' : '0');
    for (unsigned k = 0; k < nRes; ++k) Trace::detached = Trace::detached || (vResAdj[k] && !vResAdj[k]->dag());
    if (Trace::throw_adjoint) throw std::runtime_error("FFToy::adjoint");

    FFVar** pp = FFToyAdj().insert(3, vVar, vResAdj);
    for (unsigned i = 0; i < nVar; ++i) vVarAdj[i] = *pp[i];
    // y1 = sin(x0)+x2^2 does not depend on x1
    if (Trace::sparse_adjoint && !vResAdj[0] && !vResAdj[2]) vVarAdj[1] = 0.;
  }

  // Tangents: the directional derivative operation
  bool
  has_tangent(unsigned const nSweep) const
  {
    return Trace::use_tangent;
  }
  void
  tangent(unsigned const nRes, FFVar const* vRes, unsigned const nVar, FFVar const* vVar,
          FFVar const* const* vTan, FFVar* vResTan) const
  {
    ++Trace::n_tangent;
    FFVar** pp = FFToyTan().insert(3, vVar, vTan);
    for (unsigned k = 0; k < nRes; ++k) vResTan[k] = *pp[k];
  }
};

// A linear operator, y = 3 x: without an adjoint, backward mode refuses it for want of its transpose
class FFScale : public FFOp
{
 public:
  static bool use_adjoint;
  static int n_adjoint, n_apply;  // calls of adjoint() and of apply()
  FFScale() : FFOp(EXTERN) {}
  FFVar
  operator()(FFVar const& x) const
  {
    data    = nullptr;
    owndata = false;
    return *insert_external_operation(*this, 1, 1, &x)[0];
  }
  virtual void
  feval(std::type_info const& idU, unsigned const nRes, void* vRes, unsigned const nVar,
        void const* vVar, unsigned const* mVar) const
  {
    if (mc::same_type(idU, typeid(double)))
      *static_cast<double*>(vRes) = 3. * *static_cast<double const*>(vVar);
    else if (mc::same_type(idU, typeid(FFVar)))
      *static_cast<FFVar*>(vRes) = operator()(*static_cast<FFVar const*>(vVar));
    else
      throw std::runtime_error("FFScale: no evaluation in this type");
  }
  bool
  linear(unsigned const iout, unsigned const iin) const
  {
    return true;
  }
  // Coefficients for the non-linear pairs: there is none
  void
  deriv(unsigned const nRes, FFVar const* vRes, unsigned const nVar, FFVar const* vVar,
        FFVar** vDer) const
  {
    vDer[0][0] = 0.;
  }
  FFVar
  apply(unsigned const iout, unsigned const iin, FFVar const& tangent) const
  {
    ++n_apply;
    // the tangent is a variable or a constant of the DAG, here the only operand of the operation
    // that is inserted
    return FFScale()(tangent);
  }
  bool
  has_adjoint(unsigned const nSweep) const
  {
    return use_adjoint;
  }

  void
  adjoint(unsigned const nRes, FFVar const* vRes, unsigned const nVar, FFVar const* vVar,
          FFVar const* const* vResAdj, FFVar* vVarAdj) const
  {
    ++n_adjoint;
    // the operator is its own transpose; the adjoint is a variable or a constant of the DAG
    vVarAdj[0] = FFScale()(*vResAdj[0]);
  }
  std::string
  name() const
  {
    return "Scale";
  }
  bool
  commutative() const
  {
    return false;
  }
};
bool FFScale::use_adjoint = true;
int FFScale::n_adjoint     = 0;
int FFScale::n_apply       = 0;

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

// Dense Jacobian of F in Q at q, by BAD with the adjoint hook on or off
static std::vector<double>
bad(FFGraph& DAG, std::vector<FFVar> const& F, std::vector<FFVar> const& Q,
    std::vector<double> const& q, bool const hook)
{
  Trace::use_adjoint = hook;
  Trace::reset();
  std::vector<FFVar> const dF = DAG.BAD(F, Q);
  Trace::use_adjoint = true;
  return evaluate(DAG, dF, Q, q);
}

static std::string
counts()
{
  return std::to_string(Trace::n_deriv) + " deriv, " + std::to_string(Trace::n_adjoint) + " adjoint, " +
         std::to_string(Trace::n_tangent) + " tangent calls; of the derivative operations " +
         std::to_string(Trace::n_adjadj) + " adjoint-adjoint, " + std::to_string(Trace::n_adjtan) +
         " adjoint-tangent, " + std::to_string(Trace::n_tanadj) + " tangent-adjoint";
}

static std::vector<std::vector<double>> const points = {
    {0.3, -0.7, 1.1}, {-1.2, 0.4, 0.5}, {2.0, 1.5, -0.6}};

// Hessian of w^T y at x, row by row
static std::vector<double>
hessian(std::vector<double> const& x, std::vector<double> const& w)
{
  std::vector<double> H(9);
  for (int j = 0; j < 3; ++j)
  {
    double ej[3] = {0., 0., 0.}, r[3];
    ej[j] = 1.;
    toy::adj2(x.data(), w.data(), ej, r);
    for (int i = 0; i < 3; ++i) H[i * 3 + j] = r[i];
  }
  return H;
}

// -----------------------------------------------------------------------------
// Tests
// -----------------------------------------------------------------------------
int
main()
{
#if defined(MC__DAG_FADIFF) || defined(MC__DAG_BADIFF)
  // The tangent and adjoint hooks belong to the default differentiation engines: the legacy engines, in
  // FADBAD++ arithmetic, do not call them
  std::printf("DAG_adjoint: skipped with MC__DAG_FADIFF or MC__DAG_BADIFF, the hooks belong to the default engines\n");
  return 0;
#endif

  using test::near;
  using test::Result;
  using test::run;

  std::printf("Adjoint hook of external operations\n");

  run("value of the operation and of its derivative operations", []() -> Result {
    FFGraph DAG;
    std::vector<FFVar> X = DAG.add_vars(3);
    std::vector<FFVar> Y = FFToy()(X);
    std::vector<double> const x = points[0], y = evaluate(DAG, Y, X, x);
    double yref[3], J[9], h = 1e-6;
    toy::val(x.data(), yref);
    toy::jac(x.data(), J);
    double w[3] = {0.5, -1., 2.}, s[3] = {1., 0.25, -0.5}, r[3], r2[3], t[3], xp[3], xm[3], rp[3], rm[3];
    toy::adj(x.data(), w, r);
    toy::adj2(x.data(), w, s, r2);
    toy::tan(x.data(), s, t);
    bool ok = near(y, {yref[0], yref[1], yref[2]});
    // J^T w and J s against the Jacobian
    for (int i = 0; i < 3; ++i) ok = ok && near(r[i], J[i] * w[0] + J[3 + i] * w[1] + J[6 + i] * w[2]);
    for (int k = 0; k < 3; ++k) ok = ok && near(t[k], J[3 * k] * s[0] + J[3 * k + 1] * s[1] + J[3 * k + 2] * s[2]);
    // H_w s against differences of J^T w along s
    for (int i = 0; i < 3; ++i) { xp[i] = x[i] + h * s[i]; xm[i] = x[i] - h * s[i]; }
    toy::adj(xp, w, rp);
    toy::adj(xm, w, rm);
    for (int i = 0; i < 3; ++i) ok = ok && near((rp[i] - rm[i]) / (2 * h), r2[i], 1e-8);
    // H_w symmetric
    std::vector<double> const H = hessian(x, {w[0], w[1], w[2]});
    ok = ok && near(H[1], H[3]) && near(H[2], H[6]) && near(H[5], H[7]);
    return {ok, "values or analytic derivatives of the test function"};
  });

  run("BAD: same derivatives through adjoint() as through deriv()", []() -> Result {
    FFGraph DAG;
    std::vector<FFVar> Q = DAG.add_vars(3);
    std::vector<FFVar> Y = FFToy()({Q[0] * Q[1], Q[0] + 2. * Q[2], sqr(Q[1])});
    std::vector<FFVar> F = {Y[0] + Q[2] * Y[2], Y[1] * Y[1], exp(Y[0]) - Y[2]};
    for (auto const& q : points)
    {
      std::vector<double> const ref = bad(DAG, F, Q, q, false);
      int const nd = Trace::n_deriv, na0 = Trace::n_adjoint;
      std::vector<double> const der = bad(DAG, F, Q, q, true);
      if (nd != 1 || na0 != 0) return {false, "deriv() path: " + std::to_string(nd) + " deriv, " + std::to_string(na0) + " adjoint calls"};
      if (Trace::n_deriv != 0 || Trace::n_adjoint != 3) return {false, "adjoint() path: " + counts()};
      if (ref.size() != 9 || !near(der, ref)) return {false, "derivatives differ"};
    }
    return {true, ""};
  });

  run("BAD: analytic Jacobian, operands the independents themselves", []() -> Result {
    FFGraph DAG;
    std::vector<FFVar> X = DAG.add_vars(3);
    std::vector<FFVar> Y = FFToy()(X);
    for (auto const& x : points)
    {
      std::vector<double> J(9);
      toy::jac(x.data(), J.data());
      if (!near(bad(DAG, Y, X, x, true), J)) return {false, "adjoint() path"};
      // the seeds are constants: they are handed over as constants of the DAG
      if (Trace::pattern != "A00 0A0 00A" || Trace::detached) return {false, "result adjoints: " + Trace::pattern};
      if (!near(bad(DAG, Y, X, x, false), J)) return {false, "deriv() path"};
    }
    return {true, ""};
  });

  run("adjoint() receives the adjoints of all results, null for structural zeros", []() -> Result {
    FFGraph DAG;
    std::vector<FFVar> Q = DAG.add_vars(4);
    std::vector<FFVar> Y = FFToy()({Q[0] * Q[1], FFVar(2.), Q[1] + Q[3]});
    // F0 from y0; F1 from y1 and y2; F2 not from the operation; F3 from y0 and, through F0, y0 again
    std::vector<FFVar> F = {Y[0] + Q[2], Y[1] * Y[2], Q[2] * Q[3], FFVar()};
    F[3] = F[0] * Y[0];
    std::vector<double> const q = {0.3, -0.7, 1.1, 0.4};
    std::vector<double> const ref = bad(DAG, F, Q, q, false);
    std::vector<double> const der = bad(DAG, F, Q, q, true);
    // one sweep per dependent, restricted to its cone
    if (Trace::pattern != "A00 0AA A00") return {false, "result adjoints: " + Trace::pattern};
    return {ref.size() == 16 && near(der, ref), "derivatives differ"};
  });

  run("dependents that are results of the operation, one consumed by another", []() -> Result {
    FFGraph DAG;
    std::vector<FFVar> Q = DAG.add_vars(3);
    std::vector<FFVar> Y = FFToy()({Q[0] * Q[1], Q[0] + 2. * Q[2], sqr(Q[1])});
    std::vector<FFVar> F = {Y[0], Y[0] * Y[1], Q[1], Y[2]};
    for (auto const& q : points)
    {
      std::vector<double> const ref = bad(DAG, F, Q, q, false);
      std::vector<double> const der = bad(DAG, F, Q, q, true);
      if (Trace::pattern != "A00 AA0 00A") return {false, "result adjoints: " + Trace::pattern};
      if (!near(der, ref)) return {false, "derivatives differ"};
    }
    return {true, ""};
  });

  run("composition of two operations", []() -> Result {
    FFGraph DAG;
    std::vector<FFVar> Q = DAG.add_vars(3);
    std::vector<FFVar> Y1 = FFToy()({Q[0] * Q[1], Q[0] + 2. * Q[2], 0.1 * Q[1]});
    std::vector<FFVar> Y2 = FFToy()({0.5 * Y1[0], Y1[1] - Q[0], 0.2 * Y1[2] * Y1[0]});
    std::vector<FFVar> F  = {Y2[0] + Y1[1], Y2[1] * Y2[2]};
    for (auto const& q : points)
    {
      std::vector<double> const ref = bad(DAG, F, Q, q, false);
      std::vector<double> const der = bad(DAG, F, Q, q, true);
      // each sweep goes through the second operation, then through the first
      if (Trace::n_adjoint != 4) return {false, counts()};
      if (!near(der, ref)) return {false, "derivatives differ"};
    }
    return {true, ""};
  });

  run("SBAD: same sparse Jacobian through adjoint() as through deriv()", []() -> Result {
    FFGraph DAG;
    std::vector<FFVar> Q = DAG.add_vars(5);
    std::vector<FFVar> Y = FFToy()({Q[0] * Q[1], Q[0] + 2. * Q[2], sqr(Q[1])});
    std::vector<FFVar> F = {Y[0] + Q[4], Q[3] * Q[3], Y[1] * Y[2], Q[4] - Q[3]};
    std::vector<double> const q = {0.3, -0.7, 1.1, 0.4, -2.};
    Trace::use_adjoint = false;
    auto [r0, c0, v0] = DAG.SBAD(F, Q);
    Trace::use_adjoint = true;
    Trace::reset();
    auto [r1, c1, v1] = DAG.SBAD(F, Q);
    if (r0 != r1 || c0 != c1) return {false, "sparsity patterns differ"};
    if (Trace::n_adjoint != 2) return {false, std::to_string(Trace::n_adjoint) + " adjoint calls for 2 dependents computed from the operation"};
    // F0 in q0, q1, q2, q4; F1 in q3; F2 in q0, q1, q2; F3 in q3, q4
    if (r1.size() != 10) return {false, std::to_string(r1.size()) + " entries"};
    return {near(evaluate(DAG, v1, Q, q), evaluate(DAG, v0, Q, q)), "values differ"};
  });

  run("SBAD: a contribution left at zero is a structural zero", []() -> Result {
    FFGraph DAG;
    std::vector<FFVar> Q = DAG.add_vars(2);
    // only x1 depends on q1, and y1 does not depend on x1
    std::vector<FFVar> Y = FFToy()({Q[0], Q[1], sqr(Q[0])});
    std::vector<double> const q = {0.3, -0.7};
    Trace::sparse_adjoint = true;
    auto [r1, c1, v1] = DAG.SBAD(Y, Q);
    Trace::sparse_adjoint = false;
    auto [r0, c0, v0] = DAG.SBAD(Y, Q);
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

  run("DBAD: one call of adjoint() for a vector of weights", []() -> Result {
    FFGraph DAG;
    std::vector<FFVar> Q = DAG.add_vars(3), W = DAG.add_vars(2);
    std::vector<FFVar> Y = FFToy()({Q[0] * Q[1], Q[0] + 2. * Q[2], sqr(Q[1])});
    std::vector<FFVar> F = {Y[0] + Q[2] * Y[2], Y[1] * Y[1]};
    std::vector<FFVar> QW = Q;
    QW.insert(QW.end(), W.begin(), W.end());
    Trace::use_adjoint = false;
    std::vector<FFVar> const ref = DAG.DBAD(F, W, Q);
    Trace::use_adjoint = true;
    Trace::reset();
    std::vector<FFVar> const dir = DAG.DBAD(F, W, Q);
    if (Trace::n_adjoint != 1 || Trace::n_deriv != 0 || Trace::pattern != "AAA") return {false, counts() + ", result adjoints " + Trace::pattern};
    for (auto const& q : points)
    {
      std::vector<double> qw = q;
      qw.insert(qw.end(), {0.5, -2.});
      if (!near(evaluate(DAG, dir, QW, qw), evaluate(DAG, ref, QW, qw))) return {false, "weighted derivatives differ"};
    }
    // constant weights with zero components: their result adjoints are structural zeros
    std::vector<FFVar> const W0 = {FFVar(1.), FFVar(0.), FFVar(0.)};
    std::vector<FFVar> const Y0 = FFToy()(Q);
    Trace::reset();
    std::vector<FFVar> const dir0 = DAG.DBAD(Y0, W0, Q);
    std::vector<double> J(9);
    toy::jac(points[0].data(), J.data());
    return {Trace::pattern == "A00" && near(evaluate(DAG, dir0, Q, points[0]), {J[0], J[1], J[2]}),
            "constant weights: result adjoints " + Trace::pattern};
  });

  run("has_adjoint() is asked with the number of sweeps, and may decline", []() -> Result {
    FFGraph DAG;
    std::vector<FFVar> Q = DAG.add_vars(3), W = DAG.add_vars(5);
    std::vector<FFVar> Y = FFToy()({Q[0] * Q[1], Q[0] + 2. * Q[2], sqr(Q[1])});
    // F0, F2 and F4 are computed from the operation -- F4 through F0; F1 and F3 are not
    std::vector<FFVar> F = {Y[0] + Q[2], Q[0] * Q[2], Y[1] * Y[2], Q[1], FFVar()};
    F[4] = sqr(F[0]);
    std::vector<double> const q = {0.3, -0.7, 1.1};
    std::vector<double> const ref = bad(DAG, F, Q, q, false);
    std::vector<double> const der = bad(DAG, F, Q, q, true);
    // first with one sweep, then with the three dependents computed from the operation
    if (Trace::asked != "1 3 " || Trace::n_adjoint != 3 || Trace::n_deriv != 0) return {false, "asked with " + Trace::asked};
    if (!near(der, ref)) return {false, "derivatives differ"};
    // declining for three sweeps: the partial derivatives instead
    Trace::max_sweeps = 2;
    std::vector<double> const der2 = bad(DAG, F, Q, q, true);
    Trace::max_sweeps = 1000;
    if (Trace::asked != "1 3 " || Trace::n_adjoint != 0 || Trace::n_deriv != 1) return {false, "declining: asked with " + Trace::asked};
    if (!near(der2, ref)) return {false, "derivatives differ when declining"};
    // a weighted sum of the dependents is one sweep, whatever their number
    Trace::max_sweeps = 1;
    Trace::reset();
    std::vector<FFVar> const dir = DAG.DBAD(F, W, Q);
    Trace::max_sweeps = 1000;
    if (Trace::asked != "1 " || Trace::n_adjoint != 1 || Trace::n_deriv != 0) return {false, "directional: asked with " + Trace::asked};
    // an operation that a single dependent is computed from is not asked again
    Trace::reset();
    std::vector<FFVar> const dF2 = DAG.BAD({F[2], F[3]}, Q);
    return {Trace::asked == "1 " && Trace::n_adjoint == 1, "single dependent: asked with " + Trace::asked};
  });

  run("repeated operands", []() -> Result {
    FFGraph DAG;
    std::vector<FFVar> Q = DAG.add_vars(2);
    std::vector<FFVar> Y = FFToy()({Q[0], Q[0], Q[1]});
    for (auto const& p : points)
    {
      std::vector<double> const q = {p[0], p[2]};
      std::vector<double> const ref = bad(DAG, Y, Q, q, false);
      std::vector<double> const der = bad(DAG, Y, Q, q, true);
      double J[9], x[3] = {q[0], q[0], q[1]};
      toy::jac(x, J);
      if (!near(der, ref) || !near(der, {J[0] + J[1], J[2], J[3] + J[4], J[5], J[6] + J[7], J[8]}))
        return {false, "derivatives differ"};
    }
    return {true, ""};
  });

  run("operation that no dependent is computed from, or no independent reaches, is not called", []() -> Result {
    FFGraph DAG;
    std::vector<FFVar> Q = DAG.add_vars(2), P = DAG.add_vars(3);
    std::vector<FFVar> Y = FFToy()(P);
    std::vector<FFVar> QP = Q;
    QP.insert(QP.end(), P.begin(), P.end());
    std::vector<double> const qp = {0.3, -0.7, 1.1, 0.4, -2.};
    // the operation is in the cone of the dependent, but its operands do not depend on Q
    std::vector<FFVar> F = {Y[0] * Q[0] + Q[1]};
    Trace::reset();
    std::vector<FFVar> const dF = DAG.BAD(F, Q);
    if (Trace::n_adjoint != 0 || Trace::n_deriv != 0) return {false, "operation called: " + counts()};
    if (!near(evaluate(DAG, dF, QP, qp), {1.1 * 0.4, 1.})) return {false, "derivatives"};
    // one of two dependents is computed from the operation
    std::vector<FFVar> G = {Q[0] * Q[1], Y[1] + Q[0]};
    Trace::reset();
    std::vector<FFVar> const dG = DAG.BAD(G, QP);
    if (Trace::n_adjoint != 1) return {false, "second case: " + counts()};
    return {near(evaluate(DAG, dG, QP, qp), {-0.7, 0.3, 0., 0., 0., 1., 0., std::cos(1.1), 0., 2. * -2.}), "derivatives, second case"};
  });

  run("second order, backward over backward: adjoint() of the adjoint operation", []() -> Result {
    FFGraph DAG;
    std::vector<FFVar> X = DAG.add_vars(3);
    std::vector<FFVar> Y = FFToy()(X);
    std::vector<FFVar> F = {Y[0] - 2. * Y[1] + 0.5 * Y[2]};
    Trace::reset();
    std::vector<FFVar> const dF  = DAG.BAD(F, X);
    std::vector<FFVar> const d2F = DAG.BAD(dF, X);
    // first pass: one FFToyAdj operation; second pass: it is in each of the three sweeps
    if (Trace::n_adjoint != 1 || Trace::n_adjadj != 3 || Trace::n_deriv != 0) return {false, counts()};
    for (auto const& x : points)
      if (!near(evaluate(DAG, d2F, X, x), hessian(x, {1., -2., 0.5}))) return {false, "Hessian"};
    return {true, ""};
  });

  run("second order, forward over backward: tangent() of the adjoint operation", []() -> Result {
    FFGraph DAG;
    std::vector<FFVar> X = DAG.add_vars(3), D = DAG.add_vars(3);
    std::vector<FFVar> Y = FFToy()(X);
    std::vector<FFVar> F = {Y[0] - 2. * Y[1] + 0.5 * Y[2]};
    Trace::reset();
    std::vector<FFVar> const dF  = DAG.BAD(F, X);
    std::vector<FFVar> const d2F = DAG.FAD(dF, X);
    if (Trace::n_adjoint != 1 || Trace::n_adjtan != 3 || Trace::n_deriv != 0) return {false, counts()};
    // Hessian-vector product: one adjoint and one tangent call, whatever the number of variables
    Trace::reset();
    std::vector<FFVar> const Hd = DAG.DFAD(dF, X, D);
    if (Trace::n_adjtan != 1) return {false, "Hessian-vector product: " + counts()};
    std::vector<FFVar> XD = X;
    XD.insert(XD.end(), D.begin(), D.end());
    for (auto const& x : points)
    {
      std::vector<double> const H = hessian(x, {1., -2., 0.5});
      if (!near(evaluate(DAG, d2F, X, x), H)) return {false, "Hessian"};
      std::vector<double> xd = x, d = {0.5, -1., 2.}, Hdref(3, 0.);
      xd.insert(xd.end(), d.begin(), d.end());
      for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) Hdref[i] += H[i * 3 + j] * d[j];
      if (!near(evaluate(DAG, Hd, XD, xd), Hdref)) return {false, "Hessian-vector product"};
    }
    return {true, ""};
  });

  run("second order, backward over forward: adjoint() of the tangent operation", []() -> Result {
    FFGraph DAG;
    std::vector<FFVar> X = DAG.add_vars(3);
    std::vector<FFVar> Y = FFToy()(X);
    std::vector<FFVar> F = {Y[0] - 2. * Y[1] + 0.5 * Y[2]};
    Trace::use_tangent = true;
    Trace::reset();
    std::vector<FFVar> const dF = DAG.FAD(F, X);
    Trace::use_tangent = false;
    std::vector<FFVar> const d2F = DAG.BAD(dF, X);
    // first pass: one FFToyTan operation per independent; second pass: one of them in each sweep
    if (Trace::n_tangent != 3 || Trace::n_tanadj != 3 || Trace::n_deriv != 0) return {false, counts()};
    for (auto const& x : points)
      if (!near(evaluate(DAG, d2F, X, x), hessian(x, {1., -2., 0.5}))) return {false, "Hessian"};
    return {true, ""};
  });

  run("second order through expressions of the independents", []() -> Result {
    FFGraph DAG;
    std::vector<FFVar> Q = DAG.add_vars(2);
    std::vector<FFVar> Y = FFToy()({Q[0] * Q[1], Q[0] + 2. * Q[1], sqr(Q[1])});
    std::vector<FFVar> F = {Y[0] * Y[2] + Y[1]};
    std::vector<FFVar> const dF   = DAG.BAD(F, Q);
    std::vector<FFVar> const d2Fb = DAG.BAD(dF, Q), d2Ff = DAG.FAD(dF, Q);
    // against central differences of the first derivatives
    double const h = 1e-5;
    for (auto const& p : points)
    {
      std::vector<double> const q = {p[0], p[1]}, H = evaluate(DAG, d2Fb, Q, q);
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
      if (!near(evaluate(DAG, d2Ff, Q, q), H, 1e-10)) return {false, "forward over backward differs from backward over backward"};
    }
    return {true, ""};
  });

  run("FAD: forward mode still calls deriv()", []() -> Result {
    FFGraph DAG;
    std::vector<FFVar> Q = DAG.add_vars(3);
    std::vector<FFVar> Y = FFToy()({Q[0] * Q[1], Q[0] + 2. * Q[2], sqr(Q[1])});
    std::vector<FFVar> F = {Y[0] + Q[2] * Y[2], Y[1] * Y[1]};
    Trace::reset();
    std::vector<FFVar> const dFf = DAG.FAD(F, Q);
    if (Trace::n_deriv != 1 || Trace::n_adjoint != 0) return {false, counts()};
    for (auto const& q : points)
      if (!near(evaluate(DAG, dFf, Q, q), bad(DAG, F, Q, q, true))) return {false, "derivatives differ"};
    return {true, ""};
  });

  run("differentiation of an operation without the hook nor deriv() is refused, and the DAG stays usable", []() -> Result {
    FFGraph DAG;
    std::vector<FFVar> X = DAG.add_vars(3);
    std::vector<FFVar> Y = FFToy()(X);
    std::vector<FFVar> F = {Y[0] - 2. * Y[1] + 0.5 * Y[2]};
    std::vector<FFVar> const dF  = DAG.BAD(F, X);
    std::vector<FFVar> const d2F = DAG.BAD(dF, X);  // FFToyAdj2 operations, which have no derivative
    int refused = 0;
    for (int mode = 0; mode < 2; ++mode) try
      {
        std::vector<FFVar> const d3F = (mode ? DAG.BAD(d2F, X) : DAG.FAD(d2F, X));
      }
      catch (mc::FFBase::Exceptions& e)
      {
        if (e.ierr() == mc::FFBase::Exceptions::EXTERN) ++refused;
      }
    if (refused != 2) return {false, "not refused"};
    // forward mode through FFToyTan, which has an adjoint only
    Trace::use_tangent = true;
    std::vector<FFVar> const dFf = DAG.FAD(F, X);
    Trace::use_tangent = false;
    try
    {
      std::vector<FFVar> const d2Ff = DAG.FAD(dFf, X);
    }
    catch (mc::FFBase::Exceptions& e)
    {
      if (e.ierr() == mc::FFBase::Exceptions::EXTERN) ++refused;
    }
    if (refused != 3) return {false, "forward mode through an operation with adjoint() only not refused"};
    return {near(evaluate(DAG, d2F, X, points[1]), hessian(points[1], {1., -2., 0.5})) &&
                near(evaluate(DAG, DAG.BAD(dF, X), X, points[1]), hessian(points[1], {1., -2., 0.5})),
            "after the refusals"};
  });

  run("exception in adjoint(): propagated, and the DAG stays usable", []() -> Result {
    FFGraph DAG;
    std::vector<FFVar> X = DAG.add_vars(3);
    std::vector<FFVar> Y = FFToy()(X);
    bool thrown = false;
    Trace::throw_adjoint = true;
    try
    {
      std::vector<FFVar> const dY = DAG.BAD(Y, X);
    }
    catch (std::runtime_error const& e)
    {
      thrown = (std::string(e.what()) == "FFToy::adjoint");
    }
    Trace::throw_adjoint = false;
    if (!thrown) return {false, "not propagated"};
    std::vector<double> y(3), J(9);
    toy::val(points[2].data(), y.data());
    toy::jac(points[2].data(), J.data());
    return {near(evaluate(DAG, Y, X, points[2]), y) && near(bad(DAG, Y, X, points[2], true), J), "after the exception"};
  });

  run("linear operator: backward mode through its adjoint(), or through forward mode without", []() -> Result {
    std::vector<double> const q = {0.3, -0.7, 1.1};
    for (int const hook : {0, 1})
    {
      FFScale::use_adjoint = hook;
      FFScale::n_adjoint = FFScale::n_apply = 0;
      FFGraph DAG;
      std::vector<FFVar> Q = DAG.add_vars(3);
      // an operator of an operator, operands that other operations consume too, and an operator
      // of an independent itself
      FFVar const Z = Q[0] * Q[1];
      std::vector<FFVar> F = {FFScale()(Z) + Q[1] + sqr(Z), sqr(FFScale()(Q[0])), FFScale()(exp(FFScale()(Z * Q[2])) + Z)};
      std::vector<FFVar> const dFf = DAG.FAD(F, Q);
      int const napply = FFScale::n_apply;
      std::vector<FFVar> const dFb = DAG.BAD(F, Q);
      FFScale::use_adjoint = true;
      // with the adjoint the operator is not applied to tangents; without, it is not transposed
      if (hook ? (FFScale::n_adjoint == 0 || FFScale::n_apply != napply) : (FFScale::n_adjoint != 0 || FFScale::n_apply == napply))
        return {false, std::string(hook ? "with" : "without") + " adjoint(): " + std::to_string(FFScale::n_adjoint) + " adjoint, " +
                           std::to_string(FFScale::n_apply - napply) + " apply calls in backward mode"};
      double const z = q[0] * q[1], e = std::exp(3. * z * q[2]);
      std::vector<double> const ref = {3. * q[1] + 2. * z * q[1], 3. * q[0] + 1. + 2. * z * q[0], 0.,
                                       18. * q[0], 0., 0.,
                                       3. * (3. * e * q[1] * q[2] + q[1]), 3. * (3. * e * q[0] * q[2] + q[0]), 9. * e * z};
      if (!near(evaluate(DAG, dFf, Q, q), ref)) return {false, "forward derivatives"};
      if (!near(evaluate(DAG, dFb, Q, q), ref)) return {false, std::string("backward derivatives ") + (hook ? "with" : "without") + " adjoint()"};
      // weighted sum of the dependents
      std::vector<FFVar> const W = {FFVar(0.5), FFVar(-2.), FFVar(0.25)};
      std::vector<double> const dw = evaluate(DAG, DAG.DBAD(F, W, Q), Q, q);
      for (int j = 0; j < 3; ++j)
        if (!near(dw[j], 0.5 * ref[j] - 2. * ref[3 + j] + 0.25 * ref[6 + j])) return {false, "weighted backward derivatives"};
    }
    return {true, ""};
  });

  run("linear operator without adjoint() that no independent reaches: not differentiated at all", []() -> Result {
    FFGraph DAG;
    std::vector<FFVar> Q = DAG.add_vars(2), P = DAG.add_vars(1);
    std::vector<FFVar> F = {FFScale()(sqr(P[0])) * Q[0] + Q[1]};
    FFScale::use_adjoint = false;
    FFScale::n_apply = 0;
    std::vector<FFVar> const dF = DAG.BAD(F, Q);
    FFScale::use_adjoint = true;
    if (FFScale::n_apply != 0) return {false, "operator applied"};
    std::vector<FFVar> QP = Q;
    QP.push_back(P[0]);
    return {near(evaluate(DAG, dF, QP, {0.3, -0.7, 1.1}), {3. * 1.1 * 1.1, 1.}), "derivatives"};
  });

  run("default: an operation without adjoint() is differentiated through deriv()", []() -> Result {
    // FFOp::has_adjoint() is false unless overridden, and FFOp::adjoint() throws
    FFToyJac Op;
    FFVar res, var, adj;
    FFVar const* padj = &res;
    bool thrown       = false;
    try
    {
      Op.FFOp::adjoint(1, &res, 1, &var, &padj, &adj);
    }
    catch (mc::FFBase::Exceptions& e)
    {
      thrown = (e.ierr() == mc::FFBase::Exceptions::EXTERN);
    }
    return {!Op.has_adjoint(1) && thrown, "defaults of the hook"};
  });

  std::printf("\n=== Results: %d/%d passed ===\n", test::passed, test::total);
  return test::failed ? 1 : 0;
}
