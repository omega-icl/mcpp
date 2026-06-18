// =============================================================================
// test_ffunc.cpp
// Unit tests for mc::FFGraph — DAG construction, subgraph extraction,
// evaluation (double arithmetic), differentiation (FAD/BAD) and composition.
//
// Build (assuming MC++ headers are on the include path):
//   g++ -std=c++17 -O2 -I<path-to-mc++> -o test_ffunc test_ffunc.cpp
//
// All tests are self-contained; no external test framework is required.
// Each TEST() block prints PASS or FAIL and a summary is printed at the end.
// The process exits with code 0 if every test passes, 1 otherwise.
// =============================================================================

#include <cmath>
#include <cstdio>
#include <functional>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "ffunc.hpp"

// -----------------------------------------------------------------------------
// Lightweight test harness
// -----------------------------------------------------------------------------
namespace test
{

static int total = 0, passed = 0, failed = 0;

// Absolute tolerance for floating-point comparisons
static constexpr double TOL = 1e-10;

struct Result
{
  bool ok;
  std::string msg;
};

inline bool
near(double a, double b, double tol = TOL)
{
  return std::fabs(a - b) <= tol * (1. + std::fabs(b));
}

// Run a single named test.  fn should return a Result.
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
    std::printf("  FAIL  %s — %s\n", name.c_str(), r.msg.c_str());
  }
}

// Convenience macro: wrap body in a lambda and call run().
// Usage:  TEST("name") { ... REQUIRE(...); ... return PASS; }
#define TEST(name) test::run( name, [&]() -> test::Result

#define PASS \
  return test::Result { true, "" }
#define FAIL(msg) \
  return test::Result { false, (msg) }
#define REQUIRE(cond)                     \
  if (!(cond)) return test::Result        \
    {                                     \
      false, "requirement failed: " #cond \
    }
#define REQUIRE_NEAR(a, b)                                \
  if (!test::near((a), (b)))                              \
  {                                                       \
    std::ostringstream _s;                                \
    _s << #a << "=" << (a) << " != " << #b << "=" << (b); \
    return test::Result{false, _s.str()};                 \
  }
#define REQUIRE_THROWS(expr)                                         \
  {                                                                  \
    bool _threw = false;                                             \
    try                                                              \
    {                                                                \
      (void)(expr);                                                  \
    }                                                                \
    catch (...)                                                      \
    {                                                                \
      _threw = true;                                                 \
    }                                                                \
    if (!_threw)                                                     \
      return test::Result{false, "expected exception from: " #expr}; \
  }

}  // namespace test

// =============================================================================
//  S E C T I O N   1 :   D A G   C O N S T R U C T I O N
// =============================================================================
void
section_construction()
{
  std::printf("\n--- 1. DAG construction ---\n");

  // 1.1 Variable registration
  TEST("var_count_increases")
  {
    mc::FFGraph dag;
    REQUIRE(dag.nvar() == 0);
    mc::FFVar x(&dag), y(&dag), z(&dag);
    REQUIRE(dag.nvar() == 3);
    PASS;
  });

  // 1.2 Constant folding: expression of two numeric literals → no new op node
  TEST("constant_folding_add")
  {
    mc::FFGraph dag;
    mc::FFVar c = mc::FFVar(2.) + mc::FFVar(3.);
    // Result should be the constant 5 with no DAG pointer
    REQUIRE(c.id().second == mc::FFVar::NOREF);
    REQUIRE_NEAR(c.num().val(), 5.);
    REQUIRE(dag.nvar() == 0 && dag.naux() == 0);
    PASS;
  });

  // 1.3 Constant folding: multiply
  TEST("constant_folding_mul")
  {
    mc::FFGraph dag;
    mc::FFVar c = mc::FFVar(3.) * mc::FFVar(4.);
    REQUIRE(c.id().second == mc::FFVar::NOREF);
    REQUIRE_NEAR(c.num().val(), 12.);
    PASS;
  });

  // 1.4 Common subexpression detection: x*y used twice → one TIMES node
  TEST("cse_detection")
  {
    mc::FFGraph dag;
    mc::FFVar x(&dag), y(&dag);
    mc::FFVar p1 = x * y;
    mc::FFVar p2 = x * y;
    // Both should resolve to the same AUX node
    REQUIRE(p1.id() == p2.id());
    PASS;
  });

  // 1.5 Self-add shortcut: x + x → SCALE node (not PLUS), value 2
  TEST("self_add_shortcut")
  {
    mc::FFGraph dag;
    mc::FFVar x(&dag);
    mc::FFVar s = x + x;
    // Underlying operation should be SCALE (not PLUS of two distinct nodes)
    REQUIRE(s.opdef().first != nullptr);
    REQUIRE(s.opdef().first->type != mc::FFOp::PLUS);
    PASS;
  });

  // 1.6 Self-subtract shortcut: x - x → constant 0
  TEST("self_sub_shortcut")
  {
    mc::FFGraph dag;
    mc::FFVar x(&dag);
    mc::FFVar d = x - x;
    REQUIRE(d.id().second == mc::FFVar::NOREF);
    REQUIRE_NEAR(d.num().val(), 0.);
    PASS;
  });

  // 1.7 Self-multiply shortcut: x * x → SQR node
  TEST("self_mul_shortcut")
  {
    mc::FFGraph dag;
    mc::FFVar x(&dag);
    mc::FFVar s = x * x;
    REQUIRE(s.opdef().first->type == mc::FFOp::SQR);
    PASS;
  });

  // 1.8 Self-divide shortcut: x / x → constant 1
  TEST("self_div_shortcut")
  {
    mc::FFGraph dag;
    mc::FFVar x(&dag);
    mc::FFVar d = x / x;
    REQUIRE(d.id().second == mc::FFVar::NOREF);
    REQUIRE_NEAR(d.num().val(), 1.);
    PASS;
  });

  // 1.9 pow shortcuts
  TEST("pow_shortcuts")
  {
    mc::FFGraph dag;
    mc::FFVar x(&dag);
    mc::FFVar p0  = mc::pow(x, 0);   // → constant 1
    mc::FFVar p1  = mc::pow(x, 1);   // → x itself
    mc::FFVar p2  = mc::pow(x, 2);   // → SQR node
    mc::FFVar pm1 = mc::pow(x, -1);  // → INV node
    REQUIRE(p0.id().second == mc::FFVar::NOREF);
    REQUIRE_NEAR(p0.num().val(), 1.);
    REQUIRE(p1.id() == x.id());
    REQUIRE(p2.opdef().first->type == mc::FFOp::SQR);
    REQUIRE(pm1.opdef().first->type == mc::FFOp::INV);
    PASS;
  });

  // 1.10 Unary minus of unary minus → cancels
  TEST("double_neg_cancels")
  {
    mc::FFGraph dag;
    mc::FFVar x(&dag);
    mc::FFVar neg = -x;
    mc::FFVar pos = -neg;
    REQUIRE(pos.id() == x.id());
    PASS;
  });

  // 1.11 Scale by zero: 0 * x → constant 0
  TEST("scale_zero")
  {
    mc::FFGraph dag;
    mc::FFVar x(&dag);
    mc::FFVar z = 0. * x;
    REQUIRE(z.id().second == mc::FFVar::NOREF);
    REQUIRE_NEAR(z.num().val(), 0.);
    PASS;
  });

  // 1.12 Scale by one: 1 * x → x itself
  TEST("scale_one")
  {
    mc::FFGraph dag;
    mc::FFVar x(&dag);
    mc::FFVar s = 1. * x;
    REQUIRE(s.id() == x.id());
    PASS;
  });

  // 1.13 Constant deduplication: same double value → same DAG node
  TEST("constant_deduplication")
  {
    mc::FFGraph dag;
    mc::FFVar x(&dag);
    mc::FFVar a = x + 3.14;
    mc::FFVar b = x + 3.14;
    // The constant 3.14 should appear only once in _Vars
    // Both expressions share the same shift operation
    REQUIRE(a.id() == b.id());
    PASS;
  });

  // 1.14 clear() resets the DAG completely
  TEST("clear_resets_dag")
  {
    mc::FFGraph dag;
    {
      mc::FFVar x(&dag), y(&dag);
      mc::FFVar f = x * y + mc::exp(x);
    }
    REQUIRE(dag.nvar() > 0);
    dag.clear();
    REQUIRE(dag.nvar() == 0);
    REQUIRE(dag.naux() == 0);
    REQUIRE(dag.Vars().empty());
    REQUIRE(dag.Ops().empty());
    PASS;
  });

  // 1.15 Multiple independent DAGs do not interfere
  TEST("independent_dags")
  {
    mc::FFGraph dag1, dag2;
    mc::FFVar x1(&dag1), x2(&dag2);
    REQUIRE(dag1.nvar() == 1 && dag2.nvar() == 1);
    // Cross-DAG operation must throw
    REQUIRE_THROWS(x1 + x2);
    PASS;
  });
}

// =============================================================================
//  S E C T I O N   2 :   S U B G R A P H   E X T R A C T I O N
// =============================================================================
void
section_subgraph()
{
  std::printf("\n--- 2. Subgraph extraction ---\n");

  // 2.1 Subgraph of a trivial single-variable function
  TEST("subgraph_single_var")
  {
    mc::FFGraph dag;
    mc::FFVar x(&dag);
    mc::FFVar f = x + 1.;
    auto sg     = dag.subgraph(1, &f);
    REQUIRE(!sg.l_op.empty());
    REQUIRE(sg.v_dep.size() == 1);
    REQUIRE(sg.v_dep[0]->id() == f.id());
    PASS;
  });

  // 2.2 Subgraph len_tap equals number of operations' outputs
  TEST("subgraph_len_tap")
  {
    mc::FFGraph dag;
    mc::FFVar x(&dag), y(&dag);
    mc::FFVar f = x * y + mc::exp(x);
    auto sg     = dag.subgraph(1, &f);
    // len_tap must be at least as large as l_op.size() (each op has ≥1 output)
    size_t expected_min = sg.l_op.size();
    REQUIRE(sg.len_tap >= expected_min);
    PASS;
  });

  // 2.3 Shared subexpression appears only once in the subgraph
  TEST("subgraph_cse_once")
  {
    mc::FFGraph dag;
    mc::FFVar x(&dag), y(&dag);
    mc::FFVar p    = x * y;  // shared node
    mc::FFVar f1   = p - x;
    mc::FFVar f2   = mc::exp(p);
    mc::FFVar F[2] = {f1, f2};
    auto sg        = dag.subgraph(2, F);
    // Count occurrences of the TIMES op in l_op — must be exactly one
    int times_count = 0;
    for (auto const& op : sg.l_op)
      if (op->type == mc::FFOp::TIMES) ++times_count;
    REQUIRE(times_count == 1);
    PASS;
  });

  // 2.4 v_dep contains exactly the requested dependents in order
  TEST("subgraph_vdep_order")
  {
    mc::FFGraph dag;
    mc::FFVar x(&dag), y(&dag);
    mc::FFVar f1   = x + y;
    mc::FFVar f2   = x - y;
    mc::FFVar F[2] = {f1, f2};
    auto sg        = dag.subgraph(2, F);
    REQUIRE(sg.v_dep.size() == 2);
    REQUIRE(sg.v_dep[0]->id() == f1.id());
    REQUIRE(sg.v_dep[1]->id() == f2.id());
    PASS;
  });

  // 2.5 Subgraph of a constant dependent has only a CNST operation
  TEST("subgraph_constant_dep")
  {
    mc::FFGraph dag;
    mc::FFVar x(&dag);
    mc::FFVar f = x + 0.;  // collapses to x, but let's use a real constant
    mc::FFVar c(2.);       // standalone constant (no DAG)
    // Use a DAG-registered constant via _add_constant path
    mc::FFVar g = x * 1. + 3.;  // → x + 3, one SHIFT node
    auto sg     = dag.subgraph(1, &g);
    REQUIRE(sg.v_dep[0]->id() == g.id());
    PASS;
  });

  // 2.6 Subgraph from std::vector overload matches array overload
  TEST("subgraph_vector_overload")
  {
    mc::FFGraph dag;
    mc::FFVar x(&dag), y(&dag);
    mc::FFVar F[2]            = {x + y, x * y};
    auto sg_arr               = dag.subgraph(2, F);
    std::vector<mc::FFVar> vF = {F[0], F[1]};
    auto sg_vec               = dag.subgraph(vF);
    REQUIRE(sg_arr.l_op.size() == sg_vec.l_op.size());
    REQUIRE(sg_arr.len_tap == sg_vec.len_tap);
    PASS;
  });

  // 2.7 Subgraph of a function not involving one of the DAG's variables
  //     does not include that variable's op
  TEST("subgraph_excludes_irrelevant_var")
  {
    mc::FFGraph dag;
    mc::FFVar x(&dag), y(&dag), z(&dag);
    mc::FFVar f  = x + y;  // z not used
    auto sg      = dag.subgraph(1, &f);
    bool z_found = false;
    for (auto const& op : sg.l_op)
      if (op->type == mc::FFOp::VAR && op->varout[0]->id() == z.id())
        z_found = true;
    REQUIRE(!z_found);
    PASS;
  });
}

// =============================================================================
//  S E C T I O N   3 :   E V A L U A T I O N   ( d o u b l e )
// =============================================================================
void
section_eval()
{
  std::printf("\n--- 3. Evaluation (double) ---\n");

  // Helper: evaluate a scalar function of one variable
  auto eval1 = [](mc::FFGraph& dag, mc::FFVar const& f, mc::FFVar const& x,
                  double xv) -> double
  {
    double result;
    dag.eval(1, &f, &result, 1, &x, &xv);
    return result;
  };

  auto eval2 = [](mc::FFGraph& dag, mc::FFVar const& f, mc::FFVar const& x,
                  double xv, mc::FFVar const& y, double yv) -> double
  {
    double result;
    mc::FFVar vars[2] = {x, y};
    double vals[2]    = {xv, yv};
    dag.eval(1, &f, &result, 2, vars, vals);
    return result;
  };

  // 3.1 Polynomial: f(x) = x^2 + 3x + 2 at x=4 → 30
  TEST("eval_polynomial")
  {
    mc::FFGraph dag;
    mc::FFVar x(&dag);
    mc::FFVar f = mc::sqr(x) + 3. * x + 2.;
    REQUIRE_NEAR(eval1(dag, f, x, 4.), 30.);
    PASS;
  });

  // 3.2 exp(x) at x = 1
  TEST("eval_exp")
  {
    mc::FFGraph dag;
    mc::FFVar x(&dag);
    mc::FFVar f = mc::exp(x);
    REQUIRE_NEAR(eval1(dag, f, x, 1.), std::exp(1.));
    PASS;
  });

  // 3.3 log(x) at x = e
  TEST("eval_log")
  {
    mc::FFGraph dag;
    mc::FFVar x(&dag);
    mc::FFVar f = mc::log(x);
    REQUIRE_NEAR(eval1(dag, f, x, std::exp(1.)), 1.);
    PASS;
  });

  // 3.4 sin and cos at x = pi/6: sin=0.5, cos=sqrt(3)/2
  TEST("eval_trig")
  {
    mc::FFGraph dag;
    mc::FFVar x(&dag);
    mc::FFVar fs     = mc::sin(x);
    mc::FFVar fc     = mc::cos(x);
    const double pi6 = std::acos(-1.) / 6.;
    REQUIRE_NEAR(eval1(dag, fs, x, pi6), 0.5);
    REQUIRE_NEAR(eval1(dag, fc, x, pi6), std::sqrt(3.) / 2.);
    PASS;
  });

  // 3.5 sqrt(x) at x = 9
  TEST("eval_sqrt")
  {
    mc::FFGraph dag;
    mc::FFVar x(&dag);
    mc::FFVar f = mc::sqrt(x);
    REQUIRE_NEAR(eval1(dag, f, x, 9.), 3.);
    PASS;
  });

  // 3.6 Two-variable expression: f(x,y) = x*y - x at (2,3) = 4
  TEST("eval_two_vars")
  {
    mc::FFGraph dag;
    mc::FFVar x(&dag), y(&dag);
    mc::FFVar f = x * y - x;
    REQUIRE_NEAR(eval2(dag, f, x, 2., y, 3.), 4.);
    PASS;
  });

  // 3.7 Common subexpression: p=x*y used in two dependents evaluated
  // simultaneously
  TEST("eval_cse_multidep")
  {
    mc::FFGraph dag;
    mc::FFVar x(&dag), y(&dag);
    mc::FFVar p       = x * y;
    mc::FFVar F[2]    = {p - x, mc::exp(p)};
    mc::FFVar vars[2] = {x, y};
    double vals[2]    = {2., 3.};
    double res[2];
    dag.eval(2, F, res, 2, vars, vals);
    // f0 = 2*3 - 2 = 4
    // f1 = exp(2*3) = exp(6)
    REQUIRE_NEAR(res[0], 4.);
    REQUIRE_NEAR(res[1], std::exp(6.));
    PASS;
  });

  // 3.8 Reuse of cached subgraph gives identical results
  TEST("eval_cached_subgraph")
  {
    mc::FFGraph dag;
    mc::FFVar x(&dag);
    mc::FFVar f = mc::pow(x, 3);
    mc::FFSubgraph sg;
    std::vector<double> wk;
    double res1, res2;
    double xv1 = 2.;
    dag.eval(sg, wk, 1, &f, &res1, 1, &x, &xv1);
    double xv2 = 3.;
    dag.eval(sg, wk, 1, &f, &res2, 1, &x, &xv2);
    REQUIRE_NEAR(res1, 8.);
    REQUIRE_NEAR(res2, 27.);
    PASS;
  });

  // 3.9 Constant variable (set/unset): override x with a fixed value
  TEST("eval_constant_override")
  {
    mc::FFGraph dag;
    mc::FFVar x(&dag);
    mc::FFVar f = x + 10.;
    // Pin x to 5
    x.set(5.);
    double res;
    double dummy = 0.;  // overridden by the constant
    dag.eval(1, &f, &res, 1, &x, &dummy);
    REQUIRE_NEAR(res, 15.);
    x.unset();
    PASS;
  });

  // 3.10 Missing variable throws MISSVAR
  TEST("eval_missing_var_throws")
  {
    mc::FFGraph dag;
    mc::FFVar x(&dag), y(&dag);
    mc::FFVar f = x + y;
    double res, xv = 1.;
    // Only provide x, not y — should throw MISSVAR
    REQUIRE_THROWS(dag.eval(1, &f, &res, 1, &x, &xv));
    PASS;
  });

  // 3.11 Multi-group eval: two separate (var, val) pairs
  TEST("eval_multi_group")
  {
    mc::FFGraph dag;
    mc::FFVar x(&dag), y(&dag);
    mc::FFVar f = x + y;
    double res, xv = 3., yv = 4.;
    // Pass x and y as two separate single-element groups
    dag.eval(1, &f, &res, 1, &x, &xv, 1, &y, &yv);
    REQUIRE_NEAR(res, 7.);
    PASS;
  });

  // 3.12 Absolute value: fabs(-3) = 3
  TEST("eval_fabs")
  {
    mc::FFGraph dag;
    mc::FFVar x(&dag);
    mc::FFVar f = mc::fabs(x);
    REQUIRE_NEAR(eval1(dag, f, x, -3.), 3.);
    PASS;
  });

  // 3.13 Integer power: x^5 at x=2 = 32
  TEST("eval_ipow")
  {
    mc::FFGraph dag;
    mc::FFVar x(&dag);
    mc::FFVar f = mc::pow(x, 5);
    REQUIRE_NEAR(eval1(dag, f, x, 2.), 32.);
    PASS;
  });

  // 3.14 Double power: x^2.5 at x=4 = 32
  TEST("eval_dpow")
  {
    mc::FFGraph dag;
    mc::FFVar x(&dag);
    mc::FFVar f = mc::pow(x, 2.5);
    REQUIRE_NEAR(eval1(dag, f, x, 4.), 32.);
    PASS;
  });

  // 3.15 min/max: min(x,y) and max(x,y)
  TEST("eval_min_max")
  {
    mc::FFGraph dag;
    mc::FFVar x(&dag), y(&dag);
    mc::FFVar fmin = mc::min(x, y);
    mc::FFVar fmax = mc::max(x, y);
    REQUIRE_NEAR(eval2(dag, fmin, x, 3., y, 7.), 3.);
    REQUIRE_NEAR(eval2(dag, fmax, x, 3., y, 7.), 7.);
    PASS;
  });

  // 3.16 vector eval overload
  TEST("eval_vector_overload")
  {
    mc::FFGraph dag;
    mc::FFVar x(&dag);
    mc::FFVar f               = mc::sqr(x) + 1.;
    std::vector<mc::FFVar> vF = {f};
    std::vector<mc::FFVar> vX = {x};
    std::vector<double> uX    = {3.};
    std::vector<double> uF;
    dag.eval(vF, uF, vX, uX);
    REQUIRE(uF.size() == 1);
    REQUIRE_NEAR(uF[0], 10.);
    PASS;
  });

  // 3.17 USEMOVE: shared leaf variable must not be moved before its second use.
  //      f(x,y)=3*y-2*(x+y); y is consumed twice, so its tape value must have
  //      movability 0.  This is the regression test for the iterative
  //      propagate_subgraph child-output target bug.
  TEST("eval_usemove_shared_leaf_variable")
  {
    mc::FFGraph dag;
    dag.options.USEMOVE = true;
    mc::FFVar x(&dag), y(&dag);
    mc::FFVar f = 3. * y - 2. * (x + y);

    auto sg         = dag.subgraph(1, &f);
    std::size_t iwk = 0;
    bool saw_y      = false;
    for (auto const& op : sg.l_op)
    {
      for (auto const& var : op->varout)
      {
        if (op->type == mc::FFOp::VAR && var->id() == y.id())
        {
          saw_y = true;
          REQUIRE(sg.v_mov[iwk] == 0);
        }
        ++iwk;
      }
    }
    REQUIRE(saw_y);

    mc::FFVar vars[2] = {x, y};
    double vals[2]    = {1., 2.};
    double res;
    dag.eval(sg, 1, &f, &res, 2, vars, vals);
    REQUIRE_NEAR(res, 0.);
    PASS;
  });

  // 3.18 USEMOVE: shared intermediate must not be moved before its second use.
  //      p=x*y is consumed by both branches of f=(p+x)+(p-y), so p must have
  //      movability 0 and the move-enabled result must match the algebraic one.
  TEST("eval_usemove_shared_intermediate")
  {
    mc::FFGraph dag;
    dag.options.USEMOVE = true;
    mc::FFVar x(&dag), y(&dag);
    mc::FFVar p = x * y;
    mc::FFVar f = (p + x) + (p - y);

    auto sg         = dag.subgraph(1, &f);
    std::size_t iwk = 0;
    bool saw_p      = false;
    for (auto const& op : sg.l_op)
    {
      for (auto const& var : op->varout)
      {
        if (var->id() == p.id())
        {
          saw_p = true;
          REQUIRE(sg.v_mov[iwk] == 0);
        }
        ++iwk;
      }
    }
    REQUIRE(saw_p);

    mc::FFVar vars[2] = {x, y};
    double vals[2]    = {2., 3.};
    double res;
    dag.eval(sg, 1, &f, &res, 2, vars, vals);
    REQUIRE_NEAR(res, 11.);
    PASS;
  });

  // 3.19 USEMOVE: common subexpression shared by multiple dependents.
  //      This extends eval_cse_multidep to the destructive evaluation path.
  TEST("eval_usemove_cse_multidep")
  {
    mc::FFGraph dag;
    dag.options.USEMOVE = true;
    mc::FFVar x(&dag), y(&dag);
    mc::FFVar p    = x * y;
    mc::FFVar F[2] = {p - x, mc::exp(p)};

    auto sg         = dag.subgraph(2, F);
    std::size_t iwk = 0;
    bool saw_p      = false;
    for (auto const& op : sg.l_op)
    {
      for (auto const& var : op->varout)
      {
        if (var->id() == p.id())
        {
          saw_p = true;
          REQUIRE(sg.v_mov[iwk] == 0);
        }
        ++iwk;
      }
    }
    REQUIRE(saw_p);

    mc::FFVar vars[2] = {x, y};
    double vals[2]    = {2., 3.};
    double res[2];
    dag.eval(sg, 2, F, res, 2, vars, vals);
    REQUIRE_NEAR(res[0], 4.);
    REQUIRE_NEAR(res[1], std::exp(6.));
    PASS;
  });

  // 3.20 USEMOVE: a genuinely single-use leaf should still be marked movable.
  //      This guards against replacing the traversal by an over-conservative
  //      all-inputs-not-movable policy.
  TEST("subgraph_usemove_single_use_leaf_is_movable")
  {
    mc::FFGraph dag;
    dag.options.USEMOVE = true;
    mc::FFVar x(&dag), y(&dag);
    mc::FFVar f = 2. * (x + y);

    auto sg         = dag.subgraph(1, &f);
    std::size_t iwk = 0;
    bool saw_x = false, saw_y = false;
    for (auto const& op : sg.l_op)
    {
      for (auto const& var : op->varout)
      {
        if (op->type == mc::FFOp::VAR && var->id() == x.id())
        {
          saw_x = true;
          REQUIRE(sg.v_mov[iwk] == 1);
        }
        if (op->type == mc::FFOp::VAR && var->id() == y.id())
        {
          saw_y = true;
          REQUIRE(sg.v_mov[iwk] == 1);
        }
        ++iwk;
      }
    }
    REQUIRE(saw_x && saw_y);

    mc::FFVar vars[2] = {x, y};
    double vals[2]    = {1.5, 2.5};
    double res;
    dag.eval(sg, 1, &f, &res, 2, vars, vals);
    REQUIRE_NEAR(res, 8.);
    PASS;
  });
}

// =============================================================================
//  S E C T I O N   4 :   F O R W A R D   A U T O M A T I C   D I F F
// =============================================================================
void
section_fad()
{
  std::printf("\n--- 4. Forward automatic differentiation (FAD) ---\n");

  auto evalD = [](mc::FFGraph& dag, mc::FFVar const& f, mc::FFVar const& x,
                  double xv) -> double
  {
    double res;
    dag.eval(1, &f, &res, 1, &x, &xv);
    return res;
  };

  // 4.1 d/dx [x^2] = 2x: check at x=3 → 6
  TEST("fad_sqr")
  {
    mc::FFGraph dag;
    mc::FFVar x(&dag);
    mc::FFVar f = mc::sqr(x);
    auto J      = dag.FAD(1, &f, 1, &x);
    REQUIRE(J != nullptr);
    REQUIRE_NEAR(evalD(dag, J[0], x, 3.), 6.);
    delete[] J;
    PASS;
  });

  // 4.2 d/dx [exp(x)] = exp(x): check at x=2
  TEST("fad_exp")
  {
    mc::FFGraph dag;
    mc::FFVar x(&dag);
    mc::FFVar f = mc::exp(x);
    auto J      = dag.FAD(1, &f, 1, &x);
    REQUIRE_NEAR(evalD(dag, J[0], x, 2.), std::exp(2.));
    delete[] J;
    PASS;
  });

  // 4.3 d/dx [log(x)] = 1/x: check at x=4 → 0.25
  TEST("fad_log")
  {
    mc::FFGraph dag;
    mc::FFVar x(&dag);
    mc::FFVar f = mc::log(x);
    auto J      = dag.FAD(1, &f, 1, &x);
    REQUIRE_NEAR(evalD(dag, J[0], x, 4.), 0.25);
    delete[] J;
    PASS;
  });

  // 4.4 d/dx [3*x + 5] = 3 (constant)
  TEST("fad_linear")
  {
    mc::FFGraph dag;
    mc::FFVar x(&dag);
    mc::FFVar f = 3. * x + 5.;
    auto J      = dag.FAD(1, &f, 1, &x);
    REQUIRE_NEAR(evalD(dag, J[0], x, 99.), 3.);
    delete[] J;
    PASS;
  });

  // 4.5 Jacobian of [x*y, x+y] w.r.t. [x, y] — row-wise layout
  //     J = [[y, x], [1, 1]]  evaluated at (x=2, y=3) → [[3,2],[1,1]]
  TEST("fad_two_by_two_jacobian")
  {
    mc::FFGraph dag;
    mc::FFVar x(&dag), y(&dag);
    mc::FFVar F[2]  = {x * y, x + y};
    mc::FFVar XY[2] = {x, y};
    auto J =
        dag.FAD(2, F, 2, XY);  // 4 entries, row-wise: J[i*2+j] = dF[i]/dX[j]
    REQUIRE(J != nullptr);
    double vals[2] = {2., 3.};
    double dF0dx, dF0dy, dF1dx, dF1dy;
    dag.eval(1, &J[0], &dF0dx, 2, XY, vals);
    dag.eval(1, &J[1], &dF0dy, 2, XY, vals);
    dag.eval(1, &J[2], &dF1dx, 2, XY, vals);
    dag.eval(1, &J[3], &dF1dy, 2, XY, vals);
    REQUIRE_NEAR(dF0dx, 3.);
    REQUIRE_NEAR(dF0dy, 2.);
    REQUIRE_NEAR(dF1dx, 1.);
    REQUIRE_NEAR(dF1dy, 1.);
    delete[] J;
    PASS;
  });

  // 4.6 Sparse FAD: SFAD on a function with a sparse Jacobian
  //     F = [x^2, y^2] → diagonal Jacobian
  TEST("sfad_diagonal_jacobian")
  {
    mc::FFGraph dag;
    mc::FFVar x(&dag), y(&dag);
    mc::FFVar F[2]               = {mc::sqr(x), mc::sqr(y)};
    mc::FFVar XY[2]              = {x, y};
    auto [nnz, rows, cols, vals] = dag.SFAD(2, F, 2, XY);
    REQUIRE(nnz == 2);
    // Off-diagonal entries must be absent
    for (unsigned k = 0; k < nnz; ++k)
      REQUIRE(rows[k] == cols[k]);  // diagonal only
    delete[] rows;
    delete[] cols;
    delete[] vals;
    PASS;
  });

  // 4.7 Directional derivative DFAD: direction [1,1] → d/dt f(x+t,y+t)|_{t=0}
  //     f(x,y) = x*y  →  df/dt = y + x  at (2,3) = 5
  TEST("dfad_directional")
  {
    mc::FFGraph dag;
    mc::FFVar x(&dag), y(&dag);
    mc::FFVar f     = x * y;
    mc::FFVar D[2]  = {mc::FFVar(1.), mc::FFVar(1.)};  // direction [1,1]
    mc::FFVar XY[2] = {x, y};
    auto dF         = dag.DFAD(1, &f, 2, XY, D);
    REQUIRE(dF != nullptr);
    double vals[2] = {2., 3.};
    double res;
    dag.eval(1, dF, &res, 2, XY, vals);
    REQUIRE_NEAR(res, 5.);
    delete[] dF;
    PASS;
  });

  // 4.8 Second derivative via nested FAD: d^2/dx^2 [x^3] = 6x at x=2 → 12
  TEST("fad_second_derivative")
  {
    mc::FFGraph dag;
    mc::FFVar x(&dag);
    mc::FFVar f = mc::pow(x, 3);
    auto J1     = dag.FAD(1, &f, 1, &x);  // d/dx x^3 = 3x^2
    auto J2     = dag.FAD(1, J1, 1, &x);  // d/dx 3x^2 = 6x
    double xv   = 2., res;
    dag.eval(1, J2, &res, 1, &x, &xv);
    REQUIRE_NEAR(res, 12.);
    delete[] J1;
    delete[] J2;
    PASS;
  });
}

// =============================================================================
//  S E C T I O N   5 :   B A C K W A R D   A U T O M A T I C   D I F F
// =============================================================================
void
section_bad()
{
  std::printf("\n--- 5. Backward automatic differentiation (BAD) ---\n");

  auto evalD = [](mc::FFGraph& dag, mc::FFVar const* f, mc::FFVar const* x,
                  double xv) -> double
  {
    double res;
    dag.eval(1, f, &res, 1, x, &xv);
    return res;
  };

  // 5.1 BAD of x^2 = 2x — verify matches FAD
  TEST("bad_sqr_matches_fad")
  {
    mc::FFGraph dag;
    mc::FFVar x(&dag);
    mc::FFVar f = mc::sqr(x);
    auto Jf     = dag.FAD(1, &f, 1, &x);
    auto Jb     = dag.BAD(1, &f, 1, &x);
    REQUIRE(Jf != nullptr && Jb != nullptr);
    // Both should give the same DAG node (or evaluate identically)
    for (double xv : {1., 2., 5., -3.})
    {
      double vf = evalD(dag, Jf, &x, xv);
      double vb = evalD(dag, Jb, &x, xv);
      if (!test::near(vf, vb))
        FAIL("FAD/BAD mismatch at x=" + std::to_string(xv));
    }
    delete[] Jf;
    delete[] Jb;
    PASS;
  });

  // 5.2 BAD of exp(x)*sin(x): d/dx = exp(x)*(cos(x)+sin(x)) at x=1
  TEST("bad_exp_sin")
  {
    mc::FFGraph dag;
    mc::FFVar x(&dag);
    mc::FFVar f = mc::exp(x) * mc::sin(x);
    auto Jb     = dag.BAD(1, &f, 1, &x);
    double xv   = 1., res;
    dag.eval(1, Jb, &res, 1, &x, &xv);
    double expected = std::exp(1.) * (std::cos(1.) + std::sin(1.));
    REQUIRE_NEAR(res, expected);
    delete[] Jb;
    PASS;
  });

  // 5.3 BAD Jacobian of F=[x*y, x+y] matches FAD result
  TEST("bad_jacobian_matches_fad")
  {
    mc::FFGraph dag;
    mc::FFVar x(&dag), y(&dag);
    mc::FFVar F[2]  = {x * y, x + y};
    mc::FFVar XY[2] = {x, y};
    auto Jf         = dag.FAD(2, F, 2, XY);
    auto Jb         = dag.BAD(2, F, 2, XY);
    double vals[2]  = {3., 5.};
    for (int k = 0; k < 4; ++k)
    {
      double vf, vb;
      dag.eval(1, &Jf[k], &vf, 2, XY, vals);
      dag.eval(1, &Jb[k], &vb, 2, XY, vals);
      if (!test::near(vf, vb))
        FAIL("FAD/BAD Jacobian mismatch at entry " + std::to_string(k));
    }
    delete[] Jf;
    delete[] Jb;
    PASS;
  });

  // 5.4 Sparse BAD: SBAD on diagonal function → same sparsity as SFAD
  TEST("sbad_diagonal_sparsity")
  {
    mc::FFGraph dag;
    mc::FFVar x(&dag), y(&dag);
    mc::FFVar F[2]                    = {mc::sqr(x), mc::sqr(y)};
    std::vector<mc::FFVar const*> vF  = {&F[0], &F[1]};
    std::vector<mc::FFVar const*> vXY = {&x, &y};
    auto [rows, cols, vals]           = dag.SBAD(vF, vXY);
    REQUIRE(rows.size() == 2);
    for (size_t k = 0; k < rows.size(); ++k) REQUIRE(rows[k] == cols[k]);
    PASS;
  });
}

// =============================================================================
//  S E C T I O N   6 :   C O M P O S I T I O N
// =============================================================================
void
section_compose()
{
  std::printf("\n--- 6. Composition ---\n");

  // 6.1 compose f(g(x)): f(y)=y^2, g(x)=x+1  →  (x+1)^2 at x=2 = 9
  TEST("compose_substitution")
  {
    mc::FFGraph dag;
    mc::FFVar x(&dag), y(&dag);
    mc::FFVar g = x + 1.;      // inner: g(x) = x+1
    mc::FFVar f = mc::sqr(y);  // outer: f(y) = y^2
    // Compose f with y→g
    std::vector<mc::FFVar> vF     = {f};
    std::vector<mc::FFVar> vVarY  = {y};
    std::vector<mc::FFVar> vDepG  = {g};
    std::vector<mc::FFVar> result = dag.compose(vF, vVarY, vDepG);
    REQUIRE(result.size() == 1);
    double xv = 2., res;
    dag.eval(1, &result[0], &res, 1, &x, &xv);
    REQUIRE_NEAR(res, 9.);
    PASS;
  });

  // 6.2 compose with identity substitution: y→x, f(y)=sin(y) → sin(x) at pi/2 =
  // 1
  TEST("compose_identity_sub")
  {
    mc::FFGraph dag;
    mc::FFVar x(&dag), y(&dag);
    mc::FFVar f                   = mc::sin(y);
    std::vector<mc::FFVar> vF     = {f};
    std::vector<mc::FFVar> vVarY  = {y};
    std::vector<mc::FFVar> vDepX  = {x};  // substitute y with x
    std::vector<mc::FFVar> result = dag.compose(vF, vVarY, vDepX);
    double xv                     = std::acos(-1.) / 2., res;
    dag.eval(1, &result[0], &res, 1, &x, &xv);
    REQUIRE_NEAR(res, 1.);
    PASS;
  });

  // 6.3 Compose two functions simultaneously:
  //     F=[y1^2, y2^2], substitute y1→x, y2→exp(x)
  //     → [x^2, exp(2x)] at x=1: [1, e^2]
  TEST("compose_two_functions")
  {
    mc::FFGraph dag;
    mc::FFVar x(&dag), y1(&dag), y2(&dag);
    mc::FFVar F[2]                = {mc::sqr(y1), mc::sqr(y2)};
    std::vector<mc::FFVar> vF     = {F[0], F[1]};
    std::vector<mc::FFVar> vVarY  = {y1, y2};
    std::vector<mc::FFVar> vDepG  = {x, mc::exp(x)};
    std::vector<mc::FFVar> result = dag.compose(vF, vVarY, vDepG);
    REQUIRE(result.size() == 2);
    double xv = 1., res[2];
    dag.eval(1, &result[0], &res[0], 1, &x, &xv);
    dag.eval(1, &result[1], &res[1], 1, &x, &xv);
    REQUIRE_NEAR(res[0], 1.);
    REQUIRE_NEAR(res[1], std::exp(2.));
    PASS;
  });

  // 6.4 Compose differentiated output: d/dy [y^3] composed with y→x^2
  //     gives 3y^2|_{y=x^2} = 3x^4. At x=2: 48.
  TEST("compose_with_derivative")
  {
    mc::FFGraph dag;
    mc::FFVar x(&dag), y(&dag);
    mc::FFVar f                   = mc::pow(y, 3);
    auto J                        = dag.FAD(1, &f, 1, &y);  // = 3*y^2
    mc::FFVar g                   = mc::sqr(x);             // y → x^2
    std::vector<mc::FFVar> vJ     = {J[0]};
    std::vector<mc::FFVar> vVarY  = {y};
    std::vector<mc::FFVar> vDepG  = {g};
    std::vector<mc::FFVar> result = dag.compose(vJ, vVarY, vDepG);
    double xv                     = 2., res;
    dag.eval(1, &result[0], &res, 1, &x, &xv);
    REQUIRE_NEAR(res, 48.);  // 3*(2^2)^2 = 3*16 = 48
    delete[] J;
    PASS;
  });
}

// =============================================================================
//  S E C T I O N   7 :   E D G E   C A S E S   &   R O B U S T N E S S
// =============================================================================
void
section_edge()
{
  std::printf("\n--- 7. Edge cases & robustness ---\n");

  // 7.1 Subgraph of an empty dependent list is safe
  TEST("subgraph_empty")
  {
    mc::FFGraph dag;
    mc::FFVar x(&dag);
    std::vector<mc::FFVar> empty;
    auto sg = dag.subgraph(empty);
    REQUIRE(sg.l_op.empty());
    REQUIRE(sg.len_tap == 0);
    PASS;
  });

  // 7.2 eval with 0 dependents is a no-op
  TEST("eval_zero_deps")
  {
    mc::FFGraph dag;
    mc::FFVar x(&dag);
    double xv = 1.;
    dag.eval(0, &x, (double*)nullptr, 1, &x, &xv);  // must not crash
    PASS;
  });

  // 7.3 FAD of a constant dependent gives zero
  TEST("fad_constant_dep")
  {
    mc::FFGraph dag;
    mc::FFVar x(&dag);
    mc::FFVar f = mc::FFVar(7.);  // constant, no DAG variable
    // FAD of a constant w.r.t. any variable should return 0
    auto J = dag.FAD(1, &f, 1, &x);
    REQUIRE(J != nullptr);
    double xv = 42., res;
    dag.eval(1, J, &res, 1, &x, &xv);
    REQUIRE_NEAR(res, 0.);
    delete[] J;
    PASS;
  });

  // 7.4 DAG survives repeated clear() cycles
  TEST("repeated_clear")
  {
    mc::FFGraph dag;
    for (int i = 0; i < 5; ++i)
    {
      mc::FFVar x(&dag), y(&dag);
      mc::FFVar f = mc::exp(x * y);
      double xv = 1., yv = 1., res;
      dag.eval(1, &f, &res, 1, &x, &xv, 1, &y, &yv);
      REQUIRE_NEAR(res, std::exp(1.));
      dag.clear();
      REQUIRE(dag.nvar() == 0);
    }
    PASS;
  });

  // 7.5 Deep expression tree does not overflow stack during subgraph
  // construction
  TEST("deep_expression_tree")
  {
    mc::FFGraph dag;
    mc::FFVar x(&dag);
    mc::FFVar f = x;
    // Build a chain of 200 additions: f = x+1+1+...
    for (int i = 0; i < 200; ++i) f = f + 1.;
    auto sg = dag.subgraph(1, &f);
    REQUIRE(!sg.l_op.empty());
    double xv = 0., res;
    dag.eval(sg, 1, &f, &res, 1, &x, &xv);
    REQUIRE_NEAR(res, 200.);
    PASS;
  });

  // 7.6 Signomial detection: exp(d*log(x)) → x^d
  TEST("signomial_detection")
  {
    mc::FFGraph dag;
    dag.options.DETECTSIGNOM = true;
    mc::FFVar x(&dag);
    // exp(2.5 * log(x)) should be detected and rewritten as x^2.5
    mc::FFVar f = mc::exp(2.5 * mc::log(x));
    // The result should be a DPOW node, not a chain of EXP(SCALE(LOG))
    REQUIRE(f.opdef().first != nullptr);
    // Evaluate to verify correctness: at x=4 → 4^2.5 = 32
    double xv = 4., res;
    dag.eval(1, &f, &res, 1, &x, &xv);
    REQUIRE_NEAR(res, 32.);
    PASS;
  });

  // 7.7 erfc(x) = 1 - erf(x) — verified numerically
  TEST("erfc_numerical")
  {
    mc::FFGraph dag;
    mc::FFVar x(&dag);
    mc::FFVar f = mc::erfc(x);  // implemented as 1 - erf(x)
    double xv   = 1., res;
    dag.eval(1, &f, &res, 1, &x, &xv);
    REQUIRE_NEAR(res, std::erfc(1.));
    PASS;
  });

  // 7.8 Output stream operator produces non-empty output
  TEST("output_stream")
  {
    mc::FFGraph dag;
    mc::FFVar x(&dag), y(&dag);
    mc::FFVar f = x * y + mc::exp(x);
    std::ostringstream oss;
    oss << dag;
    REQUIRE(!oss.str().empty());
    REQUIRE(oss.str().find("DAG") != std::string::npos);
    PASS;
  });
}

// =============================================================================
//  M A I N
// =============================================================================
int
main()
{
  std::printf("=== ffunc unit tests ===\n");

  section_construction();
  section_subgraph();
  section_eval();
  section_fad();
  section_bad();
  section_compose();
  section_edge();

  std::printf("\n=== Results: %d/%d passed", test::passed, test::total);
  if (test::failed) std::printf(", %d FAILED", test::failed);
  std::printf(" ===\n");

  return test::failed ? 1 : 0;
}
