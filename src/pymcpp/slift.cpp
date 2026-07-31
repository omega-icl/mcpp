// Copyright (C) Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the Eclipse Public License.

#include "slift.hpp"

#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

void
mc_slift(py::module_& m)
{
  // Typedefs
  typedef mc::SLiftEnv SL;
  typedef mc::SLiftVar SV;
  typedef mc::FFGraph FG;
  typedef mc::FFVar FV;
  typedef mc::SMon<FV const*, mc::lt_FFVar> SM;
  typedef mc::SPoly<FV const*, mc::lt_FFVar> SP;

  py::class_<SL> pySLiftEnv(m, "SLift", R"doc(
Environment for recursive decomposition (lifting) of factorable expressions.

SLift decomposes factorable expressions held in an :class:`FFGraph` DAG into
a collection of sparse polynomial/rational subexpressions and transcendental
subexpressions, by introducing auxiliary variables. This reformulation is a
standard preprocessing step in global optimization: after lifting, the
original dependents are polynomial/rational in the original and auxiliary
variables, and every non-polynomial term is isolated in a simple lifting
constraint (e.g. ``x4 - exp(x2) = 0``).

Typical workflow:

1. Build a DAG of the expressions to be decomposed.
2. Create an ``SLift`` environment attached to that DAG.
3. Call :meth:`process` on the dependent variables.
4. Retrieve the results from the properties :attr:`dep` (reformulated
   dependents), :attr:`poly` (polynomial/rational lifting constraints),
   :attr:`trans` (transcendental lifting constraints), :attr:`var`
   (participating variables) and :attr:`aux` (mapping between original DAG
   intermediates and new auxiliary variables).

With ``add2dag=True`` (the default in :meth:`process`), the auxiliary
variables and the lifted constraint expressions are inserted as new nodes
into the attached DAG.

The behaviour of the decomposition (which terms get lifted, which
univariate terms get substituted) is controlled by the fields of
:attr:`options`.

The class is also exported under its C++ name ``SLiftEnv`` as an alias.

Examples
--------
>>> import pymcpp
>>> DAG = pymcpp.FFGraph()
>>> X = DAG.add_vars(2, "X")
>>> F = X[0]**3 + pymcpp.sqrt(X[0]**2 + X[1]**2)
>>> env = pymcpp.SLift(DAG)
>>> env.process([F])
>>> print(env)  # displays lifted variables and constraints
<BLANKLINE>
2 AUXILIARY VARIABLES:
 V2    <- Z4     = SQR( X0 ) + SQR( X1 )
 V3    <- Z5     = SQRT( SQR( X0 ) + SQR( X1 ) )
<BLANKLINE>
1 AUXILIARY POLYNOMIAL CONSTRAINT:
  0 = V2 - ( SQR( X0 ) + SQR( X1 ) )
<BLANKLINE>
1 AUXILIARY TRANSCENDENTAL CONSTRAINT:
  0 = SQR( V3 ) - V2
<BLANKLINE>
1 DEPENDENT EXPRESSION:
  V3 + IPOW( X0, 3 )
<BLANKLINE>
)doc");
  py::class_<SV> pySLiftVar(m, "SLiftVar", R"doc(
Sparse rational expression used as the arithmetic of DAG lifting.

An ``SLiftVar`` represents a ratio of two sparse multivariate polynomials
(numerator and denominator) in the DAG variables. Propagating ``SLiftVar``
objects through the operations of a DAG is how :class:`SLift` performs the
recursive decomposition: polynomial/rational operations are absorbed into
the numerator/denominator polynomials, whereas transcendental operations
trigger the introduction of an auxiliary variable and a lifting constraint
in the associated :class:`SLift` environment.

Objects of this class are mostly created internally by
:meth:`SLift.process`; direct construction is only needed in advanced use
cases such as implementing custom external operations that must support the
lifting arithmetic (see :meth:`SLift.lift`).

Standard arithmetic operators (``+``, ``-``, ``*``, ``/``, in-place
variants, and unary ``+``) are overloaded and combine the underlying sparse
rational expressions; mixed operations with ``float`` scalars are
supported. The module-level functions ``exp``, ``log``, ``sin``, ...
have ``SLiftVar`` overloads that record the corresponding transcendental
subexpression for lifting.
)doc");

  // Preserve the C++ class name as a Python alias, while keeping the existing
  // shorter historical name.
  m.attr("SLiftEnv") = pySLiftEnv;

  // Translate mc::SLiftEnv::Exceptions into a regular Python exception instead
  // of pybind11's generic "unknown exception" fallback.
  py::register_exception_translator(
      [](std::exception_ptr p)
      {
        try
        {
          if (p) std::rethrow_exception(p);
        }
        catch (SL::Exceptions& e)
        {
          PyErr_SetString(PyExc_RuntimeError, e.what().c_str());
        }
      });

  pySLiftVar
      .def(py::init<double const>(), py::arg("d") = 0., R"doc(
Construct a constant sparse rational expression with value `d`.

Parameters
----------
d : float, optional
    Constant value. Default is 0.
)doc")
      .def(py::init<SL*, FV const&>(), py::arg("env"), py::arg("x"),
           py::keep_alive<1, 2>(),  // SLiftVar stores a raw SLiftEnv*
           py::keep_alive<1, 3>(),  // sparse polynomial keys store FFVar const*
           R"doc(
Construct a sparse rational expression representing DAG variable `x`.

Parameters
----------
env : SLift
    Lifting environment the new expression is linked to.
x : FFVar
    DAG variable; must belong to the DAG attached to `env`.
)doc")
      .def(py::init<SV const&>(), R"doc(
Copy constructor.
)doc")
      .def("env", &SV::env, py::return_value_policy::reference_internal, R"doc(
Return the SLift environment this expression is linked to.

Returns
-------
env : SLift
    Associated lifting environment, or None for a constant expression
    created without an environment.
)doc")
      .def("numer", &SV::numer, py::return_value_policy::reference_internal,
           R"doc(
Return the numerator as a sparse polynomial in the DAG variables.

Returns
-------
numer : FFPoly
    Numerator sparse polynomial.
)doc")
      .def("denom", &SV::denom, py::return_value_policy::reference_internal,
           R"doc(
Return the denominator as a sparse polynomial in the DAG variables.

Returns
-------
denom : FFPoly
    Denominator sparse polynomial (equal to 1 for polynomial expressions).
)doc")
      .def("set", &SV::set, py::arg("env"), py::arg("x"),
           py::keep_alive<1, 2>(),  // SLiftVar stores a raw SLiftEnv*
           py::keep_alive<1, 3>(),  // sparse polynomial keys store FFVar const*
           py::return_value_policy::reference_internal,
           R"doc(
Reinitialize this expression as DAG variable `x` in environment `env`.

Parameters
----------
env : SLift
    Lifting environment the expression is linked to.
x : FFVar
    DAG variable; must belong to the DAG attached to `env`.

Returns
-------
var : SLiftVar
    This expression, after reinitialization.
)doc")
      .def("__str__",
           [](SV const& v)
           {
             std::ostringstream oss;
             oss << v;
             return oss.str();
           })
      .def("__repr__",
           [](SV const& v)
           {
             std::ostringstream oss;
             oss << v;
             return oss.str();
           })
      .def(+py::self)
      .def(py::self += py::self)
      .def(py::self += double())
      .def(py::self -= py::self)
      .def(py::self -= double())
      .def(py::self *= py::self)
      .def(py::self *= double())
      .def(py::self /= py::self)
      .def(py::self /= double())
      .def(py::self + py::self)
      .def(py::self + double())
      .def(double() + py::self)
      .def(py::self - py::self)
      .def(py::self - double())
      .def(double() - py::self)
      .def(py::self * py::self)
      .def(py::self * double())
      .def(double() * py::self)
      .def(py::self / py::self)
      .def(py::self / double())
      .def(double() / py::self);

  m.def(
      "inv", [](SV const& x) { return mc::inv(x); }, py::arg("x"), R"doc(
SLiftVar overload: reciprocal 1/x in the lifting arithmetic (rational).
)doc");
  m.def(
      "sqr", [](SV const& x) { return mc::sqr(x); }, py::arg("x"), R"doc(
SLiftVar overload: square x**2 in the lifting arithmetic (polynomial).
)doc");
  m.def(
      "sqrt", [](SV const& x) { return mc::sqrt(x); }, py::arg("x"), R"doc(
SLiftVar overload: sqrt(x) in the lifting arithmetic; lifted as a new
auxiliary (converted to a sqr constraint if option SQRT2SQR is enabled).
)doc");
  m.def(
      "exp", [](SV const& x) { return mc::exp(x); }, py::arg("x"), R"doc(
SLiftVar overload: exp(x) in the lifting arithmetic; introduces a
transcendental lifting constraint.
)doc");
  m.def(
      "log", [](SV const& x) { return mc::log(x); }, py::arg("x"), R"doc(
SLiftVar overload: log(x) in the lifting arithmetic; introduces a
transcendental lifting constraint (as exp if option LOG2EXP is enabled).
)doc");
  m.def(
      "xlog", [](SV const& x) { return mc::xlog(x); }, py::arg("x"), R"doc(
SLiftVar overload: x*log(x) in the lifting arithmetic; introduces a
transcendental lifting constraint.
)doc");
  m.def(
      "cos", [](SV const& x) { return mc::cos(x); }, py::arg("x"), R"doc(
SLiftVar overload: cos(x) in the lifting arithmetic; introduces a
transcendental lifting constraint.
)doc");
  m.def(
      "sin", [](SV const& x) { return mc::sin(x); }, py::arg("x"), R"doc(
SLiftVar overload: sin(x) in the lifting arithmetic; introduces a
transcendental lifting constraint.
)doc");
  m.def(
      "tan", [](SV const& x) { return mc::tan(x); }, py::arg("x"), R"doc(
SLiftVar overload: tan(x) in the lifting arithmetic; introduces a
transcendental lifting constraint (as atan if option TAN2ATAN is enabled).
)doc");
  m.def(
      "acos", [](SV const& x) { return mc::acos(x); }, py::arg("x"), R"doc(
SLiftVar overload: acos(x) in the lifting arithmetic; introduces a
transcendental lifting constraint (as cos if option ACOS2COS is enabled).
)doc");
  m.def(
      "asin", [](SV const& x) { return mc::asin(x); }, py::arg("x"), R"doc(
SLiftVar overload: asin(x) in the lifting arithmetic; introduces a
transcendental lifting constraint (as sin if option ASIN2SIN is enabled).
)doc");
  m.def(
      "atan", [](SV const& x) { return mc::atan(x); }, py::arg("x"), R"doc(
SLiftVar overload: atan(x) in the lifting arithmetic; introduces a
transcendental lifting constraint (as tan if option ATAN2TAN is enabled).
)doc");
  m.def(
      "cosh", [](SV const& x) { return mc::cosh(x); }, py::arg("x"), R"doc(
SLiftVar overload: cosh(x) in the lifting arithmetic; introduces a
transcendental lifting constraint.
)doc");
  m.def(
      "sinh", [](SV const& x) { return mc::sinh(x); }, py::arg("x"), R"doc(
SLiftVar overload: sinh(x) in the lifting arithmetic; introduces a
transcendental lifting constraint.
)doc");
  m.def(
      "tanh", [](SV const& x) { return mc::tanh(x); }, py::arg("x"), R"doc(
SLiftVar overload: tanh(x) in the lifting arithmetic; introduces a
transcendental lifting constraint.
)doc");
  m.def(
      "fabs", [](SV const& x) { return mc::fabs(x); }, py::arg("x"), R"doc(
SLiftVar overload: absolute value |x| in the lifting arithmetic; introduces
a lifting constraint.
)doc");
  m.def(
      "erf", [](SV const& x) { return mc::erf(x); }, py::arg("x"), R"doc(
SLiftVar overload: error function erf(x) in the lifting arithmetic;
introduces a transcendental lifting constraint.
)doc");
  m.def(
      "fstep", [](SV const& x) { return mc::fstep(x); }, py::arg("x"), R"doc(
SLiftVar overload: forward step function (1 if x>=0, else 0) in the lifting
arithmetic; introduces a lifting constraint.
)doc");
  m.def(
      "pow", [](SV const& x, int const n) { return mc::pow(x, n); },
      py::arg("x"), py::arg("n"), R"doc(
SLiftVar overload: integer power x**n in the lifting arithmetic
(polynomial/rational; lifted if option LIFTIPOW is enabled).
)doc");
  m.def(
      "pow", [](SV const& x, double const& r) { return mc::pow(x, r); },
      py::arg("x"), py::arg("r"), R"doc(
SLiftVar overload: real power x**r in the lifting arithmetic; introduces a
transcendental lifting constraint.
)doc");
  m.def(
      "cheb", [](SV const& x, unsigned const n) { return mc::cheb(x, n); },
      py::arg("x"), py::arg("n"), R"doc(
SLiftVar overload: Chebyshev polynomial of the first kind T_n(x) in the
lifting arithmetic (polynomial).
)doc");
  m.def(
      "prod", [](std::vector<SV> const& x)
      { return mc::prod(x.size(), x.data()); }, py::arg("x"), R"doc(
SLiftVar overload: product of all elements of `x` in the lifting arithmetic
(polynomial/rational).
)doc");
  m.def(
      "prod",
      [](unsigned n, std::vector<SV> const& x)
      {
        if (n > x.size()) throw py::value_error("prod: n exceeds len(x)");
        return mc::prod(n, x.data());
      },
      py::arg("n"), py::arg("x"), R"doc(
SLiftVar overload: product of the first `n` elements of `x` in the lifting
arithmetic. Raises ValueError if n exceeds len(x).
)doc");
  m.def(
      "monom",
      [](std::vector<SV> const& x, std::vector<unsigned> const& k,
         bool const chebbasis)
      {
        if (k.size() != x.size())
          throw py::value_error("monom: len(k) must match len(x)");
        return mc::monom(x.size(), x.data(), k.data(), chebbasis);
      },
      py::arg("x"), py::arg("k"), py::arg("chebbasis") = false, R"doc(
SLiftVar overload: monomial prod_i x[i]**k[i] in the lifting arithmetic; if
`chebbasis` is True, use Chebyshev basis terms T_k[i](x[i]) instead of
powers. Raises ValueError if len(k) differs from len(x).
)doc");
  m.def(
      "monom",
      [](unsigned n, std::vector<SV> const& x, std::vector<unsigned> const& k,
         bool const chebbasis)
      {
        if (n > x.size() || n > k.size())
          throw py::value_error("monom: n exceeds len(x) or len(k)");
        return mc::monom(n, x.data(), k.data(), chebbasis);
      },
      py::arg("n"), py::arg("x"), py::arg("k"), py::arg("chebbasis") = false,
      R"doc(
SLiftVar overload: monomial over the first `n` variables,
prod_{i<n} x[i]**k[i], in the lifting arithmetic; if `chebbasis` is True,
use Chebyshev basis terms instead of powers.
)doc");
  m.def(
      "max", [](SV const& x, SV const& y) { return mc::max(x, y); },
      py::arg("x"), py::arg("y"), R"doc(
SLiftVar overload: max(x, y) in the lifting arithmetic; introduces a
lifting constraint.
)doc");
  m.def(
      "min", [](SV const& x, SV const& y) { return mc::min(x, y); },
      py::arg("x"), py::arg("y"), R"doc(
SLiftVar overload: min(x, y) in the lifting arithmetic; introduces a
lifting constraint.
)doc");
  m.def(
      "lmtd", [](SV const& x, SV const& y) { return mc::lmtd(x, y); },
      py::arg("x"), py::arg("y"), R"doc(
SLiftVar overload: log-mean temperature difference (x-y)/log(x/y) in the
lifting arithmetic; introduces a transcendental lifting constraint.
)doc");
  m.def(
      "rlmtd", [](SV const& x, SV const& y) { return mc::rlmtd(x, y); },
      py::arg("x"), py::arg("y"), R"doc(
SLiftVar overload: reciprocal log-mean temperature difference
log(x/y)/(x-y) in the lifting arithmetic; introduces a transcendental
lifting constraint.
)doc");

  // --- Exceptions Nested Class ---
  py::class_<SL::Exceptions> pySLExceptions(pySLiftEnv, "Exceptions", R"doc(
Exception information raised by SLift operations.

Instances carry an error code (:meth:`ierr`) and a human-readable
description (:meth:`what`). At the Python level, SLift errors surface as
``RuntimeError`` with the corresponding description.
)doc");

  py::enum_<SL::Exceptions::TYPE>(pySLExceptions, "TYPE")
      .value("DAGERR", SL::Exceptions::DAGERR,
             "Operation involving a factorable expression linked to a "
             "different DAG")
      .value("ENVERR", SL::Exceptions::ENVERR,
             "Operation between factorable expressions linked to different "
             "environments")
      .value("EXTERNAL", SL::Exceptions::EXTERNAL, "Invalid external operation")
      .value("INTERNAL", SL::Exceptions::INTERNAL, "Internal error")
      .export_values();

  pySLExceptions
      .def("ierr", &SL::Exceptions::ierr, R"doc(
Return the error code of this exception (see SLift.Exceptions.TYPE).
)doc")
      .def("what", &SL::Exceptions::what, R"doc(
Return a human-readable description of this exception.
)doc");

  // --- Options Nested Struct ---
  py::class_<SL::Options> pySLOptions(pySLiftEnv, "Options", R"doc(
Options controlling the decomposition performed by SLift.

Set the fields on the :attr:`SLift.options` attribute of an existing
environment, e.g. ``env.options.LIFTDIV = True``, before calling
:meth:`SLift.process`.
)doc");

  pySLOptions
      .def(py::init<>(), R"doc(
Construct an option set initialized with the default values.
)doc")
      .def(py::init<SL::Options const&>(), R"doc(
Copy constructor.
)doc")
      .def("reset", &SL::Options::reset, R"doc(
Reset all options to their default values.
)doc")
      .def_readwrite("KEEPFACT", &SL::Options::KEEPFACT,
                     "bool: Whether to keep existing factorizations, e.g. "
                     "products between polynomial subexpressions; otherwise "
                     "the polynomials are expanded. Default is True.")
      .def_readwrite(
          "NOAUXIL", &SL::Options::NOAUXIL,
          "bool: Whether to reformulate expressions without introducing "
          "auxiliary variables or lifting constraints. Default is False.")
      .def_readwrite("LIFTDIV", &SL::Options::LIFTDIV,
                     "bool: Whether to lift division terms using auxiliary "
                     "variables. Default is False.")
      .def_readwrite("LIFTIPOW", &SL::Options::LIFTIPOW,
                     "bool: Whether to lift integral power terms using "
                     "auxiliary variables. Default is False.")
      .def_readwrite(
          "LOG2EXP", &SL::Options::LOG2EXP,
          "bool: Whether to convert log univariates to exp univariates in "
          "lifting constraints. Default is True.")
      .def_readwrite("SQRT2SQR", &SL::Options::SQRT2SQR,
                     "bool: Whether to convert sqrt univariates to sqr "
                     "univariates in lifting constraints. Default is True.")
      .def_readwrite("ACOS2COS", &SL::Options::ACOS2COS,
                     "bool: Whether to convert acos univariates to cos "
                     "univariates in lifting constraints. Default is True.")
      .def_readwrite("ASIN2SIN", &SL::Options::ASIN2SIN,
                     "bool: Whether to convert asin univariates to sin "
                     "univariates in lifting constraints. Default is True.")
      .def_readwrite("ATAN2TAN", &SL::Options::ATAN2TAN,
                     "bool: Whether to convert atan univariates to tan "
                     "univariates in lifting constraints. Default is False.")
      .def_readwrite("TAN2ATAN", &SL::Options::TAN2ATAN,
                     "bool: Whether to convert tan univariates to atan "
                     "univariates in lifting constraints. Default is True.")
      .def_readwrite("DISPFULL", &SL::Options::DISPFULL,
                     "bool: Whether to display a full subgraph of the lifted "
                     "expressions when printing the environment. Default is "
                     "False.");

  // --- SLiftEnv Main Class ---
  pySLiftEnv
      .def(py::init<FG*>(), py::arg("dag") = nullptr,
           py::keep_alive<1, 2>(),  // SLiftEnv stores a raw FFGraph*
           R"doc(
Construct a lifting environment, optionally attached to a DAG.

Parameters
----------
dag : FFGraph, optional
    DAG holding the expressions to be decomposed. May be omitted and set
    later with :meth:`set`. Default is None.
)doc")
      .def_readwrite("options", &SL::options, R"doc(
SLift.Options: Options controlling the decomposition (which terms are
lifted and which univariates are substituted).
)doc")
      .def_property_readonly("dag", &SL::dag,
                             py::return_value_policy::reference_internal,
                             R"doc(
FFGraph: The DAG attached to this environment (None if unset).
)doc")
      .def("set", &SL::set, py::arg("dag"),
           py::keep_alive<1, 2>(),  // SLiftEnv stores a raw FFGraph*
           R"doc(
Attach a DAG to this environment and reset all intermediate results.

Parameters
----------
dag : FFGraph
    DAG holding the expressions to be decomposed.
)doc")
      .def("reset", &SL::reset, R"doc(
Discard all intermediate expressions and lifting results.

The attached DAG itself is kept; nodes previously added to it by
:meth:`process` are not removed.
)doc")
      .def("process",
           py::overload_cast<std::vector<FV> const&, bool const>(&SL::process),
           py::arg("vDep"), py::arg("add2dag") = true, R"doc(
Decompose the dependent expressions in `vDep`.

Recursively decomposes each dependent into sparse polynomial/rational and
transcendental subexpressions, introducing auxiliary variables for the
non-polynomial parts. The results are made available through the
properties :attr:`dep`, :attr:`poly`, :attr:`trans`, :attr:`var` and
:attr:`aux`.

Parameters
----------
vDep : list of FFVar
    Dependent expressions to decompose; must belong to the attached DAG.
add2dag : bool, optional
    Whether to insert the auxiliary variables and the lifted constraint
    expressions as new nodes into the attached DAG. Default is True.

Raises
------
RuntimeError
    If a dependent is linked to a different DAG, or on internal errors.

Examples
--------
>>> import pymcpp
>>> DAG = pymcpp.FFGraph()
>>> X = DAG.add_vars(2, "X")
>>> F = X[0]**3 + pymcpp.sqrt(X[0]**2 + X[1]**2)
>>> env = pymcpp.SLift(DAG)
>>> env.process([F])
>>> env.poly, env.trans  # lifting constraints added to the DAG
([Z7], [Z9])
)doc")
      .def("process",
           py::overload_cast<std::set<unsigned> const&, std::vector<FV> const&,
                             bool const>(&SL::process),
           py::arg("ndxDep"), py::arg("vDep"), py::arg("add2dag") = true,
           R"doc(
Decompose a subset of the dependents in `vDep` selected by index.

Same as the first overload but only processes the entries of `vDep` whose
indices appear in `ndxDep`.

Parameters
----------
ndxDep : set of int
    Indices into `vDep` selecting the dependents to decompose.
vDep : list of FFVar
    Dependent expressions; must belong to the attached DAG.
add2dag : bool, optional
    Whether to insert the lifted variables and constraints into the
    attached DAG. Default is True.
)doc")
      // Accessors
      .def_property_readonly("dep", &SL::Dep,
                             py::return_value_policy::reference_internal,
                             R"doc(
list of FFVar: Reformulated expressions of the processed dependents, i.e.
each original dependent rewritten in terms of the original and auxiliary
variables. Populated by :meth:`process`.
)doc")
      .def_property_readonly("poly", &SL::Poly,
                             py::return_value_policy::reference_internal,
                             R"doc(
list of FFVar: Polynomial/rational lifting constraints, one DAG expression
per constraint, to be interpreted as equal to zero. Populated by
:meth:`process`.
)doc")
      .def_property_readonly("trans", &SL::Trans,
                             py::return_value_policy::reference_internal,
                             R"doc(
list of FFVar: Transcendental lifting constraints (e.g.
``exp(x4) - x5``), one DAG expression per constraint, to be interpreted as
equal to zero. Populated by :meth:`process`.
)doc")
      .def_property_readonly("var", &SL::Var,
                             py::return_value_policy::reference_internal,
                             R"doc(
list of FFVar: All independent DAG variables participating in the lifted
expressions, comprising the original variables followed by the new
auxiliary variables. Populated by :meth:`process`.
)doc")
      .def_property_readonly(
          "aux",
          //&SL::Aux,
          [](SL& self)
          {
            auto const& aux_set = self.Aux();
            return std::vector<std::pair<mc::FFVar const*, mc::FFVar const*>>(
                aux_set.begin(), aux_set.end());
          },
          py::return_value_policy::reference_internal,
          R"doc(
list of tuple of (FFVar, FFVar): Mapping between lifted DAG subexpressions
and the auxiliary variables substituting them. Each tuple holds:

- the original DAG intermediate that was lifted, and
- the new auxiliary DAG variable introduced for it.

Populated by :meth:`process`.
)doc")
      .def_property_readonly("OpLift", &SL::OpLift,
                             py::return_value_policy::reference_internal,
                             R"doc(
list of tuple of (FFOp, list of SLiftVar): Intermediate decomposition
results per lifted DAG operation. Each tuple holds the DAG operation and
the sparse rational expressions of its operands. Mainly useful for
debugging the decomposition.
)doc")
      .def_property_readonly("AuxLift", &SL::AuxLift,
                             py::return_value_policy::reference_internal,
                             R"doc(
list of tuple of (FFVar, SLiftVar): Intermediate decomposition results per
lifted variable. Each tuple holds a DAG variable and its sparse rational
expression. Mainly useful for debugging the decomposition.
)doc")
      .def(
          "insert_dag",
          py::overload_cast<SM const&, bool const, bool const>(&SL::insert_dag),
          py::arg("mon"), py::arg("useprod") = false, py::arg("dagaux") = false,
          R"doc(
Transcribe a sparse monomial into the attached DAG.

Creates DAG nodes evaluating the monomial and returns the resulting
dependent. New nodes are added to the attached DAG.

Parameters
----------
mon : FFMon
    Sparse monomial in the DAG variables.
useprod : bool, optional
    Whether to encode the monomial with a single n-ary product node
    instead of chained binary multiplications. Default is False.
dagaux : bool, optional
    Whether to reference the existing DAG variables directly instead of
    substituting the lifted auxiliary variables. Default is False.

Returns
-------
var : FFVar
    DAG expression of the monomial.
)doc")
      .def(
          "insert_dag",
          py::overload_cast<SP const&, bool const, bool const>(&SL::insert_dag),
          py::arg("poly"), py::arg("useprod") = false,
          py::arg("dagaux") = false, R"doc(
Transcribe a sparse polynomial into the attached DAG.

Same as the monomial overload but transcribes a complete sparse
polynomial, i.e. the weighted sum of its monomials.

Parameters
----------
poly : FFPoly
    Sparse polynomial in the DAG variables.
useprod : bool, optional
    Whether to encode monomials with n-ary product nodes. Default is
    False.
dagaux : bool, optional
    Whether to reference the existing DAG variables directly instead of
    substituting the lifted auxiliary variables. Default is False.

Returns
-------
var : FFVar
    DAG expression of the polynomial.
)doc")
      .def(
          "lift",
          [](SL& E, unsigned const nres, std::vector<SV> const& vVar)
          {
            std::vector<SV> vRes(nres);
            E.lift(nres, vRes.data(), vVar.size(), vVar.data());
            return vRes;
          },
          py::arg("nres"), py::arg("vVar"),
          R"doc(
Lift an external (black-box) operation with the given operands.

Introduces one new auxiliary variable per result of the external operation
and records the operation for later insertion in the DAG. This is the
building block used by external operations (e.g. FFCustom, FFVect) to
support the lifting arithmetic.

Parameters
----------
nres : int
    Number of results of the external operation.
vVar : list of SLiftVar
    Sparse rational expressions of the operands.

Returns
-------
vRes : list of SLiftVar
    New list of `nres` lifted result expressions, each representing one
    auxiliary variable.
)doc")
      .def("__str__",
           [](SL const& E)
           {
             std::ostringstream Ess;
             Ess << E;
             return Ess.str();
           })
      .def("__repr__",
           [](SL const& E)
           {
             std::ostringstream Ess;
             Ess << E;
             return Ess.str();
           });
}
